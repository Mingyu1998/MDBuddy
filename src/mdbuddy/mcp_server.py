"""MCP server: MDBuddy's tools for any MCP-capable agent (`mdbuddy mcp`).

Model Context Protocol over stdio (newline-delimited JSON-RPC 2.0), implemented with the
standard library only. Works with Claude Code/Desktop, OpenAI Codex, Gemini CLI, Cursor,
VS Code/Copilot and other MCP clients; see docs/agents.md for their configuration.
"""

from __future__ import annotations

import contextlib
import io
import json
import sys
import traceback
from pathlib import Path
from typing import Any, Callable

from . import __version__

PROTOCOL_VERSIONS = ("2025-06-18", "2025-03-26", "2024-11-05")

INSTRUCTIONS = """MDBuddy builds Martini coarse-grained membrane-protein systems for GROMACS.
Workflow: doctor -> write_spec / validate_spec -> build -> report. Read the skill
(skills/mdbuddy/SKILL.md in the MDBuddy repository) for the full procedure and decisions.
Never hand-edit build outputs; change the spec and rebuild. Lipid names differ between Martini
versions: use lipids_show / lipids_convert before writing lipids into a spec."""


def _schema(props: dict, required: list[str] = ()) -> dict:
    return {"type": "object", "properties": props, "required": list(required), "additionalProperties": False}


def _text(obj) -> str:
    return obj if isinstance(obj, str) else json.dumps(obj, indent=2, default=str)


# --- tool implementations -------------------------------------------------------------------

def t_doctor(network: bool = False):
    from . import doctor
    checks = doctor.run_checks(network=network)
    return json.loads(doctor.report(checks, as_json=True))


def t_validate_spec(spec_path: str):
    from . import spec as spec_mod
    try:
        s = spec_mod.load(spec_path)
    except spec_mod.SpecError as exc:
        return {"valid": False, "problems": exc.problems}
    return {"valid": True, "output_dir": str(s.output_dir()), "resolved": s.data}


def t_write_spec(spec_path: str, yaml_text: str, overwrite: bool = False):
    import yaml
    from . import spec as spec_mod
    path = Path(spec_path)
    if path.exists() and not overwrite:
        return {"written": False, "error": f"{path} exists; pass overwrite=true to replace it"}
    try:
        spec_mod.from_dict(yaml.safe_load(yaml_text) or {}, source_path=path.resolve())
    except spec_mod.SpecError as exc:
        return {"written": False, "problems": exc.problems}
    except yaml.YAMLError as exc:
        return {"written": False, "error": f"not valid YAML: {exc}"}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml_text)
    return {"written": True, "path": str(path.resolve())}


def t_build(spec_path: str, from_stage: str = "fetch", until_stage: str = "check", force: bool = False,
            grompp_check: bool = True):
    from . import spec as spec_mod
    from .pipeline import build
    log: list[str] = []
    try:
        s = spec_mod.load(spec_path)
        result = build(s, from_stage=from_stage, until_stage=until_stage, grompp_check=grompp_check,
                       force=force, log=log.append)
    except spec_mod.SpecError as exc:
        return {"ok": False, "problems": exc.problems}
    except RuntimeError as exc:
        return {"ok": False, "error": str(exc), "log": log}
    out = Path(result["output_dir"])
    return {"ok": True, "log": log, "output_dir": str(out), "gromacs_dir": str(out / "gromacs"),
            "lock": str(result["lock"]), "readme": (out / "gromacs" / "README.md").read_text()
            if (out / "gromacs" / "README.md").exists() else None}


def t_lipids_show(names: list[str]):
    from . import lipids
    out = {}
    for n in names:
        out[n] = {lib: lipids.lookup(n, lib) for lib in lipids.LIBRARIES}
        out[n]["equivalents"] = {f"{a}->{b}": lipids.equivalents(n, a, b)
                                 for a in lipids.LIBRARIES for b in lipids.LIBRARIES
                                 if a != b and lipids.lookup(n, a)}
    return out


def t_lipids_convert(names: list[str], source: str = "m2", target: str = "m3"):
    from . import lipids
    from .cli import LIB_ALIASES
    src, dst = LIB_ALIASES[source.lower()], LIB_ALIASES[target.lower()]
    out = {}
    for n in names:
        e = lipids.lookup(n, src)
        out[n] = {"source": e and lipids.describe(n, src),
                  "candidates": [{"name": c, "description": lipids.describe(c, dst),
                                  "charge": lipids.lookup(c, dst).get("charge"),
                                  "source_charge": e.get("charge") if e else None,
                                  "coby": lipids.lookup(c, dst)["coby"]["status"]}
                                 for c in lipids.equivalents(n, src, dst)]}
    return out


def t_lipids_search(text: str, ff: str = "", limit: int = 40):
    from . import lipids
    from .cli import LIB_ALIASES
    lib = LIB_ALIASES[ff.lower()] if ff else None
    return [{"library": l, "name": n, "description": e.get("common") or e.get("chemistry")}
            for l, n, e in lipids.search(text, lib)[:limit]]


def t_orient_guess(pdb_path: str, html_path: str = "", output_pdb: str = ""):
    from . import orient, viewer
    from .cli import _protein_atoms
    from .structure import write_pdb
    atoms = _protein_atoms(pdb_path)
    placed, record = orient.hydrophobic(atoms)
    normal, centre = orient.membrane_in_input_frame(record)
    if output_pdb:
        write_pdb(output_pdb, placed)
    if html_path:
        viewer.write_html(Path(html_path), atoms, record, thickness_A=record["placement"]["thickness_A"],
                          title=Path(pdb_path).name)
    return {"placement": record["placement"], "spec": {"protein": {"orientation": "manual", "manual": {
        "normal": [round(v, 4) for v in normal], "center_point": [round(v, 3) for v in centre]}}},
        "html": html_path or None, "oriented_pdb": output_pdb or None}


def t_orient_view(pdb_path: str, html_path: str, start: str = "guess"):
    from . import orient, viewer
    from .cli import _protein_atoms
    atoms = _protein_atoms(pdb_path)
    record, thick = ({"method": "input frame"}, 30.0)
    if start == "guess":
        _, record = orient.hydrophobic(atoms)
        thick = record["placement"]["thickness_A"]
    viewer.write_html(Path(html_path), atoms, record, thickness_A=thick, title=Path(pdb_path).name)
    return {"html": str(Path(html_path).resolve()),
            "next": "ask the user to open it, adjust, and paste the spec lines back"}


def t_selfassembly_check(build_dir: str, gro: str = ""):
    from . import selfassembly
    res = selfassembly.check_build(Path(build_dir), Path(gro) if gro else None)
    (Path(build_dir) / "self_assembly_result.json").write_text(json.dumps(res, indent=2) + "\n")
    return res


_STAGES = ["fetch", "clean", "orient", "forcefield", "martinize", "membrane", "gromacs", "check"]
TOOLS: dict[str, tuple[Callable, str, dict, bool]] = {
    # name: (function, description, input schema, read_only)
    "doctor": (t_doctor, "Check that martinize2, COBY, DSSP/mdtraj, PDBFixer and GROMACS are installed. "
               "Run first; report missing items and their install commands to the user.",
               _schema({"network": {"type": "boolean", "description": "also test OPM/RCSB/GitHub access"}}), True),
    "write_spec": (t_write_spec, "Validate a system spec (YAML text) and write it to a file. "
                   "Returns every problem if invalid.",
                   _schema({"spec_path": {"type": "string"}, "yaml_text": {"type": "string"},
                            "overwrite": {"type": "boolean"}}, ["spec_path", "yaml_text"]), False),
    "validate_spec": (t_validate_spec, "Validate a spec file and return it with all defaults filled in.",
                      _schema({"spec_path": {"type": "string"}}, ["spec_path"]), True),
    "build": (t_build, "Build the system described by a spec (fetch -> clean -> orient -> forcefield -> "
              "martinize -> membrane -> gromacs -> check). Returns the stage log; read and explain it.",
              _schema({"spec_path": {"type": "string"},
                       "from_stage": {"type": "string", "enum": _STAGES},
                       "until_stage": {"type": "string", "enum": _STAGES},
                       "force": {"type": "boolean", "description": "overwrite an existing output directory"},
                       "grompp_check": {"type": "boolean"}}, ["spec_path"]), False),
    "lipids_show": (t_lipids_show, "What lipid names mean in Martini 2, Martini 3 (2021) and the Martini 3 "
                    "lipidome, and their equivalents. Names can mean different molecules across versions.",
                    _schema({"names": {"type": "array", "items": {"type": "string"}}}, ["names"]), True),
    "lipids_convert": (t_lipids_convert, "Convert lipid names between force fields (m2, m3 = lipidome, m3-2021); "
                       "best chemical match first, with charges.",
                       _schema({"names": {"type": "array", "items": {"type": "string"}},
                                "source": {"type": "string"}, "target": {"type": "string"}}, ["names"]), True),
    "lipids_search": (t_lipids_search, "Search lipids by name or chemistry (e.g. 'cardiolipin', 'C18:2').",
                      _schema({"text": {"type": "string"}, "ff": {"type": "string"},
                               "limit": {"type": "integer"}}, ["text"]), True),
    "orient_guess": (t_orient_guess, "Educated guess of a protein's membrane placement from exposed "
                     "hydrophobicity (for proteins without OPM data). Returns the placement, a verdict, "
                     "and spec lines; optionally writes the oriented PDB and the interactive viewer.",
                     _schema({"pdb_path": {"type": "string"}, "html_path": {"type": "string"},
                              "output_pdb": {"type": "string"}}, ["pdb_path"]), False),
    "orient_view": (t_orient_view, "Write the interactive HTML viewer to place a protein in the membrane by "
                    "hand (CHARMM-GUI style); it prints spec lines the user pastes back.",
                    _schema({"pdb_path": {"type": "string"}, "html_path": {"type": "string"},
                             "start": {"type": "string", "enum": ["guess", "input"]}},
                            ["pdb_path", "html_path"]), False),
    "selfassembly_check": (t_selfassembly_check, "After a self-assembly run: did a bilayer form, does the "
                           "protein span it, and the orientation spec for a clean COBY rebuild.",
                           _schema({"build_dir": {"type": "string"}, "gro": {"type": "string"}},
                                   ["build_dir"]), False),
}


# --- protocol ------------------------------------------------------------------------------

def handle(msg: dict) -> dict | None:
    """Process one JSON-RPC message; return the response (None for notifications)."""
    method, mid = msg.get("method"), msg.get("id")
    if mid is None:  # notification (e.g. notifications/initialized)
        return None
    try:
        if method == "initialize":
            requested = (msg.get("params") or {}).get("protocolVersion")
            version = requested if requested in PROTOCOL_VERSIONS else PROTOCOL_VERSIONS[0]
            result: Any = {"protocolVersion": version, "capabilities": {"tools": {"listChanged": False}},
                           "serverInfo": {"name": "mdbuddy", "version": __version__},
                           "instructions": INSTRUCTIONS}
        elif method == "ping":
            result = {}
        elif method == "tools/list":
            result = {"tools": [{"name": n, "description": d, "inputSchema": s,
                                 "annotations": {"readOnlyHint": ro}} for n, (_, d, s, ro) in TOOLS.items()]}
        elif method == "tools/call":
            params = msg.get("params") or {}
            name, args = params.get("name"), params.get("arguments") or {}
            if name not in TOOLS:
                return {"jsonrpc": "2.0", "id": mid, "error": {"code": -32602, "message": f"unknown tool {name}"}}
            func = TOOLS[name][0]
            captured = io.StringIO()
            try:
                with contextlib.redirect_stdout(captured):  # stdout belongs to the protocol
                    value = func(**args)
                result = {"content": [{"type": "text", "text": _text(value)}], "isError": False}
            except TypeError as exc:
                result = {"content": [{"type": "text", "text": f"bad arguments for {name}: {exc}"}], "isError": True}
            except Exception as exc:  # noqa: BLE001 - tool errors are reported to the model
                result = {"content": [{"type": "text", "text": f"{type(exc).__name__}: {exc}"}], "isError": True}
            if captured.getvalue():
                sys.stderr.write(captured.getvalue())
        else:
            return {"jsonrpc": "2.0", "id": mid, "error": {"code": -32601, "message": f"method not found: {method}"}}
    except Exception as exc:  # noqa: BLE001
        traceback.print_exc(file=sys.stderr)
        return {"jsonrpc": "2.0", "id": mid, "error": {"code": -32603, "message": str(exc)}}
    return {"jsonrpc": "2.0", "id": mid, "result": result}


def serve(stdin=None, stdout=None) -> None:
    stdin = stdin or sys.stdin
    stdout = stdout or sys.stdout
    for line in stdin:
        line = line.strip()
        if not line:
            continue
        try:
            msg = json.loads(line)
        except json.JSONDecodeError:
            reply = {"jsonrpc": "2.0", "id": None, "error": {"code": -32700, "message": "parse error"}}
        else:
            reply = handle(msg) if isinstance(msg, dict) else {
                "jsonrpc": "2.0", "id": None, "error": {"code": -32600, "message": "batches are not supported"}}
        if reply is not None:
            stdout.write(json.dumps(reply) + "\n")
            stdout.flush()


CLIENT_CONFIGS = {
    "claude": "claude mcp add mdbuddy -- {cmd} mcp",
    "codex": "# ~/.codex/config.toml (or .codex/config.toml in a project)\n[mcp_servers.mdbuddy]\n"
             "command = \"{cmd}\"\nargs = [\"mcp\"]\n\n# or: codex mcp add mdbuddy -- {cmd} mcp",
    "gemini": "// ~/.gemini/settings.json (or .gemini/settings.json in a project)\n"
              "{{\"mcpServers\": {{\"mdbuddy\": {{\"command\": \"{cmd}\", \"args\": [\"mcp\"]}}}}}}",
    "cursor": "// ~/.cursor/mcp.json (or .cursor/mcp.json in a project)\n"
              "{{\"mcpServers\": {{\"mdbuddy\": {{\"command\": \"{cmd}\", \"args\": [\"mcp\"]}}}}}}",
    "vscode": "// .vscode/mcp.json\n{{\"servers\": {{\"mdbuddy\": {{\"type\": \"stdio\", \"command\": \"{cmd}\", "
              "\"args\": [\"mcp\"]}}}}}}",
    "claude-desktop": "// claude_desktop_config.json\n"
                      "{{\"mcpServers\": {{\"mdbuddy\": {{\"command\": \"{cmd}\", \"args\": [\"mcp\"]}}}}}}",
}


def client_config(agent: str) -> str:
    from .tools import find_program
    cmd = find_program("mdbuddy") or "mdbuddy"
    return CLIENT_CONFIGS[agent].format(cmd=cmd)
