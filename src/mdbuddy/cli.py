"""Command line: mdbuddy doctor | init | validate | build | stages."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import yaml

from . import __version__, doctor, lipids, spec as spec_mod
from .pipeline import STAGES, BuildError, build


def cmd_doctor(args) -> int:
    checks = doctor.run_checks(network=args.network)
    print(doctor.report(checks, as_json=args.json))
    return 0 if doctor.ready(checks) else 1


def cmd_init(args) -> int:
    path = Path(args.output or f"{args.name}.yaml")
    if path.exists() and not args.force:
        print(f"{path} exists; use --force to overwrite", file=sys.stderr)
        return 1
    path.write_text(spec_mod.TEMPLATE.format(name=args.name, pdb_id=args.pdb_id or "1bl8"))
    print(f"wrote {path}")
    return 0


def cmd_validate(args) -> int:
    try:
        s = spec_mod.load(args.spec)
    except spec_mod.SpecError as exc:
        if args.json:
            print(json.dumps({"valid": False, "problems": exc.problems}, indent=2))
        else:
            print(exc, file=sys.stderr)
        return 1
    if args.json:
        print(json.dumps({"valid": True, "resolved": s.data}, indent=2))
    else:
        print(s.to_yaml())
        print(f"# spec is valid; output directory: {s.output_dir()}")
    return 0


def cmd_build(args) -> int:
    try:
        s = spec_mod.load(args.spec)
    except spec_mod.SpecError as exc:
        print(exc, file=sys.stderr)
        return 1
    if args.from_stage == "fetch" and not args.skip_doctor:
        checks = doctor.run_checks()
        if not doctor.ready(checks, needs_martinize=s.has_atomistic_protein):
            print(doctor.report(checks), file=sys.stderr)
            return 1
    try:
        result = build(s, from_stage=args.from_stage, until_stage=args.until,
                       grompp_check=not args.no_grompp, force=args.force)
    except (BuildError, RuntimeError) as exc:
        print(f"\nBUILD FAILED: {exc}", file=sys.stderr)
        return 2
    out = result["output_dir"]
    print(f"\nDone. GROMACS inputs: {out / 'gromacs'}\nProvenance: {result['lock']}")
    return 0


LIB_ALIASES = {
    "m2": "martini2", "martini2": "martini2", "martini22": "martini2", "elnedyn22": "martini2",
    "m3": "martini3_lipidome", "martini3": "martini3_lipidome", "lipidome": "martini3_lipidome",
    "martini3_lipidome": "martini3_lipidome",
    "m3-2021": "martini3_2021", "release2021": "martini3_2021", "martini3_2021": "martini3_2021",
}


def _lib(name: str) -> str:
    try:
        return LIB_ALIASES[name.lower()]
    except KeyError:
        raise SystemExit(f"unknown force field '{name}'; use one of: m2, m3 (lipidome), m3-2021")


def _entry_line(lib: str, name: str) -> str:
    e = lipids.lookup(name, lib)
    if e is None:
        return f"  {lipids.LIBRARY_LABELS[lib]:<44} -"
    charge = f"{e['charge']:+g}" if e.get("charge") is not None else "?"
    coby = e["coby"]["status"]
    return (f"  {lipids.LIBRARY_LABELS[lib]:<44} {name:<6} {e.get('chemistry') or '':<18} "
            f"charge {charge:<4} COBY: {coby}")


def cmd_lipids(args) -> int:
    if args.action == "show":
        for name in args.names:
            print(f"{name}")
            for lib in lipids.LIBRARIES:
                print(_entry_line(lib, name))
            for src in lipids.LIBRARIES:
                if lipids.lookup(name, src) is None:
                    continue
                for dst in lipids.LIBRARIES:
                    if dst != src:
                        eq = lipids.equivalents(name, src, dst)
                        if eq and eq != [name]:
                            print(f"  {src} {name} -> {dst}: {', '.join(eq)}")
            for lib in lipids.LIBRARIES:
                note = (lipids.lookup(name, lib) or {}).get("note")
                if note:
                    print(f"  note ({lib}): {note}")
            print()
        return 0
    if args.action == "convert":
        src, dst = _lib(args.source), _lib(args.target)
        status = 0
        for name in args.names:
            eq = lipids.equivalents(name, src, dst)
            if lipids.lookup(name, src) is None:
                print(f"{name}: not a {lipids.LIBRARY_LABELS[src]} lipid")
                status = 1
            elif not eq:
                print(f"{name} ({lipids.describe(name, src)}) -> no equivalent in {lipids.LIBRARY_LABELS[dst]}")
                status = 1
            else:
                print(f"{lipids.describe(name, src)}")
                q_src = lipids.lookup(name, src).get("charge")
                for n in eq:
                    e = lipids.lookup(n, dst)
                    flag = ""
                    if q_src is not None and e.get("charge") is not None and abs(q_src - e["charge"]) > 1e-6:
                        flag = f"  !! CHARGE CHANGES {q_src:+g} -> {e['charge']:+g}"
                    print(f"   -> {lipids.describe(n, dst)}   [COBY: {e['coby']['status']}]{flag}")
        return status
    if args.action == "report":
        text = lipids.report_markdown()
        if args.names == ["-"]:
            print(text, end="")
        else:
            from pathlib import Path
            Path(args.names[0]).write_text(text)
            print(f"wrote {args.names[0]}")
        return 0
    if args.action == "search":
        lib = _lib(args.ff) if args.ff else None
        hits = lipids.search(" ".join(args.names), lib)
        for lname, name, e in hits[: args.limit]:
            print(f"{lname:<18} {name:<6} {e.get('common') or e.get('chemistry') or ''}")
        if len(hits) > args.limit:
            print(f"... {len(hits) - args.limit} more (use --limit)")
        return 0 if hits else 1
    return 2


def _protein_atoms(path: str):
    from .structure import DUMMY_RESNAMES, WATER_RESNAMES, read_pdb
    atoms, _ = read_pdb(path)
    return [a for a in atoms if a.record == "ATOM" and not a.is_hydrogen
            and a.resname not in DUMMY_RESNAMES | WATER_RESNAMES]


def cmd_orient(args) -> int:
    from . import orient, viewer
    from .structure import write_pdb
    atoms = _protein_atoms(args.pdb)
    if not atoms:
        print(f"no protein atoms in {args.pdb}", file=sys.stderr)
        return 1
    if args.action == "guess":
        placed, record = orient.hydrophobic(atoms, topology=args.topology, pore_correction=args.pore_correction)
        normal, centre = orient.membrane_in_input_frame(record)
        pl = record["placement"]
        summary = {"placement": pl, "membrane_normal_input_frame": normal, "membrane_centre_input_frame_A": centre,
                   "spec": {"protein": {"orientation": "manual",
                                        "manual": {"normal": [round(v, 4) for v in normal],
                                                   "center_point": [round(v, 3) for v in centre]}}}}
        if args.output:
            write_pdb(args.output, placed)
        if args.html:
            viewer.write_html(Path(args.html), atoms, record, thickness_A=pl["thickness_A"], title=Path(args.pdb).name)
        if args.json:
            print(json.dumps(summary, indent=2))
        else:
            print(f"{pl['verdict']}")
            print(f"  hydrophobic thickness {pl['thickness_A']:g} A, transfer energy {pl['energy_kcal']:g} kcal/mol, "
                  f"orientation spread {pl['uncertainty_deg']:g} deg")
            print(f"  TM segments: {', '.join(pl['tm_segments']) or '-'}")
            print(f"  cytoplasmic side (-z): {pl['kr_inside']} vs {pl['kr_outside']} Lys/Arg (positive-inside rule)")
            for n in pl["notes"]:
                print(f"  note: {n}")
            print("\nUse it in a spec (or start from it with `orientation: hydrophobic`):\n")
            print(yaml.safe_dump(summary["spec"], sort_keys=False, default_flow_style=None))
            if args.output:
                print(f"oriented structure (membrane normal = z, centre at z = 0): {args.output}")
            if args.html:
                print(f"viewer: {args.html}")
        return 0 if not pl["verdict"].startswith("not a") else 3
    # view
    if args.start == "guess":
        _, record = orient.hydrophobic(atoms, topology=args.topology, pore_correction=args.pore_correction)
        thick = record["placement"]["thickness_A"]
    else:
        record, thick = {"method": "input frame"}, 30.0
    out = Path(args.html or (Path(args.pdb).stem + "_orientation.html"))
    viewer.write_html(out, atoms, record, thickness_A=thick, title=Path(args.pdb).name)
    print(f"wrote {out} - open it in a browser, adjust, and copy the spec lines")
    return 0


def cmd_selfassembly(args) -> int:
    from . import selfassembly
    try:
        res = selfassembly.check_build(Path(args.build), Path(args.gro) if args.gro else None)
    except (selfassembly.SelfAssemblyError, FileNotFoundError) as exc:
        print(exc, file=sys.stderr)
        return 1
    (Path(args.build) / "self_assembly_result.json").write_text(json.dumps(res, indent=2) + "\n")
    if args.json:
        print(json.dumps(res, indent=2))
        return 0 if res["bilayer_formed"] else 3
    print("bilayer formed" if res["bilayer_formed"] else "NO complete bilayer")
    for a in res["axes"]:
        print(f"  normal {a['axis']}: tail ends in core {a['tail_ends_in_core']:.0%}, heads in leaflets "
              f"{a['heads_in_leaflets']:.0%}, in-plane coverage {a['in_plane_coverage']:.0%}")
    if res["bilayer_formed"]:
        print(f"  thickness (P-P) {res['thickness_PP_nm']} nm, lipids per leaflet {res['lipids_per_leaflet']}")
        if "protein_spans_bilayer" in res:
            print(f"  protein spans the bilayer: {res['protein_spans_bilayer']} "
                  f"({res['protein_backbone_beads_in_core']} backbone beads in the core), "
                  f"backbone RMSD to start {res['protein_backbone_rmsd_to_start_nm']} nm")
            print("\nRebuild with COBY in this orientation (any composition, also asymmetric):\n")
            print(yaml.safe_dump(res["rebuild_spec"], sort_keys=False, default_flow_style=None))
    if res.get("advice"):
        print(f"advice: {res['advice']}")
    print(f"details: {Path(args.build) / 'self_assembly_result.json'}")
    return 0 if res["bilayer_formed"] else 3


SKILL_DIRS = {
    # agent: (user-level, project-level) skill directories (Agent Skills / SKILL.md format)
    "claude": ("~/.claude/skills", ".claude/skills"),
    "codex": ("~/.agents/skills", ".agents/skills"),
    "gemini": ("~/.gemini/skills", ".gemini/skills"),
    "agents": ("~/.agents/skills", ".agents/skills"),  # shared location read by several agents
}


def skill_source() -> Path:
    src = Path(__file__).resolve().parents[2] / "skills" / "mdbuddy"
    if not (src / "SKILL.md").exists():
        raise SystemExit("skill files not found next to the package: install MDBuddy from a clone "
                         "(git clone ...; pip install -e .) or copy skills/mdbuddy by hand")
    return src


def cmd_skill(args) -> int:
    import shutil
    src = skill_source()
    targets = []
    agents = list(SKILL_DIRS) if args.agent == "all" else [args.agent]
    if args.dest:
        targets.append(Path(args.dest).expanduser() / "mdbuddy")
    else:
        for agent in agents:
            base = SKILL_DIRS[agent][0 if args.scope == "user" else 1]
            root = Path(base).expanduser() if args.scope == "user" else Path(args.project) / base
            targets.append(root / "mdbuddy")
    for dest in dict.fromkeys(targets):
        if dest.exists():
            if not args.force:
                print(f"{dest} exists; use --force to replace it")
                continue
            shutil.rmtree(dest)
        shutil.copytree(src, dest)
        print(f"installed skill -> {dest}")
    print("Restart your agent so it picks the skill up. For tools without skill support, use the MCP server "
          "(`mdbuddy mcp-config <agent>`) or point the agent at AGENTS.md.")
    return 0


def cmd_mcp(args) -> int:
    from . import mcp_server
    mcp_server.serve()
    return 0


def cmd_mcp_config(args) -> int:
    from . import mcp_server
    print(mcp_server.client_config(args.agent))
    return 0


def cmd_stages(args) -> int:
    for s in STAGES:
        print(s)
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="mdbuddy", description="Reproducible Martini membrane-protein "
                                     "system builder for GROMACS")
    parser.add_argument("--version", action="version", version=f"mdbuddy {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("doctor", help="check that required programs and packages are installed")
    p.add_argument("--json", action="store_true", help="machine-readable output")
    p.add_argument("--network", action="store_true", help="also test access to OPM, RCSB and GitHub")
    p.set_defaults(func=cmd_doctor)

    p = sub.add_parser("init", help="write a commented spec template")
    p.add_argument("name")
    p.add_argument("--pdb-id")
    p.add_argument("-o", "--output")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_init)

    p = sub.add_parser("validate", help="check a spec and print it with all defaults filled in")
    p.add_argument("spec")
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser("build", help="build the system described by a spec")
    p.add_argument("spec")
    p.add_argument("--from", dest="from_stage", default="fetch", choices=STAGES,
                   help="resume from this stage (earlier results are reused)")
    p.add_argument("--until", default="check", choices=STAGES, help="stop after this stage")
    p.add_argument("--no-grompp", action="store_true", help="skip the final gmx grompp check")
    p.add_argument("--skip-doctor", action="store_true", help="do not check the environment first")
    p.add_argument("--force", action="store_true", help="overwrite an existing output directory")
    p.set_defaults(func=cmd_build)

    p = sub.add_parser("lipids", help="look up lipid names across Martini 2, Martini 3 and the M3 lipidome")
    p.add_argument("action", choices=("show", "convert", "search", "report"),
                   help="show NAME...: meaning in every version; convert NAME... --from m2 --to m3; "
                        "search TEXT [--ff m3]; report FILE|-: Markdown list of every changed/clashing name")
    p.add_argument("names", nargs="*", default=["-"])
    p.add_argument("--from", dest="source", default="m2")
    p.add_argument("--to", dest="target", default="m3")
    p.add_argument("--ff", help="restrict search to one force field (m2, m3, m3-2021)")
    p.add_argument("--limit", type=int, default=40)
    p.set_defaults(func=cmd_lipids)

    p = sub.add_parser("orient", help="place a protein without OPM data: hydrophobic guess or interactive viewer")
    p.add_argument("action", choices=("guess", "view"),
                   help="guess: educated guess from exposed hydrophobicity; view: write the interactive HTML viewer")
    p.add_argument("pdb", help="atomistic protein structure (PDB)")
    p.add_argument("-o", "--output", help="guess: write the oriented structure here")
    p.add_argument("--html", help="write the viewer here (guess and view)")
    p.add_argument("--start", choices=("guess", "input"), default="guess", help="view: initial placement")
    p.add_argument("--topology", choices=("positive_inside", "as_found", "flip"), default="positive_inside")
    p.add_argument("--pore-correction", choices=("auto", "on", "off"), default="auto")
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_orient)

    p = sub.add_parser("selfassembly", help="check a finished self-assembly run and get the orientation")
    p.add_argument("action", choices=("check",))
    p.add_argument("build", help="the build directory (contains spec.resolved.yaml and gromacs/)")
    p.add_argument("--gro", help="final frame (default: gromacs/step7_self_assembly.gro)")
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_selfassembly)

    p = sub.add_parser("skill", help="install the MDBuddy agent skill (SKILL.md) for an AI coding agent")
    p.add_argument("action", choices=("install",))
    p.add_argument("--agent", choices=(*SKILL_DIRS, "all"), default="all")
    p.add_argument("--scope", choices=("user", "project"), default="user")
    p.add_argument("--project", default=".", help="project directory for --scope project")
    p.add_argument("--dest", help="install into this skills directory instead (other agents)")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_skill)

    p = sub.add_parser("mcp", help="run the MCP server on stdio (for any MCP-capable agent)")
    p.set_defaults(func=cmd_mcp)

    p = sub.add_parser("mcp-config", help="print the MCP configuration snippet for an agent")
    p.add_argument("agent", choices=("claude", "codex", "gemini", "cursor", "vscode", "claude-desktop"))
    p.set_defaults(func=cmd_mcp_config)

    p = sub.add_parser("stages", help="list pipeline stages")
    p.set_defaults(func=cmd_stages)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
