"""Reading .itp/.top files and writing the final GROMACS topology.

The parser is deliberately small: it understands sections, #include, and
#ifdef/#ifndef/#else/#endif with #define, which is all Martini files use.
It extracts per-moleculetype atom names and charges for validation.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

_SECTION = re.compile(r"^\s*\[\s*([A-Za-z_0-9]+)\s*\]")


@dataclass
class MoleculeType:
    name: str
    source: str
    atom_names: list[str] = field(default_factory=list)
    resnames: list[str] = field(default_factory=list)
    charges: list[float] = field(default_factory=list)

    @property
    def natoms(self) -> int:
        return len(self.atom_names)

    @property
    def charge(self) -> float:
        return sum(self.charges)


@dataclass
class TopologyInfo:
    moleculetypes: dict[str, MoleculeType] = field(default_factory=dict)
    duplicates: list[tuple[str, str, str]] = field(default_factory=list)
    defaults_count: int = 0
    molecules: list[tuple[str, int]] = field(default_factory=list)
    includes: list[str] = field(default_factory=list)
    missing_includes: list[str] = field(default_factory=list)


def _strip(line: str) -> str:
    return line.split(";", 1)[0].strip()


def parse(path: str | Path, defines: Optional[set[str]] = None,
          info: Optional[TopologyInfo] = None) -> TopologyInfo:
    """Parse a .top or .itp file, following #include relative to each file."""
    info = info or TopologyInfo()
    defines = set(defines or ())
    _parse_file(Path(path), defines, info, [None], [None])
    return info


def _parse_file(path: Path, defines: set[str], info: TopologyInfo,
                section: list, current: list) -> None:
    active: list[bool] = []  # stack of #if states
    with open(path) as fh:
        lines = fh.readlines()
    for raw in lines:
        line = _strip(raw)
        if not line:
            continue
        if line.startswith("#"):
            parts = line.split()
            directive = parts[0]
            if directive in ("#ifdef", "#ifndef"):
                cond = parts[1] in defines
                active.append(cond if directive == "#ifdef" else not cond)
            elif directive == "#else":
                if active:
                    active[-1] = not active[-1]
            elif directive == "#endif":
                if active:
                    active.pop()
            elif all(active):
                if directive == "#define" and len(parts) > 1:
                    defines.add(parts[1])
                elif directive == "#undef" and len(parts) > 1:
                    defines.discard(parts[1])
                elif directive == "#include":
                    target = parts[1].strip('"<>')
                    inc = path.parent / target
                    info.includes.append(str(inc))
                    if inc.exists():
                        _parse_file(inc, defines, info, section, current)
                    else:
                        info.missing_includes.append(f"{target} (from {path.name})")
            continue
        if not all(active):
            continue
        m = _SECTION.match(line)
        if m:
            section[0] = m.group(1).lower()
            if section[0] == "defaults":
                info.defaults_count += 1
            continue
        sec = section[0]
        fields = line.split()
        if sec == "moleculetype":
            name = fields[0]
            if name in info.moleculetypes:
                info.duplicates.append((name, info.moleculetypes[name].source, path.name))
            mt = MoleculeType(name=name, source=path.name)
            info.moleculetypes[name] = mt
            current[0] = mt
        elif sec == "atoms" and current[0] is not None and len(fields) >= 5:
            mt = current[0]
            mt.resnames.append(fields[3])
            mt.atom_names.append(fields[4])
            mt.charges.append(float(fields[6]) if len(fields) >= 7 else 0.0)
        elif sec == "molecules" and len(fields) >= 2:
            info.molecules.append((fields[0], int(fields[1])))


def read_molecules(top_path: str | Path) -> list[tuple[str, int]]:
    """Read only the [ molecules ] section of a .top file (no include following)."""
    molecules = []
    in_section = False
    for raw in Path(top_path).read_text().splitlines():
        line = _strip(raw)
        if not line:
            continue
        m = _SECTION.match(line)
        if m:
            in_section = m.group(1).lower() == "molecules"
            continue
        if in_section and not line.startswith("#"):
            name, count = line.split()[:2]
            molecules.append((name, int(count)))
    return molecules


def read_includes(top_path: str | Path) -> list[str]:
    """Return the #include targets written directly in a file (not recursive)."""
    out = []
    for raw in Path(top_path).read_text().splitlines():
        line = _strip(raw)
        if line.startswith("#include"):
            out.append(line.split(None, 1)[1].strip('"<>'))
    return out


def merge_consecutive(molecules: Iterable[tuple[str, int]]) -> list[tuple[str, int]]:
    out: list[tuple[str, int]] = []
    for name, count in molecules:
        if count == 0:
            continue
        if out and out[-1][0] == name:
            out[-1] = (name, out[-1][1] + count)
        else:
            out.append((name, count))
    return out


def write_top(path: Path, includes: list[str], molecules: list[tuple[str, int]],
              system_name: str, defines: Iterable[str] = (), header: Iterable[str] = ()) -> None:
    lines = [f"; {h}" for h in header]
    if lines:
        lines.append("")
    for d in defines:
        lines.append(f"#define {d}")
    if defines:
        lines.append("")
    lines += [f'#include "{inc}"' for inc in includes]
    lines += ["", "[ system ]", system_name, "", "[ molecules ]", "; name  count"]
    width = max((len(n) for n, _ in molecules), default=4)
    lines += [f"{name:<{width}}  {count}" for name, count in molecules]
    path.write_text("\n".join(lines) + "\n")
