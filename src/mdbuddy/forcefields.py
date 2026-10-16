"""Force-field parameter files: pinned sources, checksums and a local cache.

Files are downloaded from pinned commits (never a moving branch) and verified
against SHA-256 checksums, so two builds with the same MDBuddy version use
byte-identical parameters. Use ``forcefield_files.ff_dir`` to work offline from
a local copy, or ``extra_itps`` for files MDBuddy does not ship a source for
(e.g. most Martini 2 lipids, which are distributed via cgmartini.nl).
"""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .fetch import FetchError, download
from .provenance import sha256

MARTINI_FF_COMMIT = "784591ebdc91d762ed4df986c4650546c938f776"  # marrink-lab/martini-forcefields
COBY_COMMIT = "143dc533176ac1f9c04a7075ccfd5ffd330a16ff"  # MikkelDA/COBY

_M3 = (f"https://raw.githubusercontent.com/marrink-lab/martini-forcefields/{MARTINI_FF_COMMIT}"
       "/martini_forcefields/regular/v3.0.0/gmx_files/")
_COBY = f"https://raw.githubusercontent.com/MikkelDA/COBY/{COBY_COMMIT}/Tutorial/toppar/"


@dataclass(frozen=True)
class FFFile:
    name: str
    url: str
    sha256: str
    role: str  # "core" (always included) or "optional"
    key: str = ""  # short name used in forcefield_files.optional


REGISTRY: dict[str, list[FFFile]] = {
    "martini3": [
        FFFile("martini_v3.0.0.itp", _M3 + "martini_v3.0.0.itp",
               "0ecbf0ff334b2ec21dc5e3151aca42167b977bc2aed41eb2a44bd8c997754894", "core"),
        FFFile("martini_v3.0.0_ions_v1.itp", _M3 + "martini_v3.0.0_ions_v1.itp",
               "699aba132755b0ef13a0d7b41e6230a27777781bdc8075507413b89499c2be4f", "core"),
        FFFile("martini_v3.0.0_solvents_v1.itp", _M3 + "martini_v3.0.0_solvents_v1.itp",
               "e6a78414b317c38a1a816eded4aee6082a5da1ed229a126359f0cc92314dcb4e", "core"),
        FFFile("martini_v3.0.0_phospholipids_v1.itp", _M3 + "martini_v3.0.0_phospholipids_v1.itp",
               "21fb2eb6d8a1d48854d179dc71981c1210f608486efb3e6644347391ba9cf392", "core"),
        FFFile("martini_v3.0_sterols_v1.0.itp", _COBY + "martini_v3/martini_v3.0_sterols_v1.0.itp",
               "b7cf06babd2699bc54e8538a4bab4108378de9b2ab355eb13fa01a18b5f9b69f", "core"),
        FFFile("martini_v3.0.0_sugars_v1.itp", _M3 + "martini_v3.0.0_sugars_v1.itp",
               "342c8a7be7209a9d6f50c91a0adba3e16150863cd33411ca6bb18b4c2675bbbf", "optional", "sugars"),
        FFFile("martini_v3.0.0_small_molecules_v1.itp", _M3 + "martini_v3.0.0_small_molecules_v1.itp",
               "61eb12c3e70fa31af0503d91c3c4ec03c16e88414c884e694bd4092b0afe3cbc", "optional",
               "small_molecules"),
    ],
    "martini2": [
        FFFile("martini_v2.2.itp", _COBY + "M2_BacterialMembrane/Other/martini_v2.2.itp",
               "6caa464cdd21d8c4e24740080ede8dcb0bbacfd0cca700fd4f3f1067e78f975d", "core"),
        FFFile("martini_v2.0_ions.itp", _COBY + "M2_BacterialMembrane/Other/martini_v2.0_ions.itp",
               "9753494b8c26ea338c822bd9e2c5558030ad11a44c7684991971f73f680a318d", "core"),
        FFFile("martini_v2.0_POPE_02.itp", _COBY + "M2_BacterialMembrane/Lipids/martini_v2.0_POPE_02.itp",
               "d3692f096ecba50430cb9c12a7f9f10ace55657be1d99732df7c6e184fe81a61", "optional", "POPE"),
        FFFile("martini_v2.0_POPG_02.itp", _COBY + "M2_BacterialMembrane/Lipids/martini_v2.0_POPG_02.itp",
               "dce84ad8f1ae19841a119fea996f926ff67284368aacb2df66e7ba3c6143a5e8", "optional", "POPG"),
        FFFile("martini_v2.0_CDL0_01.itp", _COBY + "M2_BacterialMembrane/Lipids/martini_v2.0_CDL0_01.itp",
               "50b70b584e1d1f6cd834e6155a2b63eccf955fbca4009d41c43181db97564e0b", "optional", "CDL0"),
    ],
}

FAMILY = {"martini3": "martini3", "martini22": "martini2", "elnedyn22": "martini2"}


def cache_dir() -> Path:
    root = os.environ.get("MDBUDDY_CACHE") or os.path.join(
        os.environ.get("XDG_CACHE_HOME", os.path.expanduser("~/.cache")), "mdbuddy")
    return Path(root) / "forcefields"


class ForceFieldError(RuntimeError):
    pass


def _obtain(f: FFFile, dest_dir: Path, ff_dir: Optional[Path]) -> dict:
    dest = dest_dir / f.name
    origin = None
    if ff_dir is not None and (ff_dir / f.name).exists():
        shutil.copyfile(ff_dir / f.name, dest)
        origin = f"ff_dir:{ff_dir / f.name}"
    else:
        cached = cache_dir() / f.sha256[:16] / f.name
        if not cached.exists():
            cached.parent.mkdir(parents=True, exist_ok=True)
            tmp = cached.with_suffix(".part")
            try:
                download(f.url, tmp)
            except FetchError as exc:
                raise ForceFieldError(
                    f"{exc}\nCannot get {f.name}. Put a copy in a folder and set "
                    "forcefield_files.ff_dir to it (checksum is verified).") from exc
            tmp.rename(cached)
        shutil.copyfile(cached, dest)
        origin = f.url
    digest = sha256(dest)
    if digest != f.sha256:
        raise ForceFieldError(
            f"checksum mismatch for {f.name} from {origin}\n  expected {f.sha256}\n  got      {digest}\n"
            "Refusing to build with unverified parameters. If you deliberately use a different "
            "version, pass it through forcefield_files.extra_itps instead.")
    return {"file": f.name, "source": origin, "sha256": digest, "pinned": True}


def install_forcefield(spec, toppar: Path) -> list[dict]:
    """Copy all force-field files for ``spec`` into ``toppar``; return provenance records.

    The returned list is in include order: core parameter file first.
    """
    toppar.mkdir(parents=True, exist_ok=True)
    family = FAMILY[spec.cg["forcefield"]]
    wanted_optional = set(spec.ff_files["optional"])
    ff_dir = spec.resolve_path(spec.ff_files["ff_dir"]) if spec.ff_files["ff_dir"] else None
    known_optional = {f.key for f in REGISTRY[family] if f.role == "optional"}

    # Martini 2: pull pinned lipid files automatically when the membrane uses them
    if family == "martini2":
        wanted_optional |= set(spec.lipid_names()) & known_optional
    unknown = wanted_optional - known_optional
    if unknown:
        raise ForceFieldError(f"forcefield_files.optional: unknown {sorted(unknown)}; "
                              f"available for {family}: {sorted(known_optional)}")

    records = []
    for f in REGISTRY[family]:
        if f.role == "core" or f.key in wanted_optional:
            records.append(_obtain(f, toppar, ff_dir))

    for extra in spec.ff_files["extra_itps"]:
        src = spec.resolve_path(extra)
        if not src.exists():
            raise ForceFieldError(f"forcefield_files.extra_itps: file not found: {src}")
        dest = toppar / src.name
        if dest.exists() and any(r["file"] == src.name for r in records):
            raise ForceFieldError(f"extra itp {src.name} clashes with a pinned force-field file name")
        shutil.copyfile(src, dest)
        records.append({"file": src.name, "source": f"user:{src}", "sha256": sha256(dest), "pinned": False})
    return records


GO_ATOMTYPES_HOOK = '#ifdef GO_VIRT\n#include "go_atomtypes.itp"\n#endif\n'
GO_NBPARAMS_HOOK = '#ifdef GO_VIRT\n#include "go_nbparams.itp"\n#endif\n'


def patch_for_go(main_itp: Path) -> None:
    """Hook martinize2's Go-model files into the Martini 3 parameter file.

    martinize2 writes go_atomtypes.itp / go_nbparams.itp; they have to be
    included at the end of the [ atomtypes ] and [ nonbond_params ] sections.
    """
    text = main_itp.read_text()
    if "go_atomtypes.itp" in text:
        return
    marker = "[ nonbond_params ]"
    idx = text.find(marker)
    if idx < 0:
        raise ForceFieldError(f"{main_itp.name}: no [ nonbond_params ] section to hook the Go model into")
    text = text[:idx] + GO_ATOMTYPES_HOOK + "\n" + text[idx:]
    if not text.endswith("\n"):
        text += "\n"
    text += "\n" + GO_NBPARAMS_HOOK
    main_itp.write_text(text)


def check_molecules(spec, toppar: Path, files: list[str]) -> None:
    """Fail early if a lipid/solvent/ion in the spec has no moleculetype in the parameter files."""
    from . import topology

    info = topology.TopologyInfo()
    for f in files:
        topology.parse(toppar / f, info=info)
    s = spec.solvent
    wanted = {**{name: "membrane lipid" for name in spec.lipid_names()},
              s["water"]: "solvent.water", s["cation"]: "solvent.cation", s["anion"]: "solvent.anion"}
    missing = [f"{name} ({what})" for name, what in wanted.items() if name not in info.moleculetypes]
    if missing:
        family = FAMILY[spec.cg["forcefield"]]
        hint = ("Martini 2 lipid parameters (e.g. martini_v2.0_lipids_all_201506.itp) "
                "are distributed via cgmartini.nl; download them and list them in forcefield_files.extra_itps."
                if family == "martini2" else
                "Add the itp that defines them via forcefield_files.extra_itps, or enable an optional file "
                "(forcefield_files.optional: [sugars, small_molecules]).")
        raise ForceFieldError(f"not defined in any {family} parameter file: {', '.join(missing)}.\n{hint}")
