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
LIPIDOME_COMMIT = "2121527045696345069deb4d3df754724888a4f9"  # Martini-Force-Field-Initiative/M3-Lipid-Parameters

_M3 = (f"https://raw.githubusercontent.com/marrink-lab/martini-forcefields/{MARTINI_FF_COMMIT}"
       "/martini_forcefields/regular/v3.0.0/gmx_files/")
_COBY = f"https://raw.githubusercontent.com/MikkelDA/COBY/{COBY_COMMIT}/Tutorial/toppar/"
_LIPIDOME = (f"https://raw.githubusercontent.com/Martini-Force-Field-Initiative/M3-Lipid-Parameters/"
             f"{LIPIDOME_COMMIT}/ITPs/")


@dataclass(frozen=True)
class FFFile:
    name: str
    url: str
    sha256: str
    role: str  # "core" (always included), "optional" (by key) or "lipid" (when a membrane lipid needs it)
    key: str = ""  # short name used in forcefield_files.optional


REGISTRY: dict[str, list[FFFile]] = {
    "martini3": [
        FFFile("martini_v3.0.0.itp", _M3 + "martini_v3.0.0.itp",
               "0ecbf0ff334b2ec21dc5e3151aca42167b977bc2aed41eb2a44bd8c997754894", "core"),
        FFFile("martini_v3.0.0_ions_v1.itp", _M3 + "martini_v3.0.0_ions_v1.itp",
               "699aba132755b0ef13a0d7b41e6230a27777781bdc8075507413b89499c2be4f", "core"),
        FFFile("martini_v3.0.0_solvents_v1.itp", _M3 + "martini_v3.0.0_solvents_v1.itp",
               "e6a78414b317c38a1a816eded4aee6082a5da1ed229a126359f0cc92314dcb4e", "core"),
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

# Martini 3 lipid parameter sets (forcefield_files.martini3_lipids). File order matters for the
# lipidome: martini_v3.0.0_ffbonded_v2.itp defines the bonded macros and must precede the lipids.
_LIPIDOME_FILES = {
    "martini_v3.0.0_DOTAP_v2.itp": "57f68c5ff3abaca51fcca8d1a30ad2f6abb9eab53afcf71be7631becd7d1bd10",
    "martini_v3.0.0_ceramides_v2.itp": "aa10f614161b5722b9075a257aaa6638cde896d37b94fa5a2d1366cc6e6fd74b",
    "martini_v3.0.0_diglycerides_v2.itp": "3cd1dc52e3c8d99076d6f3da0a9c1c048500383d3a0afd7534b4956a212c4392",
    "martini_v3.0.0_etherphospholipids_PA_v2.itp": "9694b27ab08899aeeec4f0d08246b43e7d6591dec9c57aac98ab508a0f46bd13",
    "martini_v3.0.0_etherphospholipids_PC_v2.itp": "8177f9a373d91011b49183ffc2d1f738f4c5909d26ddbcdca859352799966b37",
    "martini_v3.0.0_etherphospholipids_PE_v2.itp": "b7bbe0260030cad8fe7dc85c88a853b3f96ef36c6613123e76e3b8f58c26fe92",
    "martini_v3.0.0_etherphospholipids_PG_v2.itp": "b5db1de6763a4037ada37e11a452756de52b70a97aef0b2fad52fc23bcb426f4",
    "martini_v3.0.0_etherphospholipids_PS_v2.itp": "01d6d2d3220134f72e7e1c546f044ed235cb3fa6efe22e157d8467de825aa27a",
    "martini_v3.0.0_fattyacids_v2.itp": "f390b85fc16866fa202a7ba59d69113f34a613a7f60374bf7250741428c77763",
    "martini_v3.0.0_hydrocarbons_v2.itp": "272e3a4c25b0d5efee5a2293ba36c6c63f40d107de914455919c463906c900db",
    "martini_v3.0.0_monoglycerides_v2.itp": "d13bb2c58770aa37b46b9560cde5e48aa0d53bed9d2ce1e70e470a5881ae36e3",
    "martini_v3.0.0_phospholipids_2,2-BMP_v2.itp": "d04abc391956a9d15d790e7f0a7528e8a5556f9a13b77e76cde735910b980599",
    "martini_v3.0.0_phospholipids_3,3-BMP_v2.itp": "cbfbbaff2dc722b4a51907cfa4cbe9e189dae018409acc64d7736f615d7d750c",
    "martini_v3.0.0_phospholipids_CL_v2.itp": "e98cfa3afa1b0f13f14b8474bed73bc4154c5f946bf29871bce3dfad6f3461b0",
    "martini_v3.0.0_phospholipids_PA_v2.itp": "62c83e6eabde0d6d20f5364486eae098c69ed408199ca9adb900e0edc5a235c5",
    "martini_v3.0.0_phospholipids_PC_v2.itp": "43fc3e9976ea2a11af1b3df5d841ed1d2677dddd93a73a7e86304ad875ffb150",
    "martini_v3.0.0_phospholipids_PE_v2.itp": "497c06ba81b47e28b95125e9bc56272f77b15f3df9acd53cc7e22339a0cc2bb2",
    "martini_v3.0.0_phospholipids_PG_v2.itp": "37c65865088bd3944a50a646cc3489db040475d165a87595ad5b7f38a6c7c258",
    "martini_v3.0.0_phospholipids_PI_v2.itp": "98f868d024a4f4cf4c7d83508f171bac5c2d47a79085239a894e30ea21cc4351",
    "martini_v3.0.0_phospholipids_PS_v2.itp": "b6af8cefb4bf4e91bea5caec37bac23b2fd6b2a334603902cce2915950fedd55",
    "martini_v3.0.0_phospholipids_SM_v2.itp": "9ad7ed22cdf5829f076b0a386e5c5322bd6fb042ce57434fde9ea4308414f8fd",
    "martini_v3.0.0_plasmalogens_PA_v2.itp": "52cb6f83756b1781633fcf73e3bbe4bd3fb4557539a23977d6c2884ecc2e926e",
    "martini_v3.0.0_plasmalogens_PC_v2.itp": "55dc543d11bd6a097e144cc2fedf398e54b9a2d594798d2e82b688285c93eee9",
    "martini_v3.0.0_plasmalogens_PE_v2.itp": "68133ad0128b4279d7fc231a9aacdefbba6e6f574da606d11fa85ece230b568b",
    "martini_v3.0.0_plasmalogens_PG_v2.itp": "b5fe4473eaf8cb7b4c178c82378d4b61182aabbfae71530a195df8014325e5b3",
    "martini_v3.0.0_plasmalogens_PS_v2.itp": "6495c6edaf8eb9e803ef829fab3aad540075349735daf53fffe505643977319b",
    "martini_v3.0.0_triglycerides_v2.itp": "11efbcf2bda7fd2e225e697684bead384d7718962c88d03adfba2c37c8457089",
}
LIPID_SETS: dict[str, list[FFFile]] = {
    "lipidome": [FFFile("martini_v3.0.0_ffbonded_v2.itp", _LIPIDOME + "martini_v3.0.0_ffbonded_v2.itp",
                        "149f98c48e32a9aab95e205b36c10e68332349d831aec3980bc55b4b428ab3fb", "core")]
                + [FFFile(n, _LIPIDOME + n, h, "lipid") for n, h in sorted(_LIPIDOME_FILES.items())],
    "release2021": [FFFile("martini_v3.0.0_phospholipids_v1.itp", _M3 + "martini_v3.0.0_phospholipids_v1.itp",
                           "21fb2eb6d8a1d48854d179dc71981c1210f608486efb3e6644347391ba9cf392", "core")],
}

FAMILY = {"martini3": "martini3", "martini22": "martini2", "elnedyn22": "martini2"}


def cache_dir() -> Path:
    root = os.environ.get("MDBUDDY_CACHE") or os.path.join(
        os.environ.get("XDG_CACHE_HOME", os.path.expanduser("~/.cache")), "mdbuddy")
    return Path(root) / "forcefields"


class ForceFieldError(RuntimeError):
    pass


def _obtain(f: FFFile, dest_dir: Path, ff_dir: Optional[Path],
            replace: Optional[dict] = None) -> dict:
    dest = dest_dir / f.name
    origin = None
    if replace and f.name in replace:
        # forcefield_files.replace: a deliberately different version of a pinned file, kept in
        # the same include position. Unverified by definition, so the lock records it as such
        # together with the checksum it is standing in for.
        src = replace[f.name]
        shutil.copyfile(src, dest)
        return {"file": f.name, "source": f"replace:{src}", "sha256": sha256(dest),
                "pinned": False, "replaces_pinned_sha256": f.sha256}
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


def install_forcefield(spec, toppar: Path, log=None) -> list[dict]:
    """Copy all force-field files for ``spec`` into ``toppar``; return provenance records.

    The returned list is in include order: core parameter file first.
    """
    toppar.mkdir(parents=True, exist_ok=True)
    family = FAMILY[spec.cg["forcefield"]]
    wanted_optional = set(spec.ff_files["optional"])
    ff_dir = spec.resolve_path(spec.ff_files["ff_dir"]) if spec.ff_files["ff_dir"] else None
    note = log or (lambda _m: None)

    replace = {}
    for name, path in (spec.ff_files["replace"] or {}).items():
        src = spec.resolve_path(path)
        if not src.exists():
            raise ForceFieldError(f"forcefield_files.replace[{name}]: file not found: {src}")
        if name not in {f.name for f in REGISTRY[family]} | {f.name for fs in LIPID_SETS.values() for f in fs}:
            raise ForceFieldError(
                f"forcefield_files.replace: {name} is not a pinned file name for {family}. "
                f"Pinned names: {sorted(f.name for f in REGISTRY[family])}. "
                "To add a file rather than replace one, use forcefield_files.extra_itps.")
        replace[name] = src
    known_optional = {f.key for f in REGISTRY[family] if f.role == "optional"}

    # Martini 2: pull pinned lipid files automatically when the membrane uses them
    if family == "martini2":
        wanted_optional |= set(spec.lipid_names()) & known_optional
    unknown = wanted_optional - known_optional
    if unknown:
        raise ForceFieldError(f"forcefield_files.optional: unknown {sorted(unknown)}; "
                              f"available for {family}: {sorted(known_optional)}")

    # A pinned optional lipid file would duplicate a moleculetype an extra_itps file already
    # defines (e.g. the cgmartini martini_v2.0_lipids_all file defines POPE and POPG as well as
    # LPS). The user's own file wins, since they asked for it by name.
    from . import topology
    extra_moltypes: dict[str, str] = {}
    for extra in spec.ff_files["extra_itps"]:
        src = spec.resolve_path(extra)
        if src.exists():
            for mt in topology.parse(src).moleculetypes:
                extra_moltypes.setdefault(mt, src.name)

    files = [f for f in REGISTRY[family] if f.role == "core" or f.key in wanted_optional]
    superseded = [f for f in files if f.role == "optional" and f.key in extra_moltypes]
    for f in superseded:
        note(f"{f.key} comes from your {extra_moltypes[f.key]}, not the pinned {f.name}")
    files = [f for f in files if f not in superseded]
    if family == "martini3":
        from . import lipids
        lipid_set = LIPID_SETS[spec.ff_files["martini3_lipids"]]
        lib = lipids.library_for(spec)
        needed = set(lipids.itp_files([n for n in spec.lipid_names() if lipids.lookup(n, lib)], lib))
        chosen = [f for f in lipid_set if f.role == "core" or f.name in needed]
        # main parameter file first, then the lipid set (ffbonded before any lipid), then the rest
        files = files[:1] + [f for f in chosen if f.role == "core"] + files[1:] + \
            [f for f in chosen if f.role != "core"]

    records = [_obtain(f, toppar, ff_dir, replace) for f in files]

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
