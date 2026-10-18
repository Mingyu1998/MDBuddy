"""Pre-coarse-grained proteins, membrane-only systems, per-lipid area, custom lipids and the coby section."""

import json
import textwrap

import pytest
import yaml

from mdbuddy import doctor, lipids, pipeline, spec, topology

ISLAND_ITP = """
[ moleculetype ]
ISL 1
[ atoms ]
1 P5 1 LEU BB  1 0
2 P5 1 LEU SC1 2 0
3 P5 2 ALA BB  3 0
"""


def write_island(d, box=(20.0, 20.0, 15.0)):
    """Two copies of a 3-bead CG protein in a .gro, with a .top that also includes a force field."""
    (d / "island.itp").write_text(textwrap.dedent(ISLAND_ITP))
    (d / "island.top").write_text('#include "martini_v3.0.0.itp"\n#include "island.itp"\n\n'
                                  "[ system ]\nisland\n\n[ molecules ]\nISL 2\n")
    lines = ["island", "6"]
    for m, (x, y) in enumerate(((5.0, 5.0), (12.0, 9.0))):
        for k, (rn, an) in enumerate((("LEU", "BB"), ("LEU", "SC1"), ("ALA", "BB"))):
            lines.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (2 * m + 1, rn, an, 3 * m + k + 1, x + 0.3 * k, y, 7.5))
    lines.append("%10.5f%10.5f%10.5f" % box)
    (d / "island.gro").write_text("\n".join(lines) + "\n")


def write_pdb_protein(d, cryst=False):
    (d / "small.itp").write_text(textwrap.dedent("""
        [ moleculetype ]
        SMALL 1
        [ atoms ]
        1 P5 1 LYS BB 1 1.0
        2 P5 2 GLY BB 2 0
    """))
    rows = ["CRYST1  100.000  100.000   80.000  90.00  90.00  90.00 P 1           1"] if cryst else []
    rows += ["ATOM      1  BB  LYS A   1      10.000  10.000  10.000  1.00  0.00",
             "ATOM      2  BB  GLY A   2      13.000  10.000  10.000  1.00  0.00", "END"]
    (d / "small.pdb").write_text("\n".join(rows) + "\n")


def load(tmp_path, data):
    p = tmp_path / "spec.yaml"
    p.write_text(yaml.safe_dump({"name": "sys", **data}))
    return spec.load(p)


def kwargs_of(s):
    return json.loads((s.output_dir() / "03_membrane" / "coby_kwargs.json").read_text())


def test_protein_islands_kept_in_place(tmp_path, fake_tools):
    write_island(tmp_path)
    write_pdb_protein(tmp_path)
    s = load(tmp_path, {
        "protein": {"source": "none"},
        "cg_proteins": [{"structure": "island.gro", "topology": "island.top"},
                        {"structure": "small.pdb", "itps": ["small.itp"], "moleculetypes": [["SMALL", 1]],
                         "placement": "center", "cx": 3.0, "rz": 90}],
        "membrane": {"box": "from_structure"},
    })
    logs = []
    pipeline.build(s, log=logs.append)
    kw = kwargs_of(s)
    assert kw["box"] == [20.0, 20.0, 15.0]
    assert kw["protein"][0] == "file:../02_cg/cg_protein_0.pdb moleculetypes:ISL:2 center_protein:False"
    assert kw["protein"][1] == "file:../02_cg/cg_protein_1.pdb moleculetypes:SMALL:1 cx:3 rz:90"
    cg = s.output_dir() / "02_cg"
    pdb = (cg / "cg_protein_0.pdb").read_text().splitlines()
    assert pdb[0].startswith("CRYST1  200.000  200.000  150.000") and (cg / "cg_protein_0_input.gro").exists()
    assert pdb[2][30:54] == "  53.000  50.000  75.000" and pdb[2][12:16] == " SC1"
    gmx = s.output_dir() / "gromacs"
    assert (gmx / "toppar" / "island.itp").exists() and (gmx / "toppar" / "small.itp").exists()
    assert not (gmx / "toppar" / "martini_v3.0.0.itp").read_text().startswith("[ moleculetype ]")
    info = topology.parse(gmx / "topol.top")
    assert info.molecules[:2] == [("ISL", 2), ("SMALL", 1)]
    ndx = (gmx / "index.ndx").read_text()
    protein = ndx.split("[ Protein ]")[1].split("[")[0].split()
    assert len(protein) == 8
    assert not (s.output_dir() / "02_cg" / "args.txt").exists()  # martinize2 never ran
    assert any("skipped: no atomistic protein" in line for line in logs)
    lock = json.loads((s.output_dir() / "mdbuddy.lock.json").read_text())
    rec = next(n for st in lock["stages"] if st["stage"] == "martinize" for n in st["notes"])["cg_proteins"]
    assert rec[0]["molecules"] == [["ISL", 2]] and rec[0]["sha256"]
    pipeline.build(s, from_stage="membrane", log=lambda m: None)  # resume re-reads the copied structures
    assert kwargs_of(s)["protein"] == kw["protein"]


def test_cg_protein_checks(tmp_path, fake_tools):
    write_island(tmp_path)
    write_pdb_protein(tmp_path)
    bad_count = load(tmp_path, {"protein": {"source": "none"}, "membrane": {"box": "from_structure"},
                                "cg_proteins": [{"structure": "island.gro", "itps": ["island.itp"],
                                                 "moleculetypes": [["ISL", 3]]}]})
    with pytest.raises(pipeline.BuildError, match="has 6 beads but"):
        pipeline.build(bad_count, log=lambda m: None)
    no_box = load(tmp_path, {"protein": {"source": "none"}, "membrane": {"box": [10, 10, 10]},
                             "cg_proteins": [{"structure": "small.pdb", "itps": ["small.itp"],
                                              "moleculetypes": [["SMALL", 1]]}]})
    with pytest.raises(pipeline.BuildError, match="placement 'keep' needs"):
        pipeline.build(no_box, log=lambda m: None, force=True)
    with pytest.raises(spec.SpecError) as err:
        load(tmp_path, {"protein": {"source": "none", "orientation": "self_assembly"},
                        "cg_proteins": [{"structure": "a.gro"}], "membrane": {"box": "from_structure"},
                        "coby": {"kwargs": {"box": [1, 1, 1], "protein_extra": "x"}}})
    msg = str(err.value)
    assert "self-assembly needs an atomistic protein" in msg
    assert "give 'topology'" in msg and "coby.kwargs may not set ['box', 'protein_extra']" in msg
    with pytest.raises(spec.SpecError, match="membrane-only system needs an explicit membrane.box"):
        load(tmp_path, {"protein": {"source": "none"}})
    with pytest.raises(spec.SpecError, match="more than one protein"):
        load(tmp_path, {"protein": {"source": "opm", "id": "1bl8"},
                        "cg_proteins": [{"structure": "small.pdb", "topology": "x.top"}]})


def test_membrane_only_patch_with_density_tuning(tmp_path, fake_tools):
    (tmp_path / "defs.py").write_text("lipid_defs = {}\n")
    s = load(tmp_path, {
        "protein": {"source": "none"},
        "membrane": {"lipids": {"POPC": 1}, "box": [30, 30, 12], "cz": -1.5, "lipid_apl": {"POPC": 0.68}},
        "coby": {"membrane_args": "hole:circle:cx:0:cy:0:radius:3", "solvation_args": "kick:0.05",
                 "kwargs": {"import_library": ["defs.py"], "verbose": 1}},
    })
    assert s.sim["coupling_groups"] == ["MEMBRANE", "SOLVENT"]
    pipeline.build(s, log=lambda m: None)
    kw = kwargs_of(s)
    assert "protein" not in kw
    memb = kw["membrane"]
    assert memb.startswith("cz:-1.5 hole:circle:cx:0:cy:0:radius:3 leaflet:upper lipid:POPC:1:params:")
    assert memb.count(":apl:0.68") == 2
    assert kw["solvation"].endswith(" kick:0.05")
    assert kw["import_library"] == ["user_files/defs.py"] and kw["verbose"] == 1
    assert (s.output_dir() / "03_membrane" / "user_files" / "defs.py").exists()
    ndx = (s.output_dir() / "gromacs" / "index.ndx").read_text()
    assert "[ MEMBRANE ]" in ndx
    mdp = (s.output_dir() / "gromacs" / "step7_production.mdp").read_text()
    assert "MEMBRANE SOLVENT" in mdp


def test_raw_membranes_replace_the_generated_one(tmp_path, fake_tools):
    s = load(tmp_path, {
        "protein": {"source": "none"},
        "membrane": {"box": [40, 40, 30]},
        "coby": {"membranes": ["lipid:POPC:1:params:LTF cz:-5",
                               "leaflet:upper lipid:POPC:2 lipid_extra:extra_val:3:name:POPC leaflet:lower "
                               "lipid:name:POPC cz:5"]},
    })
    assert s.lipid_names() == ["POPC"]
    pipeline.build(s, log=lambda m: None)
    kw = kwargs_of(s)
    assert kw["membrane"][0] == "lipid:POPC:1:params:LTF cz:-5"
    assert kw["membrane"][1].startswith("leaflet:upper lipid:POPC:2:params:LTF lipid_extra:extra_val:3:name:POPC:params:LTF")
    assert kw["membrane"][1].endswith("lipid:name:POPC:params:LTF cz:5")


def test_lipid_names_from_raw_strings(tmp_path):
    s = load(tmp_path, {"protein": {"source": "none"}, "membrane": {"box": [10, 10, 10]},
                        "coby": {"membranes": ["lipid:POPE:3:params:LTF lipid:name:TOCL:ratio:1 "
                                               "lipid_extra:name:CHOL:extra_val:2"]}})
    assert s.lipid_names() == ["POPE", "TOCL", "CHOL"]


def test_custom_lipids(tmp_path, fake_tools):
    (tmp_path / "lps.itp").write_text("[ moleculetype ]\nREMP 1\n[ atoms ]\n1 Q5 1 REMP PO1 1 0\n")
    assert "REMP" not in lipids.library("martini3_lipidome")
    base = {"protein": {"source": "none"}, "membrane": {"box": [10, 10, 10],
                                                        "upper": {"REMP": 1}, "lower": {"POPC": 1}},
            "forcefield_files": {"extra_itps": ["lps.itp"]}}
    s = load(tmp_path, base)
    assert pipeline.user_lipid_definitions(s) == {"REMP": "lps.itp"}
    with pytest.raises(pipeline.BuildError, match="COBY also needs its structure"):
        pipeline.build(s, log=lambda m: None)
    res = lipids.check(["REMP"], "martini3_lipidome", custom={"REMP": "lps.itp"}, coby_definitions=True)
    assert not res.errors and "custom lipid" in res.warnings[0]
    assert lipids.check(["REMP"], "martini3_lipidome", custom={"REMP": "lps.itp"}, coby_definitions={"X"}).errors
    (tmp_path / "remp.gro").write_text("x\n1\n    1REMP   PO1    1   0.000   0.000   0.000\n 1 1 1\n")
    base["coby"] = {"kwargs": {"molecule_import": "file:remp.gro moleculetype:REMP library_types:lipid"}}
    s = load(tmp_path, base)
    assert pipeline.coby_definitions(s) == {"REMP"}
    pipeline.build(s, log=lambda m: None, until_stage="membrane", force=True)
    kw = kwargs_of(s)
    assert kw["molecule_import"] == ["file:user_files/remp.gro moleculetype:REMP library_types:lipid params:default"]
    assert lipids.check(["REMP"], "martini3_lipidome").errors  # not defined anywhere


def test_doctor_does_not_need_martinize_without_protein():
    checks = [doctor.Check("martinize2", "missing", True, None, ""),
              doctor.Check("COBY", "ok", True, "1.0.15", "")]
    assert not doctor.ready(checks)
    assert doctor.ready(checks, needs_martinize=False)


def test_examples_validate():
    from pathlib import Path
    for path in sorted((Path(__file__).parent.parent / "examples").glob("*.yaml")):
        spec.load(path)  # raises SpecError on a problem
