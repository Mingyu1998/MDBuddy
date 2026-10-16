"""End-to-end plumbing test with stand-in tools (see conftest.py)."""

import json

import pytest

from mdbuddy import pipeline, spec, topology
from mdbuddy.forcefields import ForceFieldError, REGISTRY, install_forcefield


def make_spec(tmp_path, opm_pdb, **cg):
    s = {
        "name": "toy",
        "protein": {"source": "file", "path": str(opm_pdb), "orientation": "opm"},
        "coarse_grain": {"dssp": "none", **cg},
        "simulation": {"production_ns": 10},
    }
    p = tmp_path / "toy.yaml"
    import yaml
    p.write_text(yaml.safe_dump(s))
    return spec.load(p)


def test_full_build_with_stand_in_tools(tmp_path, opm_pdb, fake_tools):
    s = make_spec(tmp_path, opm_pdb)
    logs = []
    result = pipeline.build(s, log=logs.append)
    out = result["output_dir"]
    gmx = out / "gromacs"

    # martinize2 got the oligomer merged and the elastic network
    args = (out / "02_cg" / "args.txt").read_text()
    assert "-merge A,B" in args and "-elastic" in args

    # COBY call is saved as a runnable script with relative paths
    kwargs = json.loads((out / "03_membrane" / "coby_kwargs.json").read_text())
    assert kwargs["params"] == "default"
    assert kwargs["protein"].startswith("file:../02_cg/cg.pdb moleculetypes:protein_0:1")
    assert all(i.startswith("include:../gromacs/toppar/") for i in kwargs["itp_input"])
    assert "COBY.COBY(" in (out / "03_membrane" / "coby_build.py").read_text()

    # GROMACS folder
    for f in ["step5_input.gro", "topol.top", "index.ndx", "run.sh", "README.md",
              "step6.0_minimization.mdp", "step6.1_minimization.mdp", "step6.7_equilibration.mdp",
              "step7_production.mdp",
              "toppar/protein_0.itp", "toppar/martini_v3.0.0.itp"]:
        assert (gmx / f).exists(), f
    info = topology.parse(gmx / "topol.top")
    assert info.molecules == [("protein_0", 1), ("POPC", 4), ("W", 10), ("NA", 2), ("CL", 4)]
    ndx = (gmx / "index.ndx").read_text()
    for g in ("[ Protein ]", "[ MEMBRANE ]", "[ SOLU_MEMB ]", "[ SOLVENT ]", "[ BB ]"):
        assert g in ndx

    # provenance
    lock = json.loads((out / "mdbuddy.lock.json").read_text())
    assert [st["stage"] for st in lock["stages"]] == list(pipeline.STAGES)
    assert all(r["pinned"] for r in lock["forcefield_files"])
    assert lock["inputs"][0]["sha256"]
    assert (out / "spec.resolved.yaml").exists() and (out / "CITATIONS.md").exists()
    assert any("gap" in line for line in logs)
    check = lock["stages"][-1]["notes"][0]["grompp"]
    assert check.startswith("not run")


def test_resume_from_stage(tmp_path, opm_pdb, fake_tools):
    s = make_spec(tmp_path, opm_pdb)
    pipeline.build(s, log=lambda m: None)
    s2 = make_spec(tmp_path, opm_pdb, position_restraints="none")
    pipeline.build(s2, from_stage="gromacs", log=lambda m: None)
    eq = (s2.output_dir() / "gromacs" / "step6.2_equilibration.mdp").read_text()
    assert "POSRES" not in eq


def test_refuses_to_overwrite(tmp_path, opm_pdb, fake_tools):
    s = make_spec(tmp_path, opm_pdb)
    pipeline.build(s, log=lambda m: None)
    with pytest.raises(pipeline.BuildError):
        pipeline.build(s, log=lambda m: None)
    pipeline.build(s, log=lambda m: None, force=True)


def test_go_model_patches_forcefield(tmp_path, opm_pdb, fake_tools):
    s = make_spec(tmp_path, opm_pdb, model="go")
    pipeline.build(s, until_stage="forcefield", log=lambda m: None)
    text = (s.output_dir() / "gromacs" / "toppar" / "martini_v3.0.0.itp").read_text()
    assert text.index('#include "go_atomtypes.itp"') < text.index("[ nonbond_params ]")
    assert text.rstrip().endswith("#endif") and "go_nbparams.itp" in text


def test_checksum_mismatch_is_fatal(tmp_path, opm_pdb, fake_tools, monkeypatch):
    bad = [REGISTRY["martini3"][0].__class__(f.name, f.url, "0" * 64, f.role) for f in REGISTRY["martini3"]]
    monkeypatch.setitem(REGISTRY, "martini3", bad)
    s = make_spec(tmp_path, opm_pdb)
    with pytest.raises(ForceFieldError, match="checksum mismatch"):
        install_forcefield(s, tmp_path / "toppar")


def test_unknown_lipid_fails_before_martinize(tmp_path, opm_pdb, fake_tools):
    s = make_spec(tmp_path, opm_pdb)
    s.data["membrane"]["lipids"] = {"POPC": 1, "XXPC": 1}
    with pytest.raises(pipeline.BuildError, match="XXPC"):
        pipeline.build(s, log=lambda m: None)
    assert not (s.output_dir() / "00_input").exists()


def test_cross_version_name_is_explained(tmp_path, opm_pdb, fake_tools):
    s = make_spec(tmp_path, opm_pdb)
    s.data["membrane"]["lipids"] = {"POPC": 1, "CDL2": 1}
    with pytest.raises(pipeline.BuildError) as exc:
        pipeline.build(s, log=lambda m: None)
    assert "Martini 2.x name" in str(exc.value) and "TOCL" in str(exc.value)


def test_lipidome_files_included_on_demand(tmp_path, opm_pdb, fake_tools):
    s = make_spec(tmp_path, opm_pdb)
    pipeline.build(s, until_stage="forcefield", log=lambda m: None)
    top_files = [r["file"] for r in json.loads((s.output_dir() / "mdbuddy.lock.json").read_text())["forcefield_files"]]
    assert top_files[0] == "martini_v3.0.0.itp" and top_files[1] == "martini_v3.0.0_ffbonded_v2.itp"
    assert top_files[-1] == "martini_v3.0.0_phospholipids_PC_v2.itp"
