"""Lipid naming across Martini versions (data from the official lipidome naming table)."""

import textwrap

from mdbuddy import assemble, lipids, topology


def test_known_name_clashes():
    # same name, different molecule
    assert lipids.lookup("DLPC", "martini2")["chemistry"] == "di-C12:0-C14:0"
    assert lipids.lookup("DLPC", "martini3_lipidome")["chemistry"] == "di-C18:2"
    assert lipids.lookup("PIPC", "martini3_lipidome")["chemistry"] == "C16:0/20:2"
    assert "PIP2(3,4)" in lipids.lookup("POP2", "martini3_lipidome")["common"]
    assert "PIP2(4,5)" in lipids.lookup("POP6", "martini3_lipidome")["common"]


def test_equivalents_between_versions():
    m2, m3 = "martini2", "martini3_lipidome"
    assert lipids.equivalents("DIPC", m2, m3) == ["DLPC"]
    assert lipids.equivalents("DLPC", m2, m3) == ["DUPC", "DMPC"]  # M2 DLPC = C12:0 (or C14:0)
    assert lipids.equivalents("DPPC", m2, m3) == ["DPPC", "DSPC", "PSPC"]
    assert lipids.equivalents("PIPC", m2, m3) == ["PLPC", "SLPC"]
    assert lipids.equivalents("POP2", m2, m3) == ["POP6"]
    assert lipids.equivalents("CDL2", m2, m3) == ["TOCL"]
    assert lipids.equivalents("TOCL", m3, m2) == ["CDL2", "CDL1", "CDL0"]
    assert lipids.equivalents("DLPC", m3, m2) == ["DIPC"]
    assert lipids.equivalents("POPC", "martini3_2021", m2) == ["POPC"]


def test_check_errors_and_warnings():
    r = lipids.check(["POPC", "DLPC", "TOCL"], "martini3_lipidome")
    assert not r.errors
    assert any("DLPC means a different molecule" in w for w in r.warnings)
    assert any("TOCL" in w and "2-" in w for w in r.warnings)
    assert set(r.renames) == {"TOCL"}
    r = lipids.check(["CDL0", "DXPC", "XXXX"], "martini3_lipidome")
    text = "\n".join(r.errors)
    assert "Martini 2.x name" in text and "TOCL" in text  # CDL0 -> suggestion
    assert "COBY cannot build" in text                     # DXPC: bead count mismatch in COBY 1.0.15
    assert "XXXX" in text
    assert any("COBY cannot build" in e for e in lipids.check(["CDL0"], "martini2").errors)


def test_coby_params_and_files():
    assert lipids.coby_params("POPC", "martini3_lipidome") == "LTF"
    assert lipids.coby_params("CHOL", "martini3_lipidome") == "default"
    assert lipids.coby_params("POPC", "martini3_2021") == "default"
    assert lipids.itp_files(["POPC", "TOCL", "CHOL", "POPE"], "martini3_lipidome") == [
        "martini_v3.0.0_phospholipids_PC_v2.itp", "martini_v3.0.0_phospholipids_CL_v2.itp",
        "martini_v3.0_sterols_v1.0.itp", "martini_v3.0.0_phospholipids_PE_v2.itp"]


def test_fix_coby_names(tmp_path):
    (tmp_path / "l.itp").write_text(textwrap.dedent("""
        [ defaults ]
        1 1
        [ moleculetype ]
        TOCL 1
        [ atoms ]
        1 N6 1 TOCL GLC 1 0
        2 C1 1 TOCL D2A1 2 0
        3 C1 1 TOCL C3B2 3 0
    """))
    top = tmp_path / "t.top"
    topology.write_top(top, ["l.itp"], [("TOCL", 2)], "s")
    info = topology.parse(top)
    gro = tmp_path / "s.gro"
    names = ["GLC", "C2A1", "C3C2"] * 2
    gro.write_text("t\n6\n" + "".join(f"{1:>5}{'TOCL':<5}{n:>5}{i:>5}{0:8.3f}{0:8.3f}{0:8.3f}\n"
                                       for i, n in enumerate(names, 1)) + "   5.0   5.0   5.0\n")
    renames = {"TOCL": {1: ("C2A1", "D2A1"), 2: ("C3C2", "C3B2")}}
    assert assemble.fix_coby_names(gro, tmp_path / "none.pdb", [("TOCL", 2)], info, renames) == {"TOCL": 4}
    assert assemble.validate(gro, [("TOCL", 2)], info)["errors"] == []
