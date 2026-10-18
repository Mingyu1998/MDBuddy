"""Self-assembly: lipid suitability, the soup layout, and the bilayer analysis."""

import json
import textwrap

import pytest
import yaml

np = pytest.importorskip("numpy")

from mdbuddy import lipids, mdp, selfassembly, spec  # noqa: E402
from mdbuddy.structure import Atom  # noqa: E402


def test_suitability_rules():
    lib = "martini3_lipidome"
    ok = lipids.self_assembly_check({"POPC": 1}, {"POPC": 1}, lib, 310)
    assert not ok.errors
    assert lipids.self_assembly_check({"POPC": 1}, {"POPE": 1}, lib, 310).errors  # asymmetric
    lps = lipids.self_assembly_check({"RAMP": 1}, {"RAMP": 1}, "martini2", 310)
    assert any("does not self-assemble" in e for e in lps.errors)
    warm = lipids.self_assembly_check({"DPPC": 1}, {"DPPC": 1}, lib, 310)
    assert any("phase" in w for w in warm.warnings)
    assert not any("phase" in w for w in lipids.self_assembly_check({"DPPC": 1}, {"DPPC": 1}, lib, 323).warnings)


def test_self_assembly_mdp_protocol():
    s = spec.from_dict({"name": "sa", "protein": {"source": "opm", "id": "1bl8", "orientation": "self_assembly"},
                        "self_assembly": {"run_ns": 500}})
    st = mdp.stages(s)
    assert [x.name for x in st][-1] == "step7_self_assembly" and st[-1].nsteps == 25_000_000
    text = mdp.render(st[-1], s, posres_enabled=True)
    assert "pcoupltype               = isotropic" in text and "-DPOSRES" in text
    assert "compressibility          = 3e-4 " in text


def _write_system(d, prot_xyz, lip_xyz, box, n_lip):
    """topol.top + itps + gro for PROT (BB beads) and n_lip LIP lipids (PO4 C1A C2A)."""
    (d / "ff.itp").write_text(textwrap.dedent(f"""
        [ defaults ]
        1 1
        [ moleculetype ]
        PROT 1
        [ atoms ]
    """) + "".join(f"{i + 1} P5 {i + 1} {'LYS' if i < 3 else 'LEU'} BB {i + 1} 0\n" for i in range(len(prot_xyz)))
        + textwrap.dedent("""
        [ moleculetype ]
        LIP 1
        [ atoms ]
        1 Q5 1 LIP PO4 1 0
        2 C1 1 LIP C1A 2 0
        3 C1 1 LIP C2A 3 0
    """))
    (d / "topol.top").write_text(f'#include "ff.itp"\n[ system ]\nt\n[ molecules ]\nPROT 1\nLIP {n_lip}\n')


def _gro(path, prot_xyz, lip_xyz, box):
    lines = ["t", str(len(prot_xyz) + len(lip_xyz))]
    i = 0
    for p in prot_xyz:
        i += 1
        lines.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (1, "PROT", "BB", i, *p))
    for k, p in enumerate(lip_xyz):
        i += 1
        lines.append("%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (2 + k // 3, "LIP", ("PO4", "C1A", "C2A")[k % 3], i, *p))
    lines.append("%10.5f%10.5f%10.5f" % tuple(box))
    path.write_text("\n".join(lines) + "\n")


def _bilayer_along_y(L, yc, n_side=10):
    pts = []
    g = (np.arange(n_side) + 0.5) * L / n_side
    for x in g:
        for z in g:
            for sign in (1, -1):
                pts += [[x, yc + sign * 2.0, z], [x, yc + sign * 1.1, z], [x, yc + sign * 0.3, z]]
    return np.array(pts) % L


def test_analysis_detects_bilayer_and_rejects_soup(tmp_path):
    L = 10.0
    prot = np.array([[5.0, 2.0 + 0.3 * i, 5.0] for i in range(20)])  # spans y = 2..7.7
    lip = _bilayer_along_y(L, 5.0)
    _write_system(tmp_path, prot, lip, [L] * 3, len(lip) // 3)
    _gro(tmp_path / "bilayer.gro", prot, lip, [L] * 3)
    res = selfassembly.analyse(tmp_path / "bilayer.gro", tmp_path / "topol.top", ["LIP"], ["PROT"])
    assert res["bilayer_formed"] and res["normal_axis"] == "y"
    assert res["thickness_PP_nm"] == pytest.approx(4.0, abs=0.05)
    soup = np.random.default_rng(1).uniform(0, L, size=lip.shape)
    _gro(tmp_path / "soup.gro", prot, soup, [L] * 3)
    assert not selfassembly.analyse(tmp_path / "soup.gro", tmp_path / "topol.top", ["LIP"], ["PROT"])["bilayer_formed"]


def test_check_build_maps_orientation_to_input_frame(tmp_path):
    L, t0 = 10.0, np.array([5.0, 5.0, 5.0])
    prot_cg_A = np.array([[0.0, -27.0 + 3.0 * i, 0.0] for i in range(20)])  # protein along y in the cg frame
    prot = prot_cg_A / 10.0 + t0
    build = tmp_path / "b"
    (build / "gromacs").mkdir(parents=True)
    (build / "02_cg").mkdir()
    lip_final = _bilayer_along_y(L, 5.5)  # membrane centre 0.5 nm above the box centre
    _write_system(build / "gromacs", prot, lip_final, [L] * 3, len(lip_final) // 3)
    _gro(build / "gromacs" / "step5_input.gro", prot, np.random.default_rng(2).uniform(0, L, lip_final.shape), [L] * 3)
    _gro(build / "gromacs" / "step7_self_assembly.gro", prot, lip_final, [L] * 3)
    (build / "02_cg" / "cg.pdb").write_text("".join(
        Atom("ATOM", i + 1, "BB", "", "LYS" if i < 3 else "LEU", "A", i + 1, "", *p).to_pdb_line() + "\n"
        for i, p in enumerate(prot_cg_A)))
    (build / "spec.resolved.yaml").write_text(yaml.safe_dump({"membrane": {"lipids": {"LIP": 1}, "upper": None,
                                                                          "lower": None}}))
    (build / ".mdbuddy_state.json").write_text(json.dumps({
        "protein_molecules": [["PROT", 1]],
        "orientation": {"method": "manual", "rotation": np.eye(3).tolist(), "origin_A": [0, 0, 0], "z_shift_A": 0}}))
    res = selfassembly.check_build(build)
    assert res["bilayer_formed"] and res["protein_spans_bilayer"]
    n = np.array(res["membrane_normal_input_frame"])
    assert abs(n[1]) > 0.999
    assert n[1] > 0  # K/R beads (resid 1-3) sit at low y -> cytoplasm at -normal
    assert res["membrane_centre_input_frame_A"][1] == pytest.approx(5.0, abs=0.2)  # (5.5 - 5.0) nm in cg frame
    assert res["rebuild_spec"]["protein"]["orientation"] == "manual"
