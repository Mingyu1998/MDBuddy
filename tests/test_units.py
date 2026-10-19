import json
import textwrap

import pytest

from mdbuddy import assemble, cgprotein, martinize, mdp, membrane, prepare, spec, topology
from mdbuddy.structure import Atom, opm_half_thickness, read_pdb


def mkspec(**over):
    base = {"name": "toy", "protein": {"source": "opm", "id": "1BL8"}}
    for k, v in over.items():
        base.setdefault(k, {}).update(v) if isinstance(v, dict) else base.__setitem__(k, v)
    return spec.from_dict(base)


# --- spec -------------------------------------------------------------------

def test_spec_defaults_filled():
    s = mkspec()
    assert s.protein["id"] == "1bl8"
    assert s.protein["orientation"] == "opm"
    assert s.cg["elastic"] == {"force": 700.0, "lower": 0.5, "upper": 0.9}
    assert s.leaflets() == {"upper": {"POPC": 1}, "lower": {"POPC": 1}}


def test_spec_martini22_elastic_defaults_and_partial_override():
    s = mkspec(coarse_grain={"forcefield": "martini22", "elastic": {"force": 800}})
    assert s.cg["elastic"] == {"force": 800, "lower": 0.5, "upper": 0.9}
    assert s.ff_info["coby_params"] == "M2"


def test_spec_collects_all_problems():
    with pytest.raises(spec.SpecError) as exc:
        spec.from_dict({"name": "bad name", "protein": {"source": "rcsb"},
                        "coarse_grain": {"forcefield": "martini22", "model": "go"},
                        "membrane": {"lipids": {"POPC": -1}}, "typo": 1})
    text = "\n".join(exc.value.problems)
    for needle in ("unknown key 'typo'", "'name'", "protein.id", "Go model", "ratio"):
        assert needle in text


def test_spec_file_source_defaults_to_as_is():
    s = spec.from_dict({"name": "x", "protein": {"source": "file", "path": "p.pdb"}})
    assert s.protein["orientation"] == "as_is"


def test_template_is_valid_yaml_spec(tmp_path):
    p = tmp_path / "t.yaml"
    p.write_text(spec.TEMPLATE.format(name="t", pdb_id="1bl8"))
    s = spec.load(p)
    assert s.name == "t"
    assert s.output_dir() == tmp_path / "t"


# --- structure / prepare ----------------------------------------------------

def test_pdb_line_roundtrip():
    line = "ATOM     12  CA  LYS A  37      11.104  -6.500  12.002  1.00 20.00           C"
    a = Atom.from_pdb_line(line)
    assert (a.name, a.resname, a.chain, a.resid, a.x, a.element) == ("CA", "LYS", "A", 37, 11.104, "C")
    b = Atom.from_pdb_line(a.to_pdb_line())
    assert (b.name, b.resname, b.chain, b.resid, b.z) == ("CA", "LYS", "A", 37, 12.002)


def test_opm_remark_parsed():
    assert opm_half_thickness(["REMARK      1/2 of bilayer thickness:   16.2"]) == 16.2
    assert opm_half_thickness(["REMARK nothing"]) is None


def test_prepare_cleans_and_reports(opm_pdb, tmp_path):
    out = tmp_path / "clean.pdb"
    report = prepare.prepare(mkspec(), opm_pdb, out)
    atoms, _ = read_pdb(out)
    names = {(a.resname, a.name) for a in atoms}
    assert ("HOH", "O") not in names and ("DUM", "O") not in names and ("RET", "C1") not in names
    assert not any(a.name == "H" for a in atoms)
    assert sum(1 for a in atoms if a.name == "CB" and a.resname == "LYS") == 2  # one altloc per chain
    assert report["incomplete_residues"] == []
    assert report["chains"] == ["A", "B"]
    assert report["opm_half_thickness_A"] == 16.2
    assert report["dropped"]["hetatm:RET"] == 1
    assert len(report["gaps"]) == 1 and report["gaps"][0]["chain"] == "B"


def test_prepare_chain_selection_and_errors(opm_pdb, tmp_path):
    report = prepare.prepare(mkspec(protein={"chains": ["B"]}), opm_pdb, tmp_path / "b.pdb")
    assert report["chains"] == ["B"]
    with pytest.raises(prepare.PrepareError):
        prepare.prepare(mkspec(protein={"chains": ["Z"]}), opm_pdb, tmp_path / "z.pdb")


def test_prepare_keep_hetatm(opm_pdb, tmp_path):
    report = prepare.prepare(mkspec(protein={"keep_hetatm": ["RET"]}), opm_pdb, tmp_path / "k.pdb")
    assert report["kept_hetatm"] == ["RET"]


# --- martinize command --------------------------------------------------------

def test_martinize_command_elastic_merges_oligomer():
    cmd = martinize.build_command(mkspec(), "in.pdb", ["A", "B", "C", "D"], ["-dssp", "mkdssp"])
    joined = " ".join(cmd)
    assert "-ff martini3001" in joined
    assert "-elastic -ef 700 -el 0.5 -eu 0.9" in joined
    assert "-merge A,B,C,D" in joined
    assert "-p backbone -pf 1000" in joined


def test_martinize_command_go_and_no_merge():
    s = mkspec(coarse_grain={"model": "go", "merge_chains": "none"})
    joined = " ".join(martinize.build_command(s, "in.pdb", ["A", "B"], []))
    assert "-go -go-eps 9.414" in joined and "-elastic" not in joined and "-merge" not in joined


def test_martinize_command_m22_and_custom_merge():
    s = mkspec(coarse_grain={"forcefield": "elnedyn22", "merge_chains": [["A", "B"], ["C", "D"]]})
    joined = " ".join(martinize.build_command(s, "in.pdb", ["A", "B", "C", "D"], []))
    assert "-ff elnedyn22" in joined and "-el 0 " in joined
    assert "-merge A,B -merge C,D" in joined


def test_single_chain_never_merged():
    assert martinize.merge_args(mkspec(), ["A"]) == []


# --- COBY arguments -----------------------------------------------------------

def _beads(zs):
    return [Atom("ATOM", i, "BB", "", "ALA", "A", i, "", 0.0, 10.0 * (i % 2), z) for i, z in enumerate(zs)]


def test_box_auto_and_explicit():
    beads = _beads([-20.0, 20.0])  # 4 nm tall, 1 nm wide
    assert membrane.box_size(mkspec(), beads) == [7.0, 7.0, 10.0]
    assert membrane.box_size(mkspec(membrane={"box": [12, 12, 14]}), beads) == [12.0, 12.0, 14.0]


def test_protein_offset_puts_membrane_centre_at_box_centre():
    beads = _beads([-10.0, 30.0])  # COG z = 10 A above the OPM membrane centre
    assert membrane.protein_offset_nm(mkspec(), beads) == pytest.approx(1.0)
    assert membrane.protein_offset_nm(mkspec(), beads, 10.0) == pytest.approx(0.0)


def test_membrane_argument_asymmetric():
    s = mkspec(membrane={"upper": {"POPC": 7, "CHOL": 3}, "lower": {"POPE": 1}, "lower_apl": 0.55})
    arg = membrane.membrane_arg(s)
    assert arg == ("leaflet:upper lipid:POPC:7:params:LTF lipid:CHOL:3:params:default apl:0.6 "
                   "leaflet:lower lipid:POPE:1:params:LTF apl:0.55")
    s21 = mkspec(membrane={"lipids": {"POPC": 1}}, forcefield_files={"martini3_lipids": "release2021"})
    assert membrane.membrane_arg(s21) == ("leaflet:upper lipid:POPC:1:params:default apl:0.6 "
                                          "leaflet:lower lipid:POPC:1:params:default apl:0.6")
    assert s.lipid_names() == ["POPC", "CHOL", "POPE"]


def test_solvation_and_protein_arguments():
    s = mkspec(solvent={"salt_molarity": 0.3, "cation": "K"})
    assert membrane.solvation_arg(s) == "solv:W pos:K neg:CL salt_molarity:0.3"
    assert membrane.protein_arg(s, "../02_cg/cg.pdb", [("protein_0", 1)], 1.23456) == \
        "file:../02_cg/cg.pdb moleculetypes:protein_0:1 cz:1.235"


# --- topology -------------------------------------------------------------------

def test_topology_parser_follows_includes_and_ifdefs(tmp_path):
    (tmp_path / "toppar").mkdir()
    (tmp_path / "toppar" / "ff.itp").write_text("[ defaults ]\n1 1\n")
    (tmp_path / "toppar" / "mol.itp").write_text(textwrap.dedent("""
        [ moleculetype ]
        PROT 1
        [ atoms ]
        1 P5 1 LYS BB 1 1.0 ; comment
        #ifdef EXTRA
        2 P5 1 LYS SC1 2 0.0
        #else
        2 P5 1 LYS SC2 2 0.5
        #endif
    """))
    top = tmp_path / "topol.top"
    topology.write_top(top, ["toppar/ff.itp", "toppar/mol.itp"], [("PROT", 3)], "sys")
    info = topology.parse(top)
    assert info.defaults_count == 1
    assert info.molecules == [("PROT", 3)]
    assert info.moleculetypes["PROT"].atom_names == ["BB", "SC2"]
    assert info.moleculetypes["PROT"].charge == pytest.approx(1.5)
    assert topology.parse(top, defines={"EXTRA"}).moleculetypes["PROT"].atom_names == ["BB", "SC1"]


def test_merge_consecutive():
    assert topology.merge_consecutive([("W", 2), ("W", 3), ("NA", 0), ("CL", 1)]) == [("W", 5), ("CL", 1)]


# --- mdp ------------------------------------------------------------------------

def test_mdp_ladder_and_production_length():
    s = mkspec(simulation={"production_ns": 500})
    st = mdp.stages(s)
    assert [x.name for x in st][:2] == ["step6.0_minimization", "step6.1_minimization"]
    assert st[-1].name == "step7_production"
    assert st[-1].nsteps == 25_000_000
    fcs = [x.posres_fc for x in st if x.posres_fc]
    assert fcs == sorted(fcs, reverse=True)


def mdp_params(text):
    """Parse an mdp into {key: value}, dropping comments."""
    out = {}
    for line in text.splitlines():
        line = line.split(";", 1)[0].strip()
        if line:
            k, v = (x.strip() for x in line.split("=", 1))
            assert k not in out, f"duplicate key {k}"
            out[k] = v
    return out


def test_mdp_martini_nonbonded_in_every_stage():
    s = mkspec()
    for st in mdp.stages(s):
        p = mdp_params(mdp.render(st, s, posres_enabled=True))
        assert p["coulombtype"] == "reaction-field"
        assert p["epsilon_r"] == "15" and p["epsilon_rf"] == "0"
        assert p["rcoulomb"] == p["rvdw"] == "1.1"
        assert p["vdw-modifier"] == "Potential-shift" and p["cutoff-scheme"] == "Verlet"
        for obsolete in ("title", "ns-type", "ns_type"):
            assert obsolete not in p


def test_mdp_flexible_only_when_an_itp_defines_it():
    """grompp warns about a -D macro nothing uses, which fails the build's own -maxwarn 0.

    The Martini 2 force fields have no #ifdef FLEXIBLE block anywhere, so step6.0 must not
    write -DFLEXIBLE for them.
    """
    s = mkspec()
    st = {x.name: x for x in mdp.stages(s)}
    with_flex = mdp_params(mdp.render(st["step6.0_minimization"], s, True, flexible_available=True))
    without = mdp_params(mdp.render(st["step6.0_minimization"], s, True, flexible_available=False))
    assert "-DFLEXIBLE" in with_flex["define"]
    assert without["define"] == "-DPOSRES -DPOSRES_FC=1000"
    # the rest of the stage is unchanged
    assert without["integrator"] == "steep" and without["comm-mode"] == "none"
    # and the description stops claiming flexible bonds (both EM stages become the same run)
    assert mdp.title(st["step6.0_minimization"]) == "Energy minimisation, flexible bonds"
    assert mdp.title(st["step6.0_minimization"], False) == "Energy minimisation"
    assert mdp.title(st["step6.1_minimization"], False) == "Energy minimisation"
    assert mdp.title(st["step6.2_equilibration"], False) == "NPT, dt 2 fs"


def test_martinize_noscfix_only_for_martini2():
    """martinize2 >= 0.11 applies side-chain corrections by default; martini22 and elnedyn22
    do not define them and warn, which stops the run at -maxwarn 0."""
    for ff, expected in (("martini3", False), ("martini22", True), ("elnedyn22", True)):
        cmd = martinize.build_command(mkspec(coarse_grain={"forcefield": ff}), "in.pdb", ["A"], [])
        assert ("-noscfix" in cmd) is expected, ff


def test_ion_alias_follows_the_parameter_files(tmp_path):
    """Rename COBY's ion names only when the itps define the other spelling."""
    gmx = tmp_path / "g"
    (gmx / "toppar").mkdir(parents=True)
    aliases = assemble.MOLECULE_ALIASES["martini2"]

    def itp_with(*moltypes):
        (gmx / "toppar" / "ions.itp").write_text(
            "".join(f"[ moleculetype ]\n{m} 1\n[ atoms ]\n1 Qd 1 ION {m} 1 1.0\n" for m in moltypes))
        return ["toppar/ions.itp"]

    # COBY's and the pinned itp's spelling agree: nothing to do
    assert assemble.alias_renames(aliases, ["NA", "CL", "W"], gmx, itp_with("NA", "CL")) == {}
    # cgmartini.nl spelling in the itp, COBY's spelling in the topology: no rename either,
    # because the alias points the other way (NA+ -> NA) and NA is not defined
    assert assemble.alias_renames(aliases, ["NA", "CL"], gmx, itp_with("NA+", "CL-")) == {}
    # COBY naming the ions NA+/CL- against an itp that defines NA/CL: rename
    assert assemble.alias_renames(aliases, ["NA+", "CL-"], gmx,
                                  itp_with("NA", "CL")) == {"NA+": "NA", "CL-": "CL"}
    # both spellings defined: leave COBY's alone
    assert assemble.alias_renames(aliases, ["NA+"], gmx, itp_with("NA", "NA+")) == {}
    # Martini 3 has no aliases at all
    assert assemble.alias_renames({}, ["NA", "CL"], gmx, itp_with("NA", "CL")) == {}


def test_mdp_never_enables_free_energy():
    """mdrun cannot do free-energy perturbation with restricted bending angles, which Martini 3
    proteins (martinize2 backbone) and several lipids use: grompp accepts `free-energy = yes`,
    then mdrun stops with "Function type Restricted Angles does not support ... free energy".
    No generated stage may switch it on. Only step6.0 relaxes constraints into harmonic bonds.
    """
    fep_keys = ("free-energy", "free_energy", "couple-moltype", "couple_moltype",
                "init-lambda", "init_lambda", "sc-alpha", "sc_alpha")
    for s in (mkspec(), mkspec(protein={"source": "file", "path": "p.pdb",
                                        "orientation": "self_assembly"})):
        flexible = []
        for st in mdp.stages(s):
            p = mdp_params(mdp.render(st, s, posres_enabled=True))
            for key in fep_keys:
                assert key not in p, f"{st.name} sets {key}"
            if "-DFLEXIBLE" in p["define"]:
                flexible.append(st.name)
        assert flexible == ["step6.0_minimization"]


def test_mdp_stage_details():
    s = mkspec()
    st = {x.name: x for x in mdp.stages(s)}
    flex = mdp_params(mdp.render(st["step6.0_minimization"], s, posres_enabled=True))
    assert flex["define"] == "-DPOSRES -DPOSRES_FC=1000 -DFLEXIBLE"
    assert flex["comm-mode"] == "none"  # no COM removal in a minimisation
    em = mdp_params(mdp.render(st["step6.1_minimization"], s, posres_enabled=False))
    assert em["define"] == "" and "tcoupl" not in em
    eq = mdp_params(mdp.render(st["step6.2_equilibration"], s, posres_enabled=True))
    assert eq["define"] == "-DPOSRES -DPOSRES_FC=1000"
    assert eq["tc-grps"] == "Protein MEMBRANE SOLVENT" and eq["ref-t"] == "310 310 310"
    assert eq["comm-grps"] == "SOLU_MEMB SOLVENT"
    assert eq["pcoupl"] == "C-rescale" and eq["refcoord-scaling"] == "com"
    assert eq["gen-vel"] == "yes" and eq["gen-seed"] == "2026" and eq["continuation"] == "no"
    prod = mdp_params(mdp.render(st["step7_production"], s, posres_enabled=True))
    assert prod["define"] == "" and "refcoord-scaling" not in prod
    assert prod["pcoupl"] == "Parrinello-Rahman" and prod["tau-p"] == "12.0"
    assert prod["nstxout-compressed"] == "25000" and prod["gen-vel"] == "no"


def test_output_interval_validation():
    assert mdp_params(mdp.render(mdp.production_stage(100), mkspec(simulation={"output_ps": 100}),
                                 True))["nstxout-compressed"] == "5000"
    with pytest.raises(spec.SpecError, match="output_ps"):
        mkspec(simulation={"output_ps": 3})


# --- validation -------------------------------------------------------------------

def _write_system(tmp_path, gro_names, molecules):
    (tmp_path / "ff.itp").write_text(textwrap.dedent("""
        [ defaults ]
        1 1
        [ moleculetype ]
        W 1
        [ atoms ]
        1 W 1 W W 1 0
        [ moleculetype ]
        NA 1
        [ atoms ]
        1 Q 1 NA NA 1 1
    """))
    top = tmp_path / "topol.top"
    topology.write_top(top, ["ff.itp"], molecules, "s")
    gro = tmp_path / "s.gro"
    lines = ["t", str(len(gro_names))]
    lines += [f"{i:>5}{n:<5}{n:>5}{i:>5}{0:8.3f}{0:8.3f}{0:8.3f}" for i, n in enumerate(gro_names, 1)]
    lines.append("   5.00000   5.00000   5.00000")
    gro.write_text("\n".join(lines) + "\n")
    return gro, topology.parse(top)


def test_validate_ok_and_charge_warning(tmp_path):
    gro, info = _write_system(tmp_path, ["W", "W", "NA"], [("W", 2), ("NA", 1)])
    rep = assemble.validate(gro, [("W", 2), ("NA", 1)], info)
    assert rep["errors"] == [] and rep["net_charge"] == 1 and rep["warnings"]


def test_validate_catches_count_and_name_errors(tmp_path):
    gro, info = _write_system(tmp_path, ["W", "W", "NA"], [("W", 4)])
    assert "4 atoms" in assemble.validate(gro, [("W", 4)], info)["errors"][0]
    gro, info = _write_system(tmp_path, ["NA", "W"], [("W", 1), ("NA", 1)])
    assert "mismatch" in assemble.validate(gro, [("W", 1), ("NA", 1)], info)["errors"][0]
    gro, info = _write_system(tmp_path, ["W"], [("POPC", 1)])
    assert "not defined" in assemble.validate(gro, [("POPC", 1)], info)["errors"][0]


# --- CLI / doctor -------------------------------------------------------------------

def test_cli_init_and_validate(tmp_path, capsys):
    from mdbuddy.cli import main
    spec_path = tmp_path / "demo.yaml"
    assert main(["init", "demo", "--pdb-id", "2bl2", "-o", str(spec_path)]) == 0
    assert main(["validate", str(spec_path), "--json"]) == 0
    out = json.loads(capsys.readouterr().out.split("\n", 1)[1])
    assert out["valid"] and out["resolved"]["protein"]["id"] == "2bl2"
    spec_path.write_text("name: demo\nprotein: {source: opm}\n")
    assert main(["validate", str(spec_path), "--json"]) == 1


def test_doctor_reports_missing_required(monkeypatch):
    from mdbuddy import doctor, tools
    monkeypatch.setattr(tools, "module_version", lambda name, dist=None: None)
    monkeypatch.setattr(tools, "find_program", lambda name: None)
    monkeypatch.setattr(tools, "find_gmx", lambda: None)
    monkeypatch.setattr(tools, "find_dssp", lambda: None)
    monkeypatch.setattr(tools, "has_module", lambda name: False)
    checks = doctor.run_checks()
    assert not doctor.ready(checks)
    data = json.loads(doctor.report(checks, as_json=True))
    missing = {c["name"] for c in data["checks"] if c["status"] == "missing"}
    assert {"martinize2", "COBY", "vermouth", "dssp"} <= missing
    assert all(c["fix"] for c in data["checks"] if c["status"] == "missing")


# --- run.sh -----------------------------------------------------------------

RUN_SH_STAGES = ["step6.0_minimization", "step6.1_minimization",
                 "step6.2_equilibration", "step6.7_equilibration", "step7_production"]


def run_the_runner(tmp_path, *args):
    """Execute a generated run.sh against a gmx that only records its arguments.

    Returns (returncode, stages mdrun was asked for, [(stage, -c input)] grompp saw, stderr).
    """
    import os
    import stat
    import subprocess

    d = tmp_path / "gmx_dir"
    d.mkdir(exist_ok=True)
    (d / "run.sh").write_text(assemble.RUN_SH.format(
        name="toy", maxwarn=0, stages=" ".join(RUN_SH_STAGES)))
    (d / "run.sh").chmod(0o755)
    bindir = tmp_path / "runbin"
    bindir.mkdir(exist_ok=True)
    stub = bindir / "gmx"
    stub.write_text('#!/bin/bash\necho "$@" >> "$RUNLOG"\n')
    stub.chmod(stub.stat().st_mode | stat.S_IEXEC)
    log = d / "calls.txt"
    log.unlink(missing_ok=True)
    env = {**os.environ, "PATH": f"{bindir}{os.pathsep}{os.environ['PATH']}", "RUNLOG": str(log)}
    proc = subprocess.run(["bash", str(d / "run.sh"), *args], cwd=d, env=env,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    mdruns, grompps = [], []
    for line in (log.read_text().splitlines() if log.exists() else []):
        f = line.split()
        if f[0] == "mdrun":
            mdruns.append(f[f.index("-deffnm") + 1])
        else:
            grompps.append((f[f.index("-f") + 1].removesuffix(".mdp"), f[f.index("-c") + 1]))
    return proc.returncode, mdruns, grompps, proc.stderr


def test_run_sh_stage_range(tmp_path):
    # no arguments: every stage, each starting from the previous one's output
    rc, ran, grompps, _ = run_the_runner(tmp_path)
    assert rc == 0 and ran == RUN_SH_STAGES
    assert [c for _, c in grompps] == ["step5_input.gro"] + [f"{s}.gro" for s in RUN_SH_STAGES[:-1]]

    # one argument: from that stage to the end, continuing from the stage before it
    rc, ran, grompps, _ = run_the_runner(tmp_path, "step6.2")
    assert rc == 0 and ran == RUN_SH_STAGES[2:]
    assert grompps[0] == ("step6.2_equilibration", "step6.1_minimization.gro")

    # two arguments: inclusive range, so production is left out
    rc, ran, _, _ = run_the_runner(tmp_path, "step6.0", "step6.7")
    assert rc == 0 and ran == RUN_SH_STAGES[:-1]

    # empty start means "up to"
    rc, ran, _, _ = run_the_runner(tmp_path, "", "step6.1")
    assert rc == 0 and ran == RUN_SH_STAGES[:2]

    # a single stage
    rc, ran, _, _ = run_the_runner(tmp_path, "step7", "step7")
    assert rc == 0 and ran == ["step7_production"]


def test_run_sh_rejects_bad_stage_range(tmp_path):
    rc, ran, _, err = run_the_runner(tmp_path, "step9")
    assert rc != 0 and ran == [] and "no stage matches 'step9'" in err

    rc, ran, _, err = run_the_runner(tmp_path, "", "nonsense")
    assert rc != 0 and ran == [] and "no stage matches 'nonsense'" in err

    rc, ran, _, err = run_the_runner(tmp_path, "step7", "step6.0")
    assert rc != 0 and ran == [] and "comes after" in err


def test_readme_run_examples_name_real_stages(tmp_path):
    """run.sh errors on a stage name that does not exist, so the README's examples must come
    from the ladder in hand: a self-assembly build has no step6.7."""
    for data, last, before in (
            ({"name": "t", "protein": {"source": "opm", "id": "1bl8"}},
             "step7_production", "step6.7_equilibration"),
            ({"name": "t", "protein": {"source": "opm", "id": "1bl8",
                                       "orientation": "self_assembly"}},
             "step7_self_assembly", "step6.2_equilibration")):
        s = spec.from_dict(data)
        stage_list = mdp.stages(s)
        text = assemble.run_examples(stage_list)
        names = {st.name for st in stage_list}
        assert f"./run.sh {last}" in text and f"./run.sh '' {before}" in text
        # every stage name the examples mention is really in the ladder
        for line in text.splitlines():
            for word in line.split("#")[0].split():
                if word.startswith("step"):
                    assert word in names, f"{word} is not a stage of this build"


def test_run_sh_quotes_survive_the_template(tmp_path):
    """RUN_SH is a \"\"\"...\"\"\" literal, so a \\n in it would become a real newline in the
    generated script. Every line of the script must stay on its own line."""
    text = assemble.RUN_SH.format(name="t", maxwarn=0, stages="step6.0_minimization step7_production")
    assert "\\n" not in text.replace("\n", "")  # no stray escapes left either
    for line in text.splitlines():
        assert line.count("'") % 2 == 0, f"unbalanced quote: {line!r}"


# --- stiff bonds vs the time step ------------------------------------------

STIFF_ITP = """
[ moleculetype ]
STIFF 1
[ atoms ]
1 P5 1 ALA BB 1 0.0
2 P5 1 ALA SC1 2 0.0
[ bonds ]
1 2 1 0.350 150000
"""

SOFT_ITP = STIFF_ITP.replace("150000", "1250").replace("STIFF", "SOFT")


def test_fastest_bond_matches_gromacs_formula(tmp_path):
    """GROMACS estimates 2*pi*sqrt(mu/k); for two 72 amu beads at k = 150000 kJ/mol/nm^2 it
    prints 9.7e-02 ps, which is what a real OMP island itp hit at a 20 fs step."""
    p = tmp_path / "stiff.itp"
    p.write_text(STIFF_ITP)
    info = topology.parse(p)
    mt = info.moleculetypes["STIFF"]
    assert len(mt.bonds) == 1 and mt.atom_types == ["P5", "P5"]
    i, j, k, period = topology.fastest_bond(mt, {"P5": 72.0})
    assert (i, j, k) == (1, 2, 150000.0)
    assert abs(period - 0.0973) < 5e-4          # ps, as grompp reports it
    assert abs(period / 5 * 1000 - 19.5) < 0.1  # max safe time step, fs
    # a normal Martini bond is nowhere near the limit
    p2 = tmp_path / "soft.itp"
    p2.write_text(SOFT_ITP)
    _, _, _, soft_period = topology.fastest_bond(
        topology.parse(p2).moleculetypes["SOFT"], {"P5": 72.0})
    assert soft_period / 5 * 1000 > 100  # fs


def test_stiff_bond_warning(tmp_path):
    p = tmp_path / "stiff.itp"
    p.write_text(STIFF_ITP)
    info = topology.parse(p)
    mols = [("STIFF", 1)]
    sim = {"timestep_fs": 20.0, "output_ps": 500.0}
    msg = cgprotein.stiff_bond_warning(info, mols, sim, {"P5": 72.0})
    assert msg and "150000" in msg and "19.5 fs" in msg
    # it suggests a step that is both safe AND valid for this output_ps
    assert "simulation.timestep_fs: 10" in msg
    # silent once the step is small enough
    assert cgprotein.stiff_bond_warning(info, mols, {**sim, "timestep_fs": 10.0}, {"P5": 72.0}) is None
    # silent for an ordinary Martini bond
    p2 = tmp_path / "soft.itp"
    p2.write_text(SOFT_ITP)
    assert cgprotein.stiff_bond_warning(topology.parse(p2), [("SOFT", 1)], sim, {"P5": 72.0}) is None
    # and for a molecule with no numeric bonds at all
    assert cgprotein.stiff_bond_warning(info, [("NOPE", 1)], sim) is None


def test_suggest_timestep_divides_output_ps():
    # 500 ps of output: 100 steps at 10 fs = 1 ps divides it; 15 fs (1.5 ps) does not
    assert cgprotein.suggest_timestep(19.5, 500.0) == 10.0
    assert cgprotein.suggest_timestep(19.5, 300.0) == 15.0
    assert cgprotein.suggest_timestep(100.0, 500.0) == 20.0
    assert cgprotein.suggest_timestep(0.5, 500.0) is None


def test_timestep_fs_caps_the_ladder_keeping_durations():
    base = {x.name: x for x in mdp.stages(mkspec())}
    capped = {x.name: x for x in mdp.stages(mkspec(simulation={"timestep_fs": 10}))}
    for name, st in base.items():
        c = capped[name]
        if st.dt is None:
            assert c.dt is None
            continue
        assert c.dt == min(st.dt, 0.010)
        # the stage still covers the same amount of simulated time
        assert abs(c.nsteps * c.dt - st.nsteps * st.dt) < 1e-6, name
    assert capped["step7_production"].dt == 0.010
    assert capped["step6.2_equilibration"].dt == 0.002  # already below the cap, untouched
    # titles (mdp header, folder README) name the dt actually used
    assert capped["step6.6_equilibration"].title == "NPT, dt 10 fs"
    assert capped["step6.7_equilibration"].title == "NPT, dt 10 fs, 5 ns"
    assert capped["step6.3_equilibration"].title == base["step6.3_equilibration"].title
    # output interval follows the new dt, so it is still 500 ps
    prod = mdp_params(mdp.render(capped["step7_production"],
                                 mkspec(simulation={"timestep_fs": 10}), posres_enabled=False))
    assert prod["dt"] == "0.01" and prod["nstenergy"] == "50000"


def test_timestep_fs_validation():
    for bad in (0.5, 41, "20"):
        with pytest.raises(spec.SpecError, match="timestep_fs"):
            mkspec(simulation={"timestep_fs": bad})
    # output_ps must still be a whole number of 100-step intervals at the chosen dt
    with pytest.raises(spec.SpecError, match="output_ps"):
        mkspec(simulation={"timestep_fs": 15})          # 500 ps / 1.5 ps is not an integer
    assert mkspec(simulation={"timestep_fs": 15, "output_ps": 300}).sim["timestep_fs"] == 15


# --- user-supplied force-field files ---------------------------------------

def _m2_registry(tmp_path, monkeypatch):
    """A tiny stand-in for the pinned Martini 2 set, served from local files."""
    import hashlib
    from mdbuddy import forcefields
    d = tmp_path / "pinned"
    d.mkdir(exist_ok=True)
    files = []
    for name, body, role, key in (
            ("martini_v2.2.itp", "[ defaults ]\n1 1\n\n[ atomtypes ]\nP5 72.0 0.000 A 0.0 0.0\n", "core", ""),
            ("martini_v2.0_ions.itp", "[ moleculetype ]\nNA 1\n[ atoms ]\n1 P5 1 ION NA 1 1.0\n", "core", ""),
            ("martini_v2.0_POPE_02.itp",
             "[ moleculetype ]\nPOPE 1\n[ atoms ]\n1 P5 1 POPE PO4 1 0.0\n", "optional", "POPE"),
            ("martini_v2.0_POPG_02.itp",
             "[ moleculetype ]\nPOPG 1\n[ atoms ]\n1 P5 1 POPG PO4 1 0.0\n", "optional", "POPG")):
        p = d / name
        p.write_text(body)
        files.append(forcefields.FFFile(name, p.as_uri(),
                                        hashlib.sha256(p.read_bytes()).hexdigest(), role, key))
    monkeypatch.setitem(forcefields.REGISTRY, "martini2", files)
    monkeypatch.setenv("MDBUDDY_CACHE", str(tmp_path / "cache"))
    return d


def test_extra_itp_supersedes_the_pinned_lipid_file(tmp_path, monkeypatch):
    """The cgmartini martini_v2.0_lipids_all file defines POPE and POPG as well as LPS, so
    adding it used to define those moleculetypes twice and stop the build."""
    from mdbuddy import forcefields
    _m2_registry(tmp_path, monkeypatch)
    mine = tmp_path / "my_lipids.itp"
    mine.write_text("[ moleculetype ]\nPOPE 1\n[ atoms ]\n1 P5 1 POPE PO4 1 0.0\n"
                    "[ moleculetype ]\nREMP 1\n[ atoms ]\n1 P5 1 REMP PO1 1 0.0\n")
    s = mkspec(coarse_grain={"forcefield": "martini22"},
               membrane={"lipids": {"POPE": 9, "POPG": 1}},
               forcefield_files={"extra_itps": [str(mine)]})
    notes = []
    records = forcefields.install_forcefield(s, tmp_path / "toppar", log=notes.append)
    names = [r["file"] for r in records]
    assert "martini_v2.0_POPE_02.itp" not in names   # superseded by my_lipids.itp
    assert "martini_v2.0_POPG_02.itp" in names       # my file does not define POPG
    assert "my_lipids.itp" in names
    assert any("POPE comes from your my_lipids.itp" in n for n in notes)
    # no moleculetype ends up defined twice
    info = topology.TopologyInfo()
    for n in names:
        topology.parse(tmp_path / "toppar" / n, info=info)
    assert info.duplicates == []


def test_replace_a_pinned_forcefield_file(tmp_path, monkeypatch):
    """forcefield_files.replace: a modified core force field (e.g. a modified-LJ martini_v2.2
    for an accelerated LPS model) in the same include position, recorded as unpinned."""
    from mdbuddy import forcefields
    _m2_registry(tmp_path, monkeypatch)
    mod = tmp_path / "martini_v2.2_modlj2.itp"
    mod.write_text("[ defaults ]\n1 1\n\n[ atomtypes ]\nP5 72.0 0.000 A 0.0 0.0\n"
                   "RQa 72.0 0.000 A 0.0 0.0\n")
    s = mkspec(coarse_grain={"forcefield": "martini22"},
               membrane={"lipids": {"POPE": 1}},
               forcefield_files={"replace": {"martini_v2.2.itp": str(mod)}})
    records = forcefields.install_forcefield(s, tmp_path / "toppar")
    core = records[0]
    assert core["file"] == "martini_v2.2.itp"        # same name and include position
    assert core["source"].startswith("replace:")
    assert core["pinned"] is False and core["replaces_pinned_sha256"]
    assert core["sha256"] != core["replaces_pinned_sha256"]
    # the file really is the user's one
    assert "RQa" in (tmp_path / "toppar" / "martini_v2.2.itp").read_text()


def test_replace_rejects_bad_targets(tmp_path, monkeypatch):
    from mdbuddy import forcefields
    _m2_registry(tmp_path, monkeypatch)
    mod = tmp_path / "mod.itp"
    mod.write_text("[ defaults ]\n1 1\n")
    s = mkspec(coarse_grain={"forcefield": "martini22"}, membrane={"lipids": {"POPE": 1}},
               forcefield_files={"replace": {"not_a_pinned_name.itp": str(mod)}})
    with pytest.raises(forcefields.ForceFieldError, match="not a pinned file name"):
        forcefields.install_forcefield(s, tmp_path / "toppar")
    s2 = mkspec(coarse_grain={"forcefield": "martini22"}, membrane={"lipids": {"POPE": 1}},
                forcefield_files={"replace": {"martini_v2.2.itp": str(tmp_path / "nope.itp")}})
    with pytest.raises(forcefields.ForceFieldError, match="file not found"):
        forcefields.install_forcefield(s2, tmp_path / "toppar")
