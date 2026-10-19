import json
import textwrap

import pytest

from mdbuddy import assemble, martinize, mdp, membrane, prepare, spec, topology
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
