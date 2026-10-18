"""References to cite for a build. Written to CITATIONS.md next to the outputs."""

from __future__ import annotations

REFS = {
    "martini3": "Souza et al. (2021) Martini 3: a general purpose force field for coarse-grained molecular "
                "dynamics. Nat. Methods 18, 382-388. doi:10.1038/s41592-021-01098-3",
    "martini22": "de Jong et al. (2013) Improved parameters for the Martini coarse-grained protein force field. "
                 "J. Chem. Theory Comput. 9, 687-697. doi:10.1021/ct300646g",
    "martini2_lipids": "Marrink et al. (2007) The MARTINI force field: coarse grained model for biomolecular "
                       "simulations. J. Phys. Chem. B 111, 7812-7824. doi:10.1021/jp071097f",
    "elnedyn": "Periole et al. (2009) Combining an elastic network with a coarse-grained molecular force field. "
               "J. Chem. Theory Comput. 5, 2531-2543. doi:10.1021/ct9002114",
    "martinize2": "Kroon et al. Martinize2 and Vermouth: unified framework for topology generation. "
                  "eLife. doi:10.7554/eLife.90627",
    "go": "Souza et al. (2025) GoMartini 3: from large conformational changes in proteins to environmental bias "
          "corrections. Nat. Commun.",
    "coby": "Andreasen et al. (2025) Creating Coarse-Grained Systems with COBY: Toward Higher Accuracy of Complex "
            "Biological Systems. J. Chem. Inf. Model. doi:10.1021/acs.jcim.5c00069",
    "opm": "Lomize et al. (2012) OPM database and PPM web server: resources for positioning of proteins in "
           "membranes. Nucleic Acids Res. 40, D370-D376. doi:10.1093/nar/gkr703",
    "lipidome": "Pedersen et al. (2025) The Martini 3 Lipidome: Expanded and Refined Parameters Improve Lipid Phase "
                "Behavior. ACS Cent. Sci. doi:10.1021/acscentsci.5c00755",
    "sterols": "Borges-Araujo et al. (2023) Martini 3 coarse-grained force field for cholesterol. "
               "J. Chem. Theory Comput. 19, 7387-7404. doi:10.1021/acs.jctc.3c00547",
    "openmm": "Eastman et al. (2024) OpenMM 8: Molecular Dynamics Simulation with Machine Learning Potentials. "
              "J. Phys. Chem. B 128, 109-116. doi:10.1021/acs.jpcb.3c06662 (PDBFixer is part of OpenMM)",
    "wimley_white": "Wimley, Creamer & White (1996) Solvation energies of amino acid side chains and backbone in a "
                    "family of host-guest pentapeptides. Biochemistry 35, 5109-5124. doi:10.1021/bi9600153",
    "self_assembly": "Scott et al. (2008) Coarse-grained MD simulations of membrane protein-bilayer self-assembly. "
                     "Structure 16, 621-630. doi:10.1016/j.str.2008.01.014",
    "memprotmd": "Stansfeld et al. (2015) MemProtMD: automated insertion of membrane protein structures into "
                 "explicit lipid membranes. Structure 23, 1350-1361. doi:10.1016/j.str.2015.05.006",
    "gromacs": "Abraham et al. (2015) GROMACS: High performance molecular simulations through multi-level "
               "parallelism from laptops to supercomputers. SoftwareX 1-2, 19-25. doi:10.1016/j.softx.2015.06.001",
}


def for_spec(spec, used: frozenset = frozenset()) -> list[str]:
    keys = []
    ff = spec.cg["forcefield"]
    if ff == "martini3":
        keys.append("martini3")
        if spec.ff_files["martini3_lipids"] == "lipidome":
            keys.append("lipidome")
        if "CHOL" in spec.lipid_names():
            keys.append("sterols")
    else:
        keys += ["martini22", "martini2_lipids"]
        if ff == "elnedyn22":
            keys.append("elnedyn")
    keys.append("martinize2")
    if spec.cg["model"] == "go":
        keys.append("go")
    keys.append("coby")
    if spec.protein["source"] == "opm":
        keys.append("opm")
    if "pdbfixer" in used:
        keys.append("openmm")
    if spec.protein["orientation"] == "hydrophobic":
        keys.append("wimley_white")
    if spec.protein["orientation"] == "self_assembly":
        keys += ["self_assembly", "memprotmd"]
    keys.append("gromacs")
    return [REFS[k] for k in keys]


def markdown(spec, used: frozenset = frozenset()) -> str:
    lines = [f"# References for {spec.name}", "",
             "Please cite the tools and parameters this system was built with",
             "(check each entry against the publisher's page before submission):", ""]
    lines += [f"- {r}" for r in for_spec(spec, used)]
    return "\n".join(lines) + "\n"
