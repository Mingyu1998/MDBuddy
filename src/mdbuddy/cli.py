"""Command line: mdbuddy doctor | init | validate | build | stages."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__, doctor, spec as spec_mod
from .pipeline import STAGES, BuildError, build


def cmd_doctor(args) -> int:
    checks = doctor.run_checks(network=args.network)
    print(doctor.report(checks, as_json=args.json))
    return 0 if doctor.ready(checks) else 1


def cmd_init(args) -> int:
    path = Path(args.output or f"{args.name}.yaml")
    if path.exists() and not args.force:
        print(f"{path} exists; use --force to overwrite", file=sys.stderr)
        return 1
    path.write_text(spec_mod.TEMPLATE.format(name=args.name, pdb_id=args.pdb_id or "1bl8"))
    print(f"wrote {path}")
    return 0


def cmd_validate(args) -> int:
    try:
        s = spec_mod.load(args.spec)
    except spec_mod.SpecError as exc:
        if args.json:
            print(json.dumps({"valid": False, "problems": exc.problems}, indent=2))
        else:
            print(exc, file=sys.stderr)
        return 1
    if args.json:
        print(json.dumps({"valid": True, "resolved": s.data}, indent=2))
    else:
        print(s.to_yaml())
        print(f"# spec is valid; output directory: {s.output_dir()}")
    return 0


def cmd_build(args) -> int:
    try:
        s = spec_mod.load(args.spec)
    except spec_mod.SpecError as exc:
        print(exc, file=sys.stderr)
        return 1
    if args.from_stage == "fetch" and not args.skip_doctor:
        checks = doctor.run_checks()
        if not doctor.ready(checks):
            print(doctor.report(checks), file=sys.stderr)
            return 1
    try:
        result = build(s, from_stage=args.from_stage, until_stage=args.until,
                       grompp_check=not args.no_grompp, force=args.force)
    except (BuildError, RuntimeError) as exc:
        print(f"\nBUILD FAILED: {exc}", file=sys.stderr)
        return 2
    out = result["output_dir"]
    print(f"\nDone. GROMACS inputs: {out / 'gromacs'}\nProvenance: {result['lock']}")
    return 0


def cmd_stages(args) -> int:
    for s in STAGES:
        print(s)
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="mdbuddy", description="Reproducible Martini membrane-protein "
                                     "system builder for GROMACS")
    parser.add_argument("--version", action="version", version=f"mdbuddy {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("doctor", help="check that required programs and packages are installed")
    p.add_argument("--json", action="store_true", help="machine-readable output")
    p.add_argument("--network", action="store_true", help="also test access to OPM, RCSB and GitHub")
    p.set_defaults(func=cmd_doctor)

    p = sub.add_parser("init", help="write a commented spec template")
    p.add_argument("name")
    p.add_argument("--pdb-id")
    p.add_argument("-o", "--output")
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_init)

    p = sub.add_parser("validate", help="check a spec and print it with all defaults filled in")
    p.add_argument("spec")
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser("build", help="build the system described by a spec")
    p.add_argument("spec")
    p.add_argument("--from", dest="from_stage", default="fetch", choices=STAGES,
                   help="resume from this stage (earlier results are reused)")
    p.add_argument("--until", default="check", choices=STAGES, help="stop after this stage")
    p.add_argument("--no-grompp", action="store_true", help="skip the final gmx grompp check")
    p.add_argument("--skip-doctor", action="store_true", help="do not check the environment first")
    p.add_argument("--force", action="store_true", help="overwrite an existing output directory")
    p.set_defaults(func=cmd_build)

    p = sub.add_parser("stages", help="list pipeline stages")
    p.set_defaults(func=cmd_stages)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
