"""Command line: mdbuddy doctor | init | validate | build | stages."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__, doctor, lipids, spec as spec_mod
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


LIB_ALIASES = {
    "m2": "martini2", "martini2": "martini2", "martini22": "martini2", "elnedyn22": "martini2",
    "m3": "martini3_lipidome", "martini3": "martini3_lipidome", "lipidome": "martini3_lipidome",
    "martini3_lipidome": "martini3_lipidome",
    "m3-2021": "martini3_2021", "release2021": "martini3_2021", "martini3_2021": "martini3_2021",
}


def _lib(name: str) -> str:
    try:
        return LIB_ALIASES[name.lower()]
    except KeyError:
        raise SystemExit(f"unknown force field '{name}'; use one of: m2, m3 (lipidome), m3-2021")


def _entry_line(lib: str, name: str) -> str:
    e = lipids.lookup(name, lib)
    if e is None:
        return f"  {lipids.LIBRARY_LABELS[lib]:<44} -"
    charge = f"{e['charge']:+g}" if e.get("charge") is not None else "?"
    coby = e["coby"]["status"]
    return (f"  {lipids.LIBRARY_LABELS[lib]:<44} {name:<6} {e.get('chemistry') or '':<18} "
            f"charge {charge:<4} COBY: {coby}")


def cmd_lipids(args) -> int:
    if args.action == "show":
        for name in args.names:
            print(f"{name}")
            for lib in lipids.LIBRARIES:
                print(_entry_line(lib, name))
            for src in lipids.LIBRARIES:
                if lipids.lookup(name, src) is None:
                    continue
                for dst in lipids.LIBRARIES:
                    if dst != src:
                        eq = lipids.equivalents(name, src, dst)
                        if eq and eq != [name]:
                            print(f"  {src} {name} -> {dst}: {', '.join(eq)}")
            for lib in lipids.LIBRARIES:
                note = (lipids.lookup(name, lib) or {}).get("note")
                if note:
                    print(f"  note ({lib}): {note}")
            print()
        return 0
    if args.action == "convert":
        src, dst = _lib(args.source), _lib(args.target)
        status = 0
        for name in args.names:
            eq = lipids.equivalents(name, src, dst)
            if lipids.lookup(name, src) is None:
                print(f"{name}: not a {lipids.LIBRARY_LABELS[src]} lipid")
                status = 1
            elif not eq:
                print(f"{name} ({lipids.describe(name, src)}) -> no equivalent in {lipids.LIBRARY_LABELS[dst]}")
                status = 1
            else:
                print(f"{lipids.describe(name, src)}")
                q_src = lipids.lookup(name, src).get("charge")
                for n in eq:
                    e = lipids.lookup(n, dst)
                    flag = ""
                    if q_src is not None and e.get("charge") is not None and abs(q_src - e["charge"]) > 1e-6:
                        flag = f"  !! CHARGE CHANGES {q_src:+g} -> {e['charge']:+g}"
                    print(f"   -> {lipids.describe(n, dst)}   [COBY: {e['coby']['status']}]{flag}")
        return status
    if args.action == "report":
        text = lipids.report_markdown()
        if args.names == ["-"]:
            print(text, end="")
        else:
            from pathlib import Path
            Path(args.names[0]).write_text(text)
            print(f"wrote {args.names[0]}")
        return 0
    if args.action == "search":
        lib = _lib(args.ff) if args.ff else None
        hits = lipids.search(" ".join(args.names), lib)
        for lname, name, e in hits[: args.limit]:
            print(f"{lname:<18} {name:<6} {e.get('common') or e.get('chemistry') or ''}")
        if len(hits) > args.limit:
            print(f"... {len(hits) - args.limit} more (use --limit)")
        return 0 if hits else 1
    return 2


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

    p = sub.add_parser("lipids", help="look up lipid names across Martini 2, Martini 3 and the M3 lipidome")
    p.add_argument("action", choices=("show", "convert", "search", "report"),
                   help="show NAME...: meaning in every version; convert NAME... --from m2 --to m3; "
                        "search TEXT [--ff m3]; report FILE|-: Markdown list of every changed/clashing name")
    p.add_argument("names", nargs="*", default=["-"])
    p.add_argument("--from", dest="source", default="m2")
    p.add_argument("--to", dest="target", default="m3")
    p.add_argument("--ff", help="restrict search to one force field (m2, m3, m3-2021)")
    p.add_argument("--limit", type=int, default=40)
    p.set_defaults(func=cmd_lipids)

    p = sub.add_parser("stages", help="list pipeline stages")
    p.set_defaults(func=cmd_stages)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
