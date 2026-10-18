#!/usr/bin/env python3
"""Benchmark `orientation: hydrophobic` against OPM.

Downloads OPM structures (membrane normal = z, centre z = 0), applies a random rotation and
translation, runs the hydrophobic placement, and compares the recovered normal, membrane centre,
thickness and in/out side with OPM. Results are in skills/mdbuddy/references/orientation.md.

    python tools/benchmark_orientation.py [workdir]
"""
import math
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np

from mdbuddy import orient
from mdbuddy.structure import opm_half_thickness, read_pdb

IDS = ["1bl8", "1c3w", "1j4n", "1okc", "1ots", "1py6", "1u19", "2bg9", "2nwl", "2omf", "2rh1", "3emn",
       "3sn6", "4hyj"]
URL = "https://opm-assets.storage.googleapis.com/pdb/{}.pdb"


def random_rotation(rng):
    q = rng.normal(size=4)
    a, b, c, d = q / np.linalg.norm(q)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)], [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])


def main():
    work = Path(sys.argv[1] if len(sys.argv) > 1 else "opm_benchmark")
    work.mkdir(exist_ok=True)
    rng = np.random.default_rng(7)
    print("pdb   residues  tilt_err_deg  centre_err_A  D_found  D_opm  E_kcal  spread_deg  in/out=OPM  verdict")
    for pid in IDS:
        f = work / f"{pid}.pdb"
        if not f.exists():
            urllib.request.urlretrieve(URL.format(pid), f)
        atoms, header = read_pdb(f)
        atoms = [a for a in atoms if a.record == "ATOM" and not a.is_hydrogen and a.resname not in ("DUM", "HOH")]
        R, t = random_rotation(rng), rng.uniform(-20, 20, 3)
        moved = orient.apply(atoms, R, -R.T @ t)  # x' = R x + t
        true_n = R @ np.array([0.0, 0.0, 1.0])
        t0 = time.time()
        _, rec = orient.hydrophobic(moved)
        p = rec["placement"]
        n = np.array(p["normal"])
        tilt = math.degrees(math.acos(min(1.0, abs(float(n @ true_n)))))
        centre = float((np.array(p["center_A"]) - t) @ true_n)
        print(f"{pid}  {len(orient.residues_of(atoms)):8d}  {tilt:12.1f}  {centre:12.1f}  {p['thickness_A']:7.0f}  "
              f"{2 * opm_half_thickness(header):5.1f}  {p['energy_kcal']:6.1f}  {p['uncertainty_deg']:10.1f}  "
              f"{'yes' if n @ true_n > 0 else 'no':>10}  {p['verdict'][:40]}  ({time.time() - t0:.1f} s)")


if __name__ == "__main__":
    main()
