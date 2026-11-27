"""Replay the Python calculations. Add --magma to run the two independent Magma inputs."""
import argparse
import json
from pathlib import Path
import subprocess
import sys


ROOT=Path(__file__).resolve().parent
CHECKS=(
    ("verify_new_descent.py",),("verify_panel.py",),("verify_extended_quadratic.py",),
    ("prime_density.py",),("height_constants.py",),("rank_one_height.py","--panel"),
    ("tower_replay.py",),("gaussian_tower_fields.py",),("tower_divisors.py",),
    ("non_galois_rank_one.py",),("verify_classgroup_descent.py","--workers","3","--quadratic-check"),
    ("octic_verify.py",),("cover_transfer.py","verify"),("discriminant_obstruction.py",),
)


def replay(magma=None):
    if sys.flags.optimize:
        raise RuntimeError("Replay with normal Python assertions enabled.")
    for args in CHECKS:
        print("Running",*args,flush=True)
        subprocess.run([sys.executable,*args],cwd=ROOT,check=True)
    for script,archive in (("quartic_bridge.py","quartic-height.json"),
                           ("zeta24.py","zeta24-points.json"),("rank_growth.py","rank-growth.json")):
        value=subprocess.check_output([sys.executable,script],cwd=ROOT)
        if json.loads(value)!=json.loads((ROOT/archive).read_text()):
            raise ArithmeticError("Archived certificate disagreement in "+archive)
    from solve_field import classify,p,w
    for d,status,total in ((2,'complete',46),(7,'partial',130),(11,'partial',58),(17,'complete',10)):
        f=p.polcompositum(p.polcyclo(12,w),w*w-d)[0]
        answer=classify(list(map(int,f.Vecrev())))
        if (answer['status'],2+len(answer['verified_affine_points']))!=(status,total):
            raise ArithmeticError("Octic solver model transport failed")
    from cm_height import certificate,exhaustive_points
    rational=exhaustive_points([0,1],certificate([0,1]))
    if rational['status']!='complete' or set(rational['verified_affine_points'])!={('-1','0'),('1','0')}:
        raise ArithmeticError("Rational height-box enumeration failed")
    print("Python proofs and archived lists passed.",flush=True)
    if magma:
        for script in ("sqrt7.m","final_three.m"):
            subprocess.run([magma,"-b",script],cwd=ROOT,check=True)
    else:
        print("Magma proof inputs are supplied separately. Pass --magma /path/to/magma for a fresh Chabauty run.")


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--magma",metavar="EXECUTABLE")
    replay(parser.parse_args().magma)
