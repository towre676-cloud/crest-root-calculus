"""Complete fixed-degree scans with corrected/stock SymPy, Singular and GP."""

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import sympy as sp
from sympy.polys.densetools import _dup_left_decompose, _dup_right_decompose

from discover import X, coefficients, constant_degrees, polynomial

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "solver-comparison"))
from map_worker import corrected_right_decomposition


def dense(source, *, corrected=True):
    monic = source.monic().rep.to_list()
    finder = corrected_right_decomposition if corrected else _dup_right_decompose
    maps = []
    for d, m in constant_degrees(source):
        right = finder(monic, d, sp.QQ)
        left = _dup_left_decompose(monic, right, sp.QQ)
        if left is not None:
            u = sp.Poly.from_list(right, X, domain=sp.QQ)
            g = sp.Poly.from_list(left, X, domain=sp.QQ).mul_ground(source.LC())
            maps.append({"degree": d, "target_degree": m,
                         "inner": coefficients(u), "outer": coefficients(g)})
    return maps


def singular(source):
    expression = "+".join(f"({c})*x^{d}" for (d,), c in source.terms())
    library = os.environ.get("CREST_DECOMP_LIB", "decomp.lib")
    script = f'''ring r=0,x,dp;
LIB "{library}";
int MINS=1;
int DECMETH=0;
poly p={expression}; number unit=leadcoef(p); poly f=p/unit;
short=0; int d,m; int n=deg(f); ideal result;
for(d=2;3*d<=n;d++){{
    if(n%d==0){{
        m=n div d;
        result=decompunivmonic(f,m);
        if(deg(result[1])==m && deg(result[2])==d){{
            print("MAP "+string(d)+" "+string(m));
            print("INNER "+string(result[2]));
            print("OUTER "+string(unit*result[1]));
        }}
    }}
}}
print("DONE");
quit;
'''
    executable = os.environ.get("CREST_QUOTIENT_SINGULAR", "Singular")
    with tempfile.TemporaryDirectory(prefix="crest-quotient-singular-") as directory:
        path = Path(directory) / "decompose.sing"
        path.write_text(script, encoding="ascii")
        run = subprocess.run([executable, "-q", str(path)], capture_output=True, text=True, check=True)
    if "?" in run.stdout or "much slower" in run.stdout or "DONE" not in run.stdout:
        raise ArithmeticError(run.stdout[-2600:] + run.stderr[-2000:])
    lines = [line for line in run.stdout.splitlines() if line.startswith(("MAP ", "INNER ", "OUTER "))]
    maps = []
    for header, inner, outer in zip(lines[::3], lines[1::3], lines[2::3], strict=True):
        d, m = map(int, header.split()[1:])
        g = sp.Poly(outer[6:].replace("^", "**"), X, domain=sp.QQ)
        u = sp.Poly(inner[6:].replace("^", "**"), X, domain=sp.QQ)
        # The library offers zero-constant normalization. This also makes the
        # wrapper safe if a compatible version chooses a different translation.
        g = g.shift(u.nth(0))
        u -= u.nth(0)
        maps.append({"degree": d, "target_degree": m,
                     "inner": coefficients(u), "outer": coefficients(g)})
    return maps


def gp(values):
    """A native GP implementation with its own exact squarefree normalization."""
    scripts = HERE / "decompose.gp"
    value_list = ",".join(coefficients(polynomial(values)))
    request = scripts.read_text(encoding="ascii") + f"\nscan(Polrev([{value_list}]));\nquit\n"
    executable = os.environ.get("CREST_QUOTIENT_GP", "gp")
    run = subprocess.run([executable, "-q", "-f", "-s", "128M"], input=request,
                         capture_output=True, text=True, check=True)
    if "***" in run.stderr or "DONE" not in run.stdout:
        raise ArithmeticError(run.stdout[-2600:] + run.stderr[-2000:])
    lines = dict(line.split(" ", 1) for line in run.stdout.splitlines()
                 if line.startswith(("SQUARE ", "SOURCE ", "PHASES ")))

    def parse(value):
        return [str(sp.Rational(part.strip())) for part in value[1:-1].split(",")]

    maps = []
    for line in run.stdout.splitlines():
        if line.startswith("MAP "):
            degrees, inner, outer = line[4:].split("|")
            d, m = map(int, degrees.split())
            maps.append({"degree": d, "target_degree": m, "inner": parse(inner), "outer": parse(outer)})
    milliseconds = list(map(int, lines["PHASES"][1:-1].split(",")))
    return parse(lines["SQUARE"]), parse(lines["SOURCE"]), maps, [value / 1000 for value in milliseconds]
