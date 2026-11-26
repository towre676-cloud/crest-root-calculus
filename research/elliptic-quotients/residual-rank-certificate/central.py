"""Rigorous central orders for the elliptic constituents of the residual kernels."""
import argparse
import json
import sys
from pathlib import Path

from cypari2 import Pari
from flint import acb, acb_poly, arb, ctx


HERE = Path(__file__).resolve().parent
PARI = Pari()
TERMS = 1024
ENDPOINT = 2048
DIGITS = 12

# On the minimal twisted models, two independent rational points for PARI's descent.
TWIST_POINTS = {
    "338c1": [[111,1312],[303,5408]],
    "676b1": [[-5,8],[3,8]],
    "867e1": [[-1,4],[17,76]],
}


def rounded_enclosure(value):
    unit = 10**DIGITS
    scaled = value*unit
    lower = int(scaled.lower().floor().unique_fmpz())
    upper = int(scaled.upper().ceil().unique_fmpz())

    def decimal(n):
        sign = "-" if n<0 else ""
        whole,fraction = divmod(abs(n),unit)
        return f"{sign}{whole}.{fraction:0{DIGITS}d}"

    return [decimal(lower),decimal(upper)]


def mellin_interval(curve, order):
    """Enclose the derivative of the completed L-function for sign +1."""
    N = int(curve.ellglobalred()[0])
    if curve.ellrootno()!=1 or order not in (0,2):
        raise ValueError("The verifier uses even order zero or two and sign +1.")
    alpha = 2*arb.pi()/arb(N).sqrt()
    polynomial = acb_poly([0]+list(map(int,curve.ellan(TERMS))))

    def integrand(y, analytic):
        value = 2*polynomial((-alpha*y).exp())
        return value if order==0 else value*y.log(analytic=analytic)**2

    finite = acb(0)
    left = 1
    while left<ENDPOINT:
        right = min(2*left,ENDPOINT)
        finite += acb.integral(integrand,left,right,abs_tol=arb("1e-25"),rel_tol=arb("1e-25"))
        left = right
    if not finite.imag.contains(0):
        raise ArithmeticError("The integral along the real interval must be real.")

    # Uniformly for y>=1, |a_n|<=2n and (log y)^2<=y^2.
    R = TERMS+1
    exponential = (-alpha).exp()
    series_tail = 4*(-alpha*R).exp()/(1-exponential)*(
        1/alpha+2/(alpha**2*R)+2/(alpha**3*R*R))
    Y = ENDPOINT
    exp_y = (-alpha*Y).exp()
    integral_tail = 4*exp_y/(1-exp_y)**2*(Y*Y/alpha+2*Y/alpha**2+2/alpha**3)
    total_tail = (series_tail+integral_tail).upper()
    value = alpha*(finite.real+arb(0,total_tail))
    return value, total_tail


def certify(packet, check_symbols):
    label = packet["label"]
    E = PARI.ellinit(packet["ainvs"])
    twist = PARI.ellminimalmodel(PARI.ellinit(E.elltwist(packet["twist"])))
    N = int(twist.ellglobalred()[0])
    if N!=packet["twisted_conductor"] or twist.ellrootno()!=1 or packet["central_L_over_real_period"]!="0":
        raise ArithmeticError("The twist differs from its exact packet certificate.")
    if check_symbols:
        module,symbol = twist.msfromell(1)
        if PARI.mseval(module,symbol,[PARI("oo"),0])!=0:
            raise ArithmeticError("The exact central modular symbol is nonzero.")
    derivative,error = mellin_interval(twist,2)
    untwisted,_ = mellin_interval(E,0)
    if not (derivative>0 and untwisted>0 and error<arb("1e-20")):
        raise ArithmeticError("The numerical enclosures must exclude zero.")
    points = TWIST_POINTS[label]
    if not all(twist.ellisoncurve(point) for point in points):
        raise ArithmeticError("The descent points are not on the minimal twisted model.")
    rank = twist.ellrank(0,points)
    if list(map(int,rank[:2]))!=[2,2] or E.ellrank(0)[:2]!=PARI([0,0]):
        raise ArithmeticError("The independent Mordell-Weil bounds did not coincide.")
    return {
        "curve":label,"D":packet["twist"],"twisted_conductor":N,"twisted_root_number":1,
        "minimal_twist_ainvs":list(map(int,twist[:5])),"twisted_a3":int(twist.ellap(3)),
        "twisted_L_second_derivative":rounded_enclosure(derivative),
        "untwisted_L_at_one":rounded_enclosure(untwisted),
        "completed_derivative_tail_below_1e_20":True,
        "twisted_analytic_rank":2,"untwisted_analytic_rank":0,
        "twisted_mordell_weil_rank":2,"untwisted_mordell_weil_rank":0,
        "descent_points":points,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-symbols",action="store_true",
                        help="Rerun the central modular-symbol calculation (large PARI stack).")
    parser.add_argument("--verify",type=Path)
    args = parser.parse_args()
    ctx.prec = 128
    if args.check_symbols:
        PARI.allocatemem(128_000_000,2_500_000_000,silent=True)
    source = HERE.parent/"exact-residual-spectrum"/"packets.json"
    levels = json.loads(source.read_text())
    results = [certify(p,args.check_symbols) for level in levels for p in level["packets"]]
    if args.verify and results!=json.loads(args.verify.read_text()):
        raise ArithmeticError("The central-order replay differs from the saved receipt.")
    print(json.dumps(results,indent=2))


if __name__ == "__main__":
    main()
