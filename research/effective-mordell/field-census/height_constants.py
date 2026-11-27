"""Rigorous real-ball checks for the two uniform abscissa constants."""
from flint import arb,ctx

def verify():
    ctx.prec=160
    A=arb(132).log()/3
    B=arb(9).log()/3
    C=arb(36).log()
    for c,decimal in ((arb(2),"4.278300842"),(arb(3).sqrt(),"4.964624902")):
        k=(A+C+(2-c)*B)/c
        t=(A-B+((A+B)**2+4*((B+C)/c)**2).sqrt())/2
        if not (t>k and t>A and (t/2).exp()<arb(decimal)):
            raise ArithmeticError("Upper bound interval failed")
        print((t/2).exp())

if __name__=="__main__":
    verify()
