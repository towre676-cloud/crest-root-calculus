import math, itertools
from descent import Field, selmer, pari
K = Field('v^2+v+1'); r3 = K.M('2*v+1'); zeta6 = K.M('v+1'); v = K.M('v')

def search(es, H=12, den=6):
    """Points with x = a/b, a in O_K small, b in small integers; returns non-torsion candidates."""
    E = pari.ellinit([0, -(es[0]+es[1]+es[2]), 0, es[0]*es[1]+es[0]*es[2]+es[1]*es[2], -es[0]*es[1]*es[2]], K.bnf)
    found = []
    for b in range(1, den+1):
        for a0 in range(-H, H+1):
            for a1 in range(-H, H+1):
                x = (a0 + a1*v)/b**2
                fx = (x-es[0])*(x-es[1])*(x-es[2])
                if fx == 0: continue
                r = pari.nfroots(K.bnf, pari('zz^2') - pari.lift(fx))
                if len(r):
                    P = [pari.lift(x), r[0]]
                    if pari.ellorder(E, P) == 0:
                        return E, P
    return E, None

cases = {'E theta1': (3*zeta6-1, 1), 'E theta2': (3*zeta6-2, 1),
         "E' theta'": (zeta6-3, -1)}
for name, (theta, sgn3) in cases.items():
    for sign in (1, -1):
        eta = sign*theta*r3
        es = (K.M(0), -eta, sgn3*(-3)*eta)   # y^2 = x(x+eta)(x+3eta) or x(x+eta)(x-3eta)
        primes = [P for q in (2, 3, 7) for P in pari.idealprimedec(K.bnf, q)
                  if q != 7 or pari.idealval(K.bnf, pari.lift(eta), P) > 0]
        sel, _ = selmer(K, es, primes)
        E, P = search(es)
        print(name, 'sign', sign, 'rank <=', int(math.log2(len(sel)))-2, 'point', P, 'tors', pari.elltors(E)[1], flush=True)
