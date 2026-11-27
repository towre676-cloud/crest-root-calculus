import sys
from descent import Field, selmer, pari
K = Field('v^2+v+1')                     # K = Q(zeta3), v = zeta3; sqrt(-3) = 2v+1
r3 = K.M('2*v+1'); assert r3**2 == -3
zeta6 = K.M('v+1')
for name, theta in (('theta1', 3*zeta6 - 1), ('theta2', 3*zeta6 - 2)):
    for sign in (1, -1):
        eta = sign*theta*r3
        es = (K.M(0), -eta, -3*eta)
        primes = [P for q in (2, 3, 7) for P in pari.idealprimedec(K.bnf, q)
                  if q != 7 or pari.idealval(K.bnf, pari.lift(eta), P) > 0]
        sel, basis = selmer(K, es, primes)
        import math
        print(name, 'sign', sign, 'eta', pari.lift(eta), 'norm', pari.norm(eta), '#Sel', len(sel), 'rank <=', int(math.log2(len(sel))) - 2, flush=True)
