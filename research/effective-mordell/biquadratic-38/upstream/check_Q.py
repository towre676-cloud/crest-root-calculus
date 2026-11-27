import math
from descent import Field, selmer, pari
Q = Field('t')  # degree-1 field
bad = 0
for d in [1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 21, 30, -1, -2, -3, -5, -7, -15]:
    es = (Q.M(0), Q.M(-d), Q.M(-3*d))
    primes = [P for q in pari.factor(abs(6*d))[0] for P in pari.idealprimedec(Q.bnf, q)]
    sel, _ = selmer(Q, es, primes)
    r = pari.ellrank(pari.ellinit([0, 4*d, 0, 3*d*d, 0]))
    ub = int(math.log2(len(sel))) - 2
    flag = '' if r[1] <= ub and r[0] <= ub else ' MISMATCH'
    if flag: bad += 1
    print(d, 'selmer rank', ub, 'ellrank', list(r[:2]), flag, flush=True)
print('mismatches', bad)
