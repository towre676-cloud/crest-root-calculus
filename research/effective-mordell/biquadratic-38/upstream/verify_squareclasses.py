"""Exhaust the 64 squareclasses in each coordinate and certify the two final pairs."""
from sieve import *
assert int(pari.bnfcertify(L))==1 and int(L.bnf_get_no())==1
assert len(S)==4 and len(CLASSES)==64
pairs=global_filter([(d,g) for d in CLASSES for g in CLASSES])
assert len(pairs)==256
for q in (5,7,11,13,17,19,23,29,31,37,41,43):
    pairs=residue_test(pairs,q)
assert len(pairs)==48
for q,m in ((3,9),(2,8),(2,16)):
    pairs=ramified_test(pairs,q,m)
a,b=ZETA**3+2*ZETA,2*ZETA**3+ZETA
assert len(pairs)==2 and all(any(d==ZETA and g==gamma for d,g in pairs) for gamma in (a,b))
assert is_square(-tau(ZETA)/ZETA) and is_square(-tau(a)/b) and is_square(-tau(b)/a)
assert ZETA*a==3*ZETA**2-1 and ZETA*b==3*ZETA**2-2
assert tau(ZETA*a)==-ZETA*b and tau(R3)==-R3
print('4,096 pairs to 256 to 48 to 2. Symmetries and products verified.')
