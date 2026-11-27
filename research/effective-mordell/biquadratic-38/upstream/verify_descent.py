"""Certify the two field presentations, Selmer ranks, and the mod-two index."""
from descent import Field,pari,selmer
K=Field('v^2+v+1')
assert int(pari.bnfcertify(K.bnf))==1 and int(K.bnf.bnf_get_no())==1
v=K.M('v')
r=2*v+1
theta=3*v+2
for sign,point in [(1,[1-v,6-3*v]),(-1,[3*v-3,-9*v])]:
    eta=sign*r*theta
    es=(K.M(0),-eta,-3*eta)
    primes=[pr for q in (2,3,7) for pr in pari.idealprimedec(K.bnf,q)
            if q!=7 or pari.idealval(K.bnf,pari.lift(eta),pr)>0]
    group,basis=selmer(K,es,primes)
    e=pari.ellinit([0,4*eta,0,3*eta**2,0],K.bnf)
    assert len(group)==8 and e.ellisoncurve(point) and e.ellorder(point)==0
print('Both rank-one character groups certified.')
from mw_search import setup,two_saturated,L
assert int(pari.bnfcertify(L))==1 and int(L.bnf_get_no())==1
theta,eta,e,es,g1,g2,ts=setup('theta1')
assert e.elltors()[0]==4 and two_saturated(e,es,[g1,g2,*ts])==(True,None)
print('Rank two, full two-torsion, and 2-saturation certified.')
