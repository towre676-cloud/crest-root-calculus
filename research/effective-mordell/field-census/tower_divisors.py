"""Reduced divisors on the twelfth cyclotomic base for the next Gaussian tower."""
from collections import Counter, deque
import json
from pathlib import Path
from math import gcd
from gaussian_tower_fields import pari, x, w, i_base, base, gaussian_coordinates, gaussian_equation, gaussian_test, require

poly=w**4-w**2+1
nf=pari.nfinit(poly)
z=pari.Mod(w,poly)
i=z**3
one=i**0
f=((x+1)**6-x**6)*one

def degree(p):return int(p.poldegree(x)) if p else -1

def monic(p):return p/pari.pollead(p,x)

def add(D,H):
    u,v=D
    s,t=H
    a,b,d=pari.gcdext(u,s)
    c,h,e=pari.gcdext(d,v+t)
    U=monic(u*s/e**2)
    V=((c*(a*u*t+b*s*v)+h*(v*t+f))/e)%U
    while degree(U)>2:
        U=monic((f-V**2)/U)
        V=(-V)%U
    require((f-V**2)%U==0,'Cantor addition invalid')
    return U,V

def key(D):return tuple(tuple(str(pari.lift(pari.polcoef(k,j,x))) for j in range(degree(k)+1)) for k in D)

def group():
    gen=[(monic(x-1/(a-1)),0*one) for a in pari.nfroots(nf,x**6-1) if a!=1]
    gen.extend([(x*one,one),(x+one,i)])
    identity=(x**0*one,0*one)
    seen={key(identity):identity};q=deque([identity])
    while q:
        D=q.popleft()
        for g in gen:
            H=add(D,g);label=key(H)
            if label not in seen:
                q.append(H);seen[label]=H
                require(len(seen)<=144,'exceeded 144 divisor classes')
    require(len(seen)==144,'divisor group incomplete')
    return list(seen.values())

def conj(a):
    # The nonidentity automorphism fixing i=z^3 sends z to z^5.
    return pari.subst(pari.lift(a),w,z**5)

def abelian_rep(disc):
    nd=gaussian_coordinates(disc*conj(disc),i,4)
    roots=pari.nfroots(base,x*x-nd)
    if not len(roots):return None
    for s in roots:
        trace=gaussian_coordinates(disc+conj(disc),i,4)
        rep=trace+2*s
        if rep:
            # disc = (disc+s)^2 / rep; hence it is the class of rep over F0.
            return rep
    raise ArithmeticError('vanishing trace identity')

def in_all_tower(delta):
    primes,exps=pari.idealfactor(base,delta)
    gen=i_base**0
    for pr,n in zip(primes,exps):
        if not int(n)%2:continue
        norm=int(pari.idealnorm(base,pr))
        if norm==9:
            gen*=3
            continue
        if norm%3!=2 or norm==2:return False
        _,v=pari.bnfisprincipal(base,pr)
        h=pari.Mod(pari.nfbasistoalg(base,v),gaussian_equation)
        if not int(pari.polcoef(pari.lift(h),0,pari('g')))%2:h*=i_base
        gen*=h
    return bool(len(pari.nfroots(base,x*x-delta/gen)))

def verify():
    sizes=[int(pari.ellinit([0,-1],q).ellcard()*pari.ellinit([0,1],q).ellcard()) for q in (13,37)]
    require(sizes==[144,2304] and gcd(*sizes)==144, 'cyclotomic field good-reduction bounds')
    delta=1+2*i_base
    require((1-i_base)**3-delta**3==9 and not len(pari.nfroots(base,x*x-delta)),
            'nontorsion point on the parity-odd twist (only identity and rational 2-torsion)')
    ds=group()
    counts=Counter(); reps={}; fibre_counts=Counter(); receipt=[]
    h=2*z-z**3
    require(h*h==3,'square-root-of-three basis')
    fields=[('2',2),('7',7),('11',11),('13+8sqrt3',13+8*h),('13-8sqrt3',13-8*h)]
    def square_in_omega(d):
        if len(pari.nfroots(nf,x*x-d)):return True
        r=abelian_rep(d)
        return r is not None and gaussian_test(r)[1]
    for j,(name,a) in enumerate(fields):
        require(not square_in_omega(a*one),'extension is trivial')
        for _,b in fields[:j]:
            require(not square_in_omega((a*one)/b),'exceptional extensions coincide')
    for u,v in ds:
        counts['all']+=1
        if degree(u)!=2:continue
        D=pari.poldisc(u,x)
        if len(pari.nfroots(nf,x*x-D)):continue
        counts['irreducible quadratic supports']+=1
        possible=[name for name,radicand in fields if square_in_omega(D/radicand)]
        require(len(possible)==1,'unclassified quadratic field of exceptional support')
        fibre_counts[possible[0]]+=2
        receipt.append(dict(radicand=possible[0],u=key((u,))[0],v=key((v,))[0]))
        d=abelian_rep(D)
        if d is None:continue
        counts['square norm']+=1
        if not in_all_tower(d):continue
        counts['splits in full Gaussian tower']+=1
        reps[key((u,v))]=str(d)
    require(dict(counts)=={'all':144,'irreducible quadratic supports':96,'square norm':72},'quadratic divisor totals')
    require(not reps,'unexpected extra supports in the full Gaussian radical tower')
    require(dict(fibre_counts)=={'2':24,'7':72,'11':48,'13+8sqrt3':24,'13-8sqrt3':24},
            'full quadratic-extension counts')
    receipt.sort(key=lambda row:(row['radicand'],row['u'],row['v']))
    expected=Path(__file__).with_name('exceptional-divisors.json')
    if expected.exists():
        require(json.loads(expected.read_text())==json.loads(json.dumps(receipt)),'recorded divisors disagree')
    else:
        expected.write_text(json.dumps(receipt,indent=2)+'\n')
    print('Quadratic-extension point counts outside the base tower:',dict(fibre_counts))
    print(dict(counts))
    print('Extra supports',reps)
    return counts
if __name__=='__main__':verify()
