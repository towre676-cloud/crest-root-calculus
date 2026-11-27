"""An exact local and Mordell-Weil sieve certificate for the Gaussian projection."""
from fractions import Fraction
from itertools import product
from math import lcm
import json
from pathlib import Path
from cypari2 import Pari

p=Pari()
p.allocatemem(64_000_000,silent=True)
w=p('w')
field=w**4-w*w+1
z=p.Mod(w,field)
nf=p.nfinit(field)
q0=37


def check(condition,message):
    if not condition:
        raise ArithmeticError(message)


def model(z):
    v=z*z-1
    theta=3*v+2
    eta=theta*(2*v+1)
    coeff=[0,4*eta,0,3*eta**2,0]
    return theta,coeff,[1-v,6-3*v],[3-3*v,-9*v*z**3],[[0],[-eta,0],[-3*eta,0],[0,0]]


theta,coeff,G1,G2,T=model(z)
E=p.ellinit(coeff,nf)


def roots(q):
    return sorted((t for t in range(q) if (t**4-t*t+1)%q==0),
                  key=lambda t:(t**3%q,(2*t*t-1)%q))


def reduction_rows(q):
    out=[]
    for t in roots(q):
        th,co,g1,g2,tors=model(t)
        e=p.ellinit(co,q)
        check(int(e.disc())!=0 and e.ellisoncurve(g1) and e.ellisoncurve(g2),'Bad reduction data')
        out.append((t,th%q,e,g1,g2,tors))
    check(len(out)==4,'Prime must split completely')
    return out


def abscissa(pt,th,q):
    return None if len(pt)==1 else int(p.lift(pt[0]))*pow(th,-1,q)%q


def residue_table():
    q=q0
    rows=reduction_rows(q)
    out=[]
    for a,b,t in product(range(9),range(9),range(4)):
        vals=[]
        ys=[]
        for _,th,e,g1,g2,tors in rows:
            P=e.elladd(e.elladd(e.ellmul(g1,a),e.ellmul(g2,b)),tors[t])
            vals.append(abscissa(P,th,q))
            ys.append(None if len(P)==1 else int(p.lift(P[1])))
        if vals[0]==vals[1] and vals[2]==vals[3]:
            check(vals[0]==vals[2],'Unexpected nonrational residue projection')
            out.append(dict(a=a,b=b,t=t,u=vals[0],ordinates=ys))
    expected={(0,0,0,None),(0,0,3,0)}
    for a,b in ((0,1),(0,8),(2,0),(7,0)):
        expected.update(((a,b,1,28),(a,b,2,1)))
    check({(r['a'],r['b'],r['t'],r['u']) for r in out}==expected,'Full residue class enumeration')
    return out


def at(value,t,mod):
    f=p.lift(value)
    ans=0
    for j in range(4):
        c=Fraction(str(p.polcoef(f,j,w)))
        ans+=c.numerator*pow(c.denominator,-1,mod)*pow(t,j,mod)
    return ans%mod


def formal_data(q):
    lifts=[]
    for a in roots(q):
        candidates=[b for b in range(a,q*q,q) if (b**4-b*b+1)%(q*q)==0]
        check(len(candidates)==1,'Simple Hensel lift')
        lifts.append(candidates[0])
    def divide(q2value):
        check(q2value%q==0,'Expected formal point')
        return q2value//q
    basis=[]
    for G in (G1,G2):
        P=E.ellmul(G,9)
        t=-P[0]/P[1]
        basis.append([divide(at(t,h,q*q)) for h in lifts])
    check(basis==[[7,0,7,0],[24,15,13,22]],'Leading formal-log vectors')
    return lifts,basis


def tangent_tests(rows):
    q=q0
    lifts,basis=formal_data(q)
    ths=[at(theta,t,q) for t in lifts]
    check(ths==[6,32,6,32],'Theta embeddings')
    determinants=[]
    for r in rows:
        if r['u']!=1:
            continue
        slopes=[2*y*pow(th,-1,q)%q for y,th in zip(r['ordinates'],ths)]
        matrix=[[(slopes[j]*v[j]-slopes[j+1]*v[j+1])%q for v in basis] for j in (0,2)]
        det=(matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0])%q
        check(det==(13 if r['a'] else 24),'Unit-disc determinant')
        # Independently evaluate the actual finite differences at each global coset representative.
        center=E.elladd(E.elladd(E.ellmul(G1,r['a']),E.ellmul(G2,r['b'])),T[r['t']])
        for col,G in enumerate((G1,G2)):
            neighbor=E.elladd(center,E.ellmul(G,9))
            du=(neighbor[0]-center[0])/theta
            vals=[at(du,h,q*q) for h in lifts]
            check(all(c%q==0 for c in vals),'Difference not in the residue disk')
            differences=[(vals[j]//q-vals[j+1]//q)%q for j in (0,2)]
            check(differences==[line[col] for line in matrix],'Exact translation derivative')
        determinants.append(dict(a=r['a'],b=r['b'],jacobian=matrix,determinant=det))
    directions=[(i,1) for i in range(q)]+[(1,0)]
    residuals=[]
    for a,b in directions:
        t=[(a*u+b*v)%q for u,v in zip(*basis)]
        residuals.append([(ths[j]*t[j]**2-ths[j+1]*t[j+1]**2)%q for j in (0,2)])
    check([0,0] not in residuals,'Nonzero tangent direction at infinity')
    return dict(roots=roots(q),roots_mod_1369=lifts,theta=ths,formal_log_columns=basis,
                unit_disks=determinants,projective_directions_at_infinity=len(residuals),
                common_zeros_at_infinity=0)


def index_tests():
    # If a combination is ell-divisible globally, [exponent/ell] annihilates
    # its reduction in every component. Test all directions in P^1(F_ell).
    out=[]
    for prime,ell in ((37,3),(109,31)):
        rows=reduction_rows(prime)
        for _,_,e,g1,g2,_ in rows:
            check(list(map(int,e.ellgroup()))==([18,2] if prime==37 else [62,2]),'Finite group')
            check([int(e.ellorder(G)) for G in (g1,g2)]==([9,9] if prime==37 else [62,62]),'Reduced orders')
        count=0
        for a,b in [(a,1) for a in range(ell)]+[(1,0)]:
            detected=False
            for _,_,e,g1,g2,_ in rows:
                P=e.elladd(e.ellmul(g1,a),e.ellmul(g2,b))
                if len(e.ellmul(P,int(e.ellgroup()[0])//ell))!=1:
                    detected=True
                    break
            check(detected,'Odd saturation direction not detected')
            count+=1
        out.append(dict(prime=prime,saturation_prime=ell,excluded_directions=count))
    return out


def fixed_subgroup_sieve(q):
    rows=reduction_rows(q)[:2]
    period=lcm(*(int(e.ellorder(g1)) for _,_,e,g1,_,_ in rows))
    good=[]
    for a in range(period):
        u=[abscissa(e.elladd(e.ellmul(g1,a),tors[2]),th,q) for _,th,e,g1,_,tors in rows]
        if u[0]==u[1]:
            good.append(a)
    check((period,good)==((24,[4,20]) if q==13 else (62,[31])),'Fixed subgroup sieve')
    return dict(prime=q,period=period,coefficients=good)


def verify():
    check(int(E.elltors()[0])==4 and all(E.ellisoncurve(g) and E.ellorder(g)==0 for g in (G1,G2)),'Global torsion or generators')
    check(int(p.bnfcertify(p.bnfinit(field,1)))==1 and int(p.nfdisc(field))==144,'Quartic field certification')
    for g in (G2,E.ellneg(G2)):
        check(E.elladd(g,T[2])[0]==theta,'Known ordinary point')
        check(E.elladd(g,T[1])[0]==-9*theta,'Known reciprocal point')
    r=residue_table()
    return dict(odd_saturation=index_tests(),residues=r,formal=tangent_tests(r),
                fixed_subgroup=[fixed_subgroup_sieve(q) for q in (13,109)],
                projection_theta1=['infinity',0,1,-9],projection_theta2=['infinity',0,-1,9])


if __name__=='__main__':
    print(json.dumps(verify(),indent=2))
