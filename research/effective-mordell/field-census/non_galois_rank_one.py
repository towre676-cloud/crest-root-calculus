"""A positive-rank non-Galois quartic test, certified from its defining polynomial."""
import json
from rank_one_general import classify, points_low_degree, relative_radicand, contains_gaussian_square, small_gaussians
from rank_one_height import cube, gaussian_mul, pari, check


def verify():
    co = [2,0,-2,0,1]
    result = classify(co)
    check(result["status"]=="complete" and result["rank_upper_bound"]==1,"Rank bound")
    X=pari("x")
    f=pari.Polrev(co,pari("w"))
    nf=pari.nfinit(f)
    w=pari.Mod(pari("w"),f)
    i=w*w-1
    check(i*i==-1,"Field generator")
    delta=relative_radicand(nf,i,X)
    curve=pari.ellinit([0,-1],nf)
    torsion=int(curve.elltors()[0])
    witness=None
    for a,b,c in small_gaussians():
        u,v=cube((a,b))
        u-=c**3
        scaled=gaussian_mul((u*c**3,v*c**3),delta)
        if not contains_gaussian_square(scaled):
            continue
        xx=(a+b*i)/c
        for yy in pari.nfroots(nf,X**2-(xx**3-1)):
            P=pari([xx,yy])
            if curve.ellmul(P,torsion)!=pari([0]):
                witness=list(map(str,[pari.lift(xx),pari.lift(yy)]))
                break
        if witness:break
    check(witness is not None,"Missing nontorsion witness")
    check(len(result["verified_affine_points"])+2==6,"Number of C6 points")
    print(json.dumps(dict(curve_points=6,elliptic_rank=1,elliptic_torsion_order=torsion,
                          infinite_order_witness=witness),indent=2))

if __name__=="__main__":
    verify()
