"""Count geometric genus-two Jacobian isogeny targets for a source with geometric End=Z."""
import argparse
import json

import sympy as sp


def prime_counts(ell, r):
    if r==0:
        return 1,1
    lines=(ell+1)*(ell*ell+1)
    primitive=lines*ell**(3*(r-1))
    if r>=2:
        primitive+=lines*ell**(2*r-3)*(ell**(r-1)-1)//(ell-1)
    total=primitive+(prime_counts(ell,r-2)[1] if r>=2 else 0)
    return primitive,total


def counts(multiplier):
    if type(multiplier)!=int or multiplier<1:
        raise ValueError('The polarization multiplier must be a positive integer')
    primitive=total=1
    for ell,r in sp.factorint(multiplier).items():
        a,b=prime_counts(int(ell),int(r))
        primitive*=a
        total*=b
    return {'polarization_multiplier':multiplier, 'isogeny_degree':multiplier**2,
            'primitive_targets':primitive, 'all_targets_at_this_degree':total}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('multiplier',type=int,nargs='+')
    args=parser.parse_args()
    print(json.dumps([counts(n) for n in args.multiplier],indent=2))
