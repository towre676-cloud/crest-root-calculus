"""Exclude rational primitive polarized isogenies at all powers of each certified odd prime."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_gcd, gf_pow_mod, gf_sub

from frobenius import frobenius, good, polynomial


def irreducible_quartic(h, ell):
    # A reducible quartic has a factor of degree one or two. Both divide x^(ell^2)-x.
    f=[int(c)%ell for c in h]
    x=[1,0]
    difference=gf_sub(gf_pow_mod(x,ell*ell,f,ell,ZZ),x,ell,ZZ)
    return len(gf_gcd(f,difference,ell,ZZ))==1


def quadratic_tests(weil,p):
    a,b=weil[1:3]
    return a*a-4*(b-2*p),(b+2*p)**2-4*p*a*a


def prove_irreducible(pair,ell):
    return all(pow(d%ell,(ell-1)//2,ell)==ell-1 for d in pair)


def independent_characters(source,primes):
    # Record row echelon coordinates in the Q*/Q*² basis of prime squareclasses and -1.
    columns={};basis={};selected=[]
    def vector(d):
        v=0
        for q,e in sp.factorint(d).items():
            if e%2:
                if q not in columns:
                    columns[q]=len(columns)
                v^=1<<columns[q]
        return v
    for p in primes:
        if not good(source,p):
            continue
        counts, weil=frobenius(source,p);tests=quadratic_tests(weil,p)
        if 0 in tests:
            continue
        rows=basis.copy()
        for d in tests:
            v=vector(d)
            while v:
                i=v.bit_length()-1
                if i not in rows:
                    rows[i]=v
                    break
                v^=rows[i]
        if len(rows)==len(basis)+2:
            basis=rows
            selected.append({'p':p,'counts':counts,'weil_descending':weil,'quadratic_tests':tests})
    return {'independent_witness_pairs':len(selected),'independent_character_rank':len(basis),
            'possible_prime_density_upper_bound':'(3/4)^'+str(len(selected)),'witnesses':selected}


def exclude(source,bound,primes,independent=False):
    todo=set(map(int,sp.primerange(3,bound+1)))
    prime_count=len(todo)
    witnesses={}
    by_ell={}
    for p in primes:
        if not todo:
            break
        if not good(source,p):
            continue
        counts,weil=frobenius(source,p)
        pair=quadratic_tests(weil,p)
        settled=[ell for ell in sorted(todo) if ell!=p and (irreducible_quartic(weil,ell) if independent else prove_irreducible(pair,ell))]
        if settled:
            witnesses[p]={'counts':counts,'weil_descending':weil,'assigned_primes':len(settled)}
            by_ell.update({ell:p for ell in settled})
            todo.difference_update(settled)
    pairs=sorted(by_ell.items())
    digest=sha256(json.dumps(pairs,separators=(',',':')).encode()).hexdigest()
    return {'bound':bound, 'odd_primes_checked':prime_count,'all_powers_excluded':len(by_ell),
            'unresolved_primes':sorted(todo),'assignments_sha256':digest,'witnesses':witnesses}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('--bound',type=int,default=10000)
    parser.add_argument('--witness-limit',type=int,default=251)
    parser.add_argument('--density',action='store_true',help='certify a positive-density infinite set of excluded prime multipliers')
    args=parser.parse_args()
    source=polynomial(json.loads(args.input.read_text())['source'])
    primes=sp.primerange(3,args.witness_limit+1)
    result=independent_characters(source,primes) if args.density else exclude(source,args.bound,primes)
    print(json.dumps(result,indent=2))
