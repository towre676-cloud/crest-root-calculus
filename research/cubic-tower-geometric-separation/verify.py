"""Check geometric factor certificates; --recount rebuilds every local observation."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path
import sys

from affine_frobenius import signed_numerators
from geometric import fingerprint

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'cubic-tower-factors'))
from tower import check_stage, trace_from_factor, trace_sums
sys.path.insert(0,str(HERE.parent/'cubic-tower-cartier'))
from geometric_dimensions import subdegrees
from cartier import Cartier, numerators


def label(support, sign):
    return ','.join(map(str,support))+'/'+('+-'[sign])


def parse_label(name):
    support, sign = name.split('/')
    return tuple(map(int,support.split(','))), '+-'.index(sign)


@lru_cache(None)
def good(p,support):
    check_stage(p,support)
    return True


def profiles_for(group):
    n,k = group['indices'],group['size']
    return [label(s,e) for s in combinations(range(n),k) for e in (0,1)]


def recompute(job):
    p,S,expected = job
    vals = signed_numerators(p,S)
    for sign,coeff in expected:
        if vals[sign]!=coeff:
            raise ArithmeticError(f'Local Frobenius mismatch at {p} on {label(S,sign)}')
    return 1


def verify(data, recount, workers):
    all_polys=data['unit_polynomials']
    own=data['simplicity']
    groups=data['groups']
    required={label(S,e) for group in groups if group['size']>=4
              for S in combinations(range(group['indices']),group['size']) for e in (0,1)}
    if set(own)!=required:
        raise ArithmeticError('The dimension certificates do not cover the claimed supports')
    if [(g['indices'],g['size']) for g in groups] != [(10,3),(10,4),(6,5),(6,6)]:
        raise ArithmeticError('Unexpected scope of separation groups')
    used=set()
    def observe(p,name):
        S,sign=parse_label(name)
        d=len(S)*2**(len(S)-1)
        if not good(p,S) or p<=d:
            raise ArithmeticError('Wrong characteristic for the trace certificate')
        vals=all_polys[str(p)][name]
        if len(vals)!=d+1 or vals[0]!=1 or any(type(x)!=int or x<0 or x>=p for x in vals):
            raise ArithmeticError('Noncanonical local polynomial')
        used.add((p,name))
        return vals
    for name,primes in own.items():
        S,sign=parse_label(name)
        d=len(S)*2**(len(S)-1)
        possible=set(range(1,d//2+1))
        for p in primes:
            allowed=subdegrees(p,tuple(observe(p,name)))
            if allowed is None:
                raise ArithmeticError('A dimension witness is nonordinary')
            possible.intersection_update(allowed)
        if possible:
            raise ArithmeticError(f'{name} has surviving dimensions {sorted(possible)}')
    for group in groups:
        names=profiles_for(group)
        k=group['size']; d=k*2**(k-1)
        pairs=set(combinations(names,2))
        for row in group['separation']:
            p=row['p']
            needed={name for pair in pairs for name in pair}
            fps={name:fingerprint(observe(p,name),p,d) for name in needed}
            pairs={pair for pair in pairs if fps[pair[0]]==fps[pair[1]]}
            if sorted(map(list,pairs))!=row['remaining_pairs']:
                raise ArithmeticError('Geometric separation record disagrees')
        if pairs:
            raise ArithmeticError('Some equal-dimensional pairs were not separated')
    stored={(int(p),name) for p,table in all_polys.items() for name in table}
    if stored!=used:
        raise ArithmeticError('Unused local observations in the certificate')
    jobs={}
    for p,name in sorted(stored):
        S,sign=parse_label(name)
        jobs.setdefault((p,S),[]).append((sign,all_polys[str(p)][name]))
    if recount:
        work=[(p,S,pairs) for (p,S),pairs in jobs.items()]
    else:
        # Replay one local matrix in each dimension in the short mode.
        chosen=[next((p,S,pairs) for (p,S),pairs in jobs.items() if len(S)==k) for k in (3,4,5,6)]
        work=chosen
    if workers==1:
        count=sum(map(recompute,work))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            count=sum(pool.map(recompute,work,chunksize=3))
    # Independent point counts and archived integral supports at 41.
    for p,k in ((41,3),(41,4),(83,5),(193,6)):
        S=tuple(range(k))
        vals=[observe(p,label(S,e)) for e in range(2)]
        if k==4 and vals!=numerators(Cartier(p,k).matrices(S),p):
            raise ArithmeticError('Independent Cartier matrix disagrees in dimension 32')
        for extension in (1,2):
            traces=trace_sums(p,S,extension,k)[1][S]
            for coeff,trace in zip(vals,traces):
                if (trace_from_factor(coeff,extension)-trace)%p:
                    raise ArithmeticError('Direct tower trace disagrees')
    p=41
    historical=json.loads((HERE.parent/'cubic-tower-factors/receipt.json').read_text())
    exact=historical['stage_ten']['complete_factors_low_to_high']
    for k in (1,2):
        for S in combinations(range(10),k):
            for sign,vals in enumerate(signed_numerators(p,S)):
                if vals!=[v%p for v in exact[label(S,sign)][:len(vals)]]:
                    raise ArithmeticError('Earlier exact Frobenius polynomial disagrees')
    print(f'Verified {len(own)} higher-factor simplicity certificates and four geometric separation groups. '
          f'Recounted {count}/{len(jobs)} higher-support prime/support combinations.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--recount',action='store_true',help='recompute all stored local matrix numerators')
    parser.add_argument('--workers',type=int,default=1)
    parser.add_argument('--cartier-cross',action='store_true',help='compare all 420 four-index polynomials at 41 by holomorphic differentials')
    args=parser.parse_args()
    data=json.loads((HERE/'certificate.json').read_text())
    verify(data,args.recount,args.workers)
    if args.cartier_cross:
        from hashlib import sha256
        from cartier import Cartier, numerators
        engine=Cartier(41,4)
        digest=sha256()
        for S in combinations(range(10),4):
            values=numerators(engine.matrices(S),41)
            for e,poly in enumerate(values):
                if poly!=data['unit_polynomials']['41'][label(S,e)]:
                    raise ArithmeticError('Independent four-index Cartier comparison failed')
                digest.update(bytes(poly))
        print('All 420 independent degree-32 polynomials agree at 41. SHA-256',digest.hexdigest())
