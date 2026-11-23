"""Exact bounded cubic indexing from observed powder d-spacing intervals."""
from fractions import Fraction as Q
from math import sqrt


def rational(v):
    if isinstance(v,bool) or not isinstance(v,(str,float,int)) or len(str(v))>48:
        raise ValueError("use finite decimal bounds")
    try:
        x=Q(str(v))
    except (ValueError,ZeroDivisionError) as error:
        raise ValueError("use finite decimal bounds") from error
    if not 0 < x <= 10000:
        raise ValueError("distances must be positive and at most 10000 angstroms")
    return x


def interval(pair):
    if not isinstance(pair,list) or len(pair)!=2:
        raise ValueError("each interval needs two endpoints")
    lo,hi=map(rational,pair)
    if lo>hi: raise ValueError("interval endpoints out of order")
    return lo,hi


def merge(bands):
    out=[]
    for lo,hi in sorted(bands):
        if out and lo<=out[-1][1]:
            out[-1]=(out[-1][0],max(hi,out[-1][1]))
        else: out.append((lo,hi))
    return out


def intersect(a,b):
    out=[]
    i=j=0
    while i<len(a) and j<len(b):
        lo=max(a[i][0],b[j][0])
        hi=min(a[i][1],b[j][1])
        if lo<=hi: out.append((lo,hi))
        if a[i][1]<b[j][1]: i+=1
        else: j+=1
    return merge(out)


def solve(data):
    peaks=data.get('d_spacing_intervals_A')
    if not isinstance(peaks,list) or not 1<=len(peaks)<=80:
        raise ValueError("provide 1 to 80 powder d-spacing intervals")
    peaks=sorted(interval(p) for p in peaks)
    if any(a[1]>=b[0] for a,b in zip(peaks,peaks[1:])):
        raise ValueError("resolved peaks must have disjoint d-spacing intervals")
    cell=interval(data.get('cell_range_A'))
    classes=data.get('centerings',['P','I','F'])
    if not isinstance(classes,list) or not classes or any(c not in ('P','I','F') for c in classes) or len(set(classes))!=len(classes):
        raise ValueError("centerings must be a nonempty list of distinct entries from P, I, F")
    # Every relevant index satisfies h <= a_max/d_min. No arbitrary index truncation.
    limit=int(cell[1]/peaks[0][0])
    if limit>32: raise ValueError("narrow the physical cell range so each index is at most 32")
    retained={}
    for mode in classes:
        nset=set()
        for h in range(1,limit+1):
            for k in range(h+1):
                for l in range(k+1):
                    if mode=='I' and (h+k+l)%2: continue
                    if mode=='F' and not (h%2==k%2==l%2): continue
                    nset.add(h*h+k*k+l*l)
        bands=[(cell[0]**2,cell[1]**2)]
        for lo,hi in peaks:
            bands=intersect(bands,merge((n*lo**2,n*hi**2) for n in nset))
            if not bands: break
        if bands:
            retained[mode]=[{'a_squared_A2_exact':[str(lo),str(hi)],
                             'a_A_display_approx':[round(sqrt(lo),7),round(sqrt(hi),7)]} for lo,hi in bands]
    # For centering only. The report preserves even disconnected cell-length bands.
    status='incompatible' if not retained else 'centering_unique' if len(retained)==1 else 'ambiguous'
    return {'status':status,'complete':True,'model':'single cubic phase, observed line positions only',
            'classes':retained,'index_bound':limit,'peaks_used':len(peaks)}
