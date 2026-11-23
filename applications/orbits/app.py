"""Schedule transit observations under bounded uncertainty in period and reference epoch."""
import json
import sys
from fractions import Fraction as F
from pathlib import Path
from math import ceil, floor

ROOT = Path(__file__).resolve().parent
MAX_CYCLES = 10000
MAX_PIECES = 512

def num(x):
    if isinstance(x, bool) or not isinstance(x,(str,int,float)):
        raise ValueError('Times must be finite decimal numbers')
    return F(str(x))

def between(v):
    if not isinstance(v,list) or len(v)!=2:
        raise ValueError('Expected [lower,upper] bounds')
    a,b = map(num,v)
    if a>b: raise ValueError('Bounds are reversed')
    return a,b

def clip(poly, n, bound, sign):
    """Keep sign*(epoch+n*period-bound)<=0. Retain closed edges and points."""
    out=[]
    for i,cur in enumerate(poly):
        prev = poly[i-1]
        a=sign*(prev[0]+n*prev[1]-bound); b=sign*(cur[0]+n*cur[1]-bound)
        if (a>0 and b<0) or (a<0 and b>0):
            f = a/(a-b)
            out.append((prev[0]+f*(cur[0]-prev[0]),prev[1]+f*(cur[1]-prev[1])))
        if b<=0: out.append(cur)
    return list(dict.fromkeys(out))

def cycles(poly, lo, hi):
    ratios=[(t-e)/p for e,p in poly for t in (lo,hi)]
    start,end=ceil(min(ratios)),floor(max(ratios))
    if end-start+1>MAX_CYCLES: raise ValueError('Too many cycle aliases; narrow period/time bounds')
    return range(start,end+1)

def fit(model, events):
    elo,ehi=between(model['epoch_days']); plo,phi=between(model['period_days'])
    if plo<=0: raise ValueError('Periods must be positive')
    dlo,dhi=between(model['duration_hours'])
    if dlo<=0 or dhi/24>=plo: raise ValueError('Transit duration must be positive and shorter than a period')
    polys=[list(dict.fromkeys([(elo,plo),(ehi,plo),(ehi,phi),(elo,phi)]))]
    for lo,hi in events:
        pieces=[]
        for poly in polys:
            for n in cycles(poly,lo,hi):
                cut=clip(clip(poly,n,hi,1),n,lo,-1)
                if cut: pieces.append(cut)
                if len(pieces)>MAX_PIECES: raise ValueError('More than 512 distinct cycle histories; narrow bounds')
        polys=pieces
    return polys,(dlo/24,dhi/24)

def prediction(polys,duration,window):
    """Possible outcomes for continuous coverage with complete-transit detection."""
    left,right=window
    shortest,longest=duration
    outcomes=set()
    for poly in polys:
        # Existence of a complete event is easiest at the shortest transit duration.
        lo,hi=left+shortest/2,right-shortest/2
        if lo<=hi:
            # For fixed n the image of a convex polygon is an interval.
            if any(min(e+n*p for e,p in poly)<=hi and max(e+n*p for e,p in poly)>=lo
                   for n in cycles(poly,lo,hi)):
                outcomes.add(1)
        # Absence of a complete event is easiest at the longest duration.
        lo,hi=left+longest/2,right-longest/2
        if lo>hi:
            outcomes.add(0)
        else:
            indexes=[floor((lo-e)/p) for e,p in poly]
            if max(indexes)-min(indexes)+2>MAX_CYCLES:
                raise ValueError('Timing ambiguity is too broad for this observation')
            for n in range(min(indexes)-1,max(indexes)+1):
                # Strict gap: t_n < lo AND t_(n+1) > hi. Clip closed then
                # test each strict slack; convexity gives simultaneous slack.
                gap=clip(clip(poly,n,lo,1),n+1,hi,-1)
                if gap and min(e+n*p for e,p in gap)<lo and max(e+(n+1)*p for e,p in gap)>hi:
                    outcomes.add(0)
                    break
    return outcomes

def evaluate(data):
    if not data.get('time_system') or not isinstance(data['time_system'],str):
        raise ValueError('Give a common time system for hypotheses, events and windows')
    events=[between(e) for e in data['events']]
    negatives=[between(w) for w in data.get('non_detections',[])]
    if any(a>=b for a,b in negatives): raise ValueError('Non-detection windows must have positive length')
    fits={}
    rejected=[]
    ids=set()
    for h in data['hypotheses']:
        name=h['name']
        if not isinstance(name,str) or not name or name in ids: raise ValueError('Hypotheses need unique names')
        ids.add(name); polys,duration=fit(h,events)
        polys = [p for p in polys if all(0 in prediction([p],duration,w) for w in negatives)]
        if polys: fits[name]=(polys,duration)
        else: rejected.append(name)
    if not ids: raise ValueError('Provide at least one hypothesis')
    slots=[]; names=set()
    for w in data['windows']:
        name=w['name']; left,right=between(w['days'])
        if not isinstance(name,str) or not name or name in names or left>=right:
            raise ValueError('Windows need distinct names and positive duration')
        names.add(name)
        # Availability and sensitivity are declared by the user before planning.
        if not isinstance(w.get('available'), bool): raise ValueError('Each window needs boolean available')
        if not w['available']: continue
        outcomes={key:prediction(polys,duration,(left,right)) for key,(polys,duration) in fits.items()}
        pairs=[[a,b] for i,a in enumerate(outcomes) for b in list(outcomes)[i+1:] if outcomes[a].isdisjoint(outcomes[b])]
        labels={key: 'complete transit' if value=={1} else 'no complete transit' if value=={0} else 'ambiguous' for key,value in outcomes.items()}
        hours=(right-left)*24
        slots.append((F(len(pairs))/hours,len(pairs),left,name,
                      {'name':name,'hours':round(float(hours),5),'outcomes':labels,'distinguished_pairs':pairs}))
    slots.sort(key=lambda s:(-s[0],-s[1],s[2],s[3]))
    result={'time_system':data['time_system'],'surviving_hypotheses':list(fits),'excluded_hypotheses':rejected,
            'feasible_cycle_histories':{n:len(v[0]) for n,v in fits.items()},
            'recommended_window':None,'windows':[s[-1] for s in slots]}
    if not fits: decision='model_mismatch'
    elif len(fits)==1: decision='one_listed_ephemeris_remains'
    elif slots and slots[0][1]>0: decision='schedule_discriminating_window';result['recommended_window']=slots[0][3]
    else: decision='no_guaranteed_separation'
    return {'decision':decision,**result}

def main():
    if len(sys.argv)>2: sys.exit('Usage: python app.py [input.json]')
    try:
        p=Path(sys.argv[1]) if len(sys.argv)>1 else ROOT/'example.json'
        print(json.dumps(evaluate(json.loads(p.read_text(encoding='utf8'), parse_float=str)),indent=2))
    except (ValueError, TypeError, KeyError, OSError, ZeroDivisionError) as exc:
        sys.exit('Input error: '+str(exc))
if __name__=='__main__':
    main()
