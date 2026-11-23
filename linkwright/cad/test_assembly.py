"""Independent STEP reimport and contact screening at 26 poses, including both reversals."""
import json, math, itertools, base64, struct
from pathlib import Path
import cadquery as cq

cad=Path(__file__).resolve().parent
m=json.loads((cad/'mesh.js').read_text()[len('window.LinkwrightAssembly='):-2])
a=cq.Assembly.load(str(cad/'reference-assembly.step'))
shapes={p.name:p.obj.located(p.loc) for p in a.children}
boms=json.loads((cad/'receipt.json').read_text())['parts']
assert len(shapes)==len(m['parts'])==len(boms)==39
assert all(x.isValid() and len(x.Solids())==1 for x in shapes.values())
assert {p['name'] for p in boms}=={x['name'] for x in m['parts']}==set(shapes)
for p in boms:
    assert abs(shapes[p['name']].Volume()-p['volume_mm3'])<.0001
xy0=(0,0);D=(93.8,0); A0=(46.6,0)
a0,b,c,g=46.6,84.1,111.3,93.8
def pose(t):
    A=(a0*math.cos(t),a0*math.sin(t))
    dx,dy=A[0]-g,A[1]
    d=math.hypot(dx,dy)
    phi=math.atan2(dy,dx)-math.acos((d*d+c*c-b*b)/(2*d*c))
    return A,(g+c*math.cos(phi),c*math.sin(phi))
def ang(p,q): return math.degrees(math.atan2(q[1]-p[1],q[0]-p[0]))
def addz(v):return (v[0],v[1],0)
def transform(group,t):
    A,B=pose(math.radians(t));B0=pose(0)[1]
    return {'fixed':(xy0,xy0,0),'A':(A0,A,0),'B':(B0,B,0),'crank':(xy0,xy0,t),
       'rocker':(D,D,ang(D,B)-ang(D,B0)), 'coupler':(A0,A,ang(A,B)-ang(A0,B0))}[group]

def moved(sh,xf):
    p,q,t=xf
    sh=sh.translate((-p[0],-p[1],0))
    if abs(t)>1e-10:sh=sh.rotate((0,0,0),(0,0,1),t)
    return sh.translate(addz(q))
def motion_group_for(name):
    if name in ['crank','rocker','coupler']: return name
    if name.startswith('A_') or name.startswith('B_'):return name[0]
    return 'fixed'
angles=list(range(0,360,15))+[56.611681,288.077649]
worst=0; screened=0; breps=0; closest=[]
for element in m['parts']:
    group=motion_group_for(element['name'])
    assert element['motion']==group
    raw=base64.b64decode(element['xyz'])
    pts=list(struct.iter_unpack('<hhh',raw))
    facebytes=base64.b64decode(element['triangles'])
    triangles=list(struct.iter_unpack('<HHH',facebytes))
    assert max(max(row) for row in triangles)<len(pts)
    bb=shapes[element['name']].BoundingBox()
    for i,k in enumerate('xyz'):
        lo,hi=(getattr(bb,k+suffix) for suffix in ('min','max'))
        assert min(p[i] for p in pts)/100>=lo-.011 and max(p[i] for p in pts)/100<=hi+.011

def allowed_contact(l,r):
    if l=='base' or r=='base':
        other=r if l=='base' else l
        return other.startswith('foot_') or other in [n+k for n in 'OD' for k in ('_pin','_lower_collar','_base_collar')]
    if l[0] in 'ODAB' and r[0] in 'ODAB' and r[0]==l[0]: return True
    for link,joints in [('crank','OA'),('coupler','AB'),('rocker','DB')]:
        if l==link and r[0] in joints and r[1]=='_':return True
        if r==link and l[0] in joints and l[1]=='_':return True
    return False

for t in angles:
    at={name:moved(sh,transform(motion_group_for(name),t)) for name,sh in shapes.items()}
    boxes={n:s.BoundingBox() for n,s in at.items()}
    for l,r in itertools.combinations(at,2):
        L,R=boxes[l],boxes[r]
        gap=[max(getattr(L,k+'min')-getattr(R,k+'max'),getattr(R,k+'min')-getattr(L,k+'max'),0) for k in 'xyz']
        screened+=1
        if math.hypot(*gap)>0.051 or allowed_contact(l,r):continue
        dl=at[l].distance(at[r])
        vol=at[l].intersect(at[r]).Volume() if dl<.001 else 0
        breps+=1;worst=max(worst,vol)
        assert dl>0.05 and vol<1e-6,(t,l,r,dl,vol)
    if t==0:
        assert abs(boxes['A_finger_grip'].zmax-47)<1e-6
        assert abs(boxes['A_pin'].zmin-2.7)<1e-6
        assert abs(boxes['D_pin'].zmin+11.3)<1e-6
print(f'{len(shapes)} positive-volume STEP solids, matched scene attachments; {len(angles)} moving poses, {screened} pair checks, {breps} narrow-phase nonincidence checks, zero unintended contact')
