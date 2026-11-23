"""Build and reopen the reference STEP assembly, then tessellate the same solids for the offline viewer.
Run with CadQuery 2.7.0. Finished dimensions are in mm. O is the XY origin.
"""
import argparse
import array
import base64
from itertools import combinations, product
import json
import math
from pathlib import Path
import sys

import cadquery as cq

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, default=HERE, help='output directory for the STEP, mesh and receipt')
out = parser.parse_args().output
out.mkdir(parents=True, exist_ok=True)
spec = json.loads((HERE.parent/'kit/design.json').read_text())
s, j, base, link = (spec[k] for k in ('spacings','joint','base','links'))
O, D, A = (0.,0.), (s['ground'],0.), (s['crank'],0.)
d = D[0]-A[0]
along = (s['coupler']**2-s['rocker']**2+d*d)/(2*d)
B = (A[0]+along, math.sqrt(s['coupler']**2-along**2))
assert abs(math.dist(A,B)-s['coupler'])<1e-10 and abs(math.dist(D,B)-s['rocker'])<1e-10

colors = {'board':(0.81,0.76,0.64),'foot':(.16,.22,.24),'crank':(.16,.39,.59),
          'coupler':(.79,.40,.17),'rocker':(.05,.49,.43),'shaft':(.55,.64,.67),
          'collar':(.29,.37,.41),'washer':(.69,.74,.74),'spacer':(.79,.65,.40),'grip':(.78,.65,.42)}
parts = []
assembly = cq.Assembly(name=spec['revision'])
def add(name,category,body,motion='fixed',spread=0):
    shape=body.val()
    assert shape.isValid() and len(shape.Solids())==1 and shape.Volume()>0,name
    assembly.add(shape,name=name,color=cq.Color(*colors[category]))
    parts.append((name,category,shape,motion,spread))

def tube(p,z,length,od,bore=0):
    profile=cq.Workplane('XY',origin=(*p,z)).circle(od/2)
    if bore:
        profile=profile.circle(bore/2)
    return profile.extrude(length)

def cut_holes(solid,points,diameter,z,height):
    for point in points:
        solid=solid.cut(tube(point,z-1,height+2,diameter))
    return solid

def add_link(name,p,q,z,motion,extension=0):
    length=math.dist(p,q); phi=math.degrees(math.atan2(q[1]-p[1],q[0]-p[0]))
    axis=length+extension; w=link['width']; t=link['thickness']
    body=cq.Workplane('XY').center(axis/2,0).slot2D(axis+w,w).extrude(t)
    body=cut_holes(body,[(0,0),(length,0)],sum(j['moving_bore'])/2,0,t)
    body=body.rotate((0,0,0),(0,0,1),phi).translate((*p,z))
    add(name,name,body,motion,54 if name=='coupler' else 23)

x0,y0,x1,y1=(base[n] for n in ('xmin','ymin','xmax','ymax'))
bz=-base['thickness']
board=cq.Workplane('XY').box(x1-x0,y1-y0,-bz,centered=(False,False,False)).translate((x0,y0,bz))
add('base','board',cut_holes(board,[O,D],sum(j['ground_bore'])/2,bz,-bz))
inset,w,h=base['foot_inset'],base['foot_width'],base['foot_height']
for i,(x,y) in enumerate(product((x0+inset,x1-inset),(y0+inset,y1-inset)),1):
    add(f'foot_{i}','foot',cq.Workplane('XY').box(w,w,h).translate((x,y,bz-h/2)),'fixed',-12)
lower=spec['layers_nominal']['crank_rocker_z']
upper=spec['layers_nominal']['coupler_z']
add_link('crank',O,A,lower[0],'crank')
add_link('coupler',A,B,upper[0],'coupler')
add_link('rocker',D,B,lower[0],'rocker',link['pointer_extension'])
# Purchased parts use nominal envelope diameters. Screw and thread detail is suppressed.
collar=j['collar_height']; washer=j['washer_thickness']; play=sum(j['endplay'])/2
rod=max(j['shaft_diameter']); free=sum(j['free_end_projection'])/2
for letter,p in [('O',O),('D',D),('A',A),('B',B)]:
    fixed=letter in ('O','D')
    motion='fixed' if fixed else letter
    def hardware(label,category,z,h,od,bore,spread):
        add(f'{letter}_{label}',category,tube(p,z,h,od,bore),motion,spread)
    foot_z=bz-collar if fixed else lower[0]-washer-collar
    top_z=(lower if fixed else upper)[1]+washer+play
    top=top_z+collar
    grip_z=top+sum(j['grip_above_collar'])/2
    shaft_low=foot_z-free
    shaft_high=grip_z+sum(j['grip_engagement'])/2 if letter=='A' else top+free
    hardware('pin','shaft',shaft_low,shaft_high-shaft_low,rod,0,0 if fixed else 25)
    hardware('lower_collar','collar',foot_z,collar,j['collar_od_nominal'],rod,-6 if fixed else 12)
    hardware('upper_collar','collar',top_z,collar,j['collar_od_nominal'],rod,48 if fixed else 68)
    if fixed:
        hardware('base_collar','collar',0,collar,j['collar_od_nominal'],rod,3)
        hardware('spacer','spacer',collar,j['ground_spacer'],6,3.15,9)
        washers=[(lower[0]-washer,16),(lower[1],38)]
    else:
        hardware('spacer','spacer',lower[1]+washer,j['moving_spacer'],6,3.15,40)
        washers=[(lower[0]-washer,17),(lower[1],32),(upper[0]-washer,47),(upper[1],61)]
    for i,(z,gap) in enumerate(washers,1):
        hardware(f'washer_{i}','washer',z,washer,j['washer_od_nominal'],3.2,gap)
    if letter=='A':
        grip=tube(p,grip_z,j['grip_length'],j['grip_diameter'])
        # 10.00 mm shaft engagement in a 10.50 mm blind bore leaves room for adhesive.
        bore=tube(p,grip_z-1,11.5, sum(j['grip_blind_bore_diameter'])/2)
        add('A_finger_grip','grip',grip.cut(bore),motion,76)

def almost_contact(a,b):
    A,B=a.BoundingBox(),b.BoundingBox()
    return all(min(getattr(A,k+'max'),getattr(B,k+'max'))-max(getattr(A,k+'min'),getattr(B,k+'min'))>1e-6 for k in 'xyz')

assert len(parts)==39
for first,second in combinations(parts,2):
    a,b=first[2],second[2]
    if almost_contact(a,b):
        vol=a.intersect(b).Volume()
        assert vol<1e-6,(first[0],second[0],vol)
target=out/'reference-assembly.step'
assembly.export(str(target))
shapes=cq.importers.importStep(str(target)).solids().vals()
assert len(shapes)==len(parts)
def signature(shape):
    center=shape.Center()
    return tuple(round(x,4) for x in (shape.Volume(),center.x,center.y,center.z))
assert sorted(signature(p[2]) for p in parts)==sorted(signature(p) for p in shapes)
def packed(kind,values):
    values=array.array(kind,values)
    if sys.byteorder!='little': values.byteswap()
    return base64.b64encode(values.tobytes()).decode('ascii')
data=[]
for name,category,shape,motion,spread in parts:
    vertices,faces=shape.tessellate(.07,.26)
    data.append(dict(name=name,category=category,motion=motion,spread=spread,
                     xyz=packed('h',(round(v*100) for p in vertices for v in p.toTuple())),
                     triangles=packed('H',(v for face in faces for v in face))))
model={'revision':spec['revision'],'lengths':{key:s[key] for key in ['ground','crank','coupler','rocker']},
       'B0':B,'colors':colors,'parts':data}
(out/'mesh.js').write_text('window.LinkwrightAssembly='+json.dumps(model,separators=(',',':'))+';\n')
receipt={'revision':spec['revision'],'pose':'crank 0 deg, positive branch, rigid nominal',
         'solids_reimported':len(shapes),'pairs_screened_for_positive_volume':math.comb(len(parts),2),
         'maximum_allowed_intersection_mm3':1e-6,'part_comparison':'volume and centre of mass to 4 decimal places',
         'parts':[dict(name=name,volume_mm3=round(shape.Volume(),4)) for name,_,shape,_,_ in parts]}
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(f'{len(shapes)} valid STEP solids reimported; nominal zero-angle intersection and properties verified.')
