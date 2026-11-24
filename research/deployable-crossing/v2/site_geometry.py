"""Small planar geometry helpers for mapped feasibility screening.

GIS arithmetic uses binary64 coordinates with 1e-8 m predicate tolerance.
This is separate from exact structural arithmetic and does not assert source
accuracy. Rings use an even/odd fill rule, preserving holes and multipart
ArcGIS polygons. Boundary contact counts as a conflict.
"""
import math

EPS=1e-8


def cross(a,b,c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def bbox(points):
    return (min(p[0] for p in points),min(p[1] for p in points),
            max(p[0] for p in points),max(p[1] for p in points))


def overlap(a,b,pad=0):
    return not(a[2]+pad<b[0] or b[2]+pad<a[0] or a[3]+pad<b[1] or b[3]+pad<a[1])


def edges(ring):
    return zip(ring,ring[1:]+ring[:1])


def point_segment(p,a,b):
    dx,dy=b[0]-a[0],b[1]-a[1]
    den=dx*dx+dy*dy
    t=max(0,min(1,((p[0]-a[0])*dx+(p[1]-a[1])*dy)/den)) if den else 0
    return math.hypot(p[0]-a[0]-t*dx,p[1]-a[1]-t*dy)


def segment_hit(a,b,c,d):
    if not overlap(bbox([a,b]),bbox([c,d]),EPS): return False
    v=[cross(a,b,c),cross(a,b,d),cross(c,d,a),cross(c,d,b)]
    if v[0]*v[1]<0 and v[2]*v[3]<0: return True
    return min(point_segment(a,c,d),point_segment(b,c,d),
               point_segment(c,a,b),point_segment(d,a,b))<=EPS


def inside(point,rings):
    hit=False
    for ring in rings:
        for a,b in edges(ring):
            if point_segment(point,a,b)<=EPS: return True
            if (a[1]>point[1]) != (b[1]>point[1]):
                x=a[0]+(point[1]-a[1])*(b[0]-a[0])/(b[1]-a[1])
                if point[0]<x: hit=not hit
    return hit


def polygon_hit(a,b,buffer=0):
    """Two ring collections, with optional nonnegative separation sensitivity."""
    if buffer<0: raise ValueError("Buffer must be nonnegative")
    if not overlap(bbox([p for r in a for p in r]),bbox([p for r in b for p in r]),buffer):
        return False
    if any(inside(p,b) for r in a for p in r) or any(inside(p,a) for r in b for p in r):
        return True
    for ra in a:
        for rb in b:
            for p,q in edges(ra):
                for r,s in edges(rb):
                    if segment_hit(p,q,r,s): return True
                    if buffer and min(point_segment(p,r,s),point_segment(q,r,s),
                                      point_segment(r,p,q),point_segment(s,p,q))<=buffer+EPS:
                        return True
    return False


def line_hit(points,rings,radius=0):
    if radius<0: raise ValueError("Radius must be nonnegative")
    if not overlap(bbox(points),bbox([p for r in rings for p in r]),radius): return False
    if any(inside(p,rings) for p in points): return True
    for ring in rings:
        for a,b in edges(ring):
            for c,d in zip(points,points[1:]):
                if segment_hit(a,b,c,d): return True
                if min(point_segment(a,c,d),point_segment(b,c,d),
                       point_segment(c,a,b),point_segment(d,a,b))<=radius+EPS: return True
    if len(points)==1:
        return min(point_segment(points[0],a,b) for r in rings for a,b in edges(r))<=radius+EPS
    return False


def rectangle(origin,tangent,normal,s0,s1,n0,n1):
    return [[origin[0]+s*tangent[0]+n*normal[0],origin[1]+s*tangent[1]+n*normal[1]]
            for s,n in [(s0,n0),(s1,n0),(s1,n1),(s0,n1)]]


def polygon_samples(ring,spacing=2):
    """Vertices, edge samples and internal regular samples; not an extrema proof."""
    out=list(ring)
    for a,b in edges(ring):
        n=max(1,math.ceil(math.dist(a,b)/spacing))
        out.extend([[a[0]+(b[0]-a[0])*i/n,a[1]+(b[1]-a[1])*i/n] for i in range(1,n)])
    x0,y0,x1,y1=bbox(ring)
    for i in range(math.ceil((x1-x0)/spacing)+1):
        for j in range(math.ceil((y1-y0)/spacing)+1):
            p=[x0+i*spacing,y0+j*spacing]
            if inside(p,[ring]): out.append(p)
    return out
