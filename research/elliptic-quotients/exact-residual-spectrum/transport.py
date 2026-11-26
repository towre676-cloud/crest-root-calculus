"""Check the metric incidence transfer for the two residual projectors."""
from fractions import Fraction
from itertools import combinations
import json

from exact import ambient_genus, coefficient_projector, depleted_sturm, gram_matrix, kernel, pari, residue_field
from star import act, cross, fiber, normal_partition, orbits


def representatives(label_map, sizes):
    out = [None]*len(sizes)
    for item,label in label_map.items():
        if out[label] is None:
            out[label] = item
    return out


def incidence(functionals, points, point_orbits, ell):
    entries = []
    width = len(set(point_orbits.values()))
    for functional in functionals:
        row = [0]*width
        for point in points:
            if sum(a*b for a,b in zip(functional,point))%ell==0:
                row[point_orbits[point]] += 1
        entries.extend(row)
    return pari.matrix(len(functionals),width,entries)


def blocks(Q, automorphisms, field):
    ell = field.p
    vertices,internal = normal_partition(Q,field)
    actions = [lambda v,g=g: act(g,v,ell) for g in automorphisms]
    point_orbits, point_sizes = orbits(internal,actions)
    tables = [{v:a(v) for v in vertices} for a in actions]
    edge_actions = [lambda e,t=t: tuple(sorted((t[e[0]],t[e[1]]))) for t in tables]
    edge_orbits,edge_sizes = orbits(combinations(sorted(vertices),2),edge_actions)
    plane_s = [cross(*edge) for edge in representatives(edge_orbits,edge_sizes)]
    plane_n = [tuple(sum(a*b for a,b in zip(row,v)) for row in Q)
               for v in representatives(point_orbits,point_sizes)]
    B = incidence(plane_s,internal,point_orbits,ell)
    C = incidence(plane_n,internal,point_orbits,ell)
    widths = [point_sizes[i] for i in range(len(point_sizes))]
    Wn = pari.matdiagonal(widths)
    Ws = pari.matdiagonal([edge_sizes[i] for i in range(len(edge_sizes))])
    if C.mattranspose()*Wn != Wn*C:
        raise ArithmeticError("The nonsplit incidence must be self-adjoint.")
    adjoint = Wn**(-1) * B.mattranspose() * Ws
    normal = adjoint * B
    J = pari.matrix(len(widths),len(widths),widths*len(widths))
    if normal+C*C != ell*pari.matid(len(widths))+J:
        raise ArithmeticError("The projective incidence Gram identity failed.")
    left_inverse = B*normal**(-1)
    lift = left_inverse*C
    if C**(-1)*adjoint*lift != pari.matid(len(widths)):
        raise ArithmeticError("The minimum-norm lift does not split the adjoint transfer.")
    forced = pari.matid(len(edge_sizes))-left_inverse*adjoint
    return lift,forced


def block_diagonal(matrices):
    return pari.matconcat(pari.matrix(len(matrices),len(matrices),[
        left if i==j else pari.matrix(int(left.nrows()),int(right.ncols()))
        for i,left in enumerate(matrices) for j,right in enumerate(matrices)]))


def transported_projector(D,ell):
    field = residue_field(D,ell)
    records = {"split":[],"nonsplit":[]}
    maps,forced_maps = [],[]
    for Q,group in ambient_genus(D):
        fiber(Q,group,field,records)
        lift,forced = blocks(Q,group,field)
        maps.append(lift)
        forced_maps.append(forced)
    lift = block_diagonal(maps)
    _,bound = depleted_sturm(D,ell)
    s,n = ([gram_matrix(Q) for Q,_ in records[g]] for g in ("split","nonsplit"))
    ks,kn = (kernel(forms,D,ell,bound) for forms in (s,n))
    aut_s = [int(q.qfauto()[0]) for q in s]
    aut_n = [int(q.qfauto()[0]) for q in n]
    forced = [rho for _,rho in records["split"]]
    Ws = pari.matdiagonal(aut_s)
    forced_projector = Ws**(-1)*block_diagonal(forced_maps)*Ws
    if [Fraction(str(forced_projector[i,i])) for i in range(len(s))] != forced:
        raise ArithmeticError("The explicit forced projector disagrees with the star diagonal formula.")
    direct = coefficient_projector(ks,aut_s)-forced_projector
    if len(kn):
        pushed = Ws**(-1)*lift*pari.matdiagonal(aut_n)*kn
        transported = coefficient_projector(pushed,aut_s)
    else:
        transported = pari.matrix(len(s),len(s))
    if transported != direct or len(ks)-int(sum(forced)) != len(kn):
        raise ArithmeticError("Transported and direct residual projectors disagree.")
    return {"D":D,"ell":ell,"split_classes":len(s),"nonsplit_classes":len(n),
            "residual_dimension":len(kn),"direct_equals_transport":True,
            "max_residual":str(max(Fraction(str(transported[i,i])) for i in range(len(s))))}


if __name__ == "__main__":
    pari.allocatemem(32_000_000,768_000_000,silent=True)
    print(json.dumps([transported_projector(D,ell) for D,ell in [(-8,13),(-3,17),(-11,7)]],indent=2))
