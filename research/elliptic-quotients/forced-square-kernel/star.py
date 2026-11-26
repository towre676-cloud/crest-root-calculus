"""Exact forced kernels from incidence with the isotropic residues at ell."""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "hecke-observation"))
from controls import forms_of_determinant, isometric_primary, primary_discriminant_form
from neighbours import gram_matrix, pari
from observe import residue_field
from census import CLASS_NUMBER_ONE


def matrix_entries(matrix):
    return tuple(tuple(int(matrix[i,j]) for j in range(3)) for i in range(3))


def ambient_genus(D):
    t = D%2
    initial = gram_matrix([[1,0,0],[0,4,2*t],[0,2*t,t*t-D]])
    primes = [int(p) for p in pari.factor(-8*D)[0][0:]]
    data = [primary_discriminant_form(initial,p) for p in primes]
    forms = [initial]
    for gram in forms_of_determinant(-4*D):
        q = gram_matrix(gram)
        if any(q.qfisom(representative) for representative in forms):
            continue
        if all(isometric_primary(primary_discriminant_form(q,p),target) for p,target in zip(primes,data)):
            forms.append(q)
    return [(matrix_entries(q),[matrix_entries(g) for g in pari.qfauto(q)[1]]) for q in forms]


def projective(v, ell):
    leading = next(x for x in v if x%ell)
    scale = pow(leading,-1,ell)
    return tuple(x*scale%ell for x in v)


def act(g, v, ell):
    return projective(tuple(sum(a*b for a,b in zip(row,v)) for row in g),ell)


def orbits(elements, generators):
    members = {}
    sizes = Counter()
    for element in elements:
        if element in members:
            continue
        name = len(sizes)
        members[element] = name
        stack = [element]
        while stack:
            x = stack.pop()
            sizes[name] += 1
            for action in generators:
                y = action(x)
                if y not in members:
                    members[y] = name
                    stack.append(y)
    return members,sizes


def normal_partition(Q, field):
    p = field.p
    a,b,c = (Q[i][i] for i in range(3))
    d,e,f = 2*Q[0][1],2*Q[0][2],2*Q[1][2]
    points = [(1,x,y) for x,y in product(range(p),repeat=2)]
    points += [(0,1,y) for y in range(p)] + [(0,0,1)]
    isotropic,nonsplit = [],[]
    for u,x,y in points:
        value = (a*u*u+b*x*x+c*y*y+d*u*x+e*u*y+f*x*y)%p
        if not value:
            isotropic.append((u,x,y))
        elif field.legendre(value)==1:
            nonsplit.append((u,x,y))
    return isotropic,nonsplit


def cross(v, w):
    return tuple(v[(j+1)%3]*w[(j+2)%3] - v[(j+2)%3]*w[(j+1)%3] for j in range(3))


def sublattice(Q, functional, ell):
    pivot = next(j for j,t in enumerate(functional) if t%ell)
    basis = [[int(i==j) for i in range(3)] for j in range(3)]
    inverse = pow(functional[pivot],-1,ell)
    for j in range(3):
        basis[j][pivot] = ell if j==pivot else -functional[j]*inverse%ell
    return [[sum(v[i]*Q[i][j]*w[j] for i,j in product(range(3),repeat=2)) for w in basis] for v in basis]


def fiber(Q, autos, field, representatives=None):
    ell = field.p
    isotropic,nonsplit = normal_partition(Q,field)
    if len(isotropic)!=ell+1 or len(nonsplit)!=ell*(ell-1)//2:
        raise ArithmeticError("Unexpected local orbit sizes.")
    line_actions = [lambda v,g=g: act(g,v,ell) for g in autos]
    vertex,vertex_sizes = orbits(isotropic,line_actions)
    # Precompute the group action on vertices to keep the edge traversal small.
    tables = [{v:action(v) for v in isotropic} for action in line_actions]
    edge_actions = [lambda edge,t=t: tuple(sorted((t[edge[0]],t[edge[1]]))) for t in tables]
    edges,edge_sizes = orbits(combinations(sorted(isotropic),2),edge_actions)
    normals,nonsplit_sizes = orbits(nonsplit,line_actions)
    diagonals = {}
    n = len(isotropic)
    for edge,label in edges.items():
        if label in diagonals:
            continue
        endpoints = Counter(vertex[v] for v in edge)
        diagonal = edge_sizes[label]*(sum(Fraction(k*k,vertex_sizes[o]) for o,k in endpoints.items())
                                     - Fraction(4,n))/(n-2)
        diagonals[label] = diagonal
        if representatives is not None:
            linear = cross(*edge)
            representatives["split"].append((sublattice(Q,linear,ell),diagonal))
    if representatives is not None:
        seen = set()
        for normal,label in normals.items():
            if label in seen:
                continue
            seen.add(label)
            linear = tuple(sum(a*b for a,b in zip(row,normal)) for row in Q)
            representatives["nonsplit"].append((sublattice(Q,linear,ell),Fraction(0)))
    if sum(diagonals.values())!=len(vertex_sizes)-1:
        raise ArithmeticError("Incidence projector trace disagrees with its expected rank.")
    if len(edge_sizes)-len(nonsplit_sizes)!=len(vertex_sizes)-1:
        raise ArithmeticError("The split and nonsplit class numbers do not satisfy the incidence relation.")
    return {"split_classes":len(edge_sizes),"nonsplit_classes":len(nonsplit_sizes),
            "vertex_orbits":len(vertex_sizes),"forced_kernel":len(vertex_sizes)-1,
            "max_diagonal":str(max(diagonals.values()))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime-limit",type=int,default=199)
    parser.add_argument("--verify",type=Path)
    args = parser.parse_args()
    all_levels=[]
    for D in CLASS_NUMBER_ONE:
        amb = ambient_genus(D)
        for ell in range(3,args.prime_limit+1,2):
            try:
                field = residue_field(D,ell)
            except ValueError:
                continue
            blocks = [fiber(Q,autos,field) for Q,autos in amb]
            all_levels.append({"D":D,"ell":ell,"ambient_classes":len(amb),
                               "split_classes":sum(r["split_classes"] for r in blocks),
                               "nonsplit_classes":sum(r["nonsplit_classes"] for r in blocks),
                               "forced_kernel":sum(r["forced_kernel"] for r in blocks),
                               "max_diagonal":str(max(Fraction(r["max_diagonal"]) for r in blocks)),
                               "source_ambient_forced_rank":blocks[0]["forced_kernel"]})
        print(D,len(amb),"ambient classes",file=sys.stderr,flush=True)
    if args.verify and all_levels!=json.loads(args.verify.read_text()):
        raise ArithmeticError("The computed incidence census differs from the receipt.")
    print(json.dumps(all_levels,indent=2))


if __name__ == "__main__":
    main()
