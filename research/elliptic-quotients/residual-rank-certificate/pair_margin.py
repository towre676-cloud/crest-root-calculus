"""Exact projector energies of same-genus gluing-pair differences."""
import argparse
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/"exact-residual-spectrum"))
from exact import (ambient_genus, coefficient_projector, fiber, form_gram, gram_matrix,
                   kernel, pari, residue_field, scalar_orbits)
from star import act, normal_partition, orbits


def projection_energy(projector, i, j, weights):
    wi,wj = weights[i],weights[j]
    value = wi*(projector[i,i]-projector[i,j])+wj*(projector[j,j]-projector[j,i])
    return Fraction(str(value))/(wi+wj)


def star_columns(Q, group, field):
    vertices,_ = normal_partition(Q,field)
    tables = [{x:act(g,x,field.p) for x in vertices} for g in group]
    vertex,vertex_sizes = orbits(vertices,[lambda v,t=t:t[v] for t in tables])
    edges,edge_sizes = orbits(combinations(sorted(vertices),2),[
        lambda e,t=t:tuple(sorted((t[e[0]],t[e[1]]))) for t in tables])
    endpoints = {}
    for (u,v),orbit in edges.items():
        endpoints.setdefault(orbit,(vertex[u],vertex[v]))
    cols = [[edge_sizes[i]*(vertex_sizes[0]*endpoints[i].count(j)
                            -vertex_sizes[j]*endpoints[i].count(0)) for i in range(len(edge_sizes))]
            for j in range(1,len(vertex_sizes))]
    return len(edge_sizes),cols


def test_genus(D, ell, tag, records, params, forced_basis, certified_dimension, cutoff):
    if len(params)<2:
        return {"gluing_pairs":0}
    forms = [gram_matrix(Q) for Q,_ in records]
    weights = [int(Q.qfauto()[0]) for Q in forms]
    forced = coefficient_projector(forced_basis,weights)
    if len(forced_basis)==certified_dimension:
        projector = forced
    else:
        basis = kernel(forms,D,ell,cutoff)
        if certified_dimension is not None and len(basis)!=certified_dimension:
            raise ArithmeticError("The short kernel does not match the certified dimension.")
        projector = coefficient_projector(basis,weights)
        if projector*forced!=forced:
            raise ArithmeticError("The claimed forced directions are outside the computed kernel.")
    residual = projector-forced
    if any(Fraction(str(forced[i,i]))!=records[i][1] for i in range(len(forms))):
        raise ArithmeticError("The star basis disagrees with the incidence diagonals.")

    labels = {}
    field = residue_field(D,ell)
    for a in params:
        Q = gram_matrix(form_gram(field,a))
        index = next((i for i,R in enumerate(forms) if R.qfisom(Q)),None)
        if index is None or index in labels:
            raise ArithmeticError(f"Unmatched or duplicate gluing class in {tag}.")
        labels[index] = a
    largest_total = Fraction(-1)
    largest_residual = Fraction(-1)
    result = {}
    for i,j in combinations(labels,2):
        total = projection_energy(projector,i,j,weights)
        remainder = projection_energy(residual,i,j,weights)
        if not 0<=remainder<=total<1:
            raise ArithmeticError("The exact projector energies did not certify separation.")
        pair = [labels[i],labels[j]]
        if total>largest_total:
            largest_total = total
            result.update(maximum_total=str(total),minimum_visible=str(1-total),
                          maximum_total_pair=pair,forced_at_total_maximum=str(total-remainder),
                          residual_at_total_maximum=str(remainder))
        if remainder>largest_residual:
            largest_residual = remainder
            result.update(maximum_residual=str(remainder),maximum_residual_pair=pair)
    result["gluing_pairs"] = len(labels)*(len(labels)-1)//2
    return result


def run(D, ell, dimensions, ambient, cutoff):
    field = residue_field(D,ell)
    records = {"split":[],"nonsplit":[]}
    forced_blocks = []
    for Q,group in ambient:
        fiber(Q,group,field,records)
        forced_blocks.append(star_columns(Q,group,field))
    width = len(records["split"])
    columns = []
    offset = 0
    for size,block in forced_blocks:
        columns.extend(pari.Col([0]*offset+v+[0]*(width-offset-size)) for v in block)
        offset += size
    if offset!=width:
        raise ArithmeticError("The star basis uses a different list of plane orbits.")
    forced = {"split":pari.Mat(columns),"nonsplit":pari.Mat([])}
    params = {"split":[],"nonsplit":[]}
    for a in scalar_orbits(field):
        tag = "split" if field.legendre(1+field.norm(a))==-1 else "nonsplit"
        params[tag].append(a)
    return {"D":D,"ell":ell,
            **{g:test_genus(D,ell,g,records[g],params[g],forced[g],dimensions[g],cutoff)
               for g in records}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify",type=Path)
    parser.add_argument("--include-prefix",action="store_true",
                        help="Also bound the ten kernels that have only finite-count upper bounds.")
    args = parser.parse_args()
    forced = json.loads((HERE.parent/"forced-square-kernel"/"rank.json").read_text())
    residual = json.loads((HERE.parent/"exact-residual-spectrum"/"projectors.json").read_text())
    dimensions = {(r["D"],r["ell"]):{g:r[g]["actual_kernel_lower"] for g in ("split","nonsplit")}
                  for r in forced if r["residual_kernel_upper"]==0}
    dimensions.update({(r["D"],r["ell"]):{g:r[g]["exact_kernel_dimension"] for g in ("split","nonsplit")}
                       for r in residual})
    pari.allocatemem(32_000_000,768_000_000,silent=True)
    ambient = {D:ambient_genus(D) for D in {r["D"] for r in forced}}
    results = []
    for row in sorted(forced,key=lambda r:(r["D"],r["ell"])):
        D,ell = row["D"],row["ell"]
        key = (D,ell)
        if key in dimensions:
            result = run(D,ell,dimensions[key],ambient[D],300)
        elif args.include_prefix:
            cutoff = row["split"]["cutoff"]
            result = run(D,ell,{"split":None,"nonsplit":None},ambient[D],cutoff)
            result["finite_upper_bound_through_degree"] = cutoff
        else:
            continue
        results.append(result)
        print(D,ell,file=sys.stderr,flush=True)
    rendered = json.dumps(results,indent=2)
    if args.verify and json.loads(rendered)!=json.loads(args.verify.read_text()):
        raise ArithmeticError("The pair-projector replay differs from the saved receipt.")
    print(rendered)


if __name__ == "__main__":
    main()
