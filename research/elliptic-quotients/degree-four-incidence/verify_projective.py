"""Exact certificates for a quartic presentation of a map to a short Weierstrass cubic."""
import json
import sys
from pathlib import Path

import sympy as sp

from incidence import axis_unit, identities, source, evaluate

u, v, w = sp.symbols("u v w")


def monomials(degree):
    return [u**i * v**j * w**(degree-i-j)
            for i in range(degree+1) for j in range(degree-i+1)]


def homogeneous(p):
    return sp.Poly(p, u, v, domain=sp.QQ).homogenize(w).as_expr() * w**(4-sp.Poly(p, u, v).total_degree())


def compact_presentation(F, x, y):
    """Lift functions already specified on w=1 to three forms of degree four."""
    divisor = sp.Poly(F.subs(w, 1), u, v, domain=sp.QQ)
    def reduce(p):
        return sp.Poly(p, u, v, domain=sp.QQ).rem(divisor)
    basis = [m.subs(w, 1) for m in monomials(4) if m.subs(w, 1) != divisor.LM().as_expr()]
    x, y = reduce(x).as_expr(), reduce(y).as_expr()
    multiples = [[reduce(m * t) for m in basis] for t in (x, y)]
    equations = []
    for table in multiples:
        high = sorted({powers for p in table for powers in p.monoms() if sum(powers) > 4})
        equations.extend([[p.coeff_monomial(powers) for p in table] for powers in high])
    M = sp.Matrix(equations) if equations else sp.zeros(0, len(basis))
    columns = M.nullspace()
    if not columns:
        raise ValueError("No polynomial denominator in this degree-four affine chart")
    coeff = columns[0]
    c = sum(t*m for t, m in zip(coeff, basis))
    return [homogeneous(reduce(c*t).as_expr()).expand() for t in (x, y, 1)]


def verify(F, forms, a, b, field=sp.QQ):
    """Verify smooth source/target, cubic identity, map degree, and differential direction."""
    if not (4*a**3 + 27*b**2):
        raise ValueError("Singular target")
    for t in (u, v, w):
        other = [s for s in (u, v, w) if s != t]
        partials = [F.diff(s).subs(t, 1) for s in (u, v, w)]
        if list(sp.groebner(partials, *other, domain=field)) != [1]:
            raise ValueError("Singular source")
    if any(not sp.Poly(h, u, v, w).is_homogeneous or sp.Poly(h,u,v,w).total_degree()!=4
           for h in (F,*forms)):
        raise ValueError("Expected three quartic forms on a plane quartic")
    A, B, C = [p.subs(w, 1) for p in forms]
    f = sp.Poly(F.subs(w, 1), u, v, domain=field)
    def red(p):
        return sp.Poly(p, u, v, domain=field).rem(f)
    if red(B**2*C-A**3-a*A*C**2-b*C**3) or not red(B*C):
        raise ValueError("Invalid cubic identity or image contained in the coordinate axes")
    def D(p):
        return f.as_expr().diff(v)*p.diff(u) - f.as_expr().diff(u)*p.diff(v)
    derivative = red(C*D(A)-A*D(C))
    if not derivative:
        raise ValueError("Constant map")
    coeffs = sp.symbols("l0:3")
    basis = [red(2*B*C*m) for m in (u,v,1)]
    support = sorted(set(derivative.monoms()).union(*(p.monoms() for p in basis)))
    eq = [sum(c*p.coeff_monomial(m) for c,p in zip(coeffs,basis)) - derivative.coeff_monomial(m)
          for m in support]
    line, = sp.linsolve(eq, coeffs)
    scale = next(t for t in line if t)
    normalized = [str(sp.simplify(t/scale)) for t in line]
    terms = [sp.Poly(m, u, v, w).monoms()[0] for m in monomials(10)]
    rows = []
    for p in (F, *forms):
        for m in monomials(6):
            data = sp.Poly(p*m, u, v, w, domain=field).as_dict()
            rows.append([data.get(t, 0) for t in terms])
    rank = sp.polys.matrices.DomainMatrix.from_Matrix(sp.Matrix(rows), extension=True).rank()
    base = len(terms) - rank
    if (16-base) % 3:
        raise ArithmeticError("Base length does not yield integral degree")
    return {"rank":rank, "base_length":base, "degree":(16-base)//3,
            "differential":normalized, "differential_scale":str(scale)}


def fixture():
    F = -3*u**4+u**2*v**2-2*u**2*w**2+2*v**4+17*v**2*w**2+8*w**4
    r = (3*u**2+2*v**2+11)/5
    m = r*(u-v)
    xf = m**2+3-2*r
    yf = m*(r-xf)-r*v
    inverse_r = (v**2-u**2+3)/5
    a, b = -sp.Rational(7,3), -sp.Rational(20,27)
    rational_sum = compact_presentation(F,xf-sp.Rational(4,3),yf)
    double = compact_presentation(F,r-sp.Rational(4,3),r*v)
    composite = compact_presentation(F,v**2+sp.Rational(8,3),v*(r-3*inverse_r))
    k = sp.sqrt(-2)
    # Quartic coordinates for P+S over Q(k), obtained by clearing the
    # denominators in the chord addition law and solving for the three forms.
    conjugate_sum = [
        k*u*w*(1728*u**2-9558*v**2+5376*w**2) + v**2*(3442*u**2+5888*v**2)
        + w**2*(12862*u**2+39407*v**2+20504*w**2),
        3*v*(k*u*(2106*u**2+2559*v**2+4902*w**2)-w*(6378*u**2+9117*v**2+11376*w**2)),
        3*(k*u*w*(10080*v**2-1920*w**2) + u**2*(4031*v**2-4390*w**2)
           + 634*v**4-15893*v**2*w**2-3560*w**4),
    ]
    models = [
        ("sum",rational_sum,a,b,4,[1,1,0],False),
        ("rational_sign",[h.subs(v,-v) for h in rational_sum],a,b,4,[1,-1,0],False),
        ("quadratic_sum",conjugate_sum,a,b,4,[1,0,-k],False),
        ("quadratic_conjugate",[h.xreplace({sp.I:-sp.I}) for h in conjugate_sum],a,b,4,[1,0,k],False),
        ("isogeny_after_double_cover",composite,-sp.Rational(52,3),sp.Rational(736,27),
         4,[1,0,0],True),
        ("double_cover",double,a,b,2,[1,0,0],True),
    ]
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "involution-atlas"))
    from classify import classify
    fdict = sp.Poly(F,u,v,w).as_dict()
    atlas = classify({"terms":[[*m,int(c)] for m,c in fdict.items()]})
    assert atlas["geometric_involutions"] == 5
    receipts = {}
    for label, forms, a, b, deg, coords, bielliptic in models:
        field = sp.QQ.algebraic_field(k) if coords[-1] != 0 else sp.QQ
        result = verify(F,forms,a,b,field)
        axis = list(map(str,coords))
        assert result["degree"] == deg and result["differential"] == axis, (label,result)
        forms_dict = [sp.Poly(h,u,v,w).as_dict() for h in forms]
        values = identities(fdict,forms_dict,a,b,coords,sp.sympify(result["differential_scale"]))
        assert not any(sp.expand(value) for value in values)
        sep = sp.Symbol("s")
        matches = [o for o in atlas["orbits"] if axis_unit(o,coords,sep) == 0]
        assert bool(matches) == bielliptic
        result["is_bielliptic_axis"] = bool(matches)
        result.update({"forms":[str(h) for h in forms], "a":str(a), "b":str(b)})
        receipts[label] = result
    return receipts


def check_rational_input(source_path, map_path):
    """Read a source and a map of quartic forms given by rational monomial arrays."""
    source_data = json.loads(source_path.read_text())
    F = evaluate(source(source_data),(u,v,w))
    model = json.loads(map_path.read_text())
    forms = []
    for terms in model["forms"]:
        form = 0
        for i,j,k,c in terms:
            if any(type(e) != int or e < 0 for e in (i,j,k)) or i+j+k != 4:
                raise ValueError("Each map monomial must have degree four")
            form += sp.Rational(c)*u**i*v**j*w**k
        forms.append(form)
    result = verify(F,forms,sp.Rational(model["a"]),sp.Rational(model["b"]))
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "involution-atlas"))
    from classify import classify
    coords = list(map(sp.Rational,result["differential"]))
    result["is_bielliptic_axis"] = any(
        axis_unit(orbit,coords,sp.Symbol("s")) == 0 for orbit in classify(source_data)["orbits"])
    return result


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", nargs=2, type=Path, metavar=("SOURCE","MAP"),
                        help="check a rational model with three quartic monomial arrays")
    args = parser.parse_args()
    result = check_rational_input(*args.check) if args.check else fixture()
    output = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(output)
    else:
        print(output, end="")
