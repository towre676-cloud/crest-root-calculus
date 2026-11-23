"""Independently verify the descended sextic equations in each exact number field."""
import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

from frobenius import X
from kernel_fields import Y
from partitions import sextic_chart


def det(a):
    return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
            -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
            +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))


def check_all(source, data):
    f,chart=sextic_chart(source);a=int(f.LC())
    if data['source_chart']!=list(chart) or data['monic_branch_multiplier']!=a:
        raise ArithmeticError('A different branch chart was used')
    integral=sp.Poly.from_list([f.nth(i)*sp.Integer(a)**(5-i) for i in range(6,-1,-1)],X,domain=sp.ZZ)
    k=data['pair_parameter']; b=data['matching_evaluation']
    fields=data['kernel_field_polynomials_descending']; models=data['orbit_models']
    if len(models)!=len(fields) or sum(len(h)-1 for h in fields)!=15:
        raise ArithmeticError('Incomplete kernel-field support')
    for raw,model in zip(fields,models):
        F=sp.Poly.from_list(raw,Y,domain=QQ)
        if not F.is_irreducible:
            raise ArithmeticError('Kernel field polynomial is reducible')
        if F.degree()==1:
            K=QQ
            u=K.from_sympy(-F.nth(0)/F.nth(1))
        else:
            K=QQ.alg_field_from_poly(F,'u')
            u=K.unit
        def val(vals):
            v=K.zero
            for c in reversed(vals):
                v=v*u+K.from_sympy(sp.Rational(c))
            return v
        A,t=ring('theta',K)
        def poly(c):
            return A.from_dict({(i,):v for i,v in enumerate(c) if v})
        coeffs=[val(c) for c in model['pair_cubic_lower_coefficients']]
        m=t**3+poly(coeffs)
        if m.evaluate(t,b)!=u:
            raise ArithmeticError('Matching-field label does not identify the cubic')
        inverse,_,one=m.diff(t).gcdex(m)
        if one!=A.one:
            raise ArithmeticError('The cubic pairing algebra is not etale')
        prod=(poly([val(c) for c in model['pair_product_numerator']])*inverse)%m
        total=(t-k*prod)%m
        basis=[A.one,-total,prod]
        def coords(v):
            return [v.get((i,),K.zero) for i in range(3)]
        def norm(v):
            return det([coords((v*t**j)%m) for j in range(3)])
        def trace(v):
            return sum((coords((v*t**j)%m)[j] for j in range(3)),K.zero)
        B=[[trace((v*t**j)%m) for j in range(3)] for v in basis]
        determinant=det(B)
        delta2=det([coords(v) for v in basis])**2*det([[trace((t**(i+j))%m) for j in range(3)] for i in range(3)])
        if bool(delta2)==model['product_quotient']:
            raise ArithmeticError('Product classification disagrees with the trace form')
        output=[val(c) for c in model['target_sextic_ascending']]
        for x in range(7):
            if norm(A(x*x)-x*total+prod)!=K.convert(int(integral.eval(x))):
                raise ArithmeticError('The cubic norm does not recover the source')
            if not delta2:
                if output:
                    raise ArithmeticError('Unexpected genus-two model on a product quotient')
                continue
            rhs=[K.one,K.convert(-2*x),K.convert(x*x)]
            cols=[]
            for j in range(3):
                subs=[[rhs[i] if col==j else B[i][col] for col in range(3)] for i in range(3)]
                cols.append(det(subs)/determinant)
            observed=sum((v*x**i for i,v in enumerate(output)),K.zero)
            if observed!=-a*delta2*norm(poly(cols)):
                raise ArithmeticError('The target fails the trace-dual norm identity')
    return data['orbit_degrees']
