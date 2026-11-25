"""Frobenius on J[4] at p = 37 for C and C': same characteristic polynomial, different modules.

Invariant reported: the Z/4-rank of the image of (F+V-b) on J[4], i.e. its size, which is
an isomorphism invariant of the Z/4[F]-module and is computed from the matrix of F.
"""
import json
from certify import Certifier

C = [-36, 0, 49, 0, -14, 0, 1]
Cp = [-3003, -14850, 11968, 6600, -3575, -990, 154]
out = {}
for name, f in [("C", C), ("C'", Cp)]:
    c = Certifier(f, 37, rmax=6, sylow_limit=2 ** 16)
    r, M = c.frobenius_mod(2, 2)
    p = 37
    # V = p F^{-1}; on J[4] with F invertible mod 4 compute F^{-1} via adjugate
    import sympy as sp
    MF = sp.Matrix(M)
    Finv = MF.inv_mod(4)
    a, b = 2, -6
    phi = (MF + p * Finv - b * sp.eye(4)).applyfunc(lambda t: t % 4)
    # size of image of phi on (Z/4)^4 via Smith form over Z with modulus
    from sympy.matrices.normalforms import smith_normal_form
    S = smith_normal_form(sp.Matrix.vstack(phi, 4 * sp.eye(4)), domain=sp.ZZ)
    diag = [abs(S[i, i]) for i in range(4)]
    image = 1
    for dd in diag:
        image *= 4 // sp.gcd(dd, 4)
    charp = sp.Poly((MF.charpoly().as_expr()), sp.Symbol("lambda")).all_coeffs()
    out[name] = {"r": r, "F_mod_4": M, "charpoly_mod_4": [int(t) % 4 for t in charp],
                 "phi_mod_4": [[int(t) for t in row] for row in phi.tolist()], "image_of_phi_on_J4": int(image)}
    print(name, out[name])
json.dump(out, open("modules.json", "w"), indent=1)
