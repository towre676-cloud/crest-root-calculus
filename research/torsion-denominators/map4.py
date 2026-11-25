"""Explicit map C' -> E1 from the Richelot correspondence, verified exactly.

C : y^2 = (x^2-1)(x^2-4)(x^2-9),   E1 : Y^2 = (u-1)(u-4)(u-9),   C -> E1 : (x,y) -> (x^2, y).
C': w^2 = 11 (z^2-10z+13)(2z^2+2z-7)(7z^2+18z+3)  =  121 H1 H2 H3 / Delta.
For a point (z,w) of C' the Richelot correspondence gives two points (x_i, y_i) of C with
G1(x)H1(z) + G2(x)H2(z) = 0 and y_i v = G1(x_i) H1(z) (x_i - z), where w = 11 v.
The map is (z,w) -> (x1^2, y1) + (x2^2, y2) on E1.
"""
import json

import sympy as sp
from sympy.polys.polyfuncs import symmetrize

from richelot import richelot, x, z

x1, x2, e1, e2, v = sp.symbols("x1 x2 e1 e2 v")
G = [sp.expand((x - 1) * (x - 2)), sp.expand((x + 1) * (x - 3)), sp.expand((x + 2) * (x + 3))]
Delta, H = richelot(G)
fC = sp.expand((x**2 - 1) * (x**2 - 4) * (x**2 - 9))
cubic = lambda u: (u - 1) * (u - 4) * (u - 9)
R = sp.expand(H[0] * H[1] * H[2] / Delta).subs(x, z)       # v^2 = R(z)
Rp = sp.expand(11 * (z**2 - 10*z + 13) * (2*z**2 + 2*z - 7) * (7*z**2 + 18*z + 3))
assert sp.expand(121 * R - Rp) == 0

quad = sp.Poly(sp.expand(G[0] * H[0].subs(x, z) + G[1] * H[1].subs(x, z)), x)
a, b, c = quad.all_coeffs()
E1v, E2v = -b / a, c / a

def sym(expr):
    s, rem, *_ = symmetrize(sp.expand(expr), [x1, x2], formal=True)
    sub = dict((sym_, val) for sym_, val in rem) if rem else {}
    return s

def to_z(expr):
    s, rem = symmetrize(sp.expand(expr), x1, x2, formal=True)[:2]
    assert rem == 0
    return s

def S(expr):
    out, rem, defs = symmetrize(sp.expand(expr), x1, x2, formal=True)
    assert rem == 0
    return sp.cancel(out.subs({d[0]: (E1v if i == 0 else E2v) for i, d in enumerate(defs)}))

H1z = H[0].subs(x, z)
g1 = lambda t: G[0].subs(x, t)
# y_i = g1(x_i) H1(z) (x_i - z) / v ;  y1 y2 = P12 / v^2, with v^2 = R
y1y2 = S(g1(x1) * g1(x2) * (x1 - z) * (x2 - z)) * H1z**2 / R
ysq = S(fC.subs(x, x1) + fC.subs(x, x2))                  # y1^2 + y2^2 = f(x1) + f(x2)
du2 = S((x1**2 - x2**2) ** 2)
usum = S(x1**2 + x2**2)
lam2 = (ysq - 2 * y1y2) / du2
u3 = sp.factor(sp.cancel(lam2 + 14 - usum))
num, den = sp.fraction(u3)
deg = max(sp.degree(num, z), sp.degree(den, z))
print("u3(z) =", u3)
print("degree of u3 in z:", deg)

# Y3 = -(lambda (u3 - u1) + y1), symmetric up to the factor 1/v:
# lambda (u3-u1) + y1 = [(y2-y1)(u3-u1)/(u2-u1)] + y1; multiply numerator and denominator by (u2-u1)
# and symmetrize the combination (y2 - y1)(u3 - u1)(u2 - u1)... it is symmetric in (1,2):
yfac = lambda t: g1(t) * (t - z) * H1z  # y_i = yfac(x_i)/v
num_sym = S((yfac(x2) - yfac(x1)) * (u3 - x1**2) * (x2**2 - x1**2) + yfac(x1) * (x2**2 - x1**2) ** 2)
Y3_times_v = sp.factor(sp.cancel(-num_sym / du2))          # Y3 = Y3_times_v / v
chk = sp.cancel(Y3_times_v**2 / R - cubic(u3))
print("Y3^2 - cubic(u3) reduces to:", chk)

# In the coordinates of C': v = w/11
Y3_w = sp.factor(sp.cancel(Y3_times_v * 11))               # Y3 = Y3_w / w
out = {"u": str(u3), "Y_times_w": str(Y3_w), "degree_in_z": int(deg), "identity_check": str(chk)}
json.dump(out, open("map4.json", "w"), indent=1)
