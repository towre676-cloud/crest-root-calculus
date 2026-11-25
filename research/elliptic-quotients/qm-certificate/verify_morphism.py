"""Exact verification that (s1, s0, g1, g0) over L define a morphism psi : C -> Sym^2 C, and its tangent action.

Reduction modulo u = T^2 - s1 T + s0 uses T^k = a_k T + b_k with a_0 = 0, b_0 = 1,
a_{k+1} = s1 a_k + b_k, b_{k+1} = -s0 a_k, so every check is an identity in L(x).

Membership: f(x)(g1 T + g0)^2 - f(T) = 0 mod u.
Tangent:    psi^*(sum_i T_i^k dT_i / V_i) = H1 dx/y, where dT_i = (T_i ds1 - ds0)/u'(T_i),
            V_i = y (g1 T_i + g0), and H1 T + H0 = T^k (T s1' - s0') (g1 T + g0)^{-1} mod u
            (sum_i h(T_i)/u'(T_i) = H1 for h = H1 T + H0).  H1 must be c0 + c1 x.
By the Albanese property psi - psi(W0) extends to an endomorphism, whose action on the invariant
differentials is the matrix of the c's.
"""
import cypari2

pari = cypari2.Pari()
pari.allocatemem(128 * 10 ** 6)


def verify(fcoeffs, s1, s0, g1, g0, nf_poly):
    mod = f"Mod(t, {nf_poly})"
    rd = lambda s: pari(s.replace("t", mod))
    S1, S0, G1, G0 = rd(s1), rd(s0), rd(g1), rd(g0)
    if S1*S1 == 4*S0:
        raise ValueError("The generic degree-two divisor has repeated support.")
    x = pari("x")
    fx = sum((c * x ** i for i, c in enumerate(fcoeffs)), pari(0))
    deg = len(fcoeffs) - 1
    a, b = [pari(0)], [pari(1)]
    for k in range(max(deg, 3) + 1):
        a.append(S1 * a[-1] + b[-1])
        b.append(-S0 * a[-2])
    fa = sum((c * a[k] for k, c in enumerate(fcoeffs)), pari(0))
    fb = sum((c * b[k] for k, c in enumerate(fcoeffs)), pari(0))
    C1 = fx * (G1 * G1 * S1 + 2 * G1 * G0) - fa
    C0 = fx * (G0 * G0 - G1 * G1 * S0) - fb
    member = pari.simplify(pari.lift(C1)) == 0 and pari.simplify(pari.lift(C0)) == 0
    det = (G1 * S1 + G0) * G0 + G1 * G1 * S0
    if not det:
        raise ValueError("The vertical coefficient is not invertible in the generic fibre.")
    A = -G1 / det
    B = (G1 * S1 + G0) / det
    dS1, dS0 = pari.deriv(S1, "x"), pari.deriv(S0, "x")
    tangent = []
    linear = True
    for k in (0, 1):
        H1 = A * dS1 * a[k + 2] + (B * dS1 - A * dS0) * a[k + 1] - B * dS0 * a[k]
        H1 = pari.simplify(H1)
        num, den = pari.numerator(H1), pari.denominator(H1)
        if pari.poldegree(den, "x") > 0 or pari.poldegree(num, "x") > 1:
            linear = False
            tangent.append(str(H1)[:200])
        else:
            tangent.append([str(pari.polcoef(H1, 0, "x")), str(pari.polcoef(H1, 1, "x"))])
    return {"member": bool(member), "tangent_linear": linear, "tangent": tangent}
