"""Archived ball-enclosure experiment for explicit bielliptic models.

Use certify_link.py at D=-19 and verify_d43.py for p3 at D=-43.
The candidate-ideal overlap stage here assumes that an exact j-value
belongs to the supplied list. The maintained replays prove that premise.

Usage: python certify_link_any.py D [number_of_embeddings].  D = -19 reproduces certify_link.py; other D read equations_<D>.json.

For C: y^2 = x^6 - s1 x^4 + s2 x^2 - 1 (s1, s2 in L = Q(eta), eta^24 + eta^18 + 27 eta^12 + 163 eta^6 + 64 = 0):
  E1: Y^2 = g1(U) = U^3 - s1 U^2 + s2 U - 1,   E2: Y^2 = g2(U) = U^3 - s2 U^2 + s1 U - 1 (up to the twist by -1),
  roots of g2 are the inverses of the roots of g1, and (Kani; pullbacks of the two degree-2 maps)
        Jac(C) = (E1 x E2) / graph(psi),   psi(rho_i, 0) = (1/rho_i, 0),   pulled-back polarisation 2(th1 + th2).
Certified steps (Arb, 1200 bits, at a fixed embedding eta0 isolated by acb_poly roots):
  (a) j(E_k) is enclosed; among the proper ideal classes of the orders O_1, O_3 of Q(sqrt -19) exactly one has
      j(tau_a) (acb modular_j) overlapping it, all others disjoint  ->  E_k = C / lambda_k a_k;
  (b) the half-period values e(h) = P(h; tau_a) at h = 1/2, tau/2, (1+tau)/2 (acb elliptic_p) scaled by
      lambda^-2 match the shifted roots under exactly one bijection (all other bijections excluded)  ->  the
      2-torsion labels;  roots of g2 are matched to 1/rho_i by enclosure;
  (c) Lambda_J = a1 + a2 + Z (h1(i), h2(psi i)), theta = 2(H1 + H2): exact NS and refined Humbert form (rational
      arithmetic) and exact isometry test (qfisom) against the twin forms.
"""
import itertools, json, sys
sys.path.insert(0, "."); sys.path.insert(0, "../inverse-spectrum")
from fractions import Fraction as Fr
from flint import acb, arb, acb_poly, ctx
import cypari2
from polarised import ideals, ns_general
from gluing import K, E_of, q_from_ns, reduce_form, isometric

pari = cypari2.Pari()
ctx.prec = 1200
d = int(sys.argv[1]) if len(sys.argv) > 1 else -19
Kf = K(d)
wC = (acb(d) + acb(d).sqrt()) / 2                                                       # w = (d + sqrt d)/2, Im > 0


def disjoint(a, b):
    return not a.overlaps(b)


def eval_L(expr, eta):
    P = pari(expr)
    deg = int(pari.poldegree(P))
    s = acb(0)
    for k in range(deg + 1):
        c = pari.polcoef(P, k)
        if c != 0:
            s += acb(arb(int(pari.numerator(c))) / int(pari.denominator(c))) * eta ** k
    return s


def kC(el):
    p, q = el
    return acb(arb(p.numerator) / p.denominator) + acb(arb(q.numerator) / q.denominator) * wC


def classify(g2, g3, cands):
    jE = 1728 * g2 ** 3 / (g2 ** 3 - 27 * g3 ** 2)
    hits = []
    for idx, (f, a) in enumerate(cands):
        al, be = kC(a[0]), kC(a[1])
        tau = be / al
        if tau.imag < 0:
            tau = -tau
            a = (a[0], (-a[1][0], -a[1][1]))
        j = tau.modular_j()
        if jE.overlaps(j):
            hits.append((idx, f, a, tau))
        else:
            assert disjoint(jE, j)
    assert len(hits) == 1, f"j not certified uniquely ({len(hits)} candidates)"
    return hits[0], jE


def weierstrass(a, b, c):
    """Y^2 = U^3 + a U^2 + b U + c  ->  x = U + a/3,  4 g = 4x^3 - g2 x - g3."""
    g2 = -4 * (b - a * a / 3)
    g3 = -4 * (2 * a ** 3 / 27 - a * b / 3 + c)
    return g2, g3


def half_period_labels(roots_U, a, g2, g3, tau):
    """Bijection root index -> half period index (0: 1/2, 1: tau/2, 2: (1+tau)/2), certified."""
    hp = [acb(1) / 2, tau / 2, (1 + tau) / 2]
    eT = [acb.elliptic_p(h, tau) for h in hp]
    g2T, g3T = acb.elliptic_invariants(tau)
    s = (g3 / g3T) / (g2 / g2T)                         # lambda^-2
    eE = [r + a / 3 for r in roots_U]
    good = []
    for perm in itertools.permutations(range(3)):
        if all(eE[i].overlaps(s * eT[perm[i]]) for i in range(3)):
            good.append(perm)
        else:
            assert any(disjoint(eE[i], s * eT[perm[i]]) for i in range(3))
    assert len(good) == 1, f"half-period labelling not unique ({len(good)})"
    return good[0]


def certify(name, s1x, s2x, eta, cands):
    s1, s2 = eval_L(s1x, eta), eval_L(s2x, eta)
    res = {}
    # E1
    a1_, b1_, c1_ = -s1, s2, acb(-1)
    r1 = acb_poly([c1_, b1_, a1_, 1]).roots(tol=arb(2) ** -1100, maxprec=4 * ctx.prec)
    g21, g31 = weierstrass(a1_, b1_, c1_)
    (i1, f1, A1, tau1), j1 = classify(g21, g31, cands)
    p1 = half_period_labels(r1, a1_, g21, g31, tau1)
    # E2 (roots must be the inverses of r1)
    a2_, b2_, c2_ = -s2, s1, acb(-1)
    r2 = acb_poly([c2_, b2_, a2_, 1]).roots(tol=arb(2) ** -1100, maxprec=4 * ctx.prec)
    g22, g32 = weierstrass(a2_, b2_, c2_)
    (i2, f2, A2, tau2), j2 = classify(g22, g32, cands)
    p2 = half_period_labels(r2, a2_, g22, g32, tau2)
    inv = []
    for i in range(3):
        m = [k for k in range(3) if (1 / r1[i]).overlaps(r2[k])]
        assert len(m) == 1 and all(disjoint(1 / r1[i], r2[k]) for k in range(3) if k != m[0])
        inv.append(m[0])
    # lattice of J (exact, rational coordinates over (1, w))
    def hp_vec(A, k):
        al, be = A
        v = [al, be, (al[0] + be[0], al[1] + be[1])][k]
        return (v[0] / 2, v[1] / 2)
    gens = [(A1[0][0], A1[0][1], Fr(0), Fr(0)), (A1[1][0], A1[1][1], Fr(0), Fr(0)),
            (Fr(0), Fr(0), A2[0][0], A2[0][1]), (Fr(0), Fr(0), A2[1][0], A2[1][1])]
    for i in range(3):
        h1, h2 = hp_vec(A1, p1[i]), hp_vec(A2, p2[inv[i]])
        gens.append((h1[0], h1[1], h2[0], h2[1]))
    den = 1
    for g in gens:
        for x in g:
            den = den * x.denominator // __import__("math").gcd(den, x.denominator)
    M = pari.matrix(4, len(gens), [int(gens[j][i] * den) for i in range(4) for j in range(len(gens))])
    H = pari.mathnf(M)
    B = [tuple(Fr(int(H[i, j]), den) for i in range(4)) for j in range(4)]
    # principal polarisation on C/a: H = 1/covol, covol = |w-coeff of al * conj(be)|
    def princ(A):
        al, be = A
        cov = abs(Kf.mul(al, Kf.conj(be))[1])
        return Fr(1) / cov
    theta = (2 * princ(A1), 2 * princ(A2), Fr(0), Fr(0))
    Et = [[E_of(Kf, theta, B[i], B[j]) for j in range(4)] for i in range(4)]
    assert all(x.denominator == 1 for r in Et for x in r), "theta not integral on Lambda_J"
    Em = pari.matrix(4, 4, [int(x) for r in Et for x in r])
    assert abs(int(pari.matdet(Em))) == 1, "theta not principal on Lambda_J"
    ns, G = ns_general(Kf, B)
    th = pari.matsolve(pari.matrix(4, 4, [str(ns[j][i]) for i in range(4) for j in range(4)]),
                       pari.matrix(4, 1, [str(x) for x in theta]))
    tc = [Fr(str(th[i, 0])) for i in range(4)]
    assert all(x.denominator == 1 for x in tc)
    q = reduce_form(q_from_ns(G, [int(x) for x in tc]))
    res = {"E1_class": {"conductor": f1, "basis": [[str(x) for x in e] for e in A1]},
           "E2_class": {"conductor": f2, "basis": [[str(x) for x in e] for e in A2]},
           "half_period_labels_E1": list(p1), "half_period_labels_E2": list(p2), "psi_root_matching": inv,
           "humbert_form": q}
    return res


if __name__ == "__main__":
    if d == -19:
        E = json.load(open("equations_-19_reduced.json")); fieldpol = "x^24+x^18+27*x^12+163*x^6+64"; names = ["q1", "q2"]
    else:
        E = json.load(open(f"equations_{d}.json")); names = sorted(E)
        fields = {E[n]["model_field"] for n in names}
    T = json.load(open(f"twin_lattices_{d}.json"))
    cands = [(f, a) for f in (1, 3) for a in ideals(d, f)]
    # ideals() already lists every reduced form (a, b, c) with -a < b <= a: all classes, conjugates included
    nemb = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    out = {}
    for name in names:
        fieldpol = fieldpol if d == -19 else E[name]["model_field"]
        FP = pari(fieldpol)
        LPOL = [int(pari.polcoef(FP, k)) for k in range(int(pari.poldegree(FP)) + 1)]
        roots = acb_poly([acb(c) for c in LPOL]).roots(tol=arb(2) ** -1100, maxprec=4 * ctx.prec)
        for k, eta in enumerate(roots[:nemb]):
            r = certify(name, E[name]["s1"], E[name]["s2"], eta, cands)
            iso = {t: bool(isometric(r["humbert_form"], T[t]["q"])) for t in sorted(T)}
            out.setdefault(f"embedding_{k}", {})[name] = {**r, "isometric_to_twin": iso}
            print(f"curve {name} embedding {k}: E1 conductor {r['E1_class']['conductor']}, E2 conductor {r['E2_class']['conductor']}, "
                  f"Humbert form {r['humbert_form']}  isometric to: {[t for t, b in iso.items() if b]}", flush=True)
    json.dump(out, open(f"certified_link_any_{d}.json", "w"), indent=1)
