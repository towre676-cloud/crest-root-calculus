"""Generators of E_theta(L), 2-saturation over L, and a Mordell-Weil search for u = x/theta in Q(i)."""
import itertools, json, sys, time
from cypari2 import Pari
pari = Pari(); pari.allocatemem(4*10**9)
pari('zz')
POL = pari('w^4-w^2+1'); L = pari.bnfinit(POL, 1)
M = lambda a: pari.Mod(pari(a), POL)
I, R3, Z6 = M('w^3'), M('2*w^2-1'), M('w^2')
V = Z6 - 1                                          # zeta3 = zeta6 - 1 = w^2 - 1
tau = lambda a: pari.Mod(pari.subst(pari.lift(M(a)), 'w', pari('w^5')), POL)
def issq(a):
    a = pari.lift(M(a)); return a == 0 or len(pari.nfroots(L, pari('zz^2') - a)) > 0

DATA = {  # theta, point on E_eta over K, point on E_{-eta} over K (from points_K.py), K-variable v = zeta3
    'theta1': (3*Z6 - 1, (-V + 1, -3*V + 6), (3*V - 3, -9*V)),
    'theta2': (3*Z6 - 2, (V + 2, -3*V - 9), (2*V - 11, -9*V - 24)),
}

def setup(name):
    theta, (x1, y1), (X2, Y2) = DATA[name]
    eta = theta*R3
    E = pari.ellinit([0, 4*eta, 0, 3*eta**2, 0], L)
    G1 = [x1, y1]; G2 = [-X2, I*Y2]                 # (X,Y) on E_{-eta} -> (-X, iY) on E_eta
    assert pari.ellisoncurve(E, G1) and pari.ellisoncurve(E, G2)
    es = [M(0), -eta, -3*eta]
    T = [[e, M(0)] for e in es[1:]]
    return theta, eta, E, es, G1, G2, T

def kummer(P, es):
    if len(P) == 1: return (M(1), M(1))
    x = P[0]; out = []
    for e in es[:2]:
        if x == e:
            o = [f for f in es if f != e]; out.append((e - o[0])*(e - o[1]))
        else:
            out.append(x - e)
    return tuple(out)

def two_saturated(E, es, gens):
    for bits in itertools.product((0, 1), repeat=len(gens)):
        if not any(bits): continue
        P = [0]
        for b, g in zip(bits, gens):
            if b: P = pari.elladd(E, P, g)
        k = kummer(P, es)
        if issq(k[0]) and issq(k[1]):
            return False, bits
    return True, None

if __name__ == '__main__':
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    report = {}
    for name in DATA:
        theta, eta, E, es, G1, G2, T = setup(name)
        ok, bad = two_saturated(E, es, [G1, G2] + T)
        H = [float(pari.ellheight(E, G)) for G in (G1, G2)]
        print(name, 'eta', pari.lift(eta), '2-saturated:', ok, 'heights', H, flush=True)
        t0 = time.time()
        hits = {}
        multiples1 = {a: pari.ellmul(E, G1, a) for a in range(-B, B+1)}
        multiples2 = {b: pari.ellmul(E, G2, b) for b in range(-B, B+1)}
        tors = [[0]] + T + [pari.elladd(E, T[0], T[1])]
        for a, b in itertools.product(range(-B, B+1), repeat=2):
            base = pari.elladd(E, multiples1[a], multiples2[b])
            for t in tors:
                P = pari.elladd(E, base, t)
                if len(P) == 1:
                    u = 'oo'
                else:
                    u = M(P[0])/theta
                    if tau(u) != u: continue          # u not in Q(i)
                    u = str(pari.lift(u))
                hits.setdefault(u, []).append((a, b, tors.index(t)))
        print('  u in Q(i) found:', {k: len(v) for k, v in hits.items()}, f'{time.time()-t0:.0f}s', flush=True)
        report[name] = {'two_saturated': ok, 'u_values': {k: v[:4] for k, v in hits.items()}}
    json.dump(report, open('mw_search.json', 'w'), indent=1)
