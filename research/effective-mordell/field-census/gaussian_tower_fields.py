"""Recognize the proved towers from a monic defining polynomial and return exact points."""
from argparse import ArgumentParser
import json
from cypari2 import Pari

pari = Pari()
x, w, g = pari('x'), pari('w'), pari('g')
gaussian_equation = g**2 + 1
base = pari.bnfinit(gaussian_equation, 1)
i_base = pari.Mod(g, gaussian_equation)


def require(statement, message):
    if not statement:
        raise ArithmeticError(message)


def squarefree_signed(n):
    factor = pari.factorint(abs(n))
    result = -1 if n < 0 else 1
    for p, e in zip(*factor):
        require(pari.isprime(p), 'squarefree decomposition used an unproved prime')
        if int(e) % 2:
            result *= int(p)
    return result


def quadratic_subfields(f):
    if f.poldegree() == 1:
        return []
    return sorted({squarefree_signed(int(pari.poldisc(row[0]))) for row in pari.nfsubfields(f, 2)})


def rational_tower(characters):
    """Membership in the rational multiquadratic towers from squareclasses."""
    eligible, plus, minus, old = True, True, True, True
    for d in characters:
        for p in pari.factorint(abs(d))[0]:
            require(pari.isprime(p), 'unproved support prime')
            eligible &= int(p) in (2, 3) or int(p) % 12 == 5
        ev = int(d % 2 == 0)
        odd = d // (2**ev)
        old &= not ev
        plus &= odd % 4 == 1
        minus &= (odd % 4 == 3) == ev
    if eligible:
        if old:
            return 'Gaussian radical tower'
        if plus:
            return 'positive even tower'
        if minus:
            return 'negative even tower'
    return None


def gaussian_test(delta):
    """Return (allowed quartic support, in the large dyadic-kernel tower)."""
    primes, exps = pari.idealfactor(base, delta)
    generator = pari.Mod(1, gaussian_equation)
    support, kernel = True, True
    for prime, exp in zip(primes, exps):
        if not int(exp) % 2:
            continue
        norm = int(pari.idealnorm(base, prime))
        if norm in (2,9):
            if norm == 9:
                generator *= 3
            else:
                kernel = False
            continue
        if norm % 3 != 2:
            return False, False
        exponents, vector = pari.bnfisprincipal(base, prime)
        require(not any(exponents), 'a Gaussian ideal was nonprincipal')
        pi = pari.Mod(pari.nfbasistoalg(base, vector), gaussian_equation)
        if int(pari.polcoef(pari.lift(pi), 0, g)) % 2 == 0:
            pi *= i_base
        require(pari.idealhnf(base, pi) == pari.idealhnf(base, prime), 'Gaussian prime generator')
        generator *= pi
    kernel &= int(pari.polcoef(pari.lift(generator), 1, g)) % 4 == 0
    if kernel:
        kernel = bool(len(pari.nfroots(base, x*x-delta/generator)))
    return support, kernel


def gaussian_coordinates(alpha, i, degree):
    """Read the two Gaussian coefficients in the ambient absolute power basis."""
    A, I = pari.lift(alpha), pari.lift(i)
    s = next(j for j in range(1,degree) if pari.polcoef(I, j, w))
    imag = pari.polcoef(A,s,w) / pari.polcoef(I,s,w)
    real = pari.polcoef(A,0,w) - imag*pari.polcoef(I,0,w)
    require(alpha == real + i*imag, 'relative character radicand escaped the Gaussian base')
    return real + i_base*imag


def relative_character_basis(nf):
    """Recover a basis of quadratic characters of K/Q(i) by Galois projectors."""
    f = nf.nf_get_pol()
    degree = int(f.poldegree())
    roots = pari.nfroots(nf, x*x+1)
    if not len(roots):
        return None
    i, z = roots[0], pari.Mod(w, f)

    def act(z1, a):
        return pari.subst(pari.lift(a), w, z1)

    autos = [pari.Mod(t, f) for t in pari.nfgaloisconj(f)]
    autos = [t for t in autos if act(t,i) == i]
    if len(autos)*2 != degree or any(act(t,t) != z for t in autos):
        return None
    table = [z]
    for tau in autos:
        if tau not in table:
            table += [act(tau, s) for s in table]
    require(len(table) == len(autos), 'incorrect exponent-two automorphism basis')
    basis = []
    for bit in range(len(table).bit_length()-1):
        for j in range(1,degree):
            eig = sum((-1 if n & (1<<bit) else 1) * t**j for n,t in enumerate(table))
            if eig:
                basis.append(gaussian_coordinates(eig**2,i,degree))
                break
        else:
            raise ArithmeticError('quadratic Galois projector vanished on a field basis')
    return basis


def full_points(nf, mode):
    """Each abscissa polynomial comes from the specified proof."""
    equations = [x**6-1, x]
    if mode == 'negative even tower':
        equations.append(8*x**6+1)
    if mode == 'Gaussian quartic support':
        equations += [x*x+1, 2*x*x+1, x*x+2]
    if mode == 'inert rational support':
        equations += [2*x*x+1]
    points = set()
    for equation in equations:
        for a in pari.nfroots(nf, equation):
            for b in pari.nfroots(nf, x*x-a**6+1):
                require(b*b == a**6-1, 'source identity')
                points.add((str(pari.lift(a)), str(pari.lift(b))))
    return sorted(points)


def classify(coefficients):
    if not isinstance(coefficients,list) or len(coefficients)<2 or coefficients[-1]!=1 or any(type(c)!=int for c in coefficients):
        raise ValueError('Supply a monic polynomial as an ascending list of integers.')
    f = pari.Polrev(coefficients,w)
    if not pari.polisirreducible(f):
        raise ValueError('The defining polynomial must be irreducible.')
    degree = int(f.poldegree())
    nf = pari.nfinit(f)
    chars = quadratic_subfields(f)
    mode = None
    if len(chars)+1 == degree:
        mode = rational_tower(chars)
    if mode is None and degree == 2:
        if all(int(p)%3 != 1 for p in pari.factorint(abs(chars[0]))[0]):
            mode = 'inert rational support'
    if mode is None and degree & (degree-1) == 0:
        if -1 in chars:
            compositum = nf
        else:
            polys = pari.polcompositum(f,w*w+1)
            require(len(polys)==1, 'nondisjoint compositum with Gaussian field')
            compositum = pari.nfinit(polys[0])
        basis = relative_character_basis(compositum)
        if basis is not None:
            results = [gaussian_test(delta) for delta in basis]
            if all(kernel for support,kernel in results):
                mode = 'Gaussian radical tower'
            elif len(basis)==1 and all(support for support,kernel in results):
                mode = 'Gaussian quartic support'
    if mode is None:
        return {'field':str(f),'status':'outside these theorems'}
    points = full_points(nf,mode)
    return {'field':str(f),'status':'complete','method':mode,'verified_affine_points':points,
            'rational_infinities':2, 'total':len(points)+2}


def verify():
    tests = [
        ([1,0,0,0,1],18), ([1,0,-1,0,1],10), ([1,0,10,0,1],20),
        ([5,0,-2,0,1],6), ([17,0,-2,0,1],6), ([29,0,-10,0,1],6),
        ([-11,0,1],4), ([22,0,1],4), ([-19,0,1],None),
        ([-2,0,0,1],None), ([64,0,-12,0,1],None),
    ]
    for co, count in tests:
        result = classify(co)
        require(result.get('total') == count, f'field classification: {co}, {result}')
        print(str(pari.Polrev(co,w)), result['status'],result.get('method',''),count)
    # Compositum Q(i,sqrt(1+4i),sqrt(5+8i)), independently chosen over split primes 17 and 89.
    f1 = w**4 - 2*w*w + 17
    f2 = w**4 - 10*w*w + 89
    pols = pari.polcompositum(f1,f2)
    require(all(t.poldegree()==8 for t in pols),'non-Galois octic degree')
    for t in pols:
        co=list(map(int,t.Vecrev()))
        result=classify(co)
        require(result.get('total')==6 and result.get('method')=='Gaussian radical tower','octic tower input')
    print('Non-Galois octic inputs and both conjugate composita verified.')


if __name__ == '__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('ascending_polynomial',nargs='?')
    arguments=parser.parse_args()
    if arguments.ascending_polynomial:
        print(json.dumps(classify(json.loads(arguments.ascending_polynomial)),indent=2))
    else:
        verify()
