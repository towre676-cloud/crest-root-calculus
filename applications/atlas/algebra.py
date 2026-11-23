"""Exact arithmetic for integer ternary quartics and small finite fields."""
from functools import lru_cache
from itertools import product
from math import gcd as integer_gcd
import re


def mons(degree):
    return [(a, b, degree-a-b) for a in range(degree, -1, -1)
            for b in range(degree-a, -1, -1)]


def parse_form(text):
    """Read a signed sum of integer monomials without executing input."""
    if not isinstance(text, str) or len(text) > 2000:
        raise ValueError("quartic must be an expression of at most 2000 characters")
    src = text.replace(" ", "")
    pieces = re.findall(r"([+-]?)([^+-]+)", src)
    if not pieces or "".join(sign + term for sign, term in pieces) != src:
        raise ValueError("invalid quartic expression")
    out = {}
    for sign, term in pieces:
        powers = [0, 0, 0]
        c = -1 if sign == "-" else 1
        for item in term.split("*"):
            if re.fullmatch(r"[0-9]{1,12}", item):
                c *= int(item)
                continue
            var = re.fullmatch(r"([xyz])(?:\^([1-4]))?", item)
            if not var:
                raise ValueError("use x, y, z, integer coefficients, *, +, -, and ^")
            powers["xyz".index(var[1])] += int(var[2] or 1)
        if sum(powers) != 4:
            raise ValueError("every monomial must have total degree four")
        k = tuple(powers)
        out[k] = out.get(k, 0) + c
    out = {k: v for k,v in out.items() if v}
    if not out:
        raise ValueError("the zero polynomial is not a quartic")
    return out


def rank(rows, p):
    pivots = {}
    for row in rows:
        for j, value in enumerate(row):
            coeff = value % p
            if coeff == 0:
                continue
            if j in pivots:
                base = pivots[j]
                for k in range(j, len(row)):
                    row[k] = (row[k] - coeff * base[k]) % p
            else:
                inv = pow(coeff, -1, p)
                pivots[j] = [v*inv % p for v in row]
                break
    return len(pivots)


def smoothness(form, p):
    gradients = []
    for v in range(3):
        g = {}
        for key, coeff in form.items():
            if key[v]:
                mon = list(key)
                mon[v] -= 1
                g[tuple(mon)] = coeff * key[v]
        gradients.append(g)
    # Full column rank puts every monomial of this degree in the singular ideal.
    for degree in [7, 9]:
        columns = {key:i for i,key in enumerate(mons(degree))}
        rows = []
        for generator, d in [(form, 4)] + [(g, 3) for g in gradients]:
            for shift in mons(degree-d):
                row = [0]*len(columns)
                for key, coeff in generator.items():
                    row[columns[tuple(a+b for a,b in zip(key, shift))]] = coeff % p
                rows.append(row)
        r = rank(rows, p)
        if r == len(columns):
            return {"degree":degree,"rank":r,"columns":len(columns)}
    return None


@lru_cache(maxsize=9)
def field(p, degree):
    q = p**degree
    if degree == 1:
        modulus = (0,)
    else:
        modulus = next((c,b)+(0,)*(degree-2) for b,c in product(range(p),range(1,p))
                       if all((x**degree+b*x+c)%p for x in range(p)))
    digits = [tuple(x//p**k%p for k in range(degree)) for x in range(q)]
    def encode(vec):
        return sum((a%p)*p**k for k,a in enumerate(vec))
    def multiply(a,b):
        tmp = [0]*(2*degree-1)
        for j,v in enumerate(a):
            for k,w in enumerate(b):
                tmp[j+k] += v*w
        for k in range(len(tmp)-1, degree-1, -1):
            for j,c in enumerate(modulus):
                tmp[k-degree+j] -= tmp[k]*c
        return encode(tmp[:degree])
    add = [[encode([i+j for i,j in zip(a,b)]) for b in digits] for a in digits]
    mul = [[multiply(a,b) for b in digits] for a in digits]
    return add, mul


def point_count(form, p, n):
    add, mul = field(p,n)
    q = len(add)
    pw = [[1,x,mul[x][x],mul[mul[x][x]][x],mul[mul[x][x]][mul[x][x]]] for x in range(q)]
    # On z=1, evaluate a quartic in y for each x by Horner's rule.
    count = 0
    for x in range(q):
        coeff = [0]*5
        for (a,b,_), c in form.items():
            coeff[b] = add[coeff[b]][mul[c%p][pw[x][a]]]
        c0,c1,c2,c3,c4 = coeff
        for y in range(q):
            val = add[c3][mul[y][c4]]
            val = add[c2][mul[y][val]]
            val = add[c1][mul[y][val]]
            count += add[c0][mul[y][val]] == 0
    # At infinity the projective charts are [1:y:0] and [0:1:0].
    top = [(b,c%p) for (_,b,k),c in form.items() if k == 0]
    for y in range(q):
        v = 0
        for b,c in top:
            v = add[v][mul[c][pw[y][b]]]
        count += v == 0
    return count + (form.get((0,4,0),0)%p == 0)


def characteristic(p, counts):
    s1,s2,s3 = [p**i+1-a for i,a in enumerate(counts,1)]
    e2=(s1*s1-s2)//2
    e3=(s3-s1*s2+e2*s1)//3
    return [p**3,-p*p*s1,p*e2,-e3,e2,-s1,1]  # ascending powers


def trim(a):
    while a and not a[-1]:
        a.pop()
    return a


def rem(a, b, ell):
    a = trim([x%ell for x in a])
    while len(a) >= len(b):
        shift = len(a)-len(b)
        f = a[-1]*pow(b[-1],-1,ell)%ell
        for j,v in enumerate(b):
            a[shift+j] = (a[shift+j]-f*v)%ell
        trim(a)
    return a


def gcd(a, b, ell):
    while b:
        a,b = b,rem(a,b,ell)
    return a


def power(a, m, modulus, ell):
    def mul(b,c):
        tmp=[0]*(len(b)+len(c)-1)
        for j,v in enumerate(b):
            for k,w in enumerate(c):
                tmp[j+k] += v*w
        return rem(tmp,modulus,ell)
    result=[1]
    while m:
        if m%2:
            result=mul(a,result)
        a=mul(a,a)
        m//=2
    return result


def difference(a,b,ell):
    return trim([((a[i] if i<len(a) else 0)-(b[i] if i<len(b) else 0))%ell
                 for i in range(max(len(a),len(b)))])


def phi(n):
    return sum(integer_gcd(k,n)==1 for k in range(1,n+1))


@lru_cache(maxsize=1)
def possible_orders():
    # phi(n) >= sqrt(n/2) and a root ratio has at most 6*5 conjugates.
    return [m for m in range(2,1801) if phi(m)<=30]


def simplicity_witness(polynomial):
    x=[0,1]
    for ell in (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71):
        f=[i%ell for i in polynomial]
        if not f[0]:
            continue
        if any(len(gcd(f,difference(power(x,ell**k,f,ell),x,ell),ell))>1 for k in (2,3)):
            continue
        if power(x,ell**6,f,ell)!=x:
            continue
        # Irreducible mod ell. All six conjugates of x^m must remain distinct
        # to rule out a characteristic-zero root ratio of order m at this prime.
        for m in possible_orders():
            if any(power(x,m*(ell**k-1),f,ell)==[1] for k in (2,3)):
                break
        else:
            return {"irreducible_and_distinct_powers_mod":ell,"orders_checked":len(possible_orders())}
    return None
