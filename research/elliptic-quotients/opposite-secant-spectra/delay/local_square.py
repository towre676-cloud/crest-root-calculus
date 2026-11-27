"""Certify local-square exclusions for signed 3-neighbour theta sums."""

import argparse
from importlib.util import module_from_spec, spec_from_file_location
import json
from math import gcd, lcm
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

# The inverse-spectrum package also has a module named census.
spec = spec_from_file_location("delay_census", Path(__file__).with_name("census.py"))
delay_census = module_from_spec(spec)
spec.loader.exec_module(delay_census)
pari = delay_census.pari


def theta_level(gram):
    inverse = gram**-1
    return lcm(*(int((inverse[i,j] / (4 if i == j else 2)).denominator())
                 for i in range(3) for j in range(i+1)))


def zero_bound(level):
    index = level
    for p in map(int, pari.factor(level)[0]):
        index = index//p*(p+1)
    # Multiplication by the weight-1/2 theta series preserves order.
    return index//6


def unit_square(n, p):
    if p == 2:
        return n%8 == 1
    return pow(n%p, (p-1)//2, p) == 1


def local_square(n, p):
    valuation = 0
    while n%p == 0:
        n //= p
        valuation += 1
    return valuation%2 == 0 and unit_square(n,p)


class ThetaDifference:
    def __init__(self, forms):
        self.forms = forms
        self.level = theta_level(forms[0])
        if any(theta_level(q) != self.level for q in forms[1:]):
            raise ArithmeticError("The neighbour theta levels differ.")
        self.conductor = int(pari.coredisc(forms[0].matdet()))
        if self.level % self.conductor:
            raise ArithmeticError("The theta character conductor does not divide its level.")
        self.coefficients = None
        self.checked = 0

    def through(self, bound):
        if bound > self.checked:
            a,b,c,d = (q.qfrep(bound) for q in self.forms)
            self.coefficients = a+b-c-d
            self.checked = bound
        return self.coefficients

    def obstruction(self, p):
        """Return a Sturm-certified local obstruction, or an exact contrary witness."""
        level = lcm(self.level,4*p)
        bound = zero_bound(level)
        modulus = 8 if p == 2 else p
        projected_level = lcm(level,modulus**2,modulus*self.conductor)
        projected_bound = zero_bound(projected_level)
        vectors = []
        for step in range(8):
            factor = p**(2*step)
            coefficients = self.through(bound*factor)
            vectors.append(pari.Col([coefficients[m*factor-1] for m in range(1,bound+1)]))
            kernel = pari.Mat(vectors).matker()
            if len(kernel):
                relation = list(map(int,pari.content(kernel[0])**-1*kernel[0]))
                if not relation[-1]:
                    raise ArithmeticError("The previous prefix vectors should be independent.")
                if relation[-1] < 0:
                    relation = [-entry for entry in relation]
                return {"prime": p, "orbit_dimension": step, "relation": relation,
                        "orbit_level": level, "orbit_bound": bound,
                        "projection_level": projected_level, "projection_bound": projected_bound}
            coefficients = self.through(projected_bound*factor)
            witness = next((m*factor for m in range(1,projected_bound+1)
                            if unit_square(m,p) and coefficients[m*factor-1]),None)
            if witness is not None:
                return {"norm": witness, "signed_coefficient": 2*int(coefficients[witness-1])}
        raise ArithmeticError("The U(p²) orbit did not close within eight tests.")


def case_data(entry):
    parent = tuple(entry["parent"])
    edges = tuple(tuple(tuple(point) for point in edge) for edge in entry["edges"])
    det,a,b,c,u,v,w = parent
    gram = ((a,u,v),(u,b,w),(v,w,c))
    g = delay_census.gram_matrix(gram)
    primes = sorted({2,*map(int,pari.factor(det)[0])})
    if int(g.matdet()) != det or det%3 != 1 or delay_census.local_obstruction(parent,primes):
        raise ArithmeticError("Invalid realized parent.")
    field = delay_census.residue_field(-4,3)
    points,_ = delay_census.normal_partition(gram,field)
    if sorted(edges[0]+edges[1]) != sorted(points):
        raise ArithmeticError("The two labels are not opposite secants.")
    children = [delay_census.gram_matrix(delay_census.sublattice(gram,delay_census.cross(*edge),3)) for edge in edges]
    if any(child.qfrep(1)[0] for child in children) or children[0].qfisom(children[1]):
        raise ArithmeticError("Expected distinct non-product children.")
    matching = dict(zip(points,delay_census.neighbours(g,3),strict=True))
    left,right = ([matching[point] for point in edge] for edge in edges)
    if any(left[0].qfisom(right[j]) and left[1].qfisom(right[1-j]) for j in (0,1)):
        raise ArithmeticError("An ordinary neighbour isometry already settles this pair.")
    return (parent,edges,0),ThetaDifference([*left,*right]),primes


def check_central_zero_packet(theta):
    """Identify the Hecke packet and certify its pertinent central values."""
    equation = [0,1,0,-80,-304]
    curve = pari.ellinit(equation)
    level = int(pari.ellglobalred(curve)[0])
    if theta.level not in (928,1856) or level != 232 or theta.conductor != 232:
        raise ArithmeticError("Unexpected theta or elliptic-curve level.")

    primes = (3,5,17)
    eigenvalues = [int(pari.ellap(curve,p)) for p in primes]
    candidates = []
    for divisor in map(int,pari.divisors(theta.level//2)):
        space = pari.mfinit([divisor,2],0)
        for eigenform in pari.mfeigenbasis(space):
            if all(pari.mfcoef(eigenform,p) == value for p,value in zip(primes,eigenvalues)):
                candidates.append(divisor)
    if candidates != [level] or abs(eigenvalues[0]) == primes[0]+1:
        raise ArithmeticError("The good-prime fingerprint does not isolate the cuspidal packet.")

    bound = zero_bound(theta.level)
    coefficients = theta.through(max(primes)**2*bound)
    det = int(theta.forms[0].matdet())
    for p,value in zip(primes,eigenvalues):
        for n in range(1,bound+1):
            middle = pari.kronecker(-det*n,p)*coefficients[n-1]
            lower = p*coefficients[n//(p*p)-1] if n%(p*p) == 0 else 0
            if coefficients[n*p*p-1]+middle+lower != value*coefficients[n-1]:
                raise ArithmeticError("The parent theta sum failed an exact Hecke eigenvalue check.")

    t = 33
    lift_bound = zero_bound(theta.level//2)
    coefficients = theta.through(t*lift_bound**2)
    square_class = [2*int(coefficients[t*m*m-1]) for m in range(1,lift_bound+1)]
    lift = [0]*lift_bound
    for d in range(1,lift_bound+1):
        if gcd(d,theta.level) != 1:
            continue
        character = int(pari.kronecker(-det*t,d))
        for m in range(1,lift_bound//d+1):
            lift[d*m-1] += character*square_class[m-1]
    newform = list(map(int,pari.ellan(curve,lift_bound)))
    dilations = [1] if theta.level == 928 else [1,2]
    expected = [8*sum(newform[m//d-1] for d in dilations if m%d == 0)
                for m in range(1,lift_bound+1)]
    if lift[0] != 8 or lift != expected:
        raise ArithmeticError("The square-class-33 lift differs from the claimed newform.")

    symbols, minus = pari.msfromell(curve,-1)
    infinity = pari("oo")
    twists = []
    for square_class,replay_value in ((1,0),(33,128)):
        disc = int(pari.coredisc(-det*square_class))
        central = sum(
            pari.kronecker(disc,a)*pari.mseval(symbols,minus,[infinity,pari(a)/abs(disc)])
            for a in range(1,abs(disc)) if gcd(a,abs(disc)) == 1
        )
        sign = int(pari.ellrootno(pari.elltwist(curve,disc)))
        if central != replay_value or sign != 1:
            raise ArithmeticError("The exact twisted modular-symbol value or sign changed.")
        twists.append({"square_class": square_class, "twist_discriminant": disc,
                       "root_number": sign, "twisted_symbol_sum": int(central)})
    auxiliary_level = int(pari.ellglobalred(pari.elltwist(curve,8))[0])
    if auxiliary_level != 1856 or auxiliary_level%16:
        raise ArithmeticError("The auxiliary twist does not meet Waldspurger's 2-adic hypothesis.")
    return {"curve_weierstrass_coefficients": equation, "newform_conductor": level,
            "hecke_primes": primes, "eigenvalues": eigenvalues,
            "theta_eigenvalue_bound": bound, "shimura_33_sturm_bound": lift_bound,
            "shimura_33_multiplier": 8, "shimura_33_dilations": dilations,
            "auxiliary_twist_discriminant": 8,
            "auxiliary_twist_conductor": auxiliary_level, "twists": twists}


def check_progression_identity(records):
    matching = [r for r in records if r["parent"][0] == 7168 and "simultaneously_local_square" in r]
    first,second = matching
    _,a,_ = case_data(first)
    _,b,_ = case_data(second)
    level = lcm(a.level,b.level,64,8*a.conductor)
    bound = zero_bound(level)
    ca,cb = a.through(bound),b.through(bound)
    for n in range(1,bound+1):
        expected = ca[n-1]*((n%8 == 5)-(n%8 == 1))
        if cb[n-1] != expected:
            raise ArithmeticError("The determinant-7168 progression identity failed.")
    second["parent_theta_identity"] = {"type": "progression_sieve", "level": level, "bound": bound,
                                       "source_case": records.index(first)}


def certify(entries):
    pari.allocatemem(64_000_000,768_000_000,silent=True)
    records = []
    remaining = []
    field = delay_census.residue_field(-4,3)
    for entry in entries:
        case, theta, primes = case_data(entry)
        prefix = theta.through(2048)
        support = [n for n in range(1,2049) if prefix[n-1]]
        if not support:
            raise ArithmeticError("The signed ordinary theta sum should be nonzero.")
        witness = next((n for n in support if all(local_square(n,p) for p in primes)),None)
        record = {"parent": case[0], "edges": case[1],
                  "ordinary_first_difference": support[0],
                  "theta_level": theta.level, "theta_character_conductor": theta.conductor}
        if witness:
            record["simultaneously_local_square"] = {
                "norm": witness, "signed_coefficient": 2*int(prefix[witness-1]),
                "primes": primes,
            }
            remaining.append((len(records),case))
            if case[0][0] in (928,3712):
                record["waldspurger_packet"] = check_central_zero_packet(theta)
        else:
            candidates = [p for p in primes if not any(local_square(n,p) for n in support)]
            if not candidates:
                raise ArithmeticError("Expected either a simultaneous local square or a candidate prime.")
            record["local_obstructions"] = [theta.obstruction(p) for p in candidates]
            if any("norm" in result for result in record["local_obstructions"]):
                raise ArithmeticError(f"A proposed local explanation fails beyond the prefix: {record}")
        records.append(record)
    with TemporaryDirectory() as directory:
        binary = str(Path(directory)/"shells")
        subprocess.run(["g++","-std=c++17","-O3","-Wall","-Wextra","-Werror",
                        str(Path(__file__).with_name("shells.cpp")),"-o",binary],check=True)
        certificates,_ = delay_census.neighbour_comparison(binary,remaining,field)
    for index,(method,delta,_) in certificates.items():
        if method != "neighbour_sturm" or delta:
            raise ArithmeticError("The remaining case is not an all-square Sturm equality.")
        records[index]["square_equality"] = method
    check_progression_identity(records)
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--pairs", type=Path)
    source.add_argument("--verify", type=Path)
    options = parser.parse_args()
    inputs = json.loads((options.pairs or options.verify).read_text())
    keys = {(tuple(row["parent"]), tuple(tuple(tuple(p) for p in e) for e in row["edges"]))
            for row in inputs}
    if len(keys) != len(inputs):
        raise ArithmeticError("Repeated input pair.")
    receipt = json.dumps(certify(inputs),indent=2)
    if options.verify and json.loads(receipt) != inputs:
        raise ArithmeticError("The local-square replay differs from the receipt.")
    print(receipt)


if __name__ == "__main__":
    main()
