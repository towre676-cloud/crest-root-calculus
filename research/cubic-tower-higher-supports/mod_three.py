"""Signed good-reduction tower numerators via quadratic resolvents modulo 3."""
from itertools import combinations
import json
from pathlib import Path
import sys

import numpy as np

# The general tower routines live beside this research directory.
HERE = Path(__file__).resolve().parent
TOWER = HERE.parent / 'cubic-tower-factors'
sys.path.insert(0, str(TOWER))
from field import Field
from tower import check_stage, trace_from_factor, trace_sums, weil_polynomial


def character(logs):
    return np.where(logs < 0, 0, 1-2*(logs % 2))


def walsh(values):
    out = values.copy()
    width = 1
    while width < len(out):
        pairs = out.reshape(-1, 2*width)
        left = pairs[:, :width].copy()
        right = pairs[:, width:].copy()
        pairs[:, :width] = left+right
        pairs[:, width:] = left-right
        width *= 2
    return out


def hyperelliptic_traces(p, indices, degree):
    field = Field(p, degree)
    n = len(indices)
    masks = np.arange(1 << n)
    hist = np.zeros((2, len(masks)), dtype=np.int64)
    special = np.zeros_like(hist)
    sign = np.array([(-1)**int(m).bit_count() for m in masks])
    leading_negative = 0
    for i, index in enumerate(indices):
        a = -27*(12*index+1) % p
        if character(field.log[a]) == -1:
            leading_negative |= 1 << i
    for start in range(0, field.q, 65536):
        enc = np.arange(start, min(start+65536, field.q))
        z = field.log[enc]
        z2 = np.where(z < 0, -1, 2*z % field.order)
        z3 = np.where(z < 0, -1, 3*z % field.order)
        f = field.add(field.add(field.shift(z3, field.log[-4 % p]), z2), field.add(z, 0))
        chi_f = character(f)
        negative, zero = np.zeros_like(z), np.zeros_like(z)
        for j, i in enumerate(indices):
            a, b, c = (-27*(12*i+1)) % p, (-27*(12*i*i+1)) % p, (-27*(4*i**3+1)) % p
            az2 = field.shift(z2, field.log[a])
            bz = field.shift(z, field.log[b]) if b else np.full_like(z, -1)
            h = field.add(field.add(az2, bz), field.log[c])
            chi = character(h)
            negative |= (chi == -1).astype(np.int64) << j
            zero |= (chi == 0).astype(np.int64) << j
        regular = zero == 0
        hist[0] += np.bincount(negative[regular], minlength=len(masks))
        for s in (-1, 1):
            hist[1] += s*np.bincount(negative[regular & (chi_f == s)], minlength=len(masks))
        for j in np.flatnonzero(~regular):
            value = (masks & zero[j] == 0)*sign[masks & negative[j]]
            special[0] += value
            special[1] += chi_f[j]*value
    result = -np.array([walsh(hist[b])+special[b] for b in (0, 1)])
    # Even degree: two or zero rational points at infinity; odd degree: one.
    result[0] -= sign[masks & leading_negative]
    return result


def product3(factors):
    out = np.array([1], dtype=np.int64)
    for f in factors:
        out = np.convolve(out, f) % 3
    return out.tolist()


def remove_base_denominator(poly, p):
    # Divide by (1-T)(1-pT), from its constant end in F_3[T].
    quotient = []
    remainder = list(poly)
    for i in range(len(remainder)-2):
        a = remainder[i] % 3
        quotient.append(a)
        remainder[i+1] += (1+p)*a
        remainder[i+2] -= p*a
    if any(v % 3 for v in remainder[-2:]):
        raise ArithmeticError('Quadratic character product failed division')
    return quotient


def resolvent_factors(p, indices, max_size):
    check_stage(p, indices)
    traces = [hyperelliptic_traces(p, indices, d) for d in range(1, max_size+2)]
    hyper = {(0, 1): weil_polynomial([int(t[1, 0]) for t in traces], p, 1)}
    factors = {"base": [v % 3 for v in hyper[0, 1]]}
    for size in range(1, max_size+1):
        for subset in combinations(range(len(indices)), size):
            mask = sum(1 << i for i in subset)
            for b in (0, 1):
                g = size-1+2*b
                hyper[mask, b] = weil_polynomial([int(t[b, mask]) for t in traces], p, g)
            label = ','.join(str(indices[i]) for i in subset)
            pos, neg = [], [hyper[0, 1]]
            u = mask
            while u:
                parity = u.bit_count() % 2
                pos.append(hyper[u, parity])
                neg.append(hyper[u, 1-parity])
                u = (u-1) & mask
            factors[label+'/+'] = remove_base_denominator(product3(pos), p)
            factors[label+'/-'] = product3(neg)
            if any(len(factors[label+'/'+s])-1 != size*(2**size) for s in ('+', '-')):
                raise ArithmeticError('Wrong signed factor degree')
    return factors


def verify():
    source = json.loads((TOWER/'receipt.json').read_text())
    result = {}
    for p, indices, previous in ((11, [0, 1, 2], 'stage_two'), (41, list(range(10)), 'stage_ten')):
        factors = resolvent_factors(p, indices, 3)
        exact = source[previous]['complete_factors_low_to_high']
        for label, poly in exact.items():
            if factors[label] != [v % 3 for v in poly]:
                raise ArithmeticError('Integral prior factor disagrees at '+label)
        triples = {label: poly for label, poly in factors.items() if label.count(',') == 2}
        record = {'signed_triples_low_to_high_mod_3': triples}
        if p == 11:
            full = product3(factors.values())
            counts = []
            for d in range(1, 7):
                _, signed, points = trace_sums(p, indices, d, 3)
                pos, neg = signed[tuple(indices)]
                for s, t in (('+', pos), ('-', neg)):
                    if (trace_from_factor(triples['0,1,2/'+s], d)-t) % 3:
                        raise ArithmeticError('Direct signed trace disagrees')
                if (p**d+1-trace_from_factor(full, d)-points[0]) % 3:
                    raise ArithmeticError('Direct genus-55 trace disagrees')
                twisted = sum((1 if name.endswith('/+') else -1)*trace_from_factor(poly, d)
                              for name, poly in factors.items())
                if (p**d+1-twisted-points[1]) % 3:
                    raise ArithmeticError('Direct twisted genus-55 trace disagrees')
                counts.append({'degree': d, 'ordinary_and_twisted_counts': points, 'signed_triple_traces': [pos, neg]})
            record.update(full_numerator_low_to_high_mod_3=full, direct_counts=counts)
        result[str(p)] = record
    return result


if __name__ == '__main__':
    data = verify()
    if sys.argv[1:] == ['--compute']:
        print(json.dumps(data, indent=2))
    else:
        if data != json.loads((HERE/'mod-three.json').read_text()):
            raise ArithmeticError('Mod-3 receipt differs')
        print('Verified the full genus-55 numerator modulo 3 and all 240 signed triples at 41.')
