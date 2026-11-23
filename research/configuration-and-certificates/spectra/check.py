"""Check a weak-slab optical fibre and exact affine exclusions (standard library)."""

import argparse
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIELDS = ("wavelength_um", "calcite", "clay", "observed")


def physical_witnesses(source=None):
    """Rows are native USGS channels, one-based after the header; round to six places."""
    results = {}
    for sample in json.loads((HERE / "bands.json").read_text()):
        bands = sample["bands"]
        if source:
            for field, manifest in zip(
                FIELDS, sample["native_band_files"], strict=True
            ):
                data = (source / manifest["path"]).read_bytes()
                assert sha256(data).hexdigest() == manifest["sha256"]
                lines = data.decode("ascii").splitlines()
                for band in bands:
                    raw = D(lines[band["row"]])
                    assert format(raw, ".6f") == band[field]
        intervals = []
        for band in bands:
            a, b, observed = map(F, (band[name] for name in FIELDS[1:]))
            eps = F(1, 200)
            lo, hi = sorted((observed + sign * eps - b) / (a - b) for sign in (-1, 1))
            lo, hi = max(F(0), lo), min(F(1), hi)
            assert lo <= hi
            intervals.append((lo, hi))
        assert intervals[0][1] < intervals[1][0]
        results[sample["sample"]] = [[str(lo), str(hi)] for lo, hi in intervals]
    return results


def dot(x, y):
    return sum(a * b for a, b in zip(x, y, strict=True))


def profiled_information(score, du, dv):
    """Squared score length after least-squares projection off two nuisance columns."""
    aa, ab, bb = dot(du, du), dot(du, dv), dot(dv, dv)
    determinant = aa * bb - ab * ab
    assert determinant > 0
    au, av = dot(score, du), dot(score, dv)
    first, second = (bb * au - ab * av) / determinant, (aa * av - ab * au) / determinant
    residual = [
        s - first * a - second * b for s, a, b in zip(score, du, dv, strict=True)
    ]
    return dot(residual, residual)


def slab(z):
    external, internal = D(".08"), D(".43")
    transmission = (-z).exp()
    denominator = 1 - internal * transmission
    albedo = 1 - (1 - external) * (1 - transmission) / denominator
    derivative = -(1 - external) * (1 - internal) * transmission / denominator**2
    return albedo, derivative


def information(scale):
    u, v, x = D(".25"), D(".5"), D(".4")
    bands = [(".2", ".9"), (".8", ".15"), (".65", ".7"), ("1", ".4"), (".3", "1")]
    albedo_cols, reflectance_cols = [[], [], []], [[], [], []]
    for ka, kb in bands:
        ka, kb = scale * D(ka), scale * D(kb)
        za, zb = u * ka / x, v * kb / (1 - x)
        a, ap = slab(za)
        b, bp = slab(zb)
        w = x * a + (1 - x) * b
        deriv = (a - za * ap - b + zb * bp, ka * ap, kb * bp)
        root = (1 - w).sqrt()
        reflectance_derivative = 1 / (root * (1 + root) ** 2)
        for i in range(3):
            albedo_cols[i].append(deriv[i])
            reflectance_cols[i].append(deriv[i] * reflectance_derivative)
    return [profiled_information(*cols) for cols in (albedo_cols, reflectance_cols)]


def weak_fibre():
    with localcontext() as ctx:
        ctx.prec = 65
        results = [information(D(10) ** (-i)) for i in range(2, 7)]
        slopes = [
            [(left / right).log10() for left, right in zip(a, b, strict=True)]
            for a, b in zip(results[:-1], results[1:], strict=True)
        ]
        assert abs(slopes[-1][0] - 4) < D(".001")
        assert abs(slopes[-1][1] - 3) < D(".01")
        # A known density-times-reference-diameter ratio of 2 is illustrative.
        u, v = F(1, 4), F(1, 2)
        mass_fractions = []
        for x in (F(1, 4), F(2, 5), F(3, 5)):
            sa, sb = u / x, v / (1 - x)
            ratio = 2 * x * sa / ((1 - x) * sb)
            mass_fractions.append(str(ratio / (1 + ratio)))
        assert mass_fractions == ["1/2"] * 3
        return {
            "example_mass_fractions": mass_fractions,
            "decades_of_eta": ["1e-2", "1e-3", "1e-4", "1e-5", "1e-6"],
            "information_albedo_and_reflectance": [
                [format(v, ".10e") for v in row] for row in results
            ],
            "successive_log_slopes": [
                [format(v, ".6f") for v in row] for row in slopes
            ],
        }


def simplex_witness():
    lo, hi = F(35, 100), F(45, 100)
    # Each band reads just one simplex coordinate. To satisfy any two, set both at lo.
    assert 0 <= 1 - 2 * lo <= 1
    assert 2 * lo <= 1 and 3 * lo > 1 and lo < hi
    return {
        "two_bands_feasible_at": [str(lo), str(lo), str(1 - 2 * lo)],
        "three_bands_require_sum_at_least": str(3 * lo),
    }


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument(
        "--check", action="store_true", help="compare to the stored receipt"
    )
    cli.add_argument(
        "--source", type=Path, help="optional extracted assay_selector/data directory"
    )
    args = cli.parse_args()
    result = {
        "affine_band_intervals": physical_witnesses(args.source),
        "three_band_simplex": simplex_witness(),
        "weak_absorption": weak_fibre(),
    }
    if args.check:
        assert result == json.loads((HERE / "receipt.json").read_text())
        print(
            "PASS: two physical witnesses, three-band simplex, mass fibre and both information exponents"
        )
    else:
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
