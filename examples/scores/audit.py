"""Separate a score from a known ingredient, then compare the rankings."""
import argparse
import csv
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from fractions import Fraction
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1] / "docs/figures/data/lowest-zero-audit.csv"


@dataclass(frozen=True)
class Score:
    name: str
    total: Fraction
    ingredient: Fraction
    remainder: Fraction


def number(value):
    try:
        parsed = Decimal(value)
    except (InvalidOperation, TypeError, ValueError):
        raise ValueError(f"expected a decimal number, got {value!r}") from None
    if not parsed.is_finite():
        raise ValueError(f"expected a finite number, got {value!r}")
    return Fraction(parsed)


def decimal_text(value):
    """All inputs are finite decimals, so the long division always terminates."""
    if not value:
        return "0"
    sign = "-" if value < 0 else ""
    whole, rest = divmod(abs(value.numerator), value.denominator)
    digits = []
    while rest:
        digit, rest = divmod(10 * rest, value.denominator)
        digits.append(str(digit))
    return sign + str(whole) + ("." + "".join(digits) if digits else "")


def load_rows(path, columns, weight):
    rows = []
    seen = set()
    with path.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, strict=True)
        header = reader.fieldnames or []
        if len(set(header)) != len(header):
            raise ValueError("duplicate column names in the header")
        missing = [field for field in columns if field not in header]
        if missing:
            raise ValueError("missing CSV columns: " + ", ".join(missing))
        name_col, total_col, part_col = columns
        for row in reader:
            try:
                if None in row or any(v is None for v in row.values()):
                    raise ValueError("the row length differs from the header")
                name = row[name_col].strip()
                if not name or name in seen:
                    raise ValueError(f"each row needs its own nonempty identifier, got {name!r}")
                total = number(row[total_col])
                ingredient = weight * number(row[part_col])
                rows.append(Score(name, total, ingredient, total - ingredient))
                seen.add(name)
            except ValueError as e:
                raise ValueError(f"line {reader.line_num}: {e}") from e
    if not rows:
        raise ValueError("supply at least one row")
    return rows


def ranks(values):
    result = {}
    for position, value in enumerate(sorted(values, reverse=True), start=1):
        result.setdefault(value, position)
    return result


def correlation(xs, ys):
    """Pearson's r; form the bounded squared ratio before converting to float."""
    n = len(xs)
    if n < 3:
        return None
    c = n * sum(x*y for x, y in zip(xs, ys)) - sum(xs)*sum(ys)
    vx = n * sum(x*x for x in xs) - sum(xs)**2
    vy = n * sum(y*y for y in ys) - sum(ys)**2
    if not vx or not vy:
        return None
    return math.copysign(math.sqrt(float(c*c / (vx*vy))), -1 if c < 0 else 1)


def run(args):
    if args.research:
        path = RESEARCH
        columns = ("twist_d", "whole_score_scaled", "closest_pair_scaled")
        if args.file or any((args.id, args.total, args.component)) or number(args.weight) != 1:
            raise ValueError("--research uses the archived columns and a weight of 1")
    else:
        path = args.file or HERE / "players.csv"
        columns = (args.id or "player", args.total or "total", args.component or "coins")
    if len(set(columns)) != 3:
        raise ValueError("choose three different column names")
    if args.show < 1:
        raise ValueError("--show must be at least 1")
    if args.out and args.out.resolve() == path.resolve():
        raise ValueError("choose an output path different from the input CSV")
    rows = load_rows(path, columns, number(args.weight))
    totals = ranks(r.total for r in rows)
    remainders = ranks(r.remainder for r in rows)
    ordered = sorted(rows, key=lambda r: (-r.total, r.name))

    def fields(r):
        return [r.name, decimal_text(r.total), decimal_text(r.ingredient),
                decimal_text(r.remainder), totals[r.total], remainders[r.remainder]]

    header = ["id", "total", "included_points", "remaining_points", "rank_before", "rank_after"]
    print(f"{len(rows)} rows. Remaining points = total minus included points.")
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(header)
    for r in ordered[:args.show]:
        writer.writerow(fields(r))
    if len(rows) > args.show:
        print(f"Showing the top {args.show} by total. Use --out for all rows.")
    if len(rows) >= 3:
        parts = [r.ingredient for r in rows]
        for name, values in (("whole score", [r.total for r in rows]),
                             ("remainder", [r.remainder for r in rows])):
            coefficient = correlation(parts, values)
            result = "undefined (constant column)" if coefficient is None else f"{coefficient:.4f}"
            print(f"Pearson r, included points vs {name}: {result}")
    if args.out:
        # Exclusive creation keeps existing datasets and prior exports intact.
        with args.out.open("x", encoding="utf-8", newline="") as dest:
            out = csv.writer(dest)
            out.writerow(header)
            out.writerows(fields(r) for r in ordered)
        print(f"Wrote all {len(rows)} rows to {args.out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", nargs="?", type=Path, help="CSV with one observation per row")
    parser.add_argument("--research", action="store_true", help="use the 201 published curve records")
    parser.add_argument("--id", metavar="COLUMN", help="row identifier (default: player)")
    parser.add_argument("--total", metavar="COLUMN", help="whole score (default: total)")
    parser.add_argument("--component", metavar="COLUMN", help="additive ingredient (default: coins)")
    parser.add_argument("--weight", default="1", help="points per unit of the ingredient")
    parser.add_argument("--show", type=int, default=8, metavar="N", help="maximum displayed rows")
    parser.add_argument("--out", type=Path, metavar="NEW.csv", help="write exact values and ranks")
    args = parser.parse_args()
    try:
        run(args)
    except (OSError, ValueError, csv.Error) as e:
        parser.exit(2, f"{parser.prog}: {e}\n")


if __name__ == "__main__":
    main()
