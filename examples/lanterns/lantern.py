"""Count colored rows of lanterns or retrieve one by its one-based index."""

import argparse


def ending_counts(length, colors, rule):
    """Count nonempty rows by whether the final lantern is dark or lit."""
    dark, lit = 1, colors
    for _ in range(1, length):
        next_dark = dark + lit if rule == "separated" else lit
        dark, lit = next_dark, dark * colors
    return dark, lit


def count_rows(length, colors, rule):
    if length == 0:
        return 1
    return sum(ending_counts(length, colors, rule))


def row_at_index(length, colors, rule, index):
    """Unrank in order of the symbols dark, amber, blue, color-3, ... ."""
    # With zero places left, there is exactly one completion after either color.
    after_dark = [1]
    after_light = [1]
    for _ in range(length):
        light_choices = colors * after_light[-1]
        next_dark = light_choices + (after_dark[-1] if rule == "separated" else 0)
        after_light.append(after_dark[-1])
        after_dark.append(next_dark)

    rank = index - 1
    previous = None
    result = []
    for remaining in reversed(range(length)):
        dark_allowed = previous != 0 or rule == "separated"
        dark_choices = after_dark[remaining] if dark_allowed else 0
        if rank < dark_choices:
            symbol = 0
        else:
            # Every possible light color has the same number of completions.
            color, rank = divmod(rank - dark_choices, after_light[remaining])
            symbol = color + 1
        result.append(symbol)
        previous = symbol
    return result


def lamp_name(symbol):
    if symbol < 3:
        return ("dark", "amber", "blue")[symbol]
    return f"color-{symbol}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--length", type=int, default=8, help="number of lanterns (default 8)")
    parser.add_argument("--colors", type=int, default=2, help="number of light colors (default 2)")
    parser.add_argument("--rule", choices=("both", "alternating", "separated"), default="both")
    parser.add_argument("--index", type=int, help="show one row under one rule, starting at 1")
    args = parser.parse_args()
    if args.length < 0 or args.colors < 1:
        parser.error("choose a length of at least 0 and at least 1 light color")
    if args.index is not None and args.rule == "both":
        parser.error("choose --rule alternating or --rule separated when using --index")

    rules = ("alternating", "separated") if args.rule == "both" else (args.rule,)
    for rule in rules:
        total = count_rows(args.length, args.colors, rule)
        if args.index is not None and not 1 <= args.index <= total:
            parser.error(f"choose an index from 1 to {total} for this puzzle")
        print(f"{rule.capitalize()} rows: {total:,}")
        if args.index is not None:
            row = row_at_index(args.length, args.colors, rule, args.index)
            print(f"Row {args.index:,}: " + (" | ".join(map(lamp_name, row)) if row else "(empty row)"))


if __name__ == "__main__":
    main()
