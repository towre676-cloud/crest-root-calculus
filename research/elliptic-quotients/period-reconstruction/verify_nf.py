"""Exercise exact map recognition over a cubic number field."""

import mpmath as mp

from nf_certificate import certify_nf


def main():
    with mp.workdps(70):
        t = mp.root(2, 3)
        f = [-36, 0, 49, 0, -14, 0, 1]
        # Push the known elliptic target through X -> X/t^2, Y -> Y/t^3.
        args = ([-111/t**2, 0, mp.mpf(157)/(3*t**2)], [-9, 0, 1],
                mp.mpf(3556)/(3*t**4), -mp.mpf(192104)/(27*t**6), (t, mp.mpf(0)))
        result = certify_nf(f, *args)
        if not result["verified"] or result["map_degree"] != 2 or result["field_polynomial"] != "t^3 - 2":
            raise ArithmeticError(f"The cubic recognition control failed: {result}")
    print("Verified a degree-two map in coordinates over Q(cube root of 2).")


if __name__ == "__main__":
    main()
