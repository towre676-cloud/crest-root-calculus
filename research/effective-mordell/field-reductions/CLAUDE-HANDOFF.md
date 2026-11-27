# The next proof obligation over the Gaussian rationals

Determine all extra points of y² = x⁶ − 1 over B = ℚ(i,√7). The current lower bound is 38. Exactly 14 have rational abscissa when the two points at infinity are included. Exactly 24 have abscissa outside k = ℚ(i). The coordinates and the exhaustive divisor calculation for the latter are included.

After the local and elliptic tests in gaussian_descent.py, an additional point exists if and only if there are s, A, D in k with

A² = s⁴ + 3, 7D² = s⁴ + 27, s² ∉ ℚ.

The maps back are x = (3+s²)/(3−s²), v = 42sAD/(3−s²)³, and y = v/√7. Restore the reflection in x as well. A proof that every Gaussian solution has rational s² finishes the 38-point list. If other solutions exist, compute their induced points and certify that all possibilities have been enumerated.

Start from the rational elliptic Chabauty certificate on page 4 of the note. Try the descent of s² + √−3 and s² + 3√−3 over k(√−3), keeping track of its additional units, all prime ideals above 2, 3 and 7, and the relative norm conditions. Compute exact ranks and index information before applying relative elliptic Chabauty or a Mordell–Weil sieve. Over this base, the rational projection condition in sqrt7.m becomes a Gaussian projection condition and needs its own justification.

The normalized genus-five residual cover has five elliptic factors, each of rank two over k. The original biquadratic curve also has positive rank on every elliptic quotient over every finite extension of B. The attached code includes both rank checks. Increasing a coordinate search can produce points, while an exhaustive result requires a stopping argument for every surviving class.

An exact list over B would be one more solved field for the fixed-curve reduction. An arbitrary-field termination proof remains the next structural goal for hyperelliptic and superelliptic curves.
