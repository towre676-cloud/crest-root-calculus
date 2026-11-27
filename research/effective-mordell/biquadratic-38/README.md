# Exactly 38 points on the biquadratic genus-two curve

The smooth projective curve y² = x⁶ − 1 has exactly 38 points over B = ℚ(i, √7). The five-page proof combines the original Gaussian and number-field descents with a relative elliptic Chabauty calculation at 37 and an auxiliary parity sieve at 13 and 109.

| File | Contents |
|---|---|
| crest-root-biquadratic-38-proof.pdf | Complete argument for the elliptic projection theorem and its consequences |
| biquadratic-38-points.json | All 36 affine coordinates in the power basis of w⁴ − 12w² + 64, and both points at infinity |
| localization-certificate.json | The index tests at 3 and 31, ten residue classes, two formal-log columns in four embeddings, four Jacobian matrices, tangent-cone and parity checks |
| biquadratic37.py | Exact arithmetic generating the local certificate |
| points.py | The complete point classifier for the default model or any isomorphic monic integral quartic field model |
| higher-genus-points.json | Exact lists over the same field for a genus-five and a genus-55 cover |
| prior/ | Gaussian descent, quadratic divisors, cover fibres and the earlier proof of the reduction |
| upstream/ | Claude's original elliptic data, source programs, archived search receipt, and two exact check drivers |
| field-classifier.patch | Updates the previously delivered general field classifier to return complete for B |
| replay.txt | Output of the integrated proof replay |

## Replay the proof

The verified environment uses Python 3.12.14, cypari2 2.2.4 and PARI 2.17.2. In a PARI-capable Python environment, install the dependency in requirements.txt. Run with Python assertions enabled.

~~~sh
python -m pip install -r requirements.txt
python verify.py
python points.py
python higher_genus.py
~~~

All finite computations for the 38-point theorem run with PARI and Python. The proof uses exact field arithmetic for the elliptic points and evaluates their local parameters modulo 37². Each ordinary-disk matrix is also checked by translating a global representative and taking exact finite differences. The optional upstream search through a box in the Mordell–Weil subgroup is historical evidence. The verifier proves completeness without invoking that search or requiring Magma.

The genus-five cover is Y² = X¹² − 1 and has 24 points over B. The smooth plane curve (X² + X/2)⁶ − (Y⁶ + Y/2)² = 1 has genus 55 and four points over B. The exact fibres on these covers include the normalized points above infinity.

## Result for the elliptic reconstruction

The two unrestricted Gaussian projection sets are {∞, 0, 1, −9} and {∞, 0, −1, 9}. The original two quartics have exactly sixteen affine triples. Their first coordinate is ±1 or ±i, and each of the other two coordinates is independently ±2. The projection problem is stronger than the required test for points lifting to both quartics. Both are settled here.

The earlier handoff describes the squareclasses as exchanged by Gaussian conjugation. Multiplication of s by i exchanges them. On the elliptic curves the automorphism fixing ℚ(i) and negating √−3 carries the first model to the second and changes the sign of the Gaussian projection. Pages 2 and 5 give the two actions used in the proof. Descriptions in the historical documents of the field as open refer to the state before this proof.

To update the general field-census ZIP from the preceding result, run the following in the directory containing its fixed_curve.py, using the actual path to the supplied patch.

~~~sh
patch -p1 < field-classifier.patch
python solve_field.py '[64,0,-12,0,1]'
~~~

The standalone points.py already returns the proved status. This theorem resolves the indicated field and its displayed covers. The full Effective Mordell problem still requires termination for arbitrary fields and arbitrary curves of genus at least two.
