# Rational points on quartic-cover twists

Read **crest-root-twist-theorem.pdf** for the infinite-family proof, the complete classification of 122 positive squarefree parameters through 200, and the quadratic-only elliptic quotient. The editable LaTeX source is included.

## Query a twist

For a positive integer d the query refers to the smooth projective curve with affine equation y² = d(6x⁵+28x³+54x). With cypari2 installed, run

```sh
python points.py 22 55 107 118 131 166 393 2465
```

Each line reports the finite affine points and number of rational points at infinity, or that the input lies outside the proved subclasses. The program checks the prime factorization, removes square factors, applies the congruence theorem when possible and otherwise uses the certified tables. It accounts for the isomorphism between the twists of parameter n and 3n. A returned point table uses the saved arithmetic proofs. Recompute them with the commands below.

## Replay the proofs

The tested Python environment used SymPy 1.14.0, cypari2 2.2.4 with PARI 2.17.2, and python-flint 0.9.0.

```sh
python descent.py
python verify_family.py
python verify_quadratic.py
python verify_extension.py
```

The four programs check the shipped JSON receipts and compute exact identities, ranks over the rationals, torsion, local obstructions, reductions and point fibres. Two separate Magma inputs supply the number-field rank, torsion, and elliptic Chabauty certificates. Magma V2.29-10 was used with the Minkowski proof bound for class groups.

```sh
magma arithmetic.m
magma extension.m
```

The standard genus-two routine is a separate comparison. On the public Magma calculator, run a single call from comparisons.m per submission to give each curve its own 60-second budget. Recorded observations and seeds are in comparisons.txt. Successful set comparisons use completeness flags and include projective points.

The quadratic-cover input is the unmodified equation-only output of Crest–Root commit 093779c, with the original source coefficients added for reproducibility. The explicit map suffices to prove this curve’s point set, independently of geometric enumeration on any other input.
