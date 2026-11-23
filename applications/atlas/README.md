# An atlas of elliptic maps

**Does a curve hide a simpler curve inside it?** A degree-three map groups the points of a complicated curve into triples over a simpler one, except at special values. The default atlas either constructs a map, proves that no elliptic map of any degree is possible, or identifies a case for further work. It reads integer coefficients, runs offline with Python 3.9 or later and returns reproducible certificates.

Run these in this folder. On Windows, `py` can take the place of `python`.

```sh
python app.py
python app.py my-curves.json
python -m unittest -v
python select_slice.py
```

Only the last command needs internet. It verifies the pinned full-source hash and all 12 selected rows, or accepts the downloaded text file as its argument. The first command reads [example.json](example.json). The [saved report](receipt.json) contains the data for each proof. For the archived Part XXX solver and its separate certification procedure, see [the original engine](original/README.md).

## Twelve real curves, selected before scoring

Andrew Sutherland published 82,240 nonhyperelliptic genus-three curves over the rationals of bounded discriminant. Our 12 rows are the ones with the smallest SHA-256 values when hashing the ASCII text `Crest-Root atlas v1|i`, where `i` is the one-based source row. That rule sees only row numbers. The Klein quartic was added separately as a positive control. The full source SHA-256 starts with `d3f06d9a` and is recorded in the input. [Source data and citation](https://math.mit.edu/~drew/genus3curves.html)

| Population | Result of the default engine |
| --- | --- |
| 12 blind published rows | Seven certificates excluding all elliptic maps, five open |
| Klein control | Degree-three map over the rationals, with a rational origin and target invariant −3375 |

This small sample gives a reproducible demonstration. It does not measure the frequency of degree-three maps in the whole database.

## Enter a quartic

```json
{"primes":[2,3,5],"curves":[{"id":"Klein","quartic":"x^3*y+y^3*z+z^3*x"}]}
```

Use up to 32 distinct curve names. Each monomial has degree four in `x,y,z`. Write products using `*` and powers using `^`. The order of the optional primes determines the order of certificate attempts. Malformed inputs exit with status 2. An unsuccessful proof attempt returns `unresolved`. All coefficients of the target and Frobenius polynomials appear in ascending order.

### How to read a positive

The coordinate permutation $(x,y,z)\mapsto(y,z,x)$ generates a group of order three. Its affine invariants include $U=(xy+yz+zx)/(x+y+z)^2$ and $V=(x-y)(y-z)(z-x)/(x+y+z)^3$. For qualifying cyclic quartics the app gives an explicit target $Y^2=H(U)$, with $Y=2LV-S(U)$. It checks smoothness at a prime of good reduction, proves the invariant decomposition by polynomial interpolation, and exhibits squarefreeness of the target at an auxiliary prime. A rational source point supplies an origin when found. Without one, the output calls the target a genus-one curve.

For the Klein control, $L=28$, $S=4+36U-108U^2$ and $H=16-160U+432U^2-160U^3-432U^4$. The point $[0:0:1]$ maps to $(0,-4)$. To see the degree, normalize $x+y+z=1$. The curve identity determines $W=xyz$ uniquely from $U,V$. The coordinates solve $T^3-T^2+UT-W=0$. A generic point has three distinct roots and six orderings, of which precisely the three cyclic orderings have the specified sign of the Vandermonde. The rational map extends over the smooth projective models. The default constructor checks a visible cyclic permutation in the supplied coordinates and requires its coefficient of $W$ to be nonzero.

### How to read an exclusion

1. The program searches for full span of the homogeneous degree 7 or 9 piece of $(F,F_x,F_y,F_z)$ over the prime field. A full span puts every monomial of that degree in the ideal. No projective point over the algebraic closure can then be singular.
2. Exhaustive projective point counts over the first three extension fields give the three independent coefficients of the degree-six Frobenius characteristic polynomial by Newton identities and reciprocity.
3. An auxiliary prime witnesses irreducibility of that polynomial by the finite-field degree-six Rabin test. The ratios of two distinct characteristic-zero roots have at most $6\cdot5=30$ conjugates. A ratio that is a root of unity of order $m$ must therefore have $\varphi(m)\le30$. All 56 possible orders larger than one are tested. The finite enumeration to 1800 follows from $\varphi(m)\ge\sqrt{m/2}$, obtained by multiplying the elementary prime-power inequalities.
4. In the degree-six residue field, let $a$ be a root. If $a^m$ belongs to neither the degree-two nor degree-three subfield then its six conjugates are distinct. This supplies an exact modular witness that no two original roots have the same $m$-th power. Irreducibility in characteristic zero makes the roots a single Galois orbit. Their powers remain in a single six-element orbit over every finite extension, so all power characteristic polynomials remain irreducible.
5. An elliptic subvariety of the Jacobian over an algebraic closure specializes, after finite extension at good reduction, to an elliptic factor over a finite extension of the residue field. The irreducible six-dimensional Frobenius polynomial on that extension excludes such a factor. A nonconstant elliptic map induces an elliptic quotient of the Jacobian. The conclusion covers all degrees.

If the auxiliary test cannot produce a witness it makes no exclusion. A partial prime search and failure to recognize cyclic coordinates are common reasons for `unresolved`. The original §86 chart solve has different coverage. [Simplicity and extension fields](https://arxiv.org/abs/math/0002205)
