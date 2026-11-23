# A reproduced functional-decomposition failure in SymPy 1.14

The polynomial-map comparison exposed an incorrect answer from stock SymPy 1.14.0's functional-decomposition routine. An exact nontrivial composition was reported as undecomposed. This is a correctness result about the specified routine and version. It supplies no speed advantage for the incorrect output. The comparison retains stock SymPy and separately labels the corrected recurrence as `corrected-dense-decomposition`. The installed SymPy package is unchanged. [Exact evidence](sympy-regression.json) · [Reproduction](verify_sympy_regression.py)

## A small squarefree counterexample

Take

$$F(x)=x^6+4x^5+4x^4+1=(x^3+2x^2)^2+1.$$

This polynomial is squarefree: its greatest common divisor with its derivative is one. Nevertheless the following public API call, executed with SymPy 1.14.0, returns the singleton list containing $F$.

```python
import sympy as s
x = s.Symbol("x")
F = x**6 + 4*x**5 + 4*x**4 + 1
print(s.decompose(F, domain=s.QQ))
# [x**6 + 4*x**5 + 4*x**4 + 1]
```

The valid decomposition is the pair $[x^2+1,\;x^3+2x^2]$. The fixed-degree helper `_dup_right_decompose` instead produces $x^3+2x^2+2x$. Its companion left-decomposition routine cannot divide $F$ into powers of that candidate and returns `None`. The benchmark adapter originally interpreted this result as absence of a map. The exact output comparison detected that interpretation as incorrect. The public API exhibits the same missed decomposition on this example, independently of the benchmark adapter. [SymPy public API documentation](https://docs.sympy.org/latest/modules/polys/reference.html#sympy.polys.polytools.decompose) · [Version 1.14 source](https://github.com/sympy/sympy/blob/1.14.0/sympy/polys/densetools.py)

## The input that exposed the failure

The frozen panel's `g2-d3-twist` case has target genus two and degree-three horizontal map. Define

$$G(u)=u^5-4u^4-4u-4,\qquad U(x)=x^3-2x^2+2,$$

$$F(x)=\frac89 G(U(x)).$$

The panel contains the expanded coefficients of this $F$, a degree-fifteen squarefree polynomial. Its source curve has genus seven. The exact identity

$$G(U(x))=\frac98 F(x)$$

gives the two geometric maps

$$u=U(x),\qquad v=\pm\frac{3\sqrt2}{4}\,y.$$

Neither vertical multiplier is rational. These are all maps in the constant-vertical comparison class. To verify exhaustion, put

$$H=x^3-2x^2,\qquad K=x^5+6x^4+8x^3-16x^2-52x-44.$$

Then $F/\operatorname{lc}(F)=K(H)$. Any normalized degree-three inner polynomial is $H$, so every possible horizontal map is $aH+b$. The next coefficient forces $b=(6a+4)/5$. The nonzero-scale coefficient equations in $G(aX+b)-a^5K(X)$ have monic gcd $a-1$ after removing powers of $a$. Thus $a=1$, $b=2$, and the vertical equation is $c^2=9/8$. The verifier checks these identities and the gcd exactly. The correct output is **two geometric maps and zero rational maps**.

Stock SymPy's fixed-degree candidate is $x^3-2x^2+2x$, so it returns no decomposition. Both Crest–Root's independent coefficient calculation and Singular's full coefficient ideal give the correct count. The raw Singular equations reduce to the stated $U$ and $8c^2-9=0$.

## The coefficient recurrence

Let $f=K(H)$, where $H$ is monic of degree $d$, $K$ has degree $m$, and $n=md$. Normalize the constant coefficient of $H$ to zero. Write $f_k$ and $h_k$ for coefficients of the indicated powers.

The polynomial

$$f'H-mfH'=H'\bigl(K'(H)H-mK(H)\bigr)$$

has degree at most $n-1$, since the degree-$m$ terms in the bracket cancel. For $1\le i<d$, its coefficient of $x^{n+d-i-1}$ must therefore vanish. That coefficient is

$$\sum_{j=0}^{i}\bigl(-i+(m+1)j\bigr)f_{n-i+j}h_{d-j}.$$

Isolating the $j=i$ term gives the exact recurrence

$$h_{d-i}=\frac{\displaystyle\sum_{j=0}^{i-1}\bigl(i-(m+1)j\bigr)f_{n-i+j}h_{d-j}}{imf_n}.$$

The inspected SymPy helper uses the factor $i-mj$. It is missing the additional $-j$ term. The separately named corrected adapter uses the displayed recurrence, followed by SymPy's existing left-division routine and the shared exact map extraction. The independent regression verifier checks this recurrence on both witnesses. This is a correction within classical functional decomposition; the general algorithmic subject is established. [Kozen and Landau, *Polynomial decomposition algorithms*, 1989](https://www.cs.cornell.edu/kozen/Papers/poly.pdf)

## Version and source evidence, checked 23 November 2026

SymPy 1.14.0 was the latest stable release listed by both the official GitHub release endpoint and the project's PyPI release history at the time of this check. It was published on 27 April 2025. PyPI also listed **1.15.0rc3**, published on 17 November 2026, explicitly as a prerelease. [Official stable release](https://github.com/sympy/sympy/releases/tag/1.14.0) · [Latest release endpoint](https://api.github.com/repos/sympy/sympy/releases/latest) · [PyPI release history](https://pypi.org/project/sympy/) · [Prerelease and its build provenance](https://pypi.org/project/sympy/1.15.0rc3/)

The installed `sympy/polys/densetools.py` is byte-identical to the official `1.14.0` file. Its SHA-256 is

```text
034675135e9e571e7bef3c3ce26f0a6604f9a3a75b0d75e89d6f13a9f2122596
```

We also inspected two immutable upstream source snapshots: the prerelease publishing commit `3ccb2c067b651f84d35321d4f950fdb351910843`, identified by PyPI provenance, and the then-current master commit `e3aef0df27a0676abeb7dc2e13abef26e4c0cf64`, dated 23 November 2026. Their `densetools.py` files were byte-identical to each other, with SHA-256

```text
e0d5e87d3d71edd3d54a5eebe4f6c4f5c7e911d51adc1a5e6b78d836a7e977ef
```

Both source snapshots still contain the same $i-mj$ recurrence factor. This is a source inspection, not a timed or installed execution of the prerelease or development versions. The measurements and public-API counterexample in this study concern stock 1.14.0. [Prerelease publishing source](https://github.com/sympy/sympy/blob/3ccb2c067b651f84d35321d4f950fdb351910843/sympy/polys/densetools.py) · [Inspected master source](https://github.com/sympy/sympy/blob/e3aef0df27a0676abeb7dc2e13abef26e4c0cf64/sympy/polys/densetools.py)

## Replay and interpretation

Run the archived version-specific check from the repository root:

```sh
python research/solver-comparison/verify_sympy_regression.py
```

The check verifies both exact composition witnesses, squarefreeness, the complete residual affine calculation, the installed module fingerprint, and the recorded stock behavior. To inspect another installed version without requiring its behavior or fingerprint to match 1.14, run

```sh
python research/solver-comparison/verify_sympy_regression.py --observe
```

The initial map study stopped at the first output disagreement and was preserved separately. The replacement study records incorrect answers explicitly. Its timing ratios use correct completed outputs only. The corrected dense recurrence has its own method name, preserving the distinction between the published library and the independently corrected comparison implementation.
