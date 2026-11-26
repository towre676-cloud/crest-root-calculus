# Elliptic quotient spectra under complex multiplication

[Research paper (PDF)](../../../downloads/Crest-Root-Elliptic-Quotient-Spectra.pdf) · [TeX source](paper.tex) · [Exact replay](verify.py)

Let an isogeny from a power of an elliptic curve with complex multiplication carry a finite kernel $K$. A primitive CM direction $v$ is described by its coordinate ideal $\mathfrak a$. If $H$ is the pulled-back polarization and $k(v)$ is the order of the intersection of the direction with $K$, the degree is

$$
d(v)=\frac{v^*Hv}{\mathrm{N}(\mathfrak a)\,k(v)}.
$$

The paper proves an explicit asymptotic in every dimension $r\ge 2$ and for every imaginary quadratic field. Its leading constant is the product of a Hermitian volume, the primitive ideal-lattice density, and the finite average of $k(v)^r$ over admissible residues. That average factors over the rational primes in the kernel exponent.

## Two examples

**Same rational data, different spectra.** Two explicit genus-two models use the elliptic square with $j=-3375$, each with pullback polarization $2I$. Both are isogenous to the common square over a number field and therefore share Frobenius polynomials at its common good places. One has zero degree-three optimal classes and the other has eight. Their leading constants are $3/(7L(2,\chi_{-7}))$ and $4/(7L(2,\chi_{-7}))$.

**Every degree for the Gaussian gluing.** For the non-product $j=1728$ graph gluing, let $B(m)=m\prod_{p\mid m}(1-\chi_{-4}(p)/p)$ for odd $m$. The exact geometric multiplicities are:

| Degree | Classes |
|---|---:|
| Odd $m\equiv3\bmod4$ | $2B(m)$ |
| Odd $m\equiv1\bmod4$ | $0$ |
| $2m$, $m$ odd | $2(1+\mathbf1_{m\equiv1\bmod4})B(m)$ |
| $4m$, $m$ odd | $2B(m)$ |
| Multiples of 8 | $0$ |

The count through degree $X$ is $15X^2/(32L(2,\chi_{-4}))+O(X\log(2+X))$.

## An unbounded first difference

For every positive $M\equiv9\bmod24$, Theorem 5.5 constructs two genus-two curves whose optimal elliptic quotient counts agree in every degree prime to three and every degree below $M$. The counts differ by four at $M$. The two refined Humbert forms are inequivalent, and both arise from explicit two-torsion gluings. A Pell subsequence gives arbitrarily late first differences in the fixed CM field $\mathbb Q(\sqrt{-481})$ as the conductor varies. [Lattice construction and replay](../opposite-secant-spectra/bielliptic.md)

## Run the checks

From the repository root:

~~~sh
python -m pip install sympy
python -B research/elliptic-quotients/cm-spectra/verify.py
~~~

The replay checks the two covers on each of the discriminant $-7$ models, computes the two-torsion kernel tables, and enumerates primitive slopes using exact ring arithmetic. It also checks the Gaussian divisor formula against the direct enumeration and against both independent [Néron–Severi certificates](../ns-certificates/README.md) at every degree through 150.

Kani supplies the correspondence between optimal elliptic quotient classes and primitive classes on a Jacobian. Guerra gives the projection formula and estimates for elliptic degree counts. The article uses Schanuel's height-counting method with a Hermitian norm and evaluates the finite kernel correction in the leading term. [Kani](https://doi.org/10.1007/s40316-018-0105-6) · [Guerra](https://doi.org/10.13137/2464-8728/27068) · [Schanuel](https://doi.org/10.24033/bsmf.1905)
