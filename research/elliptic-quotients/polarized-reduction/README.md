# Ruling out elliptic maps with two nonsimple reductions

Every quartic obtained by adding **35 times any integral ternary quartic** to the [example](quartic.json) has an absolutely simple Jacobian over the algebraic closure of the rationals. Neither reduction used in the proof has a simple Jacobian. Each has exactly one elliptic factor over its algebraic closure. The polarization forces incompatible elliptic map degrees at the two primes. This gives an unrestricted congruence class in all fifteen quartic coefficients.

| Prime | Point counts in degrees 1, 2, 3 | Real Weil polynomial | Elliptic trace | Projector multiplier | Possible optimal degrees on the special fibre |
| --- | --- | --- | ---: | ---: | --- |
| 5 | 8, 40, 140 | $s(s^2+2s-6)$ | 0 | $-6$ | 2, 3, 6 from the multiplier; the known double cover gives 2 |
| 7 | 7, 43, 352 | $(s-5)(s^2+4s-4)$ | 5 | 41 | 41 |

At 7 the special fibre has an optimal elliptic map of degree exactly 41. Every elliptic map on any integral lift of this reduction, if one exists, has degree divisible by 41. At 5 the prime-to-5 part of any optimal degree divides 6. A lift satisfying both congruences consequently has no elliptic map in any degree. In dimension three, an abelian variety with a proper nonzero abelian subvariety has an elliptic isogeny factor by Poincaré reducibility. This proves absolute simplicity. The smoothness certificates apply to every lift because both reductions stay fixed. The reduction at 5 is also that of Sutherland row 36518, whose elliptic map was established separately. [Degree-two map](../row-36518/README.md) · [Exact computations](receipt.json)

## Why polarization constrains degree

Work over an algebraically closed field. An optimal map from a curve to an elliptic curve $E$ induces a quotient $q:J_C\to E$ with connected abelian kernel. Any nonconstant elliptic map factors, up to translation, through such a quotient followed by an isogeny. Identify the Jacobian and $E$ with their duals by their principal polarizations. The adjoint $i=q^\dagger:E\to J_C$ is a closed immersion. If $d$ is the degree of the curve map, $qi=[d]$ and the symmetric rational projector onto $i(E)$ is $e=iq/d$. Conversely, each elliptic subvariety gives such an optimal quotient, with $d$ the degree of the restricted polarization line bundle. Its class is unique up to translation and target isomorphism.

The least positive integer $m$ with $me$ integral equals $d$, in every characteristic. Indeed $me$ factors as $i h$ for a homomorphism $h:J_C\to E$, since its image lies in $i(E)$. Then $[d]h=[m]q$. On the connected abelian kernel of $q$, the image of $h$ lies in the finite group scheme $E[d]$, hence is trivial. Thus $h$ descends to an endomorphism of $E$ acting as the rational scalar $m/d$. The rational scalars in the endomorphism ring are the integers, giving $d\mid m$. Also $de=iq$ is integral. In particular $d\ge2$ for a smooth genus-three curve, since a degree-one map of smooth projective curves is an isomorphism.

On a good fibre over $\mathbb F_p$, write the Weil polynomial as $P(T)=T^3h(T+p/T)$. Let $a$ be the trace of a geometric elliptic isogeny class of multiplicity one, and $b(s)=\operatorname{rad}(h)(s)/(s-a)$, using the monic squarefree part. Frobenius $F$ and Verschiebung $V$ are integral, adjoint endomorphisms. The semisimple operator $F+V$ has eigenvalue $a$ on this elliptic block and is annihilated by $b$ on its complement. It follows that $b(F+V)=b(a)e$ on the special Jacobian. The special optimal degree therefore divides $|b(a)|$. This uses the real Weil polynomial form of the gluing-exponent bound. The geometric multiplicity check is essential when different Frobenius factors become isogenous after scalar extension. [Howe and Lauter, Proposition 2.8](https://arxiv.org/abs/1202.6308) · [Geometric component computation](../../../applications/atlas/original/vendor/crest/scan_completeness_xxx/torsion/frobenius_fields.py)

For a characteristic-zero optimal quotient of degree $n$, extend the field of definition of its projector and specialize the endomorphism at the chosen good prime. For every $\ell\ne p$, smooth proper base change identifies the integral Tate lattices and the rational projectors. On either fibre, $i$ is saturated and $q$ is surjective on $T_\ell$. Thus $v_\ell(d)$ is the least $k\ge0$ such that $\ell^k e$ preserves the lattice. The specialization has rank two and its image is an elliptic subvariety. If its geometric block has multiplicity one, the symmetric projector is the unique projector onto that block, so its special degree $d$ satisfies
$$
n/p^{v_p(n)}=d/p^{v_p(d)}.
$$

Distinct characteristic-zero optimal classes give distinct symmetric projectors and remain distinct under injective specialization of endomorphisms. Hence the number of compatible multiplicity-one blocks bounds the number of geometric optimal classes. For smooth plane quartics in characteristic zero, every optimal elliptic map is primitive as a map of curves. Riemann–Hurwitz forces any intermediate curve of genus at least two to be a genus-two unramified double quotient. Its deck involution would fix a line in the canonical plane, and that line meets the quartic. [Geometric primitivity argument](../primitive-degree-four/README.md)

## Run the screen

From this directory run the following with Python and SymPy.

```sh
python screen.py quartic.json --primes 5 7 --degree 4
python verify.py
```

Input terms have the form $[i,j,k,c]$ for $c x^iy^jz^k$ with $i+j+k=4$. The program certifies smoothness using the homogeneous Jacobian ideal, counts projective points over the first three field extensions and uses cyclotomic root ratios to identify all geometric elliptic blocks. The point counter enumerates projective points, so small good primes are practical. Primes greater than three are supported by the component verifier.

A class bound of zero excludes the requested primitive degree. A positive bound is only an upper bound. A repeated geometric elliptic factor leaves the local class bound undetermined, shown by JSON `null`. With two distinct primes giving finite bounds, the output lists every degree that survives their necessary conditions. An empty list certifies absence in all degrees. The regression verifier also checks the known degree-four example with a repeated elliptic factor, rejects singular reduction and tests the residue-characteristic allowance explicitly. [Program](screen.py) · [Verifier](verify.py)
