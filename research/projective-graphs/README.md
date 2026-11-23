# Projective graphs with poles and base points

Given an integral projective source and homogeneous formulas for a rational map, [graph.py](graph.py) computes the saturated bihomogeneous ideal of the graph closure and its bigraded Hilbert polynomial. The [worked example](row-23640.json) uses the degree-two map from Sutherland row 23640 to its elliptic quotient. Its three coordinate polynomials have greatest common divisor one, yet all three vanish at two geometric points. The graph determines their actual images. The same seven equations with a parameter give the graphs for the entire integer deformation family.

## Construction

Let the homogeneous prime source ideal be $I\subset\mathbb Q[x_0,\ldots,x_n]$, with nonempty projective zero locus. Let $h_0,\ldots,h_m$ be homogeneous forms of equal degree, at least one nonzero modulo $I$. For target coordinates $Y_i$ and one extra variable $s$, compute

$$
K=\left(I,\ Y_0-s h_0,\ldots,Y_m-s h_m\right)\cap\mathbb Q[x_0,\ldots,x_n,Y_0,\ldots,Y_m].
$$

The program uses exact lexicographic Gröbner elimination. Equivalently $K$ is the kernel of substitution into $(\mathbb Q[x]/I)[s]$. Its quotient is the Rees algebra of the coordinate ideal with its natural grading. It is a domain, and $K$ is bihomogeneous by assigning $s$ bidegree $(-\deg h_0,1)$. In the product of projective spaces, the image of the substitution away from the base locus is the graph. This open subset is dense in the irreducible image closure. Some $x_i$ and some $Y_j$ survive modulo $K$, so its prime ideal is saturated with respect to both irrelevant ideals. This proves that the elimination returns the projective graph, including limiting points. It does not require guessing a denominator. The input contract requires an integral projective source. Runtime depends on the exact elimination.

For the Hilbert polynomial, the leading monomials generate the initial ideal. Inclusion-exclusion on their multiples gives a Hilbert numerator from the least common multiples. Dividing by $(1-S)^{n+1}(1-T)^{m+1}$ and extracting the eventual bidegree coefficient gives the polynomial exactly.

## A map whose coordinate triple vanishes

In source coordinates $[t:u:v]$ let

$$
r=t^2+23k uv,\quad A=u^2+7uv+2v^2,\quad B=uv(u^2+3uv+v^2),\qquad
C_k:\ 3r^2+Ar+B=0 .
$$

At $k=0$ this is the original catalogue row after a rational linear substitution. Every integer $k$ gives the same smooth reduction at 23. Write $w=6r+A$, $G=u^2+uv+8v^2$ and $H=2(w+G)$. The projective map to $E:Y^2Z=X(X-Z)(X+15Z)$ has the cubic coordinate tuple

$$
[X:Y:Z]=[(H-15v^2)v:H(2u+v):v^3].
$$

At $v=0$ the source equation is $t^2(3t^2+u^2)=0$. There $u$ is nonzero. For either root of $3t^2+u^2=0$ all three cubics vanish, although the tuple has no common polynomial factor. Equating pairs of coordinates by minors would admit every point of the elliptic target above either base point. For example, $[0:0:1]$ passes all those minors there.

The [recorded ideal](receipt.json) contains seven bihomogeneous generators. The three variables $t,u,v$ and the three variables $X,Y,Z$ each have degree one in their own block. The full Gröbner basis, including the eliminated variable, is monic with every leading monomial independent of $k$ when $k$ is ordered last. The defining relations specialize along with their polynomial Gröbner reductions. Thus each characteristic-zero specialization retains the basis, the eliminated ideal and the same initial monomials. For the smooth integer members the resulting ideals define prime projective graphs.

The exact Hilbert polynomial is $4m+6n-2$: the source hyperplane has intersection four with the graph and a target line has intersection six. The source quartic has degree four and the target cubic has degree three. Dividing gives projection degrees $1$ and $2$. The constant term $-2=1-3$ agrees with genus three.

## The four fixed points of negation

All these fibres lie on the chart $u=1$. The verifier sets the target point in the graph ideal, computes the affine fibre ideal, and excludes points on the complementary line using both projective charts. Here $O=[0:1:0]$ and affine points refer to the chart $Z=1$.

| Target | Fibre ideal with $u=1$ | Length |
|---|---|---:|
| $O$ | $(v,t^2)$ | 2 |
| $(-15,0)$ | $(v,t^2+1/3)$ | 2 |
| $(0,0)$ | $(v+2,t^2-46k-2/3)$ | 2 |
| $(1,0)$ | $(v+2,t^2-46k-1)$ | 2 |

The two conjugate base points both map to $(-15,0)$. The rational point $[0:1:0]$ on the source maps to $O$ with ramification index two. The four displayed fibres also check a correspondence intersection in $C_k\times E$. The graph of $f$ meets the graph of $-f$ where $[2]f=O$. In characteristic zero $E[2]$ is reduced, so the intersection scheme is the disjoint union of the four fibres and has length eight. The graph’s self-intersection is zero since its normal line bundle is $f^*T_E$ and $T_E$ is trivial.

## Replay or enter another map

Run with Python 3.11 or later and SymPy 1.12–1.x.

```sh
python graph.py row-23640.json
python verify.py
```

In the generic input, `source_coordinates` is the number of homogeneous variables. Each member of `equations` defines a homogeneous source equation and each member of `map` defines one target coordinate. A polynomial is a list of monomials of the form `[i,j,k,c]` for three variables, with one exponent per variable in other dimensions. A coefficient may be a rational string. On the supplied example the order is $t,u,v$. The family verifier checks the full symbolic ideal over $\mathbb Q[k]$, each of the four fibre ideals, the Hilbert polynomial and the false image admitted by the naive minors. [Verification program](verify.py) · [Expected result](receipt.json)
