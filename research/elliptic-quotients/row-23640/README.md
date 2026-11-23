# All elliptic map degrees on Sutherland row 23640 and its deformations

On the smooth projective quartic in row 23640 of Andrew Sutherland’s table, every nonconstant map to an elliptic curve factors through one degree-two map. As the elliptic target varies, the geometric degrees are exactly the positive even integers. The degrees over the rationals are exactly $2n^2$ and $4n^2$ for positive integers $n$. In particular, the row has **zero geometric degree-three quotient classes**. This supplies the twelfth exclusion in the original twelve-row blind sample. The original Part XXX run accounted for eleven. [Source table](https://math.mit.edu/~drew/genus3curves.html) · [Original twelve-row input](../../../applications/atlas/example.json) · [Original run](../../../applications/atlas/original/catalogue.json)

The same spectra hold on infinitely many geometrically nonisomorphic smooth plane quartics. The family below keeps the elliptic quotient fixed while moving its four branch points. No input curve beyond the original twelve is included in the benchmark total.

## Recover the map from the equation

The [input file](quartic.json) gives the original ternary quartic $F(x,y,z)$. To find a linear involution, search primitive rational centers $P$ of bounded coordinate height with $F(P)\ne0$. For a smooth quartic preserved by a projective involution in characteristic zero, coordinates with center and axis as eigenspaces express $F$ as an even polynomial in the center coordinate. An odd quartic would contain the axis as a component. The center lies off the curve because an even quartic vanishing there would be singular there. Its axis is forced by the coefficient of the third power of that coordinate, namely $L(X)=\nabla F(P)\cdot X=0$. Test the cubic $D_PF$ on the axis and check $F(MX)$ coefficient by coefficient, with $M=2F(P)I-P\nabla F(P)^t$ (gradients are columns). Euler’s identity gives $M^2=(2F(P))^2I$. This finds every rational involution within the specified center height in those coordinates. The height changes under projective transformations.

The height-one search on the supplied equation finds the center $(1,1,-1)$ and axis $y=0$. The involution is $(x,y,z)\mapsto(x-2y,-y,z+2y)$. Set $x=u+t$, $y=t$, $z=v-t$. The source becomes

$$
3t^4+At^2+B=0,\qquad
A=u^2+7uv+2v^2,\qquad B=uv(u^2+3uv+v^2).
$$

The involution sends $t$ to $-t$. Its quotient $\rho:C\to E$ sends $[u:v:t]$ to $[u:v:w]$ with $w=6t^2+A$ and

$$
w^2=D=A^2-12B=u^4+2u^3v+17u^2v^2+16uv^3+4v^4 .
$$

Here $w$ has weight two. The discriminant of $D(u,1)$ is $57{,}600$. The double cover of the $[u:v]$-line is smooth of genus one and has rational points at infinity. The coordinate $t$ gives a quadratic extension of the invariant field, witnessed by the nontrivial sign involution. The formulas extend everywhere because $u,v$ cannot both vanish on $C$. In particular the quotient has degree two. Its invariant is $j=13997521/225$.

For a Weierstrass equation, on $v=1$ write $G=u^2+u+8$, so $D=G^2-60$. Set $H=2(w+G)$ and $Y=H(2u+1)$. Then $Y^2=H(H-15)(H-16)$. Conversely $u=(Y/H-1)/2$ and $w=H/2-G$. With $X=H-15$ the elliptic curve is $Y^2=X(X-1)(X+15)$.

## Why this quotient exhausts every degree

Let $J$ be the Jacobian of $C$ and $P$ the connected kernel of $\rho_*:J\to E$, an abelian surface called the Prym. The three affine coordinate charts of the gradient ideal certify smooth reduction modulo 23. Counting a squared coordinate in each fibre over $[u:v]$, including $[1:0]$, gives

| Field | Points on $C$ | Points on $E$ |
|---|---:|---:|
| $\mathbb F_{23}$ | 30 | 24 |
| $\mathbb F_{23^2}$ | 608 | 576 |
| $\mathbb F_{23^3}$ | 12,186 | 12,168 |

Newton identities and the functional equation give the Jacobian Frobenius polynomial $(T^2+23)h(T)$, where

$$
h(T)=T^4+6T^3+34T^2+138T+529 .
$$

The quadratic is the Frobenius polynomial of the exhibited elliptic quotient. To decide if the remaining surface could acquire an elliptic factor over an extension field, take the exact resultant

$$
\operatorname{Res}_T(h(T),h(UT))
=23^4(U-1)^4 r(U)^2s(U),
$$

$$
r=23U^4+12U^3+14U^2+12U+23,\qquad
s=529U^4+736U^3+558U^2+736U+529.
$$

All three quartics $h,r,s$ are primitive and irreducible modulo 11 with degree preserved, hence irreducible over the rationals. The nonmonic primitive polynomials $r,s$ cannot have roots of unity. It follows that distinct roots of $h$ have distinct $m$-th powers for every $m\ge1$. Those powers remain a single Galois orbit of size four. Thus the Frobenius polynomial on $P$ stays irreducible over every finite residue-field extension, so the reduction of $P$ is geometrically simple. Any elliptic subfactor in characteristic zero would specialize after finite extension at this good prime. Equivalently, a nontrivial rational idempotent for that subfactor would specialize injectively in the rational endomorphism algebra. This proves that $P$ is geometrically simple in characteristic zero. [Howe and Zhu on absolute simplicity](https://arxiv.org/abs/math/0002205)

For any map $f:C\to E'$ to an elliptic curve, the induced homomorphism kills $P$ since a simple abelian surface has no elliptic quotient. If $\sigma$ denotes the involution, the image of $1-\sigma_*$ lies in $P$ since $\rho\sigma=\rho$. Consequently $f-f\sigma$ is constant. Evaluate at the rational fixed point $[u:v:t]=[1:0:0]$ to see that the constant is zero. The quotient property gives $f=h'\circ\rho$, over the field of definition of $f$. Up to translation $h'$ is an isogeny. Over an algebraic closure an elliptic curve has cyclic subgroups of every positive order, so every even degree occurs. Up to target isomorphism and translation, there is exactly one degree-two map. More generally the number of geometric classes of degree $2m$ is $\sum_{d\mid m} d$, the number of index-$m$ sublattices of a rank-two lattice.

## Which degrees occur over the rationals?

The rational points of order two on $E$ supply a rational 2-isogeny. To exclude rational isogenies of odd prime degree $\ell$, Mazur’s theorem reduces the candidates to $3,5,7,11,13,17,19,37,43,67,163$. At a good prime $p\ne\ell$, a rational cyclic kernel would give a Frobenius-stable line in $E[\ell]$. A nonsquare $a_p^2-4p$ modulo $\ell$ excludes such a line. Direct counts on the Weierstrass model give the following covering set. [Mazur, Theorem 1](https://www.math.columbia.edu/~goldfeld/Mazur-Goldfeld1978.pdf)

| $p$ | $a_p$ | $a_p^2-4p$ | Odd candidates excluded at this step |
|---:|---:|---:|---|
| 23 | 0 | −92 | 5, 7, 11, 17, 19, 37, 43, 67 |
| 7 | 0 | −28 | 3, 13 |
| 17 | 2 | −64 | 163 |

A Galois-stable kernel with elementary divisors $a\mid b$ contains $E[a]$. Its image under $[a]$ is a Galois-stable cyclic group of order $b/a$. An odd prime dividing $b/a$ would give an odd prime isogeny. Thus every rational isogeny has degree $a^2 2^j$, a square or twice a square. All squares occur via multiplication by $n$. Composing this with a rational 2-isogeny realizes twice every square. Composing once more with $\rho$ proves the claimed rational spectrum. For example a degree-six map exists geometrically, while no degree-six map is defined over the rationals.

## An infinite family with identical degree spectra

For $k\in\mathbb Z$ define

$$
A_k=A+138k uv,\qquad
B_k=B+23k uv A+1587k^2u^2v^2,\qquad
C_k:\ 3t^4+A_kt^2+B_k=0 .
$$

Identically $A_k^2-12B_k=D$, so $w=6t^2+A_k$ gives every member the same elliptic quotient $E$. Every $C_k$ has the same smooth reduction at 23 and the same rational fixed point. Hence every member is a smooth plane quartic, and the entire single-prime exhaustion and both degree spectra proved above apply to every integer $k$.

These equations represent infinitely many geometric isomorphism classes. Let $\pi:E\to\mathbb P^1$ send $[u:v:w]$ to $[u:v]$. The branch divisor of $C_k\to E$ is the zero divisor of $w-A_k$, a section of $L=\pi^*\mathcal O(2)$ of degree four. The cover is smooth, so these branch points are simple. Distinct parameters give distinct divisors. Sections having the same divisor would be proportional, and the coefficient of $w$ forces the scalar to be one. An isomorphism between two source curves must identify their unique degree-two quotients up to an automorphism of $E$ and therefore carry one branch divisor to the other. Such an automorphism preserves $L$. Among translations its stabilizer is $E[4]$, of order 16; the automorphisms fixing the origin have order two because $j\notin\{0,1728\}$. Thus at most 32 parameters can give any one geometric isomorphism class.

## A second proof for the catalogue row

At 13 the complementary surface has ordinary elliptic factors with Frobenius polynomials $T^2+4T+13$ and $T^2-2T+13$. Their geometric rational endomorphism fields are $\mathbb Q(i)$ and $\mathbb Q(\sqrt{-3})$, so they are geometrically nonisogenous. This rules out a geometric square of an elliptic curve as the characteristic-zero Prym. If the Prym were geometrically decomposable, it would therefore have exactly two elliptic subvarieties, permuted by Galois. Each would be defined over an extension of degree at most two. At 7, however, the Frobenius polynomial of the Prym over $\mathbb F_{49}$ is the irreducible $T^4+46T^2+2401$. Good reduction then precludes either elliptic subvariety even after such an extension. This gives an independent exclusion on $C_0$ using two different primes. The family proof uses 23.

## Replay

Run these commands in this folder with Python 3.11 or later. On Windows, use `py` in place of `python` if needed.

```sh
python -m pip install "sympy>=1.12,<2"
python discover.py quartic.json --height 1
python verify.py
```

The discovery program uses the Python standard library. The [verifier](verify.py) checks the source identity, smoothness, point counts, root-ratio resultant, isogeny sieve and family identities. It also runs the [two-prime check](verify_involution.py). Compare the printed JSON to the [recorded output](receipt.json). The original bounded solver’s historical receipts retain their eleven exclusions. The combined twelve-row result uses that run and the independent proof given here.
