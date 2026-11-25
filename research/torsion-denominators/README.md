# Torsion denominators: recovering the degree that point counts cannot see

**24 November 2026.** Point counts determine a Jacobian only up to isogeny, and the minimal degree of a map from a curve to an elliptic curve is not an isogeny invariant. This study proves that this degree, away from the residue characteristic, is exactly a torsion invariant of Frobenius. It also turns that statement into a finite, independently replayable certificate. On a concrete pair the certificate separates degree two from degree four, and an explicit degree-four map realizes the larger value. A family sweep tests the certificate against an independent involution oracle. [Proofs](proof.md) · [Certificates](certificates.json) · [Sweep](sweep.json)

## The pair

$$
C:\ y^2=(x^2-1)(x^2-4)(x^2-9),\qquad C':\ w^2=11(z^2-10z+13)(2z^2+2z-7)(7z^2+18z+3).
$$

$C'$ is exactly the Richelot image of $C$ for the splitting $\{(x-1)(x-2),(x+1)(x-3),(x+2)(x+3)\}$. In that splitting $\Delta=-44$ and $H_1H_2H_3/\Delta=w^2/121$ ([`richelot.py`](richelot.py)). Hence the two Jacobians have identical Weil polynomials at every common good prime. Over $\overline{\mathbb Q}$, $C$ has the involution $x\mapsto-x$ with elliptic quotient of degree two. $C'$ has no involution, and its optimal elliptic degree is four.

## Theorem (exact torsion denominator)

At a good prime $p$, write the real Weil polynomial as $(s-a)(s-b)$ with geometrically non-isogenous elliptic blocks. Put $\varphi=F+V-b$ and $N=|a-b|$. The symmetric projector onto the $a$-block is $\varphi/(a-b)$. For every prime $\ell\ne p$, the optimal degree $d$ of the special elliptic quotient satisfies
$$
v_\ell(d)=v_\ell(N)-\max\{s\le v_\ell(N):\ \varphi(J[\ell^s])=0\}.
$$
The [polarized-reduction](../elliptic-quotients/polarized-reduction/README.md) screen proved only $d\mid N$. The equality comes from one observation: an endomorphism is divisible by $\ell^s$ exactly when it kills $J[\ell^s]$, because $[\ell^s]$ is separable. The certificate therefore has two parts. The first is a field $\mathbb F_{p^r}$ over which $J[\ell^s]$ is rational, together with generators that $\varphi$ annihilates. The second is one witness point of order dividing $\ell^{s+1}$ on which $\varphi$ is nonzero. With two coprime multiplicity-one primes, the characteristic-zero optimal degree, and therefore the minimal geometric degree, is determined exactly ([Corollary 3](proof.md#corollary-3-characteristic-zero-exact-degree)).

## Certified results

| Curve | $p$ | Block traces | $N$ | Kill levels (field degree $r$) | Witness level ($r$) | Prime-to-$p$ optimal degree |
|---|---:|---|---:|---|---|---:|
| $C$ | 23 | $0,\,-8$ | 8 | $J[2]$ ($r=1$), $J[4]$ ($r=2$) | closure $d\ge2$ | 2 |
| $C'$ | 23 | $0,\,-8$ | 8 | $J[2]$ ($r=2$) | $J[4]$ ($r=4$) | 4 |
| $C$ | 37 | $2,\,-6$ | 8 | $J[2]$ ($r=1$), $J[4]$ ($r=2$) | $J[8]$ ($r=8$) | 2 |
| $C'$ | 37 | $2,\,-6$ | 8 | $J[2]$ ($r=2$) | $J[4]$ ($r=6$) | 4 |
| $C$ | 61 | $6,\,-10$ | 16 | $J[2]$ ($r=1$), $J[4]$ ($r=2$), $J[8]$ ($r=4$) | $J[16]$ ($r=4$) | 2 |
| $C'$ | 61 | $6,\,-10$ | 16 | $J[2]$ ($r=1$), $J[4]$ ($r=2$) | $J[8]$ ($r=4$) | 4 |

Both elliptic blocks are certified separately in every row, and they agree each time, as the equal-exponent property requires. The primes 37 and 61 are coprime. The minimal degree of a map to an elliptic curve over $\overline{\mathbb Q}$ is therefore exactly $2$ for $C$ and exactly $4$ for $C'$. The prime 19 is excluded: its traces $\pm4$ make the blocks quadratic twists, so they are geometrically repeated. The $p=23$ rows reproduce the earlier local observation as full certificates. The $C$ row closes its top level by the degree-at-least-two remark, which is recorded in the certificate.

The integral module behind the table appears in [`modules.json`](modules.json). At $p=37$ the matrices of Frobenius on $J[4]$ have the same characteristic polynomial modulo $4$. There $\varphi=F+V+6$ vanishes identically on $J[4]$ for $C$ and has an image of order $4$ for $C'$.

## The degree-four map, explicitly

The Richelot correspondence sends $(z,w)$ to two points $(x_i,y_i)$ of $C$. Their images under $(x,y)\mapsto(x^2,y)$, summed on $E_1:Y^2=(u-1)(u-4)(u-9)$, give
$$
u=\frac{289z^4-120z^3-3446z^2+8797}{(5z+13)^2(z^2-10z+13)},\qquad
Y=-\frac{528(z+7)(3z-1)(z^2-5)(2z^2+2z-7)(7z^2+18z+3)}{w\,(5z+13)^3(z^2-10z+13)}.
$$
The identity $Y^2=(u-1)(u-4)(u-9)$ holds on $C'$ exactly. The map has degree $4$, as [Proposition 5](proof.md#proposition-5-the-degree-four-map) predicts from $q\hat\psi\hat\psi^\dagger q^\dagger=4$. The torsion certificate shows that no smaller map exists. [`map4.py`](map4.py) · [`map4.json`](map4.json)

## Family sweep

[`sweep.py`](sweep.py) takes the five sources $y^2=(x^2-\alpha^2)(x^2-\beta^2)(x^2-\gamma^2)$ for $(\alpha,\beta,\gamma)\in\{(1,2,3),(1,2,4),(1,3,4),(1,2,5),(2,3,5)\}$. For each source it runs every non-product rational Richelot neighbour. Each neighbour carries a degree-four map, so its optimal degree lies in $\{2,4\}$, and one odd multiplicity-one prime decides it. The certificate is compared with a 60-digit numerical oracle that counts Möbius involutions permuting the branch points. Genus-two curves have optimal elliptic degree two exactly when such an involution exists.

| Source $(\alpha,\beta,\gamma)$ | Neighbours | Certified degree 2 | Certified degree 4 | Undetermined in 60 s | No multiplicity-one prime $\le400$ | Oracle disagreements |
|---|---:|---:|---:|---:|---:|---:|
| $(1,2,3)$ | 14 | 6 | 5 | 3 | 0 | 0 |
| $(1,2,4)$ | 13 | 0 | 0 | 0 | 13 | — |
| $(1,2,5)$ | 14 | 6 | 8 | 0 | 0 | 0 |
| $(1,3,4)$ | 14 | 6 | 8 | 0 | 0 | 0 |
| $(2,3,5)$ | 14 | 6 | 8 | 0 | 0 | 0 |
| **Total** | **69** | **24** | **29** | **3** | **13** | **0** |

All 53 determinations agree with the oracle. Every certified degree-two neighbour has an involution, and no certified degree-four neighbour has one. In every case the two blocks were certified independently and agreed. The source $(1,2,4)$ is the instructive failure. No good prime below $400$ separates its elliptic blocks geometrically, which is the signature of $J\sim E^2$ over $\overline{\mathbb Q}$. Nine of its thirteen neighbours have involutions, and one has two. The single-projector certificate correctly abstains there. This is exactly the repeated-factor regime, where the theorem has to be replaced by the full commutant lattice. The three undetermined $(1,2,3)$ rows hit the time budget before a witness or a rational level was reached. The sweep stores the local outcomes, not the Mumford certificates. [`sweep.json`](sweep.json)

## Replay

```sh
pip install python-flint sympy mpmath
python verify.py               # replays every stored certificate and the explicit map
python build_certificates.py   # regenerates certificates.json (add --with-23 for the r = 4 witness, ~200 s)
python sweep.py                # family sweep, two processes
```

[`verify.py`](verify.py) does not trust the stored structure. It recounts points to recover the Weil polynomial and decodes each stored Mumford pair. It checks every kill generator's order and vanishing under $\varphi$, and checks by an explicit Cayley sweep that the generators span all $\ell^{4k}$ points of $J[\ell^k]$. It then checks each witness and recomputes the degree. The arithmetic is FLINT's exact $\mathbb F_q$ arithmetic with Cantor's algorithm. The certifier computes Sylow subgroups from complete Cayley graphs and a Smith normal form of their relation lattices. [`jac.py`](jac.py) · [`certify.py`](certify.py)

## What remains open

This settles the multiplicity-one genus-two case at primes away from $p$. The map is reconstructed here through Richelot structure, not from the torsion certificate. The general route runs through Frey–Kani gluing: the certificate identifies $n$, and the kernel of $E\times E'\to J$ is the graph of an anti-isometry $E[n]\to E'[n]$. Recovering that graph from the source equation, and from it the map, is the next theorem. The repeated-factor case $J\sim E^2$ needs the full lattice $\operatorname{End}(J)\otimes\mathbb Z_\ell$, the commutant of Frobenius on $T_\ell$, rather than a single projector. The $p$-part of $d$ needs the Dieudonné module. Higher genus replaces the linear $b$ by $\operatorname{rad}(h)/(s-a)$ with no change to the proof of the divisibility lemma.
