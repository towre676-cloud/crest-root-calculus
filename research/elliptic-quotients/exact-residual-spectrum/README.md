# Exact residual projectors and central zeros

At $D=-8,\ell=13$ the split-genus projector attains its maximum at the gluing class with scalar parameter $a=1$. The exact forced and residual contributions at that class are
$$
 \frac{79}{210}=\underbrace{\frac5{42}}_{\text{forced}}+
                \underbrace{\frac9{35}}_{\text{residual}}.
$$
The residual contribution $9/35$ is also its maximum over the full split genus. At $D=-3,\ell=17$, the corresponding split maximum is $1/4=1/12+1/6$, attained at a class with a norm-one vector outside the gluing classes. These values are for the true coprime kernels.

| $(D,\ell)$ | Split kernel dimension | Forced dimension | Residual dimension | Maximum nonsplit diagonal |
| --- | ---: | ---: | ---: | ---: |
| $(-8,13)$ | 5 | 3 | 2 | $11/70$ |
| $(-3,17)$ | 2 | 1 | 1 | $1/6$ |

All three residual eigensystems in the two nonsplit genera have exact twisted central zeros and root number $+1$. Thus each vanishes to order at least two at the center.

## A lower Sturm bound for the good coefficients

Let $H=\mathcal S_1\theta$ be the square-class-one Shimura lift of a combination of Humbert theta series, with theta level $M=4|D|\ell^2$. The lift has weight two, trivial character and level $N=M/2$. For every prime $r$, let $H^{(r)}$ denote the form obtained by deleting coefficients divisible by $r$, including the constant. At $r\mid N$ it is $H-V_rU_rH$. For a trivial-character form of weight two the depleted form belongs to level $\operatorname{lcm}(N,r^2)$. [Shimura lifting in weight $3/2$](https://doi.org/10.1017/S0027763000020468)

To verify the sharper level, decompose into primitive cusp forms and Eisenstein systems with their old dilations. For a primitive cusp form unramified at $r$, use $f-a_r f(rz)+r f(r^2z)$. At local conductor exponent one use $f-a_r f(rz)$. At exponent at least two, the trivial-character local Euler factor is $1$, so the sequence already vanishes at multiples of $r$. A dilation by a multiple of $r$ disappears under depletion. Other dilations commute with it. For a primitive Eisenstein system of trivial character the two character conductors coincide. When they are ramified at $r$, the coefficients at multiples of $r$ vanish. Otherwise the two local Euler roots give the same quadratic depletion. For the exceptional weight-two Eisenstein series $E_2$, the depleted combination is
$$
 E_2(z)-(r+1)E_2(rz)+rE_2(r^2z).
$$
It is modular at level $r^2$ by taking the difference of $E_2(z)-rE_2(rz)$ with its dilation. This treats the modular differences $E_2(z)-d E_2(dz)$ as well.

For the even discriminants $-4,-8$, deplete at $2$ and $\ell$. Both prime squares already divide $N$. At odd discriminants in the class-number-one list, deplete at the same two primes and use the [identity $U_{|D|}H=H$](../good-degree-spectrum/README.md#odd-discriminants-force-steinberg-at-the-ramified-prime). Depletion at $2,\ell$ commutes with the latter operator. Thus every coefficient whose index is divisible by $|D|$ agrees with the coefficient at the index after removing this prime power. The latter index is smaller and coprime to $|D|$.

Consequently the exact good-square kernel is already the kernel of the good-square coefficient matrix through the weight-two Sturm bound at
$$
 N^\sharp=\operatorname{lcm}(2|D|\ell^2,4).
$$
The triangular divisor relation for the lift changes neither vanishing through the bound at good indices nor vanishing at every good index. The bound is **728** at $(-8,13)$ and **1,224** at $(-3,17)$. The [exact projector replay](exact.py) counts through these bounds and computes the rational mass projectors. The original counts through degree 300 happen to give exactly the same kernels at both levels. [Sturm bound](https://www.lmfdb.org/knowledge/show/cmf.sturm_bound)

## Transporting the residual projector

Use function coordinates with the mass inner product, separately in each ambient class. Let $B$ be the incidence map from internal points to secants. Write $C$ for incidence with the nonsplit planes, identifying a nonsplit plane with its internal polar point. Then $C$ is self-adjoint. With $J$ denoting summation against the constant function on the internal points, the full projective plane incidence identity gives
$$
 B^*B+C^2=\ell I+J.
$$
Indeed any two distinct internal points determine one line, every internal point lies on $\ell+1$ lines, and no tangent contains an internal point. The calculation passes to invariant functions with orbit-size weights. [Finite conic incidence](https://arxiv.org/pdf/1002.1138)

Put $A=B^*B$. The count transfer from the nonsplit genus to the split genus is $T=BC^{-1}$. Its adjoint on the relation side has the minimum-norm right inverse
$$
 U=T(T^*T)^{-1}=BA^{-1}C.
$$
In particular, if $V_{\rm ns}$ is expressed in mass function coordinates, the exact split residual is $U V_{\rm ns}$. The forced projector is $I-BA^{-1}B^*$. These formulas specify both projectors directly, using only the nonsplit arithmetic kernel and local incidence matrices. To convert a function to a coefficient vector of a theta relation, multiply each coordinate by its class mass.

The [transport verifier](transport.py) compares the full rational projectors, entry by entry, with the independently computed split kernel. They agree at $(-8,13)$ and $(-3,17)$. The zero-residual control $(-11,7)$ checks the direct sum over two ambient classes. The nontrivial local normal matrix $A$ is part of the transfer. Estimates on nonsplit diagonals therefore need a transport bound before they yield split residual diagonal estimates.

## Exact central zeros in the residual

The two-dimensional nonsplit kernel at $(-8,13)$ splits into two rational neighbour eigenspaces with eigenvalues $-2$ and $1$ at $p=3$. At $(-3,17)$ there is one eigenspace, with eigenvalue $3$ at $p=5$. Their theta combinations are cuspidal, and their square-class-$5$ Shimura lifts are nonzero.

| $(D,\ell)$ | Cremona class of newform | Conductor after twisting by $D$ | Root number of twist | $L(1)/\Omega^+$ of twist |
| --- | --- | ---: | ---: | ---: |
| $(-8,13)$ | $676\mathrm b$ | 10816 | $+1$ | $0$ |
| $(-8,13)$ | $338\mathrm c$ | 10816 | $+1$ | $0$ |
| $(-3,17)$ | $867\mathrm e$ | 2601 | $+1$ | $0$ |

The [packet verifier](packets.py) forms the exact neighbour eigenspaces on the all-coprime kernel and identifies the square-class-$5$ lifts as sums of dilations of the elliptic newforms. It checks their full positive coefficient expansions through the Sturm bound at the original lift level. The Weierstrass models provide exact $a_m$ for those comparisons. The twisted root numbers are computed from the local root numbers. The modular symbol on the path from infinity to zero gives the exact rational central value divided by the real period. Thus the verified zeros occur in the even-sign sector, where they are not forced by the sign. This establishes the vanishing for these particular eigensystems. It leaves the uniform frequency and concentration of such eigensystems open. [Cremona elliptic curve data](https://github.com/JohnCremona/ecdata)

## Reproduction

~~~sh
python -m pip install cypari2 numpy
study=research/elliptic-quotients/exact-residual-spectrum
python -B "$study/exact.py" --verify "$study/projectors.json"
python -B "$study/transport.py" | diff - "$study/transport.json"
python -B "$study/packets.py" | diff - "$study/packets.json"
~~~

The exact-symbol replay at conductor 10816 can require more than a gigabyte of PARI stack. The [projector receipt](projectors.json), [transfer receipt](transport.json) and [central-value receipt](packets.json) record the exact results. The supplied comparison is commit `a055a0042ba135662940b1db6b3fc1a58c032c43`.

The [rank and pairwise certificates](../residual-rank-certificate/README.md) give rigorous nonzero second-derivative intervals for the three even-sign twists. They also bound the kernel's projection on 223,958 same-genus gluing-pair differences.
