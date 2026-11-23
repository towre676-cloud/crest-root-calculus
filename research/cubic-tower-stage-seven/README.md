# The geometric decomposition of the seventh cubic-tower stage

The Jacobian of the stage with indices $0,\ldots,6$ has genus **10,207**. Its signed subset decomposition over $\mathbb Q$ has **255 geometrically simple factors** that are pairwise nonisogenous over $\overline{\mathbb Q}$. In particular its two 448-dimensional factors both remain simple after every finite field extension of the ground field in characteristic zero. The certificate uses their good reductions and gives the complete list of excluded proper dimensions.

| Dimension | 1 | 4 | 12 | 32 | 80 | 192 | 448 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Number of geometric simple factors | 15 | 42 | 70 | 70 | 42 | 14 | 2 |

For support size $k$, there are $2\binom7k$ signed factors of dimension $k2^{k-1}$, in addition to the elliptic base. These counts give the stated dimension after summing with multiplicity. The [original projectors](../cubic-tower-factors/README.md), [triple simplicity calculation](../cubic-tower-cartier/README.md), [four-index certificates](../cubic-tower-geometric-separation/README.md) and [five-index and six-index results](../cubic-tower-size-five/README.md) supply the lower-dimensional factors and their geometric separation.

The [present certificate](certificate.json) proves absolute simplicity of the twelve further 192-dimensional factors needed beyond stage six, and the two new 448-dimensional factors on all seven indices. At each listed ordinary prime the nonzero reciprocal Frobenius roots are powered by an exponent divisible by every possible root-of-unity order arising as a ratio of the ambient Weil roots. A geometric abelian subvariety of dimension $e$ would give a Galois-invariant multiset of $e$ unit roots on every such reduction. Thus $e$ must be a subset sum of the degrees of the powered unit roots' irreducible factors. The verifier intersects the dimension sets and requires the empty intersection for $1\le e\le\lfloor d/2\rfloor$ on each signed $d$-dimensional factor. [Full argument](../cubic-tower-geometric-separation/README.md#certifying-geometric-isogeny-information)

The fourteen signed 192-dimensional factors in this stage have distinct geometric Frobenius signatures at 193. The program checks their signatures on the corresponding entries of the published six-index receipt. The two 448-dimensional signatures differ at 449. Two ordinary primes suffice to certify the positive factor, 449 and 661. For the negative factor the two primes are 449 and 607. Different dimensions distinguish the different columns of the table. These are criteria over the algebraic closure of the rationals. The characteristic-zero proof does not require that any particular local reduction stay simple.

At 449 the two signed first-power traces on the 448-dimensional factors are 2 and 10. Their second-power traces are $-314$ and $830$. Both pairs were recomputed by direct fibre enumeration. Modulo 449 they agree with the prime-field algorithm, which extracts Frobenius on the finite affine fibres and accounts for the normalized sheets at infinity. The code assembles the coefficient matrices in two groups of indices. The coefficient products are exact integer dot products with an explicit binary64 bound. A separate assembly using the parities of the two index groups agreed in all 898 local coefficients at 449. [Shared coefficient algorithm](../cubic-tower-size-five/block_frobenius.py) · [Finite-fibre formula](../cubic-tower-geometric-separation/README.md#frobenius-on-the-finite-fibres)

From this directory, with NumPy and SymPy installed, run

```sh
python verify.py --workers 4
python verify.py --recount --workers 4
```

The first command recomputes the dimension and separation checks from the stored local polynomials, reassembles the pair at 449 and checks the direct counts. The second reassembles every listed local numerator. Use `py` in place of `python` on Windows if appropriate.
