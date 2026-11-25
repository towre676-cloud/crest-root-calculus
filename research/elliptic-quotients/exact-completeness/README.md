# Exact completeness from equations and finite torsion

For four genus-two curves from the period study, two explicit rational maps and one finite kernel suffice to determine the entire geometric optimal elliptic degree spectrum. The certificates for \(T_1,T_2\) give an independent proof of the spectra in the Richelot study. They use the equations and the two maps. A separate good-reduction certificate proves that a control curve has no geometric elliptic quotient of any degree.

## Absence at one good prime

If an elliptic subvariety exists in a characteristic-zero Jacobian, it specializes to an elliptic subvariety at every good prime, after extending the residue field if necessary. It is therefore enough to certify that a good reduction is absolutely simple.

Take an irreducible quartic Weil polynomial at a certified good prime \(p\). Its distinct roots occur in two pairs under \(\alpha\mapsto p/\alpha\). The Galois group preserves this pairing, so its order divides eight. If two roots have a root-of-unity ratio of order \(n\), the splitting field contains a primitive \(n\)-th root of unity. Thus \(\varphi(n)\mid8\). The possible orders from this bound divide \(240\). If the characteristic polynomial of \(\pi^{240}\) is irreducible, no root-of-unity ratio can occur. Since the Galois action on the four original roots is transitive, all positive powers then have four distinct conjugates. If the reduction had an elliptic factor over some finite extension, its Frobenius polynomial there would factor as a product of two quadratics. This proves the certificate. See [Howe and Zhu](https://doi.org/10.1006/jnth.2001.2697) for the general absolute-simplicity criterion.

The implementation checks squarefreeness and the leading coefficient at an odd prime before doing any point count. For a sextic it counts the zero or two rational branches at infinity using the leading coefficient. As a result it can use primes without a rational branch point. The curve \(y^2=x^6-x^4+x^2+3x-2\) is certified at \(p=5\) with Weil polynomial \(U^4+2U^3+5U^2+10U+25\), although the sextic has no root over \(\mathbb F_5\). The uploaded certificate at \(p=7\) also replays. For the split control \(C\), the prime \(5\) is rejected as bad. Failure to find a good absolutely simple reduction leaves the absence question undecided.

## Exact completeness from two maps

Start with two maps \(f_i:C\to E_i\) of degrees \(d_i\) whose pulled-back differentials are linearly independent. The code rechecks both maps by the [polynomial identity](../period-reconstruction/README.md#method) and computes their degrees from their reduced horizontal coordinates. Independence implies that \(\Psi=f_1^*+f_2^*:E_1\times E_2\to J(C)\) is an isogeny.

If the curves have equal nonintegral \(j\)-invariant, identify them geometrically by an isomorphism \(\rho\). The pullback of the principal polarization has symmetric matrix
\[
M=\begin{pmatrix}d_1&m\\m&d_2\end{pmatrix},\qquad \deg\Psi=d_1d_2-m^2 .
\]
The cross term is integral because the common elliptic curve is non-CM. At a good prime not dividing \(n=\operatorname{lcm}(d_1,d_2)\), the \(n\)-torsion of the elliptic product and of the Jacobian specializes bijectively over algebraic closures. This specialization commutes with the pullbacks. Finding \(d_1d_2\) elements of the kernel on the special fibre inside \(E_1[n]\times E_2[n]\) forces \(m=0\), saturates the degree bound, and identifies the full geometric kernel in the compatible torsion coordinates. The [finite-kernel theorem](../isotypic-gluing/README.md) then gives all geometric degrees and multiplicities using a table modulo \(n\). The sign of \(\rho\) only renames the directions.

If the elliptic curves are geometrically nonisogenous, the cross term is zero and the two factors give the only elliptic subvarieties up to the isogeny. For the non-isogeny witness the code checks good reduction of both targets and compares the traces of their \(120\)-th Frobenius powers. The splitting field of two quadratic elliptic Weil polynomials has degree dividing four, so every possible root-of-unity ratio between their roots has order dividing \(120\). Unequal powered traces therefore imply geometric non-isogeny. The two optimal degrees are \(d_i/|\ker\Psi\cap E_i|\).

For the kernel calculation the reduction is checked against the same curve identity. On a suitable translate \(R\) the pullback of \(P\) is \([f_i^{-1}(P+R)]-[f_i^{-1}(R)]\). The fibre is given by the polynomial \(N-X_0D\) and its vertical coordinate by the verified map identity. Factoring the doubling equations finds all of \(E_i[2]\) or \(E_i[4]\). The code tests every addition on the enumerated groups. The search for a suitable translate exhausts a finite field and returns an unresolved result if none exists. In particular \(C\) needs a quadratic extension at \(p=7\), where one target has only two-torsion points over the prime field.

The optimal map in a given direction also has an algebraic description. If the elliptic homomorphism for that direction factors as \(\Psi h=i\circ p\), where \(p:E\to E/K_v\) is the quotient of its finite kernel and \(i:E/K_v\hookrightarrow J\) the inclusion, the optimal map is \(i^\dagger\circ\mathrm{AJ}\). The explicit sum of the input maps equals this map followed by the dual isogeny \(\widehat p\), after centering at a base point. The implementation computes the degrees and finite kernel corrections. It does not expand rational functions for all higher-degree maps.

## Certified results

| Input | Good fields for the kernel | Kernel order | Complete geometric result |
| --- | --- | ---: | --- |
| \(C\) | \(\mathbb F_{7^2},\mathbb F_{11^2}\) | 4 | Two optimal classes, both degree 2 |
| \(C'\) | \(\mathbb F_{13^2},\mathbb F_{17^2}\) | 16 | Two optimal classes, both degree 4 |
| \(T_1\) | \(\mathbb F_{17^2},\mathbb F_{29},\mathbb F_{61}\) | 16 | \(d(a,b)=2(a^2+b^2)\) for odd \(a,b\), otherwise \(4(a^2+b^2)\) |
| \(T_2\) | \(\mathbb F_{17^2},\mathbb F_{29},\mathbb F_{61},\mathbb F_{89}\) | 16 | Same formula as \(T_1\) |

For \(T_1,T_2\), the weight table has value two precisely on the primitive odd-odd residues. There are four classes of degree four and eight each at degrees \(20,52,68,100,116,148\), the only other degrees through 150. The leading asymptotic constant is \(1/\pi\). The verifier checks these weights and spectra independently at the indicated fields. It also checks rejection of the bad prime and the deterministic fibre search at \(p=7\).

## Replay and scope

Install SymPy and Python-FLINT. Run the following commands from this directory.

~~~sh
python -B verify.py
python -B certify_spectrum.py C --prime 7 --extension-degree 2 --non-isogeny-prime 13
python -B certify_spectrum.py T2 --prime 29
python -B absence.py '[-2,3,1,0,-1,0,1]' --prime 5
~~~

The tracked receipt is refreshed by running the verifier with the **--regenerate** option. The **--input** option accepts an object with the curve coefficient array **f** and the map records **maps**. With the optional numerical dependencies in the [period study](../period-reconstruction/README.md#replay), **--coefficients** instead runs map discovery first. Its candidates are accepted only if the algebraic checks and the finite-kernel certificate pass. The presence implementation handles degree-two and degree-four input maps with full torsion over the stated field. A sextic source needs a rational branch point at the kernel prime for the odd-model Jacobian arithmetic.

For the equal-\(j\) case the implementation requires a nonintegral \(j\)-invariant and a kernel attaining the product-degree bound. The cases of distinct but geometrically isogenous targets, complex multiplication, and simple Jacobians without an absolutely simple good reduction are outside these certificates. [Costa, Mascot, Sijsling and Voight](https://doi.org/10.1090/mcom/3373) give broader endomorphism-ring methods. The contribution here is a direct finite verification for these spectra, independent of a full endomorphism-ring computation and of any Richelot presentation for the input curves.
