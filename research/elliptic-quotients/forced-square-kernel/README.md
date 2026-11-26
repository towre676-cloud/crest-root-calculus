# The forced kernel of coprime square counts

The coprime-square kernel has a subspace that can be computed without any modular-form eigenvectors or $L$-values. In the split genus it has dimension of order $\ell$ and projector diagonals of order $1/\ell$ for fixed $D$. Its intersection with vectors supported on the gluing classes is zero. Thus the full kernel already has order at least $\sqrt h$, where $h$ is the split genus class number, and its dimension alone overstates the obstruction to gluing injectivity.

The split and nonsplit genera also have **identical observation rank on every common set of coprime degrees**. Their kernel dimensions differ by exactly the dimension of this forced subspace. At $D=-163,\ell=193$, the two genera have **126,488** and **125,186** classes, giving **1,302** forced split-genus directions.

## Ambient lattices and incidence relations

Let $G$ run through the integral isometry classes in the genus of $\mathbb Z\oplus O$ with the ambient norm $Q(u,c)=u^2+4N(c)$. In the local notation of the [gluing theorem](../inert-gluing/README.md), an integral lattice in either gluing genus is the inverse image of a nondegenerate plane in $G/\ell G$. Every such lattice has a unique canonical ambient overlattice of index $\ell$. Locally its discriminant group at $\ell$ is cyclic of order $\ell^2$, with just one subgroup of order $\ell$. At all other primes the ambient lattice is unchanged. The local signs classify the two possibilities. Thus in each ambient class, lattices in the split genus correspond exactly to $O(G)$-orbits of split planes, and those in the nonsplit genus to orbits of nonsplit planes.

The isotropic projective points $X$ form a conic with $n=\ell+1$ elements. Split planes are its secants, so they are indexed by unordered pairs in $X$. Let $I$ be the internal projective points, that is, points whose orthogonal planes are nonsplit. Since $(D/\ell)=-1$ and $\det G=4|D|$, a vector of norm a nonzero square lies on an internal point. For each $x\in X$ and each $v\in I$, there is a unique secant containing both. Indeed $x$ cannot be perpendicular to $v$, and their span is a nondegenerate plane containing $x$.

Consequently for every degree $m$ prime to $\ell$,
$$
 \sum_{y\in X,\ y\ne x} r_{L_{xy}}(m^2)=r_G(m^2).
$$
Here $r$ is the ordinary integral representation count and $L_{xy}$ is the inverse image of the secant $xy$. The same relation holds for the primitive counts by divisor inversion. For an $O(G)$-invariant function $f$ on $X$ of sum zero, the edge weight $Tf(xy)=f(x)+f(y)$ therefore gives a relation. On an edge orbit $\mathcal E$ the coefficient vector of that relation has entry $|\mathcal E|Tf(\mathcal E)$, up to a common constant. The map $T$ is injective for $n>2$. Taking the direct sum over ambient classes gives a forced subspace $S$ of the split kernel with
$$
 \dim S=\sum_{[G]}\bigl(\#(O(G)\backslash X)-1\bigr).
$$
For fixed $D$, the ambient groups are finite and independent of $\ell$. This gives $\dim S\asymp_D\ell$, while the number of secant orbits is $h\asymp_D\ell^2$.

## Transfer to the nonsplit genus

Let $B_{\mathrm{sp}}$ and $B_{\mathrm{ns}}$ be the incidence matrices between the respective planes and the internal points. For prime $\ell$, Wu's rank theorem gives full column rank over $\mathbb F_\ell$ for both matrices. These are the secant–internal and skew–internal blocks in that theorem, each with $\ell(\ell-1)/2$ columns. Hence both maps are injective over $\mathbb Q$, and the second is square and invertible. The maps commute with the ambient automorphisms, so these statements pass to invariant functions. [Wu, Theorem 1.3(iv)](https://arxiv.org/pdf/1002.1138)

For every $m$ prime to $\ell$, form the vector $z_G(m)$ of counts of ambient norm-$m^2$ vectors reducing to each internal point. The two observation vectors in this ambient block are $B_{\mathrm{sp}}z_G(m)$ and $B_{\mathrm{ns}}z_G(m)$. Thus the first is obtained from the second by the fixed injective map $B_{\mathrm{sp}}B_{\mathrm{ns}}^{-1}$. This works before any assumption about the actual span of the vectors $z_G(m)$.

The transpose on coefficient vectors gives the exact sequence for the kernels on any common set of degrees prime to $\ell$,
$$
 0\longrightarrow S\longrightarrow V_{\mathrm{sp}}
 \longrightarrow V_{\mathrm{ns}}\longrightarrow 0.
$$
Indeed the kernel of the transpose on the full coefficient space is exactly the incidence relations above, by injectivity of the two incidence maps and the difference of their dimensions. Passing to invariants preserves exactness. Applying the maps in all ambient blocks gives the displayed result, including for the full coprime kernels. Good-prime neighbour operators commute with the transfer. A neighbour at $p\ne\ell$ alters the ambient lattice only at $p$, where its index-$\ell$ sublattice is the same locally, and transports the plane at $\ell$ unchanged in the common rational quadratic space. The local incidence map is equivariant under these identifications.

### The forced directions cannot cause a gluing collision

The gluing classes all lie over the original ambient lattice $G=\mathbb Z\oplus O$. They correspond to planes avoiding the norm-one line $v_0=(1,0)$. The split planes through $v_0$ pair up the conic points, since $v_0$ is internal. The partner $\tau x$ of $x$ is obtained by reflecting in $v_0^\perp$, an integral ambient automorphism. Thus $f(\tau x)=f(x)$ for every invariant $f$. If a forced relation vanishes off the gluing classes, its edge value on each of these omitted pairs is $2f(x)=0$. It follows that $f=0$. The other ambient blocks also vanish, by injectivity of their incidence maps.

In particular, within either genus the residue-resolved counts already give uniform gluing separation. Here the extra label is the orbit of the internal reduction point under the original ambient automorphisms. A distinct gluing pair has different incidence on at least one such orbit by the kernel description and the zero-intersection argument. Every internal ray can be tested at a good degree $O_D(\ell^2)$. Choose a norm-one representative $(u,c)$ of the ray over $O/\ell O$, with $u\ne-1$ by changing its sign if necessary, and take $t=c/(1+u)$. Lift $t$ to $O$ with coefficients of size $O_D(\ell)$ and with $t=0\bmod 2|D|$. The integral vector $(1-4N(t),2t)$ has norm $(1+4N(t))^2$, reduces to the desired ray and has degree prime to $2D\ell$. Dividing its content preserves the ray and coprimality. For a fixed degree the ambient residue counts are constant along each orbit, so the corresponding orbit-resolved counts differ. The unmarked spectrum sums over the orbits, where cancellation remains possible.

## Concentration and the analytic question

Work in the mass inner product on functions on the lattice classes. In an ambient block with group $\Gamma=O(G)$, the mass of an edge orbit $\mathcal E$ is $s/|\Gamma|$, where $s=|\mathcal E|$. Equivalently the coefficient inner product is the inverse-mass inner product used in the supplied projector argument. Let $k_O$ count the endpoints of a representative edge in the vertex orbit $O$. The forced projector has diagonal
$$
 \rho_S(\mathcal E)=
 \frac{s}{n-2}\left(\sum_O\frac{k_O^2}{|O|}-\frac4n\right).
$$
In fact $T^*T=(n-2)I+J$ on the unquotiented vertex space. On the sum-zero space the inverse is multiplication by $1/(n-2)$. For the normalized edge-orbit indicator, $T^*$ averaged over $\Gamma$ has value $\sqrt{s}\,k_O/|O|$ on $O$. Removing its constant part gives the formula. The unrestricted diagonal on a single edge is $2/n$ by symmetry and trace, so the invariant orbit diagonal is at most $2s/n$ by Cauchy–Schwarz. Hence $\max\rho_S=O_D(1/\ell)$.

For a gluing pair put $c=e_i-e_j$ in coefficient coordinates and let $R=V_{\mathrm{sp}}\ominus S$. A sufficient condition for separation is
$$
 \frac{\|P_Sc\|^2}{\|c\|^2}+\rho_R(i)+\rho_R(j)<1.
$$
To see this, use Cauchy–Schwarz on the two normalized coordinate vectors to bound $\|P_Rc\|^2$ by $(\rho_R(i)+\rho_R(j))\|c\|^2$. The $S$ contribution is known exactly. In the nonsplit genus the same test applies with $S=0$. Thus a uniform analytic estimate could target the **residual projector** or the common residual kernel dimension $\dim V_{\mathrm{ns}}$. Counting the full split kernel as exceptional central zeros would include the explicit incidence directions.

The standard Blomer–Michel result assumes Eichler orders. At the varying prime the split stabilizer reduces to the normalizer of a split torus. Its ordered index-two subgroup fixes the length-two segment corresponding to an Eichler order of level $\ell^2$. The nonsplit stabilizer reduces to a nonsplit torus normalizer. Its torus acts transitively on the neighbours of the ambient maximal vertex, so the Eichler-level theorem does not apply there directly. Fixed bad places also need to be checked when identifying the global automorphic level. Even where the Eichler estimate applies, its squared saving of exponent less than $1/30$ in volume cannot be combined with a full kernel of size at least $h^{1/2}$ to establish the proposed rank-times-sup-norm bound. [Blomer–Michel, Theorem 1](https://www.epfl.ch/labs/tan/wp-content/uploads/2018/10/supnorms-final2.pdf)

## Exact checks

The [incidence census](star.py) enumerates the ambient genera by exact integral Gram bounds and primary discriminant forms. It then counts vertex, secant and exterior-line orbits and computes the rational diagonals. The [210-level receipt](incidence.json) records the resulting exact forced ranks. The [kernel comparison](rank.py) uses good-square counts on the full lattice genera, with upper bounds on nullity certified by ranks modulo $65\,521$ or by exact rational elimination. When the upper bound reaches the incidence lower bound, the all-degree kernel is certified exactly. The [28-level receipt](rank.json) records the sample sizes, the kernel bounds and the projector bounds where rational projectors were computed.

The supplied 12 concentration calculations at commit `b778296279fe5bc2cbb60baeb085bff83d050e94` were reproduced. The scaled experiment compares both genera at the same degree cutoff. At 16 of the 28 levels the full split kernel is exactly $S$ and the nonsplit kernel is zero. At the other 12 levels the residual dimension is bounded above by eight. For instance at $D=-8,\ell=47$ the split genus has 300 classes and 11 forced directions, with residual dimension at most two.

~~~sh
python -m pip install cypari2 numpy python-flint
study=research/elliptic-quotients/forced-square-kernel
python -B "$study/star.py" --verify "$study/incidence.json"
python -B "$study/rank.py" --verify "$study/rank.json"
~~~

Unmarked coprime-degree injectivity at arbitrary inert gluing primes remains open. The exact transfer isolates the missing arithmetic information as the span of internal-point residue-count vectors across the ambient classes. The local incidence rank and the forced relations are already determined.
