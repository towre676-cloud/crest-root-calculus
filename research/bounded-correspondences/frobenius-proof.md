# Excluding maps before coefficient reconstruction

For smooth projective curves $C,D$ of positive genus over $\mathbb Q$, any nonconstant rational morphism $f:C\to D$ of degree $d$ induces maps of Jacobians with

$$
f_*f^*=[d],\qquad e_f=\frac1d f^*f_*,\qquad e_f^2=e_f.
$$

For any prime $\ell$, the corresponding rational Tate modules split into the image of $f^*$ and a complementary Frobenius-stable space. At a prime $p$ of good reduction for both curves, taking characteristic polynomials therefore gives

$$
P_{D,p}(T)\mid P_{C,p}(T),\qquad
P_{C,p}(T)=\det(1-T\operatorname{Frob}_p\mid H^1(C)).
$$

This implication holds in every degree, including degrees divisible by $p$, by using $\ell\ne p$ and rational Tate modules. A failed divisibility test at one good prime proves that no nonconstant rational map exists. The conclusion applies even when the map's formula has bad denominators at that prime, because the Jacobian homomorphism extends between the good-reduction models. A successful test supplies a necessary condition and lets reconstruction continue. [Milne, *Jacobian Varieties*, §6](https://www.jmilne.org/math/articles/1986c.pdf) · [Kedlaya, good reduction and Frobenius, §9.2](https://kskedlaya.org/weil-cohom/chapter-9.html)

For an odd-degree hyperelliptic model $y^2=F(x)$ with squarefree reduction, the counter uses

$$
N_r=p^r+1+\sum_{x\in\mathbb F_{p^r}}\chi(F(x)).
$$

The quadratic character is zero at zero, one at a nonzero square and minus one otherwise. In even degree, add the quadratic character of the leading coefficient to include the correct infinity fibre. The first $g$ counts give the power traces $s_r=p^r+1-N_r$. Newton's identities give the first $g$ coefficients of $P$. Its functional equation supplies the remaining half. All operations in the counter use integer arithmetic and the existing [finite-field implementation](../cubic-tower-factors/field.py). A model with a vanishing leading coefficient, an excluded coefficient denominator or a repeated branch root is skipped at that prime.

The [filter](frobenius_filter.py) returns the counts, the two complete local polynomials and their exact division remainder. A witness proves exclusion of **rational** maps. Geometric map existence can behave differently under a twist and is left undetermined by this test. The five [controls](frobenius-controls.json) include a degree-two map with several poles, a cubic cover, odd- and even-degree quadratic twists, and a coefficient perturbation. The two positive maps pass all three prime checks. The odd-degree twist and perturbation are excluded at 5. The even-degree twist passes at 5 and 7 before a failed divisibility test at 11 excludes every rational map.

```sh
python research/bounded-correspondences/frobenius_filter.py my-curves.json --primes 5 7 11
```

An input gives the two increasing-power coefficient arrays `source` and `target`, as in the reconstruction programs. The optional `--field-budget` sets the maximum number of finite-field elements to enumerate. Exceeding this computational budget skips a reduction and supplies no exclusion. The final reconstruction theorem uses complete coefficient incidences after these inexpensive necessary conditions.
