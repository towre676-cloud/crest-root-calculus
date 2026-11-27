# Crest–Root visual companion

Open an exhibit's file page, select **Download raw file**, then open the downloaded HTML in a browser. Each exhibit runs offline, starts automatically and offers pause, replay and chapter controls. [Mathematical illustrations](../figures/README.md)

| Exhibit | What it shows |
| --- | --- |
| [Program overview](Crest-Root-Visual-Monograph.html) | The inverse problem, exact map reconstruction, genus-two completeness, cubic towers, CM gluings, Hecke walks and the open residual bound. |
| [Opposite secants](Crest-Root-Opposite-Secants-Autoplay.html) | The finite geometry behind equality in every degree prime to three. The illustrated degree-nine example comes from the infinite family. |
| [The discriminant −19 twins](Crest-Root-D19-Twins-Lattice-to-Curve.html) | Follow the lattice to the branch points, descend the equations to a quartic field, then recover the lattice from the equations by certified two-torsion gluing. |

![Certified identification of the discriminant −19 twins](D19-quartic-descent.webp)

The numerical branch charts show how the models were found. The closing chapter shows how the certificate identifies each model independently, from the two elliptic cubics through their half-period labels to an integral isometry. The quotient counts agree at every degree prime to three and differ at degree three. [Equations and arithmetic verifier](../../research/elliptic-quotients/twin-curves/README.md)

## Quotient spectrum comparison

The figure on the [main page](../../README.md#what-the-quotient-counts-reveal) plots these exact counts. For each certified Gram matrix $G$ and degree $d$, count primitive integer vectors $v$ with $v^{\mathsf T}Gv=d^2$, with both signs counted. To enumerate a degree range through $D$, each coordinate satisfies $|v_i|\le D\sqrt{(G^{-1})_{ii}}$. The all-degree equality for the determinant-1,440 pair follows from the neighbour-isometry theorem. [Discriminant −19 matrices](../../research/elliptic-quotients/twin-curves/twin_lattices_-19.json) · [Determinant-1,440 matrices](../../research/elliptic-quotients/full-spectrum-curves/lattices_1440.json) · [All-degree proof](../../research/elliptic-quotients/opposite-secant-spectra/delay/README.md#an-infinite-family-with-identical-full-spectra)

| Degree | −19, curve A | −19, curve B | Determinant 1,440, each curve |
| ---: | ---: | ---: | ---: |
| 2 | 2 | 2 | 0 |
| 3 | 0 | 4 | 2 |
| 4 | 0 | 0 | 2 |
| 5 | 0 | 0 | 2 |
| 6 | 4 | 4 | 2 |
| 7 | 4 | 4 | 2 |
| 8 | 4 | 4 | 0 |
| 9 | 12 | 0 | 6 |
| 10 | 4 | 4 | 2 |
