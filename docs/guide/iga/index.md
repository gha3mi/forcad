title: Isogeometric analysis
ordered_subpage: solve-poisson.md
ordered_subpage: share-interface-dofs.md
ordered_subpage: solve-two-patch-poisson.md
ordered_subpage: solve-two-volumes-poisson.md
ordered_subpage: run-surface-poisson.md
ordered_subpage: run-volume-poisson.md

In isogeometric analysis (IGA), the same basis functions \(R_i\) that define the geometry also approximate a solution field. For example, on a spline patch with parameter \(\boldsymbol{\xi}\), write

\[
\mathbf{x}(\boldsymbol{\xi})=\sum_i R_i(\boldsymbol{\xi})\,\mathbf{P}_i,
\qquad
u_h(\boldsymbol{\xi})=\sum_i R_i(\boldsymbol{\xi})\,c_i.
\]

The control points \(\mathbf{P}_i\) specify the shape. The scalar coefficients \(c_i\) are the unknowns of the field problem. ForCAD supplies the active basis indices, physical derivatives, and weighted integration measures for each element. Your Fortran program uses them to assemble the equations and impose boundary conditions.

## Worked examples

| Goal | What you can verify |
| --- | --- |
| [Solve Poisson’s equation on an interval](solve-poisson.html) | Assemble one quadratic element and recover the exact midpoint value \(1/4\). |
| [Share degrees of freedom across two patches](share-interface-dofs.html) | Two bilinear patches share two edge coefficients, leaving six global unknowns. |
| [Solve Poisson across two surface patches](solve-two-patch-poisson.html) | Assemble one system with shared edge IDs and check the field and its normal derivative from both sides. |
| [Solve 3D Poisson across two volumes](solve-two-volumes-poisson.html) | Share face coefficients, inspect the computed field, and balance the outward fluxes. |
| [Reproduce a two-dimensional Poisson solution](run-surface-poisson.html) | Run the bundled surface solver, check its error against a known solution, and inspect the exported field. |
| [Verify a three-dimensional Poisson solution](run-volume-poisson.html) | Inspect contours inside the unit cube and check the volume-integrated error. |

## Element data

Each nonzero knot span contributes an analysis element along its parameter direction. A row of `cmp_elem()` identifies the basis functions active on an element. The sampled points and cells used for visualization have their own numbering. Increasing the display resolution leaves the analysis space unchanged. For small runnable examples of these two building blocks, see [Find active controls on a knot span](../fortran/control-ids-on-a-knot-span.html) and [Integrate a field over a surface](../fortran/integrate-a-surface-field.html).

To explore the same basis interactively, [probe a volume in ForCAD Studio](../tutorials/explore-studio.html#inspect-a-volume-patch). The parameter cube and the basis plot follow the selected patch and direction.
