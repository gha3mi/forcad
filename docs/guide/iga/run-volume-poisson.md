title: Verify a three-dimensional Poisson solution

Solve Poisson’s equation on the unit cube with zero values on all six faces. The load changes sign across the middle of each coordinate direction, which makes this a useful test of a volume implementation. If you have worked through the [surface example](run-surface-poisson.html), the new ingredient is integration over physical volume.

## Run the bundled solver

From the root of a ForCAD checkout, run

```sh
pixi run -e gnu fpm run --example poisson_iga_solver_3d --profile release
```

The example starts with six controls per direction and inserts knots at \(1/3\) and \(2/3\) along each axis. The resulting degree-five patch has \(3^3=27\) analysis elements and \(8^3=512\) scalar degrees of freedom before the boundary conditions are imposed.

## See inside the cube

Open `vtk/poisson_iga_solver_3d_interp.vtk` in a VTK viewer. Select the point scalar `u` and apply a **Contour** filter at values `-0.006` and `0.006`. An isosurface joins positions where the computed field has the same value. The outer faces have prescribed value zero, so the interior contours reveal the sign changes.

<figure class="doc-figure" markdown="1">

[![Eight interior lobes of the computed Poisson field on the unit cube. Four positive and four negative isosurfaces alternate across the three coordinate directions.](../../media/iga-volume-poisson.webp)](../../media/iga-volume-poisson.webp)

<figcaption>The positive contour (ochre) and negative contour (blue) each occupy four octants. Crossing any midplane reverses the sign. The outline marks the boundary of the unit cube.</figcaption>
</figure>

For a point check, inspect the exported nodes nearest \((0.25,0.25,0.25)\) and \((0.25,0.25,0.75)\). Their values are approximately `+0.0084163` and `−0.0084163`. The export samples 50 positions per direction, so the nearest layers are at \(z=12/49\) and \(37/49\). The companion file `vtk/poisson_iga_solver_3d.vtk` contains the spline coefficients on the control net.

## Check the error

The load is \(f(x,y,z)=\sin(2\pi x)\sin(2\pi y)\sin(2\pi z)\). The exact field is \(u_{\mathrm{exact}}=f/(12\pi^2)\), since differentiating twice along each axis contributes \(4\pi^2\) to the eigenvalue of \(-\Delta\).

Compare the terminal’s `L2 error norm` line with this result from the unmodified example at default precision.

```text
L2 error norm                        =  2.4921E-06
```

The solver computes an absolute norm over the cube. Inside its quadrature loop, `T` holds the active basis values, `Xg` is the physical position, and `dV` includes the integration weight and the volume Jacobian.

```fortran
call vol%ansatz(ie, ig, T, dT, dV)
Xg = matmul(T, Xc(elem_e,:))
l2_error = l2_error + (dot_product(T, X(elem_e,1)) &
                      - exact_solution(Xg, L, m))**2 * dV
```

After summing these squared differences over all elements and quadrature points, it takes the square root.

\[
\lVert u_h-u_{\mathrm{exact}}\rVert_{L^2(\Omega)}
\;\approx\;
\left[\sum_{e,q}
\big(u_h(\mathbf{x}_{eq})-u_{\mathrm{exact}}(\mathbf{x}_{eq})\big)^2
\,dV_{eq}\right]^{1/2}.
\]

Small differences in the final digits can occur across builds. When comparing runs, check the degree, element count, and working precision first. Changing only the 50 × 50 × 50 visualization grid changes the exported samples independently of this quadrature check. [Read the Fortran assembly and error calculation](https://github.com/gha3mi/forcad/blob/dev/example/poisson_iga_solver_3d.f90).

For a volume solve with an internal interface, [assemble two connected volumes and check the two outward fluxes](solve-two-volumes-poisson.html).
