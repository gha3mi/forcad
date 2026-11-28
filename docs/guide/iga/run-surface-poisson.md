title: Reproduce a two-dimensional Poisson solution

Run the bundled Fortran solver on the unit square and verify the computed scalar field against an exact solution. The square has zero field values on all four sides and a load that changes sign in each direction.

\[
-\Delta u = \sin(2\pi x)\sin(2\pi y), \qquad
u_{\mathrm{exact}}(x,y)=\frac{\sin(2\pi x)\sin(2\pi y)}{8\pi^2}.
\]

## Run the example

From the root of a ForCAD checkout, run

```sh
pixi run -e gnu fpm run --example poisson_iga_solver_2d --profile release
```

With the default working precision, the program reports 64 analysis elements and 289 scalar degrees of freedom. The quadrature estimate of \(\lVert u_h-u_{\mathrm{exact}}\rVert_{L^2}\) is approximately \(2.67\times10^{-11}\). Compare the order of magnitude if your compiler or precision differs.

## Inspect the field

Open `vtk/poisson_iga_solver_2d_interp.vtk` in a VTK viewer and select the point scalar `u`. The field vanishes along the outer boundary and both midlines. The lower-left and upper-right quadrants are positive. The other two are negative.

<figure class="doc-figure" markdown="1">

[![Computed scalar field on the unit square. Positive peaks lie in the lower-left and upper-right quadrants, with negative peaks in the other two.](../../media/iga-poisson-field.webp)](../../media/iga-poisson-field.webp)

<figcaption>The field evaluated on a 50 × 50 grid from the example’s VTK output. The color scale is in thousandths. The analysis uses an 8 × 8 grid of knot spans.</figcaption>
</figure>

The companion file `vtk/poisson_iga_solver_2d.vtk` stores solution coefficients on the control net. To see how the program turns those coefficients into the pictured field, find the call to `basis` and the multiplication by the coefficient vector in the [example source](https://github.com/gha3mi/forcad/blob/dev/example/poisson_iga_solver_2d.f90).
