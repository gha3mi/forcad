title: Documentation guide
ordered_subpage: fortran
ordered_subpage: how-to
ordered_subpage: iga
ordered_subpage: assistant
ordered_subpage: tutorials
ordered_subpage: reference
ordered_subpage: explanation

## Get started

- **Fortran API.** [Create an fpm project and evaluate a NURBS curve](tutorials/first-curve.html). Continue with the [Fortran task guides](fortran/index.html) to construct, edit, and evaluate patches in your own program.
- **ForCAD Studio.** [Launch the application](how-to/install.html) and [explore the Chess project](tutorials/explore-studio.html). The [Studio task guides](how-to/index.html) cover patch editing, meshing, and file exchange.

## Work through an example

- **Measured data.** [Project a point onto a surface](fortran/project-a-point.html) to recover its surface parameters, or [fit a curve to measured points](fortran/fit-a-curve.html).
- **Shape sensitivity.** [Estimate how a control height or log weight changes curve length](fortran/curve-length-sensitivity.html) and check each derivative against a finite difference.
- **Multipatch analysis.** [Solve a Poisson problem on two volumes](iga/solve-two-volumes-poisson.html) with shared degrees of freedom on the joining face. The [IGA examples](iga/index.html) also cover interval and surface problems.
- **AI-assisted modeling.** [Use an agent in ForCAD Studio](assistant/index.html) to inspect the open model and review its proposed edits.

## Look up a detail

Find [public types and procedures](reference/api.html), [Studio controls and shortcuts](reference/studio.html), or [the Fortran equivalent of a Studio command](reference/studio-fortran.html). For the underlying geometry, read about [parameter space](explanation/spline-geometry.html) or [multipatch interfaces](explanation/multipatch.html).
