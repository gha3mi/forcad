title: Fortran API
ordered_subpage: build-library.md
ordered_subpage: run-curve-example.md
ordered_subpage: evaluate-rational-basis.md
ordered_subpage: evaluate-a-surface.md
ordered_subpage: project-a-point.md
ordered_subpage: locate-a-point-in-a-volume.md
ordered_subpage: measure-a-patch.md
ordered_subpage: fit-a-curve.md
ordered_subpage: fit-a-surface.md
ordered_subpage: extrude-a-surface.md
ordered_subpage: revolve-a-profile.md
ordered_subpage: sweep-a-profile.md
ordered_subpage: loft-sections.md
ordered_subpage: move-volume-control.md
ordered_subpage: change-curve-weight.md
ordered_subpage: insert-a-knot.md
ordered_subpage: remove-a-knot.md
ordered_subpage: elevate-degree.md
ordered_subpage: reduce-degree.md
ordered_subpage: refine-a-patch.md
ordered_subpage: set-domain.md
ordered_subpage: split-a-volume.md
ordered_subpage: join-volume-patches.md
ordered_subpage: control-ids-on-a-knot-span.md
ordered_subpage: integrate-a-surface-field.md
ordered_subpage: curve-length-sensitivity.md
ordered_subpage: shared-control-sensitivity.md
ordered_subpage: shared-face-volume-sensitivity.md
ordered_subpage: sample-volume-mesh.md
ordered_subpage: export-volume-vtk.md
ordered_subpage: export-curve-iges.md

## Start a program

[Evaluate your first NURBS curve](../tutorials/first-curve.html) sets up an fpm project, builds a curve, and checks its midpoint. To run a complete program from a guide below, use that project’s `app/main.f90` and run `fpm run --profile release`. To work from a checkout, [build ForCAD](build-library.html) or [run a bundled example](run-curve-example.html).

## Choose a task

| Task | Worked examples |
| --- | --- |
| Evaluate geometry | [Rational basis values on a quarter circle](evaluate-rational-basis.html), [a surface point and its tangents](evaluate-a-surface.html), [length, area, or volume](measure-a-patch.html) |
| Locate a point | [Project onto a surface or find a sampled vertex](project-a-point.html), [recover parameters inside a volume](locate-a-point-in-a-volume.html) |
| Fit measured data | [Fit a curve to a profile](fit-a-curve.html), [fit a surface to grid measurements](fit-a-surface.html) |
| Construct geometry | [Extrude a surface](extrude-a-surface.html), [revolve a line](revolve-a-profile.html), [sweep along a path](sweep-a-profile.html), [loft between two circles](loft-sections.html) |
| Edit the control net | [Move a control point](move-volume-control.html), [change a weight](change-curve-weight.html) |
| Edit knots and spans | [Insert a knot](insert-a-knot.html), [remove copies of a knot](remove-a-knot.html), [subdivide knot spans](refine-a-patch.html) |
| Change the degree | [Elevate a surface degree](elevate-degree.html), [reduce a curve degree within a tolerance](reduce-degree.html) |
| Change the domain or split a patch | [Rescale a parameter interval](set-domain.html), [split a volume](split-a-volume.html) |
| Assemble patches | [Join two volumes and check their shared controls](join-volume-patches.html) |
| Work with connectivity | [Find the active controls on a knot span](control-ids-on-a-knot-span.html), [get nodes and cells from a sampled volume](sample-volume-mesh.html) |
| Differentiate geometry | [Arc length versus control height or log weight](curve-length-sensitivity.html), [length gradient at a shared control](shared-control-sensitivity.html), [volume response on both sides of a shared face](shared-face-volume-sensitivity.html) |
| Integrate a field | [Use surface quadrature points and weights](integrate-a-surface-field.html) |
| Write a file | [Export a sampled grid to VTK](export-volume-vtk.html), [write a NURBS curve to IGES](export-curve-iges.html) |

## Isogeometric analysis

[Work through the IGA examples](../iga/index.html) to assemble Poisson systems on intervals, surfaces, and volumes. The [two-volume guide](../iga/solve-two-volumes-poisson.html) checks flux balance on an interior face.

## Look up the API

Use the [public API reference](../reference/api.html) for argument types and array shapes, [build environments](../reference/build.html) for compiler and precision options, and the [command comparison](../reference/studio-fortran.html) to find the matching Studio operation.
