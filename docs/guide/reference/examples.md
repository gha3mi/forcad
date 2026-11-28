title: Bundled examples

## Studio projects

Open **Chess** or **Head** under **Samples** on Studio Home. To select a chess piece and inspect its parameter space, follow [Explore a project in Studio](../tutorials/explore-studio.html).

<div class="doc-gallery" markdown="1">

<figure class="doc-figure" markdown="1">

![Chess board and pieces next to the Scene tree in ForCAD Studio](../../media/studio-chess.webp)

<figcaption>Chess illustrates multipatch volumes and saved meshes.</figcaption>
</figure>

<figure class="doc-figure" markdown="1">

![Two views of the Head sample with surface patch edges visible in ForCAD Studio](../../media/studio-head.webp)

<figcaption>Head uses surface patches to model the face, glasses, and cap.</figcaption>
</figure>

</div>

## Fortran examples

Run a target from the checkout root. The target names below link to their source files in the repository’s `dev` branch.

```sh
pixi run -e gnu fpm run --example example_plate_hole_4_3d --profile release
```

| Example target | Topic |
| --- | --- |
| [`demo_curve`](https://github.com/gha3mi/forcad/blob/dev/example/demo_curve.f90) | Create and sample a NURBS curve |
| [`demo_surface`](https://github.com/gha3mi/forcad/blob/dev/example/demo_surface.f90) | Refine a surface and export its parameter lines |
| [`demo_volume`](https://github.com/gha3mi/forcad/blob/dev/example/demo_volume.f90) | Create and export a volume patch |
| [`example_helix_pipe`](https://github.com/gha3mi/forcad/blob/dev/example/example_helix_pipe.f90) | Construct a helical pipe |
| [`example_plate_hole_4_3d`](https://github.com/gha3mi/forcad/blob/dev/example/example_plate_hole_4_3d.f90) | Construct a plate volume with an elliptical cutout |
| [`lsq_fit_bspline_3d`](https://github.com/gha3mi/forcad/blob/dev/example/lsq_fit_bspline_3d.f90) | Fit a spline volume to sampled data |
| [`poisson_iga_solver_2d`](https://github.com/gha3mi/forcad/blob/dev/example/poisson_iga_solver_2d.f90) | [Run the unit-square Poisson problem and inspect its field](../iga/run-surface-poisson.html) |
| [`poisson_iga_solver_3d`](https://github.com/gha3mi/forcad/blob/dev/example/poisson_iga_solver_3d.f90) | [Inspect the field inside the unit cube and verify its error](../iga/run-volume-poisson.html) |
