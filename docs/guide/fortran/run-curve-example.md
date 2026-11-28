title: Run the bundled curve example
category: Setup

From a [ForCAD checkout with Pixi](build-library.html), run the curve example.

```sh
pixi run -e gnu fpm run --example demo_curve --profile release
```

The example writes `vtk/demo_curve_Xg.vtk` for the sampled curve and `vtk/demo_curve_Xc.vtk` for its control polygon. Open the sampled curve in a VTK viewer such as [ParaView](https://www.paraview.org/).

The [example source](https://github.com/gha3mi/forcad/blob/dev/example/demo_curve.f90) shows how it creates the control polygon, samples the curve, and writes each file.
