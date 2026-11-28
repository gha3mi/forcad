title: Measure a patch in Fortran
category: Evaluation

Declare the result as `real(rk)` and call the binding for the rank of your initialized patch. All three integrate in the physical coordinates of the patch.

| Patch | Call | Runnable example |
| --- | --- | --- |
| Curve | `call curve%cmp_length(value)` | [Line segment of length 2](https://github.com/gha3mi/forcad/blob/dev/example/cmp_length.f90) |
| Surface | `call surface%cmp_area(value)` | [Square of area 4](https://github.com/gha3mi/forcad/blob/dev/example/cmp_area.f90) |
| Volume | `call volume%cmp_volume(value)` | [Cube of volume 8](https://github.com/gha3mi/forcad/blob/dev/example/cmp_volume.f90) |

Check the patch's `err%ok` before using the result. If your coordinates are in millimeters, length is in mm, area in mm², and volume in mm³.

Run the volume example from the ForCAD checkout.

```sh
pixi run -e gnu fpm run --example cmp_volume --profile release
```

It prints approximately `8`. The [public API reference](../reference/api.html) links to the optional integration arguments.

**Studio equivalent.** [Measure a selected patch](../how-to/measure-geometry.html) using **Model → Measure**.
