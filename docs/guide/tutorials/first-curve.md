title: Evaluate your first NURBS curve

Build a quadratic curve from three control points, evaluate it halfway along
its parameter interval, and see what changes when you move the middle control
point. You need [fpm](https://fpm.fortran-lang.org/), a supported Fortran
compiler, and internet access to fetch the dependencies on the first run.
See [Build environments](../reference/build.html) for compiler versions.

## Add the dependency

Create a folder named `first-curve` with an `app` subfolder. Save this manifest
as `first-curve/fpm.toml`.

```toml
name = "first-curve"
version = "0.1.0"

[dependencies]
forcad = { git = "https://github.com/gha3mi/forcad.git", branch = "dev" }
```

## Write the program

Save the following program as `first-curve/app/main.f90`. Each row of
`control` gives the X and Y coordinates of one control point.

```fortran
program first_curve
    use forcad, only: nurbs_curve, rk
    implicit none

    type(nurbs_curve) :: curve
    real(rk) :: control(3, 2)
    real(rk), parameter :: knots(6) = [0._rk, 0._rk, 0._rk, 1._rk, 1._rk, 1._rk]
    real(rk), allocatable :: point(:)

    control(1, :) = [0._rk, 0._rk]
    control(2, :) = [1._rk, 2._rk]
    control(3, :) = [2._rk, 0._rk]

    call curve%set(knots, control, degree=2)
    if (.not. curve%err%ok) error stop "Could not construct the curve"

    point = curve%cmp_Xg(0.5_rk)
    if (size(point) /= 2) error stop "Expected a 2D point"
    print '(2(f6.3,1x))', point
end program first_curve
```

## Run the program

From `first-curve`, run the program.

```sh
fpm run --profile release
```

The program prints `1.000  1.000`. At parameter `0.5`, the curve passes below
the middle control point at `(1, 2)`. Now change that point to
`[1._rk, 4._rk]` and run the program again. The result is `1.000  2.000`.
The endpoints stay in place because the knot vector clamps both ends. The
call to `set` uses unit weights because it omits `Wc`.

To sample a longer curve and write VTK files, run the
[bundled example](../fortran/run-curve-example.html).
See the [API reference](../reference/api.html) for the other patch ranks.
