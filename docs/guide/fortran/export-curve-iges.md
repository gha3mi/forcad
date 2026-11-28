title: Export a NURBS curve to IGES
category: File exchange

Export a circle of radius 1.5 inches as one rational B-spline curve. The IGES entity stores its degree, knots, weights and control points, so no sampling resolution is needed.

Use the [fpm project from the first curve tutorial](../tutorials/first-curve.html#add-the-dependency) and replace `app/main.f90` with this program.

```fortran
program export_iges_circle
    use forcad, only: nurbs_curve, rk
    implicit none
    type(nurbs_curve) :: circle

    call circle%set_circle(center=[0.0_rk, 0.0_rk, 0.0_rk], radius=1.5_rk)
    if (.not. circle%err%ok) error stop "Could not construct the circle"
    call circle%export_iges("circle.igs")
end program export_iges_circle
```

Run `fpm run --profile release` from a writable project directory. The file writer handles I/O errors separately from `circle%err`. Open `circle.igs` in a text editor to verify the export.

- The first parameter record begins with `126,6,2,`. Entity 126 is a rational B-spline curve. The upper control index 6 gives seven control points, and the degree is 2.
- The global section contains `4HINCH`. The public curve writer records coordinates as supplied and uses inches in the file header. Use inch coordinates when you want the exported physical size to match this header.

In a CAD viewer that honors IGES units, the result is a circle in the XY plane with a 3-inch diameter. In ForCAD Studio, [Export a patch as IGES](../how-to/export-iges.html) also writes the active project unit to the header.
