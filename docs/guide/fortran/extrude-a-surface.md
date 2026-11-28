title: Extrude a surface into a volume
category: Construction

Give `extrude` a surface and a displacement vector in model coordinates. The new volume retains the profile’s U and V parameterization. W runs from the original profile at zero to its translated copy at one.

Use the fpm project from [Evaluate your first NURBS curve](../tutorials/first-curve.html#add-the-dependency). Replace `app/main.f90` with this program.

```fortran
program extrude_a_surface
    use forcad, only: rk, nurbs_surface, nurbs_volume, extrude
    implicit none
    type(nurbs_surface) :: base
    type(nurbs_volume) :: block
    real(rk) :: volume

    call base%set_tetragon(L=[1.0_rk, 1.0_rk], nc=[2, 2])
    if (.not. base%err%ok) error stop "Could not create the base"
    block = extrude(base, [0.0_rk, 0.0_rk, 0.5_rk])
    if (.not. block%err%ok) error stop "Could not extrude the surface"
    call block%cmp_volume(volume)
    print '(f5.3)', volume
end program extrude_a_surface
```

Run `fpm run --profile release`. The unit square and half-unit displacement give a volume of `0.500` cubic model units. Passing a curve to the same [generic interface](../../interface/extrude.html) produces a surface.

**Studio equivalent.** [Extrude a profile](../how-to/extrude-a-profile.html) shows the same block and displacement in the Create ribbon.
