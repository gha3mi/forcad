title: Revolve a curve into a surface
category: Construction

Revolve a line parallel to the Z axis to construct a quarter-cylinder wall of radius 2 and height 3.

The [`revolve` interface](../../interface/revolve.html) accepts the angle in radians. The first surface parameter follows the line and the second follows the rotation.

```fortran
program quarter_cylinder
    use forcad, only: rk, nurbs_curve, nurbs_surface, revolve
    implicit none
    type(nurbs_curve) :: line
    type(nurbs_surface) :: wall
    real(rk) :: ends(2, 3)

    ends(1,:) = [2.0_rk, 0.0_rk, 0.0_rk]
    ends(2,:) = [2.0_rk, 0.0_rk, 3.0_rk]
    call line%set(knot=[0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xc=ends, degree=1)
    if (.not. line%err%ok) error stop "Could not create the line"
    wall = revolve(line, axis_point=[0.0_rk, 0.0_rk, 0.0_rk], &
        axis_direction=[0.0_rk, 0.0_rk, 1.0_rk], angle=acos(-1.0_rk)/2)
    if (.not. wall%err%ok) error stop 'Revolve failed'
    print '(3f8.3)', wall%cmp_Xg(Xt=[1.0_rk, 1.0_rk])
end program
```

The upper endpoint at the end of the turn is `(0.000, 2.000, 3.000)`.

**Studio equivalent.** [Revolve a profile](../how-to/revolve-a-profile.html).
