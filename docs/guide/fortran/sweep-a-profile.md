title: Sweep a curve along a path
category: Construction

Move a tilted line along a circular path. The lower endpoint follows the path and the upper endpoint retains its offset of one unit along both X and Z.

Set `origin` to the point on the profile that should coincide with the start of the path. The [sweep interface](../../interface/sweep.html) translates the profile along the path.

```fortran
program sweep_a_line
    use forcad, only: rk, nurbs_curve, nurbs_surface, sweep
    implicit none
    type(nurbs_curve) :: edge, rail
    type(nurbs_surface) :: wall
    real(rk) :: ends(2, 3)

    call rail%set_circle(center=[0.0_rk, 0.0_rk, 0.0_rk], radius=2.0_rk)
    ends(1,:) = [2.0_rk, 0.0_rk, 0.0_rk]
    ends(2,:) = [3.0_rk, 0.0_rk, 1.0_rk]
    call edge%set(knot=[0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], Xc=ends, degree=1)
    if (.not. rail%err%ok .or. .not. edge%err%ok) error stop "Could not create the curves"
    wall = sweep(edge, rail, origin=ends(1,:))
    if (.not. wall%err%ok) error stop 'Sweep failed'
    print '(3f8.3)', wall%cmp_Xg(Xt=[1.0_rk, 0.5_rk])
end program
```

Halfway around the path, the upper end is at `(-1, 0, 1)`.

**Studio equivalent.** [Sweep a profile along a curve](../how-to/sweep-a-profile.html).
