title: Move a volume control point
category: Patch editing

Raise one corner of a unit cube from 1 to 1.25 along Z. The cube has two controls in each parameter direction. Its last control, shown as **2, 2, 2** in Studio, is row 8 in the public API.

```fortran
program move_volume_control
    use forcad, only: rk, nurbs_volume
    implicit none
    type(nurbs_volume) :: block
    real(rk), allocatable :: corner(:)

    call block%set_hexahedron(L=[1.0_rk, 1.0_rk, 1.0_rk], nc=[2, 2, 2])
    if (.not. block%err%ok) error stop "Could not create the block"
    call block%modify_Xc(X=1.25_rk, num=8, dir=3)
    if (.not. block%err%ok) error stop "Could not edit the coordinate"

    corner = block%cmp_Xg([1.0_rk, 1.0_rk, 1.0_rk])
    print '(3f8.3)', corner
end program move_volume_control
```

The program prints `1.000 1.000 1.250`. The block uses implicit unit weights throughout the edit. In `modify_Xc`, `dir=3` selects the physical coordinate Z. The parameter direction W indexes the control lattice. See [control array ordering](../reference/api.html#patch-data-and-evaluation) for the row number on a larger lattice.

**Studio equivalent.** [Move a control point](../how-to/change-control-point.html) shows the same edit in the coordinate table.
