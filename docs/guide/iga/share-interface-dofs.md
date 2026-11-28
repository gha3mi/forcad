title: Share scalar degrees of freedom across an interface

Make a scalar field continuous across a matching patch edge. Two bilinear unit squares meet at \(x=1\). With \(C^0\) continuity, the two pairs of coefficients on the edge share global IDs, so eight patch-local coefficients need only six unknowns.

Save the program as `app/main.f90` in the [first-curve fpm project](../tutorials/first-curve.html#add-the-dependency) and run `fpm run --profile release`.

```fortran
program share_interface_dofs
    use forcad, only: rk, nurbs_surface, nurbs_multipatch_surface
    implicit none
    type(nurbs_surface) :: left, right
    type(nurbs_multipatch_surface) :: assembly
    integer :: left_id, right_id, element
    integer, allocatable :: dof_map(:), unshared(:,:), shared(:,:)

    call left%set_tetragon(L=[1._rk, 1._rk], nc=[2, 2])
    if (.not. left%err%ok) error stop "Could not create the left patch"
    right = left
    call right%translate_Xc([1._rk, 0._rk, 0._rk])

    call assembly%add_patch(left, id=left_id)
    call assembly%add_patch(right, id=right_id)
    call assembly%connect(left_id, 'u_max', right_id, 'u_min', continuity=0)
    if (.not. assembly%err%ok) error stop "Could not connect the patches"
    if (.not. assembly%is_valid()) error stop "Invalid interface"

    dof_map = assembly%cmp_dof_map()
    unshared = assembly%cmp_elem()
    shared = assembly%cmp_elem(shared=.true.)
    if (any(shared(1,[2,4]) /= shared(2,[1,3]))) error stop "Interface IDs differ"
    print '(a,i0)', "Scalar unknowns: ", maxval(dof_map)
    do element = 1, size(shared, 1)
        if (any(dof_map(unshared(element,:)) /= shared(element,:))) error stop "Connectivity mismatch"
        print '(a,i0,a,4i3)', "Element ", element, ": ", shared(element,:)
    end do
end program
```

| Element | Global coefficient IDs, U varying fastest | IDs on the shared edge |
| --- | --- | --- |
| 1 · left patch | 1, 2, 3, 4 | 2, 4 at U max |
| 2 · right patch | 2, 5, 4, 6 | 2, 4 at U min |

Use each row of `shared` when assembling that element into the global system. This direct sharing applies here because both patches have clamped knots in the normal (U) direction, matching linear bases along the edge (V) and unit weights. To inspect a sampled surface instead, a patch's `get_elem_Xg_vis()` indexes rows of `get_Xg()` after sampling. Changing the sampling resolution leaves the coefficient IDs above unchanged.

[Solve Poisson across the join](solve-two-patch-poisson.html) to use the shared IDs in a stiffness matrix and check the two traces of the solved field.
