title: Join two volume patches
category: Assembly and analysis

Register two adjacent blocks in a volume assembly and pair their touching faces. This example uses matching V and W directions on both sides and checks positional continuity (\(C^0\)) across the join.

```fortran
program join_volume_patches
    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume
    implicit none
    type(nurbs_volume) :: left, right
    type(nurbs_multipatch_volume) :: assembly
    integer :: a, b
    integer, allocatable :: map(:)

    call left%set_hexahedron(L=[0.5_rk, 1.0_rk, 0.5_rk], nc=[2, 2, 2])
    right = left
    call right%translate_Xc([0.5_rk, 0.0_rk, 0.0_rk])

    call assembly%add_patch(left, id=a)
    call assembly%add_patch(right, id=b)
    call assembly%connect(a, 'u_max', b, 'u_min', continuity=0)
    if (.not. assembly%err%ok) error stop "Could not record the interface"
    if (.not. assembly%is_valid()) error stop "Interface validation failed"

    map = assembly%cmp_dof_map()
    print '(a,i0)', "Patches: ", assembly%get_npatch()
    print '(a,i0)', "Interfaces: ", assembly%get_nconnection()
    print '(a,i0)', "Distinct control IDs: ", maxval(map)
end program
```

The program prints **2 patches**, **1 interface** and **12 distinct control IDs**. Each member still owns eight control points. The four pairs on the shared face receive common IDs in the assembly map.

**Studio equivalent.** [Join two volume patches](../how-to/join-volume-patches.html) shows the corresponding **U max** and **U min** faces in the interface editor. For a join with different face orientations, see [Multipatch interfaces](../explanation/multipatch.html).
