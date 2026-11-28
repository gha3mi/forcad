title: Locate a physical point in a volume
category: Evaluation

Use `nearest_point2` on a volume to recover parameters for a physical position. The search includes the interior and all boundaries. For an inverse query, check the physical distance between the input and the returned point as well as the solver diagnostic.

```fortran
program locate_a_point_in_a_volume
    use forcad, only: rk, nurbs_volume
    implicit none
    type(nurbs_volume) :: block
    real(rk) :: query(3,2), uvw(3), point(3), distance
    integer :: i

    call block%set_hexahedron(L=[4._rk, 2._rk, 3._rk], nc=[2, 2, 2])
    if (.not. block%err%ok) error stop "Could not create the block"
    query(:,1) = [1._rk, 1._rk, 2.25_rk]
    query(:,2) = [5._rk, 1._rk, 2.25_rk]

    do i = 1, size(query, 2)
        call block%nearest_point2(point_Xg=query(:,i), tol=sqrt(epsilon(1._rk)), &
            maxit=30, nearest_Xt=uvw, nearest_Xg=point)
        if (.not. block%err%ok) error stop "Volume query failed"
        distance = norm2(point - query(:,i))
        print '(a,i0,a,3f7.2,a,f7.2)', "Query ", i, "  UVW ", uvw, "  Distance (mm) ", distance
    end do
end program
```

| Query (mm) | Returned \((u,v,w)\) | Closest position (mm) | Distance (mm) |
| --- | --- | --- | ---: |
| \((1,1,2.25)\) | \((0.25,0.50,0.75)\) | \((1,1,2.25)\) | 0 |
| \((5,1,2.25)\) | \((1,0.50,0.75)\) | \((4,1,2.25)\) | 1 |

The interior query recovers its position to roundoff. The exterior query reaches the U max face. On a general curved volume the solver can converge to a local minimum. A residual within your physical tolerance verifies a match. A larger residual alone cannot establish that the point lies outside.

**Studio snapping.** On a volume, **Snap to closest point** targets the visible exterior face. For a query restricted to a known face in Fortran, extract that face with `slice`, then [project onto the resulting surface](project-a-point.html).
