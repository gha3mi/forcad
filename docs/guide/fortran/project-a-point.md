title: Project a point onto a surface
category: Evaluation

Use `nearest_point2` to project a physical position onto the continuous surface. It returns the parameters and, optionally, the physical position. If you need the ID of an existing sampled vertex, use `nearest_point` on a grid created with `create`.

The example projects \((3,0.5,2)\) onto a rectangle in the XY plane. Its sides are 4 and 2 model units, and both parameter intervals are \([0,1]\).

```fortran
program project_a_point
    use forcad, only: rk, nurbs_surface
    implicit none
    type(nurbs_surface) :: surface
    real(rk), parameter :: query(3) = [3._rk, 0.5_rk, 2._rk]
    real(rk) :: uv(2), point(3)
    integer :: sample_id

    call surface%set_tetragon(L=[4._rk, 2._rk], nc=[2, 2])
    if (.not. surface%err%ok) error stop "Could not create the rectangle"
    call surface%nearest_point2(point_Xg=query, tol=sqrt(epsilon(1._rk)), maxit=30, &
        nearest_Xt=uv, nearest_Xg=point)
    if (.not. surface%err%ok) error stop "Surface projection failed"

    print '(a,2f8.3)', "Parameters ", uv
    print '(a,3f8.3)', "Point      ", point
    print '(a,f8.3)',  "Distance   ", norm2(query - point)

    call surface%create(res1=2, res2=2)
    if (.not. surface%err%ok) error stop "Could not sample the corners"
    call surface%nearest_point(point_Xg=query, nearest_Xg=point, id=sample_id)
    if (sample_id == 0) error stop "No cached sample found"
    print '(a,i0,a,3f8.3)', "Sample ", sample_id, " at ", point
end program
```

The result is \((u,v)=(0.75,0.25)\), which gives the physical point \((3,0.5,0)\) and a distance of 2 model units. The two-by-two grid contains only the corners. Its nearest sample is vertex 2 at \((4,0,0)\). Increasing the grid resolution changes the available vertices.

On a folded surface, several local minima can exist. To start near a known feature, pass an in-domain parameter pair as `initial_Xt`. For a position inside a solid, [query the volume](locate-a-point-in-a-volume.html).

**Related Studio control.** In the viewport **Snap** menu, choose **Snap to closest point** and drag a control point near another visible patch. For a numerical query from a program, use the [surface API](../../type/nurbs_surface.html).
