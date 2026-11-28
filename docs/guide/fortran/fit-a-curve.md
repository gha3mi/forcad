title: Fit a B-spline curve to measured points
category: Construction

Fit a quadratic curve to five measured positions with known parameter values. Each row of `Xdata` holds a measured point, and the corresponding entry in `Xt` identifies its location along the curve. The least-squares fit adjusts the controls to minimize the sum of squared coordinate residuals at those parameters.

```fortran
program fit_a_curve
    use forcad, only: rk, nurbs_curve
    implicit none
    type(nurbs_curve) :: curve
    real(rk), parameter :: u(5) = [0._rk, 0.25_rk, 0.5_rk, 0.75_rk, 1._rk]
    real(rk) :: measurements(5, 2), rms
    real(rk), allocatable :: fitted(:,:)

    measurements(:,1) = [0._rk, 1._rk, 2._rk, 3._rk, 4._rk]
    measurements(:,2) = [0.125_rk, 0.5_rk, 1._rk, 1._rk, -0.125_rk]
    call curve%set(Xth_dir=[0._rk, 1._rk], degree=2, continuity=[-1, -1])
    if (.not. curve%err%ok) error stop "Could not define the spline space"
    call curve%lsq_fit_bspline(Xt=u, Xdata=measurements, ndata=size(u))
    if (.not. curve%err%ok) error stop "Curve fitting failed"

    fitted = curve%cmp_Xg(u)
    rms = sqrt(sum((fitted - measurements)**2) / size(u))
    print '(a,i0)', "Control points: ", curve%get_nc()
    print '(a,2f8.3)', "Fitted midpoint: ", fitted(3,:)
    print '(a,f9.6)', "RMS distance: ", rms
end program
```

The single quadratic span has three controls. At the middle parameter the fit gives \((2,1)\), and the RMS distance is approximately 0.176777 model units. Each residual compares positions at the supplied parameter. To estimate closest-point distances instead, project each measurement onto the curve.

For your own data, supply at least as many measurements as control points and distribute their parameters across the active spans so that the controls are independently determined. To fit both control positions and positive weights, use `lsq_fit_nurbs`. With either fitting method, evaluate the updated geometry with `cmp_Xg`, or call `create` again to refresh previously cached samples. For measurements arranged in two directions, [fit a surface](fit-a-surface.html).
