title: Fit a surface to grid measurements
category: Construction

Fit a bilinear surface to nine example height readings on a panel 100 mm wide and 60 mm high. The measurement grid supplies the parameter pairs \((u,v)=(x/100,y/60)\). Each row of `uv` corresponds to the same row of `measurements`, with U varying fastest.

```fortran
program fit_a_surface
    use forcad, only: rk, nurbs_surface
    implicit none
    type(nurbs_surface) :: panel
    real(rk), parameter :: x(3) = [0._rk, 50._rk, 100._rk], y(3) = [0._rk, 30._rk, 60._rk]
    real(rk), parameter :: height(9) = [1.10_rk, 1.40_rk, 1.70_rk, &
                                        0.90_rk, 1.58_rk, 1.90_rk, &
                                        0.70_rk, 1.40_rk, 2.10_rk]
    real(rk) :: uv(9,2), measurements(9,3)
    real(rk), allocatable :: fitted(:,:)
    integer :: i, j, row

    do j = 1, size(y)
        do i = 1, size(x)
            row = i + (j - 1) * size(x)
            uv(row,:) = [x(i)/100._rk, y(j)/60._rk]
            measurements(row,:) = [x(i), y(j), height(row)]
        end do
    end do

    call panel%set_tetragon(L=[100._rk, 60._rk], nc=[2, 2])
    if (.not. panel%err%ok) error stop "Could not define the surface"
    call panel%lsq_fit_bspline(Xt=uv, Xdata=measurements, ndata=[size(x), size(y)])
    if (.not. panel%err%ok) error stop "Surface fitting failed"

    fitted = panel%cmp_Xg(uv)
    print '(a,3f8.2)', "Fitted center (mm): ", fitted(5,:)
    print '(a,f6.2)', "Largest height residual (mm): ", maxval(abs(fitted(:,3) - height))
end program
```

The fitted center is \((50,30,1.42)\) mm. Its measured height is 1.58 mm, giving the largest absolute height residual, **0.16 mm**. The least-squares fit uses all nine positions to determine the four controls.

<figure class="doc-figure" markdown="1">

[![Height along the middle measurement row. At x equals 50 mm the measurement is 0.16 mm above the fitted line.](../../media/panel-fit-residual.webp)](../../media/panel-fit-residual.webp)

<figcaption>A section of the fitted surface at y = 30 mm, shown with that row’s three readings. All nine measurements contribute to the fit.</figcaption>
</figure>

For a larger grid, `ndata` holds the number of measurement stations in U and V. Each direction needs at least as many independent stations as controls. If you have positions along a single profile, [fit a curve](fit-a-curve.html).
