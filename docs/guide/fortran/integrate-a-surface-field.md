title: Integrate a field over a surface
category: Assembly and analysis

Integrate the scalar field \(f(x,y)=2+x\) over a rectangle spanning \([0,4]\times[0,2]\). Use `ansatz_element` when you need the physical coordinates and integration weights for your own field.

```fortran
program integrate_surface_field
    use forcad, only: rk, nurbs_surface
    implicit none
    type(nurbs_surface) :: patch
    integer, allocatable :: elements(:,:)
    real(rk), allocatable :: basis(:,:), gradients(:,:,:), weights(:), points(:,:)
    real(rk) :: integral
    integer :: element

    call patch%set_tetragon(L=[4._rk, 2._rk], nc=[2, 2])
    call patch%insert_knots(dir=1, Xth=[0.25_rk], r=[1])
    if (.not. patch%err%ok) error stop "Could not prepare the patch"

    elements = patch%cmp_elem()
    integral = 0._rk
    do element = 1, size(elements, 1)
        call patch%ansatz_element(ie=element, Tgc=basis, dTgc_dXg=gradients, &
                                  dA=weights, Xg=points)
        if (.not. patch%err%ok) error stop "Could not evaluate quadrature"
        integral = integral + sum((2._rk + points(1,:)) * weights)
    end do
    print '(a,i0,a,f6.2)', "Elements: ", size(elements, 1), "  Integral: ", integral
end program
```

The program prints two elements and an integral of `32.00`. The rectangle has area 8 and the mean value of the field is 4. Each column of `points` is one physical quadrature point. Its entry in `weights` already includes the change of area from parameter space to physical space.

For area alone, use [`cmp_area`](measure-a-patch.html). For element assembly, [find the control IDs](control-ids-on-a-knot-span.html) associated with the local basis.
