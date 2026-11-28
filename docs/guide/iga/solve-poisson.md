title: Solve a one-dimensional Poisson problem

Solve \(-y''(x)=2\) on the unit interval with prescribed endpoint values \(y(0)=y(1)=0\). You will build one quadratic spline element, assemble its stiffness matrix and load vector, and check the solution at the midpoint. Use the [fpm project from the first-curve tutorial](../tutorials/first-curve.html#add-the-dependency) to run the example.

## Set up the element

Three equally spaced geometry controls at 0, 0.5 and 1 define a straight quadratic Bézier curve. Its physical coordinate \(x\) equals its parameter. Approximate the solution as \(y_h(x)=\sum_{j=1}^3 c_j R_j(x)\), where the \(R_j\) are the same basis functions used for the geometry.

Integrating the differential equation against each test function and integrating by parts gives the stiffness and load entries

\[
K_{ij}=\int_0^1 R'_i(x) R'_j(x)\,dx,\qquad
F_i=\int_0^1 2 R_i(x)\,dx.
\]

The clamped basis interpolates the first and last coefficients at the endpoints. Setting \(c_1=c_3=0\) therefore leaves just one free equation, \(K_{22}c_2=F_2\).

## Assemble and solve

Save this program as `app/main.f90`. In `ansatz_element`, `gradient(i,1,q)` is the physical derivative of basis function `i` at quadrature point `q`. The returned `dL(q)` includes the quadrature weight and physical length factor.

```fortran
program poisson_on_an_interval
    use forcad, only: rk, nurbs_curve
    implicit none
    type(nurbs_curve) :: interval
    real(rk) :: control(3, 1), stiffness(3, 3), load(3), coeff(3)
    real(rk), allocatable :: basis(:,:), gradient(:,:,:), dL(:), midpoint_basis(:)
    integer :: i, q

    control(:,1) = [0.0_rk, 0.5_rk, 1.0_rk]
    call interval%set(Xc=control)
    if (.not. interval%err%ok) error stop "Could not build the interval"
    call interval%ansatz_element(ie=1, Tgc=basis, dTgc_dXg=gradient, dL=dL, ngauss=2)
    if (.not. interval%err%ok) error stop "Could not evaluate the element"

    stiffness = 0.0_rk
    load = 0.0_rk
    do q = 1, size(dL)
        load = load + 2.0_rk * basis(:,q) * dL(q)
        do i = 1, 3
            stiffness(i,:) = stiffness(i,:) + gradient(i,1,q) * gradient(:,1,q) * dL(q)
        end do
    end do

    coeff = 0.0_rk
    coeff(2) = load(2) / stiffness(2,2)
    call interval%basis(Xt=0.5_rk, Tgc=midpoint_basis)
    if (.not. interval%err%ok) error stop "Could not evaluate the midpoint"
    print '(a,f9.6)', "K(2,2) = ", stiffness(2,2)
    print '(a,f9.6)', "F(2)   = ", load(2)
    print '(a,3f9.6)', "Coefficients = ", coeff
    print '(a,f9.6)', "y(0.5) = ", dot_product(coeff, midpoint_basis)
end program poisson_on_an_interval
```

## Check the result

Run `fpm run --profile release`. You should get \(K_{22}=4/3\), \(F_2=2/3\), and coefficients `0.000000  0.500000  0.000000`. The program prints `y(0.5) = 0.250000`. The middle basis function has value 0.5 there, so the middle coefficient contributes \(0.5\times0.5=0.25\).

The exact solution \(y(x)=x(1-x)\) is quadratic and belongs to this spline space. Two Gauss points integrate both the stiffness and load exactly up to floating-point roundoff. For more than one span, [find the active control IDs](../fortran/control-ids-on-a-knot-span.html) to place each local contribution in the global system.
