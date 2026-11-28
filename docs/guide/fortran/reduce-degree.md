title: Reduce a curve degree within a tolerance
category: Patch editing

Find out how far a cubic curve can be simplified within a given distance in model units. The arch below has an exact quadratic representation. A straight line cannot follow it within the chosen tolerance.

```fortran
program reduce_curve_degree
    use forcad, only: rk, nurbs_curve
    implicit none
    type(nurbs_curve) :: arch
    real(rk) :: controls(4, 2), bound
    real(rk), parameter :: tolerance = sqrt(epsilon(1._rk))
    integer :: reduced

    controls(:, 1) = [0._rk, 1._rk, 2._rk, 3._rk]
    controls(:, 2) = [0._rk, 2._rk, 2._rk, 0._rk]
    call arch%set(Xc=controls)
    if (.not. arch%err%ok) error stop "Could not create the arch"

    call arch%reduce_degree(t=2, tol=tolerance, reduced=reduced, error_bound=bound)
    if (.not. arch%err%ok) error stop "Invalid degree reduction"
    print '(a,i0)', "Degrees removed: ", reduced
    print '(a,i0)', "Final degree: ", arch%get_degree()
    print '(a,es10.2)', "Error bound: ", bound
    print '(a,2f8.3)', "Midpoint: ", arch%cmp_Xg(0.5_rk)
    if (reduced /= 1 .or. bound > tolerance) error stop "Unexpected reduction"
end program
```

ForCAD accepts one of the two requested steps. The final degree is 2, and the midpoint is `(1.500, 1.500)`. The returned `error_bound` bounds the physical change from the original curve across the whole call. Check `err%ok`, then read `reduced` to distinguish zero, partial, and full acceptance.

**Studio equivalent.** [Reduce the degree of a patch](../how-to/reduce-degree.html) uses **Degrees to remove** and **Error tolerance** for the same two inputs.
