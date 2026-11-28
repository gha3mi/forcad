title: Differentiate arc length with respect to a control point and weight
category: Assembly and analysis

How much does a curve's length change if you raise one control point or increase its weight? Start with the [rational quadratic quarter circle](evaluate-rational-basis.html). Keep its endpoints fixed and differentiate its length \(L\) with respect to the middle control's height \(P_{2y}\) and log weight \(\eta_2=\log w_2\).

For one element, `ansatz_element` returns the sensitivity of each weighted length contribution: `dmeasure_dXc(i,j,q)` is the partial derivative with respect to coordinate `j` of local control `i`, and `dmeasure_dlogWc(i,q)` is the partial derivative with respect to its log weight. Sum over quadrature points `q` to obtain the derivative of the integrated length.

```fortran
program curve_length_sensitivity
    use forcad, only: rk, nurbs_curve
    implicit none
    type(nurbs_curve) :: arc
    integer, parameter :: nq = 12
    real(rk), parameter :: middle_weight = sqrt(0.5_rk)
    real(rk), parameter :: h = epsilon(1.0_rk)**(1.0_rk/3.0_rk)
    real(rk) :: control(3, 2), numerical(2)
    real(rk), allocatable :: basis(:,:), gradient(:,:,:), dL(:), dpoint(:,:,:), dlogw(:,:)

    control(:,1) = [1.0_rk, 1.0_rk, 0.0_rk]
    control(:,2) = [0.0_rk, 1.0_rk, 1.0_rk]
    call arc%set(Xc=control, Wc=[1.0_rk, middle_weight, 1.0_rk])
    if (.not. arc%err%ok) error stop "Could not build the arc"
    call arc%ansatz_element(ie=1, Tgc=basis, dTgc_dXg=gradient, dL=dL, ngauss=nq, &
                            dmeasure_dXc=dpoint, dmeasure_dlogWc=dlogw)
    if (.not. arc%err%ok) error stop "Could not evaluate the length sensitivities"

    numerical(1) = (length_at(h, 0.0_rk) - length_at(-h, 0.0_rk)) / (2*h)
    numerical(2) = (length_at(0.0_rk, h) - length_at(0.0_rk, -h)) / (2*h)
    print '(a,f10.6)',  "Length: ", sum(dL)
    print '(a,2f10.6)', "Height: ", sum(dpoint(2,2,:)), numerical(1)
    print '(a,2f10.6)', "Log weight: ", sum(dlogw(2,:)), numerical(2)

contains
    function length_at(dy, deta) result(length)
        real(rk), intent(in) :: dy, deta
        real(rk) :: length
        type(nurbs_curve) :: trial

        trial = arc
        call trial%modify_Xc(X=control(2,2)+dy, num=2, dir=2)
        call trial%modify_Wc(W=middle_weight*exp(deta), num=2)
        call trial%cmp_length(length, ngauss=nq)
        if (.not. trial%err%ok) error stop "Could not measure the perturbed arc"
    end function length_at
end program
```

The unperturbed arc length is `1.570796`, or \(\pi/2\). The two sensitivity checks give the following values with the default real kind.

| Quantity | Integrated sensitivity | Central difference |
| --- | ---: | ---: |
| Middle control height | 0.285398 | 0.285398 |
| Middle log weight | 0.141593 | 0.141593 |

The finite difference measures lengths on independent curve copies, using the same quadrature. Perturbing \(\eta_2\) by \(h\) multiplies the weight by \(e^h\). To get the derivative with respect to the weight itself, divide the log-weight derivative by \(w_2\).

For a control shared by two patches, [combine their length derivatives](shared-control-sensitivity.html).
