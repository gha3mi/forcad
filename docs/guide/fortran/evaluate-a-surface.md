title: Evaluate a surface point and its tangents
category: Evaluation

Use `cmp_Xg` to evaluate a surface at a pair of parameters. Omit `order` for the position. Pass `order=[1, 0]` for the U tangent or `order=[0, 1]` for the V tangent. Each returns one vector in physical coordinates.

This bilinear patch has three corners in the XY plane and one raised corner. Its position is \(\mathbf{S}(u,v)=(2u,3v,uv)\) for \(u,v\in[0,1]\). In the control array, U varies first.

```fortran
program evaluate_a_surface
    use forcad, only: rk, nurbs_surface
    implicit none
    type(nurbs_surface) :: patch
    real(rk) :: control(4, 3)
    real(rk), parameter :: uv(2) = [0.25_rk, 0.5_rk]

    control(1,:) = [0.0_rk, 0.0_rk, 0.0_rk]
    control(2,:) = [2.0_rk, 0.0_rk, 0.0_rk]
    control(3,:) = [0.0_rk, 3.0_rk, 0.0_rk]
    control(4,:) = [2.0_rk, 3.0_rk, 1.0_rk]
    call patch%set(nc=[2, 2], Xc=control)
    if (.not. patch%err%ok) error stop "Could not create the surface"

    print '(a,3f8.3)', "Point ", patch%cmp_Xg(uv)
    print '(a,3f8.3)', "U     ", patch%cmp_Xg(uv, order=[1, 0])
    print '(a,3f8.3)', "V     ", patch%cmp_Xg(uv, order=[0, 1])
end program
```

| At \((u,v)=(0.25,0.5)\) | X | Y | Z |
| --- | ---: | ---: | ---: |
| Position | 0.500 | 1.500 | 0.125 |
| U tangent, \(\partial\mathbf{S}/\partial u\) | 2.000 | 0.000 | 0.500 |
| V tangent, \(\partial\mathbf{S}/\partial v\) | 0.000 | 3.000 | 0.250 |

The tangent components give the change in physical position per unit of the corresponding parameter. See the [surface type reference](../../type/nurbs_surface.html) for mixed derivatives and one-sided limits at knots.

**Studio equivalent.** Select a surface and open **Inspector → Selection and probe**. Enter U and V to inspect its position and first derivatives. The Inspector uses [normalized probe coordinates](../explanation/spline-geometry.html).
