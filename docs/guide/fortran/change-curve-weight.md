title: Change a curve weight
category: Patch editing

Start with the original three-control quadratic curve in [Evaluate your first NURBS curve](../tutorials/first-curve.html). At \(u=0.5\), its position is \((1,1)\). The third control point is at \((2,0)\).

Insert these lines immediately before the tutorial's `point = curve%cmp_Xg(0.5_rk)` line.

```fortran
call curve%modify_Wc(W=2.0_rk, num=3)
if (.not. curve%err%ok) error stop "Could not change weight 3"
```

Run `fpm run --profile release` again. The program now prints `1.200  0.800`. Increasing the third weight pulls the evaluated point toward that control point. Indices start at one, and the new weight must be finite and positive.

The [curve API](../../type/nurbs_curve.html) also accepts a complete weight vector. If your program exports cached samples, call `create()` to refresh a previously sampled grid after the edit.

**Studio equivalent.** [Change a control-point weight](../how-to/change-weight.html) in **Patch → Points + weights**.
