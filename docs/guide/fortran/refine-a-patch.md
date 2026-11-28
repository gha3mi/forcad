title: Refine the knot spans of a surface
category: Patch editing

Give a rectangular surface two unequal U spans, then divide each span in half. The example also increases the U degree by one. Both operations preserve the physical surface.

```fortran
program refine_a_patch
    use forcad, only: rk, nurbs_surface
    implicit none
    type(nurbs_surface) :: patch
    real(rk), parameter :: probe(2) = [0.5_rk, 0.5_rk]
    real(rk), allocatable :: knots(:), midpoints(:)
    real(rk) :: before(3), after(3)
    integer :: first, last, i

    call patch%set_tetragon(L=[2.0_rk, 1.0_rk], nc=[2, 2])
    call patch%insert_knots(dir=1, Xth=[0.25_rk], r=[1])
    if (.not. patch%err%ok) error stop "Could not prepare the surface"
    before = patch%cmp_Xg(probe)

    call patch%elevate_degree(dir=1, t=1)
    if (.not. patch%err%ok) error stop "Could not elevate the U degree"
    knots = patch%get_knot(1)
    first = patch%get_degree(1) + 1
    last = patch%get_nc(1)
    associate (left => knots(first:last), right => knots(first+1:last+1))
        midpoints = pack(0.5_rk * left + 0.5_rk * right, right > left)
    end associate
    call patch%insert_knots(dir=1, Xth=midpoints, r=[(1, i=1, size(midpoints))])
    if (.not. patch%err%ok) error stop "Could not subdivide the U spans"

    after = patch%cmp_Xg(probe)
    if (maxval(abs(after - before)) > 100 * epsilon(1.0_rk)) error stop "Probe moved"
    print '(a,*(1x,f5.3))', "Inserted knots:", midpoints
    print '(a,2(1x,i0))', "Degrees U V:", patch%get_degree()
    print '(a,2(1x,i0))', "Controls U V:", patch%get_nc()
end program
```

The inserted knots are `0.125` and `0.625`, the midpoints of `[0, 0.25]` and `[0.25, 1]`. U now has four positive-length spans. The degrees are `2 1` and the control counts are `7 2`. The midpoint probe still evaluates to `(1, 0.5, 0)`.

<figure class="doc-figure">
  <a href="../../media/refined-u-spans.webp"><img src="../../media/refined-u-spans.webp" width="320" height="347" alt="U–V plot with knot markers along U at 0, 0.125, 0.25, 0.625 and 1. The four resulting spans have unequal widths."></a>
  <figcaption>In Studio’s <strong>Knots</strong> view, the diamonds mark distinct U knots. The original boundary at 0.25 separates the two shorter spans on the left from the two longer spans on the right.</figcaption>
</figure>

The slices select adjacent knots inside the active U interval. The test `right > left` skips repeated knots introduced by elevation. For \(d\) divisions of a span \([a,b]\), insert each interior value \(a+(b-a)j/d\) once, for \(j=1,\ldots,d-1\).

**Studio equivalent.** Create a 2-by-1 rectangle and [insert a U knot](../how-to/insert-knot.html) at 0.25. Open [Refine patch](../how-to/refine-patch.html), choose **U**, set **Divisions per span** to **2** and **Degree increase** to **1**. Studio elevates the degree before inserting the new knots.
