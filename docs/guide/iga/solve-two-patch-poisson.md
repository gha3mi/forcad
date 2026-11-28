title: Solve Poisson across two surface patches

Assemble one scalar Poisson problem on a rectangle made from two surface patches. You will use ForCAD’s shared control IDs to add both patches’ contributions to the same linear system, then evaluate the solution independently from each side of the join.

## The problem

The domain is \(\Omega=[0,3]\times[0,1]\). The left patch ends at \(x=1\), where it meets the wider right patch. Prescribe \(u=0\) on the four outer edges. The load and exact solution are

\[
-\Delta u=f,\qquad
f(x,y)=2\big[x(3-x)+y(1-y)\big],\qquad
u_{\mathrm{exact}}(x,y)=x(3-x)y(1-y).
\]

<figure class="doc-figure" markdown="1">

[![Computed scalar field on two rectangles. A vertical dashed line at x = 1 marks the interface. The field is continuous there and reaches its peak in the wider right patch.](../../media/iga-two-patch-poisson.webp)](../../media/iga-two-patch-poisson.webp)

<figcaption>The dot marks the interface midpoint. The solution is still increasing in the positive x direction there and reaches its maximum halfway across the full rectangle.</figcaption>
</figure>

## Assemble with shared IDs

1. Give each patch a \(3\times3\) control net and one biquadratic element. Connect the left patch’s `u_max` edge to the right patch’s `u_min` edge with `continuity=0`. Both clamped edges have the same quadratic trace basis, unit weights, and the same orientation.
2. Request `cmp_elem(shared=.true.)`. Each row supplies the global IDs of one element’s nine local basis functions. The three matching edge controls reuse the same IDs, giving \(9+9-3=15\) scalar unknowns. At each quadrature point, `ansatz_element` supplies the basis values, physical gradients, and weighted area. Add the local stiffness and load into those global rows and columns.
3. Fix the controls on the outer boundary. The two endpoints of the join lie on that boundary. The middle interface coefficient stays free and receives contributions from both patches. Solve for the three free coefficients together.

Each patch has exactly one element in this example, so `shared(p,:)` corresponds to local element 1 of patch `p`. For a refined assembly, `cmp_elem_patch()` and `cmp_elem_local()` identify the owning patch and local element for each row.

## Run the program

Use the [fpm setup from the first-curve tutorial](../tutorials/first-curve.html#add-the-dependency). Save the following program as `app/main.f90` and run `fpm run --profile release`.

```fortran
program two_patch_poisson
    use forcad, only: rk, nurbs_surface, nurbs_multipatch_surface, solve
    implicit none
    type(nurbs_surface) :: patch(2)
    type(nurbs_multipatch_surface) :: domain
    integer, allocatable :: shared(:,:), free(:), ids(:)
    real(rk), allocatable :: stiffness(:,:), load(:,:), coefficients(:), interior(:,:)
    real(rk), allocatable :: basis(:,:), grad(:,:,:), dA(:), points(:,:), controls(:,:)
    real(rk), allocatable :: uv_left(:,:), uv_right(:,:), edge_weights(:)
    logical, allocatable :: fixed(:)
    real(rk), parameter :: tol = 200 * epsilon(1._rk)
    real(rk) :: x, y, source, left(2), right(2), exact(2)
    integer :: left_id, right_id, n, p, i, q

    call patch(1)%set_tetragon(L=[1._rk, 1._rk], nc=[3,3])
    call patch(2)%set_tetragon(L=[2._rk, 1._rk], nc=[3,3])
    call patch(2)%translate_Xc([1._rk, 0._rk, 0._rk])
    call domain%add_patch(patch(1), id=left_id)
    call domain%add_patch(patch(2), id=right_id)
    call domain%connect(left_id, 'u_max', right_id, 'u_min', continuity=0)
    if (.not. domain%err%ok .or. .not. domain%is_valid()) error stop "Invalid join"

    shared = domain%cmp_elem(shared=.true.)
    n = maxval(shared)
    allocate(stiffness(n,n), load(n,1), coefficients(n), source=0._rk)
    allocate(fixed(n), source=.false.)
    do p = 1, 2
        ids = shared(p,:)
        call patch(p)%ansatz_element(ie=1, Tgc=basis, dTgc_dXg=grad, dA=dA, Xg=points)
        if (.not. patch(p)%err%ok) error stop "Element evaluation failed"
        do q = 1, size(dA)
            x = points(1,q)
            y = points(2,q)
            source = 2 * (x*(3-x) + y*(1-y))
            load(ids,1) = load(ids,1) + basis(:,q)*source*dA(q)
            stiffness(ids,ids) = stiffness(ids,ids) + matmul(grad(:,:,q), transpose(grad(:,:,q)))*dA(q)
        end do

        controls = patch(p)%get_Xc()
        do i = 1, size(ids)
            x = controls(i,1)
            y = controls(i,2)
            if (min(x, 3-x, y, 1-y) < tol) fixed(ids(i)) = .true.
        end do
    end do

    free = pack([(i, i=1,n)], .not. fixed)
    if (n /= 15 .or. size(free) /= 3) error stop "Unexpected global numbering"
    interior = solve(stiffness(free,free), load(free,:))
    if (size(interior,1) /= size(free)) error stop "Solve failed"
    coefficients(free) = interior(:,1)
    print '(a,i0,a,i0)', "Global DOFs: ", n, "   Free DOFs: ", size(free)

    call domain%interface_quadrature(1, uv_left, uv_right, edge_weights, ngauss=[3])
    if (size(edge_weights) /= 3) error stop "Missing interface points"
    print '(a)', "     y     u(left)    u(right)   du/dx(left)   du/dx(right)"
    do q = 1, size(edge_weights)
        left  = trace(patch(1), uv_left(q,:),  coefficients(shared(1,:)))
        right = trace(patch(2), uv_right(q,:), coefficients(shared(2,:)))
        y = uv_left(q,2)
        exact = [2*y*(1-y), y*(1-y)]
        print '(3f10.6,2f14.6)', y, left(1), right(1), left(2), right(2)
        if (.not. all(abs(left-right) < tol)) error stop "Interface jump"
        if (.not. all(abs(left-exact) < tol .and. abs(right-exact) < tol)) error stop "Wrong trace"
    end do

contains

    function trace(surface, uv, c) result(value_and_slope)
        type(nurbs_surface), intent(inout) :: surface
        real(rk), intent(in) :: uv(2), c(:)
        real(rk) :: value_and_slope(2)
        real(rk), allocatable :: values(:), derivatives(:,:), tangent(:)
        call surface%derivative(uv, dTgc=derivatives, Tgc=values)
        tangent = surface%cmp_Xg(uv, order=[1,0])
        value_and_slope = [dot_product(c, values), dot_product(c, derivatives(:,1))/tangent(1)]
    end function

end program two_patch_poisson
```

## Check the interface

`interface_quadrature` returns corresponding parameter pairs for the two edges. Here the tangential parameter equals the physical coordinate \(y\). Evaluate each patch’s basis at its own pair and use its own connectivity row to reconstruct the field. The `trace` function also converts the derivative in the local U direction to \(\partial u_h/\partial x\). Dividing by the physical U-tangent’s x component works for these two axis-aligned affine maps. That component is 1 on the left and 2 on the right. Both derivative columns use the same positive x direction. With each patch’s outward normal, the two derivatives have opposite signs and sum to zero.

| Physical \(y\) on \(x=1\) | \(u_h\), left | \(u_h\), right | \(\partial_x u_h\), left | \(\partial_x u_h\), right |
| ---: | ---: | ---: | ---: | ---: |
| 0.112702 | 0.200000 | 0.200000 | 0.100000 | 0.100000 |
| 0.500000 | 0.500000 | 0.500000 | 0.250000 | 0.250000 |
| 0.887298 | 0.200000 | 0.200000 | 0.100000 | 0.100000 |

Sharing the edge coefficients enforces equality of the field values. The derivative check is independent. For this example the exact field is biquadratic on each patch, and the default three-point Gauss rule in each direction integrates the stiffness and load exactly up to roundoff. Both values must reproduce \(2y(1-y)\) and both x-derivatives must reproduce \(y(1-y)\). The program checks these identities at all three points. A general \(C^0\) discretization can have a jump in its normal derivative.

For the numbering alone, see [Share scalar degrees of freedom across an interface](share-interface-dofs.html). To assemble across a face and integrate its outward flux, [solve the two-volume problem](solve-two-volumes-poisson.html).
