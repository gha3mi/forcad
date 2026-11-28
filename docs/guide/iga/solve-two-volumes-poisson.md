title: Solve 3D Poisson across two volumes

Join two triquadratic volume patches, assemble one Poisson system, and check that the field is continuous and the outward fluxes balance on their common face. The example writes both computed fields to VTK. For the corresponding edge assembly in two dimensions, see [Solve Poisson across two surface patches](solve-two-patch-poisson.html).

## Define the field problem

The box \([0,3]\times[0,1]\times[0,1]\) is split at \(x=1\). Set \(u=0\) on the exterior and take unit diffusivity. To verify the solution even inside the box, choose \(a=x(3-x)\), \(b=y(1-y)\), and \(c=z(1-z)\). The exact solution and its load are

\[
u_{\mathrm{exact}}=abc,\qquad -\Delta u=f=2(bc+ac+ab).
\]

<figure class="doc-figure" markdown="1">

[![A teal contour of the solved field extends across both volumes. The ochre rectangle marks their common face at x equals one.](../../media/iga-two-volumes-poisson.webp)](../../media/iga-two-volumes-poisson.webp)

<figcaption>The isosurface <i>u</i><sub>h</sub> = 0.08 from the two exported fields crosses the shared face. The field peaks in the wider volume at (1.5, 0.5, 0.5).</figcaption>
</figure>

## Assemble one system

1. Build a \(3\times3\times3\) control net on each box. Join the left volume’s `u_max` to the right volume’s `u_min` with `continuity=0`. On this face the clamped endpoints, matching V/W bases, equal weights, and aligned orientation allow direct coefficient sharing.
2. Use `cmp_elem(shared=.true.)` to obtain global indices. There is one Bézier element per patch, so each row contains 27 indices. Both rows include the same nine face IDs. This leaves \(27+27-9=45\) global scalar degrees of freedom.
3. For each element, `ansatz_element` returns basis values, physical gradients, and quadrature-weighted volume `dV`. Add its stiffness and load to the global matrix and vector at that row’s indices. Fix only controls on the exterior. The central face control and one interior control in each volume remain free.

With multiple spans, use `cmp_elem_patch()` and `cmp_elem_local()` to find the patch and local element for each global row.

## Check the field and flux

For unit diffusivity the flux is \(\mathbf j=-\nabla u_h\). On the shared face the left patch’s outward normal points along \(+x\), and the right patch’s points along \(-x\). The two field values should agree. The two **outward** fluxes should sum to zero.

| Quantity | Left volume | Right volume |
| --- | ---: | ---: |
| Field at the face center | 0.125 | 0.125 |
| Outward flux density at the face center | −0.0625 | +0.0625 |
| Outward flux integrated over the face | \(-1/36\) | \(+1/36\) |

The program checks all nine paired face quadrature points with `interface_quadrature`. It verifies their physical coordinates, reconstructs each patch’s field and flux from its own basis and coefficients, and compares each side against \(u(1,y,z)=2bc\) and \(j_n=\mp bc\). The returned physical area weights integrate each flux separately. Since \(\int_0^1y(1-y)\,dy=1/6\), the flux magnitudes integrate to \(1/36\).

<figure class="doc-figure" markdown="1">

[![Computed field on the common face in y and z. The field vanishes at the four edges and reaches 0.125 in the center.](../../media/iga-common-face-field.webp)](../../media/iga-common-face-field.webp)

<figcaption>The 25 × 25 face samples in the two VTK files coincide. The ring marks the center value in the table.</figcaption>
</figure>

The derivative returned by `volume%derivative` is with respect to the patch parameters. For these axis-aligned boxes, dividing the U derivative by \(\partial x/\partial U\) converts it to the physical x derivative. That scale is 1 on the left and 2 on the right. This exact solution belongs to both spline spaces, so the pointwise flux balance holds to roundoff. A general \(C^0\) discretization enforces continuity of the field and can have a pointwise flux jump.

## Run and reproduce the figures

Use the [fpm setup from the first-curve tutorial](../tutorials/first-curve.html#add-the-dependency). Save the complete program as `app/main.f90` and run `fpm run --profile release` from the project directory. It prints the quantities in the table and writes `poisson_left.vtk` and `poisson_right.vtk` there. Open both files together in a VTK viewer, select the point scalar `u`, and contour each at `0.08`. To inspect the common face, view the left volume’s \(x=1\) face from the positive x direction.

```fortran
program two_volume_poisson
    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume, solve
    implicit none
    type(nurbs_volume) :: patch(2)
    type(nurbs_multipatch_volume) :: domain
    integer, allocatable :: shared(:,:), ids(:), free(:)
    real(rk), allocatable :: stiffness(:,:), load(:,:), coeff(:), interior(:,:)
    real(rk), allocatable :: basis(:,:), grad(:,:,:), dV(:), points(:,:), controls(:,:)
    real(rk), allocatable :: uvw_l(:,:), uvw_r(:,:), dS(:), samples(:,:)
    logical, allocatable :: fixed(:)
    real(rk), parameter :: tol = 300*epsilon(1._rk)
    character(len=*), parameter :: files(2) = [character(len=17) :: 'poisson_left.vtk', 'poisson_right.vtk']
    integer, parameter :: nx(2) = [25,49]
    real(rk) :: x, y, z, a, b, c, source, left(2), right(2), expected, flux_sum(2), xyz(3)
    integer :: p, i, q, n, left_id, right_id

    call patch(1)%set_hexahedron(L=[1._rk,1._rk,1._rk], nc=[3,3,3])
    call patch(2)%set_hexahedron(L=[2._rk,1._rk,1._rk], nc=[3,3,3])
    call patch(2)%translate_Xc([1._rk,0._rk,0._rk])
    call domain%add_patch(patch(1), id=left_id)
    call domain%add_patch(patch(2), id=right_id)
    call domain%connect(left_id, 'u_max', right_id, 'u_min', continuity=0)
    if (.not. domain%err%ok .or. .not. domain%is_valid()) error stop "Invalid face connection"

    shared = domain%cmp_elem(shared=.true.)
    n = maxval(shared)
    allocate(stiffness(n,n), load(n,1), coeff(n), source=0._rk)
    allocate(fixed(n), source=.false.)
    do p = 1, 2
        ids = shared(p,:)
        call patch(p)%ansatz_element(ie=1, Tgc=basis, dTgc_dXg=grad, dV=dV, Xg=points)
        if (.not. patch(p)%err%ok) error stop "Element evaluation failed"
        do q = 1, size(dV)
            x = points(1,q)
            y = points(2,q)
            z = points(3,q)
            a = x*(3-x)
            b = y*(1-y)
            c = z*(1-z)
            source = 2*(b*c + a*c + a*b)
            load(ids,1) = load(ids,1) + basis(:,q)*source*dV(q)
            stiffness(ids,ids) = stiffness(ids,ids) + matmul(grad(:,:,q), transpose(grad(:,:,q)))*dV(q)
        end do
        controls = patch(p)%get_Xc()
        do i = 1, size(ids)
            x = controls(i,1)
            y = controls(i,2)
            z = controls(i,3)
            if (min(x,3-x,y,1-y,z,1-z) < tol) fixed(ids(i)) = .true.
        end do
    end do

    free = pack([(i,i=1,n)], .not. fixed)
    if (n /= 45 .or. size(free) /= 3) error stop "Unexpected global numbering"
    interior = solve(stiffness(free,free), load(free,:))
    if (size(interior,1) /= size(free)) error stop "Solve failed"
    coeff(free) = interior(:,1)
    print '(a,i0,a,i0)', "Global DOFs: ", n, "   Free DOFs: ", size(free)

    call domain%interface_quadrature(1, uvw_l, uvw_r, dS, ngauss=[3,3], physical=.true.)
    if (size(dS) /= 9) error stop "Expected nine face points"
    flux_sum = 0
    do q = 1, size(dS)
        left  = face_trace(patch(1), uvw_l(q,:), coeff(shared(1,:)),  1._rk)
        right = face_trace(patch(2), uvw_r(q,:), coeff(shared(2,:)), -1._rk)
        xyz = patch(1)%cmp_Xg(uvw_l(q,:))
        if (.not. all(abs(xyz - patch(2)%cmp_Xg(uvw_r(q,:))) < tol)) error stop "Face points differ"
        y = xyz(2)
        z = xyz(3)
        expected = y*(1-y)*z*(1-z)
        if (.not. all(abs([left(1),right(1)] - 2*expected) < tol)) error stop "Wrong face value"
        if (.not. all(abs([left(2),right(2)] - [-expected,expected]) < tol)) error stop "Wrong outward flux"
        if (abs(left(1)-right(1)) > tol .or. abs(left(2)+right(2)) > tol) error stop "Interface imbalance"
        flux_sum = flux_sum + dS(q)*[left(2),right(2)]
        if (max(abs(y-0.5_rk),abs(z-0.5_rk)) < tol) then
            print '(a,2f11.7)', "Face center, left/right field: ", left(1), right(1)
            print '(a,2f11.7)', "Face center, outward flux:     ", left(2), right(2)
        end if
    end do
    print '(a,2f11.7)', "Integrated outward flux:       ", flux_sum
    if (.not. all(abs(flux_sum - [-1._rk,1._rk]/36) < tol)) error stop "Wrong face integral"

    do p = 1, 2
        call patch(p)%create(res1=nx(p), res2=25, res3=25)
        call patch(p)%basis(Tgc=basis)
        if (.not. patch(p)%err%ok) error stop "Could not sample the patch"
        samples = matmul(basis, reshape(coeff(shared(p,:)), [size(shared,2),1]))
        call patch(p)%export_Xg(filename=trim(files(p)), point_data=samples, field_names=['u'])
    end do

contains

    function face_trace(volume, uvw, local_coeff, normal_x) result(trace)
        type(nurbs_volume), intent(inout) :: volume
        real(rk), intent(in) :: uvw(3), local_coeff(:), normal_x
        real(rk) :: trace(2)
        real(rk), allocatable :: values(:), derivative(:,:), tangent(:)
        call volume%derivative(uvw, dTgc=derivative, Tgc=values)
        tangent = volume%cmp_Xg(uvw, order=[1,0,0])
        trace(1) = dot_product(local_coeff, values)
        trace(2) = -normal_x * dot_product(local_coeff, derivative(:,1)) / tangent(1)
    end function
end program two_volume_poisson
```
