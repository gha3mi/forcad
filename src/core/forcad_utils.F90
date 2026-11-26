!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Numerical kernels and low-level geometry utilities used by ForCAD.
!!
!! This module contains the public building blocks behind the single- and
!! multipatch APIs. Real arguments use `rk`, and tensor-product indices are
!! flattened with direction 1 varying fastest. Mathematical subscripts in this
!! module are zero-based unless stated otherwise; Fortran array indices and
!! returned connectivity are one-based. A documented mathematical precondition
!! is not necessarily checked at run time: each procedure states whether an
!! invalid input is represented by a sentinel result or must be excluded by the
!! caller.
!!
!! **B-spline convention**
!!
!! Let \(n=n_c-1\). For degree \(p\), control count \(n_c\), and knot vector
!! \(U=\{u_0,\ldots,u_m\}\), a spline space requires
!! \(m=n+p+1\), equivalently \(m+1=n_c+p+1\). The active interval is
!! \([u_p,u_{n+1}]=[u_p,u_{n_c}]\). In Fortran indexing these bounds are
!! `knot(degree+1)` and `knot(nc+1)`. Basis values use the
!! Cox-de Boor recurrence and are right-continuous on interior spans; the final
!! active endpoint is included explicitly. Rational tensor bases are
!!
!! \[
!! N_{i,0}(u)=
!! \begin{cases}
!! 1,&u_i\le u<u_{i+1},\\
!! 0,&\text{otherwise},
!! \end{cases}
!! \qquad
!! N_{i,p}(u)=
!! \frac{u-u_i}{u_{i+p}-u_i}N_{i,p-1}(u)
!! +\frac{u_{i+p+1}-u}{u_{i+p+1}-u_{i+1}}N_{i+1,p-1}(u),
!! \]
!!
!! where a term with zero denominator is defined as zero. For a flattened
!! tensor index \(A\), rational normalization is
!!
!! \[
!! R_A(\boldsymbol\xi)=
!! \frac{N_A(\boldsymbol\xi)w_A}
!!      {\sum_B N_B(\boldsymbol\xi)w_B}.
!! \]
!!
!! For a derivative multi-index \(\boldsymbol\alpha\), the implementation uses
!! the quotient recurrence
!!
!! \[
!! D^{\boldsymbol\alpha}R_A
!! =\frac{1}{W}\left[
!! D^{\boldsymbol\alpha}(N_Aw_A)
!! -\sum_{\mathbf0<\boldsymbol\beta\le\boldsymbol\alpha}
!! {\boldsymbol\alpha\choose\boldsymbol\beta}
!! D^{\boldsymbol\beta}W\,
!! D^{\boldsymbol\alpha-\boldsymbol\beta}R_A
!! \right],
!! \qquad W=\sum_BN_Bw_B .
!! \]
!!
!! A polynomial B-spline derivative is zero when its order exceeds the degree
!! in that direction. Rational derivatives need not be zero above that degree:
!! the quotient recurrence combines denominator derivatives with lower-order
!! rational derivatives. Local-support variants return the first active
!! one-based control index together with only \(p+1\) values; use these variants
!! in assembly and evaluation hot paths.
!! Knot insertion and degree elevation produce a scalar operator
!! \(\mathbf S\) satisfying
!! \(\mathbf H_{\mathrm{new}}=\mathbf S\mathbf H_{\mathrm{old}}\) for
!! homogeneous controls \(\mathbf H=(w\mathbf P,w)\). Tensor-product callers
!! apply this one-dimensional operator by directional sweeps.
!! Rational sums first replace every active weight by
!! \(\widehat w_i=\operatorname{scale}(w_i,-e)\), where
!! \(e=\max_i\operatorname{exponent}(w_i)\). This exact radix-power change of
!! projective gauge leaves every rational quotient unchanged while preventing a
!! common very small or very large weight scale from exhausting the floating-
!! point exponent range.
!!
!! **Failure conventions**
!!
!! These kernels do not own a diagnostic object. Where documented, a rejected
!! shape or numerical breakdown returns an empty array, zero-sized map, NaN, or
!! false success flag. Other low-level routines deliberately assume valid
!! dimensions and metadata; violating such a precondition can cause a bounds or
!! shape error. Higher-level NURBS methods validate their public inputs and
!! translate kernel failures into their structured `err` state where supported.
!!
!! **Primary references**
!!
!! - M. G. Cox, "The Numerical Evaluation of B-Splines," *IMA Journal of
!!   Applied Mathematics* 10 (1972), 134-149.
!!   [doi:10.1093/imamat/10.2.134](https://doi.org/10.1093/imamat/10.2.134).
!! - C. de Boor, "On Calculating with B-Splines," *Journal of Approximation
!!   Theory* 6 (1972), 50-62.
!!   [doi:10.1016/0021-9045(72)90080-9](https://doi.org/10.1016/0021-9045(72)90080-9).
!! - W. Boehm, "Inserting New Knots into B-Spline Curves," *Computer-Aided
!!   Design* 12 (1980), 199-201.
!!   [doi:10.1016/0010-4485(80)90154-2](https://doi.org/10.1016/0010-4485(80)90154-2).
!! - H. Prautzsch, "Degree Elevation of B-Spline Curves," *Computer Aided
!!   Geometric Design* 1 (1984), 193-198.
!!   [doi:10.1016/0167-8396(84)90031-1](https://doi.org/10.1016/0167-8396(84)90031-1).
!! - W. Tiller, "Knot-Removal Algorithms for NURBS Curves and Surfaces,"
!!   *Computer-Aided Design* 24 (1992), 445-453.
!!   [doi:10.1016/0010-4485(92)90012-Y](https://doi.org/10.1016/0010-4485(92)90012-Y).
!! - L. Piegl and W. Tiller, *The NURBS Book*, 2nd ed., Springer, 1997.
!!   [doi:10.1007/978-3-642-59223-2](https://doi.org/10.1007/978-3-642-59223-2).
!! - IEEE, *IEEE Standard for Floating-Point Arithmetic*, IEEE Std 754-2019,
!!   2019. [IEEE 754-2019](https://standards.ieee.org/ieee/754/6210/).
!! - T. J. R. Hughes, J. A. Cottrell, and Y. Bazilevs, "Isogeometric
!!   Analysis: CAD, Finite Elements, NURBS, Exact Geometry and Mesh
!!   Refinement," *CMAME* 194 (2005), 4135-4195.
!!   [doi:10.1016/j.cma.2004.10.008](https://doi.org/10.1016/j.cma.2004.10.008).
!! - R. E. Tarjan, "Efficiency of a Good But Not Linear Set Union Algorithm,"
!!   *Journal of the ACM* 22 (1975), 215-225.
!!   [doi:10.1145/321879.321884](https://doi.org/10.1145/321879.321884).
module forcad_utils

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_quiet_nan, ieee_support_datatype, &
        ieee_value
    use, intrinsic :: iso_fortran_env, only: file_storage_size, int8, int32, int64, real64
    use forcad_kinds, only: rk
    use fordebug, only: debug

    implicit none

    private
    public :: geometry_jet, map_analysis_point, gradient_sensitivity, jacobian_log_sensitivity
    public :: axis_transfer, nested_edge_transfer, tensor_transfer
    public :: tensor_ansatz, tensor_ansatz_element, curve_frame, tensor_basis_jet, tensor_measure
    interface valid_trace_side
        module procedure valid_trace_scalar, valid_trace_vector
    end interface
    public basis_bernstein, basis_bspline, elemConn_C0, kron, ndgrid, compute_multiplicity, compute_knot_vector, &
        basis_bspline_der, basis_bspline_der_order, basis_bspline_der_order_active, basis_bspline_der_all_active, &
        geometry_derivative, tensor_basis_derivative_local, tensor_basis_derivatives2_local, tensor_basis_derivative_at, &
        insert_knot_A_5_1, insert_knot_periodic_A_5_1, findspan, &
        elevate_degree_A_5_9, hexahedron_Xc, tetragon_Xc, remove_knots_A_5_8, &
        elemConn_Cn, unique, rotation, basis_bspline_2der, det, inv, dyad, gauss_leg, export_vtk_legacy, solve, &
        solve_spd_banded, solve_least_squares_qr, sparse_left_matmul, repelem, linspace, fill_uniform, eye, kron_eye, &
        valid_knot_vector, knot_start, knot_end, active_knots, active_span_count, active_span, &
        active_knot_multiplicity, infer_knot_shape, knot_tolerance, map_parameter, valid_trace_side, periodic_topology, &
        rotate_points, transform_points, run_python_script, show_pyvista_singlepatch, show_pyvista_multipatch, boundary_index, &
        boundary_axes, boundary_layer_index, interpolatory_boundary_layer, disjoint_set_union, disjoint_set_map, structural_nonzero, &
        scaled_norm, scaled_product, scaled_difference, affine_knot_residual, jacobian_frame, map_basis_derivatives, &
        normalize_field_constraints, equilibrate_field_row, &
        measure_sensitivity, valid_weights, nonuniform_weights, slice_control_net, isoparametric_grid, restrict_control_net, &
        refine_control_lines, coarsen_control_lines, reduce_control_lines, &
        query_parameter_axis, normalize_rational_basis, refinement_operators, pack_refinement_controls, unpack_refinement_controls, &
        map_parameter_values, reparameterize_knots, reorient_patch, reverse_control_order, &
        nearest_sample_index, nearest_spline_point, &
        spline_element, spline_connectivity, spline_span_plan, bezier_element_data, transform_geometry, check_affine_map, resize_buffer, &
        refine_control_net, get_array


    !> Fine coefficients by row, coarse coefficients by column.
    type axis_transfer
        integer, allocatable :: ptr(:),col(:)
        real(rk), allocatable :: val(:)
    end type

    type :: element_rule
        integer :: nlocal=0, ncoord=0, ndim=0, nq=0
        integer, allocatable :: nodes(:)
        real(rk) :: bounds(2,3)
        real(rk), allocatable :: parameters(:,:), weights(:)
    end type

    type :: bernstein_coefficients
        real(rk), allocatable :: value(:,:)
    end type

    !> One active B-spline span expressed in its own Bernstein coordinate.
    type, public :: bernstein_span
        integer :: first
        real(rk) :: domain(2), bounds_knot(2)
        real(rk), private :: domain_width, span_width, correction(2)
        integer, private :: domain_power, span_power
        type(bernstein_coefficients), private, allocatable :: derivative(:)
    end type
    public :: prepare_bernstein_pair, eval_bernstein_pair, eval_bernstein_endpoint
    public :: interval_ratio, two_sum, add_scaled_term, weighted_tangent_measure, scaled_measure

    interface resize_buffer
        module procedure resize_buffer_1, resize_buffer_2, resize_buffer_3, resize_buffer_4, resize_buffer_5
    end interface resize_buffer

    interface get_array
        module procedure get_real_vector, get_real_matrix, get_integer_matrix
    end interface get_array

    !> Normalize polynomial values and optional gradients over a supplied weight list.
    !! Used for element-local weights by the existing low-level basis interfaces.
    interface normalize_rational_basis
        module procedure normalize_rational_basis_tensor, normalize_rational_basis_curve
    end interface

    !> Build one-based \(C^0\) tensor-product element connectivity.
    !!
    !! These constructors describe a regular positive-degree \(C^0\) layout:
    !! each directional count must satisfy `nnode=nelem*p+1` with `p>=1`.
    !! Overloads accept scalar counts for curves and directional counts for
    !! surfaces or volumes. Direction 1 varies fastest in every element row.
    !! Use [[elemConn_Cn]] or the rank-specific NURBS `cmp_elem` bindings for
    !! arbitrary knot vectors, repeated knots, or degree-zero spaces.
    interface elemConn_C0
        module procedure cmp_elemConn_C0_L
        module procedure cmp_elemConn_C0_S
        module procedure cmp_elemConn_C0_V
    end interface
    !===============================================================================


    !===============================================================================
    !> Build one-based connectivity from active knots and knot multiplicities.
    !!
    !! In each direction, `size(Xth)-1` supplies the number of active spans; the
    !! numerical values in `Xth` are not inspected. Interior entries
    !! `vecKnot_mul(2:nelem)` shift the first supported control index by the
    !! stored knot multiplicity. Entry 1 and the final active-endpoint
    !! multiplicity are not used. For a valid degree-\(p\) spline direction,
    !! the metadata must satisfy
    !!
    !! \[
    !! n_c=p+1+\sum_{i=2}^{n_e}s_i,
    !! \]
    !!
    !! where \(s_i\) is the multiplicity at the interior active knot preceding
    !! element \(i\). The routines check only basic array dimensions; callers
    !! must provide consistent positive multiplicities and control counts. At
    !! multiplicity \(s\), the spline space has guaranteed continuity
    !! \(C^{p-s}\). Returned support indices are one-based, with direction 1
    !! varying fastest.
    interface elemConn_Cn
        module procedure cmp_elemConn_Cn_L
        module procedure cmp_elemConn_Cn_S
        module procedure cmp_elemConn_Cn_V
    end interface
    !===============================================================================


    !===============================================================================
    !> Form a two- or three-dimensional tensor-product coordinate grid.
    !!
    !! Returned point rows follow Fortran ordering: the first input vector
    !! varies fastest.
    interface ndgrid
        module procedure ndgrid2
        module procedure ndgrid3
    end interface
    !===============================================================================


    !===============================================================================
    !> Count exact knot multiplicities for one knot or an array of knots.
    !!
    !! Callers that compare independently generated floating-point knots should
    !! first apply a tolerance policy; multiplicity itself follows knot-vector
    !! representation and therefore uses exact stored values.
    interface compute_multiplicity
        module procedure compute_multiplicity1
        module procedure compute_multiplicity2
    end interface
    !===============================================================================


    !===============================================================================
    !> Return stable unique integer or real values in first-occurrence order.
    interface unique
        module procedure unique_integer
        module procedure unique_real
    end interface
    !===============================================================================


    !===============================================================================
    !> Compute the dyadic product \(A_{ij}=u_i v_j\).
    interface dyad
        module procedure dyad_t1_t1
    end interface
    !===============================================================================


    !===============================================================================
    !> Generate tensor-product Gauss-Legendre points and weights.
    !!
    !! Despite its name, `degree` is the number of quadrature points minus one.
    !! Thus \(n_q=degree+1\) points integrate every univariate polynomial of
    !! degree at most \(2n_q-1=2\,degree+1\) exactly in exact arithmetic.
    !! Tensor products have this directional exactness separately in each
    !! coordinate.
    !!
    !! **Original reference:** C. F. Gauss, *Methodus nova integralium valores
    !! per approximationem inveniendi*, H. Dieterich, Gottingae, 1815.
    !! [Digitized original](https://gallica.bnf.fr/ark:/12148/bpt6k2412190).
    interface gauss_leg
        module procedure gauss_legendre_1D
        module procedure gauss_legendre_2D
        module procedure gauss_legendre_3D, gauss_legendre_tensor
    end interface
    !===============================================================================


    !===============================================================================
    !> Compute supported Kronecker products without calling `matmul`.
    !!
    !! Overloads cover vector-vector, vector-matrix, and a three-factor tensor
    !! product. Their flattening order is consistent with ForCAD control nets.
    interface kron
        module procedure kron_t1_t1
        module procedure kron_t1_t2, kron_t2_t2
        module procedure kron3
    end interface
    !===============================================================================


    !===============================================================================
    !> Evaluate first derivatives or a requested low-order local derivative table.
    !!
    !! The local overload also returns the first active control index and avoids
    !! allocating a dense `nc`-entry basis vector.
    interface basis_bspline_der
        module procedure basis_bspline_active_ders
        module procedure basis_bspline_der_A
        module procedure basis_bspline_der_B
    end interface
    !===============================================================================


    !===============================================================================
    !> Evaluate B-spline basis derivatives through second order.
    !!
    !! Optional overloads return the lower derivative orders in the same call.
    interface basis_bspline_2der
        module procedure basis_bspline_2der_A
        module procedure basis_bspline_2der_B
        module procedure basis_bspline_2der_C
    end interface
    !===============================================================================


    !===============================================================================
    !> Infer directional control counts from flattened tensor-product knot data.
    interface infer_knot_shape
        module procedure infer_knot_shape_2d
        module procedure infer_knot_shape_3d
    end interface
    !===============================================================================


    !===============================================================================
    !> Return flattened indices on a curve endpoint, surface edge, or volume face.
    !!
    !! `reverse`, `swap`, and `flip` options map interface orientation without
    !! copying the associated control net.
    interface boundary_index
        module procedure boundary_index_1d
        module procedure boundary_index_2d
        module procedure boundary_index_3d
    end interface
    !===============================================================================


    !===============================================================================
    !> Return flattened indices in a layer measured inward from a boundary.
    !!
    !! Layer zero is the boundary itself. Positive layers support arbitrary
    !! \(C^n\) and \(G^n\) interface stencils.
    interface boundary_layer_index
        module procedure boundary_layer_index_1d
        module procedure boundary_layer_index_2d
        module procedure boundary_layer_index_3d
    end interface
    !===============================================================================


contains

    !> Stream finite metric rows and compare distances as fraction*radix**power.
    !! Each coordinate pair is scaled before subtraction, so a huge common
    !! coordinate cannot erase a tiny residual in another coordinate. Only the
    !! shared query/point prefix participates; exact ties retain the first row.
    pure integer function nearest_sample_index(points, query) result(index)
        real(rk), intent(in) :: points(:,:), query(:)
        real(rk) :: difference, component, squares, distance, best_distance
        integer :: row, d, ndim, coordinate_power, component_power, distance_power
        integer :: power, best_power

        index=0
        ndim=min(size(points,2),size(query))
        if (ndim==0) return
        if (.not.all(ieee_is_finite(query(:ndim)))) return
        best_power=0
        best_distance=0.0_rk
        do row=1,size(points,1)
            if (.not.all(ieee_is_finite(points(row,:ndim)))) cycle
            squares=0.0_rk
            distance_power=0
            do d=1,ndim
                if (points(row,d)==query(d)) cycle
                if (points(row,d)==0.0_rk) then
                    coordinate_power=exponent(query(d))
                else if (query(d)==0.0_rk) then
                    coordinate_power=exponent(points(row,d))
                else
                    coordinate_power=max(exponent(points(row,d)),exponent(query(d)))
                end if
                difference=scale(points(row,d),-coordinate_power)-scale(query(d),-coordinate_power)
                component=fraction(difference)
                component_power=coordinate_power+exponent(difference)
                if (squares==0.0_rk) then
                    distance_power=component_power
                    squares=component*component
                else if (component_power>distance_power) then
                    squares=scale(squares,2*(distance_power-component_power))+component*component
                    distance_power=component_power
                else
                    squares=squares+scale(component,component_power-distance_power)**2
                end if
            end do
            if (squares==0.0_rk) then
                index=row
                return
            end if
            distance=fraction(squares)
            power=2*distance_power+exponent(squares)
            if (index==0 .or. power<best_power .or. (power==best_power .and. distance<best_distance)) then
                index=row
                best_power=power
                best_distance=distance
            end if
        end do
    end function nearest_sample_index


    !> Solve a small symmetric positive-definite Newton system on free variables.
    !! The matrix and negative gradient share one radix scale; active rows and
    !! columns become an identity block. Relative pivots below 64 epsilon, a
    !! nonfinite scaled right-hand side, or an unrepresentable solve reject this
    !! direction so the caller can try Gauss-Newton or projected descent instead.
    pure subroutine projection_direction(matrix, gradient, active, step, valid)
        real(rk), intent(in) :: matrix(:,:), gradient(:)
        logical, intent(in) :: active(:)
        real(rk), intent(out) :: step(:)
        logical, intent(out) :: valid
        real(rk) :: factor(3,3), system(3,3), rhs(3), y(3), solution(3)
        real(rk) :: matrix_scale, pivot, threshold
        integer :: n, i, j, power

        step=0.0_rk
        valid=.false.
        n=size(gradient)
        if (n<1 .or. n>3) return
        if (any(shape(matrix)/=[n,n]) .or. size(active)/=n .or. size(step)/=n) return
        if (.not.all(ieee_is_finite(matrix)) .or. .not.all(ieee_is_finite(gradient))) return
        if (all(active)) then
            valid=.true.
            return
        end if
        matrix_scale=0.0_rk
        do i=1,n
            if (active(i)) cycle
            do j=1,n
                if (.not.active(j)) matrix_scale=max(matrix_scale,abs(matrix(i,j)))
            end do
        end do
        if (matrix_scale==0.0_rk) return
        power=exponent(matrix_scale)
        system=0.0_rk
        rhs=0.0_rk
        do i=1,n
            if (active(i)) then
                system(i,i)=1.0_rk
                cycle
            end if
            if (gradient(i)/=0.0_rk) then
                if (exponent(gradient(i))-power>maxexponent(1.0_rk)) return
                rhs(i)=-scale(gradient(i),-power)
                if (rhs(i)==0.0_rk) return
            end if
            do j=1,n
                if (.not.active(j)) system(i,j)=scale(matrix(i,j),-power)
            end do
        end do
        if (.not.all(ieee_is_finite(system(:n,:n))) .or. .not.all(ieee_is_finite(rhs(:n)))) return
        threshold=64.0_rk*epsilon(1.0_rk)*maxval(abs(system(:n,:n)))
        if (any(abs(system(:n,:n)-transpose(system(:n,:n)))>threshold)) return
        factor=0.0_rk
        do i=1,n
            do j=1,i
                pivot=system(i,j)-dot_product(factor(i,:j-1),factor(j,:j-1))
                if (i==j) then
                    if (.not.ieee_is_finite(pivot) .or. pivot<=threshold) return
                    factor(i,j)=sqrt(pivot)
                else
                    factor(i,j)=pivot/factor(j,j)
                end if
            end do
            y(i)=(rhs(i)-dot_product(factor(i,:i-1),y(:i-1)))/factor(i,i)
            if (.not.ieee_is_finite(y(i))) return
        end do
        do i=n,1,-1
            solution(i)=(y(i)-dot_product(factor(i+1:n,i),solution(i+1:n)))/factor(i,i)
            if (.not.ieee_is_finite(solution(i))) return
        end do
        step=solution(:n)
        valid=.true.
    end subroutine projection_direction

    !> Stream span-aware seeds, then minimize distance with local spline jets.
    !! Work is support-sized; knots and physical offsets use independent radix gauges.
    !! Full breaks bound the selected continuous component. This is a local candidate,
    !! not a global certificate. Armijo (1966), Bertsekas (1982), Piegl--Tiller sec.6.1.
    pure subroutine nearest_spline_point(Xc,Wc,nc,degree,knot1,point,tol,maxit,nearest_Xt,nearest_Xg,err, &
        knot2,knot3,initial_Xt,wrap)
        real(rk), intent(in) :: Xc(:,:), knot1(:), point(:), tol
        real(rk), allocatable, intent(in) :: Wc(:)
        integer, intent(in) :: nc(:), degree(:), maxit
        real(rk), intent(out) :: nearest_Xt(:)
        real(rk), intent(out), optional :: nearest_Xg(:)
        type(debug), intent(inout) :: err
        real(rk), intent(in), optional :: knot2(:), knot3(:), initial_Xt(:)
        logical, intent(in), optional :: wrap(:)
        type :: parameter_axis
            real(rk), allocatable :: source(:), knots(:), seeds(:)
            real(rk) :: lo, hi
        end type
        type(parameter_axis) :: axis(3)
        type(debug) :: final_status
        real(rk), allocatable :: basis(:), gradient(:,:), hessian(:,:,:), ders1(:,:), ders2(:,:), ders3(:,:)
        real(rk), allocatable :: residual(:), trial_residual(:), jacobian(:,:), curvature(:,:,:), mapped(:)
        real(rk) :: x(size(nc)), trial(size(nc)), lower(size(nc)), upper(size(nc)), tolerance(size(nc))
        real(rk) :: g(size(nc)), H(size(nc),size(nc)), GGN(size(nc),size(nc)), step(size(nc)), projected(size(nc))
        real(rk) :: objective, trial_objective, best, alpha, slope, row_scale, width, norm_, fraction_
        real(rk) :: source_lower(size(nc)), source_upper(size(nc))
        integer :: ndim, ncoord, n(3), p(3), counts(3), index_, d, i, j, k, iteration, attempt, power, gauge, e, trial_gauge, score_power, best_power
        logical :: wrapped(3), cyclic(size(nc)), active(size(nc)), valid, accepted, converged, rational

        nearest_Xt=0.0_rk
        if (present(nearest_Xg)) nearest_Xg=0.0_rk
        if (.not.err%ok) return
        ndim=size(nc)
        ncoord=size(Xc,2)
        if (ndim<1 .or. ndim>3 .or. size(degree)/=ndim .or. size(nearest_Xt)/=ndim) then
            call projection_error(err,101,"Inconsistent spline parameter dimensions.")
            return
        end if
        if (any(nc<1) .or. any(degree<0) .or. any(degree>=nc) .or. &
            size(Xc,1)/=product(nc) .or. ncoord<1 .or. .not.all(ieee_is_finite(Xc))) then
            call projection_error(err,103,"Complete finite spline geometry is required.")
            return
        end if
        if (size(point)<ncoord) then
            call projection_error(err,101,"The query must cover every physical coordinate.")
            return
        end if
        if (present(nearest_Xg)) then
            if (size(nearest_Xg)/=ncoord) then
                call projection_error(err,101,"The output point must match the physical dimension.")
                return
            end if
        end if
        if (.not.all(ieee_is_finite(point(:ncoord))) .or. .not.ieee_is_finite(tol) .or. tol<0 .or. maxit<0) then
            call projection_error(err,100,"Use a finite query, nonnegative tolerance and iteration limit.")
            return
        end if
        rational=.false.
        if (allocated(Wc)) then
            if (size(Wc)/=size(Xc,1) .or. .not.valid_weights(Wc)) then
                call projection_error(err,103,"Spline weights must be finite and positive.")
                return
            end if
            rational=nonuniform_weights(Wc)
        end if
        n=1
        p=0
        n(:ndim)=nc
        p(:ndim)=degree
        wrapped=.false.
        if (present(wrap)) then
            if (size(wrap)/=ndim) then
                call projection_error(err,101,"Wrapping flags must match the parameter dimension.")
                return
            end if
            wrapped(:ndim)=wrap
        end if
        axis(1)%source=knot1
        if (present(knot2)) axis(2)%source=knot2
        if (present(knot3)) axis(3)%source=knot3
        cyclic=.false.
        do d=1,ndim
            if (.not.allocated(axis(d)%source)) then
                call projection_error(err,103,"Every parameter direction requires a knot vector.")
                return
            end if
            if (.not.valid_knot_vector(axis(d)%source,nc(d),degree(d))) then
                call projection_error(err,103,"Invalid spline knot vector.")
                return
            end if
            cyclic(d)=wrapped(d) .and. periodic_topology(axis(d)%source,degree(d),Xc,Wc,nc,d)
            do i=degree(d)+2,nc(d)-degree(d)
                if (axis(d)%source(i)==axis(d)%source(i+degree(d))) cyclic(d)=.false.
            end do
            axis(d)%lo=axis(d)%source(degree(d)+1)
            axis(d)%hi=axis(d)%source(nc(d)+1)
            call reparameterize_knots(axis(d)%source,nc(d),degree(d),0.0_rk,1.0_rk,axis(d)%knots,err)
            if (.not.err%ok) return
            call scaled_difference(axis(d)%lo,axis(d)%hi,width,e)
            tolerance(d)=1.0_rk
            if (tol==0) then
                tolerance(d)=0.0_rk
            else if (exponent(tol)-e-exponent(width)<=2) then
                tolerance(d)=min(1.0_rk,scale(fraction(tol)/fraction(width),exponent(tol)-e-exponent(width)))
            end if
            tolerance(d)=max(tolerance(d),64.0_rk*epsilon(1.0_rk))
        end do
        if (present(initial_Xt)) then
            if (size(initial_Xt)/=ndim) then
                call projection_error(err,101,"Initial parameters must match the parameter dimension.")
                return
            end if
            do d=1,ndim
                if (.not.ieee_is_finite(initial_Xt(d)) .or. &
                    initial_Xt(d)<axis(d)%lo .or. initial_Xt(d)>axis(d)%hi) then
                    call projection_error(err,107,"Initial parameters must be finite and inside the active domain.")
                    return
                end if
                call map_parameter_values([initial_Xt(d)],axis(d)%lo,axis(d)%hi,0.0_rk,1.0_rk,mapped,err)
                if (.not.err%ok) return
                x(d)=mapped(1)
            end do
        end if
        allocate(basis(0:product(p+1)-1),gradient(0:product(p+1)-1,ndim), &
            hessian(0:product(p+1)-1,ndim,ndim),residual(ncoord),trial_residual(ncoord), &
            jacobian(ncoord,ndim),curvature(ncoord,ndim,ndim), &
            ders1(0:2,0:p(1)),ders2(0:2,0:p(2)),ders3(0:2,0:p(3)))
        if (.not.present(initial_Xt)) then
            counts=1
            do d=1,ndim
                axis(d)%seeds=projection_seeds(axis(d)%knots,nc(d),degree(d))
                counts(d)=size(axis(d)%seeds)
            end do
            best=huge(1.0_rk)
            best_power=huge(1)
            do index_=0,product(counts)-1
                k=index_
                do d=1,ndim
                    trial(d)=axis(d)%seeds(mod(k,counts(d))+1)
                    k=k/counts(d)
                end do
                call evaluate(trial,0,.false.,trial_residual,valid,trial_gauge, &
                    basis,gradient,hessian,ders1,ders2,ders3,jacobian,curvature)
                if (.not.valid) exit
                norm_=scaled_norm(trial_residual)
                score_power=exponent(norm_)+trial_gauge
                fraction_=fraction(norm_)
                if (norm_==0.0_rk) score_power=-huge(1)
                if (score_power<best_power .or. (score_power==best_power .and. fraction_<best)) then
                    best_power=score_power
                    best=fraction_
                    x=trial
                end if
            end do
            if (.not.valid) then
                call projection_error(err,108,"A seed cannot be evaluated in the finite working range.")
                return
            end if
        end if
        do d=1,ndim
            if (wrapped(d) .and. x(d)==1.0_rk) x(d)=0.0_rk
        end do
        lower=0.0_rk
        upper=1.0_rk
        do d=1,ndim
            source_lower(d)=axis(d)%lo
            source_upper(d)=axis(d)%hi
            if (wrapped(d) .and. .not.cyclic(d)) then
                upper(d)=nearest(1.0_rk,-1.0_rk)
                source_upper(d)=nearest(axis(d)%hi,-1.0_rk)
            end if
        end do
        do d=1,ndim
            ! A full break separates independent one-sided geometry components.
            do i=degree(d)+2,nc(d)-degree(d)
                width=axis(d)%knots(i)
                if (width<=0 .or. width>=1) cycle
                if (axis(d)%knots(i+degree(d))/=width) cycle
                if (width<=x(d)) then
                    lower(d)=max(lower(d),width)
                    source_lower(d)=max(source_lower(d),axis(d)%source(i))
                else
                    upper(d)=min(upper(d),nearest(width,-1.0_rk))
                    source_upper(d)=min(source_upper(d),nearest(axis(d)%source(i),-1.0_rk))
                end if
            end do
        end do
        converged=.false.
        do iteration=0,maxit
            call evaluate(x,2,.false.,residual,valid,gauge, &
                basis,gradient,hessian,ders1,ders2,ders3,jacobian,curvature)
            if (.not.valid) then
                call projection_error(err,108,"Local spline derivatives exceed the finite working range.")
                exit
            end if
            objective=0.5_rk*dot_product(residual,residual)
            g=matmul(transpose(jacobian),residual)
            GGN=matmul(transpose(jacobian),jacobian)
            H=GGN
            do j=1,ndim
                do i=1,ndim
                    H(i,j)=H(i,j)+dot_product(residual,curvature(:,i,j))
                end do
            end do
            if (.not.all(ieee_is_finite(H)) .or. .not.all(ieee_is_finite(g))) then
                call projection_error(err,108,"Distance derivatives exceed the finite working range.")
                exit
            end if
            active=(x<=lower .and. g>=0) .or. (x>=upper .and. g<=0)
            where(cyclic) active=.false.
            projected=0.0_rk
            do d=1,ndim
                row_scale=max(maxval(abs(H(d,:))),maxval(abs(GGN(d,:))))
                if (.not.active(d) .and. row_scale>0) projected(d)=g(d)/row_scale
            end do
            if (maxit>0 .and. all(abs(projected)<=tolerance)) then
                converged=.true.
                exit
            end if
            if (iteration==maxit) exit
            accepted=.false.
            do attempt=1,3
                select case(attempt)
                case(1)
                    call projection_direction(H,g,active,step,valid)
                case(2)
                    call projection_direction(GGN,g,active,step,valid)
                case default
                    step=-projected
                    valid=all(ieee_is_finite(step))
                end select
                if (.not.valid) cycle
                trial=step
                step=max(lower,min(upper,x+step))-x
                where(cyclic) step=max(-0.5_rk,min(0.5_rk,trial))
                slope=dot_product(g,step)
                if (.not.ieee_is_finite(slope) .or. slope>=0) cycle
                alpha=1.0_rk
                do k=1,50
                    trial=max(lower,min(upper,x+alpha*step))
                    where(cyclic) trial=modulo(x+alpha*step,1.0_rk)
                    call evaluate(trial,0,.false.,trial_residual,valid,trial_gauge, &
                        basis,gradient,hessian,ders1,ders2,ders3,jacobian,curvature)
                    if (valid) then
                        norm_=scaled_norm(trial_residual)
                        power=exponent(norm_)+trial_gauge-gauge
                        if (norm_==0.0_rk) power=0
                        if (power<=maxexponent(1.0_rk)/2) then
                            norm_=scale(fraction(norm_),power)
                            trial_objective=0.5_rk*norm_*norm_
                            if (trial_objective<=objective+1.0e-4_rk*alpha*slope) then
                                accepted=.true.
                                exit
                            end if
                        end if
                    end if
                    alpha=0.5_rk*alpha
                end do
                if (accepted) exit
            end do
            if (.not.accepted) exit
            x=trial
        end do
        if (err%ok .and. .not.converged) call projection_error(err,109,"Nearest-point iteration did not converge.")
        ! Re-evaluate the literal returned parameter with the source trace convention.
        do d=1,ndim
            call map_parameter_values([x(d)],0.0_rk,1.0_rk,axis(d)%lo,axis(d)%hi,mapped,err=final_status)
            if (.not.final_status%ok) then
                err=final_status
                return
            end if
            nearest_Xt(d)=max(source_lower(d),min(source_upper(d),mapped(1)))
            call map_parameter_values([nearest_Xt(d)],axis(d)%lo,axis(d)%hi,0.0_rk,1.0_rk,mapped,final_status)
            if (err%ok .and. abs(mapped(1)-x(d))>tolerance(d)) then
              call projection_error(err,109,"The stationary parameter is not representable in the source domain.")
            end if
        end do
        call evaluate(nearest_Xt,0,.true.,residual,valid,gauge, &
            basis,gradient,hessian,ders1,ders2,ders3,jacobian,curvature)
        if (.not.valid .and. err%ok) call projection_error(err,108,"The returned parameter cannot be evaluated.")
        if (valid .and. present(nearest_Xg)) nearest_Xg=residual

    contains

        pure subroutine evaluate(parameter,order,original,offset,ok,offset_power,basis,gradient,hessian, &
            ders1,ders2,ders3,jacobian,curvature)
            real(rk), intent(in) :: parameter(:)
            integer, intent(in) :: order
            logical, intent(in) :: original
            real(rk), intent(out) :: offset(:)
            logical, intent(out) :: ok
            integer, intent(out) :: offset_power
            real(rk), intent(out), contiguous :: basis(0:), gradient(0:,:), hessian(0:,:,:), &
                ders1(0:,0:), ders2(0:,0:), ders3(0:,0:), jacobian(:,:), curvature(:,:,:)
            real(rk) :: component, delta, value
            real(rk) :: values1(0:0,0:p(1)), values2(0:0,0:p(2)), values3(0:0,0:p(3))
            integer :: first(3), a, b, local, global, anchor, l1, l2, l3, coordinate, exponent_
            ders1=0.0_rk
            ders2=0.0_rk
            ders3=0.0_rk
            ders2(0,0)=1.0_rk
            ders3(0,0)=1.0_rk
            first=1
            call axis_jet(1,parameter(1),order,original,first(1),ders1)
            if (ndim>=2) call axis_jet(2,parameter(2),order,original,first(2),ders2)
            if (ndim==3) call axis_jet(3,parameter(3),order,original,first(3),ders3)
            if (order==0) then
                values1(0,:)=ders1(0,:)
                values2(0,:)=ders2(0,:)
                values3(0,:)=ders3(0,:)
                if (rational) then
                    call tensor_basis_derivative_local(first(1),n(1),values1,first(2),n(2),values2, &
                        first(3),n(3),values3,basis,Wc)
                else
                    call tensor_basis_derivative_local(first(1),n(1),values1,first(2),n(2),values2, &
                        first(3),n(3),values3,basis)
                end if
            else if (rational) then
                call tensor_basis_derivatives2_local(first(1),n(1),ders1,first(2),n(2),ders2, &
                    first(3),n(3),ders3,basis,gradient,hessian,Wc)
            else
                call tensor_basis_derivatives2_local(first(1),n(1),ders1,first(2),n(2),ders2, &
                    first(3),n(3),ders3,basis,gradient,hessian)
            end if
            ok=all(ieee_is_finite(basis)) .and. abs(sum(basis)-1.0_rk)<=256.0_rk*epsilon(1.0_rk)*size(basis)
            offset_power=-huge(1)
            do l3=0,p(3)
                do l2=0,p(2)
                    do l1=0,p(1)
                        global=first(1)+l1+n(1)*(first(2)+l2-1+n(2)*(first(3)+l3-1))
                        do coordinate=1,ncoord
                            call scaled_difference(point(coordinate),Xc(global,coordinate),delta,exponent_)
                            if (delta/=0) offset_power=max(offset_power,exponent_+exponent(delta))
                        end do
                    end do
                end do
            end do
            if (offset_power==-huge(1)) offset_power=0
            anchor=first(1)+n(1)*(first(2)-1+n(2)*(first(3)-1))
            offset=0.0_rk
            jacobian=0.0_rk
            curvature=0.0_rk
            do l3=0,p(3)
                do l2=0,p(2)
                    do l1=0,p(1)
                        local=l1+(p(1)+1)*(l2+(p(2)+1)*l3)
                        global=first(1)+l1+n(1)*(first(2)+l2-1+n(2)*(first(3)+l3-1))
                        do coordinate=1,ncoord
                            call scaled_difference(point(coordinate),Xc(global,coordinate),delta,exponent_)
                            component=scale(delta,exponent_-offset_power)
                            if (delta/=0 .and. component==0) ok=.false.
                            offset(coordinate)=offset(coordinate)+basis(local)*component
                            if (order==0) cycle
                            call scaled_difference(Xc(anchor,coordinate),Xc(global,coordinate),delta,exponent_)
                            component=scale(delta,exponent_-offset_power)
                            if (delta/=0 .and. component==0) ok=.false.
                            do a=1,ndim
                                jacobian(coordinate,a)=jacobian(coordinate,a)+gradient(local,a)*component
                                do b=1,ndim
                                    value=hessian(local,a,b)*component
                                    curvature(coordinate,a,b)=curvature(coordinate,a,b)+value
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            if (original) then
                offset=0.0_rk
                do coordinate=1,ncoord
                    exponent_=-huge(1)
                    do l3=0,p(3)
                        do l2=0,p(2)
                            do l1=0,p(1)
                                global=first(1)+l1+n(1)*(first(2)+l2-1+n(2)*(first(3)+l3-1))
                                local=l1+(p(1)+1)*(l2+(p(2)+1)*l3)
                                if (basis(local)==0 .or. Xc(global,coordinate)==0) cycle
                                value=fraction(basis(local))*fraction(Xc(global,coordinate))
                                exponent_=max(exponent_,exponent(basis(local))+exponent(Xc(global,coordinate))+exponent(value))
                            end do
                        end do
                    end do
                    if (exponent_==-huge(1)) cycle
                    do l3=0,p(3)
                        do l2=0,p(2)
                            do l1=0,p(1)
                                local=l1+(p(1)+1)*(l2+(p(2)+1)*l3)
                                global=first(1)+l1+n(1)*(first(2)+l2-1+n(2)*(first(3)+l3-1))
                                if (basis(local)==0 .or. Xc(global,coordinate)==0) cycle
                                value=fraction(basis(local))*fraction(Xc(global,coordinate))
                                component=scale(value,exponent(basis(local))+exponent(Xc(global,coordinate))-exponent_)
                                if (component==0) ok=.false.
                                offset(coordinate)=offset(coordinate)+component
                            end do
                        end do
                    end do
                    offset(coordinate)=scale(offset(coordinate),exponent_)
                end do
            end if
            ok=ok .and. all(ieee_is_finite(offset)) .and. all(ieee_is_finite(jacobian)) &
                .and. all(ieee_is_finite(curvature))
        end subroutine evaluate

        pure subroutine axis_jet(direction,parameter,order,original,first,ders)
            integer, intent(in) :: direction, order
            real(rk), intent(in) :: parameter
            logical, intent(in) :: original
            integer, intent(out) :: first
            real(rk), intent(out), contiguous :: ders(0:,0:)
            real(rk) :: value
            if (original) then
                value=map_parameter(parameter,axis(direction)%source,n(direction),p(direction),wrapped(direction))
                call basis_bspline_der_all_active(value,axis(direction)%source,n(direction),p(direction),order,first,ders)
            else
                value=map_parameter(parameter,axis(direction)%knots,n(direction),p(direction),wrapped(direction))
                call basis_bspline_der_all_active(value,axis(direction)%knots,n(direction),p(direction),order,first,ders)
            end if
        end subroutine axis_jet
    end subroutine nearest_spline_point

    !> Seeds on each nonzero span, including right traces and adjacent left parameters.
    pure function projection_seeds(knot,nc,degree) result(seeds)
        real(rk), intent(in) :: knot(:)
        integer, intent(in) :: nc, degree
        real(rk), allocatable :: seeds(:)
        real(rk), allocatable :: work(:)
        real(rk) :: candidates(3)
        integer :: i, j, count_
        allocate(work(3*size(knot)+1))
        count_=0
        do i=degree+1,nc
            if (knot(i+1)<=knot(i)) cycle
            candidates=[knot(i),0.5_rk*knot(i)+0.5_rk*knot(i+1),nearest(knot(i+1),-1.0_rk)]
            do j=1,3
                if (count_>0) then
                    if (candidates(j)<=work(count_)) cycle
                end if
                count_=count_+1
                work(count_)=candidates(j)
            end do
        end do
        count_=count_+1
        work(count_)=knot(nc+1)
        seeds=work(:count_)
    end function projection_seeds

    pure subroutine projection_error(err,code,message)
        type(debug), intent(inout) :: err
        integer, intent(in) :: code
        character(*), intent(in) :: message
        call err%set(code=code,severity=1,category="forcad_utils",location="nearest_point2",message=message, &
            suggestion="Check geometry, dimensions, parameter bounds and numerical scale; adjust the initial point or iteration limit.")
    end subroutine projection_error


    !> Validate and apply one affine map to owned controls and any sampled cache.
    !! All input shapes and finite map entries are checked before either net changes.
    pure subroutine transform_geometry(Xc, Xg, linear, translation, center, err, category)
        real(rk), allocatable, intent(inout) :: Xc(:,:), Xg(:,:)
        real(rk), intent(in), contiguous :: linear(:,:)
        real(rk), intent(in), contiguous, optional :: translation(:), center(:)
        type(debug), intent(inout) :: err
        character(len=*), intent(in) :: category
        integer :: ncoord

        if (.not.err%ok .or. .not.allocated(Xc)) return
        ncoord=size(Xc,2)
        call check_affine_map(ncoord,linear,translation,center,err,category,101,100)
        if (.not.err%ok) return
        if (allocated(Xg)) then
            if (size(Xg,2)/=ncoord) then
                call err%set(code=101,severity=1,category=category,location="transform", &
                    message="Cached and control coordinates are inconsistent.", &
                    suggestion="Recreate or clear sampled geometry before transformation.")
                return
            end if
        end if
        call transform_points(Xc,linear,translation,center)
        if (allocated(Xg)) call transform_points(Xg,linear,translation,center)
    end subroutine transform_geometry

    !> Check an affine change of physical coordinates before updating any control net.
    pure subroutine check_affine_map(ncoord,linear,translation,center,err,category,shape_code,value_code)
        integer, intent(in) :: ncoord,shape_code,value_code
        real(rk), intent(in) :: linear(:,:)
        real(rk), optional, intent(in) :: translation(:),center(:)
        type(debug), intent(inout) :: err
        character(*), intent(in) :: category

        if (size(linear,1)/=ncoord .or. size(linear,2)/=ncoord) then
            call err%set(code=shape_code,severity=1,category=category,location="transform", &
                message="Affine matrix shape must match the coordinate dimension.", &
                suggestion="Provide linear with shape [ncoord,ncoord].")
            return
        end if
        if (present(translation)) then
            if (size(translation)/=ncoord) then
                call err%set(code=shape_code,severity=1,category=category,location="transform", &
                    message="Translation size must match the coordinate dimension.", &
                    suggestion="Provide translation with size ncoord.")
                return
            end if
        end if
        if (present(center)) then
            if (size(center)/=ncoord) then
                call err%set(code=shape_code,severity=1,category=category,location="transform", &
                    message="Transformation-center size must match coordinates.", &
                    suggestion="Provide center with size ncoord.")
                return
            end if
        end if
        if (.not.all(ieee_is_finite(linear))) then
            call err%set(code=value_code,severity=1,category=category,location="transform", &
                message="Affine matrix entries must be finite.", &
                suggestion="Remove NaN and infinite entries from linear.")
            return
        end if
        if (present(translation)) then
            if (.not.all(ieee_is_finite(translation))) then
                call err%set(code=value_code,severity=1,category=category,location="transform", &
                    message="Translation entries must be finite.", &
                    suggestion="Remove NaN and infinite entries from translation.")
                return
            end if
        end if
        if (present(center)) then
            if (.not.all(ieee_is_finite(center))) then
                call err%set(code=value_code,severity=1,category=category,location="transform", &
                    message="Transformation-center entries must be finite.", &
                    suggestion="Remove NaN and infinite entries from center.")
                return
            end if
        end if
    end subroutine


    !> Local support in literal full knots, including repeated active endpoints.
    !! Valid knots/counts/degrees and optional planned spans are caller preconditions; no source net is copied.
    pure subroutine spline_element(nc,degree,knot1,ie,ids,bounds,err,knot2,knot3,spans)
        integer, intent(in) :: nc(:), degree(:), ie
        real(rk), contiguous, intent(in) :: knot1(:)
        real(rk), contiguous, optional, intent(in) :: knot2(:), knot3(:)
        integer, optional, intent(in) :: spans(:)
        integer, allocatable, intent(out) :: ids(:)
        real(rk), intent(out) :: bounds(2,size(nc))
        type(debug), intent(inout) :: err
        integer :: ne(size(nc)), first(size(nc)), selected(size(nc)), d, e

        allocate(ids(0))
        bounds=0.0_rk
        if (.not.err%ok) return
        ne=1
        selected=0
        if (present(spans)) then
            selected=spans
        else
            ne(1)=active_span_count(knot1,nc(1),degree(1))
            if (present(knot2)) ne(2)=active_span_count(knot2,nc(2),degree(2))
            if (present(knot3)) ne(3)=active_span_count(knot3,nc(3),degree(3))
        end if
        if (.not.present(spans) .and. (ie<1 .or. ie>product(ne) .or. any(ne<1))) then
            call err%set(code=100,severity=1,category="forcad_utils",location="bezier", &
                message="The requested active element does not exist.",suggestion="Use a one-based active element index.")
            return
        end if
        e=ie-1
        do d=1,size(nc)
            select case(d)
            case(1)
                call element_axis(knot1,nc(d),degree(d),mod(e,ne(d))+1,first(d),bounds(:,d),selected(d))
            case(2)
                call element_axis(knot2,nc(d),degree(d),mod(e,ne(d))+1,first(d),bounds(:,d),selected(d))
            case default
                call element_axis(knot3,nc(d),degree(d),mod(e,ne(d))+1,first(d),bounds(:,d),selected(d))
            end select
            e=e/ne(d)
        end do
        ids=tensor_control_ids(nc,degree+1,first)
    end subroutine spline_element

    pure subroutine element_axis(knot,nc,p,element,first,bounds,stored_span)
        real(rk), contiguous, intent(in) :: knot(:)
        integer, intent(in) :: nc,p,element,stored_span
        integer, intent(out) :: first
        real(rk), intent(out) :: bounds(2)
        integer :: span
        span=stored_span
        if (span==0) span=active_span(knot,nc,p,element)
        first=span-p
        bounds=knot(span:span+1)
    end subroutine element_axis

    !> Flatten a rectangular tensor support with direction one varying fastest.
    pure function tensor_control_ids(nc,counts,first) result(ids)
        integer, intent(in) :: nc(:), counts(:), first(:)
        integer, allocatable :: ids(:)
        integer :: i, d, local, stride
        allocate(ids(product(counts)),source=1)
        do i=1,size(ids)
            local=i-1
            stride=1
            do d=1,size(nc)
                ids(i)=ids(i)+(first(d)-1+mod(local,counts(d)))*stride
                local=local/counts(d)
                stride=stride*nc(d)
            end do
        end do
    end function tensor_control_ids

    !> Plan active full-knot spans once for whole-patch decomposition.
    pure function spline_span_plan(nc,degree,knot1,knot2,knot3) result(plan)
        integer, intent(in) :: nc(:), degree(:)
        real(rk), contiguous, intent(in) :: knot1(:)
        real(rk), contiguous, optional, intent(in) :: knot2(:), knot3(:)
        integer, allocatable :: plan(:,:)
        integer :: ne(size(nc))
        ne(1)=active_span_count(knot1,nc(1),degree(1))
        if (present(knot2)) ne(2)=active_span_count(knot2,nc(2),degree(2))
        if (present(knot3)) ne(3)=active_span_count(knot3,nc(3),degree(3))
        allocate(plan(maxval(ne),size(nc)),source=0)
        plan(:ne(1),1)=positive_spans(knot1,nc(1),degree(1))
        if (present(knot2)) plan(:ne(2),2)=positive_spans(knot2,nc(2),degree(2))
        if (present(knot3)) plan(:ne(3),3)=positive_spans(knot3,nc(3),degree(3))
    end function spline_span_plan

    pure function positive_spans(knot,nc,p) result(spans)
        real(rk), intent(in) :: knot(:)
        integer, intent(in) :: nc,p
        integer, allocatable :: spans(:)
        integer :: i
        spans=pack([(i,i=p+1,nc)],knot(p+2:nc+1)>knot(p+1:nc))
    end function positive_spans

    !> Full-knot connectivity; exterior and active endpoint multiplicities stay literal.
    pure function spline_connectivity(nc,degree,knot1,knot2,knot3) result(elem)
        integer, intent(in) :: nc(:), degree(:)
        real(rk), contiguous, intent(in) :: knot1(:)
        real(rk), contiguous, optional, intent(in) :: knot2(:), knot3(:)
        integer, allocatable :: elem(:,:)
        elem=span_connectivity(knot1,nc(1),degree(1))
        if (present(knot2)) elem=connectivity_product(elem,span_connectivity(knot2,nc(2),degree(2)),nc(1))
        if (present(knot3)) elem=connectivity_product(elem,span_connectivity(knot3,nc(3),degree(3)),nc(1)*nc(2))
    end function spline_connectivity

    pure function span_connectivity(knot,nc,p) result(elem)
        real(rk), contiguous, intent(in) :: knot(:)
        integer, intent(in) :: nc,p
        integer, allocatable :: elem(:,:)
        integer :: s, row, j
        allocate(elem(active_span_count(knot,nc,p),p+1))
        row=0
        do s=p+1,nc
            if (knot(s+1)<=knot(s)) cycle
            row=row+1
            elem(row,:)=[(s-p+j,j=0,p)]
        end do
    end function span_connectivity

    !> Exact local Bernstein extraction by insertion of the two element endpoints.
    !! C is weight-independent: N=C*B, H_bezier=transpose(C)*H_active.
    !! Geometry uses successive positive directional contractions; the tensor
    !! matrix is allocated only when requested. Parameter bounds remain literal.
    pure subroutine bezier_element_data(knot1,nc,degree,ie,Xc,Wc,make_piece,Xout,Wout,ids,bounds,err,C,knot2,knot3,rational,spans)
        real(rk), contiguous, intent(in) :: knot1(:)
        real(rk), contiguous, optional, intent(in) :: knot2(:), knot3(:)
        integer, intent(in) :: nc(:), degree(:), ie
        integer, optional, intent(in) :: spans(:)
        real(rk), intent(in) :: Xc(:,:)
        real(rk), allocatable, intent(in) :: Wc(:)
        logical, intent(in) :: make_piece
        logical, optional, intent(in) :: rational
        real(rk), allocatable, intent(out) :: Xout(:,:), Wout(:)
        integer, allocatable, intent(out) :: ids(:)
        real(rk), intent(out) :: bounds(2,size(nc))
        type(debug), intent(inout) :: err
        real(rk), allocatable, optional, intent(out) :: C(:,:)
        real(rk), allocatable :: A(:,:), controls(:,:), weights(:), tensor(:,:)
        integer :: d, p, span, offset, first, counts(size(nc))
        logical :: weighted

        if (present(C)) allocate(C(0,0))
        call spline_element(nc,degree,knot1,ie,ids,bounds,err,knot2,knot3,spans)
        if (.not.err%ok) return
        if (.not.make_piece .and. .not.present(C)) return
        weighted=.false.
        if (present(rational)) then
            weighted=rational
        else if (allocated(Wc)) then
            weighted=nonuniform_weights(Wc)
        end if
        if (make_piece) then
            Xout=Xc(ids,:)
            if (weighted) Wout=Wc(ids)
        end if
        counts=degree+1
        if (present(C)) C=reshape([1.0_rk],[1,1])
        do d=1,size(nc)
            p=degree(d)
            offset=product(nc(:d-1))
            first=mod((ids(1)-1)/offset,nc(d))+1
            span=first+p
            select case(d)
            case(1)
                call bezier_axis(knot1,p,span,A,err)
            case(2)
                call bezier_axis(knot2,p,span,A,err)
            case default
                call bezier_axis(knot3,p,span,A,err)
            end select
            if (.not.err%ok) exit
            if (present(C)) then
                tensor=kron(transpose(A),C)
                call move_alloc(tensor,C)
            end if
            if (make_piece) then
                call project_control_layers(Xout,Wout,counts,d,1,A,weighted,controls,weights,err)
                if (.not.err%ok) exit
                call move_alloc(controls,Xout)
                if (allocated(Wout)) deallocate(Wout)
                if (allocated(weights)) call move_alloc(weights,Wout)
            end if
        end do
        if (.not.err%ok) then
            if (present(C)) then
                deallocate(C)
                allocate(C(0,0))
            end if
        end if
    end subroutine bezier_element_data

    !> Compute (b-a)/(hi-lo) without overflowing a difference or its reciprocal.
    pure elemental real(rk) function interval_ratio(a,b,lo,hi) result(value)
        real(rk), intent(in) :: a,b,lo,hi
        real(rk) :: part,whole
        integer :: ep,ew
        call scaled_difference(a,b,part,ep)
        call scaled_difference(lo,hi,whole,ew)
        value = scale(fraction(part)/fraction(whole),ep+exponent(part)-ew-exponent(whole))
    end function

    !> Blossom the active basis at p-k left and k right knots to obtain Bernstein row k.
    pure subroutine bezier_axis(knot,p,span,A,err)
        real(rk), contiguous, intent(in) :: knot(:)
        integer, intent(in) :: p,span
        real(rk), allocatable, intent(out) :: A(:,:)
        type(debug), intent(inout) :: err
        real(rk) :: lo,hi,alpha(2),beta(2),mass
        integer :: r,j,k,endpoint
        if (all(knot(span-p+1:span)==knot(span)) .and. all(knot(span+1:span+p)==knot(span+1))) then
            A=eye(p+1)
            return
        end if
        allocate(A(p+1,p+1),source=0.0_rk)
        A(:,p+1)=1.0_rk
        do r=p,1,-1
            do j=r,p
                lo=knot(span-p+j)
                hi=knot(span+j-r+1)
                alpha=interval_ratio(lo,knot(span:span+1),lo,hi)
                beta=interval_ratio(knot(span:span+1),hi,lo,hi)
                do k=0,p
                    endpoint=merge(2,1,r>p-k)
                    mass=A(k+1,j+1)
                    A(k+1,j)=A(k+1,j)+beta(endpoint)*mass
                    A(k+1,j+1)=alpha(endpoint)*mass
                end do
            end do
        end do
        if (any(maxval(A,dim=1)==0.0_rk)) call err%set(code=108,severity=1,category="forcad_utils",location="bezier", &
            message="The Bernstein extraction coefficients are not representable.", &
            suggestion="Rescale the knot parameter domain before extraction.")
    end subroutine bezier_axis

    pure subroutine prepare_bernstein_span(knot,p,span,axis,resolved,order)
        real(rk), contiguous, intent(in) :: knot(:)
        integer, intent(in) :: p,span,order
        type(bernstein_span), intent(out) :: axis
        logical, intent(out) :: resolved
        real(rk) :: lo,hi,denom,factor,x
        real(rk), allocatable :: base(:,:),C(:,:)
        integer :: k,d,j,i,id,power,shift
        type(debug) :: status
        lo=knot(p+1)
        hi=knot(size(knot)-p)
        axis%first=span-p
        axis%domain=[lo,hi]
        axis%bounds_knot=knot(span:span+1)
        call scaled_difference(lo,hi,axis%domain_width,axis%domain_power)
        call scaled_difference(knot(span),knot(span+1),axis%span_width,axis%span_power)
        allocate(axis%derivative(0:min(order,p)))
        resolved=.false.
        do k=0,min(order,p)
            ! The kth derivative is a degree p-k polynomial on the same span.
            call bezier_axis(knot,p-k,span,base,status)
            if (.not.status%ok) return
            if (k==0) then
                call move_alloc(base,axis%derivative(0)%value)
                cycle
            end if
            allocate(C(p-k+1,0:p),source=0.0_rk)
            C(:,k:p)=base
            ! Differentiate the spline basis via its lower-degree basis, avoiding
            ! cancellation between nearly equal extracted control values on tiny spans.
            do d=p-k+1,p
                do j=p-d+1,p
                    id=span-p+j
                    call scaled_difference(knot(id),knot(id+d),denom,power)
                    factor=d*fraction(axis%domain_width)/fraction(denom)
                    shift=axis%domain_power+exponent(axis%domain_width)-power-exponent(denom)
                    do i=1,p-k+1
                        if (C(i,j)==0.0_rk) cycle
                        x=fraction(C(i,j))*factor
                        power=shift+exponent(C(i,j))+exponent(x)
                        if (power>maxexponent(x)) return
                        x=scale(fraction(x),power)
                        C(i,j-1)=C(i,j-1)-x
                        C(i,j)=x
                    end do
                end do
            end do
            if (.not.all(ieee_is_finite(C))) return
            call move_alloc(C,axis%derivative(k)%value)
        end do
        resolved=.true.
    end subroutine

    !> Cache extraction and the affine displacement of the paired ends in B's span units.
    pure subroutine prepare_bernstein_pair(ka,kb,p,span,reverse,a,b,resolved,order)
        real(rk), contiguous, intent(in) :: ka(:),kb(:)
        integer, intent(in) :: p,span
        logical, intent(in) :: reverse
        type(bernstein_span), intent(out) :: a,b
        logical, intent(out) :: resolved
        integer, optional, intent(in) :: order
        integer :: s,j,k,n,power,highest
        real(rk) :: start(2),finish(2),d
        n=size(kb)
        s=merge(n-span,span,reverse)
        highest=p
        if (present(order)) highest=order
        call prepare_bernstein_span(ka,p,span,a,resolved,0)
        if (.not.resolved) return
        call prepare_bernstein_span(kb,p,s,b,resolved,highest)
        if (.not.resolved) return
        start=[ka(p+1),kb(merge(n-p,p+1,reverse))]
        finish=[ka(n-p),kb(merge(p+1,n-p,reverse))]
        do k=0,1
            j=s+merge(1-k,k,reverse)
            call affine_knot_residual(start,finish,[ka(span+k),kb(j)],d,power,signed_result=.true.)
            b%correction(k+1)=scale(-d/(fraction(a%domain_width)*fraction(b%span_width)), &
                power-a%domain_power-exponent(a%domain_width)-b%span_power-exponent(b%span_width))
        end do
    end subroutine

    !> Evaluate each extracted derivative with its Bernstein basis on [0,1].
    pure subroutine eval_bernstein_span(axis,xi,ders,resolved)
        type(bernstein_span), intent(in) :: axis
        real(rk), intent(in) :: xi
        real(rk), intent(out) :: ders(0:,0:)
        logical, intent(out) :: resolved
        integer :: p,n,d,k
        real(rk) :: B(0:ubound(ders,1))
        p=ubound(ders,1)
        n=min(p,ubound(ders,2))
        resolved=n<=ubound(axis%derivative,1)
        if (.not.resolved) return
        ders=0.0_rk
        B=0.0_rk
        B(0)=1.0_rk
        do d=0,p
            k=p-d
            if (k<=n) ders(:,k)=matmul(B(:d),axis%derivative(k)%value)
            if (d==p) exit
            B(1:d+1)=(1.0_rk-xi)*B(1:d+1)+xi*B(:d)
            B(0)=(1.0_rk-xi)*B(0)
        end do
        resolved=all(ieee_is_finite(ders))
    end subroutine

    pure subroutine eval_bernstein_pair(a,b,reverse,xi,da,db,resolved)
        type(bernstein_span), intent(in) :: a,b
        logical, intent(in) :: reverse
        real(rk), intent(in) :: xi
        real(rk), intent(out) :: da(:,:),db(:,:)
        logical, intent(out) :: resolved
        real(rk) :: xb
        xb=merge(1.0_rk-xi,xi,reverse)+(1.0_rk-xi)*b%correction(1)+xi*b%correction(2)
        resolved=xb>0.0_rk .and. xb<1.0_rk
        if (.not.resolved) return
        call eval_bernstein_span(a,xi,da,resolved)
        if (resolved) call eval_bernstein_span(b,xb,db,resolved)
    end subroutine

    pure subroutine eval_bernstein_endpoint(knot,p,side,first,ders,resolved)
        real(rk), contiguous, intent(in) :: knot(:)
        integer, intent(in) :: p,side
        integer, intent(out) :: first
        real(rk), intent(out) :: ders(:,:)
        logical, intent(out) :: resolved
        type(bernstein_span) :: axis
        integer :: n,span
        logical :: at_min
        n=size(knot)-p-1
        at_min=mod(side,2)==1
        span=1+findspan(n-1,p,knot(merge(p+1,n+1,at_min)),knot)
        call prepare_bernstein_span(knot,p,span,axis,resolved,size(ders,2)-1)
        first=axis%first
        if (resolved) call eval_bernstein_span(axis,merge(0.0_rk,1.0_rk,at_min),ders,resolved)
    end subroutine

    !> Affine coordinates, including exterior knots, with exact endpoint assignment.
    !! Radix-scaled differences avoid overflowing an otherwise representable map.
    !! The internal target may descend; public reparameterization requires increasing bounds.
    pure subroutine map_parameter_values(values, lo, hi, a, b, mapped, err)
        real(rk), intent(in) :: values(:), lo, hi, a, b
        real(rk), allocatable, intent(out) :: mapped(:)
        type(debug), intent(inout) :: err
        real(rk) :: source, target, delta, term
        integer :: es, et, ex, power, addition, i

        if (.not.err%ok) return
        if (.not.all(ieee_is_finite([lo,hi,a,b])) .or. .not.all(ieee_is_finite(values))) then
            call err%set(code=107,severity=1,category="forcad_utils",location="map_parameter_values", &
                message="Affine parameter mapping requires finite coordinates.", &
                suggestion="Supply finite parameter values and bounds.")
            return
        end if
        if (lo>=hi .or. a==b) then
            call err%set(code=107,severity=1,category="forcad_utils",location="map_parameter_values", &
                message="Affine parameter mapping requires nonzero intervals.", &
                suggestion="Use lo < hi and distinct target endpoints.")
            return
        end if
        mapped = values
        if (a==lo .and. b==hi) return
        if (a==hi .and. b==lo) then
            mapped=reflected_parameter(lo,hi,values)
        else
            call scaled_difference(lo,hi,source,es)
            call scaled_difference(a,b,target,et)
            do i = 1, size(values)
                if (values(i)==lo) then
                    mapped(i) = a
                else if (values(i)==hi) then
                    mapped(i) = b
                else
                    call scaled_difference(lo,values(i),delta,ex)
                    term = fraction(delta)*fraction(target)/fraction(source)
                    power = ex+et-es+exponent(delta)+exponent(target)-exponent(source)
                    addition = max(exponent(a),power+exponent(term))
                    mapped(i) = scale(scale(a,-addition)+scale(term,power-addition),addition)
                end if
            end do
        end if
        if (.not.all(ieee_is_finite(mapped))) then
            call err%set(code=108,severity=1,category="forcad_utils",location="map_parameter_values", &
                message="Affine parameter mapping exceeds the finite coordinate range.", &
                suggestion="Choose a target interval whose complete mapped coordinates are representable.")
        end if
    end subroutine map_parameter_values

    !> Reflect u as lo+hi-u. Keep the rounding residuals when the midpoint cancels.
    pure elemental real(rk) function reflected_parameter(lo,hi,u) result(value)
        real(rk), intent(in) :: lo,hi,u
        real(rk) :: outer,inner,remainder,tail,low,high,twice,candidate
        ! Pair the outer terms first: opposite extremes can cancel before either sum overflows.
        call two_sum(min(lo,-u),max(hi,-u),outer,remainder)
        value=outer
        if (.not.ieee_is_finite(outer)) return
        call two_sum(max(lo,min(hi,-u)),remainder,inner,tail)
        value=inner
        if (.not.ieee_is_finite(inner)) return
        call two_sum(outer,inner,high,low)
        value=high
        if (.not.ieee_is_finite(high)) then
            ! A tail below the overflow rounding threshold can still give the last finite float.
            low=min(abs(outer),abs(inner))-(huge(high)-max(abs(outer),abs(inner)))
            if ((low-0.5_rk*spacing(huge(high)))+sign(1.0_rk,outer)*tail<0.0_rk) value=sign(huge(high),outer)
        else if (low==0.0_rk) then
            value=high+tail
        else
            if ((tail>0.0_rk .and. low>0.0_rk) .or. (tail<0.0_rk .and. low<0.0_rk)) then
                ! A same-sign tail breaks a half-ulp tie toward the adjacent float.
                twice=2.0_rk*low
                candidate=high+twice
                if (candidate-high==twice) value=candidate
            end if
        end if
    end function reflected_parameter

    !> FastTwoSum with its larger-magnitude operand first.
    pure subroutine two_sum(a,b,total,remainder)
        real(rk), intent(in) :: a,b
        real(rk), intent(out) :: total,remainder
        total=a+b
        if (abs(a)>=abs(b)) then
            remainder=b-(total-a)
        else
            remainder=a-(total-b)
        end if
    end subroutine two_sum

    !> Represent b-a as difference*radix**power without overflowing subtraction.
    pure elemental subroutine scaled_difference(a, b, difference, power)
        real(rk), intent(in) :: a, b
        real(rk), intent(out) :: difference
        integer, intent(out) :: power
        power = max(exponent(a),exponent(b))
        difference = scale(b,-power)-scale(a,-power)
    end subroutine scaled_difference

    !> Absolute determinant of three parameter pairs as mantissa times radix**power.
    !! Its zero identifies exact affine alignment, even for a discontinuous trace.
    !! Splitting each factor at half precision gives four exact products.
    !! FastTwoSum retains the rounding residual of each addition.
    !! Separate radix exponents keep both tiny and overflowing products comparable.
    pure subroutine affine_knot_residual(start,finish,point,mantissa,power,signed_result)
        real(rk), intent(in) :: start(2),finish(2),point(2)
        logical, optional, intent(in) :: signed_result
        real(rk), intent(out) :: mantissa
        integer, intent(out) :: power
        real(rk) :: x(3),y(3),parts(24),a(2),b(2)
        integer :: powers(24),i,j,k,u,v,shift,count

        x = [start(1),finish(1),point(1)]
        y = [start(2),finish(2),point(2)]
        count = 0
        do i = 1,3
            do j = 1,2
                k = 1+mod(i-1+j,3)
                if (x(i)==0.0_rk .or. y(k)==0.0_rk) cycle
                a = split_fraction(x(i))
                b = split_fraction(merge(y(k),-y(k),j==1))
                shift = exponent(x(i))+exponent(y(k))
                do u = 1,2
                    do v = 1,2
                        call add_scaled_term(a(u)*b(v),shift,parts,powers,count)
                    end do
                end do
            end do
        end do
        mantissa = 0.0_rk
        power = 0
        if (count==0) return
        mantissa = parts(count)
        power = powers(count)
        if (count>1) mantissa = mantissa + scale(parts(count-1),powers(count-1)-power)
        power = power+exponent(mantissa)
        mantissa = fraction(mantissa)
        if (present(signed_result)) then
            if (signed_result) return
        end if
        mantissa = abs(mantissa)
    contains
        pure function split_fraction(value) result(pair)
            real(rk), intent(in) :: value
            real(rk) :: pair(2)
            integer, parameter :: half = digits(1.0_rk)/2
            pair(1) = scale(anint(scale(fraction(value),half)),-half)
            pair(2) = fraction(value)-pair(1)
        end function

    end subroutine affine_knot_residual

    !> Insert a signed term into a nonoverlapping expansion, smallest parts first.
    pure subroutine add_scaled_term(term,shift,parts,powers,n)
        real(rk), intent(in) :: term
        integer, intent(in) :: shift
        real(rk), intent(inout) :: parts(:)
        integer, intent(inout) :: powers(:),n
        real(rk) :: carry,hi,lo,total,tail
        integer :: carry_power,eh,el,i,nlow
        if (term==0.0_rk) return
        carry = fraction(term)
        carry_power = shift+exponent(term)
        nlow = 0
        do i = 1,n
            hi = carry
            eh = carry_power
            lo = parts(i)
            el = powers(i)
            if (el>eh .or. hi==0.0_rk) then
                hi = parts(i)
                eh = powers(i)
                lo = carry
                el = carry_power
            end if
            if (lo/=0.0_rk .and. eh-el<=digits(lo)+1) then
                call two_sum(hi,scale(lo,el-eh),total,tail)
                lo = fraction(tail)
                el = eh+exponent(tail)
                hi = fraction(total)
                eh = eh+exponent(total)
            end if
            if (lo/=0.0_rk) then
                nlow = nlow+1
                parts(nlow) = lo
                powers(nlow) = el
            end if
            carry = hi
            carry_power = eh
        end do
        n = nlow
        if (carry==0.0_rk) return
        n = n+1
        parts(n) = carry
        powers(n) = carry_power
    end subroutine

    !> Transform a valid full knot vector, preserving every equality and strict gap.
    pure subroutine reparameterize_knots(knot, nc, degree, a, b, mapped, err, category)
        real(rk), intent(in) :: knot(:), a, b
        integer, intent(in) :: nc, degree
        real(rk), allocatable, intent(out) :: mapped(:)
        type(debug), intent(inout) :: err
        character(*), optional, intent(in) :: category
        logical :: repeated(size(knot)-1)
        integer :: n
        if (.not.err%ok) return
        if (present(category)) then
            if (.not.all(ieee_is_finite([a,b])) .or. a>=b) then
                call err%set(code=107,severity=1,category=category,location="reparameterize", &
                    message="The new parameter interval must be increasing.",suggestion="Provide finite bounds a < b.")
                return
            end if
        end if
        call map_parameter_values(knot,knot(degree+1),knot(nc+1),a,b,mapped,err)
        if (.not.err%ok) return
        n = size(knot)
        repeated = knot(:n-1)==knot(2:)
        if (b<a) then
            mapped = mapped(n:1:-1)
            repeated = repeated(n-1:1:-1)
        end if
        if (.not.valid_knot_vector(mapped,nc,degree) .or. &
            any(repeated .neqv. (mapped(:n-1)==mapped(2:)))) then
            call err%set(code=108,severity=1,category="forcad_utils",location="reparameterize_knots", &
                message="Affine mapping cannot preserve distinct knots in this interval.", &
                suggestion="Choose a target interval that preserves every distinct knot.")
        end if
    end subroutine reparameterize_knots

    !> Source row for zero-based local indices t is origin + dot_product(stride,t).
    !! Each signed entry selects a distinct stored axis; omitted axes start at their first layer.
    pure subroutine oriented_strides(nc,axes,origin,stride)
        integer, intent(in) :: nc(:),axes(:)
        integer, intent(out) :: origin,stride(:)
        integer :: d,axis
        origin=1
        do d=1,size(axes)
            axis=abs(axes(d))
            stride(d)=product(nc(:axis-1))
            if (axes(d)>0) cycle
            origin=origin+(nc(axis)-1)*stride(d)
            stride(d)=-stride(d)
        end do
    end subroutine oriented_strides

    !> Copy a tensor patch into signed axis order, including every exterior knot and stored weight.
    !! Source geometry and directional knot inputs/outputs must match the rank (1, 2 or 3).
    pure subroutine reorient_patch(nc,degree,axes,knot1,Xc,Wc,k1out,xout,wout,err,knot2,k2out,knot3,k3out)
        integer, intent(in) :: nc(:),degree(:),axes(:)
        real(rk), intent(in) :: knot1(:),Xc(:,:)
        real(rk), allocatable, intent(in) :: Wc(:)
        real(rk), optional, intent(in) :: knot2(:),knot3(:)
        real(rk), allocatable, intent(out) :: k1out(:),xout(:,:),wout(:)
        real(rk), allocatable, optional, intent(out) :: k2out(:),k3out(:)
        type(debug), intent(inout) :: err
        type :: knot_vector
            real(rk), allocatable :: values(:)
        end type
        type(knot_vector) :: knots(size(nc))
        integer :: d,extent(3),stride(3),origin,j,k,first,last,row

        if (.not.err%ok) return
        ! Count signed matches directly so the most negative integer is safe to reject.
        if (size(axes)/=size(nc) .or. any([(count(axes==d .or. axes==-d)/=1,d=1,size(nc))])) then
            call err%set(code=100,severity=1,category="forcad_utils",location="reparameterize", &
                message="Axes must be a signed permutation of the patch directions.",suggestion="Use each direction exactly once.")
            return
        end if
        call orient_knots(1,knot1,knots(1)%values,err)
        if (present(knot2)) call orient_knots(2,knot2,knots(2)%values,err)
        if (present(knot3)) call orient_knots(3,knot3,knots(3)%values,err)
        if (.not.err%ok) return
        call move_alloc(knots(abs(axes(1)))%values,k1out)
        if (present(k2out)) call move_alloc(knots(abs(axes(2)))%values,k2out)
        if (present(k3out)) call move_alloc(knots(abs(axes(3)))%values,k3out)

        extent=1
        stride=0
        extent(:size(nc))=nc(abs(axes))
        call oriented_strides(nc,axes,origin,stride(:size(nc)))
        allocate(xout(size(Xc,1),size(Xc,2)))
        if (allocated(Wc)) allocate(wout(size(Wc)))
        do k=0,extent(3)-1
            do j=0,extent(2)-1
                first=origin+j*stride(2)+k*stride(3)
                last=first+(extent(1)-1)*stride(1)
                row=extent(1)*(j+extent(2)*k)
                xout(row+1:row+extent(1),:)=Xc(first:last:stride(1),:)
                if (allocated(Wc)) wout(row+1:row+extent(1))=Wc(first:last:stride(1))
            end do
        end do
    contains
        pure subroutine orient_knots(axis,knot,out,status)
            integer, intent(in) :: axis
            real(rk), intent(in) :: knot(:)
            real(rk), allocatable, intent(out) :: out(:)
            type(debug), intent(inout) :: status
            if (.not.status%ok) return
            if (any(axes==-axis)) then
                call reparameterize_knots(knot,nc(axis),degree(axis),knot(nc(axis)+1),knot(degree(axis)+1),out,status)
            else
                out=knot
            end if
        end subroutine orient_knots
    end subroutine reorient_patch

    !> Source row for reversing one tensor axis; direction 1 varies fastest.
    pure function reverse_control_order(nc, dir) result(order)
        integer, intent(in) :: nc(:), dir
        integer, allocatable :: order(:)
        integer :: stride, i
        allocate(order(product(nc)))
        stride = product(nc(:dir-1))
        do i = 0, size(order)-1
            order(i+1) = i+1+stride*(nc(dir)-1-2*mod(i/stride,nc(dir)))
        end do
    end function reverse_control_order
    !> Copy owned real-vector data or return an allocated empty result.
    pure subroutine get_real_vector(input,output,ok)
        real(rk), allocatable, intent(in) :: input(:)
        real(rk), allocatable, intent(out) :: output(:)
        logical, intent(in) :: ok
        if (ok.and.allocated(input)) then
            output=input
        else
            allocate(output(0))
        end if
    end subroutine get_real_vector

    !> Copy owned real-matrix data without an intermediate return-value array.
    pure subroutine get_real_matrix(input,output,ok)
        real(rk), allocatable, intent(in) :: input(:,:)
        real(rk), allocatable, intent(out) :: output(:,:)
        logical, intent(in) :: ok
        if (ok.and.allocated(input)) then
            output=input
        else
            allocate(output(0,0))
        end if
    end subroutine get_real_matrix

    !> Copy owned connectivity or return an allocated empty result.
    pure subroutine get_integer_matrix(input,output,ok)
        integer, allocatable, intent(in) :: input(:,:)
        integer, allocatable, intent(out) :: output(:,:)
        logical, intent(in) :: ok
        if (ok.and.allocated(input)) then
            output=input
        else
            allocate(output(0,0))
        end if
    end subroutine get_integer_matrix


    !> Share homogeneous insertion/elevation, decoding and optional fixed-weight maps.
    !! Valid source geometry and direction are preconditions. Supply either t,
    !! or matching Xth/r arrays. Only requested maps allocate; identity requests
    !! leave Xout unallocated so callers retain their sampling caches.
    pure subroutine refine_control_net(knot,degree,Xc,Wc,nc,dir,knots,Xout,Wout,err,B,Bs,Xth,r,t)
        real(rk), intent(in), contiguous :: knot(:), Xc(:,:)
        real(rk), allocatable, intent(in) :: Wc(:)
        integer, intent(in) :: degree, nc(:), dir
        real(rk), allocatable, intent(out) :: knots(:), Xout(:,:), Wout(:)
        type(debug), intent(inout) :: err
        real(rk), allocatable, optional, intent(out) :: B(:,:), Bs(:,:)
        real(rk), optional, intent(in) :: Xth(:)
        integer, optional, intent(in) :: r(:), t
        real(rk), allocatable :: controls(:,:), elevated(:,:), A(:,:)
        integer :: counts(size(nc)), new_degree, new_count, weight_exponent
        character(len=14) :: operation
        logical :: need_map, identity, success, cyclic

        if (.not.err%ok) return
        operation="insert_knots"
        new_degree=degree
        if (present(t)) then
            operation="elevate_degree"
            if (t<0) then
                call err%set(code=100,severity=1,category="forcad_utils",location=trim(operation), &
                    message="Invalid degree elevation request.",suggestion="Use a nonnegative degree increment.")
                return
            end if
            new_degree=degree+t
            identity=t==0
        else
            if (size(Xth)/=size(r)) then
                call err%set(code=101,severity=1,category="forcad_utils",location=trim(operation), &
                    message="Knot insertion input sizes do not match.",suggestion="Use matching Xth and r arrays.")
                return
            end if
            if (any(.not.ieee_is_finite(Xth)) .or. any(r<0)) then
                call err%set(code=100,severity=1,category="forcad_utils",location=trim(operation), &
                    message="Invalid knot insertion request.",suggestion="Use finite knots and nonnegative insertion counts.")
                return
            end if
            identity=.not.any(r>0)
        end if
        if (dir<1 .or. dir>size(nc)) then
            call err%set(code=100,severity=1,category="forcad_utils",location=trim(operation), &
                message="Invalid refinement direction.",suggestion="Use a parameter direction of the patch.")
            return
        end if
        need_map=present(B).or.present(Bs)
        if (identity) then
            if (need_map) call refinement_operators(eye(nc(dir)),nc,nc,dir,size(Xc,2),Wc,Wc,B,Bs)
            return
        end if
        cyclic=periodic_topology(knot,degree,Xc,Wc,nc,dir)
        call pack_refinement_controls(Xc,Wc,nc,dir,controls,weight_exponent,err)
        if (.not.err%ok) return
        if (present(t)) then
            if (need_map) then
                if (cyclic) then
                    call elevate_degree_periodic_A_5_9(t,knot,degree,controls,new_count,knots,elevated,A,success)
                else
                    call elevate_degree_bounded_A_5_9(t,knot,degree,controls,new_count,knots,elevated,A,success)
                end if
            else
                if (cyclic) then
                    call elevate_degree_periodic_A_5_9(t,knot,degree,controls,new_count,knots,elevated,success=success)
                else
                    call elevate_degree_bounded_A_5_9(t,knot,degree,controls,new_count,knots,elevated,success=success)
                end if
            end if
            if (.not.success) then
                call err%set(code=108,severity=1,category="forcad_utils",location=trim(operation), &
                    message="Degree elevation failed to construct a valid directional transformation.", &
                    suggestion="Check the knot vector, degree and active interval.")
                return
            end if
            call move_alloc(elevated,controls)
        else
            knots=knot
            if (need_map) then
                call refine_control_lines(knots,degree,controls,Xth,r,err,A,bounded=.not.cyclic)
            else
                call refine_control_lines(knots,degree,controls,Xth,r,err,bounded=.not.cyclic)
            end if
            if (.not.err%ok) return
        end if
        counts=nc
        counts(dir)=size(controls,1)
        call unpack_refinement_controls(controls,counts,dir,size(Xc,2),weight_exponent,Xout,Wout)
        success=valid_knot_vector(knots,counts(dir),new_degree).and.all(ieee_is_finite(Xout))
        if (allocated(Wout)) success=success.and.valid_weights(Wout)
        if (.not.success) then
            call err%set(code=100,severity=1,category="forcad_utils",location=trim(operation), &
                message="Refinement produced invalid knots, coordinates or weights.", &
                suggestion="Rescale coordinates or relative weights before refinement.")
            return
        end if
        if (need_map) call refinement_operators(A,nc,counts,dir,size(Xc,2),Wc,Wout,B,Bs)
    end subroutine refine_control_net



    !> Refine packed control lines with one shared knot-insertion policy.
    !! Source data and finite, matching requests with nonnegative counts are prevalidated.
    !! Work arrays may change before failure; callers commit only when err is OK.
    !! Zero counts need no control workspace. Dense maps allocate only on request.
    !! Bounded mode cuts the stored spline rather than its periodic extension.
    pure subroutine refine_control_lines(knot, degree, controls, Xth, r, err, Tmap, bounded)
        real(rk), allocatable, intent(inout) :: knot(:), controls(:,:)
        integer, intent(in) :: degree, r(:)
        real(rk), intent(in) :: Xth(:)
        type(debug), intent(inout) :: err
        real(rk), allocatable, optional, intent(out) :: Tmap(:,:)
        logical, intent(in), optional :: bounded
        real(rk), allocatable :: knots(:), refined(:,:), step(:,:), composed(:,:)
        real(rk) :: lo, hi, tol
        integer :: nc, i, s, k, nq
        logical :: periodic

        if (.not.err%ok) return
        nc = size(knot)-degree-1
        if (present(Tmap)) Tmap = eye(nc)
        if (.not.any(r>0)) return
        periodic = .true.
        if (present(bounded)) periodic = .not.bounded
        if (periodic) periodic = periodic_topology(knot,degree,controls)
        do i = 1, size(Xth)
            if (r(i)==0) cycle
            nc = size(controls,1)
            lo = knot_start(knot,nc,degree)
            hi = knot_end(knot,nc,degree)
            tol = knot_tolerance(knot,degree+1,nc+1)
            if (Xth(i)<lo-tol .or. Xth(i)>hi+tol) then
                call err%set(code=100,severity=1,category="forcad_utils",location="insert_knots", &
                    message="Knot insertion is outside the active interval.",suggestion="Insert inside the active domain.")
                return
            end if
            s = compute_multiplicity(knot,Xth(i))
            if (s+r(i)>degree+merge(0,1,periodic)) then
                call err%set(code=100,severity=1,category="forcad_utils",location="insert_knots", &
                    message="Knot insertion exceeds the topology's multiplicity limit.", &
                    suggestion="Use final multiplicity <= degree+1, or <= degree for periodic geometry.")
                return
            end if
            k = findspan(nc-1,degree,Xth(i),knot)
            if (periodic) then
                if (present(Tmap)) then
                    call insert_knot_periodic_A_5_1(degree,knot,controls,Xth(i),k,s,r(i),nq,knots,refined,step)
                else
                    call insert_knot_periodic_A_5_1(degree,knot,controls,Xth(i),k,s,r(i),nq,knots,refined)
                end if
            else
                if (present(Tmap)) then
                    call insert_knot_A_5_1(degree,knot,controls,Xth(i),k,s,r(i),nq,knots,refined,step)
                else
                    call insert_knot_A_5_1(degree,knot,controls,Xth(i),k,s,r(i),nq,knots,refined)
                end if
            end if
            if (nq<0) then
                call err%set(code=108,severity=1,category="forcad_utils",location="insert_knots", &
                    message="Knot insertion failed to preserve spline topology.", &
                    suggestion="Check the knot value, multiplicity and periodic representation.")
                return
            end if
            if (present(Tmap)) then
                call sparse_left_matmul(step,Tmap,composed)
                call move_alloc(composed,Tmap)
            end if
            call move_alloc(knots,knot)
            call move_alloc(refined,controls)
        end do
    end subroutine refine_control_lines


    !> Retain an analysis work buffer when its shape is unchanged.
    pure subroutine resize_buffer_1(buffer,counts)
        real(rk), allocatable, intent(inout) :: buffer(:)
        integer, intent(in) :: counts(1)
        if (allocated(buffer)) then
            if (any(shape(buffer)/=counts)) deallocate(buffer)
        end if
        if (.not.allocated(buffer)) allocate(buffer(counts(1)))
    end subroutine resize_buffer_1

    !> Retain an analysis work buffer when its shape is unchanged.
    pure subroutine resize_buffer_2(buffer,counts)
        real(rk), allocatable, intent(inout) :: buffer(:,:)
        integer, intent(in) :: counts(2)
        if (allocated(buffer)) then
            if (any(shape(buffer)/=counts)) deallocate(buffer)
        end if
        if (.not.allocated(buffer)) allocate(buffer(counts(1),counts(2)))
    end subroutine resize_buffer_2

    !> Retain an analysis work buffer when its shape is unchanged.
    pure subroutine resize_buffer_3(buffer,counts)
        real(rk), allocatable, intent(inout) :: buffer(:,:,:)
        integer, intent(in) :: counts(3)
        if (allocated(buffer)) then
            if (any(shape(buffer)/=counts)) deallocate(buffer)
        end if
        if (.not.allocated(buffer)) allocate(buffer(counts(1),counts(2),counts(3)))
    end subroutine resize_buffer_3

    !> Retain an analysis work buffer when its shape is unchanged.
    pure subroutine resize_buffer_4(buffer,counts)
        real(rk), allocatable, intent(inout) :: buffer(:,:,:,:)
        integer, intent(in) :: counts(4)
        if (allocated(buffer)) then
            if (any(shape(buffer)/=counts)) deallocate(buffer)
        end if
        if (.not.allocated(buffer)) allocate(buffer(counts(1),counts(2),counts(3),counts(4)))
    end subroutine resize_buffer_4

    !> Retain an analysis work buffer when its shape is unchanged.
    pure subroutine resize_buffer_5(buffer,counts)
        real(rk), allocatable, intent(inout) :: buffer(:,:,:,:,:)
        integer, intent(in) :: counts(5)
        if (allocated(buffer)) then
            if (any(shape(buffer)/=counts)) deallocate(buffer)
        end if
        if (.not.allocated(buffer)) allocate(buffer(counts(1),counts(2),counts(3),counts(4),counts(5)))
    end subroutine resize_buffer_5

    !> Remove knots transactionally in packed directional control lines.
    !! An explicit physical tolerance replaces A5.8's coefficient equality test.
    !! Each candidate and source are clamped for comparison, then reinserted on
    !! the call-entry interior knots, so one budget
    !! controls the complete update, including every earlier accepted request.
    !! Bounds follow nonnegative partition of unity, with no sampling or dense map.
    pure subroutine coarsen_control_lines(knot,degree,Xc,Wc,nc,dir,Xth,r,err,controls,weight_exponent,counts,tol,error_bound)
        real(rk), allocatable, intent(inout) :: knot(:)
        real(rk), contiguous, intent(in) :: Xc(:,:)
        real(rk), allocatable, intent(in) :: Wc(:)
        integer, intent(in) :: degree, nc(:), dir, r(:)
        real(rk), allocatable, intent(out) :: controls(:,:)
        integer, intent(out) :: weight_exponent
        real(rk), intent(in) :: Xth(:)
        type(debug), intent(inout) :: err
        integer, allocatable, intent(out) :: counts(:)
        real(rk), optional, intent(in) :: tol
        real(rk), optional, intent(out) :: error_bound
        real(rk), allocatable :: entry_knots(:), entry(:,:), knots(:), candidate(:,:), restored_knots(:), restored(:,:), common(:,:)
        real(rk) :: bound, trial_bound, lo, hi
        integer :: i, attempt, ncontrol, ncoord, nline, multiplicity, k, accepted, history(size(r))
        logical :: measure, certified, cyclic
        type(debug) :: status

        allocate(counts(size(r)),source=0)
        if (present(error_bound)) error_bound=0.0_rk
        if (.not.err%ok) return
        if (size(Xth)/=size(r)) then
            call err%set(code=101,severity=1,category="forcad_utils",location="remove_knots", &
                message="Knot removal input sizes do not match.",suggestion="Use matching Xth and r arrays.")
            return
        end if
        if (present(tol)) then
            if (.not.ieee_is_finite(tol) .or. tol<0.0_rk) then
                call err%set(code=100,severity=1,category="forcad_utils",location="remove_knots", &
                    message="Invalid physical knot removal tolerance.",suggestion="Use finite tol >= 0 in physical units.")
                return
            end if
        end if
        if (any(.not.ieee_is_finite(Xth)) .or. any(r<0)) then
            call err%set(code=100,severity=1,category="forcad_utils",location="remove_knots", &
                message="Invalid knot removal request.",suggestion="Use finite knots and nonnegative removal counts.")
            return
        end if
        if (.not.any(r>0)) return
        cyclic=periodic_topology(knot,degree,Xc,Wc,nc,dir)
        call pack_refinement_controls(Xc,Wc,nc,dir,controls,weight_exponent,err)
        if (.not.err%ok) return
        ncoord=size(Xc,2)
        nline=product(nc)/nc(dir)
        measure=present(tol).or.present(error_bound)
        bound=0.0_rk
        if (measure) then
            call prepare_geometry_bound(knot,degree,controls,ncoord,nline,entry_knots,entry,err)
            if (.not.err%ok) return
        end if
        do i=1,size(Xth)
            if (r(i)==0) cycle
            ncontrol=size(controls,1)
            lo=knot_start(knot,ncontrol,degree)
            hi=knot_end(knot,ncontrol,degree)
            multiplicity=compute_multiplicity(knot,Xth(i))
            if (Xth(i)<=lo .or. Xth(i)>=hi .or. multiplicity==0 .or. r(i)>multiplicity) then
                call err%set(code=100,severity=1,category="forcad_utils",location="remove_knots", &
                    message="Requested interior knot multiplicity is unavailable.", &
                    suggestion="Use exact stored interior knots and r <= their current multiplicity.")
                return
            end if
            do attempt=1,r(i)
                k=findspan(size(controls,1)-1,degree,Xth(i),knot)+1
                multiplicity=compute_multiplicity(knot,Xth(i))
                if (cyclic) then
                    call remove_knots_periodic_A_5_8(degree,knot,controls,Xth(i),k,multiplicity,merge(1,r(i),present(tol)), &
                        accepted,knots,candidate,approximate=present(tol))
                else
                    call remove_knots_bounded_A_5_8(degree,knot,controls,Xth(i),k,multiplicity,merge(1,r(i),present(tol)), &
                        accepted,knots,candidate,approximate=present(tol))
                end if
                if (accepted==0) exit
                if (.not.all(ieee_is_finite(candidate))) exit
                if (size(candidate,2)>ncoord*nline) then
                    if (any(candidate(:,ncoord*nline+1:)<=0.0_rk)) exit
                end if
                if (measure) then
                    call roundtrip_control_lines(candidate,ncoord,nline,weight_exponent)
                    if (.not.all(ieee_is_finite(candidate))) exit
                    if (size(candidate,2)>ncoord*nline) then
                        if (any(candidate(:,ncoord*nline+1:)<=0.0_rk)) exit
                    end if
                    restored=evaluated_control_lines(candidate,ncoord,nline)
                    call restrict_spline(knots,degree,restored,lo,hi,restored_knots,common,success=certified)
                    if (.not.certified) exit
                    call move_alloc(common,restored)
                    history=counts
                    history(i)=history(i)+accepted
                    call status%reset()
                    call refine_control_lines(restored_knots,degree,restored,Xth(:i),history(:i),status,bounded=.true.)
                    if (.not.status%ok) exit
                    if (size(restored_knots)/=size(entry_knots)) exit
                    if (any(restored_knots/=entry_knots)) exit
                    trial_bound=control_geometry_bound(entry,restored,ncoord,nline)
                    if (.not.ieee_is_finite(trial_bound)) exit
                    if (present(tol)) then
                        if (trial_bound>tol) exit
                    end if
                    bound=trial_bound
                end if
                counts(i)=counts(i)+accepted
                call move_alloc(knots,knot)
                call move_alloc(candidate,controls)
                if (.not.present(tol)) exit
            end do
        end do
        if (present(error_bound)) error_bound=bound
    end subroutine coarsen_control_lines

    !> Degree reduction with one physical budget against the call-entry geometry.
    !! Local inverse elevation and candidate-only knot removal construct proposals;
    !! re-elevation on a common basis bounds their complete Euclidean deviation.
    pure subroutine reduce_control_lines(knot,degree,Xc,Wc,nc,dir,t,tol,err,controls,weight_exponent,reduced,error_bound)
        real(rk), allocatable, intent(inout) :: knot(:)
        integer, intent(in) :: degree,nc(:),dir,t
        real(rk), contiguous, intent(in) :: Xc(:,:)
        real(rk), intent(in) :: tol
        real(rk), allocatable, intent(in) :: Wc(:)
        type(debug), intent(inout) :: err
        real(rk), allocatable, intent(out) :: controls(:,:)
        integer, intent(out) :: weight_exponent,reduced
        real(rk), intent(out) :: error_bound
        real(rk), allocatable :: reference(:,:), entry_knots(:), entry(:,:), candidate(:,:), knots(:)
        real(rk), allocatable :: common(:,:), common_knots(:), breaks(:)
        integer, allocatable :: missing(:)
        integer :: ncoord,nline,p,n,i
        real(rk) :: trial_bound
        logical :: success,cyclic

        reduced=0
        error_bound=0.0_rk
        if (.not.err%ok) return
        if (t<0 .or. t>degree .or. .not.ieee_is_finite(tol) .or. tol<0.0_rk) then
            call err%set(code=100,severity=1,category="forcad_utils",location="reduce_degree", &
                message="Invalid degree reduction request.",suggestion="Use 0 <= t <= degree and finite physical tol >= 0.")
            return
        end if
        if (t==0) return
        cyclic=periodic_topology(knot,degree,Xc,Wc,nc,dir)
        call pack_refinement_controls(Xc,Wc,nc,dir,controls,weight_exponent,err)
        if (.not.err%ok) return
        ncoord=size(Xc,2)
        nline=product(nc)/nc(dir)
        call prepare_geometry_bound(knot,degree,controls,ncoord,nline,entry_knots,entry,err)
        if (.not.err%ok) return
        breaks=active_knots(entry_knots,size(entry,1),degree)
        do p=degree,degree-t+1,-1
            if (cyclic) then
                if (p==1) exit
                call reduce_degree_periodic(knot,p,controls,knots,candidate,success)
            else
                call reduce_degree_bounded(knot,p,controls,knots,candidate,success)
            end if
            if (.not.success) then
                call geometry_bound_error(err)
                exit
            end if
            if (.not.all(ieee_is_finite(candidate))) then
                call geometry_bound_error(err)
                exit
            end if
            if (size(candidate,2)>ncoord*nline) then
                if (any(candidate(:,ncoord*nline+1:)<=0.0_rk)) exit
            end if
            call roundtrip_control_lines(candidate,ncoord,nline,weight_exponent)
            if (.not.all(ieee_is_finite(candidate))) then
                call geometry_bound_error(err)
                exit
            end if
            if (size(candidate,2)>ncoord*nline) then
                if (any(candidate(:,ncoord*nline+1:)<=0.0_rk)) exit
            end if
            reference=evaluated_control_lines(candidate,ncoord,nline)
            call elevate_degree_bounded_A_5_9(degree-p+1,knots,p-1,reference,n,common_knots,common,success=success)
            if (success) then
                missing=[(count(entry_knots==breaks(i))-count(common_knots==breaks(i)),i=1,size(breaks))]
                call refine_control_lines(common_knots,degree,common,breaks,missing,err,bounded=.true.)
                success=err%ok .and. size(common_knots)==size(entry_knots)
            end if
            if (success) success=all(common_knots==entry_knots)
            if (.not.success) then
                call geometry_bound_error(err)
                exit
            end if
            trial_bound=control_geometry_bound(entry,common,ncoord,nline)
            if (.not.ieee_is_finite(trial_bound)) then
                call geometry_bound_error(err)
                exit
            end if
            if (trial_bound>tol) exit
            call move_alloc(knots,knot)
            call move_alloc(candidate,controls)
            reduced=reduced+1
            error_bound=trial_bound
        end do
        if (.not.err%ok) then
            reduced=0
            error_bound=0.0_rk
        end if
    end subroutine reduce_control_lines

    !> Prepare the common bounded reference for physical coarsening error control.
    pure subroutine prepare_geometry_bound(knot,degree,controls,ncoord,nline,knots,reference,err)
        real(rk), contiguous, intent(in) :: knot(:),controls(:,:)
        integer, intent(in) :: degree,ncoord,nline
        real(rk), allocatable, intent(out) :: knots(:),reference(:,:)
        type(debug), intent(inout) :: err
        real(rk), allocatable :: evaluated(:,:)
        logical :: success
        if (any(controls/=0.0_rk .and. abs(controls)<tiny(1.0_rk))) then
            call geometry_bound_error(err)
            return
        end if
        evaluated=evaluated_control_lines(controls,ncoord,nline)
        call restrict_spline(knot,degree,evaluated,knot(degree+1),knot(size(controls,1)+1),knots,reference,success=success)
        if (.not.success) call geometry_bound_error(err)
    end subroutine prepare_geometry_bound

    pure subroutine geometry_bound_error(err)
        type(debug), intent(inout) :: err
        call err%set(code=108,severity=1,category="forcad_utils",location="control_geometry_bound", &
            message="The physical geometry-error workspace is not representable.", &
            suggestion="Rescale coordinates, relative weights or the knot parameter domain.")
    end subroutine geometry_bound_error

    !> A5.11 construction: reduce independent Bezier pieces, then restore m-1 knots.
    !! Knot removal is candidate-only; the caller must test the complete physical error.
    pure subroutine reduce_degree_bounded(knot,p,controls,knots,reduced,success)
        real(rk), contiguous, intent(in) :: knot(:),controls(:,:)
        integer, intent(in) :: p
        real(rk), allocatable, intent(out) :: knots(:),reduced(:,:)
        logical, intent(out) :: success
        real(rk), allocatable :: work(:,:),work_knots(:),breaks(:),next(:,:),next_knots(:)
        integer, allocatable :: multiplicity(:)
        integer :: e,n,s,r,removed,first
        type(debug) :: status

        n=size(controls,1)
        call restrict_spline(knot,p,controls,knot(p+1),knot(n+1),work_knots,work,success=success)
        if (.not.success) return
        breaks=active_knots(work_knots,size(work,1),p)
        multiplicity=[(count(work_knots==breaks(e)),e=1,size(breaks))]
        call refine_control_lines(work_knots,p,work,breaks,p+1-multiplicity,status,bounded=.true.)
        success=status%ok
        if (.not.success) return
        allocate(reduced((size(breaks)-1)*p,size(controls,2)))
        do e=1,size(breaks)-1
            first=lbound(work,1)+(e-1)*(p+1)
            reduced((e-1)*p+1:e*p,:)=reduce_bezier_controls(work(first:first+p,:),p)
        end do
        knots=repelem(breaks,[(p,e=1,size(breaks))])
        do e=2,size(breaks)-1
            r=p+1-multiplicity(e)
            if (r==0) cycle
            s=count(knots==breaks(e))
            n=findspan(size(reduced,1)-1,p-1,breaks(e),knots)+1
            call remove_knots_bounded_A_5_8(p-1,knots,reduced,breaks(e),n,s,r,removed,next_knots,next,approximate=.true.)
            success=removed==r
            if (.not.success) return
            call move_alloc(next_knots,knots)
            call move_alloc(next,reduced)
        end do
        success=valid_knot_vector(knots,size(reduced,1),p-1)
    end subroutine reduce_degree_bounded

    !> Construct a cyclic reduced-degree candidate from guarded bounded work.
    !! Reinsert vanished seam orbits, then repeat the q trailing control rows.
    !! The caller must bound the final stored candidate against the original.
    !! Degree zero has no native cyclic representation in periodic_topology.
    pure subroutine reduce_degree_periodic(knot,p,controls,knots,reduced,success)
        real(rk), intent(in), contiguous :: knot(:),controls(:,:)
        integer, intent(in) :: p
        real(rk), allocatable, intent(out) :: knots(:),reduced(:,:)
        logical, intent(out) :: success
        real(rk), allocatable :: extended_knots(:),extended(:,:),work_knots(:),work(:,:)
        real(rk), allocatable :: candidate_knots(:),candidate(:,:),seams(:)
        integer, allocatable :: multiplicity(:),insertions(:)
        real(rk) :: a,b,period
        integer :: q,unique,unique_new,guards,n,i,index_,base,shift,first,last,start
        type(debug) :: status

        success=.false.
        allocate(knots(0),reduced(0,size(controls,2)))
        q=p-1
        if (q<1) return
        n=size(controls,1)
        unique=n-p
        a=knot(p+1)
        b=knot(n+1)
        period=b-a
        if (.not.ieee_is_finite(period)) return
        multiplicity=active_knot_multiplicity(knot,n,p)
        unique_new=sum(multiplicity(:size(multiplicity)-1)-1)
        if (multiplicity(1)==1) unique_new=unique_new+1
        guards=(q+unique_new-1)/unique_new+1
        n=(2*guards+1)*unique+p
        allocate(extended_knots(n+p+1),extended(n,size(controls,2)))
        do i=1,n
            extended(i,:)=controls(mod(i-1,unique)+1,:)
        end do
        do i=1,size(extended_knots)
            index_=i-1-guards*unique
            base=p+modulo(index_-p,unique)
            shift=(index_-base)/unique
            extended_knots(i)=knot(base+1)+real(shift,rk)*period
            ! Keep the stored central period, including its endpoint traces, literal.
            if (index_>=0 .and. index_<size(knot)) extended_knots(i)=knot(index_+1)
        end do
        if (.not.valid_knot_vector(extended_knots,n,p)) return
        call reduce_degree_bounded(extended_knots,p,extended,work_knots,work,success)
        if (.not.success) return
        success=.false.
        if (multiplicity(1)==1) then
            seams=[(a+real(i,rk)*period,i=-guards+1,guards)]
            seams(guards)=a
            seams(guards+1)=b
            if (.not.all(ieee_is_finite(seams))) return
            insertions=[(max(0,1-count(work_knots==seams(i))),i=1,size(seams))]
            call refine_control_lines(work_knots,q,work,seams,insertions,status,bounded=.true.)
            if (.not.status%ok) return
        end if
        first=findloc(work_knots,a,dim=1)
        last=findloc(work_knots,b,dim=1)
        start=first-q
        unique_new=last-first
        if (start<1 .or. unique_new<1 .or. last+q>size(work_knots)) return
        if (start+unique_new+q-1>size(work,1)) return
        candidate_knots=work_knots(lbound(work_knots,1)+start-1:lbound(work_knots,1)+last+q-1)
        candidate=work(lbound(work,1)+start-1:lbound(work,1)+start+unique_new+q-2,:)
        do i=unique_new+1,size(candidate,1)
            candidate(i,:)=candidate(mod(i-1,unique_new)+1,:)
        end do
        if (.not.periodic_topology(candidate_knots,q,candidate)) return
        call move_alloc(candidate_knots,knots)
        call move_alloc(candidate,reduced)
        success=.true.
    end subroutine reduce_degree_periodic

    !> Invert the Bernstein elevation recurrence inward from both endpoints.
    !! Average the middle estimates for odd p; p=1 uses a constant midpoint.
    pure function reduce_bezier_controls(Pc,p) result(Qc)
        integer, intent(in) :: p
        real(rk), intent(in) :: Pc(:,:)
        real(rk) :: Qc(p,size(Pc,2)),right(size(Pc,2)),alpha
        integer :: i
        if (p==1) then
            Qc(1,:)=0.5_rk*Pc(1,:)+0.5_rk*Pc(2,:)
            return
        end if
        Qc(1,:)=Pc(1,:)
        Qc(p,:)=Pc(p+1,:)
        do i=1,(p-1)/2
            alpha=real(i,rk)/real(p,rk)
            Qc(i+1,:)=(Pc(i+1,:)-alpha*Qc(i,:))/(1.0_rk-alpha)
            right=(Pc(p-i+1,:)-alpha*Qc(p-i+1,:))/(1.0_rk-alpha)
            if (i+1==p-i) then
                Qc(i+1,:)=0.5_rk*Qc(i+1,:)+0.5_rk*right
            else
                Qc(p-i,:)=right
            end if
        end do
    end function reduce_bezier_controls

    !> Include stored-weight quantization before testing or committing a candidate.
    pure subroutine roundtrip_control_lines(controls,ncoord,nline,weight_exponent)
        real(rk), intent(inout) :: controls(:,:)
        integer, intent(in) :: ncoord, nline, weight_exponent
        real(rk) :: weight, restored, point
        integer :: i, line, c, column

        if (size(controls,2)==ncoord*nline) return
        do line=1,nline
            do i=1,size(controls,1)
                weight=controls(i,ncoord*nline+line)
                restored=scale(scale(weight,weight_exponent),-weight_exponent)
                do c=1,ncoord
                    column=line+(c-1)*nline
                    point=controls(i,column)/weight
                    controls(i,column)=scale(fraction(point)*fraction(restored),exponent(point)+exponent(restored))
                end do
                controls(i,ncoord*nline+line)=restored
            end do
        end do
    end subroutine roundtrip_control_lines

    !> Match the existing evaluator's uniform-weight polynomial classification.
    pure function evaluated_control_lines(controls,ncoord,nline) result(evaluated)
        real(rk), intent(in) :: controls(:,:)
        integer, intent(in) :: ncoord, nline
        real(rk), allocatable :: evaluated(:,:)
        integer :: c
        evaluated=controls
        if (size(controls,2)==ncoord*nline) return
        if (maxval(controls(:,ncoord*nline+1:))/=minval(controls(:,ncoord*nline+1:))) return
        do c=1,ncoord
            evaluated(:,(c-1)*nline+1:c*nline)=controls(:,(c-1)*nline+1:c*nline)/controls(:,ncoord*nline+1:)
        end do
        evaluated(:,ncoord*nline+1:)=1.0_rk
    end function evaluated_control_lines

    !> Bound physical distance for two positive-weight nets on a common basis.
    !! With anchor b, candidate radius R, original (A,w), candidate (B,v):
    !! E = max_i (norm((A_i-w_i*b)-(B_i-v_i*b)) + R*abs(w_i-v_i))/w_i.
    !! Nonnegative tensor basis functions make this a global bound on the active
    !! domain and both traces of full breaks. Each net uses a minimum-weight-one
    !! projective gauge. Centering avoids translation-sized
    !! cancellation. Floating arithmetic is used, not interval certification.
    pure function control_geometry_bound(original,candidate,ncoord,nline) result(bound)
        real(rk), intent(in) :: original(:,:), candidate(:,:)
        integer, intent(in) :: ncoord, nline
        real(rk) :: bound, anchor(ncoord), point(ncoord), residual(ncoord), radius, w, v, ratio, value, wmin, vmin
        integer :: i, line, columns(ncoord)
        logical :: weighted

        bound=ieee_value(1.0_rk,ieee_quiet_nan)
        if (any(shape(original)/=shape(candidate))) return
        weighted=size(original,2)>ncoord*nline
        if (.not.all(ieee_is_finite(original)) .or. .not.all(ieee_is_finite(candidate))) return
        if (weighted) then
            if (any(original(:,ncoord*nline+1:)<=0.0_rk) .or. any(candidate(:,ncoord*nline+1:)<=0.0_rk)) return
        end if
        wmin=1.0_rk
        vmin=1.0_rk
        if (weighted) then
            wmin=minval(original(:,ncoord*nline+1:))
            vmin=minval(candidate(:,ncoord*nline+1:))
        end if
        w=1.0_rk
        if (weighted) w=candidate(1,ncoord*nline+1)
        anchor=candidate(1,1:ncoord*nline:nline)/w
        radius=0.0_rk
        do line=1,nline
            columns=[(line+(i-1)*nline,i=1,ncoord)]
            do i=1,size(candidate,1)
                v=1.0_rk
                if (weighted) v=candidate(i,ncoord*nline+line)
                point=candidate(i,columns)/v-anchor
                if (.not.all(ieee_is_finite(point))) return
                radius=max(radius,scaled_norm(point))
            end do
        end do
        if (.not.ieee_is_finite(radius)) return
        value=0.0_rk
        do line=1,nline
            columns=[(line+(i-1)*nline,i=1,ncoord)]
            do i=1,size(original,1)
                w=1.0_rk
                v=1.0_rk
                if (weighted) then
                    w=original(i,ncoord*nline+line)
                    v=candidate(i,ncoord*nline+line)
                end if
                ! Independent projective gauges put each minimum weight at one.
                ratio=scale((fraction(v)/fraction(w))*(fraction(wmin)/fraction(vmin)), &
                    exponent(v)-exponent(w)+exponent(wmin)-exponent(vmin))
                residual=(original(i,columns)/w-anchor)-ratio*(candidate(i,columns)/v-anchor)
                if (.not.all(ieee_is_finite(residual)) .or. .not.ieee_is_finite(ratio)) return
                value=max(value,scaled_norm(residual)+radius*abs(1.0_rk-ratio))
            end do
        end do
        bound=value
    end function control_geometry_bound


    !> Restrict validated polynomial or homogeneous control lines by full-multiplicity insertion.
    !! The lower endpoint retains its right trace and the upper its left trace.
    !! Original parameter values and interior knots are retained, including C^-1 breaks.
    pure subroutine restrict_spline(knot, degree, controls, a, b, knots, restricted, Tmap, success)
        real(rk), intent(in), contiguous :: knot(:), controls(:,:)
        integer, intent(in) :: degree
        real(rk), intent(in) :: a, b
        real(rk), allocatable, intent(out) :: knots(:), restricted(:,:)
        real(rk), allocatable, optional, intent(out) :: Tmap(:,:)
        logical, intent(out) :: success
        real(rk), allocatable :: work(:,:), work_knots(:), transform(:,:)
        integer :: first, last, lower, upper, r(2)
        type(debug) :: err

        success = .false.
        allocate(knots(0),restricted(0,size(controls,2)))
        work_knots = knot
        work = controls
        r = degree+1-[compute_multiplicity(knot,a),compute_multiplicity(knot,b)]
        if (present(Tmap)) then
            call refine_control_lines(work_knots,degree,work,[a,b],r,err,transform,bounded=.true.)
        else
            call refine_control_lines(work_knots,degree,work,[a,b],r,err,bounded=.true.)
        end if
        if (.not.err%ok) return
        first = findloc(work_knots,a,dim=1)
        last = findloc(work_knots,b,dim=1,back=.true.)
        lower = findloc(work_knots,a,dim=1,back=.true.)-degree
        upper = findloc(work_knots,b,dim=1)-1
        knots = work_knots(lbound(work_knots,1)+first-1:lbound(work_knots,1)+last-1)
        restricted = work(lbound(work,1)+lower-1:lbound(work,1)+upper-1,:)
        if (present(Tmap)) Tmap = transform(lower:upper,:)
        success = valid_knot_vector(knots,size(restricted,1),degree) .and. all(ieee_is_finite(restricted))
    end subroutine restrict_spline


    !> Restrict one tensor direction without sampling or a dense extraction map.
    !! Complete source geometry and a valid direction are caller preconditions.
    !! Output weights retain a common workspace gauge to avoid rounding small
    !! weights independently when restoring the source's arbitrary projective scale.
    pure subroutine restrict_control_net(Xc, Wc, nc, dir, knot, degree, a, b, knots, Xout, Wout, err)
        real(rk), intent(in), contiguous :: Xc(:,:), knot(:)
        real(rk), allocatable, intent(in) :: Wc(:)
        integer, intent(in) :: nc(:), dir, degree
        real(rk), intent(in) :: a, b
        real(rk), allocatable, intent(out) :: knots(:), Xout(:,:), Wout(:)
        type(debug), intent(inout) :: err
        real(rk), allocatable :: controls(:,:), restricted(:,:), polynomial_weights(:)
        integer :: weight_exponent, counts(size(nc))
        logical :: weighted, success

        if (.not.err%ok) return
        if (.not.ieee_is_finite(a) .or. .not.ieee_is_finite(b) .or. &
            a<knot(degree+1) .or. b>knot(nc(dir)+1) .or. a>=b) then
            call err%set(code=107,severity=1,category="forcad_utils",location="restrict", &
                message="Restriction requires finite bounds and a positive interval inside the active domain.", &
                suggestion="Use knot_start <= a < b <= knot_end without wrapping.")
            return
        end if
        weighted = .false.
        if (allocated(Wc)) weighted = nonuniform_weights(Wc)
        if (weighted) then
            call pack_refinement_controls(Xc,Wc,nc,dir,controls,weight_exponent,err)
        else
            call pack_refinement_controls(Xc,polynomial_weights,nc,dir,controls,weight_exponent,err)
        end if
        if (.not.err%ok) return
        call restrict_spline(knot,degree,controls,a,b,knots,restricted,success=success)
        if (success) then
            counts = nc
            counts(dir) = size(restricted,1)
            call unpack_refinement_controls(restricted,counts,dir,size(Xc,2),0,Xout,Wout)
            success = all(ieee_is_finite(Xout))
            if (allocated(Wout)) success = success .and. valid_weights(Wout)
            if (success) return
        end if
        call err%set(code=108,severity=1,category="forcad_utils",location="restrict", &
            message="Restriction could not represent finite controls and positive weights.", &
            suggestion="Check the knot data or rescale coordinates and relative weights.")
    end subroutine restrict_control_net

    !> Contract one tensor direction in homogeneous coordinates, retaining the others in order.
    !! The source geometry and direction must be valid. Only the active degree+1
    !! controls on each line are read; no full homogeneous net or dense operator
    !! is allocated. Coordinate numerators and weight denominators use separate
    !! radix gauges so normalizing a small weight cannot erase a finite control.
    !! Output weights share one representable gauge; lossy subnormal scaling is rejected.
    pure subroutine slice_control_net(Xc,Wc,nc,dir,knot,degree,Xt,wrap_parameter,Xout,Wout,err)
        real(rk), intent(in), contiguous :: Xc(:,:), knot(:)
        real(rk), intent(in) :: Xt
        real(rk), allocatable, intent(in) :: Wc(:)
        integer, intent(in) :: nc(:), dir, degree
        logical, intent(in) :: wrap_parameter
        real(rk), allocatable, intent(out) :: Xout(:,:), Wout(:)
        type(debug), intent(inout) :: err
        real(rk) :: N(0:degree)
        integer :: first
        logical :: weighted

        if (.not. err%ok) return
        call basis_bspline_der_order_active(map_parameter(Xt,knot,nc(dir),degree,wrap_parameter), &
            knot,nc(dir),degree,0,first,N)
        if (.not. all(ieee_is_finite(N))) then
            call err%set(code=108,severity=1,category="forcad_utils",location="slice", &
                message="The slice basis is outside the floating-point range.", &
                suggestion="Rescale the knot parameter domain.")
            return
        end if
        if (.not. ieee_is_finite(Xt) .or. sum(N)<=0.0_rk) then
            call err%set(code=107,severity=1,category="forcad_utils",location="slice", &
                message="The slice parameter must be finite and inside the active interval unless wrapping is enabled.", &
                suggestion="Use an active knot endpoint or an interior parameter.")
            return
        end if
        weighted=.false.
        if (allocated(Wc)) weighted=nonuniform_weights(Wc)
        call project_control_layers(Xc,Wc,nc,dir,first,reshape(N,[1,degree+1]),weighted,Xout,Wout,err)
    end subroutine slice_control_net

    !> Positive control-layer contraction with independent coordinate/weight gauges.
    !! Only the supplied layer window is read; output direction one remains fastest.
    pure subroutine project_control_layers(Xc,Wc,nc,dir,first,A,weighted,Xout,Wout,err)
        real(rk), intent(in) :: Xc(:,:), A(:,:)
        real(rk), allocatable, intent(in) :: Wc(:)
        integer, intent(in) :: nc(:), dir, first
        logical, intent(in) :: weighted
        real(rk), allocatable, intent(out) :: Xout(:,:), Wout(:)
        type(debug), intent(inout) :: err
        real(rk) :: weights(size(A,2)), factors(3), total, denominator, value
        integer :: indices(size(A,2)), stride, nline, line, row, point, a_, c, pass, power, largest, weight_power
        integer :: lower, upper, gauge
        integer, allocatable :: weight_powers(:)

        stride = product(nc(:dir-1))
        nline = product(nc)/nc(dir)
        allocate(Xout(nline*size(A,1),size(Xc,2)))
        if (weighted) allocate(Wout(size(Xout,1)),weight_powers(size(Xout,1)))
        do point=0,size(Xout,1)-1
            row=mod(point/stride,size(A,1))+1
            line=mod(point,stride)+stride*(point/(stride*size(A,1)))
            indices = mod(line,stride)+nc(dir)*stride*(line/stride)+1 &
                +stride*[(first+a_-2,a_=1,size(A,2))]
            weights = 1.0_rk
            if (weighted) weights = Wc(indices)
            denominator = 1.0_rk
            weight_power = 0
            do c = merge(0,1,weighted), size(Xc,2)
                largest = -huge(0)
                total = 0.0_rk
                do pass = 1, 2
                    do a_ = 1, size(A,2)
                        factors = [A(row,a_),weights(a_),1.0_rk]
                        if (c>0) factors(3) = Xc(indices(a_),c)
                        if (any(factors==0.0_rk)) cycle
                        power = sum(exponent(factors))
                        if (pass==1) then
                            largest = max(largest,power)
                        else
                            total = total+scale(product(fraction(factors)),power-largest)
                        end if
                    end do
                end do
                if (c==0) then
                    denominator = total
                    weight_power = largest
                    Wout(point+1) = fraction(total)
                    weight_powers(point+1) = largest+exponent(total)
                else
                    value = 0.0_rk
                    if (total/=0.0_rk) value = scale(fraction(total)/fraction(denominator), &
                        exponent(total)-exponent(denominator)+largest-weight_power)
                    Xout(point+1,c) = value
                end if
            end do
        end do
        decoding: block
            if (.not. all(ieee_is_finite(Xout))) exit decoding
            if (allocated(Wout)) then
                lower = maxval(weight_powers)-maxexponent(1.0_rk)
                upper = minval(weight_powers)-minexponent(1.0_rk)
                if (lower>upper) upper = upper+digits(1.0_rk)-1
                if (lower>upper) exit decoding
                gauge = max(lower,min(maxval(weight_powers),upper))
                do line = 1, size(Xout,1)
                    power = weight_powers(line)-gauge
                    value = scale(Wout(line),power)
                    if (.not.ieee_is_finite(value) .or. value<=0.0_rk) exit decoding
                    if (scale(value,-power)/=Wout(line)) exit decoding
                    Wout(line) = value
                end do
            end if
            return
        end block decoding
        deallocate(Xout)
        if (allocated(Wout)) deallocate(Wout)
        call err%set(code=108,severity=1,category="forcad_utils",location="project_control_layers", &
            message="The contracted controls or weights are outside the floating-point range.", &
            suggestion="Rescale the source coordinates or weights.")
    end subroutine project_control_layers


    !> Sample active isoparametric grid lines in direction order with a common point count.
    !! Valid source knot vectors and counts are required. Within each line family
    !! the lowest remaining parameter direction varies fastest, matching control order.
    pure subroutine isoparametric_grid(knot1, knot2, nc, degree, Xt, elemConn, res, knot3)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
        real(rk), intent(in), contiguous, optional :: knot3(:)
        integer, intent(in) :: nc(:), degree(:)
        real(rk), allocatable, intent(out) :: Xt(:,:)
        integer, allocatable, intent(out) :: elemConn(:,:)
        integer, intent(in), optional :: res
        type :: parameter_axis
            real(rk), allocatable :: knots(:), samples(:)
        end type parameter_axis
        type(parameter_axis) :: axis(size(nc))
        integer :: counts(size(nc)), spans(size(nc)), d, j, a, b, t, intervals, resolution, n, s, r, o, line, q, row, offset

        axis(1)%knots = active_knots(knot1,nc(1),degree(1))
        axis(2)%knots = active_knots(knot2,nc(2),degree(2))
        if (present(knot3)) axis(3)%knots = active_knots(knot3,nc(3),degree(3))
        do d = 1, size(nc)
            counts(d) = size(axis(d)%knots)
        end do
        spans = counts-1
        intervals = 1
        do d = 1, size(nc)
            a = intervals
            b = spans(d)
            do
                t = mod(a,b)
                if (t==0) exit
                a = b
                b = t
            end do
            intervals = (intervals/b)*spans(d)
        end do
        resolution = 10
        if (present(res)) resolution = max(2,res)
        intervals = intervals*(resolution-1)
        n = intervals+1
        do d = 1, size(nc)
            resolution = intervals/spans(d)
            allocate(axis(d)%samples(n))
            do s = 1, spans(d)
                o = (s-1)*resolution
                do r = 1, resolution
                    axis(d)%samples(o+r) = axis(d)%knots(s) &
                        +(axis(d)%knots(s+1)-axis(d)%knots(s))*real(r-1,rk)/real(resolution,rk)
                end do
            end do
            axis(d)%samples(n) = axis(d)%knots(counts(d))
        end do
        allocate(elemConn(sum(product(counts)/counts),n))
        allocate(Xt(size(elemConn,1)*n,size(nc)))
        row = 0
        do d = 1, size(nc)
            do line = 0, product(counts)/counts(d)-1
                offset = row*n
                q = line
                do j = 1, size(nc)
                    if (j==d) then
                        Xt(offset+1:offset+n,j) = axis(j)%samples
                    else
                        Xt(offset+1:offset+n,j) = axis(j)%knots(mod(q,counts(j))+1)
                        q = q/counts(j)
                    end if
                end do
                row = row+1
                elemConn(row,:) = [(offset+r,r=1,n)]
            end do
        end do
    end subroutine isoparametric_grid

    !> NURBS weights are nonempty, finite, and strictly positive.
    pure logical function valid_weights(Wc) result(valid)
        real(rk), intent(in) :: Wc(:)
        valid = size(Wc) > 0 .and. all(ieee_is_finite(Wc) .and. Wc > 0.0_rk)
    end function valid_weights

    !> Only exactly equal weights cancel from the rational basis.
    pure logical function nonuniform_weights(Wc) result(nonuniform)
        real(rk), intent(in) :: Wc(:)
        nonuniform = .false.
        if (size(Wc) == 0) return
        nonuniform = maxval(Wc)/=minval(Wc)
    end function nonuniform_weights

    !> Apply R=Q/W and R_,a=(Q_,a-R*W_,a)/W with exact radix-power scaling.
    !! Values scale the weighted contributions, so inactive large weights cannot
    !! erase endpoint interpolation. Derivatives retain their common weight gauge.
    !! Values and gradient rows must match the finite positive weight list.
    !! A zero denominator or invalid shape/weights produces zero outputs.
    pure subroutine normalize_rational_basis_tensor(basis, weights, gradient)
        real(rk), intent(inout) :: basis(:)
        real(rk), intent(in) :: weights(:)
        real(rk), intent(inout), optional :: gradient(:,:)
        real(rk) :: denominator, derivative, weight
        integer :: i, d, weight_exponent, anchor

        if (size(basis) /= size(weights) .or. .not. valid_weights(weights) .or. .not.all(ieee_is_finite(basis))) then
            basis = 0.0_rk
            if (present(gradient)) gradient = 0.0_rk
            return
        end if
        if (present(gradient)) then
            if (size(gradient,1) /= size(basis)) then
                basis = 0.0_rk
                gradient = 0.0_rk
                return
            end if
        end if
        weight_exponent = exponent(maxval(weights))
        if (.not.present(gradient)) then
            if (all(basis==0.0_rk)) return
            weight_exponent = -huge(0)
            do i = 1, size(basis)
                if (basis(i)==0.0_rk) cycle
                weight_exponent = max(weight_exponent,exponent(basis(i))+exponent(weights(i)) &
                    +exponent(fraction(basis(i))*fraction(weights(i))))
            end do
        end if
        do i = 1, size(basis)
            if (present(gradient)) then
                weight = scale(weights(i),-weight_exponent)
                basis(i) = basis(i)*weight
                gradient(i,:) = gradient(i,:)*weight
            else
                basis(i) = scale(fraction(basis(i))*fraction(weights(i)),exponent(basis(i))+exponent(weights(i))-weight_exponent)
            end if
        end do
        denominator = sum(basis)
        if (denominator <= 0.0_rk .or. .not. ieee_is_finite(denominator)) then
            basis = 0.0_rk
            if (present(gradient)) gradient = 0.0_rk
            return
        end if
        basis = basis/denominator
        if (present(gradient)) then
            anchor=maxloc(basis,dim=1)
            do d = 1, size(gradient,2)
                derivative = sum(gradient(:,d))
                gradient(:,d) = (gradient(:,d)-basis*derivative)/denominator
                if (basis(anchor)>0.5_rk) call balance_basis_derivative(gradient(:,d),anchor)
            end do
        end if
    end subroutine normalize_rational_basis_tensor

    !> Derivatives of a partition of unity sum to zero. Recover the dominant
    !! basis derivative from the remaining rows when its quotient cancels.
    pure subroutine balance_basis_derivative(values,anchor)
        real(rk), intent(inout) :: values(:)
        integer, intent(in) :: anchor
        values(anchor)=0.0_rk
        values(anchor)=-sum(values)
    end subroutine balance_basis_derivative

    !> View a curve gradient as one column, sharing the tensor calculation without a copy.
    pure subroutine normalize_rational_basis_curve(basis, weights, gradient)
        real(rk), intent(inout) :: basis(:)
        real(rk), intent(in) :: weights(:)
        real(rk), target, intent(inout) :: gradient(:)
        real(rk), pointer :: columns(:,:)

        columns(1:size(gradient),1:1) => gradient
        call normalize_rational_basis_tensor(basis,weights,columns)
    end subroutine normalize_rational_basis_curve


    !> Pack any tensor direction as rows, with orthogonal lines then coordinates in columns.
    !! Counts, controls and optional stored weights must describe a valid geometry.
    !! Unallocated weights select polynomial controls. One radix-power gauge is
    !! shared by the complete homogeneous net and restored by unpacking. Prefer
    !! a gauge keeping every nonzero homogeneous component normal; use subnormal
    !! range only when necessary. An unrepresentable workspace returns empty
    !! Xdir and optionally sets err, before any geometry is changed.
    pure subroutine pack_refinement_controls(Xc, Wc, nc, dir, Xdir, weight_exponent, err)
        real(rk), intent(in) :: Xc(:,:)
        real(rk), allocatable, intent(in) :: Wc(:)
        integer, intent(in) :: nc(:), dir
        real(rk), allocatable, intent(out) :: Xdir(:,:)
        integer, intent(out) :: weight_exponent
        type(debug), intent(inout), optional :: err
        integer :: stride, nline, point, row, line, c, ncoord, smallest, largest, power, lower, upper
        real(rk) :: weight, value

        stride = product(nc(:dir-1))
        nline = product(nc)/nc(dir)
        ncoord = size(Xc,2)
        weight_exponent = 0
        packing: block
            if (allocated(Wc)) then
                smallest = minval(exponent(Wc))
                largest = maxval(exponent(Wc))
                do c = 1, ncoord
                    do point = 1, size(Xc,1)
                        if (Xc(point,c)==0.0_rk) cycle
                        power = exponent(Wc(point))+exponent(Xc(point,c)) &
                            +exponent(fraction(Wc(point))*fraction(Xc(point,c)))
                        smallest = min(smallest,power)
                        largest = max(largest,power)
                    end do
                end do
                lower = largest-maxexponent(1.0_rk)
                upper = smallest-minexponent(1.0_rk)
                if (lower>upper) upper = upper+digits(1.0_rk)-1
                if (lower>upper) exit packing
                weight_exponent = max(lower,min(exponent(maxval(Wc)),upper))
            end if
            allocate(Xdir(nc(dir),nline*(ncoord+merge(1,0,allocated(Wc)))))
            do point = 1, size(Xc,1)
                row = mod((point-1)/stride,nc(dir))+1
                line = mod(point-1,stride)+1+stride*((point-1)/(stride*nc(dir)))
                if (allocated(Wc)) then
                    weight = scale(Wc(point),-weight_exponent)
                    if (.not.ieee_is_finite(weight) .or. weight<=0.0_rk) exit packing
                    Xdir(row,line+ncoord*nline) = weight
                end if
                do c = 1, ncoord
                    value = Xc(point,c)
                    if (allocated(Wc)) value = scale(fraction(value)*fraction(Wc(point)), &
                        exponent(value)+exponent(Wc(point))-weight_exponent)
                    if (.not.ieee_is_finite(value) .or. (value==0.0_rk .and. Xc(point,c)/=0.0_rk)) exit packing
                    Xdir(row,line+(c-1)*nline) = value
                end do
            end do
            return
        end block packing
        if (allocated(Xdir)) deallocate(Xdir)
        allocate(Xdir(0,0))
        if (present(err)) call err%set(code=108,severity=1,category="forcad_utils",location="pack_refinement_controls", &
            message="The homogeneous refinement workspace is outside the floating-point range.", &
            suggestion="Rescale coordinates or reduce the relative weight range before refinement.")
    end subroutine pack_refinement_controls


    !> Reverse the tensor layout and dehomogenize when a final weight component is present.
    !! The layout and positive homogeneous weights must match the supplied counts.
    !! Polynomial input leaves Wc unallocated, for forwarding to optional weight arguments.
    pure subroutine unpack_refinement_controls(Xdir, nc, dir, ncoord, weight_exponent, Xc, Wc)
        real(rk), intent(in) :: Xdir(:,:)
        integer, intent(in) :: nc(:), dir, ncoord, weight_exponent
        real(rk), allocatable, intent(out) :: Xc(:,:), Wc(:)
        integer :: stride, nline, point, row, line, c
        real(rk) :: weight
        logical :: weighted

        stride = product(nc(:dir-1))
        nline = product(nc)/nc(dir)
        weighted = size(Xdir,2)>ncoord*nline
        allocate(Xc(product(nc),ncoord))
        if (weighted) allocate(Wc(product(nc)))
        do point = 1, size(Xc,1)
            row = mod((point-1)/stride,nc(dir))+1
            line = mod(point-1,stride)+1+stride*((point-1)/(stride*nc(dir)))
            weight = 1.0_rk
            if (weighted) then
                weight = Xdir(row,line+ncoord*nline)
                Wc(point) = scale(weight,weight_exponent)
            end if
            do c = 1, ncoord
                Xc(point,c) = Xdir(row,line+(c-1)*nline)/weight
            end do
        end do
    end subroutine unpack_refinement_controls


    !> Lift a directional refinement map into the stored tensor control order.
    !! Bs=D_new^(-1)*A_tensor*D_old maps Euclidean controls at fixed weights;
    !! B=Bs (x) I_ncoord uses point-major coordinate packing. Unallocated weights
    !! represent a polynomial space. Counts and weights must match the map.
    !! Only requested outputs allocate; B alone does not require a dense Bs.
    pure subroutine refinement_operators(A, nc_old, nc_new, dir, ncoord, Wc_old, Wc_new, B, Bs)
        real(rk), intent(in) :: A(:,:)
        integer, intent(in) :: nc_old(:), nc_new(:), dir, ncoord
        real(rk), allocatable, intent(in) :: Wc_old(:), Wc_new(:)
        real(rk), allocatable, optional, intent(out) :: B(:,:), Bs(:,:)
        real(rk) :: value
        integer :: stride, blocks, block_, lane, i, j, row, col, c, nrow, ncol

        stride = product(nc_old(:dir-1))
        blocks = product(nc_old(dir+1:))
        nrow = product(nc_new)
        ncol = product(nc_old)
        if (present(B)) allocate(B(nrow*ncoord,ncol*ncoord),source=0.0_rk)
        if (present(Bs)) allocate(Bs(nrow,ncol),source=0.0_rk)
        do j = 1, size(A,2)
            do i = 1, size(A,1)
                if (A(i,j)==0.0_rk) cycle
                do block_ = 0, blocks-1
                    do lane = 1, stride
                        row = lane+stride*(i-1+size(A,1)*block_)
                        col = lane+stride*(j-1+size(A,2)*block_)
                        value = A(i,j)
                        if (allocated(Wc_old)) value = scale(fraction(value)*fraction(Wc_old(col))/fraction(Wc_new(row)), &
                            exponent(value)+exponent(Wc_old(col))-exponent(Wc_new(row)))
                        if (present(Bs)) Bs(row,col) = value
                        if (present(B)) then
                            do c = 1, ncoord
                                B((row-1)*ncoord+c,(col-1)*ncoord+c) = value
                            end do
                        end if
                    end do
                end do
            end do
        end do
    end subroutine refinement_operators

    ! Sparse refinement of whole clamped edges and tensor-product faces.

    !> Oslo refinement evaluates the coarse blossom at fine(i+1), ..., fine(i+p) to obtain Q_i.
    !! Both input vectors are ascending. CSR rows and columns retain their original
    !! fine and coarse storage orders, including when the edge orientation reverses.
    !! [Lyche, A note on the Oslo algorithm, Eq. 1](https://doi.org/10.1016/0010-4485%2888%2990118-2).
    pure subroutine nested_edge_transfer(coarse,fine,degree,rowptr,col,val,ok,reverse)
        real(rk), contiguous, intent(in) :: coarse(:),fine(:)
        integer, intent(in) :: degree
        integer, allocatable, intent(out) :: rowptr(:),col(:)
        real(rk), allocatable, intent(out) :: val(:)
        logical, intent(out) :: ok
        logical, optional, intent(in) :: reverse
        real(rk) :: k(size(coarse)),a(0:degree),ends(2,2),point(2),det,lo,hi,x,left,right
        integer :: nc,nf,i,j,h,m,power,span,r,b,first,last,step
        logical :: backward,same_domain

        ok=.false.
        nc=size(coarse)-degree-1
        nf=size(fine)-degree-1
        if (nf<nc .or. .not.valid_knot_vector(coarse,nc,degree) .or. .not.valid_knot_vector(fine,nf,degree)) return
        if (coarse(1)/=coarse(degree+1) .or. coarse(nc+1)/=coarse(size(coarse))) return
        if (fine(1)/=fine(degree+1) .or. fine(nf+1)/=fine(size(fine))) return
        backward=.false.
        if (present(reverse)) backward=reverse
        first=merge(size(coarse),1,backward)
        last=size(coarse)+1-first
        step=merge(-1,1,backward)
        ends(:,1)=[coarse(first),fine(1)]
        ends(:,2)=[coarse(last),fine(size(fine))]
        same_domain=all(ends(1,:)==ends(2,:))

        ! Embed the coarse knot multiset in the fine vector. Exact determinants
        ! distinguish adjacent representable knots even on reflected or rescaled axes.
        j=1
        do i=1,size(coarse)
            h=first+(i-1)*step
            do
                if (j>size(fine)) return
                if (same_domain) then
                    if (coarse(h)==fine(j)) exit
                else
                    point=[coarse(h),fine(j)]
                    if (all(point==ends(:,1)).or.all(point==ends(:,2))) exit
                    call affine_knot_residual(ends(:,1),ends(:,2),point,det,power)
                    if (det==0.0_rk) exit
                end if
                j=j+1
            end do
            k(i)=fine(j)
            j=j+1
        end do

        if (nf==nc) then
            rowptr=[(i,i=1,nf+1)]
            col=[(merge(nc+1-i,i,backward),i=1,nf)]
            allocate(val(nf),source=1.0_rk)
            ok=.true.
            return
        end if
        if (nf>(huge(nf)-1)/(degree+1)) return
        allocate(rowptr(nf+1),col(nf*(degree+1)),val(nf*(degree+1)))
        span=degree+1
        m=0
        do i=1,nf
            do while (k(span+1)<=fine(i))
                span=span+1
            end do
            a=0.0_rk
            a(degree)=1.0_rk
            do r=1,degree
                x=fine(i+r)
                do b=degree-r+1,degree
                    if (a(b)==0.0_rk) cycle
                    h=span-degree+b
                    lo=k(h)
                    hi=k(h+r)
                    left=a(b)*interval_ratio(x,hi,lo,hi)
                    right=a(b)*interval_ratio(lo,x,lo,hi)
                    if ((x<hi.and.left==0.0_rk).or.(x>lo.and.right==0.0_rk)) then
                        deallocate(rowptr,col,val)
                        return
                    end if
                    a(b-1)=a(b-1)+left
                    a(b)=right
                end do
            end do
            rowptr(i)=m+1
            do b=0,degree
                j=merge(degree-b,b,backward)
                if (a(j)==0.0_rk) cycle
                h=span-degree+j
                m=m+1
                col(m)=merge(nc+1-h,h,backward)
                val(m)=a(j)
            end do
        end do
        rowptr(nf+1)=m+1
        col=col(:m)
        val=val(:m)
        ok=.true.
    end subroutine

    !> Form T(i,j)=T1(i1,j1)*T2(i2,j2) from successful Oslo transfers.
    !! The axis list is in common tangent order. Rows follow the fine patch's native
    !! order and columns increase in the coarse patch's native order within each row.
    !! With one tangent, the existing CSR arrays move directly to the outputs.
    pure subroutine tensor_transfer(axis,ncoarse,swap_coarse,swap_fine,rowptr,col,val,ok)
        type(axis_transfer), intent(inout) :: axis(:)
        integer, intent(in) :: ncoarse(:)
        logical, intent(in) :: swap_coarse,swap_fine
        integer, allocatable, intent(out) :: rowptr(:),col(:)
        real(rk), allocatable, intent(out) :: val(:)
        logical, intent(out) :: ok
        integer :: c,f,nfine(2),i(2),nnz,row,fast,slow,m

        ok=.false.
        if (size(axis)==1) then
            call move_alloc(axis(1)%ptr,rowptr)
            call move_alloc(axis(1)%col,col)
            call move_alloc(axis(1)%val,val)
        else
            if (size(axis)/=2 .or. size(ncoarse)/=2) return
            nfine=[size(axis(1)%ptr)-1,size(axis(2)%ptr)-1]
            if (size(axis(1)%val)>(huge(nnz)-1)/size(axis(2)%val)) return
            nnz=size(axis(1)%val)*size(axis(2)%val)
            allocate(rowptr(product(nfine)+1),col(nnz),val(nnz))
            c=merge(2,1,swap_coarse)  ! Fast native coarse tangent.
            f=merge(2,1,swap_fine)    ! Fast native fine tangent.
            m=0
            do row=1,product(nfine)
                rowptr(row)=m+1
                i(f)=mod(row-1,nfine(f))+1
                i(3-f)=(row-1)/nfine(f)+1
                do slow=axis(3-c)%ptr(i(3-c)),axis(3-c)%ptr(i(3-c)+1)-1
                    do fast=axis(c)%ptr(i(c)),axis(c)%ptr(i(c)+1)-1
                        m=m+1
                        col(m)=axis(c)%col(fast)+ncoarse(c)*(axis(3-c)%col(slow)-1)
                        val(m)=axis(c)%val(fast)*axis(3-c)%val(slow)
                        if (val(m)==0.0_rk) then
                            deallocate(rowptr,col,val)
                            return
                        end if
                    end do
                end do
            end do
            rowptr(size(rowptr))=m+1
        end if
        ok=.true.
    end subroutine

    !> Explicit parameters take precedence over a resolution, then a cached grid.
    pure subroutine query_parameter_axis(knot, nc, degree, cached, parameters, valid, resolution, explicit)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: nc, degree
        real(rk), allocatable, intent(in) :: cached(:)
        real(rk), allocatable, intent(out) :: parameters(:)
        logical, intent(out) :: valid
        integer, intent(in), optional :: resolution
        real(rk), intent(in), optional :: explicit(:)
        valid = .false.
        if (present(explicit)) then
            parameters = explicit
        else if (present(resolution)) then
            if (resolution < 1) return
            allocate(parameters(resolution))
            call fill_uniform(knot_start(knot,nc,degree),knot_end(knot,nc,degree),parameters)
        else if (allocated(cached)) then
            parameters = cached
        else
            return
        end if
        valid = size(parameters) > 0 .and. all(ieee_is_finite(parameters))
    end subroutine query_parameter_axis

    !> Euclidean norm without squaring the original scale.
    pure function scaled_norm(vector) result(value)
        real(rk), intent(in) :: vector(:)
        real(rk) :: value, largest
        value = 0.0_rk
        if (size(vector) == 0) return
        largest = maxval(abs(vector))
        if (largest > 0.0_rk) value = largest*sqrt(sum((vector/largest)**2))
    end function scaled_norm

    !> Multiply short scale factors without intermediate overflow or underflow.
    pure function scaled_product(values) result(value)
        real(rk), intent(in) :: values(:)
        real(rk), parameter :: lower = sqrt(sqrt(sqrt(tiny(1.0_rk))))
        real(rk), parameter :: upper = sqrt(sqrt(sqrt(0.5_rk*huge(1.0_rk))))
        real(rk) :: value
        ! Up to eight factors inside these bounds need no exponent scaling.
        if (size(values)<=8 .and. all(abs(values)>=lower .and. abs(values)<=upper)) then
            value = product(values)
        else if (all(ieee_is_finite(values))) then
            value = scale(product(fraction(values)),sum(exponent(values)))
        else
            value = product(values)
        end if
    end function scaled_product

    !> Contract the local spline jet with the control net in native knot coordinates.
    pure subroutine geometry_jet(Xc,param_gradient,jacobian,ids,basis,point,param_hessian,curvature)
        real(rk), intent(in) :: Xc(:,:)
        real(rk), contiguous, intent(in) :: param_gradient(:,:)
        real(rk), intent(out) :: jacobian(:,:)
        integer, optional, intent(in) :: ids(:)
        real(rk), optional, intent(in) :: basis(:), param_hessian(:,:,:)
        real(rk), optional, intent(out) :: point(:), curvature(:,:,:)
        real(rk) :: control(size(param_gradient,1)),origin
        integer :: d,i,j

        do d=1,size(Xc,2)
            if (present(ids)) then
                control=Xc(ids,d)
            else
                control=Xc(:,d)
            end if
            if (present(point)) point(d)=dot_product(basis,control)
            ! Derivatives annihilate constants. Center same-sign controls before
            ! contraction; an interval crossing zero already has a safe origin.
            origin=max(0.0_rk,minval(control))+min(0.0_rk,maxval(control))
            control=control-origin
            do j=1,size(jacobian,2)
                ! Use each contiguous basis column, as for the geometric Hessian.
                jacobian(d,j)=dot_product(control,param_gradient(:,j))
                if (.not.present(param_hessian) .or. .not.present(curvature)) cycle
                do i=1,j
                    curvature(d,i,j)=dot_product(control,param_hessian(:,i,j))
                    curvature(d,j,i)=curvature(d,i,j)
                end do
            end do
        end do
    end subroutine geometry_jet

    !> Map one tensor-product analysis point using its active controls.
    !! Widths and the Gauss weight belong to the parent-to-knot map.
    !! Parametric and physical Hessians are supplied together when requested.
    pure subroutine map_analysis_point(Xc,ids,basis,param_gradient,widths,gauss_weight,strict,gradient,measure,valid, &
        point,jacobian,param_hessian,physical_hessian,dmeasure_dXc,dmeasure_dlogWc,dgradient_dXc,dgradient_dlogWc,djacobian_dlogWc)
        real(rk), intent(in) :: Xc(:,:), basis(:), widths(:), gauss_weight
        real(rk), contiguous, intent(in) :: param_gradient(:,:)
        integer, intent(in) :: ids(:)
        logical, intent(in) :: strict
        real(rk), intent(out) :: gradient(:,:), measure
        logical, intent(out) :: valid
        real(rk), intent(out), optional :: point(:), jacobian(:,:), physical_hessian(:,:,:)
        real(rk), intent(in), optional :: param_hessian(:,:,:)
        real(rk), intent(out), optional :: dmeasure_dXc(:,:), dmeasure_dlogWc(:)
        real(rk), intent(out), optional :: dgradient_dXc(:,:,:,:), dgradient_dlogWc(:,:,:)
        real(rk), intent(out), optional :: djacobian_dlogWc(:,:,:)
        real(rk) :: q(size(Xc,2),size(widths))
        real(rk) :: tangent(size(Xc,2),size(widths)), frame(size(Xc,2),size(widths))
        real(rk) :: curvature(size(Xc,2),size(widths),size(widths)), density
        logical :: inverted, design_gradient

        call geometry_jet(Xc,param_gradient,tangent,ids,basis,point,param_hessian,curvature)
        design_gradient=present(dgradient_dXc) .or. present(dgradient_dlogWc)
        if (design_gradient) then
            call jacobian_frame(tangent,frame,density,inverted,valid,q)
        else
            call jacobian_frame(tangent,frame,density,inverted,valid)
        end if
        valid=valid .and. (.not.strict .or. .not.inverted)
        if (.not.valid) return
        if (present(jacobian)) jacobian=tangent
        call map_basis_derivatives(param_gradient,frame,gradient,param_hessian,curvature,physical_hessian)
        measure=scaled_product([density,widths,gauss_weight])
        valid=ieee_is_finite(measure) .and. measure>0.0_rk .and. all(ieee_is_finite(gradient)) .and. all(ieee_is_finite(basis))
        if (present(point)) valid=valid .and. all(ieee_is_finite(point))
        if (present(physical_hessian)) valid=valid .and. all(ieee_is_finite(physical_hessian))
        if (valid .and. (present(dmeasure_dXc) .or. present(dmeasure_dlogWc))) then
            call measure_sensitivity(measure,size(widths),basis,gradient,Xc,ids,valid,dmeasure_dXc,dmeasure_dlogWc)
        end if
        if (valid .and. design_gradient) then
            call gradient_sensitivity(basis,gradient,q,Xc,ids,valid,dgradient_dXc,dgradient_dlogWc,param_gradient)
        end if
        if (valid .and. present(djacobian_dlogWc)) then
            call jacobian_log_sensitivity(basis,param_gradient,tangent,Xc,djacobian_dlogWc,valid,ids)
        end if
    end subroutine map_analysis_point

    !> Orthonormal tangents and the factors J = Q R diag(scales).
    pure subroutine tangent_qr(jacobian,frame,r,scales,valid)
        real(rk), intent(in) :: jacobian(:,:)
        real(rk), intent(out) :: frame(:,:),r(3,3),scales(3)
        logical, intent(out) :: valid
        integer :: ndim,ncoord,d,i,pass
        real(rk) :: coefficient,length_

        frame=0.0_rk
        valid = .false.
        ndim = size(jacobian,2)
        ncoord = size(jacobian,1)
        if (ndim < 1 .or. ndim > 3 .or. ncoord < ndim) return
        if (any(shape(frame) /= shape(jacobian)) .or. any(.not. ieee_is_finite(jacobian))) return
        r = 0.0_rk
        do d = 1, ndim
            scales(d) = maxval(abs(jacobian(:,d)))
            if (.not. ieee_is_finite(scales(d)) .or. scales(d) <= 0.0_rk) return
            frame(:,d) = jacobian(:,d)/scales(d)
            do pass = 1, 2
                do i = 1, d-1
                    coefficient = dot_product(frame(:,i),frame(:,d))
                    r(i,d) = r(i,d)+coefficient
                    frame(:,d) = frame(:,d)-coefficient*frame(:,i)
                end do
            end do
            length_ = scaled_norm(frame(:,d))
            if (length_ <= 64.0_rk*epsilon(1.0_rk)) return
            r(d,d) = length_
            frame(:,d) = frame(:,d)/length_
        end do
        valid=.true.
    end subroutine

    !> A positive quadrature measure, or zero when its magnitude cannot be represented.
    pure elemental real(rk) function scaled_measure(mantissa,power) result(measure)
        real(rk), intent(in) :: mantissa
        integer, intent(in) :: power
        measure=0.0_rk
        if (power+exponent(mantissa)<=maxexponent(measure)) measure=scale(mantissa,power)
    end function

    !> Combine the factored native weight with sqrt(det(J^T J)) before rounding.
    pure real(rk) function weighted_tangent_measure(jacobian,weight,power) result(measure)
        real(rk), intent(in) :: jacobian(:,:),weight
        integer, intent(in) :: power
        real(rk) :: frame(size(jacobian,1),size(jacobian,2)),r(3,3),scales(3),factors(7)
        integer :: n,d
        logical :: valid

        call tangent_qr(jacobian,frame,r,scales,valid)
        measure=0.0_rk
        if (.not.valid) return
        n=size(jacobian,2)
        factors(:2*n+1)=[weight,scales(:n),(r(d,d),d=1,n)]
        measure=scaled_measure(product(fraction(factors(:2*n+1))),power+sum(exponent(factors(:2*n+1))))
    end function

    !> Contravariant frame J*(J^T*J)^(-1), without forming the Gram matrix.
    !! Column-scaled, reorthogonalized QR for up to three parametric directions.
    pure subroutine jacobian_frame(jacobian, contravariant, measure, inverted, valid, q)
        real(rk), intent(in) :: jacobian(:,:)
        real(rk), intent(out) :: contravariant(:,:), measure
        logical, intent(out) :: inverted, valid
        real(rk), intent(out), optional :: q(:,:)  !! Orthonormal tangent frame from the same QR factorization.
        real(rk) :: r(3,3), scales(3), z(3), factors(6), coefficient
        integer :: ndim, ncoord, d, i

        ndim=size(jacobian,2)
        ncoord=size(jacobian,1)
        measure=0.0_rk
        inverted=.false.
        call tangent_qr(jacobian,contravariant,r,scales,valid)
        if (.not.valid) return
        if (present(q)) q = contravariant
        if (ncoord == ndim) then
            select case (ndim)
            case (1)
                inverted = contravariant(1,1) < 0.0_rk
            case (2)
                inverted = contravariant(1,1)*contravariant(2,2)-contravariant(2,1)*contravariant(1,2) < 0.0_rk
            case default
                inverted = det(contravariant) < 0.0_rk
            end select
        end if
        ! R^T z=e_d gives z(1:d-1)=0, so earlier Q columns may be overwritten.
        do d = 1, ndim
            z = 0.0_rk
            do i = 1, ndim
                coefficient = 0.0_rk
                if (i == d) coefficient = 1.0_rk
                if (i > 1) coefficient = coefficient-dot_product(r(1:i-1,i),z(1:i-1))
                z(i) = coefficient/r(i,i)
            end do
            do i = 1, ncoord
                contravariant(i,d) = dot_product(contravariant(i,d:ndim),z(d:ndim))/scales(d)
            end do
        end do
        factors(:ndim) = scales(:ndim)
        do d = 1, ndim
            factors(ndim+d) = r(d,d)
        end do
        measure = scaled_product(factors(:2*ndim))
        valid = all(ieee_is_finite(contravariant)) .and. ieee_is_finite(measure) .and. measure > 0.0_rk
    end subroutine jacobian_frame

    !> Map local parametric gradients and symmetric Hessians to physical coordinates.
    !! frame=J*(J^T*J)^(-1) is supplied by jacobian_frame. The Hessian inputs
    !! are required together with physical_hessian. Shapes must match, with
    !! one through three parametric directions; arbitrary embeddings are allowed.
    pure subroutine map_basis_derivatives(param_gradient, frame, physical_gradient, &
        param_hessian, geometry_hessian, physical_hessian)
        real(rk), intent(in) :: param_gradient(:,:), frame(:,:)
        real(rk), intent(out) :: physical_gradient(:,:)
        real(rk), intent(in), optional :: param_hessian(:,:,:), geometry_hessian(:,:,:)
        real(rk), intent(out), optional :: physical_hessian(:,:,:)
        real(rk), parameter :: safe = sqrt(sqrt(sqrt(tiny(1.0_rk))))
        real(rk) :: corrected(3,3), value, fixed_frame(3,3), mapped(3,3)
        integer :: a, d, e, i, j, ndim, ncoord
        logical :: bounded

        physical_gradient = matmul(param_gradient,transpose(frame))
        if (.not. present(physical_hessian)) return
        ndim = size(param_gradient,2)
        ncoord = size(frame,1)
        fixed_frame = 0.0_rk
        if (ncoord<=3) fixed_frame(:ncoord,:ndim) = frame
        corrected = 0.0_rk
        bounded = ncoord<=3 .and. all(frame==0.0_rk .or. (abs(frame)>=safe .and. abs(frame)<=1.0_rk/safe)) &
            .and. all(param_hessian==0.0_rk .or. (abs(param_hessian)>=safe .and. abs(param_hessian)<=1.0_rk/safe)) &
            .and. all(geometry_hessian==0.0_rk .or. (abs(geometry_hessian)>=safe .and. abs(geometry_hessian)<=1.0_rk/safe)) &
            .and. all(physical_gradient==0.0_rk .or. (abs(physical_gradient)>=safe .and. abs(physical_gradient)<=1.0_rk/safe))
        ! Reuse the output prefix for the parameter-space curvature correction.
        do j = 1, ndim
            do i = 1, ndim
                physical_hessian(:,i,j) = param_hessian(:,i,j)-matmul(physical_gradient,geometry_hessian(:,i,j))
            end do
        end do
        do a = 1, size(param_gradient,1)
            corrected(:ndim,:ndim) = physical_hessian(a,:ndim,:ndim)
            if (bounded) then
                mapped = matmul(fixed_frame,matmul(corrected,transpose(fixed_frame)))
                physical_hessian(a,:,:) = mapped(:ncoord,:ncoord)
                cycle
            end if
            do e = 1, size(frame,1)
                do d = 1, e
                    value = 0.0_rk
                    do j = 1, ndim
                        do i = 1, ndim
                            value = value+scaled_product([frame(d,i),corrected(i,j),frame(e,j)])
                        end do
                    end do
                    physical_hessian(a,d,e) = value
                    if (d /= e) physical_hessian(a,e,d) = value
                end do
            end do
        end do
    end subroutine map_basis_derivatives

    !> Relative controls P_a-X expressed as offset(a,d)*radix**power(d).
    !! Center at the largest basis coefficient to preserve translations and tiny weights.
    pure subroutine relative_controls(basis,Xc,ids,offset,power,residual,residual_power,param_gradient)
        real(rk), intent(in) :: basis(:), Xc(:,:)
        integer, intent(in), optional :: ids(:)
        real(rk), intent(out) :: offset(:,:)
        integer, intent(out) :: power(:)
        real(rk), intent(out), optional :: residual(:,:)
        integer, intent(out), optional :: residual_power(:)
        real(rk), intent(in), optional :: param_gradient(:,:)
        integer :: d, anchor, largest

        anchor=maxloc(basis,dim=1)
        do d=1,size(Xc,2)
            if (present(ids)) then
                offset(:,d)=Xc(ids,d)
            else
                offset(:,d)=Xc(:,d)
            end if
            largest=exponent(maxval(abs(offset(:,d))))
            power(d)=min(largest,0)+max(0,largest-maxexponent(1.0_rk)+2)
            offset(:,d)=scale(offset(:,d),-power(d))
            offset(:,d)=offset(:,d)-offset(anchor,d)
        end do
        if (present(residual)) then
            residual=offset
            residual_power=power
            if (present(param_gradient)) call remove_tangent_offsets(residual,residual_power,param_gradient)
            do d=1,size(Xc,2)
                residual(:,d)=residual(:,d)-dot_product(basis,residual(:,d))
            end do
        end if
        do d=1,size(Xc,2)
            offset(:,d)=offset(:,d)-dot_product(basis,offset(:,d))
        end do
    end subroutine relative_controls

    !> Subtract a tangent vector from every local control displacement without
    !! changing its normal projection. A pivot control removes the affine part
    !! before differentiating, so a long flat direction cancels before scaling.
    pure subroutine remove_tangent_offsets(offset,power,param_gradient)
        real(rk), intent(inout) :: offset(:,:)
        integer, intent(inout) :: power(:)
        real(rk), intent(in) :: param_gradient(:,:)
        real(rk) :: derivative(size(offset,1)), ratio(size(offset,1)), tangent(size(offset,2)), correction
        integer :: axis, d, pivot, row, exponent_

        do axis=1,size(param_gradient,2)
            do d=1,size(offset,2)
                exponent_=exponent(maxval(abs(offset(:,d))))
                offset(:,d)=scale(offset(:,d),-exponent_)
                power(d)=power(d)+exponent_
            end do
            derivative=scale(param_gradient(:,axis),-exponent(maxval(abs(param_gradient(:,axis)))))
            tangent=matmul(derivative,offset)
            pivot=maxloc(abs(tangent),dim=1)
            if (tangent(pivot)==0.0_rk) cycle
            row=maxloc(abs(offset(:,pivot)),dim=1)
            ratio=offset(:,pivot)/offset(row,pivot)
            do d=1,size(offset,2)
                if (d==pivot) cycle
                offset(:,d)=offset(:,d)-ratio*offset(row,d)
                correction=dot_product(derivative,offset(:,d))/tangent(pivot)
                offset(:,d)=offset(:,d)-correction*offset(:,pivot)
            end do
            offset(:,pivot)=0.0_rk
        end do
    end subroutine remove_tangent_offsets

    !> Return bias + sum(x*y*radix**shifts) as value*radix**power.
    pure subroutine design_dot(x,y,shifts,bias,wide,value,power)
        real(rk), intent(in) :: x(:), y(:), bias
        integer, intent(in) :: shifts(:)
        logical, intent(in) :: wide
        real(rk), intent(out) :: value
        integer, intent(out) :: power
        real(rk) :: terms(0:size(x))
        integer :: powers(0:size(x))
        power=0
        if (.not.wide) then
            value=bias+dot_product(x,y)
        else
            terms(0)=fraction(bias)
            powers(0)=exponent(bias)
            terms(1:)=fraction(x)*fraction(y)
            powers(1:)=exponent(x)+exponent(y)+shifts
            if (any(terms/=0.0_rk)) power=maxval(powers,mask=terms/=0.0_rk)
            value=sum(scale(terms,powers-power))
        end if
    end subroutine design_dot

    !> Combine two products with independent exponents before rounding to rk.
    pure elemental real(rk) function design_sum(a,b,p,c,d,q,wide) result(v)
        real(rk), intent(in) :: a,b,c,d
        integer, intent(in) :: p,q
        logical, intent(in) :: wide
        integer :: power
        if (wide) then
            call design_dot([a,c],[b,d],[p,q],0.0_rk,.true.,v,power)
            v=scale(v,power)
        else
            v=a*b+c*d
        end if
    end function design_sum

    !> Differentiate the geometry Jacobian at fixed knot parameters.
    !! Differentiating R_b*(P_b-X) in direction alpha gives the b-th log-weight column.
    pure subroutine jacobian_log_sensitivity(basis,param_gradient,jacobian,Xc,djacobian_dlogWc,valid,ids)
        real(rk), intent(in) :: basis(:), param_gradient(:,:), jacobian(:,:), Xc(:,:)
        real(rk), intent(out) :: djacobian_dlogWc(:,:,:)
        logical, intent(out) :: valid
        integer, optional, intent(in) :: ids(:)
        real(rk) :: offset(size(basis),size(Xc,2))
        integer :: b, alpha, power(size(Xc,2))
        integer, parameter :: safe_power=maxexponent(1.0_rk)/4
        logical :: wide

        call relative_controls(basis,Xc,ids,offset,power)
        wide=any(power/=0) .or. any(abs(exponent(offset))>safe_power) .or. &
            any(abs(exponent(basis))>safe_power) .or. any(abs(exponent(param_gradient))>safe_power) .or. &
            any(abs(exponent(jacobian))>safe_power)
        do b=1,size(basis)
            do alpha=1,size(jacobian,2)
                djacobian_dlogWc(:,alpha,b)=design_sum(offset(b,:),param_gradient(b,alpha),power, &
                    -basis(b),jacobian(:,alpha),0,wide)
            end do
        end do
        valid=all(ieee_is_finite(djacobian_dlogWc))
    end subroutine jacobian_log_sensitivity

    !> Derivatives of one weighted physical measure at fixed knot parameters.
    !! For local rational basis R and physical gradient g, they are m*g_a and
    !! m*((P_a-X) dot g_a - ndim*R_a) for a regular patch of any codimension.
    pure subroutine measure_sensitivity(measure,ndim,basis,gradient,Xc,ids,valid,dmeasure_dXc,dmeasure_dlogWc)
        real(rk), intent(in) :: measure, basis(:), gradient(:,:), Xc(:,:)
        integer, intent(in) :: ndim
        integer, optional, intent(in) :: ids(:)
        logical, intent(out) :: valid
        real(rk), optional, intent(out) :: dmeasure_dXc(:,:), dmeasure_dlogWc(:)
        real(rk) :: value
        real(rk), allocatable :: offset(:,:)
        integer :: a, power, axis_power(size(Xc,2))

        valid=ieee_is_finite(measure) .and. measure>0.0_rk
        if (.not.valid) return
        if (present(dmeasure_dXc)) then
            dmeasure_dXc=measure*gradient
            valid=all(ieee_is_finite(dmeasure_dXc))
        end if
        if (.not.valid .or. .not.present(dmeasure_dlogWc)) return
        allocate(offset(size(basis),size(Xc,2)))
        call relative_controls(basis,Xc,ids,offset,axis_power)
        do a=1,size(basis)
            call design_dot(gradient(a,:),offset(a,:),axis_power,-ndim*basis(a),.true.,value,power)
            dmeasure_dlogWc(a)=scale(fraction(measure)*fraction(value),exponent(measure)+exponent(value)+power)
        end do
        valid=all(ieee_is_finite(dmeasure_dlogWc))
    end subroutine measure_sensitivity

    !> Design derivatives of the ambient physical basis gradients at fixed parameters.
    !! Write N=I-QQ^T, h_ab=g_a dot g_b, and v_b=P_b-X. Then
    !! d g_a/d P_bj = h_ab N(:,j) - g_b g_aj, including changes of the tangent plane.
    !! The log-weight derivative contracts this tensor against v_b and adds
    !! (delta_ab-R_a)g_b from differentiating the rational basis.
    pure subroutine gradient_sensitivity(basis,gradient,q,Xc,ids,valid,dgradient_dXc,dgradient_dlogWc,param_gradient)
        real(rk), intent(in) :: basis(:), gradient(:,:), q(:,:), Xc(:,:)
        integer, intent(in), optional :: ids(:)
        logical, intent(out) :: valid
        real(rk), intent(out), optional :: dgradient_dXc(:,:,:,:), dgradient_dlogWc(:,:,:)
        real(rk), intent(in), optional :: param_gradient(:,:)
        integer, parameter :: safe_power=maxexponent(1.0_rk)/8
        real(rk) :: normal(size(q,1),size(q,1)), h(size(basis)), alpha(size(basis)), normal_offset(size(q,1))
        real(rk), allocatable :: offset(:,:), residual(:,:)
        integer :: hp(size(basis)), ap(size(basis)), np(size(q,1)), axis_power(size(q,1)), zeros(size(q,1))
        integer :: a,b,i,j,normal_power(size(q,1))
        logical :: embedded, wide, log_weights

        valid=.true.
        log_weights=present(dgradient_dlogWc)
        if (log_weights .and. size(basis)==size(q,2)+1) then
            ! On an affine simplex the physical basis is barycentric for every weight set.
            dgradient_dlogWc=0.0_rk
            log_weights=.false.
        end if
        if (.not.(present(dgradient_dXc) .or. log_weights)) return
        embedded=size(q,1)>size(q,2)
        normal=0.0_rk
        if (embedded) then
            normal=-matmul(q,transpose(q))
            do i=1,size(q,1)
                normal(i,i)=1.0_rk+normal(i,i)
            end do
            ! N is an orthogonal projector. Reproject nearly tangent coordinate axes to
            ! recover small normal components lost to cancellation in I-Q*Q^T.
            do i=1,size(q,1)
                if (normal(i,i)>=0.5_rk) cycle
                normal(:,i)=matmul(normal,normal(:,i))
                normal(i,:)=normal(:,i)
            end do
        end if
        wide=any(abs(exponent(gradient))>safe_power)
        if (log_weights) then
            allocate(offset(size(basis),size(q,1)))
            if (embedded) then
                allocate(residual(size(basis),size(q,1)))
                call relative_controls(basis,Xc,ids,offset,axis_power,residual,normal_power,param_gradient)
                wide=wide .or. any(abs(normal_power)>safe_power)
            else
                call relative_controls(basis,Xc,ids,offset,axis_power)
            end if
            wide=wide .or. any(abs(axis_power)>safe_power) .or. any(abs(exponent(offset))>safe_power)
            if (.not.wide) then
                do i=1,size(q,1)
                    offset(:,i)=scale(offset(:,i),axis_power(i))
                    if (embedded) residual(:,i)=scale(residual(:,i),normal_power(i))
                end do
            end if
        end if
        zeros=0
        h=0.0_rk
        hp=0
        normal_offset=0.0_rk
        np=0
        do b=1,size(basis)
            if (embedded) then
                do a=1,size(basis)
                    call design_dot(gradient(a,:),gradient(b,:),zeros,0.0_rk,wide,h(a),hp(a))
                end do
            end if
            if (present(dgradient_dXc)) then
                do j=1,size(q,1)
                    do i=1,size(q,1)
                        dgradient_dXc(:,i,b,j)=design_sum(h,normal(i,j),hp,-gradient(:,j),gradient(b,i),0,wide)
                    end do
                end do
            end if
            if (log_weights) then
                do a=1,size(basis)
                    call design_dot(gradient(a,:),offset(b,:),axis_power,basis(a)-merge(1.0_rk,0.0_rk,a==b),wide,alpha(a),ap(a))
                end do
                do i=1,size(q,1)
                    if (embedded) call design_dot(normal(i,:),residual(b,:),normal_power,0.0_rk,wide,normal_offset(i),np(i))
                    dgradient_dlogWc(:,i,b)=design_sum(h,normal_offset(i),hp+np(i),-alpha,gradient(b,i),ap,wide)
                end do
            end if
        end do
        if (present(dgradient_dXc)) valid=all(ieee_is_finite(dgradient_dXc))
        if (log_weights) valid=valid .and. all(ieee_is_finite(dgradient_dlogWc))
    end subroutine gradient_sensitivity

    ! Internal quadrature implementation shared by single-patch curves, surfaces and volumes.

    pure subroutine analysis_error(status,ndim,code,location,message)
        type(debug), intent(inout) :: status
        integer, intent(in) :: ndim, code
        character(*), intent(in) :: location
        character(*), optional, intent(in) :: message
        character(*), parameter :: categories(3)=[character(20) :: &
            "forcad_nurbs_curve","forcad_nurbs_surface","forcad_nurbs_volume"]
        character(128) :: detail
        detail="Cannot evaluate the requested quadrature data."
        if (present(message)) detail=message
        call status%set(code=code,severity=1,category=trim(categories(ndim)), &
            location=location,message=trim(detail), &
            suggestion="Check geometry initialization, element and Gauss indices, knot-span separation, rank and orientation.")
    end subroutine

    !> Caller-owned diagnostics start fresh; object diagnostics retain an earlier failure.
    pure subroutine quadrature_status(object_err,err,status)
        type(debug), target, intent(inout) :: object_err
        type(debug), target, optional, intent(out) :: err
        type(debug), pointer, intent(out) :: status
        status=>object_err
        if (.not.present(err)) return
        call err%reset()
        status=>err
    end subroutine

    !> Resolve the caller's diagnostic and map the selected Gauss points into their knot spans.
    pure subroutine prepare_rule(Xc,degree,nc,knot1,knot2,knot3,object_err,err,status,ie,ngauss,rule,location,ig)
        real(rk), allocatable, intent(in) :: Xc(:,:), knot1(:)
        real(rk), allocatable, intent(in), optional :: knot2(:), knot3(:)
        integer, intent(in) :: degree(:), nc(:), ie
        integer, optional, intent(in) :: ngauss(:), ig
        type(debug), target, intent(inout) :: object_err
        type(debug), target, optional, intent(out) :: err
        type(debug), pointer, intent(out) :: status
        type(element_rule), intent(out) :: rule
        character(*), intent(in) :: location
        real(rk), parameter :: unit(2)=[0.0_rk,1.0_rk]
        type(debug) :: span_status
        integer :: counts(size(nc)), ndim, first, last, q, code
        integer :: span
        logical :: initialized

        call quadrature_status(object_err,err,status)
        if (.not.status%ok) return
        ndim=size(nc)
        counts=degree+1
        if (present(ngauss)) counts=ngauss
        first=1
        last=product(counts)
        if (present(ig)) first=ig
        initialized=allocated(Xc) .and. allocated(knot1)
        if (initialized .and. ndim==1) initialized=size(Xc,2)>0
        if (present(knot2)) initialized=initialized .and. allocated(knot2)
        if (present(knot3)) initialized=initialized .and. allocated(knot3)
        code=0
        if (.not.initialized) then
            code=103
        else if (ndim==3 .and. size(Xc,2)/=3) then
            code=108
        else if (any(counts<1) .or. first<1 .or. first>last) then
            code=100
        else
            if (ndim==1) then
                span=active_span(knot1,nc(1),degree(1),ie)
                if (span>0) rule%bounds(:,1)=knot1(span:span+1)
                if (span==0) code=100
            else if (ndim==2) then
                call spline_element(nc,degree,knot1,ie,rule%nodes,rule%bounds(:,:ndim),span_status,knot2)
                if (.not.span_status%ok) code=100
            else
                call spline_element(nc,degree,knot1,ie,rule%nodes,rule%bounds(:,:ndim),span_status,knot2,knot3)
                if (.not.span_status%ok) code=100
            end if
        end if
        if (code/=0) then
            call analysis_error(status,ndim,code,location)
            return
        end if
        call gauss_leg(spread(unit,2,ndim),counts-1,rule%parameters,rule%weights,ig)
        do q=1,size(rule%weights)
            rule%parameters(q,:)=rule%bounds(1,:ndim)+rule%parameters(q,:)*(rule%bounds(2,:ndim)-rule%bounds(1,:ndim))
        end do
        if (present(ig) .and. ndim==2) then
            if (.not.in_open_span(rule,1,ndim)) then
                call analysis_error(status,ndim,107,location)
                return
            end if
        end if
        rule%nlocal=product(degree+1)
        rule%ncoord=size(Xc,2)
        rule%ndim=ndim
        rule%nq=size(rule%weights)
    end subroutine

    pure logical function in_open_span(rule,q,ndim)
        type(element_rule), intent(in) :: rule
        integer, intent(in) :: q, ndim
        in_open_span=all(rule%parameters(q,:)>rule%bounds(1,:ndim) .and. rule%parameters(q,:)<rule%bounds(2,:ndim))
    end function

    pure subroutine tensor_ansatz(Xc,Wc,rational,degree,nc,knot1,knot2,object_err,err,ie,ig,basis,gradient,measure, &
        ngauss,strict,Xg,dXg_dXt,d2Tgc_dXg2,dmeasure_dXc,dmeasure_dlogWc,dgradient_dXc,dgradient_dlogWc,knot3,curve_second, &
        djacobian_dXc,djacobian_dlogWc)
        real(rk), allocatable, intent(in) :: Xc(:,:), knot1(:)
        real(rk), allocatable, target, intent(inout) :: Wc(:)
        real(rk), allocatable, optional, target, intent(inout) :: knot2(:), knot3(:)
        integer, intent(in) :: degree(:), nc(:), ie, ig
        logical, intent(in) :: rational
        integer, optional, intent(in) :: ngauss(:)
        logical, optional, intent(in) :: strict
        type(debug), target, intent(inout) :: object_err
        type(debug), target, optional, intent(out) :: err
        real(rk), allocatable, target, intent(out) :: basis(:), gradient(:,:)
        real(rk), target, intent(out) :: measure
        real(rk), allocatable, target, optional, intent(out) :: Xg(:),dXg_dXt(:,:),d2Tgc_dXg2(:,:,:)
        real(rk), allocatable, target, optional, intent(out) :: dmeasure_dXc(:,:), dmeasure_dlogWc(:)
        real(rk), allocatable, target, optional, intent(out) :: dgradient_dXc(:,:,:,:), dgradient_dlogWc(:,:,:)
        real(rk), allocatable, target, optional, intent(out) :: djacobian_dXc(:,:), djacobian_dlogWc(:,:,:)
        real(rk), pointer, contiguous :: dJ_P(:,:,:), dJ_eta(:,:,:,:)
        real(rk), pointer, contiguous :: basis_q(:,:), grad_q(:,:,:), point(:,:), jac(:,:,:), hess(:,:,:,:)
        real(rk), pointer, contiguous :: dm_P(:,:,:), dm_eta(:,:), dg_P(:,:,:,:,:), dg_eta(:,:,:,:)
        real(rk), allocatable, target, optional, intent(out) :: curve_second(:)
        type(debug), pointer :: status
        type(element_rule) :: rule
        real(rk) :: measure_q(1)
        integer :: n,d,p,failure

        call prepare_rule(Xc,degree,nc,knot1,knot2,knot3,object_err,err,status,ie,ngauss,rule,"ansatz",ig)
        n=rule%nlocal
        d=rule%ncoord
        p=rule%ndim
        allocate(basis(n),gradient(n,d))
        if (present(curve_second)) allocate(curve_second(n))
        if (present(Xg)) allocate(Xg(d))
        if (present(dXg_dXt)) allocate(dXg_dXt(d,p))
        if (present(d2Tgc_dXg2)) allocate(d2Tgc_dXg2(n,d,d))
        if (present(dmeasure_dXc)) allocate(dmeasure_dXc(n,d))
        if (present(dmeasure_dlogWc)) allocate(dmeasure_dlogWc(n))
        if (present(dgradient_dXc)) allocate(dgradient_dXc(n,d,n,d))
        if (present(dgradient_dlogWc)) allocate(dgradient_dlogWc(n,d,n))
        if (present(djacobian_dXc)) allocate(djacobian_dXc(n,p))
        if (present(djacobian_dlogWc)) allocate(djacobian_dlogWc(d,p,n))
        measure=0.0_rk
        if (rule%nq==0) return

        ! One-column views evaluate just the selected Gauss point directly in caller storage.
        basis_q(1:n,1:1)=>basis
        grad_q(1:n,1:d,1:1)=>gradient
        nullify(point,jac,hess,dm_P,dm_eta,dg_P,dg_eta,dJ_P,dJ_eta)
        if (present(curve_second)) hess(1:n,1:1,1:1,1:1)=>curve_second
        if (present(Xg)) point(1:d,1:1)=>Xg
        if (present(dXg_dXt)) jac(1:d,1:p,1:1)=>dXg_dXt
        if (present(d2Tgc_dXg2)) hess(1:n,1:d,1:d,1:1)=>d2Tgc_dXg2
        if (present(dmeasure_dXc)) dm_P(1:n,1:d,1:1)=>dmeasure_dXc
        if (present(dmeasure_dlogWc)) dm_eta(1:n,1:1)=>dmeasure_dlogWc
        if (present(dgradient_dXc)) dg_P(1:n,1:d,1:n,1:d,1:1)=>dgradient_dXc
        if (present(dgradient_dlogWc)) dg_eta(1:n,1:d,1:n,1:1)=>dgradient_dlogWc
        if (present(djacobian_dXc)) dJ_P(1:n,1:p,1:1)=>djacobian_dXc
        if (present(djacobian_dlogWc)) dJ_eta(1:d,1:p,1:n,1:1)=>djacobian_dlogWc
        call evaluate_rule(Xc,Wc,rational,degree,nc,knot1,knot2,knot3,rule,strict,basis_q,grad_q,measure_q, &
            failure,point,jac,hess,dm_P,dm_eta,dg_P,dg_eta,dJ_P,dJ_eta)
        measure=measure_q(1)
        if (failure/=0) call analysis_error(status,p,failure,"ansatz")
    end subroutine

    pure subroutine tensor_ansatz_element(Xc,Wc,rational,degree,nc,knot1,knot2,object_err,err,ie,basis,gradient,measure, &
        ngauss,strict,Xg,dXg_dXt,d2Tgc_dXg2,dmeasure_dXc,dmeasure_dlogWc,dgradient_dXc,dgradient_dlogWc,knot3,curve_second, &
        djacobian_dXc,djacobian_dlogWc)
        real(rk), allocatable, intent(in) :: Xc(:,:), knot1(:)
        real(rk), allocatable, target, intent(inout) :: Wc(:)
        real(rk), allocatable, optional, target, intent(inout) :: knot2(:), knot3(:)
        integer, intent(in) :: degree(:), nc(:), ie
        logical, intent(in) :: rational
        integer, optional, intent(in) :: ngauss(:)
        logical, optional, intent(in) :: strict
        type(debug), target, intent(inout) :: object_err
        type(debug), target, optional, intent(out) :: err
        real(rk), allocatable, intent(inout) :: basis(:,:), gradient(:,:,:), measure(:)
        real(rk), allocatable, target, optional, intent(inout) :: Xg(:,:), dXg_dXt(:,:,:), d2Tgc_dXg2(:,:,:,:)
        real(rk), allocatable, target, optional, intent(inout) :: dmeasure_dXc(:,:,:), dmeasure_dlogWc(:,:)
        real(rk), allocatable, target, optional, intent(inout) :: dgradient_dXc(:,:,:,:,:), dgradient_dlogWc(:,:,:,:)
        real(rk), allocatable, target, optional, intent(inout) :: djacobian_dXc(:,:,:), djacobian_dlogWc(:,:,:,:)
        real(rk), pointer, contiguous :: dJ_P(:,:,:), dJ_eta(:,:,:,:)
        real(rk), pointer, contiguous :: point(:,:), jac(:,:,:), hess(:,:,:,:), dm_P(:,:,:), dm_eta(:,:)
        real(rk), pointer, contiguous :: dg_P(:,:,:,:,:), dg_eta(:,:,:,:)
        real(rk), allocatable, target, optional, intent(inout) :: curve_second(:,:)
        type(debug), pointer :: status
        type(element_rule) :: rule
        integer :: n,d,p,q,failure

        call prepare_rule(Xc,degree,nc,knot1,knot2,knot3,object_err,err,status,ie,ngauss,rule,"ansatz_element")
        n=rule%nlocal
        d=rule%ncoord
        p=rule%ndim
        q=rule%nq
        call resize_buffer(basis,[n,q])
        call resize_buffer(gradient,[n,d,q])
        call resize_buffer(measure,[q])
        if (present(curve_second)) call resize_buffer(curve_second,[n,q])
        if (present(Xg)) call resize_buffer(Xg,[d,q])
        if (present(dXg_dXt)) call resize_buffer(dXg_dXt,[d,p,q])
        if (present(d2Tgc_dXg2)) call resize_buffer(d2Tgc_dXg2,[n,d,d,q])
        if (present(dmeasure_dXc)) call resize_buffer(dmeasure_dXc,[n,d,q])
        if (present(dmeasure_dlogWc)) call resize_buffer(dmeasure_dlogWc,[n,q])
        if (present(dgradient_dXc)) call resize_buffer(dgradient_dXc,[n,d,n,d,q])
        if (present(dgradient_dlogWc)) call resize_buffer(dgradient_dlogWc,[n,d,n,q])
        if (present(djacobian_dXc)) call resize_buffer(djacobian_dXc,[n,p,q])
        if (present(djacobian_dlogWc)) call resize_buffer(djacobian_dlogWc,[d,p,n,q])
        if (q==0) return
        nullify(point,jac,hess,dm_P,dm_eta,dg_P,dg_eta,dJ_P,dJ_eta)
        if (present(curve_second)) hess(1:n,1:1,1:1,1:q)=>curve_second
        if (present(Xg)) point=>Xg
        if (present(dXg_dXt)) jac=>dXg_dXt
        if (present(d2Tgc_dXg2)) hess=>d2Tgc_dXg2
        if (present(dmeasure_dXc)) dm_P=>dmeasure_dXc
        if (present(dmeasure_dlogWc)) dm_eta=>dmeasure_dlogWc
        if (present(dgradient_dXc)) dg_P=>dgradient_dXc
        if (present(dgradient_dlogWc)) dg_eta=>dgradient_dlogWc
        if (present(djacobian_dXc)) dJ_P=>djacobian_dXc
        if (present(djacobian_dlogWc)) dJ_eta=>djacobian_dlogWc
        call evaluate_rule(Xc,Wc,rational,degree,nc,knot1,knot2,knot3,rule,strict,basis,gradient,measure, &
            failure,point,jac,hess,dm_P,dm_eta,dg_P,dg_eta,dJ_P,dJ_eta)
        if (failure/=0) call analysis_error(status,p,failure,"ansatz_element")
    end subroutine

    !> Reuse the same local tensor-jet storage at every requested Gauss point.
    pure subroutine evaluate_rule(Xc,Wc,rational,degree,nc,knot1,knot2,knot3,rule,strict,basis,gradient,measure, &
        failure,Xg,dXg_dXt,d2Tgc_dXg2,dmeasure_dXc,dmeasure_dlogWc,dgradient_dXc,dgradient_dlogWc,djacobian_dXc,djacobian_dlogWc)
        real(rk), intent(in) :: Xc(:,:)
        real(rk), allocatable, target, intent(inout) :: Wc(:)
        real(rk), contiguous, intent(in) :: knot1(:)
        real(rk), allocatable, optional, target, intent(inout) :: knot2(:), knot3(:)
        integer, intent(in) :: degree(:),nc(:)
        logical, intent(in) :: rational
        type(element_rule), intent(in) :: rule
        logical, optional, intent(in) :: strict
        real(rk), contiguous, intent(out) :: basis(:,:),gradient(:,:,:),measure(:)
        integer, intent(out) :: failure
        real(rk), contiguous, target, optional, intent(out) :: Xg(:,:),dXg_dXt(:,:,:),d2Tgc_dXg2(:,:,:,:)
        real(rk), contiguous, target, optional, intent(out) :: dmeasure_dXc(:,:,:),dmeasure_dlogWc(:,:)
        real(rk), contiguous, target, optional, intent(out) :: dgradient_dXc(:,:,:,:,:),dgradient_dlogWc(:,:,:,:)
        real(rk), contiguous, target, optional, intent(out) :: djacobian_dXc(:,:,:), djacobian_dlogWc(:,:,:,:)
        real(rk), pointer, contiguous :: param_gradient(:,:), dJ_eta(:,:,:)
        real(rk), pointer, contiguous :: point(:),jac(:,:),hess(:,:,:),dm_P(:,:),dm_eta(:),dg_P(:,:,:,:),dg_eta(:,:,:)
        real(rk), allocatable :: param_hessian(:,:,:)
        real(rk), allocatable, target :: work_gradient(:,:)
        real(rk) :: widths(rule%ndim)
        real(rk) :: ders1(0:2,0:degree(1)),ders2(0:2,0:merge(degree(min(2,size(nc))),0,size(nc)>=2))
        real(rk) :: ders3(0:2,0:merge(degree(size(nc)),0,size(nc)==3))
        integer :: first(size(nc)), q, order
        real(rk), contiguous, pointer :: active_weights(:), axis2(:), axis3(:)
        logical :: strict_,valid

        failure=0
        strict_=.false.
        if (present(strict)) strict_=strict
        widths=rule%bounds(2,:rule%ndim)-rule%bounds(1,:rule%ndim)
        order=merge(2,1,present(d2Tgc_dXg2))
        if (present(d2Tgc_dXg2)) allocate(param_hessian(rule%nlocal,rule%ndim,rule%ndim))
        nullify(point,jac,hess,dm_P,dm_eta,dg_P,dg_eta,dJ_eta)
        if (.not.present(djacobian_dXc)) then
            allocate(work_gradient(rule%nlocal,rule%ndim))
            param_gradient=>work_gradient
        end if
        nullify(active_weights,axis2,axis3)
        if (present(knot2)) axis2=>knot2
        if (present(knot3)) axis3=>knot3
        if (rational .and. allocated(Wc)) active_weights=>Wc
        do q=1,size(measure)
            if (.not.in_open_span(rule,q,rule%ndim)) then
                failure=107
                exit
            end if
            if (present(Xg)) point=>Xg(:,q)
            if (present(dXg_dXt)) jac=>dXg_dXt(:,:,q)
            if (present(d2Tgc_dXg2)) hess=>d2Tgc_dXg2(:,:,:,q)
            if (present(dmeasure_dXc)) dm_P=>dmeasure_dXc(:,:,q)
            if (present(dmeasure_dlogWc)) dm_eta=>dmeasure_dlogWc(:,q)
            if (present(dgradient_dXc)) dg_P=>dgradient_dXc(:,:,:,:,q)
            if (present(dgradient_dlogWc)) dg_eta=>dgradient_dlogWc(:,:,:,q)
            if (present(djacobian_dXc)) param_gradient=>djacobian_dXc(:,:,q)
            if (present(djacobian_dlogWc)) dJ_eta=>djacobian_dlogWc(:,:,:,q)
            call tensor_basis_jet(rule%parameters(q,:),knot1,axis2,axis3,degree,nc,order,first=first, &
                basis=basis(:,q),gradient=param_gradient,hessian=param_hessian, &
                ders1=ders1,ders2=ders2,ders3=ders3,Wc=active_weights)
            if (rule%ndim==1) then
                call curve_analysis_point(Xc(first(1):first(1)+degree(1),:),basis(:,q),param_gradient,widths(1), &
                    rule%weights(q),strict_,gradient(:,:,q),measure(q),valid,point,jac,param_hessian,hess, &
                    dm_P,dm_eta,dg_P,dg_eta,dJ_eta)
            else
                call map_analysis_point(Xc,rule%nodes,basis(:,q),param_gradient,widths,rule%weights(q),strict_, &
                    gradient(:,:,q),measure(q),valid,point,jac,param_hessian,hess,dm_P,dm_eta,dg_P,dg_eta,dJ_eta)
            end if
            if (.not.valid) then
                failure=108
                exit
            end if
        end do
        if (failure==0) return
        basis=0.0_rk
        gradient=0.0_rk
        measure=0.0_rk
        if (present(Xg)) Xg=0.0_rk
        if (present(dXg_dXt)) dXg_dXt=0.0_rk
        if (present(d2Tgc_dXg2)) d2Tgc_dXg2=0.0_rk
        if (present(dmeasure_dXc)) dmeasure_dXc=0.0_rk
        if (present(dmeasure_dlogWc)) dmeasure_dlogWc=0.0_rk
        if (present(dgradient_dXc)) dgradient_dXc=0.0_rk
        if (present(dgradient_dlogWc)) dgradient_dlogWc=0.0_rk
        if (present(djacobian_dXc)) djacobian_dXc=0.0_rk
        if (present(djacobian_dlogWc)) djacobian_dlogWc=0.0_rk
    end subroutine

    !> Evaluate one local tensor jet, computing only the requested derivative order.
    !! Global weights are normalized before any selected-control gathering.
    pure subroutine tensor_basis_jet(Xt,knot1,knot2,knot3,degree,nc,order,wrap, &
        first,basis,gradient,hessian,ders1,ders2,ders3,Wc,elem, side)
        real(rk), intent(in) :: Xt(:)
        integer, intent(in), optional :: side(:)
        real(rk), intent(in), contiguous :: knot1(:)
        real(rk), optional, intent(in), contiguous :: knot2(:), knot3(:)
        integer, intent(in) :: degree(:), nc(:), order
        logical, optional, intent(in) :: wrap(:)
        integer, intent(out) :: first(:)
        real(rk), intent(out), contiguous :: basis(0:)
        real(rk), optional, intent(out), contiguous :: gradient(0:,:)
        real(rk), intent(out), contiguous, optional :: hessian(0:,:,:)
        real(rk), intent(in), contiguous, optional :: Wc(:)
        real(rk), intent(out), contiguous :: ders1(0:,0:), ders2(0:,0:), ders3(0:,0:)
        integer, intent(in), optional :: elem(:)
        integer :: trace(3), start(3), counts(3), d
        logical :: wrapping(3)

        counts=1
        counts(:size(nc))=nc
        trace=0
        if (present(side)) trace(:size(nc))=side
        wrapping=.false.
        if (present(wrap)) then
            d=min(size(wrap),size(nc))
            wrapping(:d)=wrap(:d)
        end if
        start=1
        call basis_bspline_der_all_active(Xt(1),knot1,nc(1),degree(1),order,start(1),ders1,trace(1),wrapping(1))
        if (present(knot2)) then
            call basis_bspline_der_all_active(Xt(2),knot2,nc(2),degree(2),order,start(2),ders2,trace(2),wrapping(2))
        else
            ders2=0.0_rk
            ders2(0,0)=1.0_rk
        end if
        if (present(knot3)) then
            call basis_bspline_der_all_active(Xt(3),knot3,nc(3),degree(3),order,start(3),ders3,trace(3),wrapping(3))
        else
            ders3=0.0_rk
            ders3(0,0)=1.0_rk
        end if
        first=start(:size(nc))
        select case(order)
        case(0)
            call tensor_basis_derivative_local(start(1),counts(1),ders1,start(2),counts(2),ders2, &
                start(3),counts(3),ders3,basis,Wc,elem)
        case(1)
            call tensor_basis_derivatives2_local(start(1),counts(1),ders1,start(2),counts(2),ders2, &
                start(3),counts(3),ders3,basis,gradient,Wc=Wc)
        case default
            call tensor_basis_derivatives2_local(start(1),counts(1),ders1,start(2),counts(2),ders2, &
                start(3),counts(3),ders3,basis,gradient,hessian,Wc)
        end select
    end subroutine tensor_basis_jet

    !> Integrate the metric density using the same native-knot Jacobian as ansatz.
    !! Keep valid quadrature contributions if another point in the element fails.
    pure subroutine tensor_measure(Xc,degree,nc,knot1,knot2,Wc,rational,measure,object_err,err,ngauss,strict,knot3)
        real(rk), allocatable, intent(in) :: Xc(:,:), knot1(:), knot2(:)
        real(rk), contiguous, optional, intent(in) :: knot3(:)
        real(rk), allocatable, target, intent(inout) :: Wc(:)
        integer, intent(in) :: degree(:), nc(:)
        logical, intent(in) :: rational
        type(debug), target, intent(inout) :: object_err
        type(debug), target, optional, intent(out) :: err
        integer, optional, intent(in) :: ngauss(:)
        logical, optional, intent(in) :: strict
        real(rk), intent(out) :: measure
        type(debug), pointer :: status
        real(rk), pointer, contiguous :: rational_weights(:)
        integer, allocatable :: plan(:,:)
        real(rk), allocatable :: points(:,:), weights(:)
        real(rk) :: subtotal, largest
        integer :: counts(size(nc)), ne(size(nc)), e, code, failure, ndim
        logical :: oriented, ready
        character(10) :: location

        measure=0.0_rk
        call quadrature_status(object_err,err,status)
        if (.not.status%ok) return
        ndim=size(nc)
        location=merge("cmp_area  ","cmp_volume",ndim==2)
        ready=allocated(Xc) .and. allocated(knot1) .and. allocated(knot2) .and. (ndim==2 .or. present(knot3))
        if (ready) ready=valid_knot_vector(knot1,nc(1),degree(1)) .and. valid_knot_vector(knot2,nc(2),degree(2)) &
            .and. size(Xc,1)==product(nc) .and. size(Xc,2)>0
        if (ready .and. present(knot3)) ready=valid_knot_vector(knot3,nc(3),degree(3))
        if (allocated(Wc)) ready=ready .and. size(Wc)==product(nc) .and. valid_weights(Wc)
        counts=degree+1
        if (present(ngauss)) counts=ngauss
        if (.not.ready .or. (rational .and. .not.allocated(Wc))) then
            call analysis_error(status,ndim,103,trim(location),"Set consistent control points, knots, and weights before integration.")
            return
        else if (ndim==3 .and. size(Xc,2)/=3) then
            call analysis_error(status,ndim,108,trim(location),"Volume integration requires three physical coordinates.")
            return
        else if (any(counts<1)) then
            call analysis_error(status,ndim,100,trim(location),"Every quadrature count must be positive.")
            return
        end if
        nullify(rational_weights)
        if (rational) rational_weights=>Wc

        plan=spline_span_plan(nc,degree,knot1,knot2,knot3)
        ne=count(plan>0,dim=1)
        call gauss_leg(spread([0.0_rk,1.0_rk],2,size(nc)),counts-1,points,weights)
        oriented=.false.
        if (present(strict)) oriented=strict
        largest=0.0_rk
        failure=0
        !$omp simd private(subtotal,code) reduction(+:measure) reduction(max:largest,failure)
        do e=1,product(ne)
            call integrate_element(e,subtotal,code)
            measure=measure+subtotal
            largest=max(largest,subtotal)
            failure=max(failure,code)
        end do
        if (.not.ieee_is_finite(measure)) then
            measure=largest
            failure=max(failure,108)
        end if
        if (failure/=0) call analysis_error(status,ndim,failure,trim(location),"Invalid quadrature contributions were omitted.")
    contains
        pure subroutine integrate_element(element,total,errcode)
            integer, intent(in) :: element
            real(rk), intent(out) :: total
            integer, intent(out) :: errcode
            type(debug) :: status
            integer, allocatable :: ids(:)
            integer :: spans(size(nc)),first(size(nc)),index,d,q
            real(rk) :: bounds(2,size(nc)),widths(size(nc)),parameter(size(nc))
            real(rk) :: basis(product(degree+1)),gradient(product(degree+1),size(nc))
            real(rk) :: ders1(0:1,0:degree(1)),ders2(0:1,0:degree(2)),ders3(0:1,0:merge(degree(size(nc)),0,size(nc)==3))
            real(rk) :: jac(size(Xc,2),size(nc)),frame(size(Xc,2),size(nc)),density,contribution
            logical :: inverted,valid

            index=element-1
            do d=1,size(nc)
                spans(d)=plan(mod(index,ne(d))+1,d)
                index=index/ne(d)
            end do
            call spline_element(nc,degree,knot1,element,ids,bounds,status,knot2,knot3,spans)
            widths=bounds(2,:)-bounds(1,:)
            total=0.0_rk
            errcode=0
            do q=1,size(weights)
                parameter=bounds(1,:)+points(q,:)*widths
                if (.not.all(parameter>bounds(1,:) .and. parameter<bounds(2,:))) then
                    errcode=max(errcode,107)
                    cycle
                end if
                call tensor_basis_jet(parameter,knot1,knot2,knot3,degree,nc,1,first=first,basis=basis,gradient=gradient, &
                    ders1=ders1,ders2=ders2,ders3=ders3,Wc=rational_weights)
                call geometry_jet(Xc,gradient,jac,ids)
                call jacobian_frame(jac,frame,density,inverted,valid)
                contribution=scaled_product([density,widths,weights(q)])
                if (.not.valid .or. (oriented .and. inverted) .or. .not.ieee_is_finite(contribution) .or. &
                    contribution<=0.0_rk .or. total>huge(total)-contribution) then
                    errcode=108
                    cycle
                end if
                total=total+contribution
            end do
        end subroutine integrate_element
    end subroutine tensor_measure

    !> Unnormalized curve tangent and weighted length, shared by integration and IGA.
    pure subroutine curve_frame(controls,gradient,width,weight,strict,tangent,speed,measure,valid,hessian,curvature_vector)
        real(rk), intent(in) :: controls(:,:), width, weight
        real(rk), contiguous, intent(in) :: gradient(:,:)
        logical, intent(in) :: strict
        real(rk), target, intent(out) :: tangent(:)
        real(rk), intent(out) :: speed, measure
        logical, intent(out) :: valid
        real(rk), optional, intent(in) :: hessian(:,:,:)
        real(rk), optional, target, intent(out) :: curvature_vector(:)
        real(rk), pointer :: jacobian(:,:), curvature(:,:,:)
        real(rk) :: control_scale
        integer :: a

        jacobian(1:size(tangent),1:1)=>tangent
        nullify(curvature)
        if (present(curvature_vector)) curvature(1:size(tangent),1:1,1:1)=>curvature_vector
        call geometry_jet(controls,gradient,jacobian,param_hessian=hessian,curvature=curvature)
        speed = scaled_norm(tangent)
        valid = all(ieee_is_finite(tangent)) .and. ieee_is_finite(speed)
        measure = scaled_product([speed,width,weight])
        valid = valid .and. ieee_is_finite(measure) .and. .not.(speed>0.0_rk .and. measure<=0.0_rk)
        if (strict .and. valid) then
            control_scale = 0.0_rk
            do a = 2, size(controls,1)
                control_scale = max(control_scale,maxval(abs(controls(a,:)-controls(1,:))))
            end do
            valid = speed>64.0_rk*epsilon(1.0_rk)*control_scale/abs(width) .and. &
                .not. (size(tangent)==1 .and. tangent(1)<0.0_rk)
        end if
    end subroutine curve_frame

    !> Tangential physical gradient and optional arc-length Hessian for a curve.
    pure subroutine curve_analysis_point(controls,basis,gradient,width,weight,strict,physical,measure,valid, &
        point,jacobian,hessian,second,dmeasure_dXc,dmeasure_dlogWc,dgradient_dXc,dgradient_dlogWc,djacobian_dlogWc)
        real(rk), intent(in) :: controls(:,:), basis(:), width, weight
        real(rk), contiguous, intent(in) :: gradient(:,:)
        logical, intent(in) :: strict
        real(rk), intent(out) :: physical(:,:), measure
        logical, intent(out) :: valid
        real(rk), optional, intent(in) :: hessian(:,:,:)
        real(rk), optional, intent(out) :: point(:), jacobian(:,:), second(:,:,:), dmeasure_dXc(:,:), dmeasure_dlogWc(:)
        real(rk), optional, intent(out) :: dgradient_dXc(:,:,:,:), dgradient_dlogWc(:,:,:), djacobian_dlogWc(:,:,:)
        real(rk) :: tangent(size(physical,2),1), curvature_vector(size(physical,2))
        real(rk) :: speed, curvature
        integer :: d

        call curve_frame(controls,gradient,width,weight,strict,tangent(:,1),speed,measure,valid,hessian,curvature_vector)
        if (valid .and. present(djacobian_dlogWc)) then
            call jacobian_log_sensitivity(basis,gradient,tangent,controls,djacobian_dlogWc,valid)
        end if
        if (present(jacobian)) jacobian = tangent
        if (present(point)) then
            point=matmul(basis,controls)
            valid=valid .and. all(ieee_is_finite(point))
        end if
        physical = 0.0_rk
        if (present(second)) second = 0.0_rk
        if (valid .and. speed>0.0_rk) then
            tangent = tangent/speed
            do d = 1, size(tangent,1)
                physical(:,d) = (gradient(:,1)*tangent(d,1))/speed
            end do
            if (present(second)) then
                curvature = dot_product(tangent(:,1),curvature_vector)/speed
                second(:,1,1) = ((hessian(:,1,1)-gradient(:,1)*curvature)/speed)/speed
            end if
        end if
        valid = valid .and. all(ieee_is_finite(physical))
        if (present(second)) valid = valid .and. all(ieee_is_finite(second))
        if (valid .and. (present(dmeasure_dXc) .or. present(dmeasure_dlogWc))) then
            valid=speed>0.0_rk
            if (valid) call measure_sensitivity(measure,1,basis,physical,controls, &
                valid=valid,dmeasure_dXc=dmeasure_dXc,dmeasure_dlogWc=dmeasure_dlogWc)
        end if
        if (valid .and. (present(dgradient_dXc) .or. present(dgradient_dlogWc))) then
            valid=speed>0.0_rk
            if (valid) call gradient_sensitivity(basis,physical,tangent,controls, &
                valid=valid,dgradient_dXc=dgradient_dXc,dgradient_dlogWc=dgradient_dlogWc,param_gradient=gradient)
        end if
    end subroutine curve_analysis_point

    !> Convert homogeneous CSR rows to ordinary rational field coefficients.
    !! Exponent scaling avoids overflow in A_ij*w_j; A*w=0 certifies the denominator.
    pure subroutine normalize_field_constraints(rowptr, col, val, weights, valid)
        integer, intent(in) :: rowptr(:), col(:)
        real(rk), intent(inout) :: val(:)
        real(rk), intent(in) :: weights(:)
        logical, intent(out) :: valid
        integer :: row, lo, hi

        valid = .false.
        if (size(col) /= size(val) .or. size(rowptr) < 1) return
        if (rowptr(1) /= 1 .or. rowptr(size(rowptr)) /= size(val)+1) return
        if (any(rowptr(2:) < rowptr(:size(rowptr)-1))) return
        if (any(col < 1) .or. any(col > size(weights))) return
        if (any(.not. ieee_is_finite(weights)) .or. any(weights <= 0.0_rk)) return
        do row = 1, size(rowptr)-1
            lo = rowptr(row)
            hi = rowptr(row+1)-1
            call equilibrate_field_row(val(lo:hi),weights(col(lo:hi)),valid)
            if (.not. valid) return
        end do
        valid = .true.
    end subroutine normalize_field_constraints

    !> Form (A*w)/max(abs(A*w)) with a common radix exponent before checking partition of unity.
    pure subroutine equilibrate_field_row(val, weights, valid, shift)
        real(rk), intent(inout) :: val(:)
        real(rk), intent(in) :: weights(:)
        logical, intent(out) :: valid
        integer, optional, intent(in) :: shift(:)
        integer :: k, largest, power(size(val))
        real(rk) :: product_fraction, magnitude

        valid = all(ieee_is_finite(val))
        if (.not. valid .or. size(val)==0) return
        power = -huge(1)
        do k=1,size(val)
            if (val(k)==0.0_rk) cycle
            product_fraction = fraction(val(k))*fraction(weights(k))
            power(k) = exponent(val(k))+exponent(weights(k))+exponent(product_fraction)
            if (present(shift)) power(k) = power(k)+shift(k)
            val(k) = fraction(product_fraction)
        end do
        largest = maxval(power)
        if (largest == -huge(1)) return
        magnitude = maxval(abs(val),mask=power==largest)
        do k=1,size(val)
            if (val(k)==0.0_rk) cycle
            val(k) = scale(val(k)/magnitude,power(k)-largest)
            if (val(k)==0.0_rk) valid = .false.
        end do
        valid = valid .and. abs(sum(val)) <= sqrt(epsilon(1.0_rk))*sum(abs(val))
    end subroutine equilibrate_field_row

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Test exact structural nonzero status for sparse/refinement coefficients.
    !! Returns false only for positive or negative floating-point zero. Every
    !! finite nonzero value, infinity, and NaN is structural; no tolerance is
    !! applied. Keeping NaNs structural allows invalid arithmetic to propagate.
    pure elemental logical function structural_nonzero(x) result(nonzero)
        real(rk), intent(in) :: x

        nonzero = x < 0.0_rk .or. x > 0.0_rk .or. ieee_is_nan(x)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Solve an overdetermined full-rank least-squares problem by pivoted QR.
    !!
    !! For \(A\in\mathbb R^{m\times n}\), \(m\ge n\), and one or more
    !! right-hand sides, Householder transformations with column pivoting form
    !!
    !! \[
    !! A P = Q R,\qquad
    !! X=P R^{-1}Q_1^T B.
    !! \]
    !!
    !! Columns are pivoted by recomputed trailing two-norms. Recomputing rather
    !! than downdating the norms costs additional work but avoids accumulated
    !! downdate error in this fallback path. The matrix and right-hand side are
    !! independently multiplied by exact radix powers before factorization;
    !! the solution scale is restored afterward. A diagonal of \(R\) not
    !! exceeding
    !! `64*max(m,n)*epsilon(rk)*max(initial column norm)` is classified as
    !! numerically rank deficient. Invalid shapes, nonfinite inputs, rank
    !! deficiency, or a nonfinite result return `X(0,0)`.
    !!
    !! This routine is intended as the accuracy fallback after a structured
    !! normal-equation solve fails. Its storage is \(O(mn+m\,n_{rhs})\), and
    !! it deliberately returns no arbitrary solution for a rank-deficient
    !! design.
    !!
    !! **References:** LAPACK `DGELS` uses QR for overdetermined full-rank
    !! least squares; `DGEQP3` computes QR with column pivoting; `DGELSY`
    !! extends that factorization to rank-deficient minimum-norm problems.
    !! [DGELS](https://netlib.org/lapack/explore-html/d8/d83/group__gels_gaa65298f8ef218a625e40d0da3c95803c.html),
    !! [DGEQP3](https://netlib.org/lapack/explore-html/d0/dea/group__geqp3_gae96659507d789266660af28c617ef87f.html),
    !! [DGELSY](https://netlib.org/lapack/explore-html/dc/d8b/group__gelsy_ga6d1d46ead18df76e993cd4eda6dc1bbb.html).
    pure function solve_least_squares_qr(A, B) result(X)
        real(rk), intent(in), contiguous :: A(:,:)
            !! Design matrix `[m,n]` with `m>=n>=1`.
        real(rk), intent(in) :: B(:,:)
            !! Right-hand sides `[m,nrhs]`.
        real(rk), allocatable :: X(:,:)
        real(rk), allocatable :: Awork(:,:), Bwork(:,:), solution(:,:)
        real(rk), allocatable :: column_norm(:), householder(:)
        integer, allocatable :: permutation(:)
        real(rk) :: matrix_scale, rhs_scale, first_norm, rank_tol
        real(rk) :: column_value, vector_norm, alpha, projection, factor, swap_value
        integer :: m, n, nrhs, i, j, k, r, pivot
        integer :: matrix_exponent, rhs_exponent, recovery_exponent, swap_index

        m = size(A,1)
        n = size(A,2)
        nrhs = size(B,2)
        if (m < n .or. n < 1 .or. nrhs < 1 .or. size(B,1) /= m .or. &
            .not. all(ieee_is_finite(A)) .or. .not. all(ieee_is_finite(B))) then
            allocate(X(0,0))
            return
        end if

        matrix_scale = maxval(abs(A))
        if (.not. ieee_is_finite(matrix_scale) .or. matrix_scale <= 0.0_rk) then
            allocate(X(0,0))
            return
        end if
        matrix_exponent = exponent(matrix_scale)
        rhs_scale = maxval(abs(B))
        rhs_exponent = 0
        if (rhs_scale > 0.0_rk) rhs_exponent = exponent(rhs_scale)

        allocate(Awork(m,n), Bwork(m,nrhs))
        Awork = scale(A,-matrix_exponent)
        Bwork = scale(B,-rhs_exponent)
        allocate(column_norm(n), householder(m), permutation(n))
        do j = 1, n
            column_value = 0.0_rk
            do i = 1, m
                column_value = column_value + Awork(i,j)*Awork(i,j)
            end do
            column_norm(j) = sqrt(column_value)
            permutation(j) = j
        end do
        first_norm = maxval(column_norm)
        rank_tol = 64.0_rk*real(max(m,n),rk)*epsilon(1.0_rk)*first_norm

        do k = 1, n
            pivot = k
            do j = k + 1, n
                if (column_norm(j) > column_norm(pivot)) pivot = j
            end do
            if (.not. ieee_is_finite(column_norm(pivot)) .or. &
                column_norm(pivot) <= rank_tol) then
                allocate(X(0,0))
                return
            end if
            if (pivot /= k) then
                do i = 1, m
                    swap_value = Awork(i,k)
                    Awork(i,k) = Awork(i,pivot)
                    Awork(i,pivot) = swap_value
                end do
                swap_value = column_norm(k)
                column_norm(k) = column_norm(pivot)
                column_norm(pivot) = swap_value
                swap_index = permutation(k)
                permutation(k) = permutation(pivot)
                permutation(pivot) = swap_index
            end if

            vector_norm = 0.0_rk
            do i = k, m
                vector_norm = vector_norm + Awork(i,k)*Awork(i,k)
            end do
            vector_norm = sqrt(vector_norm)
            if (.not. ieee_is_finite(vector_norm) .or. vector_norm <= rank_tol) then
                allocate(X(0,0))
                return
            end if
            alpha = -sign(vector_norm,Awork(k,k))
            householder(k:m) = Awork(k:m,k)
            householder(k) = householder(k) - alpha
            factor = 0.0_rk
            do i = k, m
                factor = factor + householder(i)*householder(i)
            end do
            if (.not. ieee_is_finite(factor) .or. factor <= 0.0_rk) then
                allocate(X(0,0))
                return
            end if
            factor = 2.0_rk/factor

            do j = k, n
                projection = 0.0_rk
                do i = k, m
                    projection = projection + householder(i)*Awork(i,j)
                end do
                projection = factor*projection
                do i = k, m
                    Awork(i,j) = Awork(i,j) - projection*householder(i)
                end do
            end do
            do r = 1, nrhs
                projection = 0.0_rk
                do i = k, m
                    projection = projection + householder(i)*Bwork(i,r)
                end do
                projection = factor*projection
                do i = k, m
                    Bwork(i,r) = Bwork(i,r) - projection*householder(i)
                end do
            end do
            Awork(k,k) = alpha
            if (k < m) Awork(k+1:m,k) = 0.0_rk

            do j = k + 1, n
                column_value = 0.0_rk
                do i = k + 1, m
                    column_value = column_value + Awork(i,j)*Awork(i,j)
                end do
                column_norm(j) = sqrt(column_value)
            end do
        end do

        allocate(solution(n,nrhs), source=0.0_rk)
        do r = 1, nrhs
            do i = n, 1, -1
                projection = Bwork(i,r)
                do j = i + 1, n
                    projection = projection - Awork(i,j)*solution(j,r)
                end do
                if (abs(Awork(i,i)) <= rank_tol) then
                    allocate(X(0,0))
                    return
                end if
                solution(i,r) = projection/Awork(i,i)
            end do
        end do

        recovery_exponent = rhs_exponent - matrix_exponent
        allocate(X(n,nrhs), source=0.0_rk)
        do r = 1, nrhs
            do i = 1, n
                X(permutation(i),r) = scale(solution(i,r),recovery_exponent)
            end do
        end do
        if (.not. all(ieee_is_finite(X))) then
            deallocate(X)
            allocate(X(0,0))
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Fill an existing array with inclusive uniformly spaced values.
    !!
    !! A one-element output receives `a`; an empty output is unchanged.
    !!
    pure subroutine fill_uniform(a, b, x)
        real(rk), intent(in) :: a
            !! First value.
        real(rk), intent(in) :: b
            !! Last value when `size(x)>1`.
        real(rk), intent(out), contiguous :: x(:)
            !! Caller-owned contiguous output array.
        integer :: i, n

        n = size(x)
        if (n < 1) return
        if (n == 1) then
            x(1) = a
        else
            do concurrent (i = 1:n)
                x(i) = a + (b - a)*real(i - 1, rk)/real(n - 1, rk)
            end do
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> Elevate a verified periodic spline without replacing its cyclic topology.
    !!
    !! The periodic coefficients are extended over guarded periods, the resulting
    !! bounded workspace is elevated by Algorithm A5.9, and one interior period
    !! with \(p+t\) exterior knots is extracted. Because the extracted period is
    !! separated from the artificial clamped boundaries by at least the elevated
    !! support width, its coefficients are the cyclic degree-elevated
    !! representation. Repeated trailing rows are assigned from the leading rows
    !! exactly after the floating-point transformation.
    pure subroutine elevate_degree_periodic_A_5_9(t, knot, degree, Xcw, nc_new, knot_new, Xcw_new, Tmap, success)
        integer, intent(in) :: t
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: degree
        real(rk), intent(in), contiguous :: Xcw(:,:)
        integer, intent(out) :: nc_new
        real(rk), allocatable, intent(out) :: knot_new(:)
        real(rk), allocatable, intent(out) :: Xcw_new(:,:)
        real(rk), allocatable, intent(out), optional :: Tmap(:,:)
        logical, intent(out), optional :: success

        real(rk), allocatable :: knot_extended(:), knot_elevated(:)
        real(rk), allocatable :: Xcw_extended(:,:), Xcw_elevated(:,:), Textended(:,:)
        real(rk) :: a, b, period, tol
        integer :: nc, elevated_degree, unique_count, guard_periods, period_count
        integer :: nc_extended, infinite_index, base_index, shift
        integer :: i, j, first_knot, first_start, first_end, unique_new
        logical :: ok

        if (present(success)) success = .false.
        nc = size(Xcw,1)
        nc_new = 0
        if (t < 1 .or. .not. periodic_topology(knot, degree, Xcw)) then
            allocate(knot_new(0), Xcw_new(0,size(Xcw,2)))
            if (present(Tmap)) allocate(Tmap(0,0))
            return
        end if

        elevated_degree = degree + t
        unique_count = nc - degree
        guard_periods = (elevated_degree + unique_count - 1)/unique_count + 1
        period_count = 2*guard_periods + 1
        nc_extended = period_count*unique_count + degree
        a = knot(degree+1)
        b = knot(nc+1)
        period = b - a

        allocate(knot_extended(nc_extended+degree+1))
        allocate(Xcw_extended(nc_extended,size(Xcw,2)))
        do concurrent (i = 1:nc_extended)
            Xcw_extended(i,:) = Xcw(mod(i-1,unique_count)+1,:)
        end do
        do concurrent (i = 1:size(knot_extended)) local(infinite_index, base_index, shift)
            infinite_index = i - 1 - guard_periods*unique_count
            base_index = degree + modulo(infinite_index-degree, unique_count)
            shift = (infinite_index-base_index)/unique_count
            knot_extended(i) = knot(base_index+1) + real(shift,rk)*period
        end do

        if (present(Tmap)) then
            call elevate_degree_bounded_A_5_9(&
                t        = t,&
                knot     = knot_extended,&
                degree   = degree,&
                Xcw      = Xcw_extended,&
                nc_new   = nc_extended,&
                knot_new = knot_elevated,&
                Xcw_new  = Xcw_elevated,&
                Tmap     = Textended,&
                success  = ok)
        else
            call elevate_degree_bounded_A_5_9(&
                t        = t,&
                knot     = knot_extended,&
                degree   = degree,&
                Xcw      = Xcw_extended,&
                nc_new   = nc_extended,&
                knot_new = knot_elevated,&
                Xcw_new  = Xcw_elevated,&
                success  = ok)
        end if
        if (.not. ok) then
            allocate(knot_new(0), Xcw_new(0,size(Xcw,2)))
            if (present(Tmap)) allocate(Tmap(0,0))
            return
        end if

        tol = knot_tolerance(knot_elevated, 1, size(knot_elevated))
        first_start = 0
        first_end = 0
        do i = 1, size(knot_elevated)
            if (first_start == 0 .and. abs(knot_elevated(i)-a) <= tol) first_start = i
            if (first_end == 0 .and. abs(knot_elevated(i)-b) <= tol) first_end = i
        end do
        first_knot = first_start - elevated_degree
        unique_new = first_end - first_start
        nc_new = unique_new + elevated_degree
        if (first_start < 1 .or. first_end <= first_start .or. first_knot < 1 .or. &
            first_end+elevated_degree > size(knot_elevated) .or. &
            first_knot+nc_new-1 > size(Xcw_elevated,1)) then
            nc_new = 0
            allocate(knot_new(0), Xcw_new(0,size(Xcw,2)))
            if (present(Tmap)) allocate(Tmap(0,0))
            return
        end if

        knot_new = knot_elevated(first_knot:first_end+elevated_degree)
        Xcw_new = Xcw_elevated(first_knot:first_knot+nc_new-1,:)
        Xcw_new(unique_new+1:nc_new,:) = Xcw_new(1:elevated_degree,:)

        if (present(Tmap)) then
            allocate(Tmap(nc_new,nc), source=0.0_rk)
            do j = 1,size(Xcw_extended,1)
                i = mod(j-1,unique_count) + 1
                Tmap(:,i) = Tmap(:,i) + Textended(first_knot:first_knot+nc_new-1,j)
            end do
            Tmap(unique_new+1:nc_new,:) = Tmap(1:elevated_degree,:)
            Xcw_new = matmul(Tmap,Xcw)
        end if

        if (.not. valid_knot_vector(knot_new, nc_new, elevated_degree) .or. &
            .not. periodic_topology(knot_new, elevated_degree, Xcw_new)) then
            nc_new = 0
            deallocate(knot_new, Xcw_new)
            allocate(knot_new(0), Xcw_new(0,size(Xcw,2)))
            if (present(Tmap)) then
                deallocate(Tmap)
                allocate(Tmap(0,0))
            end if
            return
        end if

        if (present(success)) success = .true.
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Apply one affine transformation to point rows in place.
    !!
    !! For every point \(\mathbf x\), this routine evaluates
    !!
    !! \[
    !! \widetilde{\mathbf x}
    !! =\mathbf c+\mathbf A(\mathbf x-\mathbf c)+\mathbf b ,
    !! \]
    !!
    !! where `linear` is \(\mathbf A\), optional `center` is \(\mathbf c\),
    !! and optional `translation` is \(\mathbf b\). Missing optional vectors
    !! denote zero. Any finite square matrix is accepted, including reflection,
    !! shear, anisotropic scaling, and singular maps. A singular map remains an
    !! exact affine image but can reduce geometric rank and invalidate an
    !! analysis mapping.
    !!
    !! One-, two-, and three-coordinate points use specialized scalar kernels.
    !! Higher-dimensional points use one automatic row buffer, not a temporary
    !! point array. Invalid matrix or vector shapes and nonfinite transformation
    !! data leave `points` unchanged.
    pure subroutine transform_points(points, linear, translation, center)
        real(rk), intent(inout), contiguous :: points(:,:)
            !! Point rows `[npoint,ncoord]` to transform.
        real(rk), intent(in), contiguous :: linear(:,:)
            !! Square linear map `[ncoord,ncoord]`.
        real(rk), intent(in), contiguous, optional :: translation(:)
            !! Translation vector `[ncoord]`.
        real(rk), intent(in), contiguous, optional :: center(:)
            !! Fixed center `[ncoord]` about which `linear` acts.
        integer :: d, i, j, ncoord
        real(rk) :: offset(size(points,2)), point(size(points,2))
        real(rk) :: x, y, z

        ncoord = size(points, 2)
        if (ncoord < 1 .or. size(points,1) < 1) return
        if (size(linear,1) /= ncoord .or. size(linear,2) /= ncoord) return
        if (present(translation)) then
            if (size(translation) /= ncoord) return
        end if
        if (present(center)) then
            if (size(center) /= ncoord) return
        end if
        if (.not. all(ieee_is_finite(linear))) return
        if (present(translation)) then
            if (.not. all(ieee_is_finite(translation))) return
        end if
        if (present(center)) then
            if (.not. all(ieee_is_finite(center))) return
        end if

        offset = 0.0_rk
        if (present(translation)) offset = translation
        if (present(center)) then
            do d = 1, ncoord
                offset(d) = offset(d) + center(d)
                do j = 1, ncoord
                    offset(d) = offset(d) - linear(d,j)*center(j)
                end do
            end do
        end if

        select case (ncoord)
        case (1)
            do concurrent (i = 1:size(points,1)) local(x)
                x = points(i,1)
                points(i,1) = linear(1,1)*x + offset(1)
            end do
        case (2)
            do concurrent (i = 1:size(points,1)) local(x, y)
                x = points(i,1)
                y = points(i,2)
                points(i,1) = linear(1,1)*x + linear(1,2)*y + offset(1)
                points(i,2) = linear(2,1)*x + linear(2,2)*y + offset(2)
            end do
        case (3)
            do concurrent (i = 1:size(points,1)) local(x, y, z)
                x = points(i,1)
                y = points(i,2)
                z = points(i,3)
                points(i,1) = linear(1,1)*x + linear(1,2)*y + linear(1,3)*z + offset(1)
                points(i,2) = linear(2,1)*x + linear(2,2)*y + linear(2,3)*z + offset(2)
                points(i,3) = linear(3,1)*x + linear(3,2)*y + linear(3,3)*z + offset(3)
            end do
        case default
            do i = 1, size(points,1)
                point = points(i,:)
                do d = 1, ncoord
                    points(i,d) = offset(d)
                    do j = 1, ncoord
                        points(i,d) = points(i,d) + linear(d,j)*point(j)
                    end do
                end do
            end do
        end select
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Insert a knot while retaining a verified cyclic B-spline representation.
    !!
    !! This is repeated Boehm insertion with Dierckx's periodic boundary update.
    !! After each insertion, the \(p\) exterior knots and repeated coefficient
    !! rows on the affected side are regenerated from the opposite side of the
    !! period. The result therefore satisfies [[periodic_topology]] rather than
    !! merely preserving the mapping on the active interval. `T`, when requested,
    !! includes these cyclic row identifications and satisfies
    !! `Qw=matmul(T,Pw)`.
    !!
    !! The seam multiplicity must remain at most \(p\); a full \(C^{-1}\) break
    !! is not a periodic topology. Inputs and outputs follow the zero-based
    !! bounds of [[insert_knot_A_5_1]].
    !!
    !! **References:**
    !!
    !! - W. Boehm, "Inserting New Knots into B-Spline Curves,"
    !!   *Computer-Aided Design* 12 (1980), 199--201.
    !!   [doi:10.1016/0010-4485(80)90154-2](https://doi.org/10.1016/0010-4485(80)90154-2).
    !! - P. Dierckx, *Curve and Surface Fitting with Splines*, Oxford
    !!   University Press, 1993, periodic knot insertion.
    pure subroutine insert_knot_periodic_A_5_1(p, UP, Pw, u, k, s, r, nq, UQ, Qw, T)
        integer, intent(in) :: p
            !! Polynomial degree.
        real(rk), intent(in), contiguous :: UP(0:)
            !! Verified periodic knot vector.
        real(rk), intent(in), contiguous :: Pw(0:,:)
            !! Periodic homogeneous or polynomial control variables.
        real(rk), intent(in) :: u
            !! Knot value in the active interval.
        integer, intent(in) :: k
            !! Initial zero-based span containing `u`.
        integer, intent(in) :: s
            !! Initial multiplicity of `u`.
        integer, intent(in) :: r
            !! Nonnegative insertion count.
        integer, intent(out) :: nq
            !! New highest zero-based control index.
        real(rk), allocatable, intent(out) :: UQ(:)
            !! Refined periodic knot vector.
        real(rk), allocatable, intent(out) :: Qw(:,:)
            !! Refined periodic control variables.
        real(rk), allocatable, intent(out), optional :: T(:,:)
            !! Optional periodic refinement transformation.

        real(rk), allocatable :: knot_work(:), knot_step(:)
        real(rk), allocatable :: control_work(:,:), control_step(:,:)
        real(rk), allocatable :: Twork(:,:), Tstep(:,:), Tcomposed(:,:)
        real(rk) :: period
        integer :: insertion, span, multiplicity, nc_step, unique_count, i

        nq = -1
        if (p < 1 .or. r < 0 .or. .not. ieee_is_finite(u) .or. &
            .not. periodic_topology(UP, p, Pw)) then
            allocate(UQ(0:-1), Qw(0:-1,size(Pw,2)))
            if (present(T)) allocate(T(0:-1,0:-1))
            return
        end if
        if (k /= findspan(size(Pw,1)-1, p, u, UP) .or. &
            s /= compute_multiplicity(UP, u) .or. s+r > p) then
            allocate(UQ(0:-1), Qw(0:-1,size(Pw,2)))
            if (present(T)) allocate(T(0:-1,0:-1))
            return
        end if
        if (r == 0) then
            if (present(T)) then
                call insert_knot_A_5_1(&
                    p  = p,&
                    UP = UP,&
                    Pw = Pw,&
                    u  = u,&
                    k  = k,&
                    s  = s,&
                    r  = 0,&
                    nq = nq,&
                    UQ = UQ,&
                    Qw = Qw,&
                    T  = T)
            else
                call insert_knot_A_5_1(&
                    p  = p,&
                    UP = UP,&
                    Pw = Pw,&
                    u  = u,&
                    k  = k,&
                    s  = s,&
                    r  = 0,&
                    nq = nq,&
                    UQ = UQ,&
                    Qw = Qw)
            end if
            return
        end if

        knot_work = UP
        control_work = Pw
        if (present(T)) then
            allocate(Twork(size(Pw,1),size(Pw,1)), source=0.0_rk)
            do concurrent (i = 1:size(Pw,1))
                Twork(i,i) = 1.0_rk
            end do
        end if

        do insertion = 1, r
            span = findspan(size(control_work,1)-1, p, u, knot_work)
            multiplicity = compute_multiplicity(knot_work, u)
            if (present(T)) then
                call insert_knot_A_5_1(&
                    p  = p,&
                    UP = knot_work,&
                    Pw = control_work,&
                    u  = u,&
                    k  = span,&
                    s  = multiplicity,&
                    r  = 1,&
                    nq = nq,&
                    UQ = knot_step,&
                    Qw = control_step,&
                    T  = Tstep)
            else
                call insert_knot_A_5_1(&
                    p  = p,&
                    UP = knot_work,&
                    Pw = control_work,&
                    u  = u,&
                    k  = span,&
                    s  = multiplicity,&
                    r  = 1,&
                    nq = nq,&
                    UQ = knot_step,&
                    Qw = control_step)
            end if

            nc_step = size(control_step,1)
            unique_count = nc_step - p
            period = knot_step(nc_step) - knot_step(p)

            if (span >= nc_step-p-1) then
                knot_step(0:p-1) = knot_step(nc_step-p:nc_step-1) - period
                control_step(0:p-1,:) = control_step(unique_count:unique_count+p-1,:)
                if (present(T)) Tstep(0:p-1,:) = Tstep(unique_count:unique_count+p-1,:)
            end if
            if (span <= 2*p-1) then
                knot_step(size(knot_step)-p:size(knot_step)-1) = knot_step(p+1:2*p) + period
                control_step(unique_count:unique_count+p-1,:) = control_step(0:p-1,:)
                if (present(T)) Tstep(unique_count:unique_count+p-1,:) = Tstep(0:p-1,:)
            end if

            if (present(T)) then
                call sparse_left_matmul(&
                    A = Tstep,&
                    B = Twork,&
                    C = Tcomposed)
                call move_alloc(Tcomposed, Twork)
            end if
            call move_alloc(knot_step, knot_work)
            call move_alloc(control_step, control_work)
        end do

        if (.not. periodic_topology(knot_work, p, control_work)) then
            nq = -1
            allocate(UQ(0:-1), Qw(0:-1,size(Pw,2)))
            if (present(T)) allocate(T(0:-1,0:-1))
            return
        end if

        nq = size(control_work,1) - 1
        call move_alloc(knot_work, UQ)
        call move_alloc(control_work, Qw)
        if (present(T)) call move_alloc(Twork, T)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute a tensor-product basis, gradient, and Hessian on local support.
    !! The derivative tables contain rows zero through one, or through two when
    !! `hessian` is requested. Unused directions use a matching one-column
    !! identity table. Rational derivatives apply the quotient rule through the
    !! requested total order. This is a fused
    !! tensor-product specialization of the multivariate quotient recurrence,
    !! extending the rational-surface recurrence in Piegl and Tiller, *The
    !! NURBS Book*, Algorithm A4.4 to as many as three parametric directions.
    !! Direction 1 varies fastest in `basis`. `gradient(:,i)` stores
    !! \(\partial R_A/\partial\xi_i\), while `hessian(:,i,j)` stores the mixed
    !! derivative \(\partial^2R_A/(\partial\xi_i\partial\xi_j)\).
    !! When `Wc` is present, its active entries are first multiplied by one
    !! common exact radix power selected from their largest model exponent.
    !! This preserves the projective quotient and its derivatives.
    !!
    pure subroutine tensor_basis_derivatives2_local(first1,nc1,ders1,first2,nc2,ders2,first3,nc3,ders3, &
        basis,gradient,hessian,Wc)
        integer, intent(in) :: first1, first2, first3
            !! First active control index in each parameter direction.
        integer, intent(in) :: nc1, nc2, nc3
            !! Control counts in the same direction order.
        real(rk), contiguous, intent(in) :: ders1(0:,0:), ders2(0:,0:), ders3(0:,0:)
            !! Univariate basis jets indexed by derivative order and local basis index.
        real(rk), intent(out), contiguous :: basis(0:)
            !! Local tensor-product values.
        real(rk), intent(out), contiguous :: gradient(0:,:)
            !! Local parametric gradients, shape `[nlocal,rank]`.
        real(rk), intent(out), contiguous, optional :: hessian(0:,:,:)
            !! Optional local parametric Hessians, shape `[nlocal,rank,rank]`.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Optional flattened global weight array.
        real(rk) :: denominator, denominator1(3), denominator2, wi
        integer :: p1, p2, p3, l1, l2, l3, local, global, i, j, ndim, nlocal, weight_exponent, anchor

        basis = 0.0_rk
        gradient = 0.0_rk
        if (present(hessian)) hessian = 0.0_rk
        p1 = ubound(ders1,2)
        p2 = ubound(ders2,2)
        p3 = ubound(ders3,2)
        ndim = size(gradient,2)
        if (ubound(ders1,1) < 1 .or. ubound(ders2,1) < 1 .or. ubound(ders3,1) < 1) return
        if (present(hessian)) then
            if (ubound(ders1,1) < 2 .or. ubound(ders2,1) < 2 .or. ubound(ders3,1) < 2) return
        end if
        if (nc1 < 1 .or. nc2 < 1 .or. nc3 < 1) return
        if (first1 < 1 .or. first2 < 1 .or. first3 < 1) return
        if (first1 + p1 > nc1 .or. first2 + p2 > nc2 .or. first3 + p3 > nc3) return
        nlocal = (p1+1)*(p2+1)*(p3+1)
        if (size(basis) < nlocal) return
        if (ndim < 1 .or. ndim > 3 .or. size(gradient,1) < size(basis)) return
        if (present(hessian)) then
            if (size(hessian,1) < size(basis) .or. size(hessian,2) < ndim .or. size(hessian,3) < ndim) return
        end if

        if (present(Wc)) then
            if (size(Wc) < nc1*nc2*nc3) return
            global = first1 + nc1*((first2-1) + nc2*(first3-1))
            weight_exponent = exponent(Wc(global))
            do l3 = 0, p3
                do l2 = 0, p2
                    do l1 = 0, p1
                        global = first1 + l1 + nc1*((first2+l2-1) + nc2*(first3+l3-1))
                        weight_exponent = max(weight_exponent, exponent(Wc(global)))
                    end do
                end do
            end do
        end if

        ! Form the weighted polynomial jets once, then apply the quotient rule in place.
        do l3 = 0, p3
            do l2 = 0, p2
                do l1 = 0, p1
                    local = l1 + (p1+1)*(l2 + (p2+1)*l3)
                    wi = 1.0_rk
                    if (present(Wc)) then
                        global = first1 + l1 + nc1*((first2+l2-1) + nc2*(first3+l3-1))
                        wi = scale(Wc(global), -weight_exponent)
                    end if
                    basis(local) = wi*ders1(0,l1)*ders2(0,l2)*ders3(0,l3)
                    gradient(local,1) = wi*ders1(1,l1)*ders2(0,l2)*ders3(0,l3)
                    if (ndim >= 2) gradient(local,2) = wi*ders1(0,l1)*ders2(1,l2)*ders3(0,l3)
                    if (ndim == 3) gradient(local,3) = wi*ders1(0,l1)*ders2(0,l2)*ders3(1,l3)
                    if (present(hessian)) then
                        hessian(local,1,1) = wi*ders1(2,l1)*ders2(0,l2)*ders3(0,l3)
                        if (ndim >= 2) then
                            hessian(local,1,2) = wi*ders1(1,l1)*ders2(1,l2)*ders3(0,l3)
                            hessian(local,2,2) = wi*ders1(0,l1)*ders2(2,l2)*ders3(0,l3)
                        end if
                        if (ndim == 3) then
                            hessian(local,1,3) = wi*ders1(1,l1)*ders2(0,l2)*ders3(1,l3)
                            hessian(local,2,3) = wi*ders1(0,l1)*ders2(1,l2)*ders3(1,l3)
                            hessian(local,3,3) = wi*ders1(0,l1)*ders2(0,l2)*ders3(2,l3)
                        end if
                        do j = 2, ndim
                            do i = 1, j-1
                                hessian(local,j,i) = hessian(local,i,j)
                            end do
                        end do
                    end if
                end do
            end do
        end do
        if (.not. present(Wc)) return

        denominator = sum(basis(:nlocal-1))
        if (denominator <= 0.0_rk .or. .not. ieee_is_finite(denominator)) then
            basis = 0.0_rk
            gradient = 0.0_rk
            if (present(hessian)) hessian = 0.0_rk
            return
        end if
        basis(:nlocal-1) = basis(:nlocal-1)/denominator
        anchor=maxloc(basis(:nlocal-1),dim=1)-1
        do i = 1, ndim
            denominator1(i) = sum(gradient(:nlocal-1,i))
            gradient(:nlocal-1,i) = (gradient(:nlocal-1,i)-basis(:nlocal-1)*denominator1(i))/denominator
            if (basis(anchor)>0.5_rk) call balance_basis_derivative(gradient(:nlocal-1,i),anchor+1)
        end do
        if (present(hessian)) then
            do j = 1, ndim
                do i = 1, j
                    denominator2 = sum(hessian(:nlocal-1,i,j))
                    hessian(:nlocal-1,i,j) = (hessian(:nlocal-1,i,j)-basis(:nlocal-1)*denominator2 - &
                        gradient(:nlocal-1,i)*denominator1(j)-gradient(:nlocal-1,j)*denominator1(i))/denominator
                    if (basis(anchor)>0.5_rk) call balance_basis_derivative(hessian(:nlocal-1,i,j),anchor+1)
                    if (i /= j) hessian(:nlocal-1,j,i) = hessian(:nlocal-1,i,j)
                end do
            end do
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Find a disjoint-set representative while shortening the search path.
    !!
    !! `parent` must encode a valid one-based forest and `i` must identify one
    !! of its members. Path halving updates traversed parent links in place;
    !! `r` is the root satisfying `parent(r)==r`.
    pure subroutine disjoint_set_root(parent, i, r)
        integer, intent(inout), contiguous :: parent(:)
            !! Parent forest, modified by path halving.
        integer, intent(in) :: i
            !! One-based member whose representative is requested.
        integer, intent(out) :: r
            !! One-based representative of `i`.

        r = i
        do while (parent(r) /= r)
            parent(r) = parent(parent(r))
            r = parent(r)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Merge two disjoint-set components using path halving.
    !!
    !! Roots are linked deterministically by their integer identifier. This is
    !! the union-find structure analyzed by Tarjan (1975), without a separate
    !! rank array.
    !!
    pure subroutine disjoint_set_union(parent, a, b)
        integer, intent(inout), contiguous :: parent(:)
            !! Parent forest initialized with `parent(i)=i`.
        integer, intent(in) :: a
            !! One-based member of the first component.
        integer, intent(in) :: b
            !! One-based member of the second component.
        integer :: ra, rb

        call disjoint_set_root(parent, a, ra)
        call disjoint_set_root(parent, b, rb)
        if (ra == rb) return
        if (ra < rb) then
            parent(rb) = ra
        else
            parent(ra) = rb
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compress a disjoint-set forest and return compact component identifiers.
    !!
    !! Identifiers are one-based and assigned in first-member order.
    !!
    pure subroutine disjoint_set_map(parent, map, active, include_aliases)
        integer, intent(inout), contiguous :: parent(:)
            !! Parent forest; paths are shortened in place.
        integer, allocatable, intent(out) :: map(:)
            !! Component identifier for every input member.
        logical, optional, intent(in) :: active(:)
            !! Number components in first active-member order. Entirely inactive components receive zero.
        logical, optional, intent(in) :: include_aliases
            !! Also label inactive occurrences of participating components; used for geometry design variables.
        integer :: i, r, n

        allocate(map(size(parent)), source=0)
        n = 0
        do i = 1, size(parent)
            if (present(active)) then
                if (.not.active(i)) cycle
            end if
            call disjoint_set_root(parent, i, r)
            if (map(r) == 0) then
                n = n + 1
                map(r) = n
            end if
            map(i) = map(r)
        end do
        if (.not.present(active)) return
        if (present(include_aliases)) then
            if (include_aliases) then
                do i=1,size(parent)
                    if (active(i)) cycle
                    call disjoint_set_root(parent,i,r)
                    map(i)=map(r)
                end do
                return
            end if
        end if
        ! An inactive representative may have served as temporary storage for its active members.
        where (.not.active) map = 0
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the selected endpoint index of a one-dimensional control array.
    !!
    !! Side 1 is the minimum endpoint and side 2 the maximum endpoint. Because
    !! an endpoint contains one index, `reverse` has no effect. An invalid side
    !! returns the sentinel array `[0]`.
    pure subroutine boundary_index_1d(nc, side, reverse, idx)
        integer, intent(in) :: nc
            !! Number of controls.
        integer, intent(in) :: side
            !! Endpoint identifier, 1 for minimum or 2 for maximum.
        logical, intent(in) :: reverse
            !! Orientation flag; immaterial for a single endpoint.
        integer, allocatable, intent(out) :: idx(:)
            !! One-element endpoint index, or `[0]` for an invalid side.

        call boundary_layer_index_1d(nc, side, 0, reverse, idx)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return direction-1-fastest control indices on one surface edge.
    !!
    !! Sides `1:4` denote `u_min`, `u_max`, `v_min`, and `v_max`. The output is
    !! ordered along the tangential direction and `reverse` reverses that
    !! ordering. An invalid side returns an empty array.
    pure subroutine boundary_index_2d(nc, side, reverse, idx)
        integer, intent(in), contiguous :: nc(:)
            !! Directional control counts `[nu,nv]`.
        integer, intent(in) :: side
            !! Surface-edge identifier in `1:4`.
        logical, intent(in) :: reverse
            !! Reverse tangential output ordering when true.
        integer, allocatable, intent(out) :: idx(:)
            !! Flattened one-based edge indices.

        call boundary_layer_index_2d(nc, side, 0, reverse, idx)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return direction-1-fastest control indices on one volume face.
    !!
    !! Sides `1:6` denote the minimum and maximum faces normal to directions 1,
    !! 2, and 3. `swap` exchanges the two face coordinates, `flip` reverses
    !! either mapped coordinate, and `reverse` reverses the final flattened
    !! ordering. An invalid side returns an empty array.
    pure subroutine boundary_index_3d(nc, side, reverse, swap, flip, idx)
        integer, intent(in), contiguous :: nc(:)
            !! Directional control counts `[nu,nv,nw]`.
        integer, intent(in) :: side
            !! Volume-face identifier in `1:6`.
        logical, intent(in) :: reverse
            !! Reverse the final face ordering when true.
        logical, intent(in) :: swap
            !! Exchange the two tangential face directions when true.
        logical, intent(in), contiguous :: flip(:)
            !! Tangential-axis reversal flags; up to the first two entries apply.
        integer, allocatable, intent(out) :: idx(:)
            !! Flattened one-based face indices.

        call boundary_layer_index_3d(nc, side, 0, reverse, swap, flip, idx)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the endpoint-normal control index at a specified inward layer.
    !!
    !! `layer=0` selects the endpoint itself. Side 1 advances from the minimum
    !! endpoint and side 2 from the maximum endpoint. This low-level routine
    !! assumes `0<=layer<nc`; an invalid side returns `[0]`.
    pure subroutine boundary_layer_index_1d(nc, side, layer, reverse, idx)
        integer, intent(in) :: nc
            !! Number of controls.
        integer, intent(in) :: side
            !! Endpoint identifier, 1 for minimum or 2 for maximum.
        integer, intent(in) :: layer
            !! Zero-based inward layer.
        logical, intent(in) :: reverse
            !! Orientation flag; immaterial for a one-entry layer.
        integer, allocatable, intent(out) :: idx(:)
            !! One-element layer index, or `[0]` for an invalid side.

        allocate(idx(1))
        select case(side)
        case(1)
            idx(1) = 1 + layer
        case(2)
            idx(1) = nc - layer
        case default
            idx(1) = 0
        end select
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return surface-control indices in an inward layer parallel to one edge.
    !!
    !! Sides and ordering follow [[boundary_index]]. `layer=0` is the boundary;
    !! positive values advance along its inward normal control direction. This
    !! low-level routine assumes the layer lies inside the control net.
    pure subroutine boundary_layer_index_2d(nc, side, layer, reverse, idx)
        integer, intent(in), contiguous :: nc(:)
            !! Directional control counts `[nu,nv]`.
        integer, intent(in) :: side
            !! Surface-edge identifier in `1:4`.
        integer, intent(in) :: layer
            !! Zero-based inward layer.
        logical, intent(in) :: reverse
            !! Reverse tangential output ordering when true.
        integer, allocatable, intent(out) :: idx(:)
            !! Flattened one-based layer indices; empty for an invalid side.

        call control_layer(nc, side, layer, reverse, idx)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return volume-control indices in an inward layer parallel to one face.
    !!
    !! Face identifiers and orientation options follow [[boundary_index]].
    !! `layer=0` is the boundary face; positive layers advance along the inward
    !! normal control direction. The caller must supply an in-range layer.
    pure subroutine boundary_layer_index_3d(nc, side, layer, reverse, swap, flip, idx)
        integer, intent(in), contiguous :: nc(:)
            !! Directional control counts `[nu,nv,nw]`.
        integer, intent(in) :: side
            !! Volume-face identifier in `1:6`.
        integer, intent(in) :: layer
            !! Zero-based inward layer.
        logical, intent(in) :: reverse
            !! Reverse the final face ordering when true.
        logical, intent(in) :: swap
            !! Exchange the two tangential directions when true.
        logical, intent(in), contiguous :: flip(:)
            !! Tangential-axis reversal flags; up to the first two entries apply.
        integer, allocatable, intent(out) :: idx(:)
            !! Flattened one-based layer indices; empty for an invalid side.

        call control_layer(nc, side, layer, reverse, idx, swap, flip)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> Signed tangential storage axes, followed by the normal. An edge uses the first tangent.
    pure function boundary_axes(side,swap,flip,reverse) result(axes)
        integer, intent(in) :: side
        logical, optional, intent(in) :: swap,flip(:),reverse
        integer, parameter :: directions(3)=[1,2,3]
        integer :: axes(3),d
        axes=0
        if (side<1 .or. side>6) return
        axes(3)=(side+1)/2
        axes(:2)=pack(directions,directions/=axes(3))
        if (present(flip)) then
            do d=1,min(size(flip),2)
                if (flip(d)) axes(d)=-axes(d)
            end do
        end if
        if (present(swap)) then
            if (swap) axes(:2)=axes([2,1])
        end if
        if (present(reverse)) then
            if (reverse) axes(:2)=-axes(:2)
        end if
    end function

    !> Enumerate an oriented edge or face at a fixed inward control layer.
    pure subroutine control_layer(nc, side, layer, reverse, idx, swap, flip)
        integer, intent(in) :: nc(:), side, layer
        logical, intent(in) :: reverse
        integer, allocatable, intent(out) :: idx(:)
        logical, optional, intent(in) :: swap, flip(:)
        integer :: ndim, normal, axes(3), extent(2), stride(2), origin, depth, i, j

        ndim=merge(3,2,present(swap))
        if (side<1 .or. side>2*ndim) then
            allocate(idx(0))
            return
        end if
        axes=boundary_axes(side,swap,flip,reverse)
        normal=axes(3)
        extent=1
        extent(:ndim-1)=nc(abs(axes(:ndim-1)))
        stride=0
        call oriented_strides(nc,axes(:ndim-1),origin,stride(:ndim-1))
        depth=merge(layer,nc(normal)-1-layer,mod(side,2)==1)
        origin=origin+depth*product(nc(:normal-1))
        allocate(idx(product(extent)))
        do j=0,extent(2)-1
            do i=0,extent(1)-1
                idx(1+i+j*extent(1))=origin+i*stride(1)+j*stride(2)
            end do
        end do
    end subroutine control_layer


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the inward layer interpolated by an active endpoint.
    !!
    !! The result is the zero-based normal layer whose B-spline basis value is
    !! exactly one at the selected active endpoint, provided every other active
    !! value is exactly zero. A negative result means that the endpoint trace is
    !! a linear combination of multiple layers and therefore cannot be encoded
    !! by direct control-variable identification. Odd side identifiers select
    !! the minimum endpoint and even identifiers select the maximum endpoint.
    !!
    !! The result can be positive for valid knot vectors containing exterior
    !! knots. Open-clamped vectors return layer zero. See L. Piegl and W. Tiller,
    !! *The NURBS Book*, 2nd edition, Sections 2.2 and 3.2.
    !! [doi:10.1007/978-3-642-59223-2](https://doi.org/10.1007/978-3-642-59223-2).
    pure function interpolatory_boundary_layer(knot, nc, degree, side) result(layer)
        real(rk), intent(in), contiguous :: knot(:)
            !! Valid knot vector.
        integer, intent(in) :: nc
            !! Number of basis functions.
        integer, intent(in) :: degree
            !! Polynomial degree.
        integer, intent(in) :: side
            !! Rank-specific side identifier in the range `1:6`.
        integer :: layer
        real(rk) :: B(0:max(0, degree)), endpoint
        integer :: first, i, index

        layer = -1
        if (degree < 0 .or. degree >= nc) return
        if (size(knot) < nc + degree + 1) return
        if (side < 1 .or. side > 6) return

        endpoint = merge(knot(degree+1), knot(nc+1), mod(side, 2) == 1)
        call basis_bspline_der_order_active(endpoint, knot, nc, degree, 0, first, B(0:degree))

        index = 0
        do i = 0, degree
            if (.not. structural_nonzero(B(i))) cycle
            if (.not. ieee_is_finite(B(i)) .or. abs(B(i) - 1.0_rk) > 0.0_rk .or. index /= 0) return
            index = first + i
        end do
        if (index < 1 .or. index > nc) return

        layer = merge(index - 1, nc - index, mod(side, 2) == 1)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute a scale- and span-aware tolerance for knot comparisons.
    !!
    !! The tolerance is based on endpoint spacing and interval width, then
    !! capped at one eighth of the smallest positive inspected span so that two
    !! distinct spans cannot collapse under the comparison.
    !!
    pure function knot_tolerance(knot, first, last) result(tol)
        real(rk), intent(in), contiguous :: knot(:)
            !! Knot vector; an empty vector returns zero.
        integer, intent(in) :: first
            !! First inspected one-based index, clamped to the array.
        integer, intent(in) :: last
            !! Last inspected one-based index; reversed bounds are accepted.
        real(rk) :: tol, min_span, span, width
        integer :: i, lo, hi, tmp

        tol = 0.0_rk
        if (size(knot) == 0) return

        lo = max(1, min(first, size(knot)))
        hi = max(1, min(last,  size(knot)))
        if (lo > hi) then
            tmp = lo
            lo = hi
            hi = tmp
        end if

        width = abs(knot(hi) - knot(lo))
        tol = 8.0_rk*max(&
            spacing(knot(lo)),&
            spacing(knot(hi)),&
            epsilon(1.0_rk)*width)

        min_span = huge(1.0_rk)
        do i = lo + 1, hi
            span = knot(i) - knot(i-1)
            if (span > 0.0_rk) min_span = min(min_span, span)
        end do
        if (min_span < huge(1.0_rk)) tol = min(tol, 0.125_rk*min_span)
    end function
    !===============================================================================


    !===============================================================================
    !> Accept the default trace or one selector in {-1,0,+1} per parameter.
    pure logical function valid_trace_vector(side,ndim) result(valid)
        integer, intent(in), optional :: side(:)
        integer, intent(in) :: ndim
        valid = .true.
        if (present(side)) valid = size(side)==ndim .and. all(side>=-1 .and. side<=1)
    end function

    pure logical function valid_trace_scalar(side) result(valid)
        integer, intent(in), optional :: side
        valid = .true.
        if (present(side)) valid = valid_trace_vector([side],1)
    end function

    !> Optionally wrap one spline parameter into its active interval.
    !! Disabled wrapping and invalid spline states leave the parameter unchanged.
    !! Wrapping maps into \([a,b)\), except that the left seam trace uses \(b\).
    !! It changes parameter evaluation only; cyclic geometry is not required.
    !!
    pure function map_parameter(Xt, knot, nc, degree, wrap_parameter, side) result(u)
        real(rk), intent(in) :: Xt
            !! Input parameter.
        real(rk), intent(in), contiguous :: knot(:)
            !! Knot vector.
        integer, intent(in) :: nc
            !! Number of control points.
        integer, intent(in) :: degree
            !! Polynomial degree.
        logical, intent(in) :: wrap_parameter
            !! Enable modulo mapping when true.
        integer, intent(in), optional :: side
            !! A negative side takes the left trace at the wrapped seam.
        real(rk) :: u
        real(rk) :: a, b, period
        logical :: left

        u = Xt
        if (.not. wrap_parameter .or. .not. ieee_is_finite(Xt)) return
        if (degree < 0 .or. nc < 1 .or. degree >= nc) return
        if (size(knot) < nc + degree + 1) return

        a = knot(degree+1)
        b = knot(nc+1)
        if (Xt > a .and. Xt < b) return
        left = .false.
        if (present(side)) then
            left = side < 0
            if (side/=0 .and. (Xt==a .or. Xt==b)) then
                if (b>a .and. ieee_is_finite(a) .and. ieee_is_finite(b)) u = merge(b,a,left)
                return
            end if
        end if
        if (Xt == a) return
        period = b - a
        if (.not. ieee_is_finite(period) .or. period <= 0.0_rk) return

        u = a + modulo(Xt - a, period)
        if (u >= b) u = a
        if (left .and. u == a) u = b
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Find the active 1-based span for a valid local spline space.
    !! The active parameter interval is [knot(degree+1), knot(nc+1)].
    pure subroutine active_span_at(Xt, knot, nc, degree, span, u, inside, side, wrap)
        real(rk), intent(in) :: Xt
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: nc, degree
        integer, intent(out) :: span
        real(rk), intent(out) :: u
        logical, intent(out) :: inside
        integer, intent(in), optional :: side
        logical, intent(in), optional :: wrap
        real(rk) :: a, b, endpoint

        span = 0
        u = Xt
        inside = .false.
        if (present(wrap)) u = map_parameter(Xt,knot,nc,degree,wrap,side)
        if (degree < 0 .or. nc < 1 .or. degree >= nc .or. .not. ieee_is_finite(u)) return
        if (size(knot) < nc + degree + 1) return

        a = knot(degree+1)
        b = knot(nc+1)
        if (b <= a) return
        if (u < a .or. u > b) then
            endpoint = max(a,min(b,u))
            if (abs(u-endpoint) > knot_tolerance(knot,degree+1,nc+1)) return
            u = endpoint
        end if
        span = findspan(nc-1,degree,u,knot,side)+1
        inside = .true.
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Validate the structural invariants of a one-dimensional spline space.
    !!
    !! A valid vector is finite and nondecreasing, has exactly
    !! `nc+degree+1` entries, has no run longer than `degree+1`, and has a
    !! positive active interval. Open, clamped, one-sided clamped, unclamped,
    !! uniform, nonuniform, and periodic-form vectors can all satisfy this
    !! predicate; endpoint multiplicity is not prescribed. Here
    !! "periodic-form" means only that exterior knots and unclamped active
    !! endpoints can be represented. This structural predicate does not verify
    !! periodic knot-spacing extension or repeated control data; use
    !! [[periodic_topology]] when those data are available.
    !!
    pure function valid_knot_vector(knot, nc, degree) result(valid)
        real(rk), intent(in), contiguous :: knot(:)
            !! Candidate knot vector.
        integer, intent(in) :: nc
            !! Number of control points.
        integer, intent(in) :: degree
            !! Candidate polynomial degree, satisfying `0<=degree<nc`.
        logical :: valid
        integer :: i, run

        valid = .false.
        if (degree < 0 .or. nc < 1 .or. degree >= nc) return
        if (size(knot) /= nc + degree + 1) return
        if (.not. all(ieee_is_finite(knot))) return

        run = 1
        do i = 2, size(knot)
            if (knot(i) < knot(i-1)) return
            if (knot(i) == knot(i-1)) then
                run = run + 1
                if (run > degree + 1) return
            else
                run = 1
            end if
        end do
        if (knot(nc+1) <= knot(degree+1)) return

        valid = .true.
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Verify a cyclic tensor-product spline topology in one parameter direction.
    !!
    !! Let \(p\) be `degree`, \(n_d\) the selected directional control count,
    !! \(n_u=n_d-p\), \(a=U_p\), \(b=U_{n_d}\), and \(T=b-a\). A representation
    !! is periodic only if its ordinary knot-vector invariants hold, the seam
    !! multiplicity is at most \(p\), and
    !!
    !! \[
    !! U_{i+n_u}=U_i+T,\qquad i=0,\ldots,2p,
    !! \]
    !!
    !! while every control layer with directional index \(i+n_u\) repeats layer
    !! \(i\), for \(i=0,\ldots,p-1\). Supplied rational weights must repeat in
    !! the same way. These cyclic coefficient and knot identities imply equality
    !! of the two seam traces and their derivatives through order \(p-s\), where
    !! \(s\) is the seam-knot multiplicity. They are therefore stronger than
    !! parameter wrapping or endpoint coincidence alone.
    !!
    !! `Xc` uses direction-1-fastest tensor flattening. Omit `nc` and `dir` for
    !! a curve or packed directional control lines. Every coordinate component
    !! is compared with a tolerance proportional to its control-net diameter
    !! and to the ULPs of the two values. Weights use their global positive
    !! scale and their corresponding ULPs. No absolute unit-scale floor is
    !! introduced, so uniform scaling of coordinates or weights cannot change
    !! the topology classification.
    !!
    !! **Reference:** M. W. Blake et al., *NASA Geometry Data Exchange
    !! Specification for Computational Fluid Dynamics (NASA IGES)*, NASA
    !! Reference Publication 1338, 1994, Sections 4.7.2 and 4.8.2.
    !! [NASA publication](https://ntrs.nasa.gov/citations/19950005054).
    pure function periodic_topology(knot, degree, Xc, Wc, nc, dir) result(periodic)
        real(rk), intent(in), contiguous :: knot(:)
            !! Complete candidate knot vector.
        integer, intent(in) :: degree
            !! Polynomial degree in the selected direction.
        real(rk), intent(in), contiguous :: Xc(:,:)
            !! Curve controls, packed lines, or flattened tensor-product controls.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Optional flattened positive rational weights.
        integer, intent(in), contiguous, optional :: nc(:)
            !! Optional tensor-product control counts.
        integer, intent(in), optional :: dir
            !! Selected tensor direction when `nc` is present.
        logical :: periodic
        integer :: i, idx, counterpart, layer, n_dir, stride, period_count, d
        integer :: seam_multiplicity
        real(rk) :: period, knot_tol, component_diameter, weight_scale
        real(rk) :: comparison_tol, ulp_a, ulp_b, roundoff_factor

        periodic = .false.
        if (size(Xc,1) < 1 .or. size(Xc,2) < 1) return
        if (.not. all(ieee_is_finite(Xc))) return
        if (present(nc) .neqv. present(dir)) return

        if (present(nc)) then
            if (size(nc) < 1 .or. any(nc < 1)) return
            if (dir < 1 .or. dir > size(nc)) return
            if (product(nc) /= size(Xc,1)) return
            n_dir = nc(dir)
            stride = 1
            if (dir > 1) stride = product(nc(1:dir-1))
        else
            n_dir = size(Xc,1)
            stride = 1
        end if

        if (.not. valid_knot_vector(knot, n_dir, degree) .or. degree < 1) return
        if (present(Wc)) then
            if (size(Wc) /= size(Xc,1)) return
            if (.not. all(ieee_is_finite(Wc) .and. Wc > 0.0_rk)) return
        end if

        period_count = n_dir - degree
        period = knot(n_dir+1) - knot(degree+1)
        knot_tol = knot_tolerance(knot, 1, size(knot))
        do i = 1, size(knot) - period_count
            if (abs(knot(i+period_count) - knot(i) - period) > knot_tol) return
        end do

        seam_multiplicity = compute_multiplicity(knot, knot(degree+1))
        if (seam_multiplicity /= compute_multiplicity(knot, knot(n_dir+1)) .or. &
            seam_multiplicity > degree) return

        roundoff_factor = 128.0_rk*real(degree+1, rk)
        do d = 1, size(Xc,2)
            component_diameter = maxval(Xc(:,d)) - minval(Xc(:,d))
            do idx = 0, size(Xc,1) - 1
                layer = mod(idx/stride, n_dir) + 1
                if (layer <= period_count) cycle
                counterpart = idx - period_count*stride
                ulp_a = abs(Xc(idx+1,d) - nearest(&
                    Xc(idx+1,d),&
                    merge(-1.0_rk,1.0_rk,Xc(idx+1,d) > 0.0_rk)))
                ulp_b = abs(Xc(counterpart+1,d) - nearest(&
                    Xc(counterpart+1,d),&
                    merge(-1.0_rk,1.0_rk,Xc(counterpart+1,d) > 0.0_rk)))
                comparison_tol = roundoff_factor*max(&
                    ulp_a,&
                    ulp_b,&
                    epsilon(1.0_rk)*component_diameter)
                if (abs(Xc(idx+1,d) - Xc(counterpart+1,d)) > comparison_tol) return
            end do
        end do

        if (present(Wc)) then
            weight_scale = maxval(Wc)
            do idx = 0, size(Xc,1) - 1
                layer = mod(idx/stride, n_dir) + 1
                if (layer <= period_count) cycle
                counterpart = idx - period_count*stride
                ulp_a = abs(Wc(idx+1) - nearest(Wc(idx+1),-1.0_rk))
                ulp_b = abs(Wc(counterpart+1) - nearest(Wc(counterpart+1),-1.0_rk))
                comparison_tol = roundoff_factor*max(&
                    ulp_a,&
                    ulp_b,&
                    epsilon(1.0_rk)*weight_scale)
                if (abs(Wc(idx+1) - Wc(counterpart+1)) > comparison_tol) return
            end do
        end if

        periodic = .true.
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the lower active parameter bound \(U_p\).
    !! Input is assumed to satisfy [[valid_knot_vector]].
    pure function knot_start(knot, nc, degree) result(a)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: nc, degree
        real(rk) :: a

        a = knot(degree+1)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return distinct stored knots on the closed active interval.
    !! A nondecreasing knot vector is a mathematical precondition. Invalid
    !! degree/count/size metadata returns an allocated zero-length array; other
    !! knot-vector invariants are not checked here.
    pure function active_knots(knot, nc, degree) result(Xth)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: nc, degree
        real(rk), allocatable :: Xth(:)
        integer :: nunique, i

        if (nc < 1 .or. degree < 0 .or. degree >= nc .or. size(knot) < nc + degree + 1) then
            allocate(Xth(0))
            return
        end if

        nunique = 1
        do i = degree + 2, nc + 1
            if (knot(i) > knot(i-1)) nunique = nunique + 1
        end do

        allocate(Xth(nunique))
        Xth(1) = knot(degree+1)
        nunique = 1
        do i = degree + 2, nc + 1
            if (knot(i) > Xth(nunique)) then
                nunique = nunique + 1
                Xth(nunique) = knot(i)
            end if
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Count positive-length knot spans in the active interval.
    !! Repeated knots do not create zero-measure elements. A valid
    !! nondecreasing knot vector is assumed; only degree/count/size metadata is
    !! checked.
    pure function active_span_count(knot, nc, degree) result(nspan)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: nc, degree
        integer :: nspan
        integer :: i

        nspan = 0
        if (nc < 1 .or. degree < 0 .or. degree >= nc) return
        if (size(knot) < nc + degree + 1) return

        do i = degree + 1, nc
            if (knot(i+1) > knot(i)) nspan = nspan + 1
        end do
    end function
    !===============================================================================


    !> Locate a one-based positive-length element in the stored knot vector.
    !! Return zero for an absent element or invalid degree/count/size metadata.
    pure integer function active_span(knot, nc, degree, element) result(span)
        real(rk), intent(in) :: knot(:)
        integer, intent(in) :: nc, degree, element
        integer :: i, count_
        span = 0
        if (degree < 0 .or. nc <= degree .or. element < 1) return
        if (size(knot) < nc+degree+1) return
        count_ = 0
        do i = degree+1, nc
            if (knot(i+1) <= knot(i)) cycle
            count_ = count_+1
            if (count_ /= element) cycle
            span = i
            return
        end do
    end function active_span


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return multiplicities of distinct knots on the active interval.
    !! The result includes both active endpoints and follows their order in the
    !! knot vector. A valid nondecreasing vector is assumed. Invalid
    !! degree/count/size metadata or a nonpositive active interval returns an
    !! allocated zero-length array; other invariants are not checked.
    pure function active_knot_multiplicity(knot, nc, degree) result(multiplicity)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: nc, degree
        integer, allocatable :: multiplicity(:)
        integer :: i, j, nunique, run
        real(rk) :: a, b, x

        if (nc < 1 .or. degree < 0 .or. degree >= nc .or. size(knot) < nc + degree + 1) then
            allocate(multiplicity(0))
            return
        end if

        a = knot_start(knot, nc, degree)
        b = knot_end(knot, nc, degree)
        if (b <= a) then
            allocate(multiplicity(0))
            return
        end if

        nunique = 0
        i = 1
        do while (i <= size(knot))
            x = knot(i)
            j = i + 1
            do while (j <= size(knot))
                if (knot(j) /= x) exit
                j = j + 1
            end do
            if (x >= a .and. x <= b) nunique = nunique + 1
            i = j
        end do

        allocate(multiplicity(nunique))

        nunique = 0
        i = 1
        do while (i <= size(knot))
            x = knot(i)
            j = i + 1
            do while (j <= size(knot))
                if (knot(j) /= x) exit
                j = j + 1
            end do
            if (x >= a .and. x <= b) then
                nunique = nunique + 1
                run = j - i
                multiplicity(nunique) = run
            end if
            i = j
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the upper active parameter bound \(U_{n_c}\).
    !! Input is assumed to satisfy [[valid_knot_vector]].
    pure function knot_end(knot, nc, degree) result(b)
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: nc, degree
        real(rk) :: b

        b = knot(nc+1)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Infer a two-dimensional degree and control-net shape from knot lengths.
    !!
    !! Candidates satisfy
    !! \(n_{c,d}=\operatorname{size}(U_d)-p_d-1\) and
    !! \(n_{c,1}n_{c,2}=n_{cp}\). The first candidate uses
    !! \(p_d=s_{d,0}-1\), where \(s_{d,0}\) is the exact run length of the first
    !! knot. If that candidate is valid, it takes precedence; this is the usual
    !! open-clamped inference and no ambiguity search is performed. Otherwise,
    !! all degrees are searched and `ok` is true only for exactly one valid
    !! candidate. Ignore `degree` and `nc` when `ok` is false; after an ambiguous
    !! search they can retain the first candidate found.
    pure subroutine infer_knot_shape_2d(knot1, knot2, ncp, degree, nc, ok)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:)
            !! Directional knot vectors.
        integer, intent(in) :: ncp
            !! Flattened control-point count.
        integer, intent(out) :: degree(2), nc(2)
            !! Inferred degrees and control counts; meaningful only when `ok` is true.
        logical, intent(out) :: ok
            !! True for a valid preferred candidate or a unique fallback candidate.
        integer, allocatable :: mul1(:), mul2(:)
        integer :: p1, p2, n1, n2, nsolution

        degree = 0
        nc = 0
        ok = .false.

        if (ncp <= 0 .or. size(knot1) < 2 .or. size(knot2) < 2) return

        mul1 = compute_multiplicity(knot1)
        mul2 = compute_multiplicity(knot2)
        if (size(mul1) == 0 .or. size(mul2) == 0) return
        p1 = mul1(1) - 1
        p2 = mul2(1) - 1
        n1 = size(knot1) - p1 - 1
        n2 = size(knot2) - p2 - 1
        if (n1*n2 == ncp) then
            if (valid_knot_vector(knot1, n1, p1) .and. valid_knot_vector(knot2, n2, p2)) then
                degree = [p1, p2]
                nc = [n1, n2]
                ok = .true.
                return
            end if
        end if

        nsolution = 0
        do p1 = 0, size(knot1)-2
            n1 = size(knot1) - p1 - 1
            if (.not. valid_knot_vector(knot1, n1, p1)) cycle
            if (mod(ncp, n1) /= 0) cycle
            n2 = ncp/n1
            p2 = size(knot2) - n2 - 1
            if (.not. valid_knot_vector(knot2, n2, p2)) cycle
            nsolution = nsolution + 1
            if (nsolution == 1) then
                degree = [p1, p2]
                nc = [n1, n2]
            end if
        end do
        ok = nsolution == 1
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Infer a three-dimensional degree and control-lattice shape from knot lengths.
    !!
    !! Candidates satisfy
    !! \(n_{c,d}=\operatorname{size}(U_d)-p_d-1\) and
    !! \(n_{c,1}n_{c,2}n_{c,3}=n_{cp}\). The exact first-knot run lengths first
    !! define the preferred open-clamped candidate \(p_d=s_{d,0}-1\). A valid
    !! preferred candidate is returned without an ambiguity search. If it is
    !! invalid, exhaustive search succeeds only when exactly one candidate is
    !! valid. Ignore `degree` and `nc` when `ok` is false.
    pure subroutine infer_knot_shape_3d(knot1, knot2, knot3, ncp, degree, nc, ok)
        real(rk), intent(in), contiguous :: knot1(:), knot2(:), knot3(:)
            !! Directional knot vectors.
        integer, intent(in) :: ncp
            !! Flattened control-point count.
        integer, intent(out) :: degree(3), nc(3)
            !! Inferred degrees and control counts; meaningful only when `ok` is true.
        logical, intent(out) :: ok
            !! True for a valid preferred candidate or a unique fallback candidate.
        integer, allocatable :: mul1(:), mul2(:), mul3(:)
        integer :: p1, p2, p3, n1, n2, n3, nsolution

        degree = 0
        nc = 0
        ok = .false.

        if (ncp <= 0 .or. size(knot1) < 2 .or. size(knot2) < 2 .or. size(knot3) < 2) return

        mul1 = compute_multiplicity(knot1)
        mul2 = compute_multiplicity(knot2)
        mul3 = compute_multiplicity(knot3)
        if (size(mul1) == 0 .or. size(mul2) == 0 .or. size(mul3) == 0) return
        p1 = mul1(1) - 1
        p2 = mul2(1) - 1
        p3 = mul3(1) - 1
        n1 = size(knot1) - p1 - 1
        n2 = size(knot2) - p2 - 1
        n3 = size(knot3) - p3 - 1
        if (n1*n2*n3 == ncp) then
            if (valid_knot_vector(knot1, n1, p1) .and. valid_knot_vector(knot2, n2, p2) .and. &
                valid_knot_vector(knot3, n3, p3)) then
                degree = [p1, p2, p3]
                nc = [n1, n2, n3]
                ok = .true.
                return
            end if
        end if

        nsolution = 0
        do p1 = 0, size(knot1)-2
            n1 = size(knot1) - p1 - 1
            if (.not. valid_knot_vector(knot1, n1, p1)) cycle
            do p2 = 0, size(knot2)-2
                n2 = size(knot2) - p2 - 1
                if (.not. valid_knot_vector(knot2, n2, p2)) cycle
                if (mod(ncp, n1*n2) /= 0) cycle
                n3 = ncp/(n1*n2)
                p3 = size(knot3) - n3 - 1
                if (.not. valid_knot_vector(knot3, n3, p3)) cycle
                nsolution = nsolution + 1
                if (nsolution == 1) then
                    degree = [p1, p2, p3]
                    nc = [n1, n2, n3]
                end if
            end do
        end do
        ok = nsolution == 1
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the complete degree-`degree` B-spline basis at one parameter.
    !!
    !! The result has `nc` entries but at most `degree+1` are nonzero.
    !! Parameters outside the active interval return zeros, except values within
    !! [[knot_tolerance]] of an endpoint, which are snapped to that endpoint.
    !!
    pure function basis_bspline(Xt, knot, nc, degree) result(B)
        integer, intent(in) :: degree
            !! Polynomial degree.
        real(rk), intent(in), contiguous :: knot(:)
            !! Valid knot vector.
        integer, intent(in) :: nc
            !! Number of basis functions/control points.
        real(rk), intent(in) :: Xt
            !! Parameter value.
        real(rk) :: B(nc)
        call basis_bspline_der_order(Xt,knot,nc,degree,0,B)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate all nonzero B-spline basis derivatives through `nder`.
    !! Implements Algorithm A2.3 of Piegl and Tiller, The NURBS Book,
    !! 2nd edition. Rows above `degree` are identically zero.
    !! The caller supplies `ders(0:nder,0:degree)`; on return `first` is the
    !! one-based index associated with column zero.
    !! For \(k\le p\), row \(k\) contains
    !! \(\mathrm d^kN_{i,p}/\mathrm du^k\); rows \(k>p\) are zero.
    !!
    pure subroutine basis_bspline_der_all_active(Xt, knot, nc, degree, nder, first, ders, side, wrap)
        integer, intent(in) :: degree
            !! Polynomial degree.
        integer, intent(in) :: nc
            !! Number of basis functions.
        integer, intent(in) :: nder
            !! Highest derivative order requested; must be nonnegative.
        real(rk), intent(in), contiguous :: knot(:)
            !! Valid knot vector.
        real(rk), intent(in) :: Xt
            !! Parameter value on, or within tolerance of, the active interval.
        integer, intent(out) :: first
            !! First globally active one-based basis index.
        real(rk), intent(out), contiguous :: ders(0:,0:)
            !! Derivative table; row `k` contains \(N^{(k)}\).
        integer, intent(in), optional :: side
            !! At an interior knot use the left trace (-1) or right trace (0, +1).
        logical, intent(in), optional :: wrap
            !! Enable wrapping before selecting the span.
        integer :: span, j, r, k, nder_eff
        integer :: s1, s2, rkidx, pk, j1, j2, tmp_s
        real(rk) :: left(0:degree), right(0:degree)
        real(rk) :: ndu(0:degree, 0:degree), a(0:1, 0:degree)
        real(rk) :: saved, temp, denom, d, factor, u
        logical :: inside

        first = 1
        ders = 0.0_rk
        if (degree < 0 .or. nder < 0 .or. nc < 1 .or. degree >= nc) return
        if (size(knot) < nc + degree + 1 .or. ubound(ders,1) < nder .or. ubound(ders,2) < degree) return
        call active_span_at(Xt, knot, nc, degree, span, u, inside, side, wrap)
        if (.not. inside) return
        ndu = 0.0_rk
        left = 0.0_rk
        right = 0.0_rk
        ndu(0, 0) = 1.0_rk

        do j = 1, degree
            left(j) = u - knot(span + 1 - j)
            right(j) = knot(span + j) - u
            saved = 0.0_rk
            do r = 0, j-1
                ndu(j, r) = right(r+1) + left(j-r)
                if (ndu(j, r) /= 0.0_rk) then
                    temp = ndu(r, j-1)/ndu(j, r)
                else
                    temp = 0.0_rk
                end if
                ndu(r, j) = saved + right(r+1)*temp
                saved = left(j-r)*temp
            end do
            ndu(j, j) = saved
        end do

        do j = 0, degree
            ders(0, j) = ndu(j, degree)
        end do

        nder_eff = min(nder, degree)
        if (nder_eff > 0) then
            do r = 0, degree
                a = 0.0_rk
                a(0, 0) = 1.0_rk
                s1 = 0
                s2 = 1
                do k = 1, nder_eff
                    d = 0.0_rk
                    rkidx = r - k
                    pk = degree - k

                    if (r >= k) then
                        denom = ndu(pk+1, rkidx)
                        if (denom /= 0.0_rk) then
                            a(s2, 0) = a(s1, 0)/denom
                            d = a(s2, 0)*ndu(rkidx, pk)
                        else
                            a(s2, 0) = 0.0_rk
                        end if
                    end if

                    if (rkidx >= -1) then
                        j1 = 1
                    else
                        j1 = -rkidx
                    end if

                    if (r - 1 <= pk) then
                        j2 = k - 1
                    else
                        j2 = degree - r
                    end if

                    do j = j1, j2
                        denom = ndu(pk+1, rkidx+j)
                        if (denom /= 0.0_rk) then
                            a(s2, j) = (a(s1, j) - a(s1, j-1))/denom
                            d = d + a(s2, j)*ndu(rkidx+j, pk)
                        else
                            a(s2, j) = 0.0_rk
                        end if
                    end do

                    if (r <= pk) then
                        denom = ndu(pk+1, r)
                        if (denom /= 0.0_rk) then
                            a(s2, k) = -a(s1, k-1)/denom
                            d = d + a(s2, k)*ndu(r, pk)
                        else
                            a(s2, k) = 0.0_rk
                        end if
                    end if

                    ders(k, r) = d
                    tmp_s = s1
                    s1 = s2
                    s2 = tmp_s
                end do
            end do

            factor = real(degree, rk)
            do k = 1, nder_eff
                do j = 0, degree
                    ders(k, j) = ders(k, j)*factor
                end do
                factor = factor*real(degree-k, rk)
            end do
        end if

        first = span - degree
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return local basis values and optional first or second derivatives.
    !! This convenience overload requests at most order two from
    !! [[basis_bspline_der_all_active]] and never forms a dense `nc`-vector.
    pure subroutine basis_bspline_active_ders(Xt, knot, nc, degree, nder, first, N, dN, d2N)
        integer, intent(in) :: degree, nc, nder
        real(rk), intent(in), contiguous :: knot(:)
        real(rk), intent(in) :: Xt
        integer, intent(out) :: first
        real(rk), intent(out), optional :: N(0:degree), dN(0:degree), d2N(0:degree)
        real(rk) :: ders(0:2,0:degree)

        first = 1
        if (present(N)) N = 0.0_rk
        if (present(dN)) dN = 0.0_rk
        if (present(d2N)) d2N = 0.0_rk
        if (nder < 0) return

        call basis_bspline_der_all_active(&
            Xt     = Xt,&
            knot   = knot,&
            nc     = nc,&
            degree = degree,&
            nder   = min(2, nder),&
            first  = first,&
            ders   = ders)
        if (present(N)) N = ders(0,:)
        if (present(dN) .and. nder >= 1) dN = ders(1,:)
        if (present(d2N) .and. nder >= 2) d2N = ders(2,:)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Scatter local values and derivatives into optional dense vectors.
    !! Every present output has length `nc`; entries outside local support are
    !! set to zero.
    pure subroutine basis_bspline_ders(Xt, knot, nc, degree, nder, B, dB, d2B)
        integer, intent(in) :: degree, nc, nder
        real(rk), intent(in), contiguous :: knot(:)
        real(rk), intent(in) :: Xt
        real(rk), intent(out), optional :: B(nc), dB(nc), d2B(nc)
        real(rk) :: N(0:degree), dN(0:degree), d2N(0:degree)
        integer :: first, j, idx

        if (present(B))  B  = 0.0_rk
        if (present(dB)) dB = 0.0_rk
        if (present(d2B)) d2B = 0.0_rk
        if (nc == 0) return

        call basis_bspline_active_ders(Xt, knot, nc, degree, nder, first, N, dN, d2N)
        do j = 0, degree
            idx = first + j
            if (idx >= 1 .and. idx <= nc) then
                if (present(B)) B(idx) = N(j)
                if (present(dB)) dB(idx) = dN(j)
                if (present(d2B)) d2B(idx) = d2N(j)
            end if
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate one arbitrary derivative order as a dense basis vector.
    !!
    !! `B` is zero when `nder>degree`, when metadata are invalid, or when `Xt`
    !! lies outside the active interval.
    !!
    pure subroutine basis_bspline_der_order(Xt, knot, nc, degree, nder, B)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        integer, intent(in) :: nder
            !! Nonnegative derivative order.
        real(rk), intent(in), contiguous :: knot(:)
        real(rk), intent(in) :: Xt
        real(rk), intent(out), contiguous :: B(:)
        real(rk) :: B_active(0:degree)
        integer :: j, first, idx

        B = 0.0_rk
        if (nc == 0 .or. degree < 0 .or. degree >= nc .or. nder < 0) return
        if (size(B) < nc .or. size(knot) < nc + degree + 1) return

        call basis_bspline_der_order_active(Xt, knot, nc, degree, nder, first, B_active)
        do j = 0, degree
            idx = first + j
            if (idx >= 1 .and. idx <= nc) B(idx) = B_active(j)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate one arbitrary derivative order on local support only.
    !!
    !! The fixed-size result `B(0:degree)` avoids a dense temporary. `first`
    !! maps `B(j)` to global basis index `first+j`.
    !!
    pure subroutine basis_bspline_der_order_active(Xt, knot, nc, degree, nder, first, B)
        integer, intent(in) :: degree
        integer, intent(in) :: nc
        integer, intent(in) :: nder
            !! Nonnegative derivative order.
        real(rk), intent(in), contiguous :: knot(:)
        real(rk), intent(in) :: Xt
        integer, intent(out) :: first
        real(rk), intent(out) :: B(0:degree)
        real(rk) :: ders(0:max(0, min(nder, degree)),0:degree)

        first = 1
        B = 0.0_rk
        if (nc == 0 .or. degree < 0 .or. degree >= nc .or. nder < 0) return
        if (size(knot) < nc + degree + 1) return
        if (nder > degree) return

        call basis_bspline_der_all_active(&
            Xt     = Xt,&
            knot   = knot,&
            nc     = nc,&
            degree = degree,&
            nder   = nder,&
            first  = first,&
            ders   = ders)
        B = ders(nder,:)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute one tensor-product B-spline or rational basis derivative on its local support.
    !! The three directional derivative tables include rows zero through the
    !! requested derivative order. Unused directions are represented by a
    !! \(1\times1\) table containing one. Rational derivatives use a
    !! multivariate quotient recurrence that extends the rational-surface
    !! recurrence in Piegl and Tiller, *The NURBS Book*, 2nd edition, Algorithm
    !! A4.4 to three tensor-product directions.
    !! The derivative multi-index is given by the last row present in each
    !! directional table. For example, upper row bounds `(1,0,2)` request
    !! \(\partial^3/(\partial\xi_1\partial\xi_3^2)\).
    !! `values` uses direction-1-fastest local ordering. Optional `elem` selects
    !! zero-order values in list order, retaining repeats and zeroing invalid or inactive IDs.
    !! The rational path evaluates numerator and denominator derivatives in a
    !! common exact radix-power gauge, so a uniform finite weight rescaling does
    !! not alter the recurrence.
    !!
    pure subroutine tensor_basis_derivative_local(&
        first1,&
        nc1,&
        ders1,&
        first2,&
        nc2,&
        ders2,&
        first3,&
        nc3,&
        ders3,&
        values,&
        Wc, elem, checked)
        integer, intent(in) :: first1, nc1, first2, nc2, first3, nc3
        real(rk), intent(in), contiguous :: ders1(0:,0:), ders2(0:,0:), ders3(0:,0:)
        real(rk), intent(out), contiguous :: values(0:)
            !! Supported derivative values, or selected values when `elem` is supplied.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Global rational weights, or one local weight per `elem` entry; absent for polynomial values.
        integer, intent(in), optional :: elem(:) !! Zero-order selection with matching local weights.
        logical, intent(in), optional :: checked !! Signal lost rational derivative workspace with nonfinite values.
        integer :: p1, p2, p3, o1, o2, o3, weight_exponent, nlocal, pass, weight_index
        integer :: a1, a2, a3, b1, b2, b3, l1, l2, l3, local, global
        real(rk) :: coefficient1, coefficient2, coefficient3, numerator, denominator_value, factors(4), term
        logical :: strict

        values = 0.0_rk
        strict = .false.
        if (present(checked)) strict = checked
        p1 = ubound(ders1,2)
        p2 = ubound(ders2,2)
        p3 = ubound(ders3,2)
        o1 = ubound(ders1,1)
        o2 = ubound(ders2,1)
        o3 = ubound(ders3,1)
        if (nc1 < 1 .or. nc2 < 1 .or. nc3 < 1) return
        if (first1 < 1 .or. first2 < 1 .or. first3 < 1) return
        if (first1 + p1 > nc1 .or. first2 + p2 > nc2 .or. first3 + p3 > nc3) return
        nlocal = (p1+1)*(p2+1)*(p3+1)
        if (present(elem)) nlocal = size(elem)
        if (size(values) < nlocal) return

        if (.not. present(Wc)) then
            if (present(elem)) return
            do concurrent (l3 = 0:p3, l2 = 0:p2, l1 = 0:p1)
                values(l1+(p1+1)*(l2+(p2+1)*l3)) = ders1(o1,l1)*ders2(o2,l2)*ders3(o3,l3)
            end do
            return
        end if
        if (present(elem)) then
            if (o1+o2+o3/=0 .or. size(Wc)/=nlocal) return
            if (.not.valid_weights(Wc)) return
        else
            if (size(Wc) < nc1*nc2*nc3) return
        end if

        ! Scale nonzero weighted contributions before summing Q/W.
        if (o1+o2+o3 == 0) then
            weight_exponent = -huge(0)
            do pass = 1, 2
                do local = 0, nlocal-1
                    if (present(elem)) then
                        global = elem(local+1)
                        if (global<1 .or. global>nc1*nc2*nc3) cycle
                        l1 = mod(global-1,nc1)+1-first1
                        l2 = mod((global-1)/nc1,nc2)+1-first2
                        l3 = (global-1)/(nc1*nc2)+1-first3
                        if (l1<0 .or. l1>p1 .or. l2<0 .or. l2>p2 .or. l3<0 .or. l3>p3) cycle
                        weight_index = local+1
                    else
                        l1 = mod(local,p1+1)
                        l2 = mod(local/(p1+1),p2+1)
                        l3 = local/((p1+1)*(p2+1))
                        global = first1+l1+nc1*((first2+l2-1)+nc2*(first3+l3-1))
                        weight_index = global
                    end if
                    factors = [Wc(weight_index),ders1(0,l1),ders2(0,l2),ders3(0,l3)]
                    if (pass==1) then
                        if (all(factors/=0.0_rk)) weight_exponent = max(weight_exponent, &
                            sum(exponent(factors))+exponent(product(fraction(factors))))
                    else
                        values(local) = scale(product(fraction(factors)),sum(exponent(factors))-weight_exponent)
                    end if
                end do
                if (weight_exponent==-huge(0)) return
            end do
            denominator_value = sum(values)
            if (denominator_value > 0.0_rk .and. ieee_is_finite(denominator_value)) then
                values = values/denominator_value
            else
                values = 0.0_rk
            end if
            return
        end if

        weight_exponent = -huge(0)
        do l3 = 0, p3
            do l2 = 0, p2
                global = first1 + nc1*((first2+l2-1) + nc2*(first3+l3-1))
                weight_exponent = max(weight_exponent,maxval(exponent(Wc(global:global+p1))))
            end do
        end do
        if (strict) values=ieee_value(1.0_rk,ieee_quiet_nan)
        block
            real(rk) :: denominator(0:o1,0:o2,0:o3), quotient(0:o1,0:o2,0:o3)
            denominator = 0.0_rk
            do a3 = 0, o3
                do a2 = 0, o2
                    do a1 = 0, o1
                        do l3 = 0, p3
                            do l2 = 0, p2
                                do l1 = 0, p1
                                    global = first1 + l1 + nc1*((first2+l2-1) + nc2*(first3+l3-1))
                                    denominator(a1,a2,a3) = denominator(a1,a2,a3) + &
                                        scale(Wc(global), -weight_exponent)*&
                                        ders1(a1,l1)*ders2(a2,l2)*ders3(a3,l3)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            if (denominator(0,0,0) <= 0.0_rk .or. .not. ieee_is_finite(denominator(0,0,0))) return

            do l3 = 0, p3
                do l2 = 0, p2
                    do l1 = 0, p1
                        global = first1 + l1 + nc1*((first2+l2-1) + nc2*(first3+l3-1))
                        quotient = 0.0_rk
                        do a3 = 0, o3
                            do a2 = 0, o2
                                do a1 = 0, o1
                                    numerator = scale(Wc(global), -weight_exponent)*&
                                        ders1(a1,l1)*ders2(a2,l2)*ders3(a3,l3)
                                    coefficient3 = 1.0_rk
                                    do b3 = 0, a3
                                        coefficient2 = 1.0_rk
                                        do b2 = 0, a2
                                            coefficient1 = 1.0_rk
                                            do b1 = 0, a1
                                                term = coefficient1*coefficient2*coefficient3*&
                                                    denominator(b1,b2,b3)*quotient(a1-b1,a2-b2,a3-b3)
                                                if (strict .and. (.not.ieee_is_finite(term) .or. &
                                                    (denominator(b1,b2,b3)/=0.0_rk .and. &
                                                     quotient(a1-b1,a2-b2,a3-b3)/=0.0_rk .and. abs(term)<tiny(term)))) then
                                                    return
                                                end if
                                                numerator = numerator - term
                                                coefficient1 = coefficient1*real(a1-b1,rk)/real(b1+1,rk)
                                            end do
                                            coefficient2 = coefficient2*real(a2-b2,rk)/real(b2+1,rk)
                                        end do
                                        coefficient3 = coefficient3*real(a3-b3,rk)/real(b3+1,rk)
                                    end do
                                    quotient(a1,a2,a3) = numerator/denominator(0,0,0)
                                    if (strict .and. (.not.ieee_is_finite(quotient(a1,a2,a3)) .or. &
                                        (numerator/=0.0_rk .and. abs(quotient(a1,a2,a3))<tiny(numerator)))) then
                                        return
                                    end if
                                end do
                            end do
                        end do
                        local = l1 + (p1+1)*(l2 + (p2+1)*l3)
                        values(local) = quotient(o1,o2,o3)
                    end do
                end do
            end do
        end block
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate the dense B-spline basis and its first derivative at one parameter.
    !!
    !! Both outputs have length `nc`; entries outside the local support are
    !! zero. This is the value-and-derivative overload of [[basis_bspline_der]].
    pure subroutine basis_bspline_der_A(Xt, knot, nc, degree, dB, B)
        integer, intent(in) :: degree
            !! Polynomial degree \(p\).
        real(rk), intent(in), contiguous :: knot(:)
            !! Complete knot vector of size `nc+degree+1`.
        integer, intent(in) :: nc
            !! Number of basis functions.
        real(rk), intent(in) :: Xt
            !! Evaluation parameter.
        real(rk), intent(out) :: dB(nc)
            !! First derivatives \(\mathrm dN_{i,p}/\mathrm du\).
        real(rk), intent(out) :: B(nc)
            !! Basis values \(N_{i,p}\).
        call basis_bspline_ders(Xt, knot, nc, degree, 1, B=B, dB=dB)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate only the dense first derivative of a B-spline basis.
    !!
    !! The output has length `nc`; entries outside the local support are zero.
    pure subroutine basis_bspline_der_B(Xt, knot, nc, degree, dB)
        integer, intent(in) :: degree
            !! Polynomial degree \(p\).
        real(rk), intent(in), contiguous :: knot(:)
            !! Complete knot vector of size `nc+degree+1`.
        integer, intent(in) :: nc
            !! Number of basis functions.
        real(rk), intent(in) :: Xt
            !! Evaluation parameter.
        real(rk), intent(out) :: dB(nc)
            !! First derivatives \(\mathrm dN_{i,p}/\mathrm du\).
        call basis_bspline_ders(Xt, knot, nc, degree, 1, dB=dB)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate dense B-spline values and derivatives through second order.
    !!
    !! Each output has length `nc`; unsupported derivatives, including all
    !! second derivatives when `degree<2`, are returned as zero.
    pure subroutine basis_bspline_2der_A(Xt, knot, nc, degree, d2B, dB, B)
        integer, intent(in) :: degree
            !! Polynomial degree \(p\).
        real(rk), intent(in), contiguous :: knot(:)
            !! Complete knot vector of size `nc+degree+1`.
        integer, intent(in) :: nc
            !! Number of basis functions.
        real(rk), intent(in) :: Xt
            !! Evaluation parameter.
        real(rk), intent(out) :: d2B(nc)
            !! Second derivatives \(\mathrm d^2N_{i,p}/\mathrm du^2\).
        real(rk), intent(out) :: dB(nc)
            !! First derivatives \(\mathrm dN_{i,p}/\mathrm du\).
        real(rk), intent(out) :: B(nc)
            !! Basis values \(N_{i,p}\).
        call basis_bspline_ders(Xt, knot, nc, degree, 2, B=B, dB=dB, d2B=d2B)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate dense first and second derivatives of a B-spline basis.
    !!
    !! Values are omitted; both derivative arrays have length `nc` and are
    !! zero outside the local support.
    pure subroutine basis_bspline_2der_B(Xt, knot, nc, degree, d2B, dB)
        integer, intent(in) :: degree
            !! Polynomial degree \(p\).
        real(rk), intent(in), contiguous :: knot(:)
            !! Complete knot vector of size `nc+degree+1`.
        integer, intent(in) :: nc
            !! Number of basis functions.
        real(rk), intent(in) :: Xt
            !! Evaluation parameter.
        real(rk), intent(out) :: d2B(nc)
            !! Second derivatives \(\mathrm d^2N_{i,p}/\mathrm du^2\).
        real(rk), intent(out) :: dB(nc)
            !! First derivatives \(\mathrm dN_{i,p}/\mathrm du\).
        call basis_bspline_ders(Xt, knot, nc, degree, 2, dB=dB, d2B=d2B)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate only the dense second derivative of a B-spline basis.
    !!
    !! The output has length `nc`; entries outside local support and all entries
    !! for `degree<2` are zero.
    pure subroutine basis_bspline_2der_C(Xt, knot, nc, degree, d2B)
        integer, intent(in) :: degree
            !! Polynomial degree \(p\).
        real(rk), intent(in), contiguous :: knot(:)
            !! Complete knot vector of size `nc+degree+1`.
        integer, intent(in) :: nc
            !! Number of basis functions.
        real(rk), intent(in) :: Xt
            !! Evaluation parameter.
        real(rk), intent(out) :: d2B(nc)
            !! Second derivatives \(\mathrm d^2N_{i,p}/\mathrm du^2\).
        call basis_bspline_ders(Xt, knot, nc, degree, 2, d2B=d2B)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Evaluate all Bernstein basis polynomials of degree `nc-1` on `[0,1]`.
    !!
    !! The implementation uses an adjacent-term recurrence with explicit
    !! endpoint values. `nc=1` returns the constant basis.
    !!
    !! \[
    !! B_i^p(t)={p\choose i}t^i(1-t)^{p-i},
    !! \qquad i=0,\ldots,p,\quad p=nc-1.
    !! \]
    !!
    !! `nc>=1` is an unchecked caller precondition. Values outside `[0,1]`
    !! evaluate the same polynomial formula but need not be nonnegative or form
    !! a numerically stable partition. Values within `2*epsilon(rk)` of zero or
    !! one are snapped to the corresponding exact endpoint basis vector.
    !!
    pure function basis_bernstein(Xt, nc) result(B)
        real(rk), intent(in) :: Xt
            !! Bernstein parameter, normally in `[0,1]`.
        integer, intent(in) :: nc
            !! Number of basis functions; the degree is `nc-1`.
        real(rk), allocatable :: B(:)
        integer  :: p, degree
        real(rk) :: omt, ratio

        degree = nc - 1

        allocate(B(nc), source=0.0_rk)

        if (nc == 1) then
            B(1) = 1.0_rk
            return
        end if

        if (abs(Xt) < 2.0_rk*epsilon(0.0_rk)) then
            B(1) = 1.0_rk
            return
        end if

        omt = 1.0_rk - Xt
        if (abs(omt) < 2.0_rk*epsilon(0.0_rk)) then
            B(nc) = 1.0_rk
            return
        end if

        ratio = Xt/omt
        B(1) = omt**degree
        do p = 1, degree
            B(p+1) = B(p)*real(degree - p + 1, rk)*ratio/real(p, rk)
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Form the Kronecker product of two vectors.
    !!
    !! The result is ordered in blocks of `v`:
    !! `w((i-1)*size(v)+j)=u(i)*v(j)`.
    pure function kron_t1_t1(u,v) result(w)
        real(rk), intent(in), contiguous :: u(:)
            !! Left vector.
        real(rk), intent(in), contiguous :: v(:)
            !! Right vector, varying fastest in the result.
        real(rk) :: w(size(u)*size(v))
            !! Kronecker product of length `size(u)*size(v)`.
        integer :: i, j, n

        n = size(v)

        do concurrent(i = 1:size(u), j = 1:n)
            w((i-1)*n + j) = u(i)*v(j)
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Form the Kronecker product of a vector and a matrix.
    !!
    !! Every row of `A` is scaled by each entry of `u`; the matrix columns are
    !! retained and the `A` row index varies fastest inside each `u` block.
    pure function kron_t1_t2(u,A) result(B)
        real(rk), intent(in), contiguous :: u(:)
            !! Left vector of length `m`.
        real(rk), intent(in), contiguous :: A(:,:)
            !! Right matrix of shape `[r,c]`.
        real(rk) :: B(size(u)*size(A,1), size(A,2))
            !! Kronecker product of shape `[m*r,c]`.
        integer :: i, j, k, m, r, c

        m = size(u)
        r = size(A, 1)
        c = size(A, 2)

        do concurrent (i=1:m, j=1:r, k=1:c)
            B((i-1)*r + j, k) = u(i) * A(j, k)
        end do
    end function
    !===============================================================================


    !> Matrix Kronecker product, with the right factor varying fastest.
    pure function kron_t2_t2(A,B) result(C)
        real(rk), intent(in) :: A(:,:), B(:,:)
        real(rk) :: C(size(A,1)*size(B,1),size(A,2)*size(B,2))
        integer :: i,j,a_,b_
        do concurrent (j=1:size(A,2),i=1:size(A,1),b_=1:size(B,2),a_=1:size(B,1))
            C(a_+(i-1)*size(B,1),b_+(j-1)*size(B,2))=A(i,j)*B(a_,b_)
        end do
    end function kron_t2_t2

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Form a flattened three-vector tensor product.
    !!
    !! `w` varies fastest, followed by `v` and then `u`, matching the ordering
    !! used by the ForCAD tensor-product geometry modules.
    pure function kron3(u, v, w) result(out)
        real(rk), intent(in), contiguous :: u(:)
            !! First tensor factor, varying slowest.
        real(rk), intent(in), contiguous :: v(:)
            !! Second tensor factor.
        real(rk), intent(in), contiguous :: w(:)
            !! Third tensor factor, varying fastest.
        real(rk) :: out(size(u)*size(v)*size(w))
            !! Flattened product `u(i)*v(j)*w(k)`.
        integer :: i, j, k, nv, nw

        nv = size(v)
        nw = size(w)

        do concurrent(i = 1:size(u), j = 1:nv, k = 1:nw)
            out(((i-1)*nv + j -1)*nw + k) = u(i) * v(j) * w(k)
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Expand a scalar operator independently over `ncoord` components.
    !!
    !! The result is the control-point-major representation of
    !! \(A\otimes I_{ncoord}\):
    !! `B((i-1)*ncoord+r,(j-1)*ncoord+r)=A(i,j)`. Nonpositive `ncoord` returns `B(0,0)`.
    pure function kron_eye(A, ncoord) result(B)
        real(rk), intent(in), contiguous :: A(:,:)
            !! Scalar operator `[m,n]`.
        integer, intent(in) :: ncoord
            !! Number of independent components.
        real(rk), allocatable :: B(:,:)
        integer :: i, j, r

        if (ncoord <= 0) then
            allocate(B(0,0))
            return
        end if

        allocate(B(size(A,1)*ncoord, size(A,2)*ncoord), source=0.0_rk)
        do concurrent (r = 1:ncoord, i = 1:size(A,1), j = 1:size(A,2))
            B((i-1)*ncoord + r, (j-1)*ncoord + r) = A(i,j)
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute `C=matmul(A,B)` by scanning structural nonzeros of `A`.
    !! Intended for local-support refinement maps, where A is structurally sparse
    !! but the existing API still needs a dense transformation matrix. NaN
    !! coefficients remain structural and therefore propagate rather than being
    !! silently dropped. Incompatible inner dimensions return `C(0,0)`.
    !! This is not a compressed sparse representation: every entry of dense
    !! `A` is inspected and exact zeros alone are skipped. For dimensions
    !! \(m\times k\) and \(k\times n\), work is
    !! \(O(mk+\operatorname{nnz}(A)n+mn)\), including the scan and result
    !! initialization. `A`, `B`, and the returned `C` all remain dense; result
    !! storage is \(O(mn)\).
    !!
    pure subroutine sparse_left_matmul(A, B, C)
        real(rk), intent(in), contiguous :: A(:,:)
            !! Left matrix, commonly a refinement transformation.
        real(rk), intent(in), contiguous :: B(:,:)
            !! Dense right-hand matrix.
        real(rk), allocatable, intent(out) :: C(:,:)
            !! Allocatable result of shape `[size(A,1),size(B,2)]`.
        real(rk) :: aik
        integer :: i, j, k, nrow, nmid, ncol

        if (size(A,2) /= size(B,1)) then
            allocate(C(0,0))
            return
        end if

        nrow = size(A,1)
        nmid = size(A,2)
        ncol = size(B,2)
        allocate(C(nrow, ncol))

        do concurrent (i = 1:nrow) local(aik, j, k)
            do j = 1, ncol
                C(i,j) = 0.0_rk
            end do
            do k = 1, nmid
                aik = A(i,k)
                if (structural_nonzero(aik)) then
                    do j = 1, ncol
                        C(i,j) = C(i,j) + aik*B(k,j)
                    end do
                end if
            end do
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Form all pairs from two directional coordinate vectors.
    !! Row `i+(j-1)*size(X_dir1)` contains
    !! `[X_dir1(i),X_dir2(j)]`.
    pure subroutine ndgrid2(X_dir1,X_dir2, Xt)
        real(rk), intent(in), contiguous :: X_dir1(:), X_dir2(:)
        real(rk), allocatable, intent(out) :: Xt(:,:)
        integer :: s1, s2, i, j

        s1 = size(X_dir1)
        s2 = size(X_dir2)
        allocate(Xt(s1*s2,2))
        do concurrent (j = 1:s2, i = 1:s1)
            Xt((j-1)*s1+i, 1) = X_dir1(i)
            Xt((j-1)*s1+i, 2) = X_dir2(j)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Form all triples from three directional coordinate vectors.
    !! Direction 1 varies fastest, followed by directions 2 and 3.
    pure subroutine ndgrid3(X_dir1,X_dir2,X_dir3, Xt)
        real(rk), intent(in), contiguous :: X_dir1(:), X_dir2(:), X_dir3(:)
        real(rk), allocatable, intent(out) :: Xt(:,:)
        integer :: s1, s2, s3, i, j, k, idx

        s1 = size(X_dir1)
        s2 = size(X_dir2)
        s3 = size(X_dir3)
        allocate(Xt(s1*s2*s3,3))
        do concurrent (k = 1:s3, j = 1:s2, i = 1:s1) local(idx)
            idx = ((k-1)*s2 + j - 1)*s1 + i
            Xt(idx, 1) = X_dir1(i)
            Xt(idx, 2) = X_dir2(j)
            Xt(idx, 3) = X_dir3(k)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Repeat each value according to a corresponding nonnegative count.
    !! Values retain input order. Shape mismatch or a negative count returns a
    !! zero-length array.
    pure function repelem(a, b) result(c)
        real(rk), intent(in), contiguous :: a(:)
            !! Values to repeat.
        integer, intent(in), contiguous :: b(:)
            !! Repetition counts with `size(b)==size(a)`.
        real(rk), allocatable :: c(:)
        integer :: i, l, n

        if (size(a) /= size(b) .or. any(b < 0)) then
            allocate(c(0))
            return
        end if

        allocate(c(sum(b)))
        l = 0
        do i = 1, size(a)
            n = b(i)
            if (n > 0) c(l+1:l+n) = a(i)
            l = l + n
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Allocate `n` inclusive uniformly spaced values from `a` to `b`.
    !! `n=1` returns `[a]`; `n<1` returns a zero-length array.
    pure function linspace(a, b, n) result(x)
        real(rk), intent(in) :: a, b
        integer, intent(in) :: n
        real(rk), allocatable :: x(:)

        if (n < 1) then
            allocate(x(0))
            return
        end if
        allocate(x(n))
        call fill_uniform(a, b, x)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Build curve connectivity for degree-`p` elements joined with \(C^0\).
    !! The required node count satisfies `nnode=nelem*p+1`.
    pure function cmp_elemConn_C0_L(nnode,p) result(elemConn)
      integer, intent(in) :: nnode, p
      integer, allocatable :: elemConn(:,:)
      integer :: nelem, e, j

      if (p < 1 .or. nnode < 1) then
         allocate(elemConn(0,0))
         return
      end if
      if (mod(nnode-1,p) /= 0) then
         allocate(elemConn(0,0))
         return
      end if

      nelem = (nnode-1)/p

      allocate(elemConn(nelem, p+1))
      do concurrent (e = 1:nelem, j = 0:p)
         elemConn(e, j+1) = (e-1)*p + 1 + j
      end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Build direction-1-fastest \(C^0\) tensor connectivity for a surface.
    !! The result has \((p_1+1)(p_2+1)\) local controls per element.
    pure function cmp_elemConn_C0_S(nnode1, nnode2, p1, p2) result(elemConn)
        integer, intent(in) :: nnode1, nnode2, p1, p2
        integer, allocatable :: elemConn(:,:)
        elemConn=connectivity_product(cmp_elemConn_C0_L(nnode1,p1),cmp_elemConn_C0_L(nnode2,p2),nnode1)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Build direction-1-fastest \(C^0\) tensor connectivity for a volume.
    !! The result has \((p_1+1)(p_2+1)(p_3+1)\) controls per element.
    pure function cmp_elemConn_C0_V(nnode1,nnode2,nnode3,p1,p2,p3) result(elemConn)
        integer, intent(in) :: nnode1, nnode2, nnode3, p1, p2, p3
        integer, allocatable :: elemConn(:,:)
        elemConn=connectivity_product(cmp_elemConn_C0_S(nnode1,nnode2,p1,p2), &
            cmp_elemConn_C0_L(nnode3,p3),nnode1*nnode2)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Build curve connectivity from distinct active knots and multiplicities.
    !! At an interior knot of multiplicity \(s\), adjacent degree-`p` elements
    !! share \(p-s+1\) basis functions and have guaranteed spline-space
    !! continuity \(C^{p-s}\). Special control data may produce smoother
    !! geometry than the underlying spline space. Knot values are not read;
    !! only `size(Xth)-1` determines the element count. The caller must ensure
    !! `nnode=p+1+sum(vecKnot_mul(2:nelem))` and valid multiplicities.
    !> Tensor product of element maps, with the left direction varying fastest.
    pure function connectivity_product(left,right,stride) result(elem)
        integer, intent(in) :: left(:,:), right(:,:), stride
        integer, allocatable :: elem(:,:)
        integer :: i, j, a, b
        allocate(elem(size(left,1)*size(right,1),size(left,2)*size(right,2)))
        do concurrent (j=1:size(right,1),i=1:size(left,1),b=1:size(right,2),a=1:size(left,2))
            elem(i+size(left,1)*(j-1),a+size(left,2)*(b-1))=left(i,a)+stride*(right(j,b)-1)
        end do
    end function connectivity_product

    pure subroutine cmp_elemConn_Cn_L(nnode, p, Xth, vecKnot_mul, elemConn)
        integer,  intent(in)              :: nnode, p
            !! Control count and nonnegative polynomial degree.
        real(rk), intent(in), contiguous  :: Xth(:)
            !! Active-knot metadata; only its size is used.
        integer,  intent(in), contiguous  :: vecKnot_mul(:)
            !! Run data whose entries `2:nelem` determine support shifts.
        integer,  allocatable, intent(out):: elemConn(:,:)
            !! One-based connectivity `[nelem,p+1]`, or shape `[0,0]` for invalid dimensions.

        integer :: nelem, i, j
        integer, allocatable :: pref(:)

        nelem = size(Xth)-1

        if (nnode < 1 .or. p < 0 .or. nelem < 1 .or. size(vecKnot_mul) < nelem) then
            allocate(elemConn(0,0))
            return
        end if

        allocate(elemConn(nelem, p + 1))
        allocate(pref(nelem))
        pref(1) = p + 1
        do i = 2, nelem
            pref(i) = pref(i-1) + vecKnot_mul(i)
        end do

        do concurrent (i = 1:nelem, j = 0:p)
            elemConn(i, j+1) = pref(i)-p+j
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Build surface connectivity as the tensor product of two span maps.
    !! Rows and local controls both use direction-1-fastest ordering. The knot
    !! values are not read; each directional size supplies its element count.
    !! Directional control counts and multiplicities must satisfy the scalar
    !! consistency relation documented by [[elemConn_Cn]].
    pure subroutine cmp_elemConn_Cn_S(nnode1, nnode2, p1, p2, &
        Xth1, Xth2, vecKnot_mul1, vecKnot_mul2, elemConn)
        integer,  intent(in)              :: nnode1, nnode2, p1, p2
            !! Directional control counts and nonnegative polynomial degrees.
        real(rk), intent(in), contiguous  :: Xth1(:), Xth2(:)
            !! Active-knot metadata; only the two array sizes are used.
        integer,  intent(in), contiguous  :: vecKnot_mul1(:), vecKnot_mul2(:)
            !! Stored multiplicities corresponding to `Xth1` and `Xth2`.
        integer,  allocatable, intent(out):: elemConn(:,:)
            !! Direction-1-fastest one-based connectivity, or shape `[0,0]` for invalid dimensions.

        integer, allocatable :: left(:,:), right(:,:)
        call cmp_elemConn_Cn_L(nnode1,p1,Xth1,vecKnot_mul1,left)
        call cmp_elemConn_Cn_L(nnode2,p2,Xth2,vecKnot_mul2,right)
        elemConn=connectivity_product(left,right,nnode1)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Build volume connectivity as the tensor product of three span maps.
    !! Rows and local controls both use direction-1-fastest ordering. Knot
    !! values are not read; directional array sizes supply element counts.
    !! Every directional control count and multiplicity vector must satisfy the
    !! scalar consistency relation documented by [[elemConn_Cn]].
    pure subroutine cmp_elemConn_Cn_V(nnode1,nnode2,nnode3,p1,p2,p3, &
        Xth1,Xth2,Xth3,vecKnot_mul1,vecKnot_mul2,vecKnot_mul3, elemConn)
        integer,  intent(in)              :: nnode1, nnode2, nnode3, p1, p2, p3
            !! Directional control counts and nonnegative polynomial degrees.
        real(rk), intent(in), contiguous  :: Xth1(:), Xth2(:), Xth3(:)
            !! Active-knot metadata; only the three array sizes are used.
        integer,  intent(in), contiguous  :: vecKnot_mul1(:), vecKnot_mul2(:), vecKnot_mul3(:)
            !! Stored multiplicities corresponding to the directional active knots.
        integer,  allocatable, intent(out):: elemConn(:,:)
            !! Direction-1-fastest one-based connectivity, or shape `[0,0]` for invalid dimensions.

        integer, allocatable :: left(:,:), right(:,:)
        call cmp_elemConn_Cn_S(nnode1,nnode2,p1,p2,Xth1,Xth2,vecKnot_mul1,vecKnot_mul2,left)
        call cmp_elemConn_Cn_L(nnode3,p3,Xth3,vecKnot_mul3,right)
        elemConn=connectivity_product(left,right,nnode1*nnode2)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return exact consecutive run lengths in stored knot order.
    !! For a nondecreasing knot vector, each value occupies one run and these
    !! are its mathematical knot multiplicities. For arbitrary input, equal
    !! values in separated runs produce separate entries. Equality is exact.
    pure function compute_multiplicity1(knot) result(multiplicity)
        real(rk), intent(in), contiguous :: knot(:)
        integer, allocatable :: multiplicity(:)
        integer :: i, nunique

        if (size(knot) == 0) then
            allocate(multiplicity(0))
            return
        end if

        nunique = 1
        do i = 2, size(knot)
            if (knot(i) /= knot(i-1)) nunique = nunique + 1
        end do

        allocate(multiplicity(nunique))

        multiplicity(1) = 1
        nunique = 1

        do i = 2, size(knot)
            if (knot(i) /= knot(i-1)) then
                nunique = nunique + 1
                multiplicity(nunique) = 1
            else
                multiplicity(nunique) = multiplicity(nunique) + 1
            end if
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return the longest exact consecutive run of one stored knot value.
    !! For a nondecreasing knot vector this is its knot multiplicity. For
    !! arbitrary input with separated equal runs, the result is the maximum run
    !! length rather than the total occurrence count. Zero means absent.
    pure function compute_multiplicity2(knot, Xth) result(multiplicity)
        real(rk), intent(in), contiguous :: knot(:)
        real(rk), intent(in) :: Xth
        integer :: multiplicity
        integer :: i, run_length, size_knot

        size_knot = size(knot)
        multiplicity = 0
        if (size_knot == 0) return

        i = 1
        do while (i <= size_knot)
            if (knot(i) == Xth) then
                run_length = 1
                do while ((i + run_length) <= size_knot)
                    if (knot(i + run_length) /= Xth) exit
                    run_length = run_length + 1
                end do
                if (run_length > multiplicity) then
                    multiplicity = run_length
                end if
                i = i + run_length
            else
                i = i + 1
            end if
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Construct a knot vector from breakpoints and requested continuity.
    !!
    !! Breakpoint `Xth_dir(i)` is repeated `degree-continuity(i)` times. Thus an
    !! interior continuity \(c_i\) has multiplicity \(p-c_i\); clamped
    !! endpoints use `continuity=-1` and therefore multiplicity \(p+1\).
    !! `continuity=degree` omits that breakpoint. Values are not checked for
    !! finiteness or ordering and are not sorted or normalized by this routine.
    !! To construct a standard valid knot vector, supply nondecreasing finite
    !! breakpoints and use \(-1\leq c_i\leq p-1\) for retained knots, then
    !! validate the resulting spline space. Values below `-1` produce
    !! multiplicity greater than \(p+1\) and are not rejected here.
    !!
    pure function compute_knot_vector(Xth_dir, degree, continuity) result(knot)
        real(rk), intent(in), contiguous :: Xth_dir(:)
            !! Ordered breakpoint values.
        integer, intent(in) :: degree
            !! Nonnegative polynomial degree \(p\).
        integer, intent(in), contiguous :: continuity(:)
            !! Repetition specification for every breakpoint. Shape mismatch, negative degree, or any value greater than
            !! `degree` returns an allocated empty vector.
        real(rk), allocatable :: knot(:)

        if (degree < 0 .or. size(Xth_dir) /= size(continuity) .or. any(degree - continuity < 0)) then
            allocate(knot(0))
            return
        end if

        knot = repelem(Xth_dir, (degree - continuity))
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Insert a knot without changing the represented homogeneous geometry.
    !!
    !! For `s+r<=p`, this is Algorithm A5.1 of Piegl and Tiller, equivalent to
    !! repeated Boehm knot insertion. The `s+r=p+1` extension first applies the
    !! standard insertion through multiplicity \(p\), then duplicates the
    !! interface control row to represent the same two one-sided pieces in a
    !! discontinuous \(C^{-1}\) space. Arrays use the algorithm's zero-based
    !! bounds. When requested, `T` satisfies `Qw=matmul(T,Pw)` and can transform
    !! any field in the same spline space. `r=0` returns the exact identity
    !! operation.
    !!
    pure recursive subroutine insert_knot_A_5_1(p, UP, Pw, u, k, s, r, nq, UQ, Qw, T)
        integer, intent(in) :: p
            !! Polynomial degree.
        integer, intent(in) :: k
            !! Zero-based span containing `u`.
        integer, intent(in) :: s
            !! Existing multiplicity of `u`.
        integer, intent(in) :: r
            !! Number of insertions, satisfying `s+r<=p+1`.
        real(rk), intent(in), contiguous :: UP(0:)
            !! Original valid knot vector, lower bound zero.
        real(rk), intent(in), contiguous :: Pw(0:,:)
            !! Original homogeneous or polynomial control variables by row.
        real(rk), intent(in) :: u
            !! Knot value to insert.
        real(rk), allocatable, intent(out) :: UQ(:)
            !! Refined knot vector.
        real(rk), allocatable, intent(out) :: Qw(:,:)
            !! Refined control variables.
        integer, intent(out) :: nq
            !! New highest zero-based control index; `-1` signals invalid input.
        real(rk), allocatable, intent(out), optional :: T(:,:)
            !! Optional refinement transformation.

        real(rk), allocatable :: Rw(:,:), Utemp(:), Qtemp(:,:), Ttemp(:,:)
        real(rk) :: alpha
        integer  :: i, j, L, mp, d, np, Lf, bw, nq_temp, k_full, shared, span

        d  = size(Pw, 2)
        np = size(Pw, 1)-1
        mp = np+p+1
        nq = -1

        if (p < 0 .or. np < p .or. d < 1 .or. size(UP) /= np+p+2 .or. &
            k < p .or. s < 0 .or. s > p+1 .or. k-s > np .or. r < 0 .or. s+r > p+1 .or. &
            .not. valid_knot_vector(UP, np+1, p) .or. .not. all(ieee_is_finite(Pw)) .or. &
            .not. ieee_is_finite(u)) then
            allocate(UQ(0:-1), Qw(0:-1,d))
            if (present(T)) allocate(T(0:-1,0:-1))
            return
        end if

        if (r == 0) then
            nq = np
            allocate(UQ(0:mp), Qw(0:np,d))
            UQ = UP
            Qw = Pw
            if (present(T)) then
                allocate(T(0:np,0:np), source=0.0_rk)
                do concurrent (i = 0:np)
                    T(i,i) = 1.0_rk
                end do
            end if
            return
        end if

        ! Insertion follows the last stored occurrence, while evaluation uses
        ! the interior trace at a closed endpoint with exterior or repeated knots.
        span = k
        do while (span < mp)
            if (UP(span+1)>u) exit
            span = span+1
        end do
        if (UP(span)>u .or. span-s>np) then
            allocate(UQ(0:-1),Qw(0:-1,d))
            if (present(T)) allocate(T(0:-1,0:-1))
            return
        end if

        if (s+r == p+1) then
            if (r > 1) then
                if (present(T)) then
                    call insert_knot_A_5_1(&
                        p  = p,&
                        UP = UP,&
                        Pw = Pw,&
                        u  = u,&
                        k  = span,&
                        s  = s,&
                        r  = r-1,&
                        nq = nq_temp,&
                        UQ = Utemp,&
                        Qw = Qtemp,&
                        T  = Ttemp)
                else
                    call insert_knot_A_5_1(&
                        p  = p,&
                        UP = UP,&
                        Pw = Pw,&
                        u  = u,&
                        k  = span,&
                        s  = s,&
                        r  = r-1,&
                        nq = nq_temp,&
                        UQ = Utemp,&
                        Qw = Qtemp)
                end if
            else
                nq_temp = np
                allocate(Utemp(0:mp), Qtemp(0:np,d))
                Utemp = UP
                Qtemp = Pw
                if (present(T)) then
                    allocate(Ttemp(0:np,0:np), source=0.0_rk)
                    do concurrent (i = 0:np)
                        Ttemp(i,i) = 1.0_rk
                    end do
                end if
            end if

            nq = nq_temp + 1
            k_full = span + r - 1
            shared = k_full - p
            allocate(UQ(0:mp+r), Qw(0:nq,d))
            UQ(0:k_full) = Utemp(0:k_full)
            UQ(k_full+1) = u
            UQ(k_full+2:mp+r) = Utemp(k_full+1:mp+r-1)
            Qw(0:shared,:) = Qtemp(0:shared,:)
            Qw(shared+1,:) = Qtemp(shared,:)
            Qw(shared+2:nq,:) = Qtemp(shared+1:nq-1,:)
            if (present(T)) then
                allocate(T(0:nq,0:np))
                T(0:shared,:) = Ttemp(0:shared,:)
                T(shared+1,:) = Ttemp(shared,:)
                T(shared+2:nq,:) = Ttemp(shared+1:nq-1,:)
            end if
            return
        end if

        nq = np+r

        allocate(UQ(0:mp+r))
        allocate(Qw(0:nq, 1:d))
        allocate(Rw(0:p , 1:d))

        UQ(0:span)           = UP(0:span)
        UQ(span+1:span+r)       = u
        UQ(span+1+r:mp+r)    = UP(span+1:mp)
        Qw(0:span-p, :)      = Pw(0:span-p, :)
        Qw(span-s+r:np+r, :) = Pw(span-s:np, :)
        Rw(0:p-s, :)      = Pw(span-p:span-s, :)
        if (present(T)) then
            allocate(T(0:nq, 0:np), source=0.0_rk)
            do concurrent (i = 0:span-p)
                T(i, i) = 1.0_rk
            end do
            do concurrent (i = span-s+r:nq)
                T(i, i-r) = 1.0_rk
            end do
            Lf = span-p+r
            bw = p-s
            do concurrent (i = 0:bw)
                T(Lf+i, span-p+i) = 1.0_rk
            end do
        end if
        do j = 1, r
            L = span-p+j
            do i = 0, p-j-s
                alpha = (u-UP(L+i))/(UP(i+span+1)-UP(L+i))
                Rw(i, :) = alpha*Rw(i+1, :)+(1.0_rk-alpha)*Rw(i, :)
                if (present(T)) then
                    T(Lf+i, :) = alpha*T(Lf+i+1, :)+(1.0_rk-alpha)*T(Lf+i, :)
                end if
            end do
            Qw(L, :)       = Rw(0, :)
            Qw(span+r-j-s, :) = Rw(p-j-s, :)
            if (present(T)) then
                T(L, :)       = T(Lf+0, :)
                T(span+r-j-s, :) = T(Lf+(p-j-s), :)
            end if
        end do
        Qw(L+1:span-s-1, :) = Rw(1:span-s-1-L, :)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Locate the zero-based knot span containing a parameter.
    !!
    !! Binary search follows Piegl and Tiller, Algorithm A2.1. Here `n` is the
    !! highest zero-based control index. Outside parameters saturate to the
    !! first or last positive active span. The left trace chooses the last
    !! start strictly below the knot, the right trace the last start at or below.
    !! Nonfinite parameters or invalid metadata return zero; otherwise `degree<=s<=n`.
    pure function findspan(n,degree,Xth,knot,side) result(s)
        integer, intent(in) :: n, degree
        real(rk), intent(in) :: Xth
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in), optional :: side
            !! -1 for the left trace; omitted, 0 or +1 for the right trace.
        integer :: s
        integer :: low, high, mid
        logical :: left
        real(rk) :: u

        s = 0
        if (degree < 0 .or. n < degree .or. size(knot) < n+degree+2 .or. .not. ieee_is_finite(Xth)) return

        u = max(knot(degree+1),min(Xth,knot(n+2)))
        left = u == knot(n+2)
        if (present(side)) left = left .or. side < 0
        if (u == knot(degree+1)) left = .false.
        low = degree
        high = n + 1
        do while (high-low>1)
            mid = (low+high)/2
            if (u < knot(mid+1) .or. (left .and. u == knot(mid+1))) then
                high = mid
            else
                low = mid
            end if
        end do
        s = low
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Elevate a spline degree while preserving bounded or periodic topology.
    !!
    !! Homogeneous coordinates must be supplied for rational geometry. A
    !! representation satisfying [[periodic_topology]] is elevated cyclically:
    !! each seam and interior multiplicity increases by `t`, the knot extension
    !! remains periodic, and the final elevated-degree control rows repeat the
    !! first rows. Other valid inputs are elevated on their active interval; a
    !! non-clamped bounded representation is first restricted to an equivalent
    !! open-clamped representation.
    !!
    !! `Tmap`, when requested, satisfies `Xcw_new=matmul(Tmap,Xcw)` for either
    !! topology. `t=0` is the exact identity. Invalid input returns empty arrays
    !! and `success=.false.` when present.
    !!
    !! **References:**
    !!
    !! - H. Prautzsch, "Degree Elevation of B-Spline Curves,"
    !!   *Computer Aided Geometric Design* 1 (1984), 193--198.
    !!   [doi:10.1016/0167-8396(84)90031-1](https://doi.org/10.1016/0167-8396(84)90031-1).
    !! - L. Piegl and W. Tiller, *The NURBS Book*, 2nd ed., Springer, 1997,
    !!   Algorithm A5.9.
    !!
    pure subroutine elevate_degree_A_5_9(t, knot, degree, Xcw, nc_new, knot_new, Xcw_new, Tmap, success)
        integer,  intent(in)             :: t
            !! Nonnegative number of degree elevations.
        real(rk), intent(in), contiguous :: Xcw(:,:)
            !! Original control variables, one point per row.
        real(rk), intent(in), contiguous :: knot(:)
            !! Original valid knot vector.
        integer,  intent(in)             :: degree
            !! Original polynomial degree.
        integer,  intent(out)            :: nc_new
            !! New control-point count.
        real(rk), allocatable, intent(out) :: Xcw_new(:,:)
            !! Elevated control variables.
        real(rk), allocatable, intent(out) :: knot_new(:)
            !! Elevated knot vector.
        real(rk), allocatable, intent(out), optional :: Tmap(:,:)
            !! Optional transformation satisfying `Xcw_new=matmul(Tmap,Xcw)`.
        logical, intent(out), optional :: success
            !! Optional operation status.

        if (t > 0 .and. periodic_topology(knot, degree, Xcw)) then
            call elevate_degree_periodic_A_5_9(t,knot,degree,Xcw,nc_new,knot_new,Xcw_new,Tmap,success)
        else
            call elevate_degree_bounded_A_5_9(t,knot,degree,Xcw,nc_new,knot_new,Xcw_new,Tmap,success)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> Elevate a bounded spline after periodic topology has been excluded.
    !!
    !! Open-clamped input is elevated directly. Any other valid representation
    !! is restricted to an equivalent open-clamped representation on its active
    !! interval before Algorithm A5.9 is applied.
    pure subroutine elevate_degree_bounded_A_5_9(t, knot, degree, Xcw, nc_new, knot_new, Xcw_new, Tmap, success)
        integer, intent(in) :: t, degree
        real(rk), intent(in), contiguous :: Xcw(:,:), knot(:)
        integer, intent(out) :: nc_new
        real(rk), allocatable, intent(out) :: Xcw_new(:,:), knot_new(:)
        real(rk), allocatable, intent(out), optional :: Tmap(:,:)
        logical, intent(out), optional :: success
        real(rk), allocatable :: clamped(:,:), clamped_knots(:), Telev(:,:), composed(:,:)
        integer :: nc
        logical :: open_knot, ok

        if (present(success)) success = .false.
        nc = size(Xcw,1)
        nc_new = 0
        if (t<0 .or. size(Xcw,2)<1 .or. .not.valid_knot_vector(knot,nc,degree) .or. &
            .not.all(ieee_is_finite(Xcw))) then
            allocate(knot_new(0),Xcw_new(0,size(Xcw,2)))
            if (present(Tmap)) Tmap = eye(0)
            return
        end if
        if (t==0) then
            nc_new = nc
            knot_new = knot
            Xcw_new = Xcw
            if (present(Tmap)) Tmap = eye(nc)
            if (present(success)) success = .true.
            return
        end if
        open_knot = all(knot(1:degree+1)==knot(degree+1)) .and. &
            all(knot(nc+1:nc+degree+1)==knot(nc+1))
        if (open_knot) then
            call elevate_degree_open_A_5_9(t,knot,degree,Xcw,nc_new,knot_new,Xcw_new,Tmap)
        else
            call restrict_spline(knot,degree,Xcw,knot(degree+1),knot(nc+1), &
                clamped_knots,clamped,Tmap,ok)
            if (.not.ok) then
                allocate(knot_new(0),Xcw_new(0,size(Xcw,2)))
                if (present(Tmap)) Tmap = eye(0)
                return
            end if
            if (present(Tmap)) then
                call elevate_degree_open_A_5_9(t,clamped_knots,degree,clamped,nc_new,knot_new,Xcw_new,Telev)
                call sparse_left_matmul(Telev,Tmap,composed)
                call move_alloc(composed,Tmap)
            else
                call elevate_degree_open_A_5_9(t,clamped_knots,degree,clamped,nc_new,knot_new,Xcw_new)
            end if
        end if
        if (present(success)) success = nc_new>0 .and. size(knot_new)==nc_new+degree+t+1 .and. &
            size(Xcw_new,1)==nc_new .and. all(ieee_is_finite(Xcw_new))
    end subroutine elevate_degree_bounded_A_5_9
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Elevate an open B-spline representation using Algorithm A5.9.
    !!
    !! The routine raises degree `degree` by `t`, preserves every geometric or
    !! homogeneous coordinate column of `Xcw`, and optionally returns the dense
    !! linear operator satisfying `Xcw_new=matmul(Tmap,Xcw)`. The input is the
    !! open-clamped workspace prepared by [[elevate_degree_A_5_9]].
    !!
    !! **Reference:** L. Piegl and W. Tiller, *The NURBS Book*, 2nd ed.,
    !! Springer, 1997, Algorithm A5.9.
    pure subroutine elevate_degree_open_A_5_9(t, knot, degree, Xcw, nc_new, knot_new, Xcw_new, Tmap)
        integer, intent(in) :: t
            !! Nonnegative degree increment.
        real(rk), intent(in), contiguous :: Xcw(:,:)
            !! Control or homogeneous-control rows `[nc,nvalue]`.
        real(rk), intent(in), contiguous :: knot(:)
            !! Open-clamped knot vector of degree `degree`.
        integer, intent(in) :: degree
            !! Original polynomial degree.
        integer, intent(out) :: nc_new
            !! Elevated control count.
        real(rk), allocatable, intent(out) :: Xcw_new(:,:)
            !! Elevated control or homogeneous-control rows.
        real(rk), allocatable, intent(out) :: knot_new(:)
            !! Elevated knot vector; each distinct multiplicity increases by `t`.
        real(rk), allocatable, intent(out), optional :: Tmap(:,:)
            !! Optional elevation operator `[nc_new,nc]`.

        real(rk), allocatable :: bezalfs(:,:), bpts(:,:), ebpts(:,:), Nextbpts(:,:), alfs(:)
        real(rk), allocatable :: bC(:,:), ebC(:,:), NextbC(:,:)  ! only used if Tmap present

        real(rk) :: alpha1, alpha2, Xth1, Xth2, numer, den, alpha3
        integer  :: n, lbz, rbz, sv, tr, kj, first, knoti, last, d, nc
        integer  :: i, j, q, s, m, ph, ph2, mpi, mh, r, a, b, Xcwi, oldr, mul, lo
        integer, allocatable :: mlp(:)

        nc = size(Xcw,1)
        d  = size(Xcw,2)

        ! A constant piece is a Bezier segment with every elevated control equal.
        if (degree == 0) then
            nc_new = nc*(t+1)
            allocate(knot_new(size(knot)*(t+1)), Xcw_new(nc_new,d))
            if (present(Tmap)) allocate(Tmap(nc_new,nc),source=0.0_rk)
            do i = 1, size(knot)
                knot_new((i-1)*(t+1)+1:i*(t+1)) = knot(i)
            end do
            do i = 1, nc
                do j = (i-1)*(t+1)+1, i*(t+1)
                    Xcw_new(j,:) = Xcw(i,:)
                    if (present(Tmap)) Tmap(j,i) = 1.0_rk
                end do
            end do
            return
        end if

        mlp = compute_multiplicity(knot)
        mlp = mlp + t
        nc_new = sum(mlp) - (mlp(1)-1) - 1

        allocate(Xcw_new(nc_new,d), source=0.0_rk)
        allocate(bezalfs(degree+1,degree+t+1), source=0.0_rk)
        allocate(bpts(degree+1,d), source=0.0_rk)
        allocate(ebpts(degree+t+1,d), source=0.0_rk)
        allocate(Nextbpts(degree+1,d), source=0.0_rk)
        allocate(alfs(degree), source=0.0_rk)

        if (present(Tmap)) then
            allocate(Tmap(nc_new, nc), source=0.0_rk)
            allocate(bC(degree+1,nc), ebC(degree+t+1,nc), NextbC(degree+1,nc), source=0.0_rk)
        end if

        n   = nc - 1
        m   = n + degree + 1
        ph  = degree + t
        ph2 = ph/2

        bezalfs(1,1)           = 1.0_rk
        bezalfs(degree+1,ph+1) = 1.0_rk
        do i = 1, ph2
            lo = max(0,i-t)
            mpi = min(degree,i)
            bezalfs(lo+1,i+1) = 1.0_rk
            if (i <= t) then
                do q = 0, i-1
                    bezalfs(lo+1,i+1) = bezalfs(lo+1,i+1)*real(t-q,rk)/real(ph-q,rk)
                end do
            else
                do q = 1, t
                    bezalfs(lo+1,i+1) = bezalfs(lo+1,i+1)*real(lo+q,rk)/real(degree+q,rk)
                end do
            end if
            do j = lo, mpi-1
                bezalfs(j+2,i+1) = bezalfs(j+1,i+1)*real(degree-j,rk)*real(i-j,rk) / &
                    (real(j+1,rk)*real(t-i+j+1,rk))
            end do
        end do
        do i = ph2+1, ph-1
            mpi = min(degree,i)
            do j = max(0,i-t), mpi
                bezalfs(j+1,i+1) = bezalfs(degree-j+1, ph-i+1)
            end do
        end do

        mh    = ph
        knoti = ph + 1
        r     = -1
        a     = degree
        b     = degree + 1
        Xcwi  = 1

        Xth1 = knot(1)
        Xcw_new(1,:) = Xcw(1,:)
        if (present(Tmap)) Tmap(1,1) = 1.0_rk

        allocate(knot_new(sum(mlp)), source=0.0_rk)
        knot_new(1:ph+1) = Xth1

        do i = 0, degree
            bpts(i+1,:) = Xcw(i+1,:)
            if (present(Tmap)) bC(i+1, i+1) = 1.0_rk
        end do

        do while (b < m)
            i = b
            do while (b < m .and. knot(b+1) == knot(b+2))
                b = b + 1
                if (b+2 > size(knot)) exit
            end do
            mul  = b - i + 1
            mh   = mh + mul + t
            Xth2 = knot(b+1)
            oldr = r
            r    = degree - mul

            if (oldr < 0 .and. a /= degree) then
                ! A full break shares no endpoint control point with the
                ! preceding Bezier piece, so retain this piece's first point.
                lbz = 0
            else
                lbz = merge((oldr+2)/2, 1, oldr > 0)
            end if
            rbz = merge(ph - (r+1)/2, ph, r > 0)

            if (r > 0) then
                numer = Xth2 - Xth1
                do q = degree, mul+1, -1
                    alfs(q-mul) = numer / (knot(a+q+1) - Xth1)
                end do
                do j = 1, r
                    sv = r - j;  s = mul + j
                    do q = degree, s, -1
                        bpts(q+1,:) = (1.0_rk-alfs(q-s+1))*bpts(q,:) + alfs(q-s+1)*bpts(q+1,:)
                        if (present(Tmap)) bC(q+1,:) = (1.0_rk-alfs(q-s+1))*bC(q,:) + alfs(q-s+1)*bC(q+1,:)
                    end do
                    Nextbpts(sv+1,:) = bpts(degree+1,:)
                    if (present(Tmap)) NextbC(sv+1,:) = bC(degree+1,:)
                end do
            end if

            do i = lbz, ph
                ebpts(i+1,:) = 0.0_rk
                if (present(Tmap)) ebC(i+1,:) = 0.0_rk
                mpi = min(degree,i)
                do j = max(0,i-t), mpi
                    ebpts(i+1,:) = ebpts(i+1,:) + bezalfs(j+1,i+1)*bpts(j+1,:)
                    if (present(Tmap)) ebC(i+1,:) = ebC(i+1,:) + bezalfs(j+1,i+1)*bC(j+1,:)
                end do
            end do

            if (oldr > 1) then
                first = knoti - 2
                last  = knoti
                den   = Xth2 - Xth1
                if (abs(den) <= 32.0_rk*epsilon(1.0_rk)*max(abs(Xth1), abs(Xth2))) then
                    alpha3 = 0.0_rk
                else
                    alpha3 = (Xth2 - knot_new(knoti)) / den
                end if

                do tr = 1, oldr-1
                    i  = first
                    j  = last
                    kj = j - knoti + 1
                    do while (j - i > tr)
                        if (i < Xcwi) then
                            alpha1 = (Xth2 - knot_new(i+1)) / (Xth1 - knot_new(i+1))
                            Xcw_new(i+1,:) = (1.0_rk-alpha1)*Xcw_new(i,:) + alpha1*Xcw_new(i+1,:)
                            if (present(Tmap)) Tmap(i+1,:) = (1.0_rk-alpha1)*Tmap(i,:) + alpha1*Tmap(i+1,:)
                        end if
                        if (j >= lbz) then
                            if (j-tr <= knoti - ph + oldr) then
                                alpha2 = (Xth2 - knot_new(j-tr+1)) / den
                                ebpts(kj+1,:) = alpha2*ebpts(kj+1,:) + (1.0_rk-alpha2)*ebpts(kj+2,:)
                                if (present(Tmap)) ebC(kj+1,:) = alpha2*ebC(kj+1,:) + (1.0_rk-alpha2)*ebC(kj+2,:)
                            else
                                ebpts(kj+1,:) = (1.0_rk-alpha3)*ebpts(kj+2,:) + alpha3*ebpts(kj+1,:)
                                if (present(Tmap)) ebC(kj+1,:) = (1.0_rk-alpha3)*ebC(kj+2,:) + alpha3*ebC(kj+1,:)
                            end if
                        end if
                        i  = i + 1;  j = j - 1;  kj = kj - 1
                    end do
                    first = first - 1
                    last  = last + 1
                end do
            end if

            if (a /= degree) then
                do i = 0, ph-oldr-1
                    knot_new(knoti+1) = Xth1
                    knoti = knoti + 1
                end do
            end if

            do j = lbz, rbz
                Xcw_new(Xcwi+1,:) = ebpts(j+1,:)
                if (present(Tmap)) Tmap(Xcwi+1,:) = ebC(j+1,:)
                Xcwi = Xcwi + 1
            end do

            if (b < m) then
                do j = 0, r-1
                    bpts(j+1,:) = Nextbpts(j+1,:)
                    if (present(Tmap)) bC(j+1,:) = NextbC(j+1,:)
                end do
                ! A knot of multiplicity degree+1 separates two independent
                ! spline pieces and gives r=-1.  The next piece starts at its
                ! first control point; there is no control point at index zero.
                do j = max(0, r), degree
                    bpts(j+1,:) = Xcw(b-degree+j+1,:)
                    if (present(Tmap)) then
                        bC(j+1,:)  = 0.0_rk
                        bC(j+1, b-degree+j+1) = 1.0_rk
                    end if
                end do
                a    = b
                b    = b + 1
                Xth1 = Xth2
            else
                do i = 0, ph
                    knot_new(knoti+i+1) = Xth2
                end do
            end if
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Generate a uniformly spaced axis-aligned hexahedral control net.
    !!
    !! `L` and `nc` must each have exactly three entries and every count must be
    !! greater than one; these shape preconditions are not checked. Coordinates
    !! run from zero to the corresponding signed extent, so negative extents
    !! reverse an axis and a zero extent produces a degenerate net.
    !!
    pure function hexahedron_Xc(L, nc) result(Xc)
        real(rk), intent(in), contiguous :: L(:)
            !! Signed axis extents `[Lx,Ly,Lz]`.
        integer, intent(in), contiguous :: nc(:)
            !! Directional point counts `[nx,ny,nz]`, each greater than one. The result has shape `[product(nc),3]`,
            !! with x/direction 1 fastest.
        real(rk), allocatable :: Xc(:,:)
        real(rk) :: d(3)
        integer :: i, j, k, idx

        d = L/real(nc-1, rk)

        allocate(Xc(nc(1) * nc(2) * nc(3), 3))
        do concurrent (k = 0:nc(3)-1, j = 0:nc(2)-1, i = 0:nc(1)-1) local(idx)
            idx = i+nc(1)*(j+k*nc(2))+1
            Xc(idx, 1) = real(i, rk)*d(1)
            Xc(idx, 2) = real(j, rk)*d(2)
            Xc(idx, 3) = real(k, rk)*d(3)
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Generate a uniformly spaced planar rectangular control net.
    !!
    !! `L` and `nc` must each have exactly two entries and every count must be
    !! greater than one; these shape preconditions are not checked. Coordinates
    !! run from zero to each signed extent and the third coordinate is zero.
    !!
    pure function tetragon_Xc(L, nc) result(Xc)
        real(rk), intent(in), contiguous :: L(:)
            !! Signed axis extents `[Lx,Ly]`.
        integer, intent(in), contiguous :: nc(:)
            !! Directional point counts `[nx,ny]`, each greater than one. The result has shape `[product(nc),3]`, lies
            !! in `z=0`, and stores x/direction 1 fastest.
        real(rk), allocatable :: Xc(:,:)
        real(rk) :: d(2)
        integer :: i, j, idx

        d = L/real(nc-1, rk)

        allocate(Xc(nc(1) * nc(2), 3))
        do concurrent (j = 0:nc(2)-1, i = 0:nc(1)-1) local(idx)
            idx = i+j*nc(1)+1
            Xc(idx, 1) = real(i, rk)*d(1)
            Xc(idx, 2) = real(j, rk)*d(2)
            Xc(idx, 3) = 0.0_rk
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Attempt knot removal using Algorithm A5.8 of Piegl and Tiller.
    !! `r` must be the one-based index of the final occurrence of the interior
    !! knot `u`, `s` must be its complete current multiplicity, and
    !! `1<=num<=s`. Each removal is accepted by componentwise comparisons in the
    !! supplied control-variable space using
    !!
    !! \[
    !! |a-b|\le128\,\epsilon_{rk}\max(|a|,|b|,C_c),
    !! \]
    !!
    !! where C_c is the largest magnitude among the affected control
    !! coefficients in that coordinate.
    !! Polynomial callers pass control coordinates. Rational callers must pass
    !! homogeneous controls \((w\mathbf P,w)\); the resulting test is then in
    !! homogeneous coordinate space and is not a certified Euclidean geometry
    !! error bound. The two acceptance branches are mutually exclusive:
    !! `j-i<t` compares the opposing temporary controls, whereas `j-i==t`
    !! compares the current control with its reconstructed value. If no removal
    !! is accepted, `t=0` and both outputs are copies of the inputs.
    !!
    pure subroutine remove_knots_A_5_8(p, knot, Pw, u, r, s, num, t, knot_new, Pw_new, approximate)
        integer, intent(in) :: p, r, s, num
        real(rk), intent(in) :: u
        real(rk), intent(in), contiguous :: knot(:), Pw(:,:)
        logical, optional, intent(in) :: approximate
            !! Candidate-only mode; callers must separately bound physical error.
        integer, intent(out) :: t
        real(rk), allocatable, intent(out) :: knot_new(:), Pw_new(:,:)

        if (periodic_topology(knot,p,Pw)) then
            call remove_knots_periodic_A_5_8(p,knot,Pw,u,r,s,num,t,knot_new,Pw_new,approximate)
        else
            call remove_knots_bounded_A_5_8(p,knot,Pw,u,r,s,num,t,knot_new,Pw_new,approximate)
        end if
    end subroutine remove_knots_A_5_8

    !> Remove complete cyclic knot orbits using the bounded coefficient test.
    !! Guard periods keep every extracted control support away from artificial
    !! ends. Each accepted orbit removes one copy per period; a failed orbit
    !! leaves the previously accepted prefix intact. The active endpoints stay fixed.
    pure subroutine remove_knots_periodic_A_5_8(p, knot, Pw, u, r, s, num, t, knot_new, Pw_new, approximate)
        integer, intent(in) :: p, r, s, num
        real(rk), intent(in) :: u
        real(rk), intent(in), contiguous :: knot(:), Pw(:,:)
        logical, optional, intent(in) :: approximate
        integer, intent(out) :: t
        real(rk), allocatable, intent(out) :: knot_new(:), Pw_new(:,:)
        real(rk), allocatable :: knots(:), controls(:,:), trial_knots(:), trial_controls(:,:)
        real(rk), allocatable :: next_knots(:), next_controls(:,:), candidate_knots(:), candidate_controls(:,:)
        real(rk) :: a, b, period, translated
        integer :: nc, unique_count, target_count, guards, extended_count, i, base, shift, index_, attempt, accepted
        integer :: last, multiplicity, first_start, first_end, first_knot, candidate_count

        t=0
        knot_new=knot
        Pw_new=Pw
        nc=size(Pw,1)
        if (r<1 .or. r>size(knot) .or. num<1 .or. num>s .or. s>p) return
        if (knot(r)/=u .or. count(knot==u)/=s) return
        if (r<size(knot)) then
            if (knot(r+1)==u) return
        end if
        a=knot(p+1)
        b=knot(nc+1)
        if (u<=a .or. u>=b) return
        period=b-a
        if (.not. ieee_is_finite(period)) return
        unique_count=nc-p
        target_count=unique_count-num
        if (target_count<1) return
        guards=(p+target_count-1)/target_count+1
        extended_count=(2*guards+1)*unique_count+p
        allocate(knots(extended_count+p+1),controls(extended_count,size(Pw,2)))
        do i=1,extended_count
            controls(i,:)=Pw(mod(i-1,unique_count)+1,:)
        end do
        do i=1,size(knots)
            index_=i-1-guards*unique_count
            base=p+modulo(index_-p,unique_count)
            shift=(index_-base)/unique_count
            knots(i)=knot(base+1)+real(shift,rk)*period
        end do
        if (.not. valid_knot_vector(knots,extended_count,p)) return

        removal_loop: do attempt=1,num
            trial_knots=knots
            trial_controls=controls
            do shift=-guards,guards
                translated=u+real(shift,rk)*period
                multiplicity=count(trial_knots==translated)
                last=findspan(size(trial_controls,1)-1,p,translated,trial_knots)+1
                call remove_knots_bounded_A_5_8(p,trial_knots,trial_controls,translated,last,multiplicity,1, &
                    accepted,next_knots,next_controls,approximate)
                if (accepted/=1) exit removal_loop
                call move_alloc(next_knots,trial_knots)
                call move_alloc(next_controls,trial_controls)
            end do
            first_start=findloc(trial_knots,a,dim=1)
            first_end=findloc(trial_knots,b,dim=1)
            first_knot=first_start-p
            candidate_count=first_end-first_start+p
            if (first_knot<1 .or. first_end+p>size(trial_knots) .or. candidate_count<=p) exit removal_loop
            candidate_knots=trial_knots(first_knot:first_end+p)
            candidate_controls=trial_controls(first_knot:first_knot+candidate_count-1,:)
            if (.not. periodic_topology(candidate_knots,p,candidate_controls)) exit removal_loop
            candidate_controls(candidate_count-p+1:,:)=candidate_controls(:p,:)
            call move_alloc(candidate_knots,knot_new)
            call move_alloc(candidate_controls,Pw_new)
            call move_alloc(trial_knots,knots)
            call move_alloc(trial_controls,controls)
            t=t+1
        end do removal_loop
    end subroutine remove_knots_periodic_A_5_8

    !> Shared componentwise roundoff policy for both A5.8 consistency branches.
    pure logical function removal_coefficients_equal(a,b,magnitude) result(equal)
        real(rk), intent(in) :: a(:), b(:), magnitude(:)
        equal=all(abs(a-b)<=128.0_rk*epsilon(1.0_rk)*max(abs(a),abs(b),magnitude))
    end function removal_coefficients_equal

    pure recursive subroutine remove_knots_bounded_A_5_8(p,knot,Pw,u,r,s,num,t,knot_new,Pw_new,approximate)
        integer, intent(in) :: p, r, s, num
        real(rk), intent(in) :: u
        real(rk), contiguous, intent(in) :: knot(:), Pw(:,:)
        logical, optional, intent(in) :: approximate
        integer, intent(out) :: t
        real(rk), allocatable, intent(out) :: knot_new(:), Pw_new(:,:)
        real(rk), allocatable :: Pw_copy(:,:), knot_copy(:)
        real(rk) :: alfi, alfj, denominator_i, denominator_j, candidate(size(Pw,2))
        real(rk) :: coefficient_scale(size(Pw,2))
        real(rk), allocatable :: temp(:,:)
        logical :: candidate_only
        integer :: i, j, ii, jj, remflag, off, first, last, ord, fout, m, k, n, nc, d, tt, t_more

        d  = size(Pw,2)
        nc = size(Pw,1)
        t = 0
        candidate_only=.false.
        if (present(approximate)) candidate_only=approximate
        knot_new = knot
        Pw_new = Pw

        if (p < 0 .or. nc <= p .or. size(knot) /= nc+p+1 .or. num < 1 .or. s < 1 .or. num > s) return
        if (.not. valid_knot_vector(knot, nc, p)) return
        if (r < 1 .or. r > size(knot) .or. s > p+1 .or. knot(r) /= u) return
        if (r-s+1 < 1 .or. any(knot(r-s+1:r) /= u)) return
        if (r-s >= 1) then
            if (knot(r-s) == u) return
        end if
        if (r < size(knot)) then
            if (knot(r+1) == u) return
        end if
        if (u <= knot(p+1) .or. u >= knot(nc+1)) return

        if (s == p+1) then
            first = r-p
            if (first-1 < 1 .or. first > nc) return

            if (.not.candidate_only) then
                if (.not.removal_coefficients_equal(Pw(first-1,:),Pw(first,:),abs(Pw(first,:)))) return
            end if

            deallocate(Pw_new, knot_new)
            allocate(Pw_new(nc-1,d))
            Pw_new(1:first-2,:) = Pw(1:first-2,:)
            Pw_new(first-1,:) = 0.5_rk*Pw(first-1,:) + 0.5_rk*Pw(first,:)
            Pw_new(first:nc-1,:) = Pw(first+1:nc,:)

            allocate(knot_new(size(knot)-1))
            knot_new(1:r-1) = knot(1:r-1)
            knot_new(r:size(knot)-1) = knot(r+1:size(knot))
            t = 1

            if (num > 1) then
                call remove_knots_bounded_A_5_8(p,knot_new,Pw_new,u,r-1,s-1,num-1,t_more,knot_copy,Pw_copy,approximate)
                if (t_more > 0) then
                    t = t + t_more
                    call move_alloc(knot_copy, knot_new)
                    call move_alloc(Pw_copy, Pw_new)
                end if
            end if
            return
        end if

        n = nc
        m = n+p+1
        ord = p+1
        fout = (2*r-s-p)/2
        last = r-s
        first = r-p

        if (first-1 < 1 .or. last+1 > nc) return

        Pw_copy = Pw
        knot_copy = knot

        allocate(temp(2*p+1,d), source=0.0_rk)
        removal_loop: do tt = 0, num-1
            off = first-1
            coefficient_scale=maxval(abs(Pw_copy(off:last+1,:)),dim=1)
            temp(1,:) = Pw_copy(off,:)
            temp(last+1-off+1,:) = Pw_copy(last+1,:)
            i = first
            j = last
            ii = 1
            jj = last-off
            remflag = 0
            do while (j-i>t)
                denominator_i = knot_copy(i+ord+t) - knot_copy(i)
                denominator_j = knot_copy(j+ord) - knot_copy(j-t)
                if (denominator_i == 0.0_rk .or. denominator_j == 0.0_rk) exit removal_loop
                alfi = (u-knot_copy(i))/denominator_i
                alfj = (u-knot_copy(j-t))/denominator_j
                if (alfi == 0.0_rk .or. alfj == 1.0_rk) exit removal_loop
                temp(ii+1,:) = (Pw_copy(i,:)-(1.0_rk-alfi)*temp(ii-1+1,:))/alfi
                temp(jj+1,:) = (Pw_copy(j,:)-alfj*temp(jj+1+1,:))/(1.0_rk-alfj)
                i = i+1
                ii = ii+1
                j = j-1
                jj = jj-1
            end do
            if (j-i<t) then
                remflag=merge(1,0,candidate_only.or.removal_coefficients_equal(temp(ii,:),temp(jj+2,:),coefficient_scale))
            else if (j-i==t) then
                denominator_i = knot_copy(i+ord+t) - knot_copy(i)
                if (denominator_i /= 0.0_rk) then
                    alfi = (u-knot_copy(i))/denominator_i
                    candidate=alfi*temp(ii+t+2,:) + (1.0_rk-alfi)*temp(ii,:)
                    remflag=merge(1,0,candidate_only.or.removal_coefficients_equal(Pw_copy(i,:),candidate,coefficient_scale))
                end if
            end if
            if (remflag == 0) then
                exit removal_loop
            else
                i = first
                j = last
                do while(j-i>t)
                    Pw_copy(i,:) = temp(i-off+1,:)
                    Pw_copy(j,:) = temp(j-off+1,:)
                    i = i+1
                    j = j-1
                end do
            end if
            first = first-1
            last = last+1
            t=t+1
        end do removal_loop
        if (t==0) then
            return
        end if
        do k = r+1, m
            knot_copy(k-t) = knot_copy(k)
        end do
        j = fout
        i = j
        do k = 1, t-1
            if (mod(k,2)==1) then
                i = i+1
            else
                j = j-1
            end if
        end do
        do k = i+1, n
            Pw_copy(j,:) = Pw_copy(k,:)
            j = j+1
        end do
        knot_new = knot_copy(1:size(knot_copy)-t)
        Pw_new = Pw_copy(1:size(Pw_copy,1)-t,:)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return distinct integer values in stable first-occurrence order.
    !!
    !! Equality is exact. The implementation performs \(O(n^2)\) comparisons
    !! and uses \(O(n)\) temporary storage, so it is intended for short metadata
    !! and connectivity vectors rather than large unsorted datasets.
    pure function unique_integer(vec) result(output)
        integer, intent(in), contiguous :: vec(:)
            !! Input sequence; it is not sorted or modified.
        integer, allocatable :: output(:)
            !! Stable sequence of distinct values.
        integer, allocatable :: tmp(:)
        integer :: i, j, nunique
        logical :: found

        allocate(tmp(size(vec)))
        nunique = 0
        do i = 1, size(vec)
            found = .false.
            do j = 1, nunique
                if (vec(i) == tmp(j)) then
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                nunique = nunique + 1
                tmp(nunique) = vec(i)
            end if
        end do

        allocate(output(nunique))
        if (nunique > 0) output = tmp(:nunique)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return numerically distinct real values in stable first-occurrence order.
    !!
    !! Two entries are treated as equal when their absolute difference is at
    !! most `32*epsilon(1.0_rk)*maxval(abs(vec))`. The comparison is therefore
    !! scale-relative but not transitive clustering. Complexity is \(O(n^2)\).
    pure function unique_real(vec) result(output)
        real(rk), intent(in), contiguous :: vec(:)
            !! Input sequence; an empty input returns an empty result.
        real(rk), allocatable :: output(:)
            !! Stable sequence of tolerance-distinct values.
        real(rk), allocatable :: tmp(:)
        integer :: i, j, nunique
        real(rk) :: tol
        logical :: found

        if (size(vec) == 0) then
            allocate(output(0))
            return
        end if

        allocate(tmp(size(vec)))
        tol = 32.0_rk*epsilon(1.0_rk)*maxval(abs(vec))
        nunique = 0
        do i = 1, size(vec)
            found = .false.
            do j = 1, nunique
                if (abs(vec(i) - tmp(j)) <= tol) then
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                nunique = nunique + 1
                tmp(nunique) = vec(i)
            end if
        end do

        allocate(output(nunique))
        if (nunique > 0) output = tmp(:nunique)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return an active right-handed Euler rotation matrix.
    !!
    !! Angles are in degrees and
    !! \(R=R_z(\theta)R_y(\beta)R_x(\alpha)\) acts on column vectors. Thus the
    !! active rotations are applied successively about the fixed x, y, and z
    !! axes (the extrinsic xyz convention).
    pure function rotation(alpha, beta, theta) result(R)
        real(rk), intent(in) :: alpha
            !! Rotation about the x axis in degrees.
        real(rk), intent(in) :: beta
            !! Rotation about the y axis in degrees.
        real(rk), intent(in) :: theta
            !! Rotation about the z axis in degrees.
        real(rk) :: R(3,3)

        R(1,1) = cosd(beta)*cosd(theta)
        R(2,1) = cosd(beta)*sind(theta)
        R(3,1) = -sind(beta)
        R(1,2) = sind(alpha)*sind(beta)*cosd(theta) - cosd(alpha)*sind(theta)
        R(2,2) = sind(alpha)*sind(beta)*sind(theta) + cosd(alpha)*cosd(theta)
        R(3,2) = sind(alpha)*cosd(beta)
        R(1,3) = cosd(alpha)*sind(beta)*cosd(theta) + sind(alpha)*sind(theta)
        R(2,3) = cosd(alpha)*sind(beta)*sind(theta) - sind(alpha)*cosd(theta)
        R(3,3) = cosd(alpha)*cosd(beta)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Rotate the first `n` rows of a point array in place.
    !!
    !! The same Euler convention as [[rotation]] is used. Arrays with at least
    !! three coordinates rotate their first three and leave additional
    !! coordinates unchanged. One- and two-coordinate arrays receive only the
    !! leading principal submatrix of the three-dimensional rotation; that map
    !! is not generally an orthogonal lower-dimensional rotation when omitted
    !! coordinates couple to the selected angles. `n` is clipped to the
    !! available rows; a negative `n` changes nothing.
    pure subroutine rotate_points(points, n, alpha, beta, theta)
        real(rk), intent(inout), contiguous :: points(:,:)
            !! Point rows to modify.
        integer,  intent(in) :: n
            !! Maximum number of leading rows to rotate.
        real(rk), intent(in) :: alpha
            !! Rotation about x in degrees.
        real(rk), intent(in) :: beta
            !! Rotation about y in degrees.
        real(rk), intent(in) :: theta
            !! Rotation about z in degrees.
        integer :: i, ncoord, npoint
        real(rk) :: ca, sa, cb, sb, ct, st, r11, r12, r13, r21, r22, r23, r31, r32, r33
        real(rk) :: x, y, z

        npoint = min(n, size(points, 1))
        ncoord = size(points, 2)
        if (npoint < 1 .or. ncoord < 1) return

        ca = cosd(alpha); sa = sind(alpha)
        cb = cosd(beta);  sb = sind(beta)
        ct = cosd(theta); st = sind(theta)
        r11 = cb*ct
        r21 = cb*st
        r31 = -sb
        r12 = sa*sb*ct - ca*st
        r22 = sa*sb*st + ca*ct
        r32 = sa*cb
        r13 = ca*sb*ct + sa*st
        r23 = ca*sb*st - sa*ct
        r33 = ca*cb

        select case (ncoord)
        case (1)
            do concurrent (i = 1:npoint) local(x)
                x = points(i, 1)
                points(i, 1) = r11*x
            end do
        case (2)
            do concurrent (i = 1:npoint) local(x, y)
                x = points(i, 1)
                y = points(i, 2)
                points(i, 1) = r11*x + r12*y
                points(i, 2) = r21*x + r22*y
            end do
        case default
            do concurrent (i = 1:npoint) local(x, y, z)
                x = points(i, 1)
                y = points(i, 2)
                z = points(i, 3)
                points(i, 1) = r11*x + r12*y + r13*z
                points(i, 2) = r21*x + r22*y + r23*z
                points(i, 3) = r31*x + r32*y + r33*z
            end do
        end select
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute a small determinant or a surface-Jacobian area factor.
    !!
    !! Square matrices of order one through three return their determinant. A
    !! `3x2` matrix returns \(\sqrt{\det(A^T A)}\), the magnitude of the cross
    !! product of its columns. Other shapes return quiet NaN. Finite entries
    !! whose intermediate products do not overflow are a caller precondition;
    !! this low-level routine performs no numerical-range diagnostics.
    pure function det(A) result(detA)
        real(rk), intent(in) :: A(:,:)
        real(rk) :: detA

        detA = ieee_value(1.0_rk, ieee_quiet_nan)
        if (size(A,1) == size(A,2)) then
            select case(size(A,1))
              case(1)
                detA = A(1,1)
              case(2)
                detA = A(1,1)*A(2,2) - A(1,2)*A(2,1)
              case(3)
                detA = &
                    + A(1,1)*( A(2,2)*A(3,3) - A(2,3)*A(3,2) )&
                    - A(1,2)*( A(2,1)*A(3,3) - A(2,3)*A(3,1) )&
                    + A(1,3)*( A(2,1)*A(3,2) - A(2,2)*A(3,1) )
              case default
                detA = ieee_value(1.0_rk, ieee_quiet_nan)
            end select
        else if (size(A,1) == 3 .and. size(A,2) == 2) then
            detA = sqrt( &
                (A(2,1)*A(3,2) - A(3,1)*A(2,2))**2 + &
                (A(3,1)*A(1,2) - A(1,1)*A(3,2))**2 + &
                (A(1,1)*A(2,2) - A(2,1)*A(1,2))**2)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute a small inverse or full-rank Moore-Penrose pseudoinverse.
    !!
    !! Let \(e\) be the model exponent of \(\max_{ij}|A_{ij}|\), and define the
    !! exact radix-power scaling
    !! \(s=\operatorname{radix}(A)^e\). Orders one through three use explicit
    !! formulas for \((A/s)^{-1}\); larger square matrices delegate to
    !! [[solve]] and are rescaled to the same representation. Because radix
    !! scaling introduces no rounding, the dimensionless reciprocal
    !! infinity-norm condition number
    !!
    !! \[
    !! r_\infty=
    !! \frac{1}{\|A/s\|_\infty\,\|(A/s)^{-1}\|_\infty}
    !! \]
    !!
    !! must be finite and greater than \(32\epsilon_{rk}\). Thus uniform
    !! nonzero scaling cannot change acceptance, while singular and
    !! numerically singular matrices return `A_inv(0,0)`.
    !!
    !! For a full-column-rank \(m\times n\) matrix with \(m>n\), the
    !! rectangular path computes the equivalent scaled form
    !! \(s^{-1}[(A/s)^T(A/s)]^{-1}(A/s)^T\). For full row rank with \(m<n\),
    !! it computes \(s^{-1}(A/s)^T[(A/s)(A/s)^T]^{-1}\). These are the
    !! Moore-Penrose inverse only under the stated rank assumptions and the
    !! condition test applied recursively to the normal matrix. Empty,
    !! nonfinite, or numerically rejected input returns `A_inv(0,0)`.
    !!
    !! @warning This routine is intended for geometry Jacobians, not general
    !! large least-squares problems. Normal equations square the condition
    !! number. Use a rank-revealing QR or singular-value decomposition for
    !! ill-conditioned rectangular problems.
    !! @endwarning
    !!
    !! **Reference:** R. Penrose, "A Generalized Inverse for Matrices,"
    !! *Proceedings of the Cambridge Philosophical Society* 51 (1955),
    !! 406--413.
    !! [doi:10.1017/S0305004100030401](https://doi.org/10.1017/S0305004100030401).
    !!
    !! **Numerical reference:** LAPACK `DGECON` defines the reciprocal
    !! condition number as
    !! \(1/(\|A\|\|A^{-1}\|)\):
    !! [Netlib LAPACK](https://netlib.org/lapack/explore-html/d4/daf/group__gecon_ga4f9b830e19e12c7f082ddb497a57af18.html).
    recursive pure function inv(A) result(A_inv)
        real(rk), intent(in), contiguous :: A(:,:)
        real(rk), allocatable :: A_inv(:,:)
        integer :: m, n, i, j, scale_exponent
        real(rk) :: detA, matrix_scale, matrix_norm, inverse_norm, row_norm, rcond
        real(rk) :: s11, s12, s13, s21, s22, s23, s31, s32, s33
        real(rk), allocatable :: A_transpose(:,:), normal_inv(:,:)

        m = size(A,1)
        n = size(A,2)
        if (m < 1 .or. n < 1) then
            allocate(A_inv(0,0))
            return
        end if

        if (.not. all(ieee_is_finite(A))) then
            allocate(A_inv(0,0))
            return
        end if
        matrix_scale = maxval(abs(A))
        if (.not. ieee_is_finite(matrix_scale) .or. matrix_scale <= 0.0_rk) then
            allocate(A_inv(0,0))
            return
        end if
        scale_exponent = exponent(matrix_scale)

        if (m == n) then
            select case(m)
              case(1)
                allocate(A_inv(1,1))
                A_inv(1,1) = 1.0_rk/scale(A(1,1),-scale_exponent)
              case(2)
                allocate(A_inv(m,n))
                s11 = scale(A(1,1),-scale_exponent)
                s12 = scale(A(1,2),-scale_exponent)
                s21 = scale(A(2,1),-scale_exponent)
                s22 = scale(A(2,2),-scale_exponent)
                detA = s11*s22 - s12*s21
                if (.not. ieee_is_finite(detA) .or. abs(detA) <= 0.0_rk) then
                    deallocate(A_inv)
                    allocate(A_inv(0,0))
                    return
                end if
                A_inv(1,1) =  s22
                A_inv(1,2) = -s12
                A_inv(2,1) = -s21
                A_inv(2,2) =  s11
                A_inv = A_inv/detA
              case(3)
                allocate(A_inv(m,n))
                s11 = scale(A(1,1),-scale_exponent)
                s12 = scale(A(1,2),-scale_exponent)
                s13 = scale(A(1,3),-scale_exponent)
                s21 = scale(A(2,1),-scale_exponent)
                s22 = scale(A(2,2),-scale_exponent)
                s23 = scale(A(2,3),-scale_exponent)
                s31 = scale(A(3,1),-scale_exponent)
                s32 = scale(A(3,2),-scale_exponent)
                s33 = scale(A(3,3),-scale_exponent)
                detA = s11*(s22*s33 - s23*s32) &
                    - s12*(s21*s33 - s23*s31) &
                    + s13*(s21*s32 - s22*s31)
                if (.not. ieee_is_finite(detA) .or. abs(detA) <= 0.0_rk) then
                    deallocate(A_inv)
                    allocate(A_inv(0,0))
                    return
                end if
                A_inv(1,1) = s22*s33 - s23*s32
                A_inv(1,2) = s13*s32 - s12*s33
                A_inv(1,3) = s12*s23 - s13*s22
                A_inv(2,1) = s23*s31 - s21*s33
                A_inv(2,2) = s11*s33 - s13*s31
                A_inv(2,3) = s13*s21 - s11*s23
                A_inv(3,1) = s21*s32 - s22*s31
                A_inv(3,2) = s12*s31 - s11*s32
                A_inv(3,3) = s11*s22 - s12*s21
                A_inv = A_inv/detA
              case default
                A_inv = solve(A,eye(m))
                if (size(A_inv,1) == 0) return
                A_inv = scale(A_inv,scale_exponent)
            end select

            if (.not. all(ieee_is_finite(A_inv))) then
                deallocate(A_inv)
                allocate(A_inv(0,0))
                return
            end if

            matrix_norm = 0.0_rk
            inverse_norm = 0.0_rk
            do i = 1, m
                row_norm = 0.0_rk
                do j = 1, n
                    row_norm = row_norm + abs(scale(A(i,j),-scale_exponent))
                end do
                matrix_norm = max(matrix_norm, row_norm)

                row_norm = 0.0_rk
                do j = 1, n
                    row_norm = row_norm + abs(A_inv(i,j))
                end do
                inverse_norm = max(inverse_norm, row_norm)
            end do

            if (.not. ieee_is_finite(matrix_norm) .or. &
                .not. ieee_is_finite(inverse_norm) .or. &
                matrix_norm <= 0.0_rk .or. inverse_norm <= 0.0_rk) then
                deallocate(A_inv)
                allocate(A_inv(0,0))
                return
            end if

            rcond = (1.0_rk/inverse_norm)/matrix_norm
            if (.not. ieee_is_finite(rcond) .or. rcond <= 32.0_rk*epsilon(1.0_rk)) then
                deallocate(A_inv)
                allocate(A_inv(0,0))
                return
            end if

            A_inv = scale(A_inv,-scale_exponent)
            if (.not. all(ieee_is_finite(A_inv))) then
                deallocate(A_inv)
                allocate(A_inv(0,0))
                return
            end if
        else if (m>n) then
            allocate(A_transpose(n,m))
            do i = 1, m
                do j = 1, n
                    A_transpose(j,i) = scale(A(i,j),-scale_exponent)
                end do
            end do
            normal_inv = inv(matmul(A_transpose, transpose(A_transpose)))
            if (size(normal_inv,1) == 0) then
                allocate(A_inv(0,0))
                return
            end if
            A_inv = matmul(normal_inv, A_transpose)
            A_inv = scale(A_inv,-scale_exponent)
        else if (m<n) then
            allocate(A_transpose(n,m))
            do i = 1, m
                do j = 1, n
                    A_transpose(j,i) = scale(A(i,j),-scale_exponent)
                end do
            end do
            normal_inv = inv(matmul(transpose(A_transpose), A_transpose))
            if (size(normal_inv,1) == 0) then
                allocate(A_inv(0,0))
                return
            end if
            A_inv = matmul(A_transpose, normal_inv)
            A_inv = scale(A_inv,-scale_exponent)
        end if

        if (.not. all(ieee_is_finite(A_inv))) then
            deallocate(A_inv)
            allocate(A_inv(0,0))
        end if

    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Allocate an `n` by `n` identity matrix.
    pure function eye(n) result(I)
        integer, intent(in) :: n
            !! Nonnegative matrix order.
        real(rk), allocatable :: I(:,:)

        ! local variables
        integer :: k

        allocate(I(n,n), source=0.0_rk)
        do concurrent (k = 1: n)
            I(k, k) = 1.0_rk
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Form the dyadic product of two real vectors.
    !!
    !! The allocated result has shape `[size(a),size(b)]` and entries
    !! \(c_{ij}=a_i b_j\).
    pure function dyad_t1_t1(a, b) result(c)
        real(rk), intent(in), contiguous :: a(:)
            !! Row-index factor.
        real(rk), intent(in), contiguous :: b(:)
            !! Column-index factor.
        real(rk), allocatable :: c(:,:)
            !! Allocated dyadic-product matrix.
        integer :: i, j

        allocate(c(size(a), size(b)))
        do j = 1, size(c,2)
            do i = 1, size(c,1)
                c(i,j) = a(i)*b(j)
            end do
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Generate a one-dimensional Gauss-Legendre rule on an interval.
    !!
    !! `degree+1` nodes integrate polynomials through degree `2*degree+1` in
    !! exact arithmetic. `degree>=0` and `size(interval)>=2` are unchecked shape
    !! preconditions. The first two interval entries must also be finite. When
    !! they are finite but not strictly increasing, allocated nodes and weights
    !! are filled with NaNs.
    pure subroutine gauss_legendre_1D(interval, degree, Xksi, Wksi)
        real(rk), intent(in), contiguous :: interval(:)
            !! Integration interval; the first two entries are used.
        integer, intent(in) :: degree
            !! Quadrature order minus one; must be nonnegative.
        real(rk), allocatable, intent(out) :: Xksi(:), Wksi(:)
            !! Allocated nodes and weights, each of length `degree+1`.

        allocate(Xksi(degree+1), Wksi(degree+1))

        call gauss_legendre(Xksi, Wksi, interval)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Generate a tensor-product Gauss-Legendre rule on a rectangle.
    !!
    !! Direction 1 varies fastest. `Xksi` has shape
    !! `[(degree(1)+1)*(degree(2)+1),2]`; `Wksi` contains the corresponding
    !! products of one-dimensional weights. `degree` must contain at least two
    !! nonnegative entries, and both interval arrays must contain at least two
    !! finite entries with increasing first pair; these preconditions are not
    !! shape-validated by this wrapper.
    pure subroutine gauss_legendre_2D(interval1, interval2, degree, Xksi, Wksi)
        real(rk), intent(in), contiguous :: interval1(:)
            !! Strictly increasing direction-1 interval.
        real(rk), intent(in), contiguous :: interval2(:)
            !! Strictly increasing direction-2 interval.
        integer, intent(in), contiguous :: degree(:)
            !! Quadrature orders minus one `[q1-1,q2-1]`.
        real(rk), allocatable, intent(out) :: Xksi(:,:), Wksi(:)
            !! Flattened parameter pairs and tensor-product weights.
        call gauss_legendre_tensor(reshape([interval1(:2),interval2(:2)],[2,2]),degree(:2),Xksi,Wksi)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Generate a tensor-product Gauss-Legendre rule on a cuboid.
    !!
    !! Direction 1 varies fastest, then direction 2. `Xksi` has three columns
    !! and one row per tensor-product point; `Wksi` uses the same ordering.
    !! `degree` must contain at least three nonnegative entries, and every
    !! interval must contain at least two finite entries with increasing first
    !! pair; these preconditions are not shape-validated by this wrapper.
    pure subroutine gauss_legendre_3D(interval1, interval2, interval3, degree, Xksi, Wksi)
        real(rk), intent(in), contiguous :: interval1(:)
            !! Strictly increasing direction-1 interval.
        real(rk), intent(in), contiguous :: interval2(:)
            !! Strictly increasing direction-2 interval.
        real(rk), intent(in), contiguous :: interval3(:)
            !! Strictly increasing direction-3 interval.
        integer, intent(in), contiguous :: degree(:)
            !! Quadrature orders minus one `[q1-1,q2-1,q3-1]`.
        real(rk), allocatable, intent(out) :: Xksi(:,:), Wksi(:)
            !! Flattened parameter triples and tensor-product weights.
        call gauss_legendre_tensor(reshape([interval1(:2),interval2(:2),interval3(:2)],[2,3]),degree(:3),Xksi,Wksi)
    end subroutine

    !> Tensor Gauss rule in direction-1-fastest order, or just its one-based point ig.
    !! Column d holds the interval for axis d. An invalid selection returns an empty rule.
    pure subroutine gauss_legendre_tensor(intervals,degree,Xksi,Wksi,ig)
        real(rk), contiguous, intent(in) :: intervals(:,:)
        integer, intent(in) :: degree(:)
        integer, optional, intent(in) :: ig
        real(rk), allocatable, intent(out) :: Xksi(:,:), Wksi(:)
        real(rk), allocatable :: nodes(:), weights(:)
        integer :: nq, stride, d, i, block, first, last

        nq=product(degree+1)
        if (any(degree<0) .or. size(intervals,1)<2 .or. size(intervals,2)<size(degree)) nq=0
        if (present(ig)) nq=merge(1,0,ig>=1 .and. ig<=nq)
        allocate(Xksi(nq,size(degree)),Wksi(nq))
        if (nq==0) return
        Wksi=1.0_rk
        stride=1
        do d=1,size(degree)
            call gauss_legendre_1D(intervals(:,d),degree(d),nodes,weights)
            if (present(ig)) then
                i=1+mod((ig-1)/stride,size(nodes))
                Xksi(1,d)=nodes(i)
                Wksi(1)=Wksi(1)*weights(i)
            else
                do block=0,nq-1,stride*size(nodes)
                    do i=1,size(nodes)
                        first=1+block+(i-1)*stride
                        last=first+stride-1
                        Xksi(first:last,d)=nodes(i)
                        Wksi(first:last)=Wksi(first:last)*weights(i)
                    end do
                end do
            end if
            stride=stride*size(nodes)
        end do
    end subroutine gauss_legendre_tensor
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute an \(n\)-point Gauss-Legendre rule and map it to an interval.
    !!
    !! Roots of \(P_n\) are initialized by a cosine approximation and refined
    !! with Newton iteration; symmetry halves the number of roots evaluated.
    !! `x` and `w` must have the same positive length, and `interval` must have
    !! at least two finite entries; these are caller preconditions. A finite
    !! non-increasing interval fills both arrays with quiet NaNs. On compilers other than
    !! NVFORTRAN, failure of a root iteration after 100 updates also returns
    !! NaNs. The NVFORTRAN compatibility path uses the last iterate instead of
    !! performing that post-loop failure assignment.
    pure subroutine gauss_legendre(x, w, interval)
        real(rk), intent(out) :: x(:), w(:)
            !! Nodes and weights; their common length is the quadrature order.
        real(rk), intent(in), contiguous :: interval(:)
            !! Strictly increasing interval in entries 1 and 2.
        real(rk) :: xi, delta, p_next, dp_next, p_prev, p_curr, dp_prev, dp_curr, midpoint, half_length
        integer :: i, j, k, n
        real(rk), parameter :: pi = acos(-1.0_rk)
        real(rk), parameter :: tol = 4.0_rk*epsilon(1.0_rk)
        integer, parameter :: maxit = 100
        logical :: converged

        if (size(interval) < 2 .or. interval(1) >= interval(2)) then
            x = ieee_value(1.0_rk, ieee_quiet_nan)
            w = ieee_value(1.0_rk, ieee_quiet_nan)
            return
        end if
        n = size(x)
        ! Gauss-Legendre points are symmetric, only compute half
        do i = 1, (n+1)/2
            ! Initial guess (Chebyshev approximation)
            xi = -cos(pi * (i-0.25_rk)/(n+0.5_rk))
            ! Newton iteration
            j = 0
            converged = .false.
            do while (.not. converged .and. j < maxit)
                j = j + 1
                ! Compute Legendre polynomial and derivative via recurrence
                p_prev = 1.0_rk        ! P_0(xi)
                p_curr = xi            ! P_1(xi)
                dp_prev = 0.0_rk       ! P_0d(xi)
                dp_curr = 1.0_rk       ! P_1d(xi)
                do k = 2, n
                    p_next  = ((2*k-1)*xi*p_curr-(k-1)*p_prev)/k
                    dp_next = ((2*k-1)*(xi*dp_curr+p_curr)-(k-1)*dp_prev)/k
                    p_prev  = p_curr
                    p_curr  = p_next
                    dp_prev = dp_curr
                    dp_curr = dp_next
                end do
                ! Newton correction
                delta = -p_curr / dp_curr
                xi = xi + delta
                ! Check for convergence
                converged = (abs(delta) <= tol*abs(xi))
            end do
#if !defined(__NVCOMPILER)
            if (.not. converged) then
                x = ieee_value(1.0_rk, ieee_quiet_nan)
                w = ieee_value(1.0_rk, ieee_quiet_nan)
                return
            end if
#endif
            ! Store symmetric nodes and weights
            x(i)     = xi
            x(n+1-i) = -xi
            w(i)     = 2.0_rk/((1.0_rk-xi**2)*dp_curr**2)
            w(n+1-i) = w(i)
        end do
        ! Transform from [-1,1] to [interval(1), interval(2)]
        midpoint    =0.5_rk*(interval(1)+interval(2))
        half_length =0.5_rk*(interval(2)-interval(1))
        x = midpoint+half_length*x
        w = half_length*w
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Write points, connectivity, and optional scalar fields in legacy VTK format.
    !!
    !! Connectivity is one-based on input and converted to VTK's zero-based
    !! indices; every index must lie in `1:size(points,1)`, but this is not
    !! checked. Two-dimensional points are emitted with a zero z coordinate.
    !! `vtkCellType=4` selects `POLYDATA/LINES` and permits any common row width.
    !! Other cell types require 2, 4, or 8 nodes per row and use ForCAD's
    !! line/quad/hex tensor-node reordering. The caller must pair the row width
    !! with the semantically correct VTK cell identifier.
    !!
    !! Point fields are emitted only when both `point_data` and `field_names`
    !! are present. Cell fields follow the same rule for `cell_data` and
    !! `cell_field_names`. Corresponding names are written verbatim after
    !! trimming and must already be valid VTK scalar identifiers. Supplying
    !! only one member of either data-name pair omits that association.
    !!
    !! Binary output converts coordinates and fields to IEEE `real64` and
    !! explicitly serializes big-endian 64-bit reals and 32-bit integers as
    !! required by the legacy format. Headers and binary blocks share one
    !! unformatted stream; explicit ASCII line feeds make the byte layout
    !! independent of the operating system's formatted-record convention. It
    !! requires 8-bit file storage units and the corresponding intrinsic
    !! integer and IEEE real kinds. Recognized encodings are exactly the lower-,
    !! title-, and upper-case forms of `ascii` and `binary`; binary is the
    !! default.
    !!
    !! The procedure has no status result. Rejected metadata returns before
    !! writing, and file-system I/O failures are not converted to a ForCAD
    !! diagnostic.
    !!
    !! Format reference: [VTK Legacy File Format](https://docs.vtk.org/en/latest/vtk_file_formats/vtk_legacy_file_format.html).
    impure subroutine export_vtk_legacy(&
        filename,&
        points,&
        elemConn,&
        vtkCellType,&
        point_data,&
        field_names,&
        encoding,&
        cell_data,&
        cell_field_names)
        character(len=*), intent(in) :: filename
            !! Destination path.
        real(rk), intent(in), contiguous :: points(:,:)
            !! Coordinates, shape `[npoint,2]` or `[npoint,3]`.
        integer, intent(in), contiguous :: elemConn(:,:) ! for VTK_POLY_LINE all rows have same nn
            !! One-based cell connectivity, one cell per row.
        integer, intent(in) :: vtkCellType
            !! Numeric VTK legacy cell identifier.
        real(rk), intent(in), contiguous, optional :: point_data(:,:)     ! [npoints, nfields]
            !! Optional scalar fields, shape `[npoint,nfield]`.
        character(len=*), intent(in), contiguous, optional :: field_names(:)
            !! Optional names corresponding to `point_data` columns.
        character(len=*), intent(in), optional :: encoding
            !! Optional `"ASCII"` or `"BINARY"`; binary is the default.
        real(rk), intent(in), contiguous, optional :: cell_data(:,:)
            !! Optional scalar cell fields, shape `[ncell,nfield]`.
        character(len=*), intent(in), contiguous, optional :: cell_field_names(:)
            !! Optional names corresponding to `cell_data` columns.

        integer :: byte_count, i, j, ne, np, nn, n, nunit
        integer(int8) :: byte_buffer(8192)
        character(len=128) :: header
        character(len=6) :: encoding_
        character, parameter :: lf = achar(10)
        integer, parameter :: vtk_node_order(8) = [1, 2, 4, 3, 5, 6, 8, 7]
        logical :: is_polyline

        ne = size(elemConn, 1)
        nn = size(elemConn, 2)
        np = size(points, 1)
        n = ne*(nn+1)

        is_polyline = (vtkCellType == 4)   ! VTK_POLY_LINE

        if (present(encoding)) then
            select case (trim(encoding))
                case ("ascii", "ASCII", "Ascii")
                    encoding_ = "ASCII"
                case ("binary", "BINARY", "Binary")
                    encoding_ = "BINARY"
                case default
                    return
            end select
        else
            encoding_ = "BINARY"
        end if

        if (size(points,2) /= 2 .and. size(points,2) /= 3) return
        if (.not. is_polyline .and. nn /= 2 .and. nn /= 4 .and. nn /= 8) return
        if (present(point_data) .and. present(field_names)) then
            if (size(point_data,1) /= np .or. size(field_names) < size(point_data,2)) return
        end if
        if (present(cell_data) .and. present(cell_field_names)) then
            if (size(cell_data,1) /= ne .or. size(cell_field_names) < size(cell_data,2)) return
        end if

        if (trim(encoding_) == "ASCII") then
            open(newunit=nunit, file=filename, action="write")
            write(nunit,"(a)") "# vtk DataFile Version 2.0"
            write(nunit,"(a)") "Generated by ForCAD"
            write(nunit,"(a)") "ASCII"
            if (is_polyline) then
                write(nunit,"(a)") "DATASET POLYDATA"
            else
                write(nunit,"(a)") "DATASET UNSTRUCTURED_GRID"
            end if

            write(nunit,"(a,1x,g0,1x,a)") "POINTS", np, "double"
            if (size(points,2) == 2) then
                write(nunit,"(g0,1x,g0,1x,g0)") (points(i,1), points(i,2), 0.0_rk , i = 1, np)
            else if (size(points,2) == 3) then
                write(nunit,"(g0,1x,g0,1x,g0)") (points(i,1), points(i,2), points(i,3) , i = 1, np)
            else
                return
            end if

            if (is_polyline) then
                write(nunit,"(a,1x,g0,1x,g0)") "LINES", ne, n
                do i = 1, ne
                    write(nunit,"(g0, *(1x,g0))") nn, (elemConn(i,j)-1, j=1, nn)
                end do
            else
                write(nunit,"(a,1x,g0,1x,g0)") "CELLS", ne, n
                select case (nn)
                    case (2)
                        write(nunit,"(g0,1x,g0,1x,g0)")&
                            (2, elemConn(i,1)-1,elemConn(i,2)-1, i = 1, ne)
                    case (4)
                        write(nunit,"(g0,1x,g0,1x,g0,1x,g0)")&
                            (4, elemConn(i,1)-1,elemConn(i,2)-1,elemConn(i,4)-1,elemConn(i,3)-1, i = 1, ne)
                    case (8)
                        write(nunit,"(g0,1x,g0,1x,g0,1x,g0,1x,g0,1x,g0,1x,g0,1x,g0,1x,g0)")&
                            (8, elemConn(i,1)-1,elemConn(i,2)-1,elemConn(i,4)-1,elemConn(i,3)-1,&
                            elemConn(i,5)-1,elemConn(i,6)-1,elemConn(i,8)-1,elemConn(i,7)-1, i = 1, ne)
                    case default
                        return
                end select

                write(nunit,"(a,1x,g0)") "CELL_TYPES", ne
                write(nunit,"(g0)") (vtkCellType , i = 1, ne)
            end if

            if (present(point_data) .and. present(field_names)) then
                write(nunit, "(a,1x,g0)") "POINT_DATA", size(point_data,1)
                do i = 1, size(point_data,2)
                    write(nunit, "(a,1x,a,1x,a)") "SCALARS", trim(field_names(i)), "double"
                    write(nunit, "(a)") "LOOKUP_TABLE default"
                    write(nunit, "(g0)") (point_data(j,i), j = 1, size(point_data,1))
                end do
            end if
            if (present(cell_data) .and. present(cell_field_names)) then
                write(nunit, "(a,1x,g0)") "CELL_DATA", size(cell_data,1)
                do i = 1, size(cell_data,2)
                    write(nunit, "(a,1x,a,1x,a)") "SCALARS", trim(cell_field_names(i)), "double"
                    write(nunit, "(a)") "LOOKUP_TABLE default"
                    write(nunit, "(g0)") (cell_data(j,i), j = 1, size(cell_data,1))
                end do
            end if

            close(nunit)
        end if


        if (trim(encoding_) == "BINARY") then
            if (file_storage_size /= 8 .or. bit_size(0_int8) /= 8 .or. bit_size(0_int32) /= 32 .or. &
                bit_size(0_int64) /= 64 .or. storage_size(0.0_real64) /= 64 .or. &
                .not. ieee_support_datatype(0.0_real64)) return

            open(&
                newunit = nunit,&
                file    = filename,&
                access  = "stream",&
                form    = "unformatted",&
                action  = "write",&
                status  = "replace")
            write(nunit) "# vtk DataFile Version 2.0", lf
            write(nunit) "Generated by ForCAD", lf
            write(nunit) "BINARY", lf
            if (is_polyline) then
                write(nunit) "DATASET POLYDATA", lf
            else
                write(nunit) "DATASET UNSTRUCTURED_GRID", lf
            end if

            write(header,"(a,1x,g0,1x,a)") "POINTS", np, "double"
            write(nunit) header(:len_trim(header)), lf
            byte_count = 0
            do i = 1, np
                do j = 1, size(points,2)
                    if (byte_count == size(byte_buffer)) then
                        write(nunit) byte_buffer
                        byte_count = 0
                    end if
                    call pack_big_endian_real64(&
                        value  = real(points(i,j),real64),&
                        bytes  = byte_buffer,&
                        offset = byte_count)
                    byte_count = byte_count + 8
                end do
                if (size(points,2) == 2) then
                    if (byte_count == size(byte_buffer)) then
                        write(nunit) byte_buffer
                        byte_count = 0
                    end if
                    call pack_big_endian_real64(&
                        value  = 0.0_real64,&
                        bytes  = byte_buffer,&
                        offset = byte_count)
                    byte_count = byte_count + 8
                end if
            end do
            if (byte_count > 0) write(nunit) byte_buffer(:byte_count)
            write(nunit) lf

            if (is_polyline) then
                write(header,"(a,1x,g0,1x,g0)") "LINES", ne, n
                write(nunit) header(:len_trim(header)), lf
                byte_count = 0
                do i = 1, ne
                    if (byte_count == size(byte_buffer)) then
                        write(nunit) byte_buffer
                        byte_count = 0
                    end if
                    call pack_big_endian_int32(&
                        value  = int(nn,int32),&
                        bytes  = byte_buffer,&
                        offset = byte_count)
                    byte_count = byte_count + 4
                    do j = 1, nn
                        if (byte_count == size(byte_buffer)) then
                            write(nunit) byte_buffer
                            byte_count = 0
                        end if
                        call pack_big_endian_int32(&
                            value  = int(elemConn(i,j)-1,int32),&
                            bytes  = byte_buffer,&
                            offset = byte_count)
                        byte_count = byte_count + 4
                    end do
                end do
                if (byte_count > 0) write(nunit) byte_buffer(:byte_count)
                write(nunit) lf
            else
                write(header,"(a,1x,g0,1x,g0)") "CELLS", ne, n
                write(nunit) header(:len_trim(header)), lf
                byte_count = 0
                do i = 1, ne
                    if (byte_count == size(byte_buffer)) then
                        write(nunit) byte_buffer
                        byte_count = 0
                    end if
                    call pack_big_endian_int32(&
                        value  = int(nn,int32),&
                        bytes  = byte_buffer,&
                        offset = byte_count)
                    byte_count = byte_count + 4
                    do j = 1, nn
                        if (byte_count == size(byte_buffer)) then
                            write(nunit) byte_buffer
                            byte_count = 0
                        end if
                        call pack_big_endian_int32(&
                            value  = int(elemConn(i,vtk_node_order(j))-1,int32),&
                            bytes  = byte_buffer,&
                            offset = byte_count)
                        byte_count = byte_count + 4
                    end do
                end do
                if (byte_count > 0) write(nunit) byte_buffer(:byte_count)
                write(nunit) lf

                write(header,"(a,1x,g0)") "CELL_TYPES", ne
                write(nunit) header(:len_trim(header)), lf
                byte_count = 0
                do i = 1, ne
                    if (byte_count == size(byte_buffer)) then
                        write(nunit) byte_buffer
                        byte_count = 0
                    end if
                    call pack_big_endian_int32(&
                        value  = int(vtkCellType,int32),&
                        bytes  = byte_buffer,&
                        offset = byte_count)
                    byte_count = byte_count + 4
                end do
                if (byte_count > 0) write(nunit) byte_buffer(:byte_count)
                write(nunit) lf
            end if

            if (present(point_data) .and. present(field_names)) then
                write(header,"(a,1x,g0)") "POINT_DATA", size(point_data,1)
                write(nunit) header(:len_trim(header)), lf

                do i = 1, size(point_data,2)
                    write(nunit) "SCALARS ", trim(field_names(i)), " double", lf
                    write(nunit) "LOOKUP_TABLE default", lf
                    byte_count = 0
                    do j = 1, size(point_data,1)
                        if (byte_count == size(byte_buffer)) then
                            write(nunit) byte_buffer
                            byte_count = 0
                        end if
                        call pack_big_endian_real64(&
                            value  = real(point_data(j,i),real64),&
                            bytes  = byte_buffer,&
                            offset = byte_count)
                        byte_count = byte_count + 8
                    end do
                    if (byte_count > 0) write(nunit) byte_buffer(:byte_count)
                    write(nunit) lf
                end do
            end if
            if (present(cell_data) .and. present(cell_field_names)) then
                write(header,"(a,1x,g0)") "CELL_DATA", size(cell_data,1)
                write(nunit) header(:len_trim(header)), lf

                do i = 1, size(cell_data,2)
                    write(nunit) "SCALARS ", trim(cell_field_names(i)), " double", lf
                    write(nunit) "LOOKUP_TABLE default", lf
                    byte_count = 0
                    do j = 1, size(cell_data,1)
                        if (byte_count == size(byte_buffer)) then
                            write(nunit) byte_buffer
                            byte_count = 0
                        end if
                        call pack_big_endian_real64(&
                            value  = real(cell_data(j,i),real64),&
                            bytes  = byte_buffer,&
                            offset = byte_count)
                        byte_count = byte_count + 8
                    end do
                    if (byte_count > 0) write(nunit) byte_buffer(:byte_count)
                    write(nunit) lf
                end do
            end if

            close(nunit)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Append one 32-bit integer to a byte buffer in VTK big-endian order.
    pure subroutine pack_big_endian_int32(value, bytes, offset)
        integer(int32), intent(in) :: value
        integer(int8), intent(inout), contiguous :: bytes(:)
        integer, intent(in) :: offset

        integer :: byte

        do byte = 1, 4
            bytes(offset+byte) = int(ibits(value,8*(4-byte),7),int8)
            if (btest(value,8*(4-byte)+7)) bytes(offset+byte) = ibset(bytes(offset+byte),7)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Append one IEEE 64-bit real to a byte buffer in VTK big-endian order.
    pure subroutine pack_big_endian_real64(value, bytes, offset)
        real(real64), intent(in) :: value
        integer(int8), intent(inout), contiguous :: bytes(:)
        integer, intent(in) :: offset

        integer :: byte
        integer(int64) :: bits

        bits = transfer(value, 0_int64)
        do byte = 1, 8
            bytes(offset+byte) = int(ibits(bits,8*(8-byte),7),int8)
            if (btest(bits,8*(8-byte)+7)) bytes(offset+byte) = ibset(bytes(offset+byte),7)
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Solve the square dense system \(AX=B\) for one or more right-hand sides.
    !!
    !! A matrix is classified as numerically symmetric when every compared
    !! pair satisfies the scale-invariant test
    !!
    !! \[
    !! |\widehat A_{ij}-\widehat A_{ji}|\le
    !! 32\,\epsilon_{rk},\qquad
    !! \widehat A=A/\operatorname{radix}(A)^e,
    !! \]
    !!
    !! That path first attempts unpivoted Cholesky using the lower triangle.
    !! Before accepting its solution, the full input matrix is checked with the
    !! componentwise backward-residual bound
    !!
    !! \[
    !! |(AX-B)_{ir}|\le 64n\epsilon_{rk}
    !! \left(|B_{ir}|+\sum_j|A_{ij}|\,|X_{jr}|\right).
    !! \]
    !!
    !! where \(e\) is the model exponent of the largest entry. Scaling by a
    !! radix power is exact. Cholesky and LU both solve
    !! \(\widehat A X=\widehat B\), with \(\widehat B=B/\operatorname{radix}(A)^e\),
    !! so their pivot decisions and arithmetic do not depend on the physical
    !! units of the system.
    !!
    !! If Cholesky encounters a nonpositive pivot, the residual check fails, or
    !! the symmetry test fails, LU factorization with partial pivoting solves
    !! the normalized general system. Shape mismatch, an exactly zero LU
    !! pivot, a nonfinite pivot, or a nonfinite solution returns `X(0,0)`. A
    !! small nonzero normalized pivot is accepted; its magnitude alone does not
    !! prove singularity.
    !!
    !! The routine is allocation-based and intended for small or medium dense
    !! systems. Work is \(O(n^3+n^2n_{rhs})\), and dense factor and solution
    !! storage is allocated internally. The routine does not estimate a
    !! condition number or perform iterative refinement.
    !! All entries of `A` and `B` must be finite.
    !!
    !! **Reference implementations:**
    !!
    !! - LAPACK `DPOSV`, symmetric positive-definite solution by Cholesky
    !!   factorization:
    !!   [Netlib LAPACK](https://netlib.org/lapack/explore-html/de/d6c/group__posv_ga4844053bd30fe88a17a7e08d93bbae4b.html).
    !! - LAPACK `DGESV`, general solution by LU factorization with partial
    !!   pivoting:
    !!   [Netlib LAPACK](https://netlib.org/lapack/explore-html/d8/da6/group__gesv_ga831ce6a40e7fd16295752d18aed2d541.html).
    pure function solve(A, B) result(X)
        real(rk), intent(in), contiguous :: A(:,:)
            !! Square coefficient matrix, shape `[n,n]`.
        real(rk), intent(in), contiguous :: B(:,:)
            !! Right-hand sides, shape `[n,nrhs]`.
        real(rk), allocatable :: X(:,:)

        integer :: n, m, i, j, p, r, scale_exponent
        real(rk) :: matrix_scale, sym_tol, residual, residual_scale, residual_tol
        logical :: symmetric, ok

        p = size(A,1)
        n = size(A,2)
        m = size(B,2)

        if (p /= n .or. p /= size(B,1) .or. n < 1) then
            allocate(X(0,0))
            return
        end if

        if (.not. all(ieee_is_finite(A)) .or. &
            .not. all(ieee_is_finite(B))) then
            allocate(X(0,0))
            return
        end if
        matrix_scale = maxval(abs(A))
        if (matrix_scale <= 0.0_rk) then
            allocate(X(0,0))
            return
        end if
        scale_exponent = exponent(matrix_scale)
        sym_tol = 32.0_rk*epsilon(1.0_rk)
        symmetric = .true.
        do i = 2, n
            do j = 1, i-1
                if (abs(scale(A(i,j),-scale_exponent) - &
                    scale(A(j,i),-scale_exponent)) > sym_tol) then
                    symmetric = .false.
                    exit
                end if
            end do
            if (.not. symmetric) exit
        end do

        if (symmetric) then
            call solve_cholesky(A, B, scale_exponent, X, ok)
            if (ok) then
                residual_tol = 64.0_rk*real(n,rk)*epsilon(1.0_rk)
                do r = 1, m
                    do i = 1, n
                        residual = -scale(B(i,r),-scale_exponent)
                        residual_scale = abs(residual)
                        do j = 1, n
                            residual = residual + &
                                scale(A(i,j),-scale_exponent)*X(j,r)
                            residual_scale = residual_scale + &
                                abs(scale(A(i,j),-scale_exponent))*abs(X(j,r))
                        end do
                        if (.not. ieee_is_finite(residual) .or. &
                            .not. ieee_is_finite(residual_scale) .or. &
                            abs(residual) > residual_tol*residual_scale) then
                            ok = .false.
                            exit
                        end if
                    end do
                    if (.not. ok) exit
                end do
                if (ok) return
            end if
        end if

        call solve_lu(A, B, scale_exponent, X)
        if (size(X,1) > 0 .and. .not. all(ieee_is_finite(X))) then
            deallocate(X)
            allocate(X(0,0))
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Solve a symmetric positive-definite band system by banded Cholesky factorization.
    !!
    !! `A_band(0,j)` is \(A_{jj}\), and `A_band(i-j,j)` stores the lower entry
    !! \(A_{ij}\) for `0<=i-j<=kd`; entries outside that symmetric band are
    !! assumed zero. Factorization costs \(O(n\,kd^2)\), and triangular solves
    !! cost \(O(n\,kd\,n_{rhs})\). Internal factor and solution storage is
    !! \(O(n(kd+1)+n\,n_{rhs})\), avoiding a dense \(n\times n\) matrix.
    !!
    !! Before factorization, the valid stored entries and right-hand sides are
    !! scaled by the exact radix power
    !! \(s=\operatorname{radix}(A)^e\), where \(e\) is the model exponent of
    !! the largest stored magnitude. The resulting factorization
    !! \((A/s)=LL^T\) has the same solution because it solves
    !! \((A/s)X=B/s\). A normalized Cholesky pivot not greater than
    !! \(32\epsilon_{rk}\max_{ij}|A_{ij}/s|\), or a nonfinite pivot, returns
    !! `X(0,0)`. Consequently, positive uniform rescaling cannot change the
    !! numerical-definiteness decision. The routine does not estimate a
    !! condition number. Right-hand-side entries must be finite.
    !!
    !! **Reference implementation:** LAPACK `DPBTRF` factors symmetric
    !! positive-definite band matrices and reports a failed positive-definite
    !! leading principal minor:
    !! [Netlib LAPACK](https://netlib.org/lapack/explore-html/dc/d58/group__pbtrf_gafd2c3ba375ecee6dba1fba62696993e0.html).
    !!
    pure function solve_spd_banded(A_band, B) result(X)
        real(rk), intent(in), contiguous :: A_band(0:,:)
            !! Lower symmetric band storage, lower bound zero in dimension 1.
        real(rk), intent(in), contiguous :: B(:,:)
            !! Right-hand sides, shape `[n,nrhs]`.
        real(rk), allocatable :: X(:,:)
        real(rk), allocatable :: L(:,:)
        integer :: kd, n, m, i, j, k, r, k0, i1, scale_exponent
        real(rk) :: acc, matrix_scale, pivot_tol

        kd = size(A_band, 1) - 1
        n = size(A_band, 2)
        m = size(B, 2)
        if (size(B,1) /= n .or. n < 1 .or. kd < 0) then
            allocate(X(0,0))
            return
        end if

        matrix_scale = 0.0_rk
        do j = 1, n
            do i = 0, min(kd,n-j)
                if (.not. ieee_is_finite(A_band(i,j))) then
                    allocate(X(0,0))
                    return
                end if
                matrix_scale = max(matrix_scale, abs(A_band(i,j)))
            end do
        end do
        if (matrix_scale <= 0.0_rk .or. .not. all(ieee_is_finite(B))) then
            allocate(X(0,0))
            return
        end if
        scale_exponent = exponent(matrix_scale)
        pivot_tol = 32.0_rk*epsilon(1.0_rk)*&
            scale(matrix_scale,-scale_exponent)

        allocate(L(0:kd,n), X(n,m), source=0.0_rk)
        X = scale(B,-scale_exponent)

        do j = 1, n
            acc = scale(A_band(0,j),-scale_exponent)
            k0 = max(1, j-kd)
            do k = k0, j-1
                acc = acc - L(j-k,k)*L(j-k,k)
            end do
            if (.not. ieee_is_finite(acc) .or. acc <= pivot_tol) then
                deallocate(X)
                allocate(X(0,0))
                return
            end if
            L(0,j) = sqrt(acc)

            do i = j+1, min(n, j+kd)
                acc = scale(A_band(i-j,j),-scale_exponent)
                k0 = max(1, max(j-kd, i-kd))
                do k = k0, j-1
                    acc = acc - L(i-k,k)*L(j-k,k)
                end do
                L(i-j,j) = acc/L(0,j)
            end do
        end do

        do r = 1, m
            do i = 1, n
                acc = X(i,r)
                i1 = max(1, i-kd)
                do j = i1, i-1
                    acc = acc - L(i-j,j)*X(j,r)
                end do
                X(i,r) = acc/L(0,i)
            end do

            do i = n, 1, -1
                acc = X(i,r)
                do j = i+1, min(n, i+kd)
                    acc = acc - L(j-i,i)*X(j,r)
                end do
                X(i,r) = acc/L(0,i)
            end do
        end do

        if (.not. all(ieee_is_finite(X))) then
            deallocate(X)
            allocate(X(0,0))
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Solve a dense symmetric positive-definite system by Cholesky factorization.
    !!
    !! This private path factors \(A=LL^T\) without pivoting, then performs
    !! forward and backward substitutions for every right-hand side. A
    !! normalized pivot not greater than
    !! \(32\epsilon_{rk}\max_{ij}|A_{ij}/s|\), or a nonfinite pivot, sets
    !! `ok=.false.`; `X` remains allocated but must not be used unless `ok` is
    !! true. Finite inputs are a caller precondition.
    !!
    !! **First published account:** Commandant Benoit, "Note sur une methode de
    !! resolution des equations normales provenant de l'application de la
    !! methode des moindres carres a un systeme d'equations lineaires en nombre
    !! inferieur a celui des inconnues -- Application de la methode a la
    !! resolution d'un systeme defini d'equations lineaires (procede du
    !! Commandant Cholesky)," *Bulletin Geodesique* 2 (1924), 67--77.
    !! [doi:10.1007/BF03031308](https://doi.org/10.1007/BF03031308).
    pure subroutine solve_cholesky(A, B, scale_exponent, X, ok)
        real(rk), intent(in), contiguous :: A(:,:)
            !! Symmetric coefficient matrix `[n,n]`.
        real(rk), intent(in), contiguous :: B(:,:)
            !! Right-hand sides `[n,nrhs]`.
        integer, intent(in) :: scale_exponent
            !! Radix-power exponent used to normalize `A` and `B`.
        real(rk), allocatable, intent(out) :: X(:,:)
            !! Solution matrix `[n,nrhs]` when `ok` is true.
        logical, intent(out) :: ok
            !! True only when all Cholesky pivots are positive.
        integer :: n, m, i, j, k
        real(rk), allocatable :: L(:,:)
        real(rk) :: acc, pivot_tol

        n = size(A,2)
        m = size(B,2)
        ok = .true.
        pivot_tol = 32.0_rk*epsilon(1.0_rk)*&
            scale(maxval(abs(A)),-scale_exponent)
        allocate(L(n,n), X(n,m), source=0.0_rk)
        X = scale(B,-scale_exponent)

        do i = 1, n
            do j = 1, i
                acc = scale(A(i,j),-scale_exponent)
                do k = 1, j-1
                    acc = acc - L(i,k)*L(j,k)
                end do
                if (i == j) then
                    if (.not. ieee_is_finite(acc) .or. acc <= pivot_tol) then
                        ok = .false.
                        return
                    end if
                    L(i,j) = sqrt(acc)
                else
                    L(i,j) = acc/L(j,j)
                end if
            end do
        end do

        ! Forward substitution: L*X = B
        do j = 1,m
            do i = 1, n
                acc = X(i,j)
                do k = 1, i-1
                    acc = acc - L(i,k)*X(k,j)
                end do
                X(i,j) = acc/L(i,i)
            end do
        end do

        ! Backward substitution: transpose(L)*X = X
        do j = 1,m
            do i = n, 1, -1
                acc = X(i,j)
                do k = i+1, n
                    acc = acc - L(k,i)*X(k,j)
                end do
                X(i,j) = acc/L(i,i)
            end do
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Solve a dense square system by Gaussian elimination with partial pivoting.
    !!
    !! Row interchanges use the largest available absolute pivot in each
    !! column. Invalid shapes, an exactly zero pivot, or a nonfinite pivot
    !! returns `X(0,0)`. As in LAPACK `DGESV`, no fixed magnitude threshold is
    !! used to infer singularity: uniform rescaling therefore does not change
    !! whether a finite nonzero pivot is accepted. A very small pivot can still
    !! indicate severe ill-conditioning; this routine does not estimate the
    !! condition number. Finite right-hand-side entries are a caller
    !! precondition.
    !!
    !! **Reference:** J. von Neumann and H. H. Goldstine, "Numerical Inverting
    !! of Matrices of High Order," *Bulletin of the American Mathematical
    !! Society* 53 (1947), 1021--1099.
    !! [doi:10.1090/S0002-9904-1947-08909-6](https://doi.org/10.1090/S0002-9904-1947-08909-6).
    !!
    !! **Reference implementation:** LAPACK `DGESV`,
    !! [Netlib LAPACK](https://netlib.org/lapack/explore-html/d8/da6/group__gesv_ga831ce6a40e7fd16295752d18aed2d541.html).
    pure subroutine solve_lu(A, B, scale_exponent, X)
        real(rk), intent(in), contiguous :: A(:,:)
            !! Square coefficient matrix `[n,n]`.
        real(rk), intent(in), contiguous :: B(:,:)
            !! Right-hand sides `[n,nrhs]`.
        integer, intent(in) :: scale_exponent
            !! Radix-power exponent used to normalize `A` and `B`.
        real(rk), allocatable, intent(out) :: X(:,:)
            !! Solution matrix, or an empty matrix on failure.
        integer :: n, m, i, j, k, pivot
        real(rk), allocatable :: U(:,:)
        real(rk) :: factor, pivot_abs, candidate_abs, tmp, acc

        n = size(A,2)
        m = size(B,2)
        if (n < 1 .or. size(A,1) /= n .or. size(B,1) /= n) then
            allocate(X(0,0))
            return
        end if
        allocate(U(n,n), X(n,m))
        U = scale(A,-scale_exponent)
        X = scale(B,-scale_exponent)

        do k = 1, n-1
            pivot = k
            pivot_abs = abs(U(k,k))
            do i = k+1, n
                candidate_abs = abs(U(i,k))
                if (candidate_abs > pivot_abs) then
                    pivot_abs = candidate_abs
                    pivot = i
                end if
            end do
            if (.not. ieee_is_finite(pivot_abs) .or. pivot_abs <= 0.0_rk) then
                deallocate(X)
                allocate(X(0,0))
                return
            end if

            if (pivot /= k) then
                do j = k, n
                    tmp = U(k,j)
                    U(k,j) = U(pivot,j)
                    U(pivot,j) = tmp
                end do
                do j = 1, m
                    tmp = X(k,j)
                    X(k,j) = X(pivot,j)
                    X(pivot,j) = tmp
                end do
            end if

            do i = k+1, n
                factor = U(i,k)/U(k,k)
                U(i,k) = 0.0_rk
                do j = k+1, n
                    U(i,j) = U(i,j) - factor*U(k,j)
                end do
                do j = 1, m
                    X(i,j) = X(i,j) - factor*X(k,j)
                end do
            end do
        end do

        if (.not. ieee_is_finite(U(n,n)) .or. abs(U(n,n)) <= 0.0_rk) then
            deallocate(X)
            allocate(X(0,0))
            return
        end if

        do j = 1, m
            do i = n, 1, -1
                acc = X(i,j)
                do k = i+1, n
                    acc = acc - U(i,k)*X(k,j)
                end do
                X(i,j) = acc/U(i,i)
            end do
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Run a generated Python script without making the shell parse the source.
    !! This private visualization bridge invokes `python` through
    !! `execute_command_line`; it is impure, synchronous, and intended only for
    !! the optional PyVista viewer.
    !!
    impure subroutine run_python_script(script)
        character(len=*), intent(in) :: script
            !! Complete Python source text.
        character(len=:), allocatable :: command

        command = "python - <<'FORCAD_PYVISTA_EOF'"//achar(10)//&
            trim(script)//achar(10)//&
            "FORCAD_PYVISTA_EOF"
        call execute_command_line(command)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Return a Python single-quoted string literal.
    pure function pyvista_string_literal(text) result(literal)
        character(len=*), intent(in) :: text
        character(len=:), allocatable :: literal
        integer :: i

        literal = "'"
        do i = 1, len_trim(text)
            select case (text(i:i))
            case ("'")
                literal = literal//"\'"
            case ("\")
                literal = literal//"\\"
            case default
                literal = literal//text(i:i)
            end select
        end do
        literal = literal//"'"
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Build the shared PyVista viewer body used by all single- and multi-patch ranks.
    function pyvista_viewer_script(header) result(script)
        character(len=*), intent(in) :: header
        character(len=:), allocatable :: script

        script = trim(header)//achar(10)//&
            "import numpy as np"//achar(10)//&
            "import pyvista as pv"//achar(10)//&
            "pv.global_theme.color = 'white'"//achar(10)//&
            "if len(Xc_files) == 0 or len(Xg_files) == 0:"//achar(10)//&
            "    raise RuntimeError('ForCAD PyVista viewer needs Xc and Xg VTK files')"//achar(10)//&
            "Xth_files = Xth_files if 'Xth_files' in globals() else []"//achar(10)//&
            "colormap_names = ['viridis', 'plasma', 'inferno', 'magma', 'cividis',"//&
            " 'turbo', 'coolwarm', 'RdBu', 'Spectral', 'Greys']"//achar(10)//&
            "Xc_meshes = [pv.read(path) for path in Xc_files]"//achar(10)//&
            "Xg_meshes = [pv.read(path) for path in Xg_files]"//achar(10)//&
            "Xth_meshes = [pv.read(path) for path in Xth_files]"//achar(10)//&
            "common_fields = [name for name in Xg_meshes[0].point_data.keys()"//&
            " if all(name in mesh.point_data for mesh in Xg_meshes)]"//achar(10)//&
            "field_names = []"//achar(10)//&
            "field_ranges = {}"//achar(10)//&
            "for field_name in common_fields:"//achar(10)//&
            "    field_min = None"//achar(10)//&
            "    field_max = None"//achar(10)//&
            "    scalar_field = True"//achar(10)//&
            "    for mesh in Xg_meshes:"//achar(10)//&
            "        values = np.asarray(mesh.point_data[field_name])"//achar(10)//&
            "        if values.ndim == 2 and values.shape[1] == 1:"//achar(10)//&
            "            values = values[:, 0]"//achar(10)//&
            "        if values.ndim != 1:"//achar(10)//&
            "            scalar_field = False"//achar(10)//&
            "            break"//achar(10)//&
            "        local_min = float(np.nanmin(values))"//achar(10)//&
            "        local_max = float(np.nanmax(values))"//achar(10)//&
            "        if not np.isfinite(local_min) or not np.isfinite(local_max):"//achar(10)//&
            "            scalar_field = False"//achar(10)//&
            "            break"//achar(10)//&
            "        field_min = local_min if field_min is None else min(field_min, local_min)"//achar(10)//&
            "        field_max = local_max if field_max is None else max(field_max, local_max)"//achar(10)//&
            "    if scalar_field:"//achar(10)//&
            "        field_names.append(field_name)"//achar(10)//&
            "        field_ranges[field_name] = (field_min, field_max)"//achar(10)//&
            "p = pv.Plotter(lighting='light kit')"//achar(10)//&
            "actor_Xcp = [p.add_mesh(mesh, style='points', point_size=10, color='red',"//&
            " render_points_as_spheres=True, opacity=0.5) for mesh in Xc_meshes]"//achar(10)//&
            "actor_Xcw = [p.add_mesh(mesh, show_edges=True, color='yellow', line_width=3,"//&
            " style='wireframe', opacity=0.2) for mesh in Xc_meshes]"//achar(10)//&
            "actor_Xg = [p.add_mesh(mesh, show_edges=False, color='cyan', line_width=1,"//&
            " metallic=0.6, pbr=True, smooth_shading=True, split_sharp_edges=True)"//&
            " for mesh in Xg_meshes]"//achar(10)//&
            "actor_Xth = [p.add_mesh(mesh, style='wireframe', color='magenta', line_width=2,"//&
            " opacity=0.8) for mesh in Xth_meshes]"//achar(10)//&
            "p.add_axes(interactive=False)"//achar(10)//&
            "def point_picker_callback(point):"//achar(10)//&
            "    best = None"//achar(10)//&
            "    for mesh in Xc_meshes:"//achar(10)//&
            "        point_id = mesh.find_closest_point(point)"//achar(10)//&
            "        point_coords = mesh.points[point_id]"//achar(10)//&
            "        dx = point_coords[0] - point[0]"//achar(10)//&
            "        dy = point_coords[1] - point[1]"//achar(10)//&
            "        dz = point_coords[2] - point[2]"//achar(10)//&
            "        dist2 = dx*dx + dy*dy + dz*dz"//achar(10)//&
            "        if best is None or dist2 < best[0]:"//achar(10)//&
            "            best = (dist2, point_id, point_coords)"//achar(10)//&
            "    _, point_id, point_coords = best"//achar(10)//&
            "    label = f'ID: {point_id + 1}\n({point_coords[0]:.3f}, "//&
            "{point_coords[1]:.3f}, {point_coords[2]:.3f})'"//achar(10)//&
            "    p.add_point_labels([point_coords], [label], font_size=14, text_color='black',"//&
            " show_points=False, fill_shape=False, shape=None)"//achar(10)//&
            "p.enable_point_picking(callback=point_picker_callback, show_message=False)"//achar(10)//&
            "window_size = p.window_size"//achar(10)//&
            "y_pos = window_size[1]"//achar(10)//&
            "Xg_visible = [True]"//achar(10)//&
            "active_field = [None]"//achar(10)//&
            "active_colormap = [colormap_names[0]]"//achar(10)//&
            "scalar_bar = None"//achar(10)//&
            "def Xcp_toggle_vis(flag):"//achar(10)//&
            "    for actor in actor_Xcp: actor.SetVisibility(flag)"//achar(10)//&
            "def Xcw_toggle_vis(flag):"//achar(10)//&
            "    for actor in actor_Xcw: actor.SetVisibility(flag)"//achar(10)//&
            "def Xg_toggle_vis(flag):"//achar(10)//&
            "    Xg_visible[0] = bool(flag)"//achar(10)//&
            "    for actor in actor_Xg: actor.SetVisibility(flag)"//achar(10)//&
            "    if scalar_bar is not None:"//achar(10)//&
            "        scalar_bar.SetVisibility(bool(flag) and active_field[0] is not None)"//achar(10)//&
            "def Xth_toggle_vis(flag):"//achar(10)//&
            "    for actor in actor_Xth: actor.SetVisibility(flag)"//achar(10)//&
            "p.add_checkbox_button_widget(Xcp_toggle_vis, value=True, color_on='red',"//&
            " size=25, position=(0, y_pos - 1*25))"//achar(10)//&
            "p.add_checkbox_button_widget(Xcw_toggle_vis, value=True, color_on='yellow',"//&
            " size=25, position=(0, y_pos - 2*25))"//achar(10)//&
            "p.add_checkbox_button_widget(Xg_toggle_vis, value=True, color_on='cyan',"//&
            " size=25, position=(0, y_pos - 3*25))"//achar(10)//&
            "p.add_text('Xc (Points)', position=(28, y_pos - 1*25), font_size=8,"//&
            " color='black', font='times')"//achar(10)//&
            "p.add_text('Xc (Control geometry)', position=(28, y_pos - 2*25), font_size=8,"//&
            " color='black', font='times')"//achar(10)//&
            "Xg_label = p.add_text('Xg (Geometry)', position=(28, y_pos - 3*25), font_size=8,"//&
            " color='black', font='times')"//achar(10)//&
            "if len(actor_Xth) > 0:"//achar(10)//&
            "    p.add_checkbox_button_widget(Xth_toggle_vis, value=True, color_on='magenta',"//&
            " size=25, position=(0, y_pos - 4*25))"//achar(10)//&
            "    p.add_text('Xth (Parameter)', position=(28, y_pos - 4*25), font_size=8,"//&
            " color='black', font='times')"//achar(10)//&
            "def select_Xg_field(field_name):"//achar(10)//&
            "    if field_name == 'Geometry':"//achar(10)//&
            "        active_field[0] = None"//achar(10)//&
            "        for actor in actor_Xg: actor.mapper.scalar_visibility = False"//achar(10)//&
            "        if scalar_bar is not None: scalar_bar.SetVisibility(False)"//achar(10)//&
            "        Xg_label.SetInput('Xg (Geometry)')"//achar(10)//&
            "    else:"//achar(10)//&
            "        active_field[0] = field_name"//achar(10)//&
            "        for actor in actor_Xg:"//achar(10)//&
            "            actor.mapper.set_active_scalars(field_name, preference='point')"//achar(10)//&
            "            actor.mapper.lookup_table.apply_cmap(active_colormap[0])"//achar(10)//&
            "            actor.mapper.scalar_range = field_ranges[field_name]"//achar(10)//&
            "            actor.mapper.scalar_visibility = True"//achar(10)//&
            "        if scalar_bar is not None:"//achar(10)//&
            "            scalar_bar.SetTitle(field_name)"//achar(10)//&
            "            scalar_bar.SetVisibility(Xg_visible[0])"//achar(10)//&
            "        Xg_label.SetInput(f'Xg ({field_name})')"//achar(10)//&
            "    p.render()"//achar(10)//&
            "def select_colormap(colormap_name):"//achar(10)//&
            "    active_colormap[0] = colormap_name"//achar(10)//&
            "    if active_field[0] is not None:"//achar(10)//&
            "        for actor in actor_Xg:"//achar(10)//&
            "            actor.mapper.lookup_table.apply_cmap(colormap_name)"//achar(10)//&
            "            actor.mapper.scalar_range = field_ranges[active_field[0]]"//achar(10)//&
            "    p.render()"//achar(10)//&
            "if field_names:"//achar(10)//&
            "    select_Xg_field(field_names[0])"//achar(10)//&
            "    scalar_bar = p.add_scalar_bar(title=field_names[0], mapper=actor_Xg[0].mapper)"//achar(10)//&
            "    field_menu_widgets = []"//achar(10)//&
            "    field_menu_labels = []"//achar(10)//&
            "    colormap_menu_widgets = []"//achar(10)//&
            "    colormap_menu_labels = []"//achar(10)//&
            "    field_menu_open = [False]"//achar(10)//&
            "    colormap_menu_open = [False]"//achar(10)//&
            "    field_menu_button = [None]"//achar(10)//&
            "    colormap_menu_button = [None]"//achar(10)//&
            "    def set_menu_visibility(widgets, labels, visible):"//achar(10)//&
            "        for widget in widgets:"//achar(10)//&
            "            if visible:"//achar(10)//&
            "                widget.On()"//achar(10)//&
            "                widget.GetRepresentation().VisibilityOn()"//achar(10)//&
            "            else:"//achar(10)//&
            "                widget.Off()"//achar(10)//&
            "                widget.GetRepresentation().VisibilityOff()"//achar(10)//&
            "        for label in labels: label.SetVisibility(visible)"//achar(10)//&
            "    def close_field_menu():"//achar(10)//&
            "        field_menu_open[0] = False"//achar(10)//&
            "        set_menu_visibility(field_menu_widgets, field_menu_labels, False)"//achar(10)//&
            "        if field_menu_button[0] is not None:"//achar(10)//&
            "            field_menu_button[0].GetRepresentation().SetState(0)"//achar(10)//&
            "        selected_name = 'Geometry' if active_field[0] is None else active_field[0]"//achar(10)//&
            "        field_control_label.SetInput(f'Field: {selected_name}  v')"//achar(10)//&
            "    def close_colormap_menu():"//achar(10)//&
            "        colormap_menu_open[0] = False"//achar(10)//&
            "        set_menu_visibility(colormap_menu_widgets, colormap_menu_labels, False)"//achar(10)//&
            "        if colormap_menu_button[0] is not None:"//achar(10)//&
            "            colormap_menu_button[0].GetRepresentation().SetState(0)"//achar(10)//&
            "        colormap_control_label.SetInput(f'Colormap: {active_colormap[0]}  v')"//achar(10)//&
            "    def choose_field(field_name):"//achar(10)//&
            "        select_Xg_field(field_name)"//achar(10)//&
            "        close_field_menu()"//achar(10)//&
            "        p.render()"//achar(10)//&
            "    def choose_colormap(colormap_name):"//achar(10)//&
            "        select_colormap(colormap_name)"//achar(10)//&
            "        close_colormap_menu()"//achar(10)//&
            "        p.render()"//achar(10)//&
            "    def toggle_field_menu(flag):"//achar(10)//&
            "        field_menu_open[0] = bool(flag)"//achar(10)//&
            "        if flag: close_colormap_menu()"//achar(10)//&
            "        set_menu_visibility(field_menu_widgets, field_menu_labels, bool(flag))"//achar(10)//&
            "        selected_name = 'Geometry' if active_field[0] is None else active_field[0]"//achar(10)//&
            "        marker = '^' if flag else 'v'"//achar(10)//&
            "        field_control_label.SetInput(f'Field: {selected_name}  {marker}')"//achar(10)//&
            "        p.render()"//achar(10)//&
            "    def toggle_colormap_menu(flag):"//achar(10)//&
            "        colormap_menu_open[0] = bool(flag)"//achar(10)//&
            "        if flag: close_field_menu()"//achar(10)//&
            "        set_menu_visibility(colormap_menu_widgets, colormap_menu_labels, bool(flag))"//achar(10)//&
            "        marker = '^' if flag else 'v'"//achar(10)//&
            "        colormap_control_label.SetInput(f'Colormap: {active_colormap[0]}  {marker}')"//achar(10)//&
            "        p.render()"//achar(10)//&
            "    controls_y = y_pos - (5 if actor_Xth else 4)*25"//achar(10)//&
            "    field_menu_button[0] = p.add_checkbox_button_widget(toggle_field_menu,"//&
            " value=False, position=(0, controls_y), size=16, border_size=2,"//&
            " color_on='black', color_off='grey', background_color='white')"//achar(10)//&
            "    field_control_label = p.add_text(f'Field: {field_names[0]}  v',"//&
            " position=(22, controls_y - 1), font_size=8, color='black', font='times')"//achar(10)//&
            "    colormap_menu_button[0] = p.add_checkbox_button_widget(toggle_colormap_menu,"//&
            " value=False, position=(0, controls_y - 22), size=16, border_size=2,"//&
            " color_on='black', color_off='grey', background_color='white')"//achar(10)//&
            "    colormap_control_label = p.add_text(f'Colormap: {colormap_names[0]}  v',"//&
            " position=(22, controls_y - 23), font_size=8, color='black', font='times')"//achar(10)//&
            "    menu_y = controls_y - 44"//achar(10)//&
            "    for i, field_name in enumerate(['Geometry'] + field_names):"//achar(10)//&
            "        option_y = menu_y - 18*i"//achar(10)//&
            "        widget = p.add_radio_button_widget(lambda name=field_name: choose_field(name),"//&
            " 'forcad_field', value=field_name == field_names[0], position=(4, option_y),"//&
            " size=12, border_size=3, color_on='cyan', color_off='grey', background_color='white')"//achar(10)//&
            "        label = p.add_text(field_name, position=(20, option_y - 1), font_size=7,"//&
            " color='black', font='times')"//achar(10)//&
            "        label.GetTextProperty().SetBackgroundColor(1.0, 1.0, 1.0)"//achar(10)//&
            "        label.GetTextProperty().SetBackgroundOpacity(0.9)"//achar(10)//&
            "        field_menu_widgets.append(widget)"//achar(10)//&
            "        field_menu_labels.append(label)"//achar(10)//&
            "    for i, colormap_name in enumerate(colormap_names):"//achar(10)//&
            "        option_y = menu_y - 18*i"//achar(10)//&
            "        widget = p.add_radio_button_widget("//&
            "lambda name=colormap_name: choose_colormap(name), 'forcad_colormap',"//&
            " value=i == 0, position=(4, option_y), size=12, border_size=3,"//&
            " color_on='black', color_off='grey', background_color='white')"//achar(10)//&
            "        label = p.add_text(colormap_name, position=(20, option_y - 1), font_size=7,"//&
            " color='black', font='times')"//achar(10)//&
            "        label.GetTextProperty().SetBackgroundColor(1.0, 1.0, 1.0)"//achar(10)//&
            "        label.GetTextProperty().SetBackgroundOpacity(0.9)"//achar(10)//&
            "        colormap_menu_widgets.append(widget)"//achar(10)//&
            "        colormap_menu_labels.append(label)"//achar(10)//&
            "    close_field_menu()"//achar(10)//&
            "    close_colormap_menu()"//achar(10)//&
            "p.add_text('ForCAD', position=(0.0, 10.0), font_size=14, color='black', font='times')"//achar(10)//&
            "p.add_text('https://github.com/gha3mi/forcad', position=(0.0, 0.0),"//&
            " font_size=7, color='blue', font='times')"//achar(10)//&
            "p.show(title='ForCAD', interactive=True)"//achar(10)//&
            "p.deep_clean()"//achar(10)//&
            "del p"
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Display previously exported single-patch VTK files with PyVista.
    !!
    !! This routine does not export geometry. It starts a synchronous Python
    !! process and returns after the interactive window closes. Defining
    !! `NOSHOW_PYVISTA` at compilation makes the routine a no-op.
    !!
    impure subroutine show_pyvista_singlepatch(vtkfile_Xc, vtkfile_Xg, vtkfile_Xth_in_Xg, rank_name)
        character(len=*), intent(in) :: vtkfile_Xc
            !! Existing control-geometry VTK path.
        character(len=*), intent(in) :: vtkfile_Xg
            !! Existing sampled-geometry VTK path.
        character(len=*), intent(in), optional :: vtkfile_Xth_in_Xg
            !! Optional parameter-line VTK path.
        character(len=*), intent(in), optional :: rank_name
            !! Optional viewer label.
#ifndef NOSHOW_PYVISTA
        block
        character(len=:), allocatable :: header

        header = "Xc_files = ["//pyvista_string_literal(vtkfile_Xc)//"]"//achar(10)//&
            "Xg_files = ["//pyvista_string_literal(vtkfile_Xg)//"]"//achar(10)
        if (present(vtkfile_Xth_in_Xg)) then
            header = header//"Xth_files = ["//pyvista_string_literal(vtkfile_Xth_in_Xg)//"]"//achar(10)
        else
            header = header//"Xth_files = []"//achar(10)
        end if
        if (present(rank_name)) then
            header = header//"rank_label = "//pyvista_string_literal(rank_name)//achar(10)
        else
            header = header//"rank_label = 'geometry'"//achar(10)
        end if

        call run_python_script(pyvista_viewer_script(header))
        end block
#endif
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Display previously exported multipatch VTK files using exact paths or glob patterns.
    !! Matching files are sorted before display so patch ordering is
    !! deterministic. This routine does not export geometry. Defining
    !! `NOSHOW_PYVISTA` at compilation makes it a no-op.
    !!
    impure subroutine show_pyvista_multipatch(vtkfile_Xc, vtkfile_Xg, vtkfile_Xth_in_Xg, rank_name)
        character(len=*), intent(in) :: vtkfile_Xc
            !! Control-geometry path or glob pattern.
        character(len=*), intent(in) :: vtkfile_Xg
            !! Sampled-geometry path or glob pattern.
        character(len=*), intent(in), optional :: vtkfile_Xth_in_Xg
            !! Optional parameter-line path or glob pattern.
        character(len=*), intent(in), optional :: rank_name
            !! Optional viewer label.
#ifndef NOSHOW_PYVISTA
        block
        character(len=:), allocatable :: header

        header = "import glob"//achar(10)//&
            "Xc_files = sorted(glob.glob("//pyvista_string_literal(vtkfile_Xc)//"))"//achar(10)//&
            "Xg_files = sorted(glob.glob("//pyvista_string_literal(vtkfile_Xg)//"))"//achar(10)
        if (present(vtkfile_Xth_in_Xg)) then
            header = header//"Xth_files = sorted(glob.glob("//&
                pyvista_string_literal(vtkfile_Xth_in_Xg)//"))"//achar(10)
        else
            header = header//"Xth_files = []"//achar(10)
        end if
        if (present(rank_name)) then
            header = header//"rank_label = "//pyvista_string_literal(rank_name)//achar(10)
        else
            header = header//"rank_label = 'multipatch'"//achar(10)
        end if

        call run_python_script(pyvista_viewer_script(header))
        end block
#endif
    end subroutine
    !===============================================================================
    !> Evaluate one mixed basis derivative, with local, dense or selected output.
    !! Source knots/counts/degrees and parameter/order/wrapping lengths must agree.
    !! Constructor padding avoids an ifx 2026.2 Release RESHAPE optimizer failure.
    !! Directional scratch occupies one local buffer; no global basis is allocated.
    pure subroutine tensor_basis_derivative_at(Xt,nc,degree,knot1,order,wrap,values,first_active,elem,Wc,knot2,knot3,side)
        real(rk), intent(in) :: Xt(:)
        integer, intent(in) :: nc(:), degree(:), order(:)
        real(rk), intent(in), contiguous :: knot1(:)
        real(rk), intent(in), contiguous, optional :: knot2(:), knot3(:), Wc(:)
        logical, intent(in) :: wrap(:)
        real(rk), intent(out) :: values(:)
        integer, intent(out), optional :: first_active(:)
        integer, intent(in), optional :: elem(:),side(:)
        type axis_derivatives
            real(rk), pointer, contiguous :: data(:,:)
        end type axis_derivatives
        type(axis_derivatives) :: axis(3)
        real(rk), target :: work(sum((max(0,order)+1)*(degree+1))+3-size(nc))
        real(rk) :: local_values(product(degree+1))
        integer :: counts(3), p(3), orders(3), first(3), offset(3), nd, d, start, n, i, index, global
        integer :: trace(size(nc))

        values=0.0_rk
        nd=size(nc)
        trace=0
        if (present(side)) trace=side
        counts=[nc,spread(1,1,3-nd)]
        p=[degree,spread(0,1,3-nd)]
        orders=[max(0,order),spread(0,1,3-nd)]
        first=1
        start=1
        work=1.0_rk
        do d=1,3
            n=(orders(d)+1)*(p(d)+1)
            axis(d)%data(0:orders(d),0:p(d))=>work(start:start+n-1)
            start=start+n
        end do
        call basis_bspline_der_all_active(Xt(1),knot1,nc(1),p(1),order(1),first(1),axis(1)%data,trace(1),wrap(1))
        if (nd>1) call basis_bspline_der_all_active(Xt(2),knot2,nc(2),p(2),order(2),first(2),axis(2)%data,trace(2),wrap(2))
        if (nd>2) call basis_bspline_der_all_active(Xt(3),knot3,nc(3),p(3),order(3),first(3),axis(3)%data,trace(3),wrap(3))
        call tensor_basis_derivative_local(first(1),counts(1),axis(1)%data,first(2),counts(2),axis(2)%data, &
            first(3),counts(3),axis(3)%data,local_values,Wc)
        if (present(first_active)) then
            if (size(first_active)>=nd) first_active(:nd)=first(:nd)
            n=min(size(values),size(local_values))
            values(:n)=local_values(:n)
            return
        end if
        if (present(elem)) then
            do i=1,min(size(values),size(elem))
                global=elem(i)
                if (global<1 .or. global>product(nc)) cycle
                offset=[mod(global-1,counts(1)),mod((global-1)/counts(1),counts(2)), &
                    (global-1)/(counts(1)*counts(2))]+1-first
                if (any(offset<0) .or. any(offset>p)) cycle
                index=1+offset(1)+(p(1)+1)*(offset(2)+(p(2)+1)*offset(3))
                values(i)=local_values(index)
            end do
        else
            do i=1,size(local_values)
                offset=[mod(i-1,p(1)+1),mod((i-1)/(p(1)+1),p(2)+1),(i-1)/((p(1)+1)*(p(2)+1))]+first
                if (any(offset<1) .or. any(offset>counts)) cycle
                global=offset(1)+counts(1)*(offset(2)-1+counts(2)*(offset(3)-1))
                if (global<=size(values)) values(global)=local_values(i)
            end do
        end if
    end subroutine tensor_basis_derivative_at

    !> One positive-order geometry derivative, with local span coordinates.
    !! Bernstein extraction normalizes each axis before the rational recurrence.
    !! Anchored contractions remove translations; one combined width scale
    !! restores the original parameter units, including arbitrary mixed orders.
    !! With row_power, the derivative in native units is Xg(row,:)*2**row_power(row).
    !! Any unrepresentable required workspace returns an empty result.
    !! Callers validate owned geometry, matching finite parameter axes, rank
    !! 1:3, and a nonnegative derivative multi-index with positive total order.
    pure subroutine geometry_derivative(Xt1,nc,degree,knot1,Xc,Wc,order,wrap,Xg,knot2,knot3,Xt2,Xt3,side,row_power)
        real(rk), intent(in) :: Xt1(:)
        real(rk), intent(in), optional :: Xt2(:),Xt3(:)
        integer, intent(in) :: nc(:),degree(:),order(:)
        integer, intent(in), optional :: side(:)
        real(rk), intent(in), contiguous :: knot1(:),Xc(:,:)
        real(rk), intent(in), contiguous, optional :: Wc(:),knot2(:),knot3(:)
        logical, intent(in) :: wrap(:)
        real(rk), allocatable, intent(out) :: Xg(:,:)
        integer, optional, allocatable, intent(out) :: row_power(:)
        type :: derivative_axis
            real(rk), allocatable :: knots(:),A(:,:),B(:,:),ders(:,:)
            integer :: span
        end type derivative_axis
        type(derivative_axis) :: axis(3)
        type(debug) :: status
        real(rk), allocatable :: result(:,:),basis(:),values(:),weights(:)
        integer, allocatable :: ids(:),powers(:)
        integer :: component_power(size(Xc,2))
        real(rk) :: width(3),multiplier,denominator,total,delta,term,value,factors(3)
        integer :: p(3),n(3),o(3),first(3),width_power(3),factor_power,trace(size(nc))
        integer :: ndim,nlocal,row,d,i,c,local,anchor,power,largest,pass,l1,l2,l3
        logical :: inside,weighted

        allocate(Xg(0,0))
        if (present(row_power)) then
            allocate(row_power(0),powers(size(Xt1)),source=0)
        end if
        ndim=size(nc)
        trace=0
        if (present(side)) trace=side
        axis%span=0
        p=[degree,spread(0,1,3-ndim)]
        n=[nc,spread(1,1,3-ndim)]
        o=[order,spread(0,1,3-ndim)]
        nlocal=product(p+1)
        allocate(result(size(Xt1),size(Xc,2)),source=0.0_rk)
        allocate(basis(0:nlocal-1),values(0:nlocal-1),weights(0:nlocal-1),ids(nlocal))
        do d=1,3
            axis(d)%knots=repelem([0.0_rk,1.0_rk],[p(d)+1,p(d)+1])
            allocate(axis(d)%B(0:o(d),0:p(d)),axis(d)%ders(0:o(d),0:p(d)))
            axis(d)%ders=1.0_rk
        end do
        do row=1,size(Xt1)
            first=1
            width=1.0_rk
            width_power=0
            do d=1,ndim
                select case(d)
                case(1)
                    call axis_derivatives(knot1,n(d),p(d),o(d),Xt1(row),wrap(d),trace(d),axis(d), &
                        first(d),width(d),width_power(d),inside,status)
                case(2)
                    call axis_derivatives(knot2,n(d),p(d),o(d),Xt2(row),wrap(d),trace(d),axis(d), &
                        first(d),width(d),width_power(d),inside,status)
                case default
                    call axis_derivatives(knot3,n(d),p(d),o(d),Xt3(row),wrap(d),trace(d),axis(d), &
                        first(d),width(d),width_power(d),inside,status)
                end select
                if (.not.status%ok) return
                if (.not.inside) exit
            end do
            if (.not.inside) cycle
            do local=0,nlocal-1
                l1=mod(local,p(1)+1)
                l2=mod(local/(p(1)+1),p(2)+1)
                l3=local/((p(1)+1)*(p(2)+1))
                ids(local+1)=first(1)+l1+n(1)*(first(2)+l2-1+n(2)*(first(3)+l3-1))
                basis(local)=axis(1)%ders(0,l1)*axis(2)%ders(0,l2)*axis(3)%ders(0,l3)
                factors=[minval(abs(axis(1)%ders(:,l1)),mask=axis(1)%ders(:,l1)/=0.0_rk), &
                    minval(abs(axis(2)%ders(:,l2)),mask=axis(2)%ders(:,l2)/=0.0_rk), &
                    minval(abs(axis(3)%ders(:,l3)),mask=axis(3)%ders(:,l3)/=0.0_rk)]
                values(local)=huge(1.0_rk)
                if (all(factors<huge(1.0_rk))) then
                    values(local)=min(factors(1),scaled_product(factors(:2)),scaled_product(factors))
                end if
            end do
            if (all(Xc(ids,:)==spread(Xc(ids(1),:),1,nlocal))) cycle
            weighted=.false.
            if (present(Wc)) then
                weights=Wc(ids)
                weighted=maxval(weights)/=minval(weights)
            end if
            if (.not.weighted .and. any(o>p)) cycle
            if (any(values==0.0_rk) .or. .not.all(ieee_is_finite(values))) return
            if (weighted) then
                power=maxval(exponent(weights))
                weights=scale(weights,-power)
                if (any(weights<=0.0_rk) .or. .not.all(ieee_is_finite(weights))) return
                if (any(scale(weights,power)/=Wc(ids))) return
                if (any(values/=huge(1.0_rk) .and. weights*values==0.0_rk)) return
                denominator=sum(weights*basis)
                if (denominator<=0.0_rk .or. .not.ieee_is_finite(denominator)) return
                anchor=ids(maxloc(weights*basis,dim=1))
                call tensor_basis_derivative_local(1,p(1)+1,axis(1)%ders, &
                    1,p(2)+1,axis(2)%ders,1,p(3)+1,axis(3)%ders,values,weights,checked=.true.)
            else
                anchor=ids(maxloc(basis,dim=1))
                call tensor_basis_derivative_local(1,p(1)+1,axis(1)%ders, &
                    1,p(2)+1,axis(2)%ders,1,p(3)+1,axis(3)%ders,values)
            end if
            if (.not.all(ieee_is_finite(values))) return
            multiplier=1.0_rk
            factor_power=0
            do d=1,ndim
                do i=1,o(d)
                    multiplier=multiplier/fraction(width(d))
                    factor_power=factor_power-width_power(d)-exponent(width(d))+exponent(multiplier)
                    multiplier=fraction(multiplier)
                end do
            end do
            component_power=0
            do c=1,size(Xc,2)
                largest=-huge(0)
                total=0.0_rk
                do pass=1,2
                    do i=0,nlocal-1
                        if (values(i)==0.0_rk) cycle
                        call scaled_difference(Xc(anchor,c),Xc(ids(i+1),c),delta,power)
                        if (delta==0.0_rk) cycle
                        power=power+exponent(delta)+exponent(values(i))
                        term=fraction(delta)*fraction(values(i))
                        if (pass==1) then
                            largest=max(largest,power)
                        else
                            value=scale(term,power-largest)
                            if (abs(value)<tiny(value)) return
                            total=total+value
                        end if
                    end do
                end do
                if (total==0.0_rk) cycle
                value=fraction(total)*multiplier
                component_power(c)=largest+exponent(total)+factor_power+exponent(value)
                result(row,c)=fraction(value)
            end do
            if (present(row_power)) then
                if (all(result(row,:)==0.0_rk)) cycle
                powers(row)=maxval(component_power,mask=result(row,:)/=0.0_rk)
                where (result(row,:)/=0.0_rk) component_power=component_power-powers(row)
            end if
            do c=1,size(Xc,2)
                if (result(row,c)==0.0_rk) cycle
                if (component_power(c)>maxexponent(value)) return
                result(row,c)=scale(result(row,c),component_power(c))
                if (result(row,c)==0.0_rk) return
            end do
        end do
        if (present(row_power)) call move_alloc(powers,row_power)
        call move_alloc(result,Xg)
    contains
        pure subroutine axis_derivatives(knot,n,p,o,Xt,wrapped,side,axis,first,width,power,inside,status)
            real(rk), intent(in), contiguous :: knot(:)
            integer, intent(in) :: n,p,o,side
            real(rk), intent(in) :: Xt
            logical, intent(in) :: wrapped
            type(derivative_axis), intent(inout) :: axis
            integer, intent(out) :: first,power
            real(rk), intent(out) :: width
            logical, intent(out) :: inside
            type(debug), intent(inout) :: status
            integer :: span,unused,delta_power
            real(rk) :: u,xi,delta

            call active_span_at(Xt,knot,n,p,span,u,inside,side,wrapped)
            if (.not.ieee_is_finite(u)) then
                call status%set(code=108,severity=1,category="forcad_utils",location="cmp_Xg", &
                    message="Parameter wrapping is not representable.",suggestion="Rescale the parameter domain.")
                return
            end if
            if (.not.inside) return
            first=span-p
            if (axis%span/=span) then
                call bezier_axis(knot,p,span,axis%A,status)
                if (.not.status%ok) return
                axis%span=span
            end if
            call scaled_difference(knot(span),knot(span+1),width,power)
            call scaled_difference(knot(span),u,delta,delta_power)
            xi=scale(fraction(delta)/fraction(width),delta_power-power+exponent(delta)-exponent(width))
            if (u==knot(span)) xi=0.0_rk
            if (u==knot(span+1)) xi=1.0_rk
            call basis_bspline_der_all_active(xi,axis%knots,p+1,p,o,unused,axis%B)
            axis%ders=matmul(axis%B,axis%A)
            if (.not.all(ieee_is_finite(axis%ders)) .or. &
                (u>knot(span) .and. u<knot(span+1) .and. any(axis%B(0,:)==0.0_rk))) call status%set( &
                code=108,severity=1,category="forcad_utils",location="cmp_Xg", &
                message="The local derivative workspace is not representable.",suggestion="Rescale the parameter domain.")
        end subroutine axis_derivatives
    end subroutine geometry_derivative

end module forcad_utils
