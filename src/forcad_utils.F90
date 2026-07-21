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

    implicit none

    private
    public basis_bernstein, basis_bspline, elemConn_C0, kron, ndgrid, compute_multiplicity, compute_knot_vector, &
        basis_bspline_der, basis_bspline_der_order, basis_bspline_der_order_active, basis_bspline_der_all_active, &
        tensor_basis_derivative_local, tensor_basis_derivatives2_local, &
        insert_knot_A_5_1, findspan, &
        elevate_degree_A_5_9, hexahedron_Xc, tetragon_Xc, remove_knots_A_5_8, &
        elemConn_Cn, unique, rotation, basis_bspline_2der, det, inv, dyad, gauss_leg, export_vtk_legacy, solve, &
        solve_spd_banded, sparse_left_matmul, repelem, linspace, fill_uniform, eye, kron_eye, &
        valid_knot_vector, knot_start, knot_end, active_knots, active_span_count, &
        active_knot_multiplicity, infer_knot_shape, knot_tolerance, map_parameter, &
        rotate_points, run_python_script, show_pyvista_singlepatch, show_pyvista_multipatch, boundary_index, boundary_layer_index, &
        disjoint_set_union, disjoint_set_map, structural_nonzero

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
        module procedure gauss_legendre_3D
    end interface
    !===============================================================================


    !===============================================================================
    !> Compute supported Kronecker products without calling `matmul`.
    !!
    !! Overloads cover vector-vector, vector-matrix, and a three-factor tensor
    !! product. Their flattening order is consistent with ForCAD control nets.
    interface kron
        module procedure kron_t1_t1
        module procedure kron_t1_t2
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
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Compute a tensor-product basis, gradient, and Hessian on local support.
    !! The derivative tables contain rows zero through two. Unused directions
    !! use a three-row, one-column identity table. Rational derivatives apply
    !! the quotient rule through second total order. This is a fused
    !! tensor-product specialization of the multivariate quotient recurrence,
    !! extending the rational-surface recurrence in Piegl and Tiller, *The
    !! NURBS Book*, Algorithm A4.4 to as many as three parametric directions.
    !! Direction 1 varies fastest in `basis`. `gradient(:,i)` stores
    !! \(\partial R_A/\partial\xi_i\), while `hessian(:,i,j)` stores the mixed
    !! derivative \(\partial^2R_A/(\partial\xi_i\partial\xi_j)\).
    !!
    pure subroutine tensor_basis_derivatives2_local(&
        first1,&
        nc1,&
        ders1,&
        first2,&
        nc2,&
        ders2,&
        first3,&
        nc3,&
        ders3,&
        basis,&
        gradient,&
        hessian,&
        Wc)
        integer, intent(in) :: first1
            !! First active control index in direction 1.
        integer, intent(in) :: nc1
            !! Total control count in direction 1.
        integer, intent(in) :: first2
            !! First active control index in direction 2.
        integer, intent(in) :: nc2
            !! Total control count in direction 2.
        integer, intent(in) :: first3
            !! First active control index in direction 3.
        integer, intent(in) :: nc3
            !! Total control count in direction 3.
        real(rk), intent(in), contiguous :: ders1(0:,0:)
            !! Derivative rows 0:2 for the active direction-1 basis.
        real(rk), intent(in), contiguous :: ders2(0:,0:)
            !! Derivative rows 0:2 for the active direction-2 basis.
        real(rk), intent(in), contiguous :: ders3(0:,0:)
            !! Derivative rows 0:2 for the active direction-3 basis.
        real(rk), intent(out), contiguous :: basis(0:)
            !! Local tensor-product values.
        real(rk), intent(out), contiguous :: gradient(0:,:)
            !! Local parametric gradients, shape `[nlocal,rank]`.
        real(rk), intent(out), contiguous :: hessian(0:,:,:)
            !! Local parametric Hessians, shape `[nlocal,rank,rank]`.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Optional flattened global weight array.
        real(rk) :: denominator, denominator1(3), denominator2(3,3)
        real(rk) :: value, derivative1(3), derivative2(3,3), wi
        integer :: p1, p2, p3, l1, l2, l3, local, global, i, j, ndim

        basis = 0.0_rk
        gradient = 0.0_rk
        hessian = 0.0_rk
        p1 = ubound(ders1,2)
        p2 = ubound(ders2,2)
        p3 = ubound(ders3,2)
        ndim = size(gradient,2)
        if (ubound(ders1,1) < 2 .or. ubound(ders2,1) < 2 .or. ubound(ders3,1) < 2) return
        if (nc1 < 1 .or. nc2 < 1 .or. nc3 < 1) return
        if (first1 < 1 .or. first2 < 1 .or. first3 < 1) return
        if (first1 + p1 > nc1 .or. first2 + p2 > nc2 .or. first3 + p3 > nc3) return
        if (size(basis) < (p1+1)*(p2+1)*(p3+1)) return
        if (ndim < 1 .or. ndim > 3 .or. size(gradient,1) < size(basis)) return
        if (size(hessian,1) < size(basis) .or. size(hessian,2) < ndim .or. size(hessian,3) < ndim) return

        denominator = 1.0_rk
        denominator1 = 0.0_rk
        denominator2 = 0.0_rk
        if (present(Wc)) then
            if (size(Wc) < nc1*nc2*nc3) return
            denominator = 0.0_rk
            do l3 = 0, p3
                do l2 = 0, p2
                    do l1 = 0, p1
                        global = first1 + l1 + nc1*((first2+l2-1) + nc2*(first3+l3-1))
                        wi = Wc(global)
                        value = ders1(0,l1)*ders2(0,l2)*ders3(0,l3)
                        derivative1(1) = ders1(1,l1)*ders2(0,l2)*ders3(0,l3)
                        derivative1(2) = ders1(0,l1)*ders2(1,l2)*ders3(0,l3)
                        derivative1(3) = ders1(0,l1)*ders2(0,l2)*ders3(1,l3)
                        derivative2(1,1) = ders1(2,l1)*ders2(0,l2)*ders3(0,l3)
                        derivative2(1,2) = ders1(1,l1)*ders2(1,l2)*ders3(0,l3)
                        derivative2(1,3) = ders1(1,l1)*ders2(0,l2)*ders3(1,l3)
                        derivative2(2,2) = ders1(0,l1)*ders2(2,l2)*ders3(0,l3)
                        derivative2(2,3) = ders1(0,l1)*ders2(1,l2)*ders3(1,l3)
                        derivative2(3,3) = ders1(0,l1)*ders2(0,l2)*ders3(2,l3)
                        denominator = denominator + wi*value
                        denominator1(1) = denominator1(1) + wi*derivative1(1)
                        denominator1(2) = denominator1(2) + wi*derivative1(2)
                        denominator1(3) = denominator1(3) + wi*derivative1(3)
                        denominator2(1,1) = denominator2(1,1) + wi*derivative2(1,1)
                        denominator2(1,2) = denominator2(1,2) + wi*derivative2(1,2)
                        denominator2(1,3) = denominator2(1,3) + wi*derivative2(1,3)
                        denominator2(2,2) = denominator2(2,2) + wi*derivative2(2,2)
                        denominator2(2,3) = denominator2(2,3) + wi*derivative2(2,3)
                        denominator2(3,3) = denominator2(3,3) + wi*derivative2(3,3)
                    end do
                end do
            end do
            denominator2(2,1) = denominator2(1,2)
            denominator2(3,1) = denominator2(1,3)
            denominator2(3,2) = denominator2(2,3)
            if (denominator <= 0.0_rk .or. .not. ieee_is_finite(denominator)) return
        end if

        do l3 = 0, p3
            do l2 = 0, p2
                do l1 = 0, p1
                    local = l1 + (p1+1)*(l2 + (p2+1)*l3)
                    wi = 1.0_rk
                    if (present(Wc)) then
                        global = first1 + l1 + nc1*((first2+l2-1) + nc2*(first3+l3-1))
                        wi = Wc(global)
                    end if
                    value = wi*ders1(0,l1)*ders2(0,l2)*ders3(0,l3)
                    derivative1(1) = wi*ders1(1,l1)*ders2(0,l2)*ders3(0,l3)
                    derivative1(2) = wi*ders1(0,l1)*ders2(1,l2)*ders3(0,l3)
                    derivative1(3) = wi*ders1(0,l1)*ders2(0,l2)*ders3(1,l3)
                    derivative2(1,1) = wi*ders1(2,l1)*ders2(0,l2)*ders3(0,l3)
                    derivative2(1,2) = wi*ders1(1,l1)*ders2(1,l2)*ders3(0,l3)
                    derivative2(1,3) = wi*ders1(1,l1)*ders2(0,l2)*ders3(1,l3)
                    derivative2(2,2) = wi*ders1(0,l1)*ders2(2,l2)*ders3(0,l3)
                    derivative2(2,3) = wi*ders1(0,l1)*ders2(1,l2)*ders3(1,l3)
                    derivative2(3,3) = wi*ders1(0,l1)*ders2(0,l2)*ders3(2,l3)
                    derivative2(2,1) = derivative2(1,2)
                    derivative2(3,1) = derivative2(1,3)
                    derivative2(3,2) = derivative2(2,3)
                    basis(local) = value/denominator
                    do i = 1, ndim
                        gradient(local,i) = (derivative1(i)-basis(local)*denominator1(i))/denominator
                    end do
                    do j = 1, ndim
                        do i = 1, ndim
                            hessian(local,i,j) = (derivative2(i,j)-basis(local)*denominator2(i,j) - &
                                gradient(local,i)*denominator1(j)-gradient(local,j)*denominator1(i))/denominator
                        end do
                    end do
                end do
            end do
        end do
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
    pure subroutine disjoint_set_map(parent, map)
        integer, intent(inout), contiguous :: parent(:)
            !! Parent forest; paths are shortened in place.
        integer, allocatable, intent(out) :: map(:)
            !! Component identifier for every input member.
        integer :: i, r, n

        allocate(map(size(parent)), source=0)
        n = 0
        do i = 1, size(parent)
            call disjoint_set_root(parent, i, r)
            if (map(r) == 0) then
                n = n + 1
                map(r) = n
            end if
            map(i) = map(r)
        end do
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

        allocate(idx(1))
        select case(side)
        case(1)
            idx(1) = 1
        case(2)
            idx(1) = nc
        case default
            idx(1) = 0
        end select
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
        integer :: i, p, q, n1, n2

        n1 = nc(1)
        n2 = nc(2)
        select case(side)
        case(1)
            allocate(idx(n2))
            do p = 1, n2
                q = merge(n2-p+1, p, reverse)
                idx(p) = (q-1)*n1 + 1
            end do
        case(2)
            allocate(idx(n2))
            do p = 1, n2
                q = merge(n2-p+1, p, reverse)
                idx(p) = (q-1)*n1 + n1
            end do
        case(3)
            allocate(idx(n1))
            do p = 1, n1
                i = merge(n1-p+1, p, reverse)
                idx(p) = i
            end do
        case(4)
            allocate(idx(n1))
            do p = 1, n1
                i = merge(n1-p+1, p, reverse)
                idx(p) = (n2-1)*n1 + i
            end do
        case default
            allocate(idx(0))
        end select
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
        integer :: a, b, aa, bb, i, j, k, m1, m2, out1, out2, p, q

        select case(side)
        case(1, 2)
            m1 = nc(2)
            m2 = nc(3)
            i = merge(1, nc(1), side == 1)
        case(3, 4)
            m1 = nc(1)
            m2 = nc(3)
            j = merge(1, nc(2), side == 3)
        case(5, 6)
            m1 = nc(1)
            m2 = nc(2)
            k = merge(1, nc(3), side == 5)
        case default
            allocate(idx(0))
            return
        end select

        out1 = merge(m2, m1, swap)
        out2 = merge(m1, m2, swap)
        allocate(idx(out1*out2))
        p = 0
        do b = 1, out2
            do a = 1, out1
                if (swap) then
                    aa = b
                    bb = a
                else
                    aa = a
                    bb = b
                end if
                if (size(flip) >= 1) then
                    if (flip(1)) aa = m1-aa+1
                end if
                if (size(flip) >= 2) then
                    if (flip(2)) bb = m2-bb+1
                end if
                p = p + 1
                q = merge(size(idx)-p+1, p, reverse)
                select case(side)
                case(1, 2)
                    idx(q) = (bb-1)*nc(1)*nc(2) + (aa-1)*nc(1) + i
                case(3, 4)
                    idx(q) = (bb-1)*nc(1)*nc(2) + (j-1)*nc(1) + aa
                case(5, 6)
                    idx(q) = (k-1)*nc(1)*nc(2) + (bb-1)*nc(1) + aa
                case default
                    idx(q) = 0
                end select
            end do
        end do
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
        integer :: i, j, p, q, n1, n2

        n1 = nc(1)
        n2 = nc(2)
        select case(side)
        case(1)
            allocate(idx(n2))
            i = 1 + layer
            do p = 1, n2
                q = merge(n2-p+1, p, reverse)
                idx(p) = (q-1)*n1 + i
            end do
        case(2)
            allocate(idx(n2))
            i = n1 - layer
            do p = 1, n2
                q = merge(n2-p+1, p, reverse)
                idx(p) = (q-1)*n1 + i
            end do
        case(3)
            allocate(idx(n1))
            j = 1 + layer
            do p = 1, n1
                i = merge(n1-p+1, p, reverse)
                idx(p) = (j-1)*n1 + i
            end do
        case(4)
            allocate(idx(n1))
            j = n2 - layer
            do p = 1, n1
                i = merge(n1-p+1, p, reverse)
                idx(p) = (j-1)*n1 + i
            end do
        case default
            allocate(idx(0))
        end select
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
        integer :: a, b, aa, bb, i, j, k, m1, m2, out1, out2, p, q

        select case(side)
        case(1, 2)
            m1 = nc(2)
            m2 = nc(3)
            i = merge(1 + layer, nc(1) - layer, side == 1)
        case(3, 4)
            m1 = nc(1)
            m2 = nc(3)
            j = merge(1 + layer, nc(2) - layer, side == 3)
        case(5, 6)
            m1 = nc(1)
            m2 = nc(2)
            k = merge(1 + layer, nc(3) - layer, side == 5)
        case default
            allocate(idx(0))
            return
        end select

        out1 = merge(m2, m1, swap)
        out2 = merge(m1, m2, swap)
        allocate(idx(out1*out2))
        p = 0
        do b = 1, out2
            do a = 1, out1
                if (swap) then
                    aa = b
                    bb = a
                else
                    aa = a
                    bb = b
                end if
                if (size(flip) >= 1) then
                    if (flip(1)) aa = m1-aa+1
                end if
                if (size(flip) >= 2) then
                    if (flip(2)) bb = m2-bb+1
                end if
                p = p + 1
                q = merge(size(idx)-p+1, p, reverse)
                select case(side)
                case(1, 2)
                    idx(q) = (bb-1)*nc(1)*nc(2) + (aa-1)*nc(1) + i
                case(3, 4)
                    idx(q) = (bb-1)*nc(1)*nc(2) + (j-1)*nc(1) + aa
                case(5, 6)
                    idx(q) = (k-1)*nc(1)*nc(2) + (bb-1)*nc(1) + aa
                case default
                    idx(q) = 0
                end select
            end do
        end do
    end subroutine
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
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Optionally wrap one spline parameter into its active interval.
    !! Disabled wrapping and invalid spline states leave the parameter unchanged.
    !! Wrapping maps into the half-open interval \([a,b)\), where
    !! \(a=U_p\) and \(b=U_{n_c}\). It changes parameter evaluation only; it
    !! does not alter the knot vector or assert periodic geometry closure.
    !!
    pure function map_parameter(Xt, knot, nc, degree, wrap_parameter) result(u)
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
        real(rk) :: u
        real(rk) :: a, b, period

        u = Xt
        if (.not. wrap_parameter .or. .not. ieee_is_finite(Xt)) return
        if (degree < 0 .or. nc < 1 .or. degree >= nc) return
        if (size(knot) < nc + degree + 1) return

        a = knot(degree+1)
        b = knot(nc+1)
        period = b - a
        if (.not. ieee_is_finite(period) .or. period <= 0.0_rk) return

        if (Xt >= a .and. Xt < b) return

        u = a + modulo(Xt - a, period)
        if (u >= b) u = a
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Find the active 1-based span for a valid local spline space.
    !! The active parameter interval is [knot(degree+1), knot(nc+1)].
    pure subroutine active_span_at(Xt, knot, nc, degree, span, u, inside)
        real(rk), intent(in) :: Xt
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: nc, degree
        integer, intent(out) :: span
        real(rk), intent(out) :: u
        logical, intent(out) :: inside
        integer :: low, mid, high
        real(rk) :: a, b, tol

        span = 0
        u = Xt
        inside = .false.
        if (degree < 0 .or. nc < 1 .or. degree >= nc) return
        if (size(knot) < nc + degree + 1) return

        a = knot(degree+1)
        b = knot(nc+1)
        tol = knot_tolerance(knot, degree+1, nc+1)
        if (u < a) then
            if (a - u > tol) return
            u = a
        else if (u > b) then
            if (u - b > tol) return
            u = b
        end if

        if (abs(u - b) <= tol) then
            span = nc
            inside = .true.
            return
        end if

        low = degree + 1
        high = nc
        do while (low <= high)
            mid = (low + high)/2
            if (u >= knot(mid) .and. u < knot(mid+1)) then
                span = mid
                inside = .true.
                return
            else if (u < knot(mid)) then
                high = mid - 1
            else
                low = mid + 1
            end if
        end do
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
    !! periodic knot-spacing extension, periodic control-data closure, or
    !! periodic geometry.
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
    !> Evaluate the at-most-\(p+1\) nonzero B-spline values at one parameter.
    !! `N(j)` corresponds to one-based global basis index `first+j`.
    pure subroutine basis_bspline_active(Xt, knot, nc, degree, first, N)
        real(rk), intent(in) :: Xt
        real(rk), intent(in), contiguous :: knot(:)
        integer, intent(in) :: nc, degree
        integer, intent(out) :: first
        real(rk), intent(out) :: N(0:degree)
        call basis_bspline_active_ders(Xt, knot, nc, degree, 0, first, N=N)
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
        real(rk) :: N(0:degree)
        integer :: first, i, idx

        B = 0.0_rk
        if (nc == 0) return
        call basis_bspline_active(Xt, knot, nc, degree, first, N)
        do i = 0, degree
            idx = first + i
            if (idx >= 1 .and. idx <= nc) B(idx) = N(i)
        end do
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
    pure subroutine basis_bspline_der_all_active(Xt, knot, nc, degree, nder, first, ders)
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
        call active_span_at(Xt, knot, nc, degree, span, u, inside)
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
    !! `values` uses direction-1-fastest local ordering.
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
        Wc)
        integer, intent(in) :: first1, nc1, first2, nc2, first3, nc3
        real(rk), intent(in), contiguous :: ders1(0:,0:), ders2(0:,0:), ders3(0:,0:)
        real(rk), intent(out), contiguous :: values(0:)
            !! Local derivative value for every supported basis function.
        real(rk), intent(in), contiguous, optional :: Wc(:)
            !! Optional flattened global rational weights; absence selects the polynomial tensor-product path.
        integer :: p1, p2, p3, o1, o2, o3
        integer :: a1, a2, a3, b1, b2, b3, l1, l2, l3, local, global
        real(rk) :: denominator(0:ubound(ders1,1),0:ubound(ders2,1),0:ubound(ders3,1))
        real(rk) :: quotient(0:ubound(ders1,1),0:ubound(ders2,1),0:ubound(ders3,1))
        real(rk) :: coefficient1, coefficient2, coefficient3, numerator

        values = 0.0_rk
        p1 = ubound(ders1,2)
        p2 = ubound(ders2,2)
        p3 = ubound(ders3,2)
        o1 = ubound(ders1,1)
        o2 = ubound(ders2,1)
        o3 = ubound(ders3,1)
        if (nc1 < 1 .or. nc2 < 1 .or. nc3 < 1) return
        if (first1 < 1 .or. first2 < 1 .or. first3 < 1) return
        if (first1 + p1 > nc1 .or. first2 + p2 > nc2 .or. first3 + p3 > nc3) return
        if (size(values) < (p1+1)*(p2+1)*(p3+1)) return

        if (.not. present(Wc)) then
            do concurrent (l3 = 0:p3, l2 = 0:p2, l1 = 0:p1)
                values(l1+(p1+1)*(l2+(p2+1)*l3)) = ders1(o1,l1)*ders2(o2,l2)*ders3(o3,l3)
            end do
            return
        end if
        if (size(Wc) < nc1*nc2*nc3) return

        denominator = 0.0_rk
        do a3 = 0, o3
            do a2 = 0, o2
                do a1 = 0, o1
                    do l3 = 0, p3
                        do l2 = 0, p2
                            do l1 = 0, p1
                                global = first1 + l1 + nc1*((first2+l2-1) + nc2*(first3+l3-1))
                                denominator(a1,a2,a3) = denominator(a1,a2,a3) + Wc(global)*&
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
                                numerator = Wc(global)*ders1(a1,l1)*ders2(a2,l2)*ders3(a3,l3)
                                coefficient3 = 1.0_rk
                                do b3 = 0, a3
                                    coefficient2 = 1.0_rk
                                    do b2 = 0, a2
                                        coefficient1 = 1.0_rk
                                        do b1 = 0, a1
                                            if (b1 + b2 + b3 > 0) then
                                                numerator = numerator - coefficient1*coefficient2*coefficient3*&
                                                    denominator(b1,b2,b3)*quotient(a1-b1,a2-b2,a3-b3)
                                            end if
                                            if (b1 < a1) then
                                                coefficient1 = coefficient1*real(a1-b1,rk)/real(b1+1,rk)
                                            end if
                                        end do
                                        if (b2 < a2) then
                                            coefficient2 = coefficient2*real(a2-b2,rk)/real(b2+1,rk)
                                        end if
                                    end do
                                    if (b3 < a3) then
                                        coefficient3 = coefficient3*real(a3-b3,rk)/real(b3+1,rk)
                                    end if
                                end do
                                quotient(a1,a2,a3) = numerator/denominator(0,0,0)
                            end do
                        end do
                    end do
                    local = l1 + (p1+1)*(l2 + (p2+1)*l3)
                    values(local) = quotient(o1,o2,o3)
                end do
            end do
        end do
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
        integer :: nelem1, nelem2
        integer :: e1, e2, i, j

        if (p1 < 1 .or. p2 < 1 .or. nnode1 < 1 .or. nnode2 < 1) then
            allocate(elemConn(0,0))
            return
        end if
        if (mod(nnode1-1, p1) /= 0 .or. mod(nnode2-1, p2) /= 0) then
            allocate(elemConn(0,0))
            return
        end if

        nelem1 = (nnode1-1)/p1
        nelem2 = (nnode2-1)/p2

        allocate(elemConn(nelem1*nelem2, (p1+1)*(p2+1)))

        do concurrent (e2 = 1:nelem2, e1 = 1:nelem1, j = 0:p2, i = 0:p1)
            elemConn((e2-1)*nelem1+e1, j*(p1+1)+(i+1)) = ((1+(e2-1)*p2+j)-1)*nnode1 + (1+(e1-1)*p1+i)
        end do
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
        integer :: nelem1, nelem2, nelem3, e1, e2, e3, i, j, k

        if (p1 < 1 .or. p2 < 1 .or. p3 < 1 .or. nnode1 < 1 .or. nnode2 < 1 .or. nnode3 < 1) then
            allocate(elemConn(0,0))
            return
        end if
        if (mod(nnode1-1,p1) /= 0 .or. mod(nnode2-1,p2) /= 0 .or. mod(nnode3-1,p3) /= 0) then
            allocate(elemConn(0,0))
            return
        end if

        nelem1 = (nnode1-1)/p1
        nelem2 = (nnode2-1)/p2
        nelem3 = (nnode3-1)/p3

        allocate(elemConn(nelem1*nelem2*nelem3, (p1+1)*(p2+1)*(p3+1)))

        do concurrent (e3 = 1:nelem3, e2 = 1:nelem2, e1 = 1:nelem1, k = 0:p3, j = 0:p2, i = 0:p1)
            elemConn((e3-1)*nelem1*nelem2+(e2-1)*nelem1+e1, k*(p2+1)*(p1+1)+j*(p1+1)+(i+1))&
                = ((1+(e3-1)*p3+k)-1)*nnode1*nnode2 + ((1+(e2-1)*p2+j)-1)*nnode1 + (1+(e1-1)*p1+i)
        end do
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

        integer, allocatable :: pref1(:), pref2(:)
        integer :: nelem1, nelem2, i, j, i2, j2

        nelem1 = size(Xth1)-1
        nelem2 = size(Xth2)-1

        if (nnode1 < 1 .or. nnode2 < 1 .or. p1 < 0 .or. p2 < 0 .or. &
            nelem1 < 1 .or. nelem2 < 1 .or. size(vecKnot_mul1) < nelem1 .or. &
            size(vecKnot_mul2) < nelem2) then
            allocate(elemConn(0,0))
            return
        end if

        allocate(pref1(nelem1), pref2(nelem2))
        pref1(1) = p1 + 1
        do i = 2, nelem1
            pref1(i) = pref1(i-1) + vecKnot_mul1(i)
        end do
        pref2(1) = p2 + 1
        do j = 2, nelem2
            pref2(j) = pref2(j-1) + vecKnot_mul2(j)
        end do

        allocate(elemConn(nelem1*nelem2, (p1+1)*(p2+1)))
        do concurrent (j = 1:nelem2, i = 1:nelem1, j2 = 0:p2, i2 = 0:p1)
            elemConn((j-1)*nelem1+i, j2*(p1+1)+(i2+1)) = ((pref2(j)-p2+j2)-1)*nnode1+(pref1(i)-p1+i2)
        end do
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

        integer, allocatable :: pref1(:), pref2(:), pref3(:)
        integer :: nelem1, nelem2, nelem3, i, j, k, i2, j2, k2

        nelem1 = size(Xth1)-1
        nelem2 = size(Xth2)-1
        nelem3 = size(Xth3)-1

        if (nnode1 < 1 .or. nnode2 < 1 .or. nnode3 < 1 .or. p1 < 0 .or. p2 < 0 .or. p3 < 0 .or. &
            nelem1 < 1 .or. nelem2 < 1 .or. nelem3 < 1 .or. size(vecKnot_mul1) < nelem1 .or. &
            size(vecKnot_mul2) < nelem2 .or. size(vecKnot_mul3) < nelem3) then
            allocate(elemConn(0,0))
            return
        end if

        allocate(pref1(nelem1), pref2(nelem2), pref3(nelem3))
        pref1(1) = p1 + 1
        do i = 2, nelem1
            pref1(i) = pref1(i-1) + vecKnot_mul1(i)
        end do
        pref2(1) = p2 + 1
        do j = 2, nelem2
            pref2(j) = pref2(j-1) + vecKnot_mul2(j)
        end do
        pref3(1) = p3 + 1
        do k = 2, nelem3
            pref3(k) = pref3(k-1) + vecKnot_mul3(k)
        end do

        allocate(elemConn(nelem1*nelem2*nelem3, (p1+1)*(p2+1)*(p3+1)))
        do concurrent (k = 1:nelem3, j = 1:nelem2, i = 1:nelem1, k2 = 0:p3, j2 = 0:p2, i2 = 0:p1)
            elemConn((k-1)*nelem1*nelem2+(j-1)*nelem1+i, k2*(p2+1)*(p1+1)+j2*(p1+1)+(i2+1)) =&
            ((pref3(k)-p3+k2)-1)*nnode1*nnode2+((pref2(j)-p2+j2)-1)*nnode1+(pref1(i)-p1+i2)
        end do
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
        integer  :: i, j, L, mp, d, np, Lf, bw, nq_temp, k_full, shared

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

        if (s+r == p+1) then
            if (r > 1) then
                if (present(T)) then
                    call insert_knot_A_5_1(&
                        p  = p,&
                        UP = UP,&
                        Pw = Pw,&
                        u  = u,&
                        k  = k,&
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
                        k  = k,&
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
            k_full = k + r - 1
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

        UQ(0:k)           = UP(0:k)
        UQ(k+1:k+r)       = u
        UQ(k+1+r:mp+r)    = UP(k+1:mp)
        Qw(0:k-p, :)      = Pw(0:k-p, :)
        Qw(k-s+r:np+r, :) = Pw(k-s:np, :)
        Rw(0:p-s, :)      = Pw(k-p:k-s, :)
        if (present(T)) then
            allocate(T(0:nq, 0:np), source=0.0_rk)
            do concurrent (i = 0:k-p)
                T(i, i) = 1.0_rk
            end do
            do concurrent (i = k-s+r:nq)
                T(i, i-r) = 1.0_rk
            end do
            Lf = k-p+r
            bw = p-s
            do concurrent (i = 0:bw)
                T(Lf+i, k-p+i) = 1.0_rk
            end do
        end if
        do j = 1, r
            L = k-p+j
            do i = 0, p-j-s
                alpha = (u-UP(L+i))/(UP(i+k+1)-UP(L+i))
                Rw(i, :) = alpha*Rw(i+1, :)+(1.0_rk-alpha)*Rw(i, :)
                if (present(T)) then
                    T(Lf+i, :) = alpha*T(Lf+i+1, :)+(1.0_rk-alpha)*T(Lf+i, :)
                end if
            end do
            Qw(L, :)       = Rw(0, :)
            Qw(k+r-j-s, :) = Rw(p-j-s, :)
            if (present(T)) then
                T(L, :)       = T(Lf+0, :)
                T(k+r-j-s, :) = T(Lf+(p-j-s), :)
            end if
        end do
        Qw(L+1:k-s-1, :) = Rw(1:k-s-1-L, :)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Locate the zero-based knot span containing a parameter.
    !!
    !! Implements the binary search of Piegl and Tiller, Algorithm A2.1, with
    !! endpoint saturation. Here `n` is the highest zero-based control index. A
    !! parameter at or below `knot(degree+1)+eps` returns `degree`; one at or
    !! above `knot(n+2)-eps` returns `n`, including values far outside the active
    !! interval. This routine therefore does not reject out-of-domain
    !! parameters. Invalid metadata returns the sentinel zero; otherwise
    !! `degree<=s<=n`.
    pure function findspan(n,degree,Xth,knot) result(s)
        integer, intent(in) :: n, degree
        real(rk), intent(in) :: Xth
        real(rk), intent(in), contiguous :: knot(:)
        integer :: s
        integer :: low, high, mid
        real(rk) :: eps

        s = 0
        if (degree < 0 .or. n < degree .or. size(knot) < max(n + 2, n + degree + 1)) return

        eps = knot_tolerance(knot, degree+1, n+2)
        if (Xth <= knot(degree+1) + eps) then
            s = degree
            return
        end if
        if (Xth >= knot(n+2) - eps) then
            s = n
            return
        end if
        low = degree
        high = n + 1
        mid = (low + high) / 2
        do while (Xth < knot(mid+1) .or. Xth >= knot(mid+2))
            if (Xth < knot(mid+1)) then
                high = mid
            else
                low = mid
            end if
            mid = (low + high) / 2
        end do
        s = mid
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Elevate a spline degree using Algorithm A5.9 of Piegl and Tiller.
    !! Homogeneous coordinates must be supplied for rational geometry. For
    !! `t>0`, an open-clamped input is elevated directly. Any other valid knot
    !! representation is first restricted to an equivalent open-clamped
    !! representation on its active interval; exterior knots and control
    !! extensions are therefore not retained. The returned spline represents
    !! the same mapping on that interval at degree `degree+t`, and each retained
    !! knot multiplicity is raised consistently. This operation does not impose
    !! or preserve periodic seam constraints. `t=0` returns the original arrays
    !! and identity map. A negative elevation or invalid space returns empty
    !! arrays and `success=.false.` when present.
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

        real(rk), allocatable :: Xcw_work(:,:), Xcw_inserted(:,:), Xcw_clamped(:,:)
        real(rk), allocatable :: knot_work(:), knot_inserted(:), knot_clamped(:)
        real(rk), allocatable :: Twork(:,:), Tinsert(:,:), Tcomposed(:,:), Tclamped(:,:), Telev(:,:)
        real(rk) :: u0, u1, u_boundary
        integer :: nc, i, boundary, multiplicity, add, span, nq
        integer :: first_start, last_start, first_end, last_end, start_multiplicity, end_multiplicity
        integer :: first_control, last_control, first_row, last_row, pos
        logical :: open_knot

        if (present(success)) success = .false.
        nc = size(Xcw, 1)
        nc_new = 0
        if (t < 0 .or. degree < 0 .or. nc <= degree .or. size(Xcw,2) < 1 .or. &
            .not. valid_knot_vector(knot, nc, degree) .or. .not. all(ieee_is_finite(Xcw))) then
            allocate(knot_new(0), Xcw_new(0,size(Xcw,2)))
            if (present(Tmap)) allocate(Tmap(0,0))
            return
        end if
        if (t == 0) then
            nc_new = nc
            knot_new = knot
            Xcw_new = Xcw
            if (present(Tmap)) then
                allocate(Tmap(nc, nc), source=0.0_rk)
                do concurrent (i = 1:nc)
                    Tmap(i, i) = 1.0_rk
                end do
            end if
            if (present(success)) success = .true.
            return
        end if
        open_knot = all(knot(1:degree+1) == knot(degree+1)) .and. &
            all(knot(nc+1:nc+degree+1) == knot(nc+1))

        if (open_knot) then
            call elevate_degree_open_A_5_9(&
                t        = t,&
                knot     = knot,&
                degree   = degree,&
                Xcw      = Xcw,&
                nc_new   = nc_new,&
                knot_new = knot_new,&
                Xcw_new  = Xcw_new,&
                Tmap     = Tmap)
            if (present(success)) then
                success = nc_new > 0 .and. size(knot_new) == nc_new + degree + t + 1 .and. &
                    size(Xcw_new,1) == nc_new .and. all(ieee_is_finite(Xcw_new))
            end if
            return
        end if

        u0 = knot_start(knot, nc, degree)
        u1 = knot_end(knot, nc, degree)
        knot_work = knot
        Xcw_work = Xcw
        if (present(Tmap)) then
            allocate(Twork(nc,nc), source=0.0_rk)
            do concurrent (i = 1:nc)
                Twork(i,i) = 1.0_rk
            end do
        end if

        do boundary = 1, 2
            if (boundary == 1) then
                u_boundary = u0
            else
                u_boundary = u1
            end if

            multiplicity = compute_multiplicity(knot_work, u_boundary)
            add = max(0, degree - multiplicity)
            if (add <= 0) cycle

            span = -1
            do i = lbound(knot_work,1), ubound(knot_work,1)
                if (knot_work(i) == u_boundary) span = i - lbound(knot_work,1)
            end do
            if (span < 0) then
                allocate(knot_new(0), Xcw_new(0,size(Xcw,2)))
                if (present(Tmap)) allocate(Tmap(0,0))
                return
            end if

            if (present(Tmap)) then
                call insert_knot_A_5_1(&
                    p  = degree,&
                    UP = knot_work,&
                    Pw = Xcw_work,&
                    u  = u_boundary,&
                    k  = span,&
                    s  = multiplicity,&
                    r  = add,&
                    nq = nq,&
                    UQ = knot_inserted,&
                    Qw = Xcw_inserted,&
                    T  = Tinsert)
                call sparse_left_matmul(&
                    A = Tinsert,&
                    B = Twork,&
                    C = Tcomposed)
                call move_alloc(Tcomposed, Twork)
            else
                call insert_knot_A_5_1(&
                    p  = degree,&
                    UP = knot_work,&
                    Pw = Xcw_work,&
                    u  = u_boundary,&
                    k  = span,&
                    s  = multiplicity,&
                    r  = add,&
                    nq = nq,&
                    UQ = knot_inserted,&
                    Qw = Xcw_inserted)
            end if
            call move_alloc(knot_inserted, knot_work)
            call move_alloc(Xcw_inserted, Xcw_work)
        end do

        first_start = lbound(knot_work,1) - 1
        last_start = lbound(knot_work,1) - 1
        first_end = lbound(knot_work,1) - 1
        last_end = lbound(knot_work,1) - 1
        do i = lbound(knot_work,1), ubound(knot_work,1)
            if (knot_work(i) == u0) then
                if (first_start < lbound(knot_work,1)) first_start = i
                last_start = i
            end if
            if (knot_work(i) == u1) then
                if (first_end < lbound(knot_work,1)) first_end = i
                last_end = i
            end if
        end do
        start_multiplicity = last_start - first_start + 1
        end_multiplicity = last_end - first_end + 1
        first_control = last_start - degree
        last_control = first_end - 1
        if (first_start < lbound(knot_work,1) .or. first_end < lbound(knot_work,1) .or. &
            start_multiplicity < degree .or. end_multiplicity < degree .or. &
            first_control < lbound(Xcw_work,1) .or. last_control > ubound(Xcw_work,1)) then
            allocate(knot_new(0), Xcw_new(0,size(Xcw,2)))
            if (present(Tmap)) allocate(Tmap(0,0))
            return
        end if

        pos = last_end - first_start + 1
        if (start_multiplicity == degree) pos = pos + 1
        if (end_multiplicity == degree) pos = pos + 1
        allocate(knot_clamped(pos))
        pos = 1
        if (start_multiplicity == degree) then
            knot_clamped(pos) = u0
            pos = pos + 1
        end if
        knot_clamped(pos:pos+last_end-first_start) = knot_work(first_start:last_end)
        pos = pos + last_end - first_start + 1
        if (end_multiplicity == degree) knot_clamped(pos) = u1
        Xcw_clamped = Xcw_work(first_control:last_control,:)
        if (size(knot_clamped) /= size(Xcw_clamped,1) + degree + 1) then
            nc_new = 0
            knot_new = [real(rk) ::]
            allocate(Xcw_new(0,size(Xcw,2)))
            if (present(Tmap)) allocate(Tmap(0,0))
            return
        end if
        if (present(Tmap)) then
            first_row = first_control - lbound(Xcw_work,1) + 1
            last_row = last_control - lbound(Xcw_work,1) + 1
            Tclamped = Twork(first_row:last_row,:)
            call elevate_degree_open_A_5_9(&
                t        = t,&
                knot     = knot_clamped,&
                degree   = degree,&
                Xcw      = Xcw_clamped,&
                nc_new   = nc_new,&
                knot_new = knot_new,&
                Xcw_new  = Xcw_new,&
                Tmap     = Telev)
            call sparse_left_matmul(&
                A = Telev,&
                B = Tclamped,&
                C = Tmap)
        else
            call elevate_degree_open_A_5_9(&
                t        = t,&
                knot     = knot_clamped,&
                degree   = degree,&
                Xcw      = Xcw_clamped,&
                nc_new   = nc_new,&
                knot_new = knot_new,&
                Xcw_new  = Xcw_new)
        end if
        if (present(success)) then
            success = nc_new > 0 .and. size(knot_new) == nc_new + degree + t + 1 .and. &
                size(Xcw_new,1) == nc_new .and. all(ieee_is_finite(Xcw_new))
        end if
    end subroutine
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
    !! |a-b|\le128\,\epsilon_{rk}\max(|a|,|b|).
    !! \]
    !!
    !! Polynomial callers pass control coordinates. Rational callers must pass
    !! homogeneous controls \((w\mathbf P,w)\); the resulting test is then in
    !! homogeneous coordinate space and is not a certified Euclidean geometry
    !! error bound. If no removal is accepted, `t=0` and both outputs are copies
    !! of the inputs.
    !!
    pure recursive subroutine remove_knots_A_5_8(&
        p,&
        knot,&
        Pw,&
        u,&
        r,&
        s,&
        num,&
        t,&
        knot_new,&
        Pw_new)
        real(rk), intent(in) :: u
            !! Knot value requested for removal.
        integer, intent(in) :: p
            !! Polynomial degree.
        integer, intent(in) :: r
            !! One-based index of the final current occurrence of `u`.
        integer, intent(in) :: s
            !! Complete current multiplicity of `u`.
        integer, intent(in) :: num
            !! Maximum number of removals requested.
        real(rk), intent(in), contiguous :: knot(:)
            !! Original valid knot vector.
        real(rk), intent(in), contiguous :: Pw(:,:)
            !! Original homogeneous or polynomial control variables.
        real(rk), allocatable, intent(out) :: knot_new(:)
            !! Resulting knot vector.
        real(rk), allocatable, intent(out) :: Pw_new(:,:)
            !! Resulting control variables.
        real(rk), allocatable :: Pw_copy(:,:), knot_copy(:)
        integer, intent(out) :: t
            !! Number of removals accepted.
        real(rk) :: alfi, alfj, candidate, denominator_i, denominator_j, tol
        real(rk), allocatable :: temp(:,:)
        integer :: i, j, ii, jj, remflag, off, first, last, ord, fout, m, k, n, nc, d, tt, t_more, c

        d  = size(Pw,2)
        nc = size(Pw,1)
        t = 0
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

            do c = 1, d
                tol = 128.0_rk*epsilon(1.0_rk)*max(abs(Pw(first-1,c)), abs(Pw(first,c)))
                if (abs(Pw(first-1,c) - Pw(first,c)) > tol) return
            end do

            deallocate(Pw_new, knot_new)
            allocate(Pw_new(nc-1,d))
            Pw_new(1:first-2,:) = Pw(1:first-2,:)
            Pw_new(first-1,:) = 0.5_rk*(Pw(first-1,:) + Pw(first,:))
            Pw_new(first:nc-1,:) = Pw(first+1:nc,:)

            allocate(knot_new(size(knot)-1))
            knot_new(1:r-1) = knot(1:r-1)
            knot_new(r:size(knot)-1) = knot(r+1:size(knot))
            t = 1

            if (num > 1) then
                call remove_knots_A_5_8(&
                    p        = p,&
                    knot     = knot_new,&
                    Pw       = Pw_new,&
                    u        = u,&
                    r        = r-1,&
                    s        = s-1,&
                    num      = num-1,&
                    t        = t_more,&
                    knot_new = knot_copy,&
                    Pw_new   = Pw_copy)
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
            if (j-i<=t) then
                remflag = 1
                do c = 1, d
                    tol = 128.0_rk*epsilon(1.0_rk)*max(abs(temp(ii,c)), abs(temp(jj+2,c)))
                    if (abs(temp(ii,c) - temp(jj+2,c)) > tol) then
                        remflag = 0
                        exit
                    end if
                end do
                if (remflag == 0) then
                    denominator_i = knot_copy(i+ord+t) - knot_copy(i)
                    if (denominator_i /= 0.0_rk) then
                        alfi = (u-knot_copy(i))/denominator_i
                        remflag = 1
                        do c = 1, d
                            candidate = alfi*temp(ii+t+2,c) + (1.0_rk-alfi)*temp(ii,c)
                            tol = 128.0_rk*epsilon(1.0_rk)*max(abs(Pw_copy(i,c)), abs(candidate))
                            if (abs(Pw_copy(i,c) - candidate) > tol) then
                                remflag = 0
                                exit
                            end if
                        end do
                    end if
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
        real(rk), intent(in), contiguous :: A(:,:)
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
    !! Square orders one through three use explicit formulas; larger square
    !! matrices delegate to [[solve]]. For a full-column-rank \(m\times n\)
    !! matrix with \(m>n\), the rectangular path computes
    !! \((A^TA)^{-1}A^T\). For full row rank with \(m<n\), it computes
    !! \(A^T(AA^T)^{-1}\). These are the Moore-Penrose inverse only under the
    !! stated rank assumptions and acceptance of the internal pivot thresholds.
    !! Empty or numerically rejected input returns `A_inv(0,0)`. All entries of
    !! `A` must be finite; nonfinite input is not reliably rejected.
    !!
    !! @warning This routine is intended for geometry Jacobians, not general
    !! large least-squares problems. Normal equations square the condition
    !! number. Rejection uses the heuristic
    !! `32*epsilon(rk)*max(1,maxval(abs(A)))`; the explicit 2-by-2 and 3-by-3
    !! branches compare their determinant directly with that threshold, so the
    !! test is not a rank-revealing factorization or invariant under arbitrary
    !! rescaling.
    !! @endwarning
    !!
    !! **Reference:** R. Penrose, "A Generalized Inverse for Matrices,"
    !! *Proceedings of the Cambridge Philosophical Society* 51 (1955),
    !! 406--413.
    !! [doi:10.1017/S0305004100030401](https://doi.org/10.1017/S0305004100030401).
    recursive pure function inv(A) result(A_inv)
        real(rk), intent(in), contiguous :: A(:,:)
        real(rk), allocatable :: A_inv(:,:)
        integer :: m, n
        real(rk) :: detA, tol
        real(rk), allocatable :: A_transpose(:,:), normal_inv(:,:)

        m = size(A,1)
        n = size(A,2)
        if (m < 1 .or. n < 1) then
            allocate(A_inv(0,0))
            return
        end if
        tol = 32.0_rk*epsilon(1.0_rk)*max(1.0_rk, maxval(abs(A)))

        if (m == n) then
            select case(m)
              case(1)
                allocate(A_inv(1,1))
                if (abs(A(1,1)) <= tol) then
                    deallocate(A_inv)
                    allocate(A_inv(0,0))
                    return
                end if
                A_inv(1,1) = 1.0_rk/A(1,1)
              case(2)
                allocate(A_inv(m,n))
                detA = det(A)
                if (abs(detA) <= tol) then
                    deallocate(A_inv)
                    allocate(A_inv(0,0))
                    return
                end if
                A_inv(1,1) =  A(2,2)
                A_inv(1,2) = -A(1,2)
                A_inv(2,1) = -A(2,1)
                A_inv(2,2) =  A(1,1)
                A_inv = A_inv/detA
              case(3)
                allocate(A_inv(m,n))
                detA = det(A)
                if (abs(detA) <= tol) then
                    deallocate(A_inv)
                    allocate(A_inv(0,0))
                    return
                end if
                A_inv(1,1) = A(2,2)*A(3,3) - A(2,3)*A(3,2)
                A_inv(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
                A_inv(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)
                A_inv(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
                A_inv(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
                A_inv(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)
                A_inv(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
                A_inv(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
                A_inv(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)
                A_inv = A_inv/detA
              case default
                A_inv = solve(A,eye(m))
            end select
        else if (m>n) then
            A_transpose = transpose(A)
            normal_inv = inv(matmul(A_transpose, A))
            if (size(normal_inv,1) == 0) then
                allocate(A_inv(0,0))
                return
            end if
            A_inv = matmul(normal_inv, A_transpose)
        else if (m<n) then
            A_transpose = transpose(A)
            normal_inv = inv(matmul(A, A_transpose))
            if (size(normal_inv,1) == 0) then
                allocate(A_inv(0,0))
                return
            end if
            A_inv = matmul(A_transpose, normal_inv)
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
        real(rk), allocatable :: Xksi1(:), Wksi1(:), Xksi2(:), Wksi2(:)
        integer :: i, j, n1, n2

        n1 = degree(1)+1
        n2 = degree(2)+1
        allocate(Xksi1(n1), Wksi1(n1))
        allocate(Xksi2(n2), Wksi2(n2))

        call gauss_legendre(Xksi1, Wksi1, interval1)
        call gauss_legendre(Xksi2, Wksi2, interval2)

        call ndgrid(Xksi1, Xksi2, Xksi)
        allocate(Wksi(n1*n2))
        do concurrent (j = 1:n2, i = 1:n1)
            Wksi((j-1)*n1+i) = Wksi1(i)*Wksi2(j)
        end do
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
        real(rk), allocatable :: Xksi1(:), Wksi1(:), Xksi2(:), Wksi2(:), Xksi3(:), Wksi3(:)
        integer :: i, j, k, n1, n2, n3

        n1 = degree(1)+1
        n2 = degree(2)+1
        n3 = degree(3)+1
        allocate(Xksi1(n1), Wksi1(n1))
        allocate(Xksi2(n2), Wksi2(n2))
        allocate(Xksi3(n3), Wksi3(n3))

        call gauss_legendre(Xksi1, Wksi1, interval1)
        call gauss_legendre(Xksi2, Wksi2, interval2)
        call gauss_legendre(Xksi3, Wksi3, interval3)

        call ndgrid(Xksi1, Xksi2, Xksi3, Xksi)
        allocate(Wksi(n1*n2*n3))
        do concurrent (k = 1:n3, j = 1:n2, i = 1:n1)
            Wksi(((k-1)*n2+j-1)*n1+i) = Wksi1(i)*Wksi2(j)*Wksi3(k)
        end do
    end subroutine
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
    !! Scalar fields are emitted only when both `point_data` and `field_names`
    !! are present. The first `size(point_data,2)` names are written verbatim
    !! after trimming and must already be valid VTK scalar identifiers. If only
    !! one optional argument is supplied, all point fields are omitted.
    !!
    !! Binary output converts coordinates and fields to IEEE `real64` and
    !! explicitly serializes big-endian 64-bit reals and 32-bit integers as
    !! required by the legacy format. It requires 8-bit file storage units and
    !! the corresponding intrinsic integer and IEEE real kinds. Recognized
    !! encodings are exactly the lower-, title-, and upper-case forms of
    !! `ascii` and `binary`; binary is the default.
    !!
    !! The procedure has no status result. Rejected metadata returns before
    !! writing, and file-system I/O failures are not converted to a ForCAD
    !! diagnostic.
    !!
    !! Format reference: [VTK Legacy File Format](https://docs.vtk.org/en/latest/vtk_file_formats/vtk_legacy_file_format.html).
    impure subroutine export_vtk_legacy(filename, points, elemConn, vtkCellType, point_data, field_names, encoding)
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

        integer :: byte_count, i, j, ne, np, nn, n, nunit
        integer(int8) :: byte_buffer(8192)
        character(len=6) :: encoding_
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

            close(nunit)
        end if


        if (trim(encoding_) == "BINARY") then
            if (file_storage_size /= 8 .or. bit_size(0_int8) /= 8 .or. bit_size(0_int32) /= 32 .or. &
                bit_size(0_int64) /= 64 .or. storage_size(0.0_real64) /= 64 .or. &
                .not. ieee_support_datatype(0.0_real64)) return

            open(newunit=nunit, file=filename, form="formatted", action="write")
            write(nunit,"(a)") "# vtk DataFile Version 2.0"
            write(nunit,"(a)") "Generated by ForCAD"
            write(nunit,"(a)") "BINARY"
            if (is_polyline) then
                write(nunit,"(a)") "DATASET POLYDATA"
            else
                write(nunit,"(a)") "DATASET UNSTRUCTURED_GRID"
            end if
            close(nunit)

            open(newunit=nunit, file=filename, form="formatted", action="write", position="append")
            write(nunit,"(a,1x,g0,1x,a)") "POINTS", np, "double"
            close(nunit)
            open(&
                newunit  = nunit,&
                file     = filename,&
                position = "append",&
                access   = "stream",&
                form     = "unformatted",&
                action   = "write",&
                status   = "unknown")
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
            close(nunit)

            if (is_polyline) then
                open(newunit=nunit, file=filename, form="formatted", action="write", position="append")
                write(nunit,"(a,1x,g0,1x,g0)") "LINES", ne, n
                close(nunit)
                open(&
                    newunit  = nunit,&
                    file     = filename,&
                    position = "append",&
                    access   = "stream",&
                    form     = "unformatted",&
                    action   = "write",&
                    status   = "unknown")
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
                close(nunit)
            else
                open(newunit=nunit, file=filename, form="formatted", action="write", position="append")
                write(nunit,"(a,1x,g0,1x,g0)") "CELLS", ne, n
                close(nunit)
                open(&
                    newunit  = nunit,&
                    file     = filename,&
                    position = "append",&
                    access   = "stream",&
                    form     = "unformatted",&
                    action   = "write",&
                    status   = "unknown")
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
                close(nunit)

                open(newunit=nunit, file=filename, form="formatted", action="write", position="append")
                write(nunit,"(a,1x,g0)") "CELL_TYPES", ne
                close(nunit)
                open(&
                    newunit  = nunit,&
                    file     = filename,&
                    position = "append",&
                    access   = "stream",&
                    form     = "unformatted",&
                    action   = "write",&
                    status   = "unknown")
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
                close(nunit)
            end if

            if (present(point_data) .and. present(field_names)) then
                open(newunit=nunit, file=filename, form="formatted", action="write", position="append")
                write(nunit, "(a,1x,g0)") "POINT_DATA", size(point_data,1)
                close(nunit)

                do i = 1, size(point_data,2)
                    open(newunit=nunit, file=filename, form="formatted", action="write", position="append")
                    write(nunit, "(a,1x,a,1x,a)") "SCALARS", trim(field_names(i)), "double"
                    write(nunit, "(a)") "LOOKUP_TABLE default"
                    close(nunit)

                    open(&
                        newunit  = nunit,&
                        file     = filename,&
                        position = "append",&
                        access   = "stream",&
                        form     = "unformatted",&
                        action   = "write",&
                        status   = "unknown")
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
                    close(nunit)
                end do
            end if

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
    !! A matrix is classified as symmetric when every compared pair satisfies
    !!
    !! \[
    !! |A_{ij}-A_{ji}|\le
    !! 32\,\epsilon_{rk}\max(1,\max_{kl}|A_{kl}|).
    !! \]
    !!
    !! That path first attempts unpivoted Cholesky using the lower triangle; any
    !! upper/lower discrepancies within the tolerance are therefore represented
    !! by the lower triangle. If Cholesky encounters a nonpositive pivot, or if
    !! the symmetry test fails, partial-pivoting Gaussian elimination solves the
    !! original matrix. Shape mismatch or a pivot rejected by the LU threshold
    !! returns `X(0,0)`.
    !!
    !! The routine is allocation-based and intended for small or medium dense
    !! systems. Work is \(O(n^3+n^2n_{rhs})\), and dense factor and solution
    !! storage is allocated internally. The routine does not estimate a
    !! condition number or perform iterative refinement.
    !! All entries of `A` and `B` must be finite. Nonfinite data are not
    !! reliably detected and can return a nonempty matrix containing NaNs.
    pure function solve(A, B) result(X)
        real(rk), intent(in), contiguous :: A(:,:)
            !! Square coefficient matrix, shape `[n,n]`.
        real(rk), intent(in), contiguous :: B(:,:)
            !! Right-hand sides, shape `[n,nrhs]`.
        real(rk), allocatable :: X(:,:)

        integer :: n, m, i, j, p
        real(rk) :: matrix_scale, sym_tol
        logical :: symmetric, ok

        p = size(A,1)
        n = size(A,2)
        m = size(B,2)

        if (p /= n .or. p /= size(B,1) .or. n < 1) then
            allocate(X(0,0))
            return
        end if

        matrix_scale = max(1.0_rk, maxval(abs(A)))
        sym_tol = 32.0_rk*epsilon(1.0_rk)*matrix_scale
        symmetric = .true.
        do i = 2, n
            do j = 1, i-1
                if (abs(A(i,j)-A(j,i)) > sym_tol) then
                    symmetric = .false.
                    exit
                end if
            end do
            if (.not. symmetric) exit
        end do

        if (symmetric) then
            call solve_cholesky(A, B, X, ok)
            if (ok) return
        end if

        call solve_lu(A, B, X)
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
    !! \(O(n(kd+1)+n\,n_{rhs})\), avoiding a dense \(n\times n\) matrix. Invalid
    !! shape or a computed pivot not greater than
    !! `32*epsilon(rk)*max(1,maxval(abs(A_band)))` returns `X(0,0)`.
    !! The stored band and right-hand sides must be finite; nonfinite data are
    !! not reliably rejected.
    !!
    pure function solve_spd_banded(A_band, B) result(X)
        real(rk), intent(in), contiguous :: A_band(0:,:)
            !! Lower symmetric band storage, lower bound zero in dimension 1.
        real(rk), intent(in), contiguous :: B(:,:)
            !! Right-hand sides, shape `[n,nrhs]`.
        real(rk), allocatable :: X(:,:)
        real(rk), allocatable :: L(:,:)
        integer :: kd, n, m, i, j, k, r, k0, i1
        real(rk) :: acc, tol

        kd = size(A_band, 1) - 1
        n = size(A_band, 2)
        m = size(B, 2)
        if (size(B,1) /= n .or. n < 1 .or. kd < 0) then
            allocate(X(0,0))
            return
        end if

        tol = 32.0_rk*epsilon(1.0_rk)*max(1.0_rk, maxval(abs(A_band)))
        allocate(L(0:kd,n), X(n,m), source=0.0_rk)
        X = B

        do j = 1, n
            acc = A_band(0,j)
            k0 = max(1, j-kd)
            do k = k0, j-1
                acc = acc - L(j-k,k)*L(j-k,k)
            end do
            if (acc <= tol) then
                deallocate(X)
                allocate(X(0,0))
                return
            end if
            L(0,j) = sqrt(acc)

            do i = j+1, min(n, j+kd)
                acc = A_band(i-j,j)
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
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Solve a dense symmetric positive-definite system by Cholesky factorization.
    !!
    !! This private path factors \(A=LL^T\) without pivoting, then performs
    !! forward and backward substitutions for every right-hand side. A
    !! nonpositive pivot sets `ok=.false.`; `X` remains allocated but must not
    !! be used unless `ok` is true. Finite inputs are a caller precondition.
    !!
    !! **First published account:** Commandant Benoit, "Note sur une methode de
    !! resolution des equations normales provenant de l'application de la
    !! methode des moindres carres a un systeme d'equations lineaires en nombre
    !! inferieur a celui des inconnues -- Application de la methode a la
    !! resolution d'un systeme defini d'equations lineaires (procede du
    !! Commandant Cholesky)," *Bulletin Geodesique* 2 (1924), 67--77.
    !! [doi:10.1007/BF03031308](https://doi.org/10.1007/BF03031308).
    pure subroutine solve_cholesky(A, B, X, ok)
        real(rk), intent(in), contiguous :: A(:,:)
            !! Symmetric coefficient matrix `[n,n]`.
        real(rk), intent(in), contiguous :: B(:,:)
            !! Right-hand sides `[n,nrhs]`.
        real(rk), allocatable, intent(out) :: X(:,:)
            !! Solution matrix `[n,nrhs]` when `ok` is true.
        logical, intent(out) :: ok
            !! True only when all Cholesky pivots are positive.
        integer :: n, m, i, j, k
        real(rk), allocatable :: L(:,:)
        real(rk) :: acc

        n = size(A,2)
        m = size(B,2)
        ok = .true.
        allocate(L(n,n), X(n,m), source=0.0_rk)
        X = B

        do i = 1, n
            do j = 1, i
                acc = A(i,j)
                do k = 1, j-1
                    acc = acc - L(i,k)*L(j,k)
                end do
                if (i == j) then
                    if (acc <= 0.0_rk) then
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
    !! column. Invalid shapes or a pivot not greater than
    !! `32*epsilon(rk)*max(1,maxval(abs(A)))` return `X(0,0)`. This threshold has
    !! an absolute unit-scale floor; it is an internal numerical rejection rule,
    !! not a condition-number estimate or a scale-invariant rank test. Finite
    !! matrix and right-hand-side entries are a caller precondition; NaNs are
    !! not reliably rejected by pivot comparisons.
    !!
    !! **Reference:** J. von Neumann and H. H. Goldstine, "Numerical Inverting
    !! of Matrices of High Order," *Bulletin of the American Mathematical
    !! Society* 53 (1947), 1021--1099.
    !! [doi:10.1090/S0002-9904-1947-08909-6](https://doi.org/10.1090/S0002-9904-1947-08909-6).
    pure subroutine solve_lu(A, B, X)
        real(rk), intent(in), contiguous :: A(:,:)
            !! Square coefficient matrix `[n,n]`.
        real(rk), intent(in), contiguous :: B(:,:)
            !! Right-hand sides `[n,nrhs]`.
        real(rk), allocatable, intent(out) :: X(:,:)
            !! Solution matrix, or an empty matrix on failure.
        integer :: n, m, i, j, k, pivot
        real(rk), allocatable :: U(:,:)
        real(rk) :: factor, pivot_abs, candidate_abs, tmp, tol, acc

        n = size(A,2)
        m = size(B,2)
        if (n < 1 .or. size(A,1) /= n .or. size(B,1) /= n) then
            allocate(X(0,0))
            return
        end if
        tol = 32.0_rk*epsilon(1.0_rk)*max(1.0_rk, maxval(abs(A)))
        allocate(U(n,n), X(n,m))
        U = A
        X = B

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
            if (pivot_abs <= tol) then
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

        if (abs(U(n,n)) <= tol) then
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
end module forcad_utils
