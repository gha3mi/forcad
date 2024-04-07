!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> This module defines the 'nurbs_curve' type for representing a Non-Uniform Rational B-Spline (NURBS) curve.
module forcad_nurbs_curve

    use forcad_utils, only: rk, basis_bspline, elemConn_C0, compute_multiplicity, compute_knot_vector, basis_bspline_der,&
        insert_knot_A_5_1, findspan, elevate_degree_A

    implicit none

    private
    public nurbs_curve

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    type nurbs_curve
        real(rk), allocatable, private :: Xc(:,:) !! Control points (2D array: [nc, dim])
        real(rk), allocatable, private :: Xg(:,:) !! Geometry points (2D array: [ng, dim])
        real(rk), allocatable, private :: Wc(:)   !! Weights for control points (1D array: [nc])
        real(rk), allocatable, private :: Xt(:)   !! Evaluation points (1D array: [ng])
        real(rk), allocatable, private :: knot(:) !! Knot vector (1D array)
        integer, private :: degree                !! Degree (order) of the curve
        integer, private :: nc                    !! Number of control points
        integer, private :: ng                    !! Number of geometry points
    contains
        procedure :: set1                  !!> Set knot vector, control points and weights for the NURBS curve object
        procedure :: set2                  !!> Set NURBS curve using nodes of parameter space, degree, continuity, control points and weights
        procedure :: set3                  !!> Set Bezier or Rational Bezier curve using control points and weights
        generic :: set => set1, set2, set3 !!> Set NURBS curve
        procedure :: create                !!> Generate geometry points
        procedure :: get_Xc                !!> Get control points
        procedure :: get_Xg                !!> Get geometry points
        procedure :: get_Wc                !!> Get weights
        procedure :: get_Xt                !!> Get parameter values
        procedure :: get_knot              !!> Get knot vector
        procedure :: get_ng                !!> Get number of geometry points
        procedure :: get_order             !!> Get degree of the NURBS curve
        procedure :: finalize              !!> Finalize the NURBS curve object
        procedure :: get_elem_Xc           !!> Generate connectivity for control points
        procedure :: get_elem_Xg           !!> Generate connectivity for geometry points
        procedure :: export_Xc             !!> Export control points to VTK file
        procedure :: export_Xg             !!> Export geometry points to VTK file
        procedure :: modify_Xc             !!> Modify control points
        procedure :: modify_Wc             !!> Modify weights
        procedure :: get_multiplicity      !!> Get multiplicity of the knot vector
        procedure :: get_continuity        !!> Get continuity of the curve
        procedure :: get_nc                !!> Get number of required control points
        procedure :: insert_knots          !!> Insert knots into the knot vector
        procedure :: elevate_degree        !!> Elevate the degree of the curve
        procedure :: derivative            !!> Compute the derivative of the NURBS curve
        procedure :: basis                 !!> Compute the basis functions of the NURBS curve
    end type
    !===============================================================================

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Set knot vector, control points and weights for the NURBS curve object.
    pure subroutine set1(this, knot, Xc, Wc)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: knot(:)
        real(rk), intent(in) :: Xc(:,:)
        real(rk), intent(in), optional :: Wc(:)

        if (allocated(this%knot)) deallocate(this%knot)
        if (allocated(this%Xc)) deallocate(this%Xc)

        this%knot = knot
        this%degree = this%get_order()
        this%Xc = Xc
        this%nc = size(this%Xc, 1)
        if (present(Wc)) then
            if (size(Wc) /= this%nc) then
                error stop 'Number of weights does not match the number of control points.'
            else
                if (allocated(this%Wc)) deallocate(this%Wc)
                this%Wc = Wc
            end if
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Set NURBS curve using nodes of parameter space (Xth), degree, continuity, control points and weights.
    pure subroutine set2(this, Xth_dir, degree, continuity, Xc, Wc)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: Xth_dir(:)
        integer, intent(in) :: degree
        integer, intent(in) :: continuity(:)
        real(rk), intent(in) :: Xc(:,:)
        real(rk), intent(in), optional :: Wc(:)

        this%knot = compute_knot_vector(Xth_dir, degree, continuity)
        this%degree = degree
        this%Xc = Xc
        this%nc = size(this%Xc, 1)
        if (present(Wc)) then
            if (size(Wc) /= this%nc) then
                error stop 'Number of weights does not match the number of control points.'
            else
                this%Wc = Wc
            end if
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Set Bezier or Rational Bezier curve using control points and weights.
    pure subroutine set3(this, Xc, Wc)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: Xc(:,:)
        real(rk), intent(in), optional :: Wc(:)

        if (allocated(this%knot)) deallocate(this%knot)
        if (allocated(this%Xc)) deallocate(this%Xc)

        this%Xc = Xc
        this%nc = size(this%Xc, 1)

        allocate(this%knot(2*this%nc))
        this%knot(1:this%nc) = 0.0_rk
        this%knot(this%nc+1:2*this%nc) = 1.0_rk

        this%degree = this%get_order()
        if (present(Wc)) then
            if (size(Wc) /= this%nc) then
                error stop 'Number of weights does not match the number of control points.'
            else
                if (allocated(this%Wc)) deallocate(this%Wc)
                this%Wc = Wc
            end if
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine create(this, res, Xt)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in), optional :: res
        real(rk), intent(in), optional :: Xt(:)
        real(rk), allocatable :: Tgc(:)
        integer :: i, j

        ! check
        if (.not.allocated(this%Xc)) then
            error stop 'Control points are not set.'
        end if

        ! Set parameter values
        if (present(Xt)) then
            if (allocated(this%Xt)) deallocate(this%Xt)
            this%Xt = Xt
        elseif (present(res)) then
            if (allocated(this%Xt)) deallocate(this%Xt)
            allocate(this%Xt(res))
            this%Xt = [(real(i-1, rk) / real(res-1, rk), i=1, res)]
            ! else
            ! this%Xt = this%Xt
        end if

        ! Set number of geometry points
        this%ng = size(this%Xt)

        ! Allocate memory for geometry points
        if (allocated(this%Xg)) deallocate(this%Xg)
        allocate(this%Xg(this%ng, size(this%Xc,2)))

        if (allocated(this%Wc)) then
            do i = 1, size(this%Xt, 1)
                Tgc = basis_bspline(this%Xt(i), this%knot, this%nc, this%degree)
                Tgc = Tgc*(this%Wc/(dot_product(Tgc,this%Wc)))
                do j = 1, size(this%Xc, 2)
                    this%Xg(i,j) = dot_product(Tgc,this%Xc(:,j))
                end do
            end do
        else
            do i = 1, size(this%Xt, 1)
                Tgc = basis_bspline(this%Xt(i), this%knot, this%nc, this%degree)
                do j = 1, size(this%Xc, 2)
                    this%Xg(i,j) = dot_product(Tgc,this%Xc(:,j))
                end do
            end do
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_Xc(this) result(Xc)
        class(nurbs_curve), intent(in) :: this
        real(rk), allocatable :: Xc(:,:)

        if (allocated(this%Xc)) then
            Xc = this%Xc
        else
            error stop 'Control points are not set.'
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_Xg(this) result(Xg)
        class(nurbs_curve), intent(in) :: this
        real(rk), allocatable :: Xg(:,:)

        if (allocated(this%Xg)) then
            Xg = this%Xg
        else
            error stop 'Geometry points are not set.'
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_Wc(this) result(Wc)
        class(nurbs_curve), intent(in) :: this
        real(rk), allocatable :: Wc(:)

        if (allocated(this%Wc)) then
            Wc = this%Wc
        else
            error stop 'The NURBS curve is not rational.'
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_Xt(this) result(Xt)
        class(nurbs_curve), intent(in) :: this
        real(rk), allocatable :: Xt(:)

        if (allocated(this%Xt)) then
            Xt = this%Xt
        else
            error stop 'Parameter values are not set.'
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_ng(this) result(ng)
        class(nurbs_curve), intent(in) :: this
        integer :: ng

        ng = this%ng
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_order(this) result(degree)
        class(nurbs_curve), intent(in) :: this
        integer :: degree
        integer, allocatable :: m(:)

        m = this%get_multiplicity()

        degree = m(1) - 1
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_knot(this) result(knot)
        class(nurbs_curve), intent(in) :: this
        real(rk), allocatable :: knot(:)

        if (allocated(this%knot)) then
            knot = this%knot
        else
            error stop 'Knot vector is not set.'
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine finalize(this)
        class(nurbs_curve), intent(inout) :: this
        if (allocated(this%Xc)) deallocate(this%Xc)
        if (allocated(this%Xg)) deallocate(this%Xg)
        if (allocated(this%Wc)) deallocate(this%Wc)
        if (allocated(this%Xt)) deallocate(this%Xt)
        if (allocated(this%knot)) deallocate(this%knot)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine get_elem_Xc(this, elemConn, p)
        class(nurbs_curve), intent(in) :: this
        integer, dimension(:,:), allocatable, intent(out) :: elemConn
        integer, intent(in), optional :: p

        if (present(p)) then
            elemConn = elemConn_C0(this%nc,p)
        else
            elemConn = elemConn_C0(this%nc,1)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine get_elem_Xg(this, elemConn, p)
        class(nurbs_curve), intent(in) :: this
        integer, dimension(:,:), allocatable, intent(out) :: elemConn
        integer, intent(in), optional :: p

        if (present(p)) then
            elemConn = elemConn_C0(this%ng,p)
        else
            elemConn = elemConn_C0(this%ng,1)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_Xc(this, filename)
        class(nurbs_curve), intent(in) :: this
        character(len=*), intent(in) :: filename
        integer :: i, nc, nunit
        integer, dimension(:,:), allocatable :: elemConn

        ! check
        if (.not.allocated(this%Xc)) then
            error stop 'Control points are not set.'
        end if

        call this%get_elem_Xc(elemConn)

        nc = size(this%Xc, 1)

        open(newunit=nunit, file=filename, action='write')
        write(nunit,'(a)') '# vtk DataFile Version 2.0'
        write(nunit,'(a)') 'Generated by ForCAD'
        write(nunit,'(a)') 'ASCII'
        write(nunit,'(a)') 'DATASET UNSTRUCTURED_GRID'
        write(nunit,'(a," ",g0," ",a)') 'POINTS', nc, 'double'

        if (size(this%Xc, 2) == 2) then
            write(nunit,'(g0," ",g0," ",g0)') (this%Xc(i,1), this%Xc(i,2), 0.0_rk , i = 1, nc)
        elseif (size(this%Xc, 2) == 3) then
            write(nunit,'(g0," ",g0," ",g0)') (this%Xc(i,1), this%Xc(i,2), this%Xc(i,3) , i = 1, nc)
        else
            error stop 'Invalid dimension of the control points.'
        end if

        write(nunit,'(a," ",g0," ",g0)') 'CELLS', size(elemConn,1), size(elemConn,1)*(size(elemConn,2)+1)
        write(nunit,'(g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0)')&
            (2, elemConn(i,1)-1,elemConn(i,2)-1, i = 1, size(elemConn,1))

        write(nunit,'(a," ",g0)') 'CELL_TYPES', size(elemConn,1)
        write(nunit,'(g0)') (3 , i = 1, size(elemConn,1))
        close(nunit)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_Xg(this, filename)
        class(nurbs_curve), intent(in) :: this
        character(len=*), intent(in) :: filename
        integer :: i, ng, nunit
        integer, dimension(:,:), allocatable :: elemConn

        ! check
        if (.not.allocated(this%Xg)) then
            error stop 'Geometry points are not set.'
        end if

        call this%get_elem_Xg(elemConn)

        ng = size(this%Xg, 1)

        open(newunit=nunit, file=filename, action='write')
        write(nunit,'(a)') '# vtk DataFile Version 2.0'
        write(nunit,'(a)') 'Generated by ForCAD'
        write(nunit,'(a)') 'ASCII'
        write(nunit,'(a)') 'DATASET UNSTRUCTURED_GRID'
        write(nunit,'(a," ",g0," ",a)') 'POINTS', ng, 'double'

        if (size(this%Xg, 2) == 2) then
            write(nunit,'(g0," ",g0," ",g0)') (this%Xg(i,1), this%Xg(i,2), 0.0_rk , i = 1, ng)
        elseif (size(this%Xg, 2) == 3) then
            write(nunit,'(g0," ",g0," ",g0)') (this%Xg(i,1), this%Xg(i,2), this%Xg(i,3) , i = 1, ng)
        else
            error stop 'Invalid dimension of the geometry points.'
        end if

        write(nunit,'(a," ",g0," ",g0)') 'CELLS', size(elemConn,1), size(elemConn,1)*(size(elemConn,2)+1)
        write(nunit,'(g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0)')&
            (2, elemConn(i,1)-1,elemConn(i,2)-1, i = 1, size(elemConn,1))

        write(nunit,'(a," ",g0)') 'CELL_TYPES', size(elemConn,1)
        write(nunit,'(g0)') (3 , i = 1, size(elemConn,1))
        close(nunit)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine modify_Xc(this,X,num,dir)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: X
        integer, intent(in) :: num
        integer, intent(in) :: dir

        if (allocated(this%Xc)) then
            this%Xc(num,dir) = X
            call this%set(knot = this%knot, Xc = this%Xc, Wc = this%Wc)
        else
            error stop 'Control points are not set.'
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine modify_Wc(this,W,num)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: W
        integer, intent(in) :: num

        if (allocated(this%Wc)) then
            this%Wc(num) = W
            call this%set(knot = this%knot, Xc = this%Xc, Wc = this%Wc)
        else
            error stop 'The NURBS curve is not rational.'
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_multiplicity(this) result(m)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: m(:)

        ! check
        if (.not.allocated(this%knot)) then
            error stop 'Knot vector is not set.'
        else
            m = compute_multiplicity(this%knot)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_continuity(this) result(c)
        class(nurbs_curve), intent(in) :: this
        integer, allocatable :: c(:)

        ! check
        if (.not.allocated(this%knot)) then
            error stop 'Knot vector is not set.'
        else
            c = this%degree - compute_multiplicity(this%knot)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_nc(this) result(nc)
        class(nurbs_curve), intent(in) :: this
        integer :: nc

        nc = sum(compute_multiplicity(this%knot)) - this%degree - 1
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine insert_knots(this,Xth,r)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: Xth(:)
        integer, intent(in) :: r(:)
        integer :: k, i, s, dim, j, n_new
        real(rk), allocatable :: Xcw(:,:), Xcw_new(:,:), Xc_new(:,:), Wc_new(:), knot_new(:)

        if (allocated(this%Wc)) then ! NURBS

            do i = 1, size(Xth)
                k = findspan(this%nc-1,this%degree,Xth(i),this%knot)
                if (this%knot(k+1) == Xth(i)) then
                    s = compute_multiplicity(this%knot,Xth(i))
                else
                    s = 0
                end if

                dim = size(this%Xc,2)
                allocate(Xcw(size(this%Xc,1),dim+1))
                do j = 1, size(this%Xc,1)
                    Xcw(j,1:dim) = this%Xc(j,1:dim)*this%Wc(j)
                    Xcw(j,dim+1) = this%Wc(j)
                end do

                call insert_knot_A_5_1(&
                    this%degree,&
                    this%knot,&
                    Xcw,&
                    Xth(i),&
                    k,&
                    s,&
                    r(i),&
                    n_new,&
                    knot_new,&
                    Xcw_new)

                allocate(Xc_new(1:n_new+1,1:dim))
                allocate(Wc_new(1:n_new+1))
                do j = 1, n_new+1
                    Xc_new(j,1:dim) = Xcw_new(j-1,1:dim)/Xcw_new(j-1,dim+1)
                    Wc_new(j) = Xcw_new(j-1,dim+1)
                end do

                call this%set(knot=knot_new, Xc=Xc_new, Wc=Wc_new)
                deallocate(Xcw, Xcw_new, Xc_new, Wc_new)
            end do

        else ! B-Spline

            do i = 1, size(Xth)
                k = findspan(this%nc-1,this%degree,Xth(i),this%knot)
                if (this%knot(k+1) == Xth(i)) then
                    s = compute_multiplicity(this%knot,Xth(i))
                else
                    s = 0
                end if

                call insert_knot_A_5_1(&
                    this%degree,&
                    this%knot,&
                    this%Xc,&
                    Xth(i),&
                    k,&
                    s,&
                    r(i),&
                    n_new,&
                    knot_new,&
                    Xc_new)

                deallocate(this%Xc, this%knot)
                call this%set(knot=knot_new, Xc=Xc_new)
            end do

        end if

        call this%create()

    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine elevate_degree(this, t)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in) :: t
        real(rk), allocatable :: Xcw(:,:), Xcw_new(:,:), knot_new(:), Xc_new(:,:), Wc_new(:)
        integer :: dim, j, nc_new

        if (allocated(this%Wc)) then ! NURBS

            dim = size(this%Xc,2)
            allocate(Xcw(size(this%Xc,1),dim+1))
            do j = 1, size(this%Xc,1)
                Xcw(j,1:dim) = this%Xc(j,1:dim)*this%Wc(j)
                Xcw(j,dim+1) = this%Wc(j)
            end do

            call elevate_degree_A(t, this%knot, this%degree, Xcw, nc_new, knot_new, Xcw_new)

            allocate(Xc_new(1:nc_new,1:dim))
            allocate(Wc_new(1:nc_new))
            do j = 1, nc_new
                Xc_new(j,1:dim) = Xcw_new(j,1:dim)/Xcw_new(j,dim+1)
                Wc_new(j) = Xcw_new(j,dim+1)
            end do

            call this%set(knot=knot_new, Xc=Xc_new, Wc=Wc_new)
            deallocate(Xcw, Xcw_new, Xc_new, Wc_new)

        else ! B-Spline

            dim = size(this%Xc,2)

            call elevate_degree_A(t, this%knot, this%degree, this%Xc, nc_new, knot_new, Xc_new)

            call this%set(knot=knot_new, Xc=Xc_new)
            deallocate(Xc_new)

        end if

        call this%create()

    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine derivative(this, res, Xt, dTgc)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in), optional :: res
        real(rk), intent(in), optional :: Xt(:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable :: dTgci(:)
        integer :: i

        ! Set parameter values
        if (present(Xt)) then
            if (allocated(this%Xt)) deallocate(this%Xt)
            this%Xt = Xt
        elseif (present(res)) then
            if (allocated(this%Xt)) deallocate(this%Xt)
            allocate(this%Xt(res))
            this%Xt = [(real(i-1, rk) / real(res-1, rk), i=1, res)]
            ! else
            ! this%Xt = this%Xt
        end if

        allocate(dTgc(size(this%Xt, 1), this%nc))

        if (allocated(this%Wc)) then
            do i = 1, size(this%Xt, 1)
                dTgci = basis_bspline_der(this%Xt(i), this%knot, this%nc, this%degree)
                dTgci = dTgci*(this%Wc/(dot_product(dTgci,this%Wc)))
                dTgc(i,:) = dTgci
            end do
        else
            do i = 1, size(this%Xt, 1)
                dTgci = basis_bspline_der(this%Xt(i), this%knot, this%nc, this%degree)
                dTgc(i,:) = dTgci
            end do
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine basis(this, res, Xt, Tgc)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in), optional :: res
        real(rk), intent(in), optional :: Xt(:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        real(rk), allocatable :: Tgci(:)
        integer :: i

        ! Set parameter values
        if (present(Xt)) then
            if (allocated(this%Xt)) deallocate(this%Xt)
            this%Xt = Xt
        elseif (present(res)) then
            if (allocated(this%Xt)) deallocate(this%Xt)
            allocate(this%Xt(res))
            this%Xt = [(real(i-1, rk) / real(res-1, rk), i=1, res)]
            ! else
            ! this%Xt = this%Xt
        end if

        allocate(Tgc(size(this%Xt, 1), this%nc))

        if (allocated(this%Wc)) then
            do i = 1, size(this%Xt, 1)
                Tgci = basis_bspline(this%Xt(i), this%knot, this%nc, this%degree)
                Tgci = Tgci*(this%Wc/(dot_product(Tgci,this%Wc)))
                Tgc(i,:) = Tgci
            end do
        else
            do i = 1, size(this%Xt, 1)
                Tgci = basis_bspline(this%Xt(i), this%knot, this%nc, this%degree)
                Tgc(i,:) = Tgci
            end do
        end if
    end subroutine
    !===============================================================================

end module forcad_nurbs_curve
