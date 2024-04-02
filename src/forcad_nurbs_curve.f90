module forcad_nurbs_curve

    use forcad_utils, only: rk, basis_bspline, elemConn_C0, compute_multiplicity, compute_knot_vector

    implicit none

    private
    public nurbs_curve

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    type nurbs_curve
        real(rk), allocatable, private :: Xc(:,:) !! control points
        real(rk), allocatable, private :: Xg(:,:) !! geometry points
        real(rk), allocatable, private :: Wc(:)   !! weights
        real(rk), allocatable, private :: Xt(:)   !! evaluation points
        real(rk), allocatable, private :: knot(:) !! knot vector
        integer, private :: order                 !! order of the curve
        integer, private :: nc                    !! number of control points
        integer, private :: ng                    !! number of geometry points
    contains
        procedure :: set1                !!> Set control points and weights
        procedure :: set2                !!> Set control points and weights
        generic :: set => set1, set2
        procedure :: create              !!> Generate geometry points
        procedure :: get_Xc              !!> Get control points
        procedure :: get_Xg              !!> Get geometry points
        procedure :: get_Wc              !!> Get weights
        procedure :: get_Xt              !!> Get parameter values
        procedure :: get_knot            !!> Get knot vector
        procedure :: get_ng              !!> Get number of geometry points
        procedure :: get_order           !!> Get order of the Bezier curve
        procedure :: finalize            !!> Finalize the Bezier curve object
        procedure :: get_elem_Xc         !!> Generate connectivity for control points
        procedure :: get_elem_Xg         !!> Generate connectivity for geometry points
        procedure :: export_Xc           !!> Export control points to VTK file
        procedure :: export_Xg           !!> Export geometry points to VTK file
        procedure :: modify_Xc           !!> Modify control points
        procedure :: modify_Wc           !!> Modify weights
        procedure :: get_multiplicity    !!> Get multiplicity of the knot vector
        procedure :: get_continuity      !!> Get continuity of the curve
        procedure :: get_nc              !!> Get number of required control points
        procedure :: insert_knot         !!> Insert a new knot
    end type
    !===============================================================================

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Set control points and weights for the Bezier curve object.
    pure subroutine set1(this, knot, Xc, Wc)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: knot(:)
        real(rk), intent(in) :: Xc(:,:)
        real(rk), intent(in), optional :: Wc(:)

        this%knot = knot
        this%order = this%get_order()
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
    !> Set control points and weights for the Bezier curve object.
    pure subroutine set2(this, Xth_dir, order, continuity, Xc, Wc)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: Xth_dir(:)
        integer, intent(in) :: order
        integer, intent(in) :: continuity(:)
        real(rk), intent(in) :: Xc(:,:)
        real(rk), intent(in), optional :: Wc(:)

        this%knot = compute_knot_vector(Xth_dir, order, continuity)
        this%order = order
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
                Tgc = basis_bspline(this%Xt(i), this%knot, this%nc, this%order)
                Tgc = Tgc*(this%Wc/(dot_product(Tgc,this%Wc)))
                do j = 1, size(this%Xc, 2)
                    this%Xg(i,j) = dot_product(Tgc,this%Xc(:,j))
                end do
            end do
        else
            do i = 1, size(this%Xt, 1)
                Tgc = basis_bspline(this%Xt(i), this%knot, this%nc, this%order)
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
            error stop 'The Bezier curve is not rational.'
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
    pure function get_order(this) result(order)
        class(nurbs_curve), intent(in) :: this
        integer :: order
        integer, allocatable :: m(:)

        m = this%get_multiplicity()

        order = m(1) - 1
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
            error stop 'The Bezier curve is not rational.'
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
            c = this%order - compute_multiplicity(this%knot)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_nc(this) result(nc)
        class(nurbs_curve), intent(in) :: this
        integer :: nc

        nc = sum(compute_multiplicity(this%knot)) - this%order - 1
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine insert_knot(this, Xth)
        class(nurbs_curve), intent(inout) :: this
        real(rk), intent(in) :: Xth
        real(rk) :: alpha
        real(rk), allocatable :: Xc_new(:,:), Wc_new(:), knot_new(:)
        integer, allocatable :: m(:)
        integer :: k, i, nc_new, order_new

        ! Find the index where the new knot (Xth) should be inserted
        ! knot span [k, k+1)
        k = 0
        do while (this%knot(k+1) < Xth .and. k+1 <= size(this%knot))
            k = k + 1
        end do

        ! Check if the new knot is within the range of existing knots and is not the first or the last knot
        ! Insert the new knot and update the knot vector
        if (k+1 > 1 .and. k+1 < size(this%knot)) then
            if (Xth >= this%knot(k) .and. Xth <= this%knot(k+1)) then
                allocate(knot_new(size(this%knot)+1))
                knot_new = [this%knot(1:k), Xth, this%knot(k+1:size(this%knot))]
            else
                error stop 'Invalid new knot value. It should be within the range of existing knots.'
            end if
        else
            error stop 'Invalid new knot value. It should be greater than the first knot and less than the last knot.'
        end if

        ! Compute the new control points
        m = compute_multiplicity(knot_new)
        order_new = m(1) - 1
        nc_new = sum(m) - order_new - 1


        if (allocated(this%Wc)) then ! NURBS

            allocate(Xc_new(nc_new, size(this%Xc, 2)))
            allocate(Wc_new(nc_new))

            do i = 1,nc_new
                if (i < k-this%order+1) then
                    Xc_new(i, :) = this%Xc(i, :)*this%Wc(i)
                    Wc_new(i) = this%Wc(i)
                else if (i > k) then
                    Xc_new(i, :) = this%Xc(i-1, :)*this%Wc(i-1)
                    Wc_new(i) = this%Wc(i-1)
                else
                    alpha = (Xth - this%knot(i)) / (this%knot(i+this%order) - this%knot(i))
                    Xc_new(i, :) = (1.0_rk-alpha)*this%Xc(i-1, :)*this%Wc(i-1) + alpha*this%Xc(i, :)*this%Wc(i)
                    Wc_new(i) = (1.0_rk-alpha)*this%Wc(i-1) + alpha*this%Wc(i)
                end if
            end do

            ! Normalize the new control points
            do concurrent (i = 1: size(Xc_new, 2))
                Xc_new(:, i) = Xc_new(:, i) / Wc_new(:)
            end do

            deallocate(this%Xc, this%knot, this%Xg, this%Wc)
            call this%set(knot = knot_new, Xc = Xc_new, Wc = Wc_new)
            call this%create()

        else ! B-Spline

            allocate(Xc_new(nc_new, size(this%Xc, 2)))

            do i = 1,nc_new
                if (i < k-this%order+1) then
                    Xc_new(i, :) = this%Xc(i, :)
                else if (i > k) then
                    Xc_new(i, :) = this%Xc(i-1, :)
                else
                    alpha = (Xth - this%knot(i)) / (this%knot(i+this%order) - this%knot(i))
                    Xc_new(i, :) = (1.0_rk-alpha)*this%Xc(i-1, :) + alpha*this%Xc(i, :)
                end if
            end do

            deallocate(this%Xc, this%knot, this%Xg)
            call this%set(knot = knot_new, Xc = Xc_new)
            call this%create()

        end if

    end subroutine
    !===============================================================================

end module forcad_nurbs_curve
