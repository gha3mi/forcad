module forcad_bezier_curve

    use forcad_utils, only: rk, basis_bernstein, elemConn_C0

    implicit none

    private
    public bezier_curve

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    type bezier_curve
        real(rk), allocatable, private :: Xc(:,:) !! Control points
        real(rk), allocatable, private :: Xg(:,:) !! Geometry points
        real(rk), allocatable, private :: Wc(:)   !! Weights
        real(rk), allocatable, private :: Xt(:)   !! Parameter values
        integer, private :: nc                    !! Number of control points
        integer, private :: ng                    !! Number of geometry points
    contains
        procedure :: set                 !!> Set control points and weights
        procedure :: create              !!> Generate geometry points
        procedure :: get_Xc              !!> Get control points
        procedure :: get_Xg              !!> Get geometry points
        procedure :: get_Wc              !!> Get weights
        procedure :: get_Xt              !!> Get parameter values
        procedure :: get_nc              !!> Get number of control points
        procedure :: get_ng              !!> Get number of geometry points
        procedure :: get_order           !!> Get order of the Bezier curve
        procedure :: finalize            !!> Finalize the Bezier curve object
        procedure :: get_elem_Xc         !!> Generate connectivity for control points
        procedure :: get_elem_Xg         !!> Generate connectivity for geometry points
        procedure :: export_Xc           !!> Export control points to VTK file
        procedure :: export_Xg           !!> Export geometry points to VTK file
        procedure :: modify_Xc           !!> Modify control points
        procedure :: modify_Wc           !!> Modify weights
        procedure :: elevate_degree      !!> Elevate the degree of the Bezier curve
    end type
    !===============================================================================

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Set control points and weights for the Bezier curve object.
    pure subroutine set(this, Xc, Wc)
        class(bezier_curve), intent(inout) :: this
        real(rk), intent(in) :: Xc(:,:)
        real(rk), intent(in), optional :: Wc(:)

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
    !> Generate geometry points of the Bezier curve.
    pure subroutine create(this, res, Xt)
        class(bezier_curve), intent(inout) :: this
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
        this%ng = size(this%Xt,1)

        ! Allocate memory for geometry points
        if (allocated(this%Xg)) deallocate(this%Xg)
        allocate(this%Xg(this%ng, size(this%Xc,2)))

        ! Compute geometry points
        if (allocated(this%Wc)) then ! Rational Bezier curve
            do i = 1, size(this%Xt, 1)
                Tgc = basis_bernstein(this%Xt(i), this%nc)
                Tgc = Tgc*(this%Wc/(dot_product(Tgc,this%Wc)))
                do j = 1, size(this%Xc, 2)
                    this%Xg(i,j) = dot_product(Tgc,this%Xc(:,j))
                end do
            end do
        else ! Non-rational Bezier curve
            do i = 1, size(this%Xt, 1)
                Tgc = basis_bernstein(this%Xt(i), this%nc)
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
        class(bezier_curve), intent(in) :: this
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
        class(bezier_curve), intent(in) :: this
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
        class(bezier_curve), intent(in) :: this
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
        class(bezier_curve), intent(in) :: this
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
    pure function get_nc(this) result(nc)
        class(bezier_curve), intent(in) :: this
        integer :: nc

        nc = this%nc
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_ng(this) result(ng)
        class(bezier_curve), intent(in) :: this
        integer :: ng

        ng = this%ng
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_order(this) result(order)
        class(bezier_curve), intent(in) :: this
        integer :: order

        order = this%nc - 1
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Finalize the Bezier curve object by deallocating memory.
    pure subroutine finalize(this)
        class(bezier_curve), intent(inout) :: this
        if (allocated(this%Xc)) deallocate(this%Xc)
        if (allocated(this%Xg)) deallocate(this%Xg)
        if (allocated(this%Wc)) deallocate(this%Wc)
        if (allocated(this%Xt)) deallocate(this%Xt)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Generate connectivity for control points.
    pure subroutine get_elem_Xc(this, elemConn, p)
        class(bezier_curve), intent(in) :: this
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
    !> Generate connectivity for geometry points.
    pure subroutine get_elem_Xg(this, elemConn, p)
        class(bezier_curve), intent(in) :: this
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
    !> Export control points to a VTK file.
    impure subroutine export_Xc(this, filename)
        class(bezier_curve), intent(in) :: this
        character(len=*), intent(in) :: filename
        integer :: i, j, nc, nunit
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

        if (size(this%Xc,2) == 2) then
            write(nunit,'(g0," ",g0," ",g0)') (this%Xc(i,1), this%Xc(i,2), 0.0_rk , i = 1, nc)
        elseif (size(this%Xc,2) == 3) then
            write(nunit,'(g0," ",g0," ",g0)') (this%Xc(i,1), this%Xc(i,2), this%Xc(i,3) , i = 1, nc)
        else
            error stop 'Invalid dimension of control points.'
        end if

        write(nunit,'(a," ",g0," ",g0)') 'CELLS', size(elemConn,1), size(elemConn,1)*3
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
    !> Export geometry points to a VTK file.
    impure subroutine export_Xg(this, filename)
        class(bezier_curve), intent(in) :: this
        character(len=*), intent(in) :: filename
        integer :: i, j, ng, nunit
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

        if (size(this%Xg,2) == 2) then
            write(nunit,'(g0," ",g0," ",g0)') (this%Xg(i,1), this%Xg(i,2), 0.0_rk , i = 1, ng)
        elseif (size(this%Xg,2) == 3) then
            write(nunit,'(g0," ",g0," ",g0)') (this%Xg(i,1), this%Xg(i,2), this%Xg(i,3) , i = 1, ng)
        else
            error stop 'Invalid dimension of geometry points.'
        end if

        write(nunit,'(a," ",g0," ",g0)') 'CELLS', size(elemConn,1), size(elemConn,1)*3
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
    !> Modify coordinate of a control point given its index and direction.
    pure subroutine modify_Xc(this,X,num,dir)
        class(bezier_curve), intent(inout) :: this
        real(rk), intent(in) :: X
        integer, intent(in) :: num
        integer, intent(in) :: dir

        if (allocated(this%Xc)) then
            this%Xc(num,dir) = X
            call this%set(Xc = this%Xc, Wc = this%Wc)
        else
            error stop 'Control points are not set.'
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Modify weight of a control point given its index.
    pure subroutine modify_Wc(this,W,num)
        class(bezier_curve), intent(inout) :: this
        real(rk), intent(in) :: W
        integer, intent(in) :: num

        if (allocated(this%Wc)) then
            this%Wc(num) = W
            call this%set(Xc = this%Xc, Wc = this%Wc)
        else
            error stop 'The Bezier curve is not rational.'
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Elevate the degree of the Bezier curve by one.
    pure subroutine elevate_degree(this)
        class(bezier_curve), intent(inout) :: this
        integer :: nc_new, j, i
        real(rk), allocatable :: Xc_new(:,:)
        real(rk), allocatable :: Wc_new(:)

        ! check
        if (.not.allocated(this%Xc)) then
            error stop 'Control points are not set.'
        end if

        if (allocated(this%Wc)) then ! Rational Bezier curve

            ! Calculate the new number of control points
            nc_new = this%nc + 1

            allocate(Xc_new(nc_new, size(this%Xc, 2)))
            allocate(Wc_new(nc_new))

            ! Compute new control points
            Xc_new(1,:) = this%Xc(1,:) * this%Wc(1)
            Wc_new(1) = this%Wc(1)
            Xc_new(nc_new,:) = this%Xc(this%nc,:) * this%Wc(this%nc)
            Wc_new(nc_new) = this%Wc(this%nc)
            do concurrent (j = 2: this%nc)
                do i = 1, size(this%Xc, 2)
                    Xc_new(j, i) = (j-1) / real(nc_new - 1, rk) * this%Xc(j-1, i) * this%Wc(j-1) + &
                        (1 - (j-1) / real(nc_new - 1, rk)) * this%Xc(j, i) * this%Wc(j)
                end do
                Wc_new(j) = (j-1) / real(nc_new - 1, rk) * this%Wc(j-1) + &
                    (1 - (j-1) / real(nc_new - 1, rk)) * this%Wc(j)
            end do

            ! Normalize the new control points
            do concurrent (i = 1: size(this%Xc, 2))
                Xc_new(:, i) = Xc_new(:, i) / Wc_new(:)
            end do

            ! Update geometry points
            deallocate(this%Xc, this%Wc)
            call this%set(Xc = Xc_new, Wc = Wc_new)
            call this%create(Xt = this%Xt)

        else ! Non-rational Bezier curve

            ! Calculate the new number of control points
            nc_new = this%nc + 1

            allocate(Xc_new(nc_new, size(this%Xc, 2)))

            ! Compute new control points
            Xc_new(1,:) = this%Xc(1,:)
            Xc_new(nc_new,:) = this%Xc(this%nc,:)
            do concurrent (j = 2: this%nc)
                do i = 1, size(this%Xc, 2)
                    Xc_new(j, i) = (j-1) / real(nc_new - 1, rk) * this%Xc(j-1, i) + &
                        (1 - (j-1) / real(nc_new - 1, rk)) * this%Xc(j, i)
                end do
            end do

            ! Update geometry points
            call this%set(Xc = Xc_new)
            call this%create(Xt = this%Xt)

        end if
    end subroutine
    !===============================================================================

end module forcad_bezier_curve
