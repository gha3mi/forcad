module forcad_bezier_surface

    use forcad_utils, only: rk, basis_bernstein, elemConn_C0, kron, ndgrid

    implicit none

    private
    public bezier_surface

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    type bezier_surface
        real(rk), allocatable, private :: Xc(:,:) !! control points
        real(rk), allocatable, private :: Xg(:,:) !! geometry points
        real(rk), allocatable, private :: Wc(:)   !! weights
        real(rk), allocatable, private :: Xt1(:)  !! parameter values in the first direction
        real(rk), allocatable, private :: Xt2(:)  !! parameter values in the second direction
        integer, private :: nc(2)                 !! number of control points in each direction
        integer, private :: ng(2)                 !! number of geometry points in each direction
    contains
        procedure :: set                 !!> Set control points and weights
        procedure :: create              !!> Generate geometry points
        procedure :: get_Xc              !!> Get control points
        procedure :: get_Xg              !!> Get geometry points
        procedure :: get_Wc              !!> Get weights
        procedure :: get_Xt              !!> Get parameter values
        procedure :: get_nc              !!> Get number of control points
        procedure :: get_ng              !!> Get number of geometry points
        procedure :: get_order           !!> Get order of the Bezier surface
        procedure :: finalize            !!> Finalize the Bezier surface object
        procedure :: get_elem_Xc         !!> Generate connectivity for control points
        procedure :: get_elem_Xg         !!> Generate connectivity for geometry points
        procedure :: export_Xc           !!> Export control points to VTK file
        procedure :: export_Xg           !!> Export geometry points to VTK file
        procedure :: modify_Xc           !!> Modify control points
        procedure :: modify_Wc           !!> Modify weights
    end type
    !===============================================================================

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Set control points and weights for the Bezier curve object.
    pure subroutine set(this, nc, Xc, Wc)
        class(bezier_surface), intent(inout) :: this
        integer, intent(in) :: nc(:)
        real(rk), intent(in) :: Xc(:,:)
        real(rk), intent(in), optional :: Wc(:)

        this%Xc = Xc
        this%nc = nc
        if (present(Wc)) this%Wc = Wc
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine create(this, res1, res2, Xt1, Xt2)
        class(bezier_surface), intent(inout) :: this
        integer, intent(in), optional :: res1, res2
        real(rk), intent(in), optional :: Xt1(:), Xt2(:)
        integer :: i, j
        real(rk), dimension(:), allocatable :: Tgc1, Tgc2, Tgc
        real(rk), dimension(:,:), allocatable :: Xt

        ! check
        if (.not.allocated(this%Xc)) then
            error stop 'Control points are not set.'
        end if


        ! Set parameter values
        if (present(Xt1)) then
            if (allocated(this%Xt1)) deallocate(this%Xt1)
            this%Xt1 = Xt1
        elseif (present(res1)) then
            if (allocated(this%Xt1)) deallocate(this%Xt1)
            allocate(this%Xt1(res1))
            this%Xt1 = [(real(i-1, rk) / real(res1-1, rk), i=1, res1)]
            ! else
            ! this%Xt1 = this%Xt1
        end if

        ! Set parameter values
        if (present(Xt2)) then
            if (allocated(this%Xt2)) deallocate(this%Xt2)
            this%Xt2 = Xt2
        elseif (present(res2)) then
            if (allocated(this%Xt2)) deallocate(this%Xt2)
            allocate(this%Xt2(res2))
            this%Xt2 = [(real(i-1, rk) / real(res2-1, rk), i=1, res2)]
            ! else
            ! this%Xt2 = this%Xt2
        end if


        ! Set number of geometry points
        this%ng(1) = size(this%Xt1,1)
        this%ng(2) = size(this%Xt2,1)

        call ndgrid(this%Xt1, this%Xt2, Xt)

        if (allocated(this%Xg)) deallocate(this%Xg)
        allocate(this%Xg(this%ng(1)*this%ng(2), size(this%Xc,2)))

        if (allocated(this%Wc)) then ! Rational Bezier surface
            do i = 1, size(Xt, 1)
                Tgc1 = basis_bernstein(Xt(i,1), this%nc(1))
                Tgc2 = basis_bernstein(Xt(i,2), this%nc(2))
                Tgc = kron(Tgc2, Tgc1)
                Tgc = Tgc*(this%Wc/(dot_product(Tgc,this%Wc)))
                do j = 1, size(this%Xc, 2)
                    this%Xg(i,j) = dot_product(Tgc,this%Xc(:,j))
                end do
            end do
        else ! Non-rational Bezier surface
            do i = 1, size(Xt, 1)
                Tgc1 = basis_bernstein(Xt(i,1), this%nc(1))
                Tgc2 = basis_bernstein(Xt(i,2), this%nc(2))
                Tgc = kron(Tgc2, Tgc1)
                do j = 1, size(this%Xc, 2)
                    this%Xg(i,j) = dot_product(Tgc,this%Xc(:,j))
                end do
            end do
        end if
    end subroutine create
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_Xc(this) result(Xc)
        class(bezier_surface), intent(in) :: this
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
        class(bezier_surface), intent(in) :: this
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
        class(bezier_surface), intent(in) :: this
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
    pure function get_Xt(this, dir) result(Xt)
        class(bezier_surface), intent(in) :: this
        integer, intent(in) :: dir
        real(rk), allocatable :: Xt(:)

        if (dir == 1) then
            if (allocated(this%Xt1)) then
                Xt = this%Xt1
            else
                error stop 'Parameter values are not set.'
            end if
        elseif (dir == 2) then
            if (allocated(this%Xt2)) then
                Xt = this%Xt2
            else
                error stop 'Parameter values are not set.'
            end if
        else
            error stop 'Invalid direction for parameter values.'
        end if

    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_nc(this) result(nc)
        class(bezier_surface), intent(in) :: this
        integer :: nc(2)

        nc = this%nc
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_ng(this) result(ng)
        class(bezier_surface), intent(in) :: this
        integer :: ng(2)

        ng = this%ng
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_order(this) result(order)
        class(bezier_surface), intent(in) :: this
        integer :: order(2)

        order = this%nc - 1
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine finalize(this)
        class(bezier_surface), intent(inout) :: this
        if (allocated(this%Xc)) deallocate(this%Xc)
        if (allocated(this%Xg)) deallocate(this%Xg)
        if (allocated(this%Wc)) deallocate(this%Wc)
        if (allocated(this%Xt1)) deallocate(this%Xt1)
        if (allocated(this%Xt2)) deallocate(this%Xt2)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine get_elem_Xc(this, elemConn, p)
        class(bezier_surface), intent(in) :: this
        integer, dimension(:,:), allocatable, intent(out) :: elemConn
        integer, intent(in), optional :: p(:)

        if (present(p)) then
            elemConn = elemConn_C0(this%nc(1),this%nc(2),p(1),p(2))
        else
            elemConn = elemConn_C0(this%nc(1),this%nc(2),1,1)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine get_elem_Xg(this, elemConn, p)
        class(bezier_surface), intent(in) :: this
        integer, dimension(:,:), allocatable, intent(out) :: elemConn
        integer, intent(in), optional :: p(:)

        if (present(p)) then
            elemConn = elemConn_C0(this%ng(1),this%ng(2),p(1),p(2))
        else
            elemConn = elemConn_C0(this%ng(1),this%ng(2),1,1)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_Xc(this, filename)
        class(bezier_surface), intent(in) :: this
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
        else
            write(nunit,'(g0," ",g0," ",g0)') (this%Xc(i,1), this%Xc(i,2), this%Xc(i,3) , i = 1, nc)
        end if

        write(nunit,'(a," ",g0," ",g0)') 'CELLS', size(elemConn,1), size(elemConn,1)*5
        write(nunit,'(g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0)')&
            (4, elemConn(i,1)-1,elemConn(i,2)-1,elemConn(i,4)-1,elemConn(i,3)-1, i = 1, size(elemConn,1))

        write(nunit,'(a," ",g0)') 'CELL_TYPES', size(elemConn,1)
        write(nunit,'(g0)') (9 , i = 1, size(elemConn,1))
        close(nunit)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_Xg(this, filename)
        class(bezier_surface), intent(in) :: this
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
        else
            write(nunit,'(g0," ",g0," ",g0)') (this%Xg(i,1), this%Xg(i,2), this%Xg(i,3) , i = 1, ng)
        end if

        write(nunit,'(a," ",g0," ",g0)') 'CELLS', size(elemConn,1), size(elemConn,1)*5
        write(nunit,'(g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0)')&
            (4, elemConn(i,1)-1,elemConn(i,2)-1,elemConn(i,4)-1,elemConn(i,3)-1, i = 1, size(elemConn,1))

        write(nunit,'(a," ",g0)') 'CELL_TYPES', size(elemConn,1)
        write(nunit,'(g0)') (9 , i = 1, size(elemConn,1))
        close(nunit)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Modify coordinate of a control point given its index and direction.
    pure subroutine modify_Xc(this,X,num,dir)
        class(bezier_surface), intent(inout) :: this
        real(rk), intent(in) :: X
        integer, intent(in) :: num
        integer, intent(in) :: dir

        if (allocated(this%Xc)) then
            this%Xc(num,dir) = X
            call this%set(nc = this%nc, Xc = this%Xc, Wc = this%Wc)
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
        class(bezier_surface), intent(inout) :: this
        real(rk), intent(in) :: W
        integer, intent(in) :: num

        if (allocated(this%Wc)) then
            this%Wc(num) = W
            call this%set(nc = this%nc, Xc = this%Xc, Wc = this%Wc)
        else
            error stop 'The Bezier curve is not rational.'
        end if
    end subroutine
    !===============================================================================

end module forcad_bezier_surface
