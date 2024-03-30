module forcad_nurbs_volume

    use forcad_utils, only: rk, basis_bspline, elemConn_C0, kron, ndgrid, compute_multiplicity, compute_knot_vector

    implicit none

    private
    public nurbs_volume

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    type nurbs_volume
        real(rk), allocatable, private :: Xc(:,:)  !! control points
        real(rk), allocatable, private :: Xg(:,:)  !! geometry points
        real(rk), allocatable, private :: Wc(:)    !! weights
        real(rk), allocatable, private :: Xt1(:)   !! parameter values in the first direction
        real(rk), allocatable, private :: Xt2(:)   !! parameter values in the second direction
        real(rk), allocatable, private :: Xt3(:)   !! parameter values in the third direction
        real(rk), allocatable, private :: knot1(:) !! knot vector
        real(rk), allocatable, private :: knot2(:) !! knot vector
        real(rk), allocatable, private :: knot3(:) !! knot vector
        integer, private :: order(3)               !! degree of the first direction
        integer, private :: nc(3)                  !! number of control points in each direction
        integer, private :: ng(3)                  !! number of geometry points in each direction
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
        procedure :: finalize            !!> Finalize the NURBS curve object
        procedure :: get_elem_Xc         !!> Generate connectivity for control points
        procedure :: get_elem_Xg         !!> Generate connectivity for geometry points
        procedure :: export_Xc           !!> Export control points to VTK file
        procedure :: export_Xg           !!> Export geometry points to VTK file
        procedure :: modify_Xc           !!> Modify control points
        procedure :: modify_Wc           !!> Modify weights
        procedure :: get_multiplicity    !!> Get multiplicity of the knot vector
        procedure :: get_continuity      !!> Get continuity of the curve
        procedure :: get_nc              !!> Get number of required control points
    end type
    !===============================================================================

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Set control points and weights for the NURBS curve object.
    pure subroutine set1(this, knot1, knot2, knot3, Xc, Wc)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in) :: knot1(:), knot2(:), knot3(:)
        real(rk), intent(in) :: Xc(:,:)
        real(rk), intent(in), optional :: Wc(:)

        this%knot1 = knot1
        this%knot2 = knot2
        this%knot3 = knot3
        this%order = this%get_order()
        this%nc(1) = this%get_nc(1)
        this%nc(2) = this%get_nc(2)
        this%nc(3) = this%get_nc(3)
        this%Xc = Xc
        if (present(Wc)) this%Wc = Wc
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Set control points and weights for the NURBS curve object.
    pure subroutine set2(this, Xth_dir1, Xth_dir2, Xth_dir3, order, continuity1, continuity2, continuity3, Xc, Wc)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in) :: Xth_dir1(:), Xth_dir2(:), Xth_dir3(:)
        integer, intent(in) :: order(:)
        integer, intent(in) :: continuity1(:), continuity2(:), continuity3(:)
        real(rk), intent(in) :: Xc(:,:)
        real(rk), intent(in), optional :: Wc(:)
        integer :: nc(3)

        this%knot1 = compute_knot_vector(Xth_dir1, order(1), continuity1)
        this%knot2 = compute_knot_vector(Xth_dir2, order(2), continuity2)
        this%knot3 = compute_knot_vector(Xth_dir3, order(3), continuity3)
        this%order(1) = order(1)
        this%order(2) = order(2)
        this%order(3) = order(3)
        this%nc(1) = this%get_nc(1)
        this%nc(2) = this%get_nc(2)
        this%nc(3) = this%get_nc(3)
        this%Xc = Xc
        if (present(Wc)) this%Wc = Wc
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine create(this, res1, res2, res3, Xt1, Xt2, Xt3)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in), optional :: res1, res2, res3
        real(rk), intent(in), optional :: Xt1(:), Xt2(:), Xt3(:)
        integer :: i, j
        real(rk), dimension(:), allocatable :: Tgc1, Tgc2, Tgc3, Tgc
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

        ! Set parameter values
        if (present(Xt3)) then
            if (allocated(this%Xt3)) deallocate(this%Xt3)
            this%Xt3 = Xt3
        elseif (present(res3)) then
            if (allocated(this%Xt3)) deallocate(this%Xt3)
            allocate(this%Xt3(res3))
            this%Xt3 = [(real(i-1, rk) / real(res3-1, rk), i=1, res3)]
            ! else
            ! this%Xt3 = this%Xt3
        end if

        ! Set number of geometry points
        this%ng(1) = size(this%Xt1,1)
        this%ng(2) = size(this%Xt2,1)
        this%ng(3) = size(this%Xt3,1)

        call ndgrid(this%Xt1, this%Xt2, this%Xt3, Xt)

        if (allocated(this%Xg)) deallocate(this%Xg)
        allocate(this%Xg(this%ng(1)*this%ng(2)*this%ng(3), size(this%Xc,2)))

        if (allocated(this%Wc)) then ! NURBS volume
            do i = 1, size(Xt, 1)
                Tgc1 = basis_bspline(Xt(i,1), this%knot1, this%nc(1), this%order(1))
                Tgc2 = basis_bspline(Xt(i,2), this%knot2, this%nc(2), this%order(2))
                Tgc3 = basis_bspline(Xt(i,3), this%knot3, this%nc(3), this%order(3))
                Tgc = kron(Tgc3, kron(Tgc2, Tgc1))
                Tgc = Tgc*(this%Wc/(dot_product(Tgc,this%Wc)))
                do j = 1, size(this%Xc, 2)
                    this%Xg(i,j) = dot_product(Tgc,this%Xc(:,j))
                end do
            end do
        else
            do i = 1, size(Xt, 1)
                Tgc1 = basis_bspline(Xt(i,1), this%knot1, this%nc(1), this%order(1))
                Tgc2 = basis_bspline(Xt(i,2), this%knot2, this%nc(2), this%order(2))
                Tgc3 = basis_bspline(Xt(i,3), this%knot3, this%nc(3), this%order(3))
                Tgc = kron(Tgc3, kron(Tgc2, Tgc1))
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
        class(nurbs_volume), intent(in) :: this
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
        class(nurbs_volume), intent(in) :: this
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
        class(nurbs_volume), intent(in) :: this
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
        class(nurbs_volume), intent(in) :: this
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
        elseif (dir == 3) then
            if (allocated(this%Xt3)) then
                Xt = this%Xt3
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
    pure function get_ng(this) result(ng)
        class(nurbs_volume), intent(in) :: this
        integer :: ng(3)

        ng = this%ng
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_order(this) result(order)
        class(nurbs_volume), intent(in) :: this
        integer :: order(3)
        integer, allocatable :: m1(:), m2(:), m3(:)

        m1 = this%get_multiplicity(1)
        m2 = this%get_multiplicity(2)
        m3 = this%get_multiplicity(3)

        order(1) = m1(1) - 1
        order(2) = m2(1) - 1
        order(3) = m3(1) - 1
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_knot(this, dir) result(knot)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: dir
        real(rk), allocatable :: knot(:)

        if (dir == 1) then
            if (allocated(this%knot1)) then
                knot = this%knot1
            else
                error stop 'Knot vector is not set.'
            end if
        elseif (dir == 2) then
            if (allocated(this%knot2)) then
                knot = this%knot2
            else
                error stop 'Knot vector is not set.'
            end if
        elseif (dir == 3) then
            if (allocated(this%knot3)) then
                knot = this%knot3
            else
                error stop 'Knot vector is not set.'
            end if
        else
            error stop 'Invalid direction for knot vector.'
        end if

    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine finalize(this)
        class(nurbs_volume), intent(inout) :: this
        if (allocated(this%Xc)) deallocate(this%Xc)
        if (allocated(this%Xg)) deallocate(this%Xg)
        if (allocated(this%Wc)) deallocate(this%Wc)
        if (allocated(this%Xt1)) deallocate(this%Xt1)
        if (allocated(this%Xt2)) deallocate(this%Xt2)
        if (allocated(this%Xt3)) deallocate(this%Xt3)
        if (allocated(this%knot1)) deallocate(this%knot1)
        if (allocated(this%knot2)) deallocate(this%knot2)
        if (allocated(this%knot3)) deallocate(this%knot3)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine get_elem_Xc(this, elemConn, p)
        class(nurbs_volume), intent(in) :: this
        integer, dimension(:,:), allocatable, intent(out) :: elemConn
        integer, intent(in), optional :: p(:)

        if (present(p)) then
            elemConn = elemConn_C0(this%nc(1), this%nc(2), this%nc(3), p(1), p(2), p(3))
        else
            elemConn = elemConn_C0(this%nc(1), this%nc(2), this%nc(3), 1, 1, 1)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine get_elem_Xg(this, elemConn, p)
        class(nurbs_volume), intent(in) :: this
        integer, dimension(:,:), allocatable, intent(out) :: elemConn
        integer, intent(in), optional :: p(:)

        if (present(p)) then
            elemConn = elemConn_C0(this%ng(1), this%ng(2), this%ng(3), p(1), p(2), p(3))
        else
            elemConn = elemConn_C0(this%ng(1), this%ng(2), this%ng(3), 1, 1, 1)
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_Xc(this, filename)
        class(nurbs_volume), intent(in) :: this
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
        write(nunit,'(a)') 'Xc'
        write(nunit,'(a)') 'ASCII'
        write(nunit,'(a)') 'DATASET UNSTRUCTURED_GRID'
        write(nunit,'(a," ",g0," ",a)') 'POINTS', nc+1, 'double'

        write(nunit,'(f24.18,f24.18,f24.18)') 0.0_rk, 0.0_rk, 0.0_rk
        write(nunit,'(f24.18,f24.18,f24.18)') (this%Xc(i,1), this%Xc(i,2), this%Xc(i,3) , i = 1, nc)

        write(nunit,'(a," ",g0," ",g0)') 'CELLS', size(elemConn,1), size(elemConn,1)*9
        write(nunit,'(g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0)')&
            (8, elemConn(i,1),elemConn(i,2),elemConn(i,4),elemConn(i,3),&
            elemConn(i,5),elemConn(i,6),elemConn(i,8),elemConn(i,7), i = 1, size(elemConn,1))

        write(nunit,'(a," ",g0)') 'CELL_TYPES', size(elemConn,1)
        write(nunit,'(g0)') (12 , i = 1, size(elemConn,1))
        close(nunit)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_Xg(this, filename)
        class(nurbs_volume), intent(in) :: this
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
        write(nunit,'(a)') 'Xg'
        write(nunit,'(a)') 'ASCII'
        write(nunit,'(a)') 'DATASET UNSTRUCTURED_GRID'
        write(nunit,'(a," ",g0," ",a)') 'POINTS', ng+1, 'double'

        write(nunit,'(f24.18,f24.18,f24.18)') 0.0_rk, 0.0_rk, 0.0_rk
        write(nunit,'(f24.18,f24.18,f24.18)') (this%Xg(i,1), this%Xg(i,2), this%Xg(i,3) , i = 1, ng)

        write(nunit,'(a," ",g0," ",g0)') 'CELLS', size(elemConn,1), size(elemConn,1)*9
        write(nunit,'(g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0)')&
            (8, elemConn(i,1),elemConn(i,2),elemConn(i,4),elemConn(i,3),&
            elemConn(i,5),elemConn(i,6),elemConn(i,8),elemConn(i,7), i = 1, size(elemConn,1))

        write(nunit,'(a," ",g0)') 'CELL_TYPES', size(elemConn,1)
        write(nunit,'(g0)') (12 , i = 1, size(elemConn,1))
        close(nunit)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine modify_Xc(this,X,num,dir)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in) :: X
        integer, intent(in) :: num
        integer, intent(in) :: dir

        if (allocated(this%Xc)) then
            this%Xc(num,dir) = X
            call this%set(knot1 = this%knot1, knot2 = this%knot2, knot3 = this%knot3, Xc = this%Xc, Wc = this%Wc)
        else
            error stop 'Control points are not set.'
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine modify_Wc(this,W,num)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in) :: W
        integer, intent(in) :: num

        if (allocated(this%Wc)) then
            this%Wc(num) = W
            call this%set(knot1 = this%knot1, knot2 = this%knot2, knot3 = this%knot3, Xc = this%Xc, Wc = this%Wc)
        else
            error stop 'The NURBS surface is not rational.'
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_multiplicity(this, dir) result(m)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: dir
        integer, allocatable :: m(:)

        if (dir == 1) then

            ! check
            if (.not.allocated(this%knot1)) then
                error stop 'Knot vector is not set.'
            else
                m = compute_multiplicity(this%knot1)
            end if

        elseif (dir == 2) then

            ! check
            if (.not.allocated(this%knot2)) then
                error stop 'Knot vector is not set.'
            else
                m = compute_multiplicity(this%knot2)
            end if

        elseif (dir == 3) then

            ! check
            if (.not.allocated(this%knot3)) then
                error stop 'Knot vector is not set.'
            else
                m = compute_multiplicity(this%knot3)
            end if

        else
            error stop 'Invalid direction.'
        end if

    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_continuity(this, dir) result(c)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: dir
        integer, allocatable :: c(:)

        if (dir == 1) then

            ! check
            if (.not.allocated(this%knot1)) then
                error stop 'Knot vector is not set.'
            else
                c = this%order(1) - compute_multiplicity(this%knot1)
            end if

        elseif (dir == 2) then

            ! check
            if (.not.allocated(this%knot2)) then
                error stop 'Knot vector is not set.'
            else
                c = this%order(2) - compute_multiplicity(this%knot2)
            end if

        elseif (dir == 3) then

            ! check
            if (.not.allocated(this%knot3)) then
                error stop 'Knot vector is not set.'
            else
                c = this%order(3) - compute_multiplicity(this%knot3)
            end if

        else
            error stop 'Invalid direction.'
        end if

    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_nc(this, dir) result(nc)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: dir
        integer :: nc

        if (dir == 1) then

            ! check
            if (.not.allocated(this%knot1)) then
                error stop 'Knot vector is not set.'
            else
                nc = sum(compute_multiplicity(this%knot1)) - this%order(1) - 1
            end if

        elseif (dir == 2) then

            ! check
            if (.not.allocated(this%knot2)) then
                error stop 'Knot vector is not set.'
            else
                nc = sum(compute_multiplicity(this%knot2)) - this%order(2) - 1
            end if

        elseif (dir == 3) then

            ! check
            if (.not.allocated(this%knot3)) then
                error stop 'Knot vector is not set.'
            else
                nc = sum(compute_multiplicity(this%knot3)) - this%order(3) - 1
            end if

        else
            error stop 'Invalid direction.'
        end if

    end function
    !===============================================================================

end module forcad_nurbs_volume
