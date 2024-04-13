!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> This module defines the 'nurbs_volume' type for representing a Non-Uniform Rational B-Spline (NURBS) volume.
module forcad_nurbs_volume

    use forcad_utils, only: rk, basis_bspline, elemConn_C0, kron, ndgrid, compute_multiplicity, compute_knot_vector, &
        basis_bspline_der, insert_knot_A_5_1, findspan, elevate_degree_A_5_9, hexahedron_Xc, remove_knots_A_5_8

    implicit none

    private
    public nurbs_volume

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    type nurbs_volume
        real(rk), allocatable, private :: Xc(:,:)  !! Control points (2D array: [nc(1)*nc(2)*nc(3), dim])
        real(rk), allocatable, private :: Xg(:,:)  !! Geometry points (2D array: [ng(1)*ng(2)*ng(3), dim])
        real(rk), allocatable, private :: Wc(:)    !! Weights for the control points (1D array: [nc(1)*nc(2)*nc(3)])
        real(rk), allocatable, private :: Xt1(:)   !! Evaluation parameter values in the first direction (1D array: [ng(1)])
        real(rk), allocatable, private :: Xt2(:)   !! Evaluation parameter values in the second direction (1D array: [ng(2)])
        real(rk), allocatable, private :: Xt3(:)   !! Evaluation parameter values in the third direction (1D array: [ng(3)])
        real(rk), allocatable, private :: knot1(:) !! Knot vector in the first direction (1D array)
        real(rk), allocatable, private :: knot2(:) !! Knot vector in the second direction (1D array)
        real(rk), allocatable, private :: knot3(:) !! Knot vector in the third direction (1D array)
        integer, private :: degree(3)              !! Degree (order) of the volume
        integer, private :: nc(3)                  !! Number of control points in each direction
        integer, private :: ng(3)                  !! Number of geometry points in each direction
        integer, allocatable, private :: elemConn_Xc_vis(:,:) !! Connectivity for visualization of control points
        integer, allocatable, private :: elemConn_Xg_vis(:,:) !! Connectivity for visualization of geometry points
    contains
        procedure :: set1                   !!> Set knot vectors, control points and weights for the NURBS volume object
        procedure :: set2                   !!> Set NURBS volume using nodes of parameter space, degree, continuity, control points and weights
        procedure :: set3                   !!> Set Bezier or Rational Bezier volume using control points and weights
        generic :: set => set1, set2, set3  !!> Set NURBS volume
        procedure :: create                 !!> Generate geometry points
        procedure :: get_Xc                 !!> Get control points
        procedure :: get_Xg                 !!> Get geometry points
        procedure :: get_Wc                 !!> Get weights
        procedure :: get_Xt                 !!> Get parameter values
        procedure, private :: get_knot_all  !!> Get all knot vectors
        procedure, private :: get_knoti     !!> Get i-th knot value
        generic :: get_knot => get_knoti, get_knot_all !!> Get knot vector
        procedure :: get_ng                 !!> Get number of geometry points
        procedure, private :: get_degree_all!!> Get degree of the NURBS volume in all directions
        procedure, private :: get_degree_dir!!> Get degree of the NURBS volume in a specific direction
        generic :: get_degree => get_degree_all, get_degree_dir !!> Get degree of the NURBS volume
        procedure :: finalize               !!> Finalize the NURBS volume object
        procedure :: cmp_elem_Xc_vis        !!> Generate connectivity for control points
        procedure :: cmp_elem_Xg_vis        !!> Generate connectivity for geometry points
        procedure :: get_elem_Xc_vis        !!> Get connectivity for control points
        procedure :: get_elem_Xg_vis        !!> Get connectivity for geometry points
        procedure :: set_elem_Xc_vis        !!> Set connectivity for control points
        procedure :: set_elem_Xg_vis        !!> Set connectivity for geometry points
        procedure :: export_Xc              !!> Export control points to VTK file
        procedure :: export_Xg              !!> Export geometry points to VTK file
        procedure :: modify_Xc              !!> Modify control points
        procedure :: modify_Wc              !!> Modify weights
        procedure :: get_multiplicity       !!> Get multiplicity of the knot vector
        procedure :: get_continuity         !!> Get continuity of the volume
        procedure :: get_nc                 !!> Get number of required control points
        procedure :: derivative             !!> Compute the derivative of the NURBS volume
        procedure :: basis                  !!> Compute the basis functions of the NURBS volume
        procedure :: insert_knots           !!> Insert knots into the knot vector
        procedure :: elevate_degree         !!> Elevate the degree of the NURBS volume
        procedure :: is_rational            !!> Check if the NURBS volume is rational
        procedure :: put_to_nurbs           !!> Put a shape to a NURBS volume
        procedure :: remove_knots           !!> Remove knots from the knot vector

        ! Shapes
        procedure :: set_hexahedron         !!> Set a hexahedron
    end type
    !===============================================================================

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    !> Set control points and weights for the NURBS volume object.
    pure subroutine set1(this, knot1, knot2, knot3, Xc, Wc)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in) :: knot1(:), knot2(:), knot3(:)
        real(rk), intent(in) :: Xc(:,:)
        real(rk), intent(in), optional :: Wc(:)

        this%knot1 = knot1
        this%knot2 = knot2
        this%knot3 = knot3
        this%degree = this%get_degree()
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
    !> Set control points and weights for the NURBS volume object.
    pure subroutine set2(this, Xth_dir1, Xth_dir2, Xth_dir3, degree, continuity1, continuity2, continuity3, Xc, Wc)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in) :: Xth_dir1(:), Xth_dir2(:), Xth_dir3(:)
        integer, intent(in) :: degree(:)
        integer, intent(in) :: continuity1(:), continuity2(:), continuity3(:)
        real(rk), intent(in) :: Xc(:,:)
        real(rk), intent(in), optional :: Wc(:)

        this%knot1 = compute_knot_vector(Xth_dir1, degree(1), continuity1)
        this%knot2 = compute_knot_vector(Xth_dir2, degree(2), continuity2)
        this%knot3 = compute_knot_vector(Xth_dir3, degree(3), continuity3)
        this%degree(1) = degree(1)
        this%degree(2) = degree(2)
        this%degree(3) = degree(3)
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
    !> Set Bezier or Rational Bezier volume using control points and weights.
    pure subroutine set3(this, nc, Xc, Wc)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in) :: nc(:)
        real(rk), intent(in) :: Xc(:,:)
        real(rk), intent(in), optional :: Wc(:)

        if (allocated(this%Xc)) deallocate(this%Xc)

        this%Xc = Xc
        this%nc = nc

        allocate(this%knot1(2*this%nc(1)))
        this%knot1(1:this%nc(1)) = 0.0_rk
        this%knot1(this%nc(1)+1:2*this%nc(1)) = 1.0_rk

        allocate(this%knot2(2*this%nc(2)))
        this%knot2(1:this%nc(2)) = 0.0_rk
        this%knot2(this%nc(2)+1:2*this%nc(2)) = 1.0_rk

        allocate(this%knot3(2*this%nc(3)))
        this%knot3(1:this%nc(3)) = 0.0_rk
        this%knot3(this%nc(3)+1:2*this%nc(3)) = 1.0_rk

        this%degree = this%get_degree()
        if (present(Wc)) then
            if (size(Wc) /= this%nc(1)*this%nc(2)*this%nc(3)) then
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
    pure subroutine create(this, res1, res2, res3, Xt1, Xt2, Xt3, Xt)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in), optional :: res1, res2, res3
        real(rk), intent(in), optional :: Xt1(:), Xt2(:), Xt3(:)
        real(rk), dimension(:,:), intent(in), optional :: Xt
        integer :: i, j
        real(rk), dimension(:), allocatable :: Tgc1, Tgc2, Tgc3, Tgc
        real(rk), dimension(:,:), allocatable :: Xt_

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

        if (present(Xt)) then
            Xt_ = Xt
        else

            ! Set number of geometry points
            this%ng(1) = size(this%Xt1,1)
            this%ng(2) = size(this%Xt2,1)
            this%ng(3) = size(this%Xt3,1)

            call ndgrid(this%Xt1, this%Xt2, this%Xt3, Xt_)
        end if

        if (allocated(this%Xg)) deallocate(this%Xg)
        allocate(this%Xg(size(Xt_,1), size(this%Xc,2)))

        if (allocated(this%Wc)) then ! NURBS volume
            do i = 1, size(Xt_, 1)
                Tgc1 = basis_bspline(Xt_(i,1), this%knot1, this%nc(1), this%degree(1))
                Tgc2 = basis_bspline(Xt_(i,2), this%knot2, this%nc(2), this%degree(2))
                Tgc3 = basis_bspline(Xt_(i,3), this%knot3, this%nc(3), this%degree(3))
                Tgc = kron(Tgc3, kron(Tgc2, Tgc1))
                Tgc = Tgc*(this%Wc/(dot_product(Tgc,this%Wc)))
                do j = 1, size(this%Xc, 2)
                    this%Xg(i,j) = dot_product(Tgc,this%Xc(:,j))
                end do
            end do
        else
            do i = 1, size(Xt_, 1)
                Tgc1 = basis_bspline(Xt_(i,1), this%knot1, this%nc(1), this%degree(1))
                Tgc2 = basis_bspline(Xt_(i,2), this%knot2, this%nc(2), this%degree(2))
                Tgc3 = basis_bspline(Xt_(i,3), this%knot3, this%nc(3), this%degree(3))
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
            error stop 'The NURBS volume is not rational.'
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
    pure function get_degree_all(this) result(degree)
        class(nurbs_volume), intent(in) :: this
        integer :: degree(3)
        integer, allocatable :: m1(:), m2(:), m3(:)

        m1 = this%get_multiplicity(1)
        m2 = this%get_multiplicity(2)
        m3 = this%get_multiplicity(3)

        degree(1) = m1(1) - 1
        degree(2) = m2(1) - 1
        degree(3) = m3(1) - 1
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_degree_dir(this, dir) result(degree)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: dir
        integer :: degree
        integer, allocatable :: m1(:), m2(:), m3(:)

        if (dir == 1) then
            m1 = this%get_multiplicity(1)
            degree = m1(1) - 1
        else if (dir == 2) then
            m2 = this%get_multiplicity(2)
            degree = m2(1) - 1
        else if (dir == 3) then
            m3 = this%get_multiplicity(3)
            degree = m3(1) - 1
        else
            error stop 'Invalid direction for degree.'
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_knot_all(this, dir) result(knot)
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
    pure function get_knoti(this, dir, i) result(knot)
        class(nurbs_volume), intent(in) :: this
        integer, intent(in) :: dir
        integer, intent(in) :: i
        real(rk) :: knot

        if (dir == 1) then
            if (allocated(this%knot1)) then
                if (i < 1 .or. i > size(this%knot1)) then
                    error stop 'Invalid index for knot vector.'
                else
                    knot = this%knot1(i)
                end if
            else
                error stop 'Knot vector is not set.'
            end if
        elseif (dir == 2) then
            if (allocated(this%knot2)) then
                if (i < 1 .or. i > size(this%knot2)) then
                    error stop 'Invalid index for knot vector.'
                else
                    knot = this%knot2(i)
                end if
            else
                error stop 'Knot vector is not set.'
            end if
        elseif (dir == 3) then
            if (allocated(this%knot3)) then
                if (i < 1 .or. i > size(this%knot3)) then
                    error stop 'Invalid index for knot vector.'
                else
                    knot = this%knot3(i)
                end if
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
        if (allocated(this%elemConn_Xc_vis)) deallocate(this%elemConn_Xc_vis)
        if (allocated(this%elemConn_Xg_vis)) deallocate(this%elemConn_Xg_vis)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function cmp_elem_Xc_vis(this, p) result(elemConn)
        class(nurbs_volume), intent(in) :: this
        integer, dimension(:,:), allocatable :: elemConn
        integer, intent(in), optional :: p(:)

        if (present(p)) then
            elemConn = elemConn_C0(this%nc(1), this%nc(2), this%nc(3), p(1), p(2), p(3))
        else
            elemConn = elemConn_C0(this%nc(1), this%nc(2), this%nc(3), 1, 1, 1)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function cmp_elem_Xg_vis(this, p) result(elemConn)
        class(nurbs_volume), intent(in) :: this
        integer, dimension(:,:), allocatable :: elemConn
        integer, intent(in), optional :: p(:)

        if (present(p)) then
            elemConn = elemConn_C0(this%ng(1), this%ng(2), this%ng(3), p(1), p(2), p(3))
        else
            elemConn = elemConn_C0(this%ng(1), this%ng(2), this%ng(3), 1, 1, 1)
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    impure subroutine export_Xc(this, filename)
        class(nurbs_volume), intent(in) :: this
        character(len=*), intent(in) :: filename
        integer :: i, nc, nunit
        integer, dimension(:,:), allocatable :: elemConn

        ! check
        if (.not.allocated(this%Xc)) then
            error stop 'Control points are not set.'
        end if

        if (.not.allocated(this%elemConn_Xc_vis)) then
            elemConn = this%cmp_elem_Xc_vis()
        else
            elemConn = this%elemConn_Xc_vis
        end if

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
            error stop 'Invalid dimension for control points.'
        end if

        write(nunit,'(a," ",g0," ",g0)') 'CELLS', size(elemConn,1), size(elemConn,1)*(size(elemConn,2)+1)
        write(nunit,'(g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0)')&
            (8, elemConn(i,1)-1,elemConn(i,2)-1,elemConn(i,4)-1,elemConn(i,3)-1,&
            elemConn(i,5)-1,elemConn(i,6)-1,elemConn(i,8)-1,elemConn(i,7)-1, i = 1, size(elemConn,1))

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
        integer :: i, ng, nunit
        integer, dimension(:,:), allocatable :: elemConn

        ! check
        if (.not.allocated(this%Xg)) then
            error stop 'Geometry points are not set.'
        end if

        if (.not.allocated(this%elemConn_Xg_vis)) then
            elemConn = this%cmp_elem_Xg_vis()
        else
            elemConn = this%elemConn_Xg_vis
        end if

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
            error stop 'Invalid dimension for geometry points.'
        end if

        write(nunit,'(a," ",g0," ",g0)') 'CELLS', size(elemConn,1), size(elemConn,1)*(size(elemConn,2)+1)
        write(nunit,'(g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0," ",g0)')&
            (8, elemConn(i,1)-1,elemConn(i,2)-1,elemConn(i,4)-1,elemConn(i,3)-1,&
            elemConn(i,5)-1,elemConn(i,6)-1,elemConn(i,8)-1,elemConn(i,7)-1, i = 1, size(elemConn,1))

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
                c = this%degree(1) - compute_multiplicity(this%knot1)
            end if

        elseif (dir == 2) then

            ! check
            if (.not.allocated(this%knot2)) then
                error stop 'Knot vector is not set.'
            else
                c = this%degree(2) - compute_multiplicity(this%knot2)
            end if

        elseif (dir == 3) then

            ! check
            if (.not.allocated(this%knot3)) then
                error stop 'Knot vector is not set.'
            else
                c = this%degree(3) - compute_multiplicity(this%knot3)
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
                nc = sum(compute_multiplicity(this%knot1)) - this%degree(1) - 1
            end if

        elseif (dir == 2) then

            ! check
            if (.not.allocated(this%knot2)) then
                error stop 'Knot vector is not set.'
            else
                nc = sum(compute_multiplicity(this%knot2)) - this%degree(2) - 1
            end if

        elseif (dir == 3) then

            ! check
            if (.not.allocated(this%knot3)) then
                error stop 'Knot vector is not set.'
            else
                nc = sum(compute_multiplicity(this%knot3)) - this%degree(3) - 1
            end if

        else
            error stop 'Invalid direction.'
        end if

    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine derivative(this, res1, res2, res3, Xt1, Xt2, Xt3, dTgc)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in), optional :: res1, res2, res3
        real(rk), intent(in), optional :: Xt1(:), Xt2(:), Xt3(:)
        real(rk), allocatable, intent(out) :: dTgc(:,:)
        real(rk), allocatable :: dTgci(:)
        integer :: i
        real(rk), dimension(:), allocatable :: dTgc1, dTgc2, dTgc3
        real(rk), dimension(:,:), allocatable :: Xt

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

        allocate(dTgc(this%ng(1)*this%ng(2)*this%ng(3), this%nc(1)*this%nc(2)*this%nc(3)))

        if (allocated(this%Wc)) then ! NURBS volume
            do i = 1, size(Xt, 1)
                dTgc1 = basis_bspline_der(Xt(i,1), this%knot1, this%nc(1), this%degree(1))
                dTgc2 = basis_bspline_der(Xt(i,2), this%knot2, this%nc(2), this%degree(2))
                dTgc3 = basis_bspline_der(Xt(i,3), this%knot3, this%nc(3), this%degree(3))
                dTgci = kron(dTgc3, kron(dTgc2, dTgc1))
                dTgci = dTgci*(this%Wc/(dot_product(dTgci,this%Wc)))
                dTgc(i,:) = dTgci
            end do
        else
            do i = 1, size(Xt, 1)
                dTgc1 = basis_bspline_der(Xt(i,1), this%knot1, this%nc(1), this%degree(1))
                dTgc2 = basis_bspline_der(Xt(i,2), this%knot2, this%nc(2), this%degree(2))
                dTgc3 = basis_bspline_der(Xt(i,3), this%knot3, this%nc(3), this%degree(3))
                dTgci = kron(dTgc3, kron(dTgc2, dTgc1))
                dTgc(i,:) = dTgci
            end do
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine basis(this, res1, res2, res3, Xt1, Xt2, Xt3, Tgc)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in), optional :: res1, res2, res3
        real(rk), intent(in), optional :: Xt1(:), Xt2(:), Xt3(:)
        real(rk), allocatable, intent(out) :: Tgc(:,:)
        real(rk), allocatable :: Tgci(:)
        integer :: i
        real(rk), dimension(:), allocatable :: Tgc1, Tgc2, Tgc3
        real(rk), dimension(:,:), allocatable :: Xt

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

        allocate(Tgc(this%ng(1)*this%ng(2)*this%ng(3), this%nc(1)*this%nc(2)*this%nc(3)))

        if (allocated(this%Wc)) then ! NURBS volume
            do i = 1, size(Xt, 1)
                Tgc1 = basis_bspline(Xt(i,1), this%knot1, this%nc(1), this%degree(1))
                Tgc2 = basis_bspline(Xt(i,2), this%knot2, this%nc(2), this%degree(2))
                Tgc3 = basis_bspline(Xt(i,3), this%knot3, this%nc(3), this%degree(3))
                Tgci = kron(Tgc3, kron(Tgc2, Tgc1))
                Tgci = Tgci*(this%Wc/(dot_product(Tgci,this%Wc)))
                Tgc(i,:) = Tgci
            end do
        else
            do i = 1, size(Xt, 1)
                Tgc1 = basis_bspline(Xt(i,1), this%knot1, this%nc(1), this%degree(1))
                Tgc2 = basis_bspline(Xt(i,2), this%knot2, this%nc(2), this%degree(2))
                Tgc3 = basis_bspline(Xt(i,3), this%knot3, this%nc(3), this%degree(3))
                Tgci = kron(Tgc3, kron(Tgc2, Tgc1))
                Tgc(i,:) = Tgci
            end do
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine insert_knots(this, dir ,Xth,r)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in) :: dir
        real(rk), intent(in) :: Xth(:)
        integer, intent(in) :: r(:)
        integer :: k, i, s, dim, j, n_new
        real(rk), allocatable :: Xc(:,:), Xcw(:,:), Xcw_new(:,:), Xc_new(:,:), Wc_new(:), knot_new(:)
        real(rk), allocatable :: Xc4(:,:,:,:)


        if (dir == 1) then ! direction 1

            if (allocated(this%Wc)) then ! NURBS

                do i = 1, size(Xth)
                    k = findspan(this%nc(1)-1,this%degree(1),Xth(i),this%knot1)
                    if (this%knot1(k+1) == Xth(i)) then
                        s = compute_multiplicity(this%knot1, Xth(i))
                    else
                        s = 0
                    end if

                    dim = size(this%Xc,2)
                    allocate(Xcw(size(this%Xc,1),dim+1))
                    do j = 1, size(this%Xc,1)
                        Xcw(j,1:dim) = this%Xc(j,1:dim)*this%Wc(j)
                    end do
                    Xcw(:,dim+1) = this%Wc(:)

                    Xcw = reshape(Xcw,[this%nc(1),this%nc(2)*this%nc(3)*(dim+1)])

                    call insert_knot_A_5_1(&
                        this%degree(1),&
                        this%knot1,&
                        Xcw,&
                        Xth(i),&
                        k,&
                        s,&
                        r(i),&
                        n_new,&
                        knot_new,&
                        Xcw_new)

                    Xcw_new = reshape(Xcw_new,[(n_new+1)*this%nc(2)*this%nc(3),dim+1])

                    allocate(Xc_new(1:(n_new+1)*this%nc(2)*this%nc(3),1:dim))
                    allocate(Wc_new(1:(n_new+1)*this%nc(2)*this%nc(3)))
                    do j = 1, (n_new+1)*this%nc(2)*this%nc(3)
                        Xc_new(j,1:dim) = Xcw_new(j,1:dim)/Xcw_new(j,dim+1)
                    end do
                    Wc_new(:) = Xcw_new(:,dim+1)

                    deallocate(this%Xc, this%knot1, this%Wc)
                    call this%set(knot1=knot_new, knot2=this%knot2, knot3=this%knot3, Xc=Xc_new, Wc=Wc_new)
                    deallocate(Xcw, Xcw_new, Xc_new, Wc_new)
                end do


            else ! B-Spline

                do i = 1, size(Xth)
                    k = findspan(this%nc(1)-1,this%degree(1),Xth(i),this%knot1)
                    if (this%knot1(k+1) == Xth(i)) then
                        s = compute_multiplicity(this%knot1, Xth(i))
                    else
                        s = 0
                    end if

                    dim = size(this%Xc,2)

                    Xc = reshape(this%Xc,[this%nc(1),this%nc(2)*this%nc(3)*dim])

                    call insert_knot_A_5_1(&
                        this%degree(1),&
                        this%knot1,&
                        Xc,&
                        Xth(i),&
                        k,&
                        s,&
                        r(i),&
                        n_new,&
                        knot_new,&
                        Xc_new)

                    Xc_new = reshape(Xc_new,[(n_new+1)*this%nc(2)*this%nc(3),dim])

                    deallocate(this%Xc, this%knot1)
                    call this%set(knot1=knot_new, knot2=this%knot2, knot3=this%knot3, Xc=Xc_new)
                end do

            end if

        elseif (dir == 2) then! direction 2

            if (allocated(this%Wc)) then ! NURBS

                do i = 1, size(Xth)
                    k = findspan(this%nc(2)-1,this%degree(2),Xth(i),this%knot2)
                    if (this%knot2(k+1) == Xth(i)) then
                        s = compute_multiplicity(this%knot2, Xth(i))
                    else
                        s = 0
                    end if

                    dim = size(this%Xc,2)
                    allocate(Xcw(size(this%Xc,1),dim+1))
                    do j = 1, size(this%Xc,1)
                        Xcw(j,1:dim) = this%Xc(j,1:dim)*this%Wc(j)
                    end do
                    Xcw(:,dim+1) = this%Wc(:)

                    Xc4 = reshape(Xcw, [this%nc(1),this%nc(2),this%nc(3),dim+1])
                    Xc4 = reshape(Xc4, [this%nc(2),this%nc(1),this%nc(3),dim+1], order=[2,1,3,4])
                    Xcw = reshape(Xc4,[this%nc(2),this%nc(1)*this%nc(3)*(dim+1)])

                    call insert_knot_A_5_1(&
                        this%degree(2),&
                        this%knot2,&
                        Xcw,&
                        Xth(i),&
                        k,&
                        s,&
                        r(i),&
                        n_new,&
                        knot_new,&
                        Xcw_new)

                    Xc4 = reshape(Xcw_new, [n_new+1,this%nc(1),this%nc(3),dim+1])
                    Xc4 = reshape(Xc4, [this%nc(1),n_new+1,this%nc(3),dim+1], order=[2,1,3,4])
                    Xcw_new = reshape(Xc4,[this%nc(1)*(n_new+1)*this%nc(3),dim+1])

                    allocate(Xc_new(1:this%nc(1)*(n_new+1)*this%nc(3),1:dim))
                    allocate(Wc_new(1:this%nc(1)*(n_new+1)*this%nc(3)))
                    do j = 1, this%nc(1)*(n_new+1)*this%nc(3)
                        Xc_new(j,1:dim) = Xcw_new(j,1:dim)/Xcw_new(j,dim+1)
                    end do
                    Wc_new(:) = Xcw_new(:,dim+1)

                    deallocate(this%Xc, this%knot2, this%Wc)
                    call this%set(knot1=this%knot1, knot2=knot_new, knot3=this%knot3, Xc=Xc_new, Wc=Wc_new)
                    deallocate(Xcw, Xcw_new, Xc_new, Wc_new)
                end do

            else ! B-Spline

                do i = 1, size(Xth)
                    k = findspan(this%nc(2)-1,this%degree(2),Xth(i),this%knot2)
                    if (this%knot2(k+1) == Xth(i)) then
                        s = compute_multiplicity(this%knot2, Xth(i))
                    else
                        s = 0
                    end if

                    dim = size(this%Xc,2)

                    Xc4 = reshape(this%Xc, [this%nc(1),this%nc(2),this%nc(3),dim])
                    Xc4 = reshape(Xc4, [this%nc(2),this%nc(1),this%nc(3),dim], order=[2,1,3,4])
                    Xc = reshape(Xc4,[this%nc(2),this%nc(1)*this%nc(3)*dim])

                    call insert_knot_A_5_1(&
                        this%degree(2),&
                        this%knot2,&
                        Xc,&
                        Xth(i),&
                        k,&
                        s,&
                        r(i),&
                        n_new,&
                        knot_new,&
                        Xc_new)

                    Xc4 = reshape(Xc_new, [n_new+1,this%nc(1),this%nc(3),dim])
                    Xc4 = reshape(Xc4, [this%nc(1),n_new+1,this%nc(3),dim], order=[2,1,3,4])
                    Xc_new = reshape(Xc4,[this%nc(1)*(n_new+1)*this%nc(3),dim])

                    deallocate(this%Xc, this%knot2)
                    call this%set(knot1=this%knot1, knot2=knot_new, knot3= this%knot3, Xc=Xc_new)
                end do


            end if

        elseif (dir == 3) then! direction 3

            if (allocated(this%Wc)) then ! NURBS

                do i = 1, size(Xth)
                    k = findspan(this%nc(3)-1,this%degree(3),Xth(i),this%knot3)
                    if (this%knot3(k+1) == Xth(i)) then
                        s = compute_multiplicity(this%knot3, Xth(i))
                    else
                        s = 0
                    end if

                    dim = size(this%Xc,2)
                    allocate(Xcw(size(this%Xc,1),dim+1))
                    do j = 1, size(this%Xc,1)
                        Xcw(j,1:dim) = this%Xc(j,1:dim)*this%Wc(j)
                    end do
                    Xcw(:,dim+1) = this%Wc(:)

                    Xc4 = reshape(Xcw, [this%nc(1),this%nc(2),this%nc(3),dim+1])
                    Xc4 = reshape(Xc4, [this%nc(3),this%nc(2),this%nc(1),dim+1], order=[3,2,1,4])
                    Xcw = reshape(Xc4,[this%nc(3),this%nc(2)*this%nc(1)*(dim+1)])

                    call insert_knot_A_5_1(&
                        this%degree(3),&
                        this%knot3,&
                        Xcw,&
                        Xth(i),&
                        k,&
                        s,&
                        r(i),&
                        n_new,&
                        knot_new,&
                        Xcw_new)

                    Xc4 = reshape(Xcw_new, [n_new+1,this%nc(2),this%nc(1),dim+1])
                    Xc4 = reshape(Xc4, [this%nc(1),this%nc(2),n_new+1,dim+1], order=[3,2,1,4])
                    Xcw_new = reshape(Xc4,[this%nc(1)*this%nc(2)*(n_new+1),dim+1])

                    allocate(Xc_new(1:this%nc(1)*this%nc(2)*(n_new+1),1:dim))
                    allocate(Wc_new(1:this%nc(1)*this%nc(2)*(n_new+1)))
                    do j = 1, this%nc(1)*this%nc(2)*(n_new+1)
                        Xc_new(j,1:dim) = Xcw_new(j,1:dim)/Xcw_new(j,dim+1)
                    end do
                    Wc_new(:) = Xcw_new(:,dim+1)

                    deallocate(this%Xc, this%knot3, this%Wc)
                    call this%set(knot1=this%knot1, knot2=this%knot2, knot3=knot_new, Xc=Xc_new, Wc=Wc_new)
                    deallocate(Xcw, Xcw_new, Xc_new, Wc_new)
                end do

            else ! B-Spline

                do i = 1, size(Xth)
                    k = findspan(this%nc(3)-1,this%degree(3),Xth(i),this%knot3)
                    if (this%knot3(k+1) == Xth(i)) then
                        s = compute_multiplicity(this%knot3, Xth(i))
                    else
                        s = 0
                    end if

                    dim = size(this%Xc,2)

                    Xc4 = reshape(this%Xc, [this%nc(1),this%nc(2),this%nc(3),dim])
                    Xc4 = reshape(Xc4, [this%nc(3),this%nc(2),this%nc(1),dim], order=[3,2,1,4])
                    Xc = reshape(Xc4,[this%nc(3),this%nc(2)*this%nc(1)*dim])

                    call insert_knot_A_5_1(&
                        this%degree(3),&
                        this%knot3,&
                        Xc,&
                        Xth(i),&
                        k,&
                        s,&
                        r(i),&
                        n_new,&
                        knot_new,&
                        Xc_new)

                    Xc4 = reshape(Xc_new, [n_new+1,this%nc(2),this%nc(1),dim])
                    Xc4 = reshape(Xc4, [this%nc(1),this%nc(2),n_new+1,dim], order=[3,2,1,4])
                    Xc_new = reshape(Xc4, [this%nc(1)*this%nc(2)*(n_new+1),dim])

                    deallocate(this%Xc, this%knot3)
                    call this%set(knot1=this%knot1, knot2=this%knot2, knot3=knot_new, Xc=Xc_new)
                end do

            end if

        else
            error stop 'Invalid direction.'
        end if

    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine elevate_degree(this, dir, t)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in) :: dir
        integer, intent(in) :: t
        real(rk), allocatable :: Xc(:,:), Xcw(:,:), Xcw_new(:,:), Xc_new(:,:), Wc_new(:), knot_new(:)
        integer :: nc_new, dim, j
        real(rk), allocatable:: Xc4(:,:,:,:)


        if (dir == 1) then ! direction 1

            if (allocated(this%Wc)) then ! NURBS

                dim = size(this%Xc,2)
                allocate(Xcw(size(this%Xc,1),dim+1))
                do j = 1, size(this%Xc,1)
                    Xcw(j,1:dim) = this%Xc(j,1:dim)*this%Wc(j)
                end do
                Xcw(:,dim+1) = this%Wc(:)

                Xcw = reshape(Xcw,[this%nc(1),this%nc(2)*this%nc(3)*(dim+1)],order=[1,2])

                call elevate_degree_A_5_9(t, this%knot1, this%degree(1), Xcw, nc_new, knot_new, Xcw_new)

                Xcw_new = reshape(Xcw_new,[nc_new*this%nc(2)*this%nc(3),dim+1],order=[1,2])

                allocate(Xc_new(1:nc_new*this%nc(2)*this%nc(3),1:dim))
                allocate(Wc_new(1:nc_new*this%nc(2)*this%nc(3)))
                do j = 1, nc_new*this%nc(2)*this%nc(3)
                    Xc_new(j,1:dim) = Xcw_new(j,1:dim)/Xcw_new(j,dim+1)
                end do
                Wc_new(:) = Xcw_new(:,dim+1)

                deallocate(this%Xc, this%knot1, this%Wc)
                call this%set(knot1=knot_new, knot2=this%knot2, knot3=this%knot3, Xc=Xc_new, Wc=Wc_new)
                deallocate(Xcw, Xcw_new, Xc_new, Wc_new)

            else ! B-Spline

                dim = size(this%Xc,2)

                Xc = reshape(this%Xc,[this%nc(1),this%nc(2)*this%nc(3)*dim],order=[1,2])

                call elevate_degree_A_5_9(t, this%knot1, this%degree(1), Xc, nc_new, knot_new, Xc_new)

                Xc_new = reshape(Xc_new,[nc_new*this%nc(2)*this%nc(3),dim],order=[1,2])

                deallocate(this%Xc, this%knot1)
                call this%set(knot1=knot_new, knot2=this%knot2, knot3=this%knot3, Xc=Xc_new)

            end if

        elseif (dir == 2) then ! direction 2

            if (allocated(this%Wc)) then ! NURBS


                dim = size(this%Xc,2)
                allocate(Xcw(size(this%Xc,1),dim+1))
                do j = 1, size(this%Xc,1)
                    Xcw(j,1:dim) = this%Xc(j,1:dim)*this%Wc(j)
                end do

                Xcw(:,dim+1) = this%Wc(:)

                Xc4 = reshape(Xcw, [this%nc(1),this%nc(2),this%nc(3),dim+1])
                Xc4 = reshape(Xc4, [this%nc(2),this%nc(1),this%nc(3),dim+1], order=[2,1,3,4])
                Xcw = reshape(Xc4,[this%nc(2),this%nc(1)*this%nc(3)*(dim+1)])


                call elevate_degree_A_5_9(t, this%knot2, this%degree(2), Xcw, nc_new, knot_new, Xcw_new)

                Xc4 = reshape(Xcw_new, [nc_new,this%nc(1),this%nc(3),dim+1])
                Xc4 = reshape(Xc4, [this%nc(1),nc_new,this%nc(3),dim+1], order=[2,1,3,4])
                Xcw_new = reshape(Xc4,[this%nc(1)*nc_new*this%nc(3),dim+1])

                allocate(Xc_new(1:this%nc(1)*nc_new*this%nc(3),1:dim))
                allocate(Wc_new(1:this%nc(1)*nc_new*this%nc(3)))
                do j = 1, this%nc(1)*nc_new*this%nc(3)
                    Xc_new(j,1:dim) = Xcw_new(j,1:dim)/Xcw_new(j,dim+1)
                end do

                Wc_new(:) = Xcw_new(:,dim+1)

                deallocate(this%Xc, this%knot2, this%Wc)
                call this%set(knot1=this%knot1, knot2=knot_new, knot3=this%knot3, Xc=Xc_new, Wc=Wc_new)
                deallocate(Xcw, Xcw_new, Xc_new, Wc_new)

            else ! B-Spline

                dim = size(this%Xc,2)

                Xc4 = reshape(this%Xc, [this%nc(1),this%nc(2),this%nc(3),dim])
                Xc4 = reshape(Xc4, [this%nc(2),this%nc(1),this%nc(3),dim], order=[2,1,3,4])
                Xc = reshape(Xc4,[this%nc(2),this%nc(1)*this%nc(3)*dim])

                call elevate_degree_A_5_9(t, this%knot2, this%degree(2), Xc, nc_new, knot_new, Xc_new)

                Xc4 = reshape(Xc_new, [nc_new,this%nc(1),this%nc(3),dim])
                Xc4 = reshape(Xc4, [this%nc(1),nc_new,this%nc(3),dim], order=[2,1,3,4])
                Xc_new = reshape(Xc4,[this%nc(1)*nc_new*this%nc(3),dim])

                deallocate(this%Xc, this%knot2)
                call this%set(knot1=this%knot1, knot2=knot_new, knot3= this%knot3, Xc=Xc_new)

            end if

        elseif (dir == 3) then ! direction 3

            if (allocated(this%Wc)) then ! NURBS

                dim = size(this%Xc,2)
                allocate(Xcw(size(this%Xc,1),dim+1))
                do j = 1, size(this%Xc,1)
                    Xcw(j,1:dim) = this%Xc(j,1:dim)*this%Wc(j)
                end do

                Xcw(:,dim+1) = this%Wc(:)

                Xc4 = reshape(Xcw, [this%nc(1),this%nc(2),this%nc(3),dim+1])
                Xc4 = reshape(Xc4, [this%nc(3),this%nc(2),this%nc(1),dim+1], order=[3,2,1,4])
                Xcw = reshape(Xc4,[this%nc(3),this%nc(2)*this%nc(1)*(dim+1)])

                call elevate_degree_A_5_9(t, this%knot3, this%degree(3), Xcw, nc_new, knot_new, Xcw_new)

                Xc4 = reshape(Xcw_new, [nc_new,this%nc(2),this%nc(1),dim+1])
                Xc4 = reshape(Xc4, [this%nc(1),this%nc(2),nc_new,dim+1], order=[3,2,1,4])
                Xcw_new = reshape(Xc4,[this%nc(1)*this%nc(2)*nc_new,dim+1])

                allocate(Xc_new(1:this%nc(1)*this%nc(2)*nc_new,1:dim))
                allocate(Wc_new(1:this%nc(1)*this%nc(2)*nc_new))
                do j = 1, this%nc(1)*this%nc(2)*nc_new
                    Xc_new(j,1:dim) = Xcw_new(j,1:dim)/Xcw_new(j,dim+1)
                end do

                Wc_new(:) = Xcw_new(:,dim+1)

                deallocate(this%Xc, this%knot3, this%Wc)
                call this%set(knot1=this%knot1, knot2=this%knot2, knot3=knot_new, Xc=Xc_new, Wc=Wc_new)
                deallocate(Xcw, Xcw_new, Xc_new, Wc_new)

            else ! B-Spline

                dim = size(this%Xc,2)

                Xc4 = reshape(this%Xc, [this%nc(1),this%nc(2),this%nc(3),dim])
                Xc4 = reshape(Xc4, [this%nc(3),this%nc(2),this%nc(1),dim], order=[3,2,1,4])
                Xc = reshape(Xc4,[this%nc(3),this%nc(2)*this%nc(1)*dim])

                call elevate_degree_A_5_9(t, this%knot3, this%degree(3), Xc, nc_new, knot_new, Xc_new)

                Xc4 = reshape(Xc_new, [nc_new,this%nc(2),this%nc(1),dim])
                Xc4 = reshape(Xc4, [this%nc(1),this%nc(2),nc_new,dim], order=[3,2,1,4])
                Xc_new = reshape(Xc4,[this%nc(1)*this%nc(2)*nc_new,dim])

                deallocate(this%Xc, this%knot3)
                call this%set(knot1=this%knot1, knot2=this%knot2, knot3=knot_new, Xc=Xc_new)

            end if

        else
            error stop 'Invalid direction.'
        end if

    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function is_rational(this) result(r)
        class(nurbs_volume), intent(in) :: this
        logical :: r

        r = .false.
        if (allocated(this%Wc)) then
            if (any(this%Wc /= this%Wc(1))) then
                r = .true.
            end if
        end if
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine set_elem_Xc_vis(this, elemConn)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in) :: elemConn(:,:)

        if (allocated(this%elemConn_Xc_vis)) deallocate(this%elemConn_Xc_vis)
        this%elemConn_Xc_vis = elemConn
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine set_elem_Xg_vis(this, elemConn)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in) :: elemConn(:,:)

        if (allocated(this%elemConn_Xg_vis)) deallocate(this%elemConn_Xg_vis)
        this%elemConn_Xg_vis = elemConn
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_elem_Xc_vis(this) result(elemConn)
        class(nurbs_volume), intent(in) :: this
        integer, dimension(:,:), allocatable :: elemConn

        elemConn = this%elemConn_Xc_vis
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function get_elem_Xg_vis(this) result(elemConn)
        class(nurbs_volume), intent(in) :: this
        integer, dimension(:,:), allocatable :: elemConn

        elemConn = this%elemConn_Xg_vis
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine set_hexahedron(this, L, nc, Wc)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in) :: L(:)
        integer, intent(in) :: nc(:)
        real(rk), intent(in), optional :: Wc(:)

        if (present(Wc)) then
            call this%set(nc, hexahedron_Xc(L, nc), Wc)
        else
            call this%set(nc, hexahedron_Xc(L, nc))
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine put_to_nurbs(this, X, elemConn)
        class(nurbs_volume), intent(inout) :: this
        real(rk), intent(in) :: X(:,:)
        integer, intent(in) :: elemConn(:,:)
        integer :: i, j
        real(rk), dimension(:), allocatable :: Tgc1, Tgc2, Tgc3, Tgc
        real(rk), dimension(:,:), allocatable :: Xt
        real(rk) :: min_X1, max_X1, min_X2, max_X2, min_X3, max_X3

        ! Assuming knot vectors are in the range [0,1]
        ! Normalize the X coordinates to the range [0,1]
        allocate(Xt(size(X,1), size(X,2)))
        min_X1 = minval(X(:,1))
        max_X1 = maxval(X(:,1))

        min_X2 = minval(X(:,2))
        max_X2 = maxval(X(:,2))

        min_X3 = minval(X(:,3))
        max_X3 = maxval(X(:,3))

        Xt(:,1) = (X(:,1) -  min_X1) / (max_X1 -  min_X1)
        Xt(:,2) = (X(:,2) -  min_X2) / (max_X2 -  min_X2)
        Xt(:,3) = (X(:,3) -  min_X3) / (max_X3 -  min_X3)

        allocate(this%Xg(size(Xt,1), size(this%Xc,2)))

        if (allocated(this%Wc)) then ! NURBS volume
            do i = 1, size(Xt, 1)
                Tgc1 = basis_bspline(Xt(i,1), this%knot1, this%nc(1), this%degree(1))
                Tgc2 = basis_bspline(Xt(i,2), this%knot2, this%nc(2), this%degree(2))
                Tgc3 = basis_bspline(Xt(i,3), this%knot3, this%nc(3), this%degree(3))
                Tgc = kron(Tgc3, kron(Tgc2, Tgc1))
                Tgc = Tgc*(this%Wc/(dot_product(Tgc,this%Wc)))
                do j = 1, size(this%Xc, 2)
                    this%Xg(i,j) = dot_product(Tgc,this%Xc(:,j))
                end do
            end do
        else ! B-Spline volume
            do i = 1, size(Xt, 1)
                Tgc1 = basis_bspline(Xt(i,1), this%knot1, this%nc(1), this%degree(1))
                Tgc2 = basis_bspline(Xt(i,2), this%knot2, this%nc(2), this%degree(2))
                Tgc3 = basis_bspline(Xt(i,3), this%knot3, this%nc(3), this%degree(3))
                Tgc = kron(Tgc3, kron(Tgc2, Tgc1))
                do j = 1, size(this%Xc, 2)
                    this%Xg(i,j) = dot_product(Tgc,this%Xc(:,j))
                end do
            end do
        end if

        call this%set_elem_Xg_vis(elemConn)

    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine remove_knots(this, dir ,Xth,r)
        class(nurbs_volume), intent(inout) :: this
        integer, intent(in) :: dir
        real(rk), intent(in) :: Xth(:)
        integer, intent(in) :: r(:)
        integer :: k, i, s, dim, j, nc_new, t
        real(rk), allocatable :: Xc(:,:), Xcw(:,:), Xcw_new(:,:), Xc_new(:,:), Wc_new(:), knot_new(:)
        real(rk), allocatable :: Xc4(:,:,:,:)


        if (dir == 1) then ! direction 1

            if (allocated(this%Wc)) then ! NURBS

                do i = 1, size(Xth)
                    k = findspan(this%nc(1)-1,this%degree(1),Xth(i),this%knot1)
                    if (this%knot1(k+1) == Xth(i)) then
                        s = compute_multiplicity(this%knot1, Xth(i))
                    else
                        s = 0
                    end if
                    k = k + 1

                    dim = size(this%Xc,2)
                    allocate(Xcw(size(this%Xc,1),dim+1))
                    do j = 1, size(this%Xc,1)
                        Xcw(j,1:dim) = this%Xc(j,1:dim)*this%Wc(j)
                    end do
                    Xcw(:,dim+1) = this%Wc(:)

                    Xcw = reshape(Xcw,[this%nc(1),this%nc(2)*this%nc(3)*(dim+1)],order=[1,2])

                    call remove_knots_A_5_8(&
                        this%degree(1),&
                        this%knot1,&
                        Xcw,&
                        Xth(i),&
                        k,&
                        s,&
                        r(i),&
                        t,&
                        knot_new,&
                        Xcw_new)

                    if (allocated(Xcw)) deallocate(Xcw)

                    if (t == 0) then
                        ! no change
                    else
                        nc_new = size(Xcw_new,1)
                        Xcw_new = reshape(Xcw_new,[(nc_new)*this%nc(2)*this%nc(3),dim+1],order=[1,2])

                        allocate(Xc_new(1:(nc_new)*this%nc(2)*this%nc(3),1:dim))
                        allocate(Wc_new(1:(nc_new)*this%nc(2)*this%nc(3)))
                        do j = 1, (nc_new)*this%nc(2)*this%nc(3)
                            Xc_new(j,1:dim) = Xcw_new(j,1:dim)/Xcw_new(j,dim+1)
                        end do
                        Wc_new(:) = Xcw_new(:,dim+1)

                        deallocate(this%Xc, this%knot1, this%Wc)
                        call this%set(knot1=knot_new, knot2=this%knot2, knot3=this%knot3, Xc=Xc_new, Wc=Wc_new)
                        deallocate(Xcw_new, Xc_new, Wc_new)
                    end if
                end do


            else ! B-Spline

                do i = 1, size(Xth)
                    k = findspan(this%nc(1)-1,this%degree(1),Xth(i),this%knot1)
                    if (this%knot1(k+1) == Xth(i)) then
                        s = compute_multiplicity(this%knot1, Xth(i))
                    else
                        s = 0
                    end if
                    k = k + 1

                    dim = size(this%Xc,2)

                    Xc = reshape(this%Xc,[this%nc(1),this%nc(2)*this%nc(3)*dim])

                    call remove_knots_A_5_8(&
                        this%degree(1),&
                        this%knot1,&
                        Xc,&
                        Xth(i),&
                        k,&
                        s,&
                        r(i),&
                        t,&
                        knot_new,&
                        Xc_new)

                    if (allocated(Xc)) deallocate(Xc)

                    if (t == 0) then
                        ! no change
                    else
                        nc_new = size(Xcw_new,1)
                        Xc_new = reshape(Xc_new,[(nc_new)*this%nc(2)*this%nc(3),dim])

                        deallocate(this%Xc, this%knot1)
                        call this%set(knot1=knot_new, knot2=this%knot2, knot3=this%knot3, Xc=Xc_new)
                    end if
                end do

            end if

        elseif (dir == 2) then! direction 2

            if (allocated(this%Wc)) then ! NURBS

                do i = 1, size(Xth)
                    k = findspan(this%nc(2)-1,this%degree(2),Xth(i),this%knot2)
                    if (this%knot2(k+1) == Xth(i)) then
                        s = compute_multiplicity(this%knot2, Xth(i))
                    else
                        s = 0
                    end if
                    k = k + 1

                    dim = size(this%Xc,2)
                    allocate(Xcw(size(this%Xc,1),dim+1))
                    do j = 1, size(this%Xc,1)
                        Xcw(j,1:dim) = this%Xc(j,1:dim)*this%Wc(j)
                    end do

                    Xcw(:,dim+1) = this%Wc(:)

                    Xc4 = reshape(Xcw, [this%nc(1),this%nc(2),this%nc(3),dim+1])
                    Xc4 = reshape(Xc4, [this%nc(2),this%nc(1),this%nc(3),dim+1], order=[2,1,3,4])
                    Xcw = reshape(Xc4,[this%nc(2),this%nc(1)*this%nc(3)*(dim+1)])

                    call remove_knots_A_5_8(&
                        this%degree(2),&
                        this%knot2,&
                        Xcw,&
                        Xth(i),&
                        k,&
                        s,&
                        r(i),&
                        t,&
                        knot_new,&
                        Xcw_new)

                    if (allocated(Xcw)) deallocate(Xcw)

                    if (t == 0) then
                        ! no change
                    else
                        nc_new = size(Xcw_new,1)

                        Xc4 = reshape(Xcw_new, [nc_new,this%nc(1),this%nc(3),dim+1])
                        Xc4 = reshape(Xc4, [this%nc(1),nc_new,this%nc(3),dim+1], order=[2,1,3,4])
                        Xcw_new = reshape(Xc4,[this%nc(1)*(nc_new)*this%nc(3),dim+1])

                        allocate(Xc_new(1:this%nc(1)*(nc_new)*this%nc(3),1:dim))
                        allocate(Wc_new(1:this%nc(1)*(nc_new)*this%nc(3)))
                        do j = 1, this%nc(1)*(nc_new)*this%nc(3)
                            Xc_new(j,1:dim) = Xcw_new(j,1:dim)/Xcw_new(j,dim+1)
                        end do

                        Wc_new(:) = Xcw_new(:,dim+1)

                        deallocate(this%Xc, this%knot2, this%Wc)
                        call this%set(knot1=this%knot1, knot2=knot_new, knot3=this%knot3, Xc=Xc_new, Wc=Wc_new)
                        deallocate(Xcw_new, Xc_new, Wc_new)
                    end if
                end do

            else ! B-Spline

                do i = 1, size(Xth)
                    k = findspan(this%nc(2)-1,this%degree(2),Xth(i),this%knot2)
                    if (this%knot2(k+1) == Xth(i)) then
                        s = compute_multiplicity(this%knot2, Xth(i))
                    else
                        s = 0
                    end if
                    k = k + 1

                    dim = size(this%Xc,2)

                    Xc4 = reshape(this%Xc, [this%nc(1),this%nc(2),this%nc(3),dim])
                    Xc4 = reshape(Xc4, [this%nc(2),this%nc(1),this%nc(3),dim], order=[2,1,3,4])
                    Xc = reshape(Xc4, [this%nc(2),this%nc(1)*this%nc(3)*dim])

                    call remove_knots_A_5_8(&
                        this%degree(2),&
                        this%knot2,&
                        Xc,&
                        Xth(i),&
                        k,&
                        s,&
                        r(i),&
                        t,&
                        knot_new,&
                        Xc_new)

                    if (allocated(Xc)) deallocate(Xc)

                    if (t == 0) then
                        ! no change
                    else
                        nc_new = size(Xcw_new,1)

                        Xc4 = reshape(Xc_new, [nc_new,this%nc(1),this%nc(3),dim])
                        Xc4 = reshape(Xc4, [this%nc(1),nc_new,this%nc(3),dim], order=[2,1,3,4])
                        Xc_new = reshape(Xc4, [this%nc(1)*(nc_new)*this%nc(3),dim])

                        deallocate(this%Xc, this%knot2)
                        call this%set(knot1=this%knot1, knot2=knot_new, knot3= this%knot3, Xc=Xc_new)
                    end if
                end do

            end if

        elseif (dir == 3) then! direction 3

            if (allocated(this%Wc)) then ! NURBS

                do i = 1, size(Xth)
                    k = findspan(this%nc(3)-1,this%degree(3),Xth(i),this%knot3)
                    if (this%knot3(k+1) == Xth(i)) then
                        s = compute_multiplicity(this%knot3, Xth(i))
                    else
                        s = 0
                    end if
                    k = k + 1

                    dim = size(this%Xc,2)
                    allocate(Xcw(size(this%Xc,1),dim+1))
                    do j = 1, size(this%Xc,1)
                        Xcw(j,1:dim) = this%Xc(j,1:dim)*this%Wc(j)
                    end do
                    Xcw(:,dim+1) = this%Wc(:)

                    Xc4 = reshape(Xcw, [this%nc(1),this%nc(2),this%nc(3),dim+1])
                    Xc4 = reshape(Xc4, [this%nc(3),this%nc(2),this%nc(1),dim+1], order=[3,2,1,4])
                    Xcw = reshape(Xc4, [this%nc(3),this%nc(2)*this%nc(1)*(dim+1)])

                    call remove_knots_A_5_8(&
                        this%degree(3),&
                        this%knot3,&
                        Xcw,&
                        Xth(i),&
                        k,&
                        s,&
                        r(i),&
                        t,&
                        knot_new,&
                        Xcw_new)

                    if (allocated(Xcw)) deallocate(Xcw)

                    if (t == 0) then
                        ! no change
                    else
                        nc_new = size(Xcw_new,1)

                        Xc4 = reshape(Xcw_new, [nc_new,this%nc(2),this%nc(1),dim+1])
                        Xc4 = reshape(Xc4, [this%nc(1),this%nc(2),nc_new,dim+1], order=[3,2,1,4])
                        Xcw_new = reshape(Xc4, [this%nc(1)*this%nc(2)*(nc_new),dim+1])

                        allocate(Xc_new(1:this%nc(1)*this%nc(2)*(nc_new),1:dim))
                        allocate(Wc_new(1:this%nc(1)*this%nc(2)*(nc_new)))
                        do j = 1, this%nc(1)*this%nc(2)*(nc_new)
                            Xc_new(j,1:dim) = Xcw_new(j,1:dim)/Xcw_new(j,dim+1)
                        end do
                        Wc_new(:) = Xcw_new(:,dim+1)

                        deallocate(this%Xc, this%knot3, this%Wc)
                        call this%set(knot1=this%knot1, knot2=this%knot2, knot3=knot_new, Xc=Xc_new, Wc=Wc_new)
                        deallocate(Xcw_new, Xc_new, Wc_new)
                    end if
                end do

            else ! B-Spline

                do i = 1, size(Xth)
                    k = findspan(this%nc(3)-1,this%degree(3),Xth(i),this%knot3)
                    if (this%knot3(k+1) == Xth(i)) then
                        s = compute_multiplicity(this%knot3, Xth(i))
                    else
                        s = 0
                    end if
                    k = k + 1

                    dim = size(this%Xc,2)

                    Xc4 = reshape(this%Xc, [this%nc(1),this%nc(2),this%nc(3),dim])
                    Xc4 = reshape(Xc4, [this%nc(3),this%nc(2),this%nc(1),dim], order=[3,2,1,4])
                    Xc = reshape(Xc4, [this%nc(3),this%nc(2)*this%nc(1)*dim])

                    call remove_knots_A_5_8(&
                        this%degree(3),&
                        this%knot3,&
                        Xc,&
                        Xth(i),&
                        k,&
                        s,&
                        r(i),&
                        t,&
                        knot_new,&
                        Xc_new)

                    if (allocated(Xc)) deallocate(Xc)

                    if (t == 0) then
                        ! no change
                    else
                        nc_new = size(Xcw_new,1)

                        Xc4 = reshape(Xc_new, [nc_new,this%nc(2),this%nc(1),dim])
                        Xc4 = reshape(Xc4, [this%nc(1),this%nc(2),nc_new,dim], order=[3,2,1,4])
                        Xc_new = reshape(Xc4, [this%nc(1)*this%nc(2)*(nc_new),dim])

                        deallocate(this%Xc, this%knot3)
                        call this%set(knot1=this%knot1, knot2=this%knot2, knot3=knot_new, Xc=Xc_new)
                    end if
                end do


            end if

        else
            error stop 'Invalid direction.'
        end if

    end subroutine
    !===============================================================================

end module forcad_nurbs_volume
