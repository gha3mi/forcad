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
        procedure :: elevate_degree      !!> Elevate the degree of the curve
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
        real(rk), intent(in) :: Xth(:)
        real(rk), allocatable :: Xcw(:,:), Xcw_new(:,:), knot_new(:), Xc_new(:,:), Wc_new(:)
        real(rk) :: alpha
        integer :: dim, nc, nXth, nknot, a, b, r, l, i, j, m, n, s, q, ind

        if (allocated(this%Wc)) then
            allocate(Xcw(size(this%Xc,1),size(this%Xc,2)+1))
            do i = 1, size(this%Xc,2)
                Xcw(:,i) = (this%Xc(:,i)*this%Wc)
            end do
            Xcw(:,size(this%Xc,2)+1) = this%Wc
        else
            allocate(Xcw(size(this%Xc,1),size(this%Xc,2)))
            do i = 1, size(this%Xc,2)
                Xcw(:,i) = (this%Xc(:,i))
            end do
        end if
        nc = size(Xcw, 1)
        dim = size(Xcw, 2)
        nXth = size(Xth)
        nknot = size(this%knot)
        allocate(Xcw_new(nc+nXth, dim), source=0.0_rk)
        allocate(knot_new(nknot+nXth), source=0.0_rk)
        n = nc - 1
        r = nXth - 1
        m = n + this%order + 1
        a = findspan(n, this%order, Xth(1), this%knot)
        b = findspan(n, this%order, Xth(r+1), this%knot)
        b = b + 1
        do q = 0,dim-1
            do j = 0,a-this%order
                Xcw_new(j+1,q+1) = Xcw(j+1,q+1)
            end do
            do j = b-1,n
                Xcw_new(j+r+2,q+1) = Xcw(j+1,q+1)
            end do
        end do
        do j = 0,a
            knot_new(j+1) = this%knot(j+1)
        end do
        do j = b+this%order,m
            knot_new(j+r+2) = this%knot(j+1)
        end do
        i = b + this%order - 1
        s = b + this%order + r
        do j = r,0,-1
            do while (Xth(j+1) <= this%knot(i+1) .and. i > a)
                do q = 0,dim-1
                    Xcw_new(s-this%order,q+1) = Xcw(i-this%order,q+1)
                end do
                knot_new(s+1) = this%knot(i+1)
                s = s - 1
                i = i - 1
            end do
            do q = 0,dim-1
                Xcw_new(s-this%order,q+1) = Xcw_new(s-this%order+1,q+1)
            end do
            do l = 1,this%order
                ind = s - this%order + l
                alpha = knot_new(s+l+1) - Xth(j+1)
                if (abs(alpha) == 0) then
                    do q = 0,dim-1
                        Xcw_new(ind,q+1) = Xcw_new(ind+1,q+1)
                    end do
                else
                    alpha = alpha/(knot_new(s+l+1) - this%knot(i-this%order+l+1))
                    do q = 0,dim-1
                        Xcw_new(ind,q+1) = (1.0_rk-alpha)*Xcw_new(ind+1,q+1) + alpha*Xcw_new(ind,q+1)
                    end do
                end if
            end do
            knot_new(s+1) = Xth(j+1)
            s = s - 1
        end do
        if (allocated(this%Wc)) then
            Xc_new = Xcw_new(:,1:size(this%Xc,2))
            Wc_new = Xcw_new(:,size(this%Xc,2)+1)
            do concurrent (i = 1: size(this%Xc, 2))
                Xc_new(:,i) = Xc_new(:,i) / Wc_new(:)
            end do
            deallocate(this%Xc, this%knot, this%Xg, this%Wc)
            call this%set(knot = knot_new, Xc = Xc_new, Wc = Wc_new)
            call this%create()
        else
            Xc_new = Xcw_new(1:size(this%Xc,2),:)
            deallocate(this%Xc, this%knot, this%Xg)
            call this%set(knot = knot_new, Xc = Xc_new)
            call this%create()
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function findspan(n,order,Xth,knot) result(s)
        integer, intent(in) :: n, order
        real(rk), intent(in) :: Xth
        real(rk), intent(in) :: knot(:)
        integer :: s
        integer :: low, high, mid
        if (Xth == knot(n+2)) then
            s = n
            return
        end if
        low = order
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
    pure subroutine elevate_degree(this, t)
        class(nurbs_curve), intent(inout) :: this
        integer, intent(in) :: t
        real(rk), allocatable :: Xcw(:,:), Xcw_new(:,:), knot_new(:), Xc_new(:,:), Wc_new(:)
        real(rk), allocatable :: bezalfs(:,:), bpts(:,:), ebpts(:,:), Nextbpts(:,:), alfs(:)
        real(rk) :: inv, alpha1, alpha2, Xth1, Xth2, numer, den
        integer :: n, lbz, rbz, sv, tr, kj, first, knoti, last, alpha3, ii, dim, nc, nc_new
        integer :: i, j, q, s, m, ph, ph2, mpi, mh, r, a, b, Xcwi, oldr, mul
        integer, allocatable :: mlp(:)

        if (allocated(this%Wc)) then
            allocate(Xcw(size(this%Xc,1),size(this%Xc,2)+1))
            do i = 1, size(this%Xc,2)
                Xcw(:,i) = (this%Xc(:,i)*this%Wc)
            end do
            Xcw(:,size(this%Xc,2)+1) = this%Wc
        else
            allocate(Xcw(size(this%Xc,1),size(this%Xc,2)))
            do i = 1, size(this%Xc,2)
                Xcw(:,i) = (this%Xc(:,i))
            end do
        end if
        nc = size(Xcw,1)
        dim = size(Xcw,2)
        mlp = compute_multiplicity(this%knot)
        mlp = mlp + t
        nc_new = sum(mlp) - (mlp(1)-1) - 1
        allocate(Xcw_new(nc_new,dim), source=0.0_rk)
        allocate(bezalfs(this%order+1,this%order+t+1), source=0.0_rk)
        allocate(bpts(this%order+1,dim), source=0.0_rk)
        allocate(ebpts(this%order+t+1,dim), source=0.0_rk)
        allocate(Nextbpts(this%order+1,dim), source=0.0_rk)
        allocate(alfs(this%order), source=0.0_rk)
        n = nc - 1
        m = n + this%order + 1
        ph = this%order + t
        ph2 = ph / 2
        bezalfs(1,1) = 1.0_rk
        bezalfs(this%order+1,ph+1) = 1.0_rk
        do i = 1,ph2
            inv = 1.0_rk/bincoeff(ph,i)
            mpi = min(this%order,i)
            do j = max(0,i-t),mpi
                bezalfs(j+1,i+1) = inv*bincoeff(this%order,j)*bincoeff(t,i-j)
            end do
        end do
        do i = ph2+1,ph-1
            mpi = min(this%order,i)
            do j = max(0,i-t),mpi
                bezalfs(j+1,i+1) = bezalfs(this%order-j+1,ph-i+1)
            end do
        end do
        mh = ph
        knoti = ph+1
        r = -1
        a = this%order
        b = this%order+1
        Xcwi = 1
        Xth1 = this%knot(1)
        do ii =0,dim-1
            Xcw_new(1,ii+1) = Xcw(1,ii+1)
        end do
        allocate(knot_new(sum(mlp)), source=0.0_rk)
        do i = 0,ph
            knot_new(i+1) = Xth1
        end do
        do i = 0,this%order
            do ii = 0,dim-1
                bpts(i+1,ii+1) = Xcw(i+1,ii+1)
            end do
        end do
        do while (b<m)
            i = b
            do while (b<m .and. this%knot(b+1) == this%knot(b+2))
                b = b + 1
                if (b+2 > size(this%knot)) then
                    exit
                end if
            end do
            mul = b - i + 1
            mh = mh + mul + t
            Xth2 = this%knot(b+1)
            oldr = r
            r = this%order - mul
            if (oldr > 0) then
                lbz = (oldr+2)/2
            else
                lbz = 1
            end if
            if (r > 0) then
                rbz = ph - (r+1)/2
            else
                rbz = ph
            end if
            if (r>0) then
                numer = Xth2 - Xth1
                do q = this%order,mul+1,-1
                    alfs(q-mul) = numer / (this%knot(a+q+1)-Xth1)
                end do
                do j = 1,r
                    sv = r - j
                    s = mul + j
                    do q = this%order,s,-1
                        do ii = 0,dim-1
                            bpts(q+1,ii+1) = (1.0_rk-alfs(q-s+1))*bpts(q,ii+1) + alfs(q-s+1)*bpts(q+1,ii+1)
                        end do
                    end do
                    do ii = 0,dim-1
                        Nextbpts(sv+1,ii+1) = bpts(this%order+1,ii+1)
                    end do
                end do
            end if
            do i = lbz,ph
                do ii = 0,dim-1
                    ebpts(i+1,ii+1) = 0.0_rk
                end do
                mpi = min(this%order,i)
                do j = max(0,i-t),mpi
                    do ii = 0,dim-1
                        ebpts(i+1,ii+1) = bezalfs(j+1,i+1)*bpts(j+1,ii+1) + ebpts(i+1,ii+1)
                    end do
                end do
            end do
            if (oldr > 1) then
                first = knoti - 2
                last = knoti
                den = Xth2 - Xth1
                alpha3 = floor((Xth2-this%knot(knoti)) / den)
                do tr = 1,oldr-1
                    i = first
                    j = last
                    kj = j - knoti + 1
                    do while (j-i > tr)
                        if (i < Xcwi) then
                            alpha1 = (Xth2-this%knot(i+1))/(Xth1-this%knot(i+1))
                            do ii = 0,dim-1
                                Xcw_new(i+1,ii+1) = (1-alpha1)*Xcw_new(i,ii+1) + alpha1*Xcw_new(i+1,ii+1)
                            end do
                        end if
                        if (j >= lbz) then
                            if (j-tr <= knoti-ph+oldr) then
                                alpha2 = (Xth2-knot_new(j-tr+1)) / den
                                do ii = 0,dim-1
                                    ebpts(kj+1,ii+1) = alpha2*ebpts(kj+1,ii+1) + (1-alpha2)*ebpts(kj+2,ii+1)
                                end do
                            else
                                do ii = 0,dim-1
                                    ebpts(kj+1,ii+1) = (1-alpha3)*ebpts(ii+1,kj+2) + alpha3*ebpts(kj+1,ii+1)
                                end do
                            end if
                        end if
                        i = i + 1
                        j = j - 1
                        kj = kj - 1
                    end do
                    first = first - 1
                    last = last + 1
                end do
            end if
            if (a /= this%order) then
                do i = 0,ph-oldr-1
                    knot_new(knoti+1) = Xth1
                    knoti = knoti + 1
                end do
            end if
            do j = lbz,rbz
                do ii = 0,dim-1
                    Xcw_new(Xcwi+1,ii+1) = ebpts(j+1,ii+1)
                end do
                Xcwi = Xcwi + 1
            end do
            if (b<m) then
                do j = 0,r-1
                    do ii = 0,dim-1
                        bpts(j+1,ii+1) = Nextbpts(j+1,ii+1)
                    end do
                end do
                do j = r,this%order
                    do ii = 0,dim-1
                        bpts(j+1,ii+1) = Xcw(b-this%order+j+1,ii+1)
                    end do
                end do
                a = b
                b = b+1
                Xth1 = Xth2
            else
                do i = 0,ph
                    knot_new(knoti+i+1) = Xth2
                end do
            end if
        end do
        if (allocated(this%Wc)) then
            Xc_new = Xcw_new(:,1:size(this%Xc,2))
            Wc_new = Xcw_new(:,size(this%Xc,2)+1)
            do concurrent (i = 1: size(this%Xc, 2))
                Xc_new(:,i) = Xc_new(:,i) / Wc_new(:)
            end do
            deallocate(this%Xc, this%knot, this%Xg, this%Wc)
            call this%set(knot = knot_new, Xc = Xc_new, Wc = Wc_new)
            call this%create()
        else
            Xc_new = Xcw_new(1:size(this%Xc,2),:)
            deallocate(this%Xc, this%knot, this%Xg)
            call this%set(knot = knot_new, Xc = Xc_new)
            call this%create()
        end if
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function bincoeff(n,k) result(b)
        integer, intent(in) :: n, k
        real(rk) :: b
        b = floor(0.5_rk+exp(factln(n)-factln(k)-factln(n-k)))
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function factln(n) result(f)
        integer, intent(in) :: n
        real(rk) :: f
        if (n <= 1) then
            f = 0.0_rk
            return
        end if
        f = log(gamma(real(n+1,rk)))
    end function
    !===============================================================================

end module forcad_nurbs_curve
