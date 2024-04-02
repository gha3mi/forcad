module forcad_utils

    implicit none

    private
    public :: rk, basis_bernstein, basis_bspline, elemConn_C0, kron, ndgrid, compute_multiplicity, compute_knot_vector

    integer, parameter :: rk = kind(1.0d0)

    !===============================================================================
    interface elemConn_C0
        module procedure cmp_elemConn_C0_L
        module procedure cmp_elemConn_C0_S
        module procedure cmp_elemConn_C0_V
    end interface
    !===============================================================================


    !===============================================================================
    interface ndgrid
        module procedure ndgrid2
        module procedure ndgrid3
    end interface
    !===============================================================================

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function basis_bspline(Xt, knot, nc, order) result(B)
        integer, intent(in)   :: order
        real(rk), intent(in)  :: knot(:)
        integer, intent(in)   :: nc
        real(rk), intent(in)  :: Xt
        real(rk)              :: temp, Xth_i, Xth_i1, Xth_ip, Xth_ip1
        real(rk), allocatable :: Nt(:,:)
        integer               :: i, p
        real(rk), allocatable :: B(:)

        temp = abs(Xt - knot(size(knot)))
        allocate(Nt(nc, 0:order), source=0.0_rk)

        do p = 0, order
            do i = 1, nc
                Xth_i   = knot(i)
                Xth_i1  = knot(i+1)
                Xth_ip  = knot(i+p)
                Xth_ip1 = knot(i+p+1)

                if ( temp /= tiny(0.0_rk) .and. Xt >= Xth_i .and. Xt <= Xth_i1 ) Nt(i,0) = 1.0_rk
                if ( Xth_ip /= Xth_i ) Nt(i,p) = (Xt - Xth_i)/(Xth_ip - Xth_i) * Nt(i,p-1)
                if ( Xth_ip1 /= Xth_i1 ) Nt(i,p) = Nt(i,p) + (Xth_ip1 - Xt)/(Xth_ip1  - Xth_i1) * Nt(i+1,p-1)
            end do
        end do
        B = Nt(:,order)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function basis_bernstein(Xt, nc) result(B)
        real(rk), intent(in) :: Xt
        integer, intent(in) :: nc
        real(rk), allocatable :: B(:)
        integer  :: p, order

        order = nc - 1

        allocate(B(nc), source=0.0_rk)

        do concurrent (p = 0:order)
            B(p+1) = gamma(real(nc, kind=rk))/(gamma(real(p+1, kind=rk))*gamma(real(nc-p, kind=rk)))
            if (Xt == 0.0_rk .and. p == 0) then
                B(p+1) = B(p+1)*(1.0_rk-Xt)**(order-p)
            else if (Xt == 0.0_rk .and. order-p == 0) then
                B(p+1) = B(p+1)*(Xt**p)
            else
                B(p+1) = B(p+1)*(Xt**p)*(1.0_rk-Xt)**(order-p)
            end if
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function kron(u,v) result(w)
        real(rk), dimension(:), intent(in), contiguous  :: u, v
        real(rk), dimension(size(u)*size(v)) :: w
        integer                              :: i, j, m, n

        m = size(u)
        n = size(v)

        do i = 1, m
            do j = 1, n
                w((i-1)*n + j) = u(i)*v(j)
            end do
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine ndgrid2(X_dir1,X_dir2, Xt)
        real(rk), dimension(:), intent(in), contiguous :: X_dir1, X_dir2
        real(rk), dimension(:,:), allocatable, intent(out) :: Xt
        integer :: s1, s2, i, j, n

        s1 = size(X_dir1)
        s2 = size(X_dir2)
        allocate(Xt(s1*s2,2))
        n = 0
        do j = 1, s2
            do i = 1, s1
                n = n + 1
                Xt(n,1) = X_dir1(i)
                Xt(n,2) = X_dir2(j)
            end do
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine ndgrid3(X_dir1,X_dir2,X_dir3, Xt)
        real(rk), dimension(:), intent(in), contiguous :: X_dir1, X_dir2, X_dir3
        real(rk), dimension(:,:), allocatable, intent(out) :: Xt
        integer :: s1, s2, s3, i, j, k, n

        s1 = size(X_dir1)
        s2 = size(X_dir2)
        s3 = size(X_dir3)
        allocate(Xt(s1*s2*s3,3))
        n = 0
        do k = 1, s3
            do j = 1, s2
                do i = 1, s1
                    n = n + 1
                    Xt(n,1) = X_dir1(i)
                    Xt(n,2) = X_dir2(j)
                    Xt(n,3) = X_dir3(k)
                end do
            end do
        end do
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function repelem(a, b) result(c)
        real(rk), dimension(:), intent(in), contiguous  :: a
        integer, dimension(:), intent(in), contiguous  :: b
        real(rk), dimension(sum(b)) :: c
        integer :: i, l, n

        l = 0
        do i = 1, size(a)
            n = b(i)
            c(l+1:l+n) = a(i)
            l = l + n
        end do
    end function repelem
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function cmp_elemConn_C0_L(nnode,p) result(elemConn)
        integer, intent(in) :: nnode
        integer, intent(in) :: p
        integer, dimension(:,:), allocatable :: elemConn
        integer :: i, l
        integer, dimension(:), allocatable :: nodes

        allocate(elemConn( ((nnode-p) / p) ,2))
        nodes = [(i, i=1,nnode)]
        l = 0
        do i = 1,nnode-p, p
            l = l+1
            elemConn(l,:) = reshape(nodes(i:i+p),[(p + 1)])
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function cmp_elemConn_C0_S(nnode1,nnode2,p1,p2) result(elemConn)
        integer, intent(in) :: nnode1, nnode2
        integer, intent(in) :: p1, p2
        integer, dimension(:,:), allocatable :: elemConn
        integer :: i, j, l
        integer, dimension(:,:), allocatable :: nodes

        allocate(elemConn( ((nnode1-p1) / p1) * ((nnode2-p2) / p2) ,4))
        nodes = reshape([(i, i=1,nnode1*nnode2)], [nnode1, nnode2])
        l = 0
        do j = 1,nnode2-p2, p2
            do i = 1,nnode1-p1, p1
                l = l+1
                elemConn(l,:) = reshape(nodes(i:i+p1,j:j+p2),[(p1 + 1)*(p2 + 1)])
            end do
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function cmp_elemConn_C0_V(nnode1,nnode2,nnode3,p1,p2,p3) result(elemConn)
        integer, intent(in) :: nnode1, nnode2, nnode3
        integer, intent(in) :: p1, p2, p3
        integer, dimension(:,:), allocatable :: elemConn
        integer :: i, j, k, l
        integer, dimension(:,:,:), allocatable :: nodes

        allocate(elemConn( ((nnode1-p1) / p1) * ((nnode2-p2) / p2) * ((nnode3-p3) / p3) ,8))
        nodes = reshape([(i, i=1,nnode1*nnode2*nnode3)], [nnode1, nnode2, nnode3])
        l = 0
        do k = 1,nnode3-p3, p3
            do j = 1,nnode2-p2, p2
                do i = 1,nnode1-p1, p1
                    l = l+1
                    elemConn(l,:) = reshape(nodes(i:i+p1,j:j+p2,k:k+p3),[(p1 + 1)*(p2 + 1)*(p3 + 1)])
                end do
            end do
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function compute_multiplicity(knot) result(multiplicity)
        real(rk), intent(in) :: knot(:)
        integer, dimension(:), allocatable :: multiplicity
        integer :: i, count

        count = 1
        do i = 2, size(knot)
            if (knot(i) /= knot(i-1)) count = count + 1
        end do

        allocate(multiplicity(count))

        multiplicity(1) = 1
        count = 1

        do i = 2, size(knot)
            if (knot(i) /= knot(i-1)) then
                count = count + 1
                multiplicity(count) = 1
            else
                multiplicity(count) = multiplicity(count) + 1
            end if
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function compute_knot_vector(Xth_dir, order, continuity) result(knot)
        real(rk), intent(in) :: Xth_dir(:)
        integer, intent(in) :: order
        integer, intent(in) :: continuity(:)
        real(rk), allocatable :: knot(:)

        knot = repelem(Xth_dir, (order - continuity))
    end function
    !===============================================================================

end module forcad_utils
