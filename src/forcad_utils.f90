!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> This module contains parameters, functions and subroutines that are used in the library.
module forcad_utils

    implicit none

    private
    public :: rk, basis_bernstein, basis_bspline, elemConn_C0, kron, ndgrid, compute_multiplicity, compute_knot_vector, &
        basis_bspline_der, insert_knot_A_5_1, findspan, elevate_degree_A_5_9, hexahedron_Xc, remove_knots_A_5_8

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


    !===============================================================================
    interface compute_multiplicity
        module procedure compute_multiplicity1
        module procedure compute_multiplicity2
    end interface
    !===============================================================================

contains

    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function basis_bspline(Xt, knot, nc, degree) result(B)
        integer, intent(in)   :: degree
        real(rk), intent(in)  :: knot(:)
        integer, intent(in)   :: nc
        real(rk), intent(in)  :: Xt
        real(rk)              :: temp, Xth_i, Xth_i1, Xth_ip, Xth_ip1
        real(rk), allocatable :: Nt(:,:)
        integer               :: i, p
        real(rk), allocatable :: B(:)

        temp = abs(Xt - knot(size(knot)))
        allocate(Nt(nc, 0:degree), source=0.0_rk)

        do p = 0, degree
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
        B = Nt(:,degree)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    pure function basis_bspline_der(Xt, knot, nc, degree) result(dB)
        integer, intent(in)   :: degree
        real(rk), intent(in)  :: knot(:)
        integer, intent(in)   :: nc
        real(rk), intent(in)  :: Xt
        real(rk), allocatable :: dB(:)
        real(rk), allocatable :: Nt(:,:), dNt_dXt(:,:)
        real(rk)              :: R, L, Rp, Lp, knot_i, knot_ip, knot_jk, knot_jkm, knot_end, a, b, c, d
        integer               :: i, k, n, m, jk

        k = degree + 1
        n = nc - 1
        allocate(Nt(nc+degree, degree+1))
        Nt = 0.0_rk
        do i = 1, n+k
            knot_i   = knot(i)
            knot_ip  = knot(i+1)
            knot_end = knot(size(knot))
            if ( abs(Xt - knot_end) > tiny(0.0_rk) ) then
                if ( Xt >= knot_i .and. Xt < knot_ip ) Nt(i,1) = 1.0_rk
            elseif ( abs(Xt - knot_end) < tiny(0.0_rk) ) then
                if ( Xt >= knot_i .and. Xt <= knot_ip ) Nt(i,1) = 1.0_rk
            end if
        end do
        allocate(dNt_dXt(nc+degree, degree+1))
        dNt_dXt = 0.0_rk
        m = 0
        do jk = 2, k
            m = m + 1
            do i = 1, n + k - m
                knot_i   = knot(i)
                knot_ip  = knot(i+1)
                knot_jk  = knot(i+jk)
                knot_jkm = knot(i+jk-1)
                a        = (knot_jkm - knot_i)
                b        = (knot_jk - Xt)
                c        = (knot_jk - knot_ip)
                d        = (Xt - knot_i)
                R = d/a
                if ( isnan(R) .or. isinf(R) .or. abs(R) < tiny(0.0_rk) ) R = 0.0_rk
                L = b/c
                if ( isnan(L) .or. isinf(L) .or. abs(L) < tiny(0.0_rk) ) L = 0.0_rk
                Nt(i,jk) = R*Nt(i,jk-1) + L*Nt(i+1,jk-1)
                Rp = (Nt(i,jk-1) + d*dNt_dXt(i,jk-1)) / a
                if ( isnan(Rp) .or. isinf(Rp) .or. abs(Rp) < tiny(0.0_rk) ) Rp = 0.0_rk
                Lp = (b*dNt_dXt(i+1,jk-1) - Nt(i+1,jk-1)) / c
                if ( isnan(Lp) .or. isinf(Lp) .or. abs(Lp) < tiny(0.0_rk) ) Lp = 0.0_rk
                dNt_dXt(i,jk) = Rp + Lp
            end do
        end do
        dB = dNt_dXt(1:nc,k)
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function basis_bernstein(Xt, nc) result(B)
        real(rk), intent(in) :: Xt
        integer, intent(in) :: nc
        real(rk), allocatable :: B(:)
        integer  :: p, degree

        degree = nc - 1

        allocate(B(nc), source=0.0_rk)

        do concurrent (p = 0:degree)
            B(p+1) = gamma(real(nc, kind=rk))/(gamma(real(p+1, kind=rk))*gamma(real(nc-p, kind=rk)))
            if (Xt == 0.0_rk .and. p == 0) then
                B(p+1) = B(p+1)*(1.0_rk-Xt)**(degree-p)
            else if (Xt == 0.0_rk .and. degree-p == 0) then
                B(p+1) = B(p+1)*(Xt**p)
            else
                B(p+1) = B(p+1)*(Xt**p)*(1.0_rk-Xt)**(degree-p)
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
    pure function compute_multiplicity1(knot) result(multiplicity)
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
    pure function compute_multiplicity2(knot, Xth) result(multiplicity)
        real(rk), dimension(:), intent(in) :: knot
        real(rk), intent(in) :: Xth
        integer :: multiplicity
        integer :: i, count, size_knot

        size_knot = size(knot)
        multiplicity = 0
        i = 1
        do while (i <= size_knot)
            if (knot(i) == Xth) then
                count = 1
                do while (i + count <= size_knot .and. knot(i + count) == Xth)
                    count = count + 1
                end do
                if (count > multiplicity) then
                    multiplicity = count
                end if
                i = i + count
            else
                i = i + 1
            end if
        end do
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function compute_knot_vector(Xth_dir, degree, continuity) result(knot)
        real(rk), intent(in) :: Xth_dir(:)
        integer, intent(in) :: degree
        integer, intent(in) :: continuity(:)
        real(rk), allocatable :: knot(:)

        knot = repelem(Xth_dir, (degree - continuity))
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    elemental pure function isinf(x) result(output)
        real(rk), intent(in) :: x
        logical :: output

        output=.false.
        if (x >  huge(x)) output=.true.
        if (x < -huge(x)) output=.true.
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    elemental pure function isnan(x) result(output)
        real(rk), intent(in) :: x
        logical :: output

        output =.false.
        if (x /= x) output = .true.
    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine insert_knot_A_5_1(p, UP, Pw, u, k, s, r, nq, UQ, Qw)
        integer, intent(in) :: p, k, s, r
        real(rk), intent(in) :: UP(0:), Pw(0:,:)
        real(rk), intent(in) :: u
        real(rk), allocatable, intent(out) :: UQ(:), Qw(:,:)
        integer, intent(out) :: nq
        integer :: i, j, L, mp, dim, np
        real(rk), allocatable :: Rw(:,:)
        real(rk) :: alpha

        dim = size(Pw, 2)
        np  = size(Pw, 1) - 1
        mp  = np + p + 1
        nq  = np + r

        allocate(UQ(0:mp+r))
        allocate(Qw(0:nq,1:dim))
        allocate(Rw(0:p ,1:dim))

        UQ(0:k) = UP(0:k)
        UQ(k+1:k+r) = u
        UQ(k+1+r:mp+r) = UP(k+1:mp)
        Qw(0:k-p,:) = Pw(0:k-p,:)
        Qw(k-s+r:np+r,:) = Pw(k-s:np,:)
        Rw(0:p-s,:) = Pw(k-p:k-s,:)
        do j = 1, r
            L = k-p+j
            do i = 0, p-j-s
                alpha = (u - UP(L+i)) / (UP(i+k+1) - UP(L+i))
                Rw(i,:) = alpha*Rw(i+1,:) + (1.0_rk - alpha) * Rw(i,:)
            end do
            Qw(L,:) = Rw(0,:)
            Qw(k+r-j-s,:) = Rw(p-j-s,:)
        end do
        Qw(L+1:k-s-1,:) = Rw(1:k-s-1-L,:)
    end subroutine
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function findspan(n,degree,Xth,knot) result(s)
        integer, intent(in) :: n, degree
        real(rk), intent(in) :: Xth
        real(rk), intent(in) :: knot(:)
        integer :: s
        integer :: low, high, mid
        if (Xth == knot(n+2)) then
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
    pure subroutine elevate_degree_A_5_9(t, knot, degree, Xcw, nc_new, knot_new, Xcw_new)
        integer, intent(in) :: t
        real(rk), intent(in) :: Xcw(:,:), knot(:)
        integer, intent(in) :: degree
        integer, intent(out) :: nc_new
        real(rk), allocatable, intent(out) :: Xcw_new(:,:), knot_new(:)
        real(rk), allocatable :: bezalfs(:,:), bpts(:,:), ebpts(:,:), Nextbpts(:,:), alfs(:)
        real(rk) :: inv, alpha1, alpha2, Xth1, Xth2, numer, den
        integer :: n, lbz, rbz, sv, tr, kj, first, knoti, last, alpha3, dim, nc
        integer :: i, j, q, s, m, ph, ph2, mpi, mh, r, a, b, Xcwi, oldr, mul
        integer, allocatable :: mlp(:)

        nc = size(Xcw,1)
        dim = size(Xcw,2)
        mlp = compute_multiplicity(knot)
        mlp = mlp + t
        nc_new = sum(mlp) - (mlp(1)-1) - 1
        allocate(Xcw_new(nc_new,dim), source=0.0_rk)
        allocate(bezalfs(degree+1,degree+t+1), source=0.0_rk)
        allocate(bpts(degree+1,dim), source=0.0_rk)
        allocate(ebpts(degree+t+1,dim), source=0.0_rk)
        allocate(Nextbpts(degree+1,dim), source=0.0_rk)
        allocate(alfs(degree), source=0.0_rk)
        n = nc - 1
        m = n + degree + 1
        ph = degree + t
        ph2 = ph / 2
        bezalfs(1,1) = 1.0_rk
        bezalfs(degree+1,ph+1) = 1.0_rk
        do i = 1,ph2
            inv = 1.0_rk/bincoeff(ph,i)
            mpi = min(degree,i)
            do j = max(0,i-t),mpi
                bezalfs(j+1,i+1) = inv*bincoeff(degree,j)*bincoeff(t,i-j)
            end do
        end do
        do i = ph2+1,ph-1
            mpi = min(degree,i)
            do j = max(0,i-t),mpi
                bezalfs(j+1,i+1) = bezalfs(degree-j+1,ph-i+1)
            end do
        end do
        mh = ph
        knoti = ph+1
        r = -1
        a = degree
        b = degree+1
        Xcwi = 1
        Xth1 = knot(1)
        Xcw_new(1,:) = Xcw(1,:)
        allocate(knot_new(sum(mlp)), source=0.0_rk)
        knot_new(1:ph+1) = Xth1
        do i = 0,degree
            bpts(i+1,:) = Xcw(i+1,:)
        end do
        do while (b<m)
            i = b
            do while (b<m .and. knot(b+1) == knot(b+2))
                b = b + 1
                if (b+2 > size(knot)) then
                    exit
                end if
            end do
            mul = b - i + 1
            mh = mh + mul + t
            Xth2 = knot(b+1)
            oldr = r
            r = degree - mul
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
                do q = degree,mul+1,-1
                    alfs(q-mul) = numer / (knot(a+q+1)-Xth1)
                end do
                do j = 1,r
                    sv = r - j
                    s = mul + j
                    do q = degree,s,-1
                        bpts(q+1,:) = (1.0_rk-alfs(q-s+1))*bpts(q,:) + alfs(q-s+1)*bpts(q+1,:)
                    end do
                    Nextbpts(sv+1,:) = bpts(degree+1,:)
                end do
            end if
            do i = lbz,ph
                ebpts(i+1,:) = 0.0_rk
                mpi = min(degree,i)
                do j = max(0,i-t),mpi
                    ebpts(i+1,:) = bezalfs(j+1,i+1)*bpts(j+1,:) + ebpts(i+1,:)
                end do
            end do
            if (oldr > 1) then
                first = knoti - 2
                last = knoti
                den = Xth2 - Xth1
                alpha3 = floor((Xth2-knot(knoti)) / den)
                do tr = 1,oldr-1
                    i = first
                    j = last
                    kj = j - knoti + 1
                    do while (j-i > tr)
                        if (i < Xcwi) then
                            alpha1 = (Xth2-knot(i+1))/(Xth1-knot(i+1))
                            Xcw_new(i+1,:) = (1-alpha1)*Xcw_new(i,:) + alpha1*Xcw_new(i+1,:)
                        end if
                        if (j >= lbz) then
                            if (j-tr <= knoti-ph+oldr) then
                                alpha2 = (Xth2-knot_new(j-tr+1)) / den
                                ebpts(kj+1,:) = alpha2*ebpts(kj+1,:) + (1-alpha2)*ebpts(kj+2,:)
                            else
                                ebpts(kj+1,:) = (1-alpha3)*ebpts(:,kj+2) + alpha3*ebpts(kj+1,:)
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
            if (a /= degree) then
                do i = 0,ph-oldr-1
                    knot_new(knoti+1) = Xth1
                    knoti = knoti + 1
                end do
            end if
            do j = lbz,rbz
                Xcw_new(Xcwi+1,:) = ebpts(j+1,:)
                Xcwi = Xcwi + 1
            end do
            if (b<m) then
                do j = 0,r-1
                    bpts(j+1,:) = Nextbpts(j+1,:)
                end do
                do j = r,degree
                    bpts(j+1,:) = Xcw(b-degree+j+1,:)
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


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure function hexahedron_Xc(L, nc) result(Xc)
        real(rk), intent(in) :: L(3)
        integer, intent(in) :: nc(3)
        real(rk), dimension(:,:), allocatable :: Xc
        real(rk) :: dx, dy, dz
        integer :: i, j, k, nci

        dx = L(1) / real(nc(1)-1, rk)
        dy = L(2) / real(nc(2)-1, rk)
        dz = L(3) / real(nc(3)-1, rk)

        allocate(Xc(nc(1) * nc(2) * nc(3), 3))
        nci = 1
        do k = 0, nc(3)-1
            do j = 0, nc(2)-1
                do i = 0, nc(1)-1
                    Xc(nci, 1) = real(i,rk) * dx
                    Xc(nci, 2) = real(j,rk) * dy
                    Xc(nci, 3) = real(k,rk) * dz
                    nci = nci + 1
                end do
            end do
        end do

    end function
    !===============================================================================


    !===============================================================================
    !> author: Seyed Ali Ghasemi
    !> license: BSD 3-Clause
    pure subroutine remove_knots_A_5_8(p,knot,Pw,u,r,s,num,t, knot_new, Pw_new)
        real(rk), intent(in) :: u
        integer, intent(in) :: p, r, s, num
        real(rk), intent(in) :: knot(:)
        real(rk), intent(in) :: Pw(:,:)
        real(rk), allocatable, intent(out) :: knot_new(:)
        real(rk), allocatable, intent(out) :: Pw_new(:,:)
        real(rk), allocatable :: Pw_copy(:,:), knot_copy(:)
        integer, intent(out) :: t
        real(rk) :: tol, alfi, alfj
        real(rk), allocatable :: temp(:,:)
        integer :: i, j, ii, jj, remflag, off, first, last, ord, fout, m, k, n, nc, dim, tt

        dim = size(Pw,2)
        nc = size(Pw,1)
        n = nc
        m = n+p+1
        ord = p+1
        fout = (2*r-s-p)/2
        last = r-s
        first = r-p

        Pw_copy = Pw
        knot_copy = knot

        ! TODO:
        tol = 1.0e-6_rk * minval(Pw(:,dim))/(1.0_rk + maxval(sqrt(sum(Pw**2, 2))))

        allocate(temp(2*p+1,dim), source=0.0_rk)
        t = 0
        do tt = 0, num-1
            off = first-1
            temp(1,:) = Pw_copy(off,:)
            temp(last+1-off+1,:) = Pw_copy(last+1,:)
            i = first
            j = last
            ii = 1
            jj = last-off
            remflag = 0
            do while (j-i>t)
                alfi = (u-knot_copy(i))/(knot_copy(i+ord+t)-knot_copy(i))
                alfj = (u-knot_copy(j-t))/(knot_copy(j+ord)-knot_copy(j-t))
                temp(ii+1,:) = (Pw_copy(i,:)-(1.0_rk-alfi)*temp(ii-1+1,:))/alfi
                temp(jj+1,:) = (Pw_copy(j,:)-alfj*temp(jj+1+1,:))/(1.0_rk-alfj)
                i = i+1
                ii = ii+1
                j = j-1
                jj = jj-1
            end do
            if (j-i<=t) then
                if (norm2(temp(ii-1+1,:) - temp(jj+1+1,:))<=tol) then
                    remflag = 1
                else
                    alfi = (u-knot_copy(i))/(knot_copy(i+ord+t)-knot_copy(i))
                    if (norm2(Pw_copy(i,:) - (alfi*temp(ii+t+1+1,:)+(1.0_rk-alfi)*temp(ii-1+1,:)))<=tol) then
                        remflag = 1
                    end if
                end if
            end if
            if (remflag == 0) then
                exit
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
        end do
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

end module forcad_utils
