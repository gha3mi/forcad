module forcad_utils

    implicit none

    private
    public :: rk, basis_bernstein, basis_bspline, elemConn_C0, kron, ndgrid, compute_multiplicity, compute_knot_vector, &
        basis_bspline_der, insert_knot_A_5_1, findspan, elevate_degree_A

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
    pure function basis_bspline_der(Xt, knot, nc, order) result(dB)
        integer, intent(in)   :: order
        real(rk), intent(in)  :: knot(:)
        integer, intent(in)   :: nc
        real(rk), intent(in)  :: Xt
        real(rk), allocatable :: dB(:)
        real(rk), allocatable :: Nt(:,:), dNt_dXt(:,:)
        real(rk)              :: R, L, Rp, Lp, knot_i, knot_ip, knot_jk, knot_jkm, knot_end, a, b, c, d
        integer               :: i, p, k, n, m, jk

        k = order + 1
        n = nc - 1
        allocate(Nt(nc+order, order+1))
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
        allocate(dNt_dXt(nc+order, order+1))
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
    pure subroutine elevate_degree_A(t, knot, order, Xcw, nc_new, knot_new, Xcw_new)
        integer, intent(in) :: t
        real(rk), intent(in) :: Xcw(:,:), knot(:)
        integer, intent(in) :: order
        integer, intent(out) :: nc_new
        real(rk), allocatable, intent(out) :: Xcw_new(:,:), knot_new(:)
        real(rk), allocatable :: bezalfs(:,:), bpts(:,:), ebpts(:,:), Nextbpts(:,:), alfs(:)
        real(rk) :: inv, alpha1, alpha2, Xth1, Xth2, numer, den
        integer :: n, lbz, rbz, sv, tr, kj, first, knoti, last, alpha3, ii, dim, nc
        integer :: i, j, q, s, m, ph, ph2, mpi, mh, r, a, b, Xcwi, oldr, mul
        integer, allocatable :: mlp(:)

        nc = size(Xcw,1)
        dim = size(Xcw,2)
        mlp = compute_multiplicity(knot)
        mlp = mlp + t
        nc_new = sum(mlp) - (mlp(1)-1) - 1
        allocate(Xcw_new(nc_new,dim), source=0.0_rk)
        allocate(bezalfs(order+1,order+t+1), source=0.0_rk)
        allocate(bpts(order+1,dim), source=0.0_rk)
        allocate(ebpts(order+t+1,dim), source=0.0_rk)
        allocate(Nextbpts(order+1,dim), source=0.0_rk)
        allocate(alfs(order), source=0.0_rk)
        n = nc - 1
        m = n + order + 1
        ph = order + t
        ph2 = ph / 2
        bezalfs(1,1) = 1.0_rk
        bezalfs(order+1,ph+1) = 1.0_rk
        do i = 1,ph2
            inv = 1.0_rk/bincoeff(ph,i)
            mpi = min(order,i)
            do j = max(0,i-t),mpi
                bezalfs(j+1,i+1) = inv*bincoeff(order,j)*bincoeff(t,i-j)
            end do
        end do
        do i = ph2+1,ph-1
            mpi = min(order,i)
            do j = max(0,i-t),mpi
                bezalfs(j+1,i+1) = bezalfs(order-j+1,ph-i+1)
            end do
        end do
        mh = ph
        knoti = ph+1
        r = -1
        a = order
        b = order+1
        Xcwi = 1
        Xth1 = knot(1)
        do ii =0,dim-1
            Xcw_new(1,ii+1) = Xcw(1,ii+1)
        end do
        allocate(knot_new(sum(mlp)), source=0.0_rk)
        do i = 0,ph
            knot_new(i+1) = Xth1
        end do
        do i = 0,order
            do ii = 0,dim-1
                bpts(i+1,ii+1) = Xcw(i+1,ii+1)
            end do
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
            r = order - mul
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
                do q = order,mul+1,-1
                    alfs(q-mul) = numer / (knot(a+q+1)-Xth1)
                end do
                do j = 1,r
                    sv = r - j
                    s = mul + j
                    do q = order,s,-1
                        do ii = 0,dim-1
                            bpts(q+1,ii+1) = (1.0_rk-alfs(q-s+1))*bpts(q,ii+1) + alfs(q-s+1)*bpts(q+1,ii+1)
                        end do
                    end do
                    do ii = 0,dim-1
                        Nextbpts(sv+1,ii+1) = bpts(order+1,ii+1)
                    end do
                end do
            end if
            do i = lbz,ph
                do ii = 0,dim-1
                    ebpts(i+1,ii+1) = 0.0_rk
                end do
                mpi = min(order,i)
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
                alpha3 = floor((Xth2-knot(knoti)) / den)
                do tr = 1,oldr-1
                    i = first
                    j = last
                    kj = j - knoti + 1
                    do while (j-i > tr)
                        if (i < Xcwi) then
                            alpha1 = (Xth2-knot(i+1))/(Xth1-knot(i+1))
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
            if (a /= order) then
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
                do j = r,order
                    do ii = 0,dim-1
                        bpts(j+1,ii+1) = Xcw(b-order+j+1,ii+1)
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

end module forcad_utils
