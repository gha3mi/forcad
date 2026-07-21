!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Demonstrates a radial multi-patch pavilion with an oculus and petal panels.
program surface_multipatch_pavilion

    use forcad, only: rk, nurbs_surface, nurbs_multipatch_surface

    implicit none

    integer, parameter :: ns = 16, nr = 3, ncp = 4
    real(rk), parameter :: pi = 3.1415926535897932384626433832795_rk
    real(rk), parameter :: inner = 0.72_rk, outer = 2.75_rk
    real(rk), parameter :: dtheta = 2.0_rk*pi/real(ns, rk)

    type(nurbs_surface) :: patch
    type(nurbs_multipatch_surface) :: patches
    real(rk), allocatable :: Xc(:,:)
    integer, allocatable :: dof_map(:), elem(:,:)
    integer :: ring_id(ns,nr), petal_id(ns), column_id(ns)
    integer :: i, j, n
    real(rk) :: u, v, q, theta, theta_warp, radius, z, petal_lift

    do j = 1, nr
        do i = 1, ns
            call patch%set_tetragon([1.0_rk, 1.0_rk], [ncp, ncp])
            Xc = patch%get_Xc()

            do n = 1, size(Xc, 1)
                u = Xc(n,1)
                v = Xc(n,2)
                q = (real(j - 1, rk) + v)/real(nr, rk)
                theta = dtheta*(real(i - 1, rk) + u)
                radius = inner + (outer - inner)*q + 0.10_rk*q*(1.0_rk - q)*cos(8.0_rk*theta)
                z = 0.22_rk*sin(pi*q) + 0.14_rk*cos(4.0_rk*theta)*q*(1.0_rk - q) + &
                    0.05_rk*sin(12.0_rk*theta)*q*(1.0_rk - q)
                Xc(n,:) = [radius*cos(theta), radius*sin(theta), z]
            end do

            call patch%set([ncp, ncp], Xc)
            call patches%add_patch(patch, ring_id(i,j))
        end do
    end do

    do i = 1, ns
        call patch%set_tetragon([1.0_rk, 1.0_rk], [ncp, ncp])
        Xc = patch%get_Xc()

        do n = 1, size(Xc, 1)
            u = Xc(n,1)
            q = Xc(n,2)
            theta = dtheta*(real(i - 1, rk) + u)
            theta_warp = theta + 0.08_rk*q*(1.0_rk - q)*sin(4.0_rk*theta)
            petal_lift = sin(pi*u)
            radius = outer + q*(0.82_rk + 1.05_rk*petal_lift) + &
                0.10_rk*q*(1.0_rk - q)*cos(3.0_rk*theta)
            z = 0.12_rk*q + (0.35_rk + 0.55_rk*petal_lift)*sin(pi*q) + &
                0.55_rk*q*q*petal_lift
            Xc(n,:) = [radius*cos(theta_warp), radius*sin(theta_warp), z]
        end do

        call patch%set([ncp, ncp], Xc)
        call patches%add_patch(patch, petal_id(i))
    end do

    do i = 1, ns
        call patch%set_tetragon([1.0_rk, 1.0_rk], [ncp, ncp])
        Xc = patch%get_Xc()

        do n = 1, size(Xc, 1)
            u = Xc(n,1)
            q = Xc(n,2)
            theta = dtheta*(real(i, rk) - 0.5_rk) + 0.10_rk*(u - 0.5_rk)
            radius = inner*(0.74_rk + 0.18_rk*u)
            z = -1.35_rk*q + 0.08_rk*sin(pi*q)
            Xc(n,:) = [radius*cos(theta), radius*sin(theta), z]
        end do

        call patch%set([ncp, ncp], Xc)
        call patches%add_patch(patch, column_id(i))
    end do

    do j = 1, nr
        do i = 1, ns - 1
            call patches%connect(ring_id(i,j), "u_max", ring_id(i+1,j), "u_min", 0)
        end do
        call patches%connect(ring_id(ns,j), "u_max", ring_id(1,j), "u_min", 0)
    end do

    do j = 1, nr - 1
        do i = 1, ns
            call patches%connect(ring_id(i,j), "v_max", ring_id(i,j+1), "v_min", 0)
        end do
    end do

    do i = 1, ns - 1
        call patches%connect(petal_id(i), "u_max", petal_id(i+1), "u_min", 0)
    end do
    call patches%connect(petal_id(ns), "u_max", petal_id(1), "u_min", 0)

    do i = 1, ns
        call patches%connect(ring_id(i,nr), "v_max", petal_id(i), "v_min", 0)
    end do

    dof_map = patches%cmp_dof_map()
    elem = patches%cmp_elem(shared=.true.)

    write(*,"(a,*(1x,g0))") "pavilion ring patches    :", ns*nr
    write(*,"(a,*(1x,g0))") "pavilion petal patches   :", size(petal_id)
    write(*,"(a,*(1x,g0))") "pavilion column patches  :", size(column_id)
    write(*,"(a,*(1x,g0))") "pavilion total patches   :", patches%get_npatch()
    write(*,"(a,*(1x,g0))") "pavilion connections     :", patches%get_nconnection()
    write(*,"(a,*(1x,g0))") "pavilion valid           :", patches%is_valid()
    write(*,"(a,*(1x,g0))") "local control dofs       :", size(dof_map)
    write(*,"(a,*(1x,g0))") "shared global dofs       :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "elements                 :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element        :", size(elem, 2)

    call patches%create(35, 35)
    call patches%export_Xc("vtk/surface_multipatch_pavilion")
    call patches%export_Xg("vtk/surface_multipatch_pavilion")
    call patches%export_Xth_in_Xg("vtk/surface_multipatch_pavilion", 35)
    call patches%show(&
        vtkfile_Xc        = "vtk/surface_multipatch_pavilion_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/surface_multipatch_pavilion_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_multipatch_pavilion_Xth_*.vtk")

end program surface_multipatch_pavilion
