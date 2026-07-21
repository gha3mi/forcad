!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Demonstrates a C2 multi-patch dome with an oculus, petals, and an inner drum.
program surface_multipatch_c2_dome

    use forcad, only: rk, nurbs_surface, nurbs_multipatch_surface

    implicit none

    integer, parameter :: ns = 18, nr = 4, ncp = 4
    real(rk), parameter :: pi = 3.1415926535897932384626433832795_rk
    real(rk), parameter :: inner = 0.70_rk, outer = 3.55_rk
    real(rk), parameter :: dtheta = 2.0_rk*pi/real(ns, rk)

    type(nurbs_surface) :: patch
    type(nurbs_multipatch_surface) :: patches
    real(rk), allocatable :: Xc(:,:)
    integer, allocatable :: dof_map(:), elem(:,:)
    integer :: roof_id(ns,nr), petal_id(ns), drum_id(ns)
    integer :: i, j, n
    real(rk) :: u, v, q, theta, radius, z, scallop

    do j = 1, nr
        do i = 1, ns
            call patch%set_tetragon([1.0_rk, 1.0_rk], [ncp, ncp])
            Xc = patch%get_Xc()

            do n = 1, size(Xc, 1)
                u = Xc(n,1)
                v = Xc(n,2)
                q = (real(j - 1, rk) + v)/real(nr, rk)
                theta = dtheta*(real(i - 1, rk) + u)
                scallop = q*(1.0_rk - q)
                radius = inner + (outer - inner)*q + 0.12_rk*scallop*cos(9.0_rk*theta)
                z = 1.85_rk*(1.0_rk - q)**0.70_rk + &
                    0.22_rk*scallop*cos(6.0_rk*theta) + 0.08_rk*scallop*sin(12.0_rk*theta)
                Xc(n,:) = [radius*cos(theta), radius*sin(theta), z]
            end do

            call patch%set([ncp, ncp], Xc)
            call patches%add_patch(patch, roof_id(i,j))
        end do
    end do

    do i = 1, ns
        call patch%set_tetragon([1.0_rk, 1.0_rk], [ncp, ncp])
        Xc = patch%get_Xc()

        do n = 1, size(Xc, 1)
            u = Xc(n,1)
            q = Xc(n,2)
            theta = dtheta*(real(i - 1, rk) + u)
            scallop = sin(pi*u)
            radius = outer + q*(0.70_rk + 1.25_rk*scallop)
            z = 0.10_rk*q + (0.35_rk + 0.60_rk*scallop)*sin(pi*q) + &
                0.42_rk*q*q*scallop
            Xc(n,:) = [radius*cos(theta), radius*sin(theta), z]
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
            theta = dtheta*(real(i - 1, rk) + u)
            radius = inner*(1.0_rk - 0.10_rk*q*sin(pi*u))
            z = 1.85_rk - 1.72_rk*q + 0.08_rk*sin(pi*q)*cos(3.0_rk*theta)
            Xc(n,:) = [radius*cos(theta), radius*sin(theta), z]
        end do

        call patch%set([ncp, ncp], Xc)
        call patches%add_patch(patch, drum_id(i))
    end do

    do j = 1, nr
        do i = 1, ns - 1
            call patches%connect(roof_id(i,j), "u_max", roof_id(i+1,j), "u_min", 2)
        end do
        call patches%connect(roof_id(ns,j), "u_max", roof_id(1,j), "u_min", 2)
    end do

    do j = 1, nr - 1
        do i = 1, ns
            call patches%connect(roof_id(i,j), "v_max", roof_id(i,j+1), "v_min", 2)
        end do
    end do

    do i = 1, ns - 1
        call patches%connect(petal_id(i), "u_max", petal_id(i+1), "u_min", 2)
        call patches%connect(drum_id(i), "u_max", drum_id(i+1), "u_min", 2)
    end do
    call patches%connect(petal_id(ns), "u_max", petal_id(1), "u_min", 2)
    call patches%connect(drum_id(ns), "u_max", drum_id(1), "u_min", 2)

    do i = 1, ns
        call patches%connect(roof_id(i,nr), "v_max", petal_id(i), "v_min", 2)
        call patches%connect(roof_id(i,1), "v_min", drum_id(i), "v_min", 0)
    end do

    dof_map = patches%cmp_dof_map()
    elem = patches%cmp_elem(shared=.true.)

    write(*,"(a,*(1x,g0))") "C2 dome roof patches     :", ns*nr
    write(*,"(a,*(1x,g0))") "C2 dome petal patches    :", size(petal_id)
    write(*,"(a,*(1x,g0))") "C2 dome drum patches     :", size(drum_id)
    write(*,"(a,*(1x,g0))") "C2 dome total patches    :", patches%get_npatch()
    write(*,"(a,*(1x,g0))") "C2 dome connections      :", patches%get_nconnection()
    write(*,"(a,*(1x,g0))") "C2 dome valid            :", patches%is_valid()
    write(*,"(a,*(1x,g0))") "local control dofs       :", size(dof_map)
    write(*,"(a,*(1x,g0))") "shared global dofs       :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "elements                 :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element        :", size(elem, 2)

    call patches%create(33, 33)
    call patches%export_Xc("vtk/surface_multipatch_c2_dome")
    call patches%export_Xg("vtk/surface_multipatch_c2_dome")
    call patches%export_Xth_in_Xg("vtk/surface_multipatch_c2_dome", 33)
    call patches%show(&
        vtkfile_Xc        = "vtk/surface_multipatch_c2_dome_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/surface_multipatch_c2_dome_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_multipatch_c2_dome_Xth_*.vtk")

end program surface_multipatch_c2_dome
