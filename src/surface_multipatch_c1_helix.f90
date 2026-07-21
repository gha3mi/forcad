!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Demonstrates a C1 multi-patch helical ramp with smooth rails.
program surface_multipatch_c1_helix

    use forcad, only: rk, nurbs_surface, nurbs_multipatch_surface

    implicit none

    integer, parameter :: ns = 24, nb = 2, ncp = 4
    real(rk), parameter :: pi = 3.1415926535897932384626433832795_rk
    real(rk), parameter :: turns = 1.65_rk, inner = 1.05_rk, outer = 2.65_rk
    real(rk), parameter :: height = 5.2_rk

    type(nurbs_surface) :: patch
    type(nurbs_multipatch_surface) :: patches
    real(rk), allocatable :: Xc(:,:)
    integer, allocatable :: dof_map(:), elem(:,:)
    integer :: ramp_id(ns,nb), inner_rail_id(ns), outer_rail_id(ns)
    integer :: i, j, n
    real(rk) :: u, v, q, t, w, theta, radius, z, base_z

    do j = 1, nb
        do i = 1, ns
            call patch%set_tetragon([1.0_rk, 1.0_rk], [ncp, ncp])
            Xc = patch%get_Xc()

            do n = 1, size(Xc, 1)
                u = Xc(n,1)
                v = Xc(n,2)
                t = (real(i - 1, rk) + u)/real(ns, rk)
                w = (real(j - 1, rk) + v)/real(nb, rk)
                theta = 2.0_rk*pi*turns*t
                radius = inner + (outer - inner)*w + 0.10_rk*sin(pi*w)*sin(5.0_rk*pi*t)
                z = height*t + 0.18_rk*sin(pi*w)*cos(4.0_rk*pi*t)
                Xc(n,:) = [radius*cos(theta), radius*sin(theta), z]
            end do

            call patch%set([ncp, ncp], Xc)
            call patches%add_patch(patch, ramp_id(i,j))
        end do
    end do

    do i = 1, ns
        call patch%set_tetragon([1.0_rk, 1.0_rk], [ncp, ncp])
        Xc = patch%get_Xc()

        do n = 1, size(Xc, 1)
            u = Xc(n,1)
            q = Xc(n,2)
            t = (real(i - 1, rk) + u)/real(ns, rk)
            theta = 2.0_rk*pi*turns*t
            base_z = height*t
            radius = inner - 0.10_rk*q*(1.0_rk - 0.35_rk*sin(2.0_rk*pi*t))
            z = base_z + 0.62_rk*q + 0.08_rk*sin(pi*q)*sin(3.0_rk*pi*t)
            Xc(n,:) = [radius*cos(theta), radius*sin(theta), z]
        end do

        call patch%set([ncp, ncp], Xc)
        call patches%add_patch(patch, inner_rail_id(i))

        call patch%set_tetragon([1.0_rk, 1.0_rk], [ncp, ncp])
        Xc = patch%get_Xc()

        do n = 1, size(Xc, 1)
            u = Xc(n,1)
            q = Xc(n,2)
            t = (real(i - 1, rk) + u)/real(ns, rk)
            theta = 2.0_rk*pi*turns*t
            base_z = height*t
            radius = outer + 0.14_rk*q*(1.0_rk + 0.20_rk*cos(2.0_rk*pi*t))
            z = base_z + 0.86_rk*q + 0.10_rk*sin(pi*q)*cos(3.0_rk*pi*t)
            Xc(n,:) = [radius*cos(theta), radius*sin(theta), z]
        end do

        call patch%set([ncp, ncp], Xc)
        call patches%add_patch(patch, outer_rail_id(i))
    end do

    do j = 1, nb
        do i = 1, ns - 1
            call patches%connect(ramp_id(i,j), "u_max", ramp_id(i+1,j), "u_min", 1)
        end do
    end do

    do i = 1, ns
        call patches%connect(ramp_id(i,1), "v_max", ramp_id(i,2), "v_min", 1)
        call patches%connect(ramp_id(i,1), "v_min", inner_rail_id(i), "v_min", 0)
        call patches%connect(ramp_id(i,2), "v_max", outer_rail_id(i), "v_min", 0)
    end do

    do i = 1, ns - 1
        call patches%connect(inner_rail_id(i), "u_max", inner_rail_id(i+1), "u_min", 1)
        call patches%connect(outer_rail_id(i), "u_max", outer_rail_id(i+1), "u_min", 1)
    end do

    dof_map = patches%cmp_dof_map()
    elem = patches%cmp_elem(shared=.true.)

    write(*,"(a,*(1x,g0))") "C1 helix ramp patches    :", ns*nb
    write(*,"(a,*(1x,g0))") "C1 helix rail patches    :", 2*ns
    write(*,"(a,*(1x,g0))") "C1 helix total patches   :", patches%get_npatch()
    write(*,"(a,*(1x,g0))") "C1 helix connections     :", patches%get_nconnection()
    write(*,"(a,*(1x,g0))") "C1 helix valid           :", patches%is_valid()
    write(*,"(a,*(1x,g0))") "local control dofs       :", size(dof_map)
    write(*,"(a,*(1x,g0))") "shared global dofs       :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "elements                 :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element        :", size(elem, 2)

    call patches%create(29, 31)
    call patches%export_Xc("vtk/surface_multipatch_c1_helix")
    call patches%export_Xg("vtk/surface_multipatch_c1_helix")
    call patches%export_Xth_in_Xg("vtk/surface_multipatch_c1_helix", 31)
    call patches%show(&
        vtkfile_Xc        = "vtk/surface_multipatch_c1_helix_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/surface_multipatch_c1_helix_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_multipatch_c1_helix_Xth_*.vtk")

end program surface_multipatch_c1_helix
