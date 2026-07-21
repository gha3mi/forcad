!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Demonstrates a tapered, twisted multi-patch NURBS surface tower.
program surface_multipatch_tower

    use forcad, only: rk, nurbs_surface, nurbs_multipatch_surface

    implicit none

    integer, parameter :: nt = 12, nz = 6, nf = 6, ncp = 4
    real(rk), parameter :: pi = 3.1415926535897932384626433832795_rk
    real(rk), parameter :: height = 8.0_rk
    real(rk), parameter :: crown_height = 1.65_rk
    real(rk), parameter :: dtheta = 2.0_rk*pi/real(nt, rk)

    type(nurbs_surface) :: patch
    type(nurbs_multipatch_surface) :: patches
    real(rk), allocatable :: Xc(:,:)
    integer, allocatable :: dof_map(:), elem(:,:)
    integer :: shell_id(nt,nz+1), fin_id(nf)
    integer :: i, k, n
    real(rk) :: u, v, s, q, theta, phase, twist, radius, taper, flute

    do k = 1, nz + 1
        do i = 1, nt
            call patch%set_tetragon([1.0_rk, 1.0_rk], [ncp, ncp])
            Xc = patch%get_Xc()

            do n = 1, size(Xc, 1)
                u = Xc(n,1)
                v = Xc(n,2)
                theta = dtheta*(real(i - 1, rk) + u)

                if (k <= nz) then
                    s = (real(k - 1, rk) + v)/real(nz, rk)
                    taper = 1.62_rk - 0.72_rk*s**1.18_rk
                    flute = 0.10_rk*(1.0_rk - 0.25_rk*s)*cos(6.0_rk*theta)
                    radius = taper + flute + 0.055_rk*sin(8.0_rk*pi*s)
                    twist = 0.78_rk*pi*s
                    phase = theta + twist
                    Xc(n,1) = radius*cos(phase)
                    Xc(n,2) = radius*sin(phase)
                    Xc(n,3) = height*s + 0.05_rk*sin(3.0_rk*theta)*sin(pi*s)
                else
                    q = v
                    taper = 1.62_rk - 0.72_rk
                    flute = 0.10_rk*(1.0_rk - 0.25_rk)*cos(6.0_rk*theta)
                    radius = (taper + flute)*(1.0_rk - q)**1.35_rk
                    twist = 0.78_rk*pi + 0.20_rk*pi*q
                    phase = theta + twist
                    Xc(n,1) = radius*cos(phase)
                    Xc(n,2) = radius*sin(phase)
                    Xc(n,3) = height + crown_height*(1.0_rk - cos(0.5_rk*pi*q))
                end if
            end do

            call patch%set([ncp, ncp], Xc)
            call patches%add_patch(patch, shell_id(i,k))
        end do
    end do

    do i = 1, nf
        call patch%set_tetragon([1.0_rk, 1.0_rk], [ncp, ncp])
        Xc = patch%get_Xc()

        do n = 1, size(Xc, 1)
            q = Xc(n,1)
            s = 0.04_rk + 0.68_rk*Xc(n,2)
            theta = 2.0_rk*pi*real(i - 1, rk)/real(nf, rk)
            taper = 1.62_rk - 0.72_rk*s**1.18_rk
            flute = 0.10_rk*(1.0_rk - 0.25_rk*s)*cos(6.0_rk*theta)
            radius = taper + flute + 0.055_rk*sin(8.0_rk*pi*s) + 0.42_rk*q*(1.0_rk - s)**0.70_rk
            twist = 0.78_rk*pi*s
            phase = theta + twist
            Xc(n,1) = radius*cos(phase)
            Xc(n,2) = radius*sin(phase)
            Xc(n,3) = height*s
        end do

        call patch%set([ncp, ncp], Xc)
        call patches%add_patch(patch, fin_id(i))
    end do

    do k = 1, nz + 1
        do i = 1, nt - 1
            call patches%connect(shell_id(i,k), "u_max", shell_id(i+1,k), "u_min", 0)
        end do
        call patches%connect(shell_id(nt,k), "u_max", shell_id(1,k), "u_min", 0)
    end do

    do k = 1, nz
        do i = 1, nt
            call patches%connect(shell_id(i,k), "v_max", shell_id(i,k+1), "v_min", 0)
        end do
    end do

    dof_map = patches%cmp_dof_map()
    elem = patches%cmp_elem(shared=.true.)

    write(*,"(a,*(1x,g0))") "tower shell patches      :", nt*(nz + 1)
    write(*,"(a,*(1x,g0))") "tower buttress patches   :", size(fin_id)
    write(*,"(a,*(1x,g0))") "tower total patches      :", patches%get_npatch()
    write(*,"(a,*(1x,g0))") "tower connections        :", patches%get_nconnection()
    write(*,"(a,*(1x,g0))") "tower valid              :", patches%is_valid()
    write(*,"(a,*(1x,g0))") "local control dofs       :", size(dof_map)
    write(*,"(a,*(1x,g0))") "shared global dofs       :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "elements                 :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element        :", size(elem, 2)

    call patches%create(43, 61)
    call patches%export_Xc("vtk/surface_multipatch_tower")
    call patches%export_Xg("vtk/surface_multipatch_tower")
    call patches%export_Xth_in_Xg("vtk/surface_multipatch_tower", 61)
    call patches%show(&
        vtkfile_Xc        = "vtk/surface_multipatch_tower_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/surface_multipatch_tower_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_multipatch_tower_Xth_*.vtk")

end program surface_multipatch_tower
