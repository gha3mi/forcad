!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Demonstrates a curved multi-patch NURBS surface courtyard with an opening.
program surface_multipatch_courtyard

    use forcad, only: rk, nurbs_surface, nurbs_multipatch_surface

    implicit none

    integer, parameter :: nx = 3, ny = 3
    real(rk), parameter :: pi = 3.1415926535897932384626433832795_rk

    type(nurbs_surface) :: patch
    type(nurbs_multipatch_surface) :: patches
    real(rk), allocatable :: Xc(:,:)
    real(rk) :: x, y, r2, scale
    integer, allocatable :: dof_map(:), elem(:,:)
    integer :: id(nx,ny), i, j, n

    id = 0
    do j = 1, ny
        do i = 1, nx
            if (i == 2 .and. j == 2) cycle

            call patch%set_tetragon([1.0_rk, 1.0_rk], [3, 3])
            call patch%translate_Xc([real(i - 2, rk), real(j - 2, rk), 0.0_rk])

            Xc = patch%get_Xc()
            do n = 1, size(Xc, 1)
                x = Xc(n,1) - 0.5_rk
                y = Xc(n,2) - 0.5_rk
                r2 = x*x + y*y
                scale = 1.0_rk + 0.10_rk*exp(-0.65_rk*r2)
                Xc(n,1) = 0.5_rk + scale*x
                Xc(n,2) = 0.5_rk + scale*y
                Xc(n,3) = 0.18_rk + 0.28_rk*exp(-0.50_rk*r2) + 0.10_rk*sin(pi*Xc(n,1))*cos(pi*Xc(n,2))
            end do
            call patch%set([3, 3], Xc)

            call patches%add_patch(patch, id(i,j))
        end do
    end do

    do j = 1, ny
        do i = 1, nx - 1
            if (id(i,j) > 0 .and. id(i+1,j) > 0) then
                call patches%connect(id(i,j), "u_max", id(i+1,j), "u_min", 0)
            end if
        end do
    end do
    do j = 1, ny - 1
        do i = 1, nx
            if (id(i,j) > 0 .and. id(i,j+1) > 0) then
                call patches%connect(id(i,j), "v_max", id(i,j+1), "v_min", 0)
            end if
        end do
    end do

    dof_map = patches%cmp_dof_map()
    elem = patches%cmp_elem(shared=.true.)

    write(*,"(a,*(1x,g0))") "courtyard patches        :", patches%get_npatch()
    write(*,"(a,*(1x,g0))") "courtyard connections    :", patches%get_nconnection()
    write(*,"(a,*(1x,g0))") "courtyard valid          :", patches%is_valid()
    write(*,"(a,*(1x,g0))") "local control dofs       :", size(dof_map)
    write(*,"(a,*(1x,g0))") "shared global dofs       :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "elements                 :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element        :", size(elem, 2)
    call patches%create(43, 43)
    call patches%export_Xc("vtk/surface_multipatch_courtyard")
    call patches%export_Xg("vtk/surface_multipatch_courtyard")
    call patches%export_Xth_in_Xg("vtk/surface_multipatch_courtyard", 43)
    call patches%show(&
        vtkfile_Xc        = "vtk/surface_multipatch_courtyard_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/surface_multipatch_courtyard_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_multipatch_courtyard_Xth_*.vtk")

end program surface_multipatch_courtyard
