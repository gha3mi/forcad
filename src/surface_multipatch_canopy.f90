!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Demonstrates a curved 4-by-3 multi-patch NURBS surface canopy.
program surface_multipatch_canopy

    use forcad, only: rk, nurbs_surface, nurbs_multipatch_surface

    implicit none

    integer, parameter :: nx = 4, ny = 3
    real(rk), parameter :: pi = 3.1415926535897932384626433832795_rk

    type(nurbs_surface) :: patch
    type(nurbs_multipatch_surface) :: patches
    real(rk), allocatable :: Xc(:,:)
    integer, allocatable :: dof_map(:), elem(:,:)
    integer :: id(nx,ny), i, j, n

    do j = 1, ny
        do i = 1, nx
            call patch%set_tetragon([1.0_rk, 1.0_rk], [3, 3])
            call patch%translate_Xc([real(i - 1, rk), real(j - 1, rk), 0.0_rk])

            Xc = patch%get_Xc()
            do concurrent (n = 1:size(Xc, 1))
                Xc(n,3) = 0.45_rk*sin(pi*Xc(n,1)/real(nx, rk))*sin(pi*Xc(n,2)/real(ny, rk)) + &
                    0.12_rk*cos(2.0_rk*pi*(Xc(n,1) + 0.5_rk*Xc(n,2))/real(nx + ny, rk))
            end do
            call patch%set([3, 3], Xc)

            call patches%add_patch(patch, id(i,j))
        end do
    end do

    do j = 1, ny
        do i = 1, nx - 1
            call patches%connect(id(i,j), "u_max", id(i+1,j), "u_min", 0)
        end do
    end do
    do j = 1, ny - 1
        do i = 1, nx
            call patches%connect(id(i,j), "v_max", id(i,j+1), "v_min", 0)
        end do
    end do

    dof_map = patches%cmp_dof_map()
    elem = patches%cmp_elem(shared=.true.)

    write(*,"(a,*(1x,g0))") "canopy patches           :", patches%get_npatch()
    write(*,"(a,*(1x,g0))") "canopy connections       :", patches%get_nconnection()
    write(*,"(a,*(1x,g0))") "canopy valid             :", patches%is_valid()
    write(*,"(a,*(1x,g0))") "local control dofs       :", size(dof_map)
    write(*,"(a,*(1x,g0))") "shared global dofs       :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "expected global dofs     :", (2*nx + 1)*(2*ny + 1)
    write(*,"(a,*(1x,g0))") "elements                 :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element        :", size(elem, 2)
    call patches%create(45, 39)
    call patches%export_Xc("vtk/surface_multipatch_canopy")
    call patches%export_Xg("vtk/surface_multipatch_canopy")
    call patches%export_Xth_in_Xg("vtk/surface_multipatch_canopy", 45)
    call patches%show(&
        vtkfile_Xc        = "vtk/surface_multipatch_canopy_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/surface_multipatch_canopy_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_multipatch_canopy_Xth_*.vtk")

end program surface_multipatch_canopy
