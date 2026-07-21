!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Beam-like block assembled from conforming multi-patch NURBS volumes.
program volume_multipatch_beam

    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume

    implicit none

    type(nurbs_volume) :: patch
    type(nurbs_multipatch_volume) :: beam
    integer, allocatable :: dof_map(:), elem(:,:)
    integer, parameter :: nx = 4, ny = 2, nz = 2
    integer :: id(nx, ny, nz), i, j, k
    real(rk), parameter :: length = 6.0_rk, width = 1.2_rk, height = 1.0_rk
    real(rk) :: dx, dy, dz

    dx = length/real(nx, rk)
    dy = width /real(ny, rk)
    dz = height/real(nz, rk)

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                call patch%set_hexahedron([dx, dy, dz], [2, 2, 2])
                call patch%translate_Xc([real(i - 1, rk)*dx, real(j - 1, rk)*dy, real(k - 1, rk)*dz])
                call beam%add_patch(patch, id(i, j, k))
            end do
        end do
    end do

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx - 1
                call beam%connect(id(i, j, k), "u_max", id(i + 1, j, k), "u_min", 0)
            end do
        end do
    end do
    do k = 1, nz
        do j = 1, ny - 1
            do i = 1, nx
                call beam%connect(id(i, j, k), "v_max", id(i, j + 1, k), "v_min", 0)
            end do
        end do
    end do
    do k = 1, nz - 1
        do j = 1, ny
            do i = 1, nx
                call beam%connect(id(i, j, k), "w_max", id(i, j, k + 1), "w_min", 0)
            end do
        end do
    end do

    dof_map = beam%cmp_dof_map()
    elem = beam%cmp_elem(shared=.true.)

    write(*,"(a,*(1x,g0))") "beam patches          :", beam%get_npatch()
    write(*,"(a,*(1x,g0))") "beam connections      :", beam%get_nconnection()
    write(*,"(a,*(1x,g0))") "beam valid            :", beam%is_valid()
    write(*,"(a,*(1x,g0))") "shared control nodes  :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "shared elements       :", size(elem, 1)
    call beam%create(5, 5, 5)
    call beam%export_Xc("vtk/volume_multipatch_beam_geometry")
    call beam%export_Xg("vtk/volume_multipatch_beam_geometry")
    call beam%export_Xth_in_Xg("vtk/volume_multipatch_beam_geometry", 5)
    call beam%show(&
        vtkfile_Xc        = "vtk/volume_multipatch_beam_geometry_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/volume_multipatch_beam_geometry_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_multipatch_beam_geometry_Xth_*.vtk")
    call beam%finalize()

end program volume_multipatch_beam
