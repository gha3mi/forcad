!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Loft compatible rational ring sections into a curved annular NURBS volume.
program volume_loft

    use forcad, only: rk, nurbs_surface, nurbs_volume, loft

    implicit none

    integer, parameter :: nsection = 5
    real(rk), parameter :: center_x(nsection) = [0.0_rk,0.25_rk,-0.15_rk,0.35_rk,0.10_rk]
    real(rk), parameter :: center_y(nsection) = [0.0_rk,-0.18_rk,0.22_rk,-0.10_rk,0.15_rk]
    real(rk), parameter :: center_z(nsection) = [0.0_rk,1.2_rk,2.7_rk,4.5_rk,6.2_rk]
    real(rk), parameter :: inner_radius(nsection) = [0.42_rk,0.55_rk,0.48_rk,0.62_rk,0.46_rk]
    real(rk), parameter :: outer_radius(nsection) = [0.92_rk,1.18_rk,1.05_rk,1.32_rk,0.98_rk]

    type(nurbs_surface) :: sections(nsection)
    type(nurbs_volume) :: duct
    integer :: i

    do i = 1, nsection
        call sections(i)%set_ring(&
            center  = [center_x(i),center_y(i),center_z(i)],&
            radius1 = inner_radius(i),&
            radius2 = outer_radius(i))
    end do
    duct = loft(sections=sections)

    write(*,"(a,l1)") "Valid loft             : ", duct%err%ok
    write(*,"(a,l1)") "Rational geometry      : ", duct%is_rational()
    write(*,"(a,3(i0,1x))") "Degrees                 : ", duct%get_degree()
    write(*,"(a,3(i0,1x))") "Control points          : ", duct%get_nc()

    call duct%create(61, 13, 71)
    call duct%export_Xc("vtk/volume_loft_Xc.vtk")
    call duct%export_Xg("vtk/volume_loft_Xg.vtk")
    call duct%export_Xth_in_Xg("vtk/volume_loft_Xth.vtk", 13)
    call duct%show(&
        vtkfile_Xc        = "vtk/volume_loft_Xc.vtk",&
        vtkfile_Xg        = "vtk/volume_loft_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_loft_Xth.vtk")

    do i = 1, nsection
        call sections(i)%finalize()
    end do
    call duct%finalize()

end program volume_loft
