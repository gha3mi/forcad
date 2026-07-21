!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Loft five compatible rational circle sections into a cubic NURBS surface.
program surface_loft

    use forcad, only: rk, nurbs_curve, nurbs_surface, loft

    implicit none

    integer, parameter :: nsection = 5
    real(rk), parameter :: pi = acos(-1.0_rk)
    real(rk), parameter :: parameters(nsection) = [0.0_rk,0.16_rk,0.43_rk,0.74_rk,1.0_rk]
    real(rk), parameter :: radius(nsection) = [0.65_rk,1.05_rk,0.82_rk,1.28_rk,0.72_rk]

    type(nurbs_curve) :: sections(nsection)
    type(nurbs_surface) :: shell
    real(rk) :: section_coordinate
    integer :: i

    do i = 1, nsection
        section_coordinate = real(i-1,rk)/real(nsection-1,rk)
        call sections(i)%set_circle(&
            center = [0.35_rk*sin(1.5_rk*pi*section_coordinate), &
                0.22_rk*cos(2.0_rk*pi*section_coordinate),5.0_rk*section_coordinate],&
            radius = radius(i))
    end do
    shell = loft(&
        sections   = sections,&
        parameters = parameters,&
        degree     = 3)

    write(*,"(a,l1)") "Valid loft             : ", shell%err%ok
    write(*,"(a,l1)") "Rational geometry      : ", shell%is_rational()
    write(*,"(a,2(i0,1x))") "Degrees                 : ", shell%get_degree()
    write(*,"(a,2(i0,1x))") "Control points          : ", shell%get_nc()

    call shell%create(81, 61)
    call shell%export_Xc("vtk/surface_loft_Xc.vtk")
    call shell%export_Xg("vtk/surface_loft_Xg.vtk")
    call shell%export_Xth_in_Xg("vtk/surface_loft_Xth.vtk", 17)
    call shell%show(&
        vtkfile_Xc        = "vtk/surface_loft_Xc.vtk",&
        vtkfile_Xg        = "vtk/surface_loft_Xg.vtk",&
        vtkfile_Xth_in_Xg = "vtk/surface_loft_Xth.vtk")

    do i = 1, nsection
        call sections(i)%finalize()
    end do
    call shell%finalize()

end program surface_loft
