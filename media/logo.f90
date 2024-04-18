program logo

    use forcad, only: rk, nurbs_volume

    implicit none
    type(nurbs_volume) :: shape

    !> F
    call shape%set_hexahedron(L=[1.0_rk, 1.0_rk, 7.0_rk], nc=[2,2,2])
    call shape%export_Xc('vtk/logo_Xc_1.vtk')
    call shape%create(3*3, 3*3, 21*3)
    call shape%export_Xg('vtk/logo_Xg_1.vtk')
    call shape%finalize()

    call shape%set_hexahedron(L=[2.0_rk, 1.0_rk, 1.0_rk], nc=[2,2,2])
    call shape%translate_Xc([1.0_rk, 0.0_rk, 6.0_rk])
    call shape%export_Xc('vtk/logo_Xc_2.vtk')
    call shape%create(6*3, 3*3, 3*3)
    call shape%export_Xg('vtk/logo_Xg_2.vtk')
    call shape%finalize()

    call shape%set_hexahedron(L=[1.0_rk, 1.0_rk, 1.0_rk], nc=[2,2,2])
    call shape%translate_Xc([1.0_rk, 0.0_rk, 3.0_rk])
    call shape%export_Xc('vtk/logo_Xc_3.vtk')
    call shape%create(3*3, 3*3, 3*3)
    call shape%export_Xg('vtk/logo_Xg_3.vtk')
    call shape%finalize()

    !> o
    call shape%set_ring([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 2.0_rk, 1.0_rk)
    call shape%rotate_Xc(90.0_rk, 0.0_rk, 0.0_rk)
    call shape%translate_Xc([5.0_rk, 1.0_rk, 2.0_rk])
    call shape%export_Xc('vtk/logo_Xc_4.vtk')
    call shape%create(60, 15, 10)
    call shape%export_Xg('vtk/logo_Xg_4.vtk')
    call shape%finalize()

    !> r
    call shape%set_hexahedron(L=[1.0_rk, 1.0_rk, 4.0_rk], nc=[2,2,2])
    call shape%translate_Xc([8.0_rk, 0.0_rk, 0.0_rk])
    call shape%export_Xc('vtk/logo_Xc_5.vtk')
    call shape%create(3*3, 3*3, 6*3)
    call shape%export_Xg('vtk/logo_Xg_5.vtk')
    call shape%finalize()

    call shape%set_hexahedron(L=[1.0_rk, 1.0_rk, 4.0_rk], nc=[2,2,4])
    call shape%modify_Xc(1.8_rk,13,1)
    call shape%modify_Xc(2.8_rk,14,1)
    call shape%modify_Xc(1.8_rk,15,1)
    call shape%modify_Xc(2.8_rk,16,1)
    call shape%modify_Xc(3.55_rk,13,3)
    call shape%modify_Xc(3.55_rk,15,3)
    
    call shape%translate_Xc([8.0_rk, 0.0_rk, 0.0_rk])
    call shape%export_Xc('vtk/logo_Xc_6.vtk')
    call shape%create(3*3, 3*3, 6*3)
    call shape%export_Xg('vtk/logo_Xg_6.vtk')
    call shape%finalize()

    !> C
    call shape%set_half_ring([0.0_rk, 0.0_rk, 0.0_rk], 5.0_rk, 7.0_rk, 1.0_rk)
    call shape%rotate_Xc(90.0_rk, -90.0_rk, 0.0_rk)
    call shape%translate_Xc([15.0_rk, 1.0_rk, 3.5_rk])
    call shape%export_Xc('vtk/logo_Xc_7.vtk')
    call shape%create(15*3, 3*3, 6*3)
    call shape%export_Xg('vtk/logo_Xg_7.vtk')
    call shape%finalize()


    !> A
    call shape%set_hexahedron(L=[1.0_rk, 1.0_rk, 7.0_rk], nc=[2,2,2])
    call shape%translate_Xc([16.0_rk, 0.0_rk, 0.0_rk])
    call shape%modify_Xc(17.0_rk,5,1)
    call shape%modify_Xc(18.0_rk,6,1)
    call shape%modify_Xc(17.0_rk,7,1)
    call shape%modify_Xc(18.0_rk,8,1)
    call shape%export_Xc('vtk/logo_Xc_8.vtk')
    call shape%create(3*3, 3*3, 21*3)
    call shape%export_Xg('vtk/logo_Xg_8.vtk')
    call shape%finalize()

    call shape%set_hexahedron(L=[1.0_rk, 1.0_rk, 7.0_rk], nc=[2,2,2])
    call shape%translate_Xc([17.0_rk, 0.0_rk, 0.0_rk])
    call shape%modify_Xc(19.0_rk,1,1)
    call shape%modify_Xc(20.0_rk,2,1)
    call shape%modify_Xc(19.0_rk,3,1)
    call shape%modify_Xc(20.0_rk,4,1)
    call shape%export_Xc('vtk/logo_Xc_9.vtk')
    call shape%create(3*3, 3*3, 21*3)
    call shape%export_Xg('vtk/logo_Xg_9.vtk')
    call shape%finalize()

    call shape%set_hexahedron(L=[1.5_rk, 1.0_rk, 1.0_rk], nc=[2,2,2])
    call shape%translate_Xc([17.0_rk, 0.0_rk, 3.0_rk])
    call shape%export_Xc('vtk/logo_Xc_10.vtk')
    call shape%create(3*3, 3*3, 3*3)
    call shape%export_Xg('vtk/logo_Xg_10.vtk')
    call shape%finalize()

    !> D
    call shape%set_hexahedron(L=[1.0_rk, 1.0_rk, 7.0_rk], nc=[2,2,2])
    call shape%translate_Xc([21.0_rk, 0.0_rk, 0.0_rk])
    call shape%export_Xc('vtk/logo_Xc_11.vtk')
    call shape%create(3*3, 3*3, 21*3)
    call shape%export_Xg('vtk/logo_Xg_11.vtk')
    call shape%finalize()

    call shape%set_half_ring([0.0_rk, 0.0_rk, 0.0_rk], 5.0_rk, 7.0_rk, 1.0_rk)
    call shape%rotate_Xc(90.0_rk, 90.0_rk, 0.0_rk)
    call shape%translate_Xc([22.0_rk, 1.0_rk, 3.5_rk])
    call shape%export_Xc('vtk/logo_Xc_12.vtk')
    call shape%create(15*3, 3*3, 6*3)
    call shape%export_Xg('vtk/logo_Xg_12.vtk')
    call shape%finalize()

end program
