!> Visualization of NURBS surfaces using ForCAD, ForImage and ForColorMap libraries
!> This example converts NURBS surfaces (vector) to an image (raster) using ForCAD, ForImage and ForColorMap libraries.

program example_ppm3

    use forcad, only: rk, nurbs_surface
    use forimage, only: ik, format_pnm, color
    use forcolormap, only: colormap, wp
    use fortime, only: timer

    implicit none
    type(nurbs_surface)      :: shape
    type(format_pnm)         :: image
    type(color)              :: background_color
    type(colormap)           :: cmap
    integer(ik), allocatable :: px(:, :)
    real(rk), allocatable    :: Xg(:,:), z_values(:)
    real(rk)                 :: aspect_ratio
    integer                  :: height, width, ng(2), red, green, blue, res1, res2, i
    integer, allocatable     :: idx(:,:)
    type(timer)              :: t

    !-----------------------------------------------------------------------------
    ! Set the image size and calculate the aspect ratio
    !-----------------------------------------------------------------------------
    width  = 2000
    height = 2000

    aspect_ratio = real(width,rk) / real(height,rk)
    allocate(px(height, 3*width))

    !-----------------------------------------------------------------------------
    ! Set the background color using ForColor class of ForImage
    !-----------------------------------------------------------------------------
    call t%timer_start()
    call background_color%set('white', use_library=.true.)

    do i = 1, width
        px(:, 3*(i-1)+1) = background_color%get_r()
        px(:, 3*(i-1)+2) = background_color%get_g()
        px(:, 3*(i-1)+3) = background_color%get_b()
    end do
    call t%timer_stop(message='Setting the background color')








    !-----------------------------------------------------------------------------
    ! Set a shape object using ForCAD library
    !-----------------------------------------------------------------------------
    call t%timer_start()
    !> Set the shape parameters for a tetragon
    res1         = 2000
    res2         = 2000

    call shape%set_tetragon(L=[1.0_rk, 1.0_rk], nc=[2,2])
    call shape%create(res1, res2)
    Xg = shape%get_Xg()
    ng = shape%get_ng()
    call shape%finalize()
    call t%timer_stop(message='Creating a tetragon')

    !-----------------------------------------------------------------------------
    ! Set the z-values of the geometry points
    !-----------------------------------------------------------------------------
    ! circular gradient
    z_values = ((Xg(:,1) - minval(Xg(:,1))) / (maxval(Xg(:,1)) - minval(Xg(:,1))))**2 &
        + ((Xg(:,2) - minval(Xg(:,2))) / (maxval(Xg(:,2)) - minval(Xg(:,2)))**2)

    !-----------------------------------------------------------------------------
    ! Set the colormap using ForColorMap library
    !-----------------------------------------------------------------------------
    call cmap%set('buda',real(0.0_rk, kind=wp), real(1.0_rk, kind=wp))

    !-----------------------------------------------------------------------------
    ! Set colors to the shape
    !-----------------------------------------------------------------------------
    call t%timer_start()
    Xg(:,2) = Xg(:,2) * aspect_ratio
    allocate(idx(ng(1)*ng(2),2))
    idx(:,1) = min(max(1, int(Xg(:,1) * width ) + 1), width )
    idx(:,2) = min(max(1, int(Xg(:,2) * height) + 1), height)
    do i = 1, ng(1)*ng(2)
        call cmap%compute_RGB(real(z_values(i), kind=wp), red, green, blue)
        px(idx(i,2), 3*(idx(i,1)-1)+1) = red
        px(idx(i,2), 3*(idx(i,1)-1)+2) = green
        px(idx(i,2), 3*(idx(i,1)-1)+3) = blue
    end do
    deallocate(idx)
    call t%timer_stop(message='Setting colors')









    !-----------------------------------------------------------------------------
    ! Set a shape object using ForCAD library
    !-----------------------------------------------------------------------------
    call t%timer_start()
    res1         = 2000
    res2         = 2000

    call shape%set_tetragon(L=[0.48_rk, 0.48_rk], nc=[2,2])
    call shape%translate_Xc([0.01_rk, 0.01_rk, 0.0_rk])
    call shape%create(res1, res2)
    Xg = shape%get_Xg()
    ng = shape%get_ng()
    call shape%finalize()
    call t%timer_stop(message='Creating a ring')

    !-----------------------------------------------------------------------------
    ! Set the z-values of the geometry points
    !-----------------------------------------------------------------------------
    ! circular gradient
    z_values = ((Xg(:,1) - minval(Xg(:,1))) / (maxval(Xg(:,1)) - minval(Xg(:,1))))**2 &
        + ((Xg(:,2) - minval(Xg(:,2))) / (maxval(Xg(:,2)) - minval(Xg(:,2)))**2)

    !-----------------------------------------------------------------------------
    ! Set the colormap using ForColorMap library
    !-----------------------------------------------------------------------------
    call cmap%set('managua',real(0.0_rk, kind=wp), real(2.2_rk, kind=wp))

    !-----------------------------------------------------------------------------
    ! Set colors to the shape
    !-----------------------------------------------------------------------------
    call t%timer_start()
    Xg(:,2) = Xg(:,2) * aspect_ratio
    allocate(idx(ng(1)*ng(2),2))
    idx(:,1) = min(max(1, int(Xg(:,1) * width ) + 1), width )
    idx(:,2) = min(max(1, int(Xg(:,2) * height) + 1), height)
    do i = 1, ng(1)*ng(2)
        call cmap%compute_RGB(real(z_values(i), kind=wp), red, green, blue)
        px(idx(i,2), 3*(idx(i,1)-1)+1) = red
        px(idx(i,2), 3*(idx(i,1)-1)+2) = green
        px(idx(i,2), 3*(idx(i,1)-1)+3) = blue
    end do
    deallocate(idx)
    call t%timer_stop(message='Setting colors')















    !-----------------------------------------------------------------------------
    ! Set a shape object using ForCAD library
    !-----------------------------------------------------------------------------
    call t%timer_start()
    res1         = 2000
    res2         = 2000

    call shape%set_tetragon(L=[0.48_rk, 0.48_rk], nc=[3,2])
    call shape%translate_Xc([0.51_rk, 0.01_rk, 0.0_rk])
    call shape%modify_Xc(0.24_rk,2,2)
    call shape%modify_Xc(0.26_rk,5,2)
    call shape%create(res1, res2)
    Xg = shape%get_Xg()
    ng = shape%get_ng()
    call shape%finalize()
    call t%timer_stop(message='Creating a ring')

    !-----------------------------------------------------------------------------
    ! Set the z-values of the geometry points
    !-----------------------------------------------------------------------------
    ! linear gradient in the y-direction
    z_values = (Xg(:,2) - minval(Xg(:,2))) / (maxval(Xg(:,2)) - minval(Xg(:,2)))

    !-----------------------------------------------------------------------------
    ! Set the colormap using ForColorMap library
    !-----------------------------------------------------------------------------
    call cmap%set('lipari',real(0.0_rk, kind=wp), real(1.0_rk, kind=wp))

    !-----------------------------------------------------------------------------
    ! Set colors to the shape
    !-----------------------------------------------------------------------------
    call t%timer_start()
    Xg(:,2) = Xg(:,2) * aspect_ratio
    allocate(idx(ng(1)*ng(2),2))
    idx(:,1) = min(max(1, int(Xg(:,1) * width ) + 1), width )
    idx(:,2) = min(max(1, int(Xg(:,2) * height) + 1), height)
    do i = 1, ng(1)*ng(2)
        call cmap%compute_RGB(real(z_values(i), kind=wp), red, green, blue)
        px(idx(i,2), 3*(idx(i,1)-1)+1) = red
        px(idx(i,2), 3*(idx(i,1)-1)+2) = green
        px(idx(i,2), 3*(idx(i,1)-1)+3) = blue
    end do
    deallocate(idx)
    call t%timer_stop(message='Setting colors')














    !-----------------------------------------------------------------------------
    ! Set a shape object using ForCAD library
    !-----------------------------------------------------------------------------
    call t%timer_start()
    res1         = 2000
    res2         = 2000

    call shape%set_tetragon(L=[0.48_rk, 0.48_rk], nc=[2,3])
    call shape%translate_Xc([0.01_rk, 0.51_rk, 0.0_rk])
    call shape%modify_Xc(0.26_rk,3,1)
    call shape%modify_Xc(0.24_rk,4,1)
    call shape%create(res1, res2)
    Xg = shape%get_Xg()
    ng = shape%get_ng()
    call shape%finalize()
    call t%timer_stop(message='Creating a ring')

    !-----------------------------------------------------------------------------
    ! Set the z-values of the geometry points
    !-----------------------------------------------------------------------------
    ! linear gradient in the x-direction
    z_values = (Xg(:,1) - minval(Xg(:,1))) / (maxval(Xg(:,1)) - minval(Xg(:,1)))

    !-----------------------------------------------------------------------------
    ! Set the colormap using ForColorMap library
    !-----------------------------------------------------------------------------
    call cmap%set('oslo10',real(0.0_rk, kind=wp), real(1.0_rk, kind=wp))

    !-----------------------------------------------------------------------------
    ! Set colors to the shape
    !-----------------------------------------------------------------------------
    call t%timer_start()
    Xg(:,2) = Xg(:,2) * aspect_ratio
    allocate(idx(ng(1)*ng(2),2))
    idx(:,1) = min(max(1, int(Xg(:,1) * width ) + 1), width )
    idx(:,2) = min(max(1, int(Xg(:,2) * height) + 1), height)
    do i = 1, ng(1)*ng(2)
        call cmap%compute_RGB(real(z_values(i), kind=wp), red, green, blue)
        px(idx(i,2), 3*(idx(i,1)-1)+1) = red
        px(idx(i,2), 3*(idx(i,1)-1)+2) = green
        px(idx(i,2), 3*(idx(i,1)-1)+3) = blue
    end do
    deallocate(idx)
    call t%timer_stop(message='Setting colors')















    !-----------------------------------------------------------------------------
    ! Set a shape object using ForCAD library
    !-----------------------------------------------------------------------------
    call t%timer_start()
    res1         = 2000
    res2         = 2000

    call shape%set_tetragon(L=[0.48_rk, 0.48_rk], nc=[3,3])
    call shape%translate_Xc([0.51_rk, 0.51_rk, 0.0_rk])
    call shape%modify_Xc(0.7_rk,1,2)
    call shape%modify_Xc(0.7_rk,3,2)
    call shape%modify_Xc(0.8_rk,7,2)
    call shape%modify_Xc(0.8_rk,9,2)
    call shape%create(res1, res2)
    Xg = shape%get_Xg()
    ng = shape%get_ng()
    call shape%finalize()
    call t%timer_stop(message='Creating a ring')

    !-----------------------------------------------------------------------------
    ! Set the z-values of the geometry points
    !-----------------------------------------------------------------------------
    ! linear gradient in the x-direction
    z_values = (Xg(:,1) - minval(Xg(:,1))) / (maxval(Xg(:,1)) - minval(Xg(:,1)))

    !-----------------------------------------------------------------------------
    ! Set the colormap using ForColorMap library
    !-----------------------------------------------------------------------------
    red   = 255
    green = 215
    blue  = 0

    !-----------------------------------------------------------------------------
    ! Set colors to the shape
    !-----------------------------------------------------------------------------
    call t%timer_start()
    Xg(:,2) = Xg(:,2) * aspect_ratio
    allocate(idx(ng(1)*ng(2),2))
    idx(:,1) = min(max(1, int(Xg(:,1) * width ) + 1), width )
    idx(:,2) = min(max(1, int(Xg(:,2) * height) + 1), height)
    do i = 1, ng(1)*ng(2)
        px(idx(i,2), 3*(idx(i,1)-1)+1) = red
        px(idx(i,2), 3*(idx(i,1)-1)+2) = green
        px(idx(i,2), 3*(idx(i,1)-1)+3) = blue
    end do
    deallocate(idx)
    call t%timer_stop(message='Setting colors')






    !-----------------------------------------------------------------------------
    ! Save the image to a PPM file using ForImage library
    !-----------------------------------------------------------------------------
    call t%timer_start()
    call image%set_pnm(&
        encoding    = 'binary', &
        file_format = 'ppm', &
        width       = width, &
        height      = height, &
        max_color   = 255, &
        comment     = 'example: ForCAD + ForImage + ForColor + ForColormap', &
        pixels      = px &
        )
    call image%export_pnm('ppm/example_ppm3')
    call image%finalize()
    call t%timer_stop(message='Saving the image')







    ! Clean up
    call cmap%finalize()
    deallocate(px, Xg, z_values)

end program
