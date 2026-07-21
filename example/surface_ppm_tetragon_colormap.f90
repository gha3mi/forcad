!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Render NURBS tetragon surfaces to a binary PPM image using ForCAD, ForImage and ForColorMap.
program surface_ppm_tetragon_colormap

    use forcad, only: rk, nurbs_surface
    use forimage, only: ik, format_pnm, color
    use forcolormap, only: colormap, wp
    use fortime, only: timer

    implicit none

    integer, parameter :: width = 2000
    integer, parameter :: height = 2000
    integer, parameter :: gradient_circular = 1
    integer, parameter :: gradient_x = 2
    integer, parameter :: gradient_y = 3
    real(rk), parameter :: width_r = real(width, rk)
    real(rk), parameter :: height_r = real(height, rk)
    real(rk), parameter :: aspect_ratio = width_r/height_r
    real(rk), parameter :: height_aspect = aspect_ratio*height_r

    type(nurbs_surface)      :: shape
    type(format_pnm)         :: image
    type(color)              :: background_color
    type(colormap)           :: cmap
    type(timer)              :: t
    integer(ik), allocatable :: px(:, :)
    real(rk), allocatable    :: Xg(:,:)

    allocate(px(height, 3*width))

    call t%timer_start()
    call background_color%set("white", use_library=.true.)
    call fill_background(background_color%get_r(), background_color%get_g(), background_color%get_b())
    call t%timer_stop(message="Setting the background color")

    call t%timer_start()
    call shape%set_tetragon(L=[1.0_rk, 1.0_rk], nc=[2,2])
    call shape%create(2000, 2000)
    Xg = shape%get_Xg()
    call shape%finalize()
    call t%timer_stop(message="Creating a tetragon")

    call cmap%set("buda", real(0.0_rk, kind=wp), real(1.0_rk, kind=wp))
    call t%timer_start()
    call paint_colormap(gradient_circular)
    call t%timer_stop(message="Setting colors")
    deallocate(Xg)

    call t%timer_start()
    call shape%set_tetragon(L=[0.48_rk, 0.48_rk], nc=[2,2])
    call shape%translate_Xc([0.01_rk, 0.01_rk/aspect_ratio, 0.0_rk])
    call shape%create(2000, 2000)
    Xg = shape%get_Xg()
    call shape%finalize()
    call t%timer_stop(message="Creating a tetragon")

    call cmap%set("managua", real(0.0_rk, kind=wp), real(2.2_rk, kind=wp))
    call t%timer_start()
    call paint_colormap(gradient_circular)
    call t%timer_stop(message="Setting colors")
    deallocate(Xg)

    call t%timer_start()
    call shape%set_tetragon(L=[0.48_rk, 0.48_rk], nc=[3,2])
    call shape%translate_Xc([0.51_rk, 0.01_rk/aspect_ratio, 0.0_rk])
    call shape%modify_Xc(0.24_rk, 2, 2)
    call shape%modify_Xc(0.26_rk, 5, 2)
    call shape%create(2000, 2000)
    Xg = shape%get_Xg()
    call shape%finalize()
    call t%timer_stop(message="Creating a tetragon")

    call cmap%set("lipari", real(0.0_rk, kind=wp), real(1.0_rk, kind=wp))
    call t%timer_start()
    call paint_colormap(gradient_y)
    call t%timer_stop(message="Setting colors")
    deallocate(Xg)

    call t%timer_start()
    call shape%set_tetragon(L=[0.48_rk, 0.48_rk], nc=[2,3])
    call shape%translate_Xc([0.01_rk, 0.51_rk/aspect_ratio, 0.0_rk])
    call shape%modify_Xc(0.26_rk, 3, 1)
    call shape%modify_Xc(0.24_rk, 4, 1)
    call shape%create(2000, 2000)
    Xg = shape%get_Xg()
    call shape%finalize()
    call t%timer_stop(message="Creating a tetragon")

    call cmap%set("oslo10", real(0.0_rk, kind=wp), real(1.0_rk, kind=wp))
    call t%timer_start()
    call paint_colormap(gradient_x)
    call t%timer_stop(message="Setting colors")
    deallocate(Xg)

    call t%timer_start()
    call shape%set_tetragon(L=[0.48_rk, 0.48_rk], nc=[3,3])
    call shape%translate_Xc([0.51_rk, 0.51_rk/aspect_ratio, 0.0_rk])
    call shape%modify_Xc(0.7_rk, 1, 2)
    call shape%modify_Xc(0.7_rk, 3, 2)
    call shape%modify_Xc(0.8_rk, 7, 2)
    call shape%modify_Xc(0.8_rk, 9, 2)
    call shape%create(2000, 2000)
    Xg = shape%get_Xg()
    call shape%finalize()
    call t%timer_stop(message="Creating a tetragon")

    call t%timer_start()
    call paint_solid(red=255, green=215, blue=0)
    call t%timer_stop(message="Setting colors")
    deallocate(Xg)

    call t%timer_start()
    call image%set_pnm(&
        encoding    = "binary",&
        file_format = "ppm",&
        width       = width,&
        height      = height,&
        max_color   = 255,&
        comment     = "example: ForCAD + ForImage + ForColorMap",&
        pixels      = px)
    call image%export_pnm("ppm/ppm_tetragon_colormap")
    call image%finalize()
    call t%timer_stop(message="Saving the image")

    call cmap%finalize()
    deallocate(px)

contains

    !===============================================================================
    subroutine fill_background(red, green, blue)
        integer, intent(in) :: red, green, blue
        integer :: ix, iy, offset

        do concurrent (ix = 1:width, iy = 1:height) local(offset) shared(px, red, green, blue) default(none)
            offset = 3*ix - 2
            px(iy, offset  ) = red
            px(iy, offset+1) = green
            px(iy, offset+2) = blue
        end do
    end subroutine fill_background
    !===============================================================================


    !===============================================================================
    subroutine paint_colormap(gradient)
        integer, intent(in) :: gradient
        real(rk) :: xmin, xmax, ymin, ymax, xden, yden, x, y, z
        integer :: i, ix, iy, ix0, ix1, iy0, iy1, offset, red, green, blue

        xmin = Xg(1,1)
        xmax = xmin
        ymin = Xg(1,2)
        ymax = ymin
        do i = 2, size(Xg, 1)
            xmin = min(xmin, Xg(i,1))
            xmax = max(xmax, Xg(i,1))
            ymin = min(ymin, Xg(i,2))
            ymax = max(ymax, Xg(i,2))
        end do
        xden = xmax - xmin
        yden = ymax - ymin
        ix0 = min(max(1, int(xmin*width_r) + 1), width)
        ix1 = min(max(1, int(xmax*width_r) + 1), width)
        iy0 = min(max(1, int(ymin*height_aspect) + 1), height)
        iy1 = min(max(1, int(ymax*height_aspect) + 1), height)

        do i = 1, size(Xg, 1)
            ix = min(max(1, int(Xg(i,1)*width_r) + 1), width)
            iy = min(max(1, int(Xg(i,2)*height_aspect) + 1), height)
            if (ix < ix0 .or. ix > ix1 .or. iy < iy0 .or. iy > iy1) cycle

            x = (Xg(i,1) - xmin)/xden
            y = (Xg(i,2) - ymin)/yden
            select case (gradient)
            case (gradient_x)
                z = x
            case (gradient_y)
                z = y
            case default
                z = x*x + (Xg(i,2) - ymin)/(yden**2)
            end select

            call cmap%compute_RGB(real(z, kind=wp), red, green, blue)

            offset = 3*ix - 2
            px(iy, offset  ) = red
            px(iy, offset+1) = green
            px(iy, offset+2) = blue
        end do
    end subroutine paint_colormap
    !===============================================================================


    !===============================================================================
    subroutine paint_solid(red, green, blue)
        integer, intent(in) :: red, green, blue
        real(rk) :: xmin, xmax, ymin, ymax
        integer :: i, ix, iy, ix0, ix1, iy0, iy1, offset

        xmin = Xg(1,1)
        xmax = xmin
        ymin = Xg(1,2)
        ymax = ymin
        do i = 2, size(Xg, 1)
            xmin = min(xmin, Xg(i,1))
            xmax = max(xmax, Xg(i,1))
            ymin = min(ymin, Xg(i,2))
            ymax = max(ymax, Xg(i,2))
        end do
        ix0 = min(max(1, int(xmin*width_r) + 1), width)
        ix1 = min(max(1, int(xmax*width_r) + 1), width)
        iy0 = min(max(1, int(ymin*height_aspect) + 1), height)
        iy1 = min(max(1, int(ymax*height_aspect) + 1), height)

        do i = 1, size(Xg, 1)
            ix = min(max(1, int(Xg(i,1)*width_r) + 1), width)
            iy = min(max(1, int(Xg(i,2)*height_aspect) + 1), height)
            if (ix < ix0 .or. ix > ix1 .or. iy < iy0 .or. iy > iy1) cycle

            offset = 3*ix - 2
            px(iy, offset  ) = red
            px(iy, offset+1) = green
            px(iy, offset+2) = blue
        end do
    end subroutine paint_solid
    !===============================================================================
end program surface_ppm_tetragon_colormap
