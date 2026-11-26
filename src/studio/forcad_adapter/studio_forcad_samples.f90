module studio_forcad_samples
    use forcad, only: rk, nurbs_curve, nurbs_surface, nurbs_volume
    implicit none(type, external)
    private

    public :: make_sample_patch

contains

    subroutine make_sample_patch(kind,variant,geometry)
        integer, intent(in) :: kind,variant
        class(*), allocatable, intent(out) :: geometry
        real(rk), parameter :: curve_controls(3,3)=reshape([ &
            0._rk,.5_rk,1._rk, 0._rk,1._rk,0._rk, 0._rk,0._rk,0._rk],[3,3])
        type(nurbs_curve) :: curve
        type(nurbs_surface) :: surface
        type(nurbs_volume) :: volume

        select case(kind)
        case(1)
            if(variant==2)then
                curve=make_showcase_curve()
            else
                call curve%set(Xc=curve_controls)
            end if
            allocate(geometry,source=curve)
        case(2)
            if(variant==2)then
                surface=make_showcase_surface()
            else
                call surface%set_tetragon(L=[1._rk,1._rk],nc=[2,2])
            end if
            allocate(geometry,source=surface)
        case(3)
            if(variant==2)then
                volume=make_showcase_volume()
            else
                call volume%set_hexahedron(L=[1._rk,1._rk,1._rk],nc=[2,2,2])
            end if
            allocate(geometry,source=volume)
        case default
            return
        end select
    end subroutine make_sample_patch

    function make_showcase_curve() result(curve)
        type(nurbs_curve) :: curve
        real(rk), parameter :: knots(12)=[0._rk,0._rk,0._rk,0._rk,.19_rk,.41_rk,.62_rk,.81_rk,1._rk,1._rk,1._rk,1._rk]
        ! The centerline rises into two banks and drops through the middle turn.
        real(rk), parameter :: controls(8,3)=reshape([ &
            -1.65_rk,-1.20_rk,-.72_rk,-.15_rk,.32_rk,.88_rk,1.30_rk,1.74_rk, &
             -.48_rk, -.40_rk, .88_rk,1.18_rk,-.72_rk,-1.18_rk,.26_rk,.46_rk, &
              .05_rk,  .41_rk, .46_rk, .13_rk,-.19_rk,  .13_rk,.55_rk,.59_rk],[8,3])

        call curve%set(knot=knots,Xc=controls,degree=3)
        call curve%rotate_Xc(alpha=30._rk,beta=0._rk,theta=-37._rk)
        call curve%translate_Xc(vec=[-3.3_rk,1.5_rk,.25_rk])
    end function make_showcase_curve

    function make_showcase_surface() result(surface)
        type(nurbs_surface) :: surface
        real(rk), parameter :: knots_u(9)=[0._rk,0._rk,0._rk,0._rk,.54_rk,1._rk,1._rk,1._rk,1._rk]
        real(rk), parameter :: knots_v(8)=[0._rk,0._rk,0._rk,0._rk,1._rk,1._rk,1._rk,1._rk]
        ! Crosswise rows flare at their ends. U varies fastest in the control net.
        real(rk), parameter :: controls(20,3)=reshape([ &
            -1.40_rk,-.93_rk,-.16_rk,.85_rk,1.48_rk, &
            -1.55_rk,-1.03_rk,-.08_rk,.90_rk,1.55_rk, &
            -1.58_rk,-1.08_rk,-.11_rk,.84_rk,1.56_rk, &
            -1.48_rk,-1.01_rk,-.18_rk,.77_rk,1.42_rk, &
            -1.15_rk,-1.20_rk,-1.11_rk,-1.02_rk,-.98_rk, &
             -.51_rk,-.50_rk,-.42_rk,-.34_rk,-.30_rk, &
              .40_rk,.45_rk,.47_rk,.46_rk,.49_rk, &
             1.05_rk,1.16_rk,1.22_rk,1.18_rk,1.08_rk, &
              .68_rk,.61_rk,.04_rk,-.32_rk,-.75_rk, &
              .42_rk,.24_rk,-.28_rk,-.08_rk,-.29_rk, &
             -.25_rk,-.18_rk,-.36_rk,.26_rk,.49_rk, &
             -.61_rk,-.35_rk,.06_rk,.81_rk,.90_rk],[20,3])

        call surface%set(knot1=knots_u,knot2=knots_v,Xc=controls,degree=[3,3])
        call surface%rotate_Xc(alpha=40._rk,beta=0._rk,theta=0._rk)
        call surface%translate_Xc(vec=[-.9_rk,-.6_rk,.18_rk])
    end function make_showcase_surface

    function make_showcase_volume() result(volume)
        type(nurbs_volume) :: volume

        ! ForCAD builds the rational angular section, radial thickness, and axial depth.
        call volume%set_half_ring(center=[0._rk,0._rk,-.57_rk],radius1=2.08_rk,radius2=3.25_rk,length=1.14_rk)
        call volume%rotate_Xc(alpha=8._rk,beta=-12._rk,theta=-18._rk)
        call volume%translate_Xc(vec=[2.3_rk,-2.7_rk,0._rk])
    end function make_showcase_volume


end module studio_forcad_samples
