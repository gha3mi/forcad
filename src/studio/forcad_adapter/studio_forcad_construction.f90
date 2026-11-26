module studio_forcad_construction
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use forcad, only: rk, nurbs_curve, nurbs_surface, nurbs_volume, extrude, revolve, sweep, loft
    use studio_forcad_transforms, only: studio_native_patch, native_patch_kind
    implicit none(type, external)
    private

    public :: construction_names, construction_actions, stage_studio_construction
    public :: extrude_operation, revolve_operation, sweep_operation, loft_operation

    character(len=7), parameter :: construction_names(4)=[character(len=7)::'Extrude','Revolve','Sweep','Loft']
    character(len=7), parameter :: construction_actions(4)=[character(len=7)::'extrude','revolve','sweep','loft']
    integer, parameter :: extrude_operation=1,revolve_operation=2,sweep_operation=3,loft_operation=4
    real(rk), parameter :: radians_per_degree=acos(-1._rk)/180._rk
contains
    subroutine stage_studio_construction(operation,profiles,object,error,vector,axis_point,angle_degrees,spine,degree)
        integer, intent(in) :: operation
        type(studio_native_patch), intent(in) :: profiles(:)
        type(studio_native_patch), intent(out) :: object
        character(:), allocatable, intent(out) :: error
        real(rk), intent(in), optional :: vector(3),axis_point(3),angle_degrees
        type(studio_native_patch), intent(in), optional :: spine
        integer, intent(in), optional :: degree
        type(nurbs_curve), allocatable :: curves(:)
        type(nurbs_surface), allocatable :: surfaces(:)
        type(nurbs_curve) :: path
        type(nurbs_surface) :: surface
        type(nurbs_volume) :: volume
        real(rk) :: origin(3),direction(3),angle
        real(rk), allocatable :: path_start(:,:)
        integer :: i,dimension,kind

        error='Choose Extrude, Revolve, Sweep, or Loft.'
        if(operation<extrude_operation.or.operation>loft_operation)return
        error='Choose a curve or surface profile.'
        if(size(profiles)==0)return
        if(.not.allocated(profiles(1)%geometry))return
        kind=native_patch_kind(profiles(1)%geometry)
        if(kind<1.or.kind>2)return
        if(operation==loft_operation)then
            error='Choose at least two sections for Loft.'
            if(size(profiles)<2)return
            if(present(degree))then
                error='Choose a loft degree lower than the section count.'
                if(degree<1.or.degree>=size(profiles))return
            end if
        else
            error='Choose one profile for '//trim(construction_names(operation))//'.'
            if(size(profiles)/=1)return
        end if

        origin=0._rk
        direction=0._rk
        angle=0._rk
        select case(operation)
        case(extrude_operation,revolve_operation)
            error='Enter a finite, nonzero '//merge('displacement.','axis vector. ',operation==extrude_operation)
            if(.not.present(vector))return
            if(.not.all(ieee_is_finite(vector)).or.all(vector==0._rk))return
            direction=vector
            if(operation==revolve_operation)then
                error='Enter a finite axis point and a signed angle up to 360 degrees.'
                if(.not.present(axis_point).or..not.present(angle_degrees))return
                if(.not.all(ieee_is_finite(axis_point)).or..not.ieee_is_finite(angle_degrees))return
                if(angle_degrees==0._rk.or.abs(angle_degrees)>360._rk)return
                origin=axis_point
                angle=angle_degrees*radians_per_degree
                if(angle==0._rk)return
            end if
        case(sweep_operation)
            error='Choose a valid curve as the sweep path.'
            if(.not.present(spine))return
            if(.not.load_curve(spine,path))return
            path_start=path%cmp_Xg(Xt=[path%get_knot(path%get_degree()+1)])
            if(size(path_start,1)/=1.or.size(path_start,2)>3)return
            if(.not.all(ieee_is_finite(path_start)))return
            origin(:size(path_start,2))=path_start(1,:)
        case default
        end select

        error='All profiles must be available patches of the same dimension.'
        select case(kind)
        case(1)
            allocate(curves(size(profiles)))
            do i=1,size(profiles)
                if(.not.load_curve(profiles(i),curves(i)))return
            end do
            select case(operation)
            case(extrude_operation)
                surface=extrude(curves(1),direction)
            case(revolve_operation)
                surface=revolve(curves(1),origin,direction,angle)
            case(sweep_operation)
                dimension=max(size(curves(1)%get_Xc(),2),size(path_start,2))
                if(dimension>3)return
                surface=sweep(curves(1),path,origin=origin(:dimension))
            case(loft_operation)
                surface=loft(curves,degree=degree)
            case default
                return
            end select
            if(surface%is_valid())allocate(object%geometry,source=surface)
        case(2)
            allocate(surfaces(size(profiles)))
            do i=1,size(profiles)
                if(.not.allocated(profiles(i)%geometry))return
                select type(profile=>profiles(i)%geometry)
                type is(nurbs_surface)
                    if(.not.profile%is_valid())return
                    surfaces(i)=profile
                class default
                    return
                end select
            end do
            select case(operation)
            case(extrude_operation)
                volume=extrude(surfaces(1),direction)
            case(revolve_operation)
                volume=revolve(surfaces(1),origin,direction,angle)
            case(sweep_operation)
                dimension=max(size(surfaces(1)%get_Xc(),2),size(path_start,2))
                if(dimension>3)return
                volume=sweep(surfaces(1),path,origin=origin(:dimension))
            case(loft_operation)
                volume=loft(surfaces,degree=degree)
            case default
                return
            end select
            if(volume%is_valid())allocate(object%geometry,source=volume)
        case default
            return
        end select

        error=''
        if(.not.allocated(object%geometry))then
            error='ForCAD could not construct the '//trim(construction_names(operation))//' result.'
            if(operation==loft_operation)then
                error='ForCAD rejected the loft. Check the section degrees, control counts, knot vectors, and order.'
            end if
        end if
    end subroutine stage_studio_construction

    logical function load_curve(source,curve) result(valid)
        type(studio_native_patch), intent(in) :: source
        type(nurbs_curve), intent(out) :: curve

        valid=.false.
        if(.not.allocated(source%geometry))return
        select type(geometry=>source%geometry)
        type is(nurbs_curve)
            if(.not.geometry%is_valid())return
            curve=geometry
            valid=.true.
        end select
    end function load_curve
end module studio_forcad_construction
