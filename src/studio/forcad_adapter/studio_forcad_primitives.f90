module studio_forcad_primitives
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use forcad, only: rk, nurbs_curve, nurbs_surface, nurbs_volume, revolve, extrude, rotation
    use studio_forcad_transforms, only: studio_native_patch
    implicit none(type, external)
    private

    public :: stage_studio_primitive, primitive_names, primitive_actions
    public :: point_shape, line_shape, circle_shape, arc_shape
    public :: rectangle_shape, box_shape, cylinder_shape, tube_shape

    integer, parameter :: point_shape=1,line_shape=2,circle_shape=3,arc_shape=4
    integer, parameter :: rectangle_shape=5,box_shape=6,cylinder_shape=7,tube_shape=8
    character(len=20), parameter :: primitive_names(8)=[character(len=20):: &
        'Point','Line','Circle','Arc','Rectangle','Box','Cylinder surface','Tube volume']
    character(len=9), parameter :: primitive_actions(8)=[character(len=9):: &
        'point','line','circle','arc','rectangle','box','cylinder','tube']
    real(rk), parameter :: identity(3,3)=reshape([1,0,0,0,1,0,0,0,1],[3,3])
    real(rk), parameter :: origin_at_zero(3)=0._rk,positive_z(3)=[0._rk,0._rk,1._rk]
    real(rk), parameter :: radians_per_degree=acos(-1._rk)/180._rk

contains

    subroutine stage_studio_primitive(primitive,origin,values,object,error)
        integer, intent(in) :: primitive
        real(rk), intent(in) :: origin(3),values(3)
        type(studio_native_patch), intent(out) :: object
        character(:), allocatable, intent(out) :: error
        type(nurbs_curve) :: curve,carrier
        type(nurbs_surface) :: surface,revolved
        type(nurbs_volume) :: volume
        real(rk) :: radial(3)
        logical :: valid

        error='Choose a supported primitive type.'
        if(primitive<1.or.primitive>size(primitive_names))return
        error='Enter finite coordinates and dimensions.'
        if(.not.all(ieee_is_finite(origin)).or..not.all(ieee_is_finite(values)))return

        select case(primitive)
        case(point_shape)
            call curve%set(knot=[0._rk,1._rk],Xc=reshape(origin,[1,3]),degree=0)
            allocate(object%geometry,source=curve)
        case(line_shape)
            error='The endpoints must be different.'
            if(all(origin==values))return
            call set_segment(curve,origin,values)
            allocate(object%geometry,source=curve)
        case(circle_shape)
            error='The radius must be greater than zero.'
            if(values(1)<=0._rk)return
            call curve%set_circle(origin,values(1))
            allocate(object%geometry,source=curve)
        case(arc_shape)
            error='Use a positive radius and a nonzero sweep from −360° to 360°.'
            if(values(1)<=0._rk.or.values(3)==0._rk.or.abs(values(3))>360._rk)return
            ! Reuse the exact ForCAD revolution, then extract the circular
            ! boundary of the cylindrical patch at the start of its carrier.
            radial=[values(1),0._rk,0._rk]
            call set_segment(carrier,radial,radial+positive_z)
            revolved=revolve(carrier,origin_at_zero,positive_z,values(3)*radians_per_degree)
            curve=revolved%slice(1,0._rk)
            call curve%transform(rotation(0._rk,0._rk,modulo(values(2),360._rk)),origin)
            allocate(object%geometry,source=curve)
        case(rectangle_shape)
            error='Width and height must be greater than zero.'
            if(any(values(:2)<=0._rk))return
            call surface%set_tetragon(values(:2),[2,2])
            call surface%transform(identity,origin)
            allocate(object%geometry,source=surface)
        case(box_shape)
            error='Width, depth and height must be greater than zero.'
            if(any(values<=0._rk))return
            call volume%set_hexahedron(values,[2,2,2])
            call volume%transform(identity,origin)
            allocate(object%geometry,source=volume)
        case(cylinder_shape)
            error='Radius and height must be greater than zero.'
            if(any(values(:2)<=0._rk))return
            call curve%set_circle(origin,values(1))
            surface=extrude(curve,values(2)*positive_z)
            allocate(object%geometry,source=surface)
        case(tube_shape)
            error='Use a positive inner radius, a larger outer radius and a positive height.'
            if(values(1)<=0._rk.or.values(2)<=values(1).or.values(3)<=0._rk)return
            call volume%set_ring(origin,values(1),values(2),values(3))
            allocate(object%geometry,source=volume)
        case default
            return
        end select
        select type(geometry=>object%geometry)
        type is(nurbs_curve)
            valid=geometry%is_valid()
        type is(nurbs_surface)
            valid=geometry%is_valid()
        type is(nurbs_volume)
            valid=geometry%is_valid()
        class default
            valid=.false.
        end select
        error=''
        if(.not.valid)then
            deallocate(object%geometry)
            error='ForCAD could not create the requested '//trim(primitive_names(primitive))//'.'
        end if
    end subroutine stage_studio_primitive

    subroutine set_segment(curve,first,last)
        type(nurbs_curve), intent(inout) :: curve
        real(rk), intent(in) :: first(3),last(3)
        real(rk) :: controls(2,3)

        controls(1,:)=first
        controls(2,:)=last
        call curve%set(knot=[0._rk,0._rk,1._rk,1._rk],Xc=controls,degree=1)
    end subroutine set_segment
end module studio_forcad_primitives
