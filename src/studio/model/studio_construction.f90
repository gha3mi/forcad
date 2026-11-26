! Stage and publish patches produced by the native ForCAD construction operators.
module studio_construction
    use studio_forcad, only: rk, studio_native_patch, construction_names, construction_actions, stage_studio_construction, &
        extrude_operation, loft_operation, stage_studio_primitive
    use studio_single_patch, only: studio_patch_catalog, copy_single_patch_geometry, publish_native_patches
    implicit none(type, external)
    private
    public :: construction_names, publish_construction, stage_studio_construction, construction_actions
    public :: construct_studio_primitive

contains
    subroutine publish_construction(operation,kind,sections,name,created_kind,error,vector,axis_point,angle_degrees,spine,degree,catalog)
        integer, intent(in) :: operation,kind,sections(:)
        character(*), intent(in) :: name
        integer, intent(out) :: created_kind
        character(:), allocatable, intent(out) :: error
        real(rk), intent(in), optional :: vector(3),axis_point(3),angle_degrees
        integer, intent(in), optional :: spine,degree
        type(studio_patch_catalog), intent(inout), target, optional :: catalog
        type(studio_native_patch) :: objects(1), path
        type(studio_native_patch), allocatable :: profiles(:)
        integer :: i,variant

        created_kind=0
        error='Choose Extrude, Revolve, Sweep, or Loft.'
        if(operation<extrude_operation.or.operation>loft_operation)return
        error='Choose a curve or surface profile.'
        if(kind<1.or.kind>2.or.size(sections)==0)return
        error='Enter a name for the result.'
        if(len_trim(name)==0)return

        allocate(profiles(size(sections)))
        do i=1,size(profiles)
            call copy_single_patch_geometry(kind,sections(i),profiles(i)%geometry,catalog)
        end do
        if(present(spine))call copy_single_patch_geometry(1,spine,path%geometry,catalog)
        call stage_studio_construction(operation,profiles,objects(1),error,vector,axis_point,angle_degrees,path,degree)
        if(len(error)>0)return
        call publish_native_patches(objects,[trim(name)],created_kind,variant,error,operation='construction',catalog=catalog)
    end subroutine publish_construction

    subroutine construct_studio_primitive(primitive,name,origin,values,kind,variant,error,catalog)
        integer, intent(in) :: primitive
        character(*), intent(in) :: name
        real(rk), intent(in) :: origin(3),values(3)
        integer, intent(out) :: kind,variant
        character(:), allocatable, intent(out) :: error
        type(studio_patch_catalog), intent(inout), target, optional :: catalog
        type(studio_native_patch) :: objects(1)

        kind=0; variant=0
        error='Enter a name for the new object.'
        if(len_trim(name)==0)return
        call stage_studio_primitive(primitive,origin,values,objects(1),error)
        if(len(error)>0)return
        call publish_native_patches(objects,[trim(name)],kind,variant,error,operation='primitive creation',catalog=catalog)
    end subroutine construct_studio_primitive

end module studio_construction
