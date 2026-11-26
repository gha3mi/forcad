! Geometry definitions for the editable examples in ForCAD Studio.
module studio_sample_geometry
    use studio_forcad, only: make_sample_patch
    implicit none(type, external)
    private
    public :: make_sample_patch, sample_name

contains

    pure function sample_name(kind,sample_index) result(name)
        integer, intent(in) :: kind,sample_index
        character(:), allocatable :: name
        character(len=32), parameter :: names(3,2)=reshape([character(len=32) :: &
            'Quadratic Bézier curve', 'Unit-square surface', 'Unit-cube volume', &
            'Waveguide centerline', 'Saddle canopy', 'Curved duct volume'],[3,2])

        name=''
        if (kind<1.or.kind>size(names,1).or.sample_index<1.or.sample_index>size(names,2)) return
        name=trim(names(kind,sample_index))
    end function sample_name


end module studio_sample_geometry
