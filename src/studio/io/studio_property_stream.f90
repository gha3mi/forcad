module studio_property_stream
    use, intrinsic :: iso_fortran_env, only: int64
    use studio_history_storage, only: history_blob,create_history_blob,open_history_blob,release_history_blob
    use studio_project_stream, only: remaining_file_units,integer_units
    use studio_units, only: valid_unit_id
    implicit none(type, external)
    private
    public :: studio_property_values,save_property_history,load_property_history

    type :: studio_property_values
        integer :: units(2)=-1
        integer, allocatable :: patches(:,:),meshes(:,:)
    end type studio_property_values

contains

    logical function save_property_history(values,blob,project) result(saved)
        type(studio_property_values), intent(in) :: values
        type(history_blob), intent(out) :: blob
        character(*), intent(in), optional :: project
        integer :: unit,status,closed,np,nm

        np=0
        nm=0
        if(allocated(values%patches))np=size(values%patches,2)
        if(allocated(values%meshes))nm=size(values%meshes,2)
        call create_history_blob(blob,unit,status,project)
        saved=.false.
        if(status/=0)return
        write(unit,iostat=status)values%units,np,nm
        if(status==0.and.np>0)write(unit,iostat=status)values%patches
        if(status==0.and.nm>0)write(unit,iostat=status)values%meshes
        close(unit,iostat=closed)
        saved=status==0.and.closed==0
        if(.not.saved)call release_history_blob(blob)
    end function save_property_history

    logical function load_property_history(blob,values) result(loaded)
        type(history_blob), intent(in) :: blob
        type(studio_property_values), intent(out) :: values
        integer :: unit,status,closed,np,nm

        call open_history_blob(blob,unit,status)
        loaded=.false.
        if(status/=0)return
        np=0
        nm=0
        read(unit,iostat=status)values%units,np,nm
        if(status==0)then
            if(min(np,nm)<0.or.7_int64*np+2_int64*nm/=remaining_file_units(unit)/integer_units)status=1
        end if
        if(status==0)allocate(values%patches(7,np),stat=status)
        if(status==0)allocate(values%meshes(2,nm),stat=status)
        if(status==0.and.np>0)read(unit,iostat=status)values%patches
        if(status==0.and.nm>0)read(unit,iostat=status)values%meshes
        close(unit,iostat=closed)
        if(status/=0.or.closed/=0)return
        loaded=valid_values(values)
    end function load_property_history

    logical function valid_values(values) result(valid)
        type(studio_property_values), intent(in) :: values

        valid=.false.
        if(.not.all(values%units==-1))then
            if(.not.valid_unit_id(values%units(1)).or..not.valid_unit_id(values%units(2)))return
        end if
        if(any(values%patches(1,:)<1).or.any(values%patches(1,:)>3).or.any(values%patches(2,:)<1))return
        if(any(values%patches(3,:)<0).or.any(values%patches(3,:)>1))return
        if(any(values%meshes(1,:)<1).or.any(values%meshes(2,:)<0).or.any(values%meshes(2,:)>1))return
        valid=.true.
    end function valid_values

end module studio_property_stream
