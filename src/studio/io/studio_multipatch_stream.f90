module studio_multipatch_stream
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use studio_multipatch, only: multipatch_draft
    use studio_history_storage, only: history_blob,create_history_blob,open_history_blob,release_history_blob
    use studio_project_stream, only: write_text,read_text,remaining_file_units,integer_units,real_units
    implicit none(type, external)
    private
    public :: write_group,read_group,save_group_history,load_group_history

contains

    logical function save_group_history(groups,blob,project) result(saved)
        type(multipatch_draft), intent(in) :: groups(:)
        type(history_blob), intent(out) :: blob
        character(*), intent(in), optional :: project
        integer :: unit,status,closed

        call create_history_blob(blob,unit,status,project)
        saved=.false.
        if(status/=0)return
        call write_groups(unit,groups,status)
        close(unit,iostat=closed)
        saved=status==0.and.closed==0
        if(.not.saved)call release_history_blob(blob)
    end function save_group_history

    logical function load_group_history(blob,groups) result(loaded)
        type(history_blob), intent(in) :: blob
        type(multipatch_draft), allocatable, intent(out) :: groups(:)
        integer :: unit,status,closed

        call open_history_blob(blob,unit,status)
        loaded=.false.
        if(status/=0)return
        call read_groups(unit,groups,status)
        close(unit,iostat=closed)
        loaded=status==0.and.closed==0
    end function load_group_history

    subroutine write_groups(unit,groups,status)
        integer, intent(in) :: unit
        type(multipatch_draft), intent(in) :: groups(:)
        integer, intent(out) :: status
        integer :: i

        write(unit,iostat=status)size(groups)
        if(status/=0)return
        do i=1,size(groups)
            call write_group(unit,groups(i),status)
            if(status/=0)return
        end do
    end subroutine write_groups

    subroutine read_groups(unit,groups,status)
        integer, intent(in) :: unit
        type(multipatch_draft), allocatable, intent(out) :: groups(:)
        integer, intent(out) :: status
        integer :: count,i

        read(unit,iostat=status)count
        if(status/=0)return
        status=1
        if(count<0.or.5_int64*count>remaining_file_units(unit)/integer_units)return
        allocate(groups(count),stat=status)
        if(status/=0)return
        do i=1,count
            call read_group(unit,groups(i),status)
            if(status/=0)return
        end do
    end subroutine read_groups

    subroutine write_group(unit,group,status)
        integer, intent(in) :: unit
        type(multipatch_draft), intent(in) :: group
        integer, intent(out) :: status
        integer :: join_count,index,jet_count,flags(5)

        call write_text(unit,group%name,status)
        if(status/=0)return
        join_count=0
        if(allocated(group%joins))join_count=size(group%joins)
        write(unit,iostat=status)group%kind,size(group%members),join_count,group%members
        if(status/=0)return
        do index=1,join_count
            associate(join=>group%joins(index))
                flags=merge(1,0,[join%reverse,join%swap,join%flip,join%geometric])
                jet_count=0
                if(allocated(join%reparameterization))jet_count=size(join%reparameterization)
                write(unit,iostat=status)join%patch_a,join%side_a,join%patch_b,join%side_b,join%continuity,flags,jet_count
                if(status/=0)return
                if(jet_count==0)cycle
                write(unit,iostat=status)join%reparameterization
                if(status/=0)return
            end associate
        end do
    end subroutine write_group

    subroutine read_group(unit,group,status)
        integer, intent(in) :: unit
        type(multipatch_draft), intent(out) :: group
        integer, intent(out) :: status
        integer :: member_count,join_count,index,jet_count,flags(5)

        call read_text(unit,group%name,status)
        if(status/=0)return
        read(unit,iostat=status)group%kind,member_count,join_count
        if(status/=0)return
        status=1
        if(group%kind<1.or.group%kind>3)return
        if(member_count<2.or.join_count<0)return
        if(int(member_count,int64)+11_int64*join_count>remaining_file_units(unit)/integer_units)return
        allocate(group%members(member_count),stat=status)
        if(status/=0)return
        allocate(group%joins(join_count),stat=status)
        if(status/=0)return
        read(unit,iostat=status)group%members
        if(status/=0)return
        do index=1,join_count
            associate(join=>group%joins(index))
                read(unit,iostat=status)join%patch_a,join%side_a,join%patch_b,join%side_b,join%continuity,flags,jet_count
                if(status/=0)return
                status=1
                if(.not.any(group%members==join%patch_a).or..not.any(group%members==join%patch_b))return
                if(min(join%side_a,join%side_b)<1.or.max(join%side_a,join%side_b)>2*group%kind)return
                if(join%continuity< -1)return
                if(any(flags<0).or.any(flags>1).or.jet_count<0)return
                if(jet_count>remaining_file_units(unit)/real_units)return
                join%reverse=flags(1)==1
                join%swap=flags(2)==1
                join%flip=flags(3:4)==1
                join%geometric=flags(5)==1
                status=0
                if(jet_count==0)cycle
                allocate(join%reparameterization(jet_count),stat=status)
                if(status/=0)return
                read(unit,iostat=status)join%reparameterization
                if(status/=0)return
                if(.not.all(ieee_is_finite(join%reparameterization)))then
                    status=1
                    return
                end if
            end associate
        end do
    end subroutine read_group

end module studio_multipatch_stream
