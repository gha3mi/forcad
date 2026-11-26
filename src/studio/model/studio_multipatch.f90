! Bridge the scene catalog to the native ForCAD multipatch containers.
module studio_multipatch
    use studio_forcad, only: multipatch_join, multipatch_side_label, initialize_assembly, &
        append_assembly_patch, append_assembly_join, inspect_native_assembly
    use studio_single_patch, only: copy_single_patch_geometry, single_patch_name, studio_native_patch, studio_patch_catalog
    implicit none(type, external)
    private

    public :: multipatch_draft, multipatch_join, multipatch_result, multipatch_join_result
    public :: inspect_multipatch, multipatch_side_label

    type :: multipatch_draft
        character(:), allocatable :: name
        integer :: kind=0
        integer, allocatable :: members(:)
        type(multipatch_join), allocatable :: joins(:)
    end type multipatch_draft

    type :: multipatch_join_result
        logical :: recorded=.false., valid=.false.
        character(:), allocatable :: message
    end type multipatch_join_result

    type :: multipatch_result
        logical :: valid=.false., rational=.false.
        integer :: patches=0, connections=0, controls=0, shared_controls=0
        character(:), allocatable :: message
        type(multipatch_join_result), allocatable :: joins(:)
    end type multipatch_result

contains

    subroutine inspect_multipatch(draft,result,native,patches,catalog,patch_ids)
        type(multipatch_draft), intent(in) :: draft
        type(multipatch_result), intent(out) :: result
        class(*), allocatable, intent(out), optional :: native
        type(studio_native_patch), optional, intent(in) :: patches(:)
        integer, optional, intent(in) :: patch_ids(:)
        type(studio_patch_catalog), intent(in), optional :: catalog
        class(*), allocatable :: group,geometry
        integer :: member,variant,id,i,a,b,join_count,override

        join_count=0
        if(allocated(draft%joins))join_count=size(draft%joins)
        allocate(result%joins(join_count))
        do i=1,join_count
            result%joins(i)%message='Check the group members first.'
        end do
        result%message='Choose one patch dimension for the group.'
        call initialize_assembly(draft%kind,group)
        if(.not.allocated(group))return
        result%message='Select at least two patches of the same dimension.'
        if(.not.allocated(draft%members))return
        if(size(draft%members)<2)return

        do member=1,size(draft%members)
            variant=draft%members(member)
            result%message='A group cannot contain the same patch twice.'
            if(any(draft%members(:member-1)==variant))return
            result%message='A member is missing from the project.'
            override=0
            if(present(patches))then
                override=variant
                if(present(patch_ids))override=findloc(patch_ids,variant,dim=1)
            end if
            if(override/=0)then
                if(override<1.or.override>size(patches))return
                if(.not.allocated(patches(override)%geometry))return
                call append_assembly_patch(group,patches(override)%geometry,id)
                result%message='The kernel rejected a supplied patch.'
            else
                call copy_single_patch_geometry(draft%kind,variant,geometry,catalog=catalog)
                if(.not.allocated(geometry))return
                call append_assembly_patch(group,geometry,id)
                result%message='The kernel rejected '//single_patch_name(draft%kind,variant,catalog=catalog)//'.'
            end if
            if(id/=member)return
            result%patches=id
        end do

        do i=1,join_count
            associate(spec=>draft%joins(i),check=>result%joins(i))
                a=findloc(draft%members,spec%patch_a,dim=1)
                b=findloc(draft%members,spec%patch_b,dim=1)
                check%message='Both sides must reference members of this group.'
                if(min(a,b)==0)cycle
                check%message='Choose a valid boundary on each member.'
                ! Supply valid array indices even in ForCAD builds with diagnostic logging disabled.
                if(min(spec%side_a,spec%side_b)<1.or.max(spec%side_a,spec%side_b)>2*draft%kind)cycle
                check%message='This patch dimension has no face orientation controls.'
                if(draft%kind<3.and.(spec%swap.or.any(spec%flip)))cycle
                call append_assembly_join(group,draft%kind,spec,a,b,check%recorded,check%valid)
                check%message='The kernel rejected the continuity settings or boundary bases.'
                if(.not.check%recorded)cycle
                result%connections=result%connections+1
                check%message='The boundary geometry fails the requested continuity.'
                if(.not.check%valid)cycle
                check%message='Continuity verified.'
                if(spec%continuity<0)check%message='Registered without continuity constraints.'
            end associate
        end do

        result%message='Review the interface checks and correct rejected joins.'
        if(.not.all(result%joins%valid))return
        result%message='The kernel rejected the combined multipatch model.'
        call inspect_native_assembly(group,result%valid,result%rational,result%controls,result%shared_controls)
        if(.not.result%valid)return
        result%message='Multipatch validated.'
        if(join_count==0)result%message='Members validated. No interfaces defined.'
        if(present(native))call move_alloc(group,native)
    end subroutine inspect_multipatch

end module studio_multipatch
