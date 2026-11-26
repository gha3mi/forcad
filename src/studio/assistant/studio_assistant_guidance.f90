module studio_assistant_guidance
    use qt, only: QCoreApplication_applicationFilePath,QJsonArray,QJsonObject
    use studio_assistant_transport, only: read_text
    use studio_platform_files, only: resolved_file_path
    use studio_project_folder, only: canonical_project_folder
    implicit none(type,external)
    private
    public :: studio_assistant_instructions,studio_skill_index,studio_skill_text

    character,parameter :: lf=new_line('a')
    character(*),parameter :: installed_guidance='../share/forcad-studio/assistant/'
    character(24),parameter :: skill_ids(*)=[character(24) :: 'scene-inspection','geometry-construction', &
        'patch-editing','multipatch','structured-meshing','project-exchange']
    character(*),parameter :: live_fallback= &
        'Use ForCAD Studio tools for live CAD edits. Discover workflow guides with studio_skill and read the relevant guide by ID. '// &
        'Use inspect for object IDs and units, inspect_view for viewport state, and list_commands for context-dependent GUI actions. '// &
        'Use explicit targets when supported. Preview control edits and assembly refinement with the same arguments before applying. '// &
        'Follow approvals and preserve Undo. Wait for each result before issuing a dependent edit. '// &
        'Do not change project files as a substitute for approved ForCAD Studio edits. Confirm changes from tool results.'
    integer,parameter :: max_guide_bytes=32768

contains

    function studio_assistant_instructions(live,role) result(text)
        logical,intent(in) :: live
        character(*),intent(in) :: role
        character(:),allocatable :: text,brief

        text=''
        if(live)then
            text=read_guidance('AGENTS.md')
            if(len_trim(text)==0)text=live_fallback
        end if
        select case(role)
        case('Lead')
            brief=read_guidance('agents/lead.md')
            if(len_trim(brief)==0)brief='Lead the work in ForCAD Studio. Recheck the active scene after reviewing the team proposals.'
        case('Geometry reviewer')
            brief=read_guidance('agents/geometry-reviewer.md')
            if(len_trim(brief)==0)brief='Review the geometry and reference dimensions. Identify what must be validated. Return a proposal without edits.'
        case('Modeling planner')
            brief=read_guidance('agents/modeling-planner.md')
            if(len_trim(brief)==0)brief='Propose ordered native modeling operations, their dependencies, and checks. Return a proposal without edits.'
        case default
            brief=''
        end select
        if(len_trim(brief)>0)then
            if(len(text)>0)text=text//lf//lf
            text=text//trim(brief)
        end if
    end function studio_assistant_instructions

    function studio_skill_index() result(catalog)
        type(QJsonArray) :: catalog
        type(QJsonObject) :: item
        character(:),allocatable :: markdown,title,description,id
        integer :: i

        catalog=QJsonArray()
        do i=1,size(skill_ids)
            id=trim(skill_ids(i))
            markdown=studio_skill_text(id)
            if(len_trim(markdown)==0)cycle
            title=guide_field(markdown,'# ')
            description=guide_field(markdown,'description:')
            item=QJsonObject()
            if(len(title)==0)title=id
            call item%insert('id',id)
            call item%insert('title',title)
            call item%insert('description',description)
            call catalog%append(item)
            call item%delete()
        end do
    end function studio_skill_index

    function studio_skill_text(id) result(markdown)
        character(*),intent(in) :: id
        character(:),allocatable :: markdown
        integer :: index

        markdown=''
        index=findloc(skill_ids,id,dim=1)
        if(index>0)markdown=read_guidance('skills/'//trim(skill_ids(index))//'/SKILL.md')
    end function studio_skill_text

    function read_guidance(relative) result(text)
        character(*),intent(in) :: relative
        character(:),allocatable :: text,executable,project_root
        integer :: separator,root_length,status

        text=''
        executable=resolved_file_path(QCoreApplication_applicationFilePath())
        separator=scan(executable,'/',back=.true.)
        if(separator>0)text=read_text(executable(:separator)//installed_guidance//relative,max_guide_bytes)
        if(len_trim(text)>0)return

        call get_environment_variable('PIXI_PROJECT_ROOT',length=root_length,status=status)
        if(status/=0.or.root_length<=0)return
        allocate(character(root_length)::project_root)
        call get_environment_variable('PIXI_PROJECT_ROOT',value=project_root,status=status)
        if(status/=0.or.len_trim(project_root)==0)return
        project_root=canonical_project_folder(project_root)
        if(len(project_root)==0)return
        if(index(executable,project_root//'/build/')/=1)return
        text=read_text(project_root//'/assistant/'//relative,max_guide_bytes)
    end function read_guidance

    function guide_field(markdown,prefix) result(value)
        character(*),intent(in) :: markdown,prefix
        character(:),allocatable :: value
        integer :: first,last,next

        value=''
        first=1
        do while(first<=len(markdown))
            next=index(markdown(first:),lf)
            last=len(markdown)
            if(next>0)last=first+next-2
            if(last-first+1>=len(prefix))then
                if(markdown(first:first+len(prefix)-1)==prefix)then
                    if(last>=first)then
                        if(markdown(last:last)==achar(13))last=last-1
                    end if
                    value=trim(adjustl(markdown(first+len(prefix):last)))
                    return
                end if
            end if
            first=last+2
        end do
    end function guide_field
end module studio_assistant_guidance
