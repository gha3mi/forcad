! Typed editor actions shared by the local MCP bridge and the ForCAD Studio Assistant.
module studio_mcp_actions
    use, intrinsic :: iso_fortran_env, only: int64
    use qt, only: QApplication_activeModalWidget, QJsonArray, QJsonObject, QJsonValue_String, QWidget
    use studio_state, only: studio_t
    use studio_edit_history, only: document_undo_stamp, document_can_redo
    use studio_ui, only: invoke_command, available_assistant_tools, set_workspace_units, copy_project_selection, paste_project_selection
    use studio_units, only: unitless, valid_unit_id
    use studio_mcp_arguments, only: only_fields, bounded_integer
    use studio_mcp_refinement, only: preview_refinement
    use studio_mcp_controls, only: preview_mcp_controls
    use studio_mcp_inspection, only: inspect_scene, inspect_patch, inspect_assembly, measure_geometry, describe_selection, active_rank, &
        default_scene_page, max_scene_page
    use studio_mcp_construction, only: build_model, construct_existing, build_curved_volume
    use studio_mcp_transforms, only: transform_scene, transform_selection
    use studio_mcp_geometry_edits, only: edit_patch, copy_patch, rename_patch, delete_existing, assemble_patches, mesh_selection, preview_mesh
    use studio_mcp_view, only: select_object, set_material, set_visibility, set_probe, select_mesh_cell, inspect_view, set_camera, set_parameter_view, set_mesh_view, set_workspace
    use studio_mcp_exchange, only: import_mcp_file, export_mcp_file, open_mcp_project, save_mcp_project, new_mcp_project, &
        list_mcp_project_tabs, switch_mcp_project_tab, close_mcp_project_tab
    use studio_mcp_selection, only: set_selection, select_region
    use studio_assistant_actions, only: assistant_action, assistant_action_context, available_assistant_actions, &
        resolve_assistant_action
    use studio_assistant_guidance, only: studio_skill_index,studio_skill_text
    use studio_commands, only: undo_patch_command, redo_patch_command, copy_geometry_command, paste_geometry_command
    use studio_knot_interaction, only: knot_edit_in_progress
    use studio_live_viewport, only: transform_dragging
    implicit none(type, external)
    private
    public :: execute_studio_mcp

contains

    subroutine execute_studio_mcp(s,tool,arguments,result)
        type(studio_t),intent(inout),target :: s
        character(*),intent(in) :: tool
        type(QJsonObject),intent(in) :: arguments
        type(QJsonObject),intent(out) :: result
        type(QWidget) :: modal
        character(:),allocatable :: error
        integer :: offset,limit

        result=QJsonObject()
        error=''
        select case(tool)
        case('inspect')
            offset=bounded_integer(arguments,'offset',0,huge(0),fallback=0)
            limit=bounded_integer(arguments,'limit',1,max_scene_page,fallback=default_scene_page)
            if(only_fields(arguments,'|offset|limit|').and.offset>=0.and.limit>=1)then
                call inspect_scene(s,result,offset,limit)
            else
                error='Supply a nonnegative offset and an optional page limit from 1 to 128.'
            end if
        case('inspect_patch')
            call inspect_patch(s,arguments,result,error)
        case('inspect_assembly')
            call inspect_assembly(s,arguments,result,error)
        case('inspect_view')
            call inspect_view(s,arguments,result,error)
        case('list_project_tabs')
            call list_mcp_project_tabs(s,arguments,result,error)
        case('preview_refinement')
            call preview_refinement(s,arguments,result,error)
        case('preview_controls')
            call preview_mcp_controls(s,arguments,result,error)
        case('preview_mesh')
            call preview_mesh(s,arguments,result,error)
        case('measure')
            call measure_geometry(s,arguments,result,error)
        case('list_commands')
            if(arguments%isEmpty())then
                call list_commands(s,result)
            else
                error='List commands takes an empty object.'
            end if
        case('studio_skill')
            call read_studio_skill(arguments,result,error)
        case default
            modal=QApplication_activeModalWidget()
            if(.not.modal%isNull().or.s%updating.or.s%dragging.or.s%control_dragging.or. &
                transform_dragging(s).or.knot_edit_in_progress(s))then
                error='ForCAD Studio is handling an active edit or dialog. Finish it before issuing a command.'
            else if(s%home_active.and.tool/='new_project'.and.tool/='open_project'.and.tool/='save_project'.and. &
                tool/='switch_project_tab'.and.tool/='close_project_tab')then
                error='Open or create a ForCAD Studio project before issuing a command.'
            else
                select case(tool)
                case('execute')
                    call execute_command(s,arguments,result,error)
                case('select')
                    call select_object(s,arguments,result,error)
                case('set_selection')
                    call set_selection(s,arguments,result,error)
                case('select_region')
                    call select_region(s,arguments,result,error)
                case('transform_scene')
                    call transform_scene(s,arguments,result,error)
                case('transform_selection')
                    call transform_selection(s,arguments,result,error)
                case('build_model')
                    call build_model(s,arguments,result,error)
                case('construct')
                    call construct_existing(s,arguments,result,error)
                case('delete')
                    call delete_existing(s,arguments,result,error)
                case('build_curved_volume')
                    call build_curved_volume(s,arguments,result,error)
                case('mesh')
                    call mesh_selection(s,arguments,result,error)
                case('set_material')
                    call set_material(s,arguments,result,error)
                case('set_visibility')
                    call set_visibility(s,arguments,result,error)
                case('rename_patch')
                    call rename_patch(s,arguments,result,error)
                case('set_probe')
                    call set_probe(s,arguments,result,error)
                case('set_camera')
                    call set_camera(s,arguments,result,error)
                case('set_parameter_view')
                    call set_parameter_view(s,arguments,result,error)
                case('set_workspace')
                    call set_workspace(s,arguments,result,error)
                case('set_project_units')
                    call set_project_units(s,arguments,result,error)
                case('set_mesh_view')
                    call set_mesh_view(s,arguments,result,error)
                case('pick_mesh_cell')
                    call select_mesh_cell(s,arguments,result,error)
                case('import_file')
                    call import_mcp_file(s,arguments,result,error)
                case('new_project')
                    call new_mcp_project(s,arguments,result,error)
                case('save_project')
                    call save_mcp_project(s,arguments,result,error)
                case('open_project')
                    call open_mcp_project(s,arguments,result,error)
                case('switch_project_tab')
                    call switch_mcp_project_tab(s,arguments,result,error)
                case('close_project_tab')
                    call close_mcp_project_tab(s,arguments,result,error)
                case('export_file')
                    call export_mcp_file(s,arguments,result,error)
                case('copy_patch')
                    call copy_patch(s,arguments,result,error)
                case('edit_patch')
                    call edit_patch(s,arguments,result,error)
                case('assemble')
                    call assemble_patches(s,arguments,result,error)
                case('undo')
                    call undo(s,arguments,result,error)
                case('redo')
                    call redo(s,arguments,result,error)
                case default
                    error='Unknown ForCAD Studio tool.'
                end select
            end if
        end select
        call result%insert('ok',len(error)==0)
        if(len(error)>0)call result%insert('error',error)
    end subroutine execute_studio_mcp

    subroutine read_studio_skill(args,result,error)
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(QJsonArray) :: skills
        character(:),allocatable :: id,markdown

        error=''
        if(args%isEmpty())then
            skills=studio_skill_index()
            call result%insert('skills',skills)
            call skills%delete()
            return
        end if
        error='Supply a nonempty skill id, or an empty object to list packaged skills.'
        if(.not.only_fields(args,'|id|').or.args%valueType('id')/=QJsonValue_String)return
        id=args%stringValue('id')
        if(len(id)==0.or.len_trim(id)/=len(id))return
        markdown=studio_skill_text(id)
        error='Unknown or unavailable ForCAD Studio skill ID. Call studio_skill with no arguments to list skills.'
        if(len(markdown)==0)return
        error=''
        call result%insert('id',id)
        call result%insert('markdown',markdown)
    end subroutine read_studio_skill

    subroutine list_commands(s,result)
        type(studio_t),intent(in) :: s
        type(QJsonObject),intent(inout) :: result
        type(assistant_action_context) :: context
        type(QJsonArray) :: list
        type(QJsonObject) :: command
        type(assistant_action), allocatable :: entries(:)
        integer :: i

        list=QJsonArray()
        if(.not.s%home_active)then
            context=available_assistant_tools(s)
            entries=available_assistant_actions(active_rank(s),context=context)
            do i=1,size(entries)
                if(.not.entries(i)%dialog_free)cycle
                command=QJsonObject()
                call command%insert('name',trim(entries(i)%key))
                call command%insert('description',trim(entries(i)%label))
                call list%append(command)
                call command%delete()
            end do
        end if
        call result%insert('commands',list)
        call list%delete()
    end subroutine list_commands

    subroutine execute_command(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(assistant_action_context) :: context
        character(:),allocatable :: key,label
        integer :: command

        error='Execute requires a command listed by list_commands.'
        if(.not.only_fields(args,'|command|'))return
        key=args%stringValue('command')
        context=available_assistant_tools(s)
        call resolve_assistant_action(key,active_rank(s),command,label,context=context,dialog_free_only=.true.)
        if(command<=0)then
            error='This command is unavailable for the current selection.'
            return
        end if
        select case(command)
        case(copy_geometry_command)
            call copy_project_selection(s,error)
        case(paste_geometry_command)
            call paste_project_selection(s,error,same_folder=.true.)
        case default
            call invoke_command(s,command)
            error=''
        end select
        if(len(error)>0)return
        call result%insert('command',key)
        call result%insert('description',label)
        call describe_selection(s,result)
        error=''
    end subroutine execute_command

    subroutine set_project_units(s,args,result,error)
        type(studio_t),intent(inout) :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        character(8),parameter :: model_units(0:5)=[character(8)::'unitless','mm','cm','m','in','ft']
        character(4),parameter :: display_units(0:5)=[character(4)::'same','mm','cm','m','in','ft']
        integer :: model,display

        error='Set model units, display units, or both.'
        if(args%isEmpty().or..not.only_fields(args,'|model|display|'))return
        model=s%document%model_unit
        display=s%document%display_unit
        if(args%contains('model'))then
            error='Choose unitless, mm, cm, m, in or ft for the model unit.'
            if(args%valueType('model')/=QJsonValue_String)return
            model=findloc(model_units,args%stringValue('model'),dim=1)-1
            if(.not.valid_unit_id(model))return
        end if
        if(args%contains('display'))then
            error='Choose same, mm, cm, m, in or ft for the display unit.'
            if(args%valueType('display')/=QJsonValue_String)return
            display=findloc(display_units,args%stringValue('display'),dim=1)-1
            if(.not.valid_unit_id(display))return
            if(model==unitless.and.display/=unitless)then
                error='A unitless project uses the model unit for display.'
                return
            end if
        end if
        if(model==unitless)display=unitless
        call set_workspace_units(s,model,display,error)
        if(len(error)>0)return
        call result%insert('model_unit',trim(model_units(model)))
        call result%insert('display_unit',trim(display_units(display)))
        error=''
    end subroutine set_project_units

    subroutine undo(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        integer(int64) :: previous

        error='Undo takes an empty object.'
        if(.not.args%isEmpty())return
        previous=document_undo_stamp(s%document)
        error='There is no undoable edit.'
        if(previous==0_int64)return
        call invoke_command(s,undo_patch_command)
        error='ForCAD Studio could not undo the latest edit.'
        if(document_undo_stamp(s%document)==previous)return
        call result%insert('undo_available',document_undo_stamp(s%document)>0_int64)
        call result%insert('redo_available',document_can_redo(s%document))
        call describe_selection(s,result)
        error=''
    end subroutine undo

    subroutine redo(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        integer(int64) :: previous

        error='Redo takes an empty object.'
        if(.not.args%isEmpty())return
        error='There is no redoable edit.'
        if(.not.document_can_redo(s%document))return
        previous=document_undo_stamp(s%document)
        call invoke_command(s,redo_patch_command)
        error='ForCAD Studio could not redo the latest edit.'
        if(document_undo_stamp(s%document)==previous)return
        call result%insert('undo_available',document_undo_stamp(s%document)>0_int64)
        call result%insert('redo_available',document_can_redo(s%document))
        call describe_selection(s,result)
        error=''
    end subroutine redo

end module studio_mcp_actions
