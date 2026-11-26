submodule(studio_ui) studio_ui_patch_tools
    use studio_selection_status, only: knot_status_text
    use studio_status_hints, only: publish_window_status
    use qt, only: QButtonGroup
    use studio_geometry_operations, only: commit_patch_section
    use studio_parameter_editor, only: parameter_direction_buttons
    use studio_patch_range_dialog, only: edit_studio_patch_section, edit_studio_parameter_interval, choose_studio_reverse_direction
    use studio_selection, only: pick_face
    use studio_theme, only: studio_style_primary
    implicit none(type, external)

    type :: refinement_inputs
        type(QSpinBox) :: subdivisions,degree
        type(QDialog) :: dialog
        integer :: kind, max_divisions(0:3), max_steps(0:3)
    end type refinement_inputs
contains
    module procedure on_control_open
        call invoke_command(s,edit_controls_command)
    end procedure on_control_open

    module procedure edit_workspace_controls
        type(multipatch_controls_t) :: shared
        integer :: assembly,point

        assembly=s%selected_multipatch
        point=s%selected_control
        if(assembly==0.and.s%live_scene%valid)then
            assembly=joined_patch_assembly(s%document,s%entity_kind,s%live_scene%variant)
            if(assembly>0)then
                call read_document_multipatch_controls(s%document,assembly,shared)
                point=assembly_control_id(shared,s%live_scene%variant,point)
            end if
        end if
        if(assembly>0)then
            call refresh_multipatch(s%document%multipatches,assembly,.false.,catalog=s%document%patches)
            if(.not.allocated(s%document%multipatches%groups(assembly)%native))then
                call show_patch_issue(s,'Resolve the multipatch interfaces before editing controls.',5000)
                return
            end if
            if(edit_multipatch_controls(s%window,s%document,assembly,point,active_length_unit(s)))then
                call refresh_edited_multipatch(s,assembly,point,'Multipatch control edit applied')
            end if
        else if(s%live_scene%valid)then
            if(edit_single_patch_controls(s%window,s%entity_kind,point, &
                s%document%model_unit,active_length_unit(s),catalog=s%document%patches,patch_id=s%live_scene%variant))then
                s%live_scene%valid=.false.
                call select_entity(s,s%entity_kind)
                call s%status%showMessage('Control edit applied',3500)
            end if
        end if
    end procedure edit_workspace_controls

    module procedure on_control_move
        type(single_patch_controls)::controls
        type(multipatch_controls_t)::shared
        logical::changed
        character(:),allocatable::error
        integer :: assembly,global_id,moved

        assembly=s%selected_multipatch
        global_id=point
        if(assembly>0)then
            shared=s%shared_controls
        else if(s%live_scene%valid)then
            assembly=joined_patch_assembly(s%document,s%entity_kind,s%live_scene%variant)
            if(assembly>0)then
                call read_document_multipatch_controls(s%document,assembly,shared)
                global_id=assembly_control_id(shared,s%live_scene%variant,point)
            end if
        end if
        if(assembly>0)then
            if(.not.shared%valid)then
                call show_patch_issue(s,'Resolve the multipatch interfaces before moving controls.',5000)
                return
            end if
            if(global_id<1.or.global_id>size(shared%points,2))return
            moved=1
            if(s%selected_multipatch==assembly.and.allocated(s%selected_shared_controls))then
                if(any(s%selected_shared_controls==global_id))moved=size(s%selected_shared_controls)
                call set_multipatch_control_displacement(shared,global_id,position,s%selected_shared_controls)
            else
                call set_multipatch_control_displacement(shared,global_id,position)
            end if
            call commit_document_multipatch_controls(s%document,assembly,shared,changed,error)
            if(len(error)>0)then
                call show_patch_issue(s,error,5000)
            else if(changed)then
                if(moved==1)then
                    call refresh_edited_multipatch(s,assembly,global_id,'Multipatch control point moved')
                else
                    call refresh_edited_multipatch(s,assembly,global_id,integer_text(moved)//' multipatch controls moved')
                end if
            end if
            return
        end if

        if(.not.s%live_scene%valid)return
        call read_single_patch_controls(s%entity_kind,controls,variant=s%live_scene%variant,catalog=s%document%patches)
        if(.not.controls%valid.or.point<1)return
        if(point>size(controls%weights))return
        controls%points(:,point)=position
        call commit_single_patch_controls(s%entity_kind,controls,changed,error,catalog=s%document%patches,patch_id=s%live_scene%variant)
        if(len(error)>0)then
            call show_patch_issue(s,error,5000)
        else if(changed)then
            s%live_scene%valid=.false.
            call select_entity(s,s%entity_kind)
            call s%status%showMessage('Control point moved',3500)
        end if
    end procedure on_control_move

    subroutine refresh_edited_multipatch(s,assembly,point,message)
        type(studio_t), intent(inout), target :: s
        integer, intent(in) :: assembly,point
        character(*), intent(in) :: message

        s%live_scene%valid=.false.
        if(s%selected_multipatch/=assembly)then
            call apply_workspace_selection(s,selection_target_t(selection_assembly,assembly),selection_replace)
        else
            call refresh_inspector(s)
        end if
        if(s%shared_controls%valid)then
            s%selected_control=point
            if(point>0.and..not.allocated(s%selected_shared_controls))s%selected_shared_controls=[point]
            call on_control_selection(s)
        end if
        call s%status%showMessage(message,3500)
        call repaint_all(s)
    end subroutine refresh_edited_multipatch

    module procedure on_interval_apply

        select type(s=>data)
        type is(studio_t)
            call invoke_command(s,apply_interval_command)
        end select
    end procedure on_interval_apply

    module procedure run_knot_editor
        real(dp)::selected,edited,probe(3)
        integer::accepted,requested,row,axis
        character(len=1)::axis_name
        character(len=:),allocatable::message

        if(.not.s%live_scene%valid)return
        if(removing.and..not.can_remove_knot(s,any_direction=.true.))return
        axis=s%knot_axis
        probe=[s%u,s%v,s%w]
        selected=s%live_scene%bounds(1,axis)
        row=s%knot_table%currentRow()+1
        if(row>0.and.row<=size(s%live_knot_groups))selected=s%live_knot_groups(row)

        if(.not.edit_single_patch_knots(s%window,s%entity_kind,axis,probe,removing,accepted,requested,selected,edited, &
            catalog=s%document%patches,patch_id=s%live_scene%variant))return
        s%live_scene%valid=.false.
        call select_entity(s,s%entity_kind)
        if(s%knot_axis==axis.and.allocated(s%live_knot_groups))then
            row=findloc(s%live_knot_groups,edited,dim=1)-1
            if(row>=0)call s%knot_table%setCurrentCell(row,1)
        end if
        axis_name=parameter_directions(axis:axis)

        if(removing)then
            message='Removed '
        else
            message='Inserted '
        end if
        message=message//integer_text(accepted)
        if(accepted<requested)message=message//' of '//integer_text(requested)
        message=message//' knot'
        if(requested/=1)message=message//'s'
        call s%status%showMessage(message//' in '//axis_name,5000)
    end procedure run_knot_editor

    module procedure run_patch_tool
        integer::kind,source,created,axis
        integer(int64)::revision,document_revision
        real(dp)::lower,upper,probe(3)
        logical::changed,splitting
        character(:),allocatable::name,error,message

        if(.not.s%live_scene%valid.or.s%selected_multipatch>0)return
        kind=s%entity_kind
        source=s%live_scene%variant
        if(command/=extract_patch_command.and.document_patch_in_assembly(s%document,kind,source))then
            call s%status%showMessage('Detach this patch from its multipatch before splitting or reversing it.',5000)
            return
        end if
        axis=s%knot_axis
        name=single_patch_name(kind,source,catalog=s%document%patches)
        revision=single_patch_revision(kind,source,catalog=s%document%patches)
        document_revision=s%document%document_revision
        if(command==reverse_patch_command)then
            if(kind>1)then
                if(.not.choose_studio_reverse_direction(s%window,name,kind,axis,s%live_scene%bounds))return
            end if
            if(.not.patch_tool_target_current(s,kind,source,revision,document_revision))return
            if(.not.can_restructure_patch(s))return
            call commit_single_patch_reverse(kind,axis,changed,error,catalog=s%document%patches,patch_id=source)
            if(len(error)>0)then
                call show_patch_issue(s,error,6500)
                return
            end if
            if(.not.changed)return
            message='Reversed '//parameter_directions(axis:axis)//' direction · '//name
        else
            splitting=command==split_patch_command
            probe=s%live_scene%bounds(1,:)
            if(s%live_probe%valid)probe=s%live_probe%parameters
            if(.not.edit_studio_patch_range(s%window,name,kind,axis,s%live_scene%bounds,probe,splitting,lower,upper))return
            if(.not.patch_tool_target_current(s,kind,source,revision,document_revision))return
            call apply_patch_range(s,axis,lower,upper,splitting,created,error)
            if(len(error)>0)call show_patch_issue(s,error,6500)
            return
        end if
        call finish_patch_tool(s,kind,s%knot_axis,message)
    end procedure run_patch_tool

    logical function patch_tool_target_current(s,kind,source,revision,document_revision) result(current)
        type(studio_t), intent(inout) :: s
        integer, intent(in) :: kind,source
        integer(int64), intent(in) :: revision,document_revision

        current=s%document%document_revision==document_revision.and.s%live_scene%valid.and.s%entity_kind==kind.and. &
            s%live_scene%variant==source.and.s%selected_multipatch==0.and.s%selected_mesh==0.and.s%selection%count()<=1
        if(current)current=single_patch_revision(kind,source,catalog=s%document%patches)==revision
        if(.not.current)call s%status%showMessage('The patch changed while the dialog was open. Run the command again.',5500)
    end function patch_tool_target_current

    module procedure apply_patch_range
        integer :: kind,source
        character(:),allocatable :: name,message

        created=0
        error='Select a patch and one of its parameter directions.'
        if(.not.s%live_scene%valid.or.s%selected_multipatch>0.or.s%selected_mesh>0)return
        kind=s%entity_kind
        if(axis<1.or.axis>kind)return
        source=s%live_scene%variant
        name=single_patch_name(kind,source,catalog=s%document%patches)
        call commit_patch_range(s%document,kind,source,axis,lower,upper,splitting,created,error)
        if(len(error)>0.or.created==0)return
        if(splitting)then
            message='Split '//name//' at '//parameter_directions(axis:axis)//' = '//number_text(lower)
        else
            message='Extracted '//parameter_directions(axis:axis)//' range from '//name
        end if
        call finish_patch_tool(s,kind,axis,message)
    end procedure apply_patch_range

    module procedure run_patch_section
        character(5), parameter :: sides(6)=[character(5) :: 'U min','U max','V min','V max','W min','W max']
        character(:), allocatable :: name,error
        integer :: kind,source,axis,side,created
        integer(int64) :: revision,document_revision
        real(dp) :: parameter,probe(3),bounds(2,3)

        if(.not.s%live_scene%valid.or.s%selected_multipatch>0.or.s%selected_mesh>0)return
        kind=s%entity_kind
        if(kind<2.or.(face.and.kind/=3))return
        source=s%live_scene%variant
        axis=max(1,min(kind,s%knot_axis))
        revision=single_patch_revision(kind,source,catalog=s%document%patches)
        document_revision=s%document%document_revision
        bounds=s%live_scene%bounds
        name=single_patch_name(kind,source,catalog=s%document%patches)
        if(face)then
            side=0
            if(all(s%viewport_feature(:3)==[kind,source,pick_face]))side=s%viewport_feature(4)
            if(side<1.or.side>size(sides))then
                side=choose_volume_face(s%window,name,sides,bounds,axis)
                if(side<1.or.side>size(sides))return
            end if
            axis=(side+1)/2
            parameter=bounds(2-mod(side,2),axis)
        else
            probe=0.5_dp*bounds(1,:)+0.5_dp*bounds(2,:)
            if(s%live_probe%valid)probe=s%live_probe%parameters
            if(.not.edit_studio_patch_section(s%window,name,kind,axis,bounds,probe,parameter))return
        end if
        if(s%document%document_revision/=document_revision.or.s%entity_kind/=kind.or.s%live_scene%variant/=source.or. &
            single_patch_revision(kind,source,catalog=s%document%patches)/=revision.or. &
            s%selected_multipatch>0.or.s%selected_mesh>0)then
            call s%status%showMessage('The source patch changed. Open Extract again.',5500)
            return
        end if
        call apply_patch_section(s,axis,parameter,created,error)
        if(len(error)>0)call show_patch_issue(s,error,6500)
    end procedure run_patch_section

    integer function choose_volume_face(parent,name,sides,bounds,axis) result(face)
        class(QWidget), intent(in) :: parent
        character(*), intent(in) :: name,sides(6)
        real(dp), intent(in) :: bounds(2,3)
        integer, intent(in) :: axis
        type(QDialog), target :: dialog
        type(QVBoxLayout) :: layout
        type(QGridLayout) :: boundaries
        type(QLabel) :: heading,source
        type(QPushButton) :: option
        type(QButtonGroup) :: choices
        type(QDialogButtonBox) :: actions
        integer :: i,direction
        character(:), allocatable :: caption

        dialog=QDialog(parent)
        call dialog%setObjectName('studioExtractFaceDialog')
        call dialog%setWindowTitle('Extract face')
        call dialog%setMinimumWidth(388)
        layout=QVBoxLayout(dialog)
        call layout%setContentsMargins(16,14,16,14)
        call layout%setSpacing(8)
        heading=QLabel('Extract face',dialog)
        call heading%setObjectName('studioDialogTitle')
        call layout%addWidget(heading)
        source=QLabel('Choose a boundary of '//name//'. The face becomes a new surface.',dialog)
        call source%setTextFormat(Qt_PlainText)
        call source%setWordWrap(.true.)
        call source%setProperty('dialogRole','detail')
        call layout%addWidget(source)
        boundaries=QGridLayout()
        call boundaries%setSpacing(6)
        choices=QButtonGroup(dialog)
        do i=1,size(sides)
            direction=(i+1)/2
            caption=trim(sides(i))//'  ·  '//number_text(bounds(2-mod(i,2),direction))
            option=QPushButton(caption,dialog)
            call option%setObjectName('studioExtractBoundary'//integer_text(i))
            call option%setAccessibleName('Extract '//caption)
            call option%setMinimumHeight(30)
            call option%setProperty('autoDefault',.false.)
            call choices%addButton(option,i)
            call boundaries%addWidget(option,direction-1,mod(i-1,2))
            if(i==2*axis-1)call option%setFocus()
        end do
        call layout%addLayout(boundaries)
        call choices%idClicked%connect(extract_volume_face_choice,dialog)
        actions=QDialogButtonBox(QDialogButtonBox_Cancel,dialog)
        call actions%rejected%connect(cancel_refinement,dialog)
        call layout%addWidget(actions)
        face=dialog%exec()
        call dialog%delete()
    end function choose_volume_face

    subroutine extract_volume_face_choice(face,data)
        integer, intent(in) :: face
        class(*), intent(inout), optional :: data

        select type(dialog=>data)
        type is(QDialog)
            call dialog%done(face)
        end select
    end subroutine extract_volume_face_choice

    module procedure apply_patch_section
        integer :: kind,source
        character(:), allocatable :: name,message

        created=0
        error='Select one surface or volume to extract a section.'
        if(.not.s%live_scene%valid.or.s%selected_multipatch>0.or.s%selected_mesh>0.or.s%selection%count()>1)return
        kind=s%entity_kind
        source=s%live_scene%variant
        name=single_patch_name(kind,source,catalog=s%document%patches)
        call commit_patch_section(s%document,kind,source,axis,parameter,created,error)
        if(len(error)>0.or.created==0)return
        message='Extracted '//parameter_directions(axis:axis)//' = '//number_text(parameter)//' from '//name
        call finish_patch_tool(s,kind-1,1,message)
    end procedure apply_patch_section

    subroutine show_patch_issue(s,error,duration)
        type(studio_t), intent(inout) :: s
        character(*), intent(in) :: error
        integer, intent(in) :: duration

        call record_studio_diagnostic(s,error)
        call refresh_studio_diagnostics(s)
        call s%status%showMessage(error,duration,warning=.true.)
    end subroutine show_patch_issue

    subroutine finish_patch_tool(s,kind,axis,message)
        type(studio_t),intent(inout),target :: s
        integer,intent(in) :: kind,axis
        character(*),intent(in) :: message

        call on_transform_exit(s)
        s%selected_control=0
        s%knot_axis=axis
        s%live_scene%valid=.false.
        call rebuild_project_tree(s)
        call select_patch(s,kind,reveal=.false.)
        call s%status%showMessage(message//' · Ctrl+Z to undo',6500)
    end subroutine finish_patch_tool

    module procedure run_degree_editor
        real(dp)::bound,probe(3)
        integer::axis,accepted,requested,initial_degree(3)
        character(len=1)::direction
        character(len=:),allocatable::message

        if(.not.s%live_scene%valid)return
        if(reducing.and..not.can_reduce_degree(s))return
        axis=s%knot_axis
        initial_degree=s%live_scene%degree
        probe=[s%u,s%v,s%w]
        if(.not.edit_single_patch_degree(s%window,s%entity_kind,axis,probe,reducing,accepted,requested,bound, &
            catalog=s%document%patches,patch_id=s%live_scene%variant))return
        s%live_scene%valid=.false.
        call select_entity(s,s%entity_kind)
        if(reducing)then
            message='Reduced '
        else
            message='Elevated '
        end if
        direction=parameter_directions(axis:axis)
        message=message//direction//' degree '//integer_text(initial_degree(axis))//' → '//integer_text(s%live_scene%degree(axis))
        if(accepted<requested)message=message//' · '//integer_text(accepted)//' of '//integer_text(requested)//' steps'
        if(reducing)message=message//' · Error bound '//compact_scalar(bound)//' model units'
        call s%status%showMessage(message,6500)
    end procedure run_degree_editor

    module procedure run_patch_interval
        integer :: kind,source,axis
        integer(int64) :: revision,document_revision
        real(dp) :: lower,upper
        character(:), allocatable :: name

        if(.not.s%live_scene%valid.or.s%selected_multipatch>0.or.s%selected_mesh>0)return
        kind=s%entity_kind
        if(axis_hint<1.or.axis_hint>kind)return
        axis=axis_hint
        source=s%live_scene%variant
        name=single_patch_name(kind,source,catalog=s%document%patches)
        revision=single_patch_revision(kind,source,catalog=s%document%patches)
        document_revision=s%document%document_revision
        if(.not.edit_studio_parameter_interval(s%window,name,kind,axis,s%live_scene%bounds,lower,upper))return
        if(.not.patch_tool_target_current(s,kind,source,revision,document_revision))return
        call apply_parameter_interval(s,axis,lower,upper)
    end procedure run_patch_interval

    module procedure apply_inline_interval
        real(dp)::lower,upper
        integer::axis

        if(.not.s%live_scene%valid.or.s%selected_multipatch>0)return
        axis=s%knot_axis
        if(any(s%widgets%interval_context/=[s%entity_kind,s%live_scene%variant,axis]).or. &
            s%widgets%interval_revision/=s%document%patches%patch_revision(s%widgets%interval_context(1),s%widgets%interval_context(2)))then
            call refresh_inline_interval(s)
            call s%status%showMessage('The selected interval changed. Review the current bounds.',5000)
            return
        end if
        if(.not.parse_parameter_number(s%widgets%interval_validator,s%widgets%interval_lower%text(),lower))then
            call s%widgets%interval_lower%setFocus()
            call s%status%showMessage('Enter a finite lower bound.',5500)
            return
        end if
        if(.not.parse_parameter_number(s%widgets%interval_validator,s%widgets%interval_upper%text(),upper))then
            call s%widgets%interval_upper%setFocus()
            call s%status%showMessage('Enter a finite upper bound.',5500)
            return
        end if
        call apply_parameter_interval(s,axis,lower,upper)
    end procedure apply_inline_interval

    subroutine apply_parameter_interval(s,axis,lower,upper)
        type(studio_t), intent(inout), target :: s
        integer, intent(in) :: axis
        real(dp), intent(in) :: lower,upper
        character(:), allocatable :: error
        logical :: changed

        call commit_single_patch_interval(s%entity_kind,axis,lower,upper,changed,error,catalog=s%document%patches,patch_id=s%live_scene%variant)
        if(len(error)>0)then
            call show_patch_issue(s,error,5500)
            return
        end if
        if(.not.changed)return
        s%live_scene%valid=.false.
        call select_entity(s,s%entity_kind)
        call s%status%showMessage(parameter_directions(axis:axis)//' parameter interval updated',5000)
    end subroutine apply_parameter_interval

    module procedure run_refinement
        type(QDialog),target::dialog
        type(QVBoxLayout)::layout
        type(QHBoxLayout)::direction_row,inputs
        type(QLabel)::heading,label
        type(QButtonGroup)::directions
        type(QWidget)::direction_choices
        type(refinement_inputs),target::fields
        type(QDialogButtonBox)::buttons
        type(QPushButton)::apply
        character(:),allocatable::error
        integer::axis,kind,variant,i,direction_id
        integer(int64)::revision
        logical::changed

        if(.not.s%live_scene%valid.or.s%selected_multipatch>0)return
        kind=s%entity_kind
        variant=s%live_scene%variant
        axis=max(1,min(kind,s%knot_axis))
        revision=s%document%patches%patch_revision(kind,variant)
        fields%kind=kind
        do i=0,kind
            call single_patch_refinement_limits(kind,i,fields%max_divisions(i),fields%max_steps(i), &
                catalog=s%document%patches,patch_id=variant)
        end do
        dialog=QDialog(s%window)
        fields%dialog=dialog
        call dialog%setObjectName('studioRefineDialog')
        call dialog%setWindowTitle('Refine patch')
        call dialog%setMinimumWidth(478)
        layout=QVBoxLayout(dialog)
        call layout%setContentsMargins(16,14,16,14)
        call layout%setSpacing(8)
        heading=QLabel('Refine patch',dialog)
        call heading%setObjectName('studioDialogTitle')
        call layout%addWidget(heading)
        direction_row=QHBoxLayout()
        call direction_row%setSpacing(10)
        label=QLabel('Direction',dialog)
        call direction_row%addWidget(label)
        direction_choices=parameter_direction_buttons(dialog,kind,axis,'studioRefineDirections', &
            directions,update_refinement_inputs,fields,include_all=.true.)
        call direction_row%addWidget(direction_choices,1)
        call layout%addLayout(direction_row)

        inputs=QHBoxLayout()
        call inputs%setSpacing(8)
        call inputs%addWidget(QLabel('Divisions per span',dialog))
        fields%subdivisions=QSpinBox(dialog)
        call fields%subdivisions%setObjectName('studioRefineSubdivisions')
        call fields%subdivisions%setAccessibleName('Divisions per span')
        call fields%subdivisions%setAccessibleDescription('1 keeps existing spans.')
        call fields%subdivisions%setRange(1,fields%max_divisions(axis))
        call fields%subdivisions%setValue(2)
        call fields%subdivisions%setMinimumWidth(96)
        call inputs%addWidget(fields%subdivisions,1)
        call inputs%addSpacing(6)
        call inputs%addWidget(QLabel('Degree increase',dialog))
        fields%degree=QSpinBox(dialog)
        call fields%degree%setObjectName('studioRefineDegree')
        call fields%degree%setAccessibleName('Degree increase')
        call fields%degree%setAccessibleDescription('0 keeps the current degree. Degree elevation runs before subdivision.')
        call fields%degree%setRange(0,fields%max_steps(axis))
        call fields%degree%setMinimumWidth(96)
        call inputs%addWidget(fields%degree,1)
        call layout%addLayout(inputs)
        call set_status_hint(fields%subdivisions,'1 keeps existing spans.')
        call set_status_hint(fields%degree,'0 keeps the current degree. Degree elevation runs before subdivision.')
        call fields%subdivisions%valueChanged%connect(refinement_value_changed,fields)
        call fields%degree%valueChanged%connect(refinement_value_changed,fields)
        buttons=QDialogButtonBox(QDialogButtonBox_Cancel,dialog)
        apply=buttons%addButton('Refine patch',QDialogButtonBox_AcceptRole)
        call apply%setObjectName('studioRefineApply')
        call studio_style_primary(apply)
        call apply%setDefault(.true.)
        call buttons%accepted%connect(accept_refinement,dialog)
        call buttons%rejected%connect(cancel_refinement,dialog)
        call layout%addWidget(buttons)
        do while(dialog%exec()==QDialog_Accepted)
            if(s%entity_kind/=kind.or.s%selected_multipatch/=0.or.active_patch_id(kind,catalog=s%document%patches)/=variant.or. &
                s%document%patches%patch_revision(kind,variant)/=revision)then
                call s%status%showMessage('The patch changed. Reopen refinement for the current selection.',5500)
                exit
            end if
            direction_id=directions%checkedId()
            axis=merge(0,direction_id+1,direction_id==kind)
            call commit_single_patch_refinement(kind,axis, &
                fields%subdivisions%value(),fields%degree%value(),changed,error,catalog=s%document%patches,patch_id=variant)
            if(changed)then
                s%live_scene%valid=.false.
                call select_entity(s,kind)
                call s%status%showMessage('Patch refined · Undo with Ctrl+Z',5000)
                exit
            end if
            if(len(error)==0)error='Increase subdivisions or degree to refine the patch.'
            call publish_window_status(dialog,error,warning=.true.,blocking=.true.)
        end do
        call dialog%delete()
    end procedure run_refinement

    subroutine update_refinement_inputs(index,data)
        integer, intent(in) :: index
        class(*), optional, intent(inout) :: data
        integer :: axis

        select type(fields=>data)
        type is(refinement_inputs)
            if(index<0.or.index>fields%kind)return
            axis=merge(0,index+1,index==fields%kind)
            call fields%subdivisions%setRange(1,fields%max_divisions(axis))
            call fields%degree%setRange(0,fields%max_steps(axis))
            call publish_window_status(fields%dialog,'',blocking=.true.)
        end select
    end subroutine update_refinement_inputs

    subroutine refinement_value_changed(value,data)
        integer, intent(in) :: value
        class(*), optional, intent(inout) :: data

        select type(fields=>data)
        type is(refinement_inputs)
            call publish_window_status(fields%dialog,'',blocking=.true.)
        end select
    end subroutine refinement_value_changed

    subroutine accept_refinement(data)
        class(*),intent(inout),optional::data

        select type(dialog=>data)
        type is(QDialog)
            call dialog%accept()
        end select
    end subroutine accept_refinement

    subroutine cancel_refinement(data)
        class(*),intent(inout),optional::data

        select type(dialog=>data)
        type is(QDialog)
            call dialog%reject()
        end select
    end subroutine cancel_refinement

    module procedure run_knot_values_editor
        character::direction
        integer::axis

        if(.not.can_remove_knot(s,any_direction=.true.))return
        axis=s%knot_axis
        if(.not.edit_single_patch_knot_values(s%window,s%entity_kind,axis,catalog=s%document%patches,patch_id=s%live_scene%variant))return
        s%live_scene%valid=.false.
        call select_entity(s,s%entity_kind)
        direction=parameter_directions(axis:axis)
        call s%status%showMessage(direction//' knot values updated',5000)
    end procedure run_knot_values_editor

    module procedure on_knot_interaction
        integer::row
        character(:),allocatable :: text

        if(operation=='Select knot tool')then
            call sync_controls(s)
            call refresh_domain_readout(s)
            call s%domain%update()
            return
        end if
        if(len(error)>0)then
            call show_patch_issue(s,error,5500)
            return
        end if
        if(changed)then
            s%live_scene%valid=.false.
            call select_entity(s,s%entity_kind)
        end if
        if(.not.allocated(s%live_knot_groups))return
        row=findloc(s%live_knot_groups,value,dim=1)-1
        if(row>=0)call s%knot_table%setCurrentCell(row,1)
        if(changed)then
            call s%status%showMessage(operation,4000)
        else
            text=knot_status_text(s)
            if(len(text)>0)call s%status%showMessage(text,3800,record=.false.)
        end if
    end procedure on_knot_interaction

    module procedure can_restructure_patch

        allowed=s%live_scene%valid.and.s%selected_multipatch==0
        if(allowed)allowed=.not.document_patch_in_assembly(s%document,s%entity_kind,s%live_scene%variant)
    end procedure can_restructure_patch
end submodule studio_ui_patch_tools
