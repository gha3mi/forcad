module studio_assistant_draft_context
    use studio_knot_interaction, only: knot_selection
    use studio_presentation_text, only: entity_name, integer_text, number_text
    use studio_scene, only: studio_scene_control_position
    use studio_selection, only: selection_target_t, selection_patch, selection_control, selection_assembly, selection_mesh
    use studio_single_patch, only: single_patch_name
    use studio_state, only: studio_t, dp, studio_mesh_count
    use studio_multipatch_collection, only: count_multipatches
    use studio_units, only: length_caption
    use studio_viewport_scene, only: pick_face, pick_edge, pick_vertex, viewport_feature_text
    use studio_workspace_selection, only: primary_workspace_target
    implicit none(type, external)
    private
    public :: capture_assistant_selection

    character(*), parameter :: axes='UVW', line=new_line('a')

contains

    subroutine capture_assistant_selection(s,label,details)
        type(studio_t), intent(in) :: s
        character(:), allocatable, intent(out) :: label,details
        type(selection_target_t), allocatable :: targets(:)
        character(:), allocatable :: caption,entry
        integer :: i,n
        logical :: active_patch

        if(allocated(s%selection%items))then
            targets=s%selection%items
        else
            targets=[primary_workspace_target(s)]
        end if
        label=''
        details=''
        n=0
        active_patch=.false.
        do i=1,size(targets)
            call describe_selection_item(s,targets(i),i==size(targets),caption,entry)
            if(len(entry)==0)cycle
            if(targets(i)%category==selection_patch.or.targets(i)%category==selection_control)then
                active_patch=active_patch.or.(targets(i)%rank==s%entity_kind.and.targets(i)%id==s%live_scene%variant)
            end if
            n=n+1
            label=caption
            details=details//line//'- '//entry
        end do
        call append_parameter_context(s,active_patch,label,details,n)
        if(n==0)return
        if(n>1)label=integer_text(n)//' selected items'
        details='Attached selection at capture time. Confirm targets against the current design before editing.'//details
    end subroutine capture_assistant_selection

    subroutine describe_selection_item(s,target,primary,caption,entry)
        type(studio_t), intent(in) :: s
        type(selection_target_t), intent(in) :: target
        logical, intent(in) :: primary
        character(:), allocatable, intent(out) :: caption,entry
        character(:), allocatable :: name,feature,orientation
        integer :: element,number,orientation_start
        real(dp) :: point(3)

        caption=''
        entry=''
        name=''
        feature=''
        orientation=''
        if(target%id<1)return
        select case(target%category)
        case(selection_patch,selection_control)
            name=single_patch_name(target%rank,target%id,catalog=s%document%patches)
            if(len(name)==0)return
            entry=entity_name(target%rank)//' patch '//integer_text(target%id)
            if(target%category==selection_control)then
                feature='Control point '//integer_text(target%item)
            else
                element=target%feature
                number=target%item
                if(element==0.and.primary.and.all(s%viewport_feature(:2)==[target%rank,target%id]))then
                    element=s%viewport_feature(3)
                    number=s%viewport_feature(4)
                end if
                select case(element)
                case(pick_face)
                    feature='Face '//integer_text(number)
                case(pick_edge)
                    feature='Edge '//integer_text(number)
                case(pick_vertex)
                    feature='Vertex '//integer_text(number)
                    if(target%rank==1)feature='Endpoint '//integer_text(number)
                case default
                    continue
                end select
                if(len(feature)>0.and.all(s%viewport_feature==[target%rank,target%id,element,number]))then
                    orientation=viewport_feature_text(s)
                end if
            end if
        case(selection_assembly)
            if(target%id>count_multipatches(s%document%multipatches))return
            associate(definition=>s%document%multipatches%groups(target%id)%definition)
                if(allocated(definition%name))name=definition%name
            end associate
            entry='Multipatch '//integer_text(target%id)
        case(selection_mesh)
            if(target%id>studio_mesh_count(s))return
            associate(mesh=>s%document%meshes(target%id))
                if(allocated(mesh%name))name=mesh%name
                entry='Mesh '//integer_text(target%id)
                if(target%item>0)then
                    feature='Cell '//integer_text(target%item)
                    if(target%feature>0)feature=feature//' · Face '//integer_text(target%feature)
                    entry=entry//', part '//integer_text(target%part)
                    if(allocated(mesh%parts))then
                        if(target%part>=1.and.target%part<=size(mesh%parts))then
                            entry=entry//' ('//entity_name(mesh%kind)//' patch '//integer_text(mesh%parts(target%part)%variant)//')'
                        end if
                    end if
                end if
            end associate
        case default
            return
        end select
        caption=entry
        if(len(name)>0)then
            entry=entry//' · '//name
            caption=name
        end if
        if(len(feature)>0)then
            entry=entry//' · '//feature
            caption=caption//' · '//feature
        end if
        orientation_start=index(orientation,' · ')
        if(orientation_start>0)entry=entry//' ('//orientation(orientation_start+len(' · '):)//')'
        if(target%category==selection_control)then
            if(studio_scene_control_position(s%document%scene,target%rank,target%id,target%item,point,s%document%patches))then
                entry=entry//' · '//length_caption('Position',s%document%model_unit)//' '//coordinate_tuple(point)
            end if
        end if
    end subroutine describe_selection_item

    subroutine append_parameter_context(s,active_patch,label,details,n)
        type(studio_t), intent(in) :: s
        logical, intent(in) :: active_patch
        character(:), allocatable, intent(inout) :: label,details
        integer, intent(inout) :: n
        character(:), allocatable :: patch,parameters,span
        real(dp) :: knot
        integer :: i,kind,axis,multiplicity
        logical :: selected_knot,selected_span

        if(.not.s%live_scene%valid.or.s%selected_mesh>0.or.s%selected_multipatch>0)return
        if(s%domain%isNull())return
        if(.not.s%domain%isVisible())return
        kind=s%live_scene%dimension
        if(kind<1.or.kind>len(axes))return
        patch=entity_name(kind)//' patch '//integer_text(s%live_scene%variant)
        selected_knot=.false.
        if(s%domain_layers(2).or.s%domain_edit_knots)selected_knot=knot_selection(s,knot)
        if(selected_knot)then
            axis=s%knot_axis
            multiplicity=count(s%live_knots==knot)
            label=axes(axis:axis)//' knot '//number_text(knot)
            details=details//line//'- '//patch//' · '//label//' · Multiplicity '//integer_text(multiplicity)
            n=n+1
        end if

        if(.not.s%live_probe%valid.or.s%live_probe%dimension/=kind)return
        selected_span=s%domain_layers(5)
        if(.not.(active_patch.or.selected_knot.or.selected_span))return
        parameters=''
        span=''
        do i=1,kind
            if(i>1)then
                parameters=parameters//', '
                span=span//', '
            end if
            parameters=parameters//axes(i:i)//'='//number_text(s%live_probe%parameters(i))
            if(selected_span)span=span//axes(i:i)//' '//coordinate_tuple(s%live_probe%cell(:,i))
        end do
        if(selected_span)then
            label=trim(merge('Parameter span   ','Parameter element',kind==1))
            details=details//line//'- '//patch//' · '//label//' · '//span
            n=n+1
        end if
        details=details//line//'- '//patch//' · Probe '//parameters//' · '// &
            length_caption('Position',s%document%model_unit)//' '//coordinate_tuple(s%live_probe%point)
    end subroutine append_parameter_context

    function coordinate_tuple(values) result(text)
        real(dp), intent(in) :: values(:)
        character(:), allocatable :: text
        integer :: i

        text='('
        do i=1,size(values)
            if(i>1)text=text//', '
            text=text//number_text(values(i))
        end do
        text=text//')'
    end function coordinate_tuple
end module studio_assistant_draft_context
