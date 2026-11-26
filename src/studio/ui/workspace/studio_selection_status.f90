module studio_selection_status
    use studio_document, only: document_has_patch
    use studio_knot_interaction, only: knot_selection
    use studio_multipatch_collection, only: count_multipatches
    use studio_mesh, only: studio_mesh_t
    use studio_presentation_text, only: entity_name,integer_text,number_text
    use studio_scene, only: studio_scene_is_visible
    use studio_selection, only: selection_target_t,selection_patch,selection_control,selection_assembly,selection_mesh, &
        pick_face,pick_edge,pick_vertex
    use studio_single_patch, only: single_patch_controls,read_single_patch_controls,single_patch_name
    use studio_state, only: studio_t,dp,studio_mesh_count,studio_mesh_valid,studio_mesh_cell_label
    use studio_units, only: display_length,unit_symbol,unitless
    use studio_viewport_scene, only: viewport_feature_text,corner_sample
    use studio_workspace_selection, only: primary_workspace_target
    implicit none(type,external)
    private

    public :: selection_status_text,probe_status_text,knot_status_text,mesh_counts_text

    character(*),parameter :: separator=' · ',axes='UVW'

contains

    function selection_status_text(s,compact) result(text)
        type(studio_t),intent(in) :: s
        logical,optional,intent(in) :: compact
        type(selection_target_t) :: target
        character(:),allocatable :: text,summary
        integer :: count,hidden,i
        logical :: brief

        text=''
        if(s%home_active)return
        brief=.false.
        if(present(compact))brief=compact
        count=1
        if(allocated(s%selection%items))then
            count=s%selection%count()
            if(count==0)return
            target=s%selection%items(count)
        else
            target=primary_workspace_target(s)
        end if

        if(count==1)then
            text=selected_item_text(s,target,brief)
            if(len(text)>0.and.target_hidden(s,target))text=text//separator//'Hidden'
            return
        end if
        text=integer_text(count)//' selected'
        if(.not.brief)then
            if(all(s%selection%items%category==selection_control))text=integer_text(count)//' control points selected'
            if(all(s%selection%items%category==selection_mesh.and.s%selection%items%item>0))text=integer_text(count)//' mesh cells selected'
        end if
        hidden=0
        do i=1,count
            if(target_hidden(s,s%selection%items(i)))hidden=hidden+1
        end do
        if(hidden>0)text=text//separator//integer_text(hidden)//' hidden'
        if(brief)return
        summary=selected_item_text(s,target,brief)
        if(len(summary)>0)text=text//separator//'Active: '//summary
    end function selection_status_text

    logical function target_hidden(s,target) result(hidden)
        type(studio_t),intent(in) :: s
        type(selection_target_t),intent(in) :: target
        integer :: i,member

        hidden=.false.
        if(target%id<1)return
        select case(target%category)
        case(selection_patch,selection_control)
            if(.not.document_has_patch(s%document,target%rank,target%id))return
            hidden=.not.studio_scene_is_visible(s%document%scene,target%rank,target%id,catalog=s%document%patches)
        case(selection_assembly)
            if(target%id>count_multipatches(s%document%multipatches))return
            associate(definition=>s%document%multipatches%groups(target%id)%definition)
                if(.not.allocated(definition%members))return
                do i=1,size(definition%members)
                    member=definition%members(i)
                    if(.not.document_has_patch(s%document,definition%kind,member))cycle
                    hidden=.not.studio_scene_is_visible(s%document%scene,definition%kind,member,catalog=s%document%patches)
                    if(.not.hidden)return
                end do
            end associate
        case(selection_mesh)
            if(target%id>studio_mesh_count(s))return
            hidden=.not.s%document%meshes(target%id)%visible
        case default
            continue
        end select
    end function target_hidden

    function selected_item_text(s,target,brief) result(text)
        type(studio_t),intent(in) :: s
        type(selection_target_t),intent(in) :: target
        logical,intent(in) :: brief
        character(:),allocatable :: text,name

        text=''
        if(target%id<1)return
        select case(target%category)
        case(selection_patch,selection_control)
            if(.not.document_has_patch(s%document,target%rank,target%id))return
            name=single_patch_name(target%rank,target%id,catalog=s%document%patches)
            if(target%category==selection_control)then
                text=control_status(s,target,name,brief)
            else
                text=patch_status(s,target,name,brief)
            end if
        case(selection_assembly)
            if(target%id>count_multipatches(s%document%multipatches))return
            text=multipatch_status(s,target%id,brief)
        case(selection_mesh)
            if(target%id>studio_mesh_count(s))return
            text=mesh_status(s,target,brief)
        case default
            continue
        end select
    end function selected_item_text

    function control_status(s,target,name,brief) result(text)
        type(studio_t),intent(in) :: s
        type(selection_target_t),intent(in) :: target
        character(*),intent(in) :: name
        logical,intent(in) :: brief
        character(:),allocatable :: text
        type(single_patch_controls) :: controls
        real(dp) :: point(3)

        text=''
        if(target%item<1)return
        text='Control point '//integer_text(target%item)
        if(brief)return
        call read_single_patch_controls(target%rank,controls,variant=target%id,catalog=s%document%patches)
        if(.not.controls%valid.or.target%item>size(controls%weights))then
            text=text//separator//'Control data unavailable'
            return
        end if
        point=controls%points(:,target%item)
        if(s%control_dragging.and.s%selected_control==target%item.and. &
            s%entity_kind==target%rank.and.s%live_scene%variant==target%id)point=s%control_preview
        text=text//separator//position_text(s,point)//separator//'Weight '//number_text(controls%weights(target%item))//separator//name
    end function control_status

    function patch_status(s,target,name,brief) result(text)
        type(studio_t),intent(in) :: s
        type(selection_target_t),intent(in) :: target
        character(*),intent(in) :: name
        logical,intent(in) :: brief
        character(:),allocatable :: text,kind,feature,orientation,degrees
        integer :: element,item,i,node,orientation_start
        logical :: active

        active=s%live_scene%valid.and.s%live_scene%dimension==target%rank.and.s%live_scene%variant==target%id
        kind=entity_name(target%rank)
        if(active.and.target%rank==1.and.s%live_scene%control_counts(1)==1)kind='Point'
        element=target%feature
        item=target%item
        if(element==0.and.active.and.all(s%viewport_feature(:2)==[target%rank,target%id]))then
            element=s%viewport_feature(3)
            item=s%viewport_feature(4)
        end if
        feature=''
        if(item>0.and.kind/='Point')then
            select case(element)
            case(pick_face); feature='Face'
            case(pick_edge); feature='Edge'
            case(pick_vertex)
                feature='Vertex'
                if(target%rank==1)feature='Endpoint'
            case default; continue
            end select
            if(len(feature)>0)feature=feature//' '//integer_text(item)
        end if

        text=name
        if(brief)then
            if(len(feature)>0)text=feature//separator//name
            return
        end if
        if(len(feature)>0)then
            text=feature
            if(all(s%viewport_feature==[target%rank,target%id,element,item]))then
                orientation=viewport_feature_text(s)
                orientation_start=index(orientation,separator)
                if(orientation_start>0)text=text//' ('//orientation(orientation_start+len(separator):)//')'
            end if
            if(active.and.element==pick_vertex.and.item<=2**target%rank.and.allocated(s%live_scene%points))then
                node=corner_sample(s%live_scene%sample_counts,target%rank,item)
                if(node<=size(s%live_scene%points,2))text=text//separator//position_text(s,s%live_scene%points(:,node))
            end if
            text=text//separator//kind//separator//name
        else
            text=kind//separator//name
            if(active.and.kind=='Point'.and.allocated(s%live_scene%controls))then
                if(size(s%live_scene%controls,2)==1)text=text//separator//position_text(s,s%live_scene%controls(:,1))
            else if(active)then
                degrees=integer_text(s%live_scene%degree(1))
                do i=2,target%rank
                    degrees=degrees//'×'//integer_text(s%live_scene%degree(i))
                end do
                text=text//separator//'Degree '//degrees//separator//integer_text(product(s%live_scene%control_counts))//' control points'
            end if
        end if
    end function patch_status

    function multipatch_status(s,id,brief) result(text)
        type(studio_t),intent(in) :: s
        integer,intent(in) :: id
        logical,intent(in) :: brief
        character(:),allocatable :: text
        integer :: members,interfaces,point,hover,selected

        associate(group=>s%document%multipatches%groups(id),definition=>s%document%multipatches%groups(id)%definition)
            text='Multipatch'
            if(allocated(definition%name))then
                if(brief)then
                    text=definition%name
                else
                    text=text//separator//definition%name
                end if
            end if
            if(brief)return
            point=s%selected_control
            hover=s%canvas_gesture%viewport_shared_hover
            if(hover>0.and..not.s%control_dragging)point=hover
            if(s%selected_multipatch==id.and.s%shared_controls%valid.and.point>0)then
                if(point<=size(s%shared_controls%weights))then
                    text='Control '//integer_text(point)
                    if(s%shared_controls%copies(point)>1)text='Shared control '//integer_text(point)
                    text=text//separator//shared_control_owners(s,point)
                    selected=0
                    if(allocated(s%selected_shared_controls))selected=size(s%selected_shared_controls)
                    if(selected>1)text=integer_text(selected)//' controls selected'//separator//text
                    text=text//separator//position_text(s,s%shared_controls%points(:,point))// &
                        separator//'Weight '//number_text(s%shared_controls%weights(point))
                    return
                end if
            end if
            members=0
            if(allocated(definition%members))members=size(definition%members)
            interfaces=0
            if(allocated(definition%joins))interfaces=size(definition%joins)
            text=text//separator//integer_text(members)//' patches'//separator//integer_text(interfaces)//' interface'
            if(interfaces/=1)text=text//'s'
            if(group%validation%valid)then
                if(interfaces>0)text=text//separator//'Interfaces valid'
            else
                text=text//separator//'Needs review'
            end if
        end associate
    end function multipatch_status

    function shared_control_owners(s,point) result(text)
        type(studio_t),intent(in) :: s
        integer,intent(in) :: point
        character(:),allocatable :: text,name
        integer :: i,owners

        text=''
        owners=0
        associate(controls=>s%shared_controls)
            do i=1,size(controls%members)
                if(.not.any(controls%map(controls%offsets(i)+1:controls%offsets(i+1))==point))cycle
                owners=owners+1
                if(owners>2)cycle
                name=single_patch_name(controls%kind,controls%members(i),catalog=s%document%patches)
                if(owners==2)text=text//' / '
                text=text//name
            end do
        end associate
        if(owners>2)text=integer_text(owners)//' linked patches'
    end function shared_control_owners

    function mesh_status(s,target,brief) result(text)
        type(studio_t),intent(in) :: s
        type(selection_target_t),intent(in) :: target
        logical,intent(in) :: brief
        character(:),allocatable :: text,cell,name

        cell=studio_mesh_cell_label(s,target)
        name='Mesh '//integer_text(target%id)
        associate(mesh=>s%document%meshes(target%id))
            if(allocated(mesh%name))name=mesh%name
            if(len(cell)>0)then
                if(target%feature>0)then
                    cell=cell//separator//'Face '//integer_text(target%feature)
                else if(mesh%kind==3.and.s%mesh_section_axis>0)then
                    cell=cell//separator//'Cut face'
                end if
                text=cell//separator//name
            else
                text=name
                if(.not.brief)text='Mesh'//separator//name
            end if
            if(brief)return
            text=text//separator//mesh_counts_text(mesh)
            if(.not.studio_mesh_valid(s,target%id))text=text//separator//'Needs update'
        end associate
    end function mesh_status

    function mesh_counts_text(mesh) result(text)
        type(studio_mesh_t),intent(in) :: mesh
        character(:),allocatable :: text

        text=integer_text(mesh%elements)//' cells'//separator//integer_text(mesh%nodes)//' nodes'
        if(allocated(mesh%connected_cells))text=text//separator//'Conforming hex mesh'
    end function mesh_counts_text

    function probe_status_text(s) result(text)
        type(studio_t),intent(in) :: s
        character(:),allocatable :: text
        integer :: i

        text=''
        if(s%home_active.or.s%selected_mesh>0.or.s%selected_multipatch>0.or..not.s%live_scene%valid.or..not.s%live_probe%valid)return
        if(s%live_probe%dimension/=s%live_scene%dimension.or.s%live_probe%dimension<1.or.s%live_probe%dimension>len(axes))return
        text='Probe'
        do i=1,s%live_probe%dimension
            text=text//separator//axes(i:i)//' '//number_text(s%live_probe%parameters(i))
        end do
        text=text//separator//position_text(s,s%live_probe%point)
        if(allocated(s%live_probe%name))text=text//separator//s%live_probe%name
    end function probe_status_text

    function knot_status_text(s) result(text)
        type(studio_t),intent(in) :: s
        character(:),allocatable :: text
        integer :: axis
        real(dp) :: value

        text=''
        axis=s%knot_axis
        if(s%home_active.or.s%selected_mesh>0.or.s%selected_multipatch>0.or.axis<1.or.axis>len(axes))return
        if(.not.knot_selection(s,value))return
        text=axes(axis:axis)//' knot '//number_text(value)//separator//'Multiplicity '//integer_text(count(s%live_knots==value))
        text=text//separator//single_patch_name(s%live_scene%dimension,s%live_scene%variant,catalog=s%document%patches)
    end function knot_status_text

    function position_text(s,point) result(text)
        type(studio_t),intent(in) :: s
        real(dp),intent(in) :: point(3)
        character(:),allocatable :: text,units
        real(dp) :: coordinate(3)
        integer :: unit

        unit=s%document%display_unit
        if(unit==unitless)unit=s%document%model_unit
        units=trim(unit_symbol(unit))
        coordinate=display_length(point,s%document%model_unit,unit)
        text='XYZ ('//number_text(coordinate(1))//', '//number_text(coordinate(2))//', '//number_text(coordinate(3))//')'
        if(len(units)>0)text=text//' '//units
    end function position_text

end module studio_selection_status
