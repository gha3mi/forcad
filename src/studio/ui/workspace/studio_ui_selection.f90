submodule(studio_ui) studio_ui_selection
    use qt, only: QScrollBar, QCoreApplication_sendEvent, QEvent_KeyPress, Qt_ControlModifier, Qt_Key_A
    use studio_selection, only: change_selection_batch
    use studio_workspace_widgets, only: project_tree_view
    implicit none(type, external)
contains

    module procedure forget_workspace_selection
        if(s%applying_selection)return
        if(.not.allocated(s%selection%items))return
        deallocate(s%selection%items)
        if(.not.s%widgets%project_browser%tree%isNull())call s%widgets%project_browser%tree%clearSelection()
    end procedure forget_workspace_selection

    module procedure apply_workspace_selection
        call apply_workspace_selection_batch(s,[target],mode)
    end procedure apply_workspace_selection

    module procedure apply_workspace_selection_batch
        type(selection_target_t) :: focused
        integer :: n
        logical :: mesh_pick

        if(size(targets)==0.and.mode/=selection_replace)return
        if(mode/=selection_replace.and..not.allocated(s%selection%items))then
            focused=primary_workspace_target(s)
            if(targets(1)%category/=selection_control.or.focused%category==selection_control)then
                call s%selection%change(focused,selection_replace)
            end if
        end if
        call change_selection_batch(s%selection,targets,mode)
        n=s%selection%count()
        if(n>1.and.s%transform_mode>0)call on_transform_exit(s)
        if(n==0)then
            call clear_workspace_selection(s)
            return
        end if

        focused=s%selection%items(n)
        mesh_pick=s%mesh_pick
        s%applying_selection=.true.
        select case(focused%category)
        case(selection_patch,selection_control)
            if(s%selected_mesh>0.or.s%selected_multipatch>0.or..not.s%live_scene%valid.or. &
                s%entity_kind/=focused%rank.or.s%live_scene%variant/=focused%id)then
                call select_patch(s,focused%rank,focused%id,reveal=.false.)
            end if
            s%viewport_feature=0
            if(focused%category==selection_patch.and.focused%feature>0) then
                s%viewport_feature=[focused%rank,focused%id,focused%feature,focused%item]
            end if
            s%selected_control=0
            if(focused%category==selection_control)s%selected_control=focused%item
            call on_control_selection(s)
        case(selection_assembly)
            if(s%selected_multipatch/=focused%id)call select_multipatch_group(s,focused%id,reveal=.false.)
        case(selection_mesh)
            if(focused%item>0)then
                call on_mesh_cell_selected(s,focused%id,focused%part,focused%item,focused%feature)
            else if(s%selected_mesh/=focused%id)then
                call select_project_mesh(s,focused%id)
            else
                s%mesh_cell=0
                call update_mesh_ribbon(s)
            end if
        case default
            continue
        end select
        if(s%mesh_pick.neqv.mesh_pick)then
            s%mesh_pick=mesh_pick
            call sync_controls(s)
        end if
        s%applying_selection=.false.
        call sync_selection_rows(s)
        call sync_patch_visibility_button(s)
        call refresh_edit_actions(s)
        call refresh_workspace_status(s)
        call s%surface%update()
    end procedure apply_workspace_selection_batch

    module procedure clear_workspace_selection
        logical :: refresh_context

        refresh_context=s%selected_mesh>0.or.s%selected_multipatch>0
        call s%selection%clear()
        if(s%transform_mode>0)call on_transform_exit(s)
        s%viewport_feature=0
        s%selected_control=0
        s%control_dragging=.false.
        s%control_snapped=.false.
        s%mesh_cell=0
        s%selected_mesh=0
        s%selected_multipatch=0
        if(refresh_context)then
            call refresh_inspector(s)
            call refresh_knot_table(s)
            call sync_controls(s)
        else
            call on_control_selection(s)
            call refresh_edit_actions(s)
        end if
        call sync_selection_rows(s)
        call sync_patch_visibility_button(s)
        call refresh_workspace_status(s)
        call s%status%showMessage('Selection cleared',2300)
        call s%surface%update()
    end procedure clear_workspace_selection

    module procedure tree_selection_target
        integer :: rank,id,assembly,branch

        id=mesh_project_id(s%widgets%project_browser,item)
        if(id>0)then
            target=selection_target_t(selection_mesh,id)
            return
        end if
        call project_item_identity(s%widgets%project_browser,s,item,rank,id,assembly,branch)
        if(assembly>0)then
            target=selection_target_t(selection_assembly,assembly)
        else if(rank>0)then
            target=selection_target_t(selection_patch,id,rank)
        else
            target=selection_target_t()
        end if
    end procedure tree_selection_target

    module procedure sync_selection_rows
        type(QTreeWidgetItem) :: item
        integer :: offset(0:3),i,n
        logical, allocatable :: patches(:),assemblies(:),meshes(:),selected_rows(:)
        logical :: updating

        if(.not.allocated(s%selection%items).or.s%widgets%project_browser%tree%isNull())return
        if(s%selection%count()>1.and..not.s%widgets%project_browser%tree%isVisible())return
        updating=s%updating
        s%updating=.true.
        if(s%selection%count()==0)then
            call s%widgets%project_browser%tree%clearSelection()
            s%updating=updating
            return
        end if

        offset(0)=0
        do i=1,3
            offset(i)=offset(i-1)+single_patch_count(i,catalog=s%document%patches)
        end do
        allocate(patches(offset(3)),assemblies(count_multipatches(s%document%multipatches)), &
            meshes(studio_mesh_count(s)),source=.false.)
        do i=1,s%selection%count()
            associate(target=>s%selection%items(i))
                if(target%id<1)cycle
                select case(target%category)
                case(selection_patch,selection_control)
                    if(target%rank<1.or.target%rank>3)cycle
                    n=offset(target%rank-1)+target%id
                    if(n<=offset(target%rank))patches(n)=.true.
                case(selection_assembly)
                    if(target%id<=size(assemblies))assemblies(target%id)=.true.
                case(selection_mesh)
                    if(target%id<=size(meshes))meshes(target%id)=.true.
                case default
                    continue
                end select
            end associate
        end do
        associate(browser=>s%widgets%project_browser)
            allocate(selected_rows(browser%row_count),source=.false.)
            do i=1,browser%row_count
                associate(row=>browser%rows(i))
                    select case(row%category)
                    case(selection_patch)
                        n=offset(row%rank-1)+row%id
                        if(n<=offset(row%rank))selected_rows(i)=patches(n)
                    case(selection_assembly)
                        if(row%id<=size(assemblies))selected_rows(i)=assemblies(row%id)
                    case(selection_mesh)
                        if(row%id<=size(meshes))selected_rows(i)=meshes(row%id)
                    case default
                        cycle
                    end select
                end associate
            end do
            call browser%tree%clearSelection()
            if(count(selected_rows)>sqrt(real(browser%row_count,dp)).and.browser%tree%isVisible())then
                if(select_visible_scene_rows(browser,selected_rows))then
                    s%updating=updating
                    return
                end if
            end if
            do i=1,browser%row_count
                if(.not.selected_rows(i))cycle
                item=browser%rows(i)%item
                call item%setSelected(.true.)
            end do
        end associate
        s%updating=updating
    end procedure sync_selection_rows

    logical function select_visible_scene_rows(browser,selected) result(handled)
        type(project_tree_view), intent(inout) :: browser
        logical, intent(in) :: selected(:)
        type(QTreeWidgetItem) :: item,current
        type(QKeyEvent) :: select_all
        type(QScrollBar) :: scrollbar
        logical, allocatable :: concealed(:)
        integer :: i,scroll

        allocate(concealed(size(selected)),source=.false.)
        current=browser%tree%currentItem()
        scrollbar=browser%tree%verticalScrollBar()
        scroll=scrollbar%value()
        do i=1,size(selected)
            if(selected(i).or.browser%rows(i)%category==0)cycle
            item=browser%rows(i)%item
            if(item%isHidden().or.c_associated(item%ptr,current%ptr))cycle
            if(item%childCount()>0.and.item%isExpanded())cycle
            call item%setHidden(.true.)
            concealed(i)=.true.
        end do
        select_all=QKeyEvent(QEvent_KeyPress,Qt_Key_A,Qt_ControlModifier)
        handled=QCoreApplication_sendEvent(browser%tree,select_all)
        handled=handled.and.select_all%isAccepted()
        call select_all%delete()
        do i=1,size(selected)
            if(selected(i))cycle
            item=browser%rows(i)%item
            if(concealed(i))call item%setHidden(.false.)
            if(item%isSelected())call item%setSelected(.false.)
        end do
        call scrollbar%setValue(scroll)
    end function select_visible_scene_rows

end submodule studio_ui_selection
