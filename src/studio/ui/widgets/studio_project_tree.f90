module studio_project_tree
    use studio_status_hints, only: set_status_hint, prepare_tree_item_hints
    use, intrinsic :: iso_c_binding, only: c_associated
    use, intrinsic :: iso_fortran_env, only: int64
    use qt, only: QAbstractItemView_ExtendedSelection, QAbstractItemView_SelectRows, QHeaderView, &
        QHeaderView_Fixed, QHeaderView_Stretch, QIcon, QScrollBar, QTreeWidget, QTreeWidgetItem, &
        QVBoxLayout, Qt_ItemIsEnabled, Qt_ItemIsSelectable, Qt_ScrollBarAlwaysOff
    use studio_document, only: document_has_patch
    use studio_icons, only: line_icon, icon_mesh, icon_refresh, icon_visibility_mixed, icon_scene_reference
    use studio_multipatch_collection, only: count_multipatches
    use studio_presentation_text, only: entity_name, integer_text
    use studio_scene, only: studio_scene_is_visible
    use studio_selection, only: selection_target_t, selection_patch, selection_control, selection_assembly, selection_mesh
    use studio_single_patch, only: single_patch_available, single_patch_count, single_patch_name, active_patch_id
    use studio_state, only: studio_t, studio_mesh_count, studio_mesh_valid, scene_drawn_patch_count, scene_mesh_drawn
    use studio_workspace_widgets, only: project_tree_view, project_tree_row, project_tree_state
    implicit none(type, external)
    private

    public :: project_tree_view
    public :: create_project_tree_widgets, build_project_rows, update_project_mesh_rows, mesh_project_id
    public :: sync_scene_visibility_icons, sync_scene_visibility_batch, select_project_patch_row
    public :: update_project_summary, project_item_identity, project_target_row
    public :: apply_project_filter, reveal_project_row, project_row_visibility

    integer, parameter, public :: scene_visibility_column=1
    integer, parameter :: visibility_column_width=27
    integer, parameter, public :: scene_filter_all=0, scene_filter_geometry=1, scene_filter_multipatches=2
    integer, parameter, public :: scene_filter_meshes=3, scene_filter_hidden=4, scene_filter_stale=5
    integer, parameter, public :: scene_filter_curves=6, scene_filter_surfaces=7, scene_filter_volumes=8

contains

    subroutine create_project_tree_widgets(view,layout)
        type(project_tree_view), intent(inout) :: view
        type(QVBoxLayout), intent(inout) :: layout
        type(QHeaderView) :: header
        integer :: rank
        integer, parameter :: icons(3)=[23,24,8]

        view%tree=QTreeWidget()
        call view%tree%setObjectName('studioProjectTree')
        call view%tree%setColumnCount(2)
        call view%tree%setHeaderLabels(['Scene object','Visibility  '])
        call prepare_tree_item_hints(view%tree)
        header=view%tree%header()
        call header%setStretchLastSection(.false.)
        call header%setMinimumSectionSize(visibility_column_width)
        call header%setSectionResizeMode(0,QHeaderView_Stretch)
        call header%setSectionResizeMode(scene_visibility_column,QHeaderView_Fixed)
        call header%resizeSection(scene_visibility_column,visibility_column_width)
        call view%tree%setHeaderHidden(.true.)
        call view%tree%setHorizontalScrollBarPolicy(Qt_ScrollBarAlwaysOff)
        call view%tree%setSelectionMode(QAbstractItemView_ExtendedSelection)
        call view%tree%setSelectionBehavior(QAbstractItemView_SelectRows)
        call view%tree%setAccessibleName('Scene objects and visibility')
        call set_status_hint(view%tree,'Click to select · Ctrl toggles · Shift adds · Esc clears · Right-click for commands')
        call view%tree%setIconSize(18,18)
        call view%tree%setProperty('uniformRowHeights',.true.)
        call view%tree%setProperty('expandsOnDoubleClick',.false.)
        view%viewport=view%tree%viewport()
        do rank=1,3
            view%geometry_icons(rank)=line_icon(icons(rank))
        end do
        view%visible_icon=line_icon(31)
        view%hidden_icon=line_icon(32)
        view%mixed_icon=line_icon(icon_visibility_mixed)
        view%reference_icon=line_icon(icon_scene_reference)
        view%assembly_icon=line_icon(87)
        view%folder_icon=line_icon(2)
        view%mesh_row_icon=line_icon(icon_mesh)
        view%stale_mesh_icon=line_icon(icon_refresh)
        call layout%addWidget(view%tree,1)
    end subroutine create_project_tree_widgets


    subroutine build_project_rows(view,s,restored_current)
        type(project_tree_view), intent(inout) :: view
        type(studio_t), intent(in) :: s
        logical, optional, intent(out) :: restored_current
        type(project_tree_state), allocatable :: previous(:)
        type(QTreeWidgetItem) :: branch,item
        type(QScrollBar) :: scrollbar
        integer :: rank,slot,id,scroll
        logical :: same_document,current_restored
        character(len=8), parameter :: section_names(3)=[character(len=8)::'Curves','Surfaces','Volumes']

        same_document=view%document_revision==s%document%document_revision
        scrollbar=view%tree%verticalScrollBar()
        scroll=0
        if(same_document)then
            previous=capture_project_tree(view)
            scroll=scrollbar%value()
        else
            view%document_revision=s%document%document_revision
            call reset_project_filter(view)
        end if
        call view%tree%clear()
        view%stamp=view%stamp+1_int64
        if(allocated(view%rows))deallocate(view%rows)
        allocate(view%rows(24))
        view%row_count=0
        view%visibility_dirty=.true.

        view%geometry_branch=QTreeWidgetItem(view%tree,['Geometry'])
        call view%geometry_branch%setIcon(0,view%folder_icon)
        call view%geometry_branch%setFlags(Qt_ItemIsEnabled)
        call register_project_row(view,view%geometry_branch,'geometry',branch=-4)
        do rank=1,3
            branch=QTreeWidgetItem(view%geometry_branch,[trim(section_names(rank))])
            call branch%setFlags(Qt_ItemIsEnabled)
            call branch%setIcon(0,view%geometry_icons(rank))
            call register_project_row(view,branch,'rank:'//integer_text(rank),branch=rank)
            view%entity_branches(rank)=branch
            do slot=1,single_patch_count(rank,catalog=s%document%patches)
                id=slot
                if(slot<=2)id=3-slot
                if(.not.document_has_patch(s%document,rank,id))cycle
                item=QTreeWidgetItem(branch,[project_patch_name(s,rank,id)])
                call register_project_row(view,item,patch_row_key(rank,id),category=selection_patch,rank=rank,id=id)
                call sync_patch_tree_item(view,item,s,rank,id)
            end do
            call branch%setExpanded(.true.)
        end do
        call view%geometry_branch%setExpanded(.true.)

        call build_multipatch_rows(view,s)
        view%mesh_branch=QTreeWidgetItem(view%tree,['Meshes'])
        call view%mesh_branch%setFlags(Qt_ItemIsEnabled)
        call view%mesh_branch%setIcon(0,view%mesh_row_icon)
        call set_status_hint(view%tree,view%mesh_branch,0,'Generate a mesh from a patch or multipatch in the Mesh tab.')
        call register_project_row(view,view%mesh_branch,'meshes',branch=-3)
        call view%mesh_branch%setExpanded(.true.)
        call fill_project_mesh_rows(view,s)
        current_restored=.false.
        if(same_document)call restore_project_tree(view,previous,current_restored, &
            expansion_only=allocated(s%selection%items).and.s%selection%count()==0)
        call apply_project_filter(view,s)
        call scrollbar%setValue(scroll)
        if(present(restored_current))restored_current=current_restored
    end subroutine build_project_rows

    subroutine build_multipatch_rows(view,s)
        type(project_tree_view), intent(inout) :: view
        type(studio_t), intent(in) :: s
        type(QTreeWidgetItem) :: group,item
        character(:), allocatable :: name,key
        integer :: id,member,rank,variant

        view%multipatch_branch=QTreeWidgetItem(view%tree,['Multipatches'])
        call view%multipatch_branch%setFlags(Qt_ItemIsEnabled)
        call view%multipatch_branch%setIcon(0,view%assembly_icon)
        call set_status_hint(view%tree,view%multipatch_branch,0,'Use Multipatch → Create to connect compatible patches.')
        call register_project_row(view,view%multipatch_branch,'multipatches',branch=-1)
        do id=1,count_multipatches(s%document%multipatches)
            associate(definition=>s%document%multipatches%groups(id)%definition)
                rank=definition%kind
                group=QTreeWidgetItem(view%multipatch_branch,[definition%name])
                call group%setIcon(0,view%assembly_icon)
                call set_status_hint(view%tree,group,0,definition%name//' · '//entity_name(rank)//' multipatch')
                key='multipatch:'//definition%name
                call register_project_row(view,group,key,category=selection_assembly,id=id)
                do member=1,size(definition%members)
                    variant=definition%members(member)
                    if(.not.single_patch_available(rank,variant,catalog=s%document%patches))cycle
                    name=project_patch_name(s,rank,variant)
                    item=QTreeWidgetItem(group,[name])
                    call register_project_row(view,item,key//'/'//patch_row_key(rank,variant), &
                        category=selection_patch,rank=rank,id=variant,assembly=id)
                    call sync_patch_tree_item(view,item,s,rank,variant)
                end do
                call group%setExpanded(id==s%selected_multipatch)
            end associate
        end do
        call view%multipatch_branch%setExpanded(.true.)
    end subroutine build_multipatch_rows

    subroutine update_project_mesh_rows(view,s)
        type(project_tree_view), intent(inout) :: view
        type(studio_t), intent(in) :: s
        type(QTreeWidgetItem) :: item
        integer :: id,n,existing,slot
        logical :: changed

        if(view%mesh_branch%isNull())return
        n=studio_mesh_count(s)
        existing=view%mesh_branch%childCount()
        if(existing>n)then
            call build_project_rows(view,s)
            return
        end if
        do id=1,existing
            item=view%mesh_branch%child(id-1)
            slot=project_row_index(view,item)
            if(slot>0)then
                if(view%rows(slot)%key=='mesh:'//s%document%meshes(id)%name)cycle
            end if
            call build_project_rows(view,s)
            return
        end do
        call fill_project_mesh_rows(view,s,changed)
        if(changed)call apply_project_filter(view,s)
    end subroutine update_project_mesh_rows

    subroutine fill_project_mesh_rows(view,s,changed)
        type(project_tree_view), intent(inout) :: view
        type(studio_t), intent(in) :: s
        logical, optional, intent(out) :: changed
        type(QTreeWidgetItem) :: item
        integer :: id,slot,existing,state
        logical :: new_row,stale,validity_changed

        if(present(changed))changed=.false.
        existing=view%mesh_branch%childCount()
        do id=1,studio_mesh_count(s)
            associate(mesh=>s%document%meshes(id))
                new_row=id>existing
                if(new_row)then
                    item=QTreeWidgetItem(view%mesh_branch,[mesh%name])
                    call register_project_row(view,item,'mesh:'//mesh%name,category=selection_mesh,id=id)
                    call item%setIcon(0,view%mesh_row_icon)
                else
                    item=view%mesh_branch%child(id-1)
                end if
                slot=project_row_index(view,item)
                stale=.not.studio_mesh_valid(s,id)
                validity_changed=new_row.or.(view%rows(slot)%stale.neqv.stale)
                view%rows(slot)%stale=stale
                if(stale)then
                    if(.not.validity_changed)cycle
                    call item%setText(0,mesh%name//' · Needs update')
                    call set_status_hint(view%tree,item,0,mesh%name//' · Source geometry changed or is unavailable. Regenerate from Mesh.')
                    call item%setData(scene_visibility_column,-1)
                    call item%setIcon(scene_visibility_column,view%stale_mesh_icon)
                    call set_status_hint(view%tree,item,scene_visibility_column,'Needs update')
                    view%visibility_dirty=.true.
                    if(present(changed))changed=.true.
                    cycle
                end if
                if(validity_changed)call item%setText(0,mesh%name)
                if(validity_changed.or.view%rows(slot)%mesh_cells/=mesh%elements)then
                    call set_status_hint(view%tree,item,0,mesh%name//' · '//integer_text(mesh%elements)//' cells')
                    view%rows(slot)%mesh_cells=mesh%elements
                end if
                state=merge(1,0,mesh%visible)
                if(validity_changed.or.item%data(scene_visibility_column)/=state)then
                    call set_row_visibility(view,item,state)
                    if(present(changed))changed=.true.
                end if
            end associate
        end do
    end subroutine fill_project_mesh_rows

    subroutine register_project_row(view,item,key,category,rank,id,branch,assembly)
        type(project_tree_view), intent(inout) :: view
        type(QTreeWidgetItem), intent(inout) :: item
        character(*), intent(in) :: key
        integer, optional, intent(in) :: category,rank,id,branch,assembly
        type(project_tree_row), allocatable :: grown(:)
        integer :: slot

        slot=view%row_count+1
        if(slot>size(view%rows))then
            allocate(grown(2*size(view%rows)))
            grown(:view%row_count)=view%rows
            call move_alloc(grown,view%rows)
        end if
        call item%setData(0,slot)
        view%row_count=slot
        associate(row=>view%rows(slot))
            row%item=item
            row%key=key
            if(present(category))then
                row%category=category
                call item%setFlags(ior(Qt_ItemIsEnabled,Qt_ItemIsSelectable))
            end if
            if(present(rank))row%rank=rank
            if(present(id))row%id=id
            if(present(branch))row%branch=branch
            if(present(assembly))row%assembly=assembly
        end associate
    end subroutine register_project_row

    integer function project_row_index(view,item) result(slot)
        type(project_tree_view), intent(in) :: view
        type(QTreeWidgetItem), intent(in) :: item

        slot=0
        if(item%isNull())return
        slot=item%data(0)
        if(slot<1.or.slot>view%row_count)then
            slot=0
        else if(.not.c_associated(item%ptr,view%rows(slot)%item%ptr))then
            slot=0
        end if
    end function project_row_index

    integer function project_row_visibility(view,item) result(state)
        type(project_tree_view), intent(in) :: view
        type(QTreeWidgetItem), intent(in) :: item

        state=-1
        if(project_row_index(view,item)>0)state=item%data(scene_visibility_column)
    end function project_row_visibility

    integer function mesh_project_id(view,item) result(id)
        type(project_tree_view), intent(in) :: view
        type(QTreeWidgetItem), intent(in) :: item
        integer :: slot

        id=0
        slot=project_row_index(view,item)
        if(slot==0)return
        if(view%rows(slot)%category==selection_mesh)id=view%rows(slot)%id
    end function mesh_project_id

    subroutine project_item_identity(view,s,item,kind,variant,group,branch_kind)
        type(project_tree_view), intent(in) :: view
        type(studio_t), intent(in) :: s
        type(QTreeWidgetItem), intent(in) :: item
        integer, intent(out) :: kind,variant,group,branch_kind
        integer :: slot

        kind=0
        variant=0
        group=0
        branch_kind=0
        slot=project_row_index(view,item)
        if(slot==0)return
        associate(row=>view%rows(slot))
            branch_kind=row%branch
            select case(row%category)
            case(selection_patch)
                if(.not.single_patch_available(row%rank,row%id,catalog=s%document%patches))return
                kind=row%rank
                variant=row%id
            case(selection_assembly)
                if(row%id<=count_multipatches(s%document%multipatches))group=row%id
            case default
                return
            end select
        end associate
    end subroutine project_item_identity

    function project_target_row(view,target) result(item)
        type(project_tree_view), intent(in) :: view
        type(selection_target_t), intent(in) :: target
        type(QTreeWidgetItem) :: item
        integer :: i,category

        category=target%category
        if(category<=0.or.target%id<=0)return
        if(category==selection_mesh)then
            if(.not.view%mesh_branch%isNull())item=view%mesh_branch%child(target%id-1)
            return
        end if
        if(category==selection_control)category=selection_patch
        do i=1,view%row_count
            associate(row=>view%rows(i))
                if(row%category/=category.or.row%id/=target%id.or.row%assembly/=0)cycle
                if(category==selection_patch.and.row%rank/=target%rank)cycle
                item=row%item
                return
            end associate
        end do
    end function project_target_row

    subroutine select_project_patch_row(view,s,kind)
        type(project_tree_view), intent(inout) :: view
        type(studio_t), intent(in) :: s
        integer, intent(in) :: kind
        type(QTreeWidgetItem) :: item
        integer :: id

        id=active_patch_id(kind,catalog=s%document%patches)
        item=project_target_row(view,selection_target_t(selection_patch,id,kind))
        if(item%isNull())return
        if(.not.item%isHidden())call expand_project_path(item)
        call view%tree%setCurrentItem(item)
        call sync_patch_row_visibility(view,item,s,kind,id)
    end subroutine select_project_patch_row

    subroutine sync_scene_visibility_icons(view,s,kind,variant)
        type(project_tree_view), intent(inout) :: view
        type(studio_t), intent(in) :: s
        integer, intent(in) :: kind,variant
        type(QTreeWidgetItem) :: item
        integer :: i

        do i=1,view%row_count
            associate(row=>view%rows(i))
                if(row%category/=selection_patch.or.row%rank/=kind.or.row%id/=variant)cycle
                item=row%item
                call sync_patch_row_visibility(view,item,s,kind,variant)
            end associate
        end do
    end subroutine sync_scene_visibility_icons

    subroutine sync_scene_visibility_batch(view,s,targets)
        type(project_tree_view), intent(inout) :: view
        type(studio_t), intent(in) :: s
        integer, intent(in) :: targets(:,:)
        logical, allocatable :: selected(:)
        type(QTreeWidgetItem) :: item
        integer :: offset(0:3),i,kind,id,variant

        if(size(targets,2)==0)return
        offset(0)=0
        do kind=1,3
            offset(kind)=offset(kind-1)+single_patch_count(kind,catalog=s%document%patches)
        end do
        allocate(selected(offset(3)),source=.false.)
        do i=1,size(targets,2)
            kind=targets(1,i)
            if(kind<1.or.kind>3)cycle
            id=offset(kind-1)+targets(2,i)
            if(id<=offset(kind-1).or.id>offset(kind))cycle
            selected(id)=.true.
        end do
        do i=1,view%row_count
            associate(row=>view%rows(i))
                if(row%category/=selection_patch)cycle
                kind=row%rank
                variant=row%id
                id=offset(kind-1)+variant
                if(id>offset(kind))cycle
                if(.not.selected(id))cycle
                item=row%item
                call sync_patch_row_visibility(view,item,s,kind,variant)
            end associate
        end do
    end subroutine sync_scene_visibility_batch

    subroutine sync_patch_row_visibility(view,item,s,kind,variant)
        type(project_tree_view), intent(inout) :: view
        type(QTreeWidgetItem), intent(inout) :: item
        type(studio_t), intent(in) :: s
        integer, intent(in) :: kind,variant
        integer :: state

        state=merge(1,0,studio_scene_is_visible(s%document%scene,kind,variant,catalog=s%document%patches))
        if(item%data(scene_visibility_column)/=state)call set_row_visibility(view,item,state)
    end subroutine sync_patch_row_visibility

    subroutine sync_patch_tree_item(view,item,s,kind,variant)
        type(project_tree_view), intent(inout) :: view
        type(QTreeWidgetItem), intent(inout) :: item
        type(studio_t), intent(in) :: s
        integer, intent(in) :: kind,variant
        character(:), allocatable :: tip
        integer :: slot
        logical :: visible

        slot=project_row_index(view,item)
        if(slot==0)return
        visible=studio_scene_is_visible(s%document%scene,kind,variant,catalog=s%document%patches)
        tip=single_patch_name(kind,variant,catalog=s%document%patches)//' · '//entity_name(kind)
        if(view%rows(slot)%assembly>0)then
            call item%setIcon(0,view%reference_icon)
            tip=tip//' · Multipatch member'
        else
            call item%setIcon(0,view%geometry_icons(kind))
        end if
        call set_status_hint(view%tree,item,0,tip)
        call set_row_visibility(view,item,merge(1,0,visible))
    end subroutine sync_patch_tree_item

    subroutine set_row_visibility(view,item,state)
        type(project_tree_view), intent(inout) :: view
        type(QTreeWidgetItem), intent(inout) :: item
        integer, intent(in) :: state
        type(QIcon) :: empty_icon

        view%visibility_dirty=.true.
        call item%setData(scene_visibility_column,state)
        select case(state)
        case(0)
            call item%setIcon(scene_visibility_column,view%hidden_icon)
            call set_status_hint(view%tree,item,scene_visibility_column,'Hidden · Click to show')
        case(1)
            call item%setIcon(scene_visibility_column,view%visible_icon)
            call set_status_hint(view%tree,item,scene_visibility_column,'Visible · Click to hide')
        case(2)
            call item%setIcon(scene_visibility_column,view%mixed_icon)
            call set_status_hint(view%tree,item,scene_visibility_column,'Partly visible · Click to show all')
        case default
            empty_icon=QIcon('')
            call item%setIcon(scene_visibility_column,empty_icon)
            call empty_icon%delete()
            call set_status_hint(view%tree,item,scene_visibility_column,'')
        end select
    end subroutine set_row_visibility

    subroutine update_project_visibility_columns(view)
        type(project_tree_view), intent(inout) :: view
        type(QTreeWidgetItem) :: item,child
        integer :: i,j,state,child_state

        if(.not.view%visibility_dirty)return
        do i=view%row_count,1,-1
            if(view%rows(i)%category==selection_patch.or.view%rows(i)%category==selection_mesh)cycle
            item=view%rows(i)%item
            state=-1
            do j=0,item%childCount()-1
                child=item%child(j)
                child_state=child%data(scene_visibility_column)
                if(child_state<0)cycle
                if(state<0)state=child_state
                if(state/=child_state)state=2
            end do
            if(item%data(scene_visibility_column)/=state)call set_row_visibility(view,item,state)
        end do
        view%visibility_dirty=.false.
    end subroutine update_project_visibility_columns

    subroutine apply_project_filter(view,s)
        type(project_tree_view), intent(inout) :: view
        type(studio_t), intent(in) :: s
        type(QScrollBar) :: scrollbar
        character(:), allocatable :: query
        integer :: results
        logical :: active,was_active,empty,restored

        if(view%tree%isNull())return
        call update_project_visibility_columns(view)
        query=''
        if(allocated(view%search_text))query=lowercase(trim(adjustl(view%search_text)))
        if(.not.allocated(view%applied_search))view%applied_search=''
        was_active=len(view%applied_search)>0.or.view%applied_filter/=scene_filter_all
        active=len(query)>0.or.view%filter_mode/=scene_filter_all
        if(query/=view%applied_search.or.view%filter_mode/=view%applied_filter)then
            if(allocated(view%reveal_key))deallocate(view%reveal_key)
        end if
        if(active.and..not.was_active)then
            view%expanded_before_filter=capture_project_tree(view)
            scrollbar=view%tree%verticalScrollBar()
            view%scroll_before_filter=scrollbar%value()
        end if
        view%applied_search=query
        view%applied_filter=view%filter_mode
        if(.not.view%filter_button%isNull())then
            call view%filter_button%setChecked(view%filter_mode/=scene_filter_all)
            call set_status_hint(view%filter_button,'Filter scene objects · '//view%filters(view%filter_mode+1)%text())
        end if
        if(active.or.was_active)then
            call filter_project_tree_rows(view,query,active,results)
        else
            results=count(view%rows(:view%row_count)%category>0)
        end if
        if(.not.active.and.was_active.and.allocated(view%expanded_before_filter))then
            call restore_project_tree(view,view%expanded_before_filter,restored,expansion_only=.true.)
            deallocate(view%expanded_before_filter)
            scrollbar=view%tree%verticalScrollBar()
            call scrollbar%setValue(view%scroll_before_filter)
        end if
        if(view%empty_state%isNull())return
        empty=results==0
        call view%tree%setVisible(.not.empty)
        call view%empty_state%setVisible(empty)
        if(.not.empty)return
        call view%empty_title%setText(merge('No matching objects','No geometry yet    ',active))
        call view%empty_clear%setVisible(active)
    end subroutine apply_project_filter

    subroutine filter_project_tree_rows(view,query,active,results)
        type(project_tree_view), intent(in) :: view
        character(*), intent(in) :: query
        logical, intent(in) :: active
        integer, intent(out) :: results
        type(QTreeWidgetItem) :: root
        integer :: i
        logical :: matched

        results=0
        do i=0,view%tree%topLevelItemCount()-1
            root=view%tree%topLevelItem(i)
            matched=filter_project_row(view,root,query,.false.,active,results)
        end do
    end subroutine filter_project_tree_rows

    recursive logical function filter_project_row(view,item,query,ancestor_matches,active,results) result(visible)
        type(project_tree_view), intent(in) :: view
        type(QTreeWidgetItem), intent(inout) :: item
        character(*), intent(in) :: query
        logical, intent(in) :: ancestor_matches,active
        integer, intent(inout) :: results
        type(QTreeWidgetItem) :: child
        logical :: name_matches,self_matches,has_child,child_matches,revealed
        integer :: i,slot

        slot=project_row_index(view,item)
        visible=.false.
        if(slot==0)return
        name_matches=ancestor_matches.or.len(query)==0.or.index(lowercase(item%text(0)),query)>0
        has_child=.false.
        do i=0,item%childCount()-1
            child=item%child(i)
            child_matches=filter_project_row(view,child,query,name_matches,active,results)
            has_child=has_child.or.child_matches
        end do
        associate(row=>view%rows(slot))
            revealed=.false.
            if(allocated(view%reveal_key))revealed=row%key==view%reveal_key
            self_matches=row%category>0.and.name_matches.and.row_matches_project_filter(row,item,view%filter_mode)
            if(self_matches.or.revealed)results=results+1
        end associate
        visible=.not.active.or.self_matches.or.has_child.or.revealed
        call item%setHidden(.not.visible)
        if(active.and.has_child)call item%setExpanded(.true.)
    end function filter_project_row

    logical function row_matches_project_filter(row,item,mode) result(matches)
        type(project_tree_row), intent(in) :: row
        type(QTreeWidgetItem), intent(in) :: item
        integer, intent(in) :: mode

        matches=.false.
        select case(mode)
        case(scene_filter_all)
            matches=.true.
        case(scene_filter_geometry)
            matches=row%category==selection_patch.and.row%assembly==0
        case(scene_filter_multipatches)
            matches=row%category==selection_assembly.or.row%assembly>0
        case(scene_filter_meshes)
            matches=row%category==selection_mesh
        case(scene_filter_hidden)
            matches=item%data(scene_visibility_column)==0
        case(scene_filter_stale)
            matches=row%stale
        case(scene_filter_curves:scene_filter_volumes)
            matches=row%category==selection_patch.and.row%assembly==0.and.row%rank==mode-scene_filter_curves+1
        case default
            return
        end select
    end function row_matches_project_filter

    logical function reveal_project_row(view,target) result(revealed)
        type(project_tree_view), intent(inout) :: view
        type(selection_target_t), intent(in) :: target
        type(QTreeWidgetItem) :: item
        integer :: slot,results

        item=project_target_row(view,target)
        revealed=.not.item%isNull()
        if(.not.revealed)return
        slot=project_row_index(view,item)
        view%reveal_key=view%rows(slot)%key
        if(allocated(view%applied_search))then
            if(len(view%applied_search)>0.or.view%applied_filter/=scene_filter_all)then
                call filter_project_tree_rows(view,view%applied_search,.true.,results)
            end if
        end if
        call expand_project_path(item,unhide=.true.)
        call view%tree%show()
        if(.not.view%empty_state%isNull())call view%empty_state%hide()
        call view%tree%setCurrentItem(item)
        call view%tree%scrollToItem(item)
    end function reveal_project_row

    subroutine expand_project_path(item,unhide)
        type(QTreeWidgetItem), intent(inout) :: item
        logical, optional, intent(in) :: unhide
        type(QTreeWidgetItem) :: ancestor,next
        logical :: show_path

        show_path=.false.
        if(present(unhide))show_path=unhide
        if(show_path)call item%setHidden(.false.)
        ancestor=item%parent()
        do while(.not.ancestor%isNull())
            if(show_path)call ancestor%setHidden(.false.)
            call ancestor%setExpanded(.true.)
            next=ancestor%parent()
            ancestor=next
        end do
    end subroutine expand_project_path

    function capture_project_tree(view) result(states)
        type(project_tree_view), intent(in) :: view
        type(project_tree_state), allocatable :: states(:)
        type(QTreeWidgetItem) :: item,current
        logical :: selected,focused
        integer :: i,n

        allocate(states(view%row_count))
        current=view%tree%currentItem()
        n=0
        do i=1,view%row_count
            item=view%rows(i)%item
            selected=item%isSelected()
            focused=c_associated(item%ptr,current%ptr)
            if(.not.(selected.or.focused.or.item%childCount()>0))cycle
            n=n+1
            states(n)=project_tree_state(view%rows(i)%key,item%isExpanded(),selected,focused)
        end do
        states=states(:n)
    end function capture_project_tree

    subroutine restore_project_tree(view,states,current_restored,expansion_only)
        type(project_tree_view), intent(inout) :: view
        type(project_tree_state), intent(in) :: states(:)
        logical, intent(out) :: current_restored
        logical, optional, intent(in) :: expansion_only
        type(QTreeWidgetItem) :: item,current
        logical :: keep_selection
        integer :: i,j

        current_restored=.false.
        keep_selection=.true.
        if(present(expansion_only))keep_selection=.not.expansion_only
        if(keep_selection)then
            do j=1,size(states)
                if(.not.states(j)%current)cycle
                do i=1,view%row_count
                    if(view%rows(i)%key/=states(j)%key)cycle
                    current=view%rows(i)%item
                    exit
                end do
                exit
            end do
            if(.not.current%isNull())then
                call view%tree%setCurrentItem(current)
                current_restored=.true.
            end if
        end if
        do i=1,view%row_count
            do j=1,size(states)
                if(view%rows(i)%key/=states(j)%key)cycle
                item=view%rows(i)%item
                if(keep_selection)call item%setSelected(states(j)%selected)
                if(item%childCount()>0)call item%setExpanded(states(j)%expanded)
                exit
            end do
        end do
    end subroutine restore_project_tree

    subroutine reset_project_filter(view)
        type(project_tree_view), intent(inout) :: view
        integer :: i

        view%changing_filter=.true.
        view%search_text=''
        view%applied_search=''
        view%filter_mode=scene_filter_all
        view%applied_filter=scene_filter_all
        if(allocated(view%reveal_key))deallocate(view%reveal_key)
        if(allocated(view%expanded_before_filter))deallocate(view%expanded_before_filter)
        if(.not.view%search%isNull())call view%search%clear()
        do i=1,size(view%filters)
            if(.not.view%filters(i)%isNull())call view%filters(i)%setChecked(i==1)
        end do
        view%changing_filter=.false.
    end subroutine reset_project_filter


    subroutine update_project_summary(view,s,summary)
        type(project_tree_view), intent(inout) :: view
        type(studio_t), intent(in) :: s
        character(:), allocatable, intent(out) :: summary
        integer :: i,visible_patches,visible_meshes
        character(:), allocatable :: patches,meshes

        if(view%filter_mode==scene_filter_hidden.and.view%visibility_dirty)then
            call apply_project_filter(view,s)
        else
            call update_project_visibility_columns(view)
        end if
        visible_patches=scene_drawn_patch_count(s)
        visible_meshes=count([(scene_mesh_drawn(s,i),i=1,studio_mesh_count(s))])
        patches=integer_text(visible_patches)//trim(merge(' patch  ',' patches',visible_patches==1))
        meshes=integer_text(visible_meshes)//trim(merge(' mesh  ',' meshes',visible_meshes==1))
        if(s%mesh_only)then
            summary='Meshes only · '//meshes
        else if(visible_meshes>0)then
            summary=patches//' · '//meshes//' visible'
        else
            summary=patches//' visible'
        end if
    end subroutine update_project_summary

    function project_patch_name(s,kind,id) result(name)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: kind,id
        character(:), allocatable :: name

        name=single_patch_name(kind,id,catalog=s%document%patches)
    end function project_patch_name

    function patch_row_key(rank,id) result(key)
        integer, intent(in) :: rank,id
        character(:), allocatable :: key

        key='patch:'//integer_text(rank)//'/'//integer_text(id)
    end function patch_row_key

    pure function lowercase(input) result(output)
        character(*), intent(in) :: input
        character(len(input)) :: output
        integer :: i,code

        output=input
        do i=1,len(input)
            code=iachar(input(i:i))
            if(code>=iachar('A').and.code<=iachar('Z'))output(i:i)=achar(code+32)
        end do
    end function lowercase

end module studio_project_tree
