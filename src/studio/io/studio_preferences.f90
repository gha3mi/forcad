module studio_preferences
    use, intrinsic :: iso_c_binding, only: c_int,c_null_char
    use studio_platform_files, only: file_metadata,inspect_path,is_regular_file,create_private_file, &
        current_directory,close_file,restrict_file
    use studio_commands, only: viewport_layer_count, parameter_lines_layer, default_viewport_layers, workspace_panels, &
        viewport_shaded,viewport_edges,viewport_style_count
    use, intrinsic :: iso_fortran_env, only: error_unit
    use qt_core_classes, only: QSettings, QSettings_NoError
    use studio_storage_paths, only: studio_application_directory
    implicit none(type, external)
    private
    public :: workspace_preferences_t, read_preferences, write_preferences, clear_preferences

    public :: application_settings

    character(*), parameter :: settings_group='workspace-v1'

    type :: workspace_preferences_t
        integer :: entity_kind=2, knot_axis=1, surface_tool=1, viewport_style=viewport_edges
        integer :: sample_ids(3)=0
        logical :: save_sample(3)=.false.
        logical :: basis_derivative=.false., basis_legend=.true., basis_probe_guide=.true., world_grid=.true., interactive_quality=.true.
        logical :: layers(viewport_layer_count)=default_viewport_layers
        logical :: domain_layers(5)=[.true.,.true.,.true.,.false.,.false.], visible_panels(size(workspace_panels))=.true.
        character(:), allocatable :: window_geometry
        integer :: columns(3)=[264,1041,342], center_rows(2)=[488,293]
    end type workspace_preferences_t
contains
    function application_settings() result(settings)
        type(QSettings) :: settings
        character(:),allocatable :: folder,path
        type(file_metadata) :: metadata
        integer(c_int) :: file,closed

        folder=studio_application_directory()
        settings=QSettings('')
        if(len(folder)==0)return
        path=folder//'/studio.ini'
        if(.not.inspect_path(path,metadata))then
            file=create_private_file(current_directory,path//c_null_char)
            if(file>=0)closed=close_file(file)
            if(.not.inspect_path(path,metadata))return
        end if
        if(.not.is_regular_file(metadata).or.metadata%owned/=1.or.metadata%links/=1)return
        if(metadata%is_private/=1)then
            if(restrict_file(path//c_null_char)/=0)return
        end if
        call settings%delete()
        settings=QSettings(path)
    end function application_settings

    function open_settings(file_name) result(settings)
        character(*), optional, intent(in) :: file_name
        type(QSettings) :: settings

        if(present(file_name))then
            settings=QSettings(file_name)
        else
            settings=application_settings()
        end if
        call settings%beginGroup(settings_group)
    end function open_settings

    function item_key(prefix,position) result(key)
        character(*), intent(in) :: prefix
        integer, intent(in) :: position
        character(:), allocatable :: key
        character(12) :: number

        write(number,'(i0)')position
        key=prefix//'/'//trim(number)
    end function item_key

    subroutine read_preferences(values,restored,file_name)
        type(workspace_preferences_t), intent(inout) :: values
        logical, intent(out) :: restored
        character(*), optional, intent(in) :: file_name
        type(QSettings) :: settings
        integer :: i,legacy_basis

        settings=open_settings(file_name)
        restored=settings%intValue('version',0)==1
        if(restored)then
            values%entity_kind=max(1,min(3,settings%intValue('entity',values%entity_kind)))
            do i=1,size(values%sample_ids)
                values%sample_ids(i)=settings%intValue(item_key('samples',i),values%sample_ids(i))
            end do
            legacy_basis=max(1,min(6,settings%intValue('basis',values%knot_axis+merge(3,0,values%basis_derivative))))
            values%basis_derivative=legacy_basis>3
            values%knot_axis=settings%intValue('knot-axis',mod(legacy_basis-1,3)+1)
            if(values%knot_axis<1.or.values%knot_axis>values%entity_kind)values%knot_axis=1
            values%surface_tool=max(1,min(2,settings%intValue('tool',values%surface_tool)))
            do i=1,size(values%layers)
                values%layers(i)=settings%boolValue(item_key('layers',i),values%layers(i))
            end do
            values%layers(parameter_lines_layer)=values%layers(parameter_lines_layer).or. &
                settings%boolValue('layers/6',.false.).or.settings%boolValue('layers/7',.false.)
            do i=1,size(values%domain_layers)
                values%domain_layers(i)=settings%boolValue(item_key('overlays',i),values%domain_layers(i))
            end do
            values%basis_legend=settings%boolValue('basis-legend',values%basis_legend)
            values%basis_probe_guide=settings%boolValue('basis-probe-guide',values%basis_probe_guide)
            values%world_grid=settings%boolValue('world-grid',values%world_grid)
            values%viewport_style=settings%intValue('viewport-style',viewport_edges)
            if(values%viewport_style<viewport_shaded.or.values%viewport_style>viewport_style_count)values%viewport_style=viewport_edges
            values%interactive_quality=settings%boolValue('interactive-quality',values%interactive_quality)
            do i=1,size(workspace_panels)
                values%visible_panels(i)=settings%boolValue(item_key('panels',workspace_panels(i)),.true.)
            end do
            if(settings%contains('geometry'))values%window_geometry=settings%stringValue('geometry')
            call read_sizes(settings,'columns',values%columns,[1,2,3])
            call read_sizes(settings,'center',values%center_rows,[1,2])
        end if
        call settings%delete()
    end subroutine read_preferences

    subroutine read_sizes(settings,prefix,sizes,keys)
        type(QSettings), intent(in) :: settings
        character(*), intent(in) :: prefix
        integer, intent(inout) :: sizes(:)
        integer, intent(in) :: keys(:)
        integer :: i,saved_size

        do i=1,size(keys)
            saved_size=settings%intValue(item_key(prefix,keys(i)),sizes(i))
            if(saved_size>=1.and.saved_size<=10000)sizes(i)=saved_size
        end do
    end subroutine read_sizes

    subroutine write_preferences(values,file_name)
        type(workspace_preferences_t), intent(in) :: values
        character(*), optional, intent(in) :: file_name
        type(QSettings) :: settings
        integer :: i

        settings=open_settings(file_name)
        call settings%setValue('version',1)
        call settings%setValue('entity',values%entity_kind)
        call settings%setValue('basis',values%knot_axis+merge(3,0,values%basis_derivative))
        do i=1,size(values%sample_ids)
            if(values%save_sample(i))call settings%setValue(item_key('samples',i),values%sample_ids(i))
        end do
        call settings%setValue('knot-axis',values%knot_axis)
        call settings%setValue('tool',values%surface_tool)
        if(allocated(values%window_geometry))call settings%setValue('geometry',values%window_geometry)
        do i=1,size(values%layers)
            call settings%setValue(item_key('layers',i),values%layers(i))
        end do
        call settings%remove('layers/6')
        call settings%remove('layers/7')
        do i=1,size(values%domain_layers)
            call settings%setValue(item_key('overlays',i),values%domain_layers(i))
        end do
        call settings%setValue('basis-legend',values%basis_legend)
        call settings%setValue('basis-probe-guide',values%basis_probe_guide)
        call settings%setValue('world-grid',values%world_grid)
        call settings%setValue('viewport-style',values%viewport_style)
        call settings%setValue('interactive-quality',values%interactive_quality)
        do i=1,size(workspace_panels)
            call settings%setValue(item_key('panels',workspace_panels(i)),values%visible_panels(i))
        end do
        call settings%remove('panels/3')
        call settings%remove('left')
        call write_sizes(settings,'columns',values%columns,[1,2,3])
        call write_sizes(settings,'center',values%center_rows,[1,2])
        call settings%sync()
        if(settings%status()/=QSettings_NoError)write(error_unit,'(a)')'Could not save ForCAD workspace preferences.'
        call settings%delete()
    end subroutine write_preferences

    subroutine write_sizes(settings,prefix,sizes,keys)
        type(QSettings), intent(in) :: settings
        character(*), intent(in) :: prefix
        integer, intent(in) :: sizes(:),keys(:)
        integer :: i

        do i=1,size(sizes)
            if(sizes(i)>0)call settings%setValue(item_key(prefix,keys(i)),sizes(i))
        end do
    end subroutine write_sizes

    subroutine clear_preferences()
        type(QSettings) :: settings

        settings=open_settings()
        call settings%remove('')
        call settings%sync()
        if(settings%status()/=QSettings_NoError)write(error_unit,'(a)')'Could not clear ForCAD workspace preferences.'
        call settings%delete()
    end subroutine clear_preferences
end module studio_preferences
