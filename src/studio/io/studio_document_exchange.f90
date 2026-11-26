module studio_document_exchange
    use studio_document, only: studio_document_t
    use studio_exchange, only: import_studio_iges, export_studio_patch, export_studio_multipatch
    use studio_multipatch_collection, only: count_multipatches, refresh_multipatch
    implicit none(type, external)
    private

    public :: import_document_iges, export_document_geometry

contains

    subroutine import_document_iges(document,path,kind,patch_id,count,skipped,error,source_unit,needs_source_unit)
        type(studio_document_t), intent(inout) :: document
        character(*), intent(in) :: path
        integer, intent(out) :: kind,patch_id,count,skipped
        character(:), allocatable, intent(out) :: error
        integer, optional, intent(in) :: source_unit
        logical, optional, intent(out) :: needs_source_unit

        call import_studio_iges(path,kind,patch_id,count,skipped,error,catalog=document%patches, &
            model_unit=document%model_unit,source_unit=source_unit,needs_source_unit=needs_source_unit, &
            project_directory=document%project_directory)
    end subroutine import_document_iges

    subroutine export_document_geometry(document,kind,patch_id,assembly_id,path,error)
        type(studio_document_t), intent(inout) :: document
        integer, intent(in) :: kind,patch_id,assembly_id
        character(*), intent(in) :: path
        character(:), allocatable, intent(out) :: error

        if(assembly_id>0)then
            error='Resolve the multipatch interfaces before exporting.'
            if(assembly_id>count_multipatches(document%multipatches))return
            call refresh_multipatch(document%multipatches,assembly_id,.false.,catalog=document%patches)
            associate(assembly=>document%multipatches%groups(assembly_id))
                if(.not.allocated(assembly%native))return
                call export_studio_multipatch(assembly%native,path,error)
            end associate
        else
            call export_studio_patch(kind,patch_id,path,error,model_unit=document%model_unit,catalog=document%patches)
        end if
    end subroutine export_document_geometry
end module studio_document_exchange
