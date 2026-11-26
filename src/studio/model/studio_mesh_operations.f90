module studio_mesh_operations
    use, intrinsic :: iso_fortran_env, only: int64
    use studio_system_memory, only: mesh_memory_budget
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use studio_document, only: studio_document_t, document_has_patch, document_mesh_count, document_mesh_source_group, &
        append_document_mesh, record_document_mesh_edit, remove_document_mesh
    use studio_forcad, only: rk, check_forcad_mesh_joins, estimate_forcad_mesh_samples
    use studio_geometry_operations, only: read_geometry_target
    use studio_mesh, only: studio_mesh_t, build_studio_mesh, stitch_studio_volume_mesh, move_studio_mesh, mesh_samples_fit
    use studio_single_patch, only: copy_single_patch_geometry,next_patch_revision
    use studio_multipatch_collection, only: count_multipatches, refresh_multipatch
    use studio_presentation_text, only: integer_text
    implicit none(type, external)
    private

    public :: mesh_source_t, mesh_publication_t, read_document_mesh_source, find_document_mesh, generate_document_mesh
    public :: estimate_document_mesh_samples
    integer, parameter, public :: mesh_source_missing=1, mesh_build_failed=2

    type :: mesh_source_t
        integer :: kind=0, assembly_id=0, mesh_id=0
        integer, allocatable :: members(:)
        character(:), allocatable :: name
    end type mesh_source_t

    type :: mesh_publication_t
        integer :: id=0, cells=0, failure=0
        logical :: remeshed=.false.
        integer, allocatable :: removed_ids(:)
    end type mesh_publication_t

contains

    function read_document_mesh_source(document,kind,patch_id,assembly_id,mesh_id) result(source)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: kind,patch_id,assembly_id,mesh_id
        type(mesh_source_t) :: source
        integer :: group

        if(mesh_id==0)then
            source%assembly_id=assembly_id
            source%mesh_id=find_document_mesh(document,kind,patch_id,assembly_id)
            call read_geometry_target(document,kind,patch_id,assembly_id,source%kind,source%members,source%name)
            return
        end if
        allocate(source%members(0))
        if(mesh_id<1.or.mesh_id>document_mesh_count(document))return
        source%mesh_id=mesh_id
        associate(mesh=>document%meshes(mesh_id))
            source%kind=mesh%kind
            source%name=mesh%name
            group=document_mesh_source_group(document,mesh_id)
            source%assembly_id=group
            if(group>0)then
                source%kind=document%multipatches%groups(group)%definition%kind
                source%members=document%multipatches%groups(group)%definition%members
            else if(mesh%group==0)then
                if(document_has_patch(document,source%kind,mesh%variant))source%members=[mesh%variant]
            end if
        end associate
    end function read_document_mesh_source

    integer function find_document_mesh(document,kind,patch_id,assembly_id) result(id)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: kind,patch_id,assembly_id
        integer :: candidate

        id=0
        if(assembly_id==0.and.patch_id==0)return
        do candidate=document_mesh_count(document),1,-1
            if(.not.mesh_matches_source(document,candidate,kind,assembly_id,patch_id))cycle
            id=candidate
            return
        end do
    end function find_document_mesh

    logical function mesh_matches_source(document,id,kind,group,variant) result(matches)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: id,kind,group,variant

        if(group>0)then
            matches=document_mesh_source_group(document,id)==group
        else
            associate(mesh=>document%meshes(id))
                matches=mesh%group==0.and.mesh%kind==kind.and.mesh%variant==variant
            end associate
        end if
    end function mesh_matches_source

    subroutine estimate_document_mesh_samples(document,source,target_size,synchronize,cells,error)
        type(studio_document_t), intent(in) :: document
        type(mesh_source_t), intent(in) :: source
        real(c_double), intent(in) :: target_size
        logical, intent(in) :: synchronize
        integer, allocatable, intent(out) :: cells(:,:)
        character(:), allocatable, intent(out) :: error
        class(*), allocatable :: patch
        integer, allocatable :: samples(:,:)
        integer :: kind,group,maximum
        integer(int64) :: budget,probes,cross_section
        integer(int64), parameter :: preview_work_bytes=64_int64*1024*1024

        error='Select a patch or assembly with an available mesh source.'
        if(.not.allocated(source%members))return
        if(size(source%members)==0)return
        kind=source%kind
        group=source%assembly_id
        if(kind<1.or.kind>3.or.group<0)return
        error='Enter a positive finite target size.'
        if(.not.ieee_is_finite(target_size).or.target_size<=0._c_double)return

        budget=mesh_memory_budget()/mesh_sample_bytes(kind)
        probes=min(budget,preview_work_bytes/mesh_sample_bytes(kind))
        cross_section=int(size(source%members),int64)*2_int64**(kind-1)
        maximum=int(min(probes,int(huge(0),int64)/3)/cross_section)
        error='The target-size preview exceeds the available mesh work budget.'
        if(maximum<3)return

        error='Select a patch or assembly with an available mesh source.'
        if(group==0)then
            if(size(source%members)/=1)return
            if(.not.document_has_patch(document,kind,source%members(1)))return
            call copy_single_patch_geometry(kind,source%members(1),patch,catalog=document%patches)
            if(.not.allocated(patch))return
            call estimate_forcad_mesh_samples(patch,kind,real(target_size,rk),synchronize,maximum,samples,error)
        else
            if(group>count_multipatches(document%multipatches))return
            associate(assembly=>document%multipatches%groups(group))
                if(.not.allocated(assembly%native))return
                if(size(assembly%definition%members)/=size(source%members))return
                if(any(assembly%definition%members/=source%members))return
                call estimate_forcad_mesh_samples(assembly%native,kind,real(target_size,rk),synchronize,maximum,samples,error)
            end associate
        end if
        if(len(error)>0)return
        error='The estimated mesh exceeds the supported node or cell index range. Increase the target size.'
        if(.not.mesh_samples_fit(kind,samples))return
        error='The estimated mesh exceeds the available memory budget. Increase the target size.'
        if(.not.mesh_target_fits_memory(kind,samples))return
        cells=samples(:kind,:)-1
        error=''
    end subroutine estimate_document_mesh_samples

    integer function mesh_sample_bytes(kind) result(bytes)
        integer, intent(in) :: kind

        bytes=(4*(3+kind)*storage_size(0._rk)+4*2**kind*storage_size(0))/8
    end function mesh_sample_bytes

    logical function mesh_target_fits_memory(kind,samples) result(fits)
        integer, intent(in) :: kind,samples(:,:)
        integer :: i
        integer(int64) :: total_samples

        fits=.false.
        if(.not.mesh_samples_fit(kind,samples))return
        total_samples=sum([(product(int(samples(:kind,i),int64)),i=1,size(samples,2))])
        fits=total_samples<=mesh_memory_budget()/mesh_sample_bytes(kind)
    end function mesh_target_fits_memory

    subroutine generate_document_mesh(document,source,cell_counts,conforming,selected_mesh,mesh_only,published,error,patch_cells, &
        target_size,synchronize)
        type(studio_document_t), intent(inout) :: document
        type(mesh_source_t), intent(in) :: source
        integer, intent(in) :: cell_counts(:),selected_mesh
        logical, intent(in) :: mesh_only
        logical, intent(in), optional :: conforming,synchronize
        real(c_double), intent(in), optional :: target_size
        type(mesh_publication_t), intent(out) :: published
        character(:), allocatable, intent(out) :: error
        integer, intent(in), optional :: patch_cells(:,:)
        type(studio_mesh_t) :: candidate
        integer, allocatable :: patch_samples(:,:)
        integer :: kind,group,variant,mesh_count,i
        integer(int64) :: construction_stamp
        logical :: connected,saved

        error='Select a patch, assembly, or mesh with an available source.'
        published%failure=mesh_source_missing
        if(.not.allocated(source%members).or..not.allocated(source%name))return
        if(size(source%members)==0)return
        kind=source%kind
        group=source%assembly_id
        if(kind<1.or.kind>3.or.group<0)return
        mesh_count=document_mesh_count(document)
        if(source%mesh_id<0.or.source%mesh_id>mesh_count)return

        published%failure=0
        connected=.false.
        if(source%mesh_id>0)connected=allocated(document%meshes(source%mesh_id)%connected_cells)
        if(present(conforming))connected=conforming
        error='Select a volume assembly to generate a conforming hex mesh.'
        if(connected.and.(kind/=3.or.group<=0))return
        published%failure=mesh_build_failed
        allocate(patch_samples(3,size(source%members)),source=1)
        if(present(patch_cells))then
            error='Supply one row of positive cell counts for every assembly patch.'
            if(group<=0.or.any(shape(patch_cells)/=[kind,size(source%members)]))return
            if(any(patch_cells<1).or.any(patch_cells>=huge(1)))return
            patch_samples(:kind,:)=patch_cells+1
        else
            error='Supply one positive cell count for each active direction.'
            if(size(cell_counts)/=kind.or.any(cell_counts<1).or.any(cell_counts>=huge(1)))return
            patch_samples(:kind,:)=spread(cell_counts+1,2,size(source%members))
        end if

        if(present(target_size))then
            if(target_size>0._c_double)then
                error='The requested mesh exceeds available memory or indices. Preview a larger target size.'
                if(.not.mesh_target_fits_memory(kind,patch_samples))return
            end if
        end if
        variant=source%members(1)
        construction_stamp=0
        if(source%mesh_id>0.and.group>0)then
            do i=1,mesh_count
                if(.not.mesh_matches_source(document,i,kind,group,variant))cycle
                if(document%meshes(i)%construction_stamp>0)construction_stamp=document%meshes(i)%construction_stamp
            end do
        end if
        if(group>0)then
            error='Resolve the multipatch interfaces before generating a mesh.'
            if(group>count_multipatches(document%multipatches))return
            call refresh_multipatch(document%multipatches,group,.false.,catalog=document%patches)
            associate(assembly=>document%multipatches%groups(group))
                if(.not.allocated(assembly%native))return
                error='The assembly members changed. Refresh the mesh source and try again.'
                if(size(assembly%definition%members)/=size(source%members))return
                if(any(assembly%definition%members/=source%members))return
                if(connected)then
                    error='The mesh needs a sampled volume assembly with validated joins.'
                    if(.not.allocated(assembly%definition%joins))return
                    call check_forcad_mesh_joins(patch_samples,source%members,assembly%definition%joins,error)
                    if(len(error)>0)return
                end if
                call build_studio_mesh(kind,0,group,source%members,patch_samples(:,1),candidate,error, &
                    assembly%native,catalog=document%patches,patch_samples=patch_samples)
                candidate%source_group=assembly%definition%name
            end associate
        else
            call build_studio_mesh(kind,variant,0,source%members,patch_samples(:,1),candidate,error,catalog=document%patches)
        end if
        if(len(error)>0)return
        candidate%per_patch_sizing=present(patch_cells)
        candidate%construction_stamp=construction_stamp
        if(present(target_size))candidate%target_size=real(target_size,rk)
        if(present(synchronize))candidate%synchronize=synchronize
        if(connected)then
            call stitch_studio_volume_mesh(candidate,document%multipatches%groups(group)%definition%joins,error)
            if(len(error)>0)return
        end if

        published%cells=candidate%elements
        if(source%mesh_id>0)then
            call replace_document_mesh(document,source,candidate,selected_mesh,mesh_only,published,error)
            if(len(error)>0)return
        else
            candidate%name=next_mesh_name(document,source%name)
            published%id=mesh_count+1
            call record_document_mesh_edit(document,[integer::],published%id,selected_mesh,mesh_only,saved)
            if(.not.saved)then
                error='Could not write the mesh history.'
                return
            end if
            call append_document_mesh(document,candidate)
        end if
        document%mesh_generation=document%mesh_generation+1
        published%failure=0
    end subroutine generate_document_mesh

    subroutine replace_document_mesh(document,source,candidate,selected_mesh,mesh_only,published,error)
        type(studio_document_t), intent(inout) :: document
        type(mesh_source_t), intent(in) :: source
        type(studio_mesh_t), intent(inout) :: candidate
        integer, intent(in) :: selected_mesh
        logical, intent(in) :: mesh_only
        type(mesh_publication_t), intent(inout) :: published
        character(:), allocatable, intent(out) :: error
        integer, allocatable :: prior_ids(:)
        integer :: id,i,mesh_count
        logical :: saved

        id=source%mesh_id
        mesh_count=document_mesh_count(document)
        candidate%name=document%meshes(id)%name
        if(candidate%construction_stamp==0)candidate%construction_stamp=document%meshes(id)%construction_stamp
        candidate%visible=document%meshes(id)%visible
        prior_ids=pack([(i,i=1,mesh_count)],[(i==id.or. &
            mesh_matches_source(document,i,source%kind,source%assembly_id,source%members(1)),i=1,mesh_count)])
        published%remeshed=.true.
        published%id=id-count(prior_ids<id)
        published%removed_ids=pack(prior_ids,prior_ids/=id)
        call record_document_mesh_edit(document,prior_ids,published%id,selected_mesh,mesh_only,saved)
        error='Could not write the mesh history.'
        if(.not.saved)return
        candidate%content_revision=next_patch_revision()
        call move_studio_mesh(candidate,document%meshes(id))
        do i=size(published%removed_ids),1,-1
            call remove_document_mesh(document,published%removed_ids(i))
        end do
        error=''
    end subroutine replace_document_mesh

    function next_mesh_name(document,source_name) result(name)
        type(studio_document_t), intent(in) :: document
        character(*), intent(in) :: source_name
        character(:), allocatable :: name,prefix
        integer :: mesh_count,suffix,i

        mesh_count=document_mesh_count(document)
        suffix=mesh_count+1
        do
            prefix='Mesh '//integer_text(suffix)//' · '
            if(.not.any([(index(document%meshes(i)%name,prefix)==1,i=1,mesh_count)]))exit
            suffix=suffix+1
        end do
        name=prefix//source_name
    end function next_mesh_name

end module studio_mesh_operations
