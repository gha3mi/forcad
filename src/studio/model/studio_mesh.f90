! Keep the sampled mesh returned by ForCAD separate from editable patch geometry.
module studio_mesh
    use, intrinsic :: iso_fortran_env, only: int64
    use studio_forcad, only: rk, studio_mesh_part=>forcad_mesh_samples,sample_forcad_patch,sample_forcad_assembly, &
        volume_assembly_size,join_forcad_volume_mesh,studio_native_problem
    use studio_single_patch, only: studio_patch_catalog,copy_single_patch_geometry,single_patch_available,single_patch_mesh_revision
    use studio_multipatch, only: multipatch_join
    implicit none(type, external)
    private
    public :: studio_mesh_t,studio_mesh_part,build_studio_mesh,move_studio_mesh
    public :: build_planned_volume_mesh,stitch_studio_volume_mesh
    public :: mesh_axis_cell_limit, mesh_samples_fit

    type :: studio_mesh_t
        integer :: kind=0,variant=0,group=0,samples(3)=0,nodes=0,elements=0
        integer(int64) :: construction_stamp=0_int64,content_revision=0_int64
        real(rk) :: target_size=0._rk
        logical :: visible=.true., per_patch_sizing=.false., synchronize=.true.
        character(:), allocatable :: name, source_group
        type(studio_mesh_part), allocatable :: parts(:)
        real(rk), allocatable :: connected_points(:,:)
        integer, allocatable :: connected_cells(:,:)
    end type studio_mesh_t
contains
    pure integer function mesh_axis_cell_limit(kind,patch_count) result(limit)
        integer, intent(in) :: kind,patch_count
        integer(int64) :: min_cross_section

        limit=0
        if(kind<1.or.kind>3.or.patch_count<1)return
        min_cross_section=int(patch_count,int64)*2**(kind-1)
        limit=max(0,int(int(huge(1),int64)/(3*min_cross_section))-1)
    end function mesh_axis_cell_limit

    pure logical function mesh_samples_fit(kind,samples) result(fits)
        integer, intent(in) :: kind,samples(:,:)
        integer(int64) :: node_budget,cell_budget,nodes,cells
        integer :: part,axis

        fits=.false.
        if(kind<1.or.kind>3.or.size(samples,1)/=3.or.size(samples,2)==0)return
        if(any(samples(:kind,:)<2))return
        node_budget=int(huge(1),int64)/3
        cell_budget=int(huge(1),int64)/(2**kind)
        do part=1,size(samples,2)
            nodes=1
            cells=1
            do axis=1,kind
                if(samples(axis,part)>node_budget/nodes.or.samples(axis,part)-1>cell_budget/cells)return
                nodes=nodes*samples(axis,part)
                cells=cells*(samples(axis,part)-1)
            end do
            node_budget=node_budget-nodes
            cell_budget=cell_budget-cells
        end do
        fits=.true.
    end function mesh_samples_fit

    subroutine move_studio_mesh(source,destination)

        type(studio_mesh_t),intent(inout) :: source
        type(studio_mesh_t),intent(out) :: destination

        destination%kind=source%kind
        destination%variant=source%variant
        destination%group=source%group
        destination%samples=source%samples
        destination%nodes=source%nodes
        destination%elements=source%elements
        destination%visible=source%visible
        destination%per_patch_sizing=source%per_patch_sizing
        destination%construction_stamp=source%construction_stamp
        destination%content_revision=source%content_revision
        destination%target_size=source%target_size
        destination%synchronize=source%synchronize
        call move_alloc(source%name,destination%name)
        call move_alloc(source%source_group,destination%source_group)
        call move_alloc(source%parts,destination%parts)
        call move_alloc(source%connected_points,destination%connected_points)
        call move_alloc(source%connected_cells,destination%connected_cells)
    end subroutine move_studio_mesh

    subroutine build_studio_mesh(kind,variant,group,members,samples,mesh,error,native_group,catalog,patch_samples)
        integer,intent(in) :: kind,variant,group,members(:),samples(3)
        type(studio_mesh_t),intent(out) :: mesh
        character(:),allocatable,intent(out) :: error
        class(*),intent(in),optional :: native_group
        type(studio_patch_catalog), intent(in), optional :: catalog
        integer, intent(in), optional :: patch_samples(:,:)
        type(studio_mesh_t) :: candidate
        class(*), allocatable :: geometry
        integer, allocatable :: resolution(:,:)
        integer :: i
        logical :: sampled

        error='Select a curve, surface, volume or multipatch to mesh.'
        if(kind<1.or.kind>3.or.group<0.or.size(members)==0)return
        if(group==0)then
            if(size(members)/=1.or.members(1)/=variant)return
        else if(.not.present(native_group))then
            return
        end if
        if(.not.all([(single_patch_available(kind,members(i),catalog),i=1,size(members))]))return

        error='Enter at least two samples along each active direction for every patch.'
        resolution=spread(samples,2,size(members))
        if(present(patch_samples))then
            if(any(shape(patch_samples)/=shape(resolution)))return
            resolution=patch_samples
        end if
        if(any(resolution(:kind,:)<2))return
        resolution(kind+1:,:)=1
        error='The mesh exceeds the supported node or cell index range. Reduce the cell counts.'
        if(.not.mesh_samples_fit(kind,resolution))return

        allocate(candidate%parts(size(members)))
        do i=1,size(members)
            candidate%parts(i)%variant=members(i)
            candidate%parts(i)%revision=single_patch_mesh_revision(kind,members(i),catalog)
        end do
        sampled=.false.
        if(group==0)then
            call copy_single_patch_geometry(kind,variant,geometry,catalog)
            if(allocated(geometry))sampled=sample_forcad_patch(geometry,resolution(:,1),candidate%parts(1))
        else
            allocate(geometry,source=native_group)
            sampled=sample_forcad_assembly(geometry,kind,resolution,candidate%parts)
        end if
        if(.not.sampled)then
            error='ForCAD returned unusable sample or connectivity data. Check the geometry and sample counts.'
            if(allocated(geometry))then
                if(studio_native_problem(geometry))then
                    error='ForCAD reported a native mesh sampling problem. Review the source geometry and sample counts.'
                end if
            end if
            return
        end if
        mesh%kind=kind
        mesh%variant=variant
        mesh%group=group
        mesh%samples(:kind)=maxval(resolution(:kind,:),dim=2)
        do i=1,size(members)
            mesh%nodes=mesh%nodes+size(candidate%parts(i)%points,1)
            mesh%elements=mesh%elements+size(candidate%parts(i)%cells,1)
        end do
        call move_alloc(candidate%parts,mesh%parts)
        error=''
    end subroutine build_studio_mesh

    subroutine build_planned_volume_mesh(samples,joins,native_group,mesh,error)
        integer,intent(in) :: samples(3)
        type(multipatch_join),intent(in) :: joins(:)
        class(*),intent(in) :: native_group
        type(studio_mesh_t),intent(out) :: mesh
        character(:),allocatable,intent(out) :: error
        type(studio_mesh_t) :: candidate
        class(*), allocatable :: geometry
        integer :: count,i

        error='Select a connected ForCAD volume assembly to mesh.'
        count=volume_assembly_size(native_group,size(joins))
        if(count==0)return
        error='Enter at least two samples along each direction.'
        if(any(samples<2))return
        error='The mesh exceeds the supported node or cell index range. Reduce the cell counts.'
        if(.not.mesh_samples_fit(3,spread(samples,2,count)))return

        candidate%kind=3
        candidate%samples=samples
        allocate(geometry,source=native_group)
        allocate(candidate%parts(count))
        do i=1,count
            candidate%parts(i)%variant=i
        end do
        if(.not.sample_forcad_assembly(geometry,3,spread(samples,2,count),candidate%parts))then
            error='ForCAD returned unusable sample or connectivity data for the planned volume assembly.'
            if(studio_native_problem(geometry))then
                error='ForCAD reported a native sampling problem while meshing the planned volume assembly.'
            end if
            return
        end if
        call stitch_studio_volume_mesh(candidate,joins,error)
        if(len(error)>0)return
        call move_studio_mesh(candidate,mesh)
    end subroutine build_planned_volume_mesh

    subroutine stitch_studio_volume_mesh(mesh,joins,error)
        type(studio_mesh_t),intent(inout) :: mesh
        type(multipatch_join),intent(in) :: joins(:)
        character(:),allocatable,intent(out) :: error
        real(rk),allocatable :: points(:,:)
        integer,allocatable :: cells(:,:)

        error='The mesh needs a sampled volume assembly with validated joins.'
        if(mesh%kind/=3.or..not.allocated(mesh%parts).or.size(joins)==0)return
        if(size(mesh%parts)<2)return

        call join_forcad_volume_mesh(mesh%parts%variant,mesh%parts,joins,points,cells,error)
        if(len(error)>0)return
        mesh%nodes=size(points,1)
        mesh%elements=size(cells,1)
        call move_alloc(points,mesh%connected_points)
        call move_alloc(cells,mesh%connected_cells)
    end subroutine stitch_studio_volume_mesh

end module studio_mesh
