module studio_forcad_mesh
    use, intrinsic :: iso_fortran_env, only: int64
    use forcad, only: rk,nurbs_curve,nurbs_surface,nurbs_volume, &
        nurbs_multipatch_curve,nurbs_multipatch_surface,nurbs_multipatch_volume
    use forcad_utils, only: boundary_axes,boundary_index,disjoint_set_union,disjoint_set_map
    use forcad_multipatch, only: interface_sampling
    use studio_forcad_assembly, only: multipatch_join,initialize_assembly,append_assembly_patch
    implicit none(type, external)
    private

    public :: forcad_mesh_samples, sample_forcad_patch, sample_forcad_assembly
    public :: volume_assembly_size, check_forcad_mesh_joins, join_forcad_volume_mesh, estimate_forcad_mesh_samples

    type :: forcad_mesh_samples
        integer :: variant=0, samples(3)=1
        integer(int64) :: revision=0_int64
        real(rk), allocatable :: points(:,:)   ! Native [node, coordinate] layout.
        integer, allocatable :: cells(:,:)    ! Native one-based [cell, local node] layout.
    end type forcad_mesh_samples

contains

    integer function volume_assembly_size(geometry,connections) result(count)
        class(*),intent(in) :: geometry
        integer,intent(in) :: connections

        count=0
        if(connections<1)return
        select type(geometry)
        type is(nurbs_multipatch_volume)
            if(geometry%get_npatch()<2.or.geometry%get_nconnection()/=connections)return
            count=geometry%get_npatch()
        end select
    end function volume_assembly_size

    logical function sample_forcad_patch(geometry,samples,part) result(valid)
        class(*),intent(inout) :: geometry
        integer,intent(in) :: samples(3)
        type(forcad_mesh_samples),intent(inout) :: part

        select type(geometry)
        type is(nurbs_curve)
            call geometry%create(res=samples(1))
        type is(nurbs_surface)
            call geometry%create(res1=samples(1),res2=samples(2))
        type is(nurbs_volume)
            call geometry%create(res1=samples(1),res2=samples(2),res3=samples(3))
        end select
        valid=read_mesh_part(geometry,part)
        if(valid)part%samples=samples
    end function sample_forcad_patch

    logical function sample_forcad_assembly(geometry,kind,samples,parts) result(valid)
        class(*),intent(inout) :: geometry
        integer,intent(in) :: kind
        integer,intent(in),contiguous :: samples(:,:)
        type(forcad_mesh_samples),intent(inout) :: parts(:)
        integer :: i

        valid=.false.
        if(any(shape(samples)/=[3,size(parts)]))return
        select type(geometry)
        type is(nurbs_multipatch_curve)
            if(kind/=1.or.geometry%get_npatch()/=size(parts))return
            call geometry%create(res=samples(1,:))
            if(.not.geometry%err%ok)return
            do i=1,size(parts)
                if(.not.read_mesh_part(geometry%get_patch(i),parts(i)))return
            end do
        type is(nurbs_multipatch_surface)
            if(kind/=2.or.geometry%get_npatch()/=size(parts))return
            call geometry%create(res=samples(:2,:),synchronize=.false.)
            if(.not.geometry%err%ok)return
            do i=1,size(parts)
                if(.not.read_mesh_part(geometry%get_patch(i),parts(i)))return
            end do
        type is(nurbs_multipatch_volume)
            if(kind/=3.or.geometry%get_npatch()/=size(parts))return
            call geometry%create(res=samples,synchronize=.false.)
            if(.not.geometry%err%ok)return
            do i=1,size(parts)
                if(.not.read_mesh_part(geometry%get_patch(i),parts(i)))return
            end do
        class default
            return
        end select
        do i=1,size(parts)
            parts(i)%samples=samples(:,i)
        end do
        valid=.true.
    end function sample_forcad_assembly

    subroutine estimate_forcad_mesh_samples(geometry,kind,target_size,synchronize,probe_limit,samples,error)
        class(*), intent(in) :: geometry
        integer, intent(in) :: kind,probe_limit
        real(rk), intent(in) :: target_size
        logical, intent(in) :: synchronize
        integer, allocatable, intent(out) :: samples(:,:)
        character(:), allocatable, intent(out) :: error
        class(*), allocatable :: group
        integer, allocatable :: requested(:,:),resolved(:,:)
        integer :: i,axis,id,counts(3)

        error='The source is unavailable for target-size mesh sampling.'
        if(kind<1.or.kind>3.or.probe_limit<3)return
        select type(geometry)
        type is(nurbs_multipatch_curve)
            if(kind/=1)return
            allocate(group,source=geometry)
        type is(nurbs_multipatch_surface)
            if(kind/=2)return
            allocate(group,source=geometry)
        type is(nurbs_multipatch_volume)
            if(kind/=3)return
            allocate(group,source=geometry)
        class default
            call initialize_assembly(kind,group)
            if(.not.allocated(group))return
            call append_assembly_patch(group,geometry,id)
            if(id/=1)return
        end select
        error='ForCAD could not sample the selected geometry at the requested target size.'
        select type(group)
        type is(nurbs_multipatch_curve)
            call group%create(target_size=target_size,max_res=probe_limit)
            if(.not.group%err%ok)return
            allocate(requested(1,group%get_npatch()))
            do i=1,size(requested,2)
                associate(patch=>group%get_patch(i))
                    requested(1,i)=patch%get_ng()
                end associate
            end do
            resolved=requested
        type is(nurbs_multipatch_surface)
            allocate(requested(2,group%get_npatch()))
            do axis=1,2
                counts=2
                counts(axis)=probe_limit
                call group%create(target_size=target_size,max_res=counts(:2),synchronize=.false.)
                if(.not.group%err%ok)return
                do i=1,size(requested,2)
                    associate(patch=>group%get_patch(i))
                        counts(:2)=patch%get_ng()
                    end associate
                    requested(axis,i)=counts(axis)
                end do
            end do
            call interface_sampling(2,size(requested,2),[(group%get_connection(i),i=1,group%get_nconnection())], &
                group%get_nconnection(),requested,synchronize,resolved,group%err)
            if(.not.group%err%ok)return
        type is(nurbs_multipatch_volume)
            allocate(requested(3,group%get_npatch()))
            do axis=1,3
                counts=2
                counts(axis)=probe_limit
                call group%create(target_size=target_size,max_res=counts,synchronize=.false.)
                if(.not.group%err%ok)return
                do i=1,size(requested,2)
                    associate(patch=>group%get_patch(i))
                        counts=patch%get_ng()
                    end associate
                    requested(axis,i)=counts(axis)
                end do
            end do
            call interface_sampling(3,size(requested,2),[(group%get_connection(i),i=1,group%get_nconnection())], &
                group%get_nconnection(),requested,synchronize,resolved,group%err)
            if(.not.group%err%ok)return
        class default
            return
        end select
        if(any(resolved>=probe_limit))then
            error='The target size exceeds the preview work budget. Increase the size or enter counts directly.'
            return
        end if
        allocate(samples(3,size(resolved,2)),source=1)
        samples(:kind,:)=resolved
        error=''
    end subroutine estimate_forcad_mesh_samples

    logical function read_mesh_part(geometry,part) result(valid)
        class(*),intent(in) :: geometry
        type(forcad_mesh_samples),intent(inout) :: part
        integer :: width

        valid=.false.
        select type(geometry)
        type is(nurbs_curve)
            if(.not.geometry%err%ok)return
            part%points=geometry%get_Xg()
            part%cells=geometry%get_elem_Xg_vis()
            width=2
        type is(nurbs_surface)
            if(.not.geometry%err%ok)return
            part%points=geometry%get_Xg()
            part%cells=geometry%get_elem_Xg_vis()
            width=4
        type is(nurbs_volume)
            if(.not.geometry%err%ok)return
            part%points=geometry%get_Xg()
            part%cells=geometry%get_elem_Xg_vis()
            width=8
        class default
            return
        end select
        if(size(part%points,2)<1.or.size(part%points,2)>3.or.size(part%cells,2)/=width)return
        if(size(part%cells,1)==0)return
        valid=all(part%cells>=1).and.all(part%cells<=size(part%points,1))
    end function read_mesh_part

    subroutine check_forcad_mesh_joins(samples,part_ids,joins,error)
        integer,intent(in) :: samples(:,:),part_ids(:)
        type(multipatch_join),intent(in) :: joins(:)
        character(:),allocatable,intent(out) :: error
        integer :: join,a,b,axes_a(3),axes_b(3)
        character(16) :: number

        error='The mesh needs a sampled volume assembly with validated joins.'
        if(size(part_ids)<2.or.any(shape(samples)/=[3,size(part_ids)]).or.size(joins)==0.or.any(samples<2))return
        do join=1,size(joins)
            associate(seam=>joins(join))
                a=findloc(part_ids,seam%patch_a,dim=1)
                b=findloc(part_ids,seam%patch_b,dim=1)
                error='A mesh join does not identify two valid volume faces.'
                if(min(a,b)<1.or.min(seam%side_a,seam%side_b)<1.or.max(seam%side_a,seam%side_b)>6)return
                if(a==b.and.seam%side_a==seam%side_b)return
                axes_a=boundary_axes(seam%side_a)
                axes_b=boundary_axes(seam%side_b,seam%swap,seam%flip,seam%reverse)
                if(any(samples(axes_a(:2),a)/=samples(abs(axes_b(:2)),b)))then
                    write(number,'(i0)')join
                    error='Interface '//trim(number)//' needs matching cell counts along corresponding face directions.'
                    return
                end if
            end associate
        end do
        error=''
    end subroutine check_forcad_mesh_joins

    subroutine join_forcad_volume_mesh(part_ids,parts,joins,connected_points,connected_cells,error)
        integer,intent(in) :: part_ids(:)
        type(forcad_mesh_samples),intent(in) :: parts(:)
        type(multipatch_join),intent(in) :: joins(:)
        character(:),allocatable,intent(out) :: error
        integer,allocatable,intent(out) :: connected_cells(:,:)
        real(rk),allocatable,intent(out) :: connected_points(:,:)
        integer,allocatable :: parent(:),node_map(:),face_a(:),face_b(:),cells(:,:)
        integer :: node_offsets(0:size(parts)),cell_offsets(0:size(parts))
        real(rk),allocatable :: points(:,:)
        integer :: sample_nodes,part_cells,part,i,j,join,a,b
        integer :: offset_node,offset_cell,local,used(8),filled,global
        real(rk) :: scale,delta

        call check_forcad_mesh_joins(reshape([(parts(part)%samples,part=1,size(parts))],[3,size(parts)]),part_ids,joins,error)
        if(len(error)>0)return
        error='The mesh does not contain the expected native sample data.'
        node_offsets(0)=0
        cell_offsets(0)=0
        do part=1,size(parts)
            associate(item=>parts(part))
                if(.not.allocated(item%points).or..not.allocated(item%cells))return
                sample_nodes=size(item%points,1)
                part_cells=size(item%cells,1)
                if(size(item%points,2)/=3.or.size(item%cells,2)/=8)return
                if(product(int(item%samples,int64))/=sample_nodes.or.product(int(item%samples-1,int64))/=part_cells)return
                if(any(item%cells<1).or.any(item%cells>sample_nodes))return
                if(sample_nodes>huge(1)-node_offsets(part-1).or.part_cells>huge(1)-cell_offsets(part-1))return
                node_offsets(part)=node_offsets(part-1)+sample_nodes
                cell_offsets(part)=cell_offsets(part-1)+part_cells
            end associate
        end do

        parent=[(i,i=1,node_offsets(size(parts)))]
        do join=1,size(joins)
            associate(seam=>joins(join))
                a=findloc(part_ids,seam%patch_a,dim=1)
                b=findloc(part_ids,seam%patch_b,dim=1)
                call boundary_index(parts(a)%samples,seam%side_a,.false.,.false.,[.false.,.false.],face_a)
                call boundary_index(parts(b)%samples,seam%side_b,seam%reverse,seam%swap,seam%flip,face_b)
                error='ForCAD mesh samples do not coincide across a planned join.'
                scale=max(maxval(abs(parts(a)%points)),maxval(abs(parts(b)%points)),tiny(1._rk))
                delta=maxval(abs(parts(a)%points(face_a,:)-parts(b)%points(face_b,:)))
                if(delta>512._rk*epsilon(1._rk)*scale)return
                do i=1,size(face_a)
                    call disjoint_set_union(parent,node_offsets(a-1)+face_a(i),node_offsets(b-1)+face_b(i))
                end do
            end associate
        end do
        call disjoint_set_map(parent,node_map)

        allocate(points(maxval(node_map),3),cells(cell_offsets(size(parts)),8))
        filled=0
        do part=1,size(parts)
            offset_node=node_offsets(part-1)
            offset_cell=cell_offsets(part-1)
            associate(item=>parts(part))
                do local=1,size(item%points,1)
                    global=node_map(offset_node+local)
                    if(global<=filled)cycle
                    points(global,:)=item%points(local,:)
                    filled=global
                end do
                error='A joined hexahedron has repeated node indices.'
                do i=1,size(item%cells,1)
                    used=node_map(offset_node+item%cells(i,:))
                    do j=2,8
                        if(any(used(:j-1)==used(j)))return
                    end do
                    cells(offset_cell+i,:)=used
                end do
            end associate
        end do
        call move_alloc(points,connected_points)
        call move_alloc(cells,connected_cells)
        error=''
    end subroutine join_forcad_volume_mesh

end module studio_forcad_mesh
