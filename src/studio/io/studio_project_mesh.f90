module studio_project_mesh
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use studio_mesh, only: studio_mesh_t
    use studio_forcad, only: rk
    use studio_single_patch, only: studio_patch_catalog,studio_patch_is_staged,single_patch_mesh_revision,next_patch_revision
    use studio_project_stream, only: write_text,read_text,remaining_file_units,integer_units,real_units
    implicit none(type, external)
    private
    public :: write_project_meshes,read_project_meshes,rebind_project_meshes
    public :: write_session_mesh,read_session_mesh

    integer(int64), parameter :: stale_revision=-2_int64, current_part=1_int64

contains

    subroutine write_project_meshes(unit,meshes,catalog,status)
        integer, intent(in) :: unit
        type(studio_mesh_t), allocatable, intent(in), optional :: meshes(:)
        type(studio_patch_catalog), intent(in), optional :: catalog
        integer, intent(out) :: status
        integer :: count,i

        count=0
        if(present(meshes))then
            if(allocated(meshes))count=size(meshes)
        end if
        write(unit,iostat=status)count
        if(status/=0)return
        do i=1,count
            call write_mesh(unit,meshes(i),catalog,status)
            if(status/=0)return
        end do
    end subroutine write_project_meshes

    subroutine write_session_mesh(unit,mesh,status)
        integer, intent(in) :: unit
        type(studio_mesh_t), intent(in) :: mesh
        integer, intent(out) :: status
        integer :: labels(2)

        labels=merge(1,0,[allocated(mesh%name),allocated(mesh%source_group)])
        write(unit,iostat=status)mesh%construction_stamp,mesh%content_revision,mesh%nodes,mesh%elements,labels
        if(status/=0)return
        call write_mesh(unit,mesh,status=status,labels=labels)
    end subroutine write_session_mesh

    subroutine write_mesh(unit,mesh,catalog,status,labels)
        integer, intent(in) :: unit
        type(studio_mesh_t), intent(in) :: mesh
        type(studio_patch_catalog), intent(in), optional :: catalog
        integer, intent(out) :: status
        integer, intent(in), optional :: labels(2)
        integer :: flags(4),i,current,part_count,names(2)
        integer(int64) :: revision
        logical :: connected

        status=1
        part_count=-1
        if(allocated(mesh%parts))part_count=size(mesh%parts)
        names=[1,merge(1,0,mesh%group>0)]
        if(present(labels))then
            names=labels
        else
            if(part_count<1)return
        end if
        connected=allocated(mesh%connected_cells)
        if(connected.neqv.allocated(mesh%connected_points))return
        flags=merge(1,0,[mesh%visible,mesh%per_patch_sizing,mesh%synchronize,connected])
        write(unit,iostat=status)mesh%kind,mesh%variant,mesh%group,mesh%samples,flags,part_count,mesh%target_size
        if(status/=0)return
        if(names(1)==1)then
            status=1
            if(.not.allocated(mesh%name))return
            call write_text(unit,mesh%name,status)
            if(status/=0)return
        end if
        if(names(2)==1)then
            status=1
            if(.not.allocated(mesh%source_group))return
            call write_text(unit,mesh%source_group,status)
            if(status/=0)return
        end if

        do i=1,part_count
            associate(part=>mesh%parts(i))
                status=1
                if(.not.allocated(part%points).or..not.allocated(part%cells))return
                if(present(labels))then
                    write(unit,iostat=status)part%variant,part%samples,part%revision
                else
                    revision=single_patch_mesh_revision(mesh%kind,part%variant,catalog)
                    current=merge(1,0,revision>=0_int64.and.part%revision==revision)
                    write(unit,iostat=status)part%variant,part%samples,current
                end if
                if(status/=0)return
                call write_mesh_arrays(unit,part%points,part%cells,status)
                if(status/=0)return
            end associate
        end do
        if(connected)call write_mesh_arrays(unit,mesh%connected_points,mesh%connected_cells,status)
    end subroutine write_mesh

    subroutine write_mesh_arrays(unit,points,cells,status)
        integer, intent(in) :: unit,cells(:,:)
        real(rk), intent(in) :: points(:,:)
        integer, intent(out) :: status
        integer :: io

        write(unit,iostat=io)shape(points),size(cells,1),points,cells
        status=merge(0,1,io==0)
    end subroutine write_mesh_arrays

    subroutine read_project_meshes(unit,catalog,meshes,status)
        integer, intent(in) :: unit
        type(studio_patch_catalog), intent(in) :: catalog
        type(studio_mesh_t), allocatable, intent(out) :: meshes(:)
        integer, intent(out) :: status
        integer :: count,i

        read(unit,iostat=status)count
        if(status/=0)return
        status=1
        if(count<0.or.11_int64*count>remaining_file_units(unit)/integer_units)return
        if(count==0)then
            status=0
            return
        end if
        allocate(meshes(count),stat=status)
        if(status/=0)return
        do i=1,count
            call read_mesh(unit,catalog,meshes(i),status)
            if(status/=0)return
        end do
    end subroutine read_project_meshes

    subroutine read_session_mesh(unit,mesh,status)
        integer, intent(in) :: unit
        type(studio_mesh_t), intent(out) :: mesh
        integer, intent(out) :: status
        integer(int64) :: construction_stamp,content_revision
        integer :: nodes,elements,labels(2)

        read(unit,iostat=status)construction_stamp,content_revision,nodes,elements,labels
        if(status/=0)return
        status=1
        if(any(labels<0).or.any(labels>1))return
        call read_mesh(unit,mesh=mesh,status=status,labels=labels)
        if(status/=0)return
        mesh%construction_stamp=construction_stamp
        mesh%content_revision=content_revision
        mesh%nodes=nodes
        mesh%elements=elements
    end subroutine read_session_mesh

    subroutine read_mesh(unit,catalog,mesh,status,labels)
        integer, intent(in) :: unit
        type(studio_patch_catalog), intent(in), optional :: catalog
        type(studio_mesh_t), intent(out) :: mesh
        integer, intent(out) :: status
        integer, intent(in), optional :: labels(2)
        integer :: flags(4),part_count,i,width,current,names(2)
        integer(int64) :: nodes,elements

        read(unit,iostat=status)mesh%kind,mesh%variant,mesh%group,mesh%samples,flags,part_count,mesh%target_size
        if(status/=0)return
        status=1
        if(any(flags<0).or.any(flags>1))return
        if(.not.ieee_is_finite(mesh%target_size).or.mesh%target_size<0._rk)return
        if(present(labels))then
            names=labels
            if(part_count< -1)return
        else
            if(part_count<1)return
            names=[1,merge(1,0,mesh%group>0)]
        end if
        if(part_count>0)then
            if(mesh%kind<1.or.mesh%kind>3.or.mesh%group<0.or.mesh%variant<0)return
            if(any(mesh%samples(:mesh%kind)<2).or.any(mesh%samples<0))return
            if((mesh%group==0).neqv.(mesh%variant>0))return
            if(mesh%group==0.and.part_count/=1)return
        end if
        if(flags(4)==1.and.(mesh%kind/=3.or.mesh%group==0.or.part_count<2))return
        if(8_int64*part_count>remaining_file_units(unit)/integer_units)return
        mesh%visible=flags(1)==1
        mesh%per_patch_sizing=flags(2)==1
        mesh%synchronize=flags(3)==1

        if(names(1)==1)then
            call read_text(unit,mesh%name,status)
            if(status/=0)return
        end if
        if(names(2)==1)then
            call read_text(unit,mesh%source_group,status)
            if(status/=0)return
        end if
        if(part_count<0)then
            status=0
            return
        end if
        allocate(mesh%parts(part_count),stat=status)
        if(status/=0.or.part_count==0)return
        width=2**mesh%kind
        nodes=0
        elements=0
        do i=1,part_count
            associate(part=>mesh%parts(i))
                if(present(labels))then
                    read(unit,iostat=status)part%variant,part%samples,part%revision
                else
                    read(unit,iostat=status)part%variant,part%samples,current
                end if
                if(status/=0)return
                status=1
                if(part%variant<1.or.any(part%samples(:mesh%kind)<2).or.any(part%samples(mesh%kind+1:)/=1))return
                if(mesh%group==0.and.mesh%variant/=part%variant)return
                if(.not.present(labels))then
                    if(current<0.or.current>1)return
                    part%revision=stale_revision
                    if(current==1)then
                        if(.not.studio_patch_is_staged(catalog,mesh%kind,part%variant))return
                        part%revision=current_part
                    end if
                end if
                call read_mesh_arrays(unit,width,part%points,part%cells,status,part%samples(:mesh%kind))
                if(status/=0)return
                nodes=nodes+size(part%points,1)
                elements=elements+size(part%cells,1)
            end associate
        end do

        status=1
        if(nodes>huge(mesh%nodes).or.elements>huge(mesh%elements))return
        mesh%nodes=int(nodes)
        mesh%elements=int(elements)
        if(flags(4)==1)then
            call read_mesh_arrays(unit,width,mesh%connected_points,mesh%connected_cells,status)
            if(status/=0)return
            status=1
            if(size(mesh%connected_points,2)/=3.or.size(mesh%connected_points,1)>nodes)return
            if(size(mesh%connected_cells,1)/=elements)return
            mesh%nodes=size(mesh%connected_points,1)
        end if
        status=0
    end subroutine read_mesh

    subroutine read_mesh_arrays(unit,width,points,cells,status,samples)
        integer, intent(in) :: unit,width
        real(rk), allocatable, intent(out) :: points(:,:)
        integer, allocatable, intent(out) :: cells(:,:)
        integer, intent(out) :: status
        integer, intent(in), optional :: samples(:)
        integer :: nodes,coordinates,elements
        integer(int64) :: units

        read(unit,iostat=status)nodes,coordinates,elements
        if(status/=0)return
        status=1
        if(nodes<1.or.coordinates<1.or.coordinates>3.or.elements<1)return
        if(present(samples))then
            if(.not.grid_has_size(samples,nodes).or..not.grid_has_size(samples-1,elements))return
        end if
        units=int(nodes,int64)*coordinates*real_units+int(elements,int64)*width*integer_units
        if(units>remaining_file_units(unit))return
        allocate(points(nodes,coordinates),stat=status)
        if(status/=0)return
        allocate(cells(elements,width),stat=status)
        if(status/=0)return
        read(unit,iostat=status)points,cells
        if(status/=0)return
        if(.not.all(ieee_is_finite(points)).or.any(cells<1).or.any(cells>nodes))status=1
    end subroutine read_mesh_arrays

    pure logical function grid_has_size(dimensions,entries) result(matches)
        integer, intent(in) :: dimensions(:),entries
        integer :: remaining,i

        matches=.false.
        remaining=entries
        do i=1,size(dimensions)
            if(mod(remaining,dimensions(i))/=0)return
            remaining=remaining/dimensions(i)
        end do
        matches=remaining==1
    end function grid_has_size

    subroutine rebind_project_meshes(meshes,catalog)
        type(studio_mesh_t), intent(inout) :: meshes(:)
        type(studio_patch_catalog), intent(in), optional :: catalog
        integer :: i,j

        do i=1,size(meshes)
            meshes(i)%content_revision=next_patch_revision()
            do j=1,size(meshes(i)%parts)
                associate(part=>meshes(i)%parts(j))
                    if(part%revision/=current_part)cycle
                    part%revision=single_patch_mesh_revision(meshes(i)%kind,part%variant,catalog)
                end associate
            end do
        end do
    end subroutine rebind_project_meshes

end module studio_project_mesh
