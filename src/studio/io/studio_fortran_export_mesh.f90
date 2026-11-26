module studio_fortran_export_mesh
    use studio_document, only: studio_document_t,document_mesh_count,document_mesh_valid,document_mesh_source_group
    use studio_fortran_export_text, only: source_writer
    use studio_mesh, only: studio_mesh_t
    use studio_presentation_text, only: integer_text
    use studio_forcad, only: rk
    implicit none(type,external)
    private
    public :: write_mesh_calls,write_mesh_routines

    character(8),parameter :: patch_arrays(3)=[character(8)::'curves','surfaces','volumes']
    character(18),parameter :: assembly_arrays(3)=[character(18)::'curve_assemblies','surface_assemblies','volume_assemblies']
    character(7),parameter :: geometry_ranks(3)=[character(7)::'curve','surface','volume']

contains

    subroutine write_mesh_calls(source,document)
        type(source_writer),intent(inout) :: source
        type(studio_document_t),intent(in) :: document
        integer :: i

        do i=1,document_mesh_count(document)
            if(source%status/=0)return
            call source%line('    call build_mesh_'//integer_text(i)//'()')
        end do
    end subroutine write_mesh_calls

    subroutine write_mesh_routines(source,document)
        type(source_writer),intent(inout) :: source
        type(studio_document_t),intent(in) :: document
        logical :: sampled,conforming
        integer :: i,group

        sampled=.false.
        conforming=.false.
        do i=1,document_mesh_count(document)
            if(source%status/=0)return
            associate(mesh=>document%meshes(i))
                if(document_mesh_valid(document,i))then
                    group=document_mesh_source_group(document,i)
                    call write_mesh_recipe(source,mesh,i,group)
                    sampled=.true.
                    conforming=conforming.or.allocated(mesh%connected_cells)
                else
                    call write_mesh_snapshot(source,mesh,i)
                end if
            end associate
        end do
        if(sampled)call write_sampling_helpers(source,conforming)
        if(conforming)call write_stitching_helper(source)
        if(document_mesh_count(document)>0)call write_vtk_helper(source)
    end subroutine write_mesh_routines

    subroutine write_mesh_recipe(source,mesh,id,group)
        type(source_writer),intent(inout) :: source
        type(studio_mesh_t),intent(in) :: mesh
        integer,intent(in) :: id,group
        character(:),allocatable :: number,geometry,part_count,stitch
        integer :: i

        number=integer_text(id)
        part_count=integer_text(size(mesh%parts))
        if(group>0)then
            geometry=trim(assembly_arrays(mesh%kind))//'('//integer_text(group)//')'
        else
            geometry=trim(patch_arrays(mesh%kind))//'('//integer_text(mesh%variant)//')'
        end if
        stitch=trim(merge('.true. ','.false.',allocated(mesh%connected_cells)))

        call source%line('')
        call source%line('    subroutine build_mesh_'//number//'()')
        call source%line('        integer :: samples('//integer_text(3*size(mesh%parts))//')')
        call source%integers('samples',[(mesh%parts(i)%samples,i=1,size(mesh%parts))])
        call source%line("        call sample_mesh('mesh_"//number//".vtk', "//geometry//', '// &
            integer_text(mesh%kind)//', reshape(samples,[3,'//part_count//']), '//stitch//')')
        call source%line('    end subroutine build_mesh_'//number)
    end subroutine write_mesh_recipe

    subroutine write_mesh_snapshot(source,mesh,id)
        type(source_writer),intent(inout) :: source
        type(studio_mesh_t),intent(in) :: mesh
        integer,intent(in) :: id
        character(:),allocatable :: number,part_count,point_shape,cell_shape
        real(rk),allocatable :: points(:,:)
        integer,allocatable :: cells(:,:),cell_counts(:)
        integer :: nodes,elements,width,i,node_offset,element_offset,num_nodes,num_elements

        if(.not.allocated(mesh%parts))then
            source%status=1
            return
        end if
        if(size(mesh%parts)==0.or.mesh%kind<1.or.mesh%kind>3)then
            source%status=1
            return
        end if
        if(allocated(mesh%connected_points).neqv.allocated(mesh%connected_cells))then
            source%status=1
            return
        end if
        do i=1,size(mesh%parts)
            if(.not.allocated(mesh%parts(i)%points).or..not.allocated(mesh%parts(i)%cells))then
                source%status=1
                return
            end if
        end do

        number=integer_text(id)
        part_count=integer_text(size(mesh%parts))
        width=2**mesh%kind
        cell_counts=[(size(mesh%parts(i)%cells,1),i=1,size(mesh%parts))]
        if(allocated(mesh%connected_points))then
            points=mesh%connected_points
            cells=mesh%connected_cells
        else
            nodes=sum([(size(mesh%parts(i)%points,1),i=1,size(mesh%parts))])
            elements=sum(cell_counts)
            allocate(points(nodes,3),source=0._rk)
            allocate(cells(elements,width))
            node_offset=0
            element_offset=0
            do i=1,size(mesh%parts)
                associate(part=>mesh%parts(i))
                    num_nodes=size(part%points,1)
                    num_elements=cell_counts(i)
                    points(node_offset+1:node_offset+num_nodes,:size(part%points,2))=part%points
                    cells(element_offset+1:element_offset+num_elements,:)=part%cells+node_offset
                    node_offset=node_offset+num_nodes
                    element_offset=element_offset+num_elements
                end associate
            end do
        end if

        point_shape=integer_text(size(points,1))//','//integer_text(size(points,2))
        cell_shape=integer_text(size(cells,1))//','//integer_text(width)
        call source%line('')
        call source%line('    subroutine build_mesh_'//number//'()')
        call source%line('        real(rk), allocatable :: coordinates(:)')
        call source%line('        integer, allocatable :: connectivity(:)')
        call source%line('        integer :: cell_counts('//part_count//')')
        call source%line('        allocate(coordinates('//integer_text(size(points))//'), '// &
            'connectivity('//integer_text(size(cells))//'))')
        call source%reals('coordinates',reshape(points,[size(points)]))
        call source%integers('connectivity',reshape(cells,[size(cells)]))
        call source%integers('cell_counts',cell_counts)
        call source%line("        call write_mesh_vtk('mesh_"//number//".vtk', "//integer_text(mesh%kind)//', &')
        call source%line('            reshape(coordinates,['//point_shape//']), reshape(connectivity,['//cell_shape//']), cell_counts)')
        call source%line('    end subroutine build_mesh_'//number)
    end subroutine write_mesh_snapshot

    subroutine write_sampling_helpers(source,conforming)
        type(source_writer),intent(inout) :: source
        logical,intent(in) :: conforming
        character(:),allocatable :: rank_type,kind_arg
        character(40),parameter :: assembly_resolution(3)=[character(40):: &
            'res=samples(1,:)','res=samples(:2,:), synchronize=.false.','res=samples, synchronize=.false.']
        character(6),parameter :: resolution(3)=[character(6)::'res=','res1=','res1=']
        integer :: kind

        call write_lines(source,[character(132):: &
            '', &
            '    subroutine sample_mesh(filename, geometry, kind, samples, stitch)', &
            '        character(*), intent(in) :: filename', &
            '        class(*), intent(in) :: geometry', &
            '        integer, intent(in) :: kind, samples(:,:)', &
            '        logical, intent(in) :: stitch', &
            '        real(rk), allocatable :: points(:,:)', &
            '        integer, allocatable :: cells(:,:)', &
            '        integer :: node_offsets(0:size(samples,2)), cell_offsets(0:size(samples,2)), i', &
            '', &
            '        node_offsets(0)=0', &
            '        cell_offsets(0)=0', &
            '        do i=1,size(samples,2)', &
            '            node_offsets(i)=node_offsets(i-1)+product(samples(:kind,i))', &
            '            cell_offsets(i)=cell_offsets(i-1)+product(samples(:kind,i)-1)', &
            '        end do', &
            '        allocate(points(node_offsets(size(samples,2)),3),source=0._rk)', &
            '        allocate(cells(cell_offsets(size(samples,2)),2**kind))', &
            '        select type(geometry)'])
        do kind=1,3
            rank_type=trim(geometry_ranks(kind))
            call source%line('        type is(nurbs_'//rank_type//')')
            call source%line('            block')
            call source%line('                type(nurbs_'//rank_type//') :: patch')
            call source%line('                patch=geometry')
            kind_arg=trim(resolution(kind))//'samples(1,1)'
            if(kind>=2)kind_arg=kind_arg//', res2=samples(2,1)'
            if(kind==3)kind_arg=kind_arg//', res3=samples(3,1)'
            call source%line('                call patch%create('//kind_arg//')')
            call write_lines(source,[character(132):: &
                "                if(.not.patch%err%ok)error stop 'ForCAD could not sample a patch mesh.'", &
                '                call store_mesh_part(points,cells,node_offsets,cell_offsets,1, &', &
                '                    patch%get_Xg(),patch%get_elem_Xg_vis())', &
                '            end block'])
        end do
        do kind=1,3
            rank_type=trim(geometry_ranks(kind))
            call source%line('        type is(nurbs_multipatch_'//rank_type//')')
            call source%line('            block')
            call source%line('                type(nurbs_multipatch_'//rank_type//') :: assembly')
            call source%line('                assembly=geometry')
            call source%line('                call assembly%create('//trim(assembly_resolution(kind))//')')
            call write_lines(source,[character(132):: &
                "                if(.not.assembly%err%ok)error stop 'ForCAD could not sample a multipatch mesh.'", &
                '                do i=1,assembly%get_npatch()', &
                '                    associate(patch=>assembly%get_patch(i))', &
                '                        call store_mesh_part(points,cells,node_offsets,cell_offsets,i, &', &
                '                            patch%get_Xg(),patch%get_elem_Xg_vis())', &
                '                    end associate', &
                '                end do'])
            if(kind==3.and.conforming)then
                call write_lines(source,[character(132):: &
                    '                if(stitch)call stitch_mesh(points,cells,node_offsets,samples, &', &
                    '                    [(assembly%get_connection(i),i=1,assembly%get_nconnection())])'])
            end if
            call source%line('            end block')
        end do
        call write_lines(source,[character(132):: &
            '        class default', &
            "            error stop 'Unsupported ForCAD mesh source.'", &
            '        end select', &
            '        call write_mesh_vtk(filename,kind,points,cells,cell_offsets(1:)-cell_offsets(:size(samples,2)-1))', &
            '    end subroutine sample_mesh', &
            '', &
            '    subroutine store_mesh_part(points,cells,node_offsets,cell_offsets,i,local_points,local_cells)', &
            '        real(rk), intent(inout) :: points(:,:)', &
            '        integer, intent(inout) :: cells(:,:)', &
            '        integer, intent(in) :: node_offsets(0:),cell_offsets(0:),i,local_cells(:,:)', &
            '        real(rk), intent(in) :: local_points(:,:)', &
            '', &
            '        if(size(local_points,1)/=node_offsets(i)-node_offsets(i-1).or. &', &
            '            size(local_cells,1)/=cell_offsets(i)-cell_offsets(i-1).or. &', &
            '            size(local_points,2)>size(points,2).or.size(local_cells,2)/=size(cells,2)) &', &
            "            error stop 'ForCAD returned unexpected mesh dimensions.'", &
            '        points(node_offsets(i-1)+1:node_offsets(i),:size(local_points,2))=local_points', &
            '        cells(cell_offsets(i-1)+1:cell_offsets(i),:)=local_cells+node_offsets(i-1)', &
            '    end subroutine store_mesh_part'])
    end subroutine write_sampling_helpers

    subroutine write_stitching_helper(source)
        type(source_writer),intent(inout) :: source

        call write_lines(source,[character(132):: &
            '', &
            '    subroutine stitch_mesh(points,cells,node_offsets,samples,connections)', &
            '        real(rk), allocatable, intent(inout) :: points(:,:)', &
            '        integer, intent(inout) :: cells(:,:)', &
            '        integer, intent(in) :: node_offsets(0:),samples(:,:)', &
            '        type(multipatch_connection), intent(in) :: connections(:)', &
            '        real(rk), allocatable :: shared_points(:,:)', &
            '        integer, allocatable :: parent(:),node_map(:),face_a(:),face_b(:)', &
            '        integer :: i,j,a,b,filled,node', &
            '', &
            '        parent=[(i,i=1,size(points,1))]', &
            '        do i=1,size(connections)', &
            '            associate(join=>connections(i))', &
            '                a=join%get_patch_a()', &
            '                b=join%get_patch_b()', &
            '                call boundary_index(samples(:,a),join%get_side_a(),.false.,.false.,[.false.,.false.],face_a)', &
            '                call boundary_index(samples(:,b),join%get_side_b(),join%is_reversed(),join%is_swapped(), &', &
            '                    [join%is_flipped(1),join%is_flipped(2)],face_b)', &
            "                if(size(face_a)/=size(face_b))error stop 'ForCAD returned incompatible mesh interfaces.'", &
            '                do j=1,size(face_a)', &
            '                    call disjoint_set_union(parent,node_offsets(a-1)+face_a(j),node_offsets(b-1)+face_b(j))', &
            '                end do', &
            '            end associate', &
            '        end do', &
            '        call disjoint_set_map(parent,node_map)', &
            '        allocate(shared_points(maxval(node_map),size(points,2)))', &
            '        filled=0', &
            '        do i=1,size(node_map)', &
            '            node=node_map(i)', &
            '            if(node<=filled)cycle', &
            '            shared_points(node,:)=points(i,:)', &
            '            filled=node', &
            '        end do', &
            '        call move_alloc(shared_points,points)', &
            '        do i=1,size(cells,2)', &
            '            cells(:,i)=node_map(cells(:,i))', &
            '        end do', &
            '    end subroutine stitch_mesh'])
    end subroutine write_stitching_helper

    subroutine write_vtk_helper(source)
        type(source_writer),intent(inout) :: source

        call write_lines(source,[character(132):: &
            '', &
            '    subroutine write_mesh_vtk(filename,kind,points,cells,cell_counts)', &
            '        character(*), intent(in) :: filename', &
            '        integer, intent(in) :: kind,cell_counts(:)', &
            '        real(rk), contiguous, intent(in) :: points(:,:)', &
            '        integer, contiguous, intent(in) :: cells(:,:)', &
            '        integer, parameter :: vtk_types(3)=[3,9,12]', &
            '        real(rk), allocatable :: patch_ids(:,:)', &
            '        integer :: i,offset', &
            '', &
            '        if(size(cell_counts)==1)then', &
            "            call export_vtk_legacy(filename,points,cells,vtk_types(kind),encoding='BINARY')", &
            '            return', &
            '        end if', &
            '        allocate(patch_ids(size(cells,1),1))', &
            '        offset=0', &
            '        do i=1,size(cell_counts)', &
            '            patch_ids(offset+1:offset+cell_counts(i),1)=real(i,rk)', &
            '            offset=offset+cell_counts(i)', &
            '        end do', &
            "        call export_vtk_legacy(filename,points,cells,vtk_types(kind),encoding='BINARY', &", &
            "            cell_data=patch_ids,cell_field_names=['patch_id'])", &
            '    end subroutine write_mesh_vtk'])
    end subroutine write_vtk_helper

    subroutine write_lines(source,lines)
        type(source_writer),intent(inout) :: source
        character(*),intent(in) :: lines(:)
        integer :: i

        do i=1,size(lines)
            call source%line(trim(lines(i)))
        end do
    end subroutine write_lines

end module studio_fortran_export_mesh
