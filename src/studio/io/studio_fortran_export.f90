module studio_fortran_export
    use, intrinsic :: iso_c_binding, only: c_int,c_null_char
    use qt, only: QDir_removeRecursively
    use studio_document, only: studio_document_t,document_mesh_count,document_mesh_valid
    use studio_forcad, only: rk,native_patch_data,collect_native_patch_data,multipatch_side_label
    use studio_fortran_export_text, only: source_writer,integer_list,logical_literal
    use studio_fortran_export_mesh, only: write_mesh_calls,write_mesh_routines
    use studio_multipatch, only: multipatch_draft
    use studio_multipatch_collection, only: count_multipatches
    use studio_platform_files, only: temporary_private_directory,create_private_file,current_directory,close_file, &
        replace_file,last_path_separator
    use studio_presentation_text, only: integer_text
    use studio_scene, only: project_patch_exists
    use studio_single_patch, only: single_patch_count,single_patch_available,copy_single_patch_geometry
    implicit none(type,external)
    private
    public :: export_project_fortran

    character(7),parameter :: ranks(3)=[character(7)::'curve','surface','volume']
    character(8),parameter :: patch_names(3)=[character(8)::'curves','surfaces','volumes']
    character(1),parameter :: directions(3)=['u','v','w']

contains

    subroutine export_project_fortran(document,path,error)
        type(studio_document_t),intent(in) :: document
        character(*),intent(in) :: path
        character(:),allocatable,intent(out) :: error
        type(source_writer) :: source
        character(:),allocatable :: parent,staging,temporary
        integer :: separator,io,closed
        integer(c_int) :: file
        logical :: removed

        error='Could not prepare the Fortran export destination.'
        if(len_trim(path)==0.or.index(path,c_null_char)>0)return
        separator=last_path_separator(path)
        parent='.'
        if(separator>0)parent=path(:separator)
        staging=temporary_private_directory(parent,'.forcad-source.')
        if(len(staging)==0)return
        temporary=staging//'/project.f90'
        file=create_private_file(current_directory,temporary//c_null_char)
        if(file>=0)then
            if(close_file(file)==0)then
                open(newunit=source%unit,file=temporary,form='formatted',status='old',action='write',iostat=io)
                if(io==0)then
                    call write_project(source,document,error)
                    close(source%unit,iostat=closed)
                    if(len(error)==0)then
                        error='Could not write the Fortran source.'
                        if(source%status==0.and.closed==0)then
                            error='Could not replace the destination Fortran file.'
                            if(replace_file(temporary//c_null_char,path//c_null_char)==0)error=''
                        end if
                    end if
                end if
            end if
        end if
        removed=QDir_removeRecursively(staging)
    end subroutine export_project_fortran

    subroutine write_project(source,document,error)
        type(source_writer),intent(inout) :: source
        type(studio_document_t),intent(in) :: document
        character(:),allocatable,intent(out) :: error
        logical,allocatable :: include(:,:)
        class(*),allocatable :: geometry
        type(native_patch_data) :: data
        integer :: counts(3),kind,patch,assembly,assemblies,joins,mesh,part

        counts=[(single_patch_count(kind,catalog=document%patches),kind=1,3)]
        allocate(include(maxval(counts),3),source=.false.)
        do kind=1,3
            do patch=1,counts(kind)
                include(patch,kind)=project_patch_exists(document%show_templates,kind,patch,document%patches)
            end do
        end do
        assemblies=count_multipatches(document%multipatches)
        joins=0
        do assembly=1,assemblies
            associate(group=>document%multipatches%groups(assembly))
                error='Multipatch '//integer_text(assembly)//' must be valid before exporting Fortran.'
                if(.not.group%validation%valid)return
                kind=group%definition%kind
                if(kind<1.or.kind>3.or..not.allocated(group%definition%members))return
                do part=1,size(group%definition%members)
                    patch=group%definition%members(part)
                    if(.not.single_patch_available(kind,patch,document%patches))return
                    include(patch,kind)=.true.
                end do
                if(allocated(group%definition%joins))joins=joins+size(group%definition%joins)
            end associate
        end do
        do mesh=1,document_mesh_count(document)
            if(.not.document_mesh_valid(document,mesh))cycle
            associate(item=>document%meshes(mesh))
                do part=1,size(item%parts)
                    include(item%parts(part)%variant,item%kind)=.true.
                end do
            end associate
        end do

        error=''
        call source%line('program rebuild_forcad_project')
        call source%line('    use forcad, only: rk, nurbs_curve, nurbs_surface, nurbs_volume, multipatch_connection, &')
        call source%line('        nurbs_multipatch_curve, nurbs_multipatch_surface, nurbs_multipatch_volume')
        call source%line('    use forcad_utils, only: export_vtk_legacy, boundary_index, disjoint_set_union, disjoint_set_map')
        call source%line('    implicit none')
        call source%line('')
        do kind=1,3
            call source%line('    type(nurbs_'//trim(ranks(kind))//') :: '//trim(patch_names(kind))//'('//integer_text(counts(kind))//')')
            call source%line('    type(nurbs_multipatch_'//trim(ranks(kind))//') :: '//trim(ranks(kind))// &
                '_assemblies('//integer_text(assemblies)//')')
        end do
        call source%line('')
        call source%line('    call build_patches()')
        call source%line('    call build_assemblies()')
        call write_mesh_calls(source,document)
        call source%line("    print '(a)', 'Rebuilt "//integer_text(count(include))//' patches, '//integer_text(assemblies)// &
            ' multipatches, '//integer_text(joins)//' joins and '//integer_text(document_mesh_count(document))//" meshes.'")
        call source%line('')
        call source%line('contains')
        call source%line('')
        call source%line('    subroutine build_patches()')
        do kind=1,3
            do patch=1,counts(kind)
                if(include(patch,kind))call source%line('        call '//patch_routine(kind,patch)//'('//patch_expression(kind,patch)//')')
            end do
        end do
        call source%line('    end subroutine build_patches')
        call source%line('')
        call source%line('    subroutine build_assemblies()')
        do assembly=1,assemblies
            call source%line('        call build_assembly_'//integer_text(assembly)//'()')
        end do
        call source%line('    end subroutine build_assemblies')
        call source%line('')

        do kind=1,3
            do patch=1,counts(kind)
                if(source%status/=0)return
                if(.not.include(patch,kind))cycle
                call copy_single_patch_geometry(kind,patch,geometry,catalog=document%patches)
                error='Could not read '//trim(ranks(kind))//' '//integer_text(patch)//' for Fortran export.'
                if(.not.allocated(geometry))return
                if(.not.collect_native_patch_data(geometry,data))return
                call write_patch_routine(source,kind,patch,data)
                error=''
            end do
        end do
        do assembly=1,assemblies
            call write_assembly_routine(source,assembly,document%multipatches%groups(assembly)%definition)
        end do
        call write_mesh_routines(source,document)
        call source%line('end program rebuild_forcad_project')
    end subroutine write_project

    subroutine write_patch_routine(source,kind,patch,data)
        type(source_writer),intent(inout) :: source
        integer,intent(in) :: kind,patch
        type(native_patch_data),intent(in) :: data
        character(:),allocatable :: routine,knots,weights,degree,wrap
        integer :: axis,points

        routine=patch_routine(kind,patch)
        points=size(data%points,1)
        call source%line('    subroutine '//routine//'(patch)')
        call source%line('        type(nurbs_'//trim(ranks(kind))//'), intent(out) :: patch')
        do axis=1,kind
            call source%line('        real(rk) :: knots_'//directions(axis)//'('//integer_text(size(data%axes(axis)%knots))//')')
        end do
        call source%line('        real(rk) :: controls('//integer_text(size(data%points))//')')
        weights=''
        if(allocated(data%weights))then
            if(size(data%weights)>0)then
                call source%line('        real(rk) :: weights('//integer_text(size(data%weights))//')')
                weights=', Wc=weights'
            end if
        end if
        call source%line('')
        do axis=1,kind
            call source%reals('knots_'//directions(axis),data%axes(axis)%knots)
        end do
        call source%reals('controls',reshape(data%points,[size(data%points)]))
        if(len(weights)>0)call source%reals('weights',data%weights)
        degree=integer_text(data%degree(1))
        wrap=logical_literal(data%wrap(1))
        knots='knot=knots_u'
        if(kind>1)then
            knots='knot1=knots_u, knot2=knots_v'
            if(kind==3)knots=knots//', knot3=knots_w'
            degree='['//integer_list(data%degree(:kind))//']'
            do axis=2,kind
                wrap=wrap//', '//logical_literal(data%wrap(axis))
            end do
            wrap='['//wrap//']'
        end if
        call source%line('        call patch%set('//knots//', &')
        call source%line('            Xc=reshape(controls, ['//integer_text(points)//', '// &
            integer_text(size(data%points,2))//'])'//weights//', &')
        call source%line('            degree='//degree//', wrap_parameters='//wrap//')')
        call source%line("        if (.not. patch%err%ok) error stop 'ForCAD rejected "//trim(ranks(kind))//' '//integer_text(patch)//".'")
        call source%line('    end subroutine '//routine)
        call source%line('')
    end subroutine write_patch_routine

    subroutine write_assembly_routine(source,id,definition)
        type(source_writer),intent(inout) :: source
        integer,intent(in) :: id
        type(multipatch_draft),intent(in) :: definition
        character(:),allocatable :: group,routine,arguments,checks
        integer :: kind,i,a,b,join_count

        kind=definition%kind
        join_count=0
        if(allocated(definition%joins))join_count=size(definition%joins)
        group=trim(ranks(kind))//'_assemblies('//integer_text(id)//')'
        routine='build_assembly_'//integer_text(id)
        call source%line('    subroutine '//routine//'()')
        if(allocated(definition%joins))then
            do i=1,size(definition%joins)
                if(.not.allocated(definition%joins(i)%reparameterization))cycle
                call source%line('        real(rk), allocatable :: jet(:)')
                exit
            end do
        end if
        call source%line('        associate(group => '//group//')')
        do i=1,size(definition%members)
            call source%line('        call group%add_patch('//patch_expression(kind,definition%members(i))//')')
        end do
        checks='        if (.not. group%err%ok) error stop '
        call source%line(checks//"'Could not add patches to multipatch "//integer_text(id)//".'")
        call source%line('        if (group%get_npatch() /= '//integer_text(size(definition%members))//') &')
        call source%line("            error stop 'Unexpected patch count in multipatch "//integer_text(id)//".'")
        if(allocated(definition%joins))then
            do i=1,size(definition%joins)
                associate(join=>definition%joins(i))
                    a=findloc(definition%members,join%patch_a,dim=1)
                    b=findloc(definition%members,join%patch_b,dim=1)
                    arguments=integer_text(a)//", '"//multipatch_side_label(kind,join%side_a)//"', "// &
                        integer_text(b)//", '"//multipatch_side_label(kind,join%side_b)//"', "//integer_text(join%continuity)
                    if(allocated(join%reparameterization))then
                        call source%line('        allocate(jet('//integer_text(size(join%reparameterization))//'))')
                        call source%reals('jet',real(join%reparameterization,rk))
                    end if
                    call source%line('        call group%connect('//arguments//', &')
                    arguments='            reverse='//logical_literal(join%reverse)//', geometric='//logical_literal(join%geometric)
                    if(kind==3)then
                        call source%line(arguments//', &')
                        arguments='            swap='//logical_literal(join%swap)//', flip=['// &
                            logical_literal(join%flip(1))//', '//logical_literal(join%flip(2))//']'
                    end if
                    if(allocated(join%reparameterization))arguments=arguments//', reparameterization=jet'
                    call source%line(arguments//')')
                    if(allocated(join%reparameterization))call source%line('        deallocate(jet)')
                end associate
            end do
        end if
        call source%line(checks//"'Could not reconnect multipatch "//integer_text(id)//".'")
        call source%line('        if (group%get_nconnection() /= '//integer_text(join_count)//') &')
        call source%line("            error stop 'Unexpected join count in multipatch "//integer_text(id)//".'")
        call source%line('        end associate')
        call source%line('    end subroutine '//routine)
        call source%line('')
    end subroutine write_assembly_routine

    function patch_routine(kind,patch) result(name)
        integer,intent(in) :: kind,patch
        character(:),allocatable :: name

        name='build_'//trim(ranks(kind))//'_'//integer_text(patch)
    end function patch_routine

    function patch_expression(kind,patch) result(expression)
        integer,intent(in) :: kind,patch
        character(:),allocatable :: expression

        expression=trim(patch_names(kind))//'('//integer_text(patch)//')'
    end function patch_expression

end module studio_fortran_export
