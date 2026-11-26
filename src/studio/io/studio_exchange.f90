! Connect file exchange to ForCAD and its IGES dependency. The dependency's
! reader runs in a child process because malformed input can terminate it.
module studio_exchange
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_c_binding, only: c_double,c_null_char
    use, intrinsic :: iso_fortran_env, only: int64
    use studio_single_patch, only: studio_patch_catalog,publish_native_patches,copy_single_patch_geometry
    use studio_mesh, only: studio_mesh_t
    use studio_forcad, only: rk,studio_native_patch,native_patch_kind,affine_transform_t,transform_patch, &
        decode_studio_iges=>write_iges_worker_output, read_iges_worker_output, &
        write_native_iges, write_iges_units, write_native_vtk, write_native_multipatch_vtk, write_mesh_vtk, native_multipatch_vtk_path_limit
    use qt, only: QCoreApplication_applicationFilePath, QProcess, QProcess_NormalExit, &
        QProcess_NotRunning, QProcess_nullDevice, QTemporaryDir
    use studio_project_file, only: read_patch,write_patch
    use studio_platform_files, only: replace_file,last_path_separator
    use studio_units, only: unitless,unit_m,valid_unit_id,display_length
    use studio_storage_paths, only: studio_work_directory
    use studio_presentation_text, only: entity_name,integer_text
    implicit none(type, external)
    private
    public :: import_studio_iges,export_studio_patch,export_studio_multipatch,dispatch_studio_iges_worker
    public :: studio_exchange_suffix,export_studio_mesh

    character(*), parameter :: decode_switch='--decode-iges',encode_switch='--encode-iges'
    character(*), parameter :: patch_payload='patch.bin',raw_iges_file='geometry.igs'
    character(*), parameter :: staged_iges_source='source.igs',decoded_patches='geometry'
    integer, parameter :: name_length=32

contains
    function studio_exchange_suffix(path) result(suffix)
        character(*),intent(in) :: path
        character(:),allocatable :: suffix
        integer :: i,dot,slash,code
        suffix=''
        dot=index(path,'.',back=.true.)
        slash=last_path_separator(path)
        if(dot<=slash+1)return
        suffix=path(dot:)
        do i=1,len(suffix)
            code=iachar(suffix(i:i))
            if(code>=iachar('A').and.code<=iachar('Z'))suffix(i:i)=achar(code+32)
        end do
    end function studio_exchange_suffix

    subroutine import_studio_iges(path,kind,variant,count,skipped,error,catalog,model_unit,source_unit,needs_source_unit,project_directory)
        character(*),intent(in) :: path
        character(*),intent(in),optional :: project_directory
        integer,intent(out) :: kind,variant,count,skipped
        character(:),allocatable,intent(out) :: error
        type(studio_patch_catalog),intent(inout),optional :: catalog
        integer,intent(in),optional :: model_unit,source_unit
        logical,intent(out),optional :: needs_source_unit
        type(studio_native_patch),allocatable :: geometries(:)
        character(name_length),allocatable :: names(:)
        character(:),allocatable :: directory,source,decoded
        type(QTemporaryDir) :: workspace
        integer :: destination_unit,i
        real(rk) :: file_unit_metres,file_scale

        kind=0; variant=0; count=0; skipped=0
        if(present(needs_source_unit))needs_source_unit=.false.
        destination_unit=unitless
        if(present(model_unit))destination_unit=model_unit
        directory=studio_work_directory(project_directory)
        if(len(directory)==0)then
            error='Could not prepare a private workspace for IGES import.'
            return
        end if
        workspace=QTemporaryDir(directory//'/iges-import.XXXXXX')
        if(.not.workspace%isValid())then
            call workspace%delete()
            error='Could not prepare a private workspace for IGES import.'
            return
        end if
        directory=workspace%path()
        source=workspace%filePath(staged_iges_source)
        decoded=workspace%filePath(decoded_patches)
        if(.not.copy_source(path,source))then
            error='Could not read the IGES source file.'
        else if(.not.run_iges_worker([character(len(decode_switch))::decode_switch,staged_iges_source,decoded_patches],directory))then
            error='The IGES reader could not decode this file.'
        else
            call read_iges_worker_output(decoded,geometries,skipped,error,file_unit_metres,file_scale)
            if(len(error)==0)call convert_iges_units(geometries,destination_unit,file_unit_metres,file_scale,source_unit,needs_source_unit,error)
            if(len(error)==0)then
                allocate(names(size(geometries)))
                do i=1,size(geometries)
                    names(i)='Imported '//entity_name(native_patch_kind(geometries(i)%geometry))//' '//integer_text(i)
                end do
                call publish_native_patches(geometries,names,kind,variant,error,catalog=catalog)
                if(len(error)==0)count=size(geometries)
            end if
        end if
        call workspace%delete()
    end subroutine import_studio_iges

    subroutine convert_iges_units(geometries,model_unit,file_unit_metres,file_scale,source_unit,needs_source_unit,error)
        type(studio_native_patch), intent(inout) :: geometries(:)
        integer, intent(in) :: model_unit
        real(rk), intent(in) :: file_unit_metres,file_scale
        integer, intent(in), optional :: source_unit
        logical, intent(out), optional :: needs_source_unit
        character(:), allocatable, intent(out) :: error
        type(affine_transform_t) :: transform
        real(rk) :: metres,factor
        real(c_double) :: conversion
        integer :: i

        error=''
        if(present(needs_source_unit))needs_source_unit=.false.
        if(model_unit==unitless)return
        error='The project model unit is unavailable.'
        if(.not.valid_unit_id(model_unit))return
        error='The IGES file contains an invalid model-space scale.'
        if(.not.ieee_is_finite(file_scale).or.file_scale<=0._rk)return

        metres=file_unit_metres
        if(metres==0._rk)then
            error='The IGES file does not identify its length units. Choose a source unit or keep its coordinates.'
            if(.not.present(source_unit))then
                if(present(needs_source_unit))needs_source_unit=.true.
                return
            end if
            if(.not.valid_unit_id(source_unit))return
            if(source_unit==unitless)then
                error=''
                return
            end if
            metres=real(display_length(1._c_double,source_unit,unit_m),rk)
        end if

        error='The IGES length conversion is outside the supported coordinate range.'
        factor=metres*real(display_length(1._c_double,unit_m,model_unit),rk)
        if(file_scale<1._rk)then
            if(factor>huge(factor)*file_scale)return
        end if
        factor=factor/file_scale
        if(factor<=0._rk.or.exponent(factor)>maxexponent(conversion))return
        if(exponent(factor)==maxexponent(conversion))then
            if(fraction(factor)>real(fraction(huge(conversion)),rk))return
        end if
        conversion=real(factor,c_double)
        if(.not.ieee_is_finite(conversion).or.conversion<=0._c_double)return
        if(conversion/=1._c_double)then
            transform%linear=conversion*transform%linear
            do i=1,size(geometries)
                if(.not.transform_patch(geometries(i)%geometry,transform))return
            end do
        end if
        error=''
    end subroutine convert_iges_units

    logical function copy_source(path,destination) result(copied)
        character(*),intent(in) :: path,destination
        character(65536) :: block
        integer :: input,output,status,n
        integer(int64) :: total,remaining
        copied=.false.
        inquire(file=path,size=total,iostat=status)
        if(status/=0.or.total<=0)return
        open(newunit=input,file=path,status='old',access='stream',form='unformatted',action='read',iostat=status)
        if(status/=0)return
        open(newunit=output,file=destination,status='new',access='stream',form='unformatted',action='write',iostat=status)
        if(status/=0)then
            close(input)
            return
        end if
        remaining=total
        do while(remaining>0)
            n=int(min(remaining,int(len(block),int64)))
            read(input,iostat=status)block(:n)
            if(status/=0)exit
            write(output,iostat=status)block(:n)
            if(status/=0)exit
            remaining=remaining-n
        end do
        close(input)
        close(output,iostat=status)
        copied=remaining==0.and.status==0
    end function copy_source

    logical function run_iges_worker(arguments,working_directory) result(succeeded)
        character(*),intent(in) :: arguments(:)
        character(*),intent(in) :: working_directory
        character(:),allocatable :: executable
        type(QProcess) :: worker
        logical :: finished

        succeeded=.false.
        executable=QCoreApplication_applicationFilePath()
        if(len(executable)==0)return

        worker=QProcess()
        call worker%setWorkingDirectory(working_directory)
        call worker%setStandardOutputFile(QProcess_nullDevice())
        call worker%setStandardErrorFile(QProcess_nullDevice())
        call worker%start(executable,arguments)
        if(worker%waitForStarted())then
            finished=worker%waitForFinished(-1)
            if(finished.or.worker%state()==QProcess_NotRunning)then
                succeeded=worker%exitStatus()==QProcess_NormalExit.and.worker%exitCode()==0
            end if
        end if
        call worker%delete()
    end function run_iges_worker

    subroutine dispatch_studio_iges_worker()
        character(len=len(decode_switch)) :: mode
        character(:),allocatable :: source,destination
        integer :: length,status

        call get_command_argument(1,mode,status=status)
        if(status/=0)return
        select case(mode)
        case(decode_switch)
            if(command_argument_count()/=3)stop 1,quiet=.true.
            call get_command_argument(2,length=length,status=status)
            if(status/=0)stop 1,quiet=.true.
            allocate(character(length)::source)
            call get_command_argument(2,source,status=status)
            if(status/=0)stop 1,quiet=.true.
            call get_command_argument(3,length=length,status=status)
            if(status/=0)stop 1,quiet=.true.
            allocate(character(length)::destination)
            call get_command_argument(3,destination,status=status)
            if(status/=0)stop 1,quiet=.true.
            status=decode_studio_iges(source,destination)
        case(encode_switch)
            if(command_argument_count()/=1)stop 1,quiet=.true.
            status=merge(0,1,encode_staged_iges())
        case default
            return
        end select
        stop status,quiet=.true.
    end subroutine dispatch_studio_iges_worker

    logical function encode_staged_iges() result(written)
        class(*),allocatable :: geometry
        integer :: input,kind,status

        written=.false.
        open(newunit=input,file=patch_payload,status='old',access='stream',form='unformatted',action='read',iostat=status)
        if(status/=0)return
        read(input,iostat=status)kind
        if(status==0)then
            if(kind==1.or.kind==2)call read_patch(input,kind,geometry,status)
        end if
        close(input)
        if(status/=0.or..not.allocated(geometry))return
        written=write_native_iges(geometry,raw_iges_file)
    end function encode_staged_iges

    logical function export_iges_in_worker(kind,geometry,directory) result(written)
        integer,intent(in) :: kind
        class(*),intent(in) :: geometry
        character(*),intent(in) :: directory
        integer :: output,status,closed

        written=.false.
        open(newunit=output,file=directory//'/'//patch_payload,status='new',access='stream',form='unformatted',action='write',iostat=status)
        if(status/=0)return
        write(output,iostat=status)kind
        if(status==0)call write_patch(output,kind,geometry,status)
        close(output,iostat=closed)
        if(status/=0.or.closed/=0)return
        written=run_iges_worker([encode_switch],working_directory=directory)
    end function export_iges_in_worker

    subroutine export_studio_patch(kind,variant,path,error,model_unit,catalog)
        integer,intent(in) :: kind,variant
        character(*),intent(in) :: path
        character(:),allocatable,intent(out) :: error
        integer,intent(in),optional :: model_unit
        type(studio_patch_catalog),intent(in),optional :: catalog
        class(*),allocatable :: geometry
        character(:),allocatable :: directory,output,with_units,suffix
        integer :: export_unit
        type(QTemporaryDir) :: workspace
        logical :: valid

        suffix=studio_exchange_suffix(path)
        if(suffix/='.igs'.and.suffix/='.iges'.and.suffix/='.vtk')then
            error='Choose a destination ending in .igs, .iges, or .vtk.'
            return
        end if
        if(kind==3.and.suffix/='.vtk')then
            error='Use VTK to export a volume patch.'
            return
        end if
        call copy_single_patch_geometry(kind,variant,geometry,catalog=catalog)
        if(.not.allocated(geometry))then
            error='Select a patch to export.'
            return
        end if
        call prepare_export(path,workspace,error)
        if(len(error)>0)return
        directory=workspace%path()
        output=directory//'/patch'//suffix
        if(suffix=='.vtk')then
            valid=write_native_vtk(geometry,output,generate=.true.)
        else
            output=directory//'/'//raw_iges_file
            valid=export_iges_in_worker(kind,geometry,directory)
        end if
        export_unit=unitless
        if(present(model_unit))export_unit=model_unit
        if(valid.and.suffix/='.vtk')then
            with_units=directory//'/patch-with-units'//suffix
            valid=write_iges_units(output,with_units,export_unit)
            call finish_export(valid,with_units,path,workspace,error)
            if(.not.valid)error='Could not write the selected model units into the IGES file.'
            return
        end if
        call finish_export(valid,output,path,workspace,error)
    end subroutine export_studio_patch

    subroutine export_studio_multipatch(native,path,error)
        class(*),intent(in) :: native
        character(*),intent(in) :: path
        character(:),allocatable,intent(out) :: error
        class(*),allocatable :: geometry
        character(:),allocatable :: directory,prefix
        type(QTemporaryDir) :: workspace
        logical :: valid

        if(studio_exchange_suffix(path)/='.vtk')then
            error='Choose a .vtk destination for the multipatch.'
            return
        end if
        call prepare_export(path,workspace,error)
        if(len(error)>0)return
        directory=workspace%path()
        prefix=directory//'/group'
        valid=.false.
        if(len(prefix)+len('_Xg.vtk')<=native_multipatch_vtk_path_limit)then
            allocate(geometry,source=native)
            valid=write_native_multipatch_vtk(geometry,prefix,generate=.true.)
        end if
        call finish_export(valid,prefix//'_Xg.vtk',path,workspace,error)
        if(len(prefix)+len('_Xg.vtk')>native_multipatch_vtk_path_limit)error='The destination directory is too long for the native multipatch exporter.'
    end subroutine export_studio_multipatch

    subroutine export_studio_mesh(mesh,path,error)
        type(studio_mesh_t),intent(in) :: mesh
        character(*),intent(in) :: path
        character(:),allocatable,intent(out) :: error
        character(:),allocatable :: output
        type(QTemporaryDir) :: workspace
        logical :: valid

        error='Generate a mesh before exporting it.'
        if(.not.allocated(mesh%parts).or.mesh%elements==0)return
        if(studio_exchange_suffix(path)/='.vtk')then
            error='Choose a .vtk destination for the mesh.'
            return
        end if
        call prepare_export(path,workspace,error)
        if(len(error)>0)return
        output=workspace%filePath('mesh.vtk')
        valid=write_mesh_vtk(output,mesh%kind,mesh%parts,mesh%connected_points,mesh%connected_cells)
        call finish_export(valid,output,path,workspace,error)
    end subroutine export_studio_mesh

    subroutine prepare_export(path,workspace,error)
        character(*),intent(in) :: path
        type(QTemporaryDir),intent(out) :: workspace
        character(:),allocatable,intent(out) :: error
        integer :: slash

        slash=last_path_separator(path)
        workspace=QTemporaryDir(path(:slash)//'.forcad-exchange.XXXXXX')
        error=''
        if(.not.workspace%isValid())then
            call workspace%delete()
            error='Could not write to the selected destination directory.'
        end if
    end subroutine prepare_export

    subroutine finish_export(valid,output,path,workspace,error)
        logical,intent(in) :: valid
        character(*),intent(in) :: output,path
        type(QTemporaryDir),intent(inout) :: workspace
        character(:),allocatable,intent(out) :: error
        integer(int64) :: bytes
        integer :: io

        error='ForCAD could not write the selected geometry.'
        bytes=-1
        inquire(file=output,size=bytes,iostat=io)
        if(valid.and.io==0.and.bytes>0)then
            error='Could not save the file at the selected destination.'
            if(replace_file(output//c_null_char,path//c_null_char)==0)error=''
        end if
        call workspace%delete()
    end subroutine finish_export
end module studio_exchange
