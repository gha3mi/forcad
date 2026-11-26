module studio_forcad_exchange
    use, intrinsic :: iso_fortran_env, only: int64,file_storage_size
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use forcad, only: rk,nurbs_curve,nurbs_surface,nurbs_volume, &
        nurbs_multipatch_curve,nurbs_multipatch_surface,nurbs_multipatch_volume
    use forIGES, only: IGES_t,entity124_t,entity126_t,entity128_t,makeGsection,makeTsection, &
        MILLIMETERS,CENTIMETERS,METERS,INCHES,FEET,PARAM15,MILES,KILOMETERS,MILS,MICRONS,MICROINCHES, &
        field_separator,record_terminator,FORIGES_VERSION
    use forcad_utils, only: export_vtk_legacy
    use studio_forcad_transforms, only: studio_native_patch
    use studio_forcad_mesh, only: forcad_mesh_samples
    use studio_units, only: unitless,unit_mm,unit_cm,unit_m,unit_in,unit_ft
    implicit none(type, external)
    private

    public :: write_iges_worker_output, read_iges_worker_output
    public :: write_native_iges, write_iges_units, write_native_vtk, write_native_multipatch_vtk, write_mesh_vtk
    integer, parameter, public :: native_multipatch_vtk_path_limit=512

    character(*), parameter :: worker_signature='FORCAD-IG27'
    character(*), parameter :: foriges_native_id='forIGES '//FORIGES_VERSION
    integer, parameter :: curve_resolution=81, patch_resolution=17

contains

    integer function write_iges_worker_output(path,output) result(status)
        character(*),intent(in) :: path,output
        type(IGES_t) :: model
        integer :: unit,io,i,count,skipped,entity_type
        real(rk) :: file_unit_metres,file_scale
        logical :: written

        status=1
        model%filename=path
        call model%makeFileImage()
        call model%splitFileImage()
        call read_iges_geometry_units(model,file_unit_metres,file_scale)
        call model%getD()
        call model%getP()
        call model%setPtmatForm()
        open(newunit=unit,file=output,access='stream',form='unformatted',status='new',action='write',iostat=io)
        if(io/=0)return
        write(unit,iostat=io)worker_signature,0,0,file_unit_metres,file_scale
        count=0; skipped=0
        do i=1,model%nentities
            if(io/=0)exit
            entity_type=model%D%DE(i)%entity_type
            if(entity_type==124)cycle
            if(any(model%D%DE(i)%status_numbers(:3)/=0))then
                skipped=skipped+1
                cycle
            end if
            written=.false.
            select type(entity=>model%P%PE(i)%entity)
            type is(entity126_t)
                call decode_curve(model,entity,unit,io,written)
            type is(entity128_t)
                call decode_surface(model,entity,unit,io,written)
            end select
            if(written)then
                count=count+1
            else
                skipped=skipped+1
            end if
        end do
        if(io==0)write(unit,pos=1,iostat=io)worker_signature,count,skipped,file_unit_metres,file_scale
        if(io==0)status=0
        close(unit,iostat=io)
        if(io/=0)status=1
    end function write_iges_worker_output

    subroutine read_iges_geometry_units(model,metres,scale)
        type(IGES_t), intent(in) :: model
        real(rk), intent(out) :: metres,scale
        character(:), allocatable :: globals,value,unit_name
        character :: separator,terminator
        integer :: cursor,record,field,io,unit_flag
        logical :: good

        metres=0._rk
        scale=1._rk
        separator=','
        terminator=';'
        if(.not.allocated(model%Grecords))return
        if(size(model%Grecords)>huge(record)/72)return
        allocate(character(72*size(model%Grecords))::globals)
        do record=1,size(model%Grecords)
            globals(72*(record-1)+1:72*record)=model%Grecords(record)(:72)
        end do
        cursor=1
        call read_global_field(globals,cursor,separator,terminator,value,good,declares_separator=.true.)
        if(.not.good)return
        field_separator=separator
        call read_global_field(globals,cursor,separator,terminator,value,good)
        if(.not.good.or.len(value)>1)return
        if(len(value)==1)terminator=value(1:1)
        if(separator==terminator)return
        record_terminator=terminator
        do field=3,13
            call read_global_field(globals,cursor,separator,terminator,value,good)
            if(.not.good)return
        end do
        if(len(value)>0)then
            scale=0._rk
            if(verify(value,'0123456789+-.eEdD')==0)then
                read(value,*,iostat=io)scale
                if(io/=0.or..not.ieee_is_finite(scale))scale=0._rk
            end if
        end if
        call read_global_field(globals,cursor,separator,terminator,value,good)
        if(.not.good.or.len(value)==0)return
        if(verify(value,'0123456789')/=0)return
        read(value,*,iostat=io)unit_flag
        if(io/=0)return
        unit_name=''
        if(unit_flag==PARAM15)then
            call read_global_field(globals,cursor,separator,terminator,value,good)
            if(good)unit_name=value
        end if
        metres=iges_metres_per_unit(unit_flag,unit_name)
    end subroutine read_iges_geometry_units

    subroutine read_global_field(text,cursor,separator,terminator,value,valid,declares_separator)
        character(*), intent(in) :: text
        integer, intent(inout) :: cursor
        character, intent(inout) :: separator
        character, intent(in) :: terminator
        character(:), allocatable, intent(out) :: value
        logical, intent(out) :: valid
        logical, intent(in), optional :: declares_separator
        integer :: first,last,digit,count,io

        valid=.false.
        value=''
        if(cursor<1.or.cursor>len(text))return
        first=verify(text(cursor:),' '//achar(9))
        if(first==0)return
        cursor=cursor+first-1
        first=cursor
        digit=verify(text(first:),'0123456789')
        if(digit>1)then
            last=first+digit-1
            if(index('Hh',text(last:last))>0)then
                read(text(first:last-1),*,iostat=io)count
                if(io/=0)return
                if(count<0.or.count>len(text)-last)return
                value=text(last+1:last+count)
                cursor=last+count+1
                if(present(declares_separator))then
                    if(declares_separator.and.len(value)==1)separator=value(1:1)
                end if
                if(cursor>len(text))return
                first=verify(text(cursor:),' '//achar(9))
                if(first==0)return
                cursor=cursor+first-1
                if(text(cursor:cursor)/=separator.and.text(cursor:cursor)/=terminator)then
                    if(value/=trim(foriges_native_id)//separator)return
                    value=trim(foriges_native_id)
                    cursor=last+count
                end if
                valid=.true.
            end if
        end if
        if(.not.valid)then
            last=scan(text(cursor:),separator//terminator)
            if(last==0)return
            value=trim(text(cursor:cursor+last-2))
            cursor=cursor+last-1
        end if
        valid=.true.
        if(text(cursor:cursor)==separator)then
            cursor=cursor+1
        else
            cursor=0
        end if
    end subroutine read_global_field

    real(rk) function iges_metres_per_unit(flag,name) result(metres)
        integer, intent(in) :: flag
        character(*), intent(in) :: name
        character(len(name)) :: upper
        integer :: unit,i,letter

        unit=flag
        if(unit==PARAM15)then
            upper=adjustl(name)
            do i=1,len_trim(upper)
                letter=iachar(upper(i:i))
                if(letter>=iachar('a').and.letter<=iachar('z'))upper(i:i)=achar(letter-32)
            end do
            select case(trim(upper))
            case('IN','INCH','INCHES')
                unit=INCHES
            case('MM')
                unit=MILLIMETERS
            case('FT','FOOT','FEET')
                unit=FEET
            case('MI')
                unit=MILES
            case('M')
                unit=METERS
            case('KM')
                unit=KILOMETERS
            case('MIL')
                unit=MILS
            case('UM')
                unit=MICRONS
            case('CM')
                unit=CENTIMETERS
            case('UIN')
                unit=MICROINCHES
            case default
                continue
            end select
        end if

        metres=0._rk
        select case(unit)
        case(INCHES)
            metres=.0254_rk
        case(MILLIMETERS)
            metres=.001_rk
        case(FEET)
            metres=.3048_rk
        case(MILES)
            metres=1609.344_rk
        case(METERS)
            metres=1._rk
        case(KILOMETERS)
            metres=1000._rk
        case(MILS)
            metres=.0000254_rk
        case(MICRONS)
            metres=.000001_rk
        case(CENTIMETERS)
            metres=.01_rk
        case(MICROINCHES)
            metres=.0000000254_rk
        case default
            continue
        end select
    end function iges_metres_per_unit

    subroutine decode_curve(model,entity,unit,io,written)
        type(IGES_t),intent(in) :: model
        type(entity126_t),intent(in) :: entity
        integer,intent(in) :: unit
        integer,intent(inout) :: io
        logical,intent(out) :: written
        type(nurbs_curve) :: curve
        real(rk),allocatable :: controls(:,:),weights(:),knots(:)
        real(rk) :: limits(2)
        integer :: n,degree

        written=.false.
        if(entity%form/=0)then
            if(entity%form<1.or.entity%form>5)return
        end if
        n=size(entity%X)
        controls=reshape(real([entity%X,entity%Y,entity%Z],rk),[n,3])
        call curve%set(knot=real(entity%T,rk),Xc=controls,Wc=real(entity%W,rk), &
            degree=entity%M,wrap_parameters=entity%PROP4==1)
        if(curve%get_nc()/=n)return
        limits=real(entity%V,rk)
        if(.not.match_interval(real(entity%T,rk),entity%M,n,limits))then
            curve=curve%restrict(limits(1),limits(2))
            if(curve%get_nc()==0)return
        end if
        if(.not.apply_iges_transform(model,entity%tmat,curve))return
        controls=curve%get_Xc(); weights=curve%get_Wc(); knots=curve%get_knot()
        n=curve%get_nc(); degree=curve%get_degree()
        if(size(controls,1)/=n.or.n<=0)return
        if(size(weights)==0)weights=spread(1._rk,1,n)
        write(unit,iostat=io)1,[degree,0],[n,1],curve%get_parameter_wrapping(),.false.,knots,controls,weights
        written=io==0
    end subroutine decode_curve

    subroutine decode_surface(model,entity,unit,io,written)
        type(IGES_t),intent(in) :: model
        type(entity128_t),intent(in) :: entity
        integer,intent(in) :: unit
        integer,intent(inout) :: io
        logical,intent(out) :: written
        type(nurbs_surface) :: surface
        real(rk),allocatable :: controls(:,:),weights(:),knots_u(:),knots_v(:)
        real(rk) :: limits(2)
        integer :: n,nc(2),degree(2)

        written=.false.
        ! Analytic form labels still carry the same explicit NURBS data.
        if(entity%form<0.or.entity%form>9)return
        n=size(entity%X)
        nc=shape(entity%X)
        controls=reshape(real([reshape(entity%X,[n]),reshape(entity%Y,[n]),reshape(entity%Z,[n])],rk),[n,3])
        call surface%set(knot1=real(entity%S,rk),knot2=real(entity%T,rk),Xc=controls, &
            Wc=real(reshape(entity%W,[n]),rk),degree=[entity%M1,entity%M2], &
            wrap_parameters=[entity%PROP4==1,entity%PROP5==1])
        if(any(surface%get_nc()/=nc))return
        limits=real(entity%U,rk)
        if(.not.match_interval(real(entity%S,rk),entity%M1,nc(1),limits))surface=surface%restrict(1,limits(1),limits(2))
        if(any(surface%get_nc()<=0))return
        limits=real(entity%V,rk)
        if(.not.match_interval(real(entity%T,rk),entity%M2,nc(2),limits))surface=surface%restrict(2,limits(1),limits(2))
        if(any(surface%get_nc()<=0))return
        if(.not.apply_iges_transform(model,entity%tmat,surface))return
        controls=surface%get_Xc(); weights=surface%get_Wc()
        knots_u=surface%get_knot(1); knots_v=surface%get_knot(2)
        nc=surface%get_nc(); degree=surface%get_degree()
        n=product(nc)
        if(size(controls,1)/=n.or.n<=0)return
        if(size(weights)==0)weights=spread(1._rk,1,n)
        write(unit,iostat=io)2,degree,nc,surface%get_parameter_wrapping(1),surface%get_parameter_wrapping(2), &
            knots_u,knots_v,controls,weights
        written=io==0
    end subroutine decode_surface

    logical function match_interval(knots,degree,count,limits) result(matches)
        real(rk),intent(in) :: knots(:)
        integer,intent(in) :: degree,count
        real(rk),intent(inout) :: limits(2)
        real(rk) :: lower,upper

        lower=knots(degree+1); upper=knots(count+1)
        ! Clamp file rounding at native parameter endpoints before restricting.
        if(abs(limits(1)-lower)<=32._rk*epsilon(lower)*max(1._rk,abs(lower)))limits(1)=lower
        if(abs(limits(2)-upper)<=32._rk*epsilon(upper)*max(1._rk,abs(upper)))limits(2)=upper
        matches=limits(1)==lower.and.limits(2)==upper
    end function match_interval

    logical function apply_iges_transform(model,first,geometry) result(applied)
        type(IGES_t),intent(in) :: model
        integer,intent(in) :: first
        class(*),intent(inout) :: geometry
        integer :: i,index

        index=first
        applied=.false.
        do i=1,model%nentities+1
            if(index==0)then
                applied=.true.
                return
            end if
            if(index<1.or.index>model%nentities)return
            select type(matrix=>model%P%PE(index)%entity)
            type is(entity124_t)
                if(.not.all(ieee_is_finite(matrix%R)).or..not.all(ieee_is_finite(matrix%T)))return
                select type(geometry)
                type is(nurbs_curve)
                    call geometry%transform(linear=real(matrix%R,rk),translation=real(matrix%T,rk))
                    if(.not.geometry%err%ok)return
                type is(nurbs_surface)
                    call geometry%transform(linear=real(matrix%R,rk),translation=real(matrix%T,rk))
                    if(.not.geometry%err%ok)return
                end select
                index=matrix%tmat
            class default
                return
            end select
        end do
    end function apply_iges_transform

    subroutine read_iges_worker_output(path,geometries,skipped,error,file_unit_metres,file_scale)
        character(*),intent(in) :: path
        type(studio_native_patch),allocatable,intent(out) :: geometries(:)
        integer,intent(out) :: skipped
        character(:),allocatable,intent(out) :: error
        real(rk),intent(out) :: file_unit_metres,file_scale
        type(nurbs_curve) :: curve
        type(nurbs_surface) :: surface
        real(rk),allocatable :: ku(:),kv(:),controls(:,:),weights(:)
        character(len(worker_signature)) :: signature
        integer :: unit,io,count,i,rank,degree(2),nc(2),n,nku,nkv
        integer(int64) :: file_size,nvalues,wide_count
        logical :: wrap(2)

        skipped=0; count=0; file_size=-1
        file_unit_metres=0._rk; file_scale=1._rk
        error='The IGES reader returned incomplete geometry.'
        inquire(file=path,size=file_size,iostat=io)
        if(io/=0.or.file_size<=0)return
        open(newunit=unit,file=path,status='old',access='stream',form='unformatted',action='read',iostat=io)
        if(io/=0)return
        read(unit,iostat=io)signature,count,skipped,file_unit_metres,file_scale
        if(io==0)then
            if(signature/=worker_signature.or.count<0.or.count>file_size.or.skipped<0)io=1
        end if
        if(io==0)allocate(geometries(count),stat=io)
        do i=1,merge(count,0,io==0)
            read(unit,iostat=io)rank,degree,nc,wrap
            if(io/=0)exit
            if(rank<1.or.rank>2.or.any(nc<1).or.any(degree<0))then
                io=1
                exit
            end if
            if(rank==1.and.nc(2)/=1)then
                io=1
                exit
            end if
            wide_count=int(nc(1),int64)*int(nc(2),int64)
            nvalues=sum(int(nc(:rank),int64))+sum(int(degree(:rank),int64))+rank
            if(wide_count>huge(n).or.nvalues>huge(n))then
                io=1
                exit
            end if
            if((4_int64*wide_count+nvalues)*(storage_size(0._rk)/file_storage_size)>file_size)then
                io=1
                exit
            end if
            n=int(wide_count)
            nku=nc(1)+degree(1)+1
            nkv=0
            if(rank==2)nkv=nc(2)+degree(2)+1
            allocate(ku(nku),stat=io)
            if(io==0)allocate(kv(nkv),stat=io)
            if(io==0)allocate(controls(n,3),stat=io)
            if(io==0)allocate(weights(n),stat=io)
            if(io/=0)exit
            read(unit,iostat=io)ku,kv,controls,weights
            if(io/=0)exit
            if(rank==1)then
                curve=nurbs_curve()
                call curve%set(knot=ku,Xc=controls,Wc=weights,degree=degree(1),wrap_parameters=wrap(1))
                allocate(geometries(i)%geometry,source=curve)
            else
                surface=nurbs_surface()
                call surface%set(knot1=ku,knot2=kv,Xc=controls,Wc=weights,degree=degree,wrap_parameters=wrap)
                allocate(geometries(i)%geometry,source=surface)
            end if
            deallocate(ku,kv,controls,weights)
        end do
        close(unit)
        if(io/=0)return
        error=''
        if(count==0)error='No standalone IGES NURBS curves or surfaces were found.'
    end subroutine read_iges_worker_output

    logical function write_native_iges(geometry,path) result(written)
        class(*), intent(inout) :: geometry
        character(*), intent(in) :: path

        written=.false.
        select type(geometry)
        type is(nurbs_curve)
            call geometry%export_iges(path)
            written=geometry%err%ok
        type is(nurbs_surface)
            call geometry%export_iges(path)
            written=geometry%err%ok
        end select
    end function write_native_iges

    logical function write_iges_units(source,destination,model_unit) result(written)
        character(*),intent(in) :: source,destination
        integer,intent(in) :: model_unit
        type(IGES_t) :: model
        character(80),allocatable :: global_records(:)
        character(80) :: termination
        character(11) :: name
        integer(int64) :: bytes
        integer :: flag,file,io,close_status

        written=.false.
        select case(model_unit)
        case(unitless)
            flag=PARAM15
            name='UNSPECIFIED'
        case(unit_mm)
            flag=MILLIMETERS
            name='MM'
        case(unit_cm)
            flag=CENTIMETERS
            name='CM'
        case(unit_m)
            flag=METERS
            name='M'
        case(unit_in)
            flag=INCHES
            name='INCH'
        case(unit_ft)
            flag=FEET
            name='FT'
        case default
            return
        end select
        inquire(file=source,size=bytes,iostat=io)
        if(io/=0)return
        if(bytes<=0)return
        ! Preserve ForCAD's start, directory and parameter records. The
        ! Global section describes the existing coordinates in model units.
        model%filename=source
        call model%makeFileImage()
        call model%splitFileImage()
        call model%G%init(filename='geometry.iges',native_id='ForCAD Studio',units_flag=flag,units_name=trim(name))
        call makeGsection(model%G,global_records)
        call makeTsection(model%numSrecs,size(global_records),model%numDrecs,model%numPrecs,termination)
        open(newunit=file,file=destination,status='new',action='write',iostat=io)
        if(io/=0)return
        write(file,'(a)',iostat=io)model%Srecords,global_records,model%Drecords,model%Precords,termination
        close(file,iostat=close_status)
        written=io==0.and.close_status==0
    end function write_iges_units

    logical function write_native_vtk(geometry,path,generate) result(written)
        class(*), intent(inout) :: geometry
        character(*), intent(in) :: path
        logical, intent(in) :: generate

        written=.false.
        select type(geometry)
        type is(nurbs_curve)
            if(generate)call geometry%create(res=curve_resolution)
            call geometry%export_Xg(path,encoding='ASCII')
            written=geometry%err%ok
        type is(nurbs_surface)
            if(generate)call geometry%create(res1=patch_resolution,res2=patch_resolution)
            call geometry%export_Xg(path,encoding='ASCII')
            written=geometry%err%ok
        type is(nurbs_volume)
            if(generate)call geometry%create(res1=patch_resolution,res2=patch_resolution,res3=patch_resolution)
            call geometry%export_Xg(path,encoding='ASCII')
            written=geometry%err%ok
        end select
    end function write_native_vtk

    logical function write_native_multipatch_vtk(geometry,prefix,generate) result(written)
        class(*), intent(inout) :: geometry
        character(*), intent(in) :: prefix
        logical, intent(in) :: generate

        written=.false.
        select type(geometry)
        type is(nurbs_multipatch_curve)
            if(generate)call geometry%create(res=curve_resolution)
            call geometry%export_Xg(prefix,encoding='ASCII',combined=.true.)
            written=geometry%err%ok
        type is(nurbs_multipatch_surface)
            if(generate)call geometry%create(res1=patch_resolution,res2=patch_resolution)
            call geometry%export_Xg(prefix,encoding='ASCII',combined=.true.)
            written=geometry%err%ok
        type is(nurbs_multipatch_volume)
            if(generate)call geometry%create(res1=patch_resolution,res2=patch_resolution,res3=patch_resolution)
            call geometry%export_Xg(prefix,encoding='ASCII',combined=.true.)
            written=geometry%err%ok
        end select
    end function write_native_multipatch_vtk

    logical function write_mesh_vtk(path,kind,parts,connected_points,connected_cells) result(written)
        character(*), intent(in) :: path
        integer, intent(in) :: kind
        type(forcad_mesh_samples), intent(in) :: parts(:)
        real(rk), allocatable, intent(in) :: connected_points(:,:)
        integer, allocatable, intent(in) :: connected_cells(:,:)
        integer, parameter :: vtk_cell_types(3)=[3,9,12]
        real(rk), allocatable :: points(:,:)
        integer, allocatable :: cells(:,:)
        integer :: i,width,coordinates,node_count,cell_count,node_offset,cell_offset
        integer :: patch_cell_counts(size(parts))
        integer(int64) :: total_nodes,total_cells

        written=.false.
        if(kind<1.or.kind>size(vtk_cell_types).or.size(parts)==0)return
        if(allocated(connected_points).neqv.allocated(connected_cells))return
        width=2**kind
        total_nodes=0
        do i=1,size(parts)
            associate(part=>parts(i))
                if(.not.allocated(part%points).or..not.allocated(part%cells))return
                node_count=size(part%points,1)
                coordinates=size(part%points,2)
                if(node_count==0.or.coordinates<1.or.coordinates>3)return
                if(size(part%cells,1)==0.or.size(part%cells,2)/=width)return
                if(any(part%cells<1).or.any(part%cells>node_count))return
                patch_cell_counts(i)=size(part%cells,1)
                total_nodes=total_nodes+node_count
            end associate
        end do
        total_cells=sum(int(patch_cell_counts,int64))
        if(total_nodes>huge(1).or.total_cells>huge(1)/(width+1))return

        if(allocated(connected_points))then
            if(kind/=3.or.size(connected_points,2)/=3.or.size(connected_cells,2)/=width)return
            if(size(connected_cells,1)/=total_cells.or.size(connected_points,1)==0)return
            if(any(connected_cells<1).or.any(connected_cells>size(connected_points,1)))return
            call write_vtk_with_part_ids(path,connected_points,connected_cells,vtk_cell_types(kind),patch_cell_counts)
        else if(size(parts)==1.and.size(parts(1)%points,2)>=2)then
            call export_vtk_legacy(path,parts(1)%points,parts(1)%cells,vtk_cell_types(kind),encoding='ASCII')
        else
            allocate(points(int(total_nodes),3),cells(int(total_cells),width))
            points=0._rk
            node_offset=0
            cell_offset=0
            do i=1,size(parts)
                associate(part=>parts(i))
                    node_count=size(part%points,1)
                    cell_count=patch_cell_counts(i)
                    coordinates=size(part%points,2)
                    points(node_offset+1:node_offset+node_count,:coordinates)=part%points
                    cells(cell_offset+1:cell_offset+cell_count,:)=part%cells+node_offset
                    node_offset=node_offset+node_count
                    cell_offset=cell_offset+cell_count
                end associate
            end do
            if(size(parts)==1)then
                call export_vtk_legacy(path,points,cells,vtk_cell_types(kind),encoding='ASCII')
            else
                call write_vtk_with_part_ids(path,points,cells,vtk_cell_types(kind),patch_cell_counts)
            end if
        end if
        written=.true.
    end function write_mesh_vtk

    subroutine write_vtk_with_part_ids(path,points,cells,cell_type,patch_cell_counts)
        character(*), intent(in) :: path
        real(rk), contiguous, intent(in) :: points(:,:)
        integer, contiguous, intent(in) :: cells(:,:)
        integer, intent(in) :: cell_type,patch_cell_counts(:)
        character(8), parameter :: field_names(1)=['patch_id']
        real(rk), allocatable :: patch_ids(:,:)
        integer :: offset,i,count

        allocate(patch_ids(size(cells,1),1))
        offset=0
        do i=1,size(patch_cell_counts)
            count=patch_cell_counts(i)
            patch_ids(offset+1:offset+count,1)=real(i,rk)
            offset=offset+count
        end do
        call export_vtk_legacy(path,points,cells,cell_type,encoding='ASCII',cell_data=patch_ids,cell_field_names=field_names)
    end subroutine write_vtk_with_part_ids

end module studio_forcad_exchange
