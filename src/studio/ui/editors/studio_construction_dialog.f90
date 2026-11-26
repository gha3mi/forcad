! Native Qt editors for the four ForCAD construction operations.
module studio_construction_dialog
    use studio_status_hints, only: set_status_hint, publish_window_status
    use qt, only: QButtonGroup, QComboBox, QDialog, QDialog_Accepted, QDialogButtonBox, QDialogButtonBox_AcceptRole, &
        QDialogButtonBox_Cancel, QGridLayout, QHBoxLayout, QLabel, QLineEdit, QListWidget, QPushButton, &
        QRegularExpressionValidator, QSizePolicy_Fixed, QSizePolicy_Ignored, QSpinBox, QVBoxLayout, QWidget
    use studio_icons, only: line_icon
    use studio_theme, only: studio_frame_edge, studio_selected_fill, studio_style_primary
    use studio_forcad, only: rk, extrude_operation, revolve_operation, sweep_operation, loft_operation
    use studio_single_patch, only: studio_patch_catalog, single_patch_available, single_patch_count, single_patch_name
    use studio_parameter_editor, only: parse_parameter_number,parameter_number_pattern
    use studio_construction, only: construction_names,publish_construction
    use studio_units, only: unitless,model_length,length_caption
    implicit none(type, external)
    private
    public :: create_studio_construction

    character(len=3),parameter :: axes='XYZ'

    type :: construction_editor
        type(QDialog) :: dialog
        type(QComboBox) :: source,spine
        type(QListWidget) :: sections
        type(QSpinBox) :: degree
        type(QLineEdit) :: name,vector(3),axis_point(3),angle
        type(QLabel) :: result_label
        type(QPushButton) :: kinds(2),apply,add_section,remove_section,up,down
        type(QRegularExpressionValidator) :: decimal
        type(studio_patch_catalog), pointer :: catalog=>null()
        integer,allocatable :: source_ids(:),spine_ids(:),section_ids(:)
        integer :: operation=0,profile_kind=0,created_kind=0,first_variant=1,model_unit=unitless,display_unit=unitless
        logical :: updating=.false.
    end type construction_editor
contains
    integer function create_studio_construction(parent,operation,active_kind,active_variant,model_unit,display_unit,first_variant,catalog) result(kind)
        class(QWidget),intent(in) :: parent
        integer,intent(in) :: operation,active_kind,active_variant
        integer,intent(in),optional :: model_unit,display_unit,first_variant
        type(studio_patch_catalog), intent(inout), target, optional :: catalog
        type(construction_editor),target :: editor
        type(QVBoxLayout) :: body
        type(QHBoxLayout) :: header,types
        type(QButtonGroup) :: group
        type(QGridLayout) :: form
        type(QLabel) :: heading
        type(QDialogButtonBox) :: buttons
        character(:),allocatable :: title,instructions
        integer :: start_kind,i
        character(7),parameter :: profile_types(2)=[character(7)::'Curve','Surface']

        kind=0
        if(operation<1.or.operation>size(construction_names))return
        editor%operation=operation
        if(present(catalog))editor%catalog=>catalog
        if(present(first_variant))editor%first_variant=first_variant
        if(present(model_unit))editor%model_unit=model_unit
        if(present(display_unit))editor%display_unit=display_unit
        title=trim(construction_names(operation))
        start_kind=1
        if(active_kind==2)start_kind=2
        allocate(editor%section_ids(0))
        editor%dialog=QDialog(parent)
        call editor%dialog%setObjectName('studioConstructionDialog')
        call editor%dialog%setWindowTitle(title)
        call editor%dialog%setMinimumWidth(490)
        editor%decimal=QRegularExpressionValidator(parameter_number_pattern,editor%dialog)
        body=QVBoxLayout(editor%dialog)
        call body%setContentsMargins(18,16,18,16)
        call body%setSpacing(10)

        header=QHBoxLayout()
        heading=QLabel(title,editor%dialog)
        call heading%setObjectName('studioDialogTitle')
        call header%addWidget(heading,1)
        editor%result_label=QLabel('',editor%dialog)
        call editor%result_label%setObjectName('studioConstructionResult')
        call editor%result_label%setStyleSheet('color:#d4eeeb; background:'//studio_selected_fill//'; padding:5px 9px;')
        call header%addWidget(editor%result_label)
        call body%addLayout(header)
        select case(operation)
        case(extrude_operation)
            instructions='Extend the profile by an exact XYZ displacement.'
        case(revolve_operation)
            instructions='Rotate the profile around an axis. Enter a signed angle up to one full turn.'
        case(sweep_operation)
            instructions='Translate the profile from the path start while keeping its orientation fixed.'
        case(loft_operation)
            instructions='Interpolate compatible sections in the order shown, from first to last.'
        end select
        call set_status_hint(editor%dialog,instructions)

        form=QGridLayout()
        call form%setColumnStretch(1,1)
        call form%setSpacing(8)
        group=QButtonGroup(editor%dialog)
        call group%setExclusive(.true.)
        types=QHBoxLayout()
        call types%setSpacing(4)
        do i=1,size(profile_types)
            editor%kinds(i)=QPushButton(trim(profile_types(i)),editor%dialog)
            call editor%kinds(i)%setObjectName('studioConstruction'//trim(profile_types(i)))
            call editor%kinds(i)%setCheckable(.true.)
            call editor%kinds(i)%setProperty('autoDefault',.false.)
            call group%addButton(editor%kinds(i),i)
            call types%addWidget(editor%kinds(i),1)
        end do
        call form%addWidget(QLabel('Input type',editor%dialog),0,0)
        call form%addLayout(types,0,1)
        editor%name=QLineEdit(parent=editor%dialog)
        call editor%name%setObjectName('studioConstructionName')
        call editor%name%setProperty('accessibleName','Result name')
        call editor%name%setText(next_name(title,editor%catalog))
        call form%addWidget(QLabel('Name',editor%dialog),1,0)
        call form%addWidget(editor%name,1,1)
        editor%source=QComboBox(editor%dialog)
        call editor%source%setObjectName('studioConstructionSource')
        call editor%source%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Fixed)
        if(operation==loft_operation)then
            call add_loft_fields(editor,form)
        else
            call editor%source%setProperty('accessibleName','Profile patch')
            call form%addWidget(QLabel('Profile',editor%dialog),2,0)
            call form%addWidget(editor%source,2,1)
        end if
        select case(operation)
        case(extrude_operation)
            call add_vector(editor%dialog,form,3,length_caption('Displacement',display_unit), &
                'Vector',editor%vector,['0','0','1'])
        case(revolve_operation)
            call add_vector(editor%dialog,form,3,length_caption('Axis point',display_unit), &
                'Axis',editor%axis_point,['0','0','0'])
            call add_vector(editor%dialog,form,4,'Axis direction','Vector',editor%vector,['0','0','1'])
            editor%angle=QLineEdit(parent=editor%dialog)
            call editor%angle%setObjectName('studioConstructionAngle')
            call editor%angle%setProperty('accessibleName','Revolution angle in degrees')
            call editor%angle%setMaxLength(128)
            call editor%angle%setText('90')
            call form%addWidget(QLabel('Angle (°)',editor%dialog),5,0)
            call form%addWidget(editor%angle,5,1)
        case(sweep_operation)
            editor%spine=QComboBox(editor%dialog)
            call editor%spine%setObjectName('studioConstructionSpine')
            call editor%spine%setProperty('accessibleName','Sweep path curve')
            call editor%spine%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Fixed)
            call fill_patch_picker(editor%spine,1,0,editor%spine_ids,editor%first_variant,editor%catalog)
            if(active_kind==1.and.size(editor%spine_ids)>1)then
                if(editor%spine_ids(1)==active_variant)call editor%spine%setCurrentIndex(1)
            end if
            call form%addWidget(QLabel('Path',editor%dialog),3,0)
            call form%addWidget(editor%spine,3,1)
        end select
        call body%addLayout(form)
        buttons=QDialogButtonBox(QDialogButtonBox_Cancel,editor%dialog)
        editor%apply=buttons%addButton(title,QDialogButtonBox_AcceptRole)
        call editor%apply%setIcon(line_icon(76+operation))
        call editor%apply%setIconSize(17,17)
        call editor%apply%setObjectName('studioConstructionApply')
        call editor%apply%setDefault(.true.)
        call studio_style_primary(editor%apply)
        call body%addWidget(buttons)

        call load_profile_kind(editor,start_kind,merge(active_variant,0,active_kind==start_kind))
        if(operation==loft_operation.and.active_kind==start_kind)then
            if(single_patch_available(start_kind,active_variant,editor%catalog))then
                editor%section_ids=[active_variant]
                call refresh_sections(editor,0)
            end if
        end if
        call group%idClicked%connect(on_kind_changed,editor)
        if(operation==loft_operation)then
            call editor%source%currentIndexChanged%connect(on_picker_changed,editor)
            call editor%sections%currentRowChanged%connect(on_picker_changed,editor)
            call editor%add_section%clicked%connect(on_add_section,editor)
            call editor%remove_section%clicked%connect(on_remove_section,editor)
            call editor%up%clicked%connect(on_move_up,editor)
            call editor%down%clicked%connect(on_move_down,editor)
        end if
        call editor%name%textChanged%connect(input_changed,editor)
        if(operation==extrude_operation.or.operation==revolve_operation)then
            do i=1,3
                call editor%vector(i)%textChanged%connect(input_changed,editor)
            end do
        end if
        if(operation==revolve_operation)then
            do i=1,3
                call editor%axis_point(i)%textChanged%connect(input_changed,editor)
            end do
            call editor%angle%textChanged%connect(input_changed,editor)
        end if
        call buttons%accepted%connect(on_create,editor)
        call buttons%rejected%connect(on_cancel,editor)
        if(editor%dialog%exec()==QDialog_Accepted)kind=editor%created_kind
        call editor%dialog%delete()
    end function create_studio_construction

    subroutine add_vector(parent,form,row,caption,prefix,fields,defaults)
        class(QWidget),intent(in) :: parent
        type(QGridLayout),intent(inout) :: form
        integer,intent(in) :: row
        character(*),intent(in) :: caption,prefix,defaults(3)
        type(QLineEdit),intent(out) :: fields(3)
        type(QHBoxLayout) :: values
        type(QLabel) :: axis_label
        integer :: i

        values=QHBoxLayout()
        call values%setSpacing(5)
        do i=1,3
            axis_label=QLabel(axes(i:i),parent)
            call axis_label%setProperty('dialogRole','detail')
            call values%addWidget(axis_label)
            fields(i)=QLineEdit(parent=parent)
            call fields(i)%setObjectName('studioConstruction'//prefix//axes(i:i))
            call fields(i)%setProperty('accessibleName',caption//' '//axes(i:i))
            call fields(i)%setMaxLength(128)
            call fields(i)%setMinimumHeight(27)
            call fields(i)%setText(defaults(i))
            call values%addWidget(fields(i),1)
        end do
        call form%addWidget(QLabel(caption,parent),row,0)
        call form%addLayout(values,row,1)
    end subroutine add_vector

    subroutine add_loft_fields(editor,form)
        type(construction_editor),intent(inout) :: editor
        type(QGridLayout),intent(inout) :: form
        type(QHBoxLayout) :: picker,list_row,degree_row
        type(QVBoxLayout) :: tools
        type(QLabel) :: automatic

        picker=QHBoxLayout()
        call editor%source%setProperty('accessibleName','Available section patch')
        call picker%addWidget(editor%source,1)
        editor%add_section=QPushButton('Add')
        call editor%add_section%setProperty('autoDefault',.false.)
        call editor%add_section%setObjectName('studioConstructionAddSection')
        call editor%add_section%setProperty('accessibleName','Add the selected patch to the section order')
        call picker%addWidget(editor%add_section)
        call form%addWidget(QLabel('Section',editor%dialog),2,0)
        call form%addLayout(picker,2,1)

        editor%sections=QListWidget(editor%dialog)
        call editor%sections%setObjectName('studioConstructionSections')
        call editor%sections%setStyleSheet('background:#141f27; border:1px solid '//studio_frame_edge//'; padding:5px;')
        call editor%sections%setProperty('accessibleName','Loft sections in order')
        call editor%sections%setMinimumHeight(126)
        call editor%sections%setMaximumHeight(150)
        editor%up=QPushButton('Move up')
        call editor%up%setProperty('autoDefault',.false.)
        call editor%up%setObjectName('studioConstructionSectionUp')
        editor%down=QPushButton('Move down')
        call editor%down%setProperty('autoDefault',.false.)
        call editor%down%setObjectName('studioConstructionSectionDown')
        editor%remove_section=QPushButton('Remove')
        call editor%remove_section%setProperty('autoDefault',.false.)
        call editor%remove_section%setObjectName('studioConstructionRemoveSection')
        tools=QVBoxLayout()
        call tools%setSpacing(4)
        call tools%addWidget(editor%up)
        call tools%addWidget(editor%down)
        call tools%addWidget(editor%remove_section)
        call tools%addStretch()
        list_row=QHBoxLayout()
        call list_row%addWidget(editor%sections,1)
        call list_row%addLayout(tools)
        call form%addWidget(QLabel('Order',editor%dialog),3,0)
        call form%addLayout(list_row,3,1)

        editor%degree=QSpinBox(editor%dialog)
        call editor%degree%setObjectName('studioConstructionDegree')
        call editor%degree%setProperty('accessibleName','Loft degree, zero for automatic')
        call editor%degree%setRange(0,1)
        call editor%degree%setValue(0)
        automatic=QLabel('0 = Automatic',editor%dialog)
        call automatic%setProperty('dialogRole','detail')
        degree_row=QHBoxLayout()
        call degree_row%addWidget(editor%degree)
        call degree_row%addWidget(automatic,1)
        call form%addWidget(QLabel('Loft degree',editor%dialog),4,0)
        call form%addLayout(degree_row,4,1)
    end subroutine add_loft_fields

    subroutine fill_patch_picker(picker,kind,preferred,ids,first_variant,catalog)
        type(QComboBox),intent(inout) :: picker
        integer,intent(in) :: kind,preferred,first_variant
        integer,allocatable,intent(out) :: ids(:)
        type(studio_patch_catalog), intent(in), optional :: catalog
        integer :: i,n

        n=single_patch_count(kind,catalog)
        ids=pack([(i,i=first_variant,n)],[(single_patch_available(kind,i,catalog),i=first_variant,n)])
        call picker%clear()
        do i=1,size(ids)
            call picker%addItem(single_patch_name(kind,ids(i),catalog))
        end do
        if(size(ids)>0)call picker%setCurrentIndex(max(0,findloc(ids,preferred,dim=1)-1))
    end subroutine fill_patch_picker

    integer function picked_id(picker,ids) result(id)
        type(QComboBox),intent(in) :: picker
        integer,intent(in) :: ids(:)
        integer :: row

        row=picker%currentIndex()+1
        id=0
        if(row>=1.and.row<=size(ids))id=ids(row)
    end function picked_id

    subroutine load_profile_kind(editor,kind,preferred)
        type(construction_editor),intent(inout) :: editor
        integer,intent(in) :: kind,preferred
        logical :: has_profile,has_path

        if(kind==editor%profile_kind)return
        editor%profile_kind=kind
        editor%updating=.true.
        call editor%kinds(kind)%setChecked(.true.)
        call fill_patch_picker(editor%source,kind,preferred,editor%source_ids,editor%first_variant,editor%catalog)
        call publish_window_status(editor%dialog,'',blocking=.true.)
        if(kind==1)then
            call editor%result_label%setText('Surface result')
        else
            call editor%result_label%setText('Volume result')
        end if
        select case(editor%operation)
        case(loft_operation)
            editor%section_ids=[integer::]
            call refresh_sections(editor,-1)
        case(sweep_operation)
            has_profile=size(editor%source_ids)>0
            has_path=size(editor%spine_ids)>0
            call editor%apply%setEnabled(has_profile.and.has_path)
            if(.not.has_profile.and..not.has_path)then
                call show_error(editor,'Add a profile and a curve to use as the sweep path.')
            else if(.not.has_path)then
                call show_error(editor,'Add a curve to use as the sweep path.')
            else if(.not.has_profile)then
                call show_error(editor,'Add a profile of the selected type or choose another input type.')
            end if
        case default
            call editor%apply%setEnabled(size(editor%source_ids)>0)
            if(size(editor%source_ids)==0)call show_error(editor,'Add a profile of the selected type or choose another input type.')
        end select
        editor%updating=.false.
    end subroutine load_profile_kind

    subroutine refresh_sections(editor,row)
        type(construction_editor),intent(inout) :: editor
        integer,intent(in) :: row
        integer :: i,kind
        character(len=12) :: number

        kind=editor%profile_kind
        call editor%sections%clear()
        do i=1,size(editor%section_ids)
            write(number,'(i0)')i
            call editor%sections%addItem(trim(number)//'  '//single_patch_name(kind,editor%section_ids(i),editor%catalog))
        end do
        call editor%sections%setCurrentRow(row)
        call editor%degree%setRange(0,max(0,size(editor%section_ids)-1))
        call update_section_tools(editor)
    end subroutine refresh_sections

    subroutine update_section_tools(editor)
        type(construction_editor),intent(inout) :: editor
        integer :: row,id,n

        id=picked_id(editor%source,editor%source_ids)
        row=editor%sections%currentRow()
        n=size(editor%section_ids)
        call editor%add_section%setEnabled(id>0.and..not.any(editor%section_ids==id))
        call editor%remove_section%setEnabled(row>=0.and.row<n)
        call editor%up%setEnabled(row>0.and.row<n)
        call editor%down%setEnabled(row>=0.and.row<n-1)
        call editor%apply%setEnabled(n>=2)
    end subroutine update_section_tools

    function next_name(prefix,catalog) result(name)
        character(*),intent(in) :: prefix
        type(studio_patch_catalog), intent(in), optional :: catalog
        character(:),allocatable :: name
        character(len=20) :: number
        integer :: rank,item,serial

        serial=1
        names: do
            write(number,'(i0)')serial
            name=prefix//' '//trim(number)
            do rank=1,3
                do item=1,single_patch_count(rank,catalog)
                    if(single_patch_name(rank,item,catalog)/=name)cycle
                    serial=serial+1
                    cycle names
                end do
            end do
            return
        end do names
    end function next_name

    subroutine on_kind_changed(kind,data)
        integer,intent(in) :: kind
        class(*),intent(inout),optional :: data

        select type(editor=>data)
        type is(construction_editor)
            call load_profile_kind(editor,kind,0)
        end select
    end subroutine on_kind_changed

    subroutine on_picker_changed(index,data)
        integer,intent(in) :: index
        class(*),intent(inout),optional :: data

        select type(editor=>data)
        type is(construction_editor)
            if(.not.editor%updating)call update_section_tools(editor)
        end select
    end subroutine on_picker_changed

    subroutine on_add_section(data)
        class(*),intent(inout),optional :: data
        integer :: id

        select type(editor=>data)
        type is(construction_editor)
            id=picked_id(editor%source,editor%source_ids)
            if(id<=0.or.any(editor%section_ids==id))return
            editor%section_ids=[editor%section_ids,id]
            call publish_window_status(editor%dialog,'',blocking=.true.)
            call refresh_sections(editor,size(editor%section_ids)-1)
        end select
    end subroutine on_add_section

    subroutine on_remove_section(data)
        class(*),intent(inout),optional :: data
        integer :: row,n

        select type(editor=>data)
        type is(construction_editor)
            row=editor%sections%currentRow()+1
            n=size(editor%section_ids)
            if(row<1.or.row>n)return
            editor%section_ids=[editor%section_ids(:row-1),editor%section_ids(row+1:)]
            call publish_window_status(editor%dialog,'',blocking=.true.)
            call refresh_sections(editor,min(row-1,n-2))
        end select
    end subroutine on_remove_section

    subroutine on_move_up(data)
        class(*),intent(inout),optional :: data

        select type(editor=>data)
        type is(construction_editor)
            call move_section(editor,-1)
        end select
    end subroutine on_move_up

    subroutine on_move_down(data)
        class(*),intent(inout),optional :: data

        select type(editor=>data)
        type is(construction_editor)
            call move_section(editor,1)
        end select
    end subroutine on_move_down

    subroutine move_section(editor,step)
        type(construction_editor),intent(inout) :: editor
        integer,intent(in) :: step
        integer :: row,destination,id,n

        row=editor%sections%currentRow()+1
        destination=row+step
        n=size(editor%section_ids)
        if(row<1.or.row>n.or.destination<1.or.destination>n)return
        id=editor%section_ids(row)
        editor%section_ids(row)=editor%section_ids(destination)
        editor%section_ids(destination)=id
        call publish_window_status(editor%dialog,'',blocking=.true.)
        call refresh_sections(editor,destination-1)
    end subroutine move_section

    logical function read_vector(decimal,fields,caption,values,error) result(valid)
        type(QRegularExpressionValidator),intent(in) :: decimal
        type(QLineEdit),intent(inout) :: fields(3)
        character(*),intent(in) :: caption
        real(rk),intent(out) :: values(3)
        character(:),allocatable,intent(out) :: error
        integer :: i

        valid=.false.
        error=''
        do i=1,3
            if(parse_parameter_number(decimal,fields(i)%text(),values(i)))cycle
            error='Enter a finite '//axes(i:i)//' value for the '//caption//'.'
            call fields(i)%setFocus()
            return
        end do
        valid=.true.
    end function read_vector

    subroutine input_changed(text,data)
        character(*), intent(in) :: text
        class(*), optional, intent(inout) :: data

        select type(editor=>data)
        type is(construction_editor)
            if(editor%apply%isEnabled())call publish_window_status(editor%dialog,'',blocking=.true.)
        end select
    end subroutine input_changed

    subroutine on_create(data)
        class(*),intent(inout),optional :: data
        character(:),allocatable :: error,name
        integer :: kind,profile,spine,degree
        real(rk) :: vector(3),point(3),angle

        select type(editor=>data)
        type is(construction_editor)
            call publish_window_status(editor%dialog,'',blocking=.true.)
            name=editor%name%text()
            kind=editor%profile_kind
            profile=picked_id(editor%source,editor%source_ids)
            select case(editor%operation)
            case(extrude_operation)
                if(.not.read_vector(editor%decimal,editor%vector,'displacement',vector,error))then
                    call show_error(editor,error)
                    return
                end if
                vector=model_length(vector,editor%model_unit,editor%display_unit)
                call publish_construction(editor%operation,kind,[profile],name,editor%created_kind,error,vector=vector,catalog=editor%catalog)
            case(revolve_operation)
                if(.not.read_vector(editor%decimal,editor%axis_point,'axis point',point,error))then
                    call show_error(editor,error)
                    return
                end if
                if(.not.read_vector(editor%decimal,editor%vector,'axis direction',vector,error))then
                    call show_error(editor,error)
                    return
                end if
                if(.not.parse_parameter_number(editor%decimal,editor%angle%text(),angle))then
                    call show_error(editor,'Enter a finite revolution angle in degrees.')
                    call editor%angle%setFocus()
                    return
                end if
                point=model_length(point,editor%model_unit,editor%display_unit)
                call publish_construction(editor%operation,kind,[profile],name,editor%created_kind,error, &
                    vector=vector,axis_point=point,angle_degrees=angle,catalog=editor%catalog)
            case(sweep_operation)
                spine=picked_id(editor%spine,editor%spine_ids)
                call publish_construction(editor%operation,kind,[profile],name,editor%created_kind,error,spine=spine,catalog=editor%catalog)
            case(loft_operation)
                degree=editor%degree%value()
                if(degree>0)then
                    call publish_construction(editor%operation,kind,editor%section_ids,name,editor%created_kind,error,degree=degree,catalog=editor%catalog)
                else
                    call publish_construction(editor%operation,kind,editor%section_ids,name,editor%created_kind,error,catalog=editor%catalog)
                end if
            end select
            if(editor%created_kind>0)then
                call editor%dialog%accept()
            else
                call show_error(editor,error)
            end if
        end select
    end subroutine on_create

    subroutine show_error(editor,error)
        type(construction_editor),intent(inout) :: editor
        character(*),intent(in) :: error

        call publish_window_status(editor%dialog,error,warning=.true.,blocking=.true.)
    end subroutine show_error

    subroutine on_cancel(data)
        class(*),intent(inout),optional :: data
        select type(editor=>data)
        type is(construction_editor)
            call editor%dialog%reject()
        end select
    end subroutine on_cancel
end module studio_construction_dialog
