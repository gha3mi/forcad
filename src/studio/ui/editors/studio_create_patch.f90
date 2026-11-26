! Native dialog for adding an independent patch to the current ForCAD Studio session.
module studio_create_patch
    use studio_status_hints, only: set_status_hint, publish_window_status
    use, intrinsic :: iso_c_binding, only: c_double
    use qt, only: QButtonGroup, QComboBox, QDialog, QDialog_Accepted, QDialogButtonBox, QDialogButtonBox_AcceptRole, &
        QDialogButtonBox_Cancel, QDoubleSpinBox, QGridLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, &
        QSizePolicy_Fixed, QSizePolicy_Ignored, QVBoxLayout, QWidget
    use studio_units, only: unitless,length_caption,model_length
    use studio_single_patch, only: add_single_patch, single_patch_count, single_patch_name, active_patch_id, &
        single_patch_available, studio_patch_catalog
    use studio_theme, only: studio_style_primary
    implicit none(type, external)
    private
    public :: create_studio_patch

    character(len=7), parameter :: geometry_names(3)=[character(len=7)::'Curve','Surface','Volume']

    type :: patch_creator
        type(QDialog) :: dialog
        type(QComboBox) :: source
        type(QLineEdit) :: name
        type(QPushButton) :: kinds(3),add
        type(QDoubleSpinBox) :: offset(3)
        type(studio_patch_catalog), pointer :: catalog=>null()
        character(len=:), allocatable :: suggested_name
        integer,allocatable :: source_variants(:)
        integer :: kind=0,created_kind=0,model_unit=unitless,display_unit=unitless
    end type patch_creator

contains

    integer function create_studio_patch(parent,initial_kind,model_unit,display_unit,catalog) result(kind)
        class(QWidget), intent(in) :: parent
        integer, intent(in) :: initial_kind
        integer,intent(in),optional::model_unit,display_unit
        type(studio_patch_catalog), intent(inout), target, optional :: catalog
        type(patch_creator), target :: editor
        type(QVBoxLayout) :: layout
        type(QHBoxLayout) :: placement,types
        type(QButtonGroup) :: group
        type(QGridLayout) :: form
        type(QLabel) :: label
        type(QDialogButtonBox) :: buttons
        type(QPushButton) :: cancel
        integer :: rank,axis
        character(3), parameter :: axes='XYZ'

        kind=0
        if(present(catalog))editor%catalog=>catalog
        if(present(model_unit))editor%model_unit=model_unit
        if(present(display_unit))editor%display_unit=display_unit
        editor%dialog=QDialog(parent)
        call editor%dialog%setObjectName('studioCreatePatchDialog')
        call editor%dialog%setWindowTitle('Copy patch')
        call editor%dialog%setMinimumWidth(464)
        layout=QVBoxLayout(editor%dialog)
        call layout%setContentsMargins(18,16,18,16)
        call layout%setSpacing(10)

        label=QLabel('Copy patch',editor%dialog)
        call label%setObjectName('studioDialogTitle')
        call layout%addWidget(label)
        call set_status_hint(editor%dialog,'Choose a source patch and place the copy with an optional offset.')

        form=QGridLayout()
        call form%setSpacing(8)
        call form%setColumnStretch(1,1)
        label=QLabel('Geometry type',editor%dialog)
        call form%addWidget(label,0,0)
        group=QButtonGroup(editor%dialog)
        call group%setExclusive(.true.)
        types=QHBoxLayout()
        call types%setSpacing(4)
        do rank=1,size(geometry_names)
            editor%kinds(rank)=QPushButton(trim(geometry_names(rank)),editor%dialog)
            call editor%kinds(rank)%setObjectName('studioCreatePatch'//trim(geometry_names(rank)))
            call editor%kinds(rank)%setCheckable(.true.)
            call editor%kinds(rank)%setProperty('autoDefault',.false.)
            call group%addButton(editor%kinds(rank),rank)
            call types%addWidget(editor%kinds(rank),1)
        end do
        call form%addLayout(types,0,1)

        label=QLabel('Source patch',editor%dialog)
        call form%addWidget(label,1,0)
        editor%source=QComboBox(editor%dialog)
        call editor%source%setObjectName('studioCreatePatchSource')
        call editor%source%setProperty('accessibleName','Source patch to copy')
        call editor%source%setMinimumHeight(30)
        call editor%source%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Fixed)
        call form%addWidget(editor%source,1,1)

        label=QLabel('Name',editor%dialog)
        call form%addWidget(label,2,0)
        editor%name=QLineEdit(parent=editor%dialog)
        call editor%name%setObjectName('studioCreatePatchName')
        call editor%name%setProperty('accessibleName','New patch name')
        call editor%name%setPlaceholderText('Enter a patch name')
        call editor%name%setMinimumHeight(30)
        call form%addWidget(editor%name,2,1)

        label=QLabel(length_caption('Offset',editor%display_unit),editor%dialog)
        call set_status_hint(label,'Offset from the source patch along the world axes.')
        call form%addWidget(label,3,0)
        placement=QHBoxLayout()
        call placement%setSpacing(6)
        do axis=1,size(editor%offset)
            editor%offset(axis)=QDoubleSpinBox(editor%dialog)
            call editor%offset(axis)%setObjectName('studioCreatePatchOffset'//axes(axis:axis))
            call editor%offset(axis)%setProperty('accessibleName',axes(axis:axis)//' placement offset')
            call editor%offset(axis)%setPrefix(axes(axis:axis)//'  ')
            call editor%offset(axis)%setDecimals(4)
            call editor%offset(axis)%setRange(-huge(1._c_double),huge(1._c_double))
            call editor%offset(axis)%setSingleStep(.25_c_double)
            call editor%offset(axis)%setMinimumHeight(30)
            call editor%offset(axis)%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Fixed)
            call placement%addWidget(editor%offset(axis),1)
        end do
        call form%addLayout(placement,3,1)
        call layout%addLayout(form)


        buttons=QDialogButtonBox(QDialogButtonBox_Cancel,editor%dialog)
        editor%add=buttons%addButton('Copy patch',QDialogButtonBox_AcceptRole)
        call editor%add%setObjectName('studioCreatePatchApply')
        call studio_style_primary(editor%add)
        call editor%add%setDefault(.true.)
        cancel=buttons%button(QDialogButtonBox_Cancel)
        call cancel%setObjectName('studioCreatePatchCancel')
        call cancel%setMinimumWidth(80)
        call layout%addWidget(buttons)

        call update_geometry(editor,max(1,min(size(geometry_names),initial_kind)))
        call group%idClicked%connect(geometry_changed,editor)
        call editor%source%currentIndexChanged%connect(source_changed,editor)
        call editor%name%textChanged%connect(name_changed,editor)
        call buttons%accepted%connect(add_patch,editor)
        call buttons%rejected%connect(cancel_patch,editor)
        call editor%name%setFocus()

        if (editor%dialog%exec()==QDialog_Accepted) kind=editor%created_kind
        call editor%dialog%delete()
    end function create_studio_patch

    subroutine update_geometry(editor,kind)
        type(patch_creator), intent(inout) :: editor
        integer, intent(in) :: kind
        integer :: variant,count
        logical :: use_suggestion

        if (kind==editor%kind.or.kind<1.or.kind>size(geometry_names)) return
        editor%kind=kind
        use_suggestion=.not.allocated(editor%suggested_name)
        if (allocated(editor%suggested_name)) use_suggestion=editor%name%text()==editor%suggested_name
        call editor%kinds(kind)%setChecked(.true.)
        call editor%source%clear()
        count=single_patch_count(kind,editor%catalog)
        editor%source_variants=pack([(variant,variant=1,count)], &
            [(single_patch_available(kind,variant,editor%catalog),variant=1,count)])
        do variant=1,size(editor%source_variants)
            call editor%source%addItem(single_patch_name(kind,editor%source_variants(variant),editor%catalog))
        end do
        call editor%source%setEnabled(size(editor%source_variants)>0)
        call editor%add%setEnabled(size(editor%source_variants)>0)
        if(size(editor%source_variants)>0)then
            call editor%source%setCurrentIndex(max(0,findloc(editor%source_variants,active_patch_id(kind,editor%catalog),dim=1)-1))
        end if
        editor%suggested_name=next_patch_name(kind,editor%catalog)
        if (use_suggestion) call editor%name%setText(editor%suggested_name)
        if(size(editor%source_variants)==0)then
            call publish_window_status(editor%dialog,'There are no patches of this type to copy.',warning=.true.,blocking=.true.)
        else
            call publish_window_status(editor%dialog,'',blocking=.true.)
        end if
    end subroutine update_geometry

    function next_patch_name(kind,catalog) result(name)
        integer, intent(in) :: kind
        type(studio_patch_catalog), intent(in), optional :: catalog
        character(len=:), allocatable :: name
        character(len=24) :: number
        integer :: serial,rank,variant,catalog_size

        catalog_size=sum([(single_patch_count(rank,catalog),rank=1,size(geometry_names))])
        candidate: do serial=1,catalog_size+1
            write(number,'(i0)') serial
            name=trim(geometry_names(kind))//' '//trim(number)
            do rank=1,size(geometry_names)
                do variant=1,single_patch_count(rank,catalog)
                    if (single_patch_name(rank,variant,catalog)==name) cycle candidate
                end do
            end do
            return
        end do candidate
    end function next_patch_name

    subroutine geometry_changed(kind,data)
        integer, intent(in) :: kind
        class(*), intent(inout), optional :: data

        select type(editor=>data)
        type is(patch_creator)
            call update_geometry(editor,kind)
        end select
    end subroutine geometry_changed

    subroutine source_changed(index,data)
        integer, intent(in) :: index
        class(*), intent(inout), optional :: data

        select type(editor=>data)
        type is(patch_creator)
            call publish_window_status(editor%dialog,'',blocking=.true.)
        end select
    end subroutine source_changed

    subroutine name_changed(text,data)
        character(len=*), intent(in) :: text
        class(*), intent(inout), optional :: data

        select type(editor=>data)
        type is(patch_creator)
            if(size(editor%source_variants)>0)call publish_window_status(editor%dialog,'',blocking=.true.)
        end select
    end subroutine name_changed

    subroutine add_patch(data)
        class(*), intent(inout), optional :: data
        character(len=:), allocatable :: error
        integer :: kind,source,variant,axis
        real(c_double) :: offset(3)

        select type(editor=>data)
        type is(patch_creator)
            kind=editor%kind
            source=editor%source%currentIndex()+1
            if(source<1.or.source>size(editor%source_variants))return
            source=editor%source_variants(source)
            offset=[(editor%offset(axis)%value(),axis=1,size(editor%offset))]
            call add_single_patch(kind,source,editor%name%text(),variant,error, &
                model_length(offset,editor%model_unit,editor%display_unit),catalog=editor%catalog)
            if (variant==0) then
                if (len(error)==0) error='Could not add the patch.'
                call publish_window_status(editor%dialog,error,warning=.true.,blocking=.true.)
                call editor%name%setFocus()
                return
            end if
            editor%created_kind=kind
            call editor%dialog%accept()
        end select
    end subroutine add_patch

    subroutine cancel_patch(data)
        class(*), intent(inout), optional :: data

        select type(editor=>data)
        type is(patch_creator)
            call editor%dialog%reject()
        end select
    end subroutine cancel_patch

end module studio_create_patch
