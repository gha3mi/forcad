module studio_inspector
    use studio_status_hints, only: set_status_hint
    use studio_command_presentation, only: command_visual
    use studio_commands, only: diagnose_with_ai_command
    use, intrinsic :: iso_c_binding, only: dp=>c_double
    use qt, only: QDoubleSpinBox, QGridLayout, QHBoxLayout, QLabel, QPushButton, QScrollArea, QVBoxLayout, QWidget, &
        QFrame_NoFrame, Qt_AlignHCenter, Qt_AlignLeft, Qt_AlignTop, Qt_PlainText, Qt_ScrollBarAlwaysOff, Qt_StrongFocus
    use studio_icons, only: line_icon, icon_chevron, icon_chevron_right, icon_refresh
    implicit none(type, external)
    private

    public :: studio_inspector_t, build_inspector_widgets, update_inspector_section

    character(len=24), parameter :: section_titles(2)=[character(len=24)::'Selection and probe','Diagnostics']

    type :: studio_inspector_t
        type(QPushButton) :: headers(2), reset_probe, edit_control, repair_diagnostic
        type(QWidget) :: sections(2), bodies(2), selected_control_row
        type(QLabel) :: context_icon, context_identity, context_dimensions
        type(QLabel) :: probe_v_label, probe_w_label, selected_control_label
        integer :: context_kind=0
    end type studio_inspector_t

contains

    subroutine build_context_card(parent,view)
        type(QVBoxLayout),intent(inout)::parent
        type(studio_inspector_t),intent(inout)::view
        type(QWidget)::wrapper,card
        type(QHBoxLayout)::outer,row
        type(QVBoxLayout)::identity
        wrapper=QWidget(); outer=QHBoxLayout(wrapper)
        call outer%setContentsMargins(8,6,8,2); call outer%setSpacing(0)
        card=QWidget(); call card%setObjectName('entityContext')
        row=QHBoxLayout(card); call row%setContentsMargins(10,6,10,6); call row%setSpacing(8)
        view%context_icon=QLabel(''); call view%context_icon%setFixedSize(28,28)
        call view%context_icon%setScaledContents(.true.)
        call view%context_icon%setAccessibleName('Selected geometry icon')
        call row%addWidget(view%context_icon)
        identity=QVBoxLayout(); call identity%setContentsMargins(0,0,0,0); call identity%setSpacing(1)
        view%context_identity=QLabel(''); call view%context_identity%setObjectName('contextIdentity')
        call view%context_identity%setWordWrap(.true.)
        call view%context_identity%setTextFormat(Qt_PlainText)
        call view%context_identity%setAccessibleName('Selected geometry')
        call set_status_hint(view%context_identity,'Select an object in the 3D view or Scene tree')
        call identity%addWidget(view%context_identity)
        view%context_dimensions=QLabel(''); call view%context_dimensions%setObjectName('contextDimensions')
        call view%context_dimensions%setWordWrap(.true.)
        call view%context_dimensions%setTextFormat(Qt_PlainText)
        call view%context_dimensions%setAccessibleName('Selected object details')
        call identity%addWidget(view%context_dimensions); call row%addLayout(identity,1)
        call outer%addWidget(card); call parent%addWidget(wrapper)
    end subroutine build_context_card

    subroutine build_inspector_widgets(layout,position,view,probe_u,probe_v,probe_w,probe_label,diagnostic_label)
        type(QVBoxLayout),intent(inout)::layout
        real(dp),intent(in)::position(3)
        type(studio_inspector_t),intent(out)::view
        type(QDoubleSpinBox),intent(out)::probe_u,probe_v,probe_w
        type(QLabel),intent(out)::probe_label,diagnostic_label
        type(QScrollArea)::scroll
        type(QWidget)::page
        type(QVBoxLayout)::content,section
        type(QHBoxLayout)::probe_heading
        type(QGridLayout)::probe
        type(QLabel)::label,u_axis
        call build_context_card(layout,view)
        page=QWidget(); content=QVBoxLayout(page)
        call content%setContentsMargins(8,4,8,8); call content%setSpacing(6)
        call add_section(content,section,view,1,.true.)
        probe_heading=QHBoxLayout()
        call probe_heading%setContentsMargins(10,4,10,0)
        call probe_heading%setSpacing(6)
        label=QLabel('Normalized inputs · 0–1')
        call label%setObjectName('studioNormalizedProbe')
        call probe_heading%addWidget(label,1)
        call section%addLayout(probe_heading)
        probe=QGridLayout()
        call probe%setContentsMargins(10,0,10,6)
        call probe%setSpacing(4)
        call build_probe_axis(probe,'U',0,position(1),u_axis,probe_u)
        call build_probe_axis(probe,'V',2,position(2),view%probe_v_label,probe_v)
        call build_probe_axis(probe,'W',4,position(3),view%probe_w_label,probe_w)
        call probe%setColumnStretch(6,1)
        call section%addLayout(probe)
        view%reset_probe=QPushButton('Reset')
        call view%reset_probe%setObjectName('probeReset')
        call view%reset_probe%setFixedHeight(26)
        call view%reset_probe%setIcon(line_icon(icon_refresh))
        call view%reset_probe%setIconSize(14,14)
        call set_status_hint(view%reset_probe,'Restore default selection coordinates; preserve the view')
        call view%reset_probe%setAccessibleName('Reset selection coordinates')
        call probe_heading%addWidget(view%reset_probe)
        call build_selected_control(section,view)
        probe_label=QLabel(''); call probe_label%setAlignment(ior(Qt_AlignLeft,Qt_AlignTop))
        call probe_label%setObjectName('studioProbeValues')
        call probe_label%setAccessibleName('Probe point and first derivatives')
        call probe_label%setWordWrap(.true.)
        call probe_label%setContentsMargins(10,0,10,7)
        call section%addWidget(probe_label)
        call add_section(content,section,view,2,.false.)
        diagnostic_label=QLabel(''); call diagnostic_label%setWordWrap(.true.)
        call diagnostic_label%setObjectName('studioPatchDiagnostics')
        call diagnostic_label%setTextFormat(Qt_PlainText)
        call diagnostic_label%setContentsMargins(10,4,10,7)
        call section%addWidget(diagnostic_label)
        associate(visual=>command_visual(diagnose_with_ai_command,'Ask Assistant…',29))
            view%repair_diagnostic=QPushButton(visual%caption)
            call view%repair_diagnostic%setIcon(line_icon(visual%glyph))
        end associate
        call view%repair_diagnostic%setObjectName('studioRepairDiagnostic')
        call set_status_hint(view%repair_diagnostic,'Add the current diagnostic to an Assistant draft for review.')
        call view%repair_diagnostic%hide()
        call section%addWidget(view%repair_diagnostic,alignment=Qt_AlignHCenter)

        call content%addStretch()
        scroll=QScrollArea(); call scroll%setWidgetResizable(.true.)
        call scroll%setObjectName('studioInspectorScroll')
        call scroll%setAccessibleName('Inspector details')
        call scroll%setFrameShape(QFrame_NoFrame)
        call scroll%setHorizontalScrollBarPolicy(Qt_ScrollBarAlwaysOff)
        call scroll%setWidget(page)
        call layout%addWidget(scroll,1)
    end subroutine build_inspector_widgets

    subroutine build_probe_axis(layout,axis,column,coordinate,label,field)
        type(QGridLayout),intent(inout)::layout
        character(1),intent(in)::axis
        integer,intent(in)::column
        real(dp),intent(in)::coordinate
        type(QLabel),intent(out)::label
        type(QDoubleSpinBox),intent(out)::field

        label=QLabel(axis)
        call label%setProperty('probeAxis',axis)
        call layout%addWidget(label,0,column)
        field=QDoubleSpinBox()
        call field%setRange(0._dp,1._dp)
        call field%setDecimals(4)
        call field%setSingleStep(0.01_dp)
        call field%setValue(coordinate)
        call set_status_hint(field,'Set '//axis//' relative to the patch interval.')
        call field%setAccessibleName('Normalized probe '//axis//' coordinate')
        call field%setObjectName('probe'//axis)
        call field%setFixedSize(74,28)
        call layout%addWidget(field,0,column+1)
    end subroutine build_probe_axis

    subroutine build_selected_control(layout,view)
        type(QVBoxLayout),intent(inout)::layout
        type(studio_inspector_t),intent(inout)::view
        type(QHBoxLayout)::row

        view%selected_control_row=QWidget()
        call view%selected_control_row%setObjectName('studioSelectedControl')
        row=QHBoxLayout(view%selected_control_row)
        call row%setContentsMargins(10,4,10,5)
        view%selected_control_label=QLabel('')
        call view%selected_control_label%setObjectName('studioSelectedControlReadout')
        call view%selected_control_label%setTextFormat(Qt_PlainText)
        call view%selected_control_label%setWordWrap(.true.)
        call row%addWidget(view%selected_control_label,1)
        view%edit_control=QPushButton('Edit…')
        call view%edit_control%setObjectName('studioEditSelectedControl')
        call view%edit_control%setAccessibleName('Edit selected control point coordinates and weight')
        call view%edit_control%setMinimumSize(68,27)
        call view%edit_control%setIcon(line_icon(68))
        call view%edit_control%setIconSize(15,15)
        call row%addWidget(view%edit_control)
        call layout%addWidget(view%selected_control_row)
        call view%selected_control_row%hide()
    end subroutine build_selected_control

    subroutine add_section(parent,layout,view,index,expanded)
        type(QVBoxLayout),intent(inout)::parent
        type(QVBoxLayout),intent(out)::layout
        type(studio_inspector_t),intent(inout)::view
        integer,intent(in)::index
        logical,intent(in)::expanded
        type(QWidget)::container,body
        type(QVBoxLayout)::container_layout
        container=QWidget(); call container%setObjectName('inspectorSection')
        view%sections(index)=container
        container_layout=QVBoxLayout(container)
        call container_layout%setContentsMargins(0,0,0,0); call container_layout%setSpacing(0)
        view%headers(index)=QPushButton(trim(section_titles(index)))
        call view%headers(index)%setIconSize(14,14)
        call view%headers(index)%setObjectName('sectionToggle')
        call view%headers(index)%setCheckable(.true.)
        call view%headers(index)%setChecked(expanded)
        call view%headers(index)%setFocusPolicy(Qt_StrongFocus)
        call view%headers(index)%setAccessibleName(trim(section_titles(index)))
        call container_layout%addWidget(view%headers(index))
        body=QWidget(); view%bodies(index)=body
        layout=QVBoxLayout(body)
        call layout%setContentsMargins(0,0,0,0); call layout%setSpacing(0)
        call container_layout%addWidget(body); call parent%addWidget(container)
        call update_inspector_section(view,index)
    end subroutine add_section

    subroutine update_inspector_section(view,index)
        type(studio_inspector_t),intent(inout)::view
        integer,intent(in)::index
        logical::expanded
        expanded=view%headers(index)%isChecked()
        if(expanded)then
            call view%headers(index)%setIcon(line_icon(icon_chevron))
            call set_status_hint(view%headers(index),'Collapse '//trim(section_titles(index))//' · Space')
            call view%headers(index)%setAccessibleDescription('Expanded section. Activate to collapse.')
        else
            call view%headers(index)%setIcon(line_icon(icon_chevron_right))
            call set_status_hint(view%headers(index),'Expand '//trim(section_titles(index))//' · Space')
            call view%headers(index)%setAccessibleDescription('Collapsed section. Activate to expand.')
        end if
        call view%bodies(index)%setVisible(expanded)
    end subroutine update_inspector_section
end module studio_inspector
