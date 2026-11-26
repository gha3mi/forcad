module studio_project_widgets
    use studio_status_hints, only: set_status_hint
    use studio_theme, only: studio_style_primary
    use qt, only: QWidget, QLabel, QVBoxLayout, QHBoxLayout, QFormLayout, QPushButton, QIcon, QDialog, QComboBox, &
        QDialogButtonBox, QDialogButtonBox_Cancel, QDialogButtonBox_AcceptRole, QDialog_Accepted, Qt_AlignHCenter
    use studio_state, only: studio_t
    use studio_icons, only: line_icon
    use studio_command_presentation, only: command_visual
    use studio_commands, only: new_project_command, open_project_command, resume_project_command, &
        sample_project_command, head_sample_command
    use studio_presentation_text, only: integer_text, studio_product_name
    use studio_units, only: unitless, unit_labels
    implicit none(type, external)
    private
    public :: build_project_home, edit_project_units

    abstract interface
        subroutine project_button_connector(button,s,command)
            import :: QPushButton, studio_t
            implicit none(type, external)
            type(QPushButton), intent(inout) :: button
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: command
        end subroutine project_button_connector
    end interface

contains

    subroutine build_project_home(s,connect_button)
        type(studio_t),intent(inout),target::s
        procedure(project_button_connector) :: connect_button
        type(QVBoxLayout)::outer,body,actions_layout,samples_layout
        type(QHBoxLayout)::columns
        type(QWidget)::content,card
        type(QLabel)::label

        s%home_page=QWidget()
        call s%home_page%setObjectName('studioHome')
        outer=QVBoxLayout(s%home_page)
        call outer%setContentsMargins(40,32,40,36)
        call outer%addStretch(1)
        content=QWidget(s%home_page)
        call content%setMaximumWidth(1080)
        body=QVBoxLayout(content)
        call body%setContentsMargins(0,0,0,0)
        call body%setSpacing(12)
        label=QLabel(studio_product_name)
        call label%setObjectName('studioHomeTitle')
        call body%addWidget(label)
        label=QLabel('Parametric geometry workspace')
        call label%setObjectName('studioHomeSubtitle')
        call body%addWidget(label)
        call body%addSpacing(18)
        columns=QHBoxLayout()
        call columns%setSpacing(20)
        call home_card(card,actions_layout,'Start')
        call home_button(actions_layout,s,connect_button,'New project',new_project_command,81,'Start with an empty viewport',.true.)
        call home_button(actions_layout,s,connect_button,'Open project…',open_project_command,82,'Open a ForCAD Studio project file')
        s%widgets%resume_button=home_button_widget(actions_layout,s,connect_button,'Resume project',resume_project_command,86, &
            'Return to the current workspace')
        call actions_layout%addStretch()
        call columns%addWidget(card,2)
        call home_card(card,samples_layout,'Samples')
        call home_sample_button(samples_layout,s,connect_button,'Chess',sample_project_command,85, &
            'Multipatch board and pieces with saved meshes')
        call home_sample_button(samples_layout,s,connect_button,'Head',head_sample_command,24, &
            'Portrait modeled with surface patches')
        call samples_layout%addStretch()
        call columns%addWidget(card,3)
        call body%addLayout(columns)
        call outer%addWidget(content,alignment=Qt_AlignHCenter)
        call outer%addStretch(2)
        call s%pages%addWidget(s%home_page)
    end subroutine build_project_home

    subroutine home_card(card,layout,caption)
        type(QWidget),intent(out)::card
        type(QVBoxLayout),intent(out)::layout
        character(*),intent(in)::caption
        type(QLabel)::label

        card=QWidget()
        call card%setObjectName('studioHomeCard')
        layout=QVBoxLayout(card)
        call layout%setContentsMargins(20,17,20,20)
        call layout%setSpacing(9)
        label=QLabel(caption)
        call label%setObjectName('studioHomeSection')
        call layout%addWidget(label)
        call layout%addSpacing(4)
    end subroutine home_card

    subroutine home_button(layout,s,connect_button,title,command,glyph,detail,primary)
        type(QVBoxLayout),intent(inout)::layout
        type(studio_t),intent(inout),target::s
        procedure(project_button_connector) :: connect_button
        character(*),intent(in)::title,detail
        integer,intent(in)::command,glyph
        logical,intent(in),optional::primary
        type(QPushButton)::button

        button=home_button_widget(layout,s,connect_button,title,command,glyph,detail)
        if(present(primary))call button%setProperty('primary',primary)
    end subroutine home_button

    subroutine home_sample_button(layout,s,connect_button,title,command,glyph,detail)
        type(QVBoxLayout),intent(inout)::layout
        type(studio_t),intent(inout),target::s
        procedure(project_button_connector) :: connect_button
        character(*),intent(in)::title,detail
        integer,intent(in)::command,glyph
        type(QPushButton)::button

        button=home_button_widget(layout,s,connect_button,title,command,glyph,detail)
        call button%setMinimumHeight(90)
        call button%setIconSize(46,46)
    end subroutine home_sample_button

    function home_button_widget(layout,s,connect_button,title,command,glyph,detail) result(button)
        type(QVBoxLayout),intent(inout)::layout
        type(studio_t),intent(inout),target::s
        procedure(project_button_connector) :: connect_button
        character(*),intent(in)::title,detail
        integer,intent(in)::command,glyph
        type(QPushButton)::button
        type(QIcon)::icon

        associate(visual=>command_visual(command,title,glyph))
            button=QPushButton(trim(visual%caption)//new_line('a')//detail)
            call button%setAccessibleName(trim(visual%caption))
            icon=line_icon(visual%glyph)
        end associate
        call button%setObjectName('studioHomeAction'//integer_text(command))
        call button%setProperty('welcomeAction',.true.)
        call set_status_hint(button,detail)
        call button%setIcon(icon)
        call icon%delete()
        call button%setIconSize(22,22)
        call button%setMinimumHeight(60)
        call connect_button(button,s,command)
        call layout%addWidget(button)
    end function home_button_widget

    logical function edit_project_units(parent,model_unit,display_unit) result(accepted)
        class(QWidget), intent(in) :: parent
        integer, intent(inout) :: model_unit,display_unit
        type(QDialog),target::dialog
        type(QVBoxLayout)::layout
        type(QFormLayout)::form
        type(QLabel)::header,note
        type(QComboBox)::model
        type(QComboBox),target::display
        type(QDialogButtonBox)::buttons
        type(QPushButton)::apply
        integer::i

        dialog=QDialog(parent)
        call dialog%setObjectName('studioUnitsDialog')
        call dialog%setWindowTitle('Project units')
        call dialog%setMinimumWidth(436)
        layout=QVBoxLayout(dialog)
        call layout%setContentsMargins(18,16,18,16)
        call layout%setSpacing(10)
        header=QLabel('Project units')
        call header%setObjectName('studioDialogTitle')
        call layout%addWidget(header)
        note=QLabel('Set the unit of the stored coordinates and choose a unit for measurements and numeric inputs.')
        call note%setProperty('dialogRole','detail')
        call note%setWordWrap(.true.)
        call layout%addWidget(note)
        form=QFormLayout()
        call form%setSpacing(10)
        model=QComboBox(dialog)
        display=QComboBox(dialog)
        call model%setObjectName('studioModelUnits')
        call display%setObjectName('studioDisplayUnits')
        call display%addItem('Same as model')
        do i=lbound(unit_labels,1),ubound(unit_labels,1)
            call model%addItem(trim(unit_labels(i)))
            if(i/=unitless)call display%addItem(trim(unit_labels(i)))
        end do
        call model%setCurrentIndex(model_unit)
        call display%setCurrentIndex(display_unit)
        call model_unit_changed(model_unit,display)
        call model%currentIndexChanged%connect(model_unit_changed,display)
        call model%setAccessibleName('Model unit')
        call display%setAccessibleName('Display unit')
        call set_status_hint(model,'The length unit attached to the geometry coordinates.')
        call set_status_hint(display,'The length unit shown in measurements and numeric inputs.')
        call form%addRow('Model unit',model)
        call form%addRow('Display unit',display)
        call layout%addLayout(form)
        buttons=QDialogButtonBox(QDialogButtonBox_Cancel,dialog)
        apply=buttons%addButton('Apply units',QDialogButtonBox_AcceptRole)
        call apply%setObjectName('studioUnitsApply')
        call studio_style_primary(apply)
        call apply%setDefault(.true.)
        call model%setFocus()
        call buttons%accepted%connect(accept_units,dialog)
        call buttons%rejected%connect(reject_units,dialog)
        call layout%addWidget(buttons)
        accepted=dialog%exec()==QDialog_Accepted
        if(accepted)then
            model_unit=model%currentIndex()
            display_unit=display%currentIndex()
            if(model_unit==unitless)display_unit=unitless
        end if
        call dialog%delete()
    end function edit_project_units

    subroutine model_unit_changed(index,data)
        integer,intent(in)::index
        class(*),intent(inout),optional::data
        select type(display=>data)
        type is(QComboBox)
            call display%setEnabled(index/=unitless)
            if(index==unitless)call display%setCurrentIndex(unitless)
        end select
    end subroutine model_unit_changed

    subroutine accept_units(data)
        class(*),intent(inout),optional::data
        select type(dialog=>data)
        type is(QDialog)
            call dialog%accept()
        end select
    end subroutine accept_units

    subroutine reject_units(data)
        class(*),intent(inout),optional::data
        select type(dialog=>data)
        type is(QDialog)
            call dialog%reject()
        end select
    end subroutine reject_units
end module studio_project_widgets
