module studio_assistant_team_settings
    use qt, only: QWidget, QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox, QComboBox, QLabel, QDialogButtonBox, QPushButton, &
        QDialogButtonBox_Cancel, QDialogButtonBox_AcceptRole, QDialog_Accepted, QSizePolicy_Ignored, QSizePolicy_Fixed, Qt_PlainText
    use studio_assistant_session, only: assistant_session_t, agent_model_entry
    use studio_assistant_transport, only: assistant_worker_t, assistant_worker_count
    use studio_presentation_text, only: integer_text
    use studio_theme, only: studio_style_primary
    implicit none(type, external)
    private
    public :: edit_assistant_team

    type :: worker_editor
        type(QComboBox) :: provider, model, effort
        type(assistant_worker_t) :: saved
        type(agent_model_entry), allocatable :: models(:)
        type(assistant_session_t), pointer :: session=>null()
        integer :: provider_count=0, first_model=0
        logical :: updating=.false.
    end type worker_editor

contains

    logical function edit_assistant_team(parent,session,providers,lead,roles,workers) result(accepted)
        class(QWidget),intent(in) :: parent
        type(assistant_session_t),intent(in),target :: session
        type(QComboBox),intent(in) :: providers
        character(*),intent(in) :: lead,roles(assistant_worker_count)
        type(assistant_worker_t),intent(inout) :: workers(assistant_worker_count)
        type(worker_editor),target :: editors(assistant_worker_count)
        type(QDialog),target :: dialog
        type(QVBoxLayout) :: layout
        type(QHBoxLayout) :: worker_row
        type(QLabel) :: header,summary
        type(QDialogButtonBox) :: buttons
        type(QPushButton) :: save,cancel
        integer :: i

        dialog=QDialog(parent)
        call dialog%setObjectName('studioAssistantTeamDialog')
        call dialog%setWindowTitle('Team settings')
        call dialog%setMinimumWidth(590)
        layout=QVBoxLayout(dialog)
        call layout%setContentsMargins(18,16,18,16)
        call layout%setSpacing(9)
        header=QLabel('Team settings')
        call header%setObjectName('studioDialogTitle')
        call layout%addWidget(header)
        summary=QLabel(lead)
        call summary%setTextFormat(Qt_PlainText)
        call summary%setProperty('dialogRole','context')
        call summary%setWordWrap(.true.)
        call layout%addWidget(summary)

        worker_row=QHBoxLayout()
        call worker_row%setSpacing(10)
        do i=1,size(workers)
            call build_worker_editor(worker_row,dialog,editors(i),session,providers,trim(roles(i)),i,workers(i))
        end do
        call layout%addLayout(worker_row)

        buttons=QDialogButtonBox(QDialogButtonBox_Cancel,dialog)
        cancel=buttons%button(QDialogButtonBox_Cancel)
        call cancel%setObjectName('studioAssistantTeamCancel')
        save=buttons%addButton('Save team',QDialogButtonBox_AcceptRole)
        call save%setObjectName('studioAssistantTeamSave')
        call studio_style_primary(save)
        call save%setDefault(.true.)
        call buttons%accepted%connect(accept_settings,dialog)
        call buttons%rejected%connect(cancel_settings,dialog)
        call layout%addWidget(buttons)
        accepted=dialog%exec()==QDialog_Accepted
        if(accepted)then
            do i=1,size(workers)
                workers(i)=current_worker(editors(i))
            end do
        end if
        call dialog%delete()
    end function edit_assistant_team

    subroutine build_worker_editor(row,parent,editor,session,providers,role,number,worker)
        type(QHBoxLayout),intent(inout) :: row
        type(QDialog),intent(in) :: parent
        type(worker_editor),intent(inout),target :: editor
        type(assistant_session_t),intent(in),target :: session
        type(QComboBox),intent(in) :: providers
        character(*),intent(in) :: role
        integer,intent(in) :: number
        type(assistant_worker_t),intent(in) :: worker
        type(QGroupBox) :: card
        type(QFormLayout) :: form
        character(:),allocatable :: suffix
        integer :: i,selected

        card=QGroupBox(role,parent)
        form=QFormLayout(card)
        call form%setContentsMargins(9,10,9,8)
        call form%setSpacing(6)
        suffix=integer_text(number)
        editor%saved=worker
        editor%session=>session
        editor%provider_count=providers%count()
        editor%provider=QComboBox()
        editor%model=QComboBox()
        editor%effort=QComboBox()
        call editor%provider%setObjectName('studioAssistantTeamProvider'//suffix)
        call editor%model%setObjectName('studioAssistantTeamModel'//suffix)
        call editor%effort%setObjectName('studioAssistantTeamEffort'//suffix)
        call editor%provider%setAccessibleName(role//' agent')
        call editor%model%setAccessibleName(role//' model')
        call editor%effort%setAccessibleName(role//' effort')
        call editor%provider%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Fixed)
        call editor%model%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Fixed)
        call editor%effort%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Fixed)
        call editor%provider%addItem('Same as lead')
        do i=0,providers%count()-1
            call editor%provider%addItem(providers%itemText(i))
        end do
        selected=0
        if(len(worker%provider)>0)then
            selected=editor%provider%findText(worker%provider)
            if(selected<0)then
                call editor%provider%addItem(worker%provider//' (unavailable)')
                selected=editor%provider%count()-1
            end if
        end if
        call editor%provider%setCurrentIndex(selected)
        call worker_provider_changed(selected,editor)
        call editor%provider%currentIndexChanged%connect(worker_provider_changed,editor)
        call editor%model%currentIndexChanged%connect(worker_model_changed,editor)
        call form%addRow('Agent',editor%provider)
        call form%addRow('Model',editor%model)
        call form%addRow('Effort',editor%effort)
        call row%addWidget(card,1)
    end subroutine build_worker_editor

    subroutine worker_provider_changed(index,data)
        integer,intent(in) :: index
        class(*),intent(inout),optional :: data
        character(:),allocatable :: provider,model
        integer :: i,selected

        select type(editor=>data)
        type is(worker_editor)
            editor%updating=.true.
            provider=worker_provider(editor)
            call editor%session%read_agent_models(provider,editor%models)
            editor%first_model=0
            if(len(editor%models(0)%id)==0)editor%first_model=1
            call editor%model%clear()
            if(index>0)then
                call editor%model%addItem('Provider default')
            else
                call editor%model%addItem('Same as lead')
            end if
            selected=0
            model=''
            if(index>0.and.provider==editor%saved%provider)model=editor%saved%model
            do i=editor%first_model,ubound(editor%models,1)
                call editor%model%addItem(editor%models(i)%label)
                if(editor%models(i)%id==model)selected=editor%model%count()-1
            end do
            if(selected==0.and.len(model)>0)then
                call editor%model%addItem(model//' (unavailable)')
                selected=editor%model%count()-1
            end if
            call editor%model%setCurrentIndex(selected)
            call editor%model%setEnabled(index>0)
            editor%updating=.false.
            call worker_model_changed(selected,editor)
        end select
    end subroutine worker_provider_changed

    subroutine worker_model_changed(index,data)
        integer,intent(in) :: index
        class(*),intent(inout),optional :: data
        character(:),allocatable :: effort
        integer :: i,model,selected

        select type(editor=>data)
        type is(worker_editor)
            if(editor%updating)return
            call editor%effort%clear()
            if(editor%provider%currentIndex()>0)then
                call editor%effort%addItem('Default')
            else
                call editor%effort%addItem('Same as lead')
            end if
            selected=0
            effort=''
            if(editor%provider%currentIndex()>0.and.worker_provider(editor)==editor%saved%provider.and. &
                worker_model(editor)==editor%saved%model)effort=editor%saved%effort
            if(index>0)then
                model=index+editor%first_model-1
                if(model<=ubound(editor%models,1))then
                    associate(options=>editor%models(model)%efforts)
                        do i=1,size(options)
                            call editor%effort%addItem(options(i)%label)
                            if(options(i)%id==effort)selected=i
                        end do
                    end associate
                end if
            end if
            if(selected==0.and.len(effort)>0)then
                call editor%effort%addItem(effort//' (unavailable)')
                selected=editor%effort%count()-1
            end if
            call editor%effort%setCurrentIndex(selected)
            call editor%effort%setEnabled(editor%effort%count()>1)
        end select
    end subroutine worker_model_changed

    function worker_provider(editor) result(provider)
        type(worker_editor),intent(in) :: editor
        character(:),allocatable :: provider
        integer :: index

        provider=''
        index=editor%provider%currentIndex()
        if(index<1)return
        provider=editor%provider%currentText()
        if(index>editor%provider_count)provider=editor%saved%provider
    end function worker_provider

    function worker_model(editor) result(model)
        type(worker_editor),intent(in) :: editor
        character(:),allocatable :: model
        integer :: index

        model=''
        index=editor%model%currentIndex()
        if(index<1)return
        index=index+editor%first_model-1
        if(index<=ubound(editor%models,1))then
            model=editor%models(index)%id
        else
            model=editor%saved%model
        end if
    end function worker_model

    function current_worker(editor) result(worker)
        type(worker_editor),intent(in) :: editor
        type(assistant_worker_t) :: worker
        integer :: effort,model

        worker%provider=worker_provider(editor)
        worker%model=worker_model(editor)
        worker%effort=''
        effort=editor%effort%currentIndex()
        if(effort<1)return
        model=editor%model%currentIndex()+editor%first_model-1
        if(editor%model%currentIndex()>0.and.model<=ubound(editor%models,1))then
            if(effort<=size(editor%models(model)%efforts))then
                worker%effort=editor%models(model)%efforts(effort)%id
                return
            end if
        end if
        worker%effort=editor%saved%effort
    end function current_worker

    subroutine accept_settings(data)
        class(*),intent(inout),optional :: data
        select type(dialog=>data)
        type is(QDialog)
            call dialog%accept()
        end select
    end subroutine accept_settings

    subroutine cancel_settings(data)
        class(*),intent(inout),optional :: data
        select type(dialog=>data)
        type is(QDialog)
            call dialog%reject()
        end select
    end subroutine cancel_settings

end module studio_assistant_team_settings
