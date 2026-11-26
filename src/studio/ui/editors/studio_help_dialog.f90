module studio_help_dialog
    use qt, only: QWidget, QDialog, QVBoxLayout, QFormLayout, QLabel, QFontMetrics, QTabWidget, QDialogButtonBox, QDialogButtonBox_Close
    implicit none(type, external)
    private
    public :: show_interaction_guide

contains

    subroutine show_interaction_guide(parent)
        class(QWidget), intent(in) :: parent
        type(QDialog), target :: dialog
        type(QVBoxLayout) :: layout
        type(QLabel) :: heading
        type(QTabWidget) :: tabs
        type(QDialogButtonBox) :: buttons
        integer :: result

        dialog=QDialog(parent)
        call dialog%setObjectName('studioInteractionGuide')
        call dialog%setWindowTitle('Interaction guide')
        call dialog%setMinimumWidth(494)
        call dialog%resize(556,372)
        layout=QVBoxLayout(dialog)
        call layout%setContentsMargins(16,14,16,14)
        call layout%setSpacing(8)
        heading=QLabel('Interaction guide',dialog)
        call heading%setObjectName('studioDialogTitle')
        call layout%addWidget(heading)
        tabs=QTabWidget(dialog)
        call tabs%setObjectName('studioInteractionGuideTabs')
        call tabs%setAccessibleName('Interaction guide topics')
        call add_guide_topic(tabs,'Camera', &
            [character(16)::'Orbit','Pan','Zoom','Fit scene','Fit selection','Isometric'], &
            [character(90)::'Shift+middle-drag or right-drag in the 3D viewer.','Middle-drag in the 3D viewer.', &
                'Scroll to zoom the viewer under the pointer.','Press F to frame all visible geometry.', &
                'Press Shift+F to frame your selection.','Press R to restore the isometric view.'])
        call add_guide_topic(tabs,'Selection', &
            [character(16)::'Pick','Box or Lasso','Combine','Delete','Visibility','Recall'], &
            [character(100)::'Click an item in the viewer or Scene.','Use Model → Box or Lasso, then drag in the viewer.', &
                'Shift adds, Ctrl toggles, and Esc cancels an area selection.', &
                'Press Del with the viewer or Scene focused to delete the selected objects.', &
                'Use the eye next to an item in Scene to show or hide it.', &
                'Use Ctrl+1, Ctrl+2, or Ctrl+3 to recall the last curve, surface, or volume.'])
        call add_guide_topic(tabs,'Patch editing', &
            [character(17)::'Points + weights','Parameter Space','Knots','Set domain','Undo / Redo'], &
            [character(112)::'Use Control points to drag. Ctrl+E opens exact values.', &
                'Scroll to zoom. Middle-drag to pan. Shift+middle-drag or right-drag orbits volumes.', &
                'In Knots mode, drag handles. Use + or − to insert or remove.', &
                'Use Patch → Set domain to edit parameter bounds.', &
                'Ctrl+Z undoes. Ctrl+Shift+Z or Ctrl+Y redoes.'])
        call add_guide_topic(tabs,'Workspace', &
            [character(16)::'Find','Fullscreen','Panels','Assistant','Project tabs','Save As','Screenshot'], &
            [character(108)::'Ctrl+K opens Find command. Ctrl+F searches Scene.', &
                'Press F11 to enter or exit fullscreen. Ctrl+Shift+F toggles Canvas Only.', &
                'Press Shift+F11 to toggle the parameter panels. Drag a panel divider to resize it.', &
                'Press Ctrl+Shift+A to show or hide the Assistant.', &
                'Ctrl+Tab switches to the next project. Ctrl+Shift+Tab switches to the previous one.', &
                'Ctrl+Shift+S saves the current project under a new name.', &
                'Ctrl+Alt+S or Project → Capture saves a workspace image.'])
        call layout%addWidget(tabs,1)
        buttons=QDialogButtonBox(QDialogButtonBox_Close,dialog)
        call buttons%rejected%connect(close_interaction_guide,dialog)
        call layout%addWidget(buttons)
        result=dialog%exec()
        call dialog%delete()
    end subroutine show_interaction_guide

    subroutine add_guide_topic(tabs,title,actions,instructions)
        type(QTabWidget), intent(inout) :: tabs
        character(*), intent(in) :: title,actions(:),instructions(:)
        type(QWidget) :: page
        type(QFormLayout) :: form
        type(QLabel) :: detail
        type(QFontMetrics) :: metrics
        integer :: i,row_height

        page=QWidget(tabs)
        form=QFormLayout(page)
        call form%setContentsMargins(12,10,12,10)
        call form%setSpacing(5)
        metrics=QFontMetrics(page%font())
        row_height=2*metrics%lineSpacing()+2
        call metrics%delete()
        do i=1,size(actions)
            detail=QLabel(trim(instructions(i)),page)
            call detail%setProperty('dialogRole','detail')
            call detail%setWordWrap(.true.)
            call detail%setMinimumHeight(row_height)
            call form%addRow(trim(actions(i)),detail)
        end do
        call tabs%addTab(page,title)
    end subroutine add_guide_topic

    subroutine close_interaction_guide(data)
        class(*), intent(inout), optional :: data
        select type(dialog=>data)
        type is(QDialog)
            call dialog%reject()
        end select
    end subroutine close_interaction_guide

end module studio_help_dialog
