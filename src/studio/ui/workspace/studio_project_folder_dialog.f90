module studio_project_folder_dialog
    use qt, only: QWidget, QFileDialog_getExistingDirectory, QDir_currentPath
    use studio_project_folder, only: canonical_project_folder
    implicit none(type, external)
    private
    public :: choose_project_folder
contains
    function choose_project_folder(window,start_path) result(folder)
        class(QWidget), intent(in) :: window
        character(*), optional, intent(in) :: start_path
        character(:), allocatable :: folder,start

        start=QDir_currentPath()
        if(present(start_path))start=start_path
        folder=canonical_project_folder(QFileDialog_getExistingDirectory('Choose working folder',start,parent=window))
    end function choose_project_folder
end module studio_project_folder_dialog
