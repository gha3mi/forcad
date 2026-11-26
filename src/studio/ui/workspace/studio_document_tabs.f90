module studio_document_tabs
    use, intrinsic :: iso_c_binding, only: dp=>c_double
    use, intrinsic :: iso_fortran_env, only: int64
    use qt, only: QAction,QShortcut,QTabBar,qt_string
    use studio_status_bar, only: studio_status_log
    use studio_document, only: studio_document_t
    use studio_selection, only: selection_set_t,selection_target_t
    implicit none(type,external)
    private
    public :: studio_document_tabs_t,studio_document_tab_t,studio_document_view_t

    type :: studio_document_view_t
        type(selection_set_t) :: selection,isolation
        type(qt_string),allocatable :: isolated_meshes(:)
        real(dp) :: camera(5)=0._dp,domain(5)=0._dp,probe(3)=0._dp,section_fraction=.5_dp
        integer :: ribbon=0,rank=2,axis=1,mesh_display(4)=0,feature(4)=0
        logical :: mesh_only=.false.,mesh_pick=.false.,mesh_shaded=.true.,physical_knots=.false.,derivative=.false.
    end type studio_document_view_t

    type :: studio_document_tab_t
        type(studio_document_t),allocatable :: document
        type(studio_document_view_t) :: view
        type(studio_status_log) :: status_log
    end type studio_document_tab_t

    type :: studio_document_tabs_t
        type(QTabBar) :: bar
        type(QAction) :: copy,paste,close
        type(QShortcut) :: next_shortcut,previous_shortcut
        type(studio_document_tab_t),allocatable :: items(:)
        type(studio_document_t),allocatable :: clipboard
        type(selection_target_t),allocatable :: clipboard_selection(:)
        integer(int64) :: revision=0_int64
        integer :: active=0
        logical :: changing=.false.
    end type studio_document_tabs_t
end module studio_document_tabs
