! Shared vector artwork for the native Qt workspace.
module studio_icons
    use, intrinsic :: iso_c_binding, only: dp=>c_double, c_int, c_funptr, c_funloc, c_associated
    use qt, only: QBrush, QColor, QIcon, QImage, QLabel, QLinearGradient, QPainter, QPainter_Antialiasing, &
        QPainterPath, QPen, QPixmap, QPixmap_fromImage, QSvgGenerator, QTemporaryDir, &
        Qt_DashLine, Qt_NoBrush, Qt_NoPen, Qt_transparent
    use studio_presentation_text, only: integer_text
    use studio_storage_paths, only: studio_work_directory
    implicit none(type, external)
    private

    public :: line_icon, stylesheet_icon
    integer, parameter, public :: icon_copy=11, icon_screenshot=26, icon_settings=33, icon_check=89
    integer, parameter, public :: icon_close=94, icon_minimize=95, icon_maximize=96, icon_restore=97
    integer, parameter, public :: icon_more=98, icon_send=99, icon_stop=100, icon_history=101, icon_new_chat=102
    integer, parameter, public :: icon_dock=103, icon_float=104, icon_chevron=105, icon_chevron_right=106

    integer, parameter, public :: icon_fit_selection=107, icon_fit_panels=108, icon_refresh=109, icon_export=110
    integer, parameter, public :: icon_mesh=111, icon_delete=112, icon_measure=113
    integer, parameter, public :: icon_visibility_mixed=114, icon_scene_reference=115
    integer, parameter, public :: icon_box_selection=116, icon_lasso_selection=117
    integer, parameter, public :: icon_mesh_size=118, icon_mesh_pick=119, icon_mesh_shaded=120
    integer, parameter, public :: icon_degree_elevate=121, icon_degree_reduce=122, icon_extract_section=123, icon_extract_face=124
    integer, parameter, public :: icon_camera_front=125, icon_camera_back=126, icon_camera_left=127, icon_camera_right=128
    integer, parameter, public :: icon_camera_top=129, icon_camera_bottom=130
    integer, parameter, public :: icon_bounding_face=131, icon_bounding_edge=132, icon_bounding_vertex=133, icon_generate_mesh=134
    integer, parameter, public :: icon_scene_filter=135, icon_chevron_up=136, icon_fast_shaded=137
    integer, parameter, public :: icon_fullscreen=138
    integer, parameter, public :: icon_face_grid=141, icon_isolate=142, icon_paste=143
    integer, parameter, public :: viewport_style_icons(*)=[8,icon_bounding_edge,21,icon_fast_shaded]
    integer, parameter :: face_top=1, face_left=2, face_right=3

    character(*), parameter :: paper='#e1edf9', quiet='#9caec1', dark='#152232'
    character(*), parameter :: aqua='#70ead3', blue='#8dbdff', amber='#ffd094', violet='#d9c4ff'
    character(*), parameter :: success='#a6eac0', danger='#f4a8a7'
    integer, parameter, public :: icon_scene_branch_closed=139, icon_scene_branch_open=140
    type(QTemporaryDir) :: stylesheet_icons

    interface
        integer(c_int) function register_icon_cleanup(callback) bind(C,name='atexit')
            import :: c_int,c_funptr
            implicit none(type, external)
            type(c_funptr), value :: callback
        end function register_icon_cleanup
    end interface

contains

    function line_icon(kind,label) result(icon)
        integer, intent(in) :: kind
        type(QLabel), intent(inout), optional :: label
        type(QIcon) :: icon
        type(QImage) :: image
        type(QPixmap) :: pixmap
        type(QPainter) :: p

        ! Render the 32-unit artwork with oversampling for small icons and high-DPI displays.
        image=QImage(192,192)
        call image%fill(Qt_transparent)
        p=QPainter(image)
        call p%scale(6._dp,6._dp)
        call paint_icon(p,kind)
        call p%end()
        pixmap=QPixmap_fromImage(image)
        icon=QIcon(pixmap)
        if(present(label))call label%setPixmap(pixmap)
        call pixmap%delete()
        call image%delete()
    end function line_icon

    function stylesheet_icon(kind,width,height,property) result(rule)
        integer, intent(in) :: kind,width,height
        character(*), intent(in) :: property
        character(:), allocatable :: rule,path,folder
        type(QSvgGenerator) :: image
        type(QPainter) :: painter
        integer :: length
        logical :: exists

        rule=''
        if(.not.c_associated(stylesheet_icons%ptr))then
            folder=studio_work_directory()
            if(len(folder)==0)return
            stylesheet_icons=QTemporaryDir(folder//'/icons.XXXXXX')
            if(.not.stylesheet_icons%isValid())then
                call stylesheet_icons%delete()
                return
            end if
            if(register_icon_cleanup(c_funloc(remove_stylesheet_icons))/=0_c_int)then
                call stylesheet_icons%delete()
                return
            end if
        end if
        path=stylesheet_icons%filePath(integer_text(kind)//'-'//integer_text(width)//'x'//integer_text(height)//'.svg')
        inquire(file=path,exist=exists,size=length)
        if(.not.exists.or.length<=0)then
            image=QSvgGenerator(path,width,height)
            select case(kind)
            case(icon_chevron,icon_chevron_up)
                call image%setViewBox(5._dp,9._dp,22._dp,14._dp)
            case(icon_chevron_right)
                call image%setViewBox(4._dp,4._dp,24._dp,24._dp)
            case default
                call image%setViewBox(0._dp,0._dp,32._dp,32._dp)
            end select
            painter=QPainter(image)
            call paint_icon(painter,kind)
            call painter%end()
            call image%delete()
            inquire(file=path,exist=exists,size=length)
        end if
        if(exists.and.length>0)rule=property//':url("'//path//'");'
    end function stylesheet_icon

    subroutine remove_stylesheet_icons() bind(C)
        call stylesheet_icons%delete()
    end subroutine remove_stylesheet_icons

    subroutine paint_icon(p,kind)
        type(QPainter), intent(inout) :: p
        integer, intent(in) :: kind
        type(QPainterPath) :: path
        type(QLinearGradient) :: gradient
        integer :: i

        call p%setRenderHint(QPainter_Antialiasing)
        call ink(p,paper)
        select case(kind)
        case(1)
            path=QPainterPath()
            call path%moveTo(7,6); call path%quadTo(7,4,9,4); call path%lineTo(19,4)
            call path%lineTo(25,10); call path%lineTo(25,26); call path%quadTo(25,28,23,28)
            call path%lineTo(9,28); call path%quadTo(7,28,7,26); call path%closeSubpath()
            call p%drawPath(path); call path%delete()
            call ink(p,aqua,1.9_dp)
            call p%drawPolyline([18,18,24],[5,11,11])
            call p%drawLine(12,19,20,19)
            call p%drawLine(12,23,17,23)
        case(2,82)
            call ink(p,aqua)
            path=QPainterPath()
            call path%moveTo(4,24); call path%lineTo(4,9); call path%quadTo(4,7,6,7)
            call path%lineTo(13,7); call path%lineTo(16,10); call path%lineTo(26,10)
            call path%quadTo(28,10,28,12)
            if(kind==2)then
                call path%lineTo(28,24); call path%quadTo(28,26,26,26)
                call path%lineTo(6,26); call path%quadTo(4,26,4,24)
                call p%setBrush(QColor(86,187,176,46))
                call p%drawPath(path)
                call ink(p,aqua,1.5_dp)
                call p%drawLine(5,14,27,14)
            else
                call p%drawPath(path)
                call p%setBrush(QColor(86,187,176,55))
                call p%drawPolygon([4,8,29,25],[26,15,15,26])
            end if
            call path%delete()
        case(3)
            call ink(p,aqua)
            call p%setBrush(QColor(86,187,176,42))
            call p%drawPolygon([5,22,27,27,5],[5,5,10,27,27])
            call ink(p,paper,1.8_dp)
            call p%drawPolyline([10,10,20,20],[6,12,12,6])
            call p%drawPolyline([10,10,22,22],[26,17,17,26])
            call p%drawLine(14,21,18,21)
        case(4,5)
            if(kind==5)then
                call p%translate(32,0); call p%scale(-1._dp,1._dp)
            end if
            path=QPainterPath()
            call path%moveTo(5,13); call path%lineTo(19,13)
            call path%cubicTo(24,13,27,16,27,20)
            call path%cubicTo(27,25,23,27,19,27); call path%lineTo(15,27)
            call p%drawPath(path); call path%delete()
            call ink(p,aqua)
            call p%drawPolyline([11,5,11],[7,13,19])
        case(6,7,27)
            call p%drawEllipse(4,4,17,17)
            call p%drawLine(19,19,28,28)
            if(kind/=27)then
                call ink(p,blue)
                call p%drawLine(8,12,17,12)
                if(kind==7)call p%drawLine(12,8,12,17)
            end if
        case(8,21,22,74,25,icon_camera_front:icon_bounding_vertex,icon_extract_section)
            call cube(p,kind)
        case(icon_fast_shaded)
            call p%save()
            call p%translate(8,1)
            call p%scale(.76_dp,.76_dp)
            call cube(p,8)
            call p%restore()
            call ink(p,aqua)
            call p%drawLine(3,12,11,12)
            call p%drawLine(1,17,9,17)
            call p%drawLine(5,22,11,22)
        case(icon_extract_face)
            call p%save()
            call p%translate(-1,5)
            call p%scale(.78_dp,.78_dp)
            call cube_outline(p,quiet)
            call p%restore()
            call ink(p,aqua)
            call p%setBrush(QColor(73,198,181,139))
            call p%drawPolygon([20,29,29,20],[13,8,23,28])
        case(9)
            call ink(p,amber)
            call p%drawLine(16,4,16,28)
            call p%drawLine(4,16,28,16)
            call p%drawPolyline([12,16,20],[8,4,8])
            call p%drawPolyline([12,16,20],[24,28,24])
            call p%drawPolyline([8,4,8],[12,16,20])
            call p%drawPolyline([24,28,24],[12,16,20])
            call solid(p,paper)
            call p%drawRect(14,14,4,4)
        case(10)
            call ink(p,violet)
            call p%drawRect(4,5,24,22)
            call p%drawLine(14,6,14,26)
            call p%drawLine(15,15,27,15)
            call solid(p,violet)
            call p%drawRect(7,9,4,14)
        case(11)
            call p%drawRect(4,5,15,15)
            call ink(p,blue)
            call p%setBrush(QColor(24,38,50))
            call p%drawRect(13,12,15,15)
        case(12,64,67)
            call ink(p,blue)
            call p%drawRect(5,5,22,22)
            call p%drawLine(5,16,27,16)
            select case(kind)
            case(12)
                call p%drawLine(5,10,27,10)
                call p%drawLine(5,22,27,22)
                call p%drawLine(10,5,10,27)
                call p%drawLine(16,5,16,27)
                call p%drawLine(22,5,22,27)
            case(64)
                call p%drawLine(16,5,16,27)
                call solid(p,blue)
                call p%drawEllipse(12,12,8,8)
            case(67)
                call p%drawLine(13,5,13,27)
            case default
                continue
            end select
        case(13)
            call ink(p,blue)
            call p%drawRect(4,5,24,22)
            call p%drawLine(4,11,28,11)
            call p%drawLine(14,17,23,17)
            call p%drawLine(14,22,21,22)
            call solid(p,blue)
            call p%drawRect(7,15,4,9)
        case(14)
            call p%drawEllipse(8,8,16,16)
            do i=0,7
                call p%save()
                call p%translate(16,16)
                call p%rotate(real(i*45,dp))
                call p%drawLine(0,-10,0,-13)
                call p%restore()
            end do
            call solid(p,blue)
            call p%drawEllipse(13,13,6,6)
        case(15)
            call solid(p,aqua)
            call p%drawPolygon([10,25,10],[6,16,26])
        case(16)
            call p%drawEllipse(4,4,24,24)
            path=QPainterPath()
            call path%moveTo(12,12); call path%quadTo(12,8,16,8)
            call path%quadTo(21,8,21,12); call path%quadTo(21,15,17,17)
            call path%quadTo(16,18,16,20)
            call p%drawPath(path); call path%delete()
            call solid(p,blue)
            call p%drawEllipse(15,22,2,2)
        case(17)
            gradient=QLinearGradient(2._dp,2._dp,30._dp,30._dp)
            call gradient%setColorAt(0._dp,QColor('#246c79'))
            call gradient%setColorAt(1._dp,QColor('#394b83'))
            call p%setBrush(gradient)
            call p%setPen(QPen(QColor('#8ee8e1'),1.3_dp))
            call p%drawRect(2,2,28,28)
            call solid(p,'#f1fcff')
            call p%drawPolygon([9,23,23,13,13,21,21,13,13,9],[9,9,13,13,16,16,20,20,24,24])
            call solid(p,aqua)
            call p%drawRect(22,22,3,3)
            call gradient%delete()
        case(18)
            call solid(p,paper)
            call p%drawPolygon([8,8,14,19,23,18,26],[4,25,20,28,26,18,17])
        case(icon_box_selection)
            call p%setPen(QPen(QColor(blue),2._dp,style=Qt_DashLine))
            call p%drawRect(5,6,22,20)
            call solid(p,aqua)
            call p%drawRect(3,4,5,5)
            call p%drawRect(24,23,5,5)
        case(icon_lasso_selection)
            call ink(p,blue)
            path=QPainterPath()
            call path%moveTo(9,22)
            call path%cubicTo(3,20,2,12,8,8)
            call path%cubicTo(15,3,28,5,28,13)
            call path%cubicTo(28,21,15,25,9,22)
            call path%cubicTo(6,21,6,24,8,26)
            call path%quadTo(10,28,15,28)
            call p%drawPath(path); call path%delete()
        case(19,icon_fit_selection,icon_fit_panels)
            call p%drawPolyline([5,5,12],[12,5,5])
            call p%drawPolyline([20,27,27],[5,5,12])
            call p%drawPolyline([5,5,12],[20,27,27])
            call p%drawPolyline([20,27,27],[27,27,20])
            call ink(p,blue,1.65_dp)
            select case(kind)
            case(icon_fit_selection)
                call p%drawEllipse(12,12,8,8)
            case(icon_fit_panels)
                call p%drawLine(12,10,12,22)
                call p%drawLine(21,10,21,22)
                call p%drawLine(13,16,20,16)
            case default
                call p%drawLine(16,12,16,20)
                call p%drawLine(12,16,20,16)
            end select
        case(23)
            path=QPainterPath()
            call path%moveTo(5,23)
            call path%cubicTo(11,3,13,5,18,16)
            call path%cubicTo(22,25,24,24,27,9)
            call ink(p,aqua)
            call p%drawPath(path); call path%delete()
            call solid(p,'#c3f7ee')
            call p%drawEllipse(2,20,6,6)
            call p%drawEllipse(24,6,6,6)
        case(24)
            path=QPainterPath()
            call path%moveTo(4,11); call path%quadTo(12,12,22,5)
            call path%quadTo(24,15,28,21); call path%quadTo(20,28,10,28)
            call path%quadTo(7,21,4,11); call path%closeSubpath()
            call ink(p,blue)
            call p%setBrush(QColor(78,146,212,80))
            call p%drawPath(path); call path%delete()
            call ink(p,blue,1.55_dp)
            call p%drawLine(12,11,17,27)
            call p%drawLine(18,8,23,24)
            path=QPainterPath()
            call path%moveTo(7,19); call path%quadTo(17,19,25,13)
            call p%drawPath(path); call path%delete()
        case(icon_refresh)
            call ink(p,blue)
            call p%drawArc(5,5,22,22,47._dp,290._dp)
            call p%drawPolyline([27,27,20],[4,11,11])
        case(26)
            call p%drawRect(4,10,24,16)
            call p%drawPolyline([9,12,20,23],[10,6,6,10])
            call ink(p,blue)
            call p%drawEllipse(12,13,8,8)
        case(28)
            call ink(p,blue,1.6_dp)
            call p%drawEllipse(8,8,16,16)
            call p%save()
            call p%translate(16,16); call p%rotate(-34._dp)
            call ink(p,aqua,2.1_dp)
            call p%drawArc(-14,-6,28,12,35._dp,282._dp)
            call p%drawPolyline([9,13,16],[1,4,0])
            call p%restore()
        case(29)
            call solid(p,violet)
            call spark(p,13,14,10,3)
            call solid(p,'#eee3ff')
            call spark(p,24,22,5,2)
            call solid(p,blue)
            call spark(p,24,7,4,1)
        case(30)
            call ink(p,aqua,2.7_dp)
            call p%drawLine(7,16,25,16)
            call p%drawLine(16,7,16,25)
        case(31,32,icon_visibility_mixed)
            if(kind==32)call ink(p,quiet)
            path=QPainterPath()
            call path%moveTo(3,16); call path%quadTo(16,1,29,16)
            call path%quadTo(16,31,3,16)
            call p%drawPath(path); call path%delete()
            if(kind==icon_visibility_mixed)then
                call ink(p,aqua)
                call p%drawLine(12,16,20,16)
            else
                call p%drawEllipse(12,12,8,8)
            end if
            if(kind==32)call p%drawLine(5,27,27,5)
        case(icon_scene_reference)
            call p%translate(16,16)
            call p%rotate(-39._dp)
            call ink(p,blue)
            call p%drawRoundedRect(-14,-4,14,8,4,4)
            call p%drawRoundedRect(0,-4,14,8,4,4)
            call ink(p,aqua)
            call p%drawLine(-5,0,5,0)
        case(33)
            call p%drawLine(5,8,27,8); call p%drawLine(5,16,27,16); call p%drawLine(5,24,27,24)
            call ink(p,blue,1.8_dp)
            call p%setBrush(QColor(dark))
            call p%drawEllipse(9,5,6,6); call p%drawEllipse(18,13,6,6); call p%drawEllipse(11,21,6,6)
        case(34)
            call ink(p,blue)
            call p%setBrush(QColor(67,125,180,61))
            call p%drawPolygon([3,13,21,11],[13,7,15,22])
            call ink(p,aqua)
            call p%setBrush(QColor(74,157,148,62))
            call p%drawPolygon([13,23,29,21],[7,8,20,15])
        case(35)
            call ink(p,success)
            call p%drawEllipse(4,4,24,24)
            call p%drawPolyline([10,14,23],[16,20,11])
        case(36)
            call ink(p,danger)
            call p%drawRect(5,6,22,20)
            call p%drawLine(10,16,22,16)
        case(37)
            path=QPainterPath()
            call path%moveTo(27,14); call path%cubicTo(27,7,21,4,16,4)
            call path%cubicTo(8,4,4,9,4,16); call path%cubicTo(4,23,8,28,15,28)
            call path%quadTo(20,28,20,24); call path%quadTo(17,19,24,19)
            call path%quadTo(27,19,27,14)
            call p%drawPath(path); call path%delete()
            call solid(p,aqua)
            call p%drawEllipse(8,13,5,5)
            call solid(p,amber)
            call p%drawEllipse(13,8,5,5)
            call solid(p,violet)
            call p%drawEllipse(20,10,5,5)
        case(icon_measure)
            call ink(p,aqua)
            call p%drawRect(3,10,26,12)
            do i=8,24,4
                call p%drawLine(i,11,i,14+2*mod(i/4,2))
            end do
        case(38)
            call ink(p,amber)
            call p%drawLine(10,22,27,22); call p%drawLine(10,22,10,5); call p%drawLine(10,22,23,9)
            call p%drawPolyline([6,10,14],[9,5,9])
            call p%drawPolyline([23,27,23],[18,22,26])
            call solid(p,amber)
            call p%drawEllipse(7,19,6,6)
        case(39)
            path=QPainterPath()
            call path%moveTo(7,7); call path%lineTo(7,18)
            call path%cubicTo(7,29,25,29,25,18); call path%lineTo(25,7)
            call path%lineTo(20,7); call path%lineTo(20,18)
            call path%cubicTo(20,22,12,22,12,18); call path%lineTo(12,7)
            call path%closeSubpath()
            call ink(p,aqua,1.7_dp)
            call p%setBrush(QColor(81,169,157,72))
            call p%drawPath(path); call path%delete()
            call ink(p,paper,1.65_dp)
            call p%drawLine(7,12,12,12); call p%drawLine(20,12,25,12)
        case(40,icon_dock)
            call p%drawRect(4,5,24,22)
            call p%drawLine(15,5,15,27)
            call ink(p,blue,1.85_dp)
            if(kind==40)then
                call p%drawLine(19,11,24,11); call p%drawLine(19,16,24,16); call p%drawLine(19,21,24,21)
            else
                call p%drawLine(18,16,25,16)
                call p%drawPolyline([22,18,22],[12,16,20])
            end if
        case(60)
            call ink(p,amber,2.2_dp)
            call p%drawArc(5,5,22,22,48._dp,286._dp)
            call p%drawPolyline([21,24,17],[3,8,9])
            call solid(p,paper)
            call p%drawEllipse(14,14,4,4)
        case(61)
            call ink(p,quiet,1.65_dp)
            call p%drawRect(5,14,13,13)
            call ink(p,amber)
            call p%drawPolyline([19,27,27],[5,5,13])
            call p%drawLine(13,19,26,6)
        case(62)
            call ink(p,amber)
            call p%setBrush(QColor(205,147,87,52))
            call p%drawPolygon([10,27,22,5],[6,6,26,26])
            call ink(p,amber,1.5_dp)
            call p%drawLine(18,6,13,26); call p%drawLine(8,16,25,16)
        case(63,65,66)
            call ink(p,blue)
            if(kind==63)then
                call p%drawLine(4,16,28,16)
                call p%drawLine(7,10,7,22); call p%drawLine(16,7,16,25); call p%drawLine(25,10,25,22)
            else
                call p%drawLine(4,23,28,23)
                call p%drawLine(7,18,7,28); call p%drawLine(25,18,25,28)
                call ink(p,merge(aqua,danger,kind==65))
                call p%drawLine(11,10,21,10)
                if(kind==65)call p%drawLine(16,5,16,15)
            end if
        case(68,88)
            call ink(p,blue)
            call p%drawPolygon([5,7,22,27,12],[26,19,4,9,24])
            call p%drawLine(19,7,24,12)
            call p%drawLine(7,19,12,24)
            call p%drawLine(16,27,27,27)
        case(69:73)
            call primitive(p,kind)
        case(75,76)
            call ink(p,blue)
            call p%setBrush(QColor(78,146,212,55))
            call p%drawEllipse(6,5,20,8)
            call ink(p,blue)
            call p%drawLine(6,9,6,23); call p%drawLine(26,9,26,23)
            call p%drawArc(6,19,20,8,180._dp,180._dp)
            if(kind==76)call p%drawEllipse(10,7,12,4)
        case(77)
            call ink(p,blue,1.8_dp)
            call p%drawPolygon([4,16,28,16],[20,15,20,25])
            call p%drawPolyline([4,4,16,28,28],[20,24,29,24,20])
            call ink(p,aqua)
            call p%drawLine(16,18,16,4)
            call p%drawPolyline([11,16,21],[9,4,9])
        case(78)
            call ink(p,blue,1.75_dp)
            call p%drawLine(16,4,16,28)
            call p%drawArc(4,12,24,11,49._dp,270._dp)
            call ink(p,aqua)
            call p%drawPolyline([22,27,27],[19,17,22])
        case(79)
            call ink(p,blue,1.8_dp)
            call p%drawEllipse(4,19,10,8)
            call ink(p,aqua)
            path=QPainterPath()
            call path%moveTo(10,24); call path%cubicTo(10,14,17,9,27,8)
            call p%drawPath(path); call path%delete()
            call p%drawPolyline([22,27,24],[5,8,13])
        case(80)
            call ink(p,blue,1.75_dp)
            call p%drawEllipse(10,4,12,7)
            call p%drawArc(4,19,24,9,180._dp,180._dp)
            call ink(p,aqua,1.85_dp)
            call p%drawLine(10,8,4,24); call p%drawLine(22,8,28,24)
        case(81)
            call ink(p,aqua)
            call p%setBrush(QColor(86,187,176,40))
            call p%drawRect(6,4,20,24)
            call ink(p,paper,2.15_dp)
            call p%drawLine(11,16,21,16); call p%drawLine(16,11,16,21)
        case(83)
            call ink(p,paper)
            call p%drawPolyline([5,5,27,27],[17,25,25,17])
            call ink(p,aqua)
            call p%drawLine(16,4,16,19)
            call p%drawPolyline([11,16,21],[14,19,14])
        case(84)
            call p%translate(16,16); call p%rotate(-39._dp)
            call ink(p,aqua)
            call p%setBrush(QColor(86,187,176,45))
            call p%drawRect(-13,-5,26,10)
            call ink(p,aqua,1.55_dp)
            do i=-9,9,6
                call p%drawLine(i,-5,i,0)
            end do
        case(85)
            call ink(p,blue,1.7_dp)
            call p%setBrush(QColor(84,150,198,45))
            call p%drawPolygon([7,11,11,14,14,18,18,21,21,25,25,21,20,24,24,8,8,12,11,7], &
                               [5,5,9,9,5,5,9,9,5,5,12,15,21,24,27,27,24,21,15,12])
            call ink(p,aqua,1.55_dp)
            call p%drawLine(9,24,23,24)
        case(86)
            call ink(p,aqua)
            call p%drawPolyline([4,16,28],[14,4,14])
            call p%drawPolyline([7,7,13,13,19,19,25,25],[12,27,27,19,19,27,27,12])
        case(87)
            call ink(p,blue,1.8_dp)
            call p%drawPolyline([13,5,5,13],[7,10,23,26])
            call p%drawPolyline([19,27,27,19],[7,10,23,26])
            call ink(p,aqua,2.3_dp)
            call p%drawLine(13,7,13,26); call p%drawLine(19,7,19,26)
            call p%drawLine(13,16,19,16)
        case(90)
            call ink(p,blue,1.85_dp)
            call p%drawPolygon([5,12,12,5],[11,7,22,26])
            call p%drawPolygon([20,27,27,20],[7,4,19,22])
            call p%setPen(QPen(QColor(aqua),1.8_dp,style=Qt_DashLine))
            call p%drawLine(16,5,16,26)
        case(icon_check)
            call ink(p,success,2.5_dp)
            call p%drawPolyline([6,13,27],[17,24,8])
        case(91)
            call ink(p,quiet,1.65_dp)
            call p%drawLine(4,21,28,21)
            call ink(p,aqua,2.3_dp)
            call p%drawLine(10,21,22,21)
            call p%drawPolyline([12,9,9],[8,8,25])
            call p%drawPolyline([20,23,23],[8,8,25])
        case(icon_degree_elevate,icon_degree_reduce)
            call ink(p,blue,1.85_dp)
            path=QPainterPath()
            call path%moveTo(5,25); call path%cubicTo(10,25,11,11,19,11)
            call p%drawPath(path); call path%delete()
            call ink(p,merge(aqua,danger,kind==icon_degree_elevate),2.1_dp)
            if(kind==icon_degree_elevate)then
                call p%drawLine(24,17,24,5)
                call p%drawPolyline([20,24,28],[9,5,9])
            else
                call p%drawLine(24,5,24,17)
                call p%drawPolyline([20,24,28],[13,17,13])
            end if
        case(92)
            call ink(p,blue)
            call p%drawLine(5,11,27,11)
            call p%drawPolyline([22,27,22],[6,11,16])
            call ink(p,aqua)
            call p%drawLine(5,22,27,22)
            call p%drawPolyline([10,5,10],[17,22,27])
        case(93)
            call p%drawEllipse(7,7,18,18)
            call ink(p,aqua)
            call p%drawLine(16,4,16,10); call p%drawLine(16,22,16,28)
            call p%drawLine(4,16,10,16); call p%drawLine(22,16,28,16)
            call solid(p,aqua)
            call p%drawEllipse(13,13,6,6)
        case(icon_close)
            call p%drawLine(9,9,23,23); call p%drawLine(9,23,23,9)
        case(icon_minimize)
            call p%drawLine(8,21,24,21)
        case(icon_maximize)
            call p%drawRect(8,8,16,16)
        case(icon_fullscreen)
            call p%drawPolyline([5,5,13],[13,5,5])
            call p%drawPolyline([19,27,27],[5,5,13])
            call p%drawPolyline([5,5,13],[19,27,27])
            call p%drawPolyline([19,27,27],[27,27,19])
        case(icon_restore)
            call p%drawRect(6,12,14,14)
            call p%drawPolyline([12,12,26,26,21],[10,6,6,20,20])
        case(icon_more)
            call solid(p,paper)
            do i=6,22,8
                call p%drawEllipse(i,14,4,4)
            end do
        case(icon_send)
            call solid(p,violet)
            call p%drawPolygon([4,28,20,14,10],[12,4,28,19,25])
            call ink(p,dark,1.8_dp)
            call p%drawLine(14,19,23,9)
        case(icon_stop)
            call ink(p,danger)
            call p%drawEllipse(4,4,24,24)
            call solid(p,danger)
            call p%drawRect(11,11,10,10)
        case(icon_history)
            call ink(p,violet)
            call p%drawArc(5,5,22,22,144._dp,318._dp)
            call p%drawPolyline([4,4,11],[6,13,13])
            call p%drawPolyline([16,16,21],[10,16,19])
        case(icon_new_chat)
            call ink(p,violet)
            path=QPainterPath()
            call path%moveTo(4,8); call path%quadTo(4,5,8,5); call path%lineTo(24,5)
            call path%quadTo(28,5,28,9); call path%lineTo(28,19); call path%quadTo(28,23,24,23)
            call path%lineTo(12,23); call path%lineTo(6,28); call path%lineTo(6,22)
            call path%quadTo(4,21,4,19); call path%closeSubpath()
            call p%drawPath(path); call path%delete()
            call ink(p,'#eee3ff')
            call p%drawLine(11,14,21,14); call p%drawLine(16,9,16,19)
        case(icon_float)
            call p%drawPolyline([11,5,5,22,22],[9,9,26,26,20])
            call ink(p,blue)
            call p%drawPolyline([18,27,27],[5,5,14])
            call p%drawLine(15,17,26,6)
        case(icon_chevron,icon_chevron_up)
            call ink(p,paper,2.5_dp)
            if(kind==icon_chevron_up)then
                call p%drawPolyline([8,16,24],[20,12,20])
            else
                call p%drawPolyline([8,16,24],[12,20,12])
            end if
        case(icon_chevron_right)
            call ink(p,paper,2.5_dp)
            call p%drawPolyline([12,20,12],[8,16,24])
        case(icon_scene_branch_closed,icon_scene_branch_open)
            call solid(p,quiet)
            if(kind==icon_scene_branch_closed)then
                call p%drawPolygon([9,24,9],[5,16,27])
            else
                call p%drawPolygon([5,27,16],[9,9,24])
            end if
        case(icon_face_grid)
            call p%drawPolyline([5,24,28,9,5],[11,7,21,25,11])
            call p%drawLine(11,10,15,24)
            call p%drawLine(18,8,22,22)
            call p%drawLine(7,16,25,12)
            call p%drawLine(8,20,27,16)
        case(icon_isolate)
            call p%drawPolyline([5,5,11],[11,5,5])
            call p%drawPolyline([21,27,27],[5,5,11])
            call p%drawPolyline([5,5,11],[21,27,27])
            call p%drawPolyline([21,27,27],[27,27,21])
            call p%drawPolygon([16,22,22,16,10,10],[8,12,20,24,20,12])
        case(icon_paste)
            call p%drawPolyline([11,7,7,25,25,21],[7,7,27,27,7,7])
            call p%drawRect(11,4,10,6)
            call p%drawLine(11,16,21,16)
            call p%drawLine(11,21,19,21)
        case(icon_scene_filter)
            call p%drawPolyline([4,28,19,19,13,13,4],[6,6,16,22,27,16,6])
        case(icon_export)
            call p%drawPolyline([5,5,27,27],[19,26,26,19])
            call ink(p,aqua)
            call p%drawLine(16,6,16,21)
            call p%drawPolyline([10,16,22],[12,6,12])
        case(icon_delete)
            call ink(p,danger)
            call p%drawLine(5,8,27,8)
            call p%drawPolyline([11,12,20,21],[8,4,4,8])
            call p%drawPolyline([8,9,23,24],[9,27,27,9])
            call p%drawLine(13,13,13,22)
            call p%drawLine(19,13,19,22)
        case(icon_mesh,icon_generate_mesh,icon_mesh_size,icon_mesh_pick,icon_mesh_shaded)
            call mesh_symbol(p,kind)
        case default
            continue
        end select
    end subroutine paint_icon

    subroutine ink(p,color,width)
        type(QPainter), intent(inout) :: p
        character(*), intent(in) :: color
        real(dp), intent(in), optional :: width

        call p%setBrush(QBrush(style=Qt_NoBrush))
        if(present(width))then
            call p%setPen(QPen(QColor(color),width))
        else
            call p%setPen(QPen(QColor(color),2._dp))
        end if
    end subroutine ink

    subroutine solid(p,color)
        type(QPainter), intent(inout) :: p
        character(*), intent(in) :: color

        call p%setPen(QPen(style=Qt_NoPen))
        call p%setBrush(QColor(color))
    end subroutine solid

    subroutine spark(p,x,y,size,waist)
        type(QPainter), intent(inout) :: p
        integer, intent(in) :: x,y,size,waist

        call p%drawPolygon([x,x+waist,x+size,x+waist,x,x-waist,x-size,x-waist], &
            [y-size,y-waist,y,y+waist,y+size,y+waist,y,y-waist])
    end subroutine spark

    subroutine cube(p,kind)
        type(QPainter), intent(inout) :: p
        integer, intent(in) :: kind
        integer :: highlighted
        logical :: hidden_face

        if(kind==22)then
            call control_net(p)
            return
        end if
        highlighted=0
        hidden_face=.false.
        select case(kind)
        case(8,74)
            call cube_face(p,face_top,QColor(103,175,240,188))
            call cube_face(p,face_left,QColor(52,108,187,144))
            call cube_face(p,face_right,QColor(72,143,223,164))
        case(icon_camera_front,icon_bounding_face)
            highlighted=face_right
        case(icon_camera_back)
            highlighted=face_right
            hidden_face=.true.
        case(icon_camera_left)
            highlighted=face_left
        case(icon_camera_right)
            highlighted=face_left
            hidden_face=.true.
        case(icon_camera_top)
            highlighted=face_top
        case(icon_camera_bottom)
            highlighted=face_top
            hidden_face=.true.
        case default
            continue
        end select
        if(highlighted>0.and..not.hidden_face)call cube_face(p,highlighted,QColor(73,198,181,139))
        call cube_outline(p,blue)
        if(highlighted>0)then
            if(hidden_face)then
                call p%setPen(QPen(QColor(aqua),2.05_dp,style=Qt_DashLine))
                select case(highlighted)
                case(face_top)
                    call p%drawPolyline([5,16,27],[22,16,22])
                    call p%drawPolyline([5,16,27],[22,28,22])
                case(face_left)
                    call p%drawPolyline([16,27,27,16],[4,10,22,16])
                case(face_right)
                    call p%drawPolyline([16,5,5,16],[4,10,22,16])
                case default
                    continue
                end select
            else
                call outline_cube_face(p,highlighted)
            end if
        end if

        select case(kind)
        case(25)
            call ink(p,aqua,1.85_dp)
            call p%drawPolyline([10,16,22],[13,16,13])
            call p%drawLine(16,16,16,24)
        case(icon_bounding_edge)
            call ink(p,aqua,3.15_dp)
            call p%drawLine(16,16,27,10)
        case(icon_bounding_face)
            call solid(p,aqua)
            call p%drawPolygon([19,24,24,19],[17,14,20,23])
        case(icon_bounding_vertex)
            call solid(p,aqua)
            call p%drawEllipse(12,12,8,8)
        case(icon_extract_section)
            call ink(p,aqua,1.85_dp)
            call p%setBrush(QColor(73,198,181,119))
            call p%drawPolygon([5,16,27,16],[16,10,16,22])
        case default
            continue
        end select
    end subroutine cube

    subroutine cube_outline(p,color)
        type(QPainter), intent(inout) :: p
        character(*), intent(in) :: color

        call ink(p,color,1.85_dp)
        call p%drawPolygon([5,16,27,27,16,5],[10,4,10,22,28,22])
        call p%drawPolyline([5,16,27],[10,16,10])
        call p%drawLine(16,16,16,28)
    end subroutine cube_outline

    subroutine cube_face(p,face,color)
        type(QPainter), intent(inout) :: p
        integer, intent(in) :: face
        type(QColor), intent(in) :: color
        integer, parameter :: xs(4,3)=reshape([5,16,27,16, 5,16,16,5, 16,27,27,16],[4,3])
        integer, parameter :: ys(4,3)=reshape([10,4,10,16, 10,16,28,22, 16,10,22,28],[4,3])

        call p%setPen(QPen(style=Qt_NoPen))
        call p%setBrush(color)
        call p%drawPolygon(xs(:,face),ys(:,face))
    end subroutine cube_face

    subroutine outline_cube_face(p,face)
        type(QPainter), intent(inout) :: p
        integer, intent(in) :: face

        call ink(p,aqua,2.2_dp)
        select case(face)
        case(face_top)
            call p%drawPolyline([5,16,27],[10,4,10])
        case(face_left)
            call p%drawPolyline([5,5,16],[10,22,28])
        case(face_right)
            call p%drawPolyline([27,27,16],[10,22,28])
        case default
            continue
        end select
    end subroutine outline_cube_face

    subroutine control_net(p)
        type(QPainter), intent(inout) :: p
        integer, parameter :: x(4)=[5,22,27,10], y(4)=[10,6,21,25]
        integer :: i

        call ink(p,blue,1.7_dp)
        call p%drawPolygon(x,y)
        call p%drawLine(13,8,19,23)
        call p%drawLine(8,18,25,14)
        call solid(p,aqua)
        do i=1,size(x)
            call p%drawRect(x(i)-2,y(i)-2,4,4)
        end do
    end subroutine control_net

    subroutine mesh_symbol(p,kind)
        type(QPainter), intent(inout) :: p
        integer, intent(in) :: kind

        if(kind==icon_mesh_pick.or.kind==icon_mesh_size)then
            call ink(p,blue,1.8_dp)
            call p%drawRect(5,5,22,merge(22,15,kind==icon_mesh_pick))
            call p%drawLine(16,5,16,merge(27,20,kind==icon_mesh_pick))
            call p%drawLine(5,merge(16,12,kind==icon_mesh_pick),27,merge(16,12,kind==icon_mesh_pick))
            call ink(p,aqua,2.1_dp)
            if(kind==icon_mesh_pick)then
                call p%setBrush(QColor(73,198,181,132))
                call p%drawRect(16,5,11,11)
            else
                call p%drawLine(6,25,26,25)
                call p%drawPolyline([9,6,9],[22,25,28])
                call p%drawPolyline([23,26,23],[22,25,28])
            end if
            return
        end if

        call p%save()
        if(kind==icon_generate_mesh)then
            call p%translate(-1,5)
            call p%scale(.84_dp,.84_dp)
        end if
        if(kind==icon_mesh_shaded)then
            call cube_face(p,face_top,QColor(236,160,67,179))
            call cube_face(p,face_left,QColor(181,110,42,139))
            call cube_face(p,face_right,QColor(214,134,55,167))
        end if
        call cube_outline(p,amber)
        call ink(p,amber,1.35_dp)
        call p%drawPolyline([10,22,22],[7,13,25])
        call p%drawPolyline([22,10,10],[7,13,25])
        call p%drawPolyline([5,16,27],[16,22,16])
        call p%restore()
        if(kind==icon_generate_mesh)then
            call solid(p,aqua)
            call spark(p,25,7,6,2)
        end if
    end subroutine mesh_symbol

    subroutine primitive(p,kind)
        type(QPainter), intent(inout) :: p
        integer, intent(in) :: kind
        type(QPainterPath) :: path

        call ink(p,blue)
        select case(kind)
        case(69)
            call p%drawEllipse(6,6,20,20)
        case(70)
            call p%drawLine(7,24,25,8)
        case(71)
            call p%drawEllipse(5,5,22,22)
        case(72)
            path=QPainterPath()
            call path%moveTo(5,22); call path%cubicTo(8,6,21,5,27,16)
            call p%drawPath(path); call path%delete()
        case(73)
            call p%drawRect(5,8,22,16)
        case default
            continue
        end select
        call solid(p,aqua)
        select case(kind)
        case(69)
            call p%drawEllipse(12,12,8,8)
        case(70)
            call p%drawEllipse(4,21,6,6); call p%drawEllipse(22,5,6,6)
        case(72)
            call p%drawEllipse(3,20,4,4); call p%drawEllipse(25,14,4,4)
        case(73)
            call p%drawRect(3,22,4,4); call p%drawRect(25,6,4,4)
        case default
            continue
        end select
    end subroutine primitive
end module studio_icons
