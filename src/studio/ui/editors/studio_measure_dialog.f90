module studio_measure_dialog
    use qt, only: QDialog, QDialogButtonBox, QDialogButtonBox_Close, QFrame, QFrame_HLine, QLabel, QVBoxLayout, Qt_PlainText, QWidget
    use studio_document, only: studio_document_t
    use studio_presentation_text, only: number_text
    use studio_measurements, only: measurement_t, measurement_label, measurement_unit, measurement_warning
    implicit none(type, external)
    private
    public :: show_geometry_measure
contains
    subroutine show_geometry_measure(parent,document,measures,name)
        class(QWidget), intent(in) :: parent
        type(studio_document_t), intent(in) :: document
        type(measurement_t), intent(in) :: measures(3)
        character(*), intent(in) :: name
        type(QDialog), target :: dialog
        type(QDialogButtonBox) :: buttons
        type(QVBoxLayout) :: layout
        type(QLabel) :: heading, readout
        type(QFrame) :: separator
        character(:), allocatable :: result,warning,units
        integer :: dismissed,kind,rows

        dialog=QDialog(parent)
        call dialog%setObjectName('studioMeasureDialog')
        call dialog%setWindowTitle('Measurements')
        call dialog%setMinimumWidth(336)
        layout=QVBoxLayout(dialog)
        call layout%setContentsMargins(18,16,18,16)
        call layout%setSpacing(10)
        if(name/='')then
            heading=QLabel(name,dialog)
            call heading%setObjectName('studioDialogTitle')
            call heading%setTextFormat(Qt_PlainText)
            call heading%setWordWrap(.true.)
            call layout%addWidget(heading)
        end if
        rows=0
        do kind=1,3
            associate(measure=>measures(kind))
                if(measure%count==0)cycle
                if(rows>0)then
                    separator=QFrame(dialog)
                    call separator%setFrameShape(QFrame_HLine)
                    call layout%addWidget(separator)
                end if
                if(measure%error=='')then
                    units=measurement_unit(document,kind)
                    result=measurement_label(kind)//'  '//number_text(measure%value)
                    if(units/='')result=result//' '//units
                else
                    result=measure%error
                end if
                readout=QLabel(result,dialog)
                if(rows==0)call readout%setObjectName('studioMeasureResult')
                call readout%setTextFormat(Qt_PlainText)
                call readout%setWordWrap(.true.)
                if(measure%error=='')then
                    call readout%setStyleSheet('font-size:18px; font-weight:650; color:#83ead9;')
                else
                    call readout%setProperty('dialogRole','error')
                end if
                call layout%addWidget(readout)
                rows=rows+1
                warning=measurement_warning(measure)
                if(len(warning)==0.or.len(measure%error)>0)cycle
                readout=QLabel(warning,dialog)
                call readout%setObjectName('studioMeasureNote')
                call readout%setTextFormat(Qt_PlainText)
                call readout%setWordWrap(.true.)
                call readout%setProperty('dialogRole','detail')
                if(measure%reversed)call readout%setStyleSheet('color:#f4ce8b;')
                call layout%addWidget(readout)
            end associate
        end do
        if(rows==0)then
            readout=QLabel('Select a curve, edge, surface, face, volume, or assembly to measure.',dialog)
            call readout%setWordWrap(.true.)
            call readout%setProperty('dialogRole','detail')
            call layout%addWidget(readout)
        end if
        buttons=QDialogButtonBox(QDialogButtonBox_Close,dialog)
        call buttons%rejected%connect(close_measure,dialog)
        call layout%addWidget(buttons)
        dismissed=dialog%exec()
        call dialog%delete()
    end subroutine show_geometry_measure
    subroutine close_measure(data)
        class(*), optional, intent(inout) :: data

        select type(data)
        type is(QDialog)
            call data%reject()
        end select
    end subroutine close_measure
end module studio_measure_dialog
