! Paint read-only ForCAD data. No spline evaluation is performed here.
module studio_live_plots
    use qt, only: QBrush, QColor, QFont, QFontMetrics, QPainter, QPainter_Antialiasing, QPainter_TextAntialiasing, &
        QPen, Qt_DashDotDotLine, Qt_DashDotLine, Qt_DashLine, Qt_DotLine, Qt_NoBrush, Qt_NoPen, Qt_SolidLine
    use studio_state, only: studio_t, dp, axis_rgb
    use studio_theme, only: studio_divider, studio_frame_edge, studio_heading, studio_secondary_text
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none(type, external)
    private
    public :: paint_live_basis

    integer, parameter :: key_limit=8
    real(dp), parameter :: chart_text_size=10._dp, axis_title_size=10.5_dp
    real(dp), parameter :: legend_swatch_start=13._dp, legend_swatch_end=36._dp, legend_label_offset=43._dp
    integer, parameter :: line_styles(5)=[Qt_SolidLine,Qt_DashLine,Qt_DashDotLine,Qt_DotLine,Qt_DashDotDotLine]
contains
    subroutine paint_live_basis(painter,s)
        type(QPainter), intent(inout) :: painter
        type(studio_t), intent(in) :: s
        type(QColor) :: colors(key_limit), axis_color
        type(QFontMetrics) :: metrics
        real(dp), allocatable :: xs(:),ys(:),peaks(:)
        real(dp) :: width,height,xlo,xhi,ylo,yhi,step,l,r,t,b,probe
        real(dp) :: tick_width,font_height,ascent,descent,title_x,title_half,legend_width,legend_step,legend_first
        real(dp) :: x,y,value,a,legend_text_width,horizontal_tick_width,probe_width,probe_x,badge_left,row_pitch
        integer :: axis,npoint,nfunction,nactive,nkey,keys(key_limit),i,j,k,ny,nx,legend_rows,legend_capacity,y_stride
        character :: axis_name
        character(:), allocatable :: title,probe_text,legend_summary
        logical :: derivative,show_probe

        width=real(s%basis%width(),dp)
        height=real(s%basis%height(),dp)
        call backdrop(painter,width,height)
        if (.not.s%live_scene%valid.or..not.s%live_plot%valid) then
            call unavailable(painter,width,height)
            return
        end if
        if (.not.allocated(s%live_plot%parameters).or..not.allocated(s%live_plot%values)) then
            call unavailable(painter,width,height)
            return
        end if
        npoint=size(s%live_plot%parameters)
        nfunction=size(s%live_plot%values,2)
        if (npoint<1.or.nfunction<1.or.size(s%live_plot%values,1)/=npoint) then
            call unavailable(painter,width,height)
            return
        end if
        if (.not.all(ieee_is_finite(s%live_plot%parameters)).or. &
            .not.all(ieee_is_finite(s%live_plot%values))) then
            call unavailable(painter,width,height)
            return
        end if
        if (s%live_scene%dimension<1.or.s%live_scene%dimension>3) then
            call unavailable(painter,width,height)
            return
        end if
        derivative=s%basis_derivative
        axis=s%knot_axis
        if (axis<1.or.axis>s%live_scene%dimension) then
            call unavailable(painter,width,height)
            return
        end if
        axis_name=achar(iachar('u')+axis-1)
        axis_color=QColor(axis_rgb(1,axis),axis_rgb(2,axis),axis_rgb(3,axis))
        xlo=s%live_scene%bounds(1,axis)
        xhi=s%live_scene%bounds(2,axis)
        if (.not.all(ieee_is_finite([xlo,xhi])).or.xhi<=xlo) then
            call unavailable(painter,width,height)
            return
        end if
        if (.not.ieee_is_finite(xhi-xlo)) then
            call unavailable(painter,width,height)
            return
        end if
        allocate(peaks(nfunction),xs(npoint),ys(npoint))
        do j=1,nfunction
            peaks(j)=maxval(abs(s%live_plot%values(:,j)))
        end do
        nactive=count(peaks>0._dp)
        colors=[QColor('#76d8ee'),QColor('#f6be6c'),QColor('#9de49b'),QColor('#b29cff'), &
            QColor('#f09ac5'),QColor('#e3de86'),QColor('#64d8bc'),QColor('#e8a189')]
        probe=s%live_probe%parameters(axis)
        show_probe=s%basis_probe_guide.and.s%live_probe%valid
        show_probe=show_probe.and.ieee_is_finite(probe).and.probe>=xlo.and.probe<=xhi
        probe_text=''
        if(show_probe)probe_text=achar(iachar('U')+axis-1)//' = '//tick_text(probe)

        call painter%setFont(QFont(family='Noto Sans',pointSize=chart_text_size))
        metrics=painter%fontMetrics()
        font_height=real(metrics%height(),dp)
        ascent=real(metrics%ascent(),dp)
        descent=real(metrics%descent(),dp)
        probe_width=real(metrics%horizontalAdvance(probe_text),dp)+14._dp
        t=max(12._dp,.5_dp*font_height+4._dp)
        if(show_probe)t=font_height+15._dp
        b=height-2._dp*font_height-13._dp
        if (b-t<24._dp) then
            call metrics%delete()
            call unavailable(painter,width,height)
            return
        end if
        call ordinate_limits(s%live_plot%values,max(2,min(5,int((b-t)/(font_height+6._dp)))),derivative,ylo,yhi,step)
        if (.not.all(ieee_is_finite([ylo,yhi,step])).or.yhi<=ylo.or.step<=0._dp) then
            call metrics%delete()
            call unavailable(painter,width,height)
            return
        end if
        ny=max(1,min(12,nint((yhi-ylo)/step)))
        row_pitch=(b-t)/ny
        y_stride=max(1,ceiling((font_height+3._dp)/row_pitch))
        tick_width=0._dp
        do i=0,ny
            tick_width=max(tick_width,real(metrics%horizontalAdvance(tick_text(ylo+real(i,dp)*step)),dp))
        end do
        horizontal_tick_width=0._dp
        do i=0,4
            value=xlo+(.25_dp*i)*(xhi-xlo)
            horizontal_tick_width=max(horizontal_tick_width,real(metrics%horizontalAdvance(tick_text(value)),dp))
        end do
        ! Reserve a line for the count of dimmed traces when the key cannot fit them all.
        legend_step=max(16._dp,font_height+2._dp)
        legend_capacity=max(1,int((height-12._dp)/legend_step))
        if(nactive>min(key_limit,legend_capacity))legend_capacity=max(1,legend_capacity-1)
        nkey=min(nactive,key_limit,legend_capacity)
        call choose_keys(peaks,keys,nkey)
        legend_summary=''
        if(nactive>nkey)legend_summary='+'//integer_text(nactive-nkey)//' dimmed'
        legend_rows=nkey
        if(len(legend_summary)>0)legend_rows=legend_rows+1
        legend_text_width=0._dp
        legend_width=max(9._dp,.5_dp*horizontal_tick_width+4._dp)
        if(s%basis_legend.and.nkey>0)then
            do k=1,nkey
                legend_text_width=max(legend_text_width,real(metrics%horizontalAdvance(series_text(keys(k),derivative)),dp))
            end do
            legend_width=max(legend_width,legend_text_width+legend_label_offset+4._dp, &
                real(metrics%horizontalAdvance(legend_summary),dp)+12._dp)
        end if
        call metrics%delete()
        title='Rᵢ'
        if (derivative) title='∂Rᵢ/∂'//axis_name
        call painter%setFont(QFont(family='Noto Sans',pointSize=axis_title_size,italic=.true.))
        metrics=painter%fontMetrics()
        title_x=5._dp+real(metrics%ascent(),dp)
        title_half=.5_dp*real(metrics%horizontalAdvance(title),dp)
        l=max(45._dp,title_x+real(metrics%descent(),dp)+tick_width+14._dp)
        call metrics%delete()
        r=width-legend_width-5._dp
        if (r-l<45._dp) then
            call unavailable(painter,width,height)
            return
        end if
        nx=4
        if (r-l<max(240._dp,4._dp*(horizontal_tick_width+8._dp))) nx=2
        if (r-l<2._dp*(horizontal_tick_width+8._dp)) nx=1
        call painter%fillRect(l,t,r-l,b-t,QColor('#17232c'))
        call painter%setBrush(QBrush(style=Qt_NoBrush))
        call painter%setPen(QPen(QColor(studio_divider),1._dp))
        call painter%drawRect(l,t,r-l,b-t)
        do i=1,nx-1
            x=l+real(i,dp)*(r-l)/real(nx,dp)
            call painter%drawLine(x,t,x,b)
        end do
        do i=1,ny-1
            value=ylo+real(i,dp)*step
            y=b-(value-ylo)/(yhi-ylo)*(b-t)
            if (abs(value)<step*1.e-8_dp) then
                call painter%setPen(QPen(QColor(123,150,161,165),1._dp,Qt_DashLine))
            else
                call painter%setPen(QPen(QColor(studio_divider),.8_dp))
            end if
            call painter%drawLine(l,y,r,y)
        end do
        call painter%setPen(QPen(QColor(169,189,202),1._dp))
        call painter%drawLine(l,t,l,b)
        call painter%drawLine(l,b,r,b)
        do i=0,ny
            if (i/=ny.and.mod(i,y_stride)/=0) cycle
            if (i>0.and.i<ny.and.real(ny-i,dp)*row_pitch<font_height+3._dp) cycle
            value=ylo+real(i,dp)*step
            if (abs(value)<step*1.e-8_dp) value=0._dp
            y=b-real(i,dp)*(b-t)/real(ny,dp)
            call painter%setPen(QPen(QColor(190,207,222),.8_dp))
            call painter%drawLine(l-4._dp,y,l,y)
            call label_right(painter,l-7._dp,y,tick_text(value),chart_text_size)
        end do
        do i=0,nx
            a=real(i,dp)/real(nx,dp)
            x=l+a*(r-l)
            call painter%setPen(QPen(QColor(190,207,222),.8_dp))
            call painter%drawLine(x,b,x,b+4._dp)
            call label_center(painter,x,b+font_height+4._dp,tick_text(xlo+a*(xhi-xlo)),chart_text_size)
        end do

        ! Point arrays are already evaluated by ForCAD. Only map them to pixels.
        xs=l+(s%live_plot%parameters-xlo)/(xhi-xlo)*(r-l)
        call painter%save()
        call painter%setClipRect(l-1._dp,t-2._dp,r-l+2._dp,b-t+3._dp)
        do j=1,nfunction
            if (any(keys(:nkey)==j).or.peaks(j)==0._dp) cycle
            ys=b-(s%live_plot%values(:,j)-ylo)/(yhi-ylo)*(b-t)
            call painter%setPen(QPen(QColor(151,181,191,85),1._dp,line_styles(mod(j-1,5)+1)))
            call draw_trace(painter,xs,ys)
        end do
        do k=1,nkey
            j=keys(k)
            ys=b-(s%live_plot%values(:,j)-ylo)/(yhi-ylo)*(b-t)
            call painter%setPen(QPen(colors(mod(j-1,key_limit)+1),1.9_dp,line_styles(mod(j-1,5)+1)))
            call draw_trace(painter,xs,ys)
        end do
        if (show_probe) then
            probe_x=l+(probe-xlo)/(xhi-xlo)*(r-l)
            call painter%setPen(QPen(QColor(axis_rgb(1,axis),axis_rgb(2,axis),axis_rgb(3,axis),160),1._dp,Qt_DotLine))
            call painter%drawLine(probe_x,t,probe_x,b)
        end if
        call painter%restore()
        if (show_probe) then
            call painter%setPen(QPen(style=Qt_NoPen))
            call painter%setBrush(axis_color)
            call painter%drawPolygon([probe_x-3.5_dp,probe_x+3.5_dp,probe_x],[t-3._dp,t-3._dp,t+3._dp])
            badge_left=max(7._dp,min(r-probe_width,probe_x-.5_dp*probe_width))
            call painter%setPen(QPen(QColor(axis_rgb(1,axis),axis_rgb(2,axis),axis_rgb(3,axis),88),1._dp))
            call painter%setBrush(QColor(studio_heading))
            call painter%drawRect(badge_left,5._dp,probe_width,font_height+4._dp)
            call label(painter,badge_left+7._dp,7._dp+ascent,probe_text,chart_text_size,axis_color)
        end if

        if(s%basis_legend.and.nkey>0)then
            legend_first=.5_dp*(t+b)-.5_dp*real(legend_rows-1,dp)*legend_step
            legend_first=max(.5_dp*font_height+6._dp,min(legend_first, &
                height-6._dp-.5_dp*font_height-real(legend_rows-1,dp)*legend_step))
            call painter%setPen(QPen(QColor(studio_frame_edge),.8_dp))
            call painter%drawLine(r+7._dp,max(5._dp,legend_first-.5_dp*font_height),r+7._dp, &
                min(height-5._dp,legend_first+real(legend_rows-1,dp)*legend_step+.5_dp*font_height))
            do k=1,nkey
                j=keys(k)
                y=legend_first+real(k-1,dp)*legend_step
                call painter%setPen(QPen(colors(mod(j-1,key_limit)+1),1.9_dp,line_styles(mod(j-1,5)+1)))
                call painter%drawLine(r+legend_swatch_start,y,r+legend_swatch_end,y)
                call label(painter,r+legend_label_offset,y+.5_dp*(ascent-descent),series_text(j,derivative),chart_text_size)
            end do
            if (len(legend_summary)>0) then
                call label(painter,r+9._dp,legend_first+nkey*legend_step+.5_dp*(ascent-descent), &
                           legend_summary,9.5_dp,QColor(studio_secondary_text))
            end if
        end if
        call label_center(painter,.5_dp*(l+r),height-4._dp,axis_name,axis_title_size,axis_color,italic=.true.)
        call painter%save()
        call painter%translate(title_x,.5_dp*(t+b)+title_half)
        call painter%rotate(-90._dp)
        call label(painter,0._dp,0._dp,title,axis_title_size,italic=.true.)
        call painter%restore()
    end subroutine paint_live_basis

    subroutine ordinate_limits(values,divisions,derivative,lower,upper,step)
        real(dp), intent(in) :: values(:,:)
        integer, intent(in) :: divisions
        logical, intent(in) :: derivative
        real(dp), intent(out) :: lower,upper,step
        real(dp) :: raw,magnitude,fraction,span
        integer :: exponent_value
        lower=min(0._dp,minval(values))
        upper=max(0._dp,maxval(values))
        if (upper==lower) then
            upper=1._dp
            if (derivative) lower=-1._dp
        end if
        span=upper-lower
        step=0._dp
        if (.not.ieee_is_finite(span).or.span<=0._dp) return
        raw=max(tiny(1._dp),span/real(divisions,dp))
        exponent_value=floor(log10(raw))
        magnitude=10._dp**exponent_value
        if (magnitude==0._dp.or..not.ieee_is_finite(magnitude)) return
        fraction=raw/magnitude
        if (fraction<=1._dp) then
            step=magnitude
        else if (fraction<=2._dp) then
            step=2._dp*magnitude
        else if (fraction<=2.5_dp) then
            step=2.5_dp*magnitude
        else if (fraction<=5._dp) then
            step=5._dp*magnitude
        else
            step=10._dp*magnitude
        end if
        lower=real(floor(lower/step),dp)*step
        upper=real(ceiling(upper/step),dp)*step
    end subroutine ordinate_limits

    subroutine choose_keys(peaks,keys,count)
        real(dp), intent(in) :: peaks(:)
        integer, intent(in) :: count
        integer, intent(out) :: keys(key_limit)
        real(dp), allocatable :: remaining(:)
        integer :: i,j,selected,swap
        keys=0
        remaining=peaks
        do i=1,count
            selected=maxloc(remaining,dim=1)
            keys(i)=selected
            remaining(selected)=-1._dp
        end do
        do i=2,count
            j=i
            do while (j>1)
                if (keys(j-1)<keys(j)) exit
                swap=keys(j-1); keys(j-1)=keys(j); keys(j)=swap
                j=j-1
            end do
        end do
    end subroutine choose_keys

    subroutine draw_trace(p,x,y)
        type(QPainter), intent(inout) :: p
        real(dp), intent(in) :: x(:),y(:)
        if (size(x)==1) then
            call p%drawEllipse(x(1)-1.5_dp,y(1)-1.5_dp,3._dp,3._dp)
        else
            call p%drawPolyline(x,y)
        end if
    end subroutine draw_trace

    subroutine backdrop(p,w,h)
        type(QPainter), intent(inout) :: p
        real(dp), intent(in) :: w,h
        call p%fillRect(0._dp,0._dp,w,h,QColor('#131e27'))
        call p%setRenderHint(QPainter_Antialiasing)
        call p%setRenderHint(QPainter_TextAntialiasing)
    end subroutine backdrop

    subroutine unavailable(p,w,h)
        type(QPainter), intent(inout) :: p
        real(dp), intent(in) :: w,h
        call label_center(p,.5_dp*w,.5_dp*h,'Basis functions unavailable',10._dp,color=QColor(175,195,206))
    end subroutine unavailable

    subroutine label(p,x,y,text,size,color,italic)
        type(QPainter), intent(inout) :: p
        real(dp), intent(in) :: x,y,size
        character(*), intent(in) :: text
        type(QColor), intent(in), optional :: color
        logical, intent(in), optional :: italic
        logical :: slanted
        slanted=.false.
        if (present(italic)) slanted=italic
        call p%setFont(QFont(family='Noto Sans',pointSize=size,italic=slanted))
        call p%setPen(QColor(237,243,247))
        if (present(color)) call p%setPen(color)
        call p%drawText(x,y,text)
    end subroutine label

    subroutine label_center(p,x,y,text,size,color,italic)
        type(QPainter), intent(inout) :: p
        real(dp), intent(in) :: x,y,size
        character(*), intent(in) :: text
        type(QColor), intent(in), optional :: color
        logical, intent(in), optional :: italic
        type(QFontMetrics) :: metrics
        logical :: slanted
        real(dp) :: left
        slanted=.false.
        if (present(italic)) slanted=italic
        call p%setFont(QFont(family='Noto Sans',pointSize=size,italic=slanted))
        metrics=p%fontMetrics()
        left=x-.5_dp*real(metrics%horizontalAdvance(text),dp)
        call metrics%delete()
        call label(p,left,y,text,size,color,slanted)
    end subroutine label_center

    subroutine label_right(p,x,y,text,size)
        type(QPainter), intent(inout) :: p
        real(dp), intent(in) :: x,y,size
        character(*), intent(in) :: text
        type(QFontMetrics) :: metrics
        real(dp) :: left,baseline
        call p%setFont(QFont(family='Noto Sans',pointSize=size))
        metrics=p%fontMetrics()
        left=x-real(metrics%horizontalAdvance(text),dp)
        baseline=y+.5_dp*real(metrics%ascent()-metrics%descent(),dp)
        call metrics%delete()
        call label(p,left,baseline,text,size)
    end subroutine label_right

    function integer_text(value) result(text)
        integer, intent(in) :: value
        character(:), allocatable :: text
        character(32) :: buffer
        write(buffer,'(i0)') value
        text=trim(buffer)
    end function integer_text

    function subscript(value) result(text)
        integer, intent(in) :: value
        character(:), allocatable :: text,plain
        character(3), parameter :: digits(0:9)=[character(3)::'₀','₁','₂','₃','₄','₅','₆','₇','₈','₉']
        integer :: i,digit
        plain=integer_text(value)
        text=''
        do i=1,len(plain)
            digit=iachar(plain(i:i))-iachar('0')
            if (digit>=0.and.digit<=9) text=text//digits(digit)
        end do
    end function subscript

    function series_text(column,derivative) result(text)
        integer, intent(in) :: column
        logical, intent(in) :: derivative
        character(:), allocatable :: text
        text='R'
        if (derivative) text=text//'′'
        text=text//subscript(column-1)
    end function series_text

    function tick_text(value) result(text)
        real(dp), intent(in) :: value
        character(:), allocatable :: text,mantissa,exponent_text
        character(40) :: buffer
        integer :: e,exponent_value,last
        if (value==0._dp) then
            text='0'
            return
        end if
        if (abs(value)>=10000._dp.or.abs(value)<.001_dp) then
            write(buffer,'(es18.2e3)') value
            text=trim(adjustl(buffer))
            e=index(text,'E')
            if (e==0) return
            mantissa=text(:e-1)
            read(text(e+1:),*) exponent_value
            exponent_text=integer_text(exponent_value)
            last=len(mantissa)
            do while (last>1)
                if (mantissa(last:last)/='0') exit
                last=last-1
            end do
            if (mantissa(last:last)=='.') last=last-1
            text=mantissa(:last)//'e'//exponent_text
        else
            write(buffer,'(f20.4)') value
            text=trim(adjustl(buffer))
            last=len(text)
            do while (last>1)
                if (text(last:last)/='0') exit
                last=last-1
            end do
            if (text(last:last)=='.') last=last-1
            text=text(:last)
        end if
    end function tick_text
end module studio_live_plots
