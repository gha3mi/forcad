module studio_viewport_depth
    use, intrinsic :: iso_fortran_env, only: int8
    use, intrinsic :: iso_c_binding, only: dp=>c_double,c_int32_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use qt, only: QColor,QImage,QPainter,QPainter_SmoothPixmapTransform
    implicit none(type, external)
    private

    integer, parameter :: max_samples=4
    real(dp), parameter :: sample_x(max_samples)=[.25_dp,.75_dp,.25_dp,.75_dp]
    real(dp), parameter :: sample_y(max_samples)=[.25_dp,.25_dp,.75_dp,.75_dp]
    real(dp), parameter :: empty_depth=-huge(1._dp)

    type, public :: viewport_depth_t
        private
        real(dp),allocatable :: depth(:,:,:)
        integer(c_int32_t),allocatable :: colors(:,:,:)
        integer :: samples=0,pixel_step=1
        real(dp) :: offsets(2,max_samples)=0._dp,sample_margin=0._dp
        integer :: origin(2)=[0,0], extent(2)=[0,0], first(2)=[1,1], last(2)=[0,0]
    contains
        procedure :: initialize
        procedure :: ready
        procedure :: cover_full_frame
        procedure :: triangle
        procedure :: flat_quad
        procedure :: line
        procedure :: grid_edges
        procedure :: paint
        procedure :: depth_at
    end type viewport_depth_t

contains

    subroutine initialize(self,width,height,x_origin,y_origin,multisample,pixel_step)
        class(viewport_depth_t),intent(out) :: self
        integer,intent(in) :: width,height
        integer,intent(in),optional :: x_origin,y_origin,pixel_step
        logical,intent(in),optional :: multisample
        integer :: allocation_status

        if(width<=0.or.height<=0)return
        if(present(pixel_step))self%pixel_step=max(1,pixel_step)
        self%extent=1+([width,height]-1)/self%pixel_step
        self%sample_margin=8._dp*spacing(real(maxval(self%extent),dp))
        self%samples=max_samples
        if(present(multisample))then
            if(.not.multisample)self%samples=1
        end if
        if(self%samples==1)then
            self%offsets(:,1)=.5_dp
        else
            self%offsets(1,:)=sample_x
            self%offsets(2,:)=sample_y
        end if
        allocate(self%depth(self%samples,self%extent(1),self%extent(2)),stat=allocation_status)
        if(allocation_status/=0)return
        allocate(self%colors(self%samples,self%extent(1),self%extent(2)),stat=allocation_status)
        if(allocation_status/=0)then
            deallocate(self%depth)
            return
        end if
        self%depth=empty_depth
        self%colors=0_c_int32_t
        if(present(x_origin))self%origin(1)=x_origin
        if(present(y_origin))self%origin(2)=y_origin
        self%first=self%extent
    end subroutine initialize

    logical function ready(self)
        class(viewport_depth_t),intent(in) :: self

        ready=allocated(self%depth)
    end function ready

    subroutine cover_full_frame(self)
        class(viewport_depth_t),intent(inout) :: self

        self%first=1
        self%last=self%extent
    end subroutine cover_full_frame

    subroutine triangle(self,points,colors,opacity,rows)
        class(viewport_depth_t),intent(inout) :: self
        real(dp),intent(in) :: points(3,3),colors(3,3)
        integer,intent(in) :: opacity
        integer,optional,intent(in) :: rows(2)

        if(self%samples==1)then
            call raster_triangle(self,points,opacity,1,colors=colors,rows=rows)
        else
            call raster_triangle(self,points,opacity,max_samples,colors=colors,rows=rows)
        end if
    end subroutine triangle

    subroutine flat_quad(self,points,rgb,opacity,rows)
        class(viewport_depth_t),intent(inout) :: self
        real(dp),intent(in) :: points(3,4),rgb(3)
        integer,intent(in) :: opacity
        integer,optional,intent(in) :: rows(2)
        integer(c_int32_t) :: fragment
        integer :: alpha

        alpha=min(opacity,255)
        if(.not.self%ready().or.alpha<=0.or..not.all(ieee_is_finite(rgb)))return
        fragment=fragment_color(rgb,alpha)
        call raster_triangle(self,points(:,:3),alpha,self%samples,fill=fragment,rows=rows)
        call raster_triangle(self,points(:,[1,3,4]),alpha,self%samples,fill=fragment,rows=rows)
    end subroutine flat_quad

    subroutine raster_triangle(self,points,opacity,sample_count,colors,fill,rows)
        class(viewport_depth_t),intent(inout) :: self
        real(dp),intent(in) :: points(3,3)
        integer,intent(in) :: opacity
        integer,intent(in),value :: sample_count
        real(dp),optional,intent(in) :: colors(3,3)
        integer(c_int32_t),optional,intent(in) :: fill
        integer,optional,intent(in) :: rows(2)
        integer,parameter :: edge_start(3)=[2,3,1],edge_end(3)=[3,1,2]
        real(dp) :: pos(2,3),dx(3),dy(3),area,origin_edge(3),sample_edge(3,max_samples)
        real(dp) :: base_edge(3),edge(3),depth(max_samples),inverse_area,direction,depth_offset(max_samples)
        real(dp) :: base_depth,depth_step(2),base_rgb(3),rgb_step(3,2),sample_reach(3)
        integer :: first(2),last(2),x,y,y_first,y_last,sample,alpha,row_first,row_last,span
        integer(c_int32_t) :: fragment
        logical :: owned_edge(3),visible(max_samples),flat_color

        alpha=min(opacity,255)
        if(.not.self%ready().or.alpha<=0)return
        if(.not.all(ieee_is_finite(points)))return
        if(present(colors))then
            if(.not.all(ieee_is_finite(colors)))return
        end if
        pos(1,:)=points(1,:)-self%origin(1)
        pos(2,:)=points(2,:)-self%origin(2)
        if(self%pixel_step/=1)pos=pos/self%pixel_step
        call triangle_sample_span(self,minval(pos,dim=2),maxval(pos,dim=2),first,last)
        if(any(first>last))return
        y_first=first(2)
        y_last=last(2)
        if(present(rows))then
            y_first=max(y_first,rows(1))
            y_last=min(y_last,rows(2))
            if(y_first>y_last)return
        end if

        dx=pos(1,edge_end)-pos(1,edge_start)
        dy=pos(2,edge_end)-pos(2,edge_start)
        area=dx(1)*(pos(2,1)-pos(2,2))-dy(1)*(pos(1,1)-pos(1,2))
        if(.not.ieee_is_finite(area).or.abs(area)<=tiny(area))return
        direction=sign(1._dp,area)
        inverse_area=1._dp/abs(area)
        dx=direction*dx
        dy=direction*dy
        owned_edge=dy>0._dp.or.(dy==0._dp.and.dx<0._dp)
        origin_edge=dy*pos(1,edge_start)-dx*pos(2,edge_start)
        depth_step=[-dot_product(dy,points(3,:)),dot_product(dx,points(3,:))]*inverse_area
        flat_color=.true.
        if(present(fill))then
            fragment=fill
        else
            flat_color=all(colors(:,1)==colors(:,2)).and.all(colors(:,1)==colors(:,3))
            if(flat_color)then
                fragment=fragment_color(colors(:,1),alpha)
            else
                rgb_step(:,1)=-matmul(colors,dy)*inverse_area
                rgb_step(:,2)=matmul(colors,dx)*inverse_area
            end if
        end if
        do sample=1,sample_count
            sample_edge(:,sample)=dx*self%offsets(2,sample)-dy*self%offsets(1,sample)
            depth_offset(sample)=dot_product(sample_edge(:,sample),points(3,:))*inverse_area
        end do

        span=last(1)-first(1)
        if(span>12)sample_reach=maxval(sample_edge(:,:sample_count),dim=2)+64._dp*epsilon(1._dp)* &
            (abs(origin_edge)+abs(dx)*real(last(2),dp)+abs(dy)*real(last(1),dp))
        if(.not.present(rows))then
            self%first=min(self%first,first)
            self%last=max(self%last,last)
        end if
        do y=y_first,y_last
            base_edge=origin_edge+dx*real(y-1,dp)-dy*real(first(1)-1,dp)
            row_first=0
            row_last=span
            if(span>12)call triangle_scanline(base_edge+sample_reach,dy,span,row_first,row_last)
            if(row_first>row_last)cycle
            base_edge=base_edge-dy*real(row_first,dp)
            base_depth=dot_product(base_edge,points(3,:))*inverse_area
            if(.not.flat_color)base_rgb=matmul(colors,base_edge*inverse_area)+.5_dp*(rgb_step(:,1)+rgb_step(:,2))
            do x=first(1)+row_first,first(1)+row_last
                do sample=1,sample_count
                    depth(sample)=base_depth+depth_offset(sample)
                    visible(sample)=depth(sample)>=self%depth(sample,x,y)
                    if(.not.visible(sample))cycle
                    edge=base_edge+sample_edge(:,sample)
                    visible(sample)=all(edge>0._dp.or.(edge==0._dp.and.owned_edge))
                end do
                if(any(visible(:sample_count)))then
                    if(.not.flat_color)fragment=fragment_color(base_rgb,alpha)
                    do sample=1,sample_count
                        if(.not.visible(sample))cycle
                        if(alpha==255)then
                            self%colors(sample,x,y)=fragment
                            self%depth(sample,x,y)=depth(sample)
                        else
                            self%colors(sample,x,y)=blend_sample(fragment,self%colors(sample,x,y),alpha)
                        end if
                    end do
                end if
                base_edge=base_edge-dy
                base_depth=base_depth+depth_step(1)
                if(.not.flat_color)base_rgb=base_rgb+rgb_step(:,1)
            end do
        end do
    end subroutine raster_triangle

    pure subroutine triangle_sample_span(self,lower,upper,first,last)
        class(viewport_depth_t),intent(in) :: self
        real(dp),intent(in) :: lower(2),upper(2)
        integer,intent(out) :: first(2),last(2)
        real(dp) :: lo(2),hi(2)

        lo=lower-self%offsets(:,self%samples)-self%sample_margin
        hi=upper-self%offsets(:,1)+self%sample_margin
        first=1
        last=0
        if(any(hi<0._dp).or.any(lo>real(self%extent-1,dp)))return
        first=ceiling(max(lo,0._dp))+1
        last=floor(min(hi,real(self%extent-1,dp)))+1
    end subroutine triangle_sample_span

    pure subroutine triangle_scanline(edges,slopes,span,first,last)
        real(dp),intent(in) :: edges(3),slopes(3)
        integer,intent(in) :: span
        integer,intent(out) :: first,last
        real(dp) :: far,crossing
        integer :: side

        first=0
        last=span
        do side=1,3
            far=edges(side)-slopes(side)*real(span,dp)
            if(.not.ieee_is_finite(edges(side)).or..not.ieee_is_finite(far))cycle
            if(max(edges(side),far)<0._dp)then
                last=-1
                return
            end if
            if(min(edges(side),far)>=0._dp)cycle
            crossing=max(0._dp,min(real(span,dp),edges(side)/slopes(side)))
            if(slopes(side)>0._dp)then
                last=min(last,ceiling(crossing))
            else
                first=max(first,floor(crossing))
            end if
        end do
    end subroutine triangle_scanline

    subroutine line(self,first,last,rgb,opacity,width,bias)
        class(viewport_depth_t),intent(inout) :: self
        real(dp),intent(in) :: first(3),last(3),rgb(3),width
        integer,intent(in) :: opacity
        real(dp),intent(in),optional :: bias
        real(dp) :: depth_offset
        integer :: alpha

        alpha=min(opacity,255)
        if(.not.self%ready().or.alpha<=0.or.width<=0._dp)return
        if(.not.ieee_is_finite(width).or..not.all(ieee_is_finite(rgb)))return
        depth_offset=0._dp
        if(present(bias))depth_offset=bias
        if(.not.ieee_is_finite(depth_offset))return
        call raster_line(self,first,last,alpha,fragment_color(rgb,alpha),stroke_radius(self,width),depth_offset, &
            alpha==255.and..not.present(bias))
    end subroutine line

    subroutine grid_edges(self,points,edges,strides,rgb,opacities,width,bias)
        class(viewport_depth_t),intent(inout) :: self
        real(dp),contiguous,intent(in) :: points(:,:)
        real(dp),intent(in) :: rgb(3),width,bias
        integer(int8),contiguous,intent(in) :: edges(:,:)
        integer,intent(in) :: strides(:),opacities(:)
        integer(c_int32_t) :: fragments(size(opacities))
        real(dp) :: radius
        integer :: i,node,axis,weight

        if(.not.self%ready().or.width<=0._dp)return
        if(.not.ieee_is_finite(width).or..not.ieee_is_finite(bias).or..not.all(ieee_is_finite(rgb)))return
        radius=stroke_radius(self,width)
        do i=1,size(opacities)
            fragments(i)=fragment_color(rgb,opacities(i))
        end do
        do node=1,size(points,2)
            do axis=1,size(strides)
                weight=int(edges(axis,node))
                if(weight==0)cycle
                call raster_line(self,points(:,node),points(:,node+strides(axis)), &
                    opacities(weight),fragments(weight),radius,bias,.false.)
            end do
        end do
    end subroutine grid_edges

    pure real(dp) function stroke_radius(self,width) result(radius)
        class(viewport_depth_t),intent(in) :: self
        real(dp),intent(in) :: width

        radius=.5_dp*width
        if(self%pixel_step>1)radius=max(.5_dp,radius/self%pixel_step)
    end function stroke_radius

    subroutine raster_line(self,first,last,alpha,fragment,radius,depth_offset,store_depth)
        class(viewport_depth_t),intent(inout) :: self
        real(dp),intent(in) :: first(3),last(3),radius,depth_offset
        integer,intent(in) :: alpha
        integer(c_int32_t),intent(in) :: fragment
        logical,intent(in) :: store_depth
        real(dp) :: start(2),finish(2),delta(2),v(2),length_squared,pixel_fraction,depth
        real(dp) :: t_start,t_end,low_t,high_t,lower_x,upper_x,first_sample,last_sample,margin(2),radius_squared
        integer :: x,y,sample,lo,hi,row_lo,row_hi
        logical :: full_detail,short_segment

        if(.not.all(ieee_is_finite(first)).or..not.all(ieee_is_finite(last)))return
        start=first(:2)-real(self%origin,dp)
        finish=last(:2)-real(self%origin,dp)
        if(self%pixel_step>1)then
            start=start/self%pixel_step
            finish=finish/self%pixel_step
        end if
        delta=finish-start
        length_squared=sum(delta**2)
        if(.not.ieee_is_finite(length_squared))return
        radius_squared=radius**2
        first_sample=self%offsets(1,1)
        last_sample=self%offsets(1,self%samples)
        full_detail=self%pixel_step==1
        margin=0._dp
        if(full_detail)margin=(abs(start)+abs(finish)+radius+1._dp)*(16._dp*epsilon(1._dp))
        call line_sample_span(min(start(2),finish(2))-radius,max(start(2),finish(2))+radius, &
            first_sample,last_sample,self%extent(2),margin(2),full_detail,lo,hi)

        short_segment=full_detail.and.maxval(abs(delta))<=2._dp
        if(short_segment)then
            call line_sample_span(min(start(1),finish(1))-radius,max(start(1),finish(1))+radius, &
                first_sample,last_sample,self%extent(1),margin(1),full_detail,row_lo,row_hi)
        end if
        do y=lo,hi
            if(.not.short_segment)then
                low_t=0._dp
                high_t=1._dp
                if(abs(delta(2))>tiny(1._dp))then
                    t_start=(real(y-1,dp)+first_sample-radius-margin(2)-start(2))/delta(2)
                    t_end=(real(y-1,dp)+last_sample+radius+margin(2)-start(2))/delta(2)
                    low_t=max(0._dp,min(t_start,t_end))
                    high_t=min(1._dp,max(t_start,t_end))
                    if(low_t>high_t)cycle
                end if
                lower_x=start(1)+min(low_t*delta(1),high_t*delta(1))-radius
                upper_x=start(1)+max(low_t*delta(1),high_t*delta(1))+radius
                call line_sample_span(lower_x,upper_x,first_sample,last_sample,self%extent(1),margin(1),full_detail,row_lo,row_hi)
            end if
            if(row_lo>row_hi)cycle
            self%first=min(self%first,[row_lo,y])
            self%last=max(self%last,[row_hi,y])
            do x=row_lo,row_hi
                do sample=1,self%samples
                    v=[real(x-1,dp),real(y-1,dp)]+self%offsets(:,sample)-start
                    pixel_fraction=0._dp
                    if(length_squared>tiny(1._dp))pixel_fraction=max(0._dp,min(1._dp,dot_product(v,delta)/length_squared))
                    if(sum((v-pixel_fraction*delta)**2)>radius_squared)cycle
                    depth=(1._dp-pixel_fraction)*first(3)+pixel_fraction*last(3)+depth_offset
                    if(depth<self%depth(sample,x,y))cycle
                    if(alpha==255)then
                        self%colors(sample,x,y)=fragment
                    else
                        self%colors(sample,x,y)=blend_sample(fragment,self%colors(sample,x,y),alpha)
                    end if
                    if(store_depth)self%depth(sample,x,y)=depth
                end do
            end do
        end do
    end subroutine raster_line

    pure subroutine line_sample_span(lower,upper,first_sample,last_sample,extent,margin,full_detail,first,last)
        real(dp),intent(in) :: lower,upper,first_sample,last_sample,margin
        integer,intent(in) :: extent
        logical,intent(in) :: full_detail
        integer,intent(out) :: first,last
        real(dp) :: lo,hi

        first=1
        last=0
        if(full_detail)then
            lo=lower-last_sample-margin
            hi=upper-first_sample+margin
            if(hi<0._dp.or.lo>real(extent-1,dp))return
            first=ceiling(max(lo,0._dp))+1
        else
            if(upper<0._dp.or.lower>=real(extent,dp))return
            first=floor(max(lower,0._dp))+1
            hi=upper
        end if
        last=floor(min(hi,real(extent-1,dp)))+1
    end subroutine line_sample_span

    pure integer(c_int32_t) function packed_rgba(rgb,alpha) result(color)
        integer,intent(in) :: rgb(3),alpha

        color=ior(shiftl(int(alpha,c_int32_t),24),ior(shiftl(int(rgb(1),c_int32_t),16), &
            ior(shiftl(int(rgb(2),c_int32_t),8),int(rgb(3),c_int32_t))))
    end function packed_rgba

    pure integer(c_int32_t) function fragment_color(rgb,alpha) result(color)
        real(dp),intent(in) :: rgb(3)
        integer,intent(in) :: alpha
        integer :: channels(3)
        real(dp) :: scaled(3)

        scaled=max(0._dp,min(255._dp,rgb))
        channels=int(scaled)
        where(scaled-real(channels,dp)>=.5_dp)channels=channels+1
        if(alpha<255)channels=(channels*alpha+127)/255
        color=packed_rgba(channels,alpha)
    end function fragment_color

    pure integer(c_int32_t) function blend_sample(front,back,alpha) result(color)
        integer(c_int32_t),value :: front,back
        integer,value :: alpha
        integer :: channels(3),a

        channels=[ibits(front,16,8),ibits(front,8,8),ibits(front,0,8)]+ &
            ([ibits(back,16,8),ibits(back,8,8),ibits(back,0,8)]*(255-alpha)+127)/255
        a=alpha+(ibits(back,24,8)*(255-alpha)+127)/255
        color=packed_rgba(channels,a)
    end function blend_sample

    real(dp) function depth_at(self,x,y) result(depth)
        class(viewport_depth_t),intent(in) :: self
        real(dp),intent(in) :: x,y
        real(dp) :: local(2)
        integer :: pixel(2)

        depth=empty_depth
        if(.not.self%ready())return
        local=([x,y]-real(self%origin,dp))/self%pixel_step
        if(.not.all(ieee_is_finite(local)))return
        if(any(local<0._dp).or.any(local>=real(self%extent,dp)))return
        pixel=floor(local)+1
        depth=maxval(self%depth(:,pixel(1),pixel(2)))
    end function depth_at

    subroutine paint(self,painter)
        class(viewport_depth_t),intent(in) :: self
        type(QPainter),intent(inout) :: painter
        type(QImage) :: image
        type(QColor) :: color
        integer(c_int32_t) :: pixel
        integer(c_int32_t) :: pixels(max_samples)
        integer :: x,y,alpha,rgb(3),extent(2),origin(2)

        extent=self%last-self%first+1
        if(.not.self%ready().or.any(extent<1))return
        image=QImage(extent(1),extent(2))
        if(image%isNull())then
            call image%delete()
            return
        end if
        call image%fill(QColor(0,0,0,0))
        do y=self%first(2),self%last(2)
            do x=self%first(1),self%last(1)
                pixel=self%colors(1,x,y)
                if(all(self%colors(:,x,y)==pixel))then
                    alpha=ibits(pixel,24,8)
                    if(alpha==0)cycle
                    rgb=[ibits(pixel,16,8),ibits(pixel,8,8),ibits(pixel,0,8)]
                    if(alpha<255)rgb=(255*rgb+alpha/2)/alpha
                else
                    pixels(:self%samples)=self%colors(:,x,y)
                    alpha=sum(ibits(pixels(:self%samples),24,8))
                    rgb=[sum(ibits(pixels(:self%samples),16,8)),sum(ibits(pixels(:self%samples),8,8)), &
                        sum(ibits(pixels(:self%samples),0,8))]
                    rgb=(255*rgb+alpha/2)/alpha
                    alpha=(alpha+self%samples/2)/self%samples
                end if
                color%r=rgb(1)
                color%g=rgb(2)
                color%b=rgb(3)
                color%a=alpha
                call image%setPixelColor(x-self%first(1),y-self%first(2),color)
            end do
        end do
        origin=self%origin+(self%first-1)*self%pixel_step
        if(self%pixel_step==1)then
            call painter%drawImage(origin(1),origin(2),image)
        else
            extent=extent*self%pixel_step
            call painter%save()
            call painter%setRenderHint(QPainter_SmoothPixmapTransform)
            call painter%drawImage(origin(1),origin(2),extent(1),extent(2),image)
            call painter%restore()
        end if
        call image%delete()
    end subroutine paint

end module studio_viewport_depth
