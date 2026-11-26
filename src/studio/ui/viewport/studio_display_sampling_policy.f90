module studio_display_sampling_policy
    use, intrinsic :: iso_c_binding, only: dp=>c_double
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none(type,external)
    private
    public :: display_sampling_error,display_sample_counts

    real(dp),parameter :: target_error=.35_dp,refine_error=.44_dp,coarsen_error=.24_dp
    integer,parameter :: curve_limit=1025,face_limit=257,display_budget=32768

contains

    pure function display_sampling_error(points,counts,rank) result(error)
        real(dp),intent(in) :: points(:,:)
        integer,intent(in) :: counts(3),rank
        real(dp) :: error(3)
        integer :: size_in_points,axis,face,first,second,stride(3)

        error=0._dp
        if(rank<1.or.rank>3.or.size(points,1)/=3.or.any(counts<1))return
        size_in_points=1
        do axis=1,3
            if(counts(axis)>size(points,2)/size_in_points)return
            size_in_points=size_in_points*counts(axis)
        end do
        if(size_in_points/=size(points,2).or.any(counts(rank+1:)/=1))return

        stride=[1,counts(1),counts(1)*counts(2)]
        select case(rank)
        case(1)
            call measure_line(points,1,counts(1),1,error(1))
        case(2)
            call measure_face(points,1,counts(:2),stride(:2),error(:2))
        case default
            do face=1,3
                first=mod(face,3)+1
                second=mod(face+1,3)+1
                call measure_volume_side(1,first,second)
                if(counts(face)>1)call measure_volume_side(1+(counts(face)-1)*stride(face),first,second)
            end do
        end select
    contains

        pure subroutine measure_volume_side(start,first,second)
            integer,intent(in) :: start,first,second
            real(dp) :: side_error(2)

            side_error=error([first,second])
            call measure_face(points,start,counts([first,second]),stride([first,second]),side_error)
            error([first,second])=side_error
        end subroutine measure_volume_side

    end function display_sampling_error

    pure subroutine measure_line(points,start,count,stride,error)
        real(dp),intent(in) :: points(:,:)
        integer,intent(in) :: start,count,stride
        real(dp),intent(inout) :: error
        real(dp) :: chord(3),offset(3),direction(3),length,depth,deviation
        integer :: i

        if(count<3)return
        do i=start,start+(count-3)*stride,stride
            chord=points(:,i+2*stride)-points(:,i)
            offset=points(:,i+stride)-points(:,i)
            length=norm2(chord)
            if(.not.ieee_is_finite(length))cycle
            if(length>tiny(1._dp))then
                direction=chord/length
                depth=max(0._dp,min(length,dot_product(direction,offset)))
                offset=offset-depth*direction
            end if
            deviation=.25_dp*norm2(offset)
            if(ieee_is_finite(deviation))error=max(error,deviation)
        end do
    end subroutine measure_line

    pure subroutine measure_face(points,start,counts,stride,error)
        real(dp),intent(in) :: points(:,:)
        integer,intent(in) :: start,counts(2),stride(2)
        real(dp),intent(inout) :: error(2)
        real(dp) :: tangent(3),diagonal(3),normal(3),offset(3),length,deviation
        integer :: i,j,node

        do j=0,counts(2)-1
            call measure_line(points,start+j*stride(2),counts(1),stride(1),error(1))
        end do
        do i=0,counts(1)-1
            call measure_line(points,start+i*stride(1),counts(2),stride(2),error(2))
        end do
        do j=0,counts(2)-2
            do i=0,counts(1)-2
                node=start+i*stride(1)+j*stride(2)
                tangent=unit_vector(points(:,node+stride(1))-points(:,node))
                diagonal=unit_vector(points(:,node+sum(stride))-points(:,node))
                normal=[tangent(2)*diagonal(3)-tangent(3)*diagonal(2), &
                    tangent(3)*diagonal(1)-tangent(1)*diagonal(3),tangent(1)*diagonal(2)-tangent(2)*diagonal(1)]
                length=norm2(normal)
                if(length<=tiny(1._dp))cycle
                offset=points(:,node+stride(2))-points(:,node)
                deviation=.25_dp*abs(dot_product(offset,normal/length))
                if(ieee_is_finite(deviation))error=max(error,deviation)
            end do
        end do
    end subroutine measure_face

    pure function unit_vector(vector) result(direction)
        real(dp),intent(in) :: vector(3)
        real(dp) :: direction(3),length

        direction=0._dp
        length=norm2(vector)
        if(ieee_is_finite(length).and.length>tiny(1._dp))direction=vector/length
    end function unit_vector

    pure function display_sample_counts(base,rank,error,scale,previous,max_points) result(counts)
        integer,intent(in) :: base(3),rank
        real(dp),intent(in) :: error(3),scale
        integer,intent(in),optional :: previous(3),max_points
        integer :: counts(3),limit,axis,candidate(3),best(3)
        real(dp) :: budget,pixels(3),threshold,worst,deviation

        counts=base
        if(rank<1.or.rank>3.or.any(base<1))return
        if(.not.ieee_is_finite(scale).or.scale<=0._dp)return
        limit=merge(curve_limit,face_limit,rank==1)
        budget=real(display_budget,dp)
        if(present(max_points))budget=real(max(1,max_points),dp)
        budget=max(budget,display_point_count(base,rank))
        pixels=0._dp
        do axis=1,rank
            if(base(axis)<2.or..not.ieee_is_finite(error(axis)))cycle
            pixels(axis)=max(0._dp,min(error(axis),huge(1._dp)/max(1._dp,scale)))*scale
        end do

        if(present(previous))then
            do axis=1,rank
                candidate=counts
                do while(candidate(axis)<limit.and.candidate(axis)<previous(axis))
                    if(estimated_error(pixels(axis),base(axis),candidate(axis))<coarsen_error)exit
                    candidate(axis)=min(limit,2*candidate(axis)-1)
                end do
                if(display_point_count(candidate,rank)<=budget)counts=candidate
            end do
        end if

        threshold=target_error
        if(present(previous))threshold=refine_error
        do
            best=counts
            worst=threshold
            do axis=1,rank
                if(counts(axis)<2.or.counts(axis)>=limit)cycle
                deviation=estimated_error(pixels(axis),base(axis),counts(axis))
                if(deviation<=worst)cycle
                candidate=counts
                candidate(axis)=min(limit,2*candidate(axis)-1)
                if(display_point_count(candidate,rank)>budget)cycle
                worst=deviation
                best=candidate
            end do
            if(all(best==counts))exit
            counts=best
        end do
    end function display_sample_counts

    pure real(dp) function estimated_error(pixels,base,count) result(error)
        real(dp),intent(in) :: pixels
        integer,intent(in) :: base,count

        error=pixels*(real(base-1,dp)/real(max(1,count-1),dp))**2
    end function estimated_error

    pure real(dp) function display_point_count(counts,rank) result(points)
        integer,intent(in) :: counts(3),rank
        real(dp) :: u,v,w

        u=real(counts(1),dp); v=real(counts(2),dp); w=real(counts(3),dp)
        select case(rank)
        case(1); points=u
        case(2); points=u*v
        case default; points=2._dp*(u*v+u*w+v*w)
        end select
    end function display_point_count

end module studio_display_sampling_policy
