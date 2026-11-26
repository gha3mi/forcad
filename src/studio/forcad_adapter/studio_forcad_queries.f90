module studio_forcad_queries
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use forcad, only: rk, nurbs_curve, nurbs_surface, nurbs_volume
    use forcad_utils, only: active_span_count
    implicit none(type, external)
    private

    public :: closest_point_on_patch, physical_knot_grid, sample_physical_knots

    integer, parameter :: preferred_samples=8
    integer(int64), parameter :: point_budget=60000_int64

    type :: physical_knot_grid
        real(rk), allocatable :: points(:,:) ! [coordinate, point]
        integer, allocatable :: lines(:,:)   ! Native [polyline, point on line]
        integer, allocatable :: faces(:,:)   ! For volumes, exterior faces incident to each polyline
    end type physical_knot_grid

contains

    logical function closest_point_on_patch(geometry,face,bounds,normalized,query,point) result(found)
        class(*), intent(inout) :: geometry
        integer, intent(in) :: face
        real(c_double), intent(in) :: bounds(2,3),normalized(3),query(3)
        real(c_double), intent(out) :: point(3)
        type(nurbs_surface) :: boundary
        real(rk) :: seed(3),parameters(3),native_point(3),target(3)
        real(rk), parameter :: tolerance=max(1.e-10_rk,epsilon(1._rk)*16._rk)
        integer, parameter :: iterations=32
        integer :: axis,end_point,retained(2),ncoord

        found=.false.
        native_point=0._rk
        target=real(query,rk)
        seed=real(bounds(1,:)*(1._c_double-normalized)+bounds(2,:)*normalized,rk)
        seed=max(real(bounds(1,:),rk),min(real(bounds(2,:),rk),seed))
        select type(geometry)
        type is(nurbs_curve)
            ncoord=size(geometry%get_Xc(1))
            if(ncoord<1.or.ncoord>size(native_point))return
            call geometry%nearest_point2(target,tolerance,iterations,parameters(1),native_point(:ncoord),seed(1))
            if(.not.geometry%err%ok)return
        type is(nurbs_surface)
            ncoord=size(geometry%get_Xc(1))
            if(ncoord<1.or.ncoord>size(native_point))return
            call geometry%nearest_point2(target,tolerance,iterations,parameters(:2),native_point(:ncoord),seed(:2))
            if(.not.geometry%err%ok)return
        type is(nurbs_volume)
            ! Snap to the visible exterior. Projecting onto the volume interior may leave the cursor in free space.
            if(face<1.or.face>6)return
            axis=(face+1)/2
            end_point=2-mod(face,2)
            retained=pack([1,2,3],[1,2,3]/=axis)
            boundary=geometry%slice(axis,real(bounds(end_point,axis),rk))
            if(.not.boundary%err%ok)return
            ncoord=size(boundary%get_Xc(1))
            if(ncoord<1.or.ncoord>size(native_point))return
            call boundary%nearest_point2(target,tolerance,iterations,parameters(:2),native_point(:ncoord),seed(retained))
            if(.not.boundary%err%ok)return
        class default
            return
        end select
        point=real(native_point,c_double)
        found=all(ieee_is_finite(point))
    end function closest_point_on_patch

    subroutine sample_physical_knots(native,grid)
        class(*),intent(inout) :: native
        type(physical_knot_grid),intent(out) :: grid
        real(rk),allocatable :: points(:,:)
        integer,allocatable :: lines(:,:)
        integer :: spans(3),axis,samples,ncoords
        logical :: volume

        volume=.false.
        select type(native)
        type is(nurbs_curve)
            spans(1)=active_span_count(native%get_knot(),native%get_nc(),native%get_degree())
            if(spans(1)<1.or.int(spans(1),int64)>=point_budget)return
            ! Two native samples per span return exactly the distinct knot locations.
            call native%cmp_Xth_in_Xg(points,lines,res=2)
        type is(nurbs_surface)
            do axis=1,2
                spans(axis)=active_span_count(native%get_knot(axis),native%get_nc(axis),native%get_degree(axis))
            end do
            samples=bounded_sampling(spans(:2))
            if(samples<2)return
            call native%cmp_Xth_in_Xg(points,lines,res=samples)
        type is(nurbs_volume)
            volume=.true.
            do axis=1,3
                spans(axis)=active_span_count(native%get_knot(axis),native%get_nc(axis),native%get_degree(axis))
            end do
            samples=bounded_sampling(spans)
            if(samples<2)return
            call native%cmp_Xth_in_Xg(points,lines,res=samples)
        class default
            return
        end select
        if(.not.allocated(points).or..not.allocated(lines))return
        if(size(points,1)==0.or.size(lines,1)==0.or.size(lines,2)<2)return
        ncoords=size(points,2)
        if(ncoords<2.or.ncoords>3)return
        if(any(lines<1).or.any(lines>size(points,1)).or..not.all(ieee_is_finite(points)))return
        if(volume)then
            if(size(lines,1)/=sum(product(spans+1)/(spans+1)))return
            call volume_line_faces(spans+1,grid%faces)
        end if
        allocate(grid%points(3,size(points,1)),source=0._rk)
        grid%points(:ncoords,:)=transpose(points)
        call move_alloc(lines,grid%lines)
    end subroutine sample_physical_knots

    integer function bounded_sampling(spans) result(samples)
        integer,intent(in) :: spans(:)
        integer(int64) :: axis_counts(size(spans)),max_intervals,common,a,b,remainder,factor
        integer :: axis

        ! From the native polyline layout: all lines have one common subdivision
        ! per span ratio. Restrict the allocation size before asking ForCAD to sample.
        samples=0
        if(any(spans<1).or.any(int(spans,int64)>=point_budget))return
        axis_counts=int(spans,int64)+1_int64
        max_intervals=point_budget/sum(product(axis_counts)/axis_counts)-1_int64
        if(max_intervals<1_int64)return
        common=1_int64
        do axis=1,size(spans)
            a=common
            b=int(spans(axis),int64)
            do while(b/=0_int64)
                remainder=mod(a,b)
                a=b
                b=remainder
            end do
            factor=int(spans(axis),int64)/a
            if(common>max_intervals/factor)return
            common=common*factor
        end do
        samples=min(preferred_samples,1+int(max_intervals/common))
    end function bounded_sampling

    subroutine volume_line_faces(counts,faces)
        integer,intent(in) :: counts(3)
        integer,allocatable,intent(out) :: faces(:,:)
        integer :: varying,line,row,quotient,fixed,axis,at

        ! Native line families are ordered by varying axis, then by their two
        ! fixed knot indices with the lowest fixed direction varying fastest.
        allocate(faces(2,sum(product(counts)/counts)),source=0)
        row=0
        do varying=1,3
            do line=0,product(counts)/counts(varying)-1
                row=row+1
                fixed=0
                quotient=line
                do axis=1,3
                    if(axis==varying)cycle
                    fixed=fixed+1
                    at=mod(quotient,counts(axis))+1
                    if(at==1)faces(fixed,row)=2*axis-1
                    if(at==counts(axis))faces(fixed,row)=2*axis
                    quotient=quotient/counts(axis)
                end do
            end do
        end do
    end subroutine volume_line_faces
end module studio_forcad_queries
