module studio_forcad_patch_data
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: int64
    use forcad, only: rk, nurbs_curve, nurbs_surface, nurbs_volume
    use forcad_utils, only: linspace, ndgrid
    implicit none(type, external)
    private

    type, public :: native_knot_axis
        real(rk), allocatable :: knots(:)
    end type native_knot_axis

    public :: describe_patch, copy_native_patch, get_native_controls, sample_patch_geometry, sample_volume_slice
    public :: sample_display_grid

contains

    logical function get_native_controls(geometry,count,points,weights) result(valid)
        class(*), intent(in) :: geometry
        integer, intent(in) :: count
        real(rk), allocatable, intent(out) :: points(:,:),weights(:)
        valid=.false.
        select type (geometry)
        type is(nurbs_curve)
            points=geometry%get_Xc()
            weights=geometry%get_Wc()
        type is(nurbs_surface)
            points=geometry%get_Xc()
            weights=geometry%get_Wc()
        type is(nurbs_volume)
            points=geometry%get_Xc()
            weights=geometry%get_Wc()
        class default
            return
        end select
        if (any(shape(points)/=[count,3]).or..not.all(ieee_is_finite(points))) return
        if (size(weights)==0) weights=spread(1._rk,1,count)
        valid=size(weights)==count.and.all(ieee_is_finite(weights)).and.all(weights>0._rk)
    end function get_native_controls

    logical function describe_patch(geometry,dimension,degrees,control_counts,axes,bounds) result(valid)
        class(*), intent(in) :: geometry
        integer, intent(in) :: dimension
        type(native_knot_axis), intent(inout) :: axes(3)
        integer, intent(inout) :: degrees(3),control_counts(3)
        real(rk), intent(inout) :: bounds(2,3)
        real(rk), allocatable :: knots(:)
        integer :: direction,degree,controls

        valid=.false.
        select type (geometry)
        type is(nurbs_curve)
            if (.not.geometry%err%ok.or..not.geometry%is_valid()) return
            degrees(1)=geometry%get_degree()
            control_counts(1)=geometry%get_nc()
            axes(1)%knots=geometry%get_knot()
        type is(nurbs_surface)
            if (.not.geometry%err%ok.or..not.geometry%is_valid()) return
            degrees(:2)=geometry%get_degree()
            control_counts(:2)=geometry%get_nc()
            do direction=1,2
                axes(direction)%knots=geometry%get_knot(direction)
            end do
        type is(nurbs_volume)
            if (.not.geometry%err%ok.or..not.geometry%is_valid()) return
            degrees=geometry%get_degree()
            control_counts=geometry%get_nc()
            do direction=1,3
                axes(direction)%knots=geometry%get_knot(direction)
            end do
        class default
            return
        end select

        do direction=1,dimension
            knots=axes(direction)%knots
            degree=degrees(direction)
            controls=control_counts(direction)
            if (degree<0.or.controls<=degree) return
            if (size(knots)/=controls+degree+1.or..not.all(ieee_is_finite(knots))) return
            ! Full stored knot vectors can extend beyond the active interval.
            bounds(:,direction)=[knots(degree+1),knots(controls+1)]
            if (bounds(2,direction)<=bounds(1,direction)) return
        end do
        valid=.true.
    end function describe_patch

    logical function sample_patch_geometry(geometry,counts,points,controls) result(valid)
        class(*), intent(inout) :: geometry
        integer, intent(in) :: counts(3)
        real(rk), allocatable, intent(out) :: points(:,:),controls(:,:)

        valid=.false.
        select type (geometry)
        type is(nurbs_curve)
            call geometry%create(res=counts(1))
            if (.not.geometry%err%ok) return
            points=geometry%get_Xg()
            controls=geometry%get_Xc()
        type is(nurbs_surface)
            call geometry%create(res1=counts(1),res2=counts(2))
            if (.not.geometry%err%ok) return
            points=geometry%get_Xg()
            controls=geometry%get_Xc()
        type is(nurbs_volume)
            call geometry%create(res1=counts(1),res2=counts(2),res3=counts(3))
            if (.not.geometry%err%ok) return
            points=geometry%get_Xg()
            controls=geometry%get_Xc()
        class default
            return
        end select
        valid=.true.
    end function sample_patch_geometry

    recursive logical function sample_display_grid(geometry,counts,bounds,side,points) result(valid)
        class(*), intent(in) :: geometry
        integer, intent(in) :: counts(3),side
        real(rk), intent(in) :: bounds(2,3)
        real(rk), allocatable, intent(out) :: points(:,:)
        type(nurbs_surface) :: boundary
        real(rk), allocatable :: rows(:,:),parameters(:,:),first_values(:),second_values(:)
        real(rk) :: active_bounds(2,3)
        integer :: rank,normal,endpoint,first,second,width,height,point_count,coordinates,status,i,j,direction

        valid=.false.
        select type(geometry)
        type is(nurbs_curve)
            rank=1
            active_bounds(:,1)=[geometry%get_knot(geometry%get_degree()+1),geometry%get_knot(geometry%get_nc()+1)]
        type is(nurbs_surface)
            rank=2
            do direction=1,rank
                active_bounds(:,direction)=[geometry%get_knot(direction,geometry%get_degree(direction)+1), &
                    geometry%get_knot(direction,geometry%get_nc(direction)+1)]
            end do
        type is(nurbs_volume)
            rank=3
            do direction=1,rank
                active_bounds(:,direction)=[geometry%get_knot(direction,geometry%get_degree(direction)+1), &
                    geometry%get_knot(direction,geometry%get_nc(direction)+1)]
            end do
        class default
            return
        end select

        if(side<1.or.side>merge(6,1,rank==3))return
        if(.not.all(ieee_is_finite(bounds(:,:rank))))return
        if(any(bounds(2,:rank)<=bounds(1,:rank)))return
        if(any(bounds(1,:rank)<active_bounds(1,:rank)).or.any(bounds(2,:rank)>active_bounds(2,:rank)))return
        if(any(counts(:rank)<1))return
        normal=3
        if(rank==3)normal=(side+1)/2
        first=mod(normal,3)+1
        second=mod(normal+1,3)+1
        width=counts(first)
        height=1
        if(rank>1)height=counts(second)
        if(int(width,int64)*height>int(huge(point_count),int64)/3)return
        point_count=width*height

        first_values=linspace(bounds(1,first),bounds(2,first),width)
        if(width>1)first_values(width)=bounds(2,first)
        if(rank>1)then
            second_values=linspace(bounds(1,second),bounds(2,second),height)
            if(height>1)second_values(height)=bounds(2,second)
        end if
        select type(geometry)
        type is(nurbs_curve)
            rows=geometry%cmp_Xg(first_values)
        type is(nurbs_surface)
            call ndgrid(first_values,second_values,parameters)
            rows=geometry%cmp_Xg(parameters)
        type is(nurbs_volume)
            endpoint=1+mod(side+1,2)
            boundary=geometry%slice(normal,bounds(endpoint,normal))
            if(.not.boundary%err%ok)return
            if(normal==2)then
                call boundary%create(Xt1=second_values,Xt2=first_values)
            else
                call boundary%create(Xt1=first_values,Xt2=second_values)
            end if
            if(.not.boundary%err%ok)return
            rows=boundary%get_Xg()
        end select

        if(.not.allocated(rows))return
        coordinates=size(rows,2)
        if(size(rows,1)/=point_count.or.coordinates<1.or.coordinates>3)return
        if(.not.all(ieee_is_finite(rows)))return
        allocate(points(3,point_count),source=0._rk,stat=status)
        if(status/=0)return
        if(normal==2)then
            do j=1,height
                do i=1,width
                    points(:coordinates,i+(j-1)*width)=rows(j+(i-1)*height,:)
                end do
            end do
        else
            points(:coordinates,:)=transpose(rows)
        end if
        valid=.true.
    end function sample_display_grid

    logical function sample_volume_slice(geometry,parameter,counts,points) result(valid)
        class(*), intent(in) :: geometry
        real(rk), intent(in) :: parameter
        integer, intent(in) :: counts(2)
        real(rk), allocatable, intent(out) :: points(:,:)
        type(nurbs_surface) :: slice

        valid=.false.
        select type (geometry)
        type is(nurbs_volume)
            slice=geometry%slice(dir=3,Xt=parameter)
        class default
            return
        end select
        if (.not.slice%err%ok.or..not.slice%is_valid()) return
        call slice%create(res1=counts(1),res2=counts(2))
        if (.not.slice%err%ok) return
        points=slice%get_Xg()
        valid=.true.
    end function sample_volume_slice

    subroutine copy_native_patch(geometry,copy)
        class(*), intent(in) :: geometry
        class(*), allocatable, intent(out) :: copy
        type(nurbs_curve) :: curve
        type(nurbs_surface) :: surface
        real(rk), allocatable :: points(:,:),weights(:)
        integer :: ncoord

        select type (geometry)
        type is (nurbs_curve)
            ncoord=size(geometry%get_Xc(1))
            if(ncoord>=1.and.ncoord<3.and.geometry%is_valid())then
                points=xyz_controls(geometry%get_Xc())
                weights=geometry%get_Wc()
                if(size(weights)==0)deallocate(weights)
                call curve%set(knot=geometry%get_knot(),Xc=points,Wc=weights, &
                    degree=geometry%get_degree(),wrap_parameters=geometry%get_parameter_wrapping())
                allocate(copy,source=curve)
            else
                allocate(copy,source=geometry)
            end if
        type is (nurbs_surface)
            ncoord=size(geometry%get_Xc(1))
            if(ncoord>=1.and.ncoord<3.and.geometry%is_valid())then
                points=xyz_controls(geometry%get_Xc())
                weights=geometry%get_Wc()
                if(size(weights)==0)deallocate(weights)
                call surface%set(knot1=geometry%get_knot(1),knot2=geometry%get_knot(2),Xc=points,Wc=weights, &
                    degree=geometry%get_degree(),wrap_parameters=geometry%get_parameter_wrapping([1,2]))
                allocate(copy,source=surface)
            else
                allocate(copy,source=geometry)
            end if
        type is (nurbs_volume)
            allocate(copy,source=geometry)
        end select
    end subroutine copy_native_patch

    function xyz_controls(points) result(xyz)
        real(rk), intent(in) :: points(:,:)
        real(rk) :: xyz(size(points,1),3)

        xyz=0._rk
        xyz(:,:size(points,2))=points
    end function xyz_controls
end module studio_forcad_patch_data
