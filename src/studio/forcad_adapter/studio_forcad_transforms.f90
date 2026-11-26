module studio_forcad_transforms
    use, intrinsic :: iso_c_binding, only: c_double
    use forcad, only: rk, rotation, nurbs_curve, nurbs_surface, nurbs_volume, &
        nurbs_multipatch_curve, nurbs_multipatch_surface, nurbs_multipatch_volume
    implicit none(type, external)
    private

    public :: affine_transform_t, studio_native_patch, native_patch_kind, rotation_matrix, transform_patch, transform_assembly

    real(c_double), parameter :: identity_matrix(3,3)=reshape([1._c_double,0._c_double,0._c_double, &
        0._c_double,1._c_double,0._c_double,0._c_double,0._c_double,1._c_double],[3,3])

    type :: affine_transform_t
        real(c_double) :: linear(3,3)=identity_matrix
        real(c_double) :: offset(3)=0._c_double
        real(c_double) :: pivot(3)=0._c_double
    end type affine_transform_t

    type :: studio_native_patch
        class(*), allocatable :: geometry
    end type studio_native_patch

contains

    integer function native_patch_kind(geometry) result(kind)
        class(*), intent(in) :: geometry
        kind=0
        select type(geometry)
        type is(nurbs_curve)
            kind=1
        type is(nurbs_surface)
            kind=2
        type is(nurbs_volume)
            kind=3
        end select
    end function native_patch_kind

    pure function rotation_matrix(angles) result(matrix)
        real(c_double), intent(in) :: angles(3)
        real(c_double) :: matrix(3,3)

        matrix=real(rotation(real(angles(1),rk),real(angles(2),rk),real(angles(3),rk)),c_double)
    end function rotation_matrix

    logical function transform_patch(geometry,transform) result(valid)
        class(*), intent(inout) :: geometry
        type(affine_transform_t), intent(in) :: transform
        real(rk) :: linear(3,3),offset(3),pivot(3)

        linear=real(transform%linear,rk)
        offset=real(transform%offset,rk)
        pivot=real(transform%pivot,rk)
        valid=.false.
        select type(geometry)
        type is(nurbs_curve)
            call geometry%transform(linear,offset,pivot)
            valid=geometry%err%ok
        type is(nurbs_surface)
            call geometry%transform(linear,offset,pivot)
            valid=geometry%err%ok
        type is(nurbs_volume)
            call geometry%transform(linear,offset,pivot)
            valid=geometry%err%ok
        end select
    end function transform_patch

    logical function transform_assembly(native_assembly,patch_kind,transform,patches) result(valid)
        class(*), intent(in) :: native_assembly
        integer, intent(in) :: patch_kind
        type(affine_transform_t), intent(in) :: transform
        type(studio_native_patch), intent(out) :: patches(:)
        class(*), allocatable :: assembly
        real(rk) :: linear(3,3),offset(3),pivot(3)
        integer :: i

        allocate(assembly,source=native_assembly)
        linear=real(transform%linear,rk)
        offset=real(transform%offset,rk)
        pivot=real(transform%pivot,rk)
        valid=.false.
        select type(assembly)
        type is(nurbs_multipatch_curve)
            if(assembly%get_npatch()/=size(patches).or.patch_kind/=1)return
            call assembly%transform(linear,offset,pivot)
            if(.not.assembly%err%ok.or..not.assembly%is_valid())return
            do i=1,size(patches)
                allocate(patches(i)%geometry,source=assembly%get_patch(i))
            end do
        type is(nurbs_multipatch_surface)
            if(assembly%get_npatch()/=size(patches).or.patch_kind/=2)return
            call assembly%transform(linear,offset,pivot)
            if(.not.assembly%err%ok.or..not.assembly%is_valid())return
            do i=1,size(patches)
                allocate(patches(i)%geometry,source=assembly%get_patch(i))
            end do
        type is(nurbs_multipatch_volume)
            if(assembly%get_npatch()/=size(patches).or.patch_kind/=3)return
            call assembly%transform(linear,offset,pivot)
            if(.not.assembly%err%ok.or..not.assembly%is_valid())return
            do i=1,size(patches)
                allocate(patches(i)%geometry,source=assembly%get_patch(i))
            end do
        class default
            return
        end select
        valid=.true.
    end function transform_assembly
end module studio_forcad_transforms
