module forcad_test_multipatch_curve

    use forcad_kinds, only: rk
    use forcad_multipatch, only: multipatch_connection
    use forcad_multipatch_curve, only: nurbs_multipatch_curve
    use forcad_nurbs_curve, only: nurbs_curve
    use forunittest, only: unit_tests

    implicit none

    private
    public :: run_multipatch_curve_tests

contains

    subroutine forcad_multipatch_curve_0001(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1
        type(nurbs_curve) :: c2
        type(nurbs_multipatch_curve) :: mp
        integer :: id1
        integer :: id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)

        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)
        call c1%err%print()
        call c2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve npatch",&
            res      = mp%get_npatch(),&
            expected = 2,&
            msg      = "curve multipatch patch count",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0001


    subroutine forcad_multipatch_curve_0002(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1
        type(nurbs_curve) :: c2
        type(nurbs_multipatch_curve) :: mp
        integer :: id1
        integer :: id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)

        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)
        call c1%err%print()
        call c2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve nconnection",&
            res      = mp%get_nconnection(),&
            expected = 1,&
            msg      = "curve multipatch connection count",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0002


    subroutine forcad_multipatch_curve_0003(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_tol = 1.0e-10_rk
        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1
        type(nurbs_curve) :: c2
        type(nurbs_multipatch_curve) :: mp
        integer :: id1
        integer :: id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)

        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)
        call c1%err%print()
        call c2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve valid C1",&
            res      = mp%is_valid(multipatch_curve_tol),&
            expected = .true.,&
            msg      = "curve C1 connection should be valid",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0003


    subroutine forcad_multipatch_curve_0004(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1
        type(nurbs_curve) :: c2
        type(nurbs_multipatch_curve) :: mp
        type(multipatch_connection) :: conn
        integer :: id1
        integer :: id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)

        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)

        conn = mp%get_connection(1)
        call c1%err%print()
        call c2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve connection patch a",&
            res      = conn%get_patch_a(),&
            expected = 1,&
            msg      = "curve connection first patch id mismatch",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0004


    subroutine forcad_multipatch_curve_0005(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1
        type(nurbs_curve) :: c2
        type(nurbs_multipatch_curve) :: mp
        type(multipatch_connection) :: conn
        integer :: id1
        integer :: id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)

        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)

        conn = mp%get_connection(1)
        call c1%err%print()
        call c2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve connection patch b",&
            res      = conn%get_patch_b(),&
            expected = 2,&
            msg      = "curve connection second patch id mismatch",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0005


    subroutine forcad_multipatch_curve_0006(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1
        type(nurbs_curve) :: c2
        type(nurbs_multipatch_curve) :: mp
        type(multipatch_connection) :: conn
        integer :: id1
        integer :: id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)

        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)

        conn = mp%get_connection(1)
        call c1%err%print()
        call c2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve continuity metadata",&
            res      = conn%get_continuity(),&
            expected = 1,&
            msg      = "curve connection continuity was not preserved",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0006


    subroutine forcad_multipatch_curve_0007(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1
        type(nurbs_curve) :: c2
        type(nurbs_multipatch_curve) :: mp
        type(multipatch_connection) :: conn
        integer :: id1
        integer :: id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)

        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)

        conn = mp%get_connection(1)
        call c1%err%print()
        call c2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve reverse metadata",&
            res      = conn%is_reversed(),&
            expected = .false.,&
            msg      = "curve connection reverse flag mismatch",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0007


    subroutine forcad_multipatch_curve_0008(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1
        type(nurbs_curve) :: c2
        type(nurbs_multipatch_curve) :: mp
        type(multipatch_connection) :: conn
        integer, allocatable :: offsets(:)
        integer :: id1
        integer :: id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)

        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)

        conn = mp%get_connection(1)

        offsets = mp%cmp_dof_offsets()
        call c1%err%print()
        call c2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve offsets",&
            res      = offsets,&
            expected = [0, 2, 4],&
            msg      = "curve patch dof offsets mismatch",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0008


    subroutine forcad_multipatch_curve_0009(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1
        type(nurbs_curve) :: c2
        type(nurbs_multipatch_curve) :: mp
        type(multipatch_connection) :: conn
        integer, allocatable :: map(:)
        integer, allocatable :: offsets(:)
        integer :: id1
        integer :: id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)

        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)

        conn = mp%get_connection(1)

        offsets = mp%cmp_dof_offsets()

        map = mp%cmp_dof_map()
        call c1%err%print()
        call c2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve shared map",&
            res      = map,&
            expected = [1, 2, 2, 3],&
            msg      = "curve C0+ interface dofs were not shared",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0009


    subroutine forcad_multipatch_curve_0010(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1
        type(nurbs_curve) :: c2
        type(nurbs_multipatch_curve) :: mp
        type(multipatch_connection) :: conn
        real(rk), allocatable :: val(:)
        integer, allocatable :: map(:)
        integer, allocatable :: offsets(:)
        integer, allocatable :: rowptr(:)
        integer, allocatable :: col(:)
        integer :: id1
        integer :: id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)

        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)

        conn = mp%get_connection(1)

        offsets = mp%cmp_dof_offsets()

        map = mp%cmp_dof_map()

        call mp%cmp_dof_constraint(rowptr, col, val)
        call c1%err%print()
        call c2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve C1 constraint count",&
            res      = size(rowptr)-1,&
            expected = 2,&
            msg      = "Curve C1 constraint count is incorrect.",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0010


    subroutine forcad_multipatch_curve_0011(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_tol = 1.0e-10_rk
        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1
        type(nurbs_curve) :: c2
        type(nurbs_multipatch_curve) :: mp
        type(multipatch_connection) :: conn
        real(rk), allocatable :: val(:)
        real(rk), allocatable :: dof(:)
        real(rk) :: residual, row_res
        integer, allocatable :: map(:)
        integer, allocatable :: offsets(:)
        integer, allocatable :: rowptr(:)
        integer, allocatable :: col(:)
        integer :: id1, id2, row, j

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)

        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)

        conn = mp%get_connection(1)

        offsets = mp%cmp_dof_offsets()

        map = mp%cmp_dof_map()

        call mp%cmp_dof_constraint(rowptr, col, val)
        dof = [multipatch_curve_x1(:,1), multipatch_curve_x2(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do
        call c1%err%print()
        call c2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve C1 constraint residual",&
            res      = residual,&
            expected = 0.0_rk,&
            tol      = multipatch_curve_tol,&
            msg      = "curve C1 constraints should vanish for a straight two-patch line",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0011


    subroutine forcad_multipatch_curve_0012(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1
        type(nurbs_curve) :: c2
        type(nurbs_multipatch_curve) :: mp
        type(multipatch_connection) :: conn
        real(rk), allocatable :: val(:)
        real(rk), allocatable :: dof(:)
        real(rk) :: residual, row_res
        integer, allocatable :: map(:)
        integer, allocatable :: elem(:,:)
        integer, allocatable :: offsets(:)
        integer, allocatable :: rowptr(:)
        integer, allocatable :: col(:)
        integer :: id1, id2, row, j

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)

        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)

        conn = mp%get_connection(1)

        offsets = mp%cmp_dof_offsets()

        map = mp%cmp_dof_map()

        call mp%cmp_dof_constraint(rowptr, col, val)
        dof = [multipatch_curve_x1(:,1), multipatch_curve_x2(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do

        elem = mp%get_elem(shared=.true.)
        call c1%err%print()
        call c2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve shared elem",&
            res      = elem,&
            expected = reshape([1, 2, 2, 3], [2,2]),&
            msg      = "curve shared element connectivity mismatch",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0012


    subroutine forcad_multipatch_curve_0013(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1
        type(nurbs_curve) :: c2
        type(nurbs_multipatch_curve) :: mp
        type(multipatch_connection) :: conn
        real(rk), allocatable :: val(:)
        real(rk), allocatable :: dof(:)
        real(rk) :: residual, row_res
        integer, allocatable :: map(:)
        integer, allocatable :: elem(:,:)
        integer, allocatable :: offsets(:)
        integer, allocatable :: rowptr(:)
        integer, allocatable :: col(:)
        integer :: id1, id2, row, j

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)

        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)

        conn = mp%get_connection(1)

        offsets = mp%cmp_dof_offsets()

        map = mp%cmp_dof_map()

        call mp%cmp_dof_constraint(rowptr, col, val)
        dof = [multipatch_curve_x1(:,1), multipatch_curve_x2(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do

        elem = mp%get_elem(shared=.true.)
        elem = mp%cmp_elem(shared=.false.)
        call c1%err%print()
        call c2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve unshared elem",&
            res      = elem,&
            expected = reshape([1, 3, 2, 4], [2,2]),&
            msg      = "curve unshared element connectivity mismatch",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0013


    subroutine forcad_multipatch_curve_0014(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1
        type(nurbs_curve) :: c2
        type(nurbs_multipatch_curve) :: mp
        type(multipatch_connection) :: conn
        real(rk), allocatable :: val(:)
        real(rk), allocatable :: dof(:)
        real(rk) :: residual, row_res
        integer, allocatable :: map(:), elem(:,:), patch_id(:), local_id(:), offsets(:), rowptr(:), col(:)
        integer :: id1, id2, row, j

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)

        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)

        conn = mp%get_connection(1)

        offsets = mp%cmp_dof_offsets()

        map = mp%cmp_dof_map()

        call mp%cmp_dof_constraint(rowptr, col, val)
        dof = [multipatch_curve_x1(:,1), multipatch_curve_x2(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do

        elem = mp%get_elem(shared=.true.)
        elem = mp%cmp_elem(shared=.false.)

        patch_id = mp%cmp_elem_patch()
        local_id = mp%cmp_elem_local()
        call c1%err%print()
        call c2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve elem patch ids",&
            res      = patch_id,&
            expected = [1, 2],&
            msg      = "curve element patch ids mismatch",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0014


    subroutine forcad_multipatch_curve_0015(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1
        type(nurbs_curve) :: c2
        type(nurbs_multipatch_curve) :: mp
        type(multipatch_connection) :: conn
        real(rk), allocatable :: val(:)
        real(rk), allocatable :: dof(:)
        real(rk) :: residual, row_res
        integer, allocatable :: map(:), elem(:,:), patch_id(:), local_id(:), offsets(:), rowptr(:), col(:)
        integer :: id1, id2, row, j

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)

        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)

        conn = mp%get_connection(1)

        offsets = mp%cmp_dof_offsets()

        map = mp%cmp_dof_map()

        call mp%cmp_dof_constraint(rowptr, col, val)
        dof = [multipatch_curve_x1(:,1), multipatch_curve_x2(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do

        elem = mp%get_elem(shared=.true.)
        elem = mp%cmp_elem(shared=.false.)

        patch_id = mp%cmp_elem_patch()
        local_id = mp%cmp_elem_local()
        call c1%err%print()
        call c2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve elem local ids",&
            res      = local_id,&
            expected = [1, 1],&
            msg      = "curve local element ids mismatch",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0015


    subroutine forcad_multipatch_curve_0016(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_tol = 1.0e-10_rk
        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1
        type(nurbs_curve) :: c2
        type(nurbs_curve) :: patch
        type(nurbs_multipatch_curve) :: mp
        type(multipatch_connection) :: conn
        real(rk), allocatable :: Xc(:,:), val(:), dof(:)
        real(rk) :: residual, row_res
        integer, allocatable :: map(:), elem(:,:), patch_id(:), local_id(:), offsets(:), rowptr(:), col(:)
        integer :: id1, id2, row, j

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)

        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)

        conn = mp%get_connection(1)

        offsets = mp%cmp_dof_offsets()

        map = mp%cmp_dof_map()

        call mp%cmp_dof_constraint(rowptr, col, val)
        dof = [multipatch_curve_x1(:,1), multipatch_curve_x2(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do

        elem = mp%get_elem(shared=.true.)
        elem = mp%cmp_elem(shared=.false.)

        patch_id = mp%cmp_elem_patch()
        local_id = mp%cmp_elem_local()

        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        call c1%err%print()
        call c2%err%print()
        call patch%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve get patch",&
            res      = Xc(1,:),&
            expected = [1.0_rk, 0.0_rk, 0.0_rk],&
            tol      = multipatch_curve_tol,&
            msg      = "curve get_patch returned wrong patch",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0016


    subroutine forcad_multipatch_curve_0017(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1
        type(nurbs_curve) :: c2
        type(nurbs_curve) :: patch
        type(nurbs_multipatch_curve) :: mp
        type(multipatch_connection) :: conn
        real(rk), allocatable :: Xc(:,:), val(:), dof(:)
        real(rk) :: residual, row_res
        integer, allocatable :: map(:), elem(:,:), patch_id(:), local_id(:), offsets(:), rowptr(:), col(:)
        integer :: id1, id2, row, j

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)

        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)

        conn = mp%get_connection(1)

        offsets = mp%cmp_dof_offsets()

        map = mp%cmp_dof_map()

        call mp%cmp_dof_constraint(rowptr, col, val)
        dof = [multipatch_curve_x1(:,1), multipatch_curve_x2(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do

        elem = mp%get_elem(shared=.true.)
        elem = mp%cmp_elem(shared=.false.)

        patch_id = mp%cmp_elem_patch()
        local_id = mp%cmp_elem_local()

        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        call c1%err%print()
        call c2%err%print()
        call patch%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve rational flag",&
            res      = mp%is_rational(),&
            expected = .false.,&
            msg      = "curve multipatch should report non-rational patches",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0017


    subroutine forcad_multipatch_curve_0018(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_tol = 1.0e-10_rk
        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1
        type(nurbs_curve) :: c2
        type(nurbs_curve) :: cr1
        type(nurbs_curve) :: cr2
        type(nurbs_curve) :: patch
        type(nurbs_multipatch_curve) :: mp
        type(nurbs_multipatch_curve) :: mp_rational
        type(multipatch_connection) :: conn
        real(rk) :: wr1(2)
        real(rk) :: wr2(2)
        real(rk), allocatable :: Xc(:,:), val(:), dof(:)
        real(rk) :: residual, row_res
        integer, allocatable :: map(:), elem(:,:), patch_id(:), local_id(:), offsets(:), rowptr(:), col(:)
        integer :: id1, id2, row, j

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)

        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)

        conn = mp%get_connection(1)

        offsets = mp%cmp_dof_offsets()

        map = mp%cmp_dof_map()

        call mp%cmp_dof_constraint(rowptr, col, val)
        dof = [multipatch_curve_x1(:,1), multipatch_curve_x2(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do

        elem = mp%get_elem(shared=.true.)
        elem = mp%cmp_elem(shared=.false.)

        patch_id = mp%cmp_elem_patch()
        local_id = mp%cmp_elem_local()

        patch = mp%get_patch(2)
        Xc = patch%get_Xc()

        wr1 = [1.0_rk, 2.0_rk]
        wr2 = [6.0_rk, 3.0_rk]
        call cr1%set(multipatch_curve_knot, multipatch_curve_x1, wr1)
        call cr2%set(multipatch_curve_knot, multipatch_curve_x2, wr2)
        call mp_rational%add_patch(cr1, id1)
        call mp_rational%add_patch(cr2, id2)
        call mp_rational%connect(id1, "right", id2, "left", 0)
        call c1%err%print()
        call c2%err%print()
        call cr1%err%print()
        call cr2%err%print()
        call patch%err%print()
        call mp%err%print()
        call mp_rational%err%print()

        call ut%test(ti)%check(&
            name     = "curve rational scaled interface",&
            res      = mp_rational%is_valid(multipatch_curve_tol),&
            expected = .true.,&
            msg      = "Curve rational scaled interface is incorrect.",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0018


    subroutine forcad_multipatch_curve_0019(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_tol = 1.0e-10_rk
        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1, c2, cr1, cr2, c2_scaled, patch
        type(nurbs_multipatch_curve) :: mp, mp_rational, mp_scaled_domain
        type(multipatch_connection) :: conn
        real(rk) :: knot_scaled(4), wr1(2), wr2(2)
        real(rk), allocatable :: Xc(:,:), val(:), dof(:)
        real(rk) :: residual, row_res
        integer, allocatable :: map(:), elem(:,:), patch_id(:), local_id(:), offsets(:), rowptr(:), col(:)
        integer :: id1, id2, row, j

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)

        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)

        conn = mp%get_connection(1)

        offsets = mp%cmp_dof_offsets()

        map = mp%cmp_dof_map()

        call mp%cmp_dof_constraint(rowptr, col, val)
        dof = [multipatch_curve_x1(:,1), multipatch_curve_x2(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do

        elem = mp%get_elem(shared=.true.)
        elem = mp%cmp_elem(shared=.false.)

        patch_id = mp%cmp_elem_patch()
        local_id = mp%cmp_elem_local()

        patch = mp%get_patch(2)
        Xc = patch%get_Xc()

        wr1 = [1.0_rk, 2.0_rk]
        wr2 = [6.0_rk, 3.0_rk]
        call cr1%set(multipatch_curve_knot, multipatch_curve_x1, wr1)
        call cr2%set(multipatch_curve_knot, multipatch_curve_x2, wr2)
        call mp_rational%add_patch(cr1, id1)
        call mp_rational%add_patch(cr2, id2)
        call mp_rational%connect(id1, "right", id2, "left", 0)

        knot_scaled = [0.0_rk, 0.0_rk, 2.0_rk, 2.0_rk]
        call c2_scaled%set(knot_scaled, multipatch_curve_x2)
        call mp_scaled_domain%add_patch(c1, id1)
        call mp_scaled_domain%add_patch(c2_scaled, id2)
        call mp_scaled_domain%connect(id1, "right", id2, "left", 1)
        call c1%err%print()
        call c2%err%print()
        call cr1%err%print()
        call cr2%err%print()
        call c2_scaled%err%print()
        call patch%err%print()
        call mp%err%print()
        call mp_rational%err%print()
        call mp_scaled_domain%err%print()

        call ut%test(ti)%check(&
            name     = "curve affine-scaled C1 domain",&
            res      = mp_scaled_domain%is_valid(multipatch_curve_tol),&
            expected = .true.,&
            msg      = "Curve affine-scaled C1 domain is incorrect.",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0019


    subroutine forcad_multipatch_curve_0020(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_tol = 1.0e-10_rk
        type(nurbs_curve) :: c3_b
        type(nurbs_curve) :: g3_a
        type(nurbs_curve) :: g3_b
        type(nurbs_multipatch_curve) :: mp_c3
        real(rk) :: xc3_b(4,3), xg3_a(4,3), xg3_b(4,3), knot_cubic(8)
        integer :: id1
        integer :: id2

        knot_cubic = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg3_a = 0.0_rk
        xg3_b = 0.0_rk
        xg3_a(:,1) = [0.0_rk, 1.0_rk/3.0_rk, 2.0_rk/3.0_rk, 1.0_rk]
        xg3_b(:,1) = [1.0_rk, 5.0_rk/3.0_rk, 8.0_rk/3.0_rk, 4.5_rk]
        xc3_b = 0.0_rk
        xc3_b(:,1) = [1.0_rk, 4.0_rk/3.0_rk, 5.0_rk/3.0_rk, 2.0_rk]
        call g3_a%set(knot_cubic, xg3_a)
        call g3_b%set(knot_cubic, xg3_b)
        call c3_b%set(knot_cubic, xc3_b)

        call mp_c3%add_patch(g3_a, id1)
        call mp_c3%add_patch(c3_b, id2)
        call mp_c3%connect(id1, "right", id2, "left", 3)
        call c3_b%err%print()
        call g3_a%err%print()
        call g3_b%err%print()
        call mp_c3%err%print()

        call ut%test(ti)%check(&
            name     = "curve valid C3",&
            res      = mp_c3%is_valid(multipatch_curve_tol),&
            expected = .true.,&
            msg      = "affinely parameterized collinear cubic patches should be C3",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0020


    subroutine forcad_multipatch_curve_0021(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_curve) :: c3_b
        type(nurbs_curve) :: g3_a
        type(nurbs_curve) :: g3_b
        type(nurbs_multipatch_curve) :: mp_c3
        real(rk) :: xc3_b(4,3), xg3_a(4,3), xg3_b(4,3), knot_cubic(8)
        real(rk), allocatable :: val(:)
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1
        integer :: id2

        knot_cubic = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg3_a = 0.0_rk
        xg3_b = 0.0_rk
        xg3_a(:,1) = [0.0_rk, 1.0_rk/3.0_rk, 2.0_rk/3.0_rk, 1.0_rk]
        xg3_b(:,1) = [1.0_rk, 5.0_rk/3.0_rk, 8.0_rk/3.0_rk, 4.5_rk]
        xc3_b = 0.0_rk
        xc3_b(:,1) = [1.0_rk, 4.0_rk/3.0_rk, 5.0_rk/3.0_rk, 2.0_rk]
        call g3_a%set(knot_cubic, xg3_a)
        call g3_b%set(knot_cubic, xg3_b)
        call c3_b%set(knot_cubic, xc3_b)

        call mp_c3%add_patch(g3_a, id1)
        call mp_c3%add_patch(c3_b, id2)
        call mp_c3%connect(id1, "right", id2, "left", 3)
        call mp_c3%cmp_dof_constraint(rowptr, col, val, geometric=.false.)
        call c3_b%err%print()
        call g3_a%err%print()
        call g3_b%err%print()
        call mp_c3%err%print()

        call ut%test(ti)%check(&
            name     = "curve C3 constraint count",&
            res      = size(rowptr)-1,&
            expected = 4,&
            msg      = "curve C3 should create derivative rows zero through three",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0021


    subroutine forcad_multipatch_curve_0022(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_tol = 1.0e-10_rk
        type(nurbs_curve) :: c3_b
        type(nurbs_curve) :: g3_a
        type(nurbs_curve) :: g3_b
        type(nurbs_multipatch_curve) :: mp_c3
        real(rk) :: xc3_b(4,3), xg3_a(4,3), xg3_b(4,3), knot_cubic(8)
        real(rk), allocatable :: val(:)
        real(rk), allocatable :: dof(:)
        real(rk) :: residual, row_res
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1, id2, row, j

        knot_cubic = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg3_a = 0.0_rk
        xg3_b = 0.0_rk
        xg3_a(:,1) = [0.0_rk, 1.0_rk/3.0_rk, 2.0_rk/3.0_rk, 1.0_rk]
        xg3_b(:,1) = [1.0_rk, 5.0_rk/3.0_rk, 8.0_rk/3.0_rk, 4.5_rk]
        xc3_b = 0.0_rk
        xc3_b(:,1) = [1.0_rk, 4.0_rk/3.0_rk, 5.0_rk/3.0_rk, 2.0_rk]
        call g3_a%set(knot_cubic, xg3_a)
        call g3_b%set(knot_cubic, xg3_b)
        call c3_b%set(knot_cubic, xc3_b)

        call mp_c3%add_patch(g3_a, id1)
        call mp_c3%add_patch(c3_b, id2)
        call mp_c3%connect(id1, "right", id2, "left", 3)
        call mp_c3%cmp_dof_constraint(rowptr, col, val, geometric=.false.)
        dof = [xg3_a(:,1), xc3_b(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do
        call c3_b%err%print()
        call g3_a%err%print()
        call g3_b%err%print()
        call mp_c3%err%print()

        call ut%test(ti)%check(&
            name     = "curve C3 constraint residual",&
            res      = residual,&
            expected = 0.0_rk,&
            tol      = multipatch_curve_tol,&
            msg      = "curve arbitrary-order C3 constraints should vanish",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0022


    subroutine forcad_multipatch_curve_0023(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_tol = 1.0e-10_rk
        type(nurbs_curve) :: c3_b
        type(nurbs_curve) :: g3_a
        type(nurbs_curve) :: g3_b
        type(nurbs_multipatch_curve) :: mp_c3
        type(nurbs_multipatch_curve) :: mp_g3
        real(rk) :: xc3_b(4,3), xg3_a(4,3), xg3_b(4,3), knot_cubic(8)
        real(rk), allocatable :: val(:), dof(:)
        real(rk) :: residual, row_res
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1, id2, row, j

        knot_cubic = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg3_a = 0.0_rk
        xg3_b = 0.0_rk
        xg3_a(:,1) = [0.0_rk, 1.0_rk/3.0_rk, 2.0_rk/3.0_rk, 1.0_rk]
        xg3_b(:,1) = [1.0_rk, 5.0_rk/3.0_rk, 8.0_rk/3.0_rk, 4.5_rk]
        xc3_b = 0.0_rk
        xc3_b(:,1) = [1.0_rk, 4.0_rk/3.0_rk, 5.0_rk/3.0_rk, 2.0_rk]
        call g3_a%set(knot_cubic, xg3_a)
        call g3_b%set(knot_cubic, xg3_b)
        call c3_b%set(knot_cubic, xc3_b)

        call mp_c3%add_patch(g3_a, id1)
        call mp_c3%add_patch(c3_b, id2)
        call mp_c3%connect(id1, "right", id2, "left", 3)
        call mp_c3%cmp_dof_constraint(rowptr, col, val, geometric=.false.)
        dof = [xg3_a(:,1), xc3_b(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do

        call mp_g3%add_patch(g3_a, id1)
        call mp_g3%add_patch(g3_b, id2)
        call mp_g3%connect(&
            patch_a            = id1,&
            side_a             = "right",&
            patch_b            = id2,&
            side_b             = "left",&
            continuity         = 3,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk, 0.1875_rk])
        call c3_b%err%print()
        call g3_a%err%print()
        call g3_b%err%print()
        call mp_c3%err%print()
        call mp_g3%err%print()

        call ut%test(ti)%check(&
            name     = "curve valid G3",&
            res      = mp_g3%is_valid(multipatch_curve_tol),&
            expected = .true.,&
            msg      = "nonlinearly reparameterized collinear cubic patches should be G3",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0023


    subroutine forcad_multipatch_curve_0024(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_curve) :: c3_b
        type(nurbs_curve) :: g3_a
        type(nurbs_curve) :: g3_b
        type(nurbs_multipatch_curve) :: mp_c3
        type(nurbs_multipatch_curve) :: mp_g3
        type(multipatch_connection) :: conn
        real(rk) :: xc3_b(4,3), xg3_a(4,3), xg3_b(4,3), knot_cubic(8)
        real(rk), allocatable :: val(:), dof(:)
        real(rk) :: residual, row_res
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1, id2, row, j

        knot_cubic = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg3_a = 0.0_rk
        xg3_b = 0.0_rk
        xg3_a(:,1) = [0.0_rk, 1.0_rk/3.0_rk, 2.0_rk/3.0_rk, 1.0_rk]
        xg3_b(:,1) = [1.0_rk, 5.0_rk/3.0_rk, 8.0_rk/3.0_rk, 4.5_rk]
        xc3_b = 0.0_rk
        xc3_b(:,1) = [1.0_rk, 4.0_rk/3.0_rk, 5.0_rk/3.0_rk, 2.0_rk]
        call g3_a%set(knot_cubic, xg3_a)
        call g3_b%set(knot_cubic, xg3_b)
        call c3_b%set(knot_cubic, xc3_b)

        call mp_c3%add_patch(g3_a, id1)
        call mp_c3%add_patch(c3_b, id2)
        call mp_c3%connect(id1, "right", id2, "left", 3)
        call mp_c3%cmp_dof_constraint(rowptr, col, val, geometric=.false.)
        dof = [xg3_a(:,1), xc3_b(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do

        call mp_g3%add_patch(g3_a, id1)
        call mp_g3%add_patch(g3_b, id2)
        call mp_g3%connect(&
            patch_a            = id1,&
            side_a             = "right",&
            patch_b            = id2,&
            side_b             = "left",&
            continuity         = 3,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk, 0.1875_rk])
        conn = mp_g3%get_connection(1)
        call c3_b%err%print()
        call g3_a%err%print()
        call g3_b%err%print()
        call mp_c3%err%print()
        call mp_g3%err%print()

        call ut%test(ti)%check(&
            name     = "curve G3 continuity flag",&
            res      = conn%is_geometric(),&
            expected = .true.,&
            msg      = "curve geometric continuity was not preserved",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0024


    subroutine forcad_multipatch_curve_0025(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_tol = 1.0e-10_rk
        type(nurbs_curve) :: c3_b
        type(nurbs_curve) :: g3_a
        type(nurbs_curve) :: g3_b
        type(nurbs_multipatch_curve) :: mp_c3
        type(nurbs_multipatch_curve) :: mp_g3
        type(multipatch_connection) :: conn
        real(rk) :: xc3_b(4,3), xg3_a(4,3), xg3_b(4,3), knot_cubic(8)
        real(rk), allocatable :: val(:), dof(:), reparameterization(:)
        real(rk) :: residual, row_res
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1, id2, row, j

        knot_cubic = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg3_a = 0.0_rk
        xg3_b = 0.0_rk
        xg3_a(:,1) = [0.0_rk, 1.0_rk/3.0_rk, 2.0_rk/3.0_rk, 1.0_rk]
        xg3_b(:,1) = [1.0_rk, 5.0_rk/3.0_rk, 8.0_rk/3.0_rk, 4.5_rk]
        xc3_b = 0.0_rk
        xc3_b(:,1) = [1.0_rk, 4.0_rk/3.0_rk, 5.0_rk/3.0_rk, 2.0_rk]
        call g3_a%set(knot_cubic, xg3_a)
        call g3_b%set(knot_cubic, xg3_b)
        call c3_b%set(knot_cubic, xc3_b)

        call mp_c3%add_patch(g3_a, id1)
        call mp_c3%add_patch(c3_b, id2)
        call mp_c3%connect(id1, "right", id2, "left", 3)
        call mp_c3%cmp_dof_constraint(rowptr, col, val, geometric=.false.)
        dof = [xg3_a(:,1), xc3_b(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do

        call mp_g3%add_patch(g3_a, id1)
        call mp_g3%add_patch(g3_b, id2)
        call mp_g3%connect(&
            patch_a            = id1,&
            side_a             = "right",&
            patch_b            = id2,&
            side_b             = "left",&
            continuity         = 3,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk, 0.1875_rk])
        conn = mp_g3%get_connection(1)
        reparameterization = conn%get_reparameterization()
        call c3_b%err%print()
        call g3_a%err%print()
        call g3_b%err%print()
        call mp_c3%err%print()
        call mp_g3%err%print()

        call ut%test(ti)%check(&
            name     = "curve G3 transition jet",&
            res      = reparameterization,&
            expected = [0.5_rk, -0.25_rk, 0.1875_rk],&
            tol      = multipatch_curve_tol,&
            msg      = "curve G3 transition derivatives were not preserved",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0025


    subroutine forcad_multipatch_curve_0026(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_curve) :: c3_b
        type(nurbs_curve) :: g3_a
        type(nurbs_curve) :: g3_b
        type(nurbs_multipatch_curve) :: mp_c3
        type(nurbs_multipatch_curve) :: mp_g3
        type(multipatch_connection) :: conn
        real(rk) :: xc3_b(4,3), xg3_a(4,3), xg3_b(4,3), knot_cubic(8)
        real(rk), allocatable :: val(:), dof(:), reparameterization(:)
        real(rk) :: residual, row_res
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1, id2, row, j

        knot_cubic = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg3_a = 0.0_rk
        xg3_b = 0.0_rk
        xg3_a(:,1) = [0.0_rk, 1.0_rk/3.0_rk, 2.0_rk/3.0_rk, 1.0_rk]
        xg3_b(:,1) = [1.0_rk, 5.0_rk/3.0_rk, 8.0_rk/3.0_rk, 4.5_rk]
        xc3_b = 0.0_rk
        xc3_b(:,1) = [1.0_rk, 4.0_rk/3.0_rk, 5.0_rk/3.0_rk, 2.0_rk]
        call g3_a%set(knot_cubic, xg3_a)
        call g3_b%set(knot_cubic, xg3_b)
        call c3_b%set(knot_cubic, xc3_b)

        call mp_c3%add_patch(g3_a, id1)
        call mp_c3%add_patch(c3_b, id2)
        call mp_c3%connect(id1, "right", id2, "left", 3)
        call mp_c3%cmp_dof_constraint(rowptr, col, val, geometric=.false.)
        dof = [xg3_a(:,1), xc3_b(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do

        call mp_g3%add_patch(g3_a, id1)
        call mp_g3%add_patch(g3_b, id2)
        call mp_g3%connect(&
            patch_a            = id1,&
            side_a             = "right",&
            patch_b            = id2,&
            side_b             = "left",&
            continuity         = 3,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk, 0.1875_rk])
        conn = mp_g3%get_connection(1)
        reparameterization = conn%get_reparameterization()
        call mp_g3%cmp_dof_constraint(rowptr, col, val, geometric=.true.)
        call c3_b%err%print()
        call g3_a%err%print()
        call g3_b%err%print()
        call mp_c3%err%print()
        call mp_g3%err%print()

        call ut%test(ti)%check(&
            name     = "curve G3 constraint count",&
            res      = size(rowptr)-1,&
            expected = 4,&
            msg      = "curve G3 should create derivative rows zero through three",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0026


    subroutine forcad_multipatch_curve_0027(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_tol = 1.0e-10_rk
        type(nurbs_curve) :: c3_b
        type(nurbs_curve) :: g3_a
        type(nurbs_curve) :: g3_b
        type(nurbs_multipatch_curve) :: mp_c3
        type(nurbs_multipatch_curve) :: mp_g3
        type(multipatch_connection) :: conn
        real(rk) :: xc3_b(4,3), xg3_a(4,3), xg3_b(4,3), knot_cubic(8)
        real(rk), allocatable :: val(:), dof(:), reparameterization(:)
        real(rk) :: residual, row_res
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1, id2, row, j

        knot_cubic = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg3_a = 0.0_rk
        xg3_b = 0.0_rk
        xg3_a(:,1) = [0.0_rk, 1.0_rk/3.0_rk, 2.0_rk/3.0_rk, 1.0_rk]
        xg3_b(:,1) = [1.0_rk, 5.0_rk/3.0_rk, 8.0_rk/3.0_rk, 4.5_rk]
        xc3_b = 0.0_rk
        xc3_b(:,1) = [1.0_rk, 4.0_rk/3.0_rk, 5.0_rk/3.0_rk, 2.0_rk]
        call g3_a%set(knot_cubic, xg3_a)
        call g3_b%set(knot_cubic, xg3_b)
        call c3_b%set(knot_cubic, xc3_b)

        call mp_c3%add_patch(g3_a, id1)
        call mp_c3%add_patch(c3_b, id2)
        call mp_c3%connect(id1, "right", id2, "left", 3)
        call mp_c3%cmp_dof_constraint(rowptr, col, val, geometric=.false.)
        dof = [xg3_a(:,1), xc3_b(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do

        call mp_g3%add_patch(g3_a, id1)
        call mp_g3%add_patch(g3_b, id2)
        call mp_g3%connect(&
            patch_a            = id1,&
            side_a             = "right",&
            patch_b            = id2,&
            side_b             = "left",&
            continuity         = 3,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk, 0.1875_rk])
        conn = mp_g3%get_connection(1)
        reparameterization = conn%get_reparameterization()
        call mp_g3%cmp_dof_constraint(rowptr, col, val, geometric=.true.)
        dof = [xg3_a(:,1), xg3_b(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do
        call c3_b%err%print()
        call g3_a%err%print()
        call g3_b%err%print()
        call mp_c3%err%print()
        call mp_g3%err%print()

        call ut%test(ti)%check(&
            name     = "curve G3 constraint residual",&
            res      = residual,&
            expected = 0.0_rk,&
            tol      = multipatch_curve_tol,&
            msg      = "curve G3 chain-rule constraints should vanish",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0027


    subroutine forcad_multipatch_curve_0028(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_curve) :: c3_b
        type(nurbs_curve) :: g3_a
        type(nurbs_curve) :: g3_b
        type(nurbs_multipatch_curve) :: mp_c3
        type(nurbs_multipatch_curve) :: mp_g3
        type(multipatch_connection) :: conn
        real(rk) :: xc3_b(4,3), xg3_a(4,3), xg3_b(4,3), knot_cubic(8)
        real(rk), allocatable :: val(:), dof(:), reparameterization(:)
        real(rk) :: residual, row_res
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1, id2, row, j

        knot_cubic = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg3_a = 0.0_rk
        xg3_b = 0.0_rk
        xg3_a(:,1) = [0.0_rk, 1.0_rk/3.0_rk, 2.0_rk/3.0_rk, 1.0_rk]
        xg3_b(:,1) = [1.0_rk, 5.0_rk/3.0_rk, 8.0_rk/3.0_rk, 4.5_rk]
        xc3_b = 0.0_rk
        xc3_b(:,1) = [1.0_rk, 4.0_rk/3.0_rk, 5.0_rk/3.0_rk, 2.0_rk]
        call g3_a%set(knot_cubic, xg3_a)
        call g3_b%set(knot_cubic, xg3_b)
        call c3_b%set(knot_cubic, xc3_b)

        call mp_c3%add_patch(g3_a, id1)
        call mp_c3%add_patch(c3_b, id2)
        call mp_c3%connect(id1, "right", id2, "left", 3)
        call mp_c3%cmp_dof_constraint(rowptr, col, val, geometric=.false.)
        dof = [xg3_a(:,1), xc3_b(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do

        call mp_g3%add_patch(g3_a, id1)
        call mp_g3%add_patch(g3_b, id2)
        call mp_g3%connect(&
            patch_a            = id1,&
            side_a             = "right",&
            patch_b            = id2,&
            side_b             = "left",&
            continuity         = 3,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk, 0.1875_rk])
        conn = mp_g3%get_connection(1)
        reparameterization = conn%get_reparameterization()
        call mp_g3%cmp_dof_constraint(rowptr, col, val, geometric=.true.)
        dof = [xg3_a(:,1), xg3_b(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do
        call mp_g3%cmp_dof_constraint(rowptr, col, val, geometric=.false.)
        call c3_b%err%print()
        call g3_a%err%print()
        call g3_b%err%print()
        call mp_c3%err%print()
        call mp_g3%err%print()

        call ut%test(ti)%check(&
            name     = "curve separate C constraints",&
            res      = size(rowptr)-1,&
            expected = 0,&
            msg      = "filtering a G3 connection for C constraints must return no rows",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0028


    subroutine forcad_multipatch_curve_0029(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_tol = 1.0e-10_rk
        type(nurbs_curve) :: c3_b
        type(nurbs_curve) :: g3_a
        type(nurbs_curve) :: g3_b
        type(nurbs_multipatch_curve) :: mp_c3
        type(nurbs_multipatch_curve) :: mp_g3
        type(nurbs_multipatch_curve) :: mp_not_c3
        type(multipatch_connection) :: conn
        real(rk) :: xc3_b(4,3), xg3_a(4,3), xg3_b(4,3), knot_cubic(8)
        real(rk), allocatable :: val(:), dof(:), reparameterization(:)
        real(rk) :: residual, row_res
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1, id2, row, j

        knot_cubic = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg3_a = 0.0_rk
        xg3_b = 0.0_rk
        xg3_a(:,1) = [0.0_rk, 1.0_rk/3.0_rk, 2.0_rk/3.0_rk, 1.0_rk]
        xg3_b(:,1) = [1.0_rk, 5.0_rk/3.0_rk, 8.0_rk/3.0_rk, 4.5_rk]
        xc3_b = 0.0_rk
        xc3_b(:,1) = [1.0_rk, 4.0_rk/3.0_rk, 5.0_rk/3.0_rk, 2.0_rk]
        call g3_a%set(knot_cubic, xg3_a)
        call g3_b%set(knot_cubic, xg3_b)
        call c3_b%set(knot_cubic, xc3_b)

        call mp_c3%add_patch(g3_a, id1)
        call mp_c3%add_patch(c3_b, id2)
        call mp_c3%connect(id1, "right", id2, "left", 3)
        call mp_c3%cmp_dof_constraint(rowptr, col, val, geometric=.false.)
        dof = [xg3_a(:,1), xc3_b(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do

        call mp_g3%add_patch(g3_a, id1)
        call mp_g3%add_patch(g3_b, id2)
        call mp_g3%connect(&
            patch_a            = id1,&
            side_a             = "right",&
            patch_b            = id2,&
            side_b             = "left",&
            continuity         = 3,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk, 0.1875_rk])
        conn = mp_g3%get_connection(1)
        reparameterization = conn%get_reparameterization()
        call mp_g3%cmp_dof_constraint(rowptr, col, val, geometric=.true.)
        dof = [xg3_a(:,1), xg3_b(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do
        call mp_g3%cmp_dof_constraint(rowptr, col, val, geometric=.false.)

        call mp_not_c3%add_patch(g3_a, id1)
        call mp_not_c3%add_patch(g3_b, id2)
        call mp_not_c3%connect(id1, "right", id2, "left", 3)
        call c3_b%err%print()
        call g3_a%err%print()
        call g3_b%err%print()
        call mp_c3%err%print()
        call mp_g3%err%print()
        call mp_not_c3%err%print()

        call ut%test(ti)%check(&
            name     = "curve G3 is not C3",&
            res      = mp_not_c3%is_valid(multipatch_curve_tol),&
            expected = .false.,&
            msg      = "The G3 curve interface must not satisfy C3.",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0029


    subroutine forcad_multipatch_curve_0030(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_tol = 1.0e-10_rk
        type(nurbs_curve) :: c3_b, g3_a, g3_b, rg3_a, rg3_b
        type(nurbs_multipatch_curve) :: mp_c3
        type(nurbs_multipatch_curve) :: mp_g3
        type(nurbs_multipatch_curve) :: mp_not_c3
        type(nurbs_multipatch_curve) :: mp_rg3
        type(multipatch_connection) :: conn
        real(rk) :: xc3_b(4,3), xg3_a(4,3), xg3_b(4,3), knot_cubic(8)
        real(rk), allocatable :: val(:), dof(:), reparameterization(:)
        real(rk) :: residual, row_res
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1, id2, row, j

        knot_cubic = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg3_a = 0.0_rk
        xg3_b = 0.0_rk
        xg3_a(:,1) = [0.0_rk, 1.0_rk/3.0_rk, 2.0_rk/3.0_rk, 1.0_rk]
        xg3_b(:,1) = [1.0_rk, 5.0_rk/3.0_rk, 8.0_rk/3.0_rk, 4.5_rk]
        xc3_b = 0.0_rk
        xc3_b(:,1) = [1.0_rk, 4.0_rk/3.0_rk, 5.0_rk/3.0_rk, 2.0_rk]
        call g3_a%set(knot_cubic, xg3_a)
        call g3_b%set(knot_cubic, xg3_b)
        call c3_b%set(knot_cubic, xc3_b)

        call mp_c3%add_patch(g3_a, id1)
        call mp_c3%add_patch(c3_b, id2)
        call mp_c3%connect(id1, "right", id2, "left", 3)
        call mp_c3%cmp_dof_constraint(rowptr, col, val, geometric=.false.)
        dof = [xg3_a(:,1), xc3_b(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do

        call mp_g3%add_patch(g3_a, id1)
        call mp_g3%add_patch(g3_b, id2)
        call mp_g3%connect(&
            patch_a            = id1,&
            side_a             = "right",&
            patch_b            = id2,&
            side_b             = "left",&
            continuity         = 3,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk, 0.1875_rk])
        conn = mp_g3%get_connection(1)
        reparameterization = conn%get_reparameterization()
        call mp_g3%cmp_dof_constraint(rowptr, col, val, geometric=.true.)
        dof = [xg3_a(:,1), xg3_b(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do
        call mp_g3%cmp_dof_constraint(rowptr, col, val, geometric=.false.)

        call mp_not_c3%add_patch(g3_a, id1)
        call mp_not_c3%add_patch(g3_b, id2)
        call mp_not_c3%connect(id1, "right", id2, "left", 3)

        call rg3_a%set(knot_cubic, xg3_a, [2.0_rk, 2.0_rk, 2.0_rk, 2.0_rk])
        call rg3_b%set(knot_cubic, xg3_b, [3.0_rk, 3.0_rk, 3.0_rk, 3.0_rk])
        call mp_rg3%add_patch(rg3_a, id1)
        call mp_rg3%add_patch(rg3_b, id2)
        call mp_rg3%connect(&
            patch_a            = id1,&
            side_a             = "right",&
            patch_b            = id2,&
            side_b             = "left",&
            continuity         = 3,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk, 0.1875_rk])
        call c3_b%err%print()
        call g3_a%err%print()
        call g3_b%err%print()
        call rg3_a%err%print()
        call rg3_b%err%print()
        call mp_c3%err%print()
        call mp_g3%err%print()
        call mp_not_c3%err%print()
        call mp_rg3%err%print()

        call ut%test(ti)%check(&
            name     = "curve rational G3",&
            res      = mp_rg3%is_valid(multipatch_curve_tol),&
            expected = .true.,&
            msg      = "Curve rational G3 is incorrect.",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0030


    subroutine forcad_multipatch_curve_0031(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_curve) :: c3_b, g3_a, g3_b, rg3_a, rg3_b
        type(nurbs_multipatch_curve) :: mp_c3, mp_g3, mp_not_c3, mp_rg3, mp_bad_g
        type(multipatch_connection) :: conn
        real(rk) :: xc3_b(4,3), xg3_a(4,3), xg3_b(4,3), knot_cubic(8)
        real(rk), allocatable :: val(:), dof(:), reparameterization(:)
        real(rk) :: residual, row_res
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1, id2, row, j

        knot_cubic = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg3_a = 0.0_rk
        xg3_b = 0.0_rk
        xg3_a(:,1) = [0.0_rk, 1.0_rk/3.0_rk, 2.0_rk/3.0_rk, 1.0_rk]
        xg3_b(:,1) = [1.0_rk, 5.0_rk/3.0_rk, 8.0_rk/3.0_rk, 4.5_rk]
        xc3_b = 0.0_rk
        xc3_b(:,1) = [1.0_rk, 4.0_rk/3.0_rk, 5.0_rk/3.0_rk, 2.0_rk]
        call g3_a%set(knot_cubic, xg3_a)
        call g3_b%set(knot_cubic, xg3_b)
        call c3_b%set(knot_cubic, xc3_b)

        call mp_c3%add_patch(g3_a, id1)
        call mp_c3%add_patch(c3_b, id2)
        call mp_c3%connect(id1, "right", id2, "left", 3)
        call mp_c3%cmp_dof_constraint(rowptr, col, val, geometric=.false.)
        dof = [xg3_a(:,1), xc3_b(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do

        call mp_g3%add_patch(g3_a, id1)
        call mp_g3%add_patch(g3_b, id2)
        call mp_g3%connect(&
            patch_a            = id1,&
            side_a             = "right",&
            patch_b            = id2,&
            side_b             = "left",&
            continuity         = 3,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk, 0.1875_rk])
        conn = mp_g3%get_connection(1)
        reparameterization = conn%get_reparameterization()
        call mp_g3%cmp_dof_constraint(rowptr, col, val, geometric=.true.)
        dof = [xg3_a(:,1), xg3_b(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do
        call mp_g3%cmp_dof_constraint(rowptr, col, val, geometric=.false.)

        call mp_not_c3%add_patch(g3_a, id1)
        call mp_not_c3%add_patch(g3_b, id2)
        call mp_not_c3%connect(id1, "right", id2, "left", 3)

        call rg3_a%set(knot_cubic, xg3_a, [2.0_rk, 2.0_rk, 2.0_rk, 2.0_rk])
        call rg3_b%set(knot_cubic, xg3_b, [3.0_rk, 3.0_rk, 3.0_rk, 3.0_rk])
        call mp_rg3%add_patch(rg3_a, id1)
        call mp_rg3%add_patch(rg3_b, id2)
        call mp_rg3%connect(&
            patch_a            = id1,&
            side_a             = "right",&
            patch_b            = id2,&
            side_b             = "left",&
            continuity         = 3,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk, 0.1875_rk])

        call mp_bad_g%add_patch(g3_a, id1)
        call mp_bad_g%add_patch(g3_b, id2)
        call mp_bad_g%connect(&
            patch_a         = id1,&
            side_a          = "right",&
            patch_b         = id2,&
            side_b          = "left",&
            continuity      = 2,&
            geometric          = .true.)
        call c3_b%err%print()
        call g3_a%err%print()
        call g3_b%err%print()
        call rg3_a%err%print()
        call rg3_b%err%print()
        call mp_c3%err%print()
        call mp_g3%err%print()
        call mp_not_c3%err%print()
        call mp_rg3%err%print()
        call mp_bad_g%err%print()

        call ut%test(ti)%check(&
            name     = "curve missing G2 jet err",&
            res      = mp_bad_g%err%ok,&
            expected = .false.,&
            msg      = "a G2 connection without a transition jet should set an error",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0031


    subroutine forcad_multipatch_curve_0032(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_curve) :: c3_b, g3_a, g3_b, rg3_a, rg3_b
        type(nurbs_multipatch_curve) :: mp_c3, mp_g3, mp_not_c3, mp_rg3, mp_bad_g
        type(multipatch_connection) :: conn
        real(rk) :: xc3_b(4,3), xg3_a(4,3), xg3_b(4,3), knot_cubic(8)
        real(rk), allocatable :: val(:), dof(:), reparameterization(:)
        real(rk) :: residual, row_res
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1, id2, row, j

        knot_cubic = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg3_a = 0.0_rk
        xg3_b = 0.0_rk
        xg3_a(:,1) = [0.0_rk, 1.0_rk/3.0_rk, 2.0_rk/3.0_rk, 1.0_rk]
        xg3_b(:,1) = [1.0_rk, 5.0_rk/3.0_rk, 8.0_rk/3.0_rk, 4.5_rk]
        xc3_b = 0.0_rk
        xc3_b(:,1) = [1.0_rk, 4.0_rk/3.0_rk, 5.0_rk/3.0_rk, 2.0_rk]
        call g3_a%set(knot_cubic, xg3_a)
        call g3_b%set(knot_cubic, xg3_b)
        call c3_b%set(knot_cubic, xc3_b)

        call mp_c3%add_patch(g3_a, id1)
        call mp_c3%add_patch(c3_b, id2)
        call mp_c3%connect(id1, "right", id2, "left", 3)
        call mp_c3%cmp_dof_constraint(rowptr, col, val, geometric=.false.)
        dof = [xg3_a(:,1), xc3_b(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do

        call mp_g3%add_patch(g3_a, id1)
        call mp_g3%add_patch(g3_b, id2)
        call mp_g3%connect(&
            patch_a            = id1,&
            side_a             = "right",&
            patch_b            = id2,&
            side_b             = "left",&
            continuity         = 3,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk, 0.1875_rk])
        conn = mp_g3%get_connection(1)
        reparameterization = conn%get_reparameterization()
        call mp_g3%cmp_dof_constraint(rowptr, col, val, geometric=.true.)
        dof = [xg3_a(:,1), xg3_b(:,1)]
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do
        call mp_g3%cmp_dof_constraint(rowptr, col, val, geometric=.false.)

        call mp_not_c3%add_patch(g3_a, id1)
        call mp_not_c3%add_patch(g3_b, id2)
        call mp_not_c3%connect(id1, "right", id2, "left", 3)

        call rg3_a%set(knot_cubic, xg3_a, [2.0_rk, 2.0_rk, 2.0_rk, 2.0_rk])
        call rg3_b%set(knot_cubic, xg3_b, [3.0_rk, 3.0_rk, 3.0_rk, 3.0_rk])
        call mp_rg3%add_patch(rg3_a, id1)
        call mp_rg3%add_patch(rg3_b, id2)
        call mp_rg3%connect(&
            patch_a            = id1,&
            side_a             = "right",&
            patch_b            = id2,&
            side_b             = "left",&
            continuity         = 3,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk, 0.1875_rk])

        call mp_bad_g%add_patch(g3_a, id1)
        call mp_bad_g%add_patch(g3_b, id2)
        call mp_bad_g%connect(&
            patch_a         = id1,&
            side_a          = "right",&
            patch_b         = id2,&
            side_b          = "left",&
            continuity      = 2,&
            geometric          = .true.)
        call c3_b%err%print()
        call g3_a%err%print()
        call g3_b%err%print()
        call rg3_a%err%print()
        call rg3_b%err%print()
        call mp_c3%err%print()
        call mp_g3%err%print()
        call mp_not_c3%err%print()
        call mp_rg3%err%print()
        call mp_bad_g%err%print()

        call ut%test(ti)%check(&
            name     = "curve missing G2 jet no connection",&
            res      = mp_bad_g%get_nconnection(),&
            expected = 0,&
            msg      = "an invalid G2 transition must not append a connection",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0032


    subroutine forcad_multipatch_curve_0033(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_tol = 1.0e-10_rk
        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1, c2, patch
        type(nurbs_multipatch_curve) :: mp, mp_copy
        real(rk), allocatable :: Xc(:,:), Xg(:,:)
        integer :: id1, id2, file_unit

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)
        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)
        call mp%create(res=3)
        call mp%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%translate_Xg([0.0_rk, 0.0_rk, 1.0_rk])
        call mp%rotate_Xc(0.0_rk, 0.0_rk, 0.0_rk)
        call mp%rotate_Xg(0.0_rk, 0.0_rk, 0.0_rk)
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_curve_Xc_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_curve_Xg_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_curve_Xth_1.vtk", status="replace")
        close(file_unit, status="delete")
        call mp%export_Xc("vtk/forcad_test_multipatch_curve")
        call mp%export_Xg("vtk/forcad_test_multipatch_curve")
        call mp%export_Xth_in_Xg("vtk/forcad_test_multipatch_curve", res=3)
        mp_copy = mp

        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        Xg = patch%get_Xg()
        call c1%err%print()
        call c2%err%print()
        call patch%err%print()
        call mp%err%print()
        call mp_copy%err%print()

        call ut%test(ti)%check(&
            name     = "curve translate Xc",&
            res      = Xc(1,:),&
            expected = [2.0_rk, 0.0_rk, 0.0_rk],&
            tol      = multipatch_curve_tol,&
            msg      = "curve multipatch translate_Xc failed",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0033


    subroutine forcad_multipatch_curve_0034(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_tol = 1.0e-10_rk
        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1, c2, patch
        type(nurbs_multipatch_curve) :: mp, mp_copy
        real(rk), allocatable :: Xc(:,:), Xg(:,:)
        integer :: id1, id2, file_unit

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)
        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)
        call mp%create(res=3)
        call mp%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%translate_Xg([0.0_rk, 0.0_rk, 1.0_rk])
        call mp%rotate_Xc(0.0_rk, 0.0_rk, 0.0_rk)
        call mp%rotate_Xg(0.0_rk, 0.0_rk, 0.0_rk)
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_curve_Xc_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_curve_Xg_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_curve_Xth_1.vtk", status="replace")
        close(file_unit, status="delete")
        call mp%export_Xc("vtk/forcad_test_multipatch_curve")
        call mp%export_Xg("vtk/forcad_test_multipatch_curve")
        call mp%export_Xth_in_Xg("vtk/forcad_test_multipatch_curve", res=3)
        mp_copy = mp

        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        Xg = patch%get_Xg()
        call c1%err%print()
        call c2%err%print()
        call patch%err%print()
        call mp%err%print()
        call mp_copy%err%print()

        call ut%test(ti)%check(&
            name     = "curve translate Xg",&
            res      = Xg(1,3),&
            expected = 1.0_rk,&
            tol      = multipatch_curve_tol,&
            msg      = "curve multipatch translate_Xg failed",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0034


    subroutine forcad_multipatch_curve_0035(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1, c2, patch
        type(nurbs_multipatch_curve) :: mp, mp_copy
        real(rk), allocatable :: Xc(:,:), Xg(:,:)
        integer :: id1, id2, file_unit
        logical :: file_exists

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)
        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)
        call mp%create(res=3)
        call mp%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%translate_Xg([0.0_rk, 0.0_rk, 1.0_rk])
        call mp%rotate_Xc(0.0_rk, 0.0_rk, 0.0_rk)
        call mp%rotate_Xg(0.0_rk, 0.0_rk, 0.0_rk)
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_curve_Xc_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_curve_Xg_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_curve_Xth_1.vtk", status="replace")
        close(file_unit, status="delete")
        call mp%export_Xc("vtk/forcad_test_multipatch_curve")
        call mp%export_Xg("vtk/forcad_test_multipatch_curve")
        call mp%export_Xth_in_Xg("vtk/forcad_test_multipatch_curve", res=3)
        mp_copy = mp

        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        Xg = patch%get_Xg()
        inquire(file="vtk/forcad_test_multipatch_curve_Xc_1.vtk", exist=file_exists)
        call c1%err%print()
        call c2%err%print()
        call patch%err%print()
        call mp%err%print()
        call mp_copy%err%print()

        call ut%test(ti)%check(&
            name     = "curve export Xc",&
            res      = file_exists,&
            expected = .true.,&
            msg      = "curve multipatch export_Xc did not create the first VTK file",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0035


    subroutine forcad_multipatch_curve_0036(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1, c2, patch
        type(nurbs_multipatch_curve) :: mp, mp_copy
        real(rk), allocatable :: Xc(:,:), Xg(:,:)
        integer :: id1, id2, file_unit
        logical :: file_exists

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)
        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)
        call mp%create(res=3)
        call mp%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%translate_Xg([0.0_rk, 0.0_rk, 1.0_rk])
        call mp%rotate_Xc(0.0_rk, 0.0_rk, 0.0_rk)
        call mp%rotate_Xg(0.0_rk, 0.0_rk, 0.0_rk)
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_curve_Xc_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_curve_Xg_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_curve_Xth_1.vtk", status="replace")
        close(file_unit, status="delete")
        call mp%export_Xc("vtk/forcad_test_multipatch_curve")
        call mp%export_Xg("vtk/forcad_test_multipatch_curve")
        call mp%export_Xth_in_Xg("vtk/forcad_test_multipatch_curve", res=3)
        mp_copy = mp

        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        Xg = patch%get_Xg()
        inquire(file="vtk/forcad_test_multipatch_curve_Xc_1.vtk", exist=file_exists)
        inquire(file="vtk/forcad_test_multipatch_curve_Xg_1.vtk", exist=file_exists)
        call c1%err%print()
        call c2%err%print()
        call patch%err%print()
        call mp%err%print()
        call mp_copy%err%print()

        call ut%test(ti)%check(&
            name     = "curve export Xg",&
            res      = file_exists,&
            expected = .true.,&
            msg      = "curve multipatch export_Xg did not create the first VTK file",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0036


    subroutine forcad_multipatch_curve_0037(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1, c2, patch
        type(nurbs_multipatch_curve) :: mp, mp_copy
        real(rk), allocatable :: Xc(:,:), Xg(:,:)
        integer :: id1, id2, file_unit
        logical :: file_exists

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)
        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        call mp%connect(id1, "right", id2, "left", 1)
        call mp%create(res=3)
        call mp%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%translate_Xg([0.0_rk, 0.0_rk, 1.0_rk])
        call mp%rotate_Xc(0.0_rk, 0.0_rk, 0.0_rk)
        call mp%rotate_Xg(0.0_rk, 0.0_rk, 0.0_rk)
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_curve_Xc_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_curve_Xg_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_curve_Xth_1.vtk", status="replace")
        close(file_unit, status="delete")
        call mp%export_Xc("vtk/forcad_test_multipatch_curve")
        call mp%export_Xg("vtk/forcad_test_multipatch_curve")
        call mp%export_Xth_in_Xg("vtk/forcad_test_multipatch_curve", res=3)
        mp_copy = mp

        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        Xg = patch%get_Xg()
        inquire(file="vtk/forcad_test_multipatch_curve_Xc_1.vtk", exist=file_exists)
        inquire(file="vtk/forcad_test_multipatch_curve_Xg_1.vtk", exist=file_exists)
        inquire(file="vtk/forcad_test_multipatch_curve_Xth_1.vtk", exist=file_exists)
        call c1%err%print()
        call c2%err%print()
        call patch%err%print()
        call mp%err%print()
        call mp_copy%err%print()

        call ut%test(ti)%check(&
            name     = "curve export Xth in Xg",&
            res      = file_exists,&
            expected = .true.,&
            msg      = "Curve export Xth in Xg is incorrect.",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0037


    subroutine forcad_multipatch_curve_0038(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1, c2
        type(nurbs_multipatch_curve) :: mp_discontinuous
        integer, allocatable :: map(:)
        integer :: id1, id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)
        call mp_discontinuous%add_patch(c1, id1)
        call mp_discontinuous%add_patch(c2, id2)
        call mp_discontinuous%connect(id1, "right", id2, "left", -1)
        map = mp_discontinuous%cmp_dof_map()
        call c1%err%print()
        call c2%err%print()
        call mp_discontinuous%err%print()

        call ut%test(ti)%check(&
            name     = "curve discontinuous map",&
            res      = maxval(map),&
            expected = 4,&
            msg      = "discontinuous curve interface should not share dofs",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0038


    subroutine forcad_multipatch_curve_0039(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_tol = 1.0e-10_rk
        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1, c2
        type(nurbs_multipatch_curve) :: mp_discontinuous
        type(nurbs_multipatch_curve) :: mp_bad_geom
        integer, allocatable :: map(:)
        integer :: id1, id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)
        call mp_discontinuous%add_patch(c1, id1)
        call mp_discontinuous%add_patch(c2, id2)
        call mp_discontinuous%connect(id1, "right", id2, "left", -1)
        map = mp_discontinuous%cmp_dof_map()

        call mp_bad_geom%add_patch(c1, id1)
        call mp_bad_geom%add_patch(c2, id2)
        call mp_bad_geom%connect(id1, "right", id2, "right", 0)
        call c1%err%print()
        call c2%err%print()
        call mp_discontinuous%err%print()
        call mp_bad_geom%err%print()

        call ut%test(ti)%check(&
            name     = "curve invalid geometry",&
            res      = mp_bad_geom%is_valid(multipatch_curve_tol),&
            expected = .false.,&
            msg      = "curve multipatch should reject mismatched interface geometry",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0039


    subroutine forcad_multipatch_curve_0040(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1, c2
        type(nurbs_multipatch_curve) :: mp_discontinuous
        type(nurbs_multipatch_curve) :: mp_bad_geom
        type(nurbs_multipatch_curve) :: mp_bad_continuity
        integer, allocatable :: map(:)
        integer :: id1, id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)
        call mp_discontinuous%add_patch(c1, id1)
        call mp_discontinuous%add_patch(c2, id2)
        call mp_discontinuous%connect(id1, "right", id2, "left", -1)
        map = mp_discontinuous%cmp_dof_map()

        call mp_bad_geom%add_patch(c1, id1)
        call mp_bad_geom%add_patch(c2, id2)
        call mp_bad_geom%connect(id1, "right", id2, "right", 0)

        call mp_bad_continuity%add_patch(c1, id1)
        call mp_bad_continuity%add_patch(c2, id2)
        call mp_bad_continuity%connect(id1, "right", id2, "left", 2)
        call c1%err%print()
        call c2%err%print()
        call mp_discontinuous%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_continuity%err%print()

        call ut%test(ti)%check(&
            name     = "curve invalid continuity err",&
            res      = mp_bad_continuity%err%ok,&
            expected = .false.,&
            msg      = "invalid curve continuity should set multipatch error state",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0040


    subroutine forcad_multipatch_curve_0041(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1, c2
        type(nurbs_multipatch_curve) :: mp_discontinuous
        type(nurbs_multipatch_curve) :: mp_bad_geom
        type(nurbs_multipatch_curve) :: mp_bad_continuity
        integer, allocatable :: map(:)
        integer :: id1, id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)
        call mp_discontinuous%add_patch(c1, id1)
        call mp_discontinuous%add_patch(c2, id2)
        call mp_discontinuous%connect(id1, "right", id2, "left", -1)
        map = mp_discontinuous%cmp_dof_map()

        call mp_bad_geom%add_patch(c1, id1)
        call mp_bad_geom%add_patch(c2, id2)
        call mp_bad_geom%connect(id1, "right", id2, "right", 0)

        call mp_bad_continuity%add_patch(c1, id1)
        call mp_bad_continuity%add_patch(c2, id2)
        call mp_bad_continuity%connect(id1, "right", id2, "left", 2)
        call c1%err%print()
        call c2%err%print()
        call mp_discontinuous%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_continuity%err%print()

        call ut%test(ti)%check(&
            name     = "curve invalid continuity no connection",&
            res      = mp_bad_continuity%get_nconnection(),&
            expected = 0,&
            msg      = "invalid curve continuity should not append a connection",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0041


    subroutine forcad_multipatch_curve_0042(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_tol = 1.0e-10_rk
        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1, c2
        type(nurbs_multipatch_curve) :: mp_discontinuous
        type(nurbs_multipatch_curve) :: mp_bad_geom
        type(nurbs_multipatch_curve) :: mp_bad_continuity
        integer, allocatable :: map(:)
        integer :: id1, id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)
        call mp_discontinuous%add_patch(c1, id1)
        call mp_discontinuous%add_patch(c2, id2)
        call mp_discontinuous%connect(id1, "right", id2, "left", -1)
        map = mp_discontinuous%cmp_dof_map()

        call mp_bad_geom%add_patch(c1, id1)
        call mp_bad_geom%add_patch(c2, id2)
        call mp_bad_geom%connect(id1, "right", id2, "right", 0)

        call mp_bad_continuity%add_patch(c1, id1)
        call mp_bad_continuity%add_patch(c2, id2)
        call mp_bad_continuity%connect(id1, "right", id2, "left", 2)
        call c1%err%print()
        call c2%err%print()
        call mp_discontinuous%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_continuity%err%print()

        call ut%test(ti)%check(&
            name     = "curve invalid continuity",&
            res      = mp_bad_continuity%is_valid(multipatch_curve_tol),&
            expected = .false.,&
            msg      = "curve multipatch should reject continuity above patch degree",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0042


    subroutine forcad_multipatch_curve_0043(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1, c2
        type(nurbs_multipatch_curve) :: mp_discontinuous, mp_bad_geom, mp_bad_continuity, mp_bad
        integer, allocatable :: map(:)
        integer :: id1, id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)
        call mp_discontinuous%add_patch(c1, id1)
        call mp_discontinuous%add_patch(c2, id2)
        call mp_discontinuous%connect(id1, "right", id2, "left", -1)
        map = mp_discontinuous%cmp_dof_map()

        call mp_bad_geom%add_patch(c1, id1)
        call mp_bad_geom%add_patch(c2, id2)
        call mp_bad_geom%connect(id1, "right", id2, "right", 0)

        call mp_bad_continuity%add_patch(c1, id1)
        call mp_bad_continuity%add_patch(c2, id2)
        call mp_bad_continuity%connect(id1, "right", id2, "left", 2)

        call mp_bad%add_patch(c1, id1)
        call mp_bad%connect(id1, "invalid", id1, "left", 0)
        call c1%err%print()
        call c2%err%print()
        call mp_discontinuous%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_continuity%err%print()
        call mp_bad%err%print()

        call ut%test(ti)%check(&
            name     = "curve invalid side err",&
            res      = mp_bad%err%ok,&
            expected = .false.,&
            msg      = "invalid curve side should set multipatch error state",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0043


    subroutine forcad_multipatch_curve_0044(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1, c2
        type(nurbs_multipatch_curve) :: mp_discontinuous, mp_bad_geom, mp_bad_continuity, mp_bad
        integer, allocatable :: map(:)
        integer :: id1, id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)
        call mp_discontinuous%add_patch(c1, id1)
        call mp_discontinuous%add_patch(c2, id2)
        call mp_discontinuous%connect(id1, "right", id2, "left", -1)
        map = mp_discontinuous%cmp_dof_map()

        call mp_bad_geom%add_patch(c1, id1)
        call mp_bad_geom%add_patch(c2, id2)
        call mp_bad_geom%connect(id1, "right", id2, "right", 0)

        call mp_bad_continuity%add_patch(c1, id1)
        call mp_bad_continuity%add_patch(c2, id2)
        call mp_bad_continuity%connect(id1, "right", id2, "left", 2)

        call mp_bad%add_patch(c1, id1)
        call mp_bad%connect(id1, "invalid", id1, "left", 0)
        call c1%err%print()
        call c2%err%print()
        call mp_discontinuous%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_continuity%err%print()
        call mp_bad%err%print()

        call ut%test(ti)%check(&
            name     = "curve invalid side no connection",&
            res      = mp_bad%get_nconnection(),&
            expected = 0,&
            msg      = "invalid curve side should not append a connection",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0044


    subroutine forcad_multipatch_curve_0045(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1, c2
        type(nurbs_multipatch_curve) :: mp_discontinuous, mp_bad_geom, mp_bad_continuity, mp_bad
        integer, allocatable :: map(:)
        integer :: id1, id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)
        call mp_discontinuous%add_patch(c1, id1)
        call mp_discontinuous%add_patch(c2, id2)
        call mp_discontinuous%connect(id1, "right", id2, "left", -1)
        map = mp_discontinuous%cmp_dof_map()

        call mp_bad_geom%add_patch(c1, id1)
        call mp_bad_geom%add_patch(c2, id2)
        call mp_bad_geom%connect(id1, "right", id2, "right", 0)

        call mp_bad_continuity%add_patch(c1, id1)
        call mp_bad_continuity%add_patch(c2, id2)
        call mp_bad_continuity%connect(id1, "right", id2, "left", 2)

        call mp_bad%add_patch(c1, id1)
        call mp_bad%connect(id1, "invalid", id1, "left", 0)
        call mp_bad%add_patch(c2, id2)
        call c1%err%print()
        call c2%err%print()
        call mp_discontinuous%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_continuity%err%print()
        call mp_bad%err%print()

        call ut%test(ti)%check(&
            name     = "curve add blocked on err",&
            res      = id2,&
            expected = 0,&
            msg      = "add_patch must leave id zero while multipatch error state is set",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0045


    subroutine forcad_multipatch_curve_0046(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1, c2
        type(nurbs_multipatch_curve) :: mp_discontinuous, mp_bad_geom, mp_bad_continuity, mp_bad
        integer, allocatable :: map(:)
        integer :: id1, id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)
        call mp_discontinuous%add_patch(c1, id1)
        call mp_discontinuous%add_patch(c2, id2)
        call mp_discontinuous%connect(id1, "right", id2, "left", -1)
        map = mp_discontinuous%cmp_dof_map()

        call mp_bad_geom%add_patch(c1, id1)
        call mp_bad_geom%add_patch(c2, id2)
        call mp_bad_geom%connect(id1, "right", id2, "right", 0)

        call mp_bad_continuity%add_patch(c1, id1)
        call mp_bad_continuity%add_patch(c2, id2)
        call mp_bad_continuity%connect(id1, "right", id2, "left", 2)

        call mp_bad%add_patch(c1, id1)
        call mp_bad%connect(id1, "invalid", id1, "left", 0)
        call mp_bad%add_patch(c2, id2)
        call mp_bad%err%print()
        call mp_bad%err%reset()
        call mp_bad%add_patch(c2, id2)
        call c1%err%print()
        call c2%err%print()
        call mp_discontinuous%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_continuity%err%print()
        call mp_bad%err%print()

        call ut%test(ti)%check(&
            name     = "curve add after reset",&
            res      = id2,&
            expected = 2,&
            msg      = "debug reset should allow multipatch operations again",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0046


    subroutine forcad_multipatch_curve_0047(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: invalid_patch
        type(nurbs_curve) :: c1
        type(nurbs_multipatch_curve) :: invalid_patch_mp
        integer :: id1

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call invalid_patch%set(&
            knot   = multipatch_curve_knot,&
            Xc     = multipatch_curve_x1,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        call invalid_patch_mp%add_patch(invalid_patch, id1)
        call invalid_patch%err%print()
        call c1%err%print()
        call invalid_patch_mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve rejects invalid patch",&
            res      = [invalid_patch_mp%get_npatch(), id1],&
            expected = [0, 0],&
            msg      = "An invalid curve patch must not be appended or receive an id.",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0047


    subroutine forcad_multipatch_curve_0048(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: invalid_patch
        type(nurbs_curve) :: c1
        type(nurbs_multipatch_curve) :: capacity_mp, invalid_patch_mp
        integer :: i, ids(18), id1

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call invalid_patch%set(&
            knot   = multipatch_curve_knot,&
            Xc     = multipatch_curve_x1,&
            Wc     = [1.0_rk, 0.0_rk],&
            degree = 1)
        call invalid_patch_mp%add_patch(invalid_patch, id1)

        do i = 1, size(ids)
            call capacity_mp%add_patch(c1, ids(i))
        end do
        do i = 1, size(ids)-1
            call capacity_mp%connect(ids(i), "right", ids(i+1), "right", -1)
        end do
        call invalid_patch%err%print()
        call c1%err%print()
        call capacity_mp%err%print()
        call invalid_patch_mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve dynamic storage growth",&
            res      = [capacity_mp%get_npatch(), capacity_mp%get_nconnection()],&
            expected = [18, 17],&
            msg      = "Curve dynamic storage growth is incorrect.",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0048


    subroutine forcad_multipatch_curve_0049(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1, c2
        type(nurbs_multipatch_curve) :: mp, mp_copy
        integer :: id1, id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)
        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        mp_copy = mp
        call mp%finalize()
        call c1%err%print()
        call c2%err%print()
        call mp%err%print()
        call mp_copy%err%print()

        call ut%test(ti)%check(&
            name     = "curve finalize",&
            res      = mp%get_npatch(),&
            expected = 0,&
            msg      = "curve multipatch finalize failed",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0049


    subroutine forcad_multipatch_curve_0050(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_curve_knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: multipatch_curve_x1(2,3) = reshape([&
            0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        real(rk), parameter :: multipatch_curve_x2(2,3) = reshape([&
            1.0_rk,2.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk], [2,3])
        type(nurbs_curve) :: c1, c2
        type(nurbs_multipatch_curve) :: mp, mp_copy
        integer :: id1, id2

        call c1%set(multipatch_curve_knot, multipatch_curve_x1)
        call c2%set(multipatch_curve_knot, multipatch_curve_x2)
        call mp%add_patch(c1, id1)
        call mp%add_patch(c2, id2)
        mp_copy = mp
        call mp%finalize()
        call c1%err%print()
        call c2%err%print()
        call mp%err%print()
        call mp_copy%err%print()

        call ut%test(ti)%check(&
            name     = "curve multipatch assignment",&
            res      = mp_copy%get_npatch(),&
            expected = 2,&
            msg      = "Curve multipatch assignment is incorrect.",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0050


    subroutine forcad_multipatch_curve_0051(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 1.0e-12_rk
        real(rk), parameter :: knot(7) = [0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk]
        real(rk), parameter :: Xca(4,1) = reshape([-4.0_rk,-2.0_rk,0.0_rk,2.0_rk], [4,1])
        real(rk), parameter :: Xcb(4,1) = reshape([-1.0_rk,4.0_rk,5.0_rk,6.0_rk], [4,1])
        real(rk), parameter :: Wa(4) = [2.0_rk,2.0_rk,1.0_rk,3.0_rk]
        real(rk), parameter :: Wb(4) = [2.0_rk,2.0_rk,3.0_rk,4.0_rk]
        type(nurbs_curve) :: ca, cb
        type(nurbs_multipatch_curve) :: mp
        integer :: ida, idb

        call ca%set(&
            knot   = knot,&
            Xc     = Xca,&
            Wc     = Wa,&
            degree = 2)
        call cb%set(&
            knot   = knot,&
            Xc     = Xcb,&
            Wc     = Wb,&
            degree = 2)
        call mp%add_patch(ca, ida)
        call mp%add_patch(cb, idb)
        call mp%connect(ida, "right", idb, "left", 0)
        call ca%err%print()
        call cb%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve unclamped C0 trace",&
            res      = mp%is_valid(tol),&
            expected = .true.,&
            msg      = "Basis-weighted curve endpoint traces must match.",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0051


    subroutine forcad_multipatch_curve_0052(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: knot(7) = [0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk]
        real(rk), parameter :: Xca(4,1) = reshape([-4.0_rk,-2.0_rk,0.0_rk,2.0_rk], [4,1])
        real(rk), parameter :: Xcb(4,1) = reshape([-2.0_rk,4.0_rk,5.0_rk,6.0_rk], [4,1])
        type(nurbs_curve) :: ca, cb
        type(nurbs_multipatch_curve) :: mp
        integer, allocatable :: map(:)
        integer :: ida, idb

        call ca%set(knot, Xca, degree=2)
        call cb%set(knot, Xcb, degree=2)
        call mp%add_patch(ca, ida)
        call mp%add_patch(cb, idb)
        call mp%connect(ida, "right", idb, "left", 0)
        map = mp%cmp_dof_map()
        call ca%err%print()
        call cb%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve unclamped C0 map",&
            res      = map,&
            expected = [1,2,3,4,5,6,7,8],&
            msg      = "An unclamped curve trace must remain a constraint.",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0052


    subroutine forcad_multipatch_curve_0053(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_multipatch_curve) :: mp

        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve empty multipatch invalid",&
            res      = mp%is_valid(),&
            expected = .false.,&
            msg      = "An empty curve multipatch is not geometry.",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0053


    subroutine forcad_multipatch_curve_0054(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_curve) :: curve
        type(nurbs_multipatch_curve) :: mp
        integer :: id

        call mp%add_patch(curve, id)
        call curve%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve rejects uninitialized patch",&
            res      = id,&
            expected = 0,&
            msg      = "An uninitialized curve must not be added.",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0054


    subroutine forcad_multipatch_curve_0055(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_curve) :: curve
        type(nurbs_multipatch_curve) :: mp
        integer :: id

        call curve%set(&
            Xth_dir    = [0.0_rk,1.0_rk],&
            degree     = 1,&
            continuity = [-1,-1])
        call mp%add_patch(curve, id)
        call curve%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "curve rejects space-only patch",&
            res      = id,&
            expected = 0,&
            msg      = "A curve spline space is not complete geometry.",&
            group    = "forcad_multipatch_curve")
        ti = ti + 1

    end subroutine forcad_multipatch_curve_0055


    subroutine run_multipatch_curve_tests(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        call forcad_multipatch_curve_0001(ut, ti)
        call forcad_multipatch_curve_0002(ut, ti)
        call forcad_multipatch_curve_0003(ut, ti)
        call forcad_multipatch_curve_0004(ut, ti)
        call forcad_multipatch_curve_0005(ut, ti)
        call forcad_multipatch_curve_0006(ut, ti)
        call forcad_multipatch_curve_0007(ut, ti)
        call forcad_multipatch_curve_0008(ut, ti)
        call forcad_multipatch_curve_0009(ut, ti)
        call forcad_multipatch_curve_0010(ut, ti)
        call forcad_multipatch_curve_0011(ut, ti)
        call forcad_multipatch_curve_0012(ut, ti)
        call forcad_multipatch_curve_0013(ut, ti)
        call forcad_multipatch_curve_0014(ut, ti)
        call forcad_multipatch_curve_0015(ut, ti)
        call forcad_multipatch_curve_0016(ut, ti)
        call forcad_multipatch_curve_0017(ut, ti)
        call forcad_multipatch_curve_0018(ut, ti)
        call forcad_multipatch_curve_0019(ut, ti)
        call forcad_multipatch_curve_0020(ut, ti)
        call forcad_multipatch_curve_0021(ut, ti)
        call forcad_multipatch_curve_0022(ut, ti)
        call forcad_multipatch_curve_0023(ut, ti)
        call forcad_multipatch_curve_0024(ut, ti)
        call forcad_multipatch_curve_0025(ut, ti)
        call forcad_multipatch_curve_0026(ut, ti)
        call forcad_multipatch_curve_0027(ut, ti)
        call forcad_multipatch_curve_0028(ut, ti)
        call forcad_multipatch_curve_0029(ut, ti)
        call forcad_multipatch_curve_0030(ut, ti)
        call forcad_multipatch_curve_0031(ut, ti)
        call forcad_multipatch_curve_0032(ut, ti)
        call forcad_multipatch_curve_0033(ut, ti)
        call forcad_multipatch_curve_0034(ut, ti)
        call forcad_multipatch_curve_0035(ut, ti)
        call forcad_multipatch_curve_0036(ut, ti)
        call forcad_multipatch_curve_0037(ut, ti)
        call forcad_multipatch_curve_0038(ut, ti)
        call forcad_multipatch_curve_0039(ut, ti)
        call forcad_multipatch_curve_0040(ut, ti)
        call forcad_multipatch_curve_0041(ut, ti)
        call forcad_multipatch_curve_0042(ut, ti)
        call forcad_multipatch_curve_0043(ut, ti)
        call forcad_multipatch_curve_0044(ut, ti)
        call forcad_multipatch_curve_0045(ut, ti)
        call forcad_multipatch_curve_0046(ut, ti)
        call forcad_multipatch_curve_0047(ut, ti)
        call forcad_multipatch_curve_0048(ut, ti)
        call forcad_multipatch_curve_0049(ut, ti)
        call forcad_multipatch_curve_0050(ut, ti)
        call forcad_multipatch_curve_0051(ut, ti)
        call forcad_multipatch_curve_0052(ut, ti)
        call forcad_multipatch_curve_0053(ut, ti)
        call forcad_multipatch_curve_0054(ut, ti)
        call forcad_multipatch_curve_0055(ut, ti)

    end subroutine run_multipatch_curve_tests

end module forcad_test_multipatch_curve
