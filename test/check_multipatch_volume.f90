module forcad_test_multipatch_volume

    use forcad_kinds, only: rk
    use forcad_multipatch, only: multipatch_connection
    use forcad_multipatch_volume, only: nurbs_multipatch_volume
    use forcad_nurbs_volume, only: nurbs_volume
    use forunittest, only: unit_tests

    implicit none

    private
    public :: run_multipatch_volume_tests

contains

    subroutine forcad_multipatch_volume_0001(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_volume_tol = 1.0e-10_rk
        type(nurbs_volume) :: v1, v2
        type(nurbs_multipatch_volume) :: mp
        integer :: id1
        integer :: id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])

        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        call v1%err%print()
        call v2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume valid C0",&
            res      = mp%is_valid(multipatch_volume_tol),&
            expected = .true.,&
            msg      = "volume C0 connection should be valid",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0001


    subroutine forcad_multipatch_volume_0002(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1, v2
        type(nurbs_multipatch_volume) :: mp
        integer :: id1
        integer :: id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])

        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        call v1%err%print()
        call v2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume npatch",&
            res      = mp%get_npatch(),&
            expected = 2,&
            msg      = "volume multipatch patch count mismatch",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0002


    subroutine forcad_multipatch_volume_0003(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1, v2
        type(nurbs_multipatch_volume) :: mp
        integer :: id1
        integer :: id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])

        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        call v1%err%print()
        call v2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume nconnection",&
            res      = mp%get_nconnection(),&
            expected = 1,&
            msg      = "volume multipatch connection count mismatch",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0003


    subroutine forcad_multipatch_volume_0004(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1, v2
        type(nurbs_multipatch_volume) :: mp
        integer, allocatable :: offsets(:)
        integer :: id1
        integer :: id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])

        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)

        offsets = mp%cmp_dof_offsets()
        call v1%err%print()
        call v2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume offsets",&
            res      = offsets,&
            expected = [0, 8, 16],&
            msg      = "volume patch dof offsets mismatch",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0004


    subroutine forcad_multipatch_volume_0005(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1, v2
        type(nurbs_multipatch_volume) :: mp
        integer, allocatable :: map(:)
        integer, allocatable :: offsets(:)
        integer :: id1
        integer :: id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])

        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)

        offsets = mp%cmp_dof_offsets()

        map = mp%cmp_dof_map()
        call v1%err%print()
        call v2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume shared dofs",&
            res      = maxval(map),&
            expected = 12,&
            msg      = "volume interface dofs were not shared",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0005


    subroutine forcad_multipatch_volume_0006(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1, v2
        type(nurbs_multipatch_volume) :: mp, mp_c1
        real(rk), allocatable :: val(:)
        integer, allocatable :: map(:), offsets(:), rowptr(:), col(:)
        integer :: id1
        integer :: id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])

        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)

        offsets = mp%cmp_dof_offsets()

        map = mp%cmp_dof_map()

        call mp_c1%add_patch(v1, id1)
        call mp_c1%add_patch(v2, id2)
        call mp_c1%connect(id1, "u_max", id2, "u_min", 1)
        call mp_c1%cmp_dof_constraint(rowptr, col, val)
        call v1%err%print()
        call v2%err%print()
        call mp%err%print()
        call mp_c1%err%print()

        call ut%test(ti)%check(&
            name     = "volume C1 constraint count",&
            res      = size(rowptr)-1,&
            expected = 8,&
            msg      = "Volume C1 constraint count is incorrect.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0006


    subroutine forcad_multipatch_volume_0007(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_volume_tol = 1.0e-10_rk
        type(nurbs_volume) :: v1, v2
        type(nurbs_multipatch_volume) :: mp, mp_c1
        real(rk), allocatable :: Xc(:,:), val(:), dof(:)
        real(rk) :: residual, row_res
        integer, allocatable :: map(:), offsets(:), rowptr(:), col(:)
        integer :: id1, id2, row, j

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])

        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)

        offsets = mp%cmp_dof_offsets()

        map = mp%cmp_dof_map()

        call mp_c1%add_patch(v1, id1)
        call mp_c1%add_patch(v2, id2)
        call mp_c1%connect(id1, "u_max", id2, "u_min", 1)
        call mp_c1%cmp_dof_constraint(rowptr, col, val)
        allocate(dof(16))
        Xc = v1%get_Xc()
        dof(1:8) = Xc(:,1)
        Xc = v2%get_Xc()
        dof(9:16) = Xc(:,1)
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do
        call v1%err%print()
        call v2%err%print()
        call mp%err%print()
        call mp_c1%err%print()

        call ut%test(ti)%check(&
            name     = "volume C1 constraint residual",&
            res      = residual,&
            expected = 0.0_rk,&
            tol      = multipatch_volume_tol,&
            msg      = "Volume C1 constraint residual is incorrect.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0007


    subroutine forcad_multipatch_volume_0008(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_volume_tol = 1.0e-10_rk
        type(nurbs_volume) :: vg2_a
        type(nurbs_volume) :: vg2_b
        type(nurbs_multipatch_volume) :: mp_g2
        real(rk) :: xg2_a(12,3)
        real(rk) :: xg2_b(12,3)
        real(rk) :: knot_open(4)
        real(rk) :: knot_quadratic(6)
        integer :: id1
        integer :: id2

        knot_quadratic = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot_open = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        xg2_a = 0.0_rk
        xg2_b = 0.0_rk
        xg2_a(:,1) = [0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, &
            0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk]
        xg2_a(:,2) = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, &
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_a(:,3) = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, &
            1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_b(:,1) = [1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk, &
            1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk]
        xg2_b(:,2:3) = xg2_a(:,2:3)
        call vg2_a%set(knot_quadratic, knot_open, knot_open, xg2_a, degree=[2,1,1])
        call vg2_b%set(knot_quadratic, knot_open, knot_open, xg2_b, degree=[2,1,1])
        call mp_g2%add_patch(vg2_a, id1)
        call mp_g2%add_patch(vg2_b, id2)
        call mp_g2%connect(&
            patch_a            = id1,&
            side_a             = "u_max",&
            patch_b            = id2,&
            side_b             = "u_min",&
            continuity         = 2,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk])
        call vg2_a%err%print()
        call vg2_b%err%print()
        call mp_g2%err%print()

        call ut%test(ti)%check(&
            name     = "volume valid G2",&
            res      = mp_g2%is_valid(multipatch_volume_tol),&
            expected = .true.,&
            msg      = "nonlinearly reparameterized hexahedral patches should be G2",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0008


    subroutine forcad_multipatch_volume_0009(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vg2_a
        type(nurbs_volume) :: vg2_b
        type(nurbs_multipatch_volume) :: mp_g2
        type(multipatch_connection) :: conn
        real(rk) :: xg2_a(12,3)
        real(rk) :: xg2_b(12,3)
        real(rk) :: knot_open(4)
        real(rk) :: knot_quadratic(6)
        integer :: id1
        integer :: id2

        knot_quadratic = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot_open = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        xg2_a = 0.0_rk
        xg2_b = 0.0_rk
        xg2_a(:,1) = [0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, &
            0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk]
        xg2_a(:,2) = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, &
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_a(:,3) = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, &
            1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_b(:,1) = [1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk, &
            1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk]
        xg2_b(:,2:3) = xg2_a(:,2:3)
        call vg2_a%set(knot_quadratic, knot_open, knot_open, xg2_a, degree=[2,1,1])
        call vg2_b%set(knot_quadratic, knot_open, knot_open, xg2_b, degree=[2,1,1])
        call mp_g2%add_patch(vg2_a, id1)
        call mp_g2%add_patch(vg2_b, id2)
        call mp_g2%connect(&
            patch_a            = id1,&
            side_a             = "u_max",&
            patch_b            = id2,&
            side_b             = "u_min",&
            continuity         = 2,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk])
        conn = mp_g2%get_connection(1)
        call vg2_a%err%print()
        call vg2_b%err%print()
        call mp_g2%err%print()

        call ut%test(ti)%check(&
            name     = "volume G2 continuity flag",&
            res      = conn%is_geometric(),&
            expected = .true.,&
            msg      = "volume geometric continuity was not preserved",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0009


    subroutine forcad_multipatch_volume_0010(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vg2_a
        type(nurbs_volume) :: vg2_b
        type(nurbs_multipatch_volume) :: mp_g2
        type(multipatch_connection) :: conn
        real(rk) :: xg2_a(12,3)
        real(rk) :: xg2_b(12,3)
        real(rk) :: knot_open(4)
        real(rk) :: knot_quadratic(6)
        real(rk), allocatable :: val(:)
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1
        integer :: id2

        knot_quadratic = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot_open = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        xg2_a = 0.0_rk
        xg2_b = 0.0_rk
        xg2_a(:,1) = [0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, &
            0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk]
        xg2_a(:,2) = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, &
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_a(:,3) = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, &
            1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_b(:,1) = [1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk, &
            1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk]
        xg2_b(:,2:3) = xg2_a(:,2:3)
        call vg2_a%set(knot_quadratic, knot_open, knot_open, xg2_a, degree=[2,1,1])
        call vg2_b%set(knot_quadratic, knot_open, knot_open, xg2_b, degree=[2,1,1])
        call mp_g2%add_patch(vg2_a, id1)
        call mp_g2%add_patch(vg2_b, id2)
        call mp_g2%connect(&
            patch_a            = id1,&
            side_a             = "u_max",&
            patch_b            = id2,&
            side_b             = "u_min",&
            continuity         = 2,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk])
        conn = mp_g2%get_connection(1)
        call mp_g2%cmp_dof_constraint(rowptr, col, val, geometric=.true.)
        call vg2_a%err%print()
        call vg2_b%err%print()
        call mp_g2%err%print()

        call ut%test(ti)%check(&
            name     = "volume G2 constraint count",&
            res      = size(rowptr)-1,&
            expected = 12,&
            msg      = "volume G2 should create three rows per tangential control point",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0010


    subroutine forcad_multipatch_volume_0011(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_volume_tol = 1.0e-10_rk
        type(nurbs_volume) :: vg2_a
        type(nurbs_volume) :: vg2_b
        type(nurbs_multipatch_volume) :: mp_g2
        type(multipatch_connection) :: conn
        real(rk) :: xg2_a(12,3)
        real(rk) :: xg2_b(12,3)
        real(rk) :: knot_open(4), knot_quadratic(6), residual, row_res
        real(rk), allocatable :: val(:), dof(:)
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1, id2, row, j

        knot_quadratic = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot_open = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        xg2_a = 0.0_rk
        xg2_b = 0.0_rk
        xg2_a(:,1) = [0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, &
            0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk]
        xg2_a(:,2) = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, &
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_a(:,3) = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, &
            1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_b(:,1) = [1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk, &
            1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk]
        xg2_b(:,2:3) = xg2_a(:,2:3)
        call vg2_a%set(knot_quadratic, knot_open, knot_open, xg2_a, degree=[2,1,1])
        call vg2_b%set(knot_quadratic, knot_open, knot_open, xg2_b, degree=[2,1,1])
        call mp_g2%add_patch(vg2_a, id1)
        call mp_g2%add_patch(vg2_b, id2)
        call mp_g2%connect(&
            patch_a            = id1,&
            side_a             = "u_max",&
            patch_b            = id2,&
            side_b             = "u_min",&
            continuity         = 2,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk])
        conn = mp_g2%get_connection(1)
        call mp_g2%cmp_dof_constraint(rowptr, col, val, geometric=.true.)
        allocate(dof(24))
        dof(1:12) = xg2_a(:,1)
        dof(13:24) = xg2_b(:,1)
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do
        call vg2_a%err%print()
        call vg2_b%err%print()
        call mp_g2%err%print()

        call ut%test(ti)%check(&
            name     = "volume G2 constraint residual",&
            res      = residual,&
            expected = 0.0_rk,&
            tol      = multipatch_volume_tol,&
            msg      = "volume G2 chain-rule constraints should vanish",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0011


    subroutine forcad_multipatch_volume_0012(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_volume_tol = 1.0e-10_rk
        type(nurbs_volume) :: vg2_a, vg2_b, vrg2_a, vrg2_b
        type(nurbs_multipatch_volume) :: mp_g2
        type(nurbs_multipatch_volume) :: mp_rg2
        type(multipatch_connection) :: conn
        real(rk) :: wrg2_a(12), wrg2_b(12), xg2_a(12,3), xg2_b(12,3)
        real(rk) :: knot_open(4), knot_quadratic(6), residual, row_res
        real(rk), allocatable :: val(:), dof(:)
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1, id2, row, j

        knot_quadratic = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot_open = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        xg2_a = 0.0_rk
        xg2_b = 0.0_rk
        xg2_a(:,1) = [0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, &
            0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk]
        xg2_a(:,2) = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, &
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_a(:,3) = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, &
            1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_b(:,1) = [1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk, &
            1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk]
        xg2_b(:,2:3) = xg2_a(:,2:3)
        call vg2_a%set(knot_quadratic, knot_open, knot_open, xg2_a, degree=[2,1,1])
        call vg2_b%set(knot_quadratic, knot_open, knot_open, xg2_b, degree=[2,1,1])
        call mp_g2%add_patch(vg2_a, id1)
        call mp_g2%add_patch(vg2_b, id2)
        call mp_g2%connect(&
            patch_a            = id1,&
            side_a             = "u_max",&
            patch_b            = id2,&
            side_b             = "u_min",&
            continuity         = 2,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk])
        conn = mp_g2%get_connection(1)
        call mp_g2%cmp_dof_constraint(rowptr, col, val, geometric=.true.)
        allocate(dof(24))
        dof(1:12) = xg2_a(:,1)
        dof(13:24) = xg2_b(:,1)
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do
        deallocate(dof)
        wrg2_a = 2.0_rk
        wrg2_b = 3.0_rk
        call vrg2_a%set(knot_quadratic, knot_open, knot_open, xg2_a, wrg2_a, degree=[2,1,1])
        call vrg2_b%set(knot_quadratic, knot_open, knot_open, xg2_b, wrg2_b, degree=[2,1,1])
        call mp_rg2%add_patch(vrg2_a, id1)
        call mp_rg2%add_patch(vrg2_b, id2)
        call mp_rg2%connect(&
            patch_a            = id1,&
            side_a             = "u_max",&
            patch_b            = id2,&
            side_b             = "u_min",&
            continuity         = 2,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk])
        call vg2_a%err%print()
        call vg2_b%err%print()
        call vrg2_a%err%print()
        call vrg2_b%err%print()
        call mp_g2%err%print()
        call mp_rg2%err%print()

        call ut%test(ti)%check(&
            name     = "volume rational G2",&
            res      = mp_rg2%is_valid(multipatch_volume_tol),&
            expected = .true.,&
            msg      = "Volume rational G2 is incorrect.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0012


    subroutine forcad_multipatch_volume_0013(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: vg2_a, vg2_b, vrg2_a, vrg2_b
        type(nurbs_multipatch_volume) :: mp_g2
        type(nurbs_multipatch_volume) :: mp_rg2
        type(multipatch_connection) :: conn
        real(rk) :: wrg2_a(12), wrg2_b(12), xg2_a(12,3), xg2_b(12,3)
        real(rk) :: knot_open(4), knot_quadratic(6), residual, row_res
        real(rk), allocatable :: val(:), dof(:)
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1, id2, row, j

        knot_quadratic = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot_open = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        xg2_a = 0.0_rk
        xg2_b = 0.0_rk
        xg2_a(:,1) = [0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, &
            0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk]
        xg2_a(:,2) = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, &
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_a(:,3) = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, &
            1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_b(:,1) = [1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk, &
            1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk]
        xg2_b(:,2:3) = xg2_a(:,2:3)
        call vg2_a%set(knot_quadratic, knot_open, knot_open, xg2_a, degree=[2,1,1])
        call vg2_b%set(knot_quadratic, knot_open, knot_open, xg2_b, degree=[2,1,1])
        call mp_g2%add_patch(vg2_a, id1)
        call mp_g2%add_patch(vg2_b, id2)
        call mp_g2%connect(&
            patch_a            = id1,&
            side_a             = "u_max",&
            patch_b            = id2,&
            side_b             = "u_min",&
            continuity         = 2,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk])
        conn = mp_g2%get_connection(1)
        call mp_g2%cmp_dof_constraint(rowptr, col, val, geometric=.true.)
        allocate(dof(24))
        dof(1:12) = xg2_a(:,1)
        dof(13:24) = xg2_b(:,1)
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do
        deallocate(dof)
        wrg2_a = 2.0_rk
        wrg2_b = 3.0_rk
        call vrg2_a%set(knot_quadratic, knot_open, knot_open, xg2_a, wrg2_a, degree=[2,1,1])
        call vrg2_b%set(knot_quadratic, knot_open, knot_open, xg2_b, wrg2_b, degree=[2,1,1])
        call mp_rg2%add_patch(vrg2_a, id1)
        call mp_rg2%add_patch(vrg2_b, id2)
        call mp_rg2%connect(&
            patch_a            = id1,&
            side_a             = "u_max",&
            patch_b            = id2,&
            side_b             = "u_min",&
            continuity         = 2,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk])
        call mp_g2%cmp_dof_constraint(rowptr, col, val, geometric=.false.)
        call vg2_a%err%print()
        call vg2_b%err%print()
        call vrg2_a%err%print()
        call vrg2_b%err%print()
        call mp_g2%err%print()
        call mp_rg2%err%print()

        call ut%test(ti)%check(&
            name     = "volume separate C constraints",&
            res      = size(rowptr)-1,&
            expected = 0,&
            msg      = "filtering a G2 connection for C constraints must return no rows",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0013


    subroutine forcad_multipatch_volume_0014(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_volume_tol = 1.0e-10_rk
        type(nurbs_volume) :: vg2_a, vg2_b, vrg2_a, vrg2_b
        type(nurbs_multipatch_volume) :: mp_g2, mp_rg2, mp_not_c2
        type(multipatch_connection) :: conn
        real(rk) :: wrg2_a(12), wrg2_b(12), xg2_a(12,3), xg2_b(12,3)
        real(rk) :: knot_open(4), knot_quadratic(6), residual, row_res
        real(rk), allocatable :: val(:), dof(:)
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1, id2, row, j

        knot_quadratic = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot_open = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        xg2_a = 0.0_rk
        xg2_b = 0.0_rk
        xg2_a(:,1) = [0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, &
            0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk]
        xg2_a(:,2) = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, &
            0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_a(:,3) = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, &
            1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_b(:,1) = [1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk, &
            1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk]
        xg2_b(:,2:3) = xg2_a(:,2:3)
        call vg2_a%set(knot_quadratic, knot_open, knot_open, xg2_a, degree=[2,1,1])
        call vg2_b%set(knot_quadratic, knot_open, knot_open, xg2_b, degree=[2,1,1])
        call mp_g2%add_patch(vg2_a, id1)
        call mp_g2%add_patch(vg2_b, id2)
        call mp_g2%connect(&
            patch_a            = id1,&
            side_a             = "u_max",&
            patch_b            = id2,&
            side_b             = "u_min",&
            continuity         = 2,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk])
        conn = mp_g2%get_connection(1)
        call mp_g2%cmp_dof_constraint(rowptr, col, val, geometric=.true.)
        allocate(dof(24))
        dof(1:12) = xg2_a(:,1)
        dof(13:24) = xg2_b(:,1)
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do
        deallocate(dof)
        wrg2_a = 2.0_rk
        wrg2_b = 3.0_rk
        call vrg2_a%set(knot_quadratic, knot_open, knot_open, xg2_a, wrg2_a, degree=[2,1,1])
        call vrg2_b%set(knot_quadratic, knot_open, knot_open, xg2_b, wrg2_b, degree=[2,1,1])
        call mp_rg2%add_patch(vrg2_a, id1)
        call mp_rg2%add_patch(vrg2_b, id2)
        call mp_rg2%connect(&
            patch_a            = id1,&
            side_a             = "u_max",&
            patch_b            = id2,&
            side_b             = "u_min",&
            continuity         = 2,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk])
        call mp_g2%cmp_dof_constraint(rowptr, col, val, geometric=.false.)

        call mp_not_c2%add_patch(vg2_a, id1)
        call mp_not_c2%add_patch(vg2_b, id2)
        call mp_not_c2%connect(id1, "u_max", id2, "u_min", 2)
        call vg2_a%err%print()
        call vg2_b%err%print()
        call vrg2_a%err%print()
        call vrg2_b%err%print()
        call mp_g2%err%print()
        call mp_rg2%err%print()
        call mp_not_c2%err%print()

        call ut%test(ti)%check(&
            name     = "volume G2 is not C2",&
            res      = mp_not_c2%is_valid(multipatch_volume_tol),&
            expected = .false.,&
            msg      = "The G2 volume interface must not satisfy C2.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0014


    subroutine forcad_multipatch_volume_0015(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1
        type(nurbs_volume) :: v2
        type(nurbs_multipatch_volume) :: mp
        integer, allocatable :: elem(:,:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        elem = mp%cmp_elem(shared=.true.)
        call v1%err%print()
        call v2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume elem shape",&
            res      = shape(elem),&
            expected = [2, 8],&
            msg      = "volume multipatch element shape mismatch",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0015


    subroutine forcad_multipatch_volume_0016(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1, v2
        type(nurbs_multipatch_volume) :: mp
        integer, allocatable :: elem(:,:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        elem = mp%cmp_elem(shared=.true.)
        elem = mp%cmp_elem(shared=.false.)
        call v1%err%print()
        call v2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume unshared elem max",&
            res      = maxval(elem),&
            expected = 16,&
            msg      = "Volume unshared elem max is incorrect.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0016


    subroutine forcad_multipatch_volume_0017(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1, v2
        type(nurbs_multipatch_volume) :: mp
        integer, allocatable :: elem(:,:), patch_id(:), local_id(:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        elem = mp%cmp_elem(shared=.true.)
        elem = mp%cmp_elem(shared=.false.)

        patch_id = mp%cmp_elem_patch()
        local_id = mp%cmp_elem_local()
        call v1%err%print()
        call v2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume elem patch ids",&
            res      = patch_id,&
            expected = [1, 2],&
            msg      = "volume element patch ids mismatch",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0017


    subroutine forcad_multipatch_volume_0018(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1
        type(nurbs_volume) :: v2
        type(nurbs_multipatch_volume) :: mp
        integer, allocatable :: elem(:,:), patch_id(:), local_id(:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        elem = mp%cmp_elem(shared=.true.)
        elem = mp%cmp_elem(shared=.false.)

        patch_id = mp%cmp_elem_patch()
        local_id = mp%cmp_elem_local()
        call v1%err%print()
        call v2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume elem local ids",&
            res      = local_id,&
            expected = [1, 1],&
            msg      = "volume local element ids mismatch",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0018


    subroutine forcad_multipatch_volume_0019(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1
        type(nurbs_volume) :: v2
        type(nurbs_multipatch_volume) :: mp
        type(multipatch_connection) :: conn
        integer, allocatable :: elem(:,:), patch_id(:), local_id(:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        elem = mp%cmp_elem(shared=.true.)
        elem = mp%cmp_elem(shared=.false.)

        patch_id = mp%cmp_elem_patch()
        local_id = mp%cmp_elem_local()

        conn = mp%get_connection(1)
        call v1%err%print()
        call v2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume connection continuity",&
            res      = conn%get_continuity(),&
            expected = 0,&
            msg      = "volume connection continuity metadata mismatch",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0019


    subroutine forcad_multipatch_volume_0020(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1
        type(nurbs_volume) :: v2
        type(nurbs_multipatch_volume) :: mp
        type(multipatch_connection) :: conn
        integer, allocatable :: elem(:,:), patch_id(:), local_id(:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        elem = mp%cmp_elem(shared=.true.)
        elem = mp%cmp_elem(shared=.false.)

        patch_id = mp%cmp_elem_patch()
        local_id = mp%cmp_elem_local()

        conn = mp%get_connection(1)
        call v1%err%print()
        call v2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume connection flip",&
            res      = conn%is_flipped(1),&
            expected = .false.,&
            msg      = "volume connection flip metadata mismatch",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0020


    subroutine forcad_multipatch_volume_0021(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1
        type(nurbs_volume) :: v2
        type(nurbs_multipatch_volume) :: mp
        type(multipatch_connection) :: conn
        integer, allocatable :: elem(:,:), patch_id(:), local_id(:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        elem = mp%cmp_elem(shared=.true.)
        elem = mp%cmp_elem(shared=.false.)

        patch_id = mp%cmp_elem_patch()
        local_id = mp%cmp_elem_local()

        conn = mp%get_connection(1)
        call v1%err%print()
        call v2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume connection swap",&
            res      = conn%is_swapped(),&
            expected = .false.,&
            msg      = "volume connection swap metadata mismatch",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0021


    subroutine forcad_multipatch_volume_0022(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_volume_tol = 1.0e-10_rk
        type(nurbs_volume) :: v1, v2, patch
        type(nurbs_multipatch_volume) :: mp
        type(multipatch_connection) :: conn
        real(rk), allocatable :: Xc(:,:)
        integer, allocatable :: elem(:,:), patch_id(:), local_id(:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        elem = mp%cmp_elem(shared=.true.)
        elem = mp%cmp_elem(shared=.false.)

        patch_id = mp%cmp_elem_patch()
        local_id = mp%cmp_elem_local()

        conn = mp%get_connection(1)

        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        call v1%err%print()
        call v2%err%print()
        call patch%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume get patch",&
            res      = Xc(1,:),&
            expected = [1.0_rk, 0.0_rk, 0.0_rk],&
            tol      = multipatch_volume_tol,&
            msg      = "volume get_patch returned wrong patch",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0022


    subroutine forcad_multipatch_volume_0023(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1, v2, patch
        type(nurbs_multipatch_volume) :: mp
        type(multipatch_connection) :: conn
        real(rk), allocatable :: Xc(:,:)
        integer, allocatable :: elem(:,:), patch_id(:), local_id(:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        elem = mp%cmp_elem(shared=.true.)
        elem = mp%cmp_elem(shared=.false.)

        patch_id = mp%cmp_elem_patch()
        local_id = mp%cmp_elem_local()

        conn = mp%get_connection(1)

        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        call v1%err%print()
        call v2%err%print()
        call patch%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume rational flag",&
            res      = mp%is_rational(),&
            expected = .false.,&
            msg      = "volume multipatch should report non-rational patches",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0023


    subroutine forcad_multipatch_volume_0024(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_volume_tol = 1.0e-10_rk
        type(nurbs_volume) :: v1
        type(nurbs_volume) :: v2
        type(nurbs_volume) :: vr1
        type(nurbs_volume) :: vr2
        type(nurbs_multipatch_volume) :: mp_rational
        real(rk) :: wrv1(8)
        real(rk) :: wrv2(8)
        real(rk) :: knot_open(4)
        integer, allocatable :: map(:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        knot_open = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        wrv1 = [1.0_rk, 2.0_rk, 4.0_rk, 3.0_rk, 6.0_rk, 5.0_rk, 8.0_rk, 7.0_rk]
        wrv2 = [8.0_rk, 11.0_rk, 12.0_rk, 13.0_rk, 20.0_rk, 17.0_rk, 28.0_rk, 19.0_rk]
        call vr1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2], wrv1)
        call vr2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2], wrv2)
        call vr2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_rational%add_patch(vr1, id1)
        call mp_rational%add_patch(vr2, id2)
        call mp_rational%connect(id1, "u_max", id2, "u_min", 0)
        map = mp_rational%cmp_dof_map()
        call v1%err%print()
        call v2%err%print()
        call vr1%err%print()
        call vr2%err%print()
        call mp_rational%err%print()

        call ut%test(ti)%check(&
            name     = "volume constant rational basis sharing",&
            res      = mp_rational%is_valid(multipatch_volume_tol) .and. maxval(map) == 12,&
            expected = .true.,&
            msg      = "Constant projective scaling must preserve shared face DOFs.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0024


    subroutine forcad_multipatch_volume_0025(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_volume_tol = 1.0e-10_rk
        type(nurbs_volume) :: v1
        type(nurbs_volume) :: v2
        type(nurbs_volume) :: vr1
        type(nurbs_volume) :: vr2
        type(nurbs_volume) :: v_scaled_trace
        type(nurbs_multipatch_volume) :: mp_rational
        type(nurbs_multipatch_volume) :: mp_scaled_trace
        real(rk), allocatable :: Xc(:,:)
        real(rk) :: wrv1(8), wrv2(8), knot_open(4), knot_trace_scaled(4), knot_trace_bad(4)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        knot_open = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        wrv1 = [1.0_rk, 2.0_rk, 4.0_rk, 3.0_rk, 6.0_rk, 5.0_rk, 8.0_rk, 7.0_rk]
        wrv2 = [8.0_rk, 11.0_rk, 12.0_rk, 13.0_rk, 20.0_rk, 17.0_rk, 28.0_rk, 19.0_rk]
        call vr1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2], wrv1)
        call vr2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2], wrv2)
        call vr2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_rational%add_patch(vr1, id1)
        call mp_rational%add_patch(vr2, id2)
        call mp_rational%connect(id1, "u_max", id2, "u_min", 0)

        knot_trace_scaled = [0.0_rk, 0.0_rk, 2.0_rk, 2.0_rk]
        knot_trace_bad = [0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk]
        Xc = v2%get_Xc()
        call v_scaled_trace%set(knot_open, knot_trace_scaled, knot_open, Xc, degree=[1, 1, 1])
        call mp_scaled_trace%add_patch(v1, id1)
        call mp_scaled_trace%add_patch(v_scaled_trace, id2)
        call mp_scaled_trace%connect(id1, "u_max", id2, "u_min", 0)
        call v1%err%print()
        call v2%err%print()
        call vr1%err%print()
        call vr2%err%print()
        call v_scaled_trace%err%print()
        call mp_rational%err%print()
        call mp_scaled_trace%err%print()

        call ut%test(ti)%check(&
            name     = "volume affine-scaled trace knots",&
            res      = mp_scaled_trace%is_valid(multipatch_volume_tol),&
            expected = .true.,&
            msg      = "volume C0 trace should accept affine-equivalent tangential knots",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0025


    subroutine forcad_multipatch_volume_0026(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1, v2, vr1, vr2, v_scaled_trace, v_bad_trace
        type(nurbs_multipatch_volume) :: mp_rational, mp_scaled_trace, mp_bad_trace
        real(rk), allocatable :: Xc(:,:)
        real(rk) :: wrv1(8), wrv2(8), knot_open(4), knot_trace_scaled(4), knot_trace_bad(4)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        knot_open = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        wrv1 = [1.0_rk, 2.0_rk, 4.0_rk, 3.0_rk, 6.0_rk, 5.0_rk, 8.0_rk, 7.0_rk]
        wrv2 = [8.0_rk, 11.0_rk, 12.0_rk, 13.0_rk, 20.0_rk, 17.0_rk, 28.0_rk, 19.0_rk]
        call vr1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2], wrv1)
        call vr2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2], wrv2)
        call vr2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_rational%add_patch(vr1, id1)
        call mp_rational%add_patch(vr2, id2)
        call mp_rational%connect(id1, "u_max", id2, "u_min", 0)

        knot_trace_scaled = [0.0_rk, 0.0_rk, 2.0_rk, 2.0_rk]
        knot_trace_bad = [0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk]
        Xc = v2%get_Xc()
        call v_scaled_trace%set(knot_open, knot_trace_scaled, knot_open, Xc, degree=[1, 1, 1])
        call mp_scaled_trace%add_patch(v1, id1)
        call mp_scaled_trace%add_patch(v_scaled_trace, id2)
        call mp_scaled_trace%connect(id1, "u_max", id2, "u_min", 0)

        call v_bad_trace%set(knot_open, knot_trace_bad, knot_open, Xc, degree=[1, 1, 1])
        call mp_bad_trace%add_patch(v1, id1)
        call mp_bad_trace%add_patch(v_bad_trace, id2)
        call mp_bad_trace%connect(id1, "u_max", id2, "u_min", 0)
        call v1%err%print()
        call v2%err%print()
        call vr1%err%print()
        call vr2%err%print()
        call v_scaled_trace%err%print()
        call v_bad_trace%err%print()
        call mp_rational%err%print()
        call mp_scaled_trace%err%print()
        call mp_bad_trace%err%print()

        call ut%test(ti)%check(&
            name     = "volume incompatible trace knots err",&
            res      = mp_bad_trace%err%ok,&
            expected = .false.,&
            msg      = "Volume incompatible trace knots err is incorrect.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0026


    subroutine forcad_multipatch_volume_0027(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_volume_tol = 1.0e-10_rk
        type(nurbs_volume) :: v1, v2, patch
        type(nurbs_multipatch_volume) :: mp
        real(rk), allocatable :: Xc(:,:), Xg(:,:)
        integer :: id1, id2, file_unit

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        call mp%create(res1=3, res2=3, res3=3)
        call mp%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%translate_Xg([0.0_rk, 0.0_rk, 1.0_rk])
        call mp%rotate_Xc(0.0_rk, 0.0_rk, 0.0_rk)
        call mp%rotate_Xg(0.0_rk, 0.0_rk, 0.0_rk)
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_volume_Xc_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_volume_Xg_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_volume_Xth_1.vtk", status="replace")
        close(file_unit, status="delete")
        call mp%export_Xc("vtk/forcad_test_multipatch_volume")
        call mp%export_Xg("vtk/forcad_test_multipatch_volume")
        call mp%export_Xth_in_Xg("vtk/forcad_test_multipatch_volume", res=3)
        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        Xg = patch%get_Xg()
        call v1%err%print()
        call v2%err%print()
        call patch%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume translate Xc",&
            res      = Xc(1,:),&
            expected = [2.0_rk, 0.0_rk, 0.0_rk],&
            tol      = multipatch_volume_tol,&
            msg      = "volume multipatch translate_Xc failed",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0027


    subroutine forcad_multipatch_volume_0028(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_volume_tol = 1.0e-10_rk
        type(nurbs_volume) :: v1, v2, patch
        type(nurbs_multipatch_volume) :: mp
        real(rk), allocatable :: Xc(:,:), Xg(:,:)
        integer :: id1, id2, file_unit

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        call mp%create(res1=3, res2=3, res3=3)
        call mp%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%translate_Xg([0.0_rk, 0.0_rk, 1.0_rk])
        call mp%rotate_Xc(0.0_rk, 0.0_rk, 0.0_rk)
        call mp%rotate_Xg(0.0_rk, 0.0_rk, 0.0_rk)
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_volume_Xc_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_volume_Xg_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_volume_Xth_1.vtk", status="replace")
        close(file_unit, status="delete")
        call mp%export_Xc("vtk/forcad_test_multipatch_volume")
        call mp%export_Xg("vtk/forcad_test_multipatch_volume")
        call mp%export_Xth_in_Xg("vtk/forcad_test_multipatch_volume", res=3)
        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        Xg = patch%get_Xg()
        call v1%err%print()
        call v2%err%print()
        call patch%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume translate Xg",&
            res      = Xg(1,3),&
            expected = 1.0_rk,&
            tol      = multipatch_volume_tol,&
            msg      = "volume multipatch translate_Xg failed",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0028


    subroutine forcad_multipatch_volume_0029(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1, v2, patch
        type(nurbs_multipatch_volume) :: mp
        real(rk), allocatable :: Xc(:,:), Xg(:,:)
        integer :: id1, id2, file_unit
        logical :: file_exists

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        call mp%create(res1=3, res2=3, res3=3)
        call mp%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%translate_Xg([0.0_rk, 0.0_rk, 1.0_rk])
        call mp%rotate_Xc(0.0_rk, 0.0_rk, 0.0_rk)
        call mp%rotate_Xg(0.0_rk, 0.0_rk, 0.0_rk)
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_volume_Xc_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_volume_Xg_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_volume_Xth_1.vtk", status="replace")
        close(file_unit, status="delete")
        call mp%export_Xc("vtk/forcad_test_multipatch_volume")
        call mp%export_Xg("vtk/forcad_test_multipatch_volume")
        call mp%export_Xth_in_Xg("vtk/forcad_test_multipatch_volume", res=3)
        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        Xg = patch%get_Xg()
        inquire(file="vtk/forcad_test_multipatch_volume_Xc_1.vtk", exist=file_exists)
        call v1%err%print()
        call v2%err%print()
        call patch%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume export Xc",&
            res      = file_exists,&
            expected = .true.,&
            msg      = "volume multipatch export_Xc did not create the first VTK file",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0029


    subroutine forcad_multipatch_volume_0030(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1, v2, patch
        type(nurbs_multipatch_volume) :: mp
        real(rk), allocatable :: Xc(:,:), Xg(:,:)
        integer :: id1, id2, file_unit
        logical :: file_exists

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        call mp%create(res1=3, res2=3, res3=3)
        call mp%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%translate_Xg([0.0_rk, 0.0_rk, 1.0_rk])
        call mp%rotate_Xc(0.0_rk, 0.0_rk, 0.0_rk)
        call mp%rotate_Xg(0.0_rk, 0.0_rk, 0.0_rk)
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_volume_Xc_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_volume_Xg_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_volume_Xth_1.vtk", status="replace")
        close(file_unit, status="delete")
        call mp%export_Xc("vtk/forcad_test_multipatch_volume")
        call mp%export_Xg("vtk/forcad_test_multipatch_volume")
        call mp%export_Xth_in_Xg("vtk/forcad_test_multipatch_volume", res=3)
        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        Xg = patch%get_Xg()
        inquire(file="vtk/forcad_test_multipatch_volume_Xc_1.vtk", exist=file_exists)
        inquire(file="vtk/forcad_test_multipatch_volume_Xg_1.vtk", exist=file_exists)
        call v1%err%print()
        call v2%err%print()
        call patch%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume export Xg",&
            res      = file_exists,&
            expected = .true.,&
            msg      = "volume multipatch export_Xg did not create the first VTK file",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0030


    subroutine forcad_multipatch_volume_0031(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1, v2, patch
        type(nurbs_multipatch_volume) :: mp
        real(rk), allocatable :: Xc(:,:), Xg(:,:)
        integer :: id1, id2, file_unit
        logical :: file_exists

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        call mp%create(res1=3, res2=3, res3=3)
        call mp%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%translate_Xg([0.0_rk, 0.0_rk, 1.0_rk])
        call mp%rotate_Xc(0.0_rk, 0.0_rk, 0.0_rk)
        call mp%rotate_Xg(0.0_rk, 0.0_rk, 0.0_rk)
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_volume_Xc_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_volume_Xg_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_volume_Xth_1.vtk", status="replace")
        close(file_unit, status="delete")
        call mp%export_Xc("vtk/forcad_test_multipatch_volume")
        call mp%export_Xg("vtk/forcad_test_multipatch_volume")
        call mp%export_Xth_in_Xg("vtk/forcad_test_multipatch_volume", res=3)
        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        Xg = patch%get_Xg()
        inquire(file="vtk/forcad_test_multipatch_volume_Xc_1.vtk", exist=file_exists)
        inquire(file="vtk/forcad_test_multipatch_volume_Xg_1.vtk", exist=file_exists)
        inquire(file="vtk/forcad_test_multipatch_volume_Xth_1.vtk", exist=file_exists)
        call v1%err%print()
        call v2%err%print()
        call patch%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume export Xth in Xg",&
            res      = file_exists,&
            expected = .true.,&
            msg      = "Volume export Xth in Xg is incorrect.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0031


    subroutine forcad_multipatch_volume_0032(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1
        type(nurbs_volume) :: v2
        type(nurbs_multipatch_volume) :: mp_discontinuous
        integer, allocatable :: map(:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(v1, id1)
        call mp_discontinuous%add_patch(v2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()
        call v1%err%print()
        call v2%err%print()
        call mp_discontinuous%err%print()

        call ut%test(ti)%check(&
            name     = "volume discontinuous map",&
            res      = maxval(map),&
            expected = 16,&
            msg      = "discontinuous volume interface should not share dofs",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0032


    subroutine forcad_multipatch_volume_0033(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1
        type(nurbs_volume) :: v2
        type(nurbs_volume) :: v_mismatch
        type(nurbs_multipatch_volume) :: mp_discontinuous, mp_discontinuous_mismatch
        integer, allocatable :: map(:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(v1, id1)
        call mp_discontinuous%add_patch(v2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()

        call v_mismatch%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 3, 2])
        call mp_discontinuous_mismatch%add_patch(v1, id1)
        call mp_discontinuous_mismatch%add_patch(v_mismatch, id2)
        call mp_discontinuous_mismatch%connect(id1, "u_max", id2, "u_min", -1)
        call v1%err%print()
        call v2%err%print()
        call v_mismatch%err%print()
        call mp_discontinuous%err%print()
        call mp_discontinuous_mismatch%err%print()

        call ut%test(ti)%check(&
            name     = "volume discontinuous mismatch err",&
            res      = mp_discontinuous_mismatch%err%ok,&
            expected = .true.,&
            msg      = "discontinuous volume faces must allow mismatched discretizations",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0033


    subroutine forcad_multipatch_volume_0034(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1
        type(nurbs_volume) :: v2
        type(nurbs_volume) :: v_mismatch
        type(nurbs_multipatch_volume) :: mp_discontinuous, mp_discontinuous_mismatch
        integer, allocatable :: map(:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(v1, id1)
        call mp_discontinuous%add_patch(v2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()

        call v_mismatch%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 3, 2])
        call mp_discontinuous_mismatch%add_patch(v1, id1)
        call mp_discontinuous_mismatch%add_patch(v_mismatch, id2)
        call mp_discontinuous_mismatch%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous_mismatch%cmp_dof_map()
        call v1%err%print()
        call v2%err%print()
        call v_mismatch%err%print()
        call mp_discontinuous%err%print()
        call mp_discontinuous_mismatch%err%print()

        call ut%test(ti)%check(&
            name     = "volume discontinuous mismatch map",&
            res      = maxval(map),&
            expected = 20,&
            msg      = "discontinuous mismatched volume interface should not share dofs",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0034


    subroutine forcad_multipatch_volume_0035(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_volume_tol = 1.0e-10_rk
        type(nurbs_volume) :: v1, v2, v_bad, v_mismatch
        type(nurbs_multipatch_volume) :: mp_discontinuous, mp_discontinuous_mismatch
        type(nurbs_multipatch_volume) :: mp_bad_geom
        integer, allocatable :: map(:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(v1, id1)
        call mp_discontinuous%add_patch(v2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()

        call v_mismatch%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 3, 2])
        call mp_discontinuous_mismatch%add_patch(v1, id1)
        call mp_discontinuous_mismatch%add_patch(v_mismatch, id2)
        call mp_discontinuous_mismatch%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous_mismatch%cmp_dof_map()

        call v_bad%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v_bad%translate_Xc([2.0_rk, 0.0_rk, 0.0_rk])
        call mp_bad_geom%add_patch(v1, id1)
        call mp_bad_geom%add_patch(v_bad, id2)
        call mp_bad_geom%connect(id1, "u_max", id2, "u_min", 0)
        call v1%err%print()
        call v2%err%print()
        call v_bad%err%print()
        call v_mismatch%err%print()
        call mp_discontinuous%err%print()
        call mp_discontinuous_mismatch%err%print()
        call mp_bad_geom%err%print()

        call ut%test(ti)%check(&
            name     = "volume invalid geometry",&
            res      = mp_bad_geom%is_valid(multipatch_volume_tol),&
            expected = .false.,&
            msg      = "volume multipatch should reject mismatched interface geometry",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0035


    subroutine forcad_multipatch_volume_0036(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1, v2, v_bad, v_mismatch
        type(nurbs_multipatch_volume) :: mp_discontinuous, mp_discontinuous_mismatch
        type(nurbs_multipatch_volume) :: mp_bad_geom
        type(nurbs_multipatch_volume) :: mp_bad_side
        integer, allocatable :: map(:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(v1, id1)
        call mp_discontinuous%add_patch(v2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()

        call v_mismatch%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 3, 2])
        call mp_discontinuous_mismatch%add_patch(v1, id1)
        call mp_discontinuous_mismatch%add_patch(v_mismatch, id2)
        call mp_discontinuous_mismatch%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous_mismatch%cmp_dof_map()

        call v_bad%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v_bad%translate_Xc([2.0_rk, 0.0_rk, 0.0_rk])
        call mp_bad_geom%add_patch(v1, id1)
        call mp_bad_geom%add_patch(v_bad, id2)
        call mp_bad_geom%connect(id1, "u_max", id2, "u_min", 0)

        call mp_bad_side%add_patch(v1, id1)
        call mp_bad_side%connect(id1, "bad_side", id1, "u_min", 0)
        call v1%err%print()
        call v2%err%print()
        call v_bad%err%print()
        call v_mismatch%err%print()
        call mp_discontinuous%err%print()
        call mp_discontinuous_mismatch%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_side%err%print()

        call ut%test(ti)%check(&
            name     = "volume invalid side err",&
            res      = mp_bad_side%err%ok,&
            expected = .false.,&
            msg      = "invalid volume side should set multipatch error state",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0036


    subroutine forcad_multipatch_volume_0037(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1, v2, v_bad, v_mismatch
        type(nurbs_multipatch_volume) :: mp_discontinuous, mp_discontinuous_mismatch
        type(nurbs_multipatch_volume) :: mp_bad_geom
        type(nurbs_multipatch_volume) :: mp_bad_side
        integer, allocatable :: map(:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(v1, id1)
        call mp_discontinuous%add_patch(v2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()

        call v_mismatch%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 3, 2])
        call mp_discontinuous_mismatch%add_patch(v1, id1)
        call mp_discontinuous_mismatch%add_patch(v_mismatch, id2)
        call mp_discontinuous_mismatch%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous_mismatch%cmp_dof_map()

        call v_bad%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v_bad%translate_Xc([2.0_rk, 0.0_rk, 0.0_rk])
        call mp_bad_geom%add_patch(v1, id1)
        call mp_bad_geom%add_patch(v_bad, id2)
        call mp_bad_geom%connect(id1, "u_max", id2, "u_min", 0)

        call mp_bad_side%add_patch(v1, id1)
        call mp_bad_side%connect(id1, "bad_side", id1, "u_min", 0)
        call v1%err%print()
        call v2%err%print()
        call v_bad%err%print()
        call v_mismatch%err%print()
        call mp_discontinuous%err%print()
        call mp_discontinuous_mismatch%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_side%err%print()

        call ut%test(ti)%check(&
            name     = "volume invalid side no connection",&
            res      = mp_bad_side%get_nconnection(),&
            expected = 0,&
            msg      = "invalid volume side should not append a connection",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0037


    subroutine forcad_multipatch_volume_0038(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1, v2, v_bad, v_mismatch
        type(nurbs_multipatch_volume) :: mp_discontinuous, mp_discontinuous_mismatch
        type(nurbs_multipatch_volume) :: mp_bad_geom
        type(nurbs_multipatch_volume) :: mp_bad_side
        integer, allocatable :: map(:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(v1, id1)
        call mp_discontinuous%add_patch(v2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()

        call v_mismatch%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 3, 2])
        call mp_discontinuous_mismatch%add_patch(v1, id1)
        call mp_discontinuous_mismatch%add_patch(v_mismatch, id2)
        call mp_discontinuous_mismatch%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous_mismatch%cmp_dof_map()

        call v_bad%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v_bad%translate_Xc([2.0_rk, 0.0_rk, 0.0_rk])
        call mp_bad_geom%add_patch(v1, id1)
        call mp_bad_geom%add_patch(v_bad, id2)
        call mp_bad_geom%connect(id1, "u_max", id2, "u_min", 0)

        call mp_bad_side%add_patch(v1, id1)
        call mp_bad_side%connect(id1, "bad_side", id1, "u_min", 0)
        call mp_bad_side%add_patch(v2, id2)
        call v1%err%print()
        call v2%err%print()
        call v_bad%err%print()
        call v_mismatch%err%print()
        call mp_discontinuous%err%print()
        call mp_discontinuous_mismatch%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_side%err%print()

        call ut%test(ti)%check(&
            name     = "volume add blocked on err",&
            res      = id2,&
            expected = 0,&
            msg      = "Volume add blocked on err is incorrect.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0038


    subroutine forcad_multipatch_volume_0039(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1, v2, v_bad, v_mismatch
        type(nurbs_multipatch_volume) :: mp_discontinuous, mp_discontinuous_mismatch
        type(nurbs_multipatch_volume) :: mp_bad_geom
        type(nurbs_multipatch_volume) :: mp_bad_side
        integer, allocatable :: map(:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(v1, id1)
        call mp_discontinuous%add_patch(v2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()

        call v_mismatch%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 3, 2])
        call mp_discontinuous_mismatch%add_patch(v1, id1)
        call mp_discontinuous_mismatch%add_patch(v_mismatch, id2)
        call mp_discontinuous_mismatch%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous_mismatch%cmp_dof_map()

        call v_bad%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v_bad%translate_Xc([2.0_rk, 0.0_rk, 0.0_rk])
        call mp_bad_geom%add_patch(v1, id1)
        call mp_bad_geom%add_patch(v_bad, id2)
        call mp_bad_geom%connect(id1, "u_max", id2, "u_min", 0)

        call mp_bad_side%add_patch(v1, id1)
        call mp_bad_side%connect(id1, "bad_side", id1, "u_min", 0)
        call mp_bad_side%add_patch(v2, id2)
        call mp_bad_side%err%print()
        call mp_bad_side%err%reset()
        call mp_bad_side%add_patch(v2, id2)
        call v1%err%print()
        call v2%err%print()
        call v_bad%err%print()
        call v_mismatch%err%print()
        call mp_discontinuous%err%print()
        call mp_discontinuous_mismatch%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_side%err%print()

        call ut%test(ti)%check(&
            name     = "volume add after reset",&
            res      = id2,&
            expected = 2,&
            msg      = "debug reset should allow volume multipatch operations again",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0039


    subroutine forcad_multipatch_volume_0040(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1, v2, v_bad, v_mismatch
        type(nurbs_multipatch_volume) :: mp_discontinuous, mp_discontinuous_mismatch
        type(nurbs_multipatch_volume) :: mp_bad_geom, mp_bad_side, mp_bad_continuity
        integer, allocatable :: map(:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(v1, id1)
        call mp_discontinuous%add_patch(v2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()

        call v_mismatch%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 3, 2])
        call mp_discontinuous_mismatch%add_patch(v1, id1)
        call mp_discontinuous_mismatch%add_patch(v_mismatch, id2)
        call mp_discontinuous_mismatch%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous_mismatch%cmp_dof_map()

        call v_bad%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v_bad%translate_Xc([2.0_rk, 0.0_rk, 0.0_rk])
        call mp_bad_geom%add_patch(v1, id1)
        call mp_bad_geom%add_patch(v_bad, id2)
        call mp_bad_geom%connect(id1, "u_max", id2, "u_min", 0)

        call mp_bad_side%add_patch(v1, id1)
        call mp_bad_side%connect(id1, "bad_side", id1, "u_min", 0)
        call mp_bad_side%add_patch(v2, id2)
        call mp_bad_side%err%print()
        call mp_bad_side%err%reset()
        call mp_bad_side%add_patch(v2, id2)

        call mp_bad_continuity%add_patch(v1, id1)
        call mp_bad_continuity%add_patch(v2, id2)
        call mp_bad_continuity%connect(id1, "u_max", id2, "u_min", 2)
        call v1%err%print()
        call v2%err%print()
        call v_bad%err%print()
        call v_mismatch%err%print()
        call mp_discontinuous%err%print()
        call mp_discontinuous_mismatch%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_side%err%print()
        call mp_bad_continuity%err%print()

        call ut%test(ti)%check(&
            name     = "volume invalid continuity err",&
            res      = mp_bad_continuity%err%ok,&
            expected = .false.,&
            msg      = "invalid volume continuity should set multipatch error state",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0040


    subroutine forcad_multipatch_volume_0041(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1, v2, v_bad, v_mismatch
        type(nurbs_multipatch_volume) :: mp_discontinuous, mp_discontinuous_mismatch
        type(nurbs_multipatch_volume) :: mp_bad_geom, mp_bad_side, mp_bad_continuity
        integer, allocatable :: map(:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(v1, id1)
        call mp_discontinuous%add_patch(v2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()

        call v_mismatch%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 3, 2])
        call mp_discontinuous_mismatch%add_patch(v1, id1)
        call mp_discontinuous_mismatch%add_patch(v_mismatch, id2)
        call mp_discontinuous_mismatch%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous_mismatch%cmp_dof_map()

        call v_bad%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v_bad%translate_Xc([2.0_rk, 0.0_rk, 0.0_rk])
        call mp_bad_geom%add_patch(v1, id1)
        call mp_bad_geom%add_patch(v_bad, id2)
        call mp_bad_geom%connect(id1, "u_max", id2, "u_min", 0)

        call mp_bad_side%add_patch(v1, id1)
        call mp_bad_side%connect(id1, "bad_side", id1, "u_min", 0)
        call mp_bad_side%add_patch(v2, id2)
        call mp_bad_side%err%print()
        call mp_bad_side%err%reset()
        call mp_bad_side%add_patch(v2, id2)

        call mp_bad_continuity%add_patch(v1, id1)
        call mp_bad_continuity%add_patch(v2, id2)
        call mp_bad_continuity%connect(id1, "u_max", id2, "u_min", 2)
        call v1%err%print()
        call v2%err%print()
        call v_bad%err%print()
        call v_mismatch%err%print()
        call mp_discontinuous%err%print()
        call mp_discontinuous_mismatch%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_side%err%print()
        call mp_bad_continuity%err%print()

        call ut%test(ti)%check(&
            name     = "volume invalid continuity no connection",&
            res      = mp_bad_continuity%get_nconnection(),&
            expected = 0,&
            msg      = "invalid volume continuity should not append a connection",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0041


    subroutine forcad_multipatch_volume_0042(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_volume_tol = 1.0e-10_rk
        type(nurbs_volume) :: v1, v2, v_bad, v_mismatch
        type(nurbs_multipatch_volume) :: mp_discontinuous, mp_discontinuous_mismatch
        type(nurbs_multipatch_volume) :: mp_bad_geom, mp_bad_side, mp_bad_continuity
        integer, allocatable :: map(:)
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(v1, id1)
        call mp_discontinuous%add_patch(v2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()

        call v_mismatch%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 3, 2])
        call mp_discontinuous_mismatch%add_patch(v1, id1)
        call mp_discontinuous_mismatch%add_patch(v_mismatch, id2)
        call mp_discontinuous_mismatch%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous_mismatch%cmp_dof_map()

        call v_bad%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v_bad%translate_Xc([2.0_rk, 0.0_rk, 0.0_rk])
        call mp_bad_geom%add_patch(v1, id1)
        call mp_bad_geom%add_patch(v_bad, id2)
        call mp_bad_geom%connect(id1, "u_max", id2, "u_min", 0)

        call mp_bad_side%add_patch(v1, id1)
        call mp_bad_side%connect(id1, "bad_side", id1, "u_min", 0)
        call mp_bad_side%add_patch(v2, id2)
        call mp_bad_side%err%print()
        call mp_bad_side%err%reset()
        call mp_bad_side%add_patch(v2, id2)

        call mp_bad_continuity%add_patch(v1, id1)
        call mp_bad_continuity%add_patch(v2, id2)
        call mp_bad_continuity%connect(id1, "u_max", id2, "u_min", 2)
        call v1%err%print()
        call v2%err%print()
        call v_bad%err%print()
        call v_mismatch%err%print()
        call mp_discontinuous%err%print()
        call mp_discontinuous_mismatch%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_side%err%print()
        call mp_bad_continuity%err%print()

        call ut%test(ti)%check(&
            name     = "volume invalid continuity",&
            res      = mp_bad_continuity%is_valid(multipatch_volume_tol),&
            expected = .false.,&
            msg      = "volume multipatch should reject continuity above patch degree",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0042


    subroutine forcad_multipatch_volume_0043(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: invalid_patch
        type(nurbs_volume) :: v1
        type(nurbs_multipatch_volume) :: invalid_patch_mp
        integer :: id1

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call invalid_patch%err%set(&
            code       = 901,&
            severity   = 1,&
            category   = "test_forcad_multipatch_volume",&
            message    = "Intentional invalid patch state.",&
            location   = "dynamic storage test",&
            suggestion = "Used only to verify add_patch validation.")
        call invalid_patch_mp%add_patch(invalid_patch, id1)
        call invalid_patch%err%print()
        call v1%err%print()
        call invalid_patch_mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume rejects invalid patch",&
            res      = [invalid_patch_mp%get_npatch(), id1],&
            expected = [0, 0],&
            msg      = "An invalid volume patch must not be appended or receive an id.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0043


    subroutine forcad_multipatch_volume_0044(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: invalid_patch
        type(nurbs_volume) :: v1
        type(nurbs_multipatch_volume) :: capacity_mp, invalid_patch_mp
        integer :: i, ids(18), id1

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call invalid_patch%err%set(&
            code       = 901,&
            severity   = 1,&
            category   = "test_forcad_multipatch_volume",&
            message    = "Intentional invalid patch state.",&
            location   = "dynamic storage test",&
            suggestion = "Used only to verify add_patch validation.")
        call invalid_patch_mp%add_patch(invalid_patch, id1)

        do i = 1, size(ids)
            call capacity_mp%add_patch(v1, ids(i))
        end do
        do i = 1, size(ids)-1
            call capacity_mp%connect(ids(i), "u_max", ids(i+1), "u_max", -1)
        end do
        call invalid_patch%err%print()
        call v1%err%print()
        call capacity_mp%err%print()
        call invalid_patch_mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume dynamic storage growth",&
            res      = [capacity_mp%get_npatch(), capacity_mp%get_nconnection()],&
            expected = [18, 17],&
            msg      = "Volume dynamic storage growth is incorrect.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0044


    subroutine forcad_multipatch_volume_0045(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1, v2
        type(nurbs_multipatch_volume) :: mp, mp_copy
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        mp_copy = mp
        call mp%finalize()
        call v1%err%print()
        call v2%err%print()
        call mp%err%print()
        call mp_copy%err%print()

        call ut%test(ti)%check(&
            name     = "volume finalize",&
            res      = mp%get_npatch(),&
            expected = 0,&
            msg      = "volume multipatch finalize failed",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0045


    subroutine forcad_multipatch_volume_0046(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: v1, v2
        type(nurbs_multipatch_volume) :: mp, mp_copy
        integer :: id1, id2

        call v1%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call v2%set_hexahedron([1.0_rk, 1.0_rk, 1.0_rk], [2, 2, 2])
        call mp%add_patch(v1, id1)
        call mp%add_patch(v2, id2)
        mp_copy = mp
        call mp%finalize()
        call v1%err%print()
        call v2%err%print()
        call mp%err%print()
        call mp_copy%err%print()

        call ut%test(ti)%check(&
            name     = "volume multipatch assignment",&
            res      = mp_copy%get_npatch(),&
            expected = 2,&
            msg      = "Volume multipatch assignment is incorrect.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0046


    subroutine forcad_multipatch_volume_0047(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 1.0e-12_rk
        real(rk), parameter :: knot_normal(7) = [0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk]
        real(rk), parameter :: knot_tangent(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xca(16,3) = reshape([&
            -4.0_rk,-2.0_rk,0.0_rk,2.0_rk,-4.0_rk,-2.0_rk,0.0_rk,2.0_rk,&
            -4.0_rk,-2.0_rk,0.0_rk,2.0_rk,-4.0_rk,-2.0_rk,0.0_rk,2.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,&
            1.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk], [16,3])
        real(rk), parameter :: Xcb(16,3) = reshape([&
            -1.0_rk,4.0_rk,5.0_rk,6.0_rk,-1.0_rk,4.0_rk,5.0_rk,6.0_rk,&
            -1.0_rk,4.0_rk,5.0_rk,6.0_rk,-1.0_rk,4.0_rk,5.0_rk,6.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,&
            1.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk], [16,3])
        real(rk), parameter :: Wa(16) = [2.0_rk,2.0_rk,1.0_rk,3.0_rk,2.0_rk,2.0_rk,1.0_rk,3.0_rk,&
            2.0_rk,2.0_rk,1.0_rk,3.0_rk,2.0_rk,2.0_rk,1.0_rk,3.0_rk]
        real(rk), parameter :: Wb(16) = [2.0_rk,2.0_rk,3.0_rk,4.0_rk,2.0_rk,2.0_rk,3.0_rk,4.0_rk,&
            2.0_rk,2.0_rk,3.0_rk,4.0_rk,2.0_rk,2.0_rk,3.0_rk,4.0_rk]
        integer, parameter :: degree(3) = [2,1,1]
        type(nurbs_volume) :: va, vb
        type(nurbs_multipatch_volume) :: mp
        integer :: ida, idb

        call va%set(&
            knot1  = knot_normal,&
            knot2  = knot_tangent,&
            knot3  = knot_tangent,&
            Xc     = Xca,&
            Wc     = Wa,&
            degree = degree)
        call vb%set(&
            knot1  = knot_normal,&
            knot2  = knot_tangent,&
            knot3  = knot_tangent,&
            Xc     = Xcb,&
            Wc     = Wb,&
            degree = degree)
        call mp%add_patch(va, ida)
        call mp%add_patch(vb, idb)
        call mp%connect(ida, "u_max", idb, "u_min", 0)
        call va%err%print()
        call vb%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume unclamped C0 trace",&
            res      = mp%is_valid(tol),&
            expected = .true.,&
            msg      = "Basis-weighted volume face traces must match.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0047


    subroutine forcad_multipatch_volume_0048(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: knot_normal(7) = [0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk]
        real(rk), parameter :: knot_tangent(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xca(16,1) = reshape([&
            -4.0_rk,-2.0_rk,0.0_rk,2.0_rk,-4.0_rk,-2.0_rk,0.0_rk,2.0_rk,&
            -4.0_rk,-2.0_rk,0.0_rk,2.0_rk,-4.0_rk,-2.0_rk,0.0_rk,2.0_rk], [16,1])
        real(rk), parameter :: Xcb(16,1) = reshape([&
            -2.0_rk,4.0_rk,5.0_rk,6.0_rk,-2.0_rk,4.0_rk,5.0_rk,6.0_rk,&
            -2.0_rk,4.0_rk,5.0_rk,6.0_rk,-2.0_rk,4.0_rk,5.0_rk,6.0_rk], [16,1])
        integer, parameter :: degree(3) = [2,1,1]
        type(nurbs_volume) :: va, vb
        type(nurbs_multipatch_volume) :: mp
        integer, allocatable :: map(:)
        integer :: ida, idb

        call va%set(knot_normal, knot_tangent, knot_tangent, Xca, degree=degree)
        call vb%set(knot_normal, knot_tangent, knot_tangent, Xcb, degree=degree)
        call mp%add_patch(va, ida)
        call mp%add_patch(vb, idb)
        call mp%connect(ida, "u_max", idb, "u_min", 0)
        map = mp%cmp_dof_map()
        call va%err%print()
        call vb%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume unclamped C0 map",&
            res      = map,&
            expected = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,&
                17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32],&
            msg      = "An unclamped volume trace must remain a constraint.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0048


    subroutine forcad_multipatch_volume_0049(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_multipatch_volume) :: mp

        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume empty multipatch invalid",&
            res      = mp%is_valid(),&
            expected = .false.,&
            msg      = "An empty volume multipatch is not geometry.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0049


    subroutine forcad_multipatch_volume_0050(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: volume
        type(nurbs_multipatch_volume) :: mp
        integer :: id

        call mp%add_patch(volume, id)
        call volume%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume rejects uninitialized patch",&
            res      = id,&
            expected = 0,&
            msg      = "An uninitialized volume must not be added.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0050


    subroutine forcad_multipatch_volume_0051(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_volume) :: volume
        type(nurbs_multipatch_volume) :: mp
        integer :: id

        call volume%set(&
            Xth_dir1    = [0.0_rk,1.0_rk],&
            Xth_dir2    = [0.0_rk,1.0_rk],&
            Xth_dir3    = [0.0_rk,1.0_rk],&
            degree      = [1,1,1],&
            continuity1 = [-1,-1],&
            continuity2 = [-1,-1],&
            continuity3 = [-1,-1])
        call mp%add_patch(volume, id)
        call volume%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume rejects space-only patch",&
            res      = id,&
            expected = 0,&
            msg      = "A volume spline space is not complete geometry.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0051


    subroutine forcad_multipatch_volume_0052(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 1.0e-10_rk
        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        type(nurbs_volume) :: volume_a, volume_b
        type(nurbs_multipatch_volume) :: mp
        real(rk) :: Xc_a(8,3), Xc_b(8,3), transition_jet(2,2,3,1)
        real(rk) :: drift1, drift2, normal_rate, eta
        integer :: id_a, id_b, i, j, k, index

        do k = 0, 1
            do j = 0, 1
                drift1 = 0.1_rk + 0.1_rk*real(j,rk) + 0.05_rk*real(k,rk)
                drift2 = -0.1_rk + 0.05_rk*real(j,rk) + 0.1_rk*real(k,rk)
                normal_rate = 0.5_rk + 0.1_rk*real(j+k,rk)
                transition_jet(j+1,k+1,1,1) = drift1
                transition_jet(j+1,k+1,2,1) = drift2
                transition_jet(j+1,k+1,3,1) = normal_rate
                do i = 0, 1
                    index = i + 1 + 2*j + 4*k
                    eta = real(i-1,rk)
                    Xc_a(index,:) = [&
                        real(j,rk) + drift1*eta,&
                        real(k,rk) + drift2*eta,&
                        normal_rate*eta]
                    Xc_b(index,:) = [real(j,rk),real(k,rk),real(i,rk)]
                end do
            end do
        end do
        call volume_a%set(knot, knot, knot, Xc_a, degree=[1,1,1])
        call volume_b%set(knot, knot, knot, Xc_b, degree=[1,1,1])
        call mp%add_patch(volume_a, id_a)
        call mp%add_patch(volume_b, id_b)
        call mp%connect(&
            patch_a        = id_a,&
            side_a         = "u_max",&
            patch_b        = id_b,&
            side_b         = "u_min",&
            continuity     = 1,&
            geometric      = .true.,&
            transition_jet = transition_jet)
        call volume_a%err%print()
        call volume_b%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume general G1",&
            res      = mp%is_valid(tol),&
            expected = .true.,&
            msg      = "Tangentially varying volume G1 must be valid.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0052


    subroutine forcad_multipatch_volume_0053(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        type(nurbs_volume) :: volume_a, volume_b
        type(nurbs_multipatch_volume) :: mp
        real(rk) :: Xc_a(8,3), Xc_b(8,3), transition_jet(2,2,3,1)
        real(rk) :: drift1, drift2, normal_rate, eta
        real(rk), allocatable :: val(:)
        integer, allocatable :: rowptr(:), col(:)
        integer :: id_a, id_b, i, j, k, index

        do k = 0, 1
            do j = 0, 1
                drift1 = 0.1_rk + 0.1_rk*real(j,rk) + 0.05_rk*real(k,rk)
                drift2 = -0.1_rk + 0.05_rk*real(j,rk) + 0.1_rk*real(k,rk)
                normal_rate = 0.5_rk + 0.1_rk*real(j+k,rk)
                transition_jet(j+1,k+1,1,1) = drift1
                transition_jet(j+1,k+1,2,1) = drift2
                transition_jet(j+1,k+1,3,1) = normal_rate
                do i = 0, 1
                    index = i + 1 + 2*j + 4*k
                    eta = real(i-1,rk)
                    Xc_a(index,:) = [&
                        real(j,rk) + drift1*eta,&
                        real(k,rk) + drift2*eta,&
                        normal_rate*eta]
                    Xc_b(index,:) = [real(j,rk),real(k,rk),real(i,rk)]
                end do
            end do
        end do
        call volume_a%set(knot, knot, knot, Xc_a, degree=[1,1,1])
        call volume_b%set(knot, knot, knot, Xc_b, degree=[1,1,1])
        call mp%add_patch(volume_a, id_a)
        call mp%add_patch(volume_b, id_b)
        call mp%connect(&
            patch_a        = id_a,&
            side_a         = "u_max",&
            patch_b        = id_b,&
            side_b         = "u_min",&
            continuity     = 1,&
            geometric      = .true.,&
            transition_jet = transition_jet)
        call mp%cmp_dof_constraint(rowptr, col, val, geometric=.true.)
        call volume_a%err%print()
        call volume_b%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume general G1 rows",&
            res      = size(rowptr) - 1,&
            expected = 13,&
            msg      = "Volume general G1 collocation degree is incorrect.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0053


    subroutine forcad_multipatch_volume_0054(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 1.0e-10_rk
        real(rk), parameter :: knot_linear(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: knot_quadratic(6) = [0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk]
        type(nurbs_volume) :: volume_a, volume_b
        type(nurbs_multipatch_volume) :: mp
        real(rk), allocatable :: val(:)
        real(rk) :: Xc_a(12,3), Xc_b(12,3), Wc_a(12), Wc_b(12), projective_jet(2,2,1)
        real(rk) :: coordinate_a(3), coordinate_b(3), weight
        integer, allocatable :: map(:), rowptr(:), col(:)
        integer :: id_a, id_b, i, j, k, index

        do k = 0, 1
            do j = 0, 2
                weight = 1.0_rk + 0.5_rk*real(j,rk)
                do i = 0, 1
                    index = i + 1 + 2*j + 6*k
                    coordinate_a = [&
                        merge(0.0_rk,merge(1.0_rk/3.0_rk,1.0_rk,j == 1),j == 0),&
                        real(k,rk),&
                        real(i-1,rk)]
                    coordinate_b = [0.5_rk*real(j,rk),real(k,rk),real(i,rk)]
                    Xc_a(index,:) = coordinate_a
                    Xc_b(index,:) = coordinate_b
                    Wc_a(index) = weight
                    Wc_b(index) = 1.0_rk + 0.1_rk*real(i*(1+j+3*k),rk)
                end do
            end do
        end do
        projective_jet(1,:,1) = 1.0_rk
        projective_jet(2,:,1) = 2.0_rk
        call volume_a%set(&
            knot_linear,knot_quadratic,knot_linear,Xc_a,Wc_a,degree=[1,2,1])
        call volume_b%set(&
            knot_linear,knot_quadratic,knot_linear,Xc_b,Wc_b,degree=[1,2,1])
        call mp%add_patch(volume_a, id_a)
        call mp%add_patch(volume_b, id_b)
        call mp%connect(&
            patch_a        = id_a,&
            side_a         = "u_max",&
            patch_b        = id_b,&
            side_b         = "u_min",&
            continuity     = 0,&
            geometric      = .true.,&
            projective_jet = projective_jet)
        map = mp%cmp_dof_map()
        call mp%cmp_dof_constraint(rowptr, col, val, geometric=.true.)
        call volume_a%err%print()
        call volume_b%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume variable projective constraint",&
            res      = volume_a%is_rational() .and. volume_b%is_rational() .and. &
                mp%is_valid(tol) .and. maxval(map) == size(map) .and. &
                size(rowptr) > 1 .and. size(col) == size(val) .and. size(val) > 0,&
            expected = .true.,&
            msg      = "Variable projective face bases must remain sparse constraints.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0054


    subroutine forcad_multipatch_volume_0055(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 1.0e-10_rk
        real(rk), parameter :: knot_quadratic(6) = [0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: knot_linear(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        type(nurbs_volume) :: volume_a, volume_b
        type(nurbs_multipatch_volume) :: mp
        real(rk) :: Xc_a(18,3), Xc_b(18,3), derivative0(3,3), derivative1(3,3), derivative2(3,3)
        real(rk) :: transition_jet(2,1,3,2)
        integer :: id_a, id_b, i, j, k, index

        derivative0 = 0.0_rk
        derivative1 = 0.0_rk
        derivative2 = 0.0_rk
        derivative0(:,1) = [0.0_rk,0.5_rk,1.0_rk]
        derivative1(:,1) = [0.2_rk,0.25_rk,0.3_rk]
        derivative1(:,2) = [0.5_rk,0.55_rk,0.6_rk]
        derivative1(:,3) = [-0.1_rk,0.175_rk,0.55_rk]
        derivative2(:,1) = [-0.1_rk,-0.075_rk,-0.05_rk]
        derivative2(:,2) = [0.3_rk,0.25_rk,0.2_rk]
        derivative2(:,3) = [0.4_rk,0.67_rk,0.86_rk]
        transition_jet(:,1,1,1) = [0.2_rk,0.3_rk]
        transition_jet(:,1,2,1) = [-0.1_rk,-0.05_rk]
        transition_jet(:,1,3,1) = [0.5_rk,0.6_rk]
        transition_jet(:,1,1,2) = [-0.1_rk,-0.05_rk]
        transition_jet(:,1,2,2) = [0.2_rk,0.3_rk]
        transition_jet(:,1,3,2) = [0.3_rk,0.2_rk]
        do k = 0, 1
            do j = 0, 2
                derivative0(j+1,3) = real(k,rk)
                do i = 0, 2
                    index = i + 1 + 3*j + 9*k
                    Xc_b(index,:) = [&
                        0.5_rk*real(j,rk),&
                        0.5_rk*real(i,rk),&
                        real(k,rk) + 0.25_rk*real(i*j,rk)]
                    Xc_a(index,:) = derivative0(j+1,:) - &
                        (1.0_rk - 0.5_rk*real(i,rk))*derivative1(j+1,:)
                    if (i == 0) Xc_a(index,:) = Xc_a(index,:) + 0.5_rk*derivative2(j+1,:)
                end do
            end do
        end do
        call volume_a%set(&
            knot_quadratic,knot_quadratic,knot_linear,Xc_a,degree=[2,2,1])
        call volume_b%set(&
            knot_quadratic,knot_quadratic,knot_linear,Xc_b,degree=[2,2,1])
        call mp%add_patch(volume_a, id_a)
        call mp%add_patch(volume_b, id_b)
        call mp%connect(&
            patch_a        = id_a,&
            side_a         = "u_max",&
            patch_b        = id_b,&
            side_b         = "u_min",&
            continuity     = 2,&
            geometric      = .true.,&
            transition_jet = transition_jet)
        call volume_a%err%print()
        call volume_b%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "volume general G2",&
            res      = mp%is_valid(tol),&
            expected = .true.,&
            msg      = "Volume G2 must include varying mixed transition terms.",&
            group    = "forcad_multipatch_volume")
        ti = ti + 1

    end subroutine forcad_multipatch_volume_0055


    subroutine run_multipatch_volume_tests(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        call forcad_multipatch_volume_0001(ut, ti)
        call forcad_multipatch_volume_0002(ut, ti)
        call forcad_multipatch_volume_0003(ut, ti)
        call forcad_multipatch_volume_0004(ut, ti)
        call forcad_multipatch_volume_0005(ut, ti)
        call forcad_multipatch_volume_0006(ut, ti)
        call forcad_multipatch_volume_0007(ut, ti)
        call forcad_multipatch_volume_0008(ut, ti)
        call forcad_multipatch_volume_0009(ut, ti)
        call forcad_multipatch_volume_0010(ut, ti)
        call forcad_multipatch_volume_0011(ut, ti)
        call forcad_multipatch_volume_0012(ut, ti)
        call forcad_multipatch_volume_0013(ut, ti)
        call forcad_multipatch_volume_0014(ut, ti)
        call forcad_multipatch_volume_0015(ut, ti)
        call forcad_multipatch_volume_0016(ut, ti)
        call forcad_multipatch_volume_0017(ut, ti)
        call forcad_multipatch_volume_0018(ut, ti)
        call forcad_multipatch_volume_0019(ut, ti)
        call forcad_multipatch_volume_0020(ut, ti)
        call forcad_multipatch_volume_0021(ut, ti)
        call forcad_multipatch_volume_0022(ut, ti)
        call forcad_multipatch_volume_0023(ut, ti)
        call forcad_multipatch_volume_0024(ut, ti)
        call forcad_multipatch_volume_0025(ut, ti)
        call forcad_multipatch_volume_0026(ut, ti)
        call forcad_multipatch_volume_0027(ut, ti)
        call forcad_multipatch_volume_0028(ut, ti)
        call forcad_multipatch_volume_0029(ut, ti)
        call forcad_multipatch_volume_0030(ut, ti)
        call forcad_multipatch_volume_0031(ut, ti)
        call forcad_multipatch_volume_0032(ut, ti)
        call forcad_multipatch_volume_0033(ut, ti)
        call forcad_multipatch_volume_0034(ut, ti)
        call forcad_multipatch_volume_0035(ut, ti)
        call forcad_multipatch_volume_0036(ut, ti)
        call forcad_multipatch_volume_0037(ut, ti)
        call forcad_multipatch_volume_0038(ut, ti)
        call forcad_multipatch_volume_0039(ut, ti)
        call forcad_multipatch_volume_0040(ut, ti)
        call forcad_multipatch_volume_0041(ut, ti)
        call forcad_multipatch_volume_0042(ut, ti)
        call forcad_multipatch_volume_0043(ut, ti)
        call forcad_multipatch_volume_0044(ut, ti)
        call forcad_multipatch_volume_0045(ut, ti)
        call forcad_multipatch_volume_0046(ut, ti)
        call forcad_multipatch_volume_0047(ut, ti)
        call forcad_multipatch_volume_0048(ut, ti)
        call forcad_multipatch_volume_0049(ut, ti)
        call forcad_multipatch_volume_0050(ut, ti)
        call forcad_multipatch_volume_0051(ut, ti)
        call forcad_multipatch_volume_0052(ut, ti)
        call forcad_multipatch_volume_0053(ut, ti)
        call forcad_multipatch_volume_0054(ut, ti)
        call forcad_multipatch_volume_0055(ut, ti)

    end subroutine run_multipatch_volume_tests

end module forcad_test_multipatch_volume
