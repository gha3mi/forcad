module forcad_test_multipatch_surface

    use forcad_kinds, only: rk
    use forcad_multipatch, only: multipatch_connection
    use forcad_multipatch_surface, only: nurbs_multipatch_surface
    use forcad_nurbs_surface, only: nurbs_surface
    use forunittest, only: unit_tests

    implicit none

    private
    public :: run_multipatch_surface_tests

contains

    subroutine forcad_multipatch_surface_0001(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_surface_tol = 1.0e-10_rk
        type(nurbs_surface) :: s1, s2
        type(nurbs_multipatch_surface) :: mp
        integer :: id1
        integer :: id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])

        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        call s1%err%print()
        call s2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface valid C0",&
            res      = mp%is_valid(multipatch_surface_tol),&
            expected = .true.,&
            msg      = "surface C0 connection should be valid",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0001


    subroutine forcad_multipatch_surface_0002(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2
        type(nurbs_multipatch_surface) :: mp
        integer :: id1
        integer :: id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])

        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        call s1%err%print()
        call s2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface npatch",&
            res      = mp%get_npatch(),&
            expected = 2,&
            msg      = "surface multipatch patch count mismatch",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0002


    subroutine forcad_multipatch_surface_0003(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2
        type(nurbs_multipatch_surface) :: mp
        integer :: id1
        integer :: id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])

        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        call s1%err%print()
        call s2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface nconnection",&
            res      = mp%get_nconnection(),&
            expected = 1,&
            msg      = "surface multipatch connection count mismatch",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0003


    subroutine forcad_multipatch_surface_0004(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2
        type(nurbs_multipatch_surface) :: mp
        type(multipatch_connection) :: conn
        integer :: id1
        integer :: id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])

        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)

        conn = mp%get_connection(1)
        call s1%err%print()
        call s2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface connection continuity",&
            res      = conn%get_continuity(),&
            expected = 0,&
            msg      = "surface connection continuity metadata mismatch",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0004


    subroutine forcad_multipatch_surface_0005(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2
        type(nurbs_multipatch_surface) :: mp
        type(multipatch_connection) :: conn
        integer :: id1
        integer :: id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])

        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)

        conn = mp%get_connection(1)
        call s1%err%print()
        call s2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface connection reverse",&
            res      = conn%is_reversed(),&
            expected = .false.,&
            msg      = "surface connection reverse metadata mismatch",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0005


    subroutine forcad_multipatch_surface_0006(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2
        type(nurbs_multipatch_surface) :: mp
        type(multipatch_connection) :: conn
        integer, allocatable :: offsets(:)
        integer :: id1
        integer :: id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])

        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)

        conn = mp%get_connection(1)

        offsets = mp%cmp_dof_offsets()
        call s1%err%print()
        call s2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface offsets",&
            res      = offsets,&
            expected = [0, 4, 8],&
            msg      = "surface patch dof offsets mismatch",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0006


    subroutine forcad_multipatch_surface_0007(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2
        type(nurbs_multipatch_surface) :: mp
        type(multipatch_connection) :: conn
        integer, allocatable :: map(:)
        integer, allocatable :: offsets(:)
        integer :: id1
        integer :: id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])

        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)

        conn = mp%get_connection(1)

        offsets = mp%cmp_dof_offsets()

        map = mp%cmp_dof_map()
        call s1%err%print()
        call s2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface shared dofs",&
            res      = maxval(map),&
            expected = 6,&
            msg      = "surface interface dofs were not shared",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0007


    subroutine forcad_multipatch_surface_0008(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2
        type(nurbs_multipatch_surface) :: mp, mp_c1
        type(multipatch_connection) :: conn
        real(rk), allocatable :: val(:)
        integer, allocatable :: map(:), offsets(:), rowptr(:), col(:)
        integer :: id1
        integer :: id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])

        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)

        conn = mp%get_connection(1)

        offsets = mp%cmp_dof_offsets()

        map = mp%cmp_dof_map()

        call mp_c1%add_patch(s1, id1)
        call mp_c1%add_patch(s2, id2)
        call mp_c1%connect(id1, "u_max", id2, "u_min", 1)
        call mp_c1%cmp_dof_constraint(rowptr, col, val)
        call s1%err%print()
        call s2%err%print()
        call mp%err%print()
        call mp_c1%err%print()

        call ut%test(ti)%check(&
            name     = "surface C1 constraint count",&
            res      = size(rowptr)-1,&
            expected = 4,&
            msg      = "Surface C1 constraint count is incorrect.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0008


    subroutine forcad_multipatch_surface_0009(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_surface_tol = 1.0e-10_rk
        type(nurbs_surface) :: s1, s2
        type(nurbs_multipatch_surface) :: mp, mp_c1
        type(multipatch_connection) :: conn
        real(rk), allocatable :: Xc(:,:), val(:), dof(:)
        real(rk) :: residual, row_res
        integer, allocatable :: map(:), offsets(:), rowptr(:), col(:)
        integer :: id1, id2, row, j

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])

        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)

        conn = mp%get_connection(1)

        offsets = mp%cmp_dof_offsets()

        map = mp%cmp_dof_map()

        call mp_c1%add_patch(s1, id1)
        call mp_c1%add_patch(s2, id2)
        call mp_c1%connect(id1, "u_max", id2, "u_min", 1)
        call mp_c1%cmp_dof_constraint(rowptr, col, val)
        allocate(dof(8))
        Xc = s1%get_Xc()
        dof(1:4) = Xc(:,1)
        Xc = s2%get_Xc()
        dof(5:8) = Xc(:,1)
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do
        call s1%err%print()
        call s2%err%print()
        call mp%err%print()
        call mp_c1%err%print()

        call ut%test(ti)%check(&
            name     = "surface C1 constraint residual",&
            res      = residual,&
            expected = 0.0_rk,&
            tol      = multipatch_surface_tol,&
            msg      = "Surface C1 constraint residual is incorrect.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0009


    subroutine forcad_multipatch_surface_0010(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_surface_tol = 1.0e-10_rk
        type(nurbs_surface) :: sg2_a
        type(nurbs_surface) :: sg2_b
        type(nurbs_multipatch_surface) :: mp_g2
        real(rk) :: xg2_a(6,3)
        real(rk) :: xg2_b(6,3)
        real(rk) :: knot_open(4)
        real(rk) :: knot_quadratic(6)
        integer :: id1
        integer :: id2

        knot_quadratic = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot_open = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        xg2_a = 0.0_rk
        xg2_b = 0.0_rk
        xg2_a(:,1) = [0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk]
        xg2_a(:,2) = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_b(:,1) = [1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk]
        xg2_b(:,2) = xg2_a(:,2)
        call sg2_a%set(knot_quadratic, knot_open, xg2_a, degree=[2,1])
        call sg2_b%set(knot_quadratic, knot_open, xg2_b, degree=[2,1])
        call mp_g2%add_patch(sg2_a, id1)
        call mp_g2%add_patch(sg2_b, id2)
        call mp_g2%connect(&
            patch_a            = id1,&
            side_a             = "u_max",&
            patch_b            = id2,&
            side_b             = "u_min",&
            continuity         = 2,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk])
        call sg2_a%err%print()
        call sg2_b%err%print()
        call mp_g2%err%print()

        call ut%test(ti)%check(&
            name     = "surface valid G2",&
            res      = mp_g2%is_valid(multipatch_surface_tol),&
            expected = .true.,&
            msg      = "nonlinearly reparameterized planar patches should be G2",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0010


    subroutine forcad_multipatch_surface_0011(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sg2_a
        type(nurbs_surface) :: sg2_b
        type(nurbs_multipatch_surface) :: mp_g2
        type(multipatch_connection) :: conn
        real(rk) :: xg2_a(6,3)
        real(rk) :: xg2_b(6,3)
        real(rk) :: knot_open(4)
        real(rk) :: knot_quadratic(6)
        integer :: id1
        integer :: id2

        knot_quadratic = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot_open = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        xg2_a = 0.0_rk
        xg2_b = 0.0_rk
        xg2_a(:,1) = [0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk]
        xg2_a(:,2) = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_b(:,1) = [1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk]
        xg2_b(:,2) = xg2_a(:,2)
        call sg2_a%set(knot_quadratic, knot_open, xg2_a, degree=[2,1])
        call sg2_b%set(knot_quadratic, knot_open, xg2_b, degree=[2,1])
        call mp_g2%add_patch(sg2_a, id1)
        call mp_g2%add_patch(sg2_b, id2)
        call mp_g2%connect(&
            patch_a            = id1,&
            side_a             = "u_max",&
            patch_b            = id2,&
            side_b             = "u_min",&
            continuity         = 2,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk])
        conn = mp_g2%get_connection(1)
        call sg2_a%err%print()
        call sg2_b%err%print()
        call mp_g2%err%print()

        call ut%test(ti)%check(&
            name     = "surface G2 continuity flag",&
            res      = conn%is_geometric(),&
            expected = .true.,&
            msg      = "surface geometric continuity was not preserved",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0011


    subroutine forcad_multipatch_surface_0012(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sg2_a
        type(nurbs_surface) :: sg2_b
        type(nurbs_multipatch_surface) :: mp_g2
        type(multipatch_connection) :: conn
        real(rk) :: xg2_a(6,3)
        real(rk) :: xg2_b(6,3)
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
        xg2_a(:,1) = [0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk]
        xg2_a(:,2) = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_b(:,1) = [1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk]
        xg2_b(:,2) = xg2_a(:,2)
        call sg2_a%set(knot_quadratic, knot_open, xg2_a, degree=[2,1])
        call sg2_b%set(knot_quadratic, knot_open, xg2_b, degree=[2,1])
        call mp_g2%add_patch(sg2_a, id1)
        call mp_g2%add_patch(sg2_b, id2)
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
        call sg2_a%err%print()
        call sg2_b%err%print()
        call mp_g2%err%print()

        call ut%test(ti)%check(&
            name     = "surface G2 constraint count",&
            res      = size(rowptr)-1,&
            expected = 6,&
            msg      = "surface G2 should create three rows per tangential control point",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0012


    subroutine forcad_multipatch_surface_0013(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_surface_tol = 1.0e-10_rk
        type(nurbs_surface) :: sg2_a
        type(nurbs_surface) :: sg2_b
        type(nurbs_multipatch_surface) :: mp_g2
        type(multipatch_connection) :: conn
        real(rk) :: xg2_a(6,3)
        real(rk) :: xg2_b(6,3)
        real(rk) :: knot_open(4), knot_quadratic(6), residual, row_res
        real(rk), allocatable :: val(:), dof(:)
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1, id2, row, j

        knot_quadratic = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot_open = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        xg2_a = 0.0_rk
        xg2_b = 0.0_rk
        xg2_a(:,1) = [0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk]
        xg2_a(:,2) = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_b(:,1) = [1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk]
        xg2_b(:,2) = xg2_a(:,2)
        call sg2_a%set(knot_quadratic, knot_open, xg2_a, degree=[2,1])
        call sg2_b%set(knot_quadratic, knot_open, xg2_b, degree=[2,1])
        call mp_g2%add_patch(sg2_a, id1)
        call mp_g2%add_patch(sg2_b, id2)
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
        allocate(dof(12))
        dof(1:6) = xg2_a(:,1)
        dof(7:12) = xg2_b(:,1)
        residual = 0.0_rk
        do row = 1, size(rowptr)-1
            row_res = 0.0_rk
            do j = rowptr(row), rowptr(row+1)-1
                row_res = row_res + val(j)*dof(col(j))
            end do
            residual = max(residual, abs(row_res))
        end do
        call sg2_a%err%print()
        call sg2_b%err%print()
        call mp_g2%err%print()

        call ut%test(ti)%check(&
            name     = "surface G2 constraint residual",&
            res      = residual,&
            expected = 0.0_rk,&
            tol      = multipatch_surface_tol,&
            msg      = "surface G2 chain-rule constraints should vanish",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0013


    subroutine forcad_multipatch_surface_0014(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_surface_tol = 1.0e-10_rk
        type(nurbs_surface) :: sg2_a, sg2_b, srg2_a, srg2_b
        type(nurbs_multipatch_surface) :: mp_g2
        type(nurbs_multipatch_surface) :: mp_rg2
        type(multipatch_connection) :: conn
        real(rk) :: wrg2_a(6), wrg2_b(6), xg2_a(6,3), xg2_b(6,3)
        real(rk) :: knot_open(4), knot_quadratic(6), residual, row_res
        real(rk), allocatable :: val(:), dof(:)
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1, id2, row, j

        knot_quadratic = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot_open = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        xg2_a = 0.0_rk
        xg2_b = 0.0_rk
        xg2_a(:,1) = [0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk]
        xg2_a(:,2) = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_b(:,1) = [1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk]
        xg2_b(:,2) = xg2_a(:,2)
        call sg2_a%set(knot_quadratic, knot_open, xg2_a, degree=[2,1])
        call sg2_b%set(knot_quadratic, knot_open, xg2_b, degree=[2,1])
        call mp_g2%add_patch(sg2_a, id1)
        call mp_g2%add_patch(sg2_b, id2)
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
        allocate(dof(12))
        dof(1:6) = xg2_a(:,1)
        dof(7:12) = xg2_b(:,1)
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
        call srg2_a%set(knot_quadratic, knot_open, xg2_a, wrg2_a, degree=[2,1])
        call srg2_b%set(knot_quadratic, knot_open, xg2_b, wrg2_b, degree=[2,1])
        call mp_rg2%add_patch(srg2_a, id1)
        call mp_rg2%add_patch(srg2_b, id2)
        call mp_rg2%connect(&
            patch_a            = id1,&
            side_a             = "u_max",&
            patch_b            = id2,&
            side_b             = "u_min",&
            continuity         = 2,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk])
        call sg2_a%err%print()
        call sg2_b%err%print()
        call srg2_a%err%print()
        call srg2_b%err%print()
        call mp_g2%err%print()
        call mp_rg2%err%print()

        call ut%test(ti)%check(&
            name     = "surface rational G2",&
            res      = mp_rg2%is_valid(multipatch_surface_tol),&
            expected = .true.,&
            msg      = "Surface rational G2 is incorrect.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0014


    subroutine forcad_multipatch_surface_0015(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: sg2_a, sg2_b, srg2_a, srg2_b
        type(nurbs_multipatch_surface) :: mp_g2
        type(nurbs_multipatch_surface) :: mp_rg2
        type(multipatch_connection) :: conn
        real(rk) :: wrg2_a(6), wrg2_b(6), xg2_a(6,3), xg2_b(6,3)
        real(rk) :: knot_open(4), knot_quadratic(6), residual, row_res
        real(rk), allocatable :: val(:), dof(:)
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1, id2, row, j

        knot_quadratic = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot_open = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        xg2_a = 0.0_rk
        xg2_b = 0.0_rk
        xg2_a(:,1) = [0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk]
        xg2_a(:,2) = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_b(:,1) = [1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk]
        xg2_b(:,2) = xg2_a(:,2)
        call sg2_a%set(knot_quadratic, knot_open, xg2_a, degree=[2,1])
        call sg2_b%set(knot_quadratic, knot_open, xg2_b, degree=[2,1])
        call mp_g2%add_patch(sg2_a, id1)
        call mp_g2%add_patch(sg2_b, id2)
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
        allocate(dof(12))
        dof(1:6) = xg2_a(:,1)
        dof(7:12) = xg2_b(:,1)
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
        call srg2_a%set(knot_quadratic, knot_open, xg2_a, wrg2_a, degree=[2,1])
        call srg2_b%set(knot_quadratic, knot_open, xg2_b, wrg2_b, degree=[2,1])
        call mp_rg2%add_patch(srg2_a, id1)
        call mp_rg2%add_patch(srg2_b, id2)
        call mp_rg2%connect(&
            patch_a            = id1,&
            side_a             = "u_max",&
            patch_b            = id2,&
            side_b             = "u_min",&
            continuity         = 2,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk])
        call mp_g2%cmp_dof_constraint(rowptr, col, val, geometric=.false.)
        call sg2_a%err%print()
        call sg2_b%err%print()
        call srg2_a%err%print()
        call srg2_b%err%print()
        call mp_g2%err%print()
        call mp_rg2%err%print()

        call ut%test(ti)%check(&
            name     = "surface separate C constraints",&
            res      = size(rowptr)-1,&
            expected = 0,&
            msg      = "filtering a G2 connection for C constraints must return no rows",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0015


    subroutine forcad_multipatch_surface_0016(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_surface_tol = 1.0e-10_rk
        type(nurbs_surface) :: sg2_a, sg2_b, srg2_a, srg2_b
        type(nurbs_multipatch_surface) :: mp_g2, mp_rg2, mp_not_c2
        type(multipatch_connection) :: conn
        real(rk) :: wrg2_a(6), wrg2_b(6), xg2_a(6,3), xg2_b(6,3)
        real(rk) :: knot_open(4), knot_quadratic(6), residual, row_res
        real(rk), allocatable :: val(:), dof(:)
        integer, allocatable :: rowptr(:), col(:)
        integer :: id1, id2, row, j

        knot_quadratic = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        knot_open = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        xg2_a = 0.0_rk
        xg2_b = 0.0_rk
        xg2_a(:,1) = [0.0_rk, 0.5_rk, 1.0_rk, 0.0_rk, 0.5_rk, 1.0_rk]
        xg2_a(:,2) = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk]
        xg2_b(:,1) = [1.0_rk, 2.0_rk, 4.0_rk, 1.0_rk, 2.0_rk, 4.0_rk]
        xg2_b(:,2) = xg2_a(:,2)
        call sg2_a%set(knot_quadratic, knot_open, xg2_a, degree=[2,1])
        call sg2_b%set(knot_quadratic, knot_open, xg2_b, degree=[2,1])
        call mp_g2%add_patch(sg2_a, id1)
        call mp_g2%add_patch(sg2_b, id2)
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
        allocate(dof(12))
        dof(1:6) = xg2_a(:,1)
        dof(7:12) = xg2_b(:,1)
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
        call srg2_a%set(knot_quadratic, knot_open, xg2_a, wrg2_a, degree=[2,1])
        call srg2_b%set(knot_quadratic, knot_open, xg2_b, wrg2_b, degree=[2,1])
        call mp_rg2%add_patch(srg2_a, id1)
        call mp_rg2%add_patch(srg2_b, id2)
        call mp_rg2%connect(&
            patch_a            = id1,&
            side_a             = "u_max",&
            patch_b            = id2,&
            side_b             = "u_min",&
            continuity         = 2,&
            geometric          = .true.,&
            reparameterization = [0.5_rk, -0.25_rk])
        call mp_g2%cmp_dof_constraint(rowptr, col, val, geometric=.false.)

        call mp_not_c2%add_patch(sg2_a, id1)
        call mp_not_c2%add_patch(sg2_b, id2)
        call mp_not_c2%connect(id1, "u_max", id2, "u_min", 2)
        call sg2_a%err%print()
        call sg2_b%err%print()
        call srg2_a%err%print()
        call srg2_b%err%print()
        call mp_g2%err%print()
        call mp_rg2%err%print()
        call mp_not_c2%err%print()

        call ut%test(ti)%check(&
            name     = "surface G2 is not C2",&
            res      = mp_not_c2%is_valid(multipatch_surface_tol),&
            expected = .false.,&
            msg      = "The G2 surface interface must not satisfy C2.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0016


    subroutine forcad_multipatch_surface_0017(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1
        type(nurbs_surface) :: s2
        type(nurbs_multipatch_surface) :: mp
        integer, allocatable :: elem(:,:)
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        elem = mp%cmp_elem(shared=.true.)
        call s1%err%print()
        call s2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface elem shape",&
            res      = shape(elem),&
            expected = [2, 4],&
            msg      = "surface multipatch element shape mismatch",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0017


    subroutine forcad_multipatch_surface_0018(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2
        type(nurbs_multipatch_surface) :: mp
        integer, allocatable :: elem(:,:)
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        elem = mp%cmp_elem(shared=.true.)
        elem = mp%cmp_elem(shared=.false.)
        call s1%err%print()
        call s2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface unshared elem max",&
            res      = maxval(elem),&
            expected = 8,&
            msg      = "Surface unshared elem max is incorrect.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0018


    subroutine forcad_multipatch_surface_0019(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2
        type(nurbs_multipatch_surface) :: mp
        integer, allocatable :: elem(:,:), patch_id(:), local_id(:)
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        elem = mp%cmp_elem(shared=.true.)
        elem = mp%cmp_elem(shared=.false.)

        patch_id = mp%cmp_elem_patch()
        local_id = mp%cmp_elem_local()
        call s1%err%print()
        call s2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface elem patch ids",&
            res      = patch_id,&
            expected = [1, 2],&
            msg      = "surface element patch ids mismatch",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0019


    subroutine forcad_multipatch_surface_0020(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1
        type(nurbs_surface) :: s2
        type(nurbs_multipatch_surface) :: mp
        integer, allocatable :: elem(:,:), patch_id(:), local_id(:)
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        elem = mp%cmp_elem(shared=.true.)
        elem = mp%cmp_elem(shared=.false.)

        patch_id = mp%cmp_elem_patch()
        local_id = mp%cmp_elem_local()
        call s1%err%print()
        call s2%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface elem local ids",&
            res      = local_id,&
            expected = [1, 1],&
            msg      = "surface local element ids mismatch",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0020


    subroutine forcad_multipatch_surface_0021(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_surface_tol = 1.0e-10_rk
        type(nurbs_surface) :: s1, s2, patch
        type(nurbs_multipatch_surface) :: mp
        real(rk), allocatable :: Xc(:,:)
        integer, allocatable :: elem(:,:), patch_id(:), local_id(:)
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        elem = mp%cmp_elem(shared=.true.)
        elem = mp%cmp_elem(shared=.false.)

        patch_id = mp%cmp_elem_patch()
        local_id = mp%cmp_elem_local()

        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        call s1%err%print()
        call s2%err%print()
        call patch%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface get patch",&
            res      = Xc(1,:),&
            expected = [1.0_rk, 0.0_rk, 0.0_rk],&
            tol      = multipatch_surface_tol,&
            msg      = "surface get_patch returned wrong patch",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0021


    subroutine forcad_multipatch_surface_0022(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2, patch
        type(nurbs_multipatch_surface) :: mp
        real(rk), allocatable :: Xc(:,:)
        integer, allocatable :: elem(:,:), patch_id(:), local_id(:)
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        elem = mp%cmp_elem(shared=.true.)
        elem = mp%cmp_elem(shared=.false.)

        patch_id = mp%cmp_elem_patch()
        local_id = mp%cmp_elem_local()

        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        call s1%err%print()
        call s2%err%print()
        call patch%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface rational flag",&
            res      = mp%is_rational(),&
            expected = .false.,&
            msg      = "surface multipatch should report non-rational patches",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0022


    subroutine forcad_multipatch_surface_0023(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_surface_tol = 1.0e-10_rk
        type(nurbs_surface) :: s1
        type(nurbs_surface) :: s2
        type(nurbs_surface) :: sr1
        type(nurbs_surface) :: sr2
        type(nurbs_multipatch_surface) :: mp_rational
        real(rk) :: wrs1(4)
        real(rk) :: wrs2(4)
        real(rk) :: knot_open(4)
        integer, allocatable :: map(:)
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        knot_open = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        wrs1 = [1.0_rk, 2.0_rk, 4.0_rk, 3.0_rk]
        wrs2 = [10.0_rk, 5.0_rk, 15.0_rk, 7.0_rk]
        call sr1%set_tetragon([1.0_rk, 1.0_rk], [2, 2], wrs1)
        call sr2%set_tetragon([1.0_rk, 1.0_rk], [2, 2], wrs2)
        call sr2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_rational%add_patch(sr1, id1)
        call mp_rational%add_patch(sr2, id2)
        call mp_rational%connect(id1, "u_max", id2, "u_min", 0)
        map = mp_rational%cmp_dof_map()
        call s1%err%print()
        call s2%err%print()
        call sr1%err%print()
        call sr2%err%print()
        call mp_rational%err%print()

        call ut%test(ti)%check(&
            name     = "surface constant rational basis sharing",&
            res      = mp_rational%is_valid(multipatch_surface_tol) .and. maxval(map) == 6,&
            expected = .true.,&
            msg      = "Constant projective scaling must preserve shared trace DOFs.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0023


    subroutine forcad_multipatch_surface_0024(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_surface_tol = 1.0e-10_rk
        type(nurbs_surface) :: s1
        type(nurbs_surface) :: s2
        type(nurbs_surface) :: sr1
        type(nurbs_surface) :: sr2
        type(nurbs_surface) :: s_scaled_trace
        type(nurbs_multipatch_surface) :: mp_rational
        type(nurbs_multipatch_surface) :: mp_scaled_trace
        real(rk), allocatable :: Xc(:,:)
        real(rk) :: wrs1(4), wrs2(4), knot_open(4), knot_trace_scaled(4), knot_trace_bad(4)
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        knot_open = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        wrs1 = [1.0_rk, 2.0_rk, 4.0_rk, 3.0_rk]
        wrs2 = [10.0_rk, 5.0_rk, 15.0_rk, 7.0_rk]
        call sr1%set_tetragon([1.0_rk, 1.0_rk], [2, 2], wrs1)
        call sr2%set_tetragon([1.0_rk, 1.0_rk], [2, 2], wrs2)
        call sr2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_rational%add_patch(sr1, id1)
        call mp_rational%add_patch(sr2, id2)
        call mp_rational%connect(id1, "u_max", id2, "u_min", 0)

        knot_trace_scaled = [0.0_rk, 0.0_rk, 2.0_rk, 2.0_rk]
        knot_trace_bad = [0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk]
        Xc = s2%get_Xc()
        call s_scaled_trace%set(knot_open, knot_trace_scaled, Xc, degree=[1, 1])
        call mp_scaled_trace%add_patch(s1, id1)
        call mp_scaled_trace%add_patch(s_scaled_trace, id2)
        call mp_scaled_trace%connect(id1, "u_max", id2, "u_min", 0)
        call s1%err%print()
        call s2%err%print()
        call sr1%err%print()
        call sr2%err%print()
        call s_scaled_trace%err%print()
        call mp_rational%err%print()
        call mp_scaled_trace%err%print()

        call ut%test(ti)%check(&
            name     = "surface affine-scaled trace knots",&
            res      = mp_scaled_trace%is_valid(multipatch_surface_tol),&
            expected = .true.,&
            msg      = "surface C0 trace must accept affine-equivalent tangential knots",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0024


    subroutine forcad_multipatch_surface_0025(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2, sr1, sr2, s_scaled_trace, s_bad_trace
        type(nurbs_multipatch_surface) :: mp_rational, mp_scaled_trace, mp_bad_trace
        real(rk), allocatable :: Xc(:,:)
        real(rk) :: wrs1(4), wrs2(4), knot_open(4), knot_trace_scaled(4), knot_trace_bad(4)
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        knot_open = [0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]
        wrs1 = [1.0_rk, 2.0_rk, 4.0_rk, 3.0_rk]
        wrs2 = [10.0_rk, 5.0_rk, 15.0_rk, 7.0_rk]
        call sr1%set_tetragon([1.0_rk, 1.0_rk], [2, 2], wrs1)
        call sr2%set_tetragon([1.0_rk, 1.0_rk], [2, 2], wrs2)
        call sr2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_rational%add_patch(sr1, id1)
        call mp_rational%add_patch(sr2, id2)
        call mp_rational%connect(id1, "u_max", id2, "u_min", 0)

        knot_trace_scaled = [0.0_rk, 0.0_rk, 2.0_rk, 2.0_rk]
        knot_trace_bad = [0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk]
        Xc = s2%get_Xc()
        call s_scaled_trace%set(knot_open, knot_trace_scaled, Xc, degree=[1, 1])
        call mp_scaled_trace%add_patch(s1, id1)
        call mp_scaled_trace%add_patch(s_scaled_trace, id2)
        call mp_scaled_trace%connect(id1, "u_max", id2, "u_min", 0)

        call s_bad_trace%set(knot_open, knot_trace_bad, Xc, degree=[1, 1])
        call mp_bad_trace%add_patch(s1, id1)
        call mp_bad_trace%add_patch(s_bad_trace, id2)
        call mp_bad_trace%connect(id1, "u_max", id2, "u_min", 0)
        call s1%err%print()
        call s2%err%print()
        call sr1%err%print()
        call sr2%err%print()
        call s_scaled_trace%err%print()
        call s_bad_trace%err%print()
        call mp_rational%err%print()
        call mp_scaled_trace%err%print()
        call mp_bad_trace%err%print()

        call ut%test(ti)%check(&
            name     = "surface incompatible trace knots err",&
            res      = mp_bad_trace%err%ok,&
            expected = .false.,&
            msg      = "Surface incompatible trace knots err is incorrect.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0025


    subroutine forcad_multipatch_surface_0026(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_surface_tol = 1.0e-10_rk
        type(nurbs_surface) :: s1, s2, patch
        type(nurbs_multipatch_surface) :: mp
        real(rk), allocatable :: Xc(:,:), Xg(:,:)
        integer :: id1, id2, file_unit

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        call mp%create(res1=3, res2=3)
        call mp%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%translate_Xg([0.0_rk, 0.0_rk, 1.0_rk])
        call mp%rotate_Xc(0.0_rk, 0.0_rk, 0.0_rk)
        call mp%rotate_Xg(0.0_rk, 0.0_rk, 0.0_rk)
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_surface_Xc_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_surface_Xg_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_surface_Xth_1.vtk", status="replace")
        close(file_unit, status="delete")
        call mp%export_Xc("vtk/forcad_test_multipatch_surface")
        call mp%export_Xg("vtk/forcad_test_multipatch_surface")
        call mp%export_Xth_in_Xg("vtk/forcad_test_multipatch_surface", res=3)
        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        Xg = patch%get_Xg()
        call s1%err%print()
        call s2%err%print()
        call patch%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface translate Xc",&
            res      = Xc(1,:),&
            expected = [2.0_rk, 0.0_rk, 0.0_rk],&
            tol      = multipatch_surface_tol,&
            msg      = "surface multipatch translate_Xc failed",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0026


    subroutine forcad_multipatch_surface_0027(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_surface_tol = 1.0e-10_rk
        type(nurbs_surface) :: s1, s2, patch
        type(nurbs_multipatch_surface) :: mp
        real(rk), allocatable :: Xc(:,:), Xg(:,:)
        integer :: id1, id2, file_unit

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        call mp%create(res1=3, res2=3)
        call mp%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%translate_Xg([0.0_rk, 0.0_rk, 1.0_rk])
        call mp%rotate_Xc(0.0_rk, 0.0_rk, 0.0_rk)
        call mp%rotate_Xg(0.0_rk, 0.0_rk, 0.0_rk)
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_surface_Xc_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_surface_Xg_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_surface_Xth_1.vtk", status="replace")
        close(file_unit, status="delete")
        call mp%export_Xc("vtk/forcad_test_multipatch_surface")
        call mp%export_Xg("vtk/forcad_test_multipatch_surface")
        call mp%export_Xth_in_Xg("vtk/forcad_test_multipatch_surface", res=3)
        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        Xg = patch%get_Xg()
        call s1%err%print()
        call s2%err%print()
        call patch%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface translate Xg",&
            res      = Xg(1,3),&
            expected = 1.0_rk,&
            tol      = multipatch_surface_tol,&
            msg      = "surface multipatch translate_Xg failed",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0027


    subroutine forcad_multipatch_surface_0028(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2, patch
        type(nurbs_multipatch_surface) :: mp
        real(rk), allocatable :: Xc(:,:), Xg(:,:)
        integer :: id1, id2, file_unit
        logical :: file_exists

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        call mp%create(res1=3, res2=3)
        call mp%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%translate_Xg([0.0_rk, 0.0_rk, 1.0_rk])
        call mp%rotate_Xc(0.0_rk, 0.0_rk, 0.0_rk)
        call mp%rotate_Xg(0.0_rk, 0.0_rk, 0.0_rk)
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_surface_Xc_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_surface_Xg_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_surface_Xth_1.vtk", status="replace")
        close(file_unit, status="delete")
        call mp%export_Xc("vtk/forcad_test_multipatch_surface")
        call mp%export_Xg("vtk/forcad_test_multipatch_surface")
        call mp%export_Xth_in_Xg("vtk/forcad_test_multipatch_surface", res=3)
        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        Xg = patch%get_Xg()
        inquire(file="vtk/forcad_test_multipatch_surface_Xc_1.vtk", exist=file_exists)
        call s1%err%print()
        call s2%err%print()
        call patch%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface export Xc",&
            res      = file_exists,&
            expected = .true.,&
            msg      = "surface multipatch export_Xc did not create the first VTK file",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0028


    subroutine forcad_multipatch_surface_0029(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2, patch
        type(nurbs_multipatch_surface) :: mp
        real(rk), allocatable :: Xc(:,:), Xg(:,:)
        integer :: id1, id2, file_unit
        logical :: file_exists

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        call mp%create(res1=3, res2=3)
        call mp%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%translate_Xg([0.0_rk, 0.0_rk, 1.0_rk])
        call mp%rotate_Xc(0.0_rk, 0.0_rk, 0.0_rk)
        call mp%rotate_Xg(0.0_rk, 0.0_rk, 0.0_rk)
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_surface_Xc_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_surface_Xg_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_surface_Xth_1.vtk", status="replace")
        close(file_unit, status="delete")
        call mp%export_Xc("vtk/forcad_test_multipatch_surface")
        call mp%export_Xg("vtk/forcad_test_multipatch_surface")
        call mp%export_Xth_in_Xg("vtk/forcad_test_multipatch_surface", res=3)
        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        Xg = patch%get_Xg()
        inquire(file="vtk/forcad_test_multipatch_surface_Xc_1.vtk", exist=file_exists)
        inquire(file="vtk/forcad_test_multipatch_surface_Xg_1.vtk", exist=file_exists)
        call s1%err%print()
        call s2%err%print()
        call patch%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface export Xg",&
            res      = file_exists,&
            expected = .true.,&
            msg      = "surface multipatch export_Xg did not create the first VTK file",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0029


    subroutine forcad_multipatch_surface_0030(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2, patch
        type(nurbs_multipatch_surface) :: mp
        real(rk), allocatable :: Xc(:,:), Xg(:,:)
        integer :: id1, id2, file_unit
        logical :: file_exists

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        call mp%connect(id1, "u_max", id2, "u_min", 0)
        call mp%create(res1=3, res2=3)
        call mp%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp%translate_Xg([0.0_rk, 0.0_rk, 1.0_rk])
        call mp%rotate_Xc(0.0_rk, 0.0_rk, 0.0_rk)
        call mp%rotate_Xg(0.0_rk, 0.0_rk, 0.0_rk)
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_surface_Xc_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_surface_Xg_1.vtk", status="replace")
        close(file_unit, status="delete")
        open(newunit=file_unit, file="vtk/forcad_test_multipatch_surface_Xth_1.vtk", status="replace")
        close(file_unit, status="delete")
        call mp%export_Xc("vtk/forcad_test_multipatch_surface")
        call mp%export_Xg("vtk/forcad_test_multipatch_surface")
        call mp%export_Xth_in_Xg("vtk/forcad_test_multipatch_surface", res=3)
        patch = mp%get_patch(2)
        Xc = patch%get_Xc()
        Xg = patch%get_Xg()
        inquire(file="vtk/forcad_test_multipatch_surface_Xc_1.vtk", exist=file_exists)
        inquire(file="vtk/forcad_test_multipatch_surface_Xg_1.vtk", exist=file_exists)
        inquire(file="vtk/forcad_test_multipatch_surface_Xth_1.vtk", exist=file_exists)
        call s1%err%print()
        call s2%err%print()
        call patch%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface export Xth in Xg",&
            res      = file_exists,&
            expected = .true.,&
            msg      = "Surface export Xth in Xg is incorrect.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0030


    subroutine forcad_multipatch_surface_0031(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1
        type(nurbs_surface) :: s2
        type(nurbs_multipatch_surface) :: mp_discontinuous
        integer, allocatable :: map(:)
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(s1, id1)
        call mp_discontinuous%add_patch(s2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()
        call s1%err%print()
        call s2%err%print()
        call mp_discontinuous%err%print()

        call ut%test(ti)%check(&
            name     = "surface discontinuous map",&
            res      = maxval(map),&
            expected = 8,&
            msg      = "discontinuous surface interface should not share dofs",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0031


    subroutine forcad_multipatch_surface_0032(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1
        type(nurbs_surface) :: s2
        type(nurbs_surface) :: s_mismatch
        type(nurbs_multipatch_surface) :: mp_discontinuous, mp_discontinuous_mismatch
        integer, allocatable :: map(:)
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(s1, id1)
        call mp_discontinuous%add_patch(s2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()

        call s_mismatch%set_tetragon([1.0_rk, 1.0_rk], [2, 3])
        call mp_discontinuous_mismatch%add_patch(s1, id1)
        call mp_discontinuous_mismatch%add_patch(s_mismatch, id2)
        call mp_discontinuous_mismatch%connect(id1, "u_max", id2, "u_min", -1)
        call s1%err%print()
        call s2%err%print()
        call s_mismatch%err%print()
        call mp_discontinuous%err%print()
        call mp_discontinuous_mismatch%err%print()

        call ut%test(ti)%check(&
            name     = "surface discontinuous mismatch err",&
            res      = mp_discontinuous_mismatch%err%ok,&
            expected = .true.,&
            msg      = "Surface discontinuous mismatch err is incorrect.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0032


    subroutine forcad_multipatch_surface_0033(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1
        type(nurbs_surface) :: s2
        type(nurbs_surface) :: s_mismatch
        type(nurbs_multipatch_surface) :: mp_discontinuous, mp_discontinuous_mismatch
        integer, allocatable :: map(:)
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(s1, id1)
        call mp_discontinuous%add_patch(s2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()

        call s_mismatch%set_tetragon([1.0_rk, 1.0_rk], [2, 3])
        call mp_discontinuous_mismatch%add_patch(s1, id1)
        call mp_discontinuous_mismatch%add_patch(s_mismatch, id2)
        call mp_discontinuous_mismatch%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous_mismatch%cmp_dof_map()
        call s1%err%print()
        call s2%err%print()
        call s_mismatch%err%print()
        call mp_discontinuous%err%print()
        call mp_discontinuous_mismatch%err%print()

        call ut%test(ti)%check(&
            name     = "surface discontinuous mismatch map",&
            res      = maxval(map),&
            expected = 10,&
            msg      = "discontinuous mismatched surface interface should not share dofs",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0033


    subroutine forcad_multipatch_surface_0034(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_surface_tol = 1.0e-10_rk
        type(nurbs_surface) :: s1, s2, s_bad, s_mismatch
        type(nurbs_multipatch_surface) :: mp_discontinuous, mp_discontinuous_mismatch
        type(nurbs_multipatch_surface) :: mp_bad_geom
        integer, allocatable :: map(:)
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(s1, id1)
        call mp_discontinuous%add_patch(s2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()

        call s_mismatch%set_tetragon([1.0_rk, 1.0_rk], [2, 3])
        call mp_discontinuous_mismatch%add_patch(s1, id1)
        call mp_discontinuous_mismatch%add_patch(s_mismatch, id2)
        call mp_discontinuous_mismatch%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous_mismatch%cmp_dof_map()

        call s_bad%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s_bad%translate_Xc([2.0_rk, 0.0_rk, 0.0_rk])
        call mp_bad_geom%add_patch(s1, id1)
        call mp_bad_geom%add_patch(s_bad, id2)
        call mp_bad_geom%connect(id1, "u_max", id2, "u_min", 0)
        call s1%err%print()
        call s2%err%print()
        call s_bad%err%print()
        call s_mismatch%err%print()
        call mp_discontinuous%err%print()
        call mp_discontinuous_mismatch%err%print()
        call mp_bad_geom%err%print()

        call ut%test(ti)%check(&
            name     = "surface invalid geometry",&
            res      = mp_bad_geom%is_valid(multipatch_surface_tol),&
            expected = .false.,&
            msg      = "surface multipatch should reject mismatched interface geometry",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0034


    subroutine forcad_multipatch_surface_0035(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2, s_bad, s_mismatch
        type(nurbs_multipatch_surface) :: mp_discontinuous, mp_discontinuous_mismatch
        type(nurbs_multipatch_surface) :: mp_bad_geom
        type(nurbs_multipatch_surface) :: mp_bad_continuity
        integer, allocatable :: map(:)
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(s1, id1)
        call mp_discontinuous%add_patch(s2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()

        call s_mismatch%set_tetragon([1.0_rk, 1.0_rk], [2, 3])
        call mp_discontinuous_mismatch%add_patch(s1, id1)
        call mp_discontinuous_mismatch%add_patch(s_mismatch, id2)
        call mp_discontinuous_mismatch%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous_mismatch%cmp_dof_map()

        call s_bad%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s_bad%translate_Xc([2.0_rk, 0.0_rk, 0.0_rk])
        call mp_bad_geom%add_patch(s1, id1)
        call mp_bad_geom%add_patch(s_bad, id2)
        call mp_bad_geom%connect(id1, "u_max", id2, "u_min", 0)

        call mp_bad_continuity%add_patch(s1, id1)
        call mp_bad_continuity%add_patch(s2, id2)
        call mp_bad_continuity%connect(id1, "u_max", id2, "u_min", 2)
        call s1%err%print()
        call s2%err%print()
        call s_bad%err%print()
        call s_mismatch%err%print()
        call mp_discontinuous%err%print()
        call mp_discontinuous_mismatch%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_continuity%err%print()

        call ut%test(ti)%check(&
            name     = "surface invalid continuity err",&
            res      = mp_bad_continuity%err%ok,&
            expected = .false.,&
            msg      = "invalid surface continuity should set multipatch error state",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0035


    subroutine forcad_multipatch_surface_0036(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2, s_bad, s_mismatch
        type(nurbs_multipatch_surface) :: mp_discontinuous, mp_discontinuous_mismatch
        type(nurbs_multipatch_surface) :: mp_bad_geom
        type(nurbs_multipatch_surface) :: mp_bad_continuity
        integer, allocatable :: map(:)
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(s1, id1)
        call mp_discontinuous%add_patch(s2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()

        call s_mismatch%set_tetragon([1.0_rk, 1.0_rk], [2, 3])
        call mp_discontinuous_mismatch%add_patch(s1, id1)
        call mp_discontinuous_mismatch%add_patch(s_mismatch, id2)
        call mp_discontinuous_mismatch%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous_mismatch%cmp_dof_map()

        call s_bad%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s_bad%translate_Xc([2.0_rk, 0.0_rk, 0.0_rk])
        call mp_bad_geom%add_patch(s1, id1)
        call mp_bad_geom%add_patch(s_bad, id2)
        call mp_bad_geom%connect(id1, "u_max", id2, "u_min", 0)

        call mp_bad_continuity%add_patch(s1, id1)
        call mp_bad_continuity%add_patch(s2, id2)
        call mp_bad_continuity%connect(id1, "u_max", id2, "u_min", 2)
        call s1%err%print()
        call s2%err%print()
        call s_bad%err%print()
        call s_mismatch%err%print()
        call mp_discontinuous%err%print()
        call mp_discontinuous_mismatch%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_continuity%err%print()

        call ut%test(ti)%check(&
            name     = "surface invalid continuity no connection",&
            res      = mp_bad_continuity%get_nconnection(),&
            expected = 0,&
            msg      = "invalid surface continuity should not append a connection",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0036


    subroutine forcad_multipatch_surface_0037(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: multipatch_surface_tol = 1.0e-10_rk
        type(nurbs_surface) :: s1, s2, s_bad, s_mismatch
        type(nurbs_multipatch_surface) :: mp_discontinuous, mp_discontinuous_mismatch
        type(nurbs_multipatch_surface) :: mp_bad_geom
        type(nurbs_multipatch_surface) :: mp_bad_continuity
        integer, allocatable :: map(:)
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(s1, id1)
        call mp_discontinuous%add_patch(s2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()

        call s_mismatch%set_tetragon([1.0_rk, 1.0_rk], [2, 3])
        call mp_discontinuous_mismatch%add_patch(s1, id1)
        call mp_discontinuous_mismatch%add_patch(s_mismatch, id2)
        call mp_discontinuous_mismatch%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous_mismatch%cmp_dof_map()

        call s_bad%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s_bad%translate_Xc([2.0_rk, 0.0_rk, 0.0_rk])
        call mp_bad_geom%add_patch(s1, id1)
        call mp_bad_geom%add_patch(s_bad, id2)
        call mp_bad_geom%connect(id1, "u_max", id2, "u_min", 0)

        call mp_bad_continuity%add_patch(s1, id1)
        call mp_bad_continuity%add_patch(s2, id2)
        call mp_bad_continuity%connect(id1, "u_max", id2, "u_min", 2)
        call s1%err%print()
        call s2%err%print()
        call s_bad%err%print()
        call s_mismatch%err%print()
        call mp_discontinuous%err%print()
        call mp_discontinuous_mismatch%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_continuity%err%print()

        call ut%test(ti)%check(&
            name     = "surface invalid continuity",&
            res      = mp_bad_continuity%is_valid(multipatch_surface_tol),&
            expected = .false.,&
            msg      = "surface multipatch should reject continuity above patch degree",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0037


    subroutine forcad_multipatch_surface_0038(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2, s_bad, s_mismatch
        type(nurbs_multipatch_surface) :: mp_discontinuous, mp_discontinuous_mismatch
        type(nurbs_multipatch_surface) :: mp_bad_geom, mp_bad_continuity, mp_bad_side
        integer, allocatable :: map(:)
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(s1, id1)
        call mp_discontinuous%add_patch(s2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()

        call s_mismatch%set_tetragon([1.0_rk, 1.0_rk], [2, 3])
        call mp_discontinuous_mismatch%add_patch(s1, id1)
        call mp_discontinuous_mismatch%add_patch(s_mismatch, id2)
        call mp_discontinuous_mismatch%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous_mismatch%cmp_dof_map()

        call s_bad%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s_bad%translate_Xc([2.0_rk, 0.0_rk, 0.0_rk])
        call mp_bad_geom%add_patch(s1, id1)
        call mp_bad_geom%add_patch(s_bad, id2)
        call mp_bad_geom%connect(id1, "u_max", id2, "u_min", 0)

        call mp_bad_continuity%add_patch(s1, id1)
        call mp_bad_continuity%add_patch(s2, id2)
        call mp_bad_continuity%connect(id1, "u_max", id2, "u_min", 2)

        call mp_bad_side%add_patch(s1, id1)
        call mp_bad_side%connect(id1, "bad_side", id1, "u_min", 0)
        call s1%err%print()
        call s2%err%print()
        call s_bad%err%print()
        call s_mismatch%err%print()
        call mp_discontinuous%err%print()
        call mp_discontinuous_mismatch%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_continuity%err%print()
        call mp_bad_side%err%print()

        call ut%test(ti)%check(&
            name     = "surface invalid side err",&
            res      = mp_bad_side%err%ok,&
            expected = .false.,&
            msg      = "invalid surface side should set multipatch error state",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0038


    subroutine forcad_multipatch_surface_0039(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2, s_bad, s_mismatch
        type(nurbs_multipatch_surface) :: mp_discontinuous, mp_discontinuous_mismatch
        type(nurbs_multipatch_surface) :: mp_bad_geom, mp_bad_continuity, mp_bad_side
        integer, allocatable :: map(:)
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(s1, id1)
        call mp_discontinuous%add_patch(s2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()

        call s_mismatch%set_tetragon([1.0_rk, 1.0_rk], [2, 3])
        call mp_discontinuous_mismatch%add_patch(s1, id1)
        call mp_discontinuous_mismatch%add_patch(s_mismatch, id2)
        call mp_discontinuous_mismatch%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous_mismatch%cmp_dof_map()

        call s_bad%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s_bad%translate_Xc([2.0_rk, 0.0_rk, 0.0_rk])
        call mp_bad_geom%add_patch(s1, id1)
        call mp_bad_geom%add_patch(s_bad, id2)
        call mp_bad_geom%connect(id1, "u_max", id2, "u_min", 0)

        call mp_bad_continuity%add_patch(s1, id1)
        call mp_bad_continuity%add_patch(s2, id2)
        call mp_bad_continuity%connect(id1, "u_max", id2, "u_min", 2)

        call mp_bad_side%add_patch(s1, id1)
        call mp_bad_side%connect(id1, "bad_side", id1, "u_min", 0)
        call s1%err%print()
        call s2%err%print()
        call s_bad%err%print()
        call s_mismatch%err%print()
        call mp_discontinuous%err%print()
        call mp_discontinuous_mismatch%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_continuity%err%print()
        call mp_bad_side%err%print()

        call ut%test(ti)%check(&
            name     = "surface invalid side no connection",&
            res      = mp_bad_side%get_nconnection(),&
            expected = 0,&
            msg      = "invalid surface side should not append a connection",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0039


    subroutine forcad_multipatch_surface_0040(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2, s_bad, s_mismatch
        type(nurbs_multipatch_surface) :: mp_discontinuous, mp_discontinuous_mismatch
        type(nurbs_multipatch_surface) :: mp_bad_geom, mp_bad_continuity, mp_bad_side
        integer, allocatable :: map(:)
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(s1, id1)
        call mp_discontinuous%add_patch(s2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()

        call s_mismatch%set_tetragon([1.0_rk, 1.0_rk], [2, 3])
        call mp_discontinuous_mismatch%add_patch(s1, id1)
        call mp_discontinuous_mismatch%add_patch(s_mismatch, id2)
        call mp_discontinuous_mismatch%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous_mismatch%cmp_dof_map()

        call s_bad%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s_bad%translate_Xc([2.0_rk, 0.0_rk, 0.0_rk])
        call mp_bad_geom%add_patch(s1, id1)
        call mp_bad_geom%add_patch(s_bad, id2)
        call mp_bad_geom%connect(id1, "u_max", id2, "u_min", 0)

        call mp_bad_continuity%add_patch(s1, id1)
        call mp_bad_continuity%add_patch(s2, id2)
        call mp_bad_continuity%connect(id1, "u_max", id2, "u_min", 2)

        call mp_bad_side%add_patch(s1, id1)
        call mp_bad_side%connect(id1, "bad_side", id1, "u_min", 0)
        call mp_bad_side%add_patch(s2, id2)
        call s1%err%print()
        call s2%err%print()
        call s_bad%err%print()
        call s_mismatch%err%print()
        call mp_discontinuous%err%print()
        call mp_discontinuous_mismatch%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_continuity%err%print()
        call mp_bad_side%err%print()

        call ut%test(ti)%check(&
            name     = "surface add blocked on err",&
            res      = id2,&
            expected = 0,&
            msg      = "Surface add blocked on err is incorrect.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0040


    subroutine forcad_multipatch_surface_0041(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2, s_bad, s_mismatch
        type(nurbs_multipatch_surface) :: mp_discontinuous, mp_discontinuous_mismatch
        type(nurbs_multipatch_surface) :: mp_bad_geom, mp_bad_continuity, mp_bad_side
        integer, allocatable :: map(:)
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%translate_Xc([1.0_rk, 0.0_rk, 0.0_rk])
        call mp_discontinuous%add_patch(s1, id1)
        call mp_discontinuous%add_patch(s2, id2)
        call mp_discontinuous%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous%cmp_dof_map()

        call s_mismatch%set_tetragon([1.0_rk, 1.0_rk], [2, 3])
        call mp_discontinuous_mismatch%add_patch(s1, id1)
        call mp_discontinuous_mismatch%add_patch(s_mismatch, id2)
        call mp_discontinuous_mismatch%connect(id1, "u_max", id2, "u_min", -1)
        map = mp_discontinuous_mismatch%cmp_dof_map()

        call s_bad%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s_bad%translate_Xc([2.0_rk, 0.0_rk, 0.0_rk])
        call mp_bad_geom%add_patch(s1, id1)
        call mp_bad_geom%add_patch(s_bad, id2)
        call mp_bad_geom%connect(id1, "u_max", id2, "u_min", 0)

        call mp_bad_continuity%add_patch(s1, id1)
        call mp_bad_continuity%add_patch(s2, id2)
        call mp_bad_continuity%connect(id1, "u_max", id2, "u_min", 2)

        call mp_bad_side%add_patch(s1, id1)
        call mp_bad_side%connect(id1, "bad_side", id1, "u_min", 0)
        call mp_bad_side%add_patch(s2, id2)
        call mp_bad_side%err%print()
        call mp_bad_side%err%reset()
        call mp_bad_side%add_patch(s2, id2)
        call s1%err%print()
        call s2%err%print()
        call s_bad%err%print()
        call s_mismatch%err%print()
        call mp_discontinuous%err%print()
        call mp_discontinuous_mismatch%err%print()
        call mp_bad_geom%err%print()
        call mp_bad_continuity%err%print()
        call mp_bad_side%err%print()

        call ut%test(ti)%check(&
            name     = "surface add after reset",&
            res      = id2,&
            expected = 2,&
            msg      = "debug reset should allow surface multipatch operations again",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0041


    subroutine forcad_multipatch_surface_0042(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: invalid_patch
        type(nurbs_surface) :: s1
        type(nurbs_multipatch_surface) :: invalid_patch_mp
        integer :: id1

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_patch%err%set(&
            code       = 901,&
            severity   = 1,&
            category   = "test_forcad_multipatch_surface",&
            message    = "Intentional invalid patch state.",&
            location   = "dynamic storage test",&
            suggestion = "Used only to verify add_patch validation.")
        call invalid_patch_mp%add_patch(invalid_patch, id1)
        call invalid_patch%err%print()
        call s1%err%print()
        call invalid_patch_mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface rejects invalid patch",&
            res      = [invalid_patch_mp%get_npatch(), id1],&
            expected = [0, 0],&
            msg      = "An invalid surface patch must not be appended or receive an id.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0042


    subroutine forcad_multipatch_surface_0043(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: invalid_patch
        type(nurbs_surface) :: s1
        type(nurbs_multipatch_surface) :: capacity_mp, invalid_patch_mp
        integer :: i, ids(18), id1

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call invalid_patch%err%set(&
            code       = 901,&
            severity   = 1,&
            category   = "test_forcad_multipatch_surface",&
            message    = "Intentional invalid patch state.",&
            location   = "dynamic storage test",&
            suggestion = "Used only to verify add_patch validation.")
        call invalid_patch_mp%add_patch(invalid_patch, id1)

        do i = 1, size(ids)
            call capacity_mp%add_patch(s1, ids(i))
        end do
        do i = 1, size(ids)-1
            call capacity_mp%connect(ids(i), "u_max", ids(i+1), "u_max", -1)
        end do
        call invalid_patch%err%print()
        call s1%err%print()
        call capacity_mp%err%print()
        call invalid_patch_mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface dynamic storage growth",&
            res      = [capacity_mp%get_npatch(), capacity_mp%get_nconnection()],&
            expected = [18, 17],&
            msg      = "Surface dynamic storage growth is incorrect.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0043


    subroutine forcad_multipatch_surface_0044(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2
        type(nurbs_multipatch_surface) :: mp, mp_copy
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        mp_copy = mp
        call mp%finalize()
        call s1%err%print()
        call s2%err%print()
        call mp%err%print()
        call mp_copy%err%print()

        call ut%test(ti)%check(&
            name     = "surface finalize",&
            res      = mp%get_npatch(),&
            expected = 0,&
            msg      = "surface multipatch finalize failed",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0044


    subroutine forcad_multipatch_surface_0045(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: s1, s2
        type(nurbs_multipatch_surface) :: mp, mp_copy
        integer :: id1, id2

        call s1%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call s2%set_tetragon([1.0_rk, 1.0_rk], [2, 2])
        call mp%add_patch(s1, id1)
        call mp%add_patch(s2, id2)
        mp_copy = mp
        call mp%finalize()
        call s1%err%print()
        call s2%err%print()
        call mp%err%print()
        call mp_copy%err%print()

        call ut%test(ti)%check(&
            name     = "surface multipatch assignment",&
            res      = mp_copy%get_npatch(),&
            expected = 2,&
            msg      = "Surface multipatch assignment is incorrect.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0045


    subroutine forcad_multipatch_surface_0046(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 1.0e-12_rk
        real(rk), parameter :: knot_normal(7) = [0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk]
        real(rk), parameter :: knot_tangent(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xca(8,2) = reshape([&
            -4.0_rk,-2.0_rk,0.0_rk,2.0_rk,-4.0_rk,-2.0_rk,0.0_rk,2.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk], [8,2])
        real(rk), parameter :: Xcb(8,2) = reshape([&
            -1.0_rk,4.0_rk,5.0_rk,6.0_rk,-1.0_rk,4.0_rk,5.0_rk,6.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk], [8,2])
        real(rk), parameter :: Wa(8) = [2.0_rk,2.0_rk,1.0_rk,3.0_rk,2.0_rk,2.0_rk,1.0_rk,3.0_rk]
        real(rk), parameter :: Wb(8) = [2.0_rk,2.0_rk,3.0_rk,4.0_rk,2.0_rk,2.0_rk,3.0_rk,4.0_rk]
        integer, parameter :: degree(2) = [2,1]
        type(nurbs_surface) :: sa, sb
        type(nurbs_multipatch_surface) :: mp
        integer :: ida, idb

        call sa%set(&
            knot1  = knot_normal,&
            knot2  = knot_tangent,&
            Xc     = Xca,&
            Wc     = Wa,&
            degree = degree)
        call sb%set(&
            knot1  = knot_normal,&
            knot2  = knot_tangent,&
            Xc     = Xcb,&
            Wc     = Wb,&
            degree = degree)
        call mp%add_patch(sa, ida)
        call mp%add_patch(sb, idb)
        call mp%connect(ida, "u_max", idb, "u_min", 0)
        call sa%err%print()
        call sb%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface unclamped C0 trace",&
            res      = mp%is_valid(tol),&
            expected = .true.,&
            msg      = "Basis-weighted surface edge traces must match.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0046


    subroutine forcad_multipatch_surface_0047(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: knot_normal(7) = [0.0_rk,1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk]
        real(rk), parameter :: knot_tangent(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: Xca(8,2) = reshape([&
            -4.0_rk,-2.0_rk,0.0_rk,2.0_rk,-4.0_rk,-2.0_rk,0.0_rk,2.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk], [8,2])
        real(rk), parameter :: Xcb(8,2) = reshape([&
            -2.0_rk,4.0_rk,5.0_rk,6.0_rk,-2.0_rk,4.0_rk,5.0_rk,6.0_rk,&
            0.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk,1.0_rk], [8,2])
        integer, parameter :: degree(2) = [2,1]
        type(nurbs_surface) :: sa, sb
        type(nurbs_multipatch_surface) :: mp
        integer, allocatable :: map(:)
        integer :: ida, idb

        call sa%set(knot_normal, knot_tangent, Xca, degree=degree)
        call sb%set(knot_normal, knot_tangent, Xcb, degree=degree)
        call mp%add_patch(sa, ida)
        call mp%add_patch(sb, idb)
        call mp%connect(ida, "u_max", idb, "u_min", 0)
        map = mp%cmp_dof_map()
        call sa%err%print()
        call sb%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface unclamped C0 map",&
            res      = map,&
            expected = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],&
            msg      = "An unclamped surface trace must remain a constraint.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0047


    subroutine forcad_multipatch_surface_0048(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_multipatch_surface) :: mp

        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface empty multipatch invalid",&
            res      = mp%is_valid(),&
            expected = .false.,&
            msg      = "An empty surface multipatch is not geometry.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0048


    subroutine forcad_multipatch_surface_0049(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: surface
        type(nurbs_multipatch_surface) :: mp
        integer :: id

        call mp%add_patch(surface, id)
        call surface%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface rejects uninitialized patch",&
            res      = id,&
            expected = 0,&
            msg      = "An uninitialized surface must not be added.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0049


    subroutine forcad_multipatch_surface_0050(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        type(nurbs_surface) :: surface
        type(nurbs_multipatch_surface) :: mp
        integer :: id

        call surface%set(&
            Xth_dir1    = [0.0_rk,1.0_rk],&
            Xth_dir2    = [0.0_rk,1.0_rk],&
            degree      = [1,1],&
            continuity1 = [-1,-1],&
            continuity2 = [-1,-1])
        call mp%add_patch(surface, id)
        call surface%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface rejects space-only patch",&
            res      = id,&
            expected = 0,&
            msg      = "A surface spline space is not complete geometry.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0050


    subroutine forcad_multipatch_surface_0051(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 1.0e-10_rk
        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        type(nurbs_surface) :: surface_a, surface_b
        type(nurbs_multipatch_surface) :: mp
        real(rk) :: Xc_a(4,3), Xc_b(4,3), transition_jet(2,2,1)
        integer :: id_a, id_b

        Xc_a = 0.0_rk
        Xc_b = 0.0_rk
        Xc_a(:,1) = [-0.2_rk,0.0_rk,0.5_rk,1.0_rk]
        Xc_a(:,2) = [-0.5_rk,0.0_rk,-0.7_rk,0.0_rk]
        Xc_b(:,1) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        Xc_b(:,2) = [0.0_rk,1.0_rk,0.0_rk,1.0_rk]
        transition_jet(:,1,1) = [0.2_rk,0.5_rk]
        transition_jet(:,2,1) = [0.5_rk,0.7_rk]
        call surface_a%set(knot, knot, Xc_a, degree=[1,1])
        call surface_b%set(knot, knot, Xc_b, degree=[1,1])
        call mp%add_patch(surface_a, id_a)
        call mp%add_patch(surface_b, id_b)
        call mp%connect(&
            patch_a        = id_a,&
            side_a         = "u_max",&
            patch_b        = id_b,&
            side_b         = "u_min",&
            continuity     = 1,&
            geometric      = .true.,&
            transition_jet = transition_jet)
        call surface_a%err%print()
        call surface_b%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface general G1",&
            res      = mp%is_valid(tol),&
            expected = .true.,&
            msg      = "Tangentially varying surface G1 must be valid.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0051


    subroutine forcad_multipatch_surface_0052(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: knot(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        type(nurbs_surface) :: surface_a, surface_b
        type(nurbs_multipatch_surface) :: mp
        real(rk) :: Xc_a(4,3), Xc_b(4,3), transition_jet(2,2,1)
        real(rk), allocatable :: val(:)
        integer, allocatable :: rowptr(:), col(:)
        integer :: id_a, id_b

        Xc_a = 0.0_rk
        Xc_b = 0.0_rk
        Xc_a(:,1) = [-0.2_rk,0.0_rk,0.5_rk,1.0_rk]
        Xc_a(:,2) = [-0.5_rk,0.0_rk,-0.7_rk,0.0_rk]
        Xc_b(:,1) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        Xc_b(:,2) = [0.0_rk,1.0_rk,0.0_rk,1.0_rk]
        transition_jet(:,1,1) = [0.2_rk,0.5_rk]
        transition_jet(:,2,1) = [0.5_rk,0.7_rk]
        call surface_a%set(knot, knot, Xc_a, degree=[1,1])
        call surface_b%set(knot, knot, Xc_b, degree=[1,1])
        call mp%add_patch(surface_a, id_a)
        call mp%add_patch(surface_b, id_b)
        call mp%connect(&
            patch_a        = id_a,&
            side_a         = "u_max",&
            patch_b        = id_b,&
            side_b         = "u_min",&
            continuity     = 1,&
            geometric      = .true.,&
            transition_jet = transition_jet)
        call mp%cmp_dof_constraint(rowptr, col, val, geometric=.true.)
        call surface_a%err%print()
        call surface_b%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface general G1 rows",&
            res      = size(rowptr) - 1,&
            expected = 5,&
            msg      = "Surface general G1 collocation degree is incorrect.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0052


    subroutine forcad_multipatch_surface_0053(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 1.0e-10_rk
        real(rk), parameter :: knot_linear(4) = [0.0_rk,0.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: knot_quadratic(6) = [0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk]
        type(nurbs_surface) :: surface_a, surface_b
        type(nurbs_multipatch_surface) :: mp
        real(rk), allocatable :: val(:)
        real(rk) :: Xc_a(6,3), Xc_b(6,3), Wc_a(6), Wc_b(6), projective_jet(2,1)
        integer, allocatable :: map(:), rowptr(:), col(:)
        integer :: id_a, id_b

        Xc_a = 0.0_rk
        Xc_b = 0.0_rk
        Xc_a(:,1) = [0.0_rk,0.0_rk,1.0_rk/3.0_rk,1.0_rk/3.0_rk,1.0_rk,1.0_rk]
        Xc_a(:,2) = [-1.0_rk,0.0_rk,-1.0_rk,0.0_rk,-1.0_rk,0.0_rk]
        Xc_b(:,1) = [0.0_rk,0.0_rk,0.5_rk,0.5_rk,1.0_rk,1.0_rk]
        Xc_b(:,2) = [0.0_rk,1.0_rk,0.0_rk,1.0_rk,0.0_rk,1.0_rk]
        Wc_a = [1.0_rk,1.0_rk,1.5_rk,1.5_rk,2.0_rk,2.0_rk]
        Wc_b = [1.0_rk,1.1_rk,1.0_rk,1.2_rk,1.0_rk,1.3_rk]
        projective_jet(:,1) = [1.0_rk,2.0_rk]
        call surface_a%set(knot_linear, knot_quadratic, Xc_a, Wc_a, degree=[1,2])
        call surface_b%set(knot_linear, knot_quadratic, Xc_b, Wc_b, degree=[1,2])
        call mp%add_patch(surface_a, id_a)
        call mp%add_patch(surface_b, id_b)
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
        call surface_a%err%print()
        call surface_b%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface variable projective constraint",&
            res      = surface_a%is_rational() .and. surface_b%is_rational() .and. &
                mp%is_valid(tol) .and. maxval(map) == size(map) .and. &
                size(rowptr) > 1 .and. size(col) == size(val) .and. size(val) > 0,&
            expected = .true.,&
            msg      = "Variable projective trace bases must remain sparse constraints.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0053


    subroutine forcad_multipatch_surface_0054(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 1.0e-10_rk
        real(rk), parameter :: knot(6) = [0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk]
        real(rk), parameter :: normal_weight(3) = [0.8_rk,0.9_rk,1.0_rk]
        real(rk), parameter :: normal_coordinate(3) = [-1.0_rk,-5.0_rk/9.0_rk,0.0_rk]
        real(rk), parameter :: tangent_a(3) = [0.0_rk,1.0_rk/3.0_rk,1.0_rk]
        real(rk), parameter :: tangent_b(3) = [0.0_rk,0.5_rk,1.0_rk]
        real(rk), parameter :: tangent_weight(3) = [1.0_rk,1.5_rk,2.0_rk]
        type(nurbs_surface) :: surface_a, surface_b
        type(nurbs_multipatch_surface) :: mp
        real(rk) :: Xc_a(9,3), Xc_b(9,3), Wc_a(9), Wc_b(9)
        real(rk) :: projective_jet(3,2)
        integer :: id_a, id_b, i, j, index

        Xc_a = 0.0_rk
        Xc_b = 0.0_rk
        do j = 1, 3
            do i = 1, 3
                index = i + 3*(j-1)
                Xc_a(index,1:2) = [tangent_a(j),normal_coordinate(i)]
                Xc_b(index,1:2) = [tangent_b(j),0.5_rk*real(i-1,rk)]
                Wc_a(index) = normal_weight(i)*tangent_weight(j)
                Wc_b(index) = 1.0_rk
            end do
        end do
        projective_jet(:,1) = tangent_weight
        projective_jet(:,2) = 0.2_rk*tangent_weight
        call surface_a%set(knot, knot, Xc_a, Wc_a, degree=[2,2])
        call surface_b%set(knot, knot, Xc_b, Wc_b, degree=[2,2])
        call mp%add_patch(surface_a, id_a)
        call mp%add_patch(surface_b, id_b)
        call mp%connect(&
            patch_a            = id_a,&
            side_a             = "u_max",&
            patch_b            = id_b,&
            side_b             = "u_min",&
            continuity         = 1,&
            geometric          = .true.,&
            reparameterization = [1.0_rk],&
            projective_jet     = projective_jet)
        call surface_a%err%print()
        call surface_b%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface variable projective G1",&
            res      = mp%is_valid(tol),&
            expected = .true.,&
            msg      = "Variable projective surface G1 must include product terms.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0054


    subroutine forcad_multipatch_surface_0055(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        real(rk), parameter :: tol = 1.0e-10_rk
        real(rk), parameter :: knot(6) = [0.0_rk,0.0_rk,0.0_rk,1.0_rk,1.0_rk,1.0_rk]
        type(nurbs_surface) :: surface_a, surface_b
        type(nurbs_multipatch_surface) :: mp
        real(rk) :: Xc_a(9,3), Xc_b(9,3), derivative0(3,3), derivative1(3,3), derivative2(3,3)
        real(rk) :: transition_jet(2,2,2)
        integer :: id_a, id_b, i, j, index

        derivative0 = 0.0_rk
        derivative1 = 0.0_rk
        derivative2 = 0.0_rk
        derivative0(:,1) = [0.0_rk,0.5_rk,1.0_rk]
        derivative1(:,1) = [0.2_rk,0.25_rk,0.3_rk]
        derivative1(:,2) = [0.5_rk,0.55_rk,0.6_rk]
        derivative1(:,3) = [0.0_rk,0.25_rk,0.6_rk]
        derivative2(:,1) = [-0.1_rk,-0.075_rk,-0.05_rk]
        derivative2(:,2) = [0.3_rk,0.25_rk,0.2_rk]
        derivative2(:,3) = [0.2_rk,0.42_rk,0.56_rk]
        transition_jet(:,1,1) = [0.2_rk,0.3_rk]
        transition_jet(:,2,1) = [0.5_rk,0.6_rk]
        transition_jet(:,1,2) = [-0.1_rk,-0.05_rk]
        transition_jet(:,2,2) = [0.3_rk,0.2_rk]
        do j = 0, 2
            do i = 0, 2
                index = i + 1 + 3*j
                Xc_b(index,:) = [0.5_rk*real(j,rk),0.5_rk*real(i,rk),0.25_rk*real(i*j,rk)]
                Xc_a(index,:) = derivative0(j+1,:) - &
                    (1.0_rk - 0.5_rk*real(i,rk))*derivative1(j+1,:)
                if (i == 0) Xc_a(index,:) = Xc_a(index,:) + 0.5_rk*derivative2(j+1,:)
            end do
        end do
        call surface_a%set(knot, knot, Xc_a, degree=[2,2])
        call surface_b%set(knot, knot, Xc_b, degree=[2,2])
        call mp%add_patch(surface_a, id_a)
        call mp%add_patch(surface_b, id_b)
        call mp%connect(&
            patch_a        = id_a,&
            side_a         = "u_max",&
            patch_b        = id_b,&
            side_b         = "u_min",&
            continuity     = 2,&
            geometric      = .true.,&
            transition_jet = transition_jet)
        call surface_a%err%print()
        call surface_b%err%print()
        call mp%err%print()

        call ut%test(ti)%check(&
            name     = "surface general G2",&
            res      = mp%is_valid(tol),&
            expected = .true.,&
            msg      = "Surface G2 must include varying mixed transition terms.",&
            group    = "forcad_multipatch_surface")
        ti = ti + 1

    end subroutine forcad_multipatch_surface_0055


    subroutine run_multipatch_surface_tests(ut, ti)
        type(unit_tests), intent(inout) :: ut
        integer, intent(inout) :: ti

        call forcad_multipatch_surface_0001(ut, ti)
        call forcad_multipatch_surface_0002(ut, ti)
        call forcad_multipatch_surface_0003(ut, ti)
        call forcad_multipatch_surface_0004(ut, ti)
        call forcad_multipatch_surface_0005(ut, ti)
        call forcad_multipatch_surface_0006(ut, ti)
        call forcad_multipatch_surface_0007(ut, ti)
        call forcad_multipatch_surface_0008(ut, ti)
        call forcad_multipatch_surface_0009(ut, ti)
        call forcad_multipatch_surface_0010(ut, ti)
        call forcad_multipatch_surface_0011(ut, ti)
        call forcad_multipatch_surface_0012(ut, ti)
        call forcad_multipatch_surface_0013(ut, ti)
        call forcad_multipatch_surface_0014(ut, ti)
        call forcad_multipatch_surface_0015(ut, ti)
        call forcad_multipatch_surface_0016(ut, ti)
        call forcad_multipatch_surface_0017(ut, ti)
        call forcad_multipatch_surface_0018(ut, ti)
        call forcad_multipatch_surface_0019(ut, ti)
        call forcad_multipatch_surface_0020(ut, ti)
        call forcad_multipatch_surface_0021(ut, ti)
        call forcad_multipatch_surface_0022(ut, ti)
        call forcad_multipatch_surface_0023(ut, ti)
        call forcad_multipatch_surface_0024(ut, ti)
        call forcad_multipatch_surface_0025(ut, ti)
        call forcad_multipatch_surface_0026(ut, ti)
        call forcad_multipatch_surface_0027(ut, ti)
        call forcad_multipatch_surface_0028(ut, ti)
        call forcad_multipatch_surface_0029(ut, ti)
        call forcad_multipatch_surface_0030(ut, ti)
        call forcad_multipatch_surface_0031(ut, ti)
        call forcad_multipatch_surface_0032(ut, ti)
        call forcad_multipatch_surface_0033(ut, ti)
        call forcad_multipatch_surface_0034(ut, ti)
        call forcad_multipatch_surface_0035(ut, ti)
        call forcad_multipatch_surface_0036(ut, ti)
        call forcad_multipatch_surface_0037(ut, ti)
        call forcad_multipatch_surface_0038(ut, ti)
        call forcad_multipatch_surface_0039(ut, ti)
        call forcad_multipatch_surface_0040(ut, ti)
        call forcad_multipatch_surface_0041(ut, ti)
        call forcad_multipatch_surface_0042(ut, ti)
        call forcad_multipatch_surface_0043(ut, ti)
        call forcad_multipatch_surface_0044(ut, ti)
        call forcad_multipatch_surface_0045(ut, ti)
        call forcad_multipatch_surface_0046(ut, ti)
        call forcad_multipatch_surface_0047(ut, ti)
        call forcad_multipatch_surface_0048(ut, ti)
        call forcad_multipatch_surface_0049(ut, ti)
        call forcad_multipatch_surface_0050(ut, ti)
        call forcad_multipatch_surface_0051(ut, ti)
        call forcad_multipatch_surface_0052(ut, ti)
        call forcad_multipatch_surface_0053(ut, ti)
        call forcad_multipatch_surface_0054(ut, ti)
        call forcad_multipatch_surface_0055(ut, ti)

    end subroutine run_multipatch_surface_tests

end module forcad_test_multipatch_surface
