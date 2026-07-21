!===============================================================================
!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Geometry-only six-port rounded joint with tapered curved arms.
program volume_multipatch_joint

    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume

    implicit none

    integer, parameter :: nseg = 5, ndir = 6, ncp = 3
    real(rk), parameter :: pi = acos(-1.0_rk)
    real(rk), parameter :: core = 0.55_rk
    real(rk), parameter :: arm_length = 1.75_rk

    type(nurbs_volume) :: patch
    type(nurbs_multipatch_volume) :: joint
    integer :: core_id, id(nseg, ndir), is, idir
    integer, allocatable :: dof_map(:), elem(:,:)

    call set_core_patch(patch)
    call joint%add_patch(patch, core_id)

    do idir = 1, ndir
        do is = 1, nseg
            call set_arm_patch(patch, idir, is)
            call joint%add_patch(patch, id(is, idir))
        end do
    end do

    call joint%connect(core_id, "u_max", id(1, 1), "u_min", 0)
    call joint%connect(core_id, "u_min", id(1, 2), "u_min", 0)
    call joint%connect(core_id, "v_max", id(1, 3), "u_min", 0)
    call joint%connect(core_id, "v_min", id(1, 4), "u_min", 0)
    call joint%connect(core_id, "w_max", id(1, 5), "u_min", 0)
    call joint%connect(core_id, "w_min", id(1, 6), "u_min", 0)

    do idir = 1, ndir
        do is = 1, nseg - 1
            call joint%connect(id(is, idir), "u_max", id(is + 1, idir), "u_min", 0)
        end do
    end do

    dof_map = joint%cmp_dof_map()
    elem = joint%cmp_elem(shared=.true.)

    write(*,"(a,*(1x,g0))") "joint patches         :", joint%get_npatch()
    write(*,"(a,*(1x,g0))") "joint connections     :", joint%get_nconnection()
    write(*,"(a,*(1x,g0))") "joint valid           :", joint%is_valid()
    write(*,"(a,*(1x,g0))") "shared global dofs    :", maxval(dof_map)
    write(*,"(a,*(1x,g0))") "elements              :", size(elem, 1)
    write(*,"(a,*(1x,g0))") "nodes per element     :", size(elem, 2)

    call joint%create(6, 5, 5)
    call joint%export_Xc("vtk/volume_multipatch_joint")
    call joint%export_Xg("vtk/volume_multipatch_joint")
    call joint%export_Xth_in_Xg("vtk/volume_multipatch_joint", 6)
    call joint%show(&
        vtkfile_Xc        = "vtk/volume_multipatch_joint_Xc_*.vtk",&
        vtkfile_Xg        = "vtk/volume_multipatch_joint_Xg_*.vtk",&
        vtkfile_Xth_in_Xg = "vtk/volume_multipatch_joint_Xth_*.vtk")
    call joint%finalize()

contains

    pure subroutine set_core_patch(p)
        type(nurbs_volume), intent(inout) :: p
        real(rk) :: Xc(ncp*ncp*ncp, 3)
        real(rk) :: u, v, w, x, y, z
        integer :: a, b, c, n

        do c = 1, ncp
            w = 2.0_rk*real(c - 1, rk)/real(ncp - 1, rk) - 1.0_rk
            do b = 1, ncp
                v = 2.0_rk*real(b - 1, rk)/real(ncp - 1, rk) - 1.0_rk
                do a = 1, ncp
                    u = 2.0_rk*real(a - 1, rk)/real(ncp - 1, rk) - 1.0_rk
                    n = a + (b - 1)*ncp + (c - 1)*ncp*ncp
                    x = core*u*(1.0_rk + 0.12_rk*(v*v + w*w))
                    y = core*v*(1.0_rk + 0.12_rk*(u*u + w*w))
                    z = core*w*(1.0_rk + 0.12_rk*(u*u + v*v))
                    Xc(n, :) = [x, y, z]
                end do
            end do
        end do
        call p%set([ncp, ncp, ncp], Xc)
    end subroutine
    !===============================================================================

    pure subroutine set_arm_patch(p, idir, is)
        type(nurbs_volume), intent(inout) :: p
        integer, intent(in) :: idir, is
        real(rk) :: Xc(ncp*ncp*ncp, 3)
        real(rk) :: normal(3), tangent1(3), tangent2(3)
        real(rk) :: u, v, w, s, p1, p2, axial, cross1, cross2, bend1, bend2, taper1, taper2
        integer :: a, b, c, d, n

        normal = 0.0_rk
        tangent1 = 0.0_rk
        tangent2 = 0.0_rk
        select case(idir)
        case(1)
            normal = [1.0_rk, 0.0_rk, 0.0_rk]
            tangent1 = [0.0_rk, 1.0_rk, 0.0_rk]
            tangent2 = [0.0_rk, 0.0_rk, 1.0_rk]
        case(2)
            normal = [-1.0_rk, 0.0_rk, 0.0_rk]
            tangent1 = [0.0_rk, 1.0_rk, 0.0_rk]
            tangent2 = [0.0_rk, 0.0_rk, 1.0_rk]
        case(3)
            normal = [0.0_rk, 1.0_rk, 0.0_rk]
            tangent1 = [1.0_rk, 0.0_rk, 0.0_rk]
            tangent2 = [0.0_rk, 0.0_rk, 1.0_rk]
        case(4)
            normal = [0.0_rk, -1.0_rk, 0.0_rk]
            tangent1 = [1.0_rk, 0.0_rk, 0.0_rk]
            tangent2 = [0.0_rk, 0.0_rk, 1.0_rk]
        case(5)
            normal = [0.0_rk, 0.0_rk, 1.0_rk]
            tangent1 = [1.0_rk, 0.0_rk, 0.0_rk]
            tangent2 = [0.0_rk, 1.0_rk, 0.0_rk]
        case default
            normal = [0.0_rk, 0.0_rk, -1.0_rk]
            tangent1 = [1.0_rk, 0.0_rk, 0.0_rk]
            tangent2 = [0.0_rk, 1.0_rk, 0.0_rk]
        end select

        do c = 1, ncp
            w = 2.0_rk*real(c - 1, rk)/real(ncp - 1, rk) - 1.0_rk
            do b = 1, ncp
                v = 2.0_rk*real(b - 1, rk)/real(ncp - 1, rk) - 1.0_rk
                do a = 1, ncp
                    u = real(a - 1, rk)/real(ncp - 1, rk)
                    s = (real(is - 1, rk) + u)/real(nseg, rk)
                    p1 = v
                    p2 = w
                    axial = core*(1.0_rk + 0.12_rk*(p1*p1 + p2*p2)) + arm_length*s
                    taper1 = 1.0_rk + 0.32_rk*s + 0.12_rk*sin(pi*s)
                    taper2 = 1.0_rk + 0.18_rk*s
                    cross1 = core*p1*(1.0_rk + 0.12_rk*(1.0_rk + p2*p2))*taper1
                    cross2 = core*p2*(1.0_rk + 0.12_rk*(1.0_rk + p1*p1))*taper2
                    bend1 = 0.14_rk*sin(pi*s)*sin(real(idir, rk))
                    bend2 = 0.10_rk*sin(2.0_rk*pi*s)*cos(real(idir, rk))
                    n = a + (b - 1)*ncp + (c - 1)*ncp*ncp
                    do d = 1, 3
                        Xc(n, d) = axial*normal(d) + (cross1 + bend1)*tangent1(d) + (cross2 + bend2)*tangent2(d)
                    end do
                end do
            end do
        end do
        call p%set([ncp, ncp, ncp], Xc)
    end subroutine
    !===============================================================================

end program volume_multipatch_joint
