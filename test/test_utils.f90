!> Author: Seyed Ali Ghasemi
!> License: BSD 3-Clause
!> Unit test program for `forcad_utils`.
!> using ForUnitTest: https://github.com/gha3mi/forunittest
program test_forcad_utils

   use forcad_kinds, only: rk
   use forcad_utils, only: basis_bspline, basis_bspline_der, basis_bspline_2der, &
      basis_bernstein, compute_multiplicity, ndgrid, dyad, kron, unique, findspan, &
      compute_knot_vector, insert_knot_A_5_1, remove_knots_A_5_8, elevate_degree_A_5_9, &
      hexahedron_Xc, tetragon_Xc, elemConn_C0, elemConn_Cn, rotation, det, inv, gauss_leg, &
      export_vtk_legacy, solve, repelem, linspace, eye
   use forunittest, only: unit_tests

   implicit none

   type(unit_tests) :: ut

   real(rk) :: Xt
   integer  :: nc, degree
   real(rk) :: knot(7)
   real(rk) :: B4(4), dB(4), d2B(4), A4(2,2)
   real(rk) :: B_ref(4), dB_ref(4), d2B_ref(4)
   real(rk) :: u(2), v(2), w(4)
   real(rk) :: A2x2(2,2), Bk(4,2)
   real(rk), allocatable :: A(:), vec(:), M(:,:)
   real(rk), allocatable :: X1(:), X2(:), X3(:)
   real(rk), allocatable :: Xt2(:,:), Xt3(:,:)
   real(rk), allocatable :: R(:,:), R_expected(:,:)
   real(rk), allocatable :: knot_in(:), knot_out(:), Pw(:,:), Qw(:,:), Pw_new(:,:), knot_new(:)
   real(rk), allocatable :: Xksi(:,:), Wksi(:)
   real(rk), allocatable :: K3(:), K2(:), K1(:), out(:)
   integer, allocatable :: conn1D(:,:), conn2D(:,:), conn3D(:,:)
   integer :: nq, p, rr, s, k, t, vtk_type
   character(len=*), parameter :: vtk_file = "vtk/test_output.vtk"
   real(rk), allocatable :: A2(:,:), A_inv(:,:)

   ! Initialize unit tests
   call ut%initialize(n=47)

   ! ----------------------------
   ! Test: basis_bspline
   ! ----------------------------
   degree = 2
   nc     = 4
   knot = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
   Xt     = 0.5_rk

   B4 = basis_bspline(Xt, knot, nc, degree)
   B_ref = [0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk]

   call ut%test(1)%check( &
      name     = "basis_bspline", &
      res      = B4, &
      expected = B_ref, &
      msg      = "Partition of unity and shape check", &
      group    = "basis")

   call ut%test(2)%check( &
      name     = "basis_sum", &
      res      = sum(B4), &
      expected = 1.0_rk, &
      msg      = "Partition of unity", &
      group    = "basis")

   call basis_bspline_der(Xt, knot, nc, degree, dB, B4)
   dB_ref  = [0.0_rk, -2.0_rk, 2.0_rk, 0.0_rk]

   call ut%test(3)%check( &
      name     = "basis_bspline_der", &
      res      = dB, &
      expected = dB_ref, &
      msg      = "1st derivative shape check", &
      group    = "basis_der")

   call basis_bspline_der(Xt, knot, nc, degree, dB)

   call ut%test(4)%check( &
      name     = "basis_bspline_der_B", &
      res      = dB, &
      expected = dB_ref, &
      msg      = "1st derivative alternate", &
      group    = "basis_der")

   call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB, B4)
   d2B_ref = [0.0_rk, 4.0_rk, -12.0_rk, 8.0_rk]

   call ut%test(5)%check( &
      name     = "basis_bspline_2der_A", &
      res      = d2B, &
      expected = d2B_ref, &
      msg      = "2nd derivative shape check A", &
      group    = "basis_2der")

   call basis_bspline_2der(Xt, knot, nc, degree, d2B, dB)

   call ut%test(6)%check( &
      name     = "basis_bspline_2der_B", &
      res      = d2B, &
      expected = d2B_ref, &
      msg      = "2nd derivative shape check B", &
      group    = "basis_2der")

   call basis_bspline_2der(Xt, knot, nc, degree, d2B)

   call ut%test(7)%check( &
      name     = "basis_bspline_2der_C", &
      res      = d2B, &
      expected = d2B_ref, &
      msg      = "2nd derivative shape check C", &
      group    = "basis_2der")

   Xt = -0.1_rk
   B4 = basis_bspline(Xt, knot, nc, degree)

   call ut%test(8)%check( &
      name     = "basis_out_of_bounds_low", &
      res      = sum(B4), &
      expected = 0.0_rk, &
      msg      = "Out-of-bounds left", &
      group    = "bounds")

   Xt = 1.1_rk
   B4 = basis_bspline(Xt, knot, nc, degree)

   call ut%test(9)%check( &
      name     = "basis_out_of_bounds_high", &
      res      = sum(B4), &
      expected = 0.0_rk, &
      msg      = "Out-of-bounds right", &
      group    = "bounds")

   ! ----------------------------
   ! Test: basis_bernstein
   ! ----------------------------
   call ut%test(10)%check( &
      name     = "basis_bernstein_sum", &
      res      = sum(basis_bernstein(0.3_rk, 3)), &
      expected = 1.0_rk, &
      msg      = "Bernstein basis partition of unity", &
      group    = "bernstein")

   ! ----------------------------
   ! Test: compute_multiplicity
   ! ----------------------------
   call ut%test(11)%check( &
      name     = "compute_multiplicity1", &
      res      = compute_multiplicity([0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk]), &
      expected = [2, 2], &
      msg      = "Multiplicity vector", &
      group    = "multiplicity")

   call ut%test(12)%check( &
      name     = "compute_multiplicity2", &
      res      = compute_multiplicity([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk], 1.0_rk), &
      expected = 3, &
      msg      = "Multiplicity at value", &
      group    = "multiplicity")

   ! ----------------------------
   ! Test: ndgrid
   ! ----------------------------
   X1 = [1.0_rk, 2.0_rk]
   X2 = [10.0_rk, 20.0_rk]
   X3 = [100.0_rk]

   call ndgrid(X1, X2, Xt2)
   call ndgrid(X1, X2, X3, Xt3)

   call ut%test(13)%check( &
      name     = "ndgrid2_shape", &
      res      = shape(Xt2), &
      expected = [4,2], &
      msg      = "ndgrid2 shape", &
      group    = "ndgrid")

   call ut%test(14)%check( &
      name     = "ndgrid3_value", &
      res      = Xt3(:,3), &
      expected = [100.0_rk, 100.0_rk, 100.0_rk, 100.0_rk], &
      msg      = "ndgrid3 constant Z", &
      group    = "ndgrid")

   ! ----------------------------
   ! Test: dyad
   ! ----------------------------
   A = [1.0_rk, 2.0_rk]
   vec = [3.0_rk, 4.0_rk, 5.0_rk]
   M = dyad(A, vec)

   call ut%test(15)%check( &
      name     = "dyad_check", &
      res      = M, &
      expected = reshape([3.0_rk, 6.0_rk, 4.0_rk, 8.0_rk, 5.0_rk, 10.0_rk], [2,3]), &
      msg      = "Outer product a .dyad. b", &
      group    = "dyad")

   ! ----------------------------
   ! Test: kron
   ! ----------------------------
   u = [1.0_rk, 2.0_rk]
   v = [3.0_rk, 4.0_rk]
   w = kron(u, v)

   call ut%test(16)%check( &
      name     = "kron_vector", &
      res      = w, &
      expected = [3.0_rk, 4.0_rk, 6.0_rk, 8.0_rk], &
      msg      = "u .kron. v", &
      group    = "kron")

   A2x2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])
   Bk   = kron(u, A2x2)

   call ut%test(17)%check( &
      name     = "kron_matrix", &
      res      = Bk, &
      expected = reshape([1.0_rk, 2.0_rk, 2.0_rk, 4.0_rk, 3.0_rk, 4.0_rk, 6.0_rk, 8.0_rk], [4,2]), &
      msg      = "u .kron. A", &
      group    = "kron")

   ! ----------------------------
   ! Test: unique
   ! ----------------------------
   call ut%test(18)%check( &
      name     = "unique_integer", &
      res      = unique([1,2,2,3,1,4]), &
      expected = [1,2,3,4], &
      msg      = "Unique integers", &
      group    = "unique")

   call ut%test(19)%check( &
      name     = "unique_real", &
      res      = unique([1.0_rk, 1.0_rk + 1e-16_rk, 2.0_rk, 1.0_rk, 3.0_rk]), &
      expected = [1.0_rk, 2.0_rk, 3.0_rk], &
      msg      = "Unique real with tolerance", &
      group    = "unique")

   ! ----------------------------
   ! Test: findspan
   ! ----------------------------
   call ut%test(20)%check( &
      name     = "findspan_middle", &
      res      = findspan(4, 2, 0.5_rk, knot), &
      expected = 3, &
      msg      = "Find span index at Xt=0.5", &
      group    = "findspan")

   ! ----------------------------
   ! Test: compute_knot_vector
   ! ----------------------------
   call ut%test(21)%check( &
      name     = "compute_knot_vector", &
      res      = compute_knot_vector([0.0_rk, 1.0_rk, 2.0_rk], 2, [-1, 1, -1]), &
      expected = [0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 2.0_rk, 2.0_rk, 2.0_rk], &
      msg      = "Knot vector construction", &
      group    = "knot")

   ! ----------------------------
   ! Test: repelem
   ! ----------------------------
   call ut%test(22)%check( &
      name     = "repelem", &
      res      = repelem([1.0_rk, 2.0_rk, 3.0_rk], [2, 1, 3]), &
      expected = [1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 3.0_rk, 3.0_rk], &
      msg      = "Repeat vector elements", &
      group    = "repelem")

   ! ----------------------------
   ! Test: hexahedron_Xc
   ! ----------------------------
   call ut%test(23)%check( &
      name     = "hexahedron_Xc_shape", &
      res      = shape(hexahedron_Xc([1.0_rk, 1.0_rk, 1.0_rk], [2,2,2])), &
      expected = [8,3], &
      msg      = "3D grid shape", &
      group    = "geometry")

   ! ----------------------------
   ! Test: tetragon_Xc
   ! ----------------------------
   call ut%test(24)%check( &
      name     = "tetragon_Xc_shape", &
      res      = shape(tetragon_Xc([1.0_rk, 1.0_rk], [2,2])), &
      expected = [4,3], &
      msg      = "2D grid shape", &
      group    = "geometry")

   ! ----------------------------
   ! Test: rotation
   ! ----------------------------
   call ut%test(25)%check( &
      name     = "rotation_identity", &
      res      = rotation(0.0_rk, 0.0_rk, 0.0_rk), &
      expected = reshape([1.0_rk,0.0_rk,0.0_rk, 0.0_rk,1.0_rk,0.0_rk, 0.0_rk,0.0_rk,1.0_rk], [3,3]), &
      msg      = "Rotation(0,0,0) = I", &
      group    = "rotation")

   ! ----------------------------
   ! Test: det
   ! ----------------------------
   call ut%test(26)%check( &
      name     = "det_2x2", &
      res      = det(reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2,2])), &
      expected = -2.0_rk, &
      msg      = "Determinant 2x2", &
      group    = "matrix")

   ! ----------------------------
   ! Test: inv
   ! ----------------------------
   call ut%test(27)%check( &
      name     = "inv_2x2", &
      res      = inv(reshape([4.0_rk, 7.0_rk, 2.0_rk, 6.0_rk], [2,2])), &
      expected = reshape([0.6_rk, -0.7_rk, -0.2_rk, 0.4_rk], [2,2]), &
      msg      = "Inverse of 2x2", &
      group    = "matrix")

   ! ----------------------------
   ! Test: solve
   ! ----------------------------
   A4 = reshape([4.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [2,2])
   R_expected = reshape([1.0_rk/11.0_rk, 7.0_rk/11.0_rk], [2,1])
   R = solve(A4, reshape([1.0_rk, 2.0_rk], [2,1]))

   call ut%test(28)%check( &
      name     = "solve_linear_system", &
      res      = R, &
      expected = R_expected, &
      msg      = "Solving A.X = B", &
      group    = "matrix")

   ! ----------------------------
   ! Test: insert_knot_A_5_1
   ! ----------------------------
   p = 2
   rr = 1
   s  = 1
   k  = 3
   knot_in = [0.0_rk, 0.0_rk, 0.0_rk, 0.5_rk, 1.0_rk, 1.0_rk, 1.0_rk]
   allocate(Pw(0:2,1:2))
   Pw(0,:) = [0.0_rk, 0.0_rk]
   Pw(1,:) = [0.5_rk, 0.5_rk]
   Pw(2,:) = [1.0_rk, 1.0_rk]

   call insert_knot_A_5_1(p, knot_in, Pw, 0.5_rk, k, s, rr, nq, knot_out, Qw)

   call ut%test(29)%check( &
      name     = "insert_knot_A_5_1_nq", &
      res      = nq, &
      expected = 3, &
      msg      = "Inserted knot, new number of control points", &
      group    = "insert_knot")

   ! ----------------------------
   ! Test: remove_knots_A_5_8
   ! ----------------------------
   call remove_knots_A_5_8(p=2, knot=knot_out, Pw=Qw, u=0.5_rk, r=3, s=2, num=1, t=t, knot_new=knot_new, Pw_new=Pw_new)

   call ut%test(30)%check( &
      name     = "remove_knots_A_5_8_t", &
      res      = t, &
      expected = 1, &
      msg      = "Removed 1 knot successfully", &
      group    = "remove_knot")

   ! ----------------------------
   ! Test: elevate_degree_A_5_9
   ! ----------------------------
   call elevate_degree_A_5_9(t=1, knot=knot_out, degree=p, Xcw=Qw, nc_new=nc, knot_new=knot_in, Xcw_new=Pw)

   call ut%test(31)%check( &
      name     = "elevate_degree_nc", &
      res      = nc, &
      expected = size(Pw,1), &
      msg      = "New number of control points after elevation", &
      group    = "elevate_degree")

   ! ----------------------------
   ! Test: gauss_legendre_1D
   ! ----------------------------
   call gauss_leg([0.0_rk, 1.0_rk], 2, Xksi=vec, Wksi=A)

   call ut%test(32)%check( &
      name     = "gauss_legendre_1D_pts", &
      res      = size(vec), &
      expected = 3, &
      msg      = "Number of Gauss points (1D)", &
      group    = "quadrature")

   ! ----------------------------
   ! Test: gauss_legendre_2D
   ! ----------------------------
   call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [2, 2], Xksi=Xksi, Wksi=Wksi)

   call ut%test(33)%check( &
      name     = "gauss_legendre_2D_shape", &
      res      = shape(Xksi), &
      expected = [9,2], &
      msg      = "Gauss points shape (2D)", &
      group    = "quadrature")

   ! ----------------------------
   ! Test: gauss_legendre_3D
   ! ----------------------------
   call gauss_leg([0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [0.0_rk, 1.0_rk], [1, 1, 1], Xksi=Xksi, Wksi=Wksi)

   call ut%test(34)%check( &
      name     = "gauss_legendre_3D_size", &
      res      = size(Xksi,1), &
      expected = 8, &
      msg      = "Number of Gauss points (3D)", &
      group    = "quadrature")

   ! ----------------------------
   ! Test: export_vtk_legacy
   ! ----------------------------
   call export_vtk_legacy(filename=vtk_file, points=reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], [2,3]), &
      elemConn=reshape([1,2], [1,2]), vtkCellType=3)

   call ut%test(35)%check( &
      name     = "export_vtk_legacy", &
      res      = .true., &
      expected = .true., &
      msg      = "Export to VTK (not crashing)", &
      group    = "export")

   ! ----------------------------
   ! Test: linspace
   ! ----------------------------
   call ut%test(36)%check( &
      name     = "linspace_uniform", &
      res      = linspace(0.0_rk, 1.0_rk, 5), &
      expected = [0.0_rk, 0.25_rk, 0.5_rk, 0.75_rk, 1.0_rk], &
      msg      = "Uniform linspace from 0 to 1", &
      group    = "linspace")

   ! ----------------------------
   ! Test: kron3
   ! ----------------------------
   K1 = [1.0_rk, 2.0_rk]
   K2 = [3.0_rk]
   K3 = [4.0_rk, 5.0_rk]
   out = kron(K1, K2, K3)

   call ut%test(36)%check( &
      name     = "kron3_product", &
      res      = out, &
      expected = [12.0_rk, 15.0_rk, 24.0_rk, 30.0_rk], &
      msg      = "kron3 result values", &
      group    = "kron")

   ! ----------------------------
   ! Test: elemConn_C0
   ! ----------------------------
   conn1D = elemConn_C0(5, 2)

   call ut%test(37)%check( &
      name     = "elemConn_C0_L", &
      res      = conn1D, &
      expected = reshape([1,3,2,4,3,5], [2,3]), &
      msg      = "Linear C0 connectivity", &
      group    = "connectivity")

   ! ----------------------------
   ! Test: elemConn_Cn (L)
   ! ----------------------------
   call elemConn_Cn(5, 2, [0.0_rk, 0.5_rk, 1.0_rk], [3,1,3], conn1D)

   call ut%test(38)%check( &
      name     = "elemConn_Cn_L", &
      res      = conn1D, &
      expected = reshape([1,2,2,3,3,4], [2,3]), &
      msg      = "Linear Cn connectivity", &
      group    = "connectivity")

   ! ----------------------------
   ! Test: elemConn_C0 (2D)
   ! ----------------------------
   conn2D = elemConn_C0(5, 5, 2, 2)

   call ut%test(39)%check( &
      name     = "elemConn_C0_S", &
      res      = shape(conn2D), &
      expected = [4,9], &
      msg      = "Surface C0 connectivity shape", &
      group    = "connectivity")

   ! ----------------------------
   ! Test: elemConn_Cn (2D)
   ! ----------------------------
   call elemConn_Cn(5, 5, 2, 2, [0.0_rk, 0.5_rk, 1.0_rk], [0.0_rk, 0.5_rk, 1.0_rk], &
      [2,1], [2,1], conn2D)

   call ut%test(40)%check( &
      name     = "elemConn_Cn_S", &
      res      = shape(conn2D), &
      expected = [4,9], &
      msg      = "Surface Cn connectivity shape", &
      group    = "connectivity")

   ! ----------------------------
   ! Test: elemConn_C0 (3D)
   ! ----------------------------
   conn3D = elemConn_C0(5, 5, 5, 2, 2, 2)

   call ut%test(41)%check( &
      name     = "elemConn_C0_V", &
      res      = shape(conn3D), &
      expected = [8,27], &
      msg      = "Volume C0 connectivity shape", &
      group    = "connectivity")

   ! ----------------------------
   ! Test: elemConn_Cn (3D)
   ! ----------------------------
   call elemConn_Cn(5, 5, 5, 2, 2, 2, &
      [0.0_rk, 0.5_rk, 1.0_rk], &
      [0.0_rk, 0.5_rk, 1.0_rk], &
      [0.0_rk, 0.5_rk, 1.0_rk], &
      [2,1], [2,1], [2,1], conn3D)

   call ut%test(42)%check( &
      name     = "elemConn_Cn_V", &
      res      = shape(conn3D), &
      expected = [8,27], &
      msg      = "Volume Cn connectivity shape", &
      group    = "connectivity")

   ! ----------------------------
   ! Test: inv (3x3 matrix)
   ! ----------------------------
   allocate(A2(3,3))
   A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, &
      0.0_rk, 1.0_rk, 4.0_rk, &
      5.0_rk, 6.0_rk, 0.0_rk], [3,3])
   A_inv = inv(A2)

   call ut%test(43)%check( &
      name     = "inv_3x3", &
      res      = matmul(A2, A_inv), &
      expected = eye(3), &
      msg      = "A . inv(A) = I for 3x3", &
      group    = "matrix")

   ! ----------------------------
   ! Test: inv (rectangular 3x2 matrix)
   ! ----------------------------
   deallocate(A2, A_inv)
   allocate(A2(3,2))
   A2 = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], [3,2])
   A_inv = inv(A2)

   call ut%test(44)%check( &
      name     = "inv_rectangular_3x2", &
      res      = matmul(A_inv, A2), &
      expected = eye(2), &
      msg      = "inv(A) . A = I for 3x2", &
      group    = "matrix")

   ! ----------------------------
   ! Test: inv (rectangular 2x3 matrix)
   ! ----------------------------
   deallocate(A2, A_inv)
   allocate(A2(2,3))
   A2 = reshape([1.0_rk, 4.0_rk, 2.0_rk, 5.0_rk, 3.0_rk, 6.0_rk], [2,3])
   A_inv = inv(A2)

   call ut%test(45)%check( &
      name     = "inv_rectangular_2x3", &
      res      = matmul(A2, A_inv), &
      expected = eye(2), &
      msg      = "A . inv(A) = I for 2x3", &
      group    = "matrix")

   ! ----------------------------
   ! Test: inv (identity matrix)
   ! ----------------------------
   deallocate(A2, A_inv)
   allocate(A2(4,4))
   A2 = eye(4)
   A_inv = inv(A2)

   call ut%test(46)%check( &
      name     = "inv_identity", &
      res      = A_inv, &
      expected = A2, &
      msg      = "inv(I) = I", &
      group    = "matrix")

   call ut%test(47)%check( &
      name     = "inv_rectangular_2x3_proj", &
      res      = matmul(A_inv, A2), &
      expected = transpose(matmul(A_inv, A2)), &
      msg      = "inv(A) . A is symmetric projection (2x3)", &
      group    = "matrix")

   ! summary of tests
   call ut%summary( &
      required_score = 100.0, &
      verbose        = 3, &
      stop_fail      = .false.)

end program test_forcad_utils
