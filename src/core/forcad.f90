! author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Public ForCAD API facade.
!!
!! Applications should normally import symbols from this module rather than
!! implementation modules:
!!
!! ```fortran
!! use forcad, only : rk, nurbs_curve, nurbs_surface, nurbs_volume
!! ```
!!
!! Geometry coordinates are stored as `X(i,d)`, where the first subscript is
!! the flattened control-, sample-, or element-local index and the second is
!! the physical coordinate. Tensor-product indices are flattened with the
!! first parametric direction varying fastest. Public real data use `rk`.
!!
!! **Mathematical conventions**
!!
!! A rank-\(r\) geometry maps a parameter vector
!! \(\boldsymbol{\xi}=(\xi_1,\ldots,\xi_r)\) to physical coordinates by
!!
!! \[
!! \mathbf{x}(\boldsymbol{\xi})
!!   =\sum_A R_A(\boldsymbol{\xi})\mathbf{P}_A,\qquad
!! R_A=\frac{N_Aw_A}{\sum_B N_Bw_B}.
!! \]
!!
!! Here \(N_A\) is a tensor-product B-spline basis, \(\mathbf{P}_A\) is a
!! control point, and \(w_A>0\) is its weight. Equal weights reduce \(R_A\) to
!! the polynomial B-spline basis. Curve, surface, and volume ranks are one,
!! two, and three respectively. Degrees are polynomial degrees, not orders;
!! spline order in direction \(d\) is `degree(d)+1`.
!!
!! Knot vectors are nondecreasing and use the complete NURBS representation,
!! including repeated end knots and any knots outside the active interval.
!! Parameter wrapping is an evaluation policy and is distinct from geometric
!! periodicity of the control data.
!!
!! Operations use the public structured diagnostic object `err`.
!! An operation called while `err%ok` is false normally leaves the object
!! unchanged; handle or reset that diagnostic before continuing that workflow.
!! `ansatz` and `ansatz_element` also accept a caller-owned `err` of the same
!! diagnostic type. Use it with private output arrays during parallel assembly
!! to avoid writing shared diagnostics.
!!
!! `modify_Xc` and `modify_Wc` accept scalar edits or complete arrays at fixed
!! spline topology. They invalidate sampled points while retaining parameters;
!! `create()` refreshes that grid. Invalid updates and sampling requests leave
!! prior data intact. Other setters, fitting, and refinement retain their
!! documented cache rules; rebuild connectivity after spline-space changes.
!! Dense basis queries preserve samples and reuse cached grid or paired rows;
!! `cmp_elem` computes IGA connectivity independently of visualization.
!! `cmp_Xg` evaluates one parameter point or a batch without modifying caches;
!! batch results contain one physical-coordinate row per parameter point.
!!
!! The facade exports:
!!
!! - single-patch curve, surface, and volume types;
!! - exact extrusion, revolution, and fixed-frame translational sweep, plus
!!   homogeneous global loft interpolation;
!! - curve, surface, and volume multipatch containers;
!! - selected grid, primitive-control-net, transformation, and dense-solve
!!   utilities.
!!
!! Implementation modules expose additional low-level kernels. Applications
!! should prefer this facade unless they intentionally need those contracts.
module forcad

    use forcad_kinds, only: rk
    use forcad_nurbs_curve, only: nurbs_curve
    use forcad_nurbs_surface, only: nurbs_surface
    use forcad_nurbs_volume, only: nurbs_volume
    use forcad_nurbs_construction, only: extrude, revolve, sweep, loft
    use forcad_multipatch, only: multipatch_connection
    use forcad_multipatch_curve, only: nurbs_multipatch_curve
    use forcad_multipatch_surface, only: nurbs_multipatch_surface
    use forcad_multipatch_volume, only: nurbs_multipatch_volume
    use forcad_utils, only: ndgrid, hexahedron_Xc, rotation, solve

    implicit none

    private
    public rk, nurbs_curve, nurbs_surface, nurbs_volume
    public extrude, revolve, sweep, loft
    public multipatch_connection, nurbs_multipatch_curve, nurbs_multipatch_surface, nurbs_multipatch_volume
    public ndgrid, hexahedron_Xc, rotation, solve

end module forcad
