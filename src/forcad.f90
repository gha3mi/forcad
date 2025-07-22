!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
module forcad

    use forcad_kinds, only: rk

    use forcad_nurbs_curve, only: nurbs_curve
    use forcad_nurbs_surface, only: nurbs_surface
    use forcad_nurbs_volume, only: nurbs_volume

    implicit none

    private
    public rk, nurbs_curve, nurbs_surface, nurbs_volume

end module forcad
