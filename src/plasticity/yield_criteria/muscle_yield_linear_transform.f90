! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_yield_linear_transform
    !! Module muscle_yield_linear_transform
    !! ====================================
    !! Orthotropic linear transformation S = L : sigma of the anisotropic yield criteria based on
    !! linear transformations of the stress deviator (CPB06, Yld2004-18p), with L = C P_dev:
    !! a 3x3 block on the normal components and one factor per tensorial shear. The criterion
    !! stores L in init, so no deviator is formed per call. The module gives L, S and the
    !! adjoint L^T : g that takes a gradient with respect to S to one with respect to sigma
    !! (Barlat et al., 2005, Eq. 11; Cazacu et al., 2006, Eq. 8).
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    implicit none
    private

    public :: deviatoric_block
    public :: linear_transform
    public :: pull_back

contains

    pure function deviatoric_block(C) result(L)
        !! Normal block of L = C P_dev, so that L : sigma = C : dev(sigma) (Barlat et al., 2005,
        !! Eq. 11; Cazacu et al., 2006, Eq. 8), with C the 3x3 normal block of the transformation
        !! and P_dev the deviatoric projection: L_ij = C_ij - (1/3) sum_k C_ik. Each row of L sums
        !! to zero (L annihilates a hydrostatic stress). Called by the init of CPB06 and
        !! Yld2004-18p.
        real(real64), intent(in) :: C(3,3)  !! Normal block of C (general, not only symmetric)
        real(real64) :: L(3,3)

        real(real64) :: m1, m2, m3

        ! m_i = (1/3) sum_k C_ik. The third column is taken from the row sum, so that each row
        ! sums to zero also in floating point
        m1 = (C(1,1) + C(1,2) + C(1,3))/3D0
        m2 = (C(2,1) + C(2,2) + C(2,3))/3D0
        m3 = (C(3,1) + C(3,2) + C(3,3))/3D0
        L(1,1) = C(1,1) - m1
        L(1,2) = C(1,2) - m1
        L(1,3) = -L(1,1) - L(1,2)
        L(2,1) = C(2,1) - m2
        L(2,2) = C(2,2) - m2
        L(2,3) = -L(2,1) - L(2,2)
        L(3,1) = C(3,1) - m3
        L(3,2) = C(3,2) - m3
        L(3,3) = -L(3,1) - L(3,2)
    end function deviatoric_block

    pure function linear_transform(L, c_shear, stress) result(res)
        !! S = L : stress = C : dev(stress) (Barlat et al., 2005, Eq. 11; Cazacu et al., 2006,
        !! Eq. 8): the 3x3 block L acts on the normal components and each factor of c_shear on
        !! its tensorial shear (xy, yz, xz). Called by the Yld2004-18p equivalent stress and by
        !! the CPB06 gradient.
        real(real64), intent(in) :: L(3,3)            !! Normal block, from deviatoric_block
        real(real64), intent(in) :: c_shear(3)        !! Shear factors on (xy, yz, xz)
        type(ten_3D2Osym), intent(in) :: stress
        type(ten_3D2Osym) :: res

        real(real64) :: dxx, dyy

        ! The rows of L sum to zero, so L acts on the differences of the normal stresses: the
        ! pressure cancels before the products and not by rounding after them
        dxx = stress%vals(1) - stress%vals(3)
        dyy = stress%vals(2) - stress%vals(3)
        res%vals(1) = L(1,1)*dxx + L(1,2)*dyy
        res%vals(2) = L(2,1)*dxx + L(2,2)*dyy
        res%vals(3) = L(3,1)*dxx + L(3,2)*dyy
        res%vals(4) = c_shear(1)*stress%vals(4)
        res%vals(5) = c_shear(2)*stress%vals(5)
        res%vals(6) = c_shear(3)*stress%vals(6)
    end function linear_transform

    pure function pull_back(L, c_shear, x) result(res)
        !! L^T : x, the adjoint of linear_transform: by the chain rule through S = L : stress
        !! (Barlat et al., 2005, Eq. 11; Cazacu et al., 2006, Eq. 8) it takes a derivative with
        !! respect to S to one with respect to the stress. Called by the CPB06 gradient.
        real(real64), intent(in) :: L(3,3)            !! Normal block, from deviatoric_block
        real(real64), intent(in) :: c_shear(3)        !! Shear factors on (xy, yz, xz)
        type(ten_3D2Osym), intent(in) :: x
        type(ten_3D2Osym) :: res

        real(real64) :: n(3)

        n = normal_pull_back(L, x%vals(1), x%vals(2), x%vals(3))
        call res%init(xx=n(1), yy=n(2), zz=n(3), xy=c_shear(1)*x%vals(4), &
                      yz=c_shear(2)*x%vals(5), xz=c_shear(3)*x%vals(6))
    end function pull_back

    pure function normal_pull_back(L, x, y, z) result(res)
        ! L^T (x, y, z) on the normal components, from the chain rule through S = L : stress
        ! (Barlat et al., 2005, Eq. 11; Cazacu et al., 2006, Eq. 8); used by pull_back
        real(real64), intent(in) :: L(3,3)
        real(real64), intent(in) :: x, y, z
        real(real64) :: res(3)

        res(1) = L(1,1)*x + L(2,1)*y + L(3,1)*z
        res(2) = L(1,2)*x + L(2,2)*y + L(3,2)*z
        res(3) = L(1,3)*x + L(2,3)*y + L(3,3)*z
    end function normal_pull_back

end module muscle_yield_linear_transform
