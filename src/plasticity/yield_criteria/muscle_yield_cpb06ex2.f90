! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_yield_cpb06ex2
    !! CPB06ex2: two transformations of the stress deviator (Plunkett et al., 2008, Eq. 8).
    !! Library (c44, c55, c66) on (xy, yz, xz) = paper (C66, C44, C55): xy from Eq. (13),
    !! yz and xz from the out-of-plane values of Tables 9 and 11.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_base
    implicit none
    private

    public :: CPB06ex2
    type, extends(Base_yield_critera) :: CPB06ex2
        real(real64) :: L(3,3)
            !! Normal block of C P_dev, which maps the normal stresses to those of Sigma (Eq. 9)
        real(real64) :: Lp(3,3)
            !! Normal block of Cp P_dev, for Sigmap (Eq. 9; suffix p = primed in the paper)
        real(real64) :: c_shear(3)
            !! Shear factors of C on the tensorial shears (xy, yz, xz)
        real(real64) :: cp_shear(3)
            !! Shear factors of Cp on the tensorial shears (xy, yz, xz), for Sigmap
        real(real64) :: k
            !! Strength-differential parameter k of Sigma (Eq. 8)
        real(real64) :: kp
            !! Strength-differential parameter kp of Sigmap (Eq. 8)
        real(real64) :: a
            !! Degree of homogeneity a (Eq. 8)
        real(real64) :: B
            !! Normalization: f = s in uniaxial tension s along RD (Eq. 11 with theta = 0)
    contains
        procedure :: init
        procedure :: stress_eq
    end type CPB06ex2

contains

    subroutine init(self, c11, c12, c13, c21, c22, c23, c31, c32, c33, c44, c55, c66, k, &
                    cp11, cp12, cp13, cp21, cp22, cp23, cp31, cp32, cp33, cp44, cp55, cp66, kp, a)
        !! Stores the maps of Plunkett et al. (2008), Eqs. (9)-(10). Normalizes Eq. (8)
        !! to unit uniaxial tension in RD, as in Eq. (11). Called when setting the material.
        implicit none
        class(CPB06ex2), intent(inout) :: self
        real(real64), intent(in) :: c11, c12, c13, c21, c22, c23, c31, c32, c33, c44, c55, c66, k
        real(real64), intent(in) :: cp11, cp12, cp13, cp21, cp22, cp23, cp31, cp32, cp33, &
                                    cp44, cp55, cp66, kp, a
        type(ten_3D2Osym) :: rd

        self%L = deviatoric_block(c11, c12, c13, c21, c22, c23, c31, c32, c33)
        self%Lp = deviatoric_block(cp11, cp12, cp13, cp21, cp22, cp23, cp31, cp32, cp33)
        self%c_shear = [c44, c55, c66]
        self%cp_shear = [cp44, cp55, cp66]
        self%k = k
        self%kp = kp
        self%a = a
        call rd%init(xx=1D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
        self%B = 1D0
        self%B = 1D0/self%stress_eq(rd)
    end subroutine init

    pure function stress_eq(self, stress) result(res)
        !! Equivalent stress B F**(1/a), Plunkett et al. (2008), Eq. (8).
        !! Called by the return mapping; principal values come from eigenvals.
        use muscle_math_operations, only : eigenvals
        implicit none
        class(CPB06ex2), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
        real(real64) :: res
        real(real64) :: lam(3), lamp(3)

        lam = eigenvals(linear_transform(self%L, self%c_shear, stress))
        lamp = eigenvals(linear_transform(self%Lp, self%cp_shear, stress))
        res = self%B*((abs(lam(1)) - self%k*lam(1))**self%a &
                      + (abs(lam(2)) - self%k*lam(2))**self%a &
                      + (abs(lam(3)) - self%k*lam(3))**self%a &
                      + (abs(lamp(1)) - self%kp*lamp(1))**self%a &
                      + (abs(lamp(2)) - self%kp*lamp(2))**self%a &
                      + (abs(lamp(3)) - self%kp*lamp(3))**self%a)**(1D0/self%a)
    end function stress_eq

    pure function deviatoric_block(c11, c12, c13, c21, c22, c23, c31, c32, c33) result(L)
        ! Normal block L = C P_dev (Plunkett et al., 2008, Eq. 9). Called by init.
        ! The third column comes from the zero row sum, also in floating point.
        implicit none
        real(real64), intent(in) :: c11, c12, c13, c21, c22, c23, c31, c32, c33
        real(real64) :: L(3,3), m1, m2, m3

        m1 = (c11 + c12 + c13)/3D0
        m2 = (c21 + c22 + c23)/3D0
        m3 = (c31 + c32 + c33)/3D0
        L(1,1) = c11 - m1
        L(1,2) = c12 - m1
        L(1,3) = -L(1,1) - L(1,2)
        L(2,1) = c21 - m2
        L(2,2) = c22 - m2
        L(2,3) = -L(2,1) - L(2,2)
        L(3,1) = c31 - m3
        L(3,2) = c32 - m3
        L(3,3) = -L(3,1) - L(3,2)
    end function deviatoric_block

    pure function linear_transform(L, c_shear, stress) result(res)
        ! Sigma = L : stress (Plunkett et al., 2008, Eq. 9), for either transformation.
        ! Called by stress_eq. Zero row sums cancel pressure before the products.
        implicit none
        real(real64), intent(in) :: L(3,3), c_shear(3)
        type(ten_3D2Osym), intent(in) :: stress
        type(ten_3D2Osym) :: res
        real(real64) :: dxx, dyy

        dxx = stress%vals(1) - stress%vals(3)
        dyy = stress%vals(2) - stress%vals(3)
        res%vals(1) = L(1,1)*dxx + L(1,2)*dyy
        res%vals(2) = L(2,1)*dxx + L(2,2)*dyy
        res%vals(3) = L(3,1)*dxx + L(3,2)*dyy
        res%vals(4) = c_shear(1)*stress%vals(4)
        res%vals(5) = c_shear(2)*stress%vals(5)
        res%vals(6) = c_shear(3)*stress%vals(6)
    end function linear_transform

end module muscle_yield_cpb06ex2
