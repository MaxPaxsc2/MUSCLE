! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_yield_yld2004
    !! Module muscle_yield_yld2004
    !! ===========================
    !! Yld2004-18p anisotropic yield criterion (Barlat et al., 2005).
    !!
    !! Two linear transformations of the deviator s, S' = C':s and S'' = C'':s, with
    !! principal values lam'_i and lam''_j, define
    !! \[ \phi = \sum_i \sum_j |\lambda'_i - \lambda''_j|^a = 4 \bar{\sigma}^a \]
    !! (Barlat et al., 2005, Eq. 14) and the equivalent stress \( \bar{\sigma} = (\phi/4)^{1/a} \).
    !! C' and C'' take the 9 coefficients each of Barlat et al. (2005, Eq. 15): c12, c13, c21, c23,
    !! c31, c32 on the normal components and c44 (xy), c55 (yz), c66 (xz) on the shears, in the
    !! library order; the paper orders the shears yz, zx, xy (Eq. A.15), so
    !! (c44, c55, c66) here = (c66, c44, c55) in the paper.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_base
    use muscle_yield_linear_transform, only : deviatoric_block, linear_transform
    implicit none
    private

    public :: Yld2004
    type, extends(Base_yield_critera) :: Yld2004
        !! Yld2004-18p yield criterion
        real(real64) :: Lp(3,3)
            !! Normal block of L' = C' P_dev, which maps (sxx, syy, szz) to (S'xx, S'yy, S'zz)
        real(real64) :: cp_shear(3)
            !! Shear factors of C' on (xy, yz, xz): (c'44, c'55, c'66)
        real(real64) :: Lpp(3,3)
            !! Normal block of L'' = C'' P_dev
        real(real64) :: cpp_shear(3)
            !! Shear factors of C'' on (xy, yz, xz): (c''44, c''55, c''66)
        real(real64) :: a
            !! Exponent a
    contains
        procedure :: init
        procedure :: stress_eq
    end type Yld2004

contains

    subroutine init(self, cp12, cp13, cp21, cp23, cp31, cp32, cp44, cp55, cp66, &
                    cpp12, cpp13, cpp21, cpp23, cpp31, cpp32, cpp44, cpp55, cpp66, a)
        !! Stores L' = C' P_dev and L'' = C'' P_dev (Barlat et al., 2005, Eqs. 11 and 15), so
        !! that each call applies them to the stress without forming the deviator. cpIJ = c'_IJ
        !! and cppIJ = c''_IJ; the shears follow the library order: 44 -> xy, 55 -> yz, 66 -> xz.
        implicit none
        class(Yld2004), intent(inout) :: self
        real(real64), intent(in) :: cp12, cp13, cp21, cp23, cp31, cp32, cp44, cp55, cp66
        real(real64), intent(in) :: cpp12, cpp13, cpp21, cpp23, cpp31, cpp32, cpp44, cpp55, cpp66
        real(real64), intent(in) :: a

        self%Lp = deviatoric_block(normal_block(cp12, cp13, cp21, cp23, cp31, cp32))
        self%cp_shear(1) = cp44
        self%cp_shear(2) = cp55
        self%cp_shear(3) = cp66
        self%Lpp = deviatoric_block(normal_block(cpp12, cpp13, cpp21, cpp23, cpp31, cpp32))
        self%cpp_shear(1) = cpp44
        self%cpp_shear(2) = cpp55
        self%cpp_shear(3) = cpp66
        self%a = a
    end subroutine init

    pure function stress_eq(self, stress) result(res)
        !! Equivalent stress (phi/4)**(1/a) (Barlat et al., 2005, Eq. 14), with the
        !! eigenvalues of S' and S'' from eigenvals.
        use muscle_math_operations, only : eigenvals
        implicit none
        class(Yld2004), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
        real(real64) :: res

        real(real64) :: lp(3), lpp(3)

        lp = eigenvals(linear_transform(self%Lp, self%cp_shear, stress))     ! lam'
        lpp = eigenvals(linear_transform(self%Lpp, self%cpp_shear, stress))  ! lam''
        res = (0.25D0*(abs(lp(1) - lpp(1))**self%a + abs(lp(1) - lpp(2))**self%a &
                       + abs(lp(1) - lpp(3))**self%a                             &
                       + abs(lp(2) - lpp(1))**self%a + abs(lp(2) - lpp(2))**self%a &
                       + abs(lp(2) - lpp(3))**self%a                             &
                       + abs(lp(3) - lpp(1))**self%a + abs(lp(3) - lpp(2))**self%a &
                       + abs(lp(3) - lpp(3))**self%a))**(1D0/self%a)
    end function stress_eq

    pure function normal_block(c12, c13, c21, c23, c31, c32) result(C)
        ! Normal block of C' or C'' (Barlat et al., 2005, Eq. 15): zero diagonal and -cIJ off it.
        ! Called by init.
        implicit none
        real(real64), intent(in) :: c12, c13, c21, c23, c31, c32
        real(real64) :: C(3,3)

        C(1,1) = 0D0
        C(1,2) = -c12
        C(1,3) = -c13
        C(2,1) = -c21
        C(2,2) = 0D0
        C(2,3) = -c23
        C(3,1) = -c31
        C(3,2) = -c32
        C(3,3) = 0D0
    end function normal_block

end module muscle_yield_yld2004
