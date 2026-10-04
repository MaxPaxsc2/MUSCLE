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
    !! c31, c32 on the normal components and c44 (yz), c55 (zx), c66 (xy) on the shears.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_base
    implicit none
    private

    public :: Yld2004
    type, extends(Base_yield_critera) :: Yld2004
        !! Yld2004-18p yield criterion
        real(real64) :: Ln(3,3,2)
            !! Normal block of L' = C' P_dev (k = 1) and L'' = C'' P_dev (k = 2)
        real(real64) :: Ls(3,2)
            !! Shear coefficients of C' and C'' in MUSCLE order (xy, yz, xz) = (c66, c44, c55);
            !! the paper orders the shears yz, zx, xy (Barlat et al., 2005, Eq. A.15)
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
        !! that each call applies them to the stress without forming the deviator. Arguments use
        !! the notation of the paper: cpIJ = c'_IJ, cppIJ = c''_IJ; 44, 55, 66 are the yz, zx, xy
        !! shears.
        implicit none
        class(Yld2004), intent(inout) :: self
        real(real64), intent(in) :: cp12, cp13, cp21, cp23, cp31, cp32, cp44, cp55, cp66
        real(real64), intent(in) :: cpp12, cpp13, cpp21, cpp23, cpp31, cpp32, cpp44, cpp55, cpp66
        real(real64), intent(in) :: a

        real(real64) :: C(3,3,2), P(3,3)
        integer :: k

        ! Normal blocks of C' and C'' (zero diagonal, Barlat et al., 2005, Eq. 15)
        C(:,:,1) = reshape([0D0, -cp21, -cp31, -cp12, 0D0, -cp32, -cp13, -cp23, 0D0], [3,3])
        C(:,:,2) = reshape([0D0, -cpp21, -cpp31, -cpp12, 0D0, -cpp32, -cpp13, -cpp23, 0D0], [3,3])
        P = -1D0/3D0
        do k = 1, 3
            P(k,k) = 2D0/3D0
        end do
        do k = 1, 2
            self%Ln(:,:,k) = matmul(C(:,:,k), P)
        end do
        self%Ls(:,1) = [cp66, cp44, cp55]
        self%Ls(:,2) = [cpp66, cpp44, cpp55]
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

        lp = eigenvals(transformed_stress(self, stress, 1))   ! lam'
        lpp = eigenvals(transformed_stress(self, stress, 2))  ! lam''
        res = (0.25D0*sum(abs(spread(lp, 2, 3) - spread(lpp, 1, 3))**self%a))**(1D0/self%a)
    end function stress_eq

    pure function transformed_stress(self, stress, k) result(res)
        ! S' = L':stress (k = 1) or S'' = L'':stress (k = 2) (Barlat et al., 2005, Eq. 11):
        ! 9 products on the normal components and 3 on the tensorial shears. Called by stress_eq.
        implicit none
        class(Yld2004), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
        integer, intent(in) :: k
        type(ten_3D2Osym) :: res

        call res%init([matmul(self%Ln(:,:,k), stress%vals(1:3)), self%Ls(:,k)*stress%vals(4:6)])
    end function transformed_stress

end module muscle_yield_yld2004
