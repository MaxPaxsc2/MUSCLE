! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_math_inverses
    !! Small dense matrix inverses used by the tensor engine.
    !!
    !! Algorithm by size:
    !!
    !! - 3x3 (`M33INV`) and 4x4 (`M44INV`): closed-form adjugate (cofactors / determinant),
    !!   with an LU fallback when the condition number
    !!   \( \kappa_\infty \) exceeds `ADJ_KAPPA_MAX` (the adjugate loses accuracy much faster than
    !!   LU as the matrix becomes ill-conditioned).
    !! - 5x5, 6x6 and 9x9 (`M55INV`, `M66INV`, `M99INV`): LU with partial pivoting in plain
    !!   Fortran.
    !!
    !! A singular matrix (exactly zero determinant or pivot) gives `OK_FLAG = .false.` and a zero
    !! inverse.
    use, intrinsic :: iso_fortran_env
    
    public
    private :: lu_inv3, lu_inv4

    real(real64), parameter :: ADJ_KAPPA_MAX = 1.0D3
        !! Largest \( \kappa_\infty \) accepted for the adjugate inverses (`M33INV`, `M44INV`); above
        !! it they are recomputed by LU. With 1e3 the result is as accurate as LU for every condition
        !! number, while well-conditioned matrices keep the fast adjugate path.
  
  
    contains

    PURE SUBROUTINE M66INV(A, AINV, OK_FLAG)
        !! Inverse of a general 6x6 matrix by LU with partial pivoting in plain Fortran,
        !! see `muscle_math_lu_inverse.inc`. Used by `ten_3D4O3sym` and
        !! `ten_3D4O2sym`.
        !!
        implicit none
        integer, parameter :: N = 6
        REAL(real64), INTENT(IN)  :: A(N,N)
            !! Matrix to be inverted.
        REAL(real64), INTENT(OUT) :: AINV(N,N)
            !! Inverse of `A`; zero when `A` is singular.
        LOGICAL, INTENT(OUT)      :: OK_FLAG
            !! `.TRUE.` if the inversion succeeded, `.FALSE.` if `A` is singular.
        include 'muscle_math_lu_inverse.inc'
    END SUBROUTINE M66INV
    
    PURE SUBROUTINE M99INV(A, AINV, OK_FLAG)
        !! Inverse of a general 9x9 matrix by LU with partial pivoting in plain Fortran (algorithm of
        !! LAPACK `DGETRF` + `DGETRI`, see `muscle_math_lu_inverse.inc`). Used by `ten_3D4O`.
        implicit none
        integer, parameter :: N = 9
        REAL(real64), INTENT(IN)  :: A(N,N)
            !! Matrix to be inverted.
        REAL(real64), INTENT(OUT) :: AINV(N,N)
            !! Inverse of `A`; zero when `A` is singular.
        LOGICAL, INTENT(OUT)      :: OK_FLAG
            !! `.TRUE.` if the inversion succeeded, `.FALSE.` if `A` is singular.
        include 'muscle_math_lu_inverse.inc'
    END SUBROUTINE M99INV

    PURE FUNCTION QUIET_NAN() RESULT(X)
        !! IEEE quiet NaN (`ieee_value`), used to fill the inverse of a singular tensor.
        !!
        !! It is a separate function on purpose and must only be called on the (rare) failure path:
        !! gfortran saves and restores the floating-point state on entry and exit of every procedure
        !! that uses the IEEE intrinsic modules (~150 ns per call, measured 2026-10-14). Keeping the
        !! `use ieee_arithmetic` inside this function confines that cost to the singular case; using
        !! `ieee_value` directly in the inverse routines made every inverse ~150 ns slower.
        use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
        implicit none
        REAL(real64) :: X
            !! Quiet NaN.
        X = ieee_value(1.0D0, ieee_quiet_nan)
    END FUNCTION QUIET_NAN

    PURE FUNCTION COND_INF(A, AINV) RESULT(KAPPA)
        !! Condition number in the infinity norm,
        !! \( \kappa_\infty = \|A\|_\infty \|A^{-1}\|_\infty \) (maximum absolute row sums), from a
        !! matrix and its computed inverse. Costs \( 2n^2 \) operations.
        implicit none
        REAL(real64), INTENT(IN) :: A(:,:)
            !! Matrix.
        REAL(real64), INTENT(IN) :: AINV(:,:)
            !! Its (computed) inverse, same shape as `A`.
        REAL(real64) :: KAPPA
            !! \( \kappa_\infty(A) \).
        KAPPA = maxval(sum(abs(A), dim=2)) * maxval(sum(abs(AINV), dim=2))
    END FUNCTION COND_INF

    PURE SUBROUTINE M33INV(A, AINV, OK_FLAG, KAPPA)
        !! Inverse of a general 3x3 matrix: closed-form adjugate with an LU fallback.
        !!
        !! \( A^{-1} = \mathrm{adj}(A)/\det(A) \). The condition number
        !! \( \kappa_\infty = \|A\|_\infty \|A^{-1}\|_\infty \) is computed from the result (18 extra
        !! operations). If \( \kappa_\infty > \) `ADJ_KAPPA_MAX`, or the determinant is zero or not a
        !! finite number (e.g. a well-conditioned but extremely scaled matrix), the inverse is
        !! recomputed by LU with partial pivoting.
        !!
        implicit none
        REAL(real64), INTENT(IN)  :: A(3,3)
            !! Matrix to be inverted.
        REAL(real64), INTENT(OUT) :: AINV(3,3)
            !! Inverse of `A`; zero when `A` is singular.
        LOGICAL, INTENT(OUT)      :: OK_FLAG
            !! `.TRUE.` if the inversion succeeded, `.FALSE.` if `A` is singular.
        REAL(real64), INTENT(OUT), OPTIONAL :: KAPPA
            !! Condition number \( \kappa_\infty(A) \) of `A`; `huge(1.0D0)` when `A` is singular.
        REAL(real64) :: c11, c12, c13, c21, c22, c23, c31, c32, c33, det, r, kap

        ! Cofactors (AINV = transpose of the cofactor matrix / det)
        c11 = A(2,2)*A(3,3) - A(2,3)*A(3,2)
        c12 = A(1,3)*A(3,2) - A(1,2)*A(3,3)
        c13 = A(1,2)*A(2,3) - A(1,3)*A(2,2)
        c21 = A(2,3)*A(3,1) - A(2,1)*A(3,3)
        c22 = A(1,1)*A(3,3) - A(1,3)*A(3,1)
        c23 = A(1,3)*A(2,1) - A(1,1)*A(2,3)
        c31 = A(2,1)*A(3,2) - A(2,2)*A(3,1)
        c32 = A(1,2)*A(3,1) - A(1,1)*A(3,2)
        c33 = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        det = A(1,1)*c11 + A(1,2)*c21 + A(1,3)*c31

        if (abs(det) > 0.0D0 .and. abs(det) <= huge(det)) then
            r = 1.0D0 / det
            AINV(1,1) = c11*r; AINV(1,2) = c12*r; AINV(1,3) = c13*r
            AINV(2,1) = c21*r; AINV(2,2) = c22*r; AINV(2,3) = c23*r
            AINV(3,1) = c31*r; AINV(3,2) = c32*r; AINV(3,3) = c33*r
            kap = maxval(sum(abs(A), dim=2)) * maxval(sum(abs(AINV), dim=2))
            if (kap <= ADJ_KAPPA_MAX) then      ! also .false. for NaN
                OK_FLAG = .true.
                if (present(KAPPA)) KAPPA = kap
                return
            end if
        end if

        call lu_inv3(A, AINV, OK_FLAG)
        if (present(KAPPA)) then
            if (OK_FLAG) then
                KAPPA = maxval(sum(abs(A), dim=2)) * maxval(sum(abs(AINV), dim=2))
            else
                KAPPA = huge(1.0D0)
            end if
        end if
    END SUBROUTINE M33INV

    PURE SUBROUTINE M44INV(A, AINV, OK_FLAG, KAPPA)
        !! Inverse of a general 4x4 matrix: closed-form adjugate (via 2x2 sub-determinants) with an
        !! LU fallback, same rule as `M33INV`. Used by the 2D tensor engine (4x4 Voigt matrices of
        !! `ten_2D4O3sym` / `ten_2D4O2sym`).
        !!
        implicit none
        REAL(real64), INTENT(IN)  :: A(4,4)
            !! Matrix to be inverted.
        REAL(real64), INTENT(OUT) :: AINV(4,4)
            !! Inverse of `A`; zero when `A` is singular.
        LOGICAL, INTENT(OUT)      :: OK_FLAG
            !! `.TRUE.` if the inversion succeeded, `.FALSE.` if `A` is singular.
        REAL(real64), INTENT(OUT), OPTIONAL :: KAPPA
            !! Condition number \( \kappa_\infty(A) \) of `A`; `huge(1.0D0)` when `A` is singular.
        REAL(real64) :: s0, s1, s2, s3, s4, s5, c0, c1, c2, c3, c4, c5, det, r, kap

        ! 2x2 sub-determinants of rows (1,2) and of rows (3,4)
        s0 = A(1,1)*A(2,2) - A(2,1)*A(1,2)
        s1 = A(1,1)*A(2,3) - A(2,1)*A(1,3)
        s2 = A(1,1)*A(2,4) - A(2,1)*A(1,4)
        s3 = A(1,2)*A(2,3) - A(2,2)*A(1,3)
        s4 = A(1,2)*A(2,4) - A(2,2)*A(1,4)
        s5 = A(1,3)*A(2,4) - A(2,3)*A(1,4)
        c5 = A(3,3)*A(4,4) - A(4,3)*A(3,4)
        c4 = A(3,2)*A(4,4) - A(4,2)*A(3,4)
        c3 = A(3,2)*A(4,3) - A(4,2)*A(3,3)
        c2 = A(3,1)*A(4,4) - A(4,1)*A(3,4)
        c1 = A(3,1)*A(4,3) - A(4,1)*A(3,3)
        c0 = A(3,1)*A(4,2) - A(4,1)*A(3,2)
        det = s0*c5 - s1*c4 + s2*c3 + s3*c2 - s4*c1 + s5*c0

        if (abs(det) > 0.0D0 .and. abs(det) <= huge(det)) then
            r = 1.0D0 / det
            AINV(1,1) = ( A(2,2)*c5 - A(2,3)*c4 + A(2,4)*c3) * r
            AINV(1,2) = (-A(1,2)*c5 + A(1,3)*c4 - A(1,4)*c3) * r
            AINV(1,3) = ( A(4,2)*s5 - A(4,3)*s4 + A(4,4)*s3) * r
            AINV(1,4) = (-A(3,2)*s5 + A(3,3)*s4 - A(3,4)*s3) * r
            AINV(2,1) = (-A(2,1)*c5 + A(2,3)*c2 - A(2,4)*c1) * r
            AINV(2,2) = ( A(1,1)*c5 - A(1,3)*c2 + A(1,4)*c1) * r
            AINV(2,3) = (-A(4,1)*s5 + A(4,3)*s2 - A(4,4)*s1) * r
            AINV(2,4) = ( A(3,1)*s5 - A(3,3)*s2 + A(3,4)*s1) * r
            AINV(3,1) = ( A(2,1)*c4 - A(2,2)*c2 + A(2,4)*c0) * r
            AINV(3,2) = (-A(1,1)*c4 + A(1,2)*c2 - A(1,4)*c0) * r
            AINV(3,3) = ( A(4,1)*s4 - A(4,2)*s2 + A(4,4)*s0) * r
            AINV(3,4) = (-A(3,1)*s4 + A(3,2)*s2 - A(3,4)*s0) * r
            AINV(4,1) = (-A(2,1)*c3 + A(2,2)*c1 - A(2,3)*c0) * r
            AINV(4,2) = ( A(1,1)*c3 - A(1,2)*c1 + A(1,3)*c0) * r
            AINV(4,3) = (-A(4,1)*s3 + A(4,2)*s1 - A(4,3)*s0) * r
            AINV(4,4) = ( A(3,1)*s3 - A(3,2)*s1 + A(3,3)*s0) * r
            kap = maxval(sum(abs(A), dim=2)) * maxval(sum(abs(AINV), dim=2))
            if (kap <= ADJ_KAPPA_MAX) then      ! also .false. for NaN
                OK_FLAG = .true.
                if (present(KAPPA)) KAPPA = kap
                return
            end if
        end if

        call lu_inv4(A, AINV, OK_FLAG)
        if (present(KAPPA)) then
            if (OK_FLAG) then
                KAPPA = maxval(sum(abs(A), dim=2)) * maxval(sum(abs(AINV), dim=2))
            else
                KAPPA = huge(1.0D0)
            end if
        end if
    END SUBROUTINE M44INV

    PURE SUBROUTINE lu_inv3(A, AINV, OK_FLAG)
        !! 3x3 LU inverse (fallback of `M33INV`), see `muscle_math_lu_inverse.inc`.
        implicit none
        integer, parameter :: N = 3
        REAL(real64), INTENT(IN)  :: A(N,N)
            !! Matrix to be inverted.
        REAL(real64), INTENT(OUT) :: AINV(N,N)
            !! Inverse of `A`; zero when `A` is singular.
        LOGICAL, INTENT(OUT)      :: OK_FLAG
            !! `.TRUE.` if the inversion succeeded.
        include 'muscle_math_lu_inverse.inc'
    END SUBROUTINE lu_inv3

    PURE SUBROUTINE lu_inv4(A, AINV, OK_FLAG)
        !! 4x4 LU inverse (fallback of `M44INV`), see `muscle_math_lu_inverse.inc`.
        implicit none
        integer, parameter :: N = 4
        REAL(real64), INTENT(IN)  :: A(N,N)
            !! Matrix to be inverted.
        REAL(real64), INTENT(OUT) :: AINV(N,N)
            !! Inverse of `A`; zero when `A` is singular.
        LOGICAL, INTENT(OUT)      :: OK_FLAG
            !! `.TRUE.` if the inversion succeeded.
        include 'muscle_math_lu_inverse.inc'
    END SUBROUTINE lu_inv4

    PURE SUBROUTINE M55INV(A, AINV, OK_FLAG)
        !! Inverse of a general 5x5 matrix by LU with partial pivoting in plain Fortran (algorithm of
        !! LAPACK `DGETRF` + `DGETRI`, see `muscle_math_lu_inverse.inc`). Used by `ten_2D4O`.
        implicit none
        integer, parameter :: N = 5
        REAL(real64), INTENT(IN)  :: A(N,N)
            !! Matrix to be inverted.
        REAL(real64), INTENT(OUT) :: AINV(N,N)
            !! Inverse of `A`; zero when `A` is singular.
        LOGICAL, INTENT(OUT)      :: OK_FLAG
            !! `.TRUE.` if the inversion succeeded, `.FALSE.` if `A` is singular.
        include 'muscle_math_lu_inverse.inc'
    END SUBROUTINE M55INV

end module