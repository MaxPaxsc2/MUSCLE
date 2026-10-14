! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_operator_transform
    !! Tests of the congruence transformation operator `.transform.`, one section per
    !! overload (module `muscle_tensor_ops_transform`):
    !!
    !! - `transform_3D2O_3D2Osym`:  \( \mathbf{F}\,\mathbf{S}\,\mathbf{F}^T \) (3D, 2nd order)
    !! - `transform_3D2O_3D4O2sym`: \( A_{iI}A_{jJ}A_{kK}A_{lL}C_{IJKL} \) (3D, minor symmetries)
    !! - `transform_3D2O_3D4O3sym`: same, fully symmetric 3D tensor
    !! - `transform_2D2O_2D2Osym`:  2D, 2nd order
    !! - `transform_2D2O_2D4O2sym`: 2D, minor symmetries
    !! - `transform_2D2O_2D4O3sym`: 2D, fully symmetric
    !!
    !! The 2D overloads are also checked against the 3D ones in `test_operator_2d_vs_3d`.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none

    logical :: passed

    print*, "Running tests for operator .transform. ..."

    ! --- transform_3D2O_3D2Osym ---
    call test_transform_3D2O_3D2Osym_identity(passed)
    if (.not. passed) stop 1
    call test_transform_3D2O_3D2Osym_stretch(passed)
    if (.not. passed) stop 2
    call test_transform_3D2O_3D2Osym_general(passed)
    if (.not. passed) stop 3
    call test_transform_3D2O_3D2Osym_definition(passed)
    if (.not. passed) stop 4
    call test_transform_3D2O_3D2Osym_roundtrip(passed)
    if (.not. passed) stop 5

    ! --- transform_3D2O_3D4O2sym ---
    call test_transform_3D2O_3D4O2sym_identity_stretch(passed)
    if (.not. passed) stop 6
    call test_transform_3D2O_3D4O2sym_rotation_90(passed)
    if (.not. passed) stop 7
    call test_transform_3D2O_3D4O2sym_isotropic(passed)
    if (.not. passed) stop 8
    call test_transform_3D2O_3D4O2sym_definition(passed)
    if (.not. passed) stop 9
    call test_transform_3D2O_3D4O2sym_properties(passed)
    if (.not. passed) stop 10
    call test_transform_3D2O_3D4O2sym_roundtrip(passed)
    if (.not. passed) stop 11

    ! --- transform_3D2O_3D4O3sym ---
    call test_transform_3D2O_3D4O3sym_stretch(passed)
    if (.not. passed) stop 12
    call test_transform_3D2O_3D4O3sym_isotropic(passed)
    if (.not. passed) stop 13
    call test_transform_3D2O_3D4O3sym_general(passed)
    if (.not. passed) stop 14

    ! --- transform_2D2O_2D2Osym ---
    call test_transform_2D2O_2D2Osym_stretch(passed)
    if (.not. passed) stop 15
    call test_transform_2D2O_2D2Osym_general(passed)
    if (.not. passed) stop 16

    ! --- transform_2D2O_2D4O2sym ---
    call test_transform_2D2O_2D4O2sym_rotation_90(passed)
    if (.not. passed) stop 17

    ! --- transform_2D2O_2D4O3sym ---
    call test_transform_2D2O_2D4O3sym_isotropic(passed)
    if (.not. passed) stop 18
    call test_transform_2D2O_2D4O3sym_stretch(passed)
    if (.not. passed) stop 19

    print*, "All .transform. operator tests passed successfully!"
    stop 0
end program test_operator_transform


! =========================================================================
! transform_3D2O_3D2Osym
! =========================================================================

subroutine test_transform_3D2O_3D2Osym_identity(passed)
    !! F = I gives res = S. Products of 0/1 and S, exact up to the zero terms; tol 1e-14.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2O) :: F_identity
    type(ten_3D2Osym) :: S

    passed = .true.
    call F_identity%init(xx=1.0D0, xy=0.0D0, xz=0.0D0, &
                         yx=0.0D0, yy=1.0D0, yz=0.0D0, &
                         zx=0.0D0, zy=0.0D0, zz=1.0D0)
    call S%init(xx=10.0D0, yy=-5.0D0, zz=2.0D0, xy=3.0D0, yz=-1.0D0, xz=4.0D0)

    call check_approx(passed, "I .transform. S = S", F_identity .transform. S, S, tol=1D-14)
end subroutine test_transform_3D2O_3D2Osym_identity


subroutine test_transform_3D2O_3D2Osym_stretch(passed)
    !! Diagonal stretch F = diag(2, 3, 4) (analytic): res_ij = F_ii F_jj S_ij.
    !! With S = (1, 2, 3, 4, 5, 6) in Voigt (xx, yy, zz, xy, yz, xz):
    !! R_11 = 4, R_22 = 18, R_33 = 48, R_12 = 24, R_23 = 60, R_13 = 48.
    !! Exact products of small integers; tol 1e-14.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2O) :: F_stretch
    type(ten_3D2Osym) :: S, expected

    passed = .true.
    call F_stretch%init(xx=2.0D0, xy=0.0D0, xz=0.0D0, &
                        yx=0.0D0, yy=3.0D0, yz=0.0D0, &
                        zx=0.0D0, zy=0.0D0, zz=4.0D0)
    call S%init(xx=1.0D0, yy=2.0D0, zz=3.0D0, xy=4.0D0, yz=5.0D0, xz=6.0D0)
    call expected%init(xx=4.0D0, yy=18.0D0, zz=48.0D0, xy=24.0D0, yz=60.0D0, xz=48.0D0)

    call check_approx(passed, "diag(2,3,4) .transform. S", F_stretch .transform. S, expected, tol=1D-14)
end subroutine test_transform_3D2O_3D2Osym_stretch


subroutine test_transform_3D2O_3D2Osym_general(passed)
    !! Shear + extension (analytic, computed by hand):
    !! F = [1.0, 0.5, 0.0 ; 0.0, 1.2, 0.0 ; 0.0, 0.0, 0.8],
    !! S = [100, 10, 0 ; 10, 50, 0 ; 0, 0, -20],
    !! T = F S = [105, 35, 0 ; 12, 60, 0 ; 0, 0, -16], R = T F^T:
    !! R_11 = 105 + 35*0.5 = 122.5, R_22 = 60*1.2 = 72, R_33 = -16*0.8 = -12.8,
    !! R_12 = 35*1.2 = 42, R_23 = R_13 = 0.
    !! Inputs such as 1.2 and 0.8 are not exact in binary; round-off ~1e-16 relative; tol 1e-14.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2O) :: F
    type(ten_3D2Osym) :: S, expected

    passed = .true.
    call F%init(xx=1.0D0, xy=0.5D0, xz=0.0D0, &
                yx=0.0D0, yy=1.2D0, yz=0.0D0, &
                zx=0.0D0, zy=0.0D0, zz=0.8D0)
    call S%init(xx=100.0D0, yy=50.0D0, zz=-20.0D0, xy=10.0D0, yz=0.0D0, xz=0.0D0)
    call expected%init(xx=122.5D0, yy=72.0D0, zz=-12.8D0, xy=42.0D0, yz=0.0D0, xz=0.0D0)

    call check_approx(passed, "F .transform. S (shear + extension)", F .transform. S, expected, tol=1D-14)
end subroutine test_transform_3D2O_3D2Osym_general


subroutine test_transform_3D2O_3D2Osym_definition(passed)
    !! Reference from the definition \( \mathbf{R} = \mathbf{F}\,\mathbf{S}\,\mathbf{F}^T \),
    !! evaluated with `matmul` on 3x3 matrices (independent of the unrolled library code), for
    !! a general F with all nine entries non-zero (non-symmetric, non-orthogonal) and S with all
    !! six components non-zero, so that every term of the library formula contributes.
    !! Tolerance: both evaluations sum 9 products of O(1) numbers in different orders;
    !! round-off ~1e-15 relative to max|R|; tol 1e-13.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D2O) :: F
    type(ten_3D2Osym) :: S, expected
    real(real64) :: f3(3,3), s3(3,3), r3(3,3)

    passed = .true.
    call F%init(xx=1.2D0, xy=0.4D0, xz=0.25D0, yx=-0.7D0, yy=0.9D0, yz=0.6D0, &
                zx=0.3D0, zy=-0.2D0, zz=1.1D0)
    call S%init(xx=100D0, yy=50D0, zz=-20D0, xy=10D0, yz=-15D0, xz=7.5D0)
    f3 = reshape(F%vals, (/3,3/))
    s3 = reshape((/ 100D0,  10D0,  7.5D0, &
                     10D0,  50D0, -15D0,  &
                    7.5D0, -15D0, -20D0 /), (/3,3/))
    r3 = matmul(f3, matmul(s3, transpose(f3)))
    call expected%init(xx=r3(1,1), yy=r3(2,2), zz=r3(3,3), xy=r3(1,2), yz=r3(2,3), xz=r3(1,3))

    call check_approx(passed, "F .transform. S vs matmul(F, matmul(S, F^T))", F .transform. S, &
                      expected, tol=1D-13)
end subroutine test_transform_3D2O_3D2Osym_definition


subroutine test_transform_3D2O_3D2Osym_roundtrip(passed)
    !! Push-forward followed by pull-back: \( \mathbf{F}^{-1} .transform. (\mathbf{F} .transform.
    !! \mathbf{S}) = \mathbf{S} \), since \( \mathbf{F}^{-1}\mathbf{F}\,\mathbf{S}\,\mathbf{F}^T
    !! \mathbf{F}^{-T} = \mathbf{S} \). F is general with an exact inverse (see
    !! `general_F_and_inverse`). Two chained transforms of O(1) data: round-off ~1e-15
    !! relative; tol 1e-13.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D2O) :: F, F_inv
    type(ten_3D2Osym) :: S

    passed = .true.
    call general_F_and_inverse(F, F_inv)
    call S%init(xx=100D0, yy=50D0, zz=-20D0, xy=10D0, yz=-15D0, xz=7.5D0)
    call check_approx(passed, "F^-1 .transform. (F .transform. S) = S", &
                      F_inv .transform. (F .transform. S), S, tol=1D-13)
    call check_approx(passed, "F .transform. (F^-1 .transform. S) = S", &
                      F .transform. (F_inv .transform. S), S, tol=1D-13)
end subroutine test_transform_3D2O_3D2Osym_roundtrip


! =========================================================================
! transform_3D2O_3D4O2sym
! =========================================================================

subroutine test_transform_3D2O_3D4O2sym_identity_stretch(passed)
    !! Fourth-order transformation with diagonal A (analytic):
    !! for \( A = \mathrm{diag}(a_1, a_2, a_3) \),
    !! \( res_{ijkl} = a_i a_j a_k a_l\, C_{ijkl} \) (no sum), i.e. in Voigt
    !! \( R_{\alpha\beta} = s_\alpha s_\beta C_{\alpha\beta} \) with
    !! \( s = (a_1^2, a_2^2, a_3^2, a_1 a_2, a_2 a_3, a_1 a_3) \). A = I gives R = C.
    !! C is a general minor-symmetric 6x6 matrix (no major symmetry).
    !! Tolerance: each entry is one product of exact small numbers plus a sum of
    !! zeros, so the error is a few ulp; 1e-14 relative is a safe bound.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D2O) :: A
    type(ten_3D4O2sym) :: C, expected
    real(real64) :: m(6,6), s(6)
    real(real64), parameter :: TOL = 1D-14
    integer :: i, j

    passed = .true.
    do j = 1, 6
        do i = 1, 6
            m(i,j) = 10D0*i + j + merge(50D0, 0D0, i == j)
        end do
    end do
    call C%init(m)

    call A%init(xx=1D0, xy=0D0, xz=0D0, yx=0D0, yy=1D0, yz=0D0, zx=0D0, zy=0D0, zz=1D0)
    call check_approx(passed, "I .transform. C = C", A .transform. C, C, tol=TOL)

    call A%init(xx=2D0, xy=0D0, xz=0D0, yx=0D0, yy=3D0, yz=0D0, zx=0D0, zy=0D0, zz=0.5D0)
    s = (/ 4D0, 9D0, 0.25D0, 6D0, 1.5D0, 1D0 /)
    do j = 1, 6
        do i = 1, 6
            m(i,j) = s(i)*s(j)*C%vals(i,j)
        end do
    end do
    call expected%init(m)
    call check_approx(passed, "diag(2,3,0.5) .transform. C", A .transform. C, expected, tol=TOL)
end subroutine test_transform_3D2O_3D4O2sym_identity_stretch


subroutine test_transform_3D2O_3D4O2sym_rotation_90(passed)
    !! Rotation by 90 degrees about z (analytic): e1 -> e2, e2 -> -e1, e3 -> e3, i.e.
    !! R_21 = 1, R_12 = -1, R_33 = 1. Then res_ijkl = C_{p(i)p(j)p(k)p(l)} times the signs of
    !! the R entries used, which in Voigt (11, 22, 33, 12, 23, 13) reads
    !! \( R_{\alpha\beta} = \sigma_\alpha \sigma_\beta\, C_{\pi(\alpha)\pi(\beta)} \),
    !! \( \pi = (2, 1, 3, 4, 6, 5) \), \( \sigma = (1, 1, 1, -1, 1, -1) \)
    !! (e.g. res_1112 = R_12 R_12 R_12 R_21 C_2221 = -C_2212; res_1123 = C_2213).
    !! Entries are products of +-1 and C, so the result is exact up to the zero terms;
    !! tolerance 1e-14 relative.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D2O) :: Rot
    type(ten_3D4O2sym) :: C, expected
    real(real64) :: m(6,6)
    integer, parameter :: perm(6) = (/2, 1, 3, 4, 6, 5/)
    real(real64), parameter :: sgn(6) = (/1D0, 1D0, 1D0, -1D0, 1D0, -1D0/)
    integer :: i, j

    passed = .true.
    do j = 1, 6
        do i = 1, 6
            m(i,j) = 10D0*i + j + merge(50D0, 0D0, i == j)
        end do
    end do
    call C%init(m)
    do j = 1, 6
        do i = 1, 6
            m(i,j) = sgn(i)*sgn(j)*C%vals(perm(i), perm(j))
        end do
    end do
    call expected%init(m)

    call Rot%init(xx=0D0, xy=-1D0, xz=0D0, yx=1D0, yy=0D0, yz=0D0, zx=0D0, zy=0D0, zz=1D0)
    call check_approx(passed, "Rz(90) .transform. C", Rot .transform. C, expected, tol=1D-14)
end subroutine test_transform_3D2O_3D4O2sym_rotation_90


subroutine test_transform_3D2O_3D4O2sym_isotropic(passed)
    !! Analytic invariant: the isotropic tensor \( \lambda\,\mathbf{I}\otimes\mathbf{I} +
    !! 2\mu\,\mathbb{I}^S \) is unchanged by any rotation Q.
    !! Q is the Rodrigues rotation of angle 0.7 rad about the axis (1, 2, 2)/3 (all entries
    !! non-zero). Voigt (tensorial): normal block lambda + 2 mu on the diagonal, lambda off it,
    !! shear diagonal mu.
    !! Tolerance: Q is orthogonal to ~1e-16 and each output entry sums 36 products of O(1)
    !! terms twice, so the error is ~1e-15 relative to max|C|; 1e-13 leaves a wide margin.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D2O) :: Q
    type(ten_3D4O2sym) :: C
    real(real64) :: m(6,6)
    real(real64), parameter :: LAMBDA = 121.15D0, MU = 80.77D0
    integer :: i
    interface
        function rodrigues_rotation() result(Q)
            import :: ten_3D2O
            type(ten_3D2O) :: Q
        end function rodrigues_rotation
    end interface

    passed = .true.
    m = 0D0
    m(1:3,1:3) = LAMBDA
    do i = 1, 3
        m(i,i) = LAMBDA + 2D0*MU
        m(i+3,i+3) = MU
    end do
    call C%init(m)
    Q = rodrigues_rotation()

    call check_approx(passed, "Q .transform. C_iso = C_iso", Q .transform. C, C, tol=1D-13)
end subroutine test_transform_3D2O_3D4O2sym_isotropic


subroutine test_transform_3D2O_3D4O2sym_definition(passed)
    !! Reference from the definition: res_ijkl = sum_{IJKL} A_iI A_jJ A_kK A_lL C_IJKL, evaluated
    !! here directly on the 3x3x3x3 array built from the minor-symmetric Voigt matrix
    !! (independent of the library's Voigt formulation). A is general (non-symmetric,
    !! non-orthogonal) and C has no major symmetry.
    !! Tolerance: both evaluations sum products of O(1) numbers with different orderings;
    !! the round-off difference is ~1e-15 relative to max|res|, so 1e-13 is a safe bound.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D2O) :: A
    type(ten_3D4O2sym) :: C, expected
    real(real64) :: m(6,6), a3(3,3), c4(3,3,3,3), r4(3,3,3,3)
    integer, parameter :: vi(6) = (/1, 2, 3, 1, 2, 1/), vj(6) = (/1, 2, 3, 2, 3, 3/)
    integer :: i, j, k, l, ii, jj, kk, ll

    passed = .true.
    do j = 1, 6
        do i = 1, 6
            m(i,j) = sin(1D0*i + 3D0*j) + merge(5D0, 0D0, i == j)
        end do
    end do
    call C%init(m)
    call A%init(xx=1.2D0, xy=0.4D0, xz=0.25D0, yx=-0.7D0, yy=0.9D0, yz=0.6D0, &
                zx=0.3D0, zy=-0.2D0, zz=1.1D0)
    a3 = reshape(A%vals, (/3,3/))

    do j = 1, 6
        do i = 1, 6
            c4(vi(i), vj(i), vi(j), vj(j)) = m(i,j)
            c4(vj(i), vi(i), vi(j), vj(j)) = m(i,j)
            c4(vi(i), vj(i), vj(j), vi(j)) = m(i,j)
            c4(vj(i), vi(i), vj(j), vi(j)) = m(i,j)
        end do
    end do
    r4 = 0D0
    do l = 1, 3; do k = 1, 3; do j = 1, 3; do i = 1, 3
        do ll = 1, 3; do kk = 1, 3; do jj = 1, 3; do ii = 1, 3
            r4(i,j,k,l) = r4(i,j,k,l) + a3(i,ii)*a3(j,jj)*a3(k,kk)*a3(l,ll)*c4(ii,jj,kk,ll)
        end do; end do; end do; end do
    end do; end do; end do; end do
    do j = 1, 6
        do i = 1, 6
            m(i,j) = r4(vi(i), vj(i), vi(j), vj(j))
        end do
    end do
    call expected%init(m)

    call check_approx(passed, "A .transform. C vs index definition", A .transform. C, expected, &
                      tol=1D-13)
end subroutine test_transform_3D2O_3D4O2sym_definition


subroutine test_transform_3D2O_3D4O2sym_properties(passed)
    !! Properties that follow from the definition:
    !! - composition: (A B) .transform. C = A .transform. (B .transform. C), for general A, B;
    !! - consistency with the 2nd-order transformation for a rotation Q:
    !!   (Q .transform. C) : (Q .transform. S) = Q .transform. (C : S).
    !! Tolerance: two chained transformations of O(1) data, round-off ~1e-15 relative; 1e-13.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D2O) :: A, B, Q
    type(ten_3D4O2sym) :: C
    type(ten_3D2Osym) :: S
    real(real64) :: m(6,6)
    integer :: i, j

    passed = .true.
    do j = 1, 6
        do i = 1, 6
            m(i,j) = cos(2D0*i - 1D0*j) + merge(4D0, 0D0, i == j)
        end do
    end do
    call C%init(m)
    call A%init(xx=1.2D0, xy=0.4D0, xz=0.25D0, yx=-0.7D0, yy=0.9D0, yz=0.6D0, &
                zx=0.3D0, zy=-0.2D0, zz=1.1D0)
    call B%init(xx=0.8D0, xy=-0.3D0, xz=0.5D0, yx=0.2D0, yy=1.3D0, yz=-0.4D0, &
                zx=-0.6D0, zy=0.1D0, zz=0.95D0)
    call check_approx(passed, "(A*B) .transform. C = A .transform. (B .transform. C)", &
                      (A*B) .transform. C, A .transform. (B .transform. C), tol=1D-13)

    ! Rotation e1 -> e2 -> e3 -> e1 (120 deg about (1,1,1); exact entries 0 and 1)
    call Q%init(xx=0D0, xy=0D0, xz=1D0, yx=1D0, yy=0D0, yz=0D0, zx=0D0, zy=1D0, zz=0D0)
    call S%init(xx=1D0, yy=-2D0, zz=0.5D0, xy=0.3D0, yz=-0.7D0, xz=1.1D0)
    call check_approx(passed, "(Q .transform. C) : (Q .transform. S) = Q .transform. (C : S)", &
                      (Q .transform. C) .ddot. (Q .transform. S), Q .transform. (C .ddot. S), &
                      tol=1D-13)
end subroutine test_transform_3D2O_3D4O2sym_properties


subroutine test_transform_3D2O_3D4O2sym_roundtrip(passed)
    !! Push-forward followed by pull-back: A^-1 .transform. (A .transform. C) = C, because
    !! (A^-1)_{iI} A_{Ia} = delta_{ia} in each of the four index slots. A is general with an
    !! exact inverse (see `general_F_and_inverse`); C has no major symmetry.
    !! Two chained transforms of O(1) data: round-off ~1e-15 relative; tol 1e-13.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D2O) :: A, A_inv
    type(ten_3D4O2sym) :: C
    real(real64) :: m(6,6)
    integer :: i, j

    passed = .true.
    do j = 1, 6
        do i = 1, 6
            m(i,j) = sin(1D0*i + 3D0*j) + merge(5D0, 0D0, i == j)
        end do
    end do
    call C%init(m)
    call general_F_and_inverse(A, A_inv)
    call check_approx(passed, "A^-1 .transform. (A .transform. C) = C", &
                      A_inv .transform. (A .transform. C), C, tol=1D-13)
    call check_approx(passed, "A .transform. (A^-1 .transform. C) = C", &
                      A .transform. (A_inv .transform. C), C, tol=1D-13)
end subroutine test_transform_3D2O_3D4O2sym_roundtrip


! =========================================================================
! transform_3D2O_3D4O3sym
! =========================================================================

subroutine test_transform_3D2O_3D4O3sym_stretch(passed)
    !! Diagonal A = diag(2, 3, 0.5) (analytic): \( R_{\alpha\beta} = s_\alpha s_\beta C_{\alpha\beta} \),
    !! \( s = (a_1^2, a_2^2, a_3^2, a_1 a_2, a_2 a_3, a_1 a_3) \), on a symmetric C with all 21
    !! independent entries distinct. Exact products of small numbers; tol 1e-14.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D2O) :: A
    type(ten_3D4O3sym) :: C, expected
    type(ten_3D4O2sym) :: C_full, r_full
    real(real64) :: m(6,6), sc(6)
    integer :: i, j

    passed = .true.
    do j = 1, 6
        do i = 1, j
            m(i,j) = 10D0*i + j + merge(50D0, 0D0, i == j)
            m(j,i) = m(i,j)
        end do
    end do
    call C_full%init(m)
    C = C_full%convert_3sym()

    call A%init(xx=2D0, xy=0D0, xz=0D0, yx=0D0, yy=3D0, yz=0D0, zx=0D0, zy=0D0, zz=0.5D0)
    sc = (/ 4D0, 9D0, 0.25D0, 6D0, 1.5D0, 1D0 /)
    do j = 1, 6
        do i = 1, 6
            m(i,j) = sc(i)*sc(j)*C_full%vals(i,j)
        end do
    end do
    call r_full%init(m)
    expected = r_full%convert_3sym()
    call check_approx(passed, "diag(2,3,0.5) .transform. C", A .transform. C, expected, tol=1D-14)
end subroutine test_transform_3D2O_3D4O3sym_stretch


subroutine test_transform_3D2O_3D4O3sym_isotropic(passed)
    !! Analytic invariant: the isotropic tensor is unchanged by the Rodrigues rotation of
    !! angle 0.7 rad about (1, 2, 2)/3. Round-off ~1e-15 relative; tol 1e-13.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D2O) :: Q
    type(ten_3D4O3sym) :: C
    type(ten_3D4O2sym) :: C_full
    real(real64) :: m(6,6)
    real(real64), parameter :: LAMBDA = 121.15D0, MU = 80.77D0
    integer :: i
    interface
        function rodrigues_rotation() result(Q)
            import :: ten_3D2O
            type(ten_3D2O) :: Q
        end function rodrigues_rotation
    end interface

    passed = .true.
    m = 0D0
    m(1:3,1:3) = LAMBDA
    do i = 1, 3
        m(i,i) = LAMBDA + 2D0*MU
        m(i+3,i+3) = MU
    end do
    call C_full%init(m)
    C = C_full%convert_3sym()
    Q = rodrigues_rotation()

    call check_approx(passed, "Q .transform. C_iso = C_iso", Q .transform. C, C, tol=1D-13)
end subroutine test_transform_3D2O_3D4O3sym_isotropic


subroutine test_transform_3D2O_3D4O3sym_general(passed)
    !! General A: the result equals the minor-symmetric transformation of the expanded tensor
    !! (itself checked against the index definition in
    !! `test_transform_3D2O_3D4O2sym_definition`), and that full result is symmetric (major
    !! symmetry preserved). Both evaluate the same sums, so they agree to round-off; tol 1e-13.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D2O) :: A
    type(ten_3D4O3sym) :: C, r3
    type(ten_3D4O2sym) :: C_full, r_full, r3_full
    real(real64) :: m(6,6)
    integer :: i, j

    passed = .true.
    do j = 1, 6
        do i = 1, j
            m(i,j) = 10D0*i + j + merge(50D0, 0D0, i == j)
            m(j,i) = m(i,j)
        end do
    end do
    call C_full%init(m)
    C = C_full%convert_3sym()

    call A%init(xx=1.2D0, xy=0.4D0, xz=0.25D0, yx=-0.7D0, yy=0.9D0, yz=0.6D0, &
                zx=0.3D0, zy=-0.2D0, zz=1.1D0)
    r3 = A .transform. C
    r_full = A .transform. C_full
    r3_full = r3
    call check_approx(passed, "A .transform. C (3sym) vs 2sym version", r3_full, r_full, tol=1D-13)
    call check_approx(passed, "A .transform. C (2sym) keeps major symmetry", r_full%vals, &
                      transpose(r_full%vals), tol=1D-13)
end subroutine test_transform_3D2O_3D4O3sym_general


! =========================================================================
! transform_2D2O_2D2Osym
! =========================================================================

subroutine test_transform_2D2O_2D2Osym_stretch(passed)
    !! Diagonal stretch F = diag(2, 3, 0.5) (analytic): res_ij = F_ii F_jj S_ij.
    !! With S = (xx, yy, zz, xy) = (1, 2, 3, 4): R = (4, 18, 0.75, 24).
    !! Exact products of small numbers; tol 1e-14.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2O) :: F
    type(ten_2D2Osym) :: S, expected

    passed = .true.
    call F%init(xx=2D0, xy=0D0, yx=0D0, yy=3D0, zz=0.5D0)
    call S%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0)
    call expected%init(xx=4D0, yy=18D0, zz=0.75D0, xy=24D0)
    call check_approx(passed, "diag(2,3,0.5) .transform. S", F .transform. S, expected, tol=1D-14)
end subroutine test_transform_2D2O_2D2Osym_stretch


subroutine test_transform_2D2O_2D2Osym_general(passed)
    !! General 2D F with both in-plane shears and an out-of-plane stretch (analytic, by hand):
    !! F = [[1, 0.5], [0.25, 2]], F_zz = 1.5;  S = [[2, 0.5], [0.5, -1]], S_zz = 4.
    !! T = F S = [[2.25, 0], [1.5, -1.875]];  R = T F^T:
    !! R_11 = 2.25*1 + 0*0.5 = 2.25,  R_22 = 1.5*0.25 - 1.875*2 = -3.375,
    !! R_12 = 2.25*0.25 + 0*2 = 0.5625 (= R_21 = 1.5*1 - 1.875*0.5),  R_zz = 1.5^2 * 4 = 9.
    !! All inputs are dyadic fractions, so every operation is exact; tol 1e-14.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2O) :: F
    type(ten_2D2Osym) :: S, expected

    passed = .true.
    call F%init(xx=1D0, xy=0.5D0, yx=0.25D0, yy=2D0, zz=1.5D0)
    call S%init(xx=2D0, yy=-1D0, zz=4D0, xy=0.5D0)
    call expected%init(xx=2.25D0, yy=-3.375D0, zz=9D0, xy=0.5625D0)
    call check_approx(passed, "F .transform. S (general 2D F)", F .transform. S, expected, tol=1D-14)
end subroutine test_transform_2D2O_2D2Osym_general


! =========================================================================
! transform_2D2O_2D4O2sym
! =========================================================================

subroutine test_transform_2D2O_2D4O2sym_rotation_90(passed)
    !! Analytic case: an in-plane rotation by 90 degrees (e1 -> e2, e2 -> -e1, e3 -> e3)
    !! maps C_1111 <-> C_2222, C_1133 <-> C_2233, C_1112 -> -C_2221 = -C_2212, etc.
    !! For R = [[0,-1,0],[1,0,0],[0,0,1]]: res_ijkl = R_iI R_jJ R_kK R_lL C_IJKL, so
    !! res_1111 = C_2222, res_2222 = C_1111, res_1122 = C_2211, res_3333 = C_3333,
    !! res_1133 = C_2233, res_2233 = C_1133, res_1212 = C_2121 = C_1212,
    !! res_1112 = -C_2221 = -C_2212, res_2212 = -C_1112, res_3312 = -C_3312.
    !! C has no major symmetry. Entries are products of 0/+-1, so the result is exact up to
    !! the zero terms; tol 1e-14.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2O) :: Rot
    type(ten_2D4O2sym) :: C, expected
    real(real64) :: c4(4,4)

    passed = .true.
    call C%init(reshape((/ 5D0, -1D0, 0.5D0, 1D0, &
                           1D0,  6D0, 2D0, -0.5D0, &
                           2D0,  1D0, 7D0, 0.25D0, &
                           0.5D0, 1D0, -1D0, 3D0 /), (/4,4/)))
    call Rot%init(xx=0D0, xy=-1D0, yx=1D0, yy=0D0, zz=1D0)
    c4 = C%vals
    expected%vals(1,:) = (/ c4(2,2), c4(2,1), c4(2,3), -c4(2,4) /)
    expected%vals(2,:) = (/ c4(1,2), c4(1,1), c4(1,3), -c4(1,4) /)
    expected%vals(3,:) = (/ c4(3,2), c4(3,1), c4(3,3), -c4(3,4) /)
    expected%vals(4,:) = (/ -c4(4,2), -c4(4,1), -c4(4,3), c4(4,4) /)
    call check_approx(passed, "Rz(90) .transform. C", Rot .transform. C, expected, tol=1D-14)
end subroutine test_transform_2D2O_2D4O2sym_rotation_90


! =========================================================================
! transform_2D2O_2D4O3sym
! =========================================================================

subroutine test_transform_2D2O_2D4O3sym_isotropic(passed)
    !! Analytic invariant: the isotropic tensor lambda I(x)I + 2 mu I^S (Voigt: normal block
    !! lambda + 2 mu / lambda, xyxy = mu) is unchanged by an in-plane rotation of 0.7 rad
    !! (A = [[c, -s], [s, c]], A_zz = 1). Round-off ~1e-15 relative; tol 1e-13.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2O) :: A
    type(ten_2D4O3sym) :: C
    real(real64), parameter :: LAMBDA = 121.15D0, MU = 80.77D0, THETA = 0.7D0

    passed = .true.
    call C%init(xxxx=LAMBDA + 2D0*MU, yyyy=LAMBDA + 2D0*MU, zzzz=LAMBDA + 2D0*MU, xyxy=MU, &
                xxyy=LAMBDA, yyzz=LAMBDA, zzxy=0D0, xxzz=LAMBDA, yyxy=0D0, xxxy=0D0)
    call A%init(xx=cos(THETA), xy=-sin(THETA), yx=sin(THETA), yy=cos(THETA), zz=1D0)
    call check_approx(passed, "R(0.7) .transform. C_iso = C_iso", A .transform. C, C, tol=1D-13)
end subroutine test_transform_2D2O_2D4O3sym_isotropic


subroutine test_transform_2D2O_2D4O3sym_stretch(passed)
    !! Diagonal A = diag(2, 3, 0.5) (analytic): R_ab = s_a s_b C_ab with
    !! s = (a_1^2, a_2^2, a_3^2, a_1 a_2). Exact products of small numbers; tol 1e-14.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2O) :: A
    type(ten_2D4O3sym) :: C, expected
    real(real64) :: sc(4)

    passed = .true.
    ! Storage (11,22,33,44,12,23,34,13,24,14) -> scale s_a s_b of each (a,b)
    call C%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0, 7D0, 8D0, 9D0, 10D0/))
    call A%init(xx=2D0, xy=0D0, yx=0D0, yy=3D0, zz=0.5D0)
    sc = (/4D0, 9D0, 0.25D0, 6D0/)
    call expected%init((/ sc(1)*sc(1)*1D0, sc(2)*sc(2)*2D0, sc(3)*sc(3)*3D0, sc(4)*sc(4)*4D0, &
                          sc(1)*sc(2)*5D0, sc(2)*sc(3)*6D0, sc(3)*sc(4)*7D0,                 &
                          sc(1)*sc(3)*8D0, sc(2)*sc(4)*9D0,                                  &
                          sc(1)*sc(4)*10D0 /))
    call check_approx(passed, "diag(2,3,0.5) .transform. C", A .transform. C, expected, tol=1D-14)
end subroutine test_transform_2D2O_2D4O3sym_stretch


! =========================================================================
! Shared helpers
! =========================================================================

function rodrigues_rotation() result(Q)
    !! Rotation of angle 0.7 rad about the unit axis n = (1, 2, 2)/3 (Rodrigues formula
    !! \( Q = \cos\theta\,\mathbf{I} + \sin\theta\,[\mathbf{n}]_\times + (1-\cos\theta)\,
    !! \mathbf{n}\otimes\mathbf{n} \)); all nine entries are non-zero.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    type(ten_3D2O) :: Q
        !! Rotation tensor.
    real(real64), parameter :: THETA = 0.7D0
    real(real64) :: q3(3,3), n(3), ct, st

    n = (/1D0, 2D0, 2D0/) / 3D0
    ct = cos(THETA); st = sin(THETA)
    q3(1,:) = (/ ct + n(1)**2*(1-ct),        n(1)*n(2)*(1-ct) - n(3)*st, n(1)*n(3)*(1-ct) + n(2)*st /)
    q3(2,:) = (/ n(2)*n(1)*(1-ct) + n(3)*st, ct + n(2)**2*(1-ct),        n(2)*n(3)*(1-ct) - n(1)*st /)
    q3(3,:) = (/ n(3)*n(1)*(1-ct) - n(2)*st, n(3)*n(2)*(1-ct) + n(1)*st, ct + n(3)**2*(1-ct)        /)
    call Q%init(reshape(q3, (/9/)))
end function rodrigues_rotation


subroutine general_F_and_inverse(F, F_inv)
    !! General 3D tensor with an exactly known inverse: \( \mathbf{F} = \mathbf{L}\,\mathbf{U} \)
    !! with unit lower/upper triangular factors
    !! L = [[1,0,0],[L21,1,0],[L31,L32,1]], U = [[1,U12,U13],[0,1,U23],[0,0,1]], whose inverses are
    !! L^-1 = [[1,0,0],[-L21,1,0],[L21 L32 - L31,-L32,1]],
    !! U^-1 = [[1,-U12,U12 U23 - U13],[0,1,-U23],[0,0,1]],
    !! so \( \mathbf{F}^{-1} = \mathbf{U}^{-1}\mathbf{L}^{-1} \) and det F = 1.
    !! With the dyadic values below all nine entries of F are non-zero and every product is
    !! exact in double precision.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    type(ten_3D2O), intent(out) :: F
        !! General tensor (non-symmetric, non-orthogonal).
    type(ten_3D2O), intent(out) :: F_inv
        !! Its exact inverse.
    real(real64), parameter :: L21 = 0.5D0, L31 = -0.25D0, L32 = 0.75D0
    real(real64), parameter :: U12 = 0.25D0, U13 = 0.5D0, U23 = -0.5D0
    real(real64) :: L(3,3), U(3,3), L_inv(3,3), U_inv(3,3)

    L     = reshape((/ 1D0, L21, L31,   0D0, 1D0, L32,   0D0, 0D0, 1D0 /), (/3,3/))
    U     = reshape((/ 1D0, 0D0, 0D0,   U12, 1D0, 0D0,   U13, U23, 1D0 /), (/3,3/))
    L_inv = reshape((/ 1D0, -L21, L21*L32 - L31,   0D0, 1D0, -L32,   0D0, 0D0, 1D0 /), (/3,3/))
    U_inv = reshape((/ 1D0, 0D0, 0D0,   -U12, 1D0, 0D0,   U12*U23 - U13, -U23, 1D0 /), (/3,3/))
    call F%init(reshape(matmul(L, U), (/9/)))
    call F_inv%init(reshape(matmul(U_inv, L_inv), (/9/)))
end subroutine general_F_and_inverse
