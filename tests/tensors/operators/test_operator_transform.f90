! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_operator_transform
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none

    logical :: passed

    print*, "Running tests for operator .transform. (A .transform. S)..."

    call test_transform_identity(passed)
    if (.not. passed) stop 1

    call test_transform_scaling(passed)
    if (.not. passed) stop 2

    call test_transform_general(passed)
    if (.not. passed) stop 3

    call test_transform4_identity_stretch(passed)
    if (.not. passed) stop 4

    call test_transform4_rotation_90(passed)
    if (.not. passed) stop 5

    call test_transform4_isotropic_invariance(passed)
    if (.not. passed) stop 6

    call test_transform4_definition(passed)
    if (.not. passed) stop 7

    call test_transform4_properties(passed)
    if (.not. passed) stop 8

    call test_transform4_3sym(passed)
    if (.not. passed) stop 9

    print*, "All .transform. operator tests passed successfully!"
    stop 0
end program test_operator_transform


subroutine test_transform_identity(passed)
    !! Test 1: Transformation with Identity Tensor (F = I => res = S)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2O) :: F_identity
    type(ten_3D2Osym) :: S, res

    ! Initialize F as Identity (Column-major storage)
    call F_identity%init(xx=1.0D0, xy=0.0D0, xz=0.0D0, &
                         yx=0.0D0, yy=1.0D0, yz=0.0D0, &
                         zx=0.0D0, zy=0.0D0, zz=1.0D0)

    ! Initialize S with arbitrary values
    call S%init(xx=10.0D0, yy=-5.0D0, zz=2.0D0, xy=3.0D0, yz=-1.0D0, xz=4.0D0)

    ! Perform transformation
    res = F_identity .transform. S

    passed = res .approx. S
    if (.not. passed) then
        print*, "FAIL: Identity Transformation (F = I)", new_line('A'), &
                "Expected:", S%vals, new_line('A'), &
                "Actual  :", res%vals
    end if
end subroutine test_transform_identity


subroutine test_transform_scaling(passed)
    !! Test 2: Diagonal Stretch Transformation (F = diag(2, 3, 4))
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2O) :: F_stretch
    type(ten_3D2Osym) :: S, res, expected

    ! F = diag(2, 3, 4)
    call F_stretch%init(xx=2.0D0, xy=0.0D0, xz=0.0D0, &
                        yx=0.0D0, yy=3.0D0, yz=0.0D0, &
                        zx=0.0D0, zy=0.0D0, zz=4.0D0)

    ! S = (1, 2, 3, 4, 5, 6) in Voigt (xx, yy, zz, xy, yz, xz)
    call S%init(xx=1.0D0, yy=2.0D0, zz=3.0D0, xy=4.0D0, yz=5.0D0, xz=6.0D0)

    ! Analytically:
    ! R_11 = 2^2 * 1 = 4
    ! R_22 = 3^2 * 2 = 18
    ! R_33 = 4^2 * 3 = 48
    ! R_12 = 2 * 3 * 4 = 24
    ! R_23 = 3 * 4 * 5 = 60
    ! R_13 = 2 * 4 * 6 = 48
    call expected%init(xx=4.0D0, yy=18.0D0, zz=48.0D0, xy=24.0D0, yz=60.0D0, xz=48.0D0)

    res = F_stretch .transform. S

    passed = res .approx. expected
    if (.not. passed) then
        print*, "FAIL: Diagonal Scaling Transformation", new_line('A'), &
                "Expected:", expected%vals, new_line('A'), &
                "Actual  :", res%vals
    end if
end subroutine test_transform_scaling


subroutine test_transform_general(passed)
    !! Test 3: General Shear + Extension Deformation Gradient Transformation
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2O) :: F
    type(ten_3D2Osym) :: S, res, expected

    ! F = [1.0, 0.5, 0.0 ; 0.0, 1.2, 0.0 ; 0.0, 0.0, 0.8]
    call F%init(xx=1.0D0, xy=0.5D0, xz=0.0D0, &
                yx=0.0D0, yy=1.2D0, yz=0.0D0, &
                zx=0.0D0, zy=0.0D0, zz=0.8D0)

    ! S = [100.0, 10.0, 0.0 ; 10.0, 50.0, 0.0 ; 0.0, 0.0, -20.0]
    call S%init(xx=100.0D0, yy=50.0D0, zz=-20.0D0, xy=10.0D0, yz=0.0D0, xz=0.0D0)

    ! Analytical calculation of F * S * F^T:
    ! T = F * S = [105.0, 35.0, 0.0 ; 12.0, 60.0, 0.0 ; 0.0, 0.0, -16.0]
    ! R = T * F^T:
    ! R_11 = 105.0 * 1.0 + 35.0 * 0.5 = 122.5
    ! R_22 = 60.0 * 1.2 = 72.0
    ! R_33 = -16.0 * 0.8 = -12.8
    ! R_12 = 35.0 * 1.2 = 42.0
    ! R_23 = 0.0
    ! R_13 = 0.0
    call expected%init(xx=122.5D0, yy=72.0D0, zz=-12.8D0, xy=42.0D0, yz=0.0D0, xz=0.0D0)

    res = F .transform. S

    passed = res .approx. expected
    if (.not. passed) then
        print*, "FAIL: General Transformation (F * S * F^T)", new_line('A'), &
                "Expected:", expected%vals, new_line('A'), &
                "Actual  :", res%vals, new_line('A'), &
                "Diff    :", res%vals - expected%vals
    end if
end subroutine test_transform_general


subroutine test_transform4_identity_stretch(passed)
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
end subroutine test_transform4_identity_stretch


subroutine test_transform4_rotation_90(passed)
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
end subroutine test_transform4_rotation_90


subroutine test_transform4_isotropic_invariance(passed)
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
    real(real64) :: m(6,6), q3(3,3), n(3), ct, st
    real(real64), parameter :: LAMBDA = 121.15D0, MU = 80.77D0, THETA = 0.7D0
    integer :: i

    passed = .true.
    m = 0D0
    m(1:3,1:3) = LAMBDA
    do i = 1, 3
        m(i,i) = LAMBDA + 2D0*MU
        m(i+3,i+3) = MU
    end do
    call C%init(m)

    n = (/1D0, 2D0, 2D0/) / 3D0
    ct = cos(THETA); st = sin(THETA)
    q3(1,:) = (/ ct + n(1)**2*(1-ct),        n(1)*n(2)*(1-ct) - n(3)*st, n(1)*n(3)*(1-ct) + n(2)*st /)
    q3(2,:) = (/ n(2)*n(1)*(1-ct) + n(3)*st, ct + n(2)**2*(1-ct),        n(2)*n(3)*(1-ct) - n(1)*st /)
    q3(3,:) = (/ n(3)*n(1)*(1-ct) - n(2)*st, n(3)*n(2)*(1-ct) + n(1)*st, ct + n(3)**2*(1-ct)        /)
    call Q%init(reshape(q3, (/9/)))

    call check_approx(passed, "Q .transform. C_iso = C_iso", Q .transform. C, C, tol=1D-13)
end subroutine test_transform4_isotropic_invariance


subroutine test_transform4_definition(passed)
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
end subroutine test_transform4_definition


subroutine test_transform4_properties(passed)
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
end subroutine test_transform4_properties


subroutine test_transform4_3sym(passed)
    !! `.transform.` on the fully symmetric type `ten_3D4O3sym`:
    !! - analytic, diagonal A = diag(2, 3, 0.5): \( R_{\alpha\beta} = s_\alpha s_\beta C_{\alpha\beta} \),
    !!   \( s = (a_1^2, a_2^2, a_3^2, a_1 a_2, a_2 a_3, a_1 a_3) \) (exact products; tol 1e-14);
    !! - analytic invariant: the isotropic tensor is unchanged by the Rodrigues rotation of
    !!   angle 0.7 rad about (1, 2, 2)/3 (round-off ~1e-15 relative; tol 1e-13);
    !! - general A: the result equals the minor-symmetric transformation of the expanded tensor
    !!   (itself checked against the index definition in `test_transform4_definition`), and
    !!   that full result is symmetric (major symmetry preserved). Both evaluate the same sums,
    !!   so they agree to round-off; tol 1e-13.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D2O) :: A, Q
    type(ten_3D4O3sym) :: C, expected, r3
    type(ten_3D4O2sym) :: C_full, r_full, r3_full
    real(real64) :: m(6,6), q3(3,3), n(3), sc(6), ct, st
    real(real64), parameter :: LAMBDA = 121.15D0, MU = 80.77D0, THETA = 0.7D0
    integer :: i, j

    passed = .true.
    ! Symmetric Voigt matrix with all 21 independent entries distinct
    do j = 1, 6
        do i = 1, j
            m(i,j) = 10D0*i + j + merge(50D0, 0D0, i == j)
            m(j,i) = m(i,j)
        end do
    end do
    C_full = C_full_from(m)
    C = C_full%convert_3sym()

    ! Diagonal stretch
    call A%init(xx=2D0, xy=0D0, xz=0D0, yx=0D0, yy=3D0, yz=0D0, zx=0D0, zy=0D0, zz=0.5D0)
    sc = (/ 4D0, 9D0, 0.25D0, 6D0, 1.5D0, 1D0 /)
    do j = 1, 6
        do i = 1, 6
            m(i,j) = sc(i)*sc(j)*C_full%vals(i,j)
        end do
    end do
    r_full = C_full_from(m)
    expected = r_full%convert_3sym()
    call check_approx(passed, "diag(2,3,0.5) .transform. C (3sym)", A .transform. C, expected, tol=1D-14)

    ! Isotropic invariance
    m = 0D0
    m(1:3,1:3) = LAMBDA
    do i = 1, 3
        m(i,i) = LAMBDA + 2D0*MU
        m(i+3,i+3) = MU
    end do
    r_full = C_full_from(m)
    expected = r_full%convert_3sym()
    n = (/1D0, 2D0, 2D0/) / 3D0
    ct = cos(THETA); st = sin(THETA)
    q3(1,:) = (/ ct + n(1)**2*(1-ct),        n(1)*n(2)*(1-ct) - n(3)*st, n(1)*n(3)*(1-ct) + n(2)*st /)
    q3(2,:) = (/ n(2)*n(1)*(1-ct) + n(3)*st, ct + n(2)**2*(1-ct),        n(2)*n(3)*(1-ct) - n(1)*st /)
    q3(3,:) = (/ n(3)*n(1)*(1-ct) - n(2)*st, n(3)*n(2)*(1-ct) + n(1)*st, ct + n(3)**2*(1-ct)        /)
    call Q%init(reshape(q3, (/9/)))
    call check_approx(passed, "Q .transform. C_iso = C_iso (3sym)", Q .transform. expected, expected, &
                      tol=1D-13)

    ! General A: consistency with the minor-symmetric version and preserved major symmetry
    call A%init(xx=1.2D0, xy=0.4D0, xz=0.25D0, yx=-0.7D0, yy=0.9D0, yz=0.6D0, &
                zx=0.3D0, zy=-0.2D0, zz=1.1D0)
    r3 = A .transform. C
    r_full = A .transform. C_full
    r3_full = r3
    call check_approx(passed, "A .transform. C (3sym) vs 2sym version", r3_full, r_full, tol=1D-13)
    call check_approx(passed, "A .transform. C (2sym) keeps major symmetry", r_full%vals, &
                      transpose(r_full%vals), tol=1D-13)

contains

    function C_full_from(mat) result(t)
        !! Wraps a 6x6 Voigt matrix as `ten_3D4O2sym`.
        real(real64), intent(in) :: mat(6,6)
            !! Voigt matrix.
        type(ten_3D4O2sym) :: t
            !! Tensor.
        call t%init(mat)
    end function C_full_from
end subroutine test_transform4_3sym
