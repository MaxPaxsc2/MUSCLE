! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module test_contraction_single_helpers
    !! Reference solution for the single contraction \( C_{ij} = A_{ik} B_{kj} \): the operands
    !! are written as full 3x3 matrices and multiplied with the intrinsic `matmul`.
    !!
    !! The storage maps of the tensor types are written here independently of the library, so
    !! a wrong component position in a product is detected:
    !! - `ten_3D2Osym`: (xx, yy, zz, xy, yz, xz);
    !! - `ten_3D2O`:    column-major (xx, yx, zx, xy, yy, zy, xz, yz, zz);
    !! - `ten_2D2Osym`: (xx, yy, zz, xy);
    !! - `ten_2D2O`:    (xx, yx, xy, yy, zz).
    !!
    !! All inputs are small integers, so every product and every sum of three products is
    !! exact in double precision: the checks use `tol = 0`.
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none

contains

    pure function full_3Dsym(v) result(m)
        !! Full matrix of a `ten_3D2Osym`.
        real(real64), intent(in) :: v(6)
            !! Stored components.
        real(real64) :: m(3,3)
            !! Full symmetric matrix.
        m(1,:) = [v(1), v(4), v(6)]
        m(2,:) = [v(4), v(2), v(5)]
        m(3,:) = [v(6), v(5), v(3)]
    end function full_3Dsym

    pure function full_3D(v) result(m)
        !! Full matrix of a `ten_3D2O`.
        real(real64), intent(in) :: v(9)
            !! Stored components.
        real(real64) :: m(3,3)
            !! Full matrix.
        m(1,:) = [v(1), v(4), v(7)]
        m(2,:) = [v(2), v(5), v(8)]
        m(3,:) = [v(3), v(6), v(9)]
    end function full_3D

    pure function full_2Dsym(v) result(m)
        !! Full matrix of a `ten_2D2Osym` (zero out-of-plane shear).
        real(real64), intent(in) :: v(4)
            !! Stored components.
        real(real64) :: m(3,3)
            !! Full symmetric matrix.
        m(1,:) = [v(1), v(4), 0D0]
        m(2,:) = [v(4), v(2), 0D0]
        m(3,:) = [0D0,  0D0,  v(3)]
    end function full_2Dsym

    pure function full_2D(v) result(m)
        !! Full matrix of a `ten_2D2O` (zero out-of-plane shear).
        real(real64), intent(in) :: v(5)
            !! Stored components.
        real(real64) :: m(3,3)
            !! Full matrix.
        m(1,:) = [v(1), v(3), 0D0]
        m(2,:) = [v(2), v(4), 0D0]
        m(3,:) = [0D0,  0D0,  v(5)]
    end function full_2D

end module test_contraction_single_helpers


program test_operator_contraction_single
    !! Tests of the single contraction operator `*` between second-order tensors and
    !! identities, one section per overload (module `muscle_tensor_ops_contraction_single`):
    !!
    !! - `dot_*`: \( \mathbf{A} \cdot \mathbf{B} \) for every pair of 3D (and of 2D) types,
    !!   symmetric or general; the result is always a general tensor;
    !! - `mul_*`: \( \mathbf{I} \cdot \mathbf{A} = \mathbf{A} \cdot \mathbf{I} = \mathbf{A} \) and
    !!   \( c\mathbf{I} \cdot \mathbf{A} = \mathbf{A} \cdot c\mathbf{I} = c\mathbf{A} \).
    !!
    !! The reference is `matmul` on full matrices (module `test_contraction_single_helpers`).
    !! The 2D overloads are also checked against the 3D ones in `test_operator_2d_vs_3d`.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none

    logical :: passed

    print*, "Running tests for operator * (single contraction) ..."

    ! --- dot_3D2Osym_3D2Osym ---
    call test_dot_3D2Osym_3D2Osym_definition(passed)
    if (.not. passed) stop 1
    call test_dot_3D2Osym_3D2Osym_cases(passed)
    if (.not. passed) stop 2

    ! --- dot_3D2O_3D2Osym ---
    call test_dot_3D2O_3D2Osym_definition(passed)
    if (.not. passed) stop 3

    ! --- dot_3D2Osym_3D2O ---
    call test_dot_3D2Osym_3D2O_definition(passed)
    if (.not. passed) stop 4

    ! --- dot_3D2O_3D2O ---
    call test_dot_3D2O_3D2O_definition(passed)
    if (.not. passed) stop 5

    ! --- mul_I2O_3D2O, mul_3D2O_I2O, mul_I2O_3D2Osym, mul_3D2Osym_I2O ---
    call test_mul_I2O_3D(passed)
    if (.not. passed) stop 6

    ! --- mul_I2OS_3D2O, mul_3D2O_I2OS, mul_I2OS_3D2Osym, mul_3D2Osym_I2OS ---
    call test_mul_I2OS_3D(passed)
    if (.not. passed) stop 7

    ! --- dot_2D2Osym_2D2Osym ---
    call test_dot_2D2Osym_2D2Osym_definition(passed)
    if (.not. passed) stop 8

    ! --- dot_2D2O_2D2Osym ---
    call test_dot_2D2O_2D2Osym_definition(passed)
    if (.not. passed) stop 9

    ! --- dot_2D2Osym_2D2O ---
    call test_dot_2D2Osym_2D2O_definition(passed)
    if (.not. passed) stop 10

    ! --- dot_2D2O_2D2O ---
    call test_dot_2D2O_2D2O_definition(passed)
    if (.not. passed) stop 11
    call test_dot_2D2O_2D2O_swap(passed)
    if (.not. passed) stop 12

    ! --- mul_I2O_2D2O, mul_2D2O_I2O, mul_I2O_2D2Osym, mul_2D2Osym_I2O ---
    call test_mul_I2O_2D(passed)
    if (.not. passed) stop 13

    ! --- mul_I2OS_2D2O, mul_2D2O_I2OS, mul_I2OS_2D2Osym, mul_2D2Osym_I2OS ---
    call test_mul_I2OS_2D(passed)
    if (.not. passed) stop 14

    print*, "All * (single contraction) operator tests passed successfully!"
    stop 0
end program test_operator_contraction_single


! =========================================================================
! dot_3D2Osym_3D2Osym
! =========================================================================

subroutine test_dot_3D2Osym_3D2Osym_definition(passed)
    !! a = (1..6), b = (7..12) in Voigt (xx, yy, zz, xy, yz, xz): `a * b` = matmul of the full
    !! matrices. For symmetric operands \( (\mathbf{a}\mathbf{b})^T = \mathbf{b}\mathbf{a} \),
    !! and here a and b do not commute, so a * b and b * a must differ. Exact; tol 0.
    !!
    !! Commuting case: c = 2a + I = (3, 5, 7, 8, 10, 12) commutes with a, since
    !! \( \mathbf{a}(2\mathbf{a} + \mathbf{I}) = 2\mathbf{a}^2 + \mathbf{a}
    !!   = (2\mathbf{a} + \mathbf{I})\mathbf{a} \), so a * c and c * a must be equal.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    use test_contraction_single_helpers
    implicit none
    logical, intent(out) :: passed

    real(real64), parameter :: TOL_EPS = 10D0*epsilon(1D0)
        !! ~2.2e-15 relative. With these integer inputs the results are exact, but a component
        !! of a.b is a sum of 3 products, whose rounding error is at most ~3 eps times
        !! sum |a_ik b_kj| (Higham, Accuracy and Stability of Numerical Algorithms, sec. 3.1);
        !! 10 eps covers that bound with margin. Differences below it are rounding.

    type(ten_3D2Osym) :: a, b, c
    type(ten_3D2O) :: ab, ba

    passed = .true.
    call a%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, yz=5D0, xz=6D0)
    call b%init(xx=7D0, yy=8D0, zz=9D0, xy=10D0, yz=11D0, xz=12D0)

    ab = a * b
    ba = b * a
    call check_approx(passed, "a * b (3D sym)", full_3D(ab%vals), &
        matmul(full_3Dsym(a%vals), full_3Dsym(b%vals)), tol=0D0)
    call check_approx(passed, "(a * b)^T = b * a (3D sym)", transpose(full_3D(ab%vals)), &
        full_3D(ba%vals), tol=0D0)
    call check_true(passed, "a * b /= b * a (3D sym)", &
        maxval(abs(ab%vals - ba%vals)) > &
        TOL_EPS*max(maxval(abs(ab%vals)), maxval(abs(ba%vals)), 1D-30))

    ! c = 2a + I commutes with a: a * c = c * a
    call c%init(xx=3D0, yy=5D0, zz=7D0, xy=8D0, yz=10D0, xz=12D0)
    ab = a * c
    ba = c * a
    call check_approx(passed, "a * (2a+I) = (2a+I) * a (3D sym)", ab, ba, tol=TOL_EPS)
end subroutine test_dot_3D2Osym_3D2Osym_definition


subroutine test_dot_3D2Osym_3D2Osym_cases(passed)
    !! Hand-computed cases (column-major result xx, yx, zx, xy, yy, zy, xz, yz, zz):
    !! - 0 * 0 = 0; I * I = I; (2I) * (2I) = 4I;
    !! - all-ones tensor squared: every entry 1+1+1 = 3;
    !! - A = (1..6) squared: (53, 42, 44, 42, 45, 49, 44, 49, 70);
    !! - A * 0 = 0 * A = 0;
    !! - A = (1..6), B = all ones: A*B has rows of A summed, (11, 11, 14) per column;
    !!   B*A = (A*B)^T.
    !! Exact integer arithmetic; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: a, b, zero
    type(ten_3D2O) :: expected

    passed = .true.
    call zero%init([0D0, 0D0, 0D0, 0D0, 0D0, 0D0])
    call expected%init([0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0])
    call check_approx(passed, "0 * 0", zero * zero, expected, tol=0D0)

    call a%init([1D0, 1D0, 1D0, 0D0, 0D0, 0D0])
    call expected%init([1D0, 0D0, 0D0, 0D0, 1D0, 0D0, 0D0, 0D0, 1D0])
    call check_approx(passed, "I * I", a * a, expected, tol=0D0)

    call a%init([2D0, 2D0, 2D0, 0D0, 0D0, 0D0])
    call expected%init([4D0, 0D0, 0D0, 0D0, 4D0, 0D0, 0D0, 0D0, 4D0])
    call check_approx(passed, "2I * 2I", a * a, expected, tol=0D0)

    call b%init([1D0, 1D0, 1D0, 1D0, 1D0, 1D0])
    call expected%init([3D0, 3D0, 3D0, 3D0, 3D0, 3D0, 3D0, 3D0, 3D0])
    call check_approx(passed, "ones * ones", b * b, expected, tol=0D0)

    call a%init([1D0, 2D0, 3D0, 4D0, 5D0, 6D0])
    call expected%init([53D0, 42D0, 44D0, 42D0, 45D0, 49D0, 44D0, 49D0, 70D0])
    call check_approx(passed, "A * A", a * a, expected, tol=0D0)

    call expected%init([0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0])
    call check_approx(passed, "A * 0", a * zero, expected, tol=0D0)
    call check_approx(passed, "0 * A", zero * a, expected, tol=0D0)

    call expected%init([11D0, 11D0, 14D0, 11D0, 11D0, 14D0, 11D0, 11D0, 14D0])
    call check_approx(passed, "A * ones", a * b, expected, tol=0D0)
    call expected%init([11D0, 11D0, 11D0, 11D0, 11D0, 11D0, 14D0, 14D0, 14D0])
    call check_approx(passed, "ones * A", b * a, expected, tol=0D0)
end subroutine test_dot_3D2Osym_3D2Osym_cases


! =========================================================================
! dot_3D2O_3D2Osym
! =========================================================================

subroutine test_dot_3D2O_3D2Osym_definition(passed)
    !! A general with all nine entries distinct (1..9), b = (10..15) symmetric:
    !! `A * b` = matmul of the full matrices. Exact; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    use test_contraction_single_helpers
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2O) :: a, res
    type(ten_3D2Osym) :: b

    passed = .true.
    call a%init(xx=1D0, xy=2D0, xz=3D0, yx=4D0, yy=5D0, yz=6D0, zx=7D0, zy=8D0, zz=9D0)
    call b%init(xx=10D0, yy=11D0, zz=12D0, xy=13D0, yz=14D0, xz=15D0)

    res = a * b
    call check_approx(passed, "A * b (3D gen * sym)", full_3D(res%vals), &
        matmul(reshape([1D0, 4D0, 7D0, 2D0, 5D0, 8D0, 3D0, 6D0, 9D0], [3,3]), &
               full_3Dsym(b%vals)), tol=0D0)
end subroutine test_dot_3D2O_3D2Osym_definition


! =========================================================================
! dot_3D2Osym_3D2O
! =========================================================================

subroutine test_dot_3D2Osym_3D2O_definition(passed)
    !! a = (10..15) symmetric, B general with entries 1..9: `a * B` = matmul of the full
    !! matrices. Exact; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    use test_contraction_single_helpers
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: a
    type(ten_3D2O) :: b, res

    passed = .true.
    call a%init(xx=10D0, yy=11D0, zz=12D0, xy=13D0, yz=14D0, xz=15D0)
    call b%init(xx=1D0, xy=2D0, xz=3D0, yx=4D0, yy=5D0, yz=6D0, zx=7D0, zy=8D0, zz=9D0)

    res = a * b
    call check_approx(passed, "a * B (3D sym * gen)", full_3D(res%vals), &
        matmul(full_3Dsym(a%vals), &
               reshape([1D0, 4D0, 7D0, 2D0, 5D0, 8D0, 3D0, 6D0, 9D0], [3,3])), tol=0D0)
end subroutine test_dot_3D2Osym_3D2O_definition


! =========================================================================
! dot_3D2O_3D2O
! =========================================================================

subroutine test_dot_3D2O_3D2O_definition(passed)
    !! A with entries 1..9 and B with entries 10..18 (rows written below):
    !! `A * B` = matmul of the full matrices, and `B * A` differs. Exact; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    use test_contraction_single_helpers
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2O) :: a, b, ab, ba
    real(real64) :: ma(3,3), mb(3,3)

    passed = .true.
    call a%init(xx=1D0, xy=2D0, xz=3D0, yx=4D0, yy=5D0, yz=6D0, zx=7D0, zy=8D0, zz=9D0)
    call b%init(xx=10D0, xy=11D0, xz=12D0, yx=13D0, yy=14D0, yz=15D0, zx=16D0, zy=17D0, zz=18D0)
    ma = reshape([1D0, 4D0, 7D0, 2D0, 5D0, 8D0, 3D0, 6D0, 9D0], [3,3])
    mb = reshape([10D0, 13D0, 16D0, 11D0, 14D0, 17D0, 12D0, 15D0, 18D0], [3,3])

    ab = a * b
    ba = b * a
    call check_approx(passed, "A * B (3D gen)", full_3D(ab%vals), matmul(ma, mb), tol=0D0)
    call check_approx(passed, "B * A (3D gen)", full_3D(ba%vals), matmul(mb, ma), tol=0D0)
end subroutine test_dot_3D2O_3D2O_definition


! =========================================================================
! mul_I2O_3D2O, mul_3D2O_I2O, mul_I2O_3D2Osym, mul_3D2Osym_I2O
! =========================================================================

subroutine test_mul_I2O_3D(passed)
    !! I * A = A * I = A for a general and a symmetric 3D tensor. Exact copies; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2O) :: g
    type(ten_3D2Osym) :: s
    type(iden_2O) :: I2

    passed = .true.
    call g%init(xx=1D0, xy=2D0, xz=3D0, yx=4D0, yy=5D0, yz=6D0, zx=7D0, zy=8D0, zz=9D0)
    call s%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, yz=5D0, xz=6D0)

    call check_approx(passed, "I * A (3D gen)", I2 * g, g, tol=0D0)
    call check_approx(passed, "A * I (3D gen)", g * I2, g, tol=0D0)
    call check_approx(passed, "I * a (3D sym)", I2 * s, s, tol=0D0)
    call check_approx(passed, "a * I (3D sym)", s * I2, s, tol=0D0)
end subroutine test_mul_I2O_3D


! =========================================================================
! mul_I2OS_3D2O, mul_3D2O_I2OS, mul_I2OS_3D2Osym, mul_3D2Osym_I2OS
! =========================================================================

subroutine test_mul_I2OS_3D(passed)
    !! (cI) * A = A * (cI) = c A with c = 3, written by hand:
    !! general A with entries 1..9 -> (3, 12, 21, 6, 15, 24, 9, 18, 27) column-major;
    !! symmetric a = (1..6) -> (3, 6, 9, 12, 15, 18). Exact; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2O) :: g, expected_g
    type(ten_3D2Osym) :: s, expected_s
    type(iden_2OS) :: cI

    passed = .true.
    call cI%init(3D0)
    call g%init(xx=1D0, xy=2D0, xz=3D0, yx=4D0, yy=5D0, yz=6D0, zx=7D0, zy=8D0, zz=9D0)
    call s%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, yz=5D0, xz=6D0)
    call expected_g%init([3D0, 12D0, 21D0, 6D0, 15D0, 24D0, 9D0, 18D0, 27D0])
    call expected_s%init(xx=3D0, yy=6D0, zz=9D0, xy=12D0, yz=15D0, xz=18D0)

    call check_approx(passed, "cI * A (3D gen)", cI * g, expected_g, tol=0D0)
    call check_approx(passed, "A * cI (3D gen)", g * cI, expected_g, tol=0D0)
    call check_approx(passed, "cI * a (3D sym)", cI * s, expected_s, tol=0D0)
    call check_approx(passed, "a * cI (3D sym)", s * cI, expected_s, tol=0D0)
end subroutine test_mul_I2OS_3D


! =========================================================================
! dot_2D2Osym_2D2Osym
! =========================================================================

subroutine test_dot_2D2Osym_2D2Osym_definition(passed)
    !! a = (1, 2, 3, 4), b = (5, 6, 7, 8) in Voigt (xx, yy, zz, xy): `a * b` = matmul of the
    !! full matrices (the in-plane block and zz = a_zz b_zz), and (a * b)^T = b * a. Exact; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    use test_contraction_single_helpers
    implicit none
    logical, intent(out) :: passed

    type(ten_2D2Osym) :: a, b
    type(ten_2D2O) :: ab, ba

    passed = .true.
    call a%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0)
    call b%init(xx=5D0, yy=6D0, zz=7D0, xy=8D0)

    ab = a * b
    ba = b * a
    call check_approx(passed, "a * b (2D sym)", full_2D(ab%vals), &
        matmul(full_2Dsym(a%vals), full_2Dsym(b%vals)), tol=0D0)
    call check_approx(passed, "(a * b)^T = b * a (2D sym)", transpose(full_2D(ab%vals)), &
        full_2D(ba%vals), tol=0D0)
end subroutine test_dot_2D2Osym_2D2Osym_definition


! =========================================================================
! dot_2D2O_2D2Osym
! =========================================================================

subroutine test_dot_2D2O_2D2Osym_definition(passed)
    !! A = [[1, 2], [3, 4]], A_zz = 5 (general) and b = (6, 7, 8, 9): `A * b` = matmul of the
    !! full matrices. Exact; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    use test_contraction_single_helpers
    implicit none
    logical, intent(out) :: passed

    type(ten_2D2O) :: a, res
    type(ten_2D2Osym) :: b

    passed = .true.
    call a%init(xx=1D0, xy=2D0, yx=3D0, yy=4D0, zz=5D0)
    call b%init(xx=6D0, yy=7D0, zz=8D0, xy=9D0)

    res = a * b
    call check_approx(passed, "A * b (2D gen * sym)", full_2D(res%vals), &
        matmul(reshape([1D0, 3D0, 0D0, 2D0, 4D0, 0D0, 0D0, 0D0, 5D0], [3,3]), &
               full_2Dsym(b%vals)), tol=0D0)
end subroutine test_dot_2D2O_2D2Osym_definition


! =========================================================================
! dot_2D2Osym_2D2O
! =========================================================================

subroutine test_dot_2D2Osym_2D2O_definition(passed)
    !! a = (6, 7, 8, 9) symmetric and B = [[1, 2], [3, 4]], B_zz = 5: `a * B` = matmul of the
    !! full matrices. Exact; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    use test_contraction_single_helpers
    implicit none
    logical, intent(out) :: passed

    type(ten_2D2Osym) :: a
    type(ten_2D2O) :: b, res

    passed = .true.
    call a%init(xx=6D0, yy=7D0, zz=8D0, xy=9D0)
    call b%init(xx=1D0, xy=2D0, yx=3D0, yy=4D0, zz=5D0)

    res = a * b
    call check_approx(passed, "a * B (2D sym * gen)", full_2D(res%vals), &
        matmul(full_2Dsym(a%vals), &
               reshape([1D0, 3D0, 0D0, 2D0, 4D0, 0D0, 0D0, 0D0, 5D0], [3,3])), tol=0D0)
end subroutine test_dot_2D2Osym_2D2O_definition


! =========================================================================
! dot_2D2O_2D2O
! =========================================================================

subroutine test_dot_2D2O_2D2O_definition(passed)
    !! A = [[1, 2], [3, 4]], A_zz = 5 and B = [[6, 7], [8, 9]], B_zz = 10: `A * B` and `B * A`
    !! = matmul of the full matrices. Exact; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    use test_contraction_single_helpers
    implicit none
    logical, intent(out) :: passed

    type(ten_2D2O) :: a, b, ab, ba
    real(real64) :: ma(3,3), mb(3,3)

    passed = .true.
    call a%init(xx=1D0, xy=2D0, yx=3D0, yy=4D0, zz=5D0)
    call b%init(xx=6D0, xy=7D0, yx=8D0, yy=9D0, zz=10D0)
    ma = reshape([1D0, 3D0, 0D0, 2D0, 4D0, 0D0, 0D0, 0D0, 5D0], [3,3])
    mb = reshape([6D0, 8D0, 0D0, 7D0, 9D0, 0D0, 0D0, 0D0, 10D0], [3,3])

    ab = a * b
    ba = b * a
    call check_approx(passed, "A * B (2D gen)", full_2D(ab%vals), matmul(ma, mb), tol=0D0)
    call check_approx(passed, "B * A (2D gen)", full_2D(ba%vals), matmul(mb, ma), tol=0D0)
end subroutine test_dot_2D2O_2D2O_definition


subroutine test_dot_2D2O_2D2O_swap(passed)
    !! Hand-computed: [[1,2],[3,4]] . [[0,1],[1,0]] = [[2,1],[4,3]] (right-multiplying by the
    !! swap matrix exchanges the columns), zz: 5*2 = 10; and
    !! det(A.B) = det(A) det(B) = (-10)(-2) = 20.
    !! Exact integer arithmetic (the determinants are sums of integer products); tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed

    type(ten_2D2O) :: a, b, ab, expected

    passed = .true.
    call a%init(xx=1D0, xy=2D0, yx=3D0, yy=4D0, zz=5D0)
    call b%init(xx=0D0, xy=1D0, yx=1D0, yy=0D0, zz=2D0)
    call expected%init(xx=2D0, xy=1D0, yx=4D0, yy=3D0, zz=10D0)

    ab = a * b
    call check_approx(passed, "A * swap (2D gen)", ab, expected, tol=0D0)
    call check_approx(passed, "det(A*B) = det(A) det(B)", ab%det(), a%det()*b%det(), tol=0D0)
end subroutine test_dot_2D2O_2D2O_swap


! =========================================================================
! mul_I2O_2D2O, mul_2D2O_I2O, mul_I2O_2D2Osym, mul_2D2Osym_I2O
! =========================================================================

subroutine test_mul_I2O_2D(passed)
    !! I * A = A * I = A for a general and a symmetric 2D tensor. Exact copies; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed

    type(ten_2D2O) :: g
    type(ten_2D2Osym) :: s
    type(iden_2O) :: I2

    passed = .true.
    call g%init(xx=1D0, xy=2D0, yx=3D0, yy=4D0, zz=5D0)
    call s%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0)

    call check_approx(passed, "I * A (2D gen)", I2 * g, g, tol=0D0)
    call check_approx(passed, "A * I (2D gen)", g * I2, g, tol=0D0)
    call check_approx(passed, "I * a (2D sym)", I2 * s, s, tol=0D0)
    call check_approx(passed, "a * I (2D sym)", s * I2, s, tol=0D0)
end subroutine test_mul_I2O_2D


! =========================================================================
! mul_I2OS_2D2O, mul_2D2O_I2OS, mul_I2OS_2D2Osym, mul_2D2Osym_I2OS
! =========================================================================

subroutine test_mul_I2OS_2D(passed)
    !! (cI) * A = A * (cI) = c A with c = 3, written by hand:
    !! general A = (xx, yx, xy, yy, zz) = (1, 3, 2, 4, 5) -> (3, 9, 6, 12, 15);
    !! symmetric a = (1, 2, 3, 4) -> (3, 6, 9, 12). Exact; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed

    type(ten_2D2O) :: g, expected_g
    type(ten_2D2Osym) :: s, expected_s
    type(iden_2OS) :: cI

    passed = .true.
    call cI%init(3D0)
    call g%init(xx=1D0, xy=2D0, yx=3D0, yy=4D0, zz=5D0)
    call s%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0)
    call expected_g%init(xx=3D0, xy=6D0, yx=9D0, yy=12D0, zz=15D0)
    call expected_s%init(xx=3D0, yy=6D0, zz=9D0, xy=12D0)

    call check_approx(passed, "cI * A (2D gen)", cI * g, expected_g, tol=0D0)
    call check_approx(passed, "A * cI (2D gen)", g * cI, expected_g, tol=0D0)
    call check_approx(passed, "cI * a (2D sym)", cI * s, expected_s, tol=0D0)
    call check_approx(passed, "a * cI (2D sym)", s * cI, expected_s, tol=0D0)
end subroutine test_mul_I2OS_2D
