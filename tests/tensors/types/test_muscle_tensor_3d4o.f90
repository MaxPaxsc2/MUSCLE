! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_tensor_3d4o
    !! Unit tests for `ten_3D4O` (general 3D fourth-order tensor, 9x9 storage with
    !! row \( 3(j-1)+i \) and column \( 3(l-1)+k \)), used only through `muscle_tensors`
    !! (so the program also checks that the type is exported).
    use muscle_tensors
    implicit none
    logical :: passed

    call test_3D4O_init_get_set(passed)
    if (.not. passed) STOP 1

    call test_3D4O_arithmetic_norm(passed)
    if (.not. passed) STOP 2

    call test_3D4O_approx(passed)
    if (.not. passed) STOP 3

    call test_3D4O_convert_2sym(passed)
    if (.not. passed) STOP 4

    call test_3D4O_inverse(passed)
    if (.not. passed) STOP 5

    print *, "Passed!"
    STOP 0
end program test_muscle_tensor_3d4o


subroutine test_3D4O_init_get_set(passed)
    !! init from a 3x3x3x3 array with C4_ijkl = 1000 i + 100 j + 10 k + l: every one of the 81
    !! components must be found at row 3(j-1)+i, column 3(l-1)+k, and returned by get(i,j,k,l).
    !! init from a 9x9 matrix must give the same tensor. set/get on a single component.
    !! All values are exact copies of small integers, so the comparisons are exact (tol = 0).
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D4O) :: a, b
    real(real64) :: C4(3,3,3,3), expected(9,9), got_get(9,9)
    integer :: i, j, k, l

    passed = .true.
    do l = 1, 3; do k = 1, 3; do j = 1, 3; do i = 1, 3
        C4(i,j,k,l) = 1000D0*i + 100D0*j + 10D0*k + l
        expected(3*(j-1)+i, 3*(l-1)+k) = C4(i,j,k,l)
    end do; end do; end do; end do

    call a%init(C4)
    call check_approx(passed, "init(C4): vals", a%vals, expected, tol=0D0)
    do l = 1, 3; do k = 1, 3; do j = 1, 3; do i = 1, 3
        got_get(3*(j-1)+i, 3*(l-1)+k) = a%get(i, j, k, l)
    end do; end do; end do; end do
    call check_approx(passed, "init(C4): get(i,j,k,l)", got_get, expected, tol=0D0)

    call b%init(expected)
    call check_approx(passed, "init(9x9) = init(C4)", b, a, tol=0D0)

    b = 0D0
    call b%set(2, 3, 1, 2, -7.5D0)
    call check_approx(passed, "set(2,3,1,2) -> vals(8,4)", b%vals(8, 4), -7.5D0, tol=0D0)
    call check_approx(passed, "get(2,3,1,2)", b%get(2, 3, 1, 2), -7.5D0, tol=0D0)
    call check_approx(passed, "set touches a single component", sum(abs(b%vals)), 7.5D0, tol=0D0)
end subroutine test_3D4O_init_get_set


subroutine test_3D4O_arithmetic_norm(passed)
    !! +, -, unary -, scalar *, / and scalar assignment on integer-valued entries scaled by
    !! powers of two (exact in binary, tol = 0); L1 norm sum |A_IJ| of integers (exact).
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D4O) :: a, b, r
    real(real64) :: m(9,9)
    integer :: i, j

    passed = .true.
    do j = 1, 9
        do i = 1, 9
            m(i,j) = 10D0*i - j
        end do
    end do
    call a%init(m)
    call b%init(2D0*m)

    r = a + b;   call check_approx(passed, "a + b", r%vals, 3D0*m, tol=0D0)
    r = b - a;   call check_approx(passed, "b - a", r%vals, m, tol=0D0)
    r = -a;      call check_approx(passed, "-a", r%vals, -m, tol=0D0)
    r = 2D0*a;   call check_approx(passed, "2*a", r%vals, 2D0*m, tol=0D0)
    r = a*2D0;   call check_approx(passed, "a*2", r%vals, 2D0*m, tol=0D0)
    r = a/2D0;   call check_approx(passed, "a/2", r%vals, 0.5D0*m, tol=0D0)
    r = 1.5D0;   call check_approx(passed, "r = 1.5 (scalar assignment)", r%vals, spread(spread(1.5D0, 1, 9), 1, 9), tol=0D0)
    call check_approx(passed, "a%norm() = sum |m|", a%norm(), sum(abs(m)), tol=0D0)
end subroutine test_3D4O_arithmetic_norm


subroutine test_3D4O_approx(passed)
    !! Relative tolerance of is_approx / .approx.: a 1e-9 relative perturbation of one component
    !! is rejected with the default 1e-12 and accepted with tol = 1e-8.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D4O) :: a, b

    passed = .true.
    a = 3D0
    b = a
    call check_true(passed, ".approx. of two identical tensors", a .approx. b)
    b%vals(4, 7) = b%vals(4, 7) * (1D0 + 1D-9)
    call check_true(passed, ".approx. rejects a 1e-9 relative perturbation", .not. (a .approx. b))
    call check_true(passed, "is_approx(tol=1e-8) accepts a 1e-9 relative perturbation", &
                    a%is_approx(b, tol=1D-8))
end subroutine test_3D4O_approx


subroutine test_3D4O_convert_2sym(passed)
    !! convert_2sym samples rows/columns (xx, yy, zz, xy, yz, xz) = 9D indices (1, 5, 9, 4, 8, 7)
    !! of the 9x9 matrix into the 6x6 Voigt matrix (exact copies, tol = 0).
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D4O) :: a
    type(ten_3D4O2sym) :: c
    real(real64) :: m(9,9), expected(6,6)
    integer, parameter :: map6(6) = (/1, 5, 9, 4, 8, 7/)
    integer :: i, j

    passed = .true.
    do j = 1, 9
        do i = 1, 9
            m(i,j) = 10D0*i + j
        end do
    end do
    call a%init(m)
    do j = 1, 6
        do i = 1, 6
            expected(i,j) = m(map6(i), map6(j))
        end do
    end do
    c = a%convert_2sym()
    call check_approx(passed, "a%convert_2sym()", c%vals, expected, tol=0D0)
end subroutine test_3D4O_convert_2sym


subroutine test_3D4O_inverse(passed)
    !! .inv.: A * inv(A) = inv(A) * A = I_9 (products of the 9x9 representations) for a
    !! diagonally dominant A (diagonal 10, off-diagonal |a| <= 1, so cond(A) ~ 2): the LU
    !! round-off is ~1e-15, tol 1e-12. A singular A (all entries equal, rank 1) returns quiet
    !! NaN in every component.
    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D4O) :: a, ainv
    real(real64) :: m(9,9), id9(9,9)
    integer :: i, j

    passed = .true.
    id9 = 0D0
    do j = 1, 9
        id9(j,j) = 1D0
        do i = 1, 9
            m(i,j) = sin(1D0*i + 3D0*j) + merge(10D0, 0D0, i == j)
        end do
    end do
    call a%init(m)
    ainv = .inv. a
    call check_approx(passed, "A * inv(A) = I", matmul(a%vals, ainv%vals), id9, tol=1D-12)
    call check_approx(passed, "inv(A) * A = I", matmul(ainv%vals, a%vals), id9, tol=1D-12)

    a = 1D0
    ainv = .inv. a
    call check_true(passed, ".inv. (singular) = NaN", all(ieee_is_nan(ainv%vals)))
end subroutine test_3D4O_inverse
