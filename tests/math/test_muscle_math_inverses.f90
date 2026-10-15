! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module test_inverses_helpers
    !! Dispatch over the matrix sizes of `muscle_math_inverses` and exact reference matrices.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_math_inverses, only : M33INV, M44INV, M55INV, M66INV, M99INV
    implicit none

    integer, parameter :: SIZES(5) = (/3, 4, 5, 6, 9/)
        !! Sizes with a dedicated inverse.

contains

    subroutine invert(n, a, x, ok, kappa)
        !! Calls the inverse routine of size n (`kappa` only for 3x3 and 4x4, else -1).
        integer, intent(in) :: n
            !! Matrix size (3, 4, 5, 6 or 9).
        real(real64), intent(in) :: a(n,n)
            !! Matrix to invert.
        real(real64), intent(out) :: x(n,n)
            !! Computed inverse.
        logical, intent(out) :: ok
            !! Success flag returned by the routine.
        real(real64), intent(out) :: kappa
            !! Condition number returned by M33INV/M44INV; -1 for the other sizes.
        kappa = -1.0D0
        select case (n)
        case (3); call M33INV(a, x, ok, kappa)
        case (4); call M44INV(a, x, ok, kappa)
        case (5); call M55INV(a, x, ok)
        case (6); call M66INV(a, x, ok)
        case (9); call M99INV(a, x, ok)
        end select
    end subroutine invert

    pure function hilbert(n) result(h)
        !! Hilbert matrix H_ij = 1 / (i + j - 1).
        integer, intent(in) :: n
            !! Size.
        real(real64) :: h(n,n)
            !! Hilbert matrix.
        integer :: i, j
        do j = 1, n
            do i = 1, n
                h(i,j) = 1.0D0 / dble(i + j - 1)
            end do
        end do
    end function hilbert

    pure function hilbert_inverse(n) result(h)
        !! Exact inverse of the Hilbert matrix (integer entries, exactly representable for n <= 9):
        !! (H^-1)_ij = (-1)^(i+j) (i+j-1) C(n+i-1, n-j) C(n+j-1, n-i) C(i+j-2, i-1)^2.
        integer, intent(in) :: n
            !! Size.
        real(real64) :: h(n,n)
            !! Exact inverse.
        integer :: i, j
        do j = 1, n
            do i = 1, n
                h(i,j) = (-1)**(i+j) * dble(i+j-1) * binom(n+i-1, n-j) * binom(n+j-1, n-i) &
                         * binom(i+j-2, i-1)**2
            end do
        end do
    end function hilbert_inverse

    pure real(real64) function binom(n, k)
        !! Binomial coefficient C(n, k) (exact in double precision for the sizes used here).
        integer, intent(in) :: n
            !! Upper index.
        integer, intent(in) :: k
            !! Lower index.
        integer :: m
        binom = 1.0D0
        do m = 1, k
            binom = binom * dble(n - k + m) / dble(m)
        end do
    end function binom

    pure real(real64) function norm_inf(a)
        !! Infinity norm (maximum absolute row sum).
        real(real64), intent(in) :: a(:,:)
            !! Matrix.
        norm_inf = maxval(sum(abs(a), dim=2))
    end function norm_inf

end module test_inverses_helpers


program test_muscle_math_inverses
    !! Tests of the small dense inverses of `muscle_math_inverses` (M33INV, M44INV, M55INV,
    !! M66INV, M99INV): exact cases, Hilbert matrices with their exact integer inverse (covers the
    !! adjugate and the LU-fallback paths of M33INV/M44INV), general matrices, singular input,
    !! NaN input and extreme scaling.
    use, intrinsic :: ieee_exceptions, only : ieee_set_flag, ieee_all
    use muscle_test_utils
    implicit none
    logical :: passed

    call test_exact_cases(passed)
    if (.not. passed) stop 1

    call test_hilbert(passed)
    if (.not. passed) stop 2

    call test_general(passed)
    if (.not. passed) stop 3

    call test_singular_nan(passed)
    if (.not. passed) stop 4

    call test_extreme_scaling(passed)
    if (.not. passed) stop 5

    ! The NaN and extreme-scaling cases raise IEEE invalid/overflow/underflow flags on purpose;
    ! clear them so that the runtime does not report them as unexpected at STOP.
    call ieee_set_flag(ieee_all, .false.)
    print *, "Passed!"
    stop 0
end program test_muscle_math_inverses


subroutine test_exact_cases(passed)
    !! Inverses that are exact in binary, compared exactly (tol = 0):
    !! - identity -> identity;
    !! - diagonal with powers of two -> reciprocals (also powers of two);
    !! - permutation matrix (needs pivoting / has zero leading pivots) -> its transpose.
    !! All operations involve only 0, +-1 and powers of two, so no rounding occurs.
    !! For 3x3/4x4 the returned kappa_inf must be exact as well (1, max/min diagonal, 1).
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_test_utils
    use test_inverses_helpers
    implicit none
    logical, intent(out) :: passed
    real(real64), allocatable :: a(:,:), x(:,:), expected(:,:)
    real(real64) :: kap
    integer :: s, n, i
    logical :: ok
    character(len=32) :: tag

    passed = .true.
    do s = 1, size(SIZES)
        n = SIZES(s)
        write(tag, '(i0, "x", i0)') n, n
        allocate(a(n,n), x(n,n), expected(n,n))

        a = 0.0D0; expected = 0.0D0
        do i = 1, n
            a(i,i) = 1.0D0; expected(i,i) = 1.0D0
        end do
        call invert(n, a, x, ok, kap)
        call check_true(passed, trim(tag) // " identity: ok", ok)
        call check_approx(passed, trim(tag) // " inverse of identity", x, expected, tol=0.0D0)
        if (n <= 4) call check_approx(passed, trim(tag) // " kappa(identity) = 1", kap, 1.0D0, tol=0.0D0)

        a = 0.0D0; expected = 0.0D0
        do i = 1, n
            a(i,i) = 2.0D0**(i - 3)
            expected(i,i) = 2.0D0**(3 - i)
        end do
        call invert(n, a, x, ok, kap)
        call check_approx(passed, trim(tag) // " inverse of diag(2^(i-3))", x, expected, tol=0.0D0)
        if (n <= 4) call check_approx(passed, trim(tag) // " kappa(diag) = max/min", kap, 2.0D0**(n-1), tol=0.0D0)

        ! Cyclic permutation: e_i -> e_(i+1); zero diagonal, so LU must pivot
        a = 0.0D0
        do i = 1, n
            a(1 + mod(i, n), i) = 1.0D0
        end do
        expected = transpose(a)
        call invert(n, a, x, ok, kap)
        call check_approx(passed, trim(tag) // " inverse of a permutation = transpose", x, expected, tol=0.0D0)

        deallocate(a, x, expected)
    end do
end subroutine test_exact_cases


subroutine test_hilbert(passed)
    !! Hilbert matrices against their exact integer inverse (n = 3, 4, 5, 6, 9).
    !! kappa_inf(H) grows from 748 (n = 3) to 1.1e12 (n = 9): H3 takes the adjugate path of M33INV
    !! (kappa <= 1e3), H4 (kappa 2.8e4) the LU fallback of M44INV.
    !! Tolerance: the relative error of a backward-stable inverse is O(kappa * eps); the benchmark of
    !! 2026-10-14 measured at most 0.22 kappa * eps, so tol = 10 kappa_inf * eps leaves a wide margin.
    !! For 3x3/4x4 the returned kappa_inf is compared with the exact one (from the exact inverse):
    !! it is computed from the numerical inverse, so it carries the same relative error.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_test_utils
    use test_inverses_helpers
    implicit none
    logical, intent(out) :: passed
    real(real64), allocatable :: h(:,:), x(:,:), hinv(:,:)
    real(real64) :: kap, kexact
    integer :: s, n
    logical :: ok
    character(len=32) :: tag

    passed = .true.
    do s = 1, size(SIZES)
        n = SIZES(s)
        write(tag, '("Hilbert ", i0, "x", i0)') n, n
        allocate(h(n,n), x(n,n), hinv(n,n))
        h = hilbert(n)
        hinv = hilbert_inverse(n)
        kexact = norm_inf(h) * norm_inf(hinv)
        call invert(n, h, x, ok, kap)
        call check_true(passed, trim(tag) // ": ok", ok)
        call check_approx(passed, trim(tag) // " inverse vs exact", x, hinv, &
                          tol=10.0D0 * kexact * epsilon(1.0D0))
        if (n <= 4) call check_approx(passed, trim(tag) // " kappa_inf vs exact", kap, kexact, &
                                      tol=10.0D0 * kexact * epsilon(1.0D0))
        deallocate(h, x, hinv)
    end do
end subroutine test_hilbert


subroutine test_general(passed)
    !! General non-symmetric, diagonally dominant matrices (A_ij = sin(i + 3j) + 5 delta_ij,
    !! kappa_inf ~ 2-4): A * inv(A) = inv(A) * A = I. The residual is a few ulp times kappa;
    !! tol 1e-13 (absolute on entries of order 1).
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_test_utils
    use test_inverses_helpers
    implicit none
    logical, intent(out) :: passed
    real(real64), allocatable :: a(:,:), x(:,:), id(:,:)
    real(real64) :: kap
    integer :: s, n, i, j
    logical :: ok
    character(len=32) :: tag

    passed = .true.
    do s = 1, size(SIZES)
        n = SIZES(s)
        write(tag, '(i0, "x", i0)') n, n
        allocate(a(n,n), x(n,n), id(n,n))
        id = 0.0D0
        do j = 1, n
            id(j,j) = 1.0D0
            do i = 1, n
                a(i,j) = sin(dble(i) + 3.0D0*dble(j)) + merge(5.0D0, 0.0D0, i == j)
            end do
        end do
        call invert(n, a, x, ok, kap)
        call check_true(passed, trim(tag) // " general: ok", ok)
        call check_approx(passed, trim(tag) // " A * inv(A) = I", matmul(a, x), id, tol=1.0D-13)
        call check_approx(passed, trim(tag) // " inv(A) * A = I", matmul(x, a), id, tol=1.0D-13)
        if (n <= 4) call check_approx(passed, trim(tag) // " kappa_inf = |A| |inv(A)|", kap, &
                                      norm_inf(a) * norm_inf(x), tol=1.0D-14)
        deallocate(a, x, id)
    end do
end subroutine test_general


subroutine test_singular_nan(passed)
    !! A singular matrix (rank 1, all entries equal) and a matrix containing a NaN must give
    !! ok = .false. and a zero inverse (kappa = huge for 3x3/4x4).
    use, intrinsic :: iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use muscle_test_utils
    use test_inverses_helpers
    implicit none
    logical, intent(out) :: passed
    real(real64), allocatable :: a(:,:), x(:,:), zero(:,:)
    real(real64) :: kap
    integer :: s, n
    logical :: ok
    character(len=32) :: tag

    passed = .true.
    do s = 1, size(SIZES)
        n = SIZES(s)
        write(tag, '(i0, "x", i0)') n, n
        allocate(a(n,n), x(n,n), zero(n,n))
        zero = 0.0D0

        a = 1.0D0
        call invert(n, a, x, ok, kap)
        call check_true(passed, trim(tag) // " singular: ok = .false.", .not. ok)
        call check_approx(passed, trim(tag) // " singular: inverse = 0", x, zero, tol=0.0D0)
        if (n <= 4) call check_approx(passed, trim(tag) // " singular: kappa = huge", kap, huge(1.0D0), &
                                      tol=0.0D0)

        a = 0.0D0
        a(1,1) = ieee_value(1.0D0, ieee_quiet_nan)
        call invert(n, a, x, ok, kap)
        call check_true(passed, trim(tag) // " NaN input: ok = .false.", .not. ok)

        deallocate(a, x, zero)
    end do
end subroutine test_singular_nan


subroutine test_extreme_scaling(passed)
    !! Well-conditioned matrices with extreme scale: c * B with c = 1e-100 and 1e+100 and B the
    !! general matrix of `test_general`. For 4x4 the adjugate determinant would underflow
    !! (1e-400) or overflow (1e+400), so M33INV/M44INV must fall back to LU instead of failing.
    !! Expected: inv(c B) = inv(B) / c. The scaling c is not a power of two, so c*B and inv(B)/c
    !! carry a few ulp of rounding; with kappa ~ 4, tol 1e-13 relative is a wide margin.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_test_utils
    use test_inverses_helpers
    implicit none
    logical, intent(out) :: passed
    real(real64), allocatable :: b(:,:), xb(:,:), x(:,:)
    real(real64) :: kap, c
    integer :: s, n, i, j, k
    logical :: ok
    character(len=48) :: tag

    passed = .true.
    do s = 1, size(SIZES)
        n = SIZES(s)
        allocate(b(n,n), xb(n,n), x(n,n))
        do j = 1, n
            do i = 1, n
                b(i,j) = sin(dble(i) + 3.0D0*dble(j)) + merge(5.0D0, 0.0D0, i == j)
            end do
        end do
        call invert(n, b, xb, ok, kap)
        do k = 1, 2
            c = merge(1.0D-100, 1.0D+100, k == 1)
            write(tag, '(i0, "x", i0, " scaled by ", es8.1)') n, n, c
            call invert(n, c*b, x, ok, kap)
            call check_true(passed, trim(tag) // ": ok", ok)
            call check_approx(passed, trim(tag) // ": inv(cB) = inv(B)/c", x, xb/c, tol=1.0D-13)
        end do
        deallocate(b, xb, x)
    end do
end subroutine test_extreme_scaling
