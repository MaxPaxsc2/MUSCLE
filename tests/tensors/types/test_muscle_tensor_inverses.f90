! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_tensor_inverses
    !! Inverse API of every tensor type (`call inverse(a, a_inv, ok, cond)`, `.inv.`, `%inv()`):
    !!
    !! - a regular case with an inverse that is exact in binary (entries 0, +-1 and powers of two),
    !!   so the three entry points must agree exactly and match the expected tensor (tol = 0), and
    !!   `cond` must equal the exact \( \kappa_\infty \) of the matrix that is inverted;
    !! - a singular case: `ok = .false.`, quiet NaN in every component, `cond = huge`;
    !! - 2nd order symmetric types: a well-conditioned tensor of small scale (1e-11), which the
    !!   former absolute test |det| < 1e-30 wrongly rejected;
    !! - `ten_3D4O2sym`: isotropic compliance (analytic) and a tensor without major symmetry.
    use muscle_test_utils
    implicit none
    logical :: passed

    call test_inverse_2D2O(passed)
    if (.not. passed) stop 1
    call test_inverse_2D2Osym(passed)
    if (.not. passed) stop 2
    call test_inverse_3D2O(passed)
    if (.not. passed) stop 3
    call test_inverse_3D2Osym(passed)
    if (.not. passed) stop 4
    call test_inverse_2D4O2sym(passed)
    if (.not. passed) stop 5
    call test_inverse_2D4O3sym(passed)
    if (.not. passed) stop 6
    call test_inverse_2D4O(passed)
    if (.not. passed) stop 7
    call test_inverse_3D4O2sym(passed)
    if (.not. passed) stop 8
    call test_inverse_3D4O3sym(passed)
    if (.not. passed) stop 9
    call test_inverse_3D4O(passed)
    if (.not. passed) stop 10

    print *, "Passed!"
    stop 0
end program test_muscle_tensor_inverses


! =========================================================================
! 2nd order
! =========================================================================

subroutine test_inverse_2D2O(passed)
    !! F = [[2, 1], [0, 4]], F_zz = 0.5 -> inverse [[1/2, -1/8], [0, 1/4]], 2 (exact).
    !! Inverted as the 3x3 [[2,1,0],[0,4,0],[0,0,0.5]]: ||A||_inf = 4, ||A^-1||_inf = 2 -> cond 8.
    use, intrinsic :: iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2O) :: a, x, expected
    real(real64) :: cond
    logical :: ok

    passed = .true.
    call a%init(xx=2D0, xy=1D0, yx=0D0, yy=4D0, zz=0.5D0)
    call expected%init(xx=0.5D0, xy=-0.125D0, yx=0D0, yy=0.25D0, zz=2D0)
    call inverse(a, x, ok, cond)
    call check_true(passed, "2D2O: ok", ok)
    call check_approx(passed, "2D2O: inverse", x, expected, tol=0D0)
    call check_approx(passed, "2D2O: .inv.", .inv. a, expected, tol=0D0)
    call check_approx(passed, "2D2O: %inv()", a%inv(), expected, tol=0D0)
    call check_approx(passed, "2D2O: cond", cond, 8D0, tol=0D0)

    call a%init(xx=1D0, xy=2D0, yx=2D0, yy=4D0, zz=1D0)      ! singular in-plane block
    call inverse(a, x, ok, cond)
    call check_true(passed, "2D2O singular: ok = .false.", .not. ok)
    call check_true(passed, "2D2O singular: NaN", all(ieee_is_nan(x%vals)))
    call check_approx(passed, "2D2O singular: cond = huge", cond, huge(1D0), tol=0D0)
end subroutine test_inverse_2D2O


subroutine test_inverse_2D2Osym(passed)
    !! A = [[4, 2], [2, 2]], A_zz = 0.5 -> inverse [[1/2, -1/2], [-1/2, 1]], 2 (exact);
    !! 3x3 [[4,2,0],[2,2,0],[0,0,0.5]]: ||A||_inf = 6, ||A^-1||_inf = 2 -> cond 12.
    !! Small scale: 1e-11 * diag(1, 2, 4) (cond 4) -> 1e11 * diag(1, 1/2, 1/4); the scaling by
    !! 1e-11 is not exact in binary, so tol 1e-14 relative (a few ulp).
    use, intrinsic :: iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2Osym) :: a, x, expected
    real(real64) :: cond
    logical :: ok

    passed = .true.
    call a%init(xx=4D0, yy=2D0, zz=0.5D0, xy=2D0)
    call expected%init(xx=0.5D0, yy=1D0, zz=2D0, xy=-0.5D0)
    call inverse(a, x, ok, cond)
    call check_true(passed, "2D2Osym: ok", ok)
    call check_approx(passed, "2D2Osym: inverse", x, expected, tol=0D0)
    call check_approx(passed, "2D2Osym: .inv.", .inv. a, expected, tol=0D0)
    call check_approx(passed, "2D2Osym: %inv()", a%inv(), expected, tol=0D0)
    call check_approx(passed, "2D2Osym: cond", cond, 12D0, tol=0D0)

    call a%init(xx=1D-11, yy=2D-11, zz=4D-11, xy=0D0)
    call expected%init(xx=1D11, yy=0.5D11, zz=0.25D11, xy=0D0)
    call inverse(a, x, ok, cond)
    call check_true(passed, "2D2Osym small scale: ok", ok)
    call check_approx(passed, "2D2Osym small scale: inverse", x, expected, tol=1D-14)

    call a%init(xx=1D0, yy=1D0, zz=1D0, xy=1D0)              ! singular in-plane block
    call inverse(a, x, ok, cond)
    call check_true(passed, "2D2Osym singular: ok = .false.", .not. ok)
    call check_true(passed, "2D2Osym singular: NaN", all(ieee_is_nan(x%vals)))
    call check_approx(passed, "2D2Osym singular: cond = huge", cond, huge(1D0), tol=0D0)
end subroutine test_inverse_2D2Osym


subroutine test_inverse_3D2O(passed)
    !! Upper triangular A = [[2,1,0],[0,4,2],[0,0,1/2]] -> inverse [[1/2,-1/8,1/2],[0,1/4,-1],[0,0,2]]
    !! (exact); ||A||_inf = 6, ||A^-1||_inf = 2 -> cond 12.
    use, intrinsic :: iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D2O) :: a, x, expected
    real(real64) :: cond
    logical :: ok

    passed = .true.
    call a%init(xx=2D0, xy=1D0, xz=0D0, yx=0D0, yy=4D0, yz=2D0, zx=0D0, zy=0D0, zz=0.5D0)
    call expected%init(xx=0.5D0, xy=-0.125D0, xz=0.5D0, yx=0D0, yy=0.25D0, yz=-1D0, &
                       zx=0D0, zy=0D0, zz=2D0)
    call inverse(a, x, ok, cond)
    call check_true(passed, "3D2O: ok", ok)
    call check_approx(passed, "3D2O: inverse", x, expected, tol=0D0)
    call check_approx(passed, "3D2O: .inv.", .inv. a, expected, tol=0D0)
    call check_approx(passed, "3D2O: %inv()", a%inv(), expected, tol=0D0)
    call check_approx(passed, "3D2O: cond", cond, 12D0, tol=0D0)

    a = 1D0                                                   ! rank 1
    call inverse(a, x, ok, cond)
    call check_true(passed, "3D2O singular: ok = .false.", .not. ok)
    call check_true(passed, "3D2O singular: NaN", all(ieee_is_nan(x%vals)))
    call check_approx(passed, "3D2O singular: cond = huge", cond, huge(1D0), tol=0D0)
end subroutine test_inverse_3D2O


subroutine test_inverse_3D2Osym(passed)
    !! A = [[2,0,0],[0,4,2],[0,2,2]] -> inverse [[1/2,0,0],[0,1/2,-1/2],[0,-1/2,1]] (exact);
    !! ||A||_inf = 6, ||A^-1||_inf = 3/2 -> cond 9.
    !! Small scale: 1e-11 * diag(1, 2, 4) -> 1e11 * diag(1, 1/2, 1/4), tol 1e-14 relative.
    use, intrinsic :: iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D2Osym) :: a, x, expected
    real(real64) :: cond
    logical :: ok

    passed = .true.
    call a%init(xx=2D0, yy=4D0, zz=2D0, xy=0D0, yz=2D0, xz=0D0)
    call expected%init(xx=0.5D0, yy=0.5D0, zz=1D0, xy=0D0, yz=-0.5D0, xz=0D0)
    call inverse(a, x, ok, cond)
    call check_true(passed, "3D2Osym: ok", ok)
    call check_approx(passed, "3D2Osym: inverse", x, expected, tol=0D0)
    call check_approx(passed, "3D2Osym: .inv.", .inv. a, expected, tol=0D0)
    call check_approx(passed, "3D2Osym: %inv()", a%inv(), expected, tol=0D0)
    call check_approx(passed, "3D2Osym: cond", cond, 9D0, tol=0D0)

    call a%init(xx=1D-11, yy=2D-11, zz=4D-11, xy=0D0, yz=0D0, xz=0D0)
    call expected%init(xx=1D11, yy=0.5D11, zz=0.25D11, xy=0D0, yz=0D0, xz=0D0)
    call inverse(a, x, ok, cond)
    call check_true(passed, "3D2Osym small scale: ok", ok)
    call check_approx(passed, "3D2Osym small scale: inverse", x, expected, tol=1D-14)

    a = 1D0                                                   ! rank 1
    call inverse(a, x, ok, cond)
    call check_true(passed, "3D2Osym singular: ok = .false.", .not. ok)
    call check_true(passed, "3D2Osym singular: NaN", all(ieee_is_nan(x%vals)))
    call check_approx(passed, "3D2Osym singular: cond = huge", cond, huge(1D0), tol=0D0)
end subroutine test_inverse_3D2Osym


! =========================================================================
! 4th order: diagonal Voigt matrices with powers of two (exact inverses);
! shear entries of the result are additionally scaled by 1/4 (tensorial shear)
! =========================================================================

subroutine test_inverse_2D4O2sym(passed)
    !! Voigt diag(1, 2, 4, 8) -> diag(1, 1/2, 1/4, 1/8 / 4); cond = 8.
    use, intrinsic :: iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D4O2sym) :: a, x, expected
    real(real64) :: m(4,4), cond
    logical :: ok
    integer :: i

    passed = .true.
    m = 0D0
    do i = 1, 4
        m(i,i) = 2D0**(i-1)
    end do
    call a%init(m)
    m = 0D0
    m(1,1) = 1D0; m(2,2) = 0.5D0; m(3,3) = 0.25D0; m(4,4) = 0.125D0/4D0
    call expected%init(m)
    call inverse(a, x, ok, cond)
    call check_true(passed, "2D4O2sym: ok", ok)
    call check_approx(passed, "2D4O2sym: inverse", x, expected, tol=0D0)
    call check_approx(passed, "2D4O2sym: .inv.", .inv. a, expected, tol=0D0)
    call check_approx(passed, "2D4O2sym: %inv()", a%inv(), expected, tol=0D0)
    call check_approx(passed, "2D4O2sym: cond", cond, 8D0, tol=0D0)

    a%vals = 1D0
    call inverse(a, x, ok, cond)
    call check_true(passed, "2D4O2sym singular: ok = .false.", .not. ok)
    call check_true(passed, "2D4O2sym singular: NaN", all(ieee_is_nan(x%vals)))
    call check_approx(passed, "2D4O2sym singular: cond = huge", cond, huge(1D0), tol=0D0)
end subroutine test_inverse_2D4O2sym


subroutine test_inverse_2D4O3sym(passed)
    !! Voigt diag(1, 2, 4, 8) -> diag(1, 1/2, 1/4, 1/8 / 4); cond = 8.
    use, intrinsic :: iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D4O3sym) :: a, x, expected
    real(real64) :: cond
    logical :: ok

    passed = .true.
    call a%init((/1D0, 2D0, 4D0, 8D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    call expected%init((/1D0, 0.5D0, 0.25D0, 0.125D0/4D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    call inverse(a, x, ok, cond)
    call check_true(passed, "2D4O3sym: ok", ok)
    call check_approx(passed, "2D4O3sym: inverse", x, expected, tol=0D0)
    call check_approx(passed, "2D4O3sym: .inv.", .inv. a, expected, tol=0D0)
    call check_approx(passed, "2D4O3sym: %inv()", a%inv(), expected, tol=0D0)
    call check_approx(passed, "2D4O3sym: cond", cond, 8D0, tol=0D0)

    a = 1D0
    call inverse(a, x, ok, cond)
    call check_true(passed, "2D4O3sym singular: ok = .false.", .not. ok)
    call check_true(passed, "2D4O3sym singular: NaN", all(ieee_is_nan(x%vals)))
    call check_approx(passed, "2D4O3sym singular: cond = huge", cond, huge(1D0), tol=0D0)
end subroutine test_inverse_2D4O3sym


subroutine test_inverse_2D4O(passed)
    !! 5x5 diag(1, 2, 4, 8, 16) -> diag(1, 1/2, 1/4, 1/8, 1/16) (no shear scaling); cond = 16.
    use, intrinsic :: iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D4O) :: a, x, expected
    real(real64) :: cond
    logical :: ok
    integer :: i

    passed = .true.
    a = 0D0
    expected = 0D0
    do i = 1, 5
        a%vals(i,i) = 2D0**(i-1)
        expected%vals(i,i) = 2D0**(1-i)
    end do
    call inverse(a, x, ok, cond)
    call check_true(passed, "2D4O: ok", ok)
    call check_approx(passed, "2D4O: inverse", x, expected, tol=0D0)
    call check_approx(passed, "2D4O: .inv.", .inv. a, expected, tol=0D0)
    call check_approx(passed, "2D4O: %inv()", a%inv(), expected, tol=0D0)
    call check_approx(passed, "2D4O: cond", cond, 16D0, tol=0D0)

    a = 1D0
    call inverse(a, x, ok, cond)
    call check_true(passed, "2D4O singular: ok = .false.", .not. ok)
    call check_true(passed, "2D4O singular: NaN", all(ieee_is_nan(x%vals)))
    call check_approx(passed, "2D4O singular: cond = huge", cond, huge(1D0), tol=0D0)
end subroutine test_inverse_2D4O


subroutine test_inverse_3D4O2sym(passed)
    !! - Voigt diag(1, ..., 32) -> normal 1, 1/2, 1/4 and shear (1/8, 1/16, 1/32) / 4; cond = 32.
    !! - Isotropic stiffness lambda I(x)I + 2 mu I^S (E = 210, nu = 0.3) -> isotropic compliance
    !!   S_1111 = 1/E, S_1122 = -nu/E, S_1212 = (1+nu)/(2E) (tensorial). One 6x6 LU of a matrix
    !!   with cond ~ 10: error ~1e-15 relative, tol 1e-12.
    !! - No major symmetry: with W = diag(1,1,1,2,2,2), the tensorial-shear inverse satisfies
    !!   S W C W = C W S W = I (6x6). Diagonally dominant C (cond ~ 2): tol 1e-13.
    use, intrinsic :: iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D4O2sym) :: a, x, expected
    real(real64) :: m(6,6), w(6,6), id6(6,6), cond, lambda, mu
    real(real64), parameter :: E = 210D0, NU = 0.3D0
    logical :: ok
    integer :: i, j

    passed = .true.
    m = 0D0
    do i = 1, 6
        m(i,i) = 2D0**(i-1)
    end do
    call a%init(m)
    m = 0D0
    m(1,1) = 1D0; m(2,2) = 0.5D0; m(3,3) = 0.25D0
    m(4,4) = 0.125D0/4D0; m(5,5) = 0.0625D0/4D0; m(6,6) = 0.03125D0/4D0
    call expected%init(m)
    call inverse(a, x, ok, cond)
    call check_true(passed, "3D4O2sym: ok", ok)
    call check_approx(passed, "3D4O2sym: inverse", x, expected, tol=0D0)
    call check_approx(passed, "3D4O2sym: .inv.", .inv. a, expected, tol=0D0)
    call check_approx(passed, "3D4O2sym: %inv()", a%inv(), expected, tol=0D0)
    call check_approx(passed, "3D4O2sym: cond", cond, 32D0, tol=0D0)

    ! Isotropic compliance
    lambda = E*NU/((1D0 + NU)*(1D0 - 2D0*NU))
    mu = E/(2D0*(1D0 + NU))
    m = 0D0
    m(1:3,1:3) = lambda
    do i = 1, 3
        m(i,i) = lambda + 2D0*mu
        m(i+3,i+3) = mu
    end do
    call a%init(m)
    m = 0D0
    m(1:3,1:3) = -NU/E
    do i = 1, 3
        m(i,i) = 1D0/E
        m(i+3,i+3) = (1D0 + NU)/(2D0*E)
    end do
    call expected%init(m)
    call check_approx(passed, "3D4O2sym: isotropic compliance", .inv. a, expected, tol=1D-12)

    ! Without major symmetry
    w = 0D0; id6 = 0D0
    do i = 1, 6
        w(i,i) = merge(1D0, 2D0, i <= 3)
        id6(i,i) = 1D0
    end do
    do j = 1, 6
        do i = 1, 6
            m(i,j) = sin(1D0*i + 3D0*j) + merge(8D0, 0D0, i == j)
        end do
    end do
    call a%init(m)
    x = .inv. a
    call check_approx(passed, "3D4O2sym (no major symmetry): S W C W = I", &
                      matmul(x%vals, matmul(w, matmul(a%vals, w))), id6, tol=1D-13)
    call check_approx(passed, "3D4O2sym (no major symmetry): C W S W = I", &
                      matmul(a%vals, matmul(w, matmul(x%vals, w))), id6, tol=1D-13)

    a%vals = 1D0
    call inverse(a, x, ok, cond)
    call check_true(passed, "3D4O2sym singular: ok = .false.", .not. ok)
    call check_true(passed, "3D4O2sym singular: NaN", all(ieee_is_nan(x%vals)))
    call check_approx(passed, "3D4O2sym singular: cond = huge", cond, huge(1D0), tol=0D0)
end subroutine test_inverse_3D4O2sym


subroutine test_inverse_3D4O3sym(passed)
    !! Voigt diag(1, ..., 32) -> normal 1, 1/2, 1/4 and shear (1/8, 1/16, 1/32) / 4; cond = 32.
    use, intrinsic :: iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D4O3sym) :: a, x, expected
    real(real64) :: v(21), cond
    logical :: ok

    passed = .true.
    v = 0D0
    v(1:6) = (/1D0, 2D0, 4D0, 8D0, 16D0, 32D0/)
    call a%init(v)
    v(1:6) = (/1D0, 0.5D0, 0.25D0, 0.125D0/4D0, 0.0625D0/4D0, 0.03125D0/4D0/)
    call expected%init(v)
    call inverse(a, x, ok, cond)
    call check_true(passed, "3D4O3sym: ok", ok)
    call check_approx(passed, "3D4O3sym: inverse", x, expected, tol=0D0)
    call check_approx(passed, "3D4O3sym: .inv.", .inv. a, expected, tol=0D0)
    call check_approx(passed, "3D4O3sym: %inv()", a%inv(), expected, tol=0D0)
    call check_approx(passed, "3D4O3sym: cond", cond, 32D0, tol=0D0)

    a = 1D0
    call inverse(a, x, ok, cond)
    call check_true(passed, "3D4O3sym singular: ok = .false.", .not. ok)
    call check_true(passed, "3D4O3sym singular: NaN", all(ieee_is_nan(x%vals)))
    call check_approx(passed, "3D4O3sym singular: cond = huge", cond, huge(1D0), tol=0D0)
end subroutine test_inverse_3D4O3sym


subroutine test_inverse_3D4O(passed)
    !! 9x9 diag(2^0, ..., 2^8) -> diag(2^0, ..., 2^-8) (no shear scaling); cond = 256.
    use, intrinsic :: iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_3D4O) :: a, x, expected
    real(real64) :: cond
    logical :: ok
    integer :: i

    passed = .true.
    a = 0D0
    expected = 0D0
    do i = 1, 9
        a%vals(i,i) = 2D0**(i-1)
        expected%vals(i,i) = 2D0**(1-i)
    end do
    call inverse(a, x, ok, cond)
    call check_true(passed, "3D4O: ok", ok)
    call check_approx(passed, "3D4O: inverse", x, expected, tol=0D0)
    call check_approx(passed, "3D4O: .inv.", .inv. a, expected, tol=0D0)
    call check_approx(passed, "3D4O: %inv()", a%inv(), expected, tol=0D0)
    call check_approx(passed, "3D4O: cond", cond, 256D0, tol=0D0)

    a = 1D0
    call inverse(a, x, ok, cond)
    call check_true(passed, "3D4O singular: ok = .false.", .not. ok)
    call check_true(passed, "3D4O singular: NaN", all(ieee_is_nan(x%vals)))
    call check_approx(passed, "3D4O singular: cond = huge", cond, huge(1D0), tol=0D0)
end subroutine test_inverse_3D4O
