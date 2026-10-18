! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module test_dyadic_helpers
    !! Reference solution for the dyadic products, built from the full index definition.
    !!
    !! A symmetric second-order tensor is written by hand as a full 3x3 matrix, the fourth-order
    !! result is evaluated component by component (\( C_{ijkl} = A_{ij} B_{kl} \) and its
    !! variants), and every \( C_{ijkl} \) is written to its Voigt position \( (I, J) \) with
    !! I = V(i,j), J = V(k,l). The Voigt map V is written here independently of the library,
    !! so a wrong packing order or a wrong shear position in the library is detected.
    !!
    !! All inputs are small integers, so every product, every sum of two products and every
    !! halving of such a sum is exact in double precision: the checks use `tol = 0`.
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none

    integer, parameter :: V3(3,3) = reshape([1, 4, 6, &
                                             4, 2, 5, &
                                             6, 5, 3], [3,3])
        !! 3D Voigt map (xx, yy, zz, xy, yz, xz) of the index pair (i,j).
    integer, parameter :: V2(3,3) = reshape([1, 4, 0, &
                                             4, 2, 0, &
                                             0, 0, 3], [3,3])
        !! 2D Voigt map (xx, yy, zz, xy); 0 marks the out-of-plane shear pairs, which vanish.
    real(real64), parameter :: DELTA(3,3) = reshape([1D0, 0D0, 0D0, &
                                                     0D0, 1D0, 0D0, &
                                                     0D0, 0D0, 1D0], [3,3])
        !! Kronecker delta (second-order identity).

contains

    pure function full_3D(v) result(m)
        !! Full 3x3 matrix of a 3D symmetric tensor given in Voigt order (xx, yy, zz, xy, yz, xz).
        real(real64), intent(in) :: v(6)
            !! Voigt components.
        real(real64) :: m(3,3)
            !! Full symmetric matrix.
        m(1,:) = [v(1), v(4), v(6)]
        m(2,:) = [v(4), v(2), v(5)]
        m(3,:) = [v(6), v(5), v(3)]
    end function full_3D

    pure function full_2D(v) result(m)
        !! Full 3x3 matrix of a 2D symmetric tensor given in Voigt order (xx, yy, zz, xy).
        real(real64), intent(in) :: v(4)
            !! Voigt components.
        real(real64) :: m(3,3)
            !! Full symmetric matrix (zero out-of-plane shear).
        m(1,:) = [v(1), v(4), 0D0]
        m(2,:) = [v(4), v(2), 0D0]
        m(3,:) = [0D0,  0D0,  v(3)]
    end function full_2D

    pure function dyad_voigt(A, B, symmetrize, n, V) result(C)
        !! Voigt matrix of \( C_{ijkl} = A_{ij} B_{kl} \), or of
        !! \( \frac{1}{2}(A_{ij} B_{kl} + B_{ij} A_{kl}) \) if `symmetrize`.
        real(real64), intent(in) :: A(3,3)
            !! First full second-order tensor.
        real(real64), intent(in) :: B(3,3)
            !! Second full second-order tensor.
        logical, intent(in) :: symmetrize
            !! Return the symmetrized product.
        integer, intent(in) :: n
            !! Size of the Voigt matrix (6 in 3D, 4 in 2D).
        integer, intent(in) :: V(3,3)
            !! Voigt map; pairs mapped to 0 are not stored.
        real(real64) :: C(n,n)
            !! Voigt matrix of the product.
        integer :: i, j, k, l

        C = 0D0
        do l = 1, 3
            do k = 1, 3
                do j = 1, 3
                    do i = 1, 3
                        if (V(i,j) == 0 .or. V(k,l) == 0) cycle
                        if (symmetrize) then
                            C(V(i,j), V(k,l)) = 0.5D0*(A(i,j)*B(k,l) + B(i,j)*A(k,l))
                        else
                            C(V(i,j), V(k,l)) = A(i,j)*B(k,l)
                        end if
                    end do
                end do
            end do
        end do
    end function dyad_voigt

end module test_dyadic_helpers


program test_operator_dyadic
    !! Tests of the dyadic product operators, one section per overload
    !! (module `muscle_tensor_ops_dyadic`):
    !!
    !! - `.tdot.`: \( \mathbf{a} \otimes \mathbf{b} \), \( \mathbf{a} \otimes \mathbf{I} \),
    !!   \( \mathbf{I} \otimes \mathbf{a} \), in 3D and 2D (minor symmetries only);
    !! - `.tdotsym.`: \( \frac{1}{2}(\mathbf{a} \otimes \mathbf{b} + \mathbf{b} \otimes \mathbf{a}) \),
    !!   the unary \( \mathbf{a} \otimes \mathbf{a} \) and the products with \( \mathbf{I} \),
    !!   in 3D and 2D (fully symmetric).
    !!
    !! The reference is the full index definition (module `test_dyadic_helpers`). The 2D
    !! overloads are also checked against the 3D ones in `test_operator_2d_vs_3d`.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none

    logical :: passed

    print*, "Running tests for operators .tdot. and .tdotsym. ..."

    ! --- tdot_3D2Osym_3D2Osym ---
    call test_tdot_3D2Osym_3D2Osym_definition(passed)
    if (.not. passed) stop 1

    ! --- tdot_3D2Osym_I2O ---
    call test_tdot_3D2Osym_I2O_definition(passed)
    if (.not. passed) stop 2

    ! --- tdot_I2O_3D2Osym ---
    call test_tdot_I2O_3D2Osym_definition(passed)
    if (.not. passed) stop 3

    ! --- tdot_2D2Osym_2D2Osym ---
    call test_tdot_2D2Osym_2D2Osym_definition(passed)
    if (.not. passed) stop 4

    ! --- tdot_2D2Osym_I2O, tdot_I2O_2D2Osym ---
    call test_tdot_2D2Osym_I2O_definition(passed)
    if (.not. passed) stop 5

    ! --- tdotsym_3D2Osym_3D2Osym ---
    call test_tdotsym_3D2Osym_3D2Osym_definition(passed)
    if (.not. passed) stop 6
    call test_tdotsym_3D2Osym_3D2Osym_sparse(passed)
    if (.not. passed) stop 7

    ! --- tdotsym_3D2Osym (unary) ---
    call test_tdotsym_3D2Osym_unary(passed)
    if (.not. passed) stop 8

    ! --- tdotsym_3D2Osym_I2O, tdotsym_I2O_3D2Osym ---
    call test_tdotsym_3D2Osym_I2O_definition(passed)
    if (.not. passed) stop 9

    ! --- tdotsym_2D2Osym_2D2Osym ---
    call test_tdotsym_2D2Osym_2D2Osym_definition(passed)
    if (.not. passed) stop 10

    ! --- tdotsym_2D2Osym (unary) ---
    call test_tdotsym_2D2Osym_unary(passed)
    if (.not. passed) stop 11

    ! --- tdotsym_2D2Osym_I2O, tdotsym_I2O_2D2Osym ---
    call test_tdotsym_2D2Osym_I2O_definition(passed)
    if (.not. passed) stop 12

    print*, "All .tdot. / .tdotsym. operator tests passed successfully!"
    stop 0
end program test_operator_dyadic


! =========================================================================
! tdot_3D2Osym_3D2Osym
! =========================================================================

subroutine test_tdot_3D2Osym_3D2Osym_definition(passed)
    !! a = (1..6), b = (7..12) in Voigt (xx, yy, zz, xy, yz, xz):
    !! \( C_{ijkl} = a_{ij} b_{kl} \) from the full index definition. Exact products of
    !! integers; tol 0. The result has no major symmetry (a and b are not proportional).
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    use test_dyadic_helpers
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: a, b
    type(ten_3D4O2sym) :: res

    passed = .true.
    call a%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, yz=5D0, xz=6D0)
    call b%init(xx=7D0, yy=8D0, zz=9D0, xy=10D0, yz=11D0, xz=12D0)

    res = a .tdot. b
    call check_approx(passed, "a .tdot. b (3D)", res%vals, &
        dyad_voigt(full_3D(a%vals), full_3D(b%vals), .false., 6, V3), tol=0D0)
    call check_true(passed, "a .tdot. b (3D) has no major symmetry", &
        maxval(abs(res%vals - transpose(res%vals))) > 0D0)
end subroutine test_tdot_3D2Osym_3D2Osym_definition


! =========================================================================
! tdot_3D2Osym_I2O
! =========================================================================

subroutine test_tdot_3D2Osym_I2O_definition(passed)
    !! \( C_{ijkl} = a_{ij} \delta_{kl} \): the three normal columns hold a, the shear columns
    !! are zero. a = (1..6); exact copies; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    use test_dyadic_helpers
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: a
    type(iden_2O) :: I2
    type(ten_3D4O2sym) :: res

    passed = .true.
    call a%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, yz=5D0, xz=6D0)

    res = a .tdot. I2
    call check_approx(passed, "a .tdot. I (3D)", res%vals, &
        dyad_voigt(full_3D(a%vals), DELTA, .false., 6, V3), tol=0D0)
end subroutine test_tdot_3D2Osym_I2O_definition


! =========================================================================
! tdot_I2O_3D2Osym
! =========================================================================

subroutine test_tdot_I2O_3D2Osym_definition(passed)
    !! \( C_{ijkl} = \delta_{ij} a_{kl} \): the three normal rows hold a, i.e. the transpose of
    !! `a .tdot. I`. a = (1..6); exact copies; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    use test_dyadic_helpers
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: a
    type(iden_2O) :: I2
    type(ten_3D4O2sym) :: res, res_t

    passed = .true.
    call a%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, yz=5D0, xz=6D0)

    res = I2 .tdot. a
    call check_approx(passed, "I .tdot. a (3D)", res%vals, &
        dyad_voigt(DELTA, full_3D(a%vals), .false., 6, V3), tol=0D0)
    res_t = a .tdot. I2
    call check_approx(passed, "I .tdot. a = (a .tdot. I)^T (3D)", res%vals, &
        transpose(res_t%vals), tol=0D0)
end subroutine test_tdot_I2O_3D2Osym_definition


! =========================================================================
! tdot_2D2Osym_2D2Osym
! =========================================================================

subroutine test_tdot_2D2Osym_2D2Osym_definition(passed)
    !! a = (1, 2, 3, 4), b = (5, 6, 7, 8) in Voigt (xx, yy, zz, xy):
    !! \( C_{ijkl} = a_{ij} b_{kl} \), the 4x4 matrix \( C_{IJ} = a_I b_J \) (no major symmetry).
    !! Exact products of integers; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    use test_dyadic_helpers
    implicit none
    logical, intent(out) :: passed

    type(ten_2D2Osym) :: a, b
    type(ten_2D4O2sym) :: res

    passed = .true.
    call a%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0)
    call b%init(xx=5D0, yy=6D0, zz=7D0, xy=8D0)

    res = a .tdot. b
    call check_approx(passed, "a .tdot. b (2D)", res%vals, &
        dyad_voigt(full_2D(a%vals), full_2D(b%vals), .false., 4, V2), tol=0D0)
    call check_true(passed, "a .tdot. b (2D) has no major symmetry", &
        maxval(abs(res%vals - transpose(res%vals))) > 0D0)
end subroutine test_tdot_2D2Osym_2D2Osym_definition


! =========================================================================
! tdot_2D2Osym_I2O, tdot_I2O_2D2Osym
! =========================================================================

subroutine test_tdot_2D2Osym_I2O_definition(passed)
    !! \( a_{ij} \delta_{kl} \) fills the three normal columns with a = (1, 2, 3, 4);
    !! \( \delta_{ij} a_{kl} \) is its transpose. Exact copies; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    use test_dyadic_helpers
    implicit none
    logical, intent(out) :: passed

    type(ten_2D2Osym) :: a
    type(iden_2O) :: I2
    type(ten_2D4O2sym) :: res

    passed = .true.
    call a%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0)

    res = a .tdot. I2
    call check_approx(passed, "a .tdot. I (2D)", res%vals, &
        dyad_voigt(full_2D(a%vals), DELTA, .false., 4, V2), tol=0D0)
    res = I2 .tdot. a
    call check_approx(passed, "I .tdot. a (2D)", res%vals, &
        dyad_voigt(DELTA, full_2D(a%vals), .false., 4, V2), tol=0D0)
end subroutine test_tdot_2D2Osym_I2O_definition


! =========================================================================
! tdotsym_3D2Osym_3D2Osym
! =========================================================================

subroutine test_tdotsym_3D2Osym_3D2Osym_definition(passed)
    !! a = (1..6), b = (7..12): \( \frac{1}{2}(a_{ij} b_{kl} + b_{ij} a_{kl}) \) from the full
    !! index definition, the symmetric part of `a .tdot. b`, and commutativity.
    !! Halves of sums of integer products are exact; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    use test_dyadic_helpers
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: a, b
    type(ten_3D4O2sym) :: res, ab

    passed = .true.
    call a%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, yz=5D0, xz=6D0)
    call b%init(xx=7D0, yy=8D0, zz=9D0, xy=10D0, yz=11D0, xz=12D0)

    res = a .tdotsym. b
    call check_approx(passed, "a .tdotsym. b (3D)", res%vals, &
        dyad_voigt(full_3D(a%vals), full_3D(b%vals), .true., 6, V3), tol=0D0)
    ab = a .tdot. b
    call check_approx(passed, "a .tdotsym. b = sym(a .tdot. b) (3D)", res%vals, &
        0.5D0*(ab%vals + transpose(ab%vals)), tol=0D0)
    call check_approx(passed, "a .tdotsym. b = b .tdotsym. a (3D)", a .tdotsym. b, &
        b .tdotsym. a, tol=0D0)
end subroutine test_tdotsym_3D2Osym_3D2Osym_definition


subroutine test_tdotsym_3D2Osym_3D2Osym_sparse(passed)
    !! a = diag(1, 2, 0), b with only b_12 = 3 (computed by hand): the only nonzero components
    !! are \( C_{1112} = \frac{1}{2}(a_{11} b_{12} + b_{11} a_{12}) = 1.5 \) and
    !! \( C_{2212} = 3 \), stored at positions 16 (11-12) and 13 (22-12). Exact; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: a, b
    type(ten_3D4O3sym) :: expected

    passed = .true.
    call a%init(xx=1D0, yy=2D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call b%init(xx=0D0, yy=0D0, zz=0D0, xy=3D0, yz=0D0, xz=0D0)

    expected%vals = 0D0
    expected%vals(16) = 1.5D0
    expected%vals(13) = 3.0D0
    call check_approx(passed, "a .tdotsym. b (3D, sparse)", a .tdotsym. b, expected, tol=0D0)
end subroutine test_tdotsym_3D2Osym_3D2Osym_sparse


! =========================================================================
! tdotsym_3D2Osym (unary)
! =========================================================================

subroutine test_tdotsym_3D2Osym_unary(passed)
    !! a = (1..6): `.tdotsym. a` and `a .tdotsym. a` give \( a_{ij} a_{kl} \). The 21 stored
    !! components (11, 22, 33, 44, 55, 66, 12, 23, 34, 45, 56, 13, 24, 35, 46, 14, 25, 36,
    !! 15, 26, 16) are the products a_I a_J written below by hand, and also the full index
    !! definition. Exact; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    use test_dyadic_helpers
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: a
    type(ten_3D4O3sym) :: expected
    type(ten_3D4O2sym) :: res

    passed = .true.
    call a%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, yz=5D0, xz=6D0)

    call expected%init([1D0, 4D0, 9D0, 16D0, 25D0, 36D0, &
                        2D0, 6D0, 12D0, 20D0, 30D0,       &
                        3D0, 8D0, 15D0, 24D0,             &
                        4D0, 10D0, 18D0,                  &
                        5D0, 12D0,                        &
                        6D0])
    call check_approx(passed, ".tdotsym. a (3D)", .tdotsym. a, expected, tol=0D0)
    call check_approx(passed, "a .tdotsym. a (3D)", a .tdotsym. a, expected, tol=0D0)

    res = .tdotsym. a
    call check_approx(passed, ".tdotsym. a (3D), full index definition", res%vals, &
        dyad_voigt(full_3D(a%vals), full_3D(a%vals), .false., 6, V3), tol=0D0)
end subroutine test_tdotsym_3D2Osym_unary


! =========================================================================
! tdotsym_3D2Osym_I2O, tdotsym_I2O_3D2Osym
! =========================================================================

subroutine test_tdotsym_3D2Osym_I2O_definition(passed)
    !! a = (1..6): \( \frac{1}{2}(a_{ij} \delta_{kl} + \delta_{ij} a_{kl}) \) from the full index
    !! definition: normal-normal block (a_I + a_J)/2, normal-shear block a_J/2, shear-shear
    !! block zero. Both operand orders give the same tensor. Exact; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    use test_dyadic_helpers
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: a
    type(iden_2O) :: I2
    type(ten_3D4O2sym) :: res

    passed = .true.
    call a%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, yz=5D0, xz=6D0)

    res = a .tdotsym. I2
    call check_approx(passed, "a .tdotsym. I (3D)", res%vals, &
        dyad_voigt(full_3D(a%vals), DELTA, .true., 6, V3), tol=0D0)
    res = I2 .tdotsym. a
    call check_approx(passed, "I .tdotsym. a (3D)", res%vals, &
        dyad_voigt(full_3D(a%vals), DELTA, .true., 6, V3), tol=0D0)
end subroutine test_tdotsym_3D2Osym_I2O_definition


! =========================================================================
! tdotsym_2D2Osym_2D2Osym
! =========================================================================

subroutine test_tdotsym_2D2Osym_2D2Osym_definition(passed)
    !! a = (1, 2, 3, 4), b = (5, 6, 7, 8): \( \frac{1}{2}(a_{ij} b_{kl} + b_{ij} a_{kl}) \) from
    !! the full index definition, the symmetric part of `a .tdot. b`, and commutativity.
    !! Exact; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    use test_dyadic_helpers
    implicit none
    logical, intent(out) :: passed

    type(ten_2D2Osym) :: a, b
    type(ten_2D4O2sym) :: res, ab

    passed = .true.
    call a%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0)
    call b%init(xx=5D0, yy=6D0, zz=7D0, xy=8D0)

    res = a .tdotsym. b
    call check_approx(passed, "a .tdotsym. b (2D)", res%vals, &
        dyad_voigt(full_2D(a%vals), full_2D(b%vals), .true., 4, V2), tol=0D0)
    ab = a .tdot. b
    call check_approx(passed, "a .tdotsym. b = sym(a .tdot. b) (2D)", res%vals, &
        0.5D0*(ab%vals + transpose(ab%vals)), tol=0D0)
    call check_approx(passed, "a .tdotsym. b = b .tdotsym. a (2D)", a .tdotsym. b, &
        b .tdotsym. a, tol=0D0)
end subroutine test_tdotsym_2D2Osym_2D2Osym_definition


! =========================================================================
! tdotsym_2D2Osym (unary)
! =========================================================================

subroutine test_tdotsym_2D2Osym_unary(passed)
    !! a = (1, 2, 3, 4): `.tdotsym. a` and `a .tdotsym. a` give \( a_{ij} a_{kl} \). Stored
    !! components (11, 22, 33, 44, 12, 23, 34, 13, 24, 14) = products a_I a_J written by hand.
    !! Exact; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed

    type(ten_2D2Osym) :: a
    type(ten_2D4O3sym) :: expected

    passed = .true.
    call a%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0)

    call expected%init([1D0, 4D0, 9D0, 16D0, 2D0, 6D0, 12D0, 3D0, 8D0, 4D0])
    call check_approx(passed, ".tdotsym. a (2D)", .tdotsym. a, expected, tol=0D0)
    call check_approx(passed, "a .tdotsym. a (2D)", a .tdotsym. a, expected, tol=0D0)
end subroutine test_tdotsym_2D2Osym_unary


! =========================================================================
! tdotsym_2D2Osym_I2O, tdotsym_I2O_2D2Osym
! =========================================================================

subroutine test_tdotsym_2D2Osym_I2O_definition(passed)
    !! a = (1, 2, 3, 4): \( \frac{1}{2}(a_{ij} \delta_{kl} + \delta_{ij} a_{kl}) \). Stored
    !! components (11, 22, 33, 44, 12, 23, 34, 13, 24, 14) written by hand: normal-normal
    !! (a_I + a_J)/2, normal-shear a_xy/2 = 2, shear-shear 0; also the full index definition.
    !! Both operand orders give the same tensor. Exact; tol 0.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    use test_dyadic_helpers
    implicit none
    logical, intent(out) :: passed

    type(ten_2D2Osym) :: a
    type(iden_2O) :: I2
    type(ten_2D4O3sym) :: expected
    type(ten_2D4O2sym) :: res

    passed = .true.
    call a%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0)

    call expected%init([1D0, 2D0, 3D0, 0D0, 1.5D0, 2.5D0, 2D0, 2D0, 2D0, 2D0])
    call check_approx(passed, "a .tdotsym. I (2D)", a .tdotsym. I2, expected, tol=0D0)
    call check_approx(passed, "I .tdotsym. a (2D)", I2 .tdotsym. a, expected, tol=0D0)

    res = a .tdotsym. I2
    call check_approx(passed, "a .tdotsym. I (2D), full index definition", res%vals, &
        dyad_voigt(full_2D(a%vals), DELTA, .true., 4, V2), tol=0D0)
end subroutine test_tdotsym_2D2Osym_I2O_definition
