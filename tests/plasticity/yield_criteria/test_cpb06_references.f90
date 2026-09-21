! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module test_cpb06_references
    ! Independent references for the regression levels 5-7 (no spectral decomposition):
    !  - closed-form J2 (value, gradient, full Hessian),
    !  - CPB06 with k = 0 and even exponent a: Phi = tr(Sigma**a), a polynomial invariant
    !    that stays smooth when transformed eigenvalues repeat.
    !
    ! Independence from the implementation
    ! ------------------------------------
    ! The criterion under test computes its value and derivatives through an
    ! eigen-decomposition, divided differences and a Daleckii-Krein assembly. These
    ! references use none of that: no call to eigen_sym3, no divided difference, no use
    ! of the yield module. A sign error in the eigen-solver or a wrong limit in a divided
    ! difference therefore cannot cancel between the two sides.
    !
    ! The independence is of the algorithm, not of everything. Both sides use the tensor
    ! representation and contraction operators of muscle_tensors, the same transformation
    ! conventions, and the same formula for gamma and the normalisation constant B. An
    ! error in that formula would appear on both sides and pass unnoticed.
    !
    ! The polynomial identity
    ! -----------------------
    ! For k = 0 the yield function is Phi = sum |lambda_i|**a, and for an even integer
    ! a the absolute value is redundant: Phi = sum lambda_i**a = tr(Sigma**a), a
    ! polynomial invariant computable by repeated matrix multiplication with no
    ! eigenvalue extracted.
    !
    ! The restriction to k = 0 and even a limits coverage, but the polynomial route is
    ! smooth where the spectral one is most fragile: repeated eigenvalues, which force
    ! the divided differences into their limits, are an ordinary point for
    ! tr(Sigma**a).
    !
    ! Matrix derivatives
    ! ------------------
    ! Two identities drive the derivative references:
    !   d/dX tr(X**a) = a X**(a-1)
    !   d/dX (X**n) : E = sum_{r=0}^{n-1} X**r E X**(n-1-r)
    ! The sum in the second is not an artifact: matrices do not commute, so E has to be
    ! inserted in every position along the product instead of collecting into a single
    ! power.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    private
    public :: action_matrix, rel_err, j2_reference, tracepower_reference

contains

    function action_matrix(h) result(act)
        ! Column j = H : E_j, where E_j has vals e_j (tensorial shear)
        !
        ! Flattens a fourth-order tensor into the 6x6 matrix of its action, so that a
        ! Hessian stored in the 21-component format and one assembled column by column
        ! can be compared entry by entry without depending on the storage convention.
        type(ten_3D4O3sym), intent(in) :: h
        real(real64) :: act(6,6)
        type(ten_3D2Osym) :: e
        integer :: j
        do j = 1, 6
            e = 0.0D0
            e%vals(j) = 1.0D0
            e = h .ddot. e
            act(:, j) = e%vals
        end do
    end function action_matrix

    function rel_err(got, expected) result(res)
        ! max|got - expected| / max|expected|; huge() if got is not finite
        !
        ! Returning huge() rather than NaN keeps a non-finite result from passing a
        ! threshold comparison by accident: any test of the form err < tol then fails,
        ! whereas a NaN comparison would be false for both < and > and could be read
        ! either way depending on how the check is written.
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        real(real64), intent(in) :: got(:), expected(:)
        real(real64) :: res
        if (.not. all(ieee_is_finite(got))) then
            res = huge(1.0D0)
        else
            res = maxval(abs(got - expected))/max(maxval(abs(expected)), tiny(1.0D0))
        end if
    end function rel_err

    subroutine j2_reference(s, q, g, h21)
        ! q = sqrt(3/2 d:d); g = 3/2 d/q;
        ! H = 3/(2q) (I_sym - I(x)I/3) - 9/(4q^3) d(x)d, in ten_3D4O3sym storage order
        !
        ! The isotropic baseline: von Mises in closed form. CPB06 reproduces it when C
        ! is the identity, k = 0 and a = 2, which is the configuration the calling tests
        ! use; the agreement is not a property of every exponent. It catches errors in
        ! the normalisation constant B and in the pull-back that an anisotropic
        ! comparison alone would hide.
        !
        ! pi and pj map each of the 21 stored components back to the pair of Voigt
        ! indices it represents, which is what allows the Hessian to be written as a
        ! single loop over the storage order.
        real(real64), intent(in) :: s(6)
        real(real64), intent(out) :: q, g(6), h21(21)
        integer, parameter :: pi(21) = [1,2,3,4,5,6,1,2,3,4,5,1,2,3,4,1,2,3,1,2,1]
        integer, parameter :: pj(21) = [1,2,3,4,5,6,2,3,4,5,6,3,4,5,6,4,5,6,5,6,6]
        real(real64) :: d(6), isym
        integer :: n
        d = s
        d(1:3) = d(1:3) - sum(s(1:3))/3.0D0
        q = sqrt(1.5D0*(sum(d(1:3)**2) + 2.0D0*sum(d(4:6)**2)))
        g = 1.5D0*d/q
        do n = 1, 21
            ! Symmetric identity minus its spherical part; the 0.5 on the shear
            ! diagonal and the -1/3 on the normal block are the deviatoric projector.
            isym = 0.0D0
            if (pi(n) == pj(n)) isym = merge(1.0D0, 0.5D0, pi(n) <= 3)
            if (pi(n) <= 3 .and. pj(n) <= 3) isym = isym - 1.0D0/3.0D0
            h21(n) = 1.5D0/q*isym - 2.25D0/q**3*d(pi(n))*d(pj(n))
        end do
    end subroutine j2_reference

    function sym3(v) result(m)
        ! Voigt storage to full 3x3 matrix.
        real(real64), intent(in) :: v(6)
        real(real64) :: m(3,3)
        m = reshape([v(1), v(4), v(6), v(4), v(2), v(5), v(6), v(5), v(3)], [3,3])
    end function sym3

    function vec6(m) result(v)
        ! Full 3x3 matrix back to Voigt storage, symmetrising on the way.
        real(real64), intent(in) :: m(3,3)
        real(real64) :: v(6)
        v = [m(1,1), m(2,2), m(3,3), 0.5D0*(m(1,2) + m(2,1)), 0.5D0*(m(2,3) + m(3,2)), &
             0.5D0*(m(1,3) + m(3,1))]
    end function vec6

    function map_l(c1, c2, s) result(t)
        ! Sigma = C : dev(s)  (normal block C1, diagonal shear factors C2)
        real(real64), intent(in) :: c1(3,3), c2(3), s(6)
        real(real64) :: t(6), d(6)
        d = s
        d(1:3) = d(1:3) - sum(s(1:3))/3.0D0
        t(1:3) = matmul(c1, d(1:3))
        t(4:6) = c2*d(4:6)
    end function map_l

    function map_lt(c1, c2, w) result(t)
        ! Adjoint of map_l under the tensorial double contraction: dev(C^T : w)
        !
        ! No Voigt weights appear here: for the block structure CPB06 uses, a 3x3 normal
        ! block plus three diagonal shear factors, C commutes with the weights and its
        ! adjoint is represented by C^T. The Yld2004 reference, which accepts a general
        ! 6x6, has to keep them.
        real(real64), intent(in) :: c1(3,3), c2(3), w(6)
        real(real64) :: t(6)
        t(1:3) = matmul(transpose(c1), w(1:3))
        t(4:6) = c2*w(4:6)
        t(1:3) = t(1:3) - sum(t(1:3))/3.0D0
    end function map_lt

    subroutine tracepower_reference(c1, c2, a, s, q, g, act)
        ! CPB06 with k = 0 and even integer a: q = B*tr(Sigma**a)**(1/a),
        ! B = (sum |gamma_i|**a)**(-1/a), gamma_i = (2 C_i1 - C_i2 - C_i3)/3.
        !
        ! p(:,:,n) holds Sigma**n, built once and reused for the value, the gradient
        ! and every column of the Hessian. Its upper bound of 8 caps the exponent this
        ! reference supports.
        !
        ! The gradient follows from d tr(Sigma**a) = a Sigma**(a-1), scaled by kappa
        ! from the outer power and pulled back with the adjoint. The Hessian adds the
        ! derivative of kappa, a rank-one term, to the derivative of Sigma**(a-1),
        ! which is the non-commutative sum over r.
        real(real64), intent(in) :: c1(3,3), c2(3), s(6)
        integer, intent(in) :: a
        real(real64), intent(out) :: q, g(6), act(6,6)
        real(real64) :: p(3,3,0:8), gm(3,3), em(3,3), dg(3,3), gam(3), e(6)
        real(real64) :: b, phi, kappa, dphi, dkappa
        integer :: n, r, j

        gam = (2.0D0*c1(:,1) - c1(:,2) - c1(:,3))/3.0D0
        b = sum(abs(gam)**a)**(-1.0D0/a)
        p(:,:,0) = reshape([1.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0], [3,3])
        do n = 1, a
            p(:,:,n) = matmul(p(:,:,n-1), sym3(map_l(c1, c2, s)))
        end do
        phi = p(1,1,a) + p(2,2,a) + p(3,3,a)
        q = b*phi**(1.0D0/a)
        gm = a*p(:,:,a-1)
        kappa = q/(a*phi)
        g = map_lt(c1, c2, kappa*vec6(gm))
        do j = 1, 6
            e = 0.0D0
            e(j) = 1.0D0
            em = sym3(map_l(c1, c2, e))
            dphi = sum(gm*em)
            ! d(Sigma**(a-1)) : E = sum_r Sigma**r E Sigma**(a-2-r); every insertion
            ! position contributes because the factors do not commute.
            dg = 0.0D0
            do r = 0, a - 2
                dg = dg + matmul(matmul(p(:,:,r), em), p(:,:,a-2-r))
            end do
            dg = a*dg
            dkappa = q*(1.0D0 - a)/(a*a*phi*phi)*dphi
            act(:, j) = map_lt(c1, c2, vec6(dkappa*gm + kappa*dg))
        end do
    end subroutine tracepower_reference

end module test_cpb06_references
