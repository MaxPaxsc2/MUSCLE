! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module test_yld2004_references
    ! Independent references for the regression levels 5-7 (no spectral decomposition):
    !  - closed-form J2 (value, gradient, full Hessian),
    !  - Yld2004 with even exponent a through the 9x9 Kronecker difference
    !    K = S' (x) I - I (x) S'', whose eigenvalues are lambda'_i - lambda''_j:
    !    Phi = tr(K**a), a polynomial invariant that stays smooth when S' or S''
    !    has repeated eigenvalues. Any 6x6 linear map (tensorial Voigt storage).
    !
    ! Independence from the implementation
    ! ------------------------------------
    ! The criterion under test computes its value and derivatives through two
    ! eigen-decompositions, divided differences taken within each spectrum, and a
    ! three-block Daleckii-Krein assembly. These references use none of that: no call
    ! to eigen_sym3, no divided difference, no use of the yield module. A sign error in
    ! a cross coupling or a wrong limit in a divided difference therefore cannot cancel
    ! between the two sides.
    !
    ! The independence is of the algorithm, not of everything. Both sides use the tensor
    ! representation and contraction operators of muscle_tensors, the same transformation
    ! conventions, and the same fixed normalisation factor.
    !
    ! The Kronecker identity
    ! ----------------------
    ! The difficulty specific to Yld2004 is that its nine terms mix two spectra, so no
    ! single matrix has the cross differences as eigenvalues. The Kronecker difference
    ! supplies one: for K = S' (x) I - I (x) S'' acting on the nine-dimensional tensor
    ! product space, the eigenvalues are exactly the nine lambda'_i - lambda''_j, with
    ! eigenvectors the products of the individual eigenvectors.
    !
    ! For an even integer a the absolute value in Phi is redundant, so
    !   Phi = sum_ij (lambda'_i - lambda''_j)**a = tr(K**a),
    ! a polynomial invariant of a 9x9 matrix, computable by repeated multiplication
    ! with no eigenvalue ever extracted.
    !
    ! The restriction to even a limits coverage, but the polynomial route is smooth
    ! where the spectral one is most fragile: repeated eigenvalues in either spectrum,
    ! which force the divided differences into their limits, are an ordinary point for
    ! tr(K**a).
    !
    ! Matrix derivatives
    ! ------------------
    ! Three identities drive the derivative references:
    !   d/dX tr(X**a) = a X**(a-1)
    !   d/dX (X**n) : E = sum_{r=0}^{n-1} X**r E X**(n-1-r)
    !   dK = dS' (x) I - I (x) dS''
    ! The sum in the second is not an artifact: matrices do not commute, so E has to be
    ! inserted in every position along the product. The third turns a derivative with
    ! respect to K into derivatives with respect to each transformed tensor, which is
    ! why the partial traces appear: tracing over one factor of the product space
    ! extracts the part that belongs to the other.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    private
    public :: action_matrix, rel_err, j2_reference, kron_reference

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
        ! The isotropic baseline: von Mises in closed form. Yld2004 reproduces it when
        ! both transformations are the identity and a = 2, which is the configuration
        ! the calling tests use; the agreement is not a property of every exponent. It
        ! checks the factor 1/4 in the normalisation and the sign of the second
        ! contribution to the gradient.
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

    function map_l(c, s) result(t)
        ! S = C dev(s), tensorial Voigt components
        real(real64), intent(in) :: c(6,6), s(6)
        real(real64) :: t(6), d(6)
        d = s
        d(1:3) = d(1:3) - sum(s(1:3))/3.0D0
        t = matmul(c, d)
    end function map_l

    function map_lt(c, w) result(t)
        ! Adjoint of map_l under the tensorial double contraction: dev(W^-1 C^T W w)
        !
        ! Unlike the CPB06 reference, this one accepts a general 6x6 map, so the Voigt
        ! weights cannot be dropped: the shear components stand for two equal entries
        ! each, and the adjoint has to respect that inner product.
        real(real64), intent(in) :: c(6,6), w(6)
        real(real64) :: t(6)
        real(real64), parameter :: wv(6) = [1.0D0, 1.0D0, 1.0D0, 2.0D0, 2.0D0, 2.0D0]
        t = matmul(transpose(c), wv*w)/wv
        t(1:3) = t(1:3) - sum(t(1:3))/3.0D0
    end function map_lt

    function kron_diff(a, b) result(k)
        ! a (x) I - I (x) b, row/column index 3*(i-1)+j
        !
        ! The 9x9 matrix whose eigenvalues are the nine differences a_i - b_j. Each
        ! index of the product space is the pair (i,j) flattened, and the two terms act
        ! on one member of the pair each, which is why the Kronecker deltas appear.
        real(real64), intent(in) :: a(3,3), b(3,3)
        real(real64) :: k(9,9)
        integer :: i, j, m, n
        do i = 1, 3
            do j = 1, 3
                do m = 1, 3
                    do n = 1, 3
                        k(3*(i-1)+j, 3*(m-1)+n) = merge(a(i,m), 0.0D0, j == n) &
                                                 - merge(b(j,n), 0.0D0, i == m)
                    end do
                end do
            end do
        end do
    end function kron_diff

    function ptrace2(k) result(p)
        ! Partial trace over the second factor
        !
        ! Collapses the 9x9 back to a 3x3 acting on the first spectrum alone, which is
        ! how a derivative with respect to K becomes a derivative with respect to S'.
        real(real64), intent(in) :: k(9,9)
        real(real64) :: p(3,3)
        integer :: i, m, j
        p = 0.0D0
        do i = 1, 3
            do m = 1, 3
                do j = 1, 3
                    p(i,m) = p(i,m) + k(3*(i-1)+j, 3*(m-1)+j)
                end do
            end do
        end do
    end function ptrace2

    function ptrace1(k) result(p)
        ! Partial trace over the first factor
        !
        ! The counterpart of ptrace2 for the second spectrum. Its contribution enters
        ! with a minus sign, inherited from the minus in the Kronecker difference.
        real(real64), intent(in) :: k(9,9)
        real(real64) :: p(3,3)
        integer :: i, j, n
        p = 0.0D0
        do j = 1, 3
            do n = 1, 3
                do i = 1, 3
                    p(j,n) = p(j,n) + k(3*(i-1)+j, 3*(i-1)+n)
                end do
            end do
        end do
    end function ptrace1

    subroutine kron_reference(cp, cd, a, s, q, g, act)
        ! q = (tr(K**a)/4)**(1/a); dPhi = a tr(K**(a-1) dK),
        ! dK = dS' (x) I - I (x) dS''  =>  dPhi/dS' = a Tr_2(K**(a-1)),
        !                                  dPhi/dS'' = -a Tr_1(K**(a-1))
        !
        ! p(:,:,n) holds K**n, built once and reused for the value, the gradient and
        ! every column of the Hessian. Its upper bound of 8 caps the exponent this
        ! reference supports.
        !
        ! The Hessian adds the derivative of kappa, a rank-one term, to the derivative
        ! of K**(a-1), which is the non-commutative sum over r, and then splits the
        ! result between the two transformations with the partial traces.
        real(real64), intent(in) :: cp(6,6), cd(6,6), s(6)
        integer, intent(in) :: a
        real(real64), intent(out) :: q, g(6), act(6,6)
        real(real64) :: p(9,9,0:8), m(9,9), dk(9,9), dm(9,9), gp(3,3), gd(3,3), e(6)
        real(real64) :: phi, kappa, dphi, dkappa
        integer :: n, r, j

        p(:,:,0) = 0.0D0
        do n = 1, 9
            p(n,n,0) = 1.0D0
        end do
        p(:,:,1) = kron_diff(sym3(map_l(cp, s)), sym3(map_l(cd, s)))
        do n = 2, a
            p(:,:,n) = matmul(p(:,:,n-1), p(:,:,1))
        end do
        phi = 0.0D0
        do n = 1, 9
            phi = phi + p(n,n,a)
        end do
        q = (0.25D0*phi)**(1.0D0/a)
        m = p(:,:,a-1)
        gp = a*ptrace2(m)
        gd = -a*ptrace1(m)
        kappa = q/(a*phi)
        g = map_lt(cp, vec6(kappa*gp)) + map_lt(cd, vec6(kappa*gd))
        do j = 1, 6
            e = 0.0D0
            e(j) = 1.0D0
            dk = kron_diff(sym3(map_l(cp, e)), sym3(map_l(cd, e)))
            dphi = a*sum(m*dk)
            ! d(K**(a-1)) : E = sum_r K**r E K**(a-2-r); every insertion position
            ! contributes because the factors do not commute.
            dm = 0.0D0
            do r = 0, a - 2
                dm = dm + matmul(matmul(p(:,:,r), dk), p(:,:,a-2-r))
            end do
            dkappa = q*(1.0D0 - a)/(a*a*phi*phi)*dphi
            act(:, j) = map_lt(cp, vec6(dkappa*gp + kappa*a*ptrace2(dm))) &
                      + map_lt(cd, vec6(dkappa*gd - kappa*a*ptrace1(dm)))
        end do
    end subroutine kron_reference

end module test_yld2004_references
