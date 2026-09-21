! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_yield_cpb06
    !! Module muscle_yield_cpb06
    !! =========================
    !! Implements the Cazacu-Plunkett-Barlat (2006) orthotropic yield criterion, which
    !! captures the tension-compression asymmetry of hexagonal close-packed metals,
    !! together with its analytical first and second derivatives with respect to stress.
    !!
    !! Mathematical Formulation
    !! ------------------------
    !! The stress deviator is mapped through a fourth-order orthotropic tensor,
    !! \( \boldsymbol{\Sigma} = \mathbf{C} : \mathbf{s} \), and the criterion is built
    !! from the three principal values \( \Sigma_i \) of that transformed tensor:
    !! \[ \Phi = \sum_{i=1}^{3} \left( |\Sigma_i| - k\,\Sigma_i \right)^a \]
    !! \[ \bar{\sigma} = B\,\Phi^{1/a} \]
    !! The parameter \(k\) shifts the weight between tensile and compressive principal
    !! values and is what produces the asymmetry: with \(k = 0\) the criterion is
    !! symmetric, and the sign of \(k\) selects which of the two is favoured.
    !!
    !! \(B\) normalises the result so that \( \bar{\sigma} \) equals the tensile yield
    !! stress along the reference axis \(x\). It is obtained by evaluating \(\Phi\) for a
    !! unit uniaxial stress in that direction, whose transformed tensor is diagonal with
    !! entries \( \gamma_i = (2C_{i1} - C_{i2} - C_{i3})/3 \):
    !! \[ B = \left[ \sum_{i=1}^{3} \left( |\gamma_i| - k\,\gamma_i \right)^a \right]^{-1/a} \]
    !!
    !! Since \( \mathbf{C} \) acts on the deviator, the criterion is insensitive to
    !! hydrostatic pressure. It is homogeneous of degree one in the stress.
    !!
    !! Parameters
    !! ----------
    !! - \( C_{ij} \): nine coefficients acting on the normal components plus three
    !!   acting on the shear components. The identity recovers an isotropic criterion.
    !! - \( k \in [-1, 1] \): tension-compression asymmetry. Outside that range the
    !!   quantity \( |\Sigma_i| - k\Sigma_i \) becomes negative for one sign of
    !!   \( \Sigma_i \) and the criterion loses meaning.
    !! - \( a \ge 1 \): exponent controlling the curvature of the surface. Values
    !!   \( a \le 2 \) make the yield function non-smooth at \( \Sigma_i = 0 \)
    !!   (see *Numerical treatment*).
    !!
    !! Numerical treatment
    !! -------------------
    !! This implementation evaluates \( \bar{\sigma} \) and its derivatives spectrally.
    !! A polynomial route exists for special parameter sets — with \( k = 0 \) and even
    !! integer \(a\), \( \Phi = \mathrm{tr}(\boldsymbol{\Sigma}^a) \), which the
    !! regression references exploit — but none covers the general case. Three
    !! difficulties are handled explicitly:
    !!
    !! 1. **Eigen-decomposition.** `eigen_sym3` (cyclic Jacobi) is used instead of an
    !!    analytical root formula, which would divide by quantities that vanish on
    !!    repeated roots. The same solver serves the value and both derivatives, so the
    !!    three are always consistent with each other.
    !! 2. **Repeated roots.** The second derivative of a spectral function involves the
    !!    divided differences \( (g_i - g_j)/(\Sigma_i - \Sigma_j) \) of the
    !!    Daleckii-Krein formula, which are indeterminate when two principal values
    !!    coincide. They are replaced by their limit, the derivative itself, and by a
    !!    Taylor expansion when the values are merely close, where the direct quotient
    !!    would lose most of its significant digits to cancellation.
    !! 3. **The kink at \( \Sigma_i = 0 \).** For \( a \le 2 \) and \( k \ne 0 \) the
    !!    scalar term \( (|\Sigma| - k\Sigma)^a \) is not twice differentiable at the
    !!    origin. At \( a = 2 \) the implementation returns the average of the two
    !!    one-sided slopes, which is also the limit of a central difference. For
    !!    \( 1 < a < 2 \) those slopes diverge and the returned value is a convention
    !!    that keeps the Hessian finite rather than the limit of anything; the same
    !!    holds at the contract limit \( a = 1 \).
    !!
    !! Limitations
    !! -----------
    !! The normalisation needs a strictly positive reference value. If every
    !! \( \gamma_i \) vanishes, which happens when all entries of \( C_1 \) are equal,
    !! \(B\) is not finite; the coefficients must be chosen to avoid that.
    !!
    !! Both derivatives return zero when \( \Phi \le 0 \). This occurs only in
    !! degenerate configurations: a vanishing transformed deviator, or \( k = \pm 1 \)
    !! with every principal value of the favoured sign. That zero is a safeguard, not
    !! a computed derivative, and callers should not read it as one.
    !!
    !! Reference
    !! ---------
    !! Cazacu, O., Plunkett, B., & Barlat, F. (2006). *Orthotropic yield criterion for
    !! hexagonal closed packed metals.* International Journal of Plasticity, 22(7),
    !! 1171-1194.
    !!
    !! For more information see [[muscle_yield_base]] and [[muscle_math_spectral_derivs]]

    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_base
    implicit none
    private

    real(real64), parameter :: CPB_ZERO_REL = 1.0D-12
        !! |lambda| <= CPB_ZERO_REL*max|lambda| is treated as a zero transformed eigenvalue
        !! (second-derivative convention for a <= 2, k /= 0; see cpb_dh). Relative, not absolute.
        !!
        !! The threshold is relative because the principal values scale with the stress
        !! level. A fixed absolute tolerance would change which eigenvalues count as
        !! zero when the same physical state is expressed in different units, since the
        !! numbers shrink by six orders of magnitude on passing from pascals to
        !! megapascals.
    real(real64), parameter :: CPB_DD_REL = 1.0D-3
        !! Relative gap below which divided differences use their Taylor expansion
        !!
        !! Above this gap the direct quotient is accurate; below it, cancellation in the
        !! numerator generally costs more digits than the expansion's truncation error.
        !! The threshold is a switching point, not a uniform accuracy bound: see cpb_dpow.

    public :: CPB06
    type, extends(Base_yield_critera) :: CPB06
        !! Cazacu-Plunkett-Barlat (2006) orthotropic yield criterion with
        !! tension-compression asymmetry.
        !!
        !! Coefficients are stored as supplied; no admissibility check is performed on
        !! them, and in particular \(k\) is not verified to lie in \([-1, 1]\).
        real(real64), dimension(3,3) :: C1
            !! Orthotropic coefficients acting on the normal components of the deviator.
        real(real64), dimension(3) :: C2
            !! Orthotropic coefficients acting on the \(xy\), \(yz\) and \(xz\) shear
            !! components, in that order.
        real(real64) :: k
            !! Tension-compression asymmetry parameter, expected in \([-1, 1]\).
        real(real64) :: a
            !! Exponent of the yield function, expected \( \ge 1 \).
    contains
        procedure :: init
        procedure :: stress_eq
        procedure :: dstressEq_dstress => dstressEq_dstress_cpb
        procedure :: ddstressEq_ddstress => ddstressEq_ddstress_cpb
    end type CPB06

contains

    subroutine init(self,           &
                    c11, c12, c13,  &
                    c21, c22, c23,  &
                    c31, c32, c33,  &
                    c44, c55, c66,  &
                    k, a)
        !! Stores the twelve orthotropic coefficients and the two scalar parameters.
        !!
        !! The nine normal coefficients are laid out row by row into the \(3\times3\)
        !! matrix \( C_1 \); note that `reshape` fills by columns, so the argument list
        !! is transposed on the way in. Setting \( C_1 = \mathbf{I} \), \( C_2 = 1 \)
        !! and \( k = 0 \) recovers an isotropic criterion.
        implicit none
        class(CPB06), intent(inout) :: self
        real(real64), intent(in) :: c11, c12, c13, c21, c22, c23, c31, c32, c33, c44, c55, c66
            !! Orthotropic coefficients: nine for the normal components, then the
            !! \(xy\), \(yz\) and \(xz\) shear coefficients.
        real(real64), intent(in) :: k, a
            !! Asymmetry parameter and exponent.

        self%C1 = reshape(  &
                          [c11, c21, c31, &
                           c12, c22, c32, &
                           c13, c23, c33], &
                          [3,3] )

        self%C2 = [c44, c55, c66]
        self%k = k
        self%a = a

    end subroutine init

    pure function stress_eq(self, stress) result(res)
        !! Returns the CPB06 equivalent stress of a given stress state.
        !!
        !! The deviator is transformed by \( \mathbf{C} \), its principal values are
        !! extracted, and the yield function is assembled and normalised by \(B\).
        use muscle_tensors
        use muscle_math_spectral_derivs, only: eigen_sym3
        implicit none
        class(CPB06), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
            !! Cauchy stress, expressed in the orthotropic axes of the material.
        real(real64) :: res
            !! Equivalent stress, normalised to the tensile yield stress along \(x\).

        real(real64) :: k, a
        type(ten_3D2Osym) :: sigma
        type(ten_3D2Osym) :: dev
        real(real64), dimension(3) :: sigma_eig
        real(real64), dimension(3,3) :: sigma_vec

        real(real64) :: sxx, syy, szz, sxy, syz, sxz
        real(real64), dimension(3) :: gamma
        real(real64) :: b

        k = self%k
        a = self%a

        dev = .dev. stress
        sxx = dev%xx()*self%C1(1,1) + dev%yy()*self%C1(1,2) + dev%zz()*self%C1(1,3)
        syy = dev%xx()*self%C1(2,1) + dev%yy()*self%C1(2,2) + dev%zz()*self%C1(2,3)
        szz = dev%xx()*self%C1(3,1) + dev%yy()*self%C1(3,2) + dev%zz()*self%C1(3,3)
        sxy = dev%xy()*self%C2(1)
        syz = dev%yz()*self%C2(2)
        sxz = dev%xz()*self%C2(3)

        call sigma%init(xx=sxx, yy=syy, zz=szz, xy=sxy, yz=syz, xz=sxz)
        ! Same eigen-solver as the derivatives; no division by a vanishing quantity.
        call eigen_sym3(sigma, sigma_eig, sigma_vec)

        res = sum((abs(sigma_eig) - k*sigma_eig)**a)**(1D0/a)

        ! gamma = principal values of the transformed deviator of a unit uniaxial
        ! stress along x, whose deviator is diag(2/3, -1/3, -1/3). The transformed
        ! tensor is diagonal for that state, so its entries are already its eigenvalues.
        gamma(1) = (2D0*self%C1(1,1) - self%C1(1,2) - self%C1(1,3))/3D0
        gamma(2) = (2D0*self%C1(2,1) - self%C1(2,2) - self%C1(2,3))/3D0
        gamma(3) = (2D0*self%C1(3,1) - self%C1(3,2) - self%C1(3,3))/3D0

        b = sum((abs(gamma) - k*gamma)**a)**(-1D0/a)

        res = res*b

    end function stress_eq

    ! ----------------------------------------------------------------------------
    ! Linear maps and their transposes.
    !
    ! L = C : dev is the composition applied to the stress; its transpose L^T is
    ! needed to pull derivatives taken with respect to Sigma back to derivatives
    ! with respect to the stress. The deviatoric projector is self-adjoint, and for
    ! the block structure used here (a 3x3 normal block plus three diagonal shear
    ! factors) the adjoint of C is represented by C^T. So L* amounts to applying C^T
    ! and then taking the deviator, in the reverse order of L.
    ! ----------------------------------------------------------------------------

    pure function apply_C(self, s) result(res)
        ! Applies the orthotropic transformation C to a symmetric tensor.
        implicit none
        class(CPB06), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: s
        type(ten_3D2Osym) :: res

        res%vals(1) = s%xx()*self%C1(1,1) + s%yy()*self%C1(1,2) + s%zz()*self%C1(1,3)
        res%vals(2) = s%xx()*self%C1(2,1) + s%yy()*self%C1(2,2) + s%zz()*self%C1(2,3)
        res%vals(3) = s%xx()*self%C1(3,1) + s%yy()*self%C1(3,2) + s%zz()*self%C1(3,3)
        res%vals(4) = s%xy()*self%C2(1)
        res%vals(5) = s%yz()*self%C2(2)
        res%vals(6) = s%xz()*self%C2(3)
    end function apply_C

    pure function apply_CT(self, w) result(res)
        ! Applies the transpose of C. Only the normal block differs from apply_C:
        ! the shear coefficients act diagonally and are therefore unchanged.
        implicit none
        class(CPB06), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: w
        type(ten_3D2Osym) :: res

        res%vals(1) = w%xx()*self%C1(1,1) + w%yy()*self%C1(2,1) + w%zz()*self%C1(3,1)
        res%vals(2) = w%xx()*self%C1(1,2) + w%yy()*self%C1(2,2) + w%zz()*self%C1(3,2)
        res%vals(3) = w%xx()*self%C1(1,3) + w%yy()*self%C1(2,3) + w%zz()*self%C1(3,3)
        res%vals(4) = w%xy()*self%C2(1)
        res%vals(5) = w%yz()*self%C2(2)
        res%vals(6) = w%xz()*self%C2(3)
    end function apply_CT

    pure function apply_L(self, stress) result(res)
        ! Sigma = C : dev(stress).
        implicit none
        class(CPB06), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
        type(ten_3D2Osym) :: res

        res = apply_C(self, .dev. stress)
    end function apply_L

    pure function apply_LT(self, w) result(res)
        ! Transpose of apply_L: the two operations are applied in reverse order.
        implicit none
        class(CPB06), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: w
        type(ten_3D2Osym) :: res

        res = .dev. apply_CT(self, w)
    end function apply_LT

    pure function cpb_B(self) result(b)
        ! Normalisation constant B. See the module header for the derivation of gamma.
        implicit none
        class(CPB06), intent(in) :: self
        real(real64) :: b
        real(real64), dimension(3) :: gamma
        real(real64) :: k, a

        k = self%k
        a = self%a
        gamma(1) = (2.0D0*self%C1(1,1) - self%C1(1,2) - self%C1(1,3))/3.0D0
        gamma(2) = (2.0D0*self%C1(2,1) - self%C1(2,2) - self%C1(2,3))/3.0D0
        gamma(3) = (2.0D0*self%C1(3,1) - self%C1(3,2) - self%C1(3,3))/3.0D0
        b = sum((abs(gamma) - k*gamma)**a)**(-1.0D0/a)
    end function cpb_B

    ! ----------------------------------------------------------------------------
    ! Scalar building blocks of the yield function.
    !
    ! Writing psi(lam) = |lam| - k*lam, the yield function is sum(psi**a) and its
    ! derivatives with respect to the eigenvalues reduce to derivatives of psi**a.
    ! These helpers isolate that scalar algebra, including the two places where a
    ! naive evaluation breaks down: the kink of |lam| at the origin, and the
    ! cancellation in divided differences of nearly equal arguments.
    ! ----------------------------------------------------------------------------

    pure function cpb_sgnkm(lam, k) result(s)
        ! d(psi)/dlam = sgn(lam) - k, with the convention sgn(0) = 0. The value at
        ! lam = 0 is the average of the two one-sided limits.
        implicit none
        real(real64), intent(in) :: lam, k
        real(real64) :: s

        if (lam > 0.0D0) then
            s = 1.0D0 - k
        else if (lam < 0.0D0) then
            s = -1.0D0 - k
        else
            s = -k
        end if
    end function cpb_sgnkm

    pure function cpb_h(lam, k, a) result(h)
        ! h(lam) = (|lam| - k*lam)**(a-1) * (sgn(lam) - k);  d/dlam (|lam| - k*lam)**a = a*h
        ! The guard psi > 0 covers k = +-1, where psi vanishes identically on one
        ! branch and the power would otherwise raise zero to a negative exponent.
        implicit none
        real(real64), intent(in) :: lam, k, a
        real(real64) :: h, psi

        psi = abs(lam) - k*lam
        if (psi > 0.0D0) then
            h = psi**(a - 1.0D0)*cpb_sgnkm(lam, k)
        else
            h = 0.0D0
        end if
    end function cpb_h

    pure function cpb_dh(lam, k, a, zero_tol) result(dh)
        ! h'(lam) = (a-1)*(|lam| - k*lam)**(a-2)*(sgn(lam) - k)**2 for lam /= 0.
        ! lam = 0 and a <= 2: h is not differentiable there; keep the previous
        ! convention (the regimes below say what that value is and is not).
        !
        ! The regimes differ. For a > 2 the limit at the origin is zero, so the branch
        ! is exact there; within the tolerance, small nonzero arguments are approximated
        ! by that limit. For a = 2 the function h is piecewise linear with slopes
        ! (1-k)**2 and (1+k)**2, and the value returned is both the average of the
        ! one-sided derivatives and the limit of a central difference. For 1 < a < 2 the
        ! derivative diverges on each branch that carries nonzero weight; at k = +-1 the
        ! zero-weight branch is flat instead, but a central difference still grows
        ! without bound. The value returned there is a convention that keeps the Hessian
        ! finite, not a limit, and the same holds at the contract limit a = 1, where h
        ! is discontinuous at the origin.
        implicit none
        real(real64), intent(in) :: lam, k, a, zero_tol
        real(real64) :: dh, psi, s

        if (abs(lam) <= zero_tol) then
            if (a > 2.0D0) then
                dh = 0.0D0
            else
                dh = (a - 1.0D0)*0.5D0*((1.0D0 - k)**2 + (1.0D0 + k)**2)
            end if
            return
        end if
        psi = abs(lam) - k*lam
        if (psi > 0.0D0) then
            s = cpb_sgnkm(lam, k)
            dh = (a - 1.0D0)*psi**(a - 2.0D0)*s*s
        else
            dh = 0.0D0
        end if
    end function cpb_dh

    pure function cpb_dpow(u, v, p) result(res)
        ! (u**p - v**p)/(u - v) for u, v > 0; Taylor expansion about the midpoint
        ! when the relative gap is small (the direct quotient loses digits there).
        !
        ! The expansion is written in the squared relative half-gap r, so only even
        ! orders survive; it retains terms through r**2. Its truncation error depends on
        ! both p and the gap, and not monotonically in p: it vanishes outright at small
        ! positive integer p, where the series terminates. CPB_DD_REL is a switching
        ! point, not a uniform accuracy bound.
        implicit none
        real(real64), intent(in) :: u, v, p
        real(real64) :: res, m, r

        if (abs(u - v) > CPB_DD_REL*max(u, v)) then
            res = (u**p - v**p)/(u - v)
        else
            m = 0.5D0*(u + v)
            r = (0.5D0*(u - v)/m)**2
            res = p*m**(p - 1.0D0)*(1.0D0 + (p - 1.0D0)*(p - 2.0D0)*r/6.0D0 &
                  *(1.0D0 + (p - 3.0D0)*(p - 4.0D0)*r/20.0D0))
        end if
    end function cpb_dpow

    pure function cpb_divdiff(x, y, k, a, zero_tol) result(res)
        ! (h(x) - h(y))/(x - y), with its limit h'(x) at x = y.
        !
        ! Three cases, chosen so that no branch ever divides by a small difference:
        ! both arguments at the origin, both on the same branch of psi (where the
        ! problem reduces to a difference of powers and cpb_dpow handles the near-equal
        ! case), and opposite signs, where |x - y| is at least max(|x|,|y|) and the
        ! direct quotient is well conditioned.
        implicit none
        real(real64), intent(in) :: x, y, k, a, zero_tol
        real(real64) :: res, s, u, v

        if (abs(x) <= zero_tol .and. abs(y) <= zero_tol) then
            res = cpb_dh(0.0D0, k, a, zero_tol)
        else if ((x > 0.0D0 .and. y > 0.0D0) .or. (x < 0.0D0 .and. y < 0.0D0)) then
            ! Same branch of |lam| - k*lam = s*lam: h = s*(s*lam)**(a-1)
            s = cpb_sgnkm(x, k)
            u = s*x
            v = s*y
            if (u > 0.0D0 .and. v > 0.0D0) then
                res = s*s*cpb_dpow(u, v, a - 1.0D0)
            else
                res = 0.0D0
            end if
        else
            ! Opposite signs (or one zero): |x - y| >= max(|x|,|y|), well conditioned
            res = (cpb_h(x, k, a) - cpb_h(y, k, a))/(x - y)
        end if
    end function cpb_divdiff

    pure function pack_hess_from_cols(col) result(res)
        ! Assembles the 21 stored components of the Hessian from its six columns.
        !
        ! Each column is the Hessian contracted with one Voigt basis tensor. The shear
        ! columns are halved because the basis tensor e_J carries the off-diagonal pair
        ! twice, while the stored component represents it once.
        implicit none
        type(ten_3D2Osym), intent(in) :: col(6)
        type(ten_3D4O3sym) :: res

        ! col(J) = H : e_J. Shear columns include the tensorial factor 2.
        res%vals(1)  = col(1)%vals(1)
        res%vals(2)  = col(2)%vals(2)
        res%vals(3)  = col(3)%vals(3)
        res%vals(4)  = 0.5D0*col(4)%vals(4)
        res%vals(5)  = 0.5D0*col(5)%vals(5)
        res%vals(6)  = 0.5D0*col(6)%vals(6)
        res%vals(7)  = col(1)%vals(2)
        res%vals(8)  = col(2)%vals(3)
        res%vals(9)  = 0.5D0*col(4)%vals(3)
        res%vals(10) = 0.5D0*col(5)%vals(4)
        res%vals(11) = 0.5D0*col(6)%vals(5)
        res%vals(12) = col(1)%vals(3)
        res%vals(13) = 0.5D0*col(4)%vals(2)
        res%vals(14) = 0.5D0*col(5)%vals(3)
        res%vals(15) = 0.5D0*col(6)%vals(4)
        res%vals(16) = 0.5D0*col(4)%vals(1)
        res%vals(17) = 0.5D0*col(5)%vals(2)
        res%vals(18) = 0.5D0*col(6)%vals(3)
        res%vals(19) = 0.5D0*col(5)%vals(1)
        res%vals(20) = 0.5D0*col(6)%vals(2)
        res%vals(21) = 0.5D0*col(6)%vals(1)
    end function pack_hess_from_cols

    pure function dstressEq_dstress_cpb(self, stress) result(res)
        !! Returns the gradient of the equivalent stress, which under an associated
        !! flow rule gives the direction of plastic flow.
        !!
        !! The derivative of a spectral function shares the eigenbasis of its argument,
        !! so the gradient with respect to \( \boldsymbol{\Sigma} \) is diagonal in that
        !! basis and only its three eigenvalues have to be computed. It is then pulled
        !! back to the stress with the transpose of the transformation.
        !!
        !! Repeated principal values need no special treatment here: the diagonal
        !! entries coincide on a repeated root, so the result does not depend on which
        !! eigenbasis the solver happens to return for the cluster.
        use muscle_tensors
        use muscle_math_spectral_derivs, only: eigen_sym3, spectral_compose
        implicit none
        class(CPB06), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
            !! Cauchy stress, expressed in the orthotropic axes of the material.
        type(ten_3D2Osym) :: res
            !! Gradient of the equivalent stress. Zero in the degenerate case
            !! \( \Phi \le 0 \); that value is a safeguard, not a computed derivative.

        type(ten_3D2Osym) :: Sigma, dPhi
        real(real64) :: eigs(3), vecs(3,3), gdiag(3,3)
        real(real64) :: k, a, b, Phi, seq
        integer :: i

        k = self%k
        a = self%a
        Sigma = apply_L(self, stress)
        call eigen_sym3(Sigma, eigs, vecs)
        Phi = sum((abs(eigs) - k*eigs)**a)
        if (Phi <= 0.0D0) then
            res = 0.0D0
            return
        end if

        b = cpb_B(self)
        seq = b*Phi**(1.0D0/a)

        ! dseq/dSigma = sum_i g_i v_i (x) v_i,  g_i = (seq/Phi)*h(lambda_i).
        ! g_i = g_j on a repeated root, so any eigenbasis of the cluster works.
        gdiag = 0.0D0
        do i = 1, 3
            gdiag(i,i) = (seq/Phi)*cpb_h(eigs(i), k, a)
        end do
        dPhi = spectral_compose(vecs, gdiag)

        res = apply_LT(self, dPhi)
    end function dstressEq_dstress_cpb

    pure function ddstressEq_ddstress_cpb(self, stress) result(res)
        !! Returns the second derivative of the equivalent stress, required by the
        !! consistent tangent operator of implicit return-mapping schemes.
        !!
        !! Unlike the gradient, the second derivative of a spectral function is not
        !! diagonal in the eigenbasis. The Daleckii-Krein formula gives it as two
        !! distinct blocks: the diagonal entries follow from the derivatives of the
        !! \( g_i \) with respect to the principal values, while the off-diagonal
        !! entries are the divided differences \( (g_i - g_j)/(\Sigma_i - \Sigma_j) \),
        !! which measure how the eigenvectors themselves rotate.
        !!
        !! The Hessian is assembled one column at a time: each Voigt basis tensor is
        !! pushed forward through the transformation, rotated into the eigenbasis,
        !! multiplied blockwise, and pulled back.
        use muscle_tensors
        use muscle_math_spectral_derivs, only: eigen_sym3, spectral_rotate, spectral_compose
        implicit none
        class(CPB06), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
            !! Cauchy stress, expressed in the orthotropic axes of the material.
        type(ten_3D4O3sym) :: res
            !! Second derivative of the equivalent stress. Zero in the degenerate case
            !! \( \Phi \le 0 \).

        type(ten_3D2Osym) :: Sigma, e
        type(ten_3D2Osym) :: col(6)
        real(real64) :: eigs(3), vecs(3,3), g(3), dgdl(3,3), theta(3,3), z(3,3), dm(3,3)
        real(real64) :: k, a, b, Phi, seq, c, zero_tol
        integer :: i, j

        k = self%k
        a = self%a
        Sigma = apply_L(self, stress)
        call eigen_sym3(Sigma, eigs, vecs)
        Phi = sum((abs(eigs) - k*eigs)**a)
        if (Phi <= 0.0D0) then
            res = 0.0D0
            return
        end if

        b = cpb_B(self)
        seq = b*Phi**(1.0D0/a)
        c = seq/Phi
        zero_tol = CPB_ZERO_REL*maxval(abs(eigs))

        do i = 1, 3
            g(i) = c*cpb_h(eigs(i), k, a)
        end do

        ! dg_i/dlambda_j: rank-one term from seq plus the diagonal h' term
        do j = 1, 3
            do i = 1, 3
                dgdl(i,j) = (1.0D0 - a)*g(i)*g(j)/seq
            end do
            dgdl(j,j) = dgdl(j,j) + c*cpb_dh(eigs(j), k, a, zero_tol)
        end do

        ! Daleckii-Krein: theta_ij = (g_i - g_j)/(lambda_i - lambda_j) = c*h[lambda_i, lambda_j],
        ! evaluated without cancellation; its limit on a repeated root is c*h'(lambda).
        theta = 0.0D0
        do j = 2, 3
            do i = 1, j - 1
                theta(i,j) = c*cpb_divdiff(eigs(i), eigs(j), k, a, zero_tol)
                theta(j,i) = theta(i,j)
            end do
        end do

        ! d2(seq)/dSigma2 : E = V M V^T, z = V^T E V,
        !   M_ii = sum_j dgdl_ij z_jj,  M_ij = theta_ij z_ij (i /= j)
        ! Pull back column by column: P_dev : C^T : Hsig : C : P_dev
        do j = 1, 6
            e = 0.0D0
            e%vals(j) = 1.0D0
            z = spectral_rotate(vecs, apply_L(self, e))
            dm = theta*z
            do i = 1, 3
                dm(i,i) = dgdl(i,1)*z(1,1) + dgdl(i,2)*z(2,2) + dgdl(i,3)*z(3,3)
            end do
            col(j) = apply_LT(self, spectral_compose(vecs, dm))
        end do
        res = pack_hess_from_cols(col)
    end function ddstressEq_ddstress_cpb

end module muscle_yield_cpb06
