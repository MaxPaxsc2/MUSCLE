! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_yield_yld2004
    !! Module muscle_yield_yld2004
    !! ===========================
    !! Implements the Barlat Yld2004-18p anisotropic yield criterion, built on two
    !! independent linear transformations of the stress deviator, together with its
    !! analytical first and second derivatives with respect to stress.
    !!
    !! Mathematical Formulation
    !! ------------------------
    !! The deviator is mapped through two fourth-order orthotropic tensors,
    !! \( \mathbf{S}' = \mathbf{C}' : \mathbf{s} \) and
    !! \( \mathbf{S}'' = \mathbf{C}'' : \mathbf{s} \), and the criterion is assembled
    !! from every difference between a principal value of the first and a principal
    !! value of the second:
    !! \[ \Phi = \sum_{i=1}^{3}\sum_{j=1}^{3}
    !!     \left| \lambda'_i - \lambda''_j \right|^a \]
    !! \[ \bar{\sigma} = \left( \tfrac{1}{4}\Phi \right)^{1/a} \]
    !! The nine cross differences are what give the criterion its flexibility: with a
    !! single transformation only the three differences within one spectrum would be
    !! available, which is the Yld91 family.
    !!
    !! The factor \(1/4\) is a fixed constant. It gives unit equivalent stress under
    !! unit uniaxial tension when both transformations are the identity: the principal
    !! values are then those of the deviator itself, a unit uniaxial stress gives
    !! \( \mathbf{s} = \mathrm{diag}(2/3, -1/3, -1/3) \), and its nine cross differences
    !! are four values of magnitude one and five zeros, so \( \Phi = 4 \). For
    !! anisotropic coefficients the factor alone normalises nothing: the transformations
    !! themselves must be calibrated to the chosen reference yield stress. Scaling both
    !! by two, for instance, doubles \( \bar{\sigma} \).
    !!
    !! Since both transformations act on the deviator, the criterion is insensitive to
    !! hydrostatic pressure. It is homogeneous of degree one in the stress.
    !!
    !! Parameters
    !! ----------
    !! - \( \mathbf{C}' \), \( \mathbf{C}'' \): the two transformations, supplied as
    !!   \(6\times6\) matrices acting on the Voigt components of the deviator. The
    !!   conventional orthotropic parametrisation uses eighteen coefficients, which give
    !!   the criterion its usual name; they are not all independently identifiable.
    !!   Setting both transformations to the identity recovers an isotropic criterion.
    !! - \( a \): exponent, expected \( \ge 1 \), conventionally 6 for body-centred and
    !!   8 for face-centred cubic metals. Values \( 1 < a < 2 \) make the yield function
    !!   non-smooth where a cross difference vanishes, and at the lower limit
    !!   \( a = 1 \) the scalar term is not even once differentiable there (see
    !!   *Numerical treatment*).
    !!
    !! Numerical treatment
    !! -------------------
    !! The criterion is a spectral function of two tensors rather than one, so its
    !! derivatives need the same machinery as [[muscle_yield_cpb06]] applied twice,
    !! plus a coupling block that has no counterpart there:
    !!
    !! 1. **Eigen-decomposition.** `eigen_sym3` (cyclic Jacobi) is used for both
    !!    transformed tensors, avoiding the analytical root formula and its division by
    !!    quantities that vanish on repeated roots. The same solver serves the value and
    !!    both derivatives.
    !! 2. **Repeated roots.** The Daleckii-Krein coefficients are divided differences
    !!    taken *within* each spectrum. They are computed from the cross differences,
    !!    which works because subtracting two of them with a common index leaves exactly
    !!    the gap between two principal values of the same tensor:
    !!    \( \delta_{ij} - \delta_{kj} = \lambda'_i - \lambda'_k \). Nearly equal
    !!    arguments fall back to a Taylor expansion, since the direct quotient loses
    !!    most of its significant digits to cancellation there.
    !! 3. **The kink at \( \lambda'_i = \lambda''_j \).** Where a cross difference
    !!    vanishes, \( |\delta|^a \) fails to be twice differentiable for
    !!    \( 1 < a < 2 \); at \( a = 2 \) it is the smooth function \( \delta^2 \). The
    !!    implementation keeps the previous convention, which agrees with the general
    !!    formula at \( a = 2 \).
    !!
    !! Limitations
    !! -----------
    !! Both derivatives return zero when \( \Phi \le 0 \). That requires every one of
    !! the nine cross differences to vanish, which happens only if the two transformed
    !! tensors are one and the same multiple of the identity. That zero is a safeguard,
    !! not a computed derivative, and callers should not read it as one.
    !!
    !! Reference
    !! ---------
    !! Barlat, F., Aretz, H., Yoon, J. W., Karabin, M. E., Brem, J. C., & Dick, R. E.
    !! (2005). *Linear transfomation-based anisotropic yield functions.* International
    !! Journal of Plasticity, 21(5), 1009-1039.
    !!
    !! For more information see [[muscle_yield_base]] and [[muscle_math_spectral_derivs]]

    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_base
    implicit none
    private

    real(real64), parameter :: YLD_ZERO_REL = 1.0D-12
        !! |lambda'_i - lambda''_j| <= YLD_ZERO_REL*max|lambda| is treated as zero
        !! (previous second-derivative convention for a <= 2). Relative, not absolute.
        !!
        !! The threshold is relative because the principal values scale with the stress
        !! level. A fixed absolute tolerance would change which cross differences count
        !! as zero when the same physical state is expressed in different units, since
        !! the numbers shrink by six orders of magnitude on passing from pascals to
        !! megapascals.
    real(real64), parameter :: YLD_DD_REL = 1.0D-3
        !! Relative gap below which divided differences use their Taylor expansion
        !!
        !! Above this gap the direct quotient is accurate; below it, cancellation in the
        !! numerator generally costs more digits than the expansion's truncation error.
        !! The threshold is a switching point, not a uniform accuracy bound: see yld_dpow.
    real(real64), parameter :: VOIGT_W(6) = [1.0D0, 1.0D0, 1.0D0, 2.0D0, 2.0D0, 2.0D0]
        !! Tensorial double-contraction weights of the (xx,yy,zz,xy,yz,xz) storage
        !!
        !! The three shear components stand for two equal off-diagonal entries each, so
        !! the double contraction of two stored tensors counts them twice. The weights
        !! make that inner product explicit, which is what the adjoint of a
        !! transformation has to respect.

    public :: Yld2004

    type, extends(Base_yield_critera) :: Yld2004
        !! Barlat Yld2004-18p anisotropic yield criterion.
        !!
        !! Coefficients are stored as supplied; the type performs no admissibility check
        !! and does not verify that the two transformations are orthotropic.
        real(real64) :: a
            !! Exponent of the yield function.
        real(real64), dimension(6,6) :: C_prime
            !! First linear transformation, acting on the Voigt components of the deviator.
        real(real64), dimension(6,6) :: C_dprime
            !! Second linear transformation, acting on the Voigt components of the deviator.
    contains
        procedure :: init
        procedure :: stress_eq
        procedure :: dstressEq_dstress => dstressEq_dstress_yld
        procedure :: ddstressEq_ddstress => ddstressEq_ddstress_yld
    end type Yld2004

contains

    subroutine init(self, C_prime, C_dprime, a)
        !! Stores the two linear transformations and the exponent.
        !!
        !! The matrices are taken as given, with no structural assumption: a full
        !! \(6\times6\) is accepted even though the orthotropic case fills only the
        !! normal block and the shear diagonal.
        implicit none
        class(Yld2004), intent(inout) :: self
        real(real64), intent(in), dimension(6,6) :: C_prime
            !! First transformation.
        real(real64), intent(in), dimension(6,6) :: C_dprime
            !! Second transformation.
        real(real64), intent(in) :: a
            !! Exponent.

        self%C_prime = C_prime
        self%C_dprime = C_dprime
        self%a = a
    end subroutine init

    pure function stress_eq(self, stress) result(res)
        !! Returns the Yld2004-18p equivalent stress of a given stress state.
        !!
        !! The deviator is transformed twice, both spectra are extracted, and the
        !! double sum over the nine cross differences is assembled and normalised.
        use muscle_tensors
        use muscle_math_spectral_derivs, only: eigen_sym3
        implicit none
        class(Yld2004), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
            !! Cauchy stress, expressed in the orthotropic axes of the material.
        real(real64) :: res
            !! Equivalent stress. It equals the uniaxial yield stress only for
            !! transformations calibrated to that reference; see the module header.

        type(ten_3D2Osym) :: s
        type(ten_3D2Osym) :: s_prime
        type(ten_3D2Osym) :: s_dprime
        real(real64), dimension(6) :: s_prime_vals
        real(real64), dimension(6) :: s_dprime_vals
        real(real64), dimension(3) :: eig_prime
        real(real64), dimension(3) :: eig_dprime
        real(real64), dimension(3,3) :: vec_prime, vec_dprime
        real(real64) :: sum_val
        integer :: i, j

        ! 1. Calculate deviatoric stress tensor
        s = .dev. stress

        ! 2. Apply linear transformations to Voigt components
        s_prime_vals = matmul(self%C_prime, s%vals)
        s_dprime_vals = matmul(self%C_dprime, s%vals)

        call s_prime%init(s_prime_vals)
        call s_dprime%init(s_dprime_vals)

        ! 3. Compute eigenvalues (same eigen-solver as the derivatives; no division
        !    by a vanishing quantity and no absolute tolerance)
        call eigen_sym3(s_prime, eig_prime, vec_prime)
        call eigen_sym3(s_dprime, eig_dprime, vec_dprime)

        ! 4. Evaluate Barlat Yld2004-18D double summation
        sum_val = 0.0D0
        do i = 1, 3
            do j = 1, 3
                sum_val = sum_val + (abs(eig_prime(i) - eig_dprime(j)))**self%a
            end do
        end do

        res = (0.25D0 * sum_val)**(1.0D0 / self%a)
    end function stress_eq

    ! ----------------------------------------------------------------------------
    ! Linear maps and their adjoints.
    !
    ! L = C : dev is the composition applied to the stress; its adjoint L* is needed
    ! to pull derivatives taken with respect to S' or S'' back to derivatives with
    ! respect to the stress. Both transformations use the same pair of routines,
    ! which therefore take the matrix as an argument rather than reading it from the
    ! type.
    ! ----------------------------------------------------------------------------

    pure function apply_C(C, s) result(res)
        ! Applies a transformation to the Voigt components of a symmetric tensor.
        implicit none
        real(real64), intent(in), dimension(6,6) :: C
        type(ten_3D2Osym), intent(in) :: s
        type(ten_3D2Osym) :: res

        res%vals = matmul(C, s%vals)
    end function apply_C

    pure function apply_CT(C, w) result(res)
        ! Adjoint of s -> C s under the tensorial double contraction:
        ! C* = W^-1 C^T W. Equal to C^T for the Barlat block form (no normal-shear
        ! or cross-shear coupling); init accepts any 6x6, so the weights are kept.
        !
        ! The weights cancel whenever C commutes with W, which covers the orthotropic
        ! block form and, since the three shear weights are equal, any map that keeps
        ! normal and shear components separate. Dropping them can give a silently wrong
        ! adjoint once normal and shear are coupled.
        implicit none
        real(real64), intent(in), dimension(6,6) :: C
        type(ten_3D2Osym), intent(in) :: w
        type(ten_3D2Osym) :: res
        integer :: i

        do i = 1, 6
            res%vals(i) = sum(C(:,i)*VOIGT_W*w%vals)/VOIGT_W(i)
        end do
    end function apply_CT

    pure function apply_L(C, stress) result(res)
        ! S = C : dev(stress).
        implicit none
        real(real64), intent(in), dimension(6,6) :: C
        type(ten_3D2Osym), intent(in) :: stress
        type(ten_3D2Osym) :: res

        res = apply_C(C, .dev. stress)
    end function apply_L

    pure function apply_LT(C, w) result(res)
        ! Adjoint of apply_L: the two operations are applied in reverse order.
        implicit none
        real(real64), intent(in), dimension(6,6) :: C
        type(ten_3D2Osym), intent(in) :: w
        type(ten_3D2Osym) :: res

        res = .dev. apply_CT(C, w)
    end function apply_LT

    ! ----------------------------------------------------------------------------
    ! Scalar building blocks of the yield function.
    !
    ! Every term of Phi is |delta|**a in the cross difference delta = lambda'_i -
    ! lambda''_j, so all derivatives reduce to derivatives of that single scalar
    ! function. These helpers isolate it, including the two places where a naive
    ! evaluation breaks down: the kink of |delta| at the origin, and the
    ! cancellation in divided differences of nearly equal arguments.
    ! ----------------------------------------------------------------------------

    pure function yld_sgn(x) result(s)
        ! Sign with the convention sgn(0) = 0, which keeps t odd at the origin.
        implicit none
        real(real64), intent(in) :: x
        real(real64) :: s

        if (x > 0.0D0) then
            s = 1.0D0
        else if (x < 0.0D0) then
            s = -1.0D0
        else
            s = 0.0D0
        end if
    end function yld_sgn

    pure function yld_t(delta, a) result(res)
        ! t(delta) = |delta|^(a-1) * sgn(delta);  d/d(delta) |delta|^a = a*t
        implicit none
        real(real64), intent(in) :: delta, a
        real(real64) :: res

        res = abs(delta)**(a - 1.0D0)*yld_sgn(delta)
    end function yld_t

    pure function yld_dabs(delta, a, zero_tol) result(res)
        implicit none
        real(real64), intent(in) :: delta, a, zero_tol
        real(real64) :: res

        ! d/d(delta) of |delta|^(a-1) * sgn(delta) = (a-1) |delta|^(a-2).
        ! Within |delta| <= zero_tol (relative to the spectra) keep the previous
        ! convention: 0 for a > 2, 1 otherwise (exact for a = 2).
        !
        ! For a > 2 the limit at the origin really is zero, so the branch is exact
        ! there; within the tolerance, small nonzero arguments are approximated by that
        ! limit. For a = 2 the expression collapses to 1 for any delta, so the branch
        ! agrees with the general formula. For 1 < a < 2 no finite value is correct,
        ! since the derivative diverges; the 1 returned is a convention, as it is at the
        ! contract limit a = 1, where sign(delta) has no derivative at the origin.
        if (abs(delta) > zero_tol) then
            res = (a - 1.0D0)*abs(delta)**(a - 2.0D0)
        else if (a > 2.0D0) then
            res = 0.0D0
        else
            res = 1.0D0
        end if
    end function yld_dabs

    pure function yld_dpow(u, v, p) result(res)
        ! (u**p - v**p)/(u - v) for u, v > 0; Taylor expansion about the midpoint
        ! when the relative gap is small (the direct quotient loses digits there).
        !
        ! The expansion is written in the squared relative half-gap r, so only even
        ! orders survive; it retains terms through r**2. Its truncation error depends on
        ! both p and the gap, and not monotonically in p: it vanishes outright at small
        ! positive integer p, where the series terminates. YLD_DD_REL is a switching
        ! point, not a uniform accuracy bound.
        implicit none
        real(real64), intent(in) :: u, v, p
        real(real64) :: res, m, r

        if (abs(u - v) > YLD_DD_REL*max(u, v)) then
            res = (u**p - v**p)/(u - v)
        else
            m = 0.5D0*(u + v)
            r = (0.5D0*(u - v)/m)**2
            res = p*m**(p - 1.0D0)*(1.0D0 + (p - 1.0D0)*(p - 2.0D0)*r/6.0D0 &
                  *(1.0D0 + (p - 3.0D0)*(p - 4.0D0)*r/20.0D0))
        end if
    end function yld_dpow

    pure function yld_divdiff(x, y, a, zero_tol) result(res)
        ! (t(x) - t(y))/(x - y), with its limit t'(x) at x = y.
        !
        ! Four cases, chosen so that no branch ever divides by a small difference:
        ! both arguments at the origin, both positive, both negative (where the
        ! oddness of t reduces the problem to the positive case), and opposite signs,
        ! where |x - y| is at least max(|x|,|y|) and the quotient is well conditioned.
        implicit none
        real(real64), intent(in) :: x, y, a, zero_tol
        real(real64) :: res

        if (abs(x) <= zero_tol .and. abs(y) <= zero_tol) then
            res = yld_dabs(0.0D0, a, zero_tol)
        else if (x > 0.0D0 .and. y > 0.0D0) then
            res = yld_dpow(x, y, a - 1.0D0)
        else if (x < 0.0D0 .and. y < 0.0D0) then
            ! t is odd: t(x) = -(-x)**(a-1)
            res = yld_dpow(-x, -y, a - 1.0D0)
        else
            ! Opposite signs (or one zero): |x - y| >= max(|x|,|y|), well conditioned
            res = (yld_t(x, a) - yld_t(y, a))/(x - y)
        end if
    end function yld_divdiff

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

    pure function dstressEq_dstress_yld(self, stress) result(res)
        !! Returns the gradient of the equivalent stress, which under an associated
        !! flow rule gives the direction of plastic flow.
        !!
        !! The derivative is taken with respect to each transformed tensor separately
        !! and the two contributions are pulled back and added. Each one is diagonal in
        !! its own eigenbasis, so only three numbers per spectrum are needed.
        !!
        !! The second contribution carries the opposite sign, because a cross difference
        !! enters with a minus in front of the second principal value.
        !!
        !! Repeated principal values need no special treatment here: the diagonal
        !! entries coincide on a repeated root, so the result does not depend on which
        !! eigenbasis the solver returns for the cluster.
        use muscle_tensors
        use muscle_math_spectral_derivs, only: eigen_sym3, spectral_compose
        implicit none
        class(Yld2004), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
            !! Cauchy stress, expressed in the orthotropic axes of the material.
        type(ten_3D2Osym) :: res
            !! Gradient of the equivalent stress. Zero in the degenerate case
            !! \( \Phi \le 0 \); that value is a safeguard, not a computed derivative.

        type(ten_3D2Osym) :: sp, sd, dphi_p, dphi_d
        real(real64) :: ep(3), ed(3), vp(3,3), vd(3,3), gp(3,3), gd(3,3)
        real(real64) :: a, Phi, seq, psi
        integer :: i, j

        a = self%a
        sp = apply_L(self%C_prime, stress)
        sd = apply_L(self%C_dprime, stress)
        call eigen_sym3(sp, ep, vp)
        call eigen_sym3(sd, ed, vd)

        Phi = 0.0D0
        do i = 1, 3
            do j = 1, 3
                Phi = Phi + abs(ep(i) - ed(j))**a
            end do
        end do
        if (Phi <= 0.0D0) then
            res = 0.0D0
            return
        end if
        seq = (0.25D0*Phi)**(1.0D0/a)

        ! dseq/dS' = sum_i gp_i v'_i (x) v'_i,  dseq/dS'' = sum_j gd_j v''_j (x) v''_j.
        ! gp (gd) is symmetric in the repeated roots of S' (S''), so any
        ! eigenbasis of a cluster gives the same tensor.
        gp = 0.0D0
        do i = 1, 3
            psi = 0.0D0
            do j = 1, 3
                psi = psi + yld_t(ep(i) - ed(j), a)
            end do
            gp(i,i) = (seq/Phi)*psi
        end do

        gd = 0.0D0
        do j = 1, 3
            psi = 0.0D0
            do i = 1, 3
                psi = psi + yld_t(ep(i) - ed(j), a)
            end do
            gd(j,j) = -(seq/Phi)*psi
        end do

        dphi_p = spectral_compose(vp, gp)
        dphi_d = spectral_compose(vd, gd)
        res = apply_LT(self%C_prime, dphi_p) + apply_LT(self%C_dprime, dphi_d)
    end function dstressEq_dstress_yld

    pure function ddstressEq_ddstress_yld(self, stress) result(res)
        !! Returns the second derivative of the equivalent stress, required by the
        !! consistent tangent operator of implicit return-mapping schemes.
        !!
        !! With two transformations the Daleckii-Krein structure appears twice and a
        !! third block joins it. The diagonal entries split into three couplings:
        !! within the first spectrum, within the second, and between the two. That last
        !! block is what distinguishes this criterion from a single-transformation one,
        !! and it is where the direct \( |\delta|^{a-2} \) term enters with a minus
        !! sign, since the two principal values sit on opposite sides of the difference.
        !!
        !! The off-diagonal entries remain divided differences taken within one
        !! spectrum, computed from the cross differences so that the denominator is
        !! exactly the gap between two principal values of the same tensor.
        !!
        !! The Hessian is assembled one column at a time: each Voigt basis tensor is
        !! pushed forward through both transformations, rotated into each eigenbasis,
        !! multiplied blockwise, and pulled back.
        use muscle_tensors
        use muscle_math_spectral_derivs, only: eigen_sym3, spectral_rotate, spectral_compose
        implicit none
        class(Yld2004), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
            !! Cauchy stress, expressed in the orthotropic axes of the material.
        type(ten_3D4O3sym) :: res
            !! Second derivative of the equivalent stress. Zero in the degenerate case
            !! \( \Phi \le 0 \).

        type(ten_3D2Osym) :: sp, sd, e
        type(ten_3D2Osym) :: col(6)
        real(real64) :: ep(3), ed(3), vp(3,3), vd(3,3)
        real(real64) :: gp(3), gd(3), delta(3,3), tt(3,3), dt(3,3)
        real(real64) :: kpp(3,3), kdd(3,3), kpd(3,3), thp(3,3), thd(3,3)
        real(real64) :: zp(3,3), zd(3,3), mp(3,3), md(3,3)
        real(real64) :: a, Phi, seq, c, zero_tol
        integer :: i, j, k, l

        a = self%a
        sp = apply_L(self%C_prime, stress)
        sd = apply_L(self%C_dprime, stress)
        call eigen_sym3(sp, ep, vp)
        call eigen_sym3(sd, ed, vd)

        Phi = 0.0D0
        do i = 1, 3
            do j = 1, 3
                Phi = Phi + abs(ep(i) - ed(j))**a
            end do
        end do
        if (Phi <= 0.0D0) then
            res = 0.0D0
            return
        end if
        seq = (0.25D0*Phi)**(1.0D0/a)
        c = seq/Phi
        zero_tol = YLD_ZERO_REL*max(maxval(abs(ep)), maxval(abs(ed)))

        ! delta holds the nine cross differences; tt and dt are the first and second
        ! derivatives of |delta|^a evaluated on them, up to the factor a.
        do j = 1, 3
            do i = 1, 3
                delta(i,j) = ep(i) - ed(j)
                tt(i,j) = yld_t(delta(i,j), a)
                dt(i,j) = yld_dabs(delta(i,j), a, zero_tol)
            end do
        end do
        do i = 1, 3
            gp(i) = c*sum(tt(i,:))
            gd(i) = -c*sum(tt(:,i))
        end do

        ! dg/dlambda: rank-one terms (1-a) g g / seq plus |delta|^(a-2) couplings
        do k = 1, 3
            do i = 1, 3
                kpp(i,k) = (1.0D0 - a)*gp(i)*gp(k)/seq
                kdd(i,k) = (1.0D0 - a)*gd(i)*gd(k)/seq
                kpd(i,k) = (1.0D0 - a)*gp(i)*gd(k)/seq - c*dt(i,k)
            end do
        end do
        do i = 1, 3
            kpp(i,i) = kpp(i,i) + c*sum(dt(i,:))
            kdd(i,i) = kdd(i,i) + c*sum(dt(:,i))
        end do

        ! Daleckii-Krein coefficients, evaluated without cancellation:
        !   (gp_i - gp_k)/(lambda'_i - lambda'_k)   = c*sum_j t[delta_ij, delta_kj]
        !   (gd_j - gd_l)/(lambda''_j - lambda''_l) = c*sum_i t[delta_ij, delta_il]
        ! On a repeated root they tend to the diagonal h' terms above.
        !
        ! The identity that makes this work is delta_ij - delta_kj = lambda'_i -
        ! lambda'_k: differencing the cross differences over a common second index
        ! leaves precisely the gap within the first spectrum, and symmetrically for
        ! the second. No subtraction of principal values is ever formed explicitly.
        thp = 0.0D0
        thd = 0.0D0
        do k = 2, 3
            do i = 1, k - 1
                do j = 1, 3
                    thp(i,k) = thp(i,k) + yld_divdiff(delta(i,j), delta(k,j), a, zero_tol)
                    thd(i,k) = thd(i,k) + yld_divdiff(delta(j,i), delta(j,k), a, zero_tol)
                end do
                thp(i,k) = c*thp(i,k)
                thd(i,k) = c*thd(i,k)
                thp(k,i) = thp(i,k)
                thd(k,i) = thd(i,k)
            end do
        end do

        ! Pull back column by column: H_sigma : e_J
        !   d(dseq/dS')  : E = V' Mp V'^T,  d(dseq/dS'') : E = V'' Md V''^T
        do l = 1, 6
            e = 0.0D0
            e%vals(l) = 1.0D0
            zp = spectral_rotate(vp, apply_L(self%C_prime, e))
            zd = spectral_rotate(vd, apply_L(self%C_dprime, e))
            mp = thp*zp
            md = thd*zd
            do i = 1, 3
                mp(i,i) = 0.0D0
                md(i,i) = 0.0D0
                do k = 1, 3
                    mp(i,i) = mp(i,i) + kpp(i,k)*zp(k,k) + kpd(i,k)*zd(k,k)
                    md(i,i) = md(i,i) + kdd(i,k)*zd(k,k) + kpd(k,i)*zp(k,k)
                end do
            end do
            col(l) = apply_LT(self%C_prime, spectral_compose(vp, mp)) &
                   + apply_LT(self%C_dprime, spectral_compose(vd, md))
        end do
        res = pack_hess_from_cols(col)
    end function ddstressEq_ddstress_yld

end module muscle_yield_yld2004
