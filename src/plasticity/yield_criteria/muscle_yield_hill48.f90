! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_yield_hill48
    !! Module muscle_yield_hill48
    !! ==========================
    !! Implements the Hill (1948) quadratic yield criterion for orthotropic metals,
    !! together with its analytical first and second derivatives with respect to stress.
    !!
    !! Mathematical Formulation
    !! ------------------------
    !! The equivalent stress is the square root of a quadratic form in the stress components:
    !! \[ \bar{\sigma} = \left[ F(\sigma_{yy}-\sigma_{zz})^2 + G(\sigma_{zz}-\sigma_{xx})^2
    !!     + H(\sigma_{xx}-\sigma_{yy})^2 + 2L\sigma_{yz}^2 + 2M\sigma_{xz}^2
    !!     + 2N\sigma_{xy}^2 \right]^{1/2} \]
    !! where \(F, G, H, L, M, N\) are the six anisotropy coefficients and the coordinate
    !! axes are assumed aligned with the orthotropic axes of the material.
    !!
    !! Because the criterion depends on stress only through differences of normal
    !! components and through shear components, it is insensitive to hydrostatic
    !! pressure: \( \bar{\sigma}(\boldsymbol{\sigma} + p\mathbf{I}) = \bar{\sigma}(\boldsymbol{\sigma}) \).
    !!
    !! Setting \( F = G = H = 1/2 \) and \( L = M = N = 3/2 \) recovers the isotropic
    !! von Mises criterion, which is the usual check when validating a parameter set.
    !!
    !! Derivatives
    !! -----------
    !! Writing \( A = \bar{\sigma}^2 \) for the quadratic form, both derivatives follow
    !! from the chain rule applied to \( \bar{\sigma} = A^{1/2} \):
    !! \[ \frac{\partial \bar{\sigma}}{\partial \boldsymbol{\sigma}}
    !!    = \frac{1}{2\bar{\sigma}} \frac{\partial A}{\partial \boldsymbol{\sigma}} \]
    !! \[ \frac{\partial^2 \bar{\sigma}}{\partial \boldsymbol{\sigma}^2}
    !!    = \frac{1}{2\bar{\sigma}} \frac{\partial^2 A}{\partial \boldsymbol{\sigma}^2}
    !!    - \frac{1}{\bar{\sigma}}
    !!      \frac{\partial \bar{\sigma}}{\partial \boldsymbol{\sigma}} \otimes
    !!      \frac{\partial \bar{\sigma}}{\partial \boldsymbol{\sigma}} \]
    !! Since \(A\) is quadratic, \( \partial^2 A / \partial \boldsymbol{\sigma}^2 \) is a
    !! constant tensor, assembled once from the six coefficients.
    !!
    !! Limitations
    !! -----------
    !! Both derivatives divide by \( \bar{\sigma} \) and require it to be strictly
    !! positive. It vanishes on the hydrostatic axis for every parameter set, since the
    !! quadratic form is singular in that direction whatever the coefficients. Further
    !! null modes appear only when the form restricted to deviatoric stresses is itself
    !! singular: with \( L = 0 \), for instance, pure \(yz\) shear also gives
    !! \( \bar{\sigma} = 0 \) although the state is not hydrostatic. Callers must not
    !! evaluate the derivatives at such a state.
    !! Return-mapping algorithms do not normally reach one, because the yield surface
    !! is not crossed there.
    !!
    !! Reference
    !! ---------
    !! Hill, R. (1948). *A theory of the yielding and plastic flow of anisotropic metals.*
    !! Proceedings of the Royal Society A, 193(1033), 281-297.
    !!
    !! For more information see [[muscle_yield_base]]

    use, intrinsic :: iso_fortran_env
    use muscle_yield_base
    implicit none
    private

    public :: Hill48
    type, extends(Base_yield_critera) :: Hill48
        !! Hill (1948) orthotropic quadratic yield criterion.
        !!
        !! The six coefficients are dimensionless and are usually identified either from
        !! directional yield stresses or from Lankford coefficients measured on sheet
        !! specimens. They are stored exactly as supplied: the type applies no
        !! normalisation and performs no admissibility check.
        real(real64) :: f, g, h, l, m, n
            !! Anisotropy coefficients \(F, G, H\) (normal terms) and \(L, M, N\) (shear
            !! terms in the \(yz\), \(xz\) and \(xy\) planes respectively).
    contains
        procedure :: init
        procedure :: stress_eq
        procedure :: dstressEq_dstress => dstressEq_dstress_h48
        procedure :: ddstressEq_ddstress => ddstressEq_ddstress_h48
    end type Hill48

contains

    subroutine init(self, f, g, h, l, m, n)
        !! Stores the six anisotropy coefficients in the criterion.
        !!
        !! No consistency check is performed. A physically meaningful set keeps the
        !! quadratic form positive semi-definite; otherwise `stress_eq` may take the
        !! square root of a negative number.
        implicit none
        class(Hill48), intent(inout) :: self
        real(real64), intent(in) :: f, g, h, l, m, n
            !! Anisotropy coefficients, in the order \(F, G, H, L, M, N\).

        self%f = f
        self%g = g
        self%h = h
        self%l = l
        self%m = m
        self%n = n
    end subroutine init

    pure function stress_eq(self, stress) result(res)
        !! Returns the Hill (1948) equivalent stress of a given stress state.
        use muscle_tensors
        implicit none
        class(Hill48), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
            !! Cauchy stress, expressed in the orthotropic axes of the material.
        real(real64) :: res
            !! Equivalent stress. Zero on the hydrostatic axis.

        res = (self%f*(stress%yy()-stress%zz())**2 + &
               self%g*(stress%zz()-stress%xx())**2 + &
               self%h*(stress%xx()-stress%yy())**2 + &
               2.0D0*self%l*stress%yz()**2 + &
               2.0D0*self%m*stress%xz()**2 + &
               2.0D0*self%n*stress%xy()**2   &
               )**0.5D0

    end function stress_eq

    pure function dstressEq_dstress_h48(self, stress) result(res)
        !! Returns the gradient of the equivalent stress, which under an associated
        !! flow rule gives the direction of plastic flow.
        !!
        !! The quadratic form is differentiated first and the result scaled by
        !! \( 1/(2\bar{\sigma}) \). The gradient is homogeneous of degree zero, so
        !! scaling the stress by a positive factor leaves it unchanged; reversing the
        !! sign of the stress reverses the gradient.
        use muscle_tensors
        implicit none
        class(Hill48), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
            !! Cauchy stress. Must not lie on the hydrostatic axis.
        type(ten_3D2Osym) :: res
            !! Gradient of the equivalent stress, in tensorial (not engineering) shear.

        real(real64) :: seq
        real(real64) :: sxx, syy, szz, sxy, syz, sxz
        type(ten_3D2Osym) :: df

        sxx = stress%xx()
        syy = stress%yy()
        szz = stress%zz()
        sxy = stress%xy()
        syz = stress%yz()
        sxz = stress%xz()

        ! df = d(seq^2)/dstress, tensorial shear (not engineering).
        ! Each normal component collects the two quadratic terms it appears in. The
        ! shear components carry a factor 2 because the off-diagonal pair contributes
        ! twice to the quadratic form while being stored once.
        df%vals(1) = 2.0D0*(self%g*(sxx - szz) + self%h*(sxx - syy))
        df%vals(2) = 2.0D0*(self%f*(syy - szz) + self%h*(syy - sxx))
        df%vals(3) = 2.0D0*(self%f*(szz - syy) + self%g*(szz - sxx))
        df%vals(4) = 2.0D0*self%n*sxy
        df%vals(5) = 2.0D0*self%l*syz
        df%vals(6) = 2.0D0*self%m*sxz

        seq = self%stress_eq(stress)
        res = df/(2.0D0*seq)
    end function dstressEq_dstress_h48

    pure function ddstressEq_ddstress_h48(self, stress) result(res)
        !! Returns the second derivative of the equivalent stress, required by the
        !! consistent tangent operator of implicit return-mapping schemes.
        !!
        !! The subtracted outer product makes the result singular along the stress
        !! direction, a consequence of the equivalent stress being homogeneous of
        !! degree one. Pressure insensitivity supplies a second null direction, the
        !! spherical tensor.
        use muscle_tensors
        implicit none
        class(Hill48), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
            !! Cauchy stress. Must not lie on the hydrostatic axis.
        type(ten_3D4O3sym) :: res
            !! Second derivative of the equivalent stress.

        real(real64) :: seq
        type(ten_3D2Osym) :: grad
        type(ten_3D4O3sym) :: ddf

        seq = self%stress_eq(stress)
        grad = self%dstressEq_dstress(stress)

        ! ddf = d2(seq^2)/dstress2. Voigt 21: (11,22,33,44,55,66,12,23,34,45,56,13,...)
        ! Constant, because seq^2 is quadratic in the stress. At most nine of the
        ! twenty-one components can be non-zero: three normal-normal diagonal terms,
        ! three shear-shear diagonal terms, and three normal-normal couplings. The
        ! couplings are -2H, -2F and -2G, hence negative whenever the corresponding
        ! coefficients are positive, which a general admissible set does not require.
        ! The remaining twelve always vanish, since the criterion couples neither
        ! normal with shear nor one shear plane with another.
        ddf = 0.0D0
        ddf%vals(1)  = 2.0D0*(self%g + self%h)
        ddf%vals(2)  = 2.0D0*(self%f + self%h)
        ddf%vals(3)  = 2.0D0*(self%f + self%g)
        ddf%vals(4)  = self%n
        ddf%vals(5)  = self%l
        ddf%vals(6)  = self%m
        ddf%vals(7)  = -2.0D0*self%h
        ddf%vals(8)  = -2.0D0*self%f
        ddf%vals(12) = -2.0D0*self%g

        ! d2seq/dstress2 = ddf/(2*seq) - (dseq/dstress)⊗(dseq/dstress)/seq
        res = ddf/(2.0D0*seq) - (.tdotsym. grad)/seq
    end function ddstressEq_ddstress_h48

end module muscle_yield_hill48
