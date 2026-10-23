! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_visco_frozen
    !! Module muscle_visco_frozen
    !! ==========================
    !! Algorithmic ("frozen") elasticity of a linear viscoelastic model over one time step.
    !!
    !! Integrated over a step, every linear viscoelastic model is affine in the elastic strain:
    !!
    !! \[
    !!   \boldsymbol{\sigma}_{n+1} = \mathbb{C}^{alg} : \boldsymbol{\varepsilon}^e_{n+1}
    !!                             + \boldsymbol{\sigma}^{hist}_n
    !! \]
    !!
    !! where \( \mathbb{C}^{alg} \) depends only on the time step and \( \boldsymbol{\sigma}^{hist}_n \)
    !! only on the converged state t_n. Both are constant during the Newton iterations, so this
    !! object is an ordinary `Base_elasticity` and can be handed to any return mapping solver
    !! without changes (Delta sigma = -C_alg : Delta eps_p still holds).
    !!
    !! The isotropic case is stored with two scalars (K_alg, G_alg): stress and tangent are
    !! evaluated without fourth-order products.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_elasticity_base
    implicit none
    private

    public :: Elasticity_frozen

    type, extends(Base_elasticity) :: Elasticity_frozen
        real(real64), private      :: k_alg = 0D0  !! Algorithmic bulk modulus
        real(real64), private      :: g_alg = 0D0  !! Algorithmic shear modulus
        type(ten_3D2Osym), private :: stress_hist = ten_3D2Osym((/0d0, 0d0, 0d0, 0d0, 0d0, 0d0/))
            !! History stress (constant over the step)
    contains
        procedure :: set
        procedure :: bulk
        procedure :: shear
        procedure :: history_stress
        procedure :: stress_3D              => stress_frozen_3D
        procedure :: stress_2D              => stress_frozen_2D
        procedure :: stress_dev_3D          => stress_dev_frozen_3D
        procedure :: dstress_dstrain_3D     => dstress_dstrain_frozen_3D
        procedure :: dstress_dstrain_2D     => dstress_dstrain_frozen_2D
        procedure :: dstress_dstrain_dev_3D => dstress_dstrain_dev_frozen_3D
        procedure, public :: to_string => to_string_elasticity_frozen
    end type Elasticity_frozen

contains

    pure subroutine set(self, k_alg, g_alg, stress_hist)
        class(Elasticity_frozen), intent(inout) :: self
        real(real64), intent(in)                :: k_alg, g_alg
        type(ten_3D2Osym), intent(in)           :: stress_hist
        self%k_alg = k_alg
        self%g_alg = g_alg
        self%stress_hist = stress_hist
    end subroutine set

    pure function bulk(self) result(res)
        class(Elasticity_frozen), intent(in) :: self
        real(real64) :: res
        res = self%k_alg
    end function bulk

    pure function shear(self) result(res)
        class(Elasticity_frozen), intent(in) :: self
        real(real64) :: res
        res = self%g_alg
    end function shear

    pure function history_stress(self) result(res)
        class(Elasticity_frozen), intent(in) :: self
        type(ten_3D2Osym) :: res
        res = self%stress_hist
    end function history_stress

    ! =========================================================================
    ! Stress
    ! =========================================================================
    pure function stress_frozen_3D(self, strain) result(res)
        class(Elasticity_frozen), intent(in) :: self
        class(ten_3D2Osym), intent(in)       :: strain
        type(ten_3D2Osym) :: res
        real(real64) :: two_g, p

        two_g = 2D0*self%g_alg
        p = (self%k_alg - two_g/3D0)*(strain%xx() + strain%yy() + strain%zz())
        res%vals(1:3) = two_g*strain%vals(1:3) + p + self%stress_hist%vals(1:3)
        res%vals(4:6) = two_g*strain%vals(4:6) + self%stress_hist%vals(4:6)
    end function stress_frozen_3D

    pure function stress_frozen_2D(self, strain) result(res)
        class(Elasticity_frozen), intent(in) :: self
        class(ten_2D2Osym), intent(in)       :: strain
        type(ten_2D2Osym) :: res
        real(real64) :: two_g, p

        two_g = 2D0*self%g_alg
        p = (self%k_alg - two_g/3D0)*(strain%xx() + strain%yy() + strain%zz())
        call res%init(xx = two_g*strain%xx() + p + self%stress_hist%xx(), &
                      yy = two_g*strain%yy() + p + self%stress_hist%yy(), &
                      zz = two_g*strain%zz() + p + self%stress_hist%zz(), &
                      xy = two_g*strain%xy() + self%stress_hist%xy())
    end function stress_frozen_2D

    pure function stress_dev_frozen_3D(self, strain) result(res)
        class(Elasticity_frozen), intent(in) :: self
        class(ten_3D2Osym), intent(in)       :: strain
        type(ten_3D2Osym) :: res
        res = 2D0*self%g_alg*(.dev. strain) + (.dev. self%stress_hist)
    end function stress_dev_frozen_3D

    ! =========================================================================
    ! Tangent (constant over the step)
    ! =========================================================================
    pure function dstress_dstrain_frozen_3D(self, strain) result(res)
        class(Elasticity_frozen), intent(in) :: self
        class(ten_3D2Osym), intent(in)       :: strain
        type(ten_3D4O3sym) :: res
        real(real64) :: e1, lam, mu

        mu  = self%g_alg
        lam = self%k_alg - 2D0*mu/3D0
        e1  = lam + 2D0*mu
        call res%init( xxxx=e1,    yyyy=e1,    zzzz=e1,     &
                       xxyy=lam,   yyzz=lam,   xxzz=lam,    &
                       xxxy=0.0D0, xxyz=0.0D0, xxxz=0.0D0,  &
                       yyxy=0.0D0, yyyz=0.0D0, yyxz=0.0D0,  &
                       zzxy=0.0D0, zzyz=0.0D0, zzxz=0.0D0,  &
                       xyxy=mu,    yzyz=mu,    xzxz=mu,     &
                       xyyz=0.0D0, yzxz=0.0D0, xyxz=0.0D0   )
    end function dstress_dstrain_frozen_3D

    pure function dstress_dstrain_frozen_2D(self, strain) result(res)
        class(Elasticity_frozen), intent(in) :: self
        class(ten_2D2Osym), intent(in)       :: strain
        type(ten_2D4O3sym) :: res
        real(real64) :: e1, lam, mu

        mu  = self%g_alg
        lam = self%k_alg - 2D0*mu/3D0
        e1  = lam + 2D0*mu
        call res%init( xxxx=e1,    yyyy=e1,    zzzz=e1,    &
                       xxyy=lam,   yyzz=lam,   xxzz=lam,   &
                       xxxy=0.0D0, yyxy=0.0D0, zzxy=0.0D0, &
                       xyxy=mu                             )
    end function dstress_dstrain_frozen_2D

    pure function dstress_dstrain_dev_frozen_3D(self, strain) result(res)
        class(Elasticity_frozen), intent(in) :: self
        class(ten_3D2Osym), intent(in)       :: strain
        type(ten_3D4O3sym) :: res
        real(real64) :: d, o, mu

        mu = self%g_alg
        d  =  4D0*mu/3D0
        o  = -2D0*mu/3D0
        call res%init( xxxx=d,     yyyy=d,     zzzz=d,      &
                       xxyy=o,     yyzz=o,     xxzz=o,      &
                       xxxy=0.0D0, xxyz=0.0D0, xxxz=0.0D0,  &
                       yyxy=0.0D0, yyyz=0.0D0, yyxz=0.0D0,  &
                       zzxy=0.0D0, zzyz=0.0D0, zzxz=0.0D0,  &
                       xyxy=mu,    yzyz=mu,    xzxz=mu,     &
                       xyyz=0.0D0, yzxz=0.0D0, xyxz=0.0D0   )
    end function dstress_dstrain_dev_frozen_3D

    pure function to_string_elasticity_frozen(self) result(str)
        !! Generates a formatted string representation of Elasticity_frozen's internal state.
        class(Elasticity_frozen), intent(in) :: self
        character(len=275)                   :: str
        character(len=1), parameter          :: nl = new_line('A')

        write(str, "(A, A, A, ES14.6, A, A, ES14.6, A, A, 6ES11.3, A)") &
            "=== Elasticity_frozen Internal State ===", nl, &
            "  K_alg (Bulk)          : ", self%k_alg, nl, &
            "  G_alg (Shear)         : ", self%g_alg, nl, &
            "  Stress hist           : ", self%stress_hist%vals, nl
    end function to_string_elasticity_frozen

end module muscle_visco_frozen
