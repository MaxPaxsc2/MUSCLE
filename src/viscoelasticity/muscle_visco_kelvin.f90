! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_visco_kelvin
    !! Module muscle_visco_kelvin
    !! ==========================
    !! Isotropic Kelvin-Voigt solid: Hookean spring and Newtonian dashpot in parallel,
    !!
    !! \[ \sigma = K\,\mathrm{tr}(\varepsilon)\,\mathbf{I} + 2G\,\mathrm{dev}(\varepsilon)
    !!          + \zeta\,\mathrm{tr}(\dot\varepsilon)\,\mathbf{I}
    !!          + 2\mu_v\,\mathrm{dev}(\dot\varepsilon) \]
    !! (Picard and Fafard (2011), Eqs. (47)-(48), isotropic case.)
    !!
    !! With the strain varying linearly over the step, \(\dot\varepsilon = \Delta\varepsilon/\Delta t\),
    !! the stress is affine in the elastic strain: \(K^{alg} = K + \zeta/\Delta t\),
    !! \(G^{alg} = G + \mu_v/\Delta t\) and \(\sigma^{hist}_n\) depends on \(\varepsilon_n\) alone
    !! (see Theory/viscoelasticity/generalized_maxwell.tex, pending Kelvin-Voigt item).
    !!
    !! Coefficients: `weight(1) = zeta/dt`, `weight(2) = 2 mu_v/dt`; `decay` is unused.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_base, only : Base_viscoelasticity, Visco_coefficients
    implicit none
    private

    public :: Kelvin_voigt

    type, extends(Base_viscoelasticity) :: Kelvin_voigt
        real(real64), private :: k = 0D0      !! Elastic bulk modulus K
        real(real64), private :: g = 0D0      !! Elastic shear modulus G
        real(real64), private :: zeta = 0D0   !! Volumetric viscosity of the dashpot zeta
        real(real64), private :: mu_v = 0D0   !! Deviatoric viscosity of the dashpot mu_v
    contains
        procedure :: set_parameters
        procedure :: n_branches => n_branches_kelvin
        procedure :: coefficients => coefficients_kelvin
        procedure :: freeze => freeze_kelvin
        procedure :: update => update_kelvin
    end type Kelvin_voigt

contains

    pure subroutine set_parameters(self, bulk, shear, bulk_viscosity, shear_viscosity)
        !! Sets the elastic moduli and the viscosities of the solid.
        class(Kelvin_voigt), intent(inout) :: self
        real(real64), intent(in) :: bulk             !! Bulk modulus K
        real(real64), intent(in) :: shear            !! Shear modulus G
        real(real64), intent(in) :: bulk_viscosity   !! Volumetric viscosity zeta
        real(real64), intent(in) :: shear_viscosity  !! Deviatoric viscosity mu_v

        if (bulk < 0D0 .or. shear < 0D0) then
            error stop "ERROR FATAL [Kelvin_voigt%set_parameters]: " // &
                       "Los modulos 'bulk' y 'shear' no pueden ser negativos."
        end if
        if (bulk_viscosity < 0D0 .or. shear_viscosity < 0D0) then
            error stop "ERROR FATAL [Kelvin_voigt%set_parameters]: " // &
                       "Las viscosidades no pueden ser negativas."
        end if

        self%k = bulk
        self%g = shear
        self%zeta = bulk_viscosity
        self%mu_v = shear_viscosity
    end subroutine set_parameters

    pure function n_branches_kelvin(self) result(res)
        !! No history tensors: the elastic strain at t_n is the whole history of the model.
        class(Kelvin_voigt), intent(in) :: self
        integer :: res

        res = 0
    end function n_branches_kelvin

    pure function coefficients_kelvin(self, dt) result(coefs)
        !! Algorithmic moduli over the step, \(K^{alg} = K + \zeta/\Delta t\) and
        !! \(G^{alg} = G + \mu_v/\Delta t\), and the two viscous factors used by `freeze`.
        !! Evaluated once per step.
        class(Kelvin_voigt), intent(in) :: self
        real(real64), intent(in)        :: dt
        type(Visco_coefficients)        :: coefs

        if (dt <= 0D0) then
            error stop "ERROR FATAL [Kelvin_voigt%coefficients]: " // &
                       "El paso de tiempo 'dt' debe ser positivo: el modelo es " // &
                       "infinitamente rigido para dt <= 0."
        end if

        coefs%dt = dt
        coefs%k_alg = self%k + self%zeta/dt
        coefs%g_alg = self%g + self%mu_v/dt
        coefs%weight(1) = self%zeta/dt
        coefs%weight(2) = 2D0*self%mu_v/dt
    end function coefficients_kelvin

    pure function freeze_kelvin(self, coefs, history) result(elas)
        !! History stress \(\sigma^{hist}_n = -(\zeta/\Delta t)\,\mathrm{tr}(\varepsilon_n)\,\mathbf{I}
        !! - (2\mu_v/\Delta t)\,\mathrm{dev}(\varepsilon_n)\), constant over the step.
        class(Kelvin_voigt), intent(in)        :: self
        type(Visco_coefficients), intent(in)   :: coefs
        type(Viscoelastic_history), intent(in) :: history
        type(Elasticity_frozen)                :: elas
        type(ten_3D2Osym)                      :: e_n, hist

        e_n = history%state_n%strain_e
        hist = -coefs%weight(2)*(.dev. e_n) - (coefs%weight(1)*(e_n .ddot. iden_2O()))*iden_2O()
        call elas%set(coefs%k_alg, coefs%g_alg, hist)
    end function freeze_kelvin

    pure subroutine update_kelvin(self, coefs, strain_e, history)
        !! Stores the elastic strain of the converged step; there are no branch tensors.
        class(Kelvin_voigt), intent(in)           :: self
        type(Visco_coefficients), intent(in)      :: coefs
        type(ten_3D2Osym), intent(in)             :: strain_e
        type(Viscoelastic_history), intent(inout) :: history

        history%state_np1%n_branches = self%n_branches()
        history%state_np1%strain_e = strain_e
    end subroutine update_kelvin

end module muscle_visco_kelvin
