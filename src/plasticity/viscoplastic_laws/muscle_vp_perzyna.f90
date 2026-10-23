! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_vp_perzyna
    !! Module muscle_vp_perzyna
    !! =================================
    !!
    !! Perzyna overstress law:
    !!
    !!     g = (1/eta) * <psi>^N,     psi = sigma_eq/Y(ep) - 1
    !!
    !! where `Y(ep)` is the static hardening law (`hard_law`), `eta` the viscosity
    !! (relaxation time) and `N` the rate sensitivity exponent.
    use, intrinsic :: iso_fortran_env, only: real64
    use muscle_vp_overstress_base, only: Base_overstress_law
    use muscle_hard_base, only: Base_hardening_laws
    implicit none
    private
    public :: Perzyna_overstress

    type, extends(Base_overstress_law) :: Perzyna_overstress
        class(Base_hardening_laws), allocatable :: hard_law
        real(real64) :: eta = 1D0  !! Viscosity (time)
        real(real64) :: n = 1D0    !! Rate sensitivity exponent
    contains
        procedure :: flow_rate => flow_rate_perzyna
        procedure :: dflow_dsigma => dflow_dsigma_perzyna
        procedure :: dflow_dep => dflow_dep_perzyna
    end type Perzyna_overstress

    contains

        pure function dflow_dpsi(self, psi) result(res)
        !! Derivative of `g` with respect to `psi`.
            class(Perzyna_overstress), intent(in) :: self
            real(real64), intent(in) :: psi
            real(real64) :: res
            if (psi .le. 0D0) then
                res = 0D0
            else
                res = self%n*psi**(self%n - 1D0)/self%eta
            end if
        end function dflow_dpsi

        pure function flow_rate_perzyna(self, sigma_eq, ep) result(res)
            class(Perzyna_overstress), intent(in) :: self
            real(real64), intent(in) :: sigma_eq, ep
            real(real64) :: res
            real(real64) :: psi

            psi = sigma_eq/self%hard_law%stress(ep) - 1D0
            if (psi .le. 0D0) then
                res = 0D0
            else
                res = psi**self%n/self%eta
            end if
        end function flow_rate_perzyna

        pure function dflow_dsigma_perzyna(self, sigma_eq, ep) result(res)
            class(Perzyna_overstress), intent(in) :: self
            real(real64), intent(in) :: sigma_eq, ep
            real(real64) :: res
            real(real64) :: y

            y = self%hard_law%stress(ep)
            res = dflow_dpsi(self, sigma_eq/y - 1D0)/y
        end function dflow_dsigma_perzyna

        pure function dflow_dep_perzyna(self, sigma_eq, ep) result(res)
            class(Perzyna_overstress), intent(in) :: self
            real(real64), intent(in) :: sigma_eq, ep
            real(real64) :: res
            real(real64) :: y

            y = self%hard_law%stress(ep)
            res = -dflow_dpsi(self, sigma_eq/y - 1D0)*sigma_eq*self%hard_law%dstress_dep(ep)/y**2
        end function dflow_dep_perzyna

end module muscle_vp_perzyna
