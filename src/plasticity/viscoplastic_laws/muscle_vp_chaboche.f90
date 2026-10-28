! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_vp_chaboche
    !! Additive power overstress: Chaboche (2008), Eq. (15), p. 1646.
    !! g = <(sigma_eq - Y(ep))/K>**n; K is independent of Y.
    !! g and its derivatives are 0 for sigma_eq <= Y by convention; the right-hand
    !! derivative at the threshold is not 0 for n <= 1.
    use, intrinsic :: iso_fortran_env, only: real64
    use muscle_vp_overstress_base, only: Base_overstress_law
    use muscle_hard_base, only: Base_hardening_laws
    implicit none
    private
    public :: Chaboche_overstress

    type, extends(Base_overstress_law) :: Chaboche_overstress
        class(Base_hardening_laws), allocatable :: hard_law
        real(real64) :: k = 1D0  !! Viscous resistance K, stress*time**(1/n)
        real(real64) :: n = 1D0  !! Positive rate exponent n
    contains
        procedure :: flow_rate => flow_rate_chaboche
        procedure :: dflow_dsigma => dflow_dsigma_chaboche
        procedure :: dflow_dep => dflow_dep_chaboche
    end type Chaboche_overstress

contains
    pure function dflow_dz(self, z) result(res)
        ! Slope in z=(sigma_eq-Y)/K; Chaboche (2008), Eq. (15), chain rule.
        ! Called by the two analytical derivative bindings.
        implicit none
        class(Chaboche_overstress), intent(in) :: self
        real(real64), intent(in) :: z
        real(real64) :: res
        if (z <= 0D0) then
            res = 0D0
        else
            res = self%n*z**(self%n - 1D0)
        end if
    end function dflow_dz

    pure function flow_rate_chaboche(self, sigma_eq, ep) result(res)
        !! Flow rate for Closest_point; Chaboche (2008), Eq. (15).
        implicit none
        class(Chaboche_overstress), intent(in) :: self
        real(real64), intent(in) :: sigma_eq, ep
        real(real64) :: res, z  ! z: normalized additive overstress
        z = (sigma_eq - self%hard_law%stress(ep))/self%k
        if (z <= 0D0) then
            res = 0D0
        else
            res = z**self%n
        end if
    end function flow_rate_chaboche

    pure function dflow_dsigma_chaboche(self, sigma_eq, ep) result(res)
        !! Stress sensitivity for Closest_point; Chaboche (2008), Eq. (15), chain rule.
        implicit none
        class(Chaboche_overstress), intent(in) :: self
        real(real64), intent(in) :: sigma_eq, ep
        real(real64) :: res
        res = dflow_dz(self, (sigma_eq - self%hard_law%stress(ep))/self%k)/self%k
    end function dflow_dsigma_chaboche

    pure function dflow_dep_chaboche(self, sigma_eq, ep) result(res)
        !! Hardening sensitivity for Closest_point; Chaboche (2008), Eq. (15), chain rule.
        implicit none
        class(Chaboche_overstress), intent(in) :: self
        real(real64), intent(in) :: sigma_eq, ep
        real(real64) :: res
        res = -dflow_dz(self, (sigma_eq - self%hard_law%stress(ep))/self%k) &
              *self%hard_law%dstress_dep(ep)/self%k
    end function dflow_dep_chaboche
end module muscle_vp_chaboche
