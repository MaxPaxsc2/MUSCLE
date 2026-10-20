! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_vp_peric
    !! Module muscle_vp_peric
    !! ======================
    !! Implements the Peric overstress viscoplastic law, de Souza Neto et al. (2008, Eq. 11.86):
    !! \[ \dot{\bar\epsilon}_p = \left[ (\bar\sigma/\sigma_s(\bar\epsilon_p))^{1/\epsilon} - 1 \right] / \mu \]
    !! written in the inverted form solved by `closest_point_time` (Eq. 11.94):
    !! \[ \sigma_{flow} = \sigma_s(\bar\epsilon_p)\,(1 + \mu\,\dot{\bar\epsilon}_p)^{\epsilon} \]
    !! where \( \sigma_s \) is the rate-independent hardening law `hard_law`. At zero rate the flow
    !! stress is \( \sigma_s \), and with \( \mu \to 0 \) or \( \epsilon \to 0 \) the rate-independent
    !! return is recovered (Remark 11.4). Called by `closest_point_time` through `flow_stress` and
    !! `dstress_dep`.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_vp_base, only : Base_viscoplastic_law
    use muscle_hard_base, only : Base_hardening_laws
    implicit none
    private

    public :: Peric_viscoplastic

    type, extends(Base_viscoplastic_law) :: Peric_viscoplastic
        !! Peric law over a rate-independent hardening law.
        class(Base_hardening_laws), allocatable :: hard_law
            !! Static hardening \( \sigma_s(\bar\epsilon_p) \).
        real(real64) :: mu
            !! Viscosity parameter \( \mu \) [time], \( \mu \ge 0 \).
        real(real64) :: eps_rate
            !! Rate sensitivity \( \epsilon > 0 \) (named `eps_rate` to keep the intrinsic `epsilon`).
    contains
        procedure :: flow_stress => flow_peric
        procedure :: dstress_dep => dstress_dep_peric
    end type Peric_viscoplastic

contains

    pure function flow_peric(self, ep, epd) result(res)
        !! Flow stress \( \sigma_s(\bar\epsilon_p)(1 + \mu\dot{\bar\epsilon}_p)^{\epsilon} \) (Eq. 11.94).
        class(Peric_viscoplastic), intent(in) :: self
        real(real64), intent(in) :: ep
            !! Equivalent plastic strain.
        real(real64), intent(in) :: epd
            !! Equivalent plastic strain rate.
        real(real64) :: res

        res = self%hard_law%stress(ep)*(1D0 + self%mu*epd)**self%eps_rate
    end function flow_peric

    pure function dstress_dep_peric(self, ep, epd, dt) result(res)
        !! Derivative of the flow stress with respect to the plastic strain increment:
        !! \( \partial_{\bar\epsilon_p}\sigma_{flow} + \partial_{\dot{\bar\epsilon}_p}\sigma_{flow} / \Delta t \).
        class(Peric_viscoplastic), intent(in) :: self
        real(real64), intent(in) :: ep
            !! Equivalent plastic strain.
        real(real64), intent(in) :: epd
            !! Equivalent plastic strain rate.
        real(real64), intent(in) :: dt
            !! Time increment.
        real(real64) :: res
        real(real64) :: base, factor_slope
            ! base = 1 + mu*epd; factor_slope = base**(eps - 1), in (0, 1] for epd >= 0 and eps <= 1

        base = 1D0 + self%mu*epd
        factor_slope = base**(self%eps_rate - 1D0)
        res = self%hard_law%dstress_dep(ep)*factor_slope*base &
            + self%hard_law%stress(ep)*self%eps_rate*self%mu*factor_slope/dt
    end function dstress_dep_peric

end module muscle_vp_peric
