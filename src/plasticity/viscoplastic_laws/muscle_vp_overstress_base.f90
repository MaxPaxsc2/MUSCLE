! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_vp_overstress_base
    !! Module muscle_vp_overstress_base
    !! =================================
    !!
    !! Defines the abstract base type for **overstress** viscoplastic laws (Perzyna,
    !! Chaboche, Bodner-Partom, ...).
    !!
    !! In these models the static yield surface `f = sigma_eq - Y(ep)` can be exceeded and the
    !! plastic multiplier rate is an explicit constitutive law
    !!
    !!     lambda_dot = g(sigma_eq, ep)     (g = 0 if f <= 0 for models with static yield)
    !!
    !! instead of a Lagrange multiplier fixed by the consistency condition `f = 0`
    !! (see `Base_viscoplastic_law` in `muscle_vp_base` for that formulation).
    !!
    !! A concrete law only has to provide `flow_rate`. The derivatives `dflow_dsigma` and
    !! `dflow_dep` are computed numerically by default and should be overridden when an
    !! analytical expression is available.
    !!
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private
    public :: Base_overstress_law

    type, abstract :: Base_overstress_law
    contains
        procedure(flow_rate_interface), deferred :: flow_rate
        procedure :: dflow_dsigma => dflow_dsigma_numeric
        procedure :: dflow_dep => dflow_dep_numeric
    end type Base_overstress_law

    abstract interface
        pure function flow_rate_interface(self, sigma_eq, ep) result(res)
        !! Plastic multiplier rate `g(sigma_eq, ep)` (equivalent plastic strain rate).
            use, intrinsic :: iso_fortran_env, only: real64
            import Base_overstress_law
            class(Base_overstress_law), intent(in) :: self
            real(real64), intent(in) :: sigma_eq  !! Equivalent stress
            real(real64), intent(in) :: ep        !! Accumulated equivalent plastic strain
            real(real64) :: res
        end function flow_rate_interface
    end interface

    contains

        pure function dflow_dsigma_numeric(self, sigma_eq, ep) result(res)
        !! Numerical derivative of `flow_rate` with respect to the equivalent stress.
            use muscle_math_derivatives
            class(Base_overstress_law), intent(in) :: self
            real(real64), intent(in) :: sigma_eq, ep
            real(real64) :: res

            res = derivative(wrapper, sigma_eq)

            contains
            pure function wrapper(x) result(res_w)
                real(real64), intent(in) :: x
                real(real64) :: res_w
                res_w = self%flow_rate(x, ep)
            end function wrapper
        end function dflow_dsigma_numeric

        pure function dflow_dep_numeric(self, sigma_eq, ep) result(res)
        !! Numerical derivative of `flow_rate` with respect to the equivalent plastic strain.
            use muscle_math_derivatives
            class(Base_overstress_law), intent(in) :: self
            real(real64), intent(in) :: sigma_eq, ep
            real(real64) :: res

            res = derivative(wrapper, ep)

            contains
            pure function wrapper(x) result(res_w)
                real(real64), intent(in) :: x
                real(real64) :: res_w
                res_w = self%flow_rate(sigma_eq, x)
            end function wrapper
        end function dflow_dep_numeric

end module muscle_vp_overstress_base
