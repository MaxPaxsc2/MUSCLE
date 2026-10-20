! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_vp_peric
    implicit none
    logical :: passed

    call test_peric_zero_rate_flow_and_slope(passed)
    if (.not. passed) STOP 1

    call test_peric_flow_stress_values(passed)
    if (.not. passed) STOP 2

    call test_peric_dstress_dep_analytic_and_numeric(passed)
    if (.not. passed) STOP 3

    print*, "Passed!", passed
    STOP 0
end program test_muscle_vp_peric


subroutine test_peric_zero_rate_flow_and_slope(passed)
    !! At zero rate the flow stress is the static hardening and the slope stays finite.
    ! Called by program test_muscle_vp_peric (STOP 1)
    ! Zero rate is where a power-law overstress has an infinite slope; here the slope is
    ! sigma_s' + sigma_s*eps*mu/dt. References: mpmath (50 digits), de Souza Neto et al. (2008, Eq. 11.94).
    ! Copper C11000 hardening: Navarrete et al. (2025), Table 2.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_vocemod, only : Voce_modified_hardening
    use muscle_vp_peric, only : Peric_viscoplastic
    use muscle_test_utils, only : check_approx
    implicit none
    logical, intent(out) :: passed
    ! Relative error of closed-form values (round-off only). Measured: 1.3e-16 or less.
    real(real64), parameter :: TOL = 1D-15
    type(Voce_modified_hardening) :: voce
    type(Peric_viscoplastic) :: relax, moderate

    passed = .true.
    voce = Voce_modified_hardening(Sy=154.6D0, k=345.8D0, q=122.1D0, n=620.7D0)
    relax = Peric_viscoplastic(mu=4D7, eps_rate=0.02D0)
    relax%hard_law = voce
    moderate = Peric_viscoplastic(mu=1D0, eps_rate=0.2D0)
    moderate%hard_law = voce
    call check_approx(passed, "Flow stress at zero rate equals Voce", relax%flow_stress(0.015D0, 0D0), &
                      2.8187595369856579D+2, tol=TOL)
    call check_approx(passed, "Slope at zero rate, mu = 4e7 s, eps = 0.02", relax%dstress_dep(0.015D0, 0D0, 0.5D0), &
                      4.5100187857414457D+8, tol=TOL)
    call check_approx(passed, "Slope at zero rate, mu = 1 s, eps = 0.2", moderate%dstress_dep(0.015D0, 0D0, 0.5D0), &
                      4.654068207796387D+2, tol=TOL)
end subroutine test_peric_zero_rate_flow_and_slope


subroutine test_peric_flow_stress_values(passed)
    !! Flow stress from relaxation rates (1e-9 1/s) to fast tests (1e3 1/s), two parameter sets.
    ! Called by program test_muscle_vp_peric (STOP 2)
    ! Fixes the formula and its composition with the hardening law; mpmath (50 digits).
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_vocemod, only : Voce_modified_hardening
    use muscle_vp_peric, only : Peric_viscoplastic
    use muscle_test_utils, only : check_approx
    implicit none
    logical, intent(out) :: passed
    ! Relative error of closed-form values (round-off of one real power). Measured: 1.9e-16 or less.
    real(real64), parameter :: TOL = 1D-15
    type(Voce_modified_hardening) :: voce
    type(Peric_viscoplastic) :: relax, moderate

    passed = .true.
    voce = Voce_modified_hardening(Sy=154.6D0, k=345.8D0, q=122.1D0, n=620.7D0)
    relax = Peric_viscoplastic(mu=4D7, eps_rate=0.02D0)
    relax%hard_law = voce
    moderate = Peric_viscoplastic(mu=1D0, eps_rate=0.2D0)
    moderate%hard_law = voce
    call check_approx(passed, "Flow stress, mu = 4e7 s, eps = 0.02", &
                      [relax%flow_stress(0D0, 1D-9), relax%flow_stress(0.015D0, 1D-4), relax%flow_stress(0.05D0, 1D3)], &
                      [1.5472131802064219D+2, 3.3273723214919666D+2, 4.7904218334723294D+2], tol=TOL)
    call check_approx(passed, "Flow stress, mu = 1 s, eps = 0.2", &
                      [moderate%flow_stress(0D0, 1D-9), moderate%flow_stress(0.015D0, 1D-4), &
                       moderate%flow_stress(0.05D0, 1D3)], &
                      [1.5460000003092D+2, 2.8188159099215253D+2, 1.1706292561888705D+3], tol=TOL)
end subroutine test_peric_flow_stress_values


subroutine test_peric_dstress_dep_analytic_and_numeric(passed)
    !! Slope against mpmath and against centered differences of the flow stress.
    ! Called by program test_muscle_vp_peric (STOP 3)
    ! With mu/dt large the rate term dominates, so an exponent error (eps instead of eps - 1) shows here.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_vocemod, only : Voce_modified_hardening
    use muscle_vp_peric, only : Peric_viscoplastic
    use muscle_math_derivatives, only : derivative
    use muscle_test_utils, only : check_approx
    implicit none
    logical, intent(out) :: passed
    ! Relative error of the closed form. Measured: 1.2e-16 or less.
    real(real64), parameter :: TOL = 1D-15
    ! Relative error of centered differences. Steps: 1e-6 in ep (the Voce exponent n*h stays below 1e-3)
    ! and 1e-5*(epd + 1/mu) in the rate (the default relative step, 1e-7*epd, loses the rate slope to
    ! round-off: 1.4e-6 measured). Measured: 2.6e-10 or less.
    real(real64), parameter :: TOL_NUM = 1D-9
    type(Voce_modified_hardening) :: voce
    type(Peric_viscoplastic) :: law
    real(real64) :: ep, epd, dt, numeric

    passed = .true.
    voce = Voce_modified_hardening(Sy=154.6D0, k=345.8D0, q=122.1D0, n=620.7D0)

    law = Peric_viscoplastic(mu=1D0, eps_rate=0.2D0)
    law%hard_law = voce
    ep = 0.02D0; epd = 1D-3; dt = 0.5D0
    numeric = derivative(flow_of_ep, ep, 1D-6) + derivative(flow_of_epd, epd, 1D-5*(epd + 1D0/law%mu))/dt
    call check_approx(passed, "Slope, mu = 1 s, eps = 0.2", law%dstress_dep(ep, epd, dt), &
                      4.5953251811745998D+2, tol=TOL)
    call check_approx(passed, "Slope against differences, mu = 1 s", law%dstress_dep(ep, epd, dt), numeric, tol=TOL_NUM)

    law = Peric_viscoplastic(mu=4D7, eps_rate=0.02D0)
    law%hard_law = voce
    ep = 0.02D0; epd = 1D-6; dt = 10D0
    numeric = derivative(flow_of_ep, ep, 1D-6) + derivative(flow_of_epd, epd, 1D-5*(epd + 1D0/law%mu))/dt
    call check_approx(passed, "Slope, mu = 4e7 s, eps = 0.02", law%dstress_dep(ep, epd, dt), &
                      5.9643525907136781D+5, tol=TOL)
    call check_approx(passed, "Slope against differences, mu = 4e7 s", law%dstress_dep(ep, epd, dt), numeric, tol=TOL_NUM)

contains

    pure function flow_of_ep(x) result(res)
        ! Flow stress as a function of the plastic strain at the current rate; called by derivative.
        real(real64), intent(in) :: x
        real(real64) :: res
        res = law%flow_stress(x, epd)
    end function flow_of_ep

    pure function flow_of_epd(x) result(res)
        ! Flow stress as a function of the rate at the current plastic strain; called by derivative.
        real(real64), intent(in) :: x
        real(real64) :: res
        res = law%flow_stress(ep, x)
    end function flow_of_epd

end subroutine test_peric_dstress_dep_analytic_and_numeric
