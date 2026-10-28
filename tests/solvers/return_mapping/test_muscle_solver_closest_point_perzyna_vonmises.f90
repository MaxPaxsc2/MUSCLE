! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_solver_closest_point_perzyna_vonmises
    implicit none
    logical :: passed
    call test_perzyna_vonmises_exact_relaxation_commit(passed)
    if (.not. passed) STOP 1
    call test_perzyna_vonmises_rates_unload_reload_rollback(passed)
    if (.not. passed) STOP 2
    call test_perzyna_vonmises_extreme_exponents_near_yield(passed)
    if (.not. passed) STOP 3
    print*, "Passed!", passed
    STOP 0
end program test_muscle_solver_closest_point_perzyna_vonmises

subroutine test_perzyna_vonmises_exact_relaxation_commit(passed)
    ! Called by test_muscle_solver_closest_point_perzyna_vonmises, STOP 1.
    ! Fixed total strain tests accumulated history against the exact discrete radial solution.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_bilinear, only : Bilinear_hardening
    use muscle_vp_perzyna, only : Perzyna_overstress
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point_perzyna, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    use muscle_yield_vonmises, only : VonMises
    implicit none
    logical, intent(out) :: passed
    type(Elasticity_linear) :: elas
    type(Perzyna_overstress) :: vp
    type(VonMises) :: vm
    type(Closest_point) :: solver
    type(Plastic_material_history) :: history
    type(ten_3D2Osym) :: strain
    real(real64) :: q, expected, previous, dt, factor ! Discrete decay denominator and previous exact q
    integer :: status, i
    ! Relative stress error in the exact implicit-Euler recurrence. Measured: <= 2.58e-16.
    real(real64), parameter :: TOL_Q = 1D-10
    passed = .true.
    call elas%set_parameters(young=210000D0, poisson=0.3D0)
    vp%hard_law = Bilinear_hardening(y0=300D0, K=0D0)
    vp%n = 1D0
    vp%eta = 1D0
    call solver%init(elas, vp, vm, tol=1D-13)
    call history%init()
    call strain%init(xx=0.004D0, yy=-0.001D0, zz=-0.001D0, xy=0.001D0, yz=-0.0007D0, xz=0.0004D0)
    call solver%solve(strain, history, 0.001D0, status)
    call check_true(passed, "Initial relaxation state converged", status == STATUS_CONVERGED)
    previous = vm%stress_eq(history%state_np1%stress)
    do i = 1, 4
        call history%commit()
        dt = 0.0002D0*i
        factor = 1D0 + 3D0*(210000D0/2.6D0)*dt/300D0
        expected = 300D0 + (previous-300D0)/factor
        call solver%solve(strain, history, dt, status)
        call check_true(passed, "Relaxation step converged", status == STATUS_CONVERGED)
        q = vm%stress_eq(history%state_np1%stress)
        call check_approx(passed, "Exact discrete relaxation", q, expected, tol=TOL_Q)
        previous = expected
    end do
end subroutine test_perzyna_vonmises_exact_relaxation_commit

subroutine test_perzyna_vonmises_rates_unload_reload_rollback(passed)
    ! Called by test_muscle_solver_closest_point_perzyna_vonmises, STOP 2.
    ! A rate change and elastic unload separate commit/rollback from candidate state updates.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_bilinear, only : Bilinear_hardening
    use muscle_vp_perzyna, only : Perzyna_overstress
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point_perzyna, only : Closest_point, STATUS_CONVERGED, STATUS_ELASTIC_CASE
    use muscle_test_utils, only : check_approx, check_true
    use muscle_yield_vonmises, only : VonMises
    implicit none
    logical, intent(out) :: passed
    type(Elasticity_linear) :: elas
    type(Perzyna_overstress) :: vp
    type(VonMises) :: vm
    type(Closest_point) :: solver
    type(Plastic_material_history) :: history, slow, saved
    type(ten_3D2Osym) :: strain, direction, trial_strain
    real(real64) :: totals(2), dl, dt ! Sums of Delta lambda for slow and fast loading
    integer :: status, rate, i
    ! Transactional copies. Measured: 0 for all these checks.
    real(real64), parameter :: TOL_HISTORY = 1D-14
    ! Accumulated implicit flow balance in strain units. Measured: <= 8.68e-18.
    real(real64), parameter :: TOL_SUM = 1D-14
    passed = .true.
    call elas%set_parameters(young=210000D0, poisson=0.3D0)
    vp%hard_law = Bilinear_hardening(y0=300D0, K=0D0)
    vp%n = 1D0
    vp%eta = 1D0
    call solver%init(elas, vp, vm, tol=1D-13)
    call direction%init(xx=0.004D0, yy=-0.001D0, zz=-0.001D0, xy=0.001D0, yz=-0.0007D0, xz=0.0004D0)
    do rate = 1, 2
        history = Plastic_material_history()
        dt = 0.01D0
        if (rate == 2) dt = 0.00001D0
        totals(rate) = 0D0
        do i = 1, 4
            saved = history
            strain = (0.5D0*i)*direction
            if (i == 3) strain = history%state_n%strain_p + 0.1D0*direction
            call solver%solve(strain, history, dt, status)
            if (i == 3) then
                call check_true(passed, "Elastic unload", status == STATUS_ELASTIC_CASE)
            else
                call check_true(passed, "Load/reload converged", status == STATUS_CONVERGED)
            end if
            dl = dt*vp%flow_rate(vm%stress_eq(history%state_np1%stress), history%state_np1%strain_pf)
            totals(rate) = totals(rate) + dl
            call check_approx(passed, "Sum of multiplier increments", history%state_np1%strain_pf, totals(rate), &
                              tol=TOL_SUM)
            call check_approx(passed, "Solve preserves committed stress", history%state_n%stress, &
                              saved%state_n%stress, tol=TOL_HISTORY)
            call check_approx(passed, "Solve preserves committed plastic tensor", history%state_n%strain_p, &
                              saved%state_n%strain_p, tol=TOL_HISTORY)
            call check_approx(passed, "Solve preserves committed multiplier", history%state_n%strain_pf, &
                              saved%state_n%strain_pf, tol=TOL_HISTORY)
            call history%commit()
            call check_approx(passed, "Commit stress", history%state_n%stress, history%state_np1%stress, &
                              tol=TOL_HISTORY)
            call check_approx(passed, "Commit plastic tensor", history%state_n%strain_p, &
                              history%state_np1%strain_p, tol=TOL_HISTORY)
            call check_approx(passed, "Commit multiplier", history%state_n%strain_pf, &
                              history%state_np1%strain_pf, tol=TOL_HISTORY)
            saved = history
            trial_strain = strain + direction
            call solver%solve(trial_strain, history, dt, status)
            call history%rollback()
            call check_approx(passed, "Rollback committed stress", history%state_n%stress, &
                              saved%state_n%stress, tol=TOL_HISTORY)
            call check_approx(passed, "Rollback committed plastic tensor", history%state_n%strain_p, &
                              saved%state_n%strain_p, tol=TOL_HISTORY)
            call check_approx(passed, "Rollback committed multiplier", history%state_n%strain_pf, &
                              saved%state_n%strain_pf, tol=TOL_HISTORY)
            call check_approx(passed, "Rollback candidate stress", history%state_np1%stress, &
                              saved%state_n%stress, tol=TOL_HISTORY)
            call check_approx(passed, "Rollback candidate plastic tensor", history%state_np1%strain_p, &
                              saved%state_n%strain_p, tol=TOL_HISTORY)
            call check_approx(passed, "Rollback candidate multiplier", history%state_np1%strain_pf, &
                              saved%state_n%strain_pf, tol=TOL_HISTORY)
        end do
        if (rate == 1) slow = history
    end do
    call check_true(passed, "Faster load retains higher stress", &
                    vm%stress_eq(history%state_n%stress) > vm%stress_eq(slow%state_n%stress))
    call check_true(passed, "Slower load accumulates more plastic strain", totals(1) > totals(2))
end subroutine test_perzyna_vonmises_rates_unload_reload_rollback

subroutine test_perzyna_vonmises_extreme_exponents_near_yield(passed)
    ! Called by test_muscle_solver_closest_point_perzyna_vonmises, STOP 3.
    ! A 1 percent trial overstress probes the power law close to its activation surface, with
    ! N = 0.5 (infinite slope at yield) and N = 10 (flat slope at yield).
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_bilinear, only : Bilinear_hardening
    use muscle_vp_perzyna, only : Perzyna_overstress
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point_perzyna, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    use muscle_yield_vonmises, only : VonMises
    implicit none
    logical, intent(out) :: passed
    type(Elasticity_linear) :: elas
    type(Perzyna_overstress) :: vp
    type(VonMises) :: vm
    type(Closest_point) :: solver
    type(Plastic_material_history) :: history
    type(ten_3D2Osym) :: strain
    type(ten_3D4O2sym) :: tangent, numerical
    real(real64) :: q, psi, dl
    integer :: status, i
    real(real64), parameter :: EXPONENTS(2) = [0.5D0, 10D0], ETAS(2) = [10D0, 1D-20], DTS(2) = [1D-4, 1D-3]
    ! Numerical tangent truncation at the default step. Measured: <= 4.71e-7 (N = 0.5), <= 1.18e-5 (N = 10).
    real(real64), parameter :: TOL_TAN(2) = [1D-6, 2D-5]
    ! Absolute multiplier residual in strain units. Measured: <= 5.22e-16.
    real(real64), parameter :: TOL_FLOW = 1D-11
    passed = .true.
    call elas%set_parameters(young=210000D0, poisson=0.3D0)
    vp%hard_law = Bilinear_hardening(y0=300D0, K=0D0)
    call strain%init(xx=0.004D0, yy=-0.001D0, zz=-0.001D0, xy=0.001D0, yz=-0.0007D0, xz=0.0004D0)
    strain = (303D0/vm%stress_eq(elas%stress(strain)))*strain
    do i = 1, 2
        vp%n = EXPONENTS(i)
        vp%eta = ETAS(i)
        call solver%init(elas, vp, vm, tol=1D-13)
        history = Plastic_material_history()
        call solver%solve(strain, history, DTS(i), status)
        call check_true(passed, "Extreme exponent converged", status == STATUS_CONVERGED)
        q = vm%stress_eq(history%state_np1%stress)
        psi = q/300D0 - 1D0
        call check_true(passed, "Small positive final overstress", psi > 0D0 .and. psi < 0.01D0)
        dl = history%state_np1%strain_pf
        call check_approx(passed, "Extreme exponent flow balance", dl, DTS(i)*vp%flow_rate(q, dl), tol=TOL_FLOW)
        call solver%tangent(strain, history, DTS(i), tangent)
        call solver%tangent_numerical(strain, history, DTS(i), numerical)
        call check_approx(passed, "Extreme exponent tangent", tangent, numerical, tol=TOL_TAN(i))
    end do
end subroutine test_perzyna_vonmises_extreme_exponents_near_yield
