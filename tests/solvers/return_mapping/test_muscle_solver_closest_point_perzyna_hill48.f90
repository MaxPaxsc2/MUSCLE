! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_solver_closest_point_perzyna_hill48
    implicit none
    logical :: passed
    call test_perzyna_hill48_mixed_residual_tangent_static_limit(passed)
    if (.not. passed) STOP 1
    print*, "Passed!", passed
    STOP 0
end program test_muscle_solver_closest_point_perzyna_hill48

subroutine test_perzyna_hill48_mixed_residual_tangent_static_limit(passed)
    ! Called by test_muscle_solver_closest_point_perzyna_hill48, STOP 1.
    ! All six components couple anisotropy, flow, curvature and the static limit in one plastic step.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_bilinear, only : Bilinear_hardening
    use muscle_vp_perzyna, only : Perzyna_overstress
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point_perzyna, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    use muscle_yield_hill48, only : Hill48
    use muscle_solver_closest_point, only : Static_point => Closest_point, STATIC_OK => STATUS_CONVERGED
    implicit none
    logical, intent(out) :: passed
    ! Absolute flow residuals. Measured: tensor norm <= 7.47e-18, multiplier <= 6.25e-17.
    real(real64), parameter :: TOL_RES = 1D-9
    ! Relative tangent error against the library numerical tangent. Measured: <= 2.78e-9.
    real(real64), parameter :: TOL_TAN = 1D-7
    ! Static-limit stress/plastic tensor relative error. Measured: <= 4.21e-7.
    ! Scalar multiplier absolute error. Measured: <= 1.52e-8.
    real(real64), parameter :: TOL_LIMIT = 1D-5
    type(Hill48) :: yld
    type(Elasticity_linear) :: elas
    type(Bilinear_hardening) :: hard
    type(Perzyna_overstress) :: vp
    type(Closest_point) :: solver
    type(Static_point) :: reference
    type(Plastic_material_history) :: history, h_static
    type(ten_3D2Osym) :: strain, flow, plastic_increment, residual ! Delta lambda*n, plastic increment, difference
    type(ten_3D4O2sym) :: tangent, numerical
    real(real64) :: dl, q ! Plastic multiplier increment and equivalent stress
    integer :: status, static_status, iters

    passed = .true.
    yld = Hill48(f=0.2115D0/(0.6923D0*1.2115D0), g=1D0/1.2115D0, h=0.2115D0/1.2115D0, &
                 l=1.25D0, m=1.75D0, n=(0.2115D0+0.6923D0)*4.1538D0/(2D0*0.6923D0*1.2115D0))
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    hard = Bilinear_hardening(y0=50D0, K=100D0)
    vp%hard_law = hard
    vp%eta = 0.02D0
    vp%n = 1D0
    call strain%init(xx=0.09D0, yy=-0.025D0, zz=-0.015D0, xy=0.025D0, yz=-0.015D0, xz=0.01D0)
    call history%init()
    call solver%init(elas, vp, yld, tol=1D-12)
    call solver%solve(strain, history, dt=0.1D0, status=status, iters=iters)
    call check_true(passed, "Mixed step converged", status == STATUS_CONVERGED)
    q = yld%stress_eq(history%state_np1%stress)
    dl = history%state_np1%strain_pf - history%state_n%strain_pf
    flow = dl*yld%dstressEq_dstress(history%state_np1%stress)
    plastic_increment = history%state_np1%strain_p - history%state_n%strain_p
    residual = plastic_increment - flow
    call check_approx(passed, "Flow tensor residual norm", residual%norm(), 0D0, tol=TOL_RES)
    call check_approx(passed, "Multiplier balance", dl, 0.1D0*vp%flow_rate(q,dl), tol=TOL_RES)
    call check_true(passed, "Positive overstress", q > hard%stress(dl))
    call solver%tangent(strain, history, 0.1D0, tangent)
    call solver%tangent_numerical(strain, history, 0.1D0, numerical)
    call check_approx(passed, "Mixed consistent tangent", tangent, numerical, tol=TOL_TAN)
    call reference%init(elas, hard, yld)
    call h_static%init()
    call reference%solve(strain, h_static, static_status)
    vp%eta = 1D-5
    history = Plastic_material_history()
    call solver%init(elas, vp, yld, tol=1D-9)
    call solver%solve(strain, history, 1D0, status)
    call check_true(passed, "Small-viscosity step converged", status == STATUS_CONVERGED)
    call check_true(passed, "Static reference converged", static_status == STATIC_OK)
    call check_approx(passed, "Static-limit stress", history%state_np1%stress, h_static%state_np1%stress, tol=TOL_LIMIT)
    call check_approx(passed, "Static-limit plastic tensor", history%state_np1%strain_p, &
                      h_static%state_np1%strain_p, tol=TOL_LIMIT)
    call check_approx(passed, "Static-limit multiplier", history%state_np1%strain_pf, h_static%state_np1%strain_pf, tol=TOL_LIMIT)
end subroutine test_perzyna_hill48_mixed_residual_tangent_static_limit
