! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_vp_chaboche
    implicit none
    logical :: passed
    call test_chaboche_voce_value(passed)
    if (.not. passed) STOP 1
    call test_chaboche_derivatives(passed)
    if (.not. passed) STOP 2
    call test_chaboche_perzyna_equivalence(passed)
    if (.not. passed) STOP 3
    call test_chaboche_inactive_threshold(passed)
    if (.not. passed) STOP 4
    call test_chaboche_solver_mixed_tangent(passed)
    if (.not. passed) STOP 5
    print *, "Passed!", passed
end program test_muscle_vp_chaboche

subroutine test_chaboche_voce_value(passed)
    ! Nonlinear C11000 hardening distinguishes viscous K from the current yield stress.
    ! Voce parameters: Navarrete et al. (2025), Table 2.
    use, intrinsic :: iso_fortran_env, only: real64
    use muscle_vp_chaboche, only: Chaboche_overstress
    use muscle_hard_vocemod, only: Voce_modified_hardening
    use muscle_test_utils, only: check_approx
    implicit none
    logical, intent(out) :: passed
    type(Chaboche_overstress) :: law
    ! Scalar reference evaluated independently with 60-digit decimal arithmetic.
    real(real64), parameter :: TOL = 1D-13 ! Round-off; Measured: 3.8e-16 relative.
    passed = .true.
    law%hard_law = Voce_modified_hardening(Sy=154.6D0, k=345.8D0, q=122.1D0, n=620.7D0)
    law%k = 80D0
    law%n = 2.3D0
    call check_approx(passed, "Voce flow", law%flow_rate(400D0, .002D0), 4.776681187143114305D0, tol=TOL)
end subroutine test_chaboche_voce_value

subroutine test_chaboche_derivatives(passed)
    ! Fractional exponent and steep hardening; states near (1.05 Y) and far (2 Y)
    ! from the threshold test the chain rule against central differences.
    use, intrinsic :: iso_fortran_env, only: real64
    use muscle_vp_chaboche, only: Chaboche_overstress
    use muscle_hard_vocemod, only: Voce_modified_hardening
    use muscle_math_derivatives, only: derivative
    use muscle_test_utils, only: check_approx
    implicit none
    logical, intent(out) :: passed
    type(Chaboche_overstress) :: law
    real(real64) :: sigma, y
    real(real64), parameter :: ep = .002D0
    real(real64), parameter :: TOL = 1D-8 ! Central differences; Measured: 1.8e-9 relative or less.
    real(real64), parameter :: factors(2) = [1.05D0, 2D0]
    integer :: i
    passed = .true.
    law%hard_law = Voce_modified_hardening(Sy=154.6D0, k=345.8D0, q=122.1D0, n=620.7D0)
    law%k = 80D0
    law%n = 2.3D0
    y = law%hard_law%stress(ep)
    do i = 1, 2
        sigma = factors(i)*y
        call check_approx(passed, "dsigma", law%dflow_dsigma(sigma, ep)/derivative(in_sigma, sigma), 1D0, tol=TOL)
        call check_approx(passed, "dep", law%dflow_dep(sigma, ep)/derivative(in_ep, ep), 1D0, tol=TOL)
    end do
contains
    pure function in_sigma(x) result(res)
        ! Scalar flow wrapper for the library derivative in stress.
        real(real64), intent(in) :: x
        real(real64) :: res
        res = law%flow_rate(x, ep)
    end function in_sigma
    pure function in_ep(x) result(res)
        ! Scalar flow wrapper for the library derivative in accumulated strain.
        real(real64), intent(in) :: x
        real(real64) :: res
        res = law%flow_rate(sigma, x)
    end function in_ep
end subroutine test_chaboche_derivatives

subroutine test_chaboche_perzyna_equivalence(passed)
    ! K = Y(ep) and eta = 1 give the Perzyna flow and stress slope at an active state
    ! with g /= 1. With hardening, the strain slopes differ by Y/sigma_eq because K stays
    ! fixed while Y moves; Perzyna divides by the current Y.
    use, intrinsic :: iso_fortran_env, only: real64
    use muscle_vp_chaboche, only: Chaboche_overstress
    use muscle_vp_perzyna, only: Perzyna_overstress
    use muscle_hard_bilinear, only: Bilinear_hardening
    use muscle_test_utils, only: check_approx
    implicit none
    logical, intent(out) :: passed
    type(Chaboche_overstress) :: law
    type(Perzyna_overstress) :: ref
    real(real64), parameter :: sigma = 200D0, ep = .002D0
    real(real64), parameter :: TOL = 1D-13 ! Algebraic equivalence; Measured: 3.5e-16 relative.
    real(real64) :: y
    passed = .true.
    law%hard_law = Bilinear_hardening(y0=80D0, k=100D0)
    y = law%hard_law%stress(ep)
    law%k = y
    law%n = 2.3D0
    ref%hard_law = law%hard_law
    ref%eta = 1D0
    ref%n = law%n
    call check_approx(passed, "Perzyna flow", law%flow_rate(sigma, ep), ref%flow_rate(sigma, ep), tol=TOL)
    call check_approx(passed, "Perzyna dsigma", law%dflow_dsigma(sigma, ep), ref%dflow_dsigma(sigma, ep), tol=TOL)
    call check_approx(passed, "Perzyna dep", law%dflow_dep(sigma, ep), ref%dflow_dep(sigma, ep)*y/sigma, tol=TOL)
end subroutine test_chaboche_perzyna_equivalence

subroutine test_chaboche_inactive_threshold(passed)
    ! Fractional n<1 exposes accidental powers of negative/zero overstress in derivatives.
    use muscle_vp_chaboche, only: Chaboche_overstress
    use muscle_hard_bilinear, only: Bilinear_hardening
    use muscle_test_utils, only: check_approx
    implicit none
    logical, intent(out) :: passed
    type(Chaboche_overstress) :: law
    passed = .true.
    law%hard_law = Bilinear_hardening(y0=80D0, k=100D0)
    law%k = 30D0
    law%n = .5D0
    ! Inactive scalar identities are exact; Measured: 0.
    call check_approx(passed, "below flow", law%flow_rate(70D0, .2D0), 0D0, tol=0D0)
    call check_approx(passed, "below dsigma", law%dflow_dsigma(70D0, .2D0), 0D0, tol=0D0)
    call check_approx(passed, "below dep", law%dflow_dep(70D0, .2D0), 0D0, tol=0D0)
    call check_approx(passed, "threshold flow", law%flow_rate(100D0, .2D0), 0D0, tol=0D0)
    call check_approx(passed, "threshold dsigma", law%dflow_dsigma(100D0, .2D0), 0D0, tol=0D0)
    call check_approx(passed, "threshold dep", law%dflow_dep(100D0, .2D0), 0D0, tol=0D0)
end subroutine test_chaboche_inactive_threshold

subroutine test_chaboche_solver_mixed_tangent(passed)
    ! All normal/shear components act together; the solver sees only Base_overstress_law.
    use, intrinsic :: iso_fortran_env, only: real64
    use muscle_tensors
    use muscle_vp_chaboche, only: Chaboche_overstress
    use muscle_hard_swift, only: Swift_hardening
    use muscle_yield_vonmises, only: VonMises
    use muscle_elasticity_linear, only: Elasticity_linear
    use muscle_plastic_history, only: Plastic_material_history
    use muscle_solver_closest_point_perzyna, only: Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only: check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    type(Chaboche_overstress) :: law
    type(Closest_point) :: solver
    type(Plastic_material_history) :: history
    type(Elasticity_linear) :: elas
    type(VonMises) :: vm
    type(ten_3D2Osym) :: strain, zero, residual
    type(ten_3D4O2sym) :: tangent, numerical
    real(real64) :: p, seq
    real(real64), parameter :: dt = .1D0
    real(real64), parameter :: TOL = 1D-7 ! Tangent; Measured: 8.5e-10 relative.
    integer :: status, iters
    passed = .true.
    call strain%init(xx=.08D0, yy=-.024D0, zz=-.01D0, xy=.012D0, yz=-.007D0, xz=.004D0)
    zero = 0D0
    call elas%set_parameters(young=1000D0, poisson=.3D0)
    law%hard_law = Swift_hardening(k=100D0, n=.1D0, e0=1D-4)
    law%k = 80D0
    law%n = 2.3D0
    call history%init(strain_p=zero, strain_pf=0D0)
    call solver%init(elasticity=elas, viscoplastic=law, yield=vm)
    call solver%solve(strain, history, dt, status, iters)
    p = history%state_np1%strain_pf
    seq = vm%stress_eq(history%state_np1%stress)
    residual = zero - history%state_np1%strain_p + p*vm%dstressEq_dstress(history%state_np1%stress)
    call check_true(passed, "mixed converged", status == STATUS_CONVERGED)
    call check_true(passed, "mixed plastic", p > 0D0)
    ! Absolute flow balance; Measured: 2.6e-13.
    call check_approx(passed, "mixed flow residual", p, dt*law%flow_rate(seq, p), tol=1D-8)
    ! Absolute strain balance; Measured: 2.1e-18.
    call check_approx(passed, "mixed strain residual", residual%norm(), 0D0, tol=1D-8)
    call solver%tangent(strain, history, dt, tangent)
    call solver%tangent_numerical(strain, history, dt, numerical)
    call check_approx(passed, "mixed tangent", tangent, numerical, tol=TOL)
end subroutine test_chaboche_solver_mixed_tangent
