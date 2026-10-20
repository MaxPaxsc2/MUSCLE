! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_solver_closest_point_time_peric
    implicit none
    logical :: passed

    call test_closest_point_time_peric_vonmises_uniaxial_tensile(passed)
    if (.not. passed) STOP 1

    call test_closest_point_time_peric_cpb06_three_shears(passed)
    if (.not. passed) STOP 2

    call test_closest_point_time_peric_long_step_rate_independent(passed)
    if (.not. passed) STOP 3

    print*, "Passed!", passed
    STOP 0
end program test_muscle_solver_closest_point_time_peric


subroutine test_closest_point_time_peric_vonmises_uniaxial_tensile(passed)
    !! Uniaxial tension from the virgin state: the first Newton iteration starts at zero rate.
    ! Called by program test_muscle_solver_closest_point_time_peric (STOP 1)
    ! With a power-law overstress this step stalls at zero rate; the Peric slope there is finite.
    ! Exact uniaxial state for ep = 0.01, dt = 10 s: sigma = sigma_s(ep)*(1 + mu*ep/dt)**eps
    ! de Souza Neto et al. (2008, Eq. 11.94), mpmath (50 digits). Voce: Navarrete et al. (2025), Table 2.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_vocemod, only : Voce_modified_hardening
    use muscle_vp_peric, only : Peric_viscoplastic
    use muscle_yield_vonmises, only : VonMises
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_solver_closest_point_time, only : Closest_point, Closest_point_data, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Relative stress and plastic strain error; the last Newton step lands at round-off. Measured: 4.6e-15 or less.
    real(real64), parameter :: TOL_STRESS = 1D-12
    ! Equivalent plastic strain error (scalar scale = 1). Measured: 1.9e-17.
    real(real64), parameter :: TOL_PF = 1D-14
    ! Relative tangent error against centered differences. Measured: 6.3e-10.
    real(real64), parameter :: TOL_TAN = 1D-8
    type(Voce_modified_hardening) :: voce
    type(Peric_viscoplastic) :: peric
    type(VonMises) :: vm
    type(Elasticity_linear) :: elas
    type(Closest_point) :: solver
    type(Closest_point_data) :: data
    type(ten_3D2Osym) :: strain, strain_p, stress, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: strain_pf
    integer :: status, iters

    passed = .true.
    voce = Voce_modified_hardening(Sy=154.6D0, k=345.8D0, q=122.1D0, n=620.7D0)
    peric = Peric_viscoplastic(mu=4D7, eps_rate=0.02D0)
    peric%hard_law = voce
    call elas%set_parameters(young=1D5, poisson=0.34D0)
    call expected_stress%init(xx=3.4598943176413219D+2, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call expected_strain_p%init(xx=0.01D0, yy=-5D-3, zz=-5D-3, xy=0D0, yz=0D0, xz=0D0)
    call strain%init(xx=1.3459894317641322D-2, yy=-6.1763640679980494D-3, zz=-6.1763640679980494D-3, &
                     xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call data%init(strain_pf=0D0, strain_p=strain_p, dt=10D0)
    call solver%init(elasticity=elas, hardening=peric, yield=vm)
    call solver%solve(strain=strain, data=data)
    call data%get(stress=stress, strain_pf=strain_pf, strain_p=strain_p, status=status, iters=iters)
    call check_true(passed, "Converged status", status == STATUS_CONVERGED)
    call check_approx(passed, "Stress", stress, expected_stress, tol=TOL_STRESS)
    call check_approx(passed, "Equivalent plastic strain", strain_pf, 0.01D0, tol=TOL_PF)
    call check_approx(passed, "Plastic strain", strain_p, expected_strain_p, tol=TOL_STRESS)
    ! Newton iterations from zero rate. Measured: 6 (Release and Debug).
    call check_true(passed, "Iterations <= 6", iters <= 6)
    call solver%tangent_numerical(strain=strain, data=data, tangent=numerical_tangent)
    call solver%tangent(strain=strain, data=data, tangent=tangent)
    call check_approx(passed, "Consistent tangent", tangent, numerical_tangent, tol=TOL_TAN)
end subroutine test_closest_point_time_peric_vonmises_uniaxial_tensile


subroutine test_closest_point_time_peric_cpb06_three_shears(passed)
    !! CPB06 of rolled copper with xy, yz and xz strains loaded at once: the law must not depend on the criterion.
    ! Called by program test_muscle_solver_closest_point_time_peric (STOP 2)
    ! Reference: mpmath solution (50 digits) of the discrete system eps_p = dlambda*grad(f), f = sigma_flow,
    ! with CPB06 (Cazacu et al., 2006, Eqs. 8, 9, 12) and the C11000 parameters of Navarrete et al. (2025), Table 2.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_vocemod, only : Voce_modified_hardening
    use muscle_vp_peric, only : Peric_viscoplastic
    use muscle_yield_cpb06, only : CPB06
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_solver_closest_point_time, only : Closest_point, Closest_point_data, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Relative stress and plastic strain error at solver tolerance 1e-5. Measured: 1.8e-8 or less.
    real(real64), parameter :: TOL_STRESS = 1D-7
    ! Equivalent plastic strain error (scalar scale = 1). Measured: 1.4e-10.
    real(real64), parameter :: TOL_PF = 1D-9
    ! Relative tangent error against centered differences. Measured: 6.8e-9.
    real(real64), parameter :: TOL_TAN = 1D-7
    type(Voce_modified_hardening) :: voce
    type(Peric_viscoplastic) :: peric
    type(CPB06) :: cpb
    type(Elasticity_linear) :: elas
    type(Closest_point) :: solver
    type(Closest_point_data) :: data
    type(ten_3D2Osym) :: strain, strain_p, stress, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: strain_pf
    integer :: status

    passed = .true.
    voce = Voce_modified_hardening(Sy=154.6D0, k=345.8D0, q=122.1D0, n=620.7D0)
    peric = Peric_viscoplastic(mu=4D7, eps_rate=0.02D0)
    peric%hard_law = voce
    call cpb%init(c11=1D0, c12=0.091D0, c13=0.002D0, &
                  c21=0.091D0, c22=1.07D0, c23=0.028D0, &
                  c31=0.002D0, c32=0.028D0, c33=1.04D0, &
                  c44=0.991D0, c55=0.991D0, c66=0.991D0, k=0D0, a=4D0)
    call elas%set_parameters(young=1D5, poisson=0.34D0)
    call expected_stress%init(xx=5.3216545964977764D+2, yy=1.9758493723796365D+2, zz=2.0774960311225871D+2, &
                              xy=5.7062556347203245D+1, yz=-3.5915016614278371D+1, xz=4.7405633798477728D+1)
    call expected_strain_p%init(xx=8.0564828406929796D-3, yy=-3.4601381589887129D-3, zz=-4.5963446817042667D-3, &
                                xy=2.2353617449474765D-3, yz=-1.5187387773686698D-3, xz=1.8647645071003984D-3)
    call strain%init(xx=0.012D0, yy=-0.004D0, zz=-0.005D0, xy=0.003D0, yz=-0.002D0, xz=0.0025D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call data%init(strain_pf=0D0, strain_p=strain_p, dt=1D0)
    call solver%init(elasticity=elas, hardening=peric, yield=cpb)
    call solver%solve(strain=strain, data=data)
    call data%get(stress=stress, strain_pf=strain_pf, strain_p=strain_p, status=status)
    call check_true(passed, "Converged status", status == STATUS_CONVERGED)
    call check_approx(passed, "Stress", stress, expected_stress, tol=TOL_STRESS)
    call check_approx(passed, "Equivalent plastic strain", strain_pf, 8.8468275180343146D-3, tol=TOL_PF)
    call check_approx(passed, "Plastic strain", strain_p, expected_strain_p, tol=TOL_STRESS)
    call solver%tangent_numerical(strain=strain, data=data, tangent=numerical_tangent)
    call solver%tangent(strain=strain, data=data, tangent=tangent)
    call check_approx(passed, "Consistent tangent", tangent, numerical_tangent, tol=TOL_TAN)
end subroutine test_closest_point_time_peric_cpb06_three_shears


subroutine test_closest_point_time_peric_long_step_rate_independent(passed)
    !! With dt = 1e8*mu the rate factor is 1 + 1e-12 and the step equals the rate-independent return.
    ! Called by program test_muscle_solver_closest_point_time_peric (STOP 3)
    ! Limit of de Souza Neto et al. (2008), Remark 11.4; reference: Closest_point without rate, same Voce.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_vocemod, only : Voce_modified_hardening
    use muscle_vp_peric, only : Peric_viscoplastic
    use muscle_yield_vonmises, only : VonMises
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point_rate_free => Closest_point, &
                                           STATUS_CONVERGED_RATE_FREE => STATUS_CONVERGED
    use muscle_solver_closest_point_time, only : Closest_point, Closest_point_data, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Relative difference between the two solvers; the rate factor adds 2e-12. Measured: 1.3e-12 or less.
    real(real64), parameter :: TOL_STRESS = 1D-10
    ! Equivalent plastic strain difference (scalar scale = 1). Measured: 5.3e-15.
    real(real64), parameter :: TOL_PF = 1D-13
    type(Voce_modified_hardening) :: voce
    type(Peric_viscoplastic) :: peric
    type(VonMises) :: vm
    type(Elasticity_linear) :: elas
    type(Closest_point) :: solver
    type(Closest_point_data) :: data
    type(Closest_point_rate_free) :: solver_rate_free
    type(Plastic_material_history) :: history
    type(ten_3D2Osym) :: strain, zero, strain_p, stress
    real(real64) :: strain_pf
    integer :: status, status_rate_free, iters

    passed = .true.
    voce = Voce_modified_hardening(Sy=154.6D0, k=345.8D0, q=122.1D0, n=620.7D0)
    peric = Peric_viscoplastic(mu=4D7, eps_rate=0.02D0)
    peric%hard_law = voce
    call elas%set_parameters(young=1D5, poisson=0.34D0)
    call strain%init(xx=1.3459894317641322D-2, yy=-6.1763640679980494D-3, zz=-6.1763640679980494D-3, &
                     xy=0D0, yz=0D0, xz=0D0)
    call zero%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call data%init(strain_pf=0D0, strain_p=zero, dt=4D15)
    call solver%init(elasticity=elas, hardening=peric, yield=vm)
    call solver%solve(strain=strain, data=data)
    call data%get(stress=stress, strain_pf=strain_pf, strain_p=strain_p, status=status)
    call history%init(strain_p=zero, strain_pf=0D0)
    call solver_rate_free%init(elasticity=elas, hardening=voce, yield=vm)
    call solver_rate_free%solve(strain=strain, history=history, status=status_rate_free, iters=iters)
    call check_true(passed, "Both converged", &
                    status == STATUS_CONVERGED .and. status_rate_free == STATUS_CONVERGED_RATE_FREE)
    call check_approx(passed, "Stress against rate-independent", stress, history%state_np1%stress, tol=TOL_STRESS)
    call check_approx(passed, "Equivalent plastic strain against rate-independent", strain_pf, &
                      history%state_np1%strain_pf, tol=TOL_PF)
    call check_approx(passed, "Plastic strain against rate-independent", strain_p, &
                      history%state_np1%strain_p, tol=TOL_STRESS)
end subroutine test_closest_point_time_peric_long_step_rate_independent
