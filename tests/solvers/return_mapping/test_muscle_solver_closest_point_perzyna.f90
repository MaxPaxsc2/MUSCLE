program test_muscle_solver_closest_point_perzyna
    implicit none
    logical :: passed

    call test_perzyna_uniaxial_tangent(passed)
    if (.not. passed) STOP 1

    call test_perzyna_rate_independent_limit(passed)
    if (.not. passed) STOP 2

    call test_perzyna_elastic(passed)
    if (.not. passed) STOP 3

    call test_perzyna_non_converged(passed)
    if (.not. passed) STOP 4

    call test_perzyna_reinit_other_types(passed)
    if (.not. passed) STOP 5

    print*, "Passed!", passed
end program test_muscle_solver_closest_point_perzyna

subroutine test_perzyna_uniaxial_tangent(passed)
    !! Overstress must remain after the step (f > 0), residual must vanish and the
    !! analytical tangent must match the numerical one.
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_vp_perzyna, only : Perzyna_overstress
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point_perzyna
    implicit none

    logical, intent(out) :: passed
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(VonMises) :: vm
    type(ten_3D2Osym) :: strain, strain_p, stress
    type(Perzyna_overstress) :: pz
    type(Swift_hardening) :: sw
    type(Elasticity_linear) :: elas
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: strain_pf, dt, y, f
    integer :: status, iters

    passed = .False.
    call strain%init(xx=0.05D0, yy=-0.015D0, zz=-0.015D0, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)

    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    pz%hard_law = sw  ! Error in gfortran-12 if not assigned here
    pz%eta = 1D0
    pz%n = 2D0
    dt = 1D-1

    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, viscoplastic=pz, yield=vm)
    call solver%solve(strain=strain, history=history, dt=dt, status=status, iters=iters)
    stress    = history%state_np1%stress
    strain_pf = history%state_np1%strain_pf

    passed = (status .eq. STATUS_CONVERGED)
    if (.not. passed) print*, "Case 1 Not converged. status:", status, " iters:", iters
    if (.not. passed) return

    ! Flow law: dlambda = dt*g(sigma_eq, p)
    y = sw%stress(strain_pf)
    f = vm%stress_eq(stress) - y
    passed = f .gt. 0D0
    if (.not. passed) print*, "Case 2 Overstress must be positive. f:", f
    if (.not. passed) return

    passed = abs(strain_pf - dt*pz%flow_rate(vm%stress_eq(stress), strain_pf)) < 1D-8
    if (.not. passed) print*, "Case 3 Flow law not satisfied", strain_pf, &
                              dt*pz%flow_rate(vm%stress_eq(stress), strain_pf)
    if (.not. passed) return

    call solver%tangent_numerical(strain=strain, history=history, dt=dt, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, dt=dt, tangent=tangent)

    passed = tangent%is_approx(numerical_tangent, tol=1.0D-7)
    if (.not. passed) print*, "Case 4 Tangent is no equal", new_line('A'),    &
                              "Analitical:", tangent, new_line('A'),         &
                              "Numerical:", numerical_tangent, new_line('A'), &
                              "Difference", tangent - numerical_tangent
end subroutine

subroutine test_perzyna_rate_independent_limit(passed)
    !! With eta -> 0 the solution must tend to the rate independent return mapping (f -> 0).
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_vp_perzyna, only : Perzyna_overstress
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point_perzyna
    implicit none

    logical, intent(out) :: passed
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(VonMises) :: vm
    type(ten_3D2Osym) :: strain, strain_p, stress
    type(Perzyna_overstress) :: pz
    type(Swift_hardening) :: sw
    type(Elasticity_linear) :: elas
    real(real64) :: strain_pf, f
    integer :: status

    passed = .False.
    call strain%init(xx=0.05D0, yy=-0.015D0, zz=-0.015D0, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)

    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    pz%hard_law = sw
    pz%eta = 1D-6
    pz%n = 1D0

    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, viscoplastic=pz, yield=vm)
    call solver%solve(strain=strain, history=history, dt=1D0, status=status)
    stress    = history%state_np1%stress
    strain_pf = history%state_np1%strain_pf

    passed = (status .eq. STATUS_CONVERGED)
    if (.not. passed) print*, "Case 5 Not converged. status:", status
    if (.not. passed) return

    f = vm%stress_eq(stress) - sw%stress(strain_pf)
    passed = abs(f)/sw%stress(strain_pf) < 1D-5
    if (.not. passed) print*, "Case 6 Rate independent limit not recovered. f:", f
end subroutine

subroutine test_perzyna_elastic(passed)
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_vp_perzyna, only : Perzyna_overstress
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point_perzyna
    implicit none

    logical, intent(out) :: passed
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(VonMises) :: vm
    type(ten_3D2Osym) :: strain, strain_p, expected_stress, stress
    type(Perzyna_overstress) :: pz
    type(Swift_hardening) :: sw
    type(Elasticity_linear) :: elas
    type(ten_3D4O2sym) :: tangent, expected_tangent
    integer :: status

    passed = .False.
    call strain%init(xx=1D-4, yy=-3D-5, zz=-3D-5, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    expected_stress = elas%stress(strain)

    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    pz%hard_law = sw
    pz%eta = 1D0
    pz%n = 2D0

    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, viscoplastic=pz, yield=vm)
    call solver%solve(strain=strain, history=history, dt=1D-1, status=status)
    stress = history%state_np1%stress

    passed = (status .eq. STATUS_ELASTIC_CASE) .and. stress%is_approx(expected_stress, tol=1.0D-8)
    if (.not. passed) print*, "Case 7 Elastic case failed. status:", status
    if (.not. passed) return

    call solver%tangent(strain=strain, history=history, dt=1D-1, tangent=tangent)
    expected_tangent = elas%dstress_dstrain(strain)
    passed = tangent%is_approx(expected_tangent, tol=1.0D-8)
    if (.not. passed) print*, "Case 8 Elastic tangent failed"
end subroutine

subroutine test_perzyna_non_converged(passed)
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_vp_perzyna, only : Perzyna_overstress
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point_perzyna
    implicit none

    logical, intent(out) :: passed
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(VonMises) :: vm
    type(ten_3D2Osym) :: strain, strain_p
    type(Perzyna_overstress) :: pz
    type(Swift_hardening) :: sw
    type(Elasticity_linear) :: elas
    integer :: status

    passed = .False.
    call strain%init(xx=0.05D0, yy=-0.015D0, zz=-0.015D0, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    pz%hard_law = sw
    pz%eta = 1D0
    pz%n = 2D0

    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, viscoplastic=pz, yield=vm, iter_nw=1)
    call solver%solve(strain=strain, history=history, dt=1D-1, status=status)

    passed = (status .eq. STATUS_NONCONVERGED)
    if (.not. passed) print*, "Case 9 Converged when non converged must occur. status:", status
end subroutine

subroutine test_perzyna_reinit_other_types(passed)
    ! A second init with components of larger dynamic types (Elasticity_frozen -> Elasticity_linear,
    ! VonMises -> Hill48, and a larger overstress law) must give the same step as a solver initialized
    ! with the second set only. Without reallocating each component, the new data overflow the old
    ! storage. Perzyna_overstress is the only concrete overstress law, so the larger law is a local
    ! extension of it with extra storage and the same behaviour.
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_vp_perzyna, only : Perzyna_overstress
    use muscle_hard_bilinear, only : Bilinear_hardening
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_yield_hill48, only : Hill48
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point_perzyna
    use muscle_test_utils, only : check_approx, check_true
    implicit none

    type, extends(Perzyna_overstress) :: Larger_overstress
        real(real64) :: extra(32) = 0D0
    end type Larger_overstress

    logical, intent(out) :: passed
    type(Plastic_material_history) :: history, history_ref
    type(Closest_point) :: solver, fresh
    type(VonMises) :: vm
    type(Hill48) :: hill
    type(ten_3D2Osym) :: strain, strain_p
    type(Perzyna_overstress) :: pz
    type(Larger_overstress) :: larger
    type(Elasticity_linear) :: lin
    type(Elasticity_frozen) :: frozen
    type(ten_3D4O2sym) :: tangent, tangent_ref
    real(real64) :: dt
    integer :: status, status_ref, iters

    passed = .true.
    call strain%init(xx=0.08D0, yy=-0.024D0, zz=-0.01D0, xy=0.012D0, yz=-0.007D0, xz=0.004D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call lin%set_parameters(young=1000D0, poisson=0.3D0)
    call frozen%set(k_alg=1000D0/1.2D0, g_alg=1000D0/2.6D0, stress_hist=strain_p)
    hill = Hill48(f=0.25D0, g=0.35D0, h=0.65D0, l=1.3D0, m=1.7D0, n=1.45D0)
    pz%hard_law = Bilinear_hardening(y0=50D0, k=100D0)
    larger%hard_law = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    larger%eta = 1D0
    larger%n = 2D0
    dt = 1D-1

    call solver%init(elasticity=frozen, viscoplastic=pz, yield=vm)
    call solver%init(elasticity=lin, viscoplastic=larger, yield=hill)
    call fresh%init(elasticity=lin, viscoplastic=larger, yield=hill)
    call history%init(strain_p=strain_p, strain_pf=0D0)
    call history_ref%init(strain_p=strain_p, strain_pf=0D0)
    call solver%solve(strain=strain, history=history, dt=dt, status=status, iters=iters)
    call fresh%solve(strain=strain, history=history_ref, dt=dt, status=status_ref, iters=iters)
    call check_true(passed, "reinit: converged", status == STATUS_CONVERGED)
    call check_true(passed, "fresh: converged", status_ref == STATUS_CONVERGED)
    call check_true(passed, "reinit: plastic step", history%state_np1%strain_pf > 0D0)

    ! Same operations on the same data; Measured: 0
    call solver%tangent(strain=strain, history=history, dt=dt, tangent=tangent)
    call fresh%tangent(strain=strain, history=history_ref, dt=dt, tangent=tangent_ref)
    call check_approx(passed, "reinit: stress", history%state_np1%stress, history_ref%state_np1%stress, tol=0D0)
    call check_approx(passed, "reinit: plastic strain", history%state_np1%strain_p, &
                      history_ref%state_np1%strain_p, tol=0D0)
    call check_approx(passed, "reinit: equivalent plastic strain", history%state_np1%strain_pf, &
                      history_ref%state_np1%strain_pf, tol=0D0)
    call check_approx(passed, "reinit: tangent", tangent, tangent_ref, tol=0D0)
end subroutine
