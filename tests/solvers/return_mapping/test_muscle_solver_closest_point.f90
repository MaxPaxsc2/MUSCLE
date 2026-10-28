program test_muscle_solver_closest_point
    implicit none
    logical :: passed

    call test_closest_point_vonmises_uniaxial_tensile(passed)
    if (.not. passed) STOP 1

    call test_closest_point_vonmises_zero_strain(passed)
    if (.not. passed) STOP 2

    call test_closest_point_vonmises_elastic_strain(passed)
    if (.not. passed) STOP 3

    ! call test_closest_point_druckerPrager_uniaxial_tensile(passed)
    ! if (.not. passed) STOP 4
 
    call test_closest_point_packed_shear(passed)
    if (.not. passed) STOP 5

    call test_closest_point_keeps_state_n(passed)
    if (.not. passed) STOP 6

    call test_closest_point_vonmises_elastic_status(passed)
    if (.not. passed) STOP 7

    call test_closest_point_vonmises_first_update(passed)
    if (.not. passed) STOP 8

    call test_closest_point_vonmises_iteration_limit(passed)
    if (.not. passed) STOP 9

    call test_closest_point_vonmises_rollback_retry(passed)
    if (.not. passed) STOP 10

    call test_closest_point_reinit_other_types(passed)
    if (.not. passed) STOP 11

    print*, "Passed!", passed
end program test_muscle_solver_closest_point

subroutine test_closest_point_vonmises_uniaxial_tensile(passed)
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point
    implicit none

    real(real64), parameter :: EPS=1e-5
    logical, intent(out) :: passed
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(VonMises) :: vm
    type(ten_3D2Osym) :: strain, strain_p, stress
    type(Swift_hardening) :: sw
    type(Elasticity_linear) :: elas
    real(real64) :: strain_pf
    logical :: error

    type(ten_3D2Osym) :: expected_stress, expected_strain_plastic
    real(real64) :: expected_strain_effective
    integer :: status, iters
    type(ten_3D4O2sym) :: tangent, numerical_tangent

    passed = .False.
    strain_pf = 0D0
    error = .False.
    

    call expected_stress%init(xx=90D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    expected_strain_effective = 0.3485784401D0
    call expected_strain_plastic%init(xx=expected_strain_effective,        &
                                      yy=-0.5D0*expected_strain_effective, &
                                      zz=-0.5D0*expected_strain_effective, &
                                      xy=0D0, yz=0D0, xz=0D0)
                                      
    call strain%init(xx=0.43857844D0, yy=-0.20128922D0, zz=-0.20128922D0, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)

    call history%init(strain_p=strain_p, strain_pf=strain_pf)
    call solver%init(elasticity=elas, hardening=sw, yield=vm)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    
    stress    = history%state_np1%stress
    strain_p  = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf

    passed = stress%is_approx(expected_stress, tol=EPS)
    if (.not. passed) print*, "Error: Stress is no equal", new_line('A'),          &
                              "Expected:", expected_stress, new_line('A'),  &
                              "Actual Value:", stress, new_line('A'),       &
                              "Difference", stress - expected_stress
    if (.not. passed) return
    
    passed = abs(strain_pf - expected_strain_effective) < EPS
    if (.not. passed) print*, "Effective plastic Strain is not equal", new_line('A'), &
                              "Expected:", expected_strain_effective, new_line('A'),  &
                              "Actual Value:", strain_pf, new_line('A'),              &
                              "Difference", strain_pf - expected_strain_effective
    if (.not. passed) return

    passed = strain_p%is_approx(expected_strain_plastic, tol=EPS)
    if (.not. passed) print*, "Plastic Strain is not equal", new_line('A'), &
                              "Expected:", expected_strain_plastic, new_line('A'), &
                              "Actual Value:", strain_p, new_line('A'),            & 
                              "Difference:", strain_p - expected_strain_plastic
    if (.not. passed) return
    
    passed = iters .le. 6
    if(.not. passed) print*, "Iterations are greater than expected", new_line('A'), &
                             "Iters:", iters, " Maximum:", 6 , new_line('A'), &
                             "Status code: ", status
    if(.not. passed) return


    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, tangent=tangent)

    passed = tangent%is_approx(numerical_tangent, tol=1.0D-7)
    if (.not. passed) print*, "Tangent is no equal", new_line('A'),          &
                              "Analitical:", tangent, new_line('A'),  &
                              "Numerical:", numerical_tangent, new_line('A'),       &
                              "Difference", tangent - numerical_tangent
    if (.not. passed) return


    return
end subroutine


subroutine test_closest_point_vonmises_zero_strain(passed)
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point
    implicit none

    real(real64), parameter :: EPS=1e-8
    logical, intent(out) :: passed

    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(VonMises) :: vm
    type(ten_3D2Osym) :: strain, strain_p, stress
    type(Swift_hardening) :: sw
    type(Elasticity_linear) :: elas
    real(real64) :: strain_pf
    logical :: error

    type(ten_3D2Osym) :: expected_stress, expected_strain_plastic
    real(real64) :: expected_strain_effective
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    integer :: status, iters


    passed = .False.
    strain_pf = 0D0
    error = .False.

    expected_strain_effective = 0.0D0
    call expected_stress%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call expected_strain_plastic%init(xx=0D0, yy=0D0, zz=0D0, &
                                      xy=0D0, yz=0D0, xz=0D0)
                                      
    call strain%init(xx=0.0D0, yy=0.0D0, zz=0.0D0, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)

    call history%init(strain_p=strain_p, strain_pf=strain_pf)
    call solver%init(elasticity=elas, hardening=sw, yield=vm)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    
    stress    = history%state_np1%stress
    strain_p  = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf

    passed = stress .approx. expected_stress
    if (.not. passed) print*, "Stress is no equal", new_line('A'),          &
                              "Expected:", expected_stress, new_line('A'),  &
                              "Actual Value:", stress, new_line('A'),       &
                              "Difference", stress - expected_stress
    if (.not. passed) return
    
    passed = abs(strain_pf - expected_strain_effective) < EPS
    if (.not. passed) print*, "Effective plastic Strain is not equal", new_line('A'), &
                              "Expected:", expected_strain_effective, new_line('A'),  &
                              "Actual Value:", strain_pf, new_line('A'),              &
                              "Difference", strain_pf - expected_strain_effective
    if (.not. passed) return

    passed = strain_p .approx. expected_strain_plastic
    if (.not. passed) print*, "Plastic Strain is not equal", new_line('A'), &
                              "Expected:", expected_strain_plastic, new_line('A'), &
                              "Actual Value:", strain_p, new_line('A'),            & 
                              "Difference:", strain_p - expected_strain_plastic
    if (.not. passed) return


    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, tangent=tangent)

    passed = tangent .approx. numerical_tangent
    if (.not. passed) print*, "Tangent is no equal", new_line('A'),          &
                              "Analitical:", tangent, new_line('A'),  &
                              "Numerical:", numerical_tangent, new_line('A'),       &
                              "Difference", tangent - numerical_tangent
    if (.not. passed) return
    

    return
end subroutine



subroutine test_closest_point_vonmises_elastic_strain(passed)
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point
    implicit none

    real(real64), parameter :: EPS=1e-8
    logical, intent(out) :: passed

    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(VonMises) :: vm
    type(ten_3D2Osym) :: strain, strain_p, stress
    type(Swift_hardening) :: sw
    type(Elasticity_linear) :: elas
    real(real64) :: strain_pf
    logical :: error

    type(ten_3D2Osym) :: expected_stress, expected_strain_plastic
    real(real64) :: expected_strain_effective
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    integer :: status, iters


    passed = .False.
    strain_pf = 0D0
    error = .False.

    expected_strain_effective = 0.0D0
    call expected_stress%init(xx=1D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call expected_strain_plastic%init(xx=0D0, yy=0D0, zz=0D0, &
                                      xy=0D0, yz=0D0, xz=0D0)
                                      
    call strain%init(xx=0.001D0, yy=-0.0003D0, zz=-0.0003D0, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)

    call history%init(strain_p=strain_p, strain_pf=strain_pf)
    call solver%init(elasticity=elas, hardening=sw, yield=vm)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    
    stress    = history%state_np1%stress
    strain_p  = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf

    passed = stress .approx. expected_stress
    if (.not. passed) print*, "Stress is no equal", new_line('A'),          &
                              "Expected:", expected_stress, new_line('A'),  &
                              "Actual Value:", stress, new_line('A'),       &
                              "Difference", stress - expected_stress
    if (.not. passed) return
    
    passed = abs(strain_pf - expected_strain_effective) < EPS
    if (.not. passed) print*, "Effective plastic Strain is not equal", new_line('A'), &
                              "Expected:", expected_strain_effective, new_line('A'),  &
                              "Actual Value:", strain_pf, new_line('A'),              &
                              "Difference", strain_pf - expected_strain_effective
    if (.not. passed) return

    passed = strain_p .approx. expected_strain_plastic
    if (.not. passed) print*, "Plastic Strain is not equal", new_line('A'), &
                              "Expected:", expected_strain_plastic, new_line('A'), &
                              "Actual Value:", strain_p, new_line('A'),            & 
                              "Difference:", strain_p - expected_strain_plastic
    if (.not. passed) return


    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, tangent=tangent)

    passed = tangent .approx. numerical_tangent
    if (.not. passed) print*, "Tangent is no equal", new_line('A'),          &
                              "Analitical:", tangent, new_line('A'),  &
                              "Numerical:", numerical_tangent, new_line('A'),       &
                              "Difference", tangent - numerical_tangent
    if (.not. passed) return

    return
end subroutine

subroutine test_closest_point_druckerPrager_uniaxial_tensile(passed)
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_bilinear
    use muscle_yield_druckerprager
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point
    implicit none

    real(real64), parameter :: EPS=1e-5
    logical, intent(out) :: passed
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: strain, strain_p, stress
    type(Bilinear_hardening) :: bilinear
    type(Elasticity_linear) :: elas
    real(real64) :: strain_pf
    real(real64) :: beta, K_ratio

    type(ten_3D2Osym) :: expected_stress, expected_strain_plastic
    real(real64) :: expected_strain_effective
    integer :: status, iters
    real(real64) :: strain_xx, strain_yy

    passed = .False.
    strain_pf = 0D0

    call elas%set_parameters(young=210D3, poisson=0.25D0)
    bilinear = Bilinear_hardening(y0=300D0, k=1000D0)
    beta = 16.0D0
    K_ratio = 0.85D0

    ! ---------------------------------------------------
    ! --------  Tensile type  ---------------------------
    ! ---------------------------------------------------

    strain_xx = 0.09531D0
    strain_yy = -0.03666D0
    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_TENSION)


    call expected_stress%init(xx=393.4D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call expected_strain_plastic%init(xx=0.09344D0,  &
                                      yy=-0.03619D0, &
                                      zz=-0.03619D0, &
                                      xy=0D0, yz=0D0, xz=0D0)
    expected_strain_effective = 0.09344D0

    call strain%init(xx=strain_xx, yy=strain_yy, zz=strain_yy, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)

    call history%init(strain_p=strain_p, strain_pf=strain_pf)
    call solver%init(elasticity=elas, hardening=bilinear, yield=dp)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)

    stress    = history%state_np1%stress
    strain_p  = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf

    ! passed = stress%is_approx(expected_stress, tol=EPS)
    ! if (.not. passed) print*, "Error: Stress is no equal", iters, new_line('A'),          &
    !                           "Expected:", expected_stress, new_line('A'),  &
    !                           "Actual Value:", stress, new_line('A'),       &
    !                           "Difference", stress - expected_stress
    ! if (.not. passed) return

    passed = abs(strain_pf - expected_strain_effective) < EPS
    if (.not. passed) print*, "Effective plastic Strain is not equal", new_line('A'), &
                              "Expected:", expected_strain_effective, new_line('A'),  &
                              "Actual Value:", strain_pf, new_line('A'),              &
                              "Difference", strain_pf - expected_strain_effective
    if (.not. passed) return

    passed = strain_p%is_approx(expected_strain_plastic, tol=EPS)
    if (.not. passed) print*, "Plastic Strain is not equal", new_line('A'), &
                              "Expected:", expected_strain_plastic, new_line('A'), &
                              "Actual Value:", strain_p, new_line('A'),            & 
                              "Difference:", strain_p - expected_strain_plastic
    if (.not. passed) return

    passed = iters .le. 6
    if(.not. passed) print*, "Iterations are greater than expected", new_line('A'), &
                             "Iters:", iters, " Maximum:", 6 , new_line('A'), &
                             "Status code: ", status
    if(.not. passed) return

    return
end subroutine test_closest_point_druckerPrager_uniaxial_tensile

! FE boundary: the FE program reads the result through pack_to_fea and commits only after global
! equilibrium. Both cases use von Mises with Swift hardening and a strain with all six components
! nonzero, written as a multiple of the yield strain eps_y along a fixed direction.

subroutine test_closest_point_packed_shear(passed)
    ! The FE program reads the plastic strain from slots 7-12. Its shears must be the tensorial
    ! components (eps_xy, not gamma_xy = 2*eps_xy): with them the packed stress is the elastic
    ! law of the total minus the plastic strain.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_hard_swift
    use muscle_yield_vonmises
    use muscle_elasticity_linear
    use muscle_plastic_history
    use muscle_solver_closest_point
    implicit none
    logical, intent(out) :: passed

    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: hardening
    type(VonMises) :: vm
    type(Closest_point) :: solver
    type(Plastic_material_history) :: history
    type(ten_3D2Osym) :: strain, packed_stress, packed_strain_p
    real(real64) :: hsv(13), eps_y
    integer :: status
    real(real64), parameter :: TOL = 1.0D-10
    ! Newton tolerance of the solver
    real(real64), parameter :: TOL_NW = 1.0D-5

    call elas%set_parameters(young=70000.0D0, poisson=0.33D0)
    hardening = Swift_hardening(k=500.0D0, n=0.2D0, e0=0.01D0)
    call solver%init(elasticity=elas, hardening=hardening, yield=vm)
    eps_y = hardening%stress(0.0D0)/70000.0D0

    ! Five times the yield strain: plastic, with large plastic shears
    call strain%init(xx=4.0D0*eps_y, yy=-1.5D0*eps_y, zz=-2.5D0*eps_y, &
                     xy=2.0D0*eps_y, yz=-1.0D0*eps_y, xz=1.5D0*eps_y)
    call solver%solve(strain=strain, history=history, status=status)
    if (status /= STATUS_CONVERGED) then
        print*, "FAIL: Packed shear case did not converge, status:", status
        passed = .false.
        return
    end if

    hsv = 0.0D0
    call history%pack_to_fea(hsv)
    call packed_stress%init(hsv(1:6))
    call packed_strain_p%init(hsv(7:12))

    passed = packed_stress%is_approx(elas%stress(strain - packed_strain_p), tol=TOL)
    if (.not. passed) then
        print*, "FAIL: Packed stress is not the elastic law of the packed (tensorial) plastic strain"
        return
    end if

    ! One step from the virgin state: eps_p = dgamma*n with the von Mises normal n of unit
    ! equivalent norm, so slot 13 is sqrt(2/3 eps_p:eps_p)
    passed = abs(hsv(13) - sqrt(2.0D0/3.0D0*(packed_strain_p .ddot. packed_strain_p))) < TOL_NW
    if (.not. passed) then
        print*, "FAIL: Slot 13 is not the equivalent plastic strain of slots 7-12"
        return
    end if
end subroutine test_closest_point_packed_shear


subroutine test_closest_point_keeps_state_n(passed)
    ! The FE program commits only after global equilibrium. Before that, neither a failed attempt
    ! (rolled back and retried) nor a converged one may write the committed state_n.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_hard_swift
    use muscle_yield_vonmises
    use muscle_elasticity_linear
    use muscle_plastic_history
    use muscle_solver_closest_point
    implicit none
    logical, intent(out) :: passed

    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: hardening
    type(VonMises) :: vm
    type(Closest_point) :: solver, one_iteration
    type(Plastic_material_history) :: history
    type(Plastic_state) :: state_before
    type(ten_3D2Osym) :: strain_1, strain_2
    real(real64) :: eps_y
    integer :: status
    real(real64), parameter :: TOL = 1.0D-12

    call elas%set_parameters(young=70000.0D0, poisson=0.33D0)
    hardening = Swift_hardening(k=500.0D0, n=0.2D0, e0=0.01D0)
    call solver%init(elasticity=elas, hardening=hardening, yield=vm)
    ! One Newton iteration is not enough for step 2: forces a failed attempt
    call one_iteration%init(elasticity=elas, hardening=hardening, yield=vm, iter_nw=1)
    eps_y = hardening%stress(0.0D0)/70000.0D0

    ! Step 1, converged and committed, so state_n is a plastic state and not the virgin one
    call strain_1%init(xx=1.6D0*eps_y, yy=-0.6D0*eps_y, zz=-1.0D0*eps_y, &
                       xy=0.8D0*eps_y, yz=-0.4D0*eps_y, xz=0.6D0*eps_y)
    call solver%solve(strain=strain_1, history=history, status=status)
    if (status /= STATUS_CONVERGED) then
        print*, "FAIL: Step 1 did not converge, status:", status
        passed = .false.
        return
    end if
    call history%commit()
    state_before = history%state_n

    ! Step 2: failed attempt with a single Newton iteration, rollback, converged retry
    call strain_2%init(xx=3.2D0*eps_y, yy=-1.2D0*eps_y, zz=-2.0D0*eps_y, &
                       xy=1.6D0*eps_y, yz=-0.8D0*eps_y, xz=1.2D0*eps_y)
    call one_iteration%solve(strain=strain_2, history=history, status=status)
    if (status /= STATUS_NONCONVERGED) then
        print*, "FAIL: Forced attempt was expected not to converge, status:", status
        passed = .false.
        return
    end if
    call history%rollback()
    call solver%solve(strain=strain_2, history=history, status=status)
    if (status /= STATUS_CONVERGED) then
        print*, "FAIL: Retry of step 2 did not converge, status:", status
        passed = .false.
        return
    end if

    passed = history%state_n%stress%is_approx(state_before%stress, tol=TOL) .and. &
             history%state_n%strain_p%is_approx(state_before%strain_p, tol=TOL) .and. &
             abs(history%state_n%strain_pf - state_before%strain_pf) <= TOL*state_before%strain_pf
    if (.not. passed) then
        print*, "FAIL: An attempt before commit changed state_n"
        return
    end if
end subroutine test_closest_point_keeps_state_n


subroutine test_closest_point_vonmises_elastic_status(passed)
    !! An elastic step returns STATUS_ELASTIC_CASE and zero Newton iterations.
    ! All six strain components are nonzero, with three distinct tensorial shears; the trial
    ! von Mises stress (17.8) is below the initial yield stress y0 = 200.
    use muscle_tensors, only : ten_3D2Osym
    use muscle_hard_bilinear, only : Bilinear_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_ELASTIC_CASE
    use muscle_test_utils, only : check_true
    implicit none
    logical, intent(out) :: passed
    type(Elasticity_linear) :: elas
    type(Bilinear_hardening) :: linear
    type(VonMises) :: vm
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain
    integer :: status, iters

    passed = .true.
    call elas%set_parameters(young=70000D0, poisson=0.33D0)
    linear = Bilinear_hardening(y0=200D0, k=1000D0)
    call solver%init(elasticity=elas, hardening=linear, yield=vm)
    call strain%init(xx=2D-4, yy=-6D-5, zz=-1D-4, xy=8D-5, yz=-4D-5, xz=6D-5)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    call check_true(passed, "Elastic status", status == STATUS_ELASTIC_CASE)
    call check_true(passed, "Elastic iterations", iters == 0)
end subroutine test_closest_point_vonmises_elastic_status


subroutine test_closest_point_vonmises_first_update(passed)
    !! One Newton update (iter_nw = 1) with linear hardening equals the closed-form radial return.
    ! Coupled normal and shear strains include pressure in a plastic trial state.
    ! Simo and Hughes (1998), Box 3.1 and Remark 3.3.2: with linear hardening of slope H,
    ! Newton reaches the radial return in one update; dg = (q_trial - y0)/(3G + H).
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors, only : ten_3D2Osym, iden_2O, operator(*), operator(+), &
                              operator(.ddot.), operator(.dev.)
    use muscle_hard_bilinear, only : Bilinear_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_NONCONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Relative plastic strain error against radial return. Measured: 1.35e-16 or less.
    real(real64), parameter :: TOL_STRAIN_P = 1D-12
    ! Relative equivalent strain error (tol = TOL_PF*dg, as dg < 1). Measured: 1.12e-16 or less.
    real(real64), parameter :: TOL_PF = 1D-12
    type(Elasticity_linear) :: elas
    type(Bilinear_hardening) :: linear
    type(VonMises) :: vm
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, trial, deviator
    real(real64) :: g, lame ! Lame moduli G and lambda.
    real(real64) :: q_trial, dg ! Trial equivalent stress and equivalent plastic strain increment.
    integer :: status

    passed = .true.
    call elas%set_parameters(young=70000D0, poisson=0.33D0)
    linear = Bilinear_hardening(y0=200D0, k=1000D0)
    call solver%init(elasticity=elas, hardening=linear, yield=vm, iter_nw=1)
    call strain%init(xx=2D-4, yy=-6D-5, zz=-1D-4, xy=8D-5, yz=-4D-5, xz=6D-5)
    strain = 20D0*strain
    g = 70000D0/(2D0*(1D0 + 0.33D0))
    lame = 70000D0*0.33D0/((1D0 + 0.33D0)*(1D0 - 2D0*0.33D0))
    trial = 2D0*g*strain + lame*(strain .ddot. iden_2O())*iden_2O()
    deviator = .dev. trial
    q_trial = sqrt(1.5D0*(deviator .ddot. deviator))
    dg = (q_trial - 200D0)/(3D0*g + 1000D0)
    call solver%solve(strain=strain, history=history, status=status)
    ! With iter_nw = 1 the solver stops before the iteration that checks convergence.
    call check_true(passed, "One update status", status == STATUS_NONCONVERGED)
    call check_approx(passed, "Plastic strain after one update", history%state_np1%strain_p, &
                      (1.5D0*dg/q_trial)*deviator, tol=TOL_STRAIN_P)
    call check_approx(passed, "Equivalent plastic strain after one update", &
                      history%state_np1%strain_pf, dg, tol=TOL_PF*dg)
end subroutine test_closest_point_vonmises_first_update


subroutine test_closest_point_vonmises_iteration_limit(passed)
    !! iter_nw = k - 1 does not converge; iter_nw = k converges with the default history.
    ! Nonlinear Swift hardening needs several iterations, including the convergence check.
    ! Measured: k = 4 (Release and Debug).
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors, only : ten_3D2Osym
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED, STATUS_NONCONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Same operations in the same order give identical packed histories. Measured: 0.
    real(real64), parameter :: TOL_SAME_HISTORY = 0D0
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: hardening
    type(VonMises) :: vm
    type(Plastic_material_history) :: base, short, exact
    type(Closest_point) :: solver, limited
    type(ten_3D2Osym) :: strain
    real(real64) :: eps_y, base_hsv(13), exact_hsv(13)
    integer :: status, iters, k

    passed = .true.
    call elas%set_parameters(young=70000D0, poisson=0.33D0)
    hardening = Swift_hardening(k=500D0, n=0.2D0, e0=0.01D0)
    eps_y = hardening%stress(0D0)/70000D0
    call strain%init(xx=2D0*eps_y, yy=-0.75D0*eps_y, zz=-1.25D0*eps_y, &
                     xy=eps_y, yz=-0.5D0*eps_y, xz=0.75D0*eps_y)
    call solver%init(elasticity=elas, hardening=hardening, yield=vm)
    call solver%solve(strain=strain, history=base, status=status, iters=k)
    call check_true(passed, "Default convergence", status == STATUS_CONVERGED)
    call check_true(passed, "Multiple Newton iterations", k >= 2)
    call limited%init(elasticity=elas, hardening=hardening, yield=vm, iter_nw=k - 1)
    call limited%solve(strain=strain, history=short, status=status)
    call check_true(passed, "Below limit status", status == STATUS_NONCONVERGED)
    call limited%init(elasticity=elas, hardening=hardening, yield=vm, iter_nw=k)
    call limited%solve(strain=strain, history=exact, status=status, iters=iters)
    call check_true(passed, "At limit status", status == STATUS_CONVERGED)
    call check_true(passed, "At limit iterations", iters == k)
    call base%pack_to_fea(base_hsv)
    call exact%pack_to_fea(exact_hsv)
    call check_approx(passed, "Same history at the limit", exact_hsv, base_hsv, &
                      tol=TOL_SAME_HISTORY)
end subroutine test_closest_point_vonmises_iteration_limit


subroutine test_closest_point_vonmises_rollback_retry(passed)
    !! A failed step, rollback and a smaller retry end with the history of a clean path.
    ! A failed large step precedes a smaller retry from a nonvirgin plastic state.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors, only : ten_3D2Osym, operator(*)
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED, STATUS_NONCONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Same operations in the same order give identical packed histories. Measured: 0.
    real(real64), parameter :: TOL_SAME_HISTORY = 0D0
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: hardening
    type(VonMises) :: vm
    type(Plastic_material_history) :: base, clean, retry
    type(Closest_point) :: solver, limited
    type(ten_3D2Osym) :: strain, strain_3, strain_4
    real(real64) :: eps_y, clean_hsv(13), retry_hsv(13)
    integer :: status

    passed = .true.
    call elas%set_parameters(young=70000D0, poisson=0.33D0)
    hardening = Swift_hardening(k=500D0, n=0.2D0, e0=0.01D0)
    eps_y = hardening%stress(0D0)/70000D0
    call strain%init(xx=2D0*eps_y, yy=-0.75D0*eps_y, zz=-1.25D0*eps_y, &
                     xy=eps_y, yz=-0.5D0*eps_y, xz=0.75D0*eps_y)
    call solver%init(elasticity=elas, hardening=hardening, yield=vm)
    call limited%init(elasticity=elas, hardening=hardening, yield=vm, iter_nw=1)
    call solver%solve(strain=strain, history=base, status=status)
    call check_true(passed, "Starting step converged", status == STATUS_CONVERGED)
    call base%commit()
    clean = base
    retry = base
    strain_3 = 1.5D0*strain
    strain_4 = 2D0*strain
    call solver%solve(strain=strain_3, history=clean, status=status)
    call check_true(passed, "Clean smaller step converged", status == STATUS_CONVERGED)
    call clean%commit()
    call solver%solve(strain=strain_4, history=clean, status=status)
    call check_true(passed, "Clean final step converged", status == STATUS_CONVERGED)
    call clean%commit()
    call limited%solve(strain=strain_4, history=retry, status=status)
    call check_true(passed, "Forced attempt failed", status == STATUS_NONCONVERGED)
    call retry%rollback()
    call base%pack_to_fea(clean_hsv)
    call retry%pack_to_fea(retry_hsv)
    call check_approx(passed, "History after rollback", retry_hsv, clean_hsv, tol=TOL_SAME_HISTORY)
    call solver%solve(strain=strain_3, history=retry, status=status)
    call check_true(passed, "Smaller retry converged", status == STATUS_CONVERGED)
    call retry%commit()
    call solver%solve(strain=strain_4, history=retry, status=status)
    call check_true(passed, "Final retry converged", status == STATUS_CONVERGED)
    call retry%commit()
    call clean%pack_to_fea(clean_hsv)
    call retry%pack_to_fea(retry_hsv)
    call check_approx(passed, "Same final history", retry_hsv, clean_hsv, tol=TOL_SAME_HISTORY)
    clean%state_np1 = clean%state_n
    retry%state_np1 = retry%state_n
    call clean%pack_to_fea(clean_hsv)
    call retry%pack_to_fea(retry_hsv)
    call check_approx(passed, "Same committed history", retry_hsv, clean_hsv, tol=TOL_SAME_HISTORY)
end subroutine test_closest_point_vonmises_rollback_retry


subroutine test_closest_point_reinit_other_types(passed)
    ! A second init with components of larger dynamic types (Elasticity_frozen -> Elasticity_linear,
    ! Swift -> a larger extension of Swift, VonMises -> Hill48) must give the same step as a solver
    ! initialized with the second set only. Without reallocating each component, the new data
    ! overflow the old storage.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_hard_swift
    use muscle_yield_vonmises
    use muscle_yield_hill48
    use muscle_elasticity_linear
    use muscle_visco_frozen
    use muscle_plastic_history
    use muscle_solver_closest_point
    use muscle_test_utils, only: check_approx, check_true
    implicit none

    type, extends(Swift_hardening) :: Larger_hardening
        real(real64) :: extra(32) = 0.0D0
    end type Larger_hardening

    logical, intent(out) :: passed
    type(Elasticity_frozen) :: frozen
    type(Elasticity_linear) :: lin
    type(Swift_hardening) :: swift
    type(Larger_hardening) :: larger
    type(VonMises) :: vm
    type(Hill48) :: hill
    type(Closest_point) :: solver, fresh
    type(Plastic_material_history) :: history, history_ref
    type(ten_3D2Osym) :: strain
    type(ten_3D4O2sym) :: tangent, tangent_ref
    integer :: status, status_ref

    passed = .true.
    strain = 0.0D0
    call frozen%set(k_alg=1000.0D0/1.2D0, g_alg=1000.0D0/2.6D0, stress_hist=strain)
    call lin%set_parameters(young=1000.0D0, poisson=0.3D0)
    swift = Swift_hardening(k=100.0D0, n=0.1D0, e0=1.0D-4)
    larger%k = 100.0D0
    larger%n = 0.1D0
    larger%e0 = 1.0D-4
    hill = Hill48(f=0.25D0, g=0.35D0, h=0.65D0, l=1.3D0, m=1.7D0, n=1.45D0)
    call solver%init(elasticity=frozen, hardening=swift, yield=vm)
    call solver%init(elasticity=lin, hardening=larger, yield=hill)
    call fresh%init(elasticity=lin, hardening=larger, yield=hill)

    call strain%init(xx=0.08D0, yy=-0.024D0, zz=-0.01D0, xy=0.012D0, yz=-0.007D0, xz=0.004D0)
    call solver%solve(strain=strain, history=history, status=status)
    call fresh%solve(strain=strain, history=history_ref, status=status_ref)
    call check_true(passed, "reinit: converged", status == STATUS_CONVERGED)
    call check_true(passed, "fresh: converged", status_ref == STATUS_CONVERGED)
    call check_true(passed, "reinit: plastic step", history%state_np1%strain_pf > 0.0D0)

    ! Same operations on the same data; Measured: 0
    call solver%tangent(strain=strain, history=history, tangent=tangent)
    call fresh%tangent(strain=strain, history=history_ref, tangent=tangent_ref)
    call check_approx(passed, "reinit: stress", history%state_np1%stress, history_ref%state_np1%stress, tol=0.0D0)
    call check_approx(passed, "reinit: plastic strain", history%state_np1%strain_p, &
                      history_ref%state_np1%strain_p, tol=0.0D0)
    call check_approx(passed, "reinit: equivalent plastic strain", history%state_np1%strain_pf, &
                      history_ref%state_np1%strain_pf, tol=0.0D0)
    call check_approx(passed, "reinit: tangent", tangent, tangent_ref, tol=0.0D0)
end subroutine test_closest_point_reinit_other_types
