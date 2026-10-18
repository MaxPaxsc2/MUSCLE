! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_solver_closest_point_cpb06
    implicit none
    logical :: passed

    call test_closest_point_cpb06_tension_rd(passed)
    if (.not. passed) STOP 1

    call test_closest_point_cpb06_tension_45(passed)
    if (.not. passed) STOP 2

    call test_closest_point_cpb06_compression_rd(passed)
    if (.not. passed) STOP 3

    call test_closest_point_cpb06_shear_xy(passed)
    if (.not. passed) STOP 4

    call test_closest_point_cpb06_shear_yz(passed)
    if (.not. passed) STOP 5

    call test_closest_point_cpb06_shear_xz(passed)
    if (.not. passed) STOP 6

    print*, "Passed!", passed
    STOP 0
end program test_muscle_solver_closest_point_cpb06


subroutine test_closest_point_cpb06_tension_rd(passed)
    !! RD tension checks the normal block and associated-flow r0.
    ! Cazacu et al. (2006, Eqs. 8, 9, 12); associated flow, mpmath (50 digits).
    ! Ti-6Al-4V: Tuninetti et al. (2013), Table 2, Wp = 1.857.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_cpb06, only : CPB06
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Relative stress and plastic strain at solver tolerance 1e-5. Measured: 3.97e-6 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 1.24e-8 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! r or asymmetry error / max(1, abs(got), abs(expected)). Measured: 6.22e-8 or less.
    real(real64), parameter :: TOL_R = 1D-6
    ! Relative tangent error against centered differences. Measured: 2.66e-8 or less.
    real(real64), parameter :: TOL_TAN = 1D-7
    type(CPB06) :: cpb
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: strain_pf
    integer :: status, iters
    type(ten_3D2Osym) :: width_projection ! t (x) t, transverse to the tensile axis
    real(real64) :: r_theta, r0

    passed = .true.
    r0 = 1.0180631085055389D+0
    call cpb%init(c11=1D0, c12=-2.373D0, c13=-2.364D0, &
                  c21=-2.373D0, c22=-1.838D0, c23=1.196D0, &
                  c31=-2.364D0, c32=1.196D0, c33=-2.444D0, &
                  c44=3.607D0, c55=3.607D0, c66=3.607D0, k=-0.136D0, a=2D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=7.4128254276130228D+1, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call expected_strain_p%init(xx=0.05D0, yy=-2.5223767884530076D-2, zz=-2.4776232115469924D-2, &
                                xy=0D0, yz=0D0, xz=0D0)
    call strain%init(xx=1.2412825427613023D-1, yy=-4.7462244167369144D-2, &
                     zz=-4.7014708398308993D-2, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=cpb)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    stress = history%state_np1%stress
    strain_p = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf
    call check_true(passed, "Converged status", status == STATUS_CONVERGED)
    call check_approx(passed, "Stress", stress, expected_stress, tol=TOL_STRESS)
    call check_approx(passed, "Equivalent plastic strain", strain_pf, 0.05D0, tol=TOL_PF)
    call check_approx(passed, "Plastic strain", strain_p, expected_strain_p, tol=TOL_STRESS)
    ! Newton iterations from the virgin state. Measured: 6 (Release and Debug).
    call check_true(passed, "Iterations <= 6", iters <= 6)
    call width_projection%init(xx=0D0, yy=1D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    r_theta = (width_projection .ddot. strain_p)/strain_p%zz()
    call check_approx(passed, "Lankford r-value", r_theta, r0, tol=TOL_R)
    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, tangent=tangent)
    call check_approx(passed, "Consistent tangent", tangent, numerical_tangent, tol=TOL_TAN)
end subroutine test_closest_point_cpb06_tension_rd


subroutine test_closest_point_cpb06_tension_45(passed)
    !! 45-degree tension couples normal anisotropy and tensorial xy flow in r45.
    ! Cazacu et al. (2006, Eqs. 8, 9, 12); associated flow, mpmath (50 digits).
    ! Ti-6Al-4V: Tuninetti et al. (2013), Table 2, Wp = 1.857.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_cpb06, only : CPB06
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Relative stress and plastic strain at solver tolerance 1e-5. Measured: 3.97e-6 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 1.24e-8 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! r or asymmetry error / max(1, abs(got), abs(expected)). Measured: 6.22e-8 or less.
    real(real64), parameter :: TOL_R = 1D-6
    ! Relative tangent error against centered differences. Measured: 2.66e-8 or less.
    real(real64), parameter :: TOL_TAN = 1D-7
    type(CPB06) :: cpb
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: strain_pf
    integer :: status, iters
    type(ten_3D2Osym) :: width_projection ! t (x) t, transverse to the tensile axis
    real(real64) :: r_theta, r45

    passed = .true.
    r45 = 1.3726367893107158D+0
    call cpb%init(c11=1D0, c12=-2.373D0, c13=-2.364D0, &
                  c21=-2.373D0, c22=-1.838D0, c23=1.196D0, &
                  c31=-2.364D0, c32=1.196D0, c33=-2.444D0, &
                  c44=3.607D0, c55=3.607D0, c66=3.607D0, k=-0.136D0, a=2D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=3.8898652987499279D+1, yy=3.8898652987499279D+1, zz=0D0, &
                              xy=3.8898652987499279D+1, yz=0D0, xz=0D0)
    call expected_strain_p%init(xx=1.3233458141441924D-2, yy=6.8462760103472110D-3, &
                                zz=-2.0079734151789135D-2, xy=3.7602048892219136D-2, yz=0D0, &
                                xz=0D0)
    call strain%init(xx=4.0462515232691419D-2, yy=3.4075333101596706D-2, &
                     zz=-4.3418925944288702D-2, xy=8.8170297775968199D-2, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=cpb)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    stress = history%state_np1%stress
    strain_p = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf
    call check_true(passed, "Converged status", status == STATUS_CONVERGED)
    call check_approx(passed, "Stress", stress, expected_stress, tol=TOL_STRESS)
    call check_approx(passed, "Equivalent plastic strain", strain_pf, 0.05D0, tol=TOL_PF)
    call check_approx(passed, "Plastic strain", strain_p, expected_strain_p, tol=TOL_STRESS)
    ! Newton iterations from the virgin state. Measured: 6 (Release and Debug).
    call check_true(passed, "Iterations <= 6", iters <= 6)
    call width_projection%init(xx=0.5D0, yy=0.5D0, zz=0D0, xy=-0.5D0, yz=0D0, xz=0D0)
    r_theta = (width_projection .ddot. strain_p)/strain_p%zz()
    call check_approx(passed, "Lankford r-value", r_theta, r45, tol=TOL_R)
    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, tangent=tangent)
    call check_approx(passed, "Consistent tangent", tangent, numerical_tangent, tol=TOL_TAN)
end subroutine test_closest_point_cpb06_tension_45


subroutine test_closest_point_cpb06_compression_rd(passed)
    !! RD compression reverses the principal values and checks strength asymmetry.
    ! Cazacu et al. (2006, Eqs. 8, 9, 12); associated flow, mpmath (50 digits).
    ! Ti-6Al-4V: Tuninetti et al. (2013), Table 2, Wp = 1.857.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_cpb06, only : CPB06
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Relative stress and plastic strain at solver tolerance 1e-5. Measured: 3.97e-6 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 1.24e-8 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! r or asymmetry error / max(1, abs(got), abs(expected)). Measured: 6.22e-8 or less.
    real(real64), parameter :: TOL_R = 1D-6
    ! Relative tangent error against centered differences. Measured: 2.66e-8 or less.
    real(real64), parameter :: TOL_TAN = 1D-7
    ! Expected axial stress from tension_rd at the same equivalent plastic strain.
    real(real64), parameter :: EXPECTED_TENSION_RD = 7.4128254276130228D+1
    type(CPB06) :: cpb
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: strain_pf
    integer :: status, iters

    passed = .true.
    call cpb%init(c11=1D0, c12=-2.373D0, c13=-2.364D0, &
                  c21=-2.373D0, c22=-1.838D0, c23=1.196D0, &
                  c31=-2.364D0, c32=1.196D0, c33=-2.444D0, &
                  c44=3.607D0, c55=3.607D0, c66=3.607D0, k=-0.136D0, a=2D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=-7.8641758456969284D+1, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call expected_strain_p%init(xx=-4.7130338722455750D-2, yy=2.3886247961267747D-2, &
                                zz=2.3244090761188003D-2, xy=0D0, yz=0D0, xz=0D0)
    call strain%init(xx=-1.2577209717942503D-1, yy=4.7478775498358532D-2, &
                     zz=4.6836618298278789D-2, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=cpb)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    stress = history%state_np1%stress
    strain_p = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf
    call check_true(passed, "Converged status", status == STATUS_CONVERGED)
    call check_approx(passed, "Stress", stress, expected_stress, tol=TOL_STRESS)
    call check_approx(passed, "Equivalent plastic strain", strain_pf, 0.05D0, tol=TOL_PF)
    call check_approx(passed, "Plastic strain", strain_p, expected_strain_p, tol=TOL_STRESS)
    ! Newton iterations from the virgin state. Measured: 6 (Release and Debug).
    call check_true(passed, "Iterations <= 6", iters <= 6)
    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, tangent=tangent)
    call check_approx(passed, "Consistent tangent", tangent, numerical_tangent, tol=TOL_TAN)
    call check_approx(passed, "Tension/compression ratio", &
                      EXPECTED_TENSION_RD/abs(stress%xx()), 0.94260677444911500D0, tol=TOL_R)
end subroutine test_closest_point_cpb06_compression_rd


subroutine test_closest_point_cpb06_shear_xy(passed)
    !! Pure xy shear isolates c44 and retains the normal plastic flow for nonzero k.
    ! Cazacu et al. (2006, Eqs. 8, 9, 12); associated flow, mpmath (50 digits).
    ! Ti-6Al-4V: Tuninetti et al. (2013), Table 2, Wp = 1.857.
    ! Shear variant: c44/c55/c66 = xy/yz/xz (paper c66/c44/c55); B is unchanged.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_cpb06, only : CPB06
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Relative stress and plastic strain at solver tolerance 1e-5. Measured: 3.97e-6 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 1.24e-8 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    type(CPB06) :: cpb
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    real(real64) :: strain_pf
    integer :: status, iters

    passed = .true.
    call cpb%init(c11=1D0, c12=-2.373D0, c13=-2.364D0, &
                  c21=-2.373D0, c22=-1.838D0, c23=1.196D0, &
                  c31=-2.364D0, c32=1.196D0, c33=-2.444D0, &
                  c44=3.607D0, c55=2.900D0, c66=4.100D0, k=-0.136D0, a=2D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=0D0, yy=0D0, zz=0D0, xy=4.2980297569717324D+1, yz=0D0, xz=0D0)
    call expected_strain_p%init(xx=2.8018662612699186D-3, yy=-6.2581752687156823D-3, &
                                zz=3.4563090074457636D-3, xy=4.3117578557877909D-2, yz=0D0, &
                                xz=0D0)
    call strain%init(xx=2.8018662612699186D-3, yy=-6.2581752687156823D-3, &
                     zz=3.4563090074457636D-3, xy=9.8991965398510431D-2, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=cpb)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    stress = history%state_np1%stress
    strain_p = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf
    call check_true(passed, "Converged status", status == STATUS_CONVERGED)
    call check_approx(passed, "Stress", stress, expected_stress, tol=TOL_STRESS)
    call check_approx(passed, "Equivalent plastic strain", strain_pf, 0.05D0, tol=TOL_PF)
    call check_approx(passed, "Plastic strain", strain_p, expected_strain_p, tol=TOL_STRESS)
    ! Newton iterations from the virgin state. Measured: 6 (Release and Debug).
    call check_true(passed, "Iterations <= 6", iters <= 6)
    ! Sigma has a zero principal value. For k /= 0 the one-sided curvatures of
    ! (|lambda| - k*lambda)**2 are 2*(1-k)**2 and 2*(1+k)**2: no unique tangent.
    ! Tangent comparison is omitted; centered differences differed by up to 2e-2.
end subroutine test_closest_point_cpb06_shear_xy


subroutine test_closest_point_cpb06_shear_yz(passed)
    !! Pure yz shear isolates c55 and retains the normal plastic flow for nonzero k.
    ! Cazacu et al. (2006, Eqs. 8, 9, 12); associated flow, mpmath (50 digits).
    ! Ti-6Al-4V: Tuninetti et al. (2013), Table 2, Wp = 1.857.
    ! Shear variant: c44/c55/c66 = xy/yz/xz (paper c66/c44/c55); B is unchanged.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_cpb06, only : CPB06
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Relative stress and plastic strain at solver tolerance 1e-5. Measured: 3.97e-6 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 1.24e-8 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    type(CPB06) :: cpb
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    real(real64) :: strain_pf
    integer :: status, iters

    passed = .true.
    call cpb%init(c11=1D0, c12=-2.373D0, c13=-2.364D0, &
                  c21=-2.373D0, c22=-1.838D0, c23=1.196D0, &
                  c31=-2.364D0, c32=1.196D0, c33=-2.444D0, &
                  c44=3.607D0, c55=2.900D0, c66=4.100D0, k=-0.136D0, a=2D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=5.3458597701369099D+1, xz=0D0)
    call expected_strain_p%init(xx=-8.0703964016221279D-3, yy=5.0024965036953617D-3, &
                                zz=3.0678998979267662D-3, xy=0D0, yz=3.4666198452410850D-2, &
                                xz=0D0)
    call strain%init(xx=-8.0703964016221279D-3, yy=5.0024965036953617D-3, &
                     zz=3.0678998979267662D-3, xy=0D0, yz=1.0416237546419068D-1, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=cpb)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    stress = history%state_np1%stress
    strain_p = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf
    call check_true(passed, "Converged status", status == STATUS_CONVERGED)
    call check_approx(passed, "Stress", stress, expected_stress, tol=TOL_STRESS)
    call check_approx(passed, "Equivalent plastic strain", strain_pf, 0.05D0, tol=TOL_PF)
    call check_approx(passed, "Plastic strain", strain_p, expected_strain_p, tol=TOL_STRESS)
    ! Newton iterations from the virgin state. Measured: 6 (Release and Debug).
    call check_true(passed, "Iterations <= 6", iters <= 6)
    ! Sigma has a zero principal value. For k /= 0 the one-sided curvatures of
    ! (|lambda| - k*lambda)**2 are 2*(1-k)**2 and 2*(1+k)**2: no unique tangent.
    ! Tangent comparison is omitted; centered differences differed by up to 2e-2.
end subroutine test_closest_point_cpb06_shear_yz


subroutine test_closest_point_cpb06_shear_xz(passed)
    !! Pure xz shear isolates c66 and retains the normal plastic flow for nonzero k.
    ! Cazacu et al. (2006, Eqs. 8, 9, 12); associated flow, mpmath (50 digits).
    ! Ti-6Al-4V: Tuninetti et al. (2013), Table 2, Wp = 1.857.
    ! Shear variant: c44/c55/c66 = xy/yz/xz (paper c66/c44/c55); B is unchanged.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_cpb06, only : CPB06
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Relative stress and plastic strain at solver tolerance 1e-5. Measured: 3.97e-6 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 1.24e-8 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    type(CPB06) :: cpb
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    real(real64) :: strain_pf
    integer :: status, iters

    passed = .true.
    call cpb%init(c11=1D0, c12=-2.373D0, c13=-2.364D0, &
                  c21=-2.373D0, c22=-1.838D0, c23=1.196D0, &
                  c31=-2.364D0, c32=1.196D0, c33=-2.444D0, &
                  c44=3.607D0, c55=2.900D0, c66=4.100D0, k=-0.136D0, a=2D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=3.7812178861943997D+1)
    call expected_strain_p%init(xx=3.4658862183654102D-3, yy=4.0628656990233761D-3, &
                                zz=-7.5287519173887863D-3, xy=0D0, yz=0D0, &
                                xz=4.9010832294787754D-2)
    call strain%init(xx=3.4658862183654102D-3, yy=4.0628656990233761D-3, &
                     zz=-7.5287519173887863D-3, xy=0D0, yz=0D0, xz=9.8166664815314950D-2)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=cpb)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    stress = history%state_np1%stress
    strain_p = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf
    call check_true(passed, "Converged status", status == STATUS_CONVERGED)
    call check_approx(passed, "Stress", stress, expected_stress, tol=TOL_STRESS)
    call check_approx(passed, "Equivalent plastic strain", strain_pf, 0.05D0, tol=TOL_PF)
    call check_approx(passed, "Plastic strain", strain_p, expected_strain_p, tol=TOL_STRESS)
    ! Newton iterations from the virgin state. Measured: 6 (Release and Debug).
    call check_true(passed, "Iterations <= 6", iters <= 6)
    ! Sigma has a zero principal value. For k /= 0 the one-sided curvatures of
    ! (|lambda| - k*lambda)**2 are 2*(1-k)**2 and 2*(1+k)**2: no unique tangent.
    ! Tangent comparison is omitted; centered differences differed by up to 2e-2.
end subroutine test_closest_point_cpb06_shear_xz
