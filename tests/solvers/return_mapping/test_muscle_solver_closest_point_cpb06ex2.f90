! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_solver_closest_point_cpb06ex2
    implicit none
    logical :: passed

    call test_closest_point_cpb06ex2_tension_rd(passed)
    if (.not. passed) STOP 1

    call test_closest_point_cpb06ex2_tension_td(passed)
    if (.not. passed) STOP 2

    call test_closest_point_cpb06ex2_compression_rd(passed)
    if (.not. passed) STOP 3

    call test_closest_point_cpb06ex2_shear_xy(passed)
    if (.not. passed) STOP 4

    call test_closest_point_cpb06ex2_shear_yz(passed)
    if (.not. passed) STOP 5

    call test_closest_point_cpb06ex2_shear_xz(passed)
    if (.not. passed) STOP 6

    print*, "Passed!", passed
    STOP 0
end program test_muscle_solver_closest_point_cpb06ex2


subroutine test_closest_point_cpb06ex2_tension_rd(passed)
    !! RD tension checks normalization, normal anisotropy and r0 (0.67 in Table 4).
    ! Plunkett et al. (2008), Eqs. (8)-(10), (22); associated flow, mpmath (70 digits).
    ! Steel, Table 2; library c44/c55/c66 = xy/yz/xz = paper C66/C44/C55.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_cpb06ex2, only : CPB06ex2
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Relative stress and plastic strain at solver tolerance 1e-5. Measured: 3.86e-8 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 1.95e-9 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! r or asymmetry error / max(1, abs(got), abs(expected)). Measured: 8.36e-8 or less.
    real(real64), parameter :: TOL_R = 1D-6
    ! Relative tangent error against centered differences. Measured: 3.65e-8 or less.
    real(real64), parameter :: TOL_TAN = 1D-6
    type(CPB06ex2) :: cpb
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
    r0 = 6.6907115270309783D-1
    call cpb%init(1.283D0,0.208D0,-0.421D0,0.208D0,0.535D0,-0.791D0, &
                  -0.421D0,-0.791D0,0.746D0,1.156D0,1.418D0,1.329D0,0.030D0, &
                  1.283D0,0.380D0,-0.001D0,0.380D0,2.145D0,0.792D0, &
                  -0.001D0,0.792D0,0.218D0,1.721D0,1.418D0,1.329D0,-0.027D0,5D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=7.4128254276130228D+1, yy=0.0D0, zz=0.0D0, xy=0.0D0, &
                              yz=0.0D0, xz=0.0D0)
    call expected_strain_p%init(xx=5.0000000000000000D-2, yy=-2.0043218397835294D-2, &
                                zz=-2.9956781602164706D-2, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    call strain%init(xx=1.2412825427613023D-1, yy=-4.2281694680674362D-2, &
                     zz=-5.2195257885003775D-2, xy=0.0D0, yz=0.0D0, xz=0.0D0)
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
end subroutine test_closest_point_cpb06ex2_tension_rd


subroutine test_closest_point_cpb06ex2_tension_td(passed)
    !! TD tension separates the normal coefficients and r90 (1.4 in Table 4) from RD.
    ! Plunkett et al. (2008), Eqs. (8)-(10), (22); associated flow, mpmath (70 digits).
    ! Steel, Table 2; library c44/c55/c66 = xy/yz/xz = paper C66/C44/C55.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_cpb06ex2, only : CPB06ex2
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Relative stress and plastic strain at solver tolerance 1e-5. Measured: 3.86e-8 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 1.95e-9 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! r or asymmetry error / max(1, abs(got), abs(expected)). Measured: 8.36e-8 or less.
    real(real64), parameter :: TOL_R = 1D-6
    ! Relative tangent error against centered differences. Measured: 3.65e-8 or less.
    real(real64), parameter :: TOL_TAN = 1D-6
    type(CPB06ex2) :: cpb
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: strain_pf
    integer :: status, iters
    type(ten_3D2Osym) :: width_projection ! t (x) t, transverse to the tensile axis
    real(real64) :: r_theta, r90

    passed = .true.
    r90 = 1.4000672140658833D+0
    call cpb%init(1.283D0,0.208D0,-0.421D0,0.208D0,0.535D0,-0.791D0, &
                  -0.421D0,-0.791D0,0.746D0,1.156D0,1.418D0,1.329D0,0.030D0, &
                  1.283D0,0.380D0,-0.001D0,0.380D0,2.145D0,0.792D0, &
                  -0.001D0,0.792D0,0.218D0,1.721D0,1.418D0,1.329D0,-0.027D0,5D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=0.0D0, yy=6.7745654672007246D+1, zz=0.0D0, xy=0.0D0, &
                              yz=0.0D0, xz=0.0D0)
    call expected_strain_p%init(xx=-3.1915217925684258D-2, yy=5.4710707745776864D-2, &
                                zz=-2.2795489820092606D-2, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    call strain%init(xx=-5.2238914327286432D-2, yy=1.2245636241778411D-1, &
                     zz=-4.3119186221694780D-2, xy=0.0D0, yz=0.0D0, xz=0.0D0)
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
    call width_projection%init(xx=1D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    r_theta = (width_projection .ddot. strain_p)/strain_p%zz()
    call check_approx(passed, "Lankford r-value", r_theta, r90, tol=TOL_R)
    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, tangent=tangent)
    call check_approx(passed, "Consistent tangent", tangent, numerical_tangent, tol=TOL_TAN)
end subroutine test_closest_point_cpb06ex2_tension_td


subroutine test_closest_point_cpb06ex2_compression_rd(passed)
    !! RD compression checks negative principal values, r0 in compression and the
    !! tension/compression asymmetry given by kp = 0.9.
    ! Plunkett et al. (2008), Eqs. (8)-(10), (22); associated flow, mpmath (70 digits).
    ! AZ31B, Table 1; unit out-of-plane shear factors are synthetic and inactive in RD.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_cpb06ex2, only : CPB06ex2
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Relative stress and plastic strain at solver tolerance 1e-5. Measured: 3.86e-8 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 1.95e-9 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! r or asymmetry error / max(1, abs(got), abs(expected)). Measured: 8.36e-8 or less.
    real(real64), parameter :: TOL_R = 1D-6
    ! Relative tangent error against centered differences. Measured: 3.65e-8 or less.
    real(real64), parameter :: TOL_TAN = 1D-6
    ! Swift stress at ep = 0.05: RD tensile normalization gives the tension_rd value.
    real(real64), parameter :: EXPECTED_TENSION_RD = 7.4128254276130228D+1
    type(CPB06ex2) :: cpb
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
    r0 = 1.9994508168136877D-1
    call cpb%init(0.771D0,-0.914D0,-1.312D0,-0.914D0,0.654D0,-1.209D0, &
                  -1.312D0,-1.209D0,-0.293D0,1.617D0,1D0,1D0,0.295D0, &
                  0.771D0,0.337D0,-0.764D0,0.337D0,0.646D0,-0.401D0, &
                  -0.764D0,-0.401D0,-0.022D0,0.499D0,1D0,1D0,0.900D0,4D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=-4.2304398571060402D+1, yy=0.0D0, zz=0.0D0, xy=0.0D0, &
                              yz=0.0D0, xz=0.0D0)
    call expected_strain_p%init(xx=-8.7612939528751382D-2, yy=1.4598815077332561D-2, &
                                zz=7.3014124451418822D-2, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    call strain%init(xx=-1.2991733809981178D-1, yy=2.7290134648650681D-2, &
                     zz=8.5705444022736942D-2, xy=0.0D0, yz=0.0D0, xz=0.0D0)
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
    call check_approx(passed, "Tension/compression ratio", &
                      EXPECTED_TENSION_RD/abs(stress%xx()), 1.7522587905750276D0, tol=TOL_R)
end subroutine test_closest_point_cpb06ex2_compression_rd


subroutine test_closest_point_cpb06ex2_shear_xy(passed)
    !! Pure xy shear isolates c44 and cp44 with nonzero associated normal flow.
    ! Plunkett et al. (2008), Eqs. (8)-(10); associated flow, mpmath (70 digits).
    ! Steel, Table 2; library c44/c55/c66 = xy/yz/xz = paper C66/C44/C55.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_cpb06ex2, only : CPB06ex2
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Relative stress and plastic strain at solver tolerance 1e-5. Measured: 3.86e-8 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 1.95e-9 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! Relative tangent error against centered differences. Measured: 3.65e-8 or less.
    real(real64), parameter :: TOL_TAN = 1D-6
    type(CPB06ex2) :: cpb
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: strain_pf
    integer :: status, iters

    passed = .true.
    call cpb%init(1.283D0,0.208D0,-0.421D0,0.208D0,0.535D0,-0.791D0, &
                  -0.421D0,-0.791D0,0.746D0,1.156D0,1.418D0,1.329D0,0.030D0, &
                  1.283D0,0.380D0,-0.001D0,0.380D0,2.145D0,0.792D0, &
                  -0.001D0,0.792D0,0.218D0,1.721D0,1.418D0,1.329D0,-0.027D0,5D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=0.0D0, yy=0.0D0, zz=0.0D0, xy=3.6510114242481117D+1, &
                              yz=0.0D0, xz=0.0D0)
    call expected_strain_p%init(xx=-8.9460540145381232D-4, yy=2.6958849818483326D-3, &
                                zz=-1.8012795803945203D-3, xy=5.0758711533884189D-2, yz=0.0D0, &
                                xz=0.0D0)
    call strain%init(xx=-8.9460540145381232D-4, yy=2.6958849818483326D-3, &
                     zz=-1.8012795803945203D-3, xy=9.8221860049109642D-2, yz=0.0D0, xz=0.0D0)
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
    ! a > 2: the zero principal value has zero classical curvature and a unique tangent.
    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, tangent=tangent)
    call check_approx(passed, "Consistent tangent", tangent, numerical_tangent, tol=TOL_TAN)
end subroutine test_closest_point_cpb06ex2_shear_xy


subroutine test_closest_point_cpb06ex2_shear_yz(passed)
    !! Pure yz shear isolates c55 and cp55, distinct from the xz factors.
    ! Plunkett et al. (2008), Eqs. (8)-(10); associated flow, mpmath (70 digits).
    ! Steel, Table 2; library c44/c55/c66 = xy/yz/xz = paper C66/C44/C55.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_cpb06ex2, only : CPB06ex2
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Relative stress and plastic strain at solver tolerance 1e-5. Measured: 3.86e-8 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 1.95e-9 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! Relative tangent error against centered differences. Measured: 3.65e-8 or less.
    real(real64), parameter :: TOL_TAN = 1D-6
    type(CPB06ex2) :: cpb
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: strain_pf
    integer :: status, iters

    passed = .true.
    call cpb%init(1.283D0,0.208D0,-0.421D0,0.208D0,0.535D0,-0.791D0, &
                  -0.421D0,-0.791D0,0.746D0,1.156D0,1.418D0,1.329D0,0.030D0, &
                  1.283D0,0.380D0,-0.001D0,0.380D0,2.145D0,0.792D0, &
                  -0.001D0,0.792D0,0.218D0,1.721D0,1.418D0,1.329D0,-0.027D0,5D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=0.0D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, &
                              yz=3.9571981385624044D+1, xz=0.0D0)
    call expected_strain_p%init(xx=-2.2514355494241961D-3, yy=3.5190447685552339D-3, &
                                zz=-1.2676092191310378D-3, xy=0.0D0, yz=4.6831275362332503D-2, &
                                xz=0.0D0)
    call strain%init(xx=-2.2514355494241961D-3, yy=3.5190447685552339D-3, &
                     zz=-1.2676092191310378D-3, xy=0.0D0, yz=9.8274851163643760D-2, xz=0.0D0)
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
    ! a > 2: the zero principal value has zero classical curvature and a unique tangent.
    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, tangent=tangent)
    call check_approx(passed, "Consistent tangent", tangent, numerical_tangent, tol=TOL_TAN)
end subroutine test_closest_point_cpb06ex2_shear_yz


subroutine test_closest_point_cpb06ex2_shear_xz(passed)
    !! Pure xz shear isolates c66 and cp66 with associated normal plastic flow.
    ! Plunkett et al. (2008), Eqs. (8)-(10); associated flow, mpmath (70 digits).
    ! Steel, Table 2; library c44/c55/c66 = xy/yz/xz = paper C66/C44/C55.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_cpb06ex2, only : CPB06ex2
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Relative stress and plastic strain at solver tolerance 1e-5. Measured: 3.86e-8 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 1.95e-9 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! Relative tangent error against centered differences. Measured: 3.65e-8 or less.
    real(real64), parameter :: TOL_TAN = 1D-6
    type(CPB06ex2) :: cpb
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: strain_pf
    integer :: status, iters

    passed = .true.
    call cpb%init(1.283D0,0.208D0,-0.421D0,0.208D0,0.535D0,-0.791D0, &
                  -0.421D0,-0.791D0,0.746D0,1.156D0,1.418D0,1.329D0,0.030D0, &
                  1.283D0,0.380D0,-0.001D0,0.380D0,2.145D0,0.792D0, &
                  -0.001D0,0.792D0,0.218D0,1.721D0,1.418D0,1.329D0,-0.027D0,5D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=0.0D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, &
                              xz=4.2222023780899093D+1)
    call expected_strain_p%init(xx=-7.5880977411385639D-4, yy=2.5542832986116338D-3, &
                                zz=-1.7954735244977774D-3, xy=0.0D0, yz=0.0D0, &
                                xz=4.3891935794456909D-2)
    call strain%init(xx=-7.5880977411385639D-4, yy=2.5542832986116338D-3, &
                     zz=-1.7954735244977774D-3, xy=0.0D0, yz=0.0D0, xz=9.8780566709625729D-2)
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
    ! a > 2: the zero principal value has zero classical curvature and a unique tangent.
    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, tangent=tangent)
    call check_approx(passed, "Consistent tangent", tangent, numerical_tangent, tol=TOL_TAN)
end subroutine test_closest_point_cpb06ex2_shear_xz
