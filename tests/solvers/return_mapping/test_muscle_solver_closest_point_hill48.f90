! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_solver_closest_point_hill48
    implicit none
    logical :: passed

    call test_closest_point_hill48_tension_rd(passed)
    if (.not. passed) STOP 1

    call test_closest_point_hill48_tension_45(passed)
    if (.not. passed) STOP 2

    call test_closest_point_hill48_tension_td(passed)
    if (.not. passed) STOP 3

    call test_closest_point_hill48_shear_xy(passed)
    if (.not. passed) STOP 4

    call test_closest_point_hill48_shear_yz(passed)
    if (.not. passed) STOP 5

    call test_closest_point_hill48_shear_xz(passed)
    if (.not. passed) STOP 6

    print*, "Passed!", passed
    STOP 0
end program test_muscle_solver_closest_point_hill48

subroutine test_closest_point_hill48_tension_rd(passed)
    !! RD tension checks the normal block and r0 = H/G.
    ! Hill (1948, Eq. 5), associated flow at equivalent plastic strain 0.05.
    ! Expected values evaluated independently with mpmath (50 digits).
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_hill48, only : Hill48
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Stress and plastic strain, relative, at solver tolerance 1e-5. Measured: 5.98e-09 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 5.37e-10 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! r error / max(1, abs(r), abs(r_in)); absolute below 1. Measured: 6.86e-09 or less.
    real(real64), parameter :: TOL_R = 1D-6
    ! Relative tangent error against centered differences. Measured: 4.33e-08 or less.
    real(real64), parameter :: TOL_TAN = 1D-7
    type(Hill48) :: h48
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: r0, r45, r90, strain_pf
    integer :: status, iters
    type(ten_3D2Osym) :: width_projection ! t (x) t, transverse to the tensile axis
    real(real64) :: r_theta

    passed = .true.
    r0 = 0.2115D0
    r45 = 1.5769D0
    r90 = 0.6923D0
    h48 = Hill48(f=r0/(r90*(1D0 + r0)), g=1D0/(1D0 + r0), h=r0/(1D0 + r0), &
                 l=1.25D0, m=1.75D0, &
                 n=(r0 + r90)*(1D0 + 2D0*r45)/(2D0*r90*(1D0 + r0)))
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=7.4128254276130235D+01, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call expected_strain_p%init(xx=0.05D0, yy=-8.7288485348741233D-03, zz=-4.1271151465125874D-02, &
                                xy=0D0, yz=0D0, xz=0D0)
    call strain%init(xx=1.2412825427613022D-01, yy=-3.0967324817713193D-02, &
                     zz=-6.3509627747964939D-02, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=h48)
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
end subroutine test_closest_point_hill48_tension_rd


subroutine test_closest_point_hill48_tension_45(passed)
    !! 45-degree tension activates N and tensorial xy flow in r45.
    ! Hill (1948, Eq. 5), associated flow at equivalent plastic strain 0.05.
    ! Expected values evaluated independently with mpmath (50 digits).
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_hill48, only : Hill48
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Stress and plastic strain, relative, at solver tolerance 1e-5. Measured: 5.98e-09 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 5.37e-10 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! r error / max(1, abs(r), abs(r_in)); absolute below 1. Measured: 6.86e-09 or less.
    real(real64), parameter :: TOL_R = 1D-6
    ! Relative tangent error against centered differences. Measured: 4.33e-08 or less.
    real(real64), parameter :: TOL_TAN = 1D-7
    type(Hill48) :: h48
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: r0, r45, r90, strain_pf
    integer :: status, iters
    type(ten_3D2Osym) :: width_projection ! t (x) t, transverse to the tensile axis
    real(real64) :: r_theta

    passed = .true.
    r0 = 0.2115D0
    r45 = 1.5769D0
    r90 = 0.6923D0
    h48 = Hill48(f=r0/(r90*(1D0 + r0)), g=1D0/(1D0 + r0), h=r0/(1D0 + r0), &
                 l=1.25D0, m=1.75D0, &
                 n=(r0 + r90)*(1D0 + 2D0*r45)/(2D0*r90*(1D0 + r0)))
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=3.1455217937819118D+01, yy=3.1455217937819118D+01, zz=0D0, &
                              xy=3.1455217937819118D+01, yz=0D0, xz=0D0)
    call expected_strain_p%init(xx=1.7512796929554885D-02, yy=5.3502189088557826D-03, &
                                zz=-2.2863015838410668D-02, xy=4.7484197594795115D-02, yz=0D0, &
                                xz=0D0)
    call strain%init(xx=3.9531449486028269D-02, yy=2.7368871465329167D-02, &
                     zz=-4.1736146601102142D-02, xy=8.8375980913959970D-02, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=h48)
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
end subroutine test_closest_point_hill48_tension_45


subroutine test_closest_point_hill48_tension_td(passed)
    !! TD tension activates F and checks r90 = H/F.
    ! Hill (1948, Eq. 5), associated flow at equivalent plastic strain 0.05.
    ! Expected values evaluated independently with mpmath (50 digits).
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_hill48, only : Hill48
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Stress and plastic strain, relative, at solver tolerance 1e-5. Measured: 5.98e-09 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 5.37e-10 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! r error / max(1, abs(r), abs(r_in)); absolute below 1. Measured: 6.86e-09 or less.
    real(real64), parameter :: TOL_R = 1D-6
    ! Relative tangent error against centered differences. Measured: 4.33e-08 or less.
    real(real64), parameter :: TOL_TAN = 1D-7
    type(Hill48) :: h48
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: r0, r45, r90, strain_pf
    integer :: status, iters
    type(ten_3D2Osym) :: width_projection ! t (x) t, transverse to the tensile axis
    real(real64) :: r_theta

    passed = .true.
    r0 = 0.2115D0
    r45 = 1.5769D0
    r90 = 0.6923D0
    h48 = Hill48(f=r0/(r90*(1D0 + r0)), g=1D0/(1D0 + r0), h=r0/(1D0 + r0), &
                 l=1.25D0, m=1.75D0, &
                 n=(r0 + r90)*(1D0 + 2D0*r45)/(2D0*r90*(1D0 + r0)))
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=0D0, yy=1.1347463240701646D+02, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call expected_strain_p%init(xx=-1.3362015448815358D-02, yy=3.2662918884920164D-02, &
                                zz=-1.9300903436104806D-02, xy=0D0, yz=0D0, xz=0D0)
    call strain%init(xx=-4.7404405170920295D-02, yy=1.4613755129193662D-01, &
                     zz=-5.3343293158209744D-02, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=h48)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    stress = history%state_np1%stress
    strain_p = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf
    call check_true(passed, "Converged status", status == STATUS_CONVERGED)
    call check_approx(passed, "Stress", stress, expected_stress, tol=TOL_STRESS)
    call check_approx(passed, "Equivalent plastic strain", strain_pf, 0.05D0, tol=TOL_PF)
    call check_approx(passed, "Plastic strain", strain_p, expected_strain_p, tol=TOL_STRESS)
    ! Newton iterations from the virgin state; TD needs one more. Measured: 7 (Release, Debug).
    call check_true(passed, "Iterations <= 7", iters <= 7)
    call width_projection%init(xx=1D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    r_theta = (width_projection .ddot. strain_p)/strain_p%zz()
    call check_approx(passed, "Lankford r-value", r_theta, r90, tol=TOL_R)
    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, tangent=tangent)
    call check_approx(passed, "Consistent tangent", tangent, numerical_tangent, tol=TOL_TAN)
end subroutine test_closest_point_hill48_tension_td


subroutine test_closest_point_hill48_shear_xy(passed)
    !! Pure xy shear isolates N, distinct from L and M.
    ! Hill (1948, Eq. 5), associated flow at equivalent plastic strain 0.05.
    ! Expected values evaluated independently with mpmath (50 digits).
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_hill48, only : Hill48
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Stress and plastic strain, relative, at solver tolerance 1e-5. Measured: 5.98e-09 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 5.37e-10 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! Relative tangent error against centered differences. Measured: 4.33e-08 or less.
    real(real64), parameter :: TOL_TAN = 1D-7
    type(Hill48) :: h48
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: r0, r45, r90, strain_pf
    integer :: status, iters

    passed = .true.
    r0 = 0.2115D0
    r45 = 1.5769D0
    r90 = 0.6923D0
    h48 = Hill48(f=r0/(r90*(1D0 + r0)), g=1D0/(1D0 + r0), h=r0/(1D0 + r0), &
                 l=1.25D0, m=1.75D0, &
                 n=(r0 + r90)*(1D0 + 2D0*r45)/(2D0*r90*(1D0 + r0)))
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=0D0, yy=0D0, zz=0D0, xy=3.5037546580065133D+01, yz=0D0, xz=0D0)
    call expected_strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=5.2892012649014948D-02, yz=0D0, xz=0D0)
    call strain%init(xx=0D0, yy=0D0, zz=0D0, xy=9.8440823203099631D-02, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=h48)
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
end subroutine test_closest_point_hill48_shear_xy


subroutine test_closest_point_hill48_shear_yz(passed)
    !! Pure yz shear isolates L, distinct from N and M.
    ! Hill (1948, Eq. 5), associated flow at equivalent plastic strain 0.05.
    ! Expected values evaluated independently with mpmath (50 digits).
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_hill48, only : Hill48
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Stress and plastic strain, relative, at solver tolerance 1e-5. Measured: 5.98e-09 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 5.37e-10 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! Relative tangent error against centered differences. Measured: 4.33e-08 or less.
    real(real64), parameter :: TOL_TAN = 1D-7
    type(Hill48) :: h48
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: r0, r45, r90, strain_pf
    integer :: status, iters

    passed = .true.
    r0 = 0.2115D0
    r45 = 1.5769D0
    r90 = 0.6923D0
    h48 = Hill48(f=r0/(r90*(1D0 + r0)), g=1D0/(1D0 + r0), h=r0/(1D0 + r0), &
                 l=1.25D0, m=1.75D0, &
                 n=(r0 + r90)*(1D0 + 2D0*r45)/(2D0*r90*(1D0 + r0)))
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=4.6882824496937552D+01, xz=0D0)
    call expected_strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=3.9528470752104743D-02, xz=0D0)
    call strain%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=1.0047614259812356D-01, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=h48)
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
end subroutine test_closest_point_hill48_shear_yz


subroutine test_closest_point_hill48_shear_xz(passed)
    !! Pure xz shear isolates M, distinct from N and L.
    ! Hill (1948, Eq. 5), associated flow at equivalent plastic strain 0.05.
    ! Expected values evaluated independently with mpmath (50 digits).
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_hill48, only : Hill48
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Stress and plastic strain, relative, at solver tolerance 1e-5. Measured: 5.98e-09 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 5.37e-10 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! Relative tangent error against centered differences. Measured: 4.33e-08 or less.
    real(real64), parameter :: TOL_TAN = 1D-7
    type(Hill48) :: h48
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: r0, r45, r90, strain_pf
    integer :: status, iters

    passed = .true.
    r0 = 0.2115D0
    r45 = 1.5769D0
    r90 = 0.6923D0
    h48 = Hill48(f=r0/(r90*(1D0 + r0)), g=1D0/(1D0 + r0), h=r0/(1D0 + r0), &
                 l=1.25D0, m=1.75D0, &
                 n=(r0 + r90)*(1D0 + 2D0*r45)/(2D0*r90*(1D0 + r0)))
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=3.9623218597277095D+01)
    call expected_strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=4.6770717334674271D-02)
    call strain%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=9.8280901511134500D-02)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=h48)
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
end subroutine test_closest_point_hill48_shear_xz
