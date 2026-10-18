! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_solver_closest_point_yld2004
    implicit none
    logical :: passed

    call test_closest_point_yld2004_tension_rd(passed)
    if (.not. passed) STOP 1

    call test_closest_point_yld2004_tension_45(passed)
    if (.not. passed) STOP 2

    call test_closest_point_yld2004_tension_td(passed)
    if (.not. passed) STOP 3

    call test_closest_point_yld2004_shear_xy(passed)
    if (.not. passed) STOP 4

    call test_closest_point_yld2004_shear_yz(passed)
    if (.not. passed) STOP 5

    call test_closest_point_yld2004_shear_xz(passed)
    if (.not. passed) STOP 6

    print*, "Passed!", passed
    STOP 0
end program test_muscle_solver_closest_point_yld2004

subroutine test_closest_point_yld2004_tension_rd(passed)
    !! RD tension checks the normal blocks of C' and C'' and the associated r0.
    ! Barlat et al. (2005), Table 2, Eqs. 14-15; AA2090-T3, a = 8.
    ! Library (xy, yz, xz) = table (c66, c44, c55); tensorial shears.
    ! Independent associated-flow reference, mpmath (50 digits), outside this solver.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_yld2004, only : Yld2004
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Stress and plastic strain, relative, at solver tolerance 1e-5. Measured: 2.69e-8 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 9.58e-10 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! r error / max(1, abs(r), abs(r_in)); absolute below 1. Measured: 1.85e-8 or less.
    real(real64), parameter :: TOL_R = 1D-6
    ! Tangent vs centered differences of solves stopped at solver tolerance (a = 8).
    ! Finite differences need 1e-6 even with exact Hessian. Measured: 2.61e-7 or less.
    real(real64), parameter :: TOL_TAN = 1D-6
    type(Yld2004) :: yld
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: r0, strain_pf
    integer :: status, iters
    type(ten_3D2Osym) :: width_projection ! t (x) t, transverse to the tensile axis
    real(real64) :: r_theta

    passed = .true.
    r0 = 2.4506026150375929D-01
    call yld%init(cp12=-0.069888D0, cp13=0.936408D0, cp21=0.079143D0, cp23=1.003060D0, &
                  cp31=0.524741D0, cp32=1.363180D0, &
                  cp44=0.954322D0, cp55=1.023770D0, cp66=1.069060D0, &
                  cpp12=0.981171D0, cpp13=0.476741D0, cpp21=0.575316D0, cpp23=0.866827D0, &
                  cpp31=1.145010D0, cpp32=-0.079294D0, &
                  cpp44=1.404620D0, cpp55=1.051660D0, cpp66=1.147100D0, a=8D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=7.4178581528627859D+01, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call expected_strain_p%init(xx=4.9966076964899762D-02, yy=-9.8346242876199350D-03, &
                                zz=-4.0131452677279827D-02, xy=0D0, yz=0D0, xz=0D0)
    call strain%init(xx=1.2414465849352763D-01, yy=-3.2088198746208295D-02, &
                     zz=-6.2385027135868187D-02, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=yld)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    stress = history%state_np1%stress
    strain_p = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf
    call check_true(passed, "Converged status", status == STATUS_CONVERGED)
    call check_approx(passed, "Stress", stress, expected_stress, tol=TOL_STRESS)
    call check_approx(passed, "Equivalent plastic strain", strain_pf, 0.05D0, tol=TOL_PF)
    call check_approx(passed, "Plastic strain", strain_p, expected_strain_p, tol=TOL_STRESS)
    ! Newton iterations from the virgin state. Measured: 6 (Release, Debug, both Hessians).
    call check_true(passed, "Iterations <= 6", iters <= 6)
    call width_projection%init(xx=0D0, yy=1D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    r_theta = (width_projection .ddot. strain_p)/strain_p%zz()
    call check_approx(passed, "Lankford r-value", r_theta, r0, tol=TOL_R)
    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, tangent=tangent)
    call check_approx(passed, "Consistent tangent", tangent, numerical_tangent, tol=TOL_TAN)
end subroutine test_closest_point_yld2004_tension_rd


subroutine test_closest_point_yld2004_tension_45(passed)
    !! 45-degree tension activates library c'44/c''44 and tensorial xy flow in r45.
    ! Barlat et al. (2005), Table 2, Eqs. 14-15; AA2090-T3, a = 8.
    ! Library (xy, yz, xz) = table (c66, c44, c55); tensorial shears.
    ! Independent associated-flow reference, mpmath (50 digits), outside this solver.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_yld2004, only : Yld2004
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Stress and plastic strain, relative, at solver tolerance 1e-5. Measured: 2.69e-8 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 9.58e-10 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! r error / max(1, abs(r), abs(r_in)); absolute below 1. Measured: 1.85e-8 or less.
    real(real64), parameter :: TOL_R = 1D-6
    ! Tangent vs centered differences of solves stopped at solver tolerance (a = 8).
    ! Finite differences need 1e-6 even with exact Hessian. Measured: 2.61e-7 or less.
    real(real64), parameter :: TOL_TAN = 1D-6
    type(Yld2004) :: yld
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: r45, strain_pf
    integer :: status, iters
    type(ten_3D2Osym) :: width_projection ! t (x) t, transverse to the tensile axis
    real(real64) :: r_theta

    passed = .true.
    r45 = 1.5445248728743324D+00
    call yld%init(cp12=-0.069888D0, cp13=0.936408D0, cp21=0.079143D0, cp23=1.003060D0, &
                  cp31=0.524741D0, cp32=1.363180D0, &
                  cp44=0.954322D0, cp55=1.023770D0, cp66=1.069060D0, &
                  cpp12=0.981171D0, cpp13=0.476741D0, cpp21=0.575316D0, cpp23=0.866827D0, &
                  cpp31=1.145010D0, cpp32=-0.079294D0, &
                  cpp44=1.404620D0, cpp55=1.051660D0, cpp66=1.147100D0, a=8D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=3.0304826648644006D+01, yy=3.0304826648644006D+01, zz=0D0, &
                              xy=3.0304826648644006D+01, yz=0D0, xz=0D0)
    call expected_strain_p%init(xx=4.5381962634522014D-03, yy=1.9494653784340569D-02, &
                                zz=-2.4032850047792771D-02, xy=4.9135759688771409D-02, yz=0D0, &
                                xz=0D0)
    call strain%init(xx=2.5751574917503007D-02, yy=4.0708032438391377D-02, &
                     zz=-4.2215746036979175D-02, xy=8.8532034332008624D-02, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=yld)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    stress = history%state_np1%stress
    strain_p = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf
    call check_true(passed, "Converged status", status == STATUS_CONVERGED)
    call check_approx(passed, "Stress", stress, expected_stress, tol=TOL_STRESS)
    call check_approx(passed, "Equivalent plastic strain", strain_pf, 0.05D0, tol=TOL_PF)
    call check_approx(passed, "Plastic strain", strain_p, expected_strain_p, tol=TOL_STRESS)
    ! Newton iterations from the virgin state. Measured: 6 (Release, Debug, both Hessians).
    call check_true(passed, "Iterations <= 6", iters <= 6)
    call width_projection%init(xx=0.5D0, yy=0.5D0, zz=0D0, xy=-0.5D0, yz=0D0, xz=0D0)
    r_theta = (width_projection .ddot. strain_p)/strain_p%zz()
    call check_approx(passed, "Lankford r-value", r_theta, r45, tol=TOL_R)
    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, tangent=tangent)
    call check_approx(passed, "Consistent tangent", tangent, numerical_tangent, tol=TOL_TAN)
end subroutine test_closest_point_yld2004_tension_45


subroutine test_closest_point_yld2004_tension_td(passed)
    !! TD tension checks the other rows of C' and C'' and the associated r90.
    ! Barlat et al. (2005), Table 2, Eqs. 14-15; AA2090-T3, a = 8.
    ! Library (xy, yz, xz) = table (c66, c44, c55); tensorial shears.
    ! Independent associated-flow reference, mpmath (50 digits), outside this solver.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_yld2004, only : Yld2004
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Stress and plastic strain, relative, at solver tolerance 1e-5. Measured: 2.69e-8 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 9.58e-10 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! r error / max(1, abs(r), abs(r_in)); absolute below 1. Measured: 1.85e-8 or less.
    real(real64), parameter :: TOL_R = 1D-6
    ! Tangent vs centered differences of solves stopped at solver tolerance (a = 8).
    ! Finite differences need 1e-6 even with exact Hessian. Measured: 2.61e-7 or less.
    real(real64), parameter :: TOL_TAN = 1D-6
    type(Yld2004) :: yld
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: r90, strain_pf
    integer :: status, iters
    type(ten_3D2Osym) :: width_projection ! t (x) t, transverse to the tensile axis
    real(real64) :: r_theta

    passed = .true.
    r90 = 6.8494596476625125D-01
    call yld%init(cp12=-0.069888D0, cp13=0.936408D0, cp21=0.079143D0, cp23=1.003060D0, &
                  cp31=0.524741D0, cp32=1.363180D0, &
                  cp44=0.954322D0, cp55=1.023770D0, cp66=1.069060D0, &
                  cpp12=0.981171D0, cpp13=0.476741D0, cpp21=0.575316D0, cpp23=0.866827D0, &
                  cpp31=1.145010D0, cpp32=-0.079294D0, &
                  cpp44=1.404620D0, cpp55=1.051660D0, cpp66=1.147100D0, a=8D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=0D0, yy=6.7189605089885617D+01, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call expected_strain_p%init(xx=-2.2424461463235235D-02, yy=5.5163484125976145D-02, &
                                zz=-3.2739022662740910D-02, xy=0D0, yz=0D0, xz=0D0)
    call strain%init(xx=-4.2581342990200922D-02, yy=1.2235308921586177D-01, &
                     zz=-5.2895904189706598D-02, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=yld)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    stress = history%state_np1%stress
    strain_p = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf
    call check_true(passed, "Converged status", status == STATUS_CONVERGED)
    call check_approx(passed, "Stress", stress, expected_stress, tol=TOL_STRESS)
    call check_approx(passed, "Equivalent plastic strain", strain_pf, 0.05D0, tol=TOL_PF)
    call check_approx(passed, "Plastic strain", strain_p, expected_strain_p, tol=TOL_STRESS)
    ! Newton iterations from the virgin state. Measured: 6 (Release, Debug, both Hessians).
    call check_true(passed, "Iterations <= 6", iters <= 6)
    call width_projection%init(xx=1D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    r_theta = (width_projection .ddot. strain_p)/strain_p%zz()
    call check_approx(passed, "Lankford r-value", r_theta, r90, tol=TOL_R)
    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, tangent=tangent)
    call check_approx(passed, "Consistent tangent", tangent, numerical_tangent, tol=TOL_TAN)
end subroutine test_closest_point_yld2004_tension_td


subroutine test_closest_point_yld2004_shear_xy(passed)
    !! Pure xy shear isolates library c'44/c''44, distinct from yz and xz.
    ! Barlat et al. (2005), Table 2, Eqs. 14-15; AA2090-T3, a = 8.
    ! Library (xy, yz, xz) = table (c66, c44, c55); tensorial shears.
    ! Independent associated-flow reference, mpmath (50 digits), outside this solver.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_yld2004, only : Yld2004
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Stress and plastic strain, relative, at solver tolerance 1e-5. Measured: 2.69e-8 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 9.58e-10 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! Tangent vs centered differences of solves stopped at solver tolerance (a = 8).
    ! Finite differences need 1e-6 even with exact Hessian. Measured: 2.61e-7 or less.
    real(real64), parameter :: TOL_TAN = 1D-6
    type(Yld2004) :: yld
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: strain_pf
    integer :: status, iters

    passed = .true.
    call yld%init(cp12=-0.069888D0, cp13=0.936408D0, cp21=0.079143D0, cp23=1.003060D0, &
                  cp31=0.524741D0, cp32=1.363180D0, &
                  cp44=0.954322D0, cp55=1.023770D0, cp66=1.069060D0, &
                  cpp12=0.981171D0, cpp13=0.476741D0, cpp21=0.575316D0, cpp23=0.866827D0, &
                  cpp31=1.145010D0, cpp32=-0.079294D0, &
                  cpp44=1.404620D0, cpp55=1.051660D0, cpp66=1.147100D0, a=8D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=0D0, yy=0D0, zz=0D0, xy=3.4198389998627306D+01, yz=0D0, xz=0D0)
    call expected_strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=5.4189871423117930D-02, yz=0D0, &
                                xz=0D0)
    call strain%init(xx=0D0, yy=0D0, zz=0D0, xy=9.8647778421333418D-02, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=yld)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    stress = history%state_np1%stress
    strain_p = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf
    call check_true(passed, "Converged status", status == STATUS_CONVERGED)
    call check_approx(passed, "Stress", stress, expected_stress, tol=TOL_STRESS)
    call check_approx(passed, "Equivalent plastic strain", strain_pf, 0.05D0, tol=TOL_PF)
    call check_approx(passed, "Plastic strain", strain_p, expected_strain_p, tol=TOL_STRESS)
    ! Newton iterations from the virgin state. Measured: 6 (Release, Debug, both Hessians).
    call check_true(passed, "Iterations <= 6", iters <= 6)
    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, tangent=tangent)
    call check_approx(passed, "Consistent tangent", tangent, numerical_tangent, tol=TOL_TAN)
end subroutine test_closest_point_yld2004_shear_xy


subroutine test_closest_point_yld2004_shear_yz(passed)
    !! Pure yz shear isolates library c'55/c''55, distinct from xy and xz.
    ! Barlat et al. (2005), Table 2, Eqs. 14-15; AA2090-T3, a = 8.
    ! Library (xy, yz, xz) = table (c66, c44, c55); tensorial shears.
    ! Independent associated-flow reference, mpmath (50 digits), outside this solver.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_yld2004, only : Yld2004
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Stress and plastic strain, relative, at solver tolerance 1e-5. Measured: 2.69e-8 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 9.58e-10 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! Tangent vs centered differences of solves stopped at solver tolerance (a = 8).
    ! Finite differences need 1e-6 even with exact Hessian. Measured: 2.61e-7 or less.
    real(real64), parameter :: TOL_TAN = 1D-6
    type(Yld2004) :: yld
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: strain_pf
    integer :: status, iters

    passed = .true.
    call yld%init(cp12=-0.069888D0, cp13=0.936408D0, cp21=0.079143D0, cp23=1.003060D0, &
                  cp31=0.524741D0, cp32=1.363180D0, &
                  cp44=0.954322D0, cp55=1.023770D0, cp66=1.069060D0, &
                  cpp12=0.981171D0, cpp13=0.476741D0, cpp21=0.575316D0, cpp23=0.866827D0, &
                  cpp31=1.145010D0, cpp32=-0.079294D0, &
                  cpp44=1.404620D0, cpp55=1.051660D0, cpp66=1.147100D0, a=8D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=3.8911666904622642D+01, xz=0D0)
    call expected_strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=4.7625982239354998D-02, &
                                xz=0D0)
    call strain%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=9.8211149215364432D-02, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=yld)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    stress = history%state_np1%stress
    strain_p = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf
    call check_true(passed, "Converged status", status == STATUS_CONVERGED)
    call check_approx(passed, "Stress", stress, expected_stress, tol=TOL_STRESS)
    call check_approx(passed, "Equivalent plastic strain", strain_pf, 0.05D0, tol=TOL_PF)
    call check_approx(passed, "Plastic strain", strain_p, expected_strain_p, tol=TOL_STRESS)
    ! Newton iterations from the virgin state. Measured: 6 (Release, Debug, both Hessians).
    call check_true(passed, "Iterations <= 6", iters <= 6)
    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, tangent=tangent)
    call check_approx(passed, "Consistent tangent", tangent, numerical_tangent, tol=TOL_TAN)
end subroutine test_closest_point_yld2004_shear_yz


subroutine test_closest_point_yld2004_shear_xz(passed)
    !! Pure xz shear isolates library c'66/c''66, distinct from xy and yz.
    ! Barlat et al. (2005), Table 2, Eqs. 14-15; AA2090-T3, a = 8.
    ! Library (xy, yz, xz) = table (c66, c44, c55); tensorial shears.
    ! Independent associated-flow reference, mpmath (50 digits), outside this solver.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_yld2004, only : Yld2004
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    implicit none
    logical, intent(out) :: passed
    ! Stress and plastic strain, relative, at solver tolerance 1e-5. Measured: 2.69e-8 or less.
    real(real64), parameter :: TOL_STRESS = 1D-5
    ! Absolute equivalent plastic strain error (scalar scale = 1). Measured: 9.58e-10 or less.
    real(real64), parameter :: TOL_PF = 1D-7
    ! Tangent vs centered differences of solves stopped at solver tolerance (a = 8).
    ! Finite differences need 1e-6 even with exact Hessian. Measured: 2.61e-7 or less.
    real(real64), parameter :: TOL_TAN = 1D-6
    type(Yld2004) :: yld
    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: sw
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: strain_pf
    integer :: status, iters

    passed = .true.
    call yld%init(cp12=-0.069888D0, cp13=0.936408D0, cp21=0.079143D0, cp23=1.003060D0, &
                  cp31=0.524741D0, cp32=1.363180D0, &
                  cp44=0.954322D0, cp55=1.023770D0, cp66=1.069060D0, &
                  cpp12=0.981171D0, cpp13=0.476741D0, cpp21=0.575316D0, cpp23=0.866827D0, &
                  cpp31=1.145010D0, cpp32=-0.079294D0, &
                  cpp44=1.404620D0, cpp55=1.051660D0, cpp66=1.147100D0, a=8D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call expected_stress%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=3.6439657613105283D+01)
    call expected_strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, &
                                xz=5.0856854270682338D-02)
    call strain%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=9.8228409167719216D-02)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call history%init(strain_p=strain_p, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=yld)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    stress = history%state_np1%stress
    strain_p = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf
    call check_true(passed, "Converged status", status == STATUS_CONVERGED)
    call check_approx(passed, "Stress", stress, expected_stress, tol=TOL_STRESS)
    call check_approx(passed, "Equivalent plastic strain", strain_pf, 0.05D0, tol=TOL_PF)
    call check_approx(passed, "Plastic strain", strain_p, expected_strain_p, tol=TOL_STRESS)
    ! Newton iterations from the virgin state. Measured: 6 (Release, Debug, both Hessians).
    call check_true(passed, "Iterations <= 6", iters <= 6)
    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, tangent=tangent)
    call check_approx(passed, "Consistent tangent", tangent, numerical_tangent, tol=TOL_TAN)
end subroutine test_closest_point_yld2004_shear_xz
