! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_solver_closest_point_apex
    implicit none
    logical :: passed

    call test_closest_point_dp_apex_linear(passed)
    if (.not. passed) STOP 1

    call test_closest_point_dp_apex_swift(passed)
    if (.not. passed) STOP 2

    call test_closest_point_dp_apex_threshold(passed)
    if (.not. passed) STOP 3

    call test_closest_point_no_apex(passed)
    if (.not. passed) STOP 4

    print*, "Passed!", passed
end program test_muscle_solver_closest_point_apex

subroutine check_apex_case(solver, strain, p, iters_expected, tol, history, passed)
    ! Solve and compare with the apex p*I (status, scalar Newton iterations and stress).
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    implicit none
    type(Closest_point), intent(in) :: solver
    type(ten_3D2Osym), intent(in) :: strain
    real(real64), intent(in) :: p, tol
    integer, intent(in) :: iters_expected
    type(Plastic_material_history), intent(out) :: history
    logical, intent(out) :: passed
    type(ten_3D2Osym) :: expected_stress
    type(iden_2O) :: I2O
    integer :: status, iters

    expected_stress = p * I2O
    call history%init()
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    passed = status == STATUS_CONVERGED .and. iters == iters_expected
    passed = passed .and. history%state_np1%stress%is_approx(expected_stress, tol=tol)
    if (.not. passed) print *, "FAIL: apex case", status, iters, history%state_np1%stress%vals
end subroutine check_apex_case

subroutine check_smooth_case(solver, strain, history, passed)
    ! Solve and compare bit by bit with the Newton of solve driven through the public iter
    ! (states that converge before the first relaxation, iteration 10): the smooth return is
    ! unchanged by the apex branch.
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED, STATUS_ITER_CONVERGED
    implicit none
    type(Closest_point), intent(in) :: solver
    type(ten_3D2Osym), intent(in) :: strain
    type(Plastic_material_history), intent(out) :: history
    logical, intent(out) :: passed
    type(Plastic_material_history) :: reference
    real(real64) :: dgamma, omega, ddgamma
    integer :: status, iters, iters_ref, iter_status

    call history%init()
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    call reference%init()
    reference%state_np1 = reference%state_n
    dgamma = 0D0
    omega  = 1D0
    do iters_ref = 1, 9
        call solver%iter(strain, reference, dgamma, omega, ddgamma, iter_status)
        if (iter_status == STATUS_ITER_CONVERGED) exit
    end do
    passed = status == STATUS_CONVERGED .and. iters == iters_ref
    passed = passed .and. history%state_np1%stress%is_approx(reference%state_np1%stress, tol=0D0)
    passed = passed .and. history%state_np1%strain_p%is_approx(reference%state_np1%strain_p, tol=0D0)
    passed = passed .and. abs(history%state_np1%strain_pf - reference%state_np1%strain_pf) <= 0D0
    if (.not. passed) print *, "FAIL: smooth case", status, iters, iters_ref
end subroutine check_smooth_case

subroutine test_closest_point_dp_apex_linear(passed)
    ! Drucker-Prager (K = 1, TENSION, beta = 16) with linear hardening: closed form of the apex
    ! return (Sysala et al., 2016, Eq. 3.27-3.28) and of its tangent (Eq. 3.32), also from a
    ! previous plastic state.
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_bilinear, only : Bilinear_hardening
    use muscle_yield_druckerprager, only : DruckerPrager, DP_HARDENING_TENSION
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    implicit none

    real(real64), parameter :: EPS = 1D-10, EPS_NUM = 1D-7
    real(real64), parameter :: YOUNG = 200000D0, POISSON = 0.3D0, Y0 = 100D0, H = 1000D0
    logical, intent(out) :: passed
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(DruckerPrager) :: dp
    type(Bilinear_hardening) :: bl
    type(Elasticity_linear) :: elas
    type(ten_3D2Osym) :: strain, expected_strain_p, strain_p_n, apex_stress
    type(ten_3D4O2sym) :: tangent, numerical_tangent, expected_tangent
    type(ten_3D4O3sym) :: tangent_3s
    type(iden_2O) :: I2O
    type(iden_4O3T) :: I4O3T
    real(real64) :: kb, d, tanbeta, dgamma, p, strain_pf_n
    integer :: k, status, iters

    tanbeta = tan(16D0 * acos(-1D0) / 180D0)
    d  = tanbeta / (1D0 + tanbeta / 3D0)
    kb = YOUNG / (3D0 * (1D0 - 2D0 * POISSON))
    tangent_3s       = (kb * H / (kb * d * d + H)) * I4O3T
    expected_tangent = tangent_3s

    call elas%set_parameters(young=YOUNG, poisson=POISSON)
    bl = Bilinear_hardening(y0=Y0, K=H)
    call dp%init(beta_deg=16D0, K=1D0, hardening_mode=DP_HARDENING_TENSION)
    call solver%init(elasticity=elas, hardening=bl, yield=dp)

    do k = 1, 3
        call strain_p_n%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
        strain_pf_n = 0D0
        if (k == 1) then
            ! Triaxial tension with all shears: trial beyond the tip
            call strain%init(xx=2.0D-3, yy=1.6D-3, zz=1.4D-3, xy=1.0D-4, yz=-5.0D-5, xz=8.0D-5)
        else if (k == 2) then
            ! Pure volumetric trial: zero deviatoric flow, no Lode angle
            call strain%init(xx=1.0D-3, yy=1.0D-3, zz=1.0D-3, xy=0D0, yz=0D0, xz=0D0)
        else
            ! Previous plastic state with normals and shears (trace d*strain_pf_n, as an apex
            ! history), plus a trial elastic strain with all shears
            strain_pf_n = 0.027D0
            call strain_p_n%init(xx=2.0D-3, yy=-3.0D-4, zz=7.0D-4, xy=2.0D-4, yz=-4.0D-4, xz=1.0D-4)
            strain_p_n = strain_p_n + ((d * strain_pf_n - (strain_p_n .ddot. I2O)) / 3D0) * I2O
            call strain%init(xx=3.0D-3, yy=2.7D-3, zz=2.3D-3, xy=1.5D-4, yz=-1.0D-4, xz=1.2D-4)
            strain = strain_p_n + strain
        end if
        dgamma = (d * kb * ((strain - strain_p_n) .ddot. I2O) - Y0 - H * strain_pf_n) / (kb * d * d + H)
        p      = kb * ((strain - strain_p_n) .ddot. I2O) - kb * d * dgamma
        expected_strain_p = strain - (p / (3D0 * kb)) * I2O

        ! Linear hardening: the scalar Newton ends in one iteration
        call history%init(strain_p=strain_p_n, strain_pf=strain_pf_n)
        call solver%solve(strain=strain, history=history, status=status, iters=iters)
        passed = status == STATUS_CONVERGED .and. iters == 1
        apex_stress = p * I2O
        passed = passed .and. history%state_np1%stress%is_approx(apex_stress, tol=EPS)
        call solver%tangent(strain=strain, history=history, tangent=tangent)
        call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
        passed = passed .and. history%state_np1%strain_p%is_approx(expected_strain_p, tol=EPS)
        passed = passed .and. abs(history%state_np1%strain_pf - strain_pf_n - dgamma) < EPS * dgamma
        passed = passed .and. tangent%is_approx(expected_tangent, tol=EPS)
        passed = passed .and. tangent%is_approx(numerical_tangent, tol=EPS_NUM)
        if (.not. passed) then
            print *, "FAIL: apex return, linear hardening, state", k
            return
        end if
    end do
end subroutine test_closest_point_dp_apex_linear

subroutine test_closest_point_dp_apex_swift(passed)
    ! Drucker-Prager (K = 0.85, TENSION, beta = 16) with Swift hardening: Lode-dependent cone
    ! and nonlinear scalar equation. References from mpmath (40 digits) of Sysala et al. (2016),
    ! Eq. 3.27-3.28 and 3.32.
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_druckerprager, only : DruckerPrager, DP_HARDENING_TENSION
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point
    implicit none

    real(real64), parameter :: EPS = 1D-10, EPS_NUM = 1D-7
    real(real64), parameter :: P_REF = 444.05540549227005948D0, DGAMMA_REF = 0.0099178100105184492369D0
    real(real64), parameter :: CEP_REF = 194.82815340995678924D0
    logical, intent(out) :: passed
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(DruckerPrager) :: dp
    type(Swift_hardening) :: sw
    type(Elasticity_linear) :: elas
    type(ten_3D2Osym) :: strain, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent, expected_tangent
    type(ten_3D4O3sym) :: tangent_3s
    type(iden_4O3T) :: I4O3T

    call elas%set_parameters(young=200000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D0)
    call dp%init(beta_deg=16D0, K=0.85D0, hardening_mode=DP_HARDENING_TENSION)
    call solver%init(elasticity=elas, hardening=sw, yield=dp)

    ! Triaxial tension with all shears, beyond the tip
    call strain%init(xx=2.2D-3, yy=1.5D-3, zz=1.2D-3, xy=2.0D-4, yz=1.0D-4, xz=-1.5D-4)
    call expected_strain_p%init(xx=0.001311889189015459881D0, yy=0.00061188918901545988104D0, &
                                zz=0.00031188918901545988104D0, xy=2.0D-4, yz=1.0D-4, xz=-1.5D-4)
    tangent_3s       = CEP_REF * I4O3T
    expected_tangent = tangent_3s

    ! Two Newton iterations from dgamma = 0 reach |g| < 1D-5
    call check_apex_case(solver, strain, P_REF, 2, EPS, history, passed)
    call solver%tangent(strain=strain, history=history, tangent=tangent)
    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    passed = passed .and. history%state_np1%strain_p%is_approx(expected_strain_p, tol=EPS)
    passed = passed .and. abs(history%state_np1%strain_pf - DGAMMA_REF) < EPS * DGAMMA_REF
    passed = passed .and. tangent%is_approx(expected_tangent, tol=EPS)
    passed = passed .and. tangent%is_approx(numerical_tangent, tol=EPS_NUM)
    if (.not. passed) print *, "FAIL: apex return, Swift hardening:", tangent%vals(1, 1:3)
end subroutine test_closest_point_dp_apex_swift

subroutine test_closest_point_dp_apex_threshold(passed)
    ! Trial strains (v/3) I + lambda e0 just inside (0.999 lambda*) and just outside
    ! (1.001 lambda*) the apex cone, with apex_dev_gauge(lambda* e0) = dgamma at the apex.
    ! Inside: apex. Outside: the unchanged smooth Newton, compared with the mpmath smooth return
    ! (q = 0.099 for K = 1, 0.072 for K = 0.85), which is 7D-4 below the apex pressure.
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_bilinear, only : Bilinear_hardening
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_druckerprager, only : DruckerPrager, DP_HARDENING_TENSION
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point
    implicit none

    real(real64), parameter :: EPS = 1D-10, EPS_SMOOTH = 1D-6
    logical, intent(out) :: passed
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(DruckerPrager) :: dp
    type(Elasticity_linear) :: elas
    type(ten_3D2Osym) :: strain, expected_stress
    logical :: same

    call elas%set_parameters(young=200000D0, poisson=0.3D0)

    ! K = 1, linear hardening, deviatoric direction with Lode angle 0.6
    call dp%init(beta_deg=16D0, K=1D0, hardening_mode=DP_HARDENING_TENSION)
    call solver%init(elasticity=elas, hardening=Bilinear_hardening(y0=100D0, K=1000D0), yield=dp)
    call strain%init(xx=0.00716986109250104D0, yy=-0.00038995369750034667D0, zz=-0.0022799073950006933D0, &
                     xy=0.00283493054625052D0, yz=-0.0018899536975003467D0, xz=0.00094497684875017334D0)
    call check_apex_case(solver, strain, 411.70554313931258012D0, 1, EPS, history, passed)
    if (.not. passed) return
    call strain%init(xx=0.0071812121657593004D0, yy=-0.00039373738858643345D0, zz=-0.0022874747771728669D0, &
                     xy=0.0028406060828796502D0, yz=-0.0018937373885864335D0, xz=0.00094686869429321673D0)
    call expected_stress%init(xx=411.44074323731654221D0, yy=411.37011464205624344D0, &
                              zz=411.35245749324116875D0, xy=0.026485723222612038731D0, &
                              yz=-0.017657148815074692798D0, xz=0.0088285744075373462127D0)
    call check_smooth_case(solver, strain, history, same)
    passed = same .and. history%state_np1%stress%is_approx(expected_stress, tol=EPS)
    if (.not. passed) return

    ! K = 0.85, Swift, on the tension meridian (Lode 0): off the meridians the smooth Newton is
    ! not reliable this close to the tip. It stops at TOL = 1D-5 (error 4.3D-7 here)
    call dp%init(beta_deg=16D0, K=0.85D0, hardening_mode=DP_HARDENING_TENSION)
    call solver%init(elasticity=elas, hardening=Swift_hardening(k=100D0, n=0.1D0, e0=1D0), yield=dp)
    call strain%init(xx=-0.00081026270417024298D0, yy=0.00076137966791118435D0, zz=0.0045488830362571772D0, &
                     xy=-0.0020960119722120298D0, yz=0.0045388239817337429D0, xz=-0.0031456647209452598D0)
    call check_apex_case(solver, strain, 443.97741268870125493D0, 2, EPS, history, passed)
    if (.not. passed) return
    call strain%init(xx=-0.00081488785472914237D0, yy=0.00075990094852762315D0, zz=0.004554986906199634D0, &
                     xy=-0.0021002081923766184D0, yz=0.0045479107164319086D0, xz=-0.0031519623480142193D0)
    call expected_stress%init(xx=443.66934078098584742D0, yy=443.67930905122067293D0, &
                              zz=443.70333160098221622D0, xy=-0.013294127293479743085D0, &
                              yz=0.028787862176289878213D0, xz=-0.019951635666858500001D0)
    call check_smooth_case(solver, strain, history, same)
    passed = same .and. history%state_np1%stress%is_approx(expected_stress, tol=EPS_SMOOTH)
    if (.not. passed) print *, "FAIL: just outside the apex cone, K = 0.85:", history%state_np1%stress%vals
end subroutine test_closest_point_dp_apex_threshold

subroutine test_closest_point_no_apex(passed)
    ! States where the apex return must not act: a criterion without apex (VonMises, slope 0),
    ! DP trials with g(0) <= 0 in compression and with a positive trace (the apex root would be
    ! negative, and Swift with e0 = 1D-4 has no value below -1D-4), and a DP apex trial with
    ! softening (Sysala et al., 2016, assume non-decreasing hardening).
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use, intrinsic :: ieee_exceptions
    use muscle_hard_swift, only : Swift_hardening
    use muscle_hard_bilinear, only : Bilinear_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_yield_druckerprager, only : DruckerPrager, DP_HARDENING_TENSION
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point
    implicit none

    logical, intent(out) :: passed
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(VonMises) :: vm
    type(DruckerPrager) :: dp
    type(Elasticity_linear) :: elas
    type(ten_3D2Osym) :: strain, apex_stress
    type(iden_2O) :: I2O
    real(real64) :: kb, d, p
    logical :: invalid

    call elas%set_parameters(young=200000D0, poisson=0.3D0)

    ! VonMises, coupled plastic trial with mean strain 3D-3
    call solver%init(elasticity=elas, hardening=Swift_hardening(k=100D0, n=0.1D0, e0=1D-4), yield=vm)
    call strain%init(xx=5.0D-3, yy=2.6D-3, zz=1.4D-3, xy=6.0D-4, yz=-3.0D-4, xz=4.0D-4)
    call check_smooth_case(solver, strain, history, passed)
    if (.not. passed) return

    ! DP compression with all shears: no invalid operation in the apex branch
    call dp%init(beta_deg=16D0, K=0.85D0, hardening_mode=DP_HARDENING_TENSION)
    call solver%init(elasticity=elas, hardening=Swift_hardening(k=100D0, n=0.1D0, e0=1D-4), yield=dp)
    call strain%init(xx=-3.0D-3, yy=-1.2D-3, zz=-0.8D-3, xy=4.0D-4, yz=-2.0D-4, xz=3.0D-4)
    call ieee_set_flag(ieee_invalid, .false.)
    call check_smooth_case(solver, strain, history, passed)
    call ieee_get_flag(ieee_invalid, invalid)
    passed = passed .and. .not. invalid
    if (.not. passed) return

    ! DP, positive trace, g(0) < 0
    call strain%init(xx=5.0D-3, yy=-2.0D-3, zz=-2.99D-3, xy=4.0D-4, yz=-2.0D-4, xz=3.0D-4)
    call ieee_set_flag(ieee_invalid, .false.)
    call check_smooth_case(solver, strain, history, passed)
    call ieee_get_flag(ieee_invalid, invalid)
    passed = passed .and. .not. invalid
    if (.not. passed) return

    ! DP, K = 1, softening (sigma_y' = -500): the apex root exists but is not used
    call dp%init(beta_deg=16D0, K=1D0, hardening_mode=DP_HARDENING_TENSION)
    call solver%init(elasticity=elas, hardening=Bilinear_hardening(y0=100D0, K=-500D0), yield=dp)
    call strain%init(xx=2.0D-3, yy=1.6D-3, zz=1.4D-3, xy=1.0D-4, yz=-5.0D-5, xz=8.0D-5)
    kb = 200000D0 / (3D0 * (1D0 - 2D0 * 0.3D0))
    d  = tan(16D0 * acos(-1D0) / 180D0)
    d  = d / (1D0 + d / 3D0)
    p  = kb * (strain .ddot. I2O) - kb * d * (d * kb * (strain .ddot. I2O) - 100D0) / (kb * d * d - 500D0)
    apex_stress = p * I2O
    call history%init()
    call solver%solve(strain=strain, history=history)
    passed = .not. history%state_np1%stress%is_approx(apex_stress, tol=1D-6)
    if (.not. passed) print *, "FAIL: softening, the apex return was used"
end subroutine test_closest_point_no_apex
