program test_muscle_visco_kelvin
    implicit none
    logical :: passed

    call test_constant_rate(passed)
    if (.not. passed) STOP 1

    call test_step_and_hold(passed)
    if (.not. passed) STOP 2

    call test_time_step_limit(passed)
    if (.not. passed) STOP 3

    call test_frozen_stress_tangent_2d(passed)
    if (.not. passed) STOP 4

    call test_history(passed)
    if (.not. passed) STOP 5

    call test_closest_point_coupling(passed)
    if (.not. passed) STOP 6

    print*, "Passed!", passed
end program test_muscle_visco_kelvin

subroutine test_constant_rate(passed)
    ! Constant strain rate D: sigma(t) = t (K tr(D) I + 2 G dev(D)) + zeta tr(D) I + 2 mu_v dev(D),
    ! exact for any dt.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_kelvin, only : Kelvin_voigt
    implicit none
    logical, intent(out) :: passed

    type(Kelvin_voigt) :: model
    type(Viscoelastic_history) :: hist
    type(Elasticity_frozen) :: elas
    type(Visco_coefficients) :: coefs
    type(ten_3D2Osym) :: rate, strain, stress, expected
    real(real64), parameter :: K = 2000D0, G = 800D0, ZETA = 300D0, MU_V = 150D0
    real(real64), parameter :: DTS(5) = [0.03D0, 0.2D0, 0.5D0, 1.5D0, 3D0]
    real(real64) :: t, trd
    integer :: i
    character(len=44) :: label

    passed = .true.
    call model%set_parameters(bulk=K, shear=G, bulk_viscosity=ZETA, shear_viscosity=MU_V)
    call model%init_history(hist)

    ! Measured: 0
    coefs = model%coefficients(0.25D0)
    call check_approx(passed, "coefficients: K_alg = K + zeta/dt", coefs%k_alg, K + ZETA/0.25D0, tol=1D-14)
    call check_approx(passed, "coefficients: G_alg = G + mu_v/dt", coefs%g_alg, G + MU_V/0.25D0, tol=1D-14)

    call rate%init(xx=2D-3, yy=-1D-3, zz=5D-4, xy=1D-3, yz=3D-4, xz=-2D-4)
    trd = rate .ddot. iden_2O()

    ! Measured: 1.2e-16 or less
    t = 0D0
    do i = 1, size(DTS)
        t = t + DTS(i)
        strain = t*rate
        coefs = model%coefficients(DTS(i))
        elas = model%freeze(coefs, hist)
        stress = elas%stress(strain)
        call model%update(coefs, strain, hist)
        call hist%commit()

        expected = t*(K*trd*(iden_2O()) + 2D0*G*(.dev. rate)) + (ZETA*trd*(iden_2O()) + 2D0*MU_V*(.dev. rate))
        write(label, "(A, ES9.2)") "constant rate: stress at t =", t
        call check_approx(passed, trim(label), stress, expected, tol=1D-12)
    end do
end subroutine test_constant_rate

subroutine test_step_and_hold(passed)
    ! Step strain in a very short step, then held: with zero strain rate the stress is C:eps
    ! with any dt. A wrong sign or factor in sigma_hist passes the constant rate case, not this one.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_kelvin, only : Kelvin_voigt
    implicit none
    logical, intent(out) :: passed

    type(Kelvin_voigt) :: model
    type(Viscoelastic_history) :: hist
    type(Elasticity_frozen) :: elas
    type(ten_3D2Osym) :: strain, stress, expected
    real(real64), parameter :: K = 2000D0, G = 800D0, ZETA = 300D0, MU_V = 150D0
    real(real64), parameter :: DT_STEP = 1D-6
    real(real64), parameter :: DTS(3) = [0.5D0, 1D-6, 3D4]
    real(real64) :: tr
    integer :: i
    character(len=44) :: label

    passed = .true.
    call model%set_parameters(bulk=K, shear=G, bulk_viscosity=ZETA, shear_viscosity=MU_V)
    call model%init_history(hist)

    call strain%init(xx=1D-3, yy=-2D-4, zz=3D-4, xy=5D-4, yz=-1D-4, xz=2D-4)
    tr = strain .ddot. iden_2O()

    ! Step from eps_n = 0: elastic plus viscous response. Measured: 1.1e-16
    elas = model%freeze_dt(DT_STEP, hist)
    stress = elas%stress(strain)
    expected = (K*tr + ZETA*tr/DT_STEP)*(iden_2O()) + 2D0*(G + MU_V/DT_STEP)*(.dev. strain)
    call check_approx(passed, "step: elastic plus viscous response", stress, expected, tol=1D-12)
    call model%update_dt(DT_STEP, strain, hist)
    call hist%commit()

    ! Measured: 2.8e-16 (dt = 0.5 and 3e4), 1.5e-11 (dt = 1e-6: cancellation of the zeta/dt terms)
    do i = 1, size(DTS)
        elas = model%freeze_dt(DTS(i), hist)
        stress = elas%stress(strain)
        expected = K*tr*(iden_2O()) + 2D0*G*(.dev. strain)
        write(label, "(A, ES9.2)") "hold: stress with dt =", DTS(i)
        call check_approx(passed, trim(label), stress, expected, tol=1D-9)
        call model%update_dt(DTS(i), strain, hist)
        call hist%commit()
    end do
end subroutine test_step_and_hold

subroutine test_time_step_limit(passed)
    ! dt -> inf from a non-zero state, and zero viscosity at any dt: both give Elasticity_linear
    ! with the same K and G.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_kelvin, only : Kelvin_voigt
    implicit none
    logical, intent(out) :: passed

    type(Kelvin_voigt) :: model, elastic
    type(Elasticity_linear) :: lin
    type(Viscoelastic_history) :: hist, hist0
    type(Elasticity_frozen) :: elas
    type(Visco_coefficients) :: coefs
    type(ten_3D2Osym) :: strain
    real(real64), parameter :: K = 2000D0, G = 800D0, ZETA = 300D0, MU_V = 150D0
    real(real64), parameter :: DT_INF = 1D14
    real(real64) :: young, poisson

    passed = .true.
    call model%set_parameters(bulk=K, shear=G, bulk_viscosity=ZETA, shear_viscosity=MU_V)
    call model%init_history(hist)
    call strain%init(xx=1D-3, yy=-2D-4, zz=3D-4, xy=5D-4, yz=-1D-4, xz=2D-4)
    call model%update_dt(0.5D0, strain, hist)
    call hist%commit()

    young = 9D0*K*G/(3D0*K + G)
    poisson = (3D0*K - 2D0*G)/(2D0*(3D0*K + G))
    call lin%set_parameters(young=young, poisson=poisson)

    ! Measured: 1.9e-15 or less
    coefs = model%coefficients(DT_INF)
    call check_approx(passed, "dt -> inf: K_alg = K", coefs%k_alg, K, tol=1D-12)
    call check_approx(passed, "dt -> inf: G_alg = G", coefs%g_alg, G, tol=1D-12)

    ! Measured: 1.5e-15 or less
    call strain%init(xx=7D-4, yy=-1D-4, zz=2D-4, xy=2D-4, yz=1D-4, xz=-3D-4)
    elas = model%freeze(coefs, hist)
    call check_approx(passed, "dt -> inf: stress = linear elastic", elas%stress(strain), lin%stress(strain), tol=1D-12)
    call check_approx(passed, "dt -> inf: tangent = linear elastic", elas%dstress_dstrain(strain), &
                      lin%dstress_dstrain(strain), tol=1D-12)

    ! Measured: 1.9e-16 or less
    call elastic%set_parameters(bulk=K, shear=G, bulk_viscosity=0D0, shear_viscosity=0D0)
    call elastic%init_history(hist0)
    coefs = elastic%coefficients(0.5D0)
    elas = elastic%freeze(coefs, hist0)
    call check_approx(passed, "no viscosity: stress = linear elastic", elas%stress(strain), lin%stress(strain), &
                      tol=1D-12)
    call check_approx(passed, "no viscosity: tangent = linear elastic", elas%dstress_dstrain(strain), &
                      lin%dstress_dstrain(strain), tol=1D-12)
end subroutine test_time_step_limit

subroutine test_frozen_stress_tangent_2d(passed)
    ! Frozen stress at eps_n+1 (law with eps_dot = de/dt) and at eps_n (elastic), analytical vs
    ! numerical tangent, deviatoric stress and plane strain; normals and three distinct shears.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_math_derivatives, only : derivative
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_kelvin, only : Kelvin_voigt
    implicit none
    logical, intent(out) :: passed

    type(Kelvin_voigt) :: model
    type(Viscoelastic_history) :: hist
    type(Elasticity_frozen) :: elas
    type(Visco_coefficients) :: coefs
    type(ten_3D2Osym) :: strain, strain_n, de, stress, expected
    type(ten_2D2Osym) :: strain_2d, stress_2d, expected_2d
    type(ten_3D4O2sym) :: tan_an, tan_num
    real(real64), parameter :: K = 2000D0, G = 800D0, ZETA = 300D0, MU_V = 150D0
    real(real64), parameter :: DT = 0.3D0
    real(real64) :: tr_n, tr_de

    passed = .true.
    call model%set_parameters(bulk=K, shear=G, bulk_viscosity=ZETA, shear_viscosity=MU_V)
    call model%init_history(hist)

    ! Two converged steps, the first one in plane strain
    call strain%init(xx=1D-3, yy=-4D-4, zz=0D0, xy=6D-4, yz=0D0, xz=0D0)
    call model%update_dt(0.2D0, strain, hist)
    call hist%commit()
    call strain%init(xx=4D-4, yy=2D-4, zz=0D0, xy=-3D-4, yz=0D0, xz=0D0)
    call model%update_dt(0.5D0, strain, hist)
    call hist%commit()

    strain_n = hist%state_n%strain_e
    call strain%init(xx=7D-4, yy=-1D-4, zz=2D-4, xy=2D-4, yz=1D-4, xz=-3D-4)
    coefs = model%coefficients(DT)
    elas = model%freeze(coefs, hist)

    ! Measured: 4.1e-17
    de = strain - strain_n
    tr_n = strain_n .ddot. iden_2O()
    tr_de = de .ddot. iden_2O()
    expected = K*(strain .ddot. iden_2O())*(iden_2O()) + 2D0*G*(.dev. strain) &
             + (ZETA*tr_de/DT)*(iden_2O()) + (2D0*MU_V/DT)*(.dev. de)
    call check_approx(passed, "frozen stress = C:eps_n+1 + eta:de/dt", elas%stress(strain), expected, tol=1D-12)

    ! Measured: 3.7e-17
    expected = K*tr_n*(iden_2O()) + 2D0*G*(.dev. strain_n)
    call check_approx(passed, "frozen law at eps_n = elastic stress of eps_n", elas%stress(strain_n), expected, &
                      tol=1D-12)

    ! Measured: 2.0e-16
    call check_approx(passed, "stress_dev = dev(stress)", elas%stress_dev(strain), .dev. elas%stress(strain), &
                      tol=1D-12)

    ! Measured: 1.2e-13 (finite differences)
    tan_an = elas%dstress_dstrain(strain)
    tan_num = derivative(wrapper, strain)
    call check_approx(passed, "frozen tangent: analytical vs numerical", tan_an, tan_num, tol=1D-7)

    ! Measured: 0
    stress = elas%stress(strain)
    call strain_2d%init(xx=strain%xx(), yy=strain%yy(), zz=strain%zz(), xy=strain%xy())
    stress_2d = elas%stress(strain_2d)
    call expected_2d%init(xx=stress%xx(), yy=stress%yy(), zz=stress%zz(), xy=stress%xy())
    call check_approx(passed, "plane strain stress = 3D components", stress_2d, expected_2d, tol=1D-12)

    contains
    pure function wrapper(x) result(res)
        ! Stress of the frozen elasticity, to differentiate it numerically
        type(ten_3D2Osym), intent(in) :: x
        type(ten_3D2Osym) :: res
        res = elas%stress(x)
    end function wrapper
end subroutine test_frozen_stress_tangent_2d

subroutine test_history(passed)
    ! commit / rollback, freeze between update and commit, and the FEA pack / unpack round trip
    ! (6 components: the elastic strain, no branch tensors).
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_kelvin, only : Kelvin_voigt
    implicit none
    logical, intent(out) :: passed

    type(Kelvin_voigt) :: model
    type(Viscoelastic_history) :: hist, hist2
    type(Elasticity_frozen) :: elas
    type(ten_3D2Osym) :: strain, other, expected
    real(real64) :: hsv(6)

    passed = .true.
    call model%set_parameters(bulk=2000D0, shear=800D0, bulk_viscosity=300D0, shear_viscosity=150D0)
    call model%init_history(hist)
    call check_true(passed, "init: zero branch tensors", hist%state_n%n_branches == 0)
    call check_true(passed, "n_branches = 0", model%n_branches() == 0)

    call strain%init(xx=1D-3, yy=-2D-4, zz=3D-4, xy=5D-4, yz=-1D-4, xz=2D-4)
    call model%update_dt(0.1D0, strain, hist)
    call hist%commit()
    call other%init(xx=-2D-3, yy=1D-4, zz=1D-4, xy=0D0, yz=2D-4, xz=0D0)
    call model%update_dt(0.1D0, other, hist)

    ! Before the rollback, freeze reads state_n, not the unconverged state_np1. Measured: 0
    elas = model%freeze_dt(0.1D0, hist)
    expected = -(300D0/0.1D0)*(strain .ddot. iden_2O())*iden_2O() - (2D0*150D0/0.1D0)*(.dev. strain)
    call check_approx(passed, "freeze after update: history stress of state_n", elas%history_stress(), expected, &
                      tol=1D-12)

    ! Measured: 0
    call hist%rollback()
    call check_approx(passed, "rollback: elastic strain of state_np1", hist%state_np1%strain_e, strain, tol=0D0)
    call check_approx(passed, "rollback: elastic strain of state_n", hist%state_n%strain_e, strain, tol=0D0)

    ! Measured: 0
    call hist%pack_to_fea(hsv)
    call model%init_history(hist2)
    call hist2%unpack_from_fea(hsv)
    call check_approx(passed, "pack/unpack: 6 FEA history components", hist2%state_n%strain_e, strain, tol=0D0)
end subroutine test_history

subroutine test_closest_point_coupling(passed)
    ! The frozen elasticity with a non-zero history stress in the rate-independent closest point:
    ! f = 0, consistent tangent, and the update at eps - eps_p gives the returned stress.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_kelvin, only : Kelvin_voigt
    implicit none
    logical, intent(out) :: passed

    type(Kelvin_voigt) :: model
    type(Viscoelastic_history) :: vhist
    type(Plastic_material_history) :: phist
    type(Elasticity_frozen) :: elas
    type(Visco_coefficients) :: coefs
    type(Closest_point) :: solver
    type(VonMises) :: vm
    type(Swift_hardening) :: sw
    type(ten_3D2Osym) :: strain, strain_e, strain_p0, stress_hist
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64), parameter :: K = 800D0, G = 400D0, ZETA = 50D0, MU_V = 20D0, DT = 0.05D0
    real(real64) :: f
    integer :: status, iters

    passed = .true.
    call model%set_parameters(bulk=K, shear=G, bulk_viscosity=ZETA, shear_viscosity=MU_V)
    call model%init_history(vhist)
    call strain%init(xx=2D-2, yy=-5D-3, zz=-5D-3, xy=4D-3, yz=2D-3, xz=1D-3)
    call model%update_dt(DT, strain, vhist)
    call vhist%commit()

    coefs = model%coefficients(DT)
    elas = model%freeze(coefs, vhist)
    stress_hist = elas%history_stress()
    call check_true(passed, "coupling: history stress is non-zero", stress_hist%norm() > 0D0)

    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call strain_p0%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call phist%init(strain_p=strain_p0, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=vm)

    call strain%init(xx=0.15D0, yy=-0.05D0, zz=-0.05D0, xy=1D-2, yz=-6D-3, xz=4D-3)
    call solver%solve(strain=strain, history=phist, status=status, iters=iters)
    call check_true(passed, "coupling: return mapping converged", status == STATUS_CONVERGED)
    call check_true(passed, "coupling: plastic step", phist%state_np1%strain_pf > 0D0)

    ! Measured: 2.6e-6 (solver residual, TOL = 1D-5)
    f = vm%stress_eq(phist%state_np1%stress) - sw%stress(phist%state_np1%strain_pf)
    call check_approx(passed, "coupling: f = 0", f, 0D0, tol=1D-5)

    ! Measured: 2.3e-8
    call solver%tangent(strain=strain, history=phist, tangent=tangent)
    call solver%tangent_numerical(strain=strain, history=phist, tangent=numerical_tangent)
    call check_approx(passed, "coupling: consistent tangent", tangent, numerical_tangent, tol=1D-6)

    ! Measured: 0
    strain_e = strain - phist%state_np1%strain_p
    call model%update(coefs, strain_e, vhist)
    call check_approx(passed, "coupling: frozen stress = returned stress", elas%stress(strain_e), &
                      phist%state_np1%stress, tol=1D-12)
    call check_approx(passed, "coupling: update stores eps_e", vhist%state_np1%strain_e, strain_e, tol=0D0)
end subroutine test_closest_point_coupling
