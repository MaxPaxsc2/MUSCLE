program test_muscle_visco_maxwell
    implicit none
    logical :: passed

    call test_relaxation(passed)
    if (.not. passed) STOP 1

    call test_constant_rate(passed)
    if (.not. passed) STOP 2

    call test_time_step_limits(passed)
    if (.not. passed) STOP 3

    call test_freeze_update_consistency(passed)
    if (.not. passed) STOP 4

    call test_history(passed)
    if (.not. passed) STOP 5

    call test_closest_point_coupling(passed)
    if (.not. passed) STOP 6

    print*, "Passed!", passed
end program test_muscle_visco_maxwell

subroutine test_relaxation(passed)
    !! Step strain and hold: sigma(t) = K tr(eps) I + 2 G(t) dev(eps), exact for any dt
    !! (only the first step, which applies the strain, has to be short).
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: hist
    type(Elasticity_frozen) :: elas
    type(ten_3D2Osym) :: strain, stress, expected, vol
    real(real64), parameter :: K = 2000D0, DT1 = 1D-12
    real(real64), parameter :: DTS(5) = [0.05D0, 0.3D0, 0.7D0, 2D0, 5D0]
    real(real64) :: t, tr
    integer :: i
    character(len=40) :: label

    passed = .true.
    call model%set_parameters(bulk=K, shear_inf=100D0, shear=[200D0, 300D0], tau=[0.1D0, 1D0])
    call model%init_history(hist)

    call strain%init(xx=1D-3, yy=-2D-4, zz=3D-4, xy=5D-4, yz=-1D-4, xz=2D-4)
    tr = strain%xx() + strain%yy() + strain%zz()
    call vol%init(xx=K*tr, yy=K*tr, zz=K*tr, xy=0D0, yz=0D0, xz=0D0)

    ! Apply the strain in a very short step
    t = DT1
    elas = model%freeze_dt(DT1, hist)
    stress = elas%stress(strain)
    call model%update_dt(DT1, strain, hist)
    call hist%commit()
    expected = vol + 2D0*model%shear_instantaneous()*(.dev. strain)
    call check_approx(passed, "relaxation: instantaneous response", stress, expected, tol=1D-9)

    ! Hold the strain with coarse steps
    do i = 1, size(DTS)
        t = t + DTS(i)
        elas = model%freeze_dt(DTS(i), hist)
        stress = elas%stress(strain)
        call model%update_dt(DTS(i), strain, hist)
        call hist%commit()
        expected = vol + 2D0*model%relaxation_shear(t)*(.dev. strain)
        write(label, "(A, ES9.2)") "relaxation: stress at t =", t
        call check_approx(passed, trim(label), stress, expected, tol=1D-9)
    end do
end subroutine test_relaxation

subroutine test_constant_rate(passed)
    !! Constant strain rate D: the update is exact for any dt,
    !! sigma(t) = K tr(D) t I + 2 dev(D) [G_inf t + sum G_k tau_k (1 - exp(-t/tau_k))].
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: hist
    type(Elasticity_frozen) :: elas
    type(Visco_coefficients) :: coefs
    type(ten_3D2Osym) :: rate, strain, stress, expected, vol
    real(real64), parameter :: K = 2000D0, G_INF = 100D0
    real(real64), parameter :: G(2) = [200D0, 300D0], TAU(2) = [0.1D0, 1D0]
    real(real64), parameter :: DTS(5) = [0.03D0, 0.2D0, 0.5D0, 1.5D0, 3D0]
    real(real64) :: t, tr, gt
    integer :: i
    character(len=40) :: label

    passed = .true.
    call model%set_parameters(bulk=K, shear_inf=G_INF, shear=G, tau=TAU)
    call model%init_history(hist)

    call rate%init(xx=2D-3, yy=-1D-3, zz=5D-4, xy=1D-3, yz=3D-4, xz=-2D-4)
    tr = rate%xx() + rate%yy() + rate%zz()

    t = 0D0
    do i = 1, size(DTS)
        t = t + DTS(i)
        strain = t*rate
        coefs = model%coefficients(DTS(i))
        elas = model%freeze(coefs, hist)
        stress = elas%stress(strain)
        call model%update(coefs, strain, hist)
        call hist%commit()

        gt = G_INF*t + sum(G*TAU*(1D0 - exp(-t/TAU)))
        call vol%init(xx=K*tr*t, yy=K*tr*t, zz=K*tr*t, xy=0D0, yz=0D0, xz=0D0)
        expected = vol + 2D0*gt*(.dev. rate)
        write(label, "(A, ES9.2)") "constant rate: stress at t =", t
        call check_approx(passed, trim(label), stress, expected, tol=1D-10)
    end do
end subroutine test_constant_rate

subroutine test_time_step_limits(passed)
    !! dt -> 0 gives the instantaneous (linear elastic) response, dt -> inf the long-term one.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: hist
    type(Elasticity_frozen) :: elas
    type(Elasticity_linear) :: lin
    type(Visco_coefficients) :: coefs
    type(ten_3D2Osym) :: strain
    type(ten_3D4O3sym) :: tan_frozen, tan_lin
    real(real64), parameter :: E = 1000D0, NU = 0.3D0
    real(real64) :: g0

    passed = .true.
    call model%set_parameters_normalized(young=E, poisson=NU, g=[0.2D0, 0.3D0], tau=[0.1D0, 1D0])
    call model%init_history(hist)
    call lin%set_parameters(young=E, poisson=NU)
    g0 = 0.5D0*E/(1D0 + NU)

    call check_approx(passed, "normalized: G_0", model%shear_instantaneous(), g0, tol=1D-12)
    call check_approx(passed, "normalized: G(inf)", model%relaxation_shear(1D12), 0.5D0*g0, tol=1D-12)

    coefs = model%coefficients(0D0)
    call check_approx(passed, "dt = 0: G_alg = G_0", coefs%g_alg, g0, tol=1D-12)
    call check_approx(passed, "dt = 0: K_alg", coefs%k_alg, E/(3D0*(1D0 - 2D0*NU)), tol=1D-12)

    ! dt = 0 without history must be exactly the linear elastic law
    call strain%init(xx=1D-3, yy=-2D-4, zz=3D-4, xy=5D-4, yz=-1D-4, xz=2D-4)
    elas = model%freeze(coefs, hist)
    call check_approx(passed, "dt = 0: stress = linear elastic", elas%stress(strain), &
                      lin%stress(strain), tol=1D-12)
    tan_frozen = elas%dstress_dstrain(strain)
    tan_lin = lin%dstress_dstrain(strain)
    call check_approx(passed, "dt = 0: tangent = linear elastic", tan_frozen, tan_lin, tol=1D-12)

    coefs = model%coefficients(1D12)
    call check_approx(passed, "dt -> inf: G_alg = G_inf", coefs%g_alg, 0.5D0*g0, tol=1D-9)
end subroutine test_time_step_limits

subroutine test_freeze_update_consistency(passed)
    !! The frozen stress at eps_e_n+1 equals K tr(eps) I + 2 G_inf dev(eps) + sum s_k^(n+1)
    !! after the update. Also: analytical vs numerical tangent, deviatoric part and plane strain.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_math_derivatives
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: hist
    type(Elasticity_frozen) :: elas
    type(Visco_coefficients) :: coefs
    type(ten_3D2Osym) :: strain, branches, expected, vol
    type(ten_2D2Osym) :: strain_2d, stress_2d, expected_2d
    type(ten_3D4O2sym) :: tan_num, tan_an
    real(real64), parameter :: K = 2000D0, G_INF = 100D0
    real(real64) :: tr
    integer :: i

    passed = .true.
    call model%set_parameters(bulk=K, shear_inf=G_INF, shear=[200D0, 300D0, 50D0], &
                              tau=[0.1D0, 1D0, 10D0])
    call model%init_history(hist)

    ! Build a non-trivial (plane strain) history
    call strain%init(xx=1D-3, yy=-4D-4, zz=0D0, xy=6D-4, yz=0D0, xz=0D0)
    call model%update_dt(0.2D0, strain, hist)
    call hist%commit()
    call strain%init(xx=4D-4, yy=2D-4, zz=0D0, xy=-3D-4, yz=0D0, xz=0D0)
    call model%update_dt(0.5D0, strain, hist)
    call hist%commit()

    call strain%init(xx=7D-4, yy=-1D-4, zz=0D0, xy=2D-4, yz=0D0, xz=0D0)
    coefs = model%coefficients(0.3D0)
    elas = model%freeze(coefs, hist)
    call model%update(coefs, strain, hist)

    branches%vals = 0D0
    do i = 1, model%n_branches()
        branches%vals = branches%vals + hist%state_np1%h(:, i)
    end do
    tr = strain%xx() + strain%yy() + strain%zz()
    call vol%init(xx=K*tr, yy=K*tr, zz=K*tr, xy=0D0, yz=0D0, xz=0D0)
    expected = vol + 2D0*G_INF*(.dev. strain) + branches
    call check_approx(passed, "freeze(eps_n+1) = K tr I + 2 G_inf e + sum s_k^(n+1)", &
                      elas%stress(strain), expected, tol=1D-12)

    call check_approx(passed, "stress_dev = dev(stress)", elas%stress_dev(strain), &
                      .dev. elas%stress(strain), tol=1D-12)

    tan_an = elas%dstress_dstrain(strain)
    tan_num = derivative(wrapper, strain)
    call check_approx(passed, "frozen tangent: analytical vs numerical", tan_an, tan_num, tol=1D-7)

    call strain_2d%init(xx=strain%xx(), yy=strain%yy(), zz=strain%zz(), xy=strain%xy())
    stress_2d = elas%stress(strain_2d)
    call expected_2d%init(xx=expected%xx(), yy=expected%yy(), zz=expected%zz(), xy=expected%xy())
    call check_approx(passed, "plane strain stress = 3D components", stress_2d, expected_2d, tol=1D-12)

    contains
    pure function wrapper(x) result(res)
        type(ten_3D2Osym), intent(in) :: x
        type(ten_3D2Osym) :: res
        res = elas%stress(x)
    end function wrapper
end subroutine test_freeze_update_consistency

subroutine test_history(passed)
    !! commit / rollback and FEA pack / unpack round trip.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: hist, hist2
    type(ten_3D2Osym) :: strain
    real(real64) :: hsv(18), h_n(6, 2)

    passed = .true.
    call model%set_parameters(bulk=2000D0, shear_inf=100D0, shear=[200D0, 300D0], tau=[0.1D0, 1D0])
    call model%init_history(hist)
    call check_true(passed, "init: n_branches", hist%state_n%n_branches == 2)

    call strain%init(xx=1D-3, yy=-2D-4, zz=3D-4, xy=5D-4, yz=-1D-4, xz=2D-4)
    call model%update_dt(0.1D0, strain, hist)
    call hist%commit()
    h_n = hist%state_n%h(:, 1:2)

    ! A failed global iteration must not touch state_n
    call model%update_dt(0.1D0, 2D0*strain, hist)
    call hist%rollback()
    call check_approx(passed, "rollback: history tensors", hist%state_np1%h(:, 1:2), h_n, tol=0D0)
    call check_approx(passed, "rollback: elastic strain", hist%state_np1%strain_e, strain, tol=0D0)

    call hist%pack_to_fea(hsv)
    call model%init_history(hist2)
    call hist2%unpack_from_fea(hsv)
    call check_approx(passed, "pack/unpack: history tensors", hist2%state_n%h(:, 1:2), h_n, tol=0D0)
    call check_approx(passed, "pack/unpack: elastic strain", hist2%state_n%strain_e, strain, tol=0D0)
end subroutine test_history

subroutine test_closest_point_coupling(passed)
    !! The frozen elasticity (with a non-zero history stress) plugged into the rate-independent
    !! closest point: the return mapping converges, f = 0, its consistent tangent matches the
    !! numerical one, and the viscoelastic update at eps - eps_p reproduces the returned stress.
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
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: vhist
    type(Plastic_material_history) :: phist
    type(Elasticity_frozen) :: elas
    type(Visco_coefficients) :: coefs
    type(Closest_point) :: solver
    type(VonMises) :: vm
    type(Swift_hardening) :: sw
    type(ten_3D2Osym) :: strain, strain_e, strain_p0, stress_hist
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: f
    integer :: status, iters

    passed = .true.
    call model%set_parameters_normalized(young=1000D0, poisson=0.3D0, g=[0.3D0, 0.2D0], &
                                         tau=[0.1D0, 1D0])
    call model%init_history(vhist)

    ! Previous step leaves a non-zero viscous history
    call strain%init(xx=2D-2, yy=-5D-3, zz=-5D-3, xy=4D-3, yz=0D0, xz=1D-3)
    call model%update_dt(0.05D0, strain, vhist)
    call vhist%commit()

    coefs = model%coefficients(0.05D0)
    elas = model%freeze(coefs, vhist)
    stress_hist = elas%history_stress()
    call check_true(passed, "coupling: history stress is non-zero", stress_hist%norm() > 0D0)

    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call strain_p0%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call phist%init(strain_p=strain_p0, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=vm)

    call strain%init(xx=0.15D0, yy=-0.05D0, zz=-0.05D0, xy=1D-2, yz=0D0, xz=0D0)
    call solver%solve(strain=strain, history=phist, status=status, iters=iters)
    call check_true(passed, "coupling: return mapping converged", status == STATUS_CONVERGED)
    call check_true(passed, "coupling: plastic step", phist%state_np1%strain_pf > 0D0)

    f = vm%stress_eq(phist%state_np1%stress) - sw%stress(phist%state_np1%strain_pf)
    call check_approx(passed, "coupling: f = 0", f, 0D0, tol=1D-5)

    call solver%tangent(strain=strain, history=phist, tangent=tangent)
    call solver%tangent_numerical(strain=strain, history=phist, tangent=numerical_tangent)
    call check_approx(passed, "coupling: consistent tangent", tangent, numerical_tangent, tol=1D-6)

    ! Viscoelastic update with the converged elastic strain
    strain_e = strain - phist%state_np1%strain_p
    call model%update(coefs, strain_e, vhist)
    call check_approx(passed, "coupling: frozen stress = returned stress", elas%stress(strain_e), &
                      phist%state_np1%stress, tol=1D-12)
    call check_true(passed, "coupling: update stores eps_e", &
                    vhist%state_np1%strain_e%is_approx(strain_e, tol=0D0))
end subroutine test_closest_point_coupling
