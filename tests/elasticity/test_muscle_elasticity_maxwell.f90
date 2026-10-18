! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_elasticity_maxwell
    implicit none
    logical :: passed
    call test_maxwell_equilibrium_response(passed)
    if (.not. passed) STOP 1
    print *, 'Passed!', passed
    STOP 0
end program test_muscle_elasticity_maxwell

subroutine test_maxwell_equilibrium_response(passed)
    ! History-free calls return C_inf:strain with branches present; coupled Hooke, all shears.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_elasticity_linear
    use muscle_elasticity_maxwell
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(GeneralizedMaxwell3D) :: law
    type(MaxwellInternalVars) :: iv
    type(Elasticity_linear) :: el
    type(ten_3D2Osym) :: eps, expected, zero
    type(ten_2D2Osym) :: eps2
    type(ten_3D4O3sym) :: c
    integer :: k
    ! Coupled Hooke round-off. Measured: 1.55e-16 or less.
    real(real64), parameter :: TOL = 2D-14
    passed = .true.
    call law%set_parameters(120D0, 0.25D0, [0.2D0,0.7D0,1.1D0], [0.5D0,5D0,40D0])
    call law%init_internal(iv)
    call el%set_parameters(120D0, 0.25D0)
    call eps%init([0.01D0,0.02D0,0.03D0,0.04D0,0.05D0,0.06D0])
    call expected%init([3.84D0,4.8D0,5.76D0,3.84D0,4.8D0,5.76D0])
    call check_approx(passed, 'Hooke coupled independent values', el%stress(eps), expected, tol=TOL)
    c = el%dstress_dstrain(eps)
    call check_approx(passed, 'Hooke tangent coupled action', c .ddot. eps, expected, tol=TOL)
    call check_approx(passed, 'three branches equilibrium stress', &
                      law%stress(eps), expected, tol=TOL)
    call check_approx(passed, 'three branches equilibrium tangent', &
                      law%dstress_dstrain(eps), c, tol=TOL)
    call eps2%init([0.01D0,0.02D0,0.03D0,0.04D0])
    call check_approx(passed, '2D equilibrium stress', law%stress(eps2), el%stress(eps2), tol=TOL)
    call check_approx(passed, '2D equilibrium tangent', &
                      law%dstress_dstrain(eps2), el%dstress_dstrain(eps2), tol=TOL)
    call check_true(passed, 'one branch stress per branch', size(iv%q) == 3)
    zero = 0D0
    do k = 1, size(iv%q)
        call check_approx(passed, 'branch stress starts at zero', iv%q(k), zero, tol=TOL)
    end do
end subroutine test_maxwell_equilibrium_response
