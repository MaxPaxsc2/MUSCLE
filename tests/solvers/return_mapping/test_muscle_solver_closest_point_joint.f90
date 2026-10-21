! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module test_closest_point_joint_helpers
    ! Von Mises criteria whose stressEq_and_derivatives returns the equivalent stress plus
    ! SHIFT, so a solution shows whether closest_point_iter called it. Shifted_vonmises_joint
    ! declares joint_evaluation(); Shifted_vonmises does not.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_vonmises, only : VonMises
    implicit none
    private
    public :: Shifted_vonmises, Shifted_vonmises_joint, SHIFT

    real(real64), parameter :: SHIFT = 20D0  ! Added to the equivalent stress of the joint call

    type, extends(VonMises) :: Shifted_vonmises
    contains
        procedure :: stressEq_and_derivatives => shifted_stressEq_and_derivatives
    end type Shifted_vonmises

    type, extends(Shifted_vonmises) :: Shifted_vonmises_joint
    contains
        procedure, nopass :: joint_evaluation => joint_true
    end type Shifted_vonmises_joint

contains

    pure subroutine shifted_stressEq_and_derivatives(self, stress, f, df, ddf)
        ! Von Mises stress_eq + SHIFT, with the von Mises gradient and Hessian
        class(Shifted_vonmises), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
        real(real64), intent(out) :: f
        type(ten_3D2Osym), intent(out) :: df
        type(ten_3D4O3sym), intent(out) :: ddf

        f = self%stress_eq(stress) + SHIFT
        df = self%dstressEq_dstress(stress)
        ddf = self%ddstressEq_ddstress(stress)
    end subroutine shifted_stressEq_and_derivatives

    pure function joint_true() result(res)
        logical :: res

        res = .true.
    end function joint_true

end module test_closest_point_joint_helpers

program test_muscle_solver_closest_point_joint
    implicit none
    logical :: passed

    call test_closest_point_joint_evaluation_used(passed)
    if (.not. passed) STOP 1

    call test_closest_point_joint_evaluation_ignored(passed)
    if (.not. passed) STOP 2

    print*, "Passed!", passed
    STOP 0
end program test_muscle_solver_closest_point_joint

subroutine test_closest_point_joint_evaluation_used(passed)
    ! A criterion that declares joint_evaluation() is evaluated through stressEq_and_derivatives
    ! in the Newton iterations: the shifted von Mises must give the von Mises solution with the
    ! initial yield stress lowered by SHIFT (the elastic trial check still uses stress_eq, and the
    ! step is plastic in both). Plastic step with the six strain components.
    ! Called by program test_muscle_solver_closest_point_joint (STOP 1)
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_vonmises, only : VonMises
    use muscle_hard_bilinear, only : Bilinear_hardening
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    use test_closest_point_joint_helpers, only : Shifted_vonmises_joint, SHIFT
    implicit none
    logical, intent(out) :: passed

    ! Same Newton iterates up to the rounding of (f + SHIFT) - hard against f - (hard - SHIFT).
    ! Measured: 0 in this step; 1e-12 (~5000 eps) covers that rounding in general
    real(real64), parameter :: TOL = 1D-12
    type(Shifted_vonmises_joint) :: shifted
    type(VonMises) :: vm
    type(Elasticity_linear) :: elas
    type(Plastic_material_history) :: h_shifted, h_ref
    type(Closest_point) :: solver, reference
    type(ten_3D2Osym) :: strain
    integer :: status, status_ref

    passed = .true.
    call elas%set_parameters(young=200000D0, poisson=0.3D0)
    call strain%init(xx=4D-3, yy=-1.5D-3, zz=-1D-3, xy=1.2D-3, yz=-0.8D-3, xz=0.6D-3)
    call solver%init(elasticity=elas, hardening=Bilinear_hardening(y0=300D0, K=1000D0), yield=shifted)
    call reference%init(elasticity=elas, hardening=Bilinear_hardening(y0=300D0 - SHIFT, K=1000D0), &
                        yield=vm)
    call h_shifted%init()
    call h_ref%init()
    call solver%solve(strain, h_shifted, status)
    call reference%solve(strain, h_ref, status_ref)

    call check_true(passed, "Joint evaluation used: shifted step converged", status == STATUS_CONVERGED)
    call check_true(passed, "Joint evaluation used: reference step converged", status_ref == STATUS_CONVERGED)
    call check_approx(passed, "Joint evaluation used: stress", h_shifted%state_np1%stress, &
                      h_ref%state_np1%stress, tol=TOL)
    call check_approx(passed, "Joint evaluation used: plastic strain", h_shifted%state_np1%strain_p, &
                      h_ref%state_np1%strain_p, tol=TOL)
    call check_approx(passed, "Joint evaluation used: equivalent plastic strain", &
                      h_shifted%state_np1%strain_pf, h_ref%state_np1%strain_pf, tol=TOL)
end subroutine test_closest_point_joint_evaluation_used

subroutine test_closest_point_joint_evaluation_ignored(passed)
    ! A criterion that overrides stressEq_and_derivatives without declaring joint_evaluation()
    ! is evaluated through stress_eq, dstressEq_dstress and ddstressEq_ddstress: the shifted
    ! von Mises must give the plain von Mises solution, bit for bit. Same step as STOP 1.
    ! Called by program test_muscle_solver_closest_point_joint (STOP 2)
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_vonmises, only : VonMises
    use muscle_hard_bilinear, only : Bilinear_hardening
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    use test_closest_point_joint_helpers, only : Shifted_vonmises
    implicit none
    logical, intent(out) :: passed

    ! Same calls as von Mises, so the same results bit for bit. Measured: 0
    real(real64), parameter :: TOL = 0D0
    type(Shifted_vonmises) :: shifted
    type(VonMises) :: vm
    type(Elasticity_linear) :: elas
    type(Plastic_material_history) :: h_shifted, h_ref
    type(Closest_point) :: solver, reference
    type(ten_3D2Osym) :: strain
    integer :: status, status_ref

    passed = .true.
    call elas%set_parameters(young=200000D0, poisson=0.3D0)
    call strain%init(xx=4D-3, yy=-1.5D-3, zz=-1D-3, xy=1.2D-3, yz=-0.8D-3, xz=0.6D-3)
    call solver%init(elasticity=elas, hardening=Bilinear_hardening(y0=300D0, K=1000D0), yield=shifted)
    call reference%init(elasticity=elas, hardening=Bilinear_hardening(y0=300D0, K=1000D0), yield=vm)
    call h_shifted%init()
    call h_ref%init()
    call solver%solve(strain, h_shifted, status)
    call reference%solve(strain, h_ref, status_ref)

    call check_true(passed, "Joint evaluation ignored: shifted step converged", status == STATUS_CONVERGED)
    call check_true(passed, "Joint evaluation ignored: reference step converged", status_ref == STATUS_CONVERGED)
    call check_approx(passed, "Joint evaluation ignored: stress", h_shifted%state_np1%stress, &
                      h_ref%state_np1%stress, tol=TOL)
    call check_approx(passed, "Joint evaluation ignored: plastic strain", h_shifted%state_np1%strain_p, &
                      h_ref%state_np1%strain_p, tol=TOL)
    call check_approx(passed, "Joint evaluation ignored: equivalent plastic strain", &
                      h_shifted%state_np1%strain_pf, h_ref%state_np1%strain_pf, tol=TOL)
end subroutine test_closest_point_joint_evaluation_ignored
