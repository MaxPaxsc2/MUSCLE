! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_elasticity_maxwell
    !! Generalized Maxwell model, small strains and a common Poisson ratio.
    !! M1: dq_k/dt + q_k/tau_k = alpha_k d(C_inf:strain)/dt.
    !! History-free calls return the equilibrium response C_inf:strain.
    !! Initialize material and state before stepping.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_elasticity_base
    use muscle_elasticity_linear, only : Elasticity_linear
    implicit none
    private

    public :: MaxwellInternalVars
    type, extends(Base_internal_elasticity_variables) :: MaxwellInternalVars
        type(ten_3D2Osym), allocatable :: q(:) !! Branch stresses q_k.
    end type MaxwellInternalVars

    public :: GeneralizedMaxwell3D
    type, extends(Base_elasticity) :: GeneralizedMaxwell3D
        real(real64), private :: E_inf, nu !! Equilibrium Young modulus and Poisson ratio.
        real(real64), allocatable, private :: alpha(:), tau(:) !! E_k/E_inf and relaxation times.
        type(Elasticity_linear), private :: equilibrium !! Cached equilibrium elasticity C_inf.
    contains
        procedure :: set_parameters
        procedure :: init_internal
        procedure :: stress_3D => stress_maxwell_3D
        procedure :: stress_2D => stress_maxwell_2D
        procedure :: dstress_dstrain_3D => tangent_maxwell_3D
        procedure :: dstress_dstrain_2D => tangent_maxwell_2D
        procedure :: to_string => to_string_maxwell
    end type GeneralizedMaxwell3D

contains

    pure subroutine set_parameters(self, young_inf, poisson, alpha, tau)
        !! Initialize M1 parameters and Hooke cache; called by the material owner.
        implicit none
        class(GeneralizedMaxwell3D), intent(inout) :: self
        real(real64), intent(in) :: young_inf, poisson, alpha(:), tau(:)
        if (size(alpha) /= size(tau)) error stop 'Maxwell: alpha/tau size mismatch'
        if (young_inf <= 0D0 .or. poisson <= -1D0 .or. poisson >= 0.5D0) &
            error stop 'Maxwell: invalid equilibrium parameters'
        if (any(alpha < 0D0) .or. any(tau <= 0D0)) error stop 'Maxwell: invalid branch parameters'
        self%E_inf = young_inf
        self%nu = poisson
        self%alpha = alpha
        self%tau = tau
        call self%equilibrium%set_parameters(young=young_inf, poisson=poisson)
    end subroutine set_parameters

    pure subroutine init_internal(self, internal)
        !! Allocate and reset M1 branch stresses; called by the material owner once per point.
        implicit none
        class(GeneralizedMaxwell3D), intent(in) :: self
        type(MaxwellInternalVars), intent(out) :: internal
        integer :: k
        if (.not. allocated(self%alpha)) error stop 'Maxwell: parameters not set'
        allocate(internal%q(size(self%alpha)))
        do k = 1, size(internal%q)
            internal%q(k) = 0D0
        end do
    end subroutine init_internal

    pure function stress_maxwell_3D(self, strain) result(res)
        !! Equilibrium stress from M1 at infinite time; called by generic stress.
        implicit none
        class(GeneralizedMaxwell3D), intent(in) :: self
        class(ten_3D2Osym), intent(in) :: strain
        type(ten_3D2Osym) :: res
        res = self%equilibrium%stress(strain)
    end function stress_maxwell_3D

    pure function stress_maxwell_2D(self, strain) result(res)
        !! 2D equilibrium stress from M1; called by generic stress (no 2D history).
        implicit none
        class(GeneralizedMaxwell3D), intent(in) :: self
        class(ten_2D2Osym), intent(in) :: strain
        type(ten_2D2Osym) :: res
        res = self%equilibrium%stress(strain)
    end function stress_maxwell_2D

    pure function tangent_maxwell_3D(self, strain) result(res)
        !! Equilibrium tangent C_inf; called by generic dstress_dstrain.
        implicit none
        class(GeneralizedMaxwell3D), intent(in) :: self
        class(ten_3D2Osym), intent(in) :: strain
        type(ten_3D4O3sym) :: res
        res = self%equilibrium%dstress_dstrain(strain)
    end function tangent_maxwell_3D

    pure function tangent_maxwell_2D(self, strain) result(res)
        !! 2D equilibrium tangent C_inf; called by generic dstress_dstrain.
        implicit none
        class(GeneralizedMaxwell3D), intent(in) :: self
        class(ten_2D2Osym), intent(in) :: strain
        type(ten_2D4O3sym) :: res
        res = self%equilibrium%dstress_dstrain(strain)
    end function tangent_maxwell_2D

    pure function to_string_maxwell(self) result(res)
        !! Describe M1 parameters; called by the material owner for diagnostics.
        implicit none
        class(GeneralizedMaxwell3D), intent(in) :: self
        character(len=275) :: res
        if (.not. allocated(self%alpha)) then
            res = 'GeneralizedMaxwell3D: parameters not set'
            return
        end if
        write(res, '(a,es14.6,a,es14.6,a,i0)') 'GeneralizedMaxwell3D: E_inf=', &
            self%E_inf, ', nu=', self%nu, ', branches=', size(self%alpha)
    end function to_string_maxwell

end module muscle_elasticity_maxwell
