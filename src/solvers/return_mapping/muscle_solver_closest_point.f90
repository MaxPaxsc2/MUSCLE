! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_solver_closest_point
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_base, only : Base_hardening_laws
    use muscle_yield_base, only : Base_yield_critera
    use muscle_elasticity_base, only : Base_elasticity
    use muscle_plastic_history, only : Plastic_material_history, Plastic_state

    private

    integer, parameter, public :: STATUS_UNDEFINED=-1
    integer, parameter, public :: STATUS_ELASTIC_CASE=0
    integer, parameter, public :: STATUS_CONVERGED=1
    integer, parameter, public :: STATUS_NONCONVERGED=2
    integer, parameter, public :: STATUS_ITER_CONVERGED    = 10
    integer, parameter, public :: STATUS_ITER_NONCONVERGED = 11

    ! Relative deviatoric part below which a stress is taken as the apex p*I of a cone criterion
    real(real64), parameter :: TOL_HYDROSTATIC = 1D-9


    type, public :: Closest_point
        class(Base_elasticity), allocatable     :: elasticity
        class(Base_hardening_laws), allocatable :: hardening
        class(Base_yield_critera), allocatable  :: yield
        integer, private                        :: iter_nw = 200
        
        contains
        procedure, public :: init => closest_point_init
        procedure, public :: solve => closest_point_solve
        procedure, public :: iter => closest_point_iter
        procedure, public :: tangent => closest_point_tangent
        procedure, public :: tangent_numerical => closest_point_tangent_numerical
    end type Closest_point

    public :: closest_point2

    contains

        pure subroutine closest_point_init(self, elasticity, hardening, yield, iter_nw)
            implicit none
            class(Closest_point), intent(inout) :: self
            class(Base_elasticity), intent(in) :: elasticity
            class(Base_hardening_laws), intent(in) :: hardening
            class(Base_yield_critera), intent(in) :: yield
            integer, optional, intent(in) :: iter_nw
            self%elasticity = elasticity
            self%hardening = hardening
            self%yield = yield

            self%iter_nw = 200
            if (present(iter_nw)) self%iter_nw = iter_nw
        end subroutine closest_point_init

        pure subroutine closest_point_iter(self, strain, history, dgamma, omega, ddgamma, status)
            !! Performs a single Newton-Raphson iteration step with adaptive under-relaxation omega.
            implicit none
            class(Closest_point), intent(in)              :: self
            type(ten_3D2Osym), intent(in)                 :: strain   !! Input total strain tensor
            type(Plastic_material_history), intent(inout) :: history  !! Material history container
            real(real64), intent(inout)                   :: dgamma   !! Accumulated plastic multiplier increment
            real(real64), intent(inout)                   :: omega    !! Relaxation factor (0 < omega <= 1)
            real(real64), intent(out)                     :: ddgamma  !! Newton-Raphson update step (delta^2 gamma)
            integer, intent(out)                          :: status   !! Iteration status

            real(real64), parameter :: TOL = 1D-5

            ! Local variables for the current iteration
            type(ten_3D2Osym)  :: strain_p_n, strain_p
            real(real64)       :: strain_pf_n, strain_pf
            type(ten_3D2Osym)  :: df, residual1, dstrain_p, stress
            type(ten_3D4O3sym) :: elas_tan, hess, ddf
            real(real64)       :: hard, dhard, f, norm_res

            ! 1. Read state t_n
            strain_p_n  = history%state_n%strain_p
            strain_pf_n = history%state_n%strain_pf

            ! 2. Read state t_n+1 (iter k)
            strain_p  = history%state_np1%strain_p
            strain_pf = history%state_np1%strain_pf

            ! 3. Elastic Tangent & Current Stress
            elas_tan = self%elasticity%dstress_dstrain(strain - strain_p)
            stress   = self%elasticity%stress(strain - strain_p)
            hard     = self%hardening%stress(strain_pf)
            f        = self%yield%stress_eq(stress) - hard

            ! 4. Derivatives & Plastic Strain Residual
            dhard     = self%hardening%dstress_dep(strain_pf)
            df        = self%yield%dstressEq_dstress(stress)
            residual1 = strain_p_n - strain_p + (dgamma * df)
            norm_res  = sqrt(sum(residual1%vals**2))

            ! 5. Check Convergence
            if (abs(f) < TOL .and. norm_res < TOL) then
                status = STATUS_ITER_CONVERGED
                history%state_np1%stress = stress
                return
            end if

            ! 6. Compute Hessian & Raw Newton-Raphson correction (ddgamma)
            ddf  = self%yield%ddstressEq_ddstress(stress)
            hess = .inv. ((.inv. elas_tan) + dgamma * ddf)

            ddgamma = (f - (df .ddot. hess .ddot. residual1)) / ((df .ddot. hess .ddot. df) + dhard)

            ! 7. Apply Under-Relaxed Step: dgamma = dgamma + omega * ddgamma
            dgamma = dgamma + omega * ddgamma

            ! Enforce non-negativity constraint on plastic multiplier
            if (dgamma < 0.0D0) then
                dgamma = 0.0D0
                omega  = omega * 0.75D0  ! Reduce relaxation factor if hitting non-negative boundary
            end if

            strain_pf = strain_pf_n + dgamma

            ! The plastic strain takes the same relaxed step as dgamma, so the iterate
            ! stays on the Newton direction when omega < 1 (except in an iteration where
            ! dgamma is clipped at zero)
            dstrain_p = ((.inv. elas_tan) .ddot. hess) .ddot. (residual1 + ddgamma * df)
            strain_p  = strain_p + omega * dstrain_p

            ! 8. Update candidate state t_n+1 (iter k+1)
            history%state_np1%strain_p  = strain_p
            history%state_np1%strain_pf = strain_pf
            history%state_np1%stress    = stress
            status = STATUS_ITER_NONCONVERGED

        end subroutine closest_point_iter

        pure subroutine closest_point_solve(self, strain, history, status, iters)
            implicit none
            class(Closest_point), intent(in)              :: self
            type(ten_3D2Osym), intent(in)                 :: strain
            type(Plastic_material_history), intent(inout) :: history
            integer, intent(out), optional                :: status
            integer, intent(out), optional                :: iters

            real(real64), parameter :: TOL2 = 1D-8
            type(ten_3D2Osym)       :: stress_trial
            real(real64)            :: hard_n, f_trial
            real(real64)            :: dgamma, ddgamma, omega
            integer                 :: i, iter_status, local_status
            real(real64)            :: d_apex
            logical                 :: apex_accepted


           ! *** Step 1: Elastic Trial Check ***
            stress_trial = self%elasticity%stress(strain - history%state_n%strain_p)
            hard_n       = self%hardening%stress(history%state_n%strain_pf)
            f_trial      = self%yield%stress_eq(stress_trial) - hard_n

            if (f_trial / max(hard_n, 1.0D-10) <= TOL2) then  ! Elastic Case
                history%state_np1%stress    = stress_trial
                history%state_np1%strain_p  = history%state_n%strain_p
                history%state_np1%strain_pf = history%state_n%strain_pf

                if (present(status)) status = STATUS_ELASTIC_CASE
                if (present(iters))  iters  = 0
                return
            end if

            ! *** Step 1b: Return to the Apex of a Cone Criterion ***
            d_apex = self%yield%apex_slope()
            if (d_apex > 0.0D0) then
                call closest_point_apex(self, strain, history, d_apex, hard_n, apex_accepted, i)
                if (apex_accepted) then
                    if (present(status)) status = STATUS_CONVERGED
                    if (present(iters))  iters  = i
                    return
                end if
            end if

            ! *** Step 2: Initialize Candidate State and Relaxation Factor ***
            history%state_np1 = history%state_n
            dgamma       = 0.0D0
            ddgamma      = 0.0D0
            omega        = 1.0D0  ! Initial unrelaxed factor
            local_status = STATUS_NONCONVERGED

            ! *** Step 3: Newton-Raphson Loop with Adaptive Relaxation ***
            do i = 1, self%iter_nw
                ! Stagnation protection: reduce omega if taking too many iterations
                if (i == 10 .or. i == 20 .or. i == 30 .or. i == 40) then
                    omega = omega * 0.75D0
                end if

                call self%iter(strain, history, dgamma, omega, ddgamma, iter_status)

                if (iter_status == STATUS_ITER_CONVERGED) then
                    local_status = STATUS_CONVERGED
                    exit
                end if
            end do

            if (present(status)) status = local_status
            if (present(iters))  iters  = i
        end subroutine closest_point_solve

        pure subroutine closest_point_apex(self, strain, history, d, hard_n, accepted, iters)
            !! Return to the apex p*I of a cone criterion with slope d (Sysala et al., 2016, Eq. 3.27
            !! and 3.28, with the yield function of MUSCLE scaled by a = 1/(sqrt(3) f): eta = eta_bar
            !! = a*d, xi = a, Delta lambda = dgamma/a, K = 1/w): Newton from dgamma = 0 on
            !! g = d*p - sigma_y(eps_pf_n + dgamma), p = (tr(eps_e_trial) - d*dgamma)/w, w = I:C^-1:I.
            !! Accepted only if dev(Delta eps_p) flows from the apex (apex_dev_gauge <= dgamma),
            !! sigma_y' >= 0 and the stress from elasticity%stress is p*I on the yield surface.
            !! Called by solve; history is written only when accepted.
            implicit none
            class(Closest_point), intent(in)              :: self
            type(ten_3D2Osym), intent(in)                 :: strain   !! Input total strain tensor
            type(Plastic_material_history), intent(inout) :: history  !! Material history container
            real(real64), intent(in)                      :: d        !! Slope of the cone (apex_slope)
            real(real64), intent(in)                      :: hard_n   !! Yield stress at t_n
            logical, intent(out)                          :: accepted !! Apex return accepted
            integer, intent(out)                          :: iters    !! Scalar Newton iterations

            real(real64), parameter :: TOL = 1D-5  ! same tolerance as closest_point_iter
            type(ten_3D2Osym) :: strain_e, cinv_i, dstrain_p, strain_p, stress, eye
            type(iden_2O)     :: I2O
            real(real64)      :: w, tr_e, dgamma, hard, dhard, g, p
            integer           :: i

            accepted = .false.
            iters    = 0
            strain_e = strain - history%state_n%strain_p
            eye      = I2O
            tr_e     = strain_e .ddot. I2O

            ! w > 0, so tr(eps_e_trial) <= 0 gives g(0) <= 0: no apex, and no inverse is needed
            if (tr_e <= 0.0D0) return

            ! C^-1 : I and w = I : C^-1 : I (w = 1/K_b for isotropic elasticity)
            cinv_i = (.inv. self%elasticity%dstress_dstrain(strain_e)) .ddot. eye
            w      = cinv_i .ddot. I2O

            ! g decreases when sigma_y' >= 0: g(0) <= 0 means no return to the apex
            dgamma = 0.0D0
            hard   = hard_n
            g      = d * tr_e / w - hard
            if (g <= 0.0D0) return

            do i = 1, self%iter_nw
                dhard = self%hardening%dstress_dep(history%state_n%strain_pf + dgamma)
                if (dhard < 0.0D0) return
                dgamma = dgamma + g / (d * d / w + dhard)
                hard   = self%hardening%stress(history%state_n%strain_pf + dgamma)
                g      = d * (tr_e - d * dgamma) / w - hard
                iters  = i
                if (abs(g) < TOL) exit
            end do
            if (abs(g) >= TOL) return

            ! Delta eps_p = eps_e_trial - p C^-1 : I; its deviatoric part must flow from the apex
            p         = (tr_e - d * dgamma) / w
            dstrain_p = strain_e - p * cinv_i
            if (self%yield%apex_dev_gauge(.dev. dstrain_p) > dgamma) return

            strain_p = history%state_n%strain_p + dstrain_p
            stress   = self%elasticity%stress(strain - strain_p)
            if (.not. stress%is_approx(p * eye, tol=TOL_HYDROSTATIC)) return
            if (abs(self%yield%stress_eq(stress) - hard) >= TOL) return

            history%state_np1%strain_p  = strain_p
            history%state_np1%strain_pf = history%state_n%strain_pf + dgamma
            history%state_np1%stress    = stress
            accepted = .true.
        end subroutine closest_point_apex



        pure subroutine closest_point_tangent(self, strain, history, tangent)
            !! Computes the consistent algorithmic tangent stiffness tensor C_mat = dS/dE.
            !! Pure read-only function that evaluates the tangent at the converged state t_n+1.
            implicit none
            class(Closest_point), intent(in)           :: self
            type(ten_3D2Osym), intent(in)              :: strain   !! Input total strain tensor
            type(Plastic_material_history), intent(in) :: history  !! Material history container (Read-Only)
            type(ten_3D4O2sym), intent(out)            :: tangent  !! Algorithmic tangent tensor C_mat

            type(ten_3D2Osym)  :: df, stress, strain_p
            type(ten_3D4O3sym) :: elas_tan, hess, ddf
            real(real64)       :: dhard, dgamma, strain_pf
            type(ten_3D2Osym)  :: eye
            type(ten_3D4O3sym) :: apex_tan
            real(real64)       :: d, w
            type(iden_2O)      :: I2O
            type(iden_4O3T)    :: I4O3T

            ! 1. Read state variables at t_n+1 from history
            strain_p  = history%state_np1%strain_p
            strain_pf = history%state_np1%strain_pf
            stress    = history%state_np1%stress

            ! 2. Compute step plastic multiplier increment Delta gamma = bar_eps_p_{n+1} - bar_eps_p_n
            dgamma = history%state_np1%strain_pf - history%state_n%strain_pf

            ! 3. Compute elastic tangent at current elastic strain
            elas_tan = self%elasticity%dstress_dstrain(strain - strain_p)

            ! 4. Check for elastic step (Delta gamma = 0)
            if (dgamma <= 1.0D-12) then
                tangent = elas_tan
                return
            end if

            ! 4b. Apex of a cone criterion: C_ep = sigma_y' / (d**2 + w*sigma_y') I x I, with
            ! w = I : C^-1 : I (Sysala et al., 2016, Eq. 3.32 with K = 1/w, scaled as in closest_point_apex)
            d = self%yield%apex_slope()
            if (d > 0.0D0) then
                if (stress%is_approx(stress - (.dev. stress), tol=TOL_HYDROSTATIC)) then
                    dhard    = self%hardening%dstress_dep(strain_pf)
                    eye      = I2O
                    w        = ((.inv. elas_tan) .ddot. eye) .ddot. I2O
                    apex_tan = (dhard / (d * d + w * dhard)) * I4O3T
                    tangent  = apex_tan
                    return
                end if
            end if

            ! 5. Compute derivatives at converged stress state
            dhard = self%hardening%dstress_dep(strain_pf)
            df    = self%yield%dstressEq_dstress(stress)
            ddf   = self%yield%ddstressEq_ddstress(stress)

            ! 6. Algorithmic Hessian matrix: H_alg = [ C^-1 + dgamma * d2f/dsigma2 ]^-1
            hess = .inv. ((.inv. elas_tan) + dgamma * ddf)

            ! 7. Consistent Elastoplastic Tangent Modulus:
            ! C_ep = H_alg - (H_alg : N x N : H_alg) / (N : H_alg : N + H')
            tangent = hess - (.tdotsym. (df .ddot. hess)) / ((df .ddot. hess .ddot. df) + dhard)

        end subroutine closest_point_tangent

        subroutine closest_point_tangent_numerical(self, strain, history, tangent)
            use muscle_math_derivatives
            implicit none
            class(Closest_point), intent(in)           :: self
            type(ten_3D2Osym), intent(in)              :: strain
            type(Plastic_material_history), intent(in) :: history
            type(ten_3D4O2sym), intent(out)            :: tangent

            tangent = derivative(wrapper, strain)

            contains
            pure function wrapper(x) result(stress_out)
                implicit none
                type(ten_3D2Osym), intent(in) :: x
                type(ten_3D2Osym)              :: stress_out
                type(Plastic_material_history) :: local_history

                local_history = history
                call self%solve(x, local_history)
                stress_out = local_history%state_np1%stress
            end function wrapper
        end subroutine closest_point_tangent_numerical

end module