! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_math_spectral_derivs
    !! # Module mod_muscle_math_spectral_derivs
    !!
    !! This module provides two routes to derivatives of spectral functions of a symmetric
    !! second-order tensor:
    !!
    !! * [[dEigenvalues_dTensor]] / [[d2Eigenvalues_dTensor2]]: closed-form, eigenvector-free
    !!   first and second derivatives of the eigenvalues (legacy; see the warning below).
    !! * [[eigen_sym3]] with [[spectral_rotate]] / [[spectral_compose]]: a cyclic Jacobi
    !!   eigen-decomposition and the Daleckii-Krein divided differences built on it. This
    !!   route does use eigenvectors, and is the one used by the anisotropic yield criteria.
    !!
    !! ## Mathematical Formulation (eigenvector-free route)
    !!
    !! The eigenvalues \(\Sigma_a\) (\(a=1,2,3\)) are roots of the characteristic polynomial:
    !!
    !! \[ P(\Sigma_a) = \Sigma_a^3 - I_1 \Sigma_a^2 + I_2 \Sigma_a - I_3 = 0 \]
    !!
    !! By applying implicit differentiation and the Cayley-Hamilton theorem, the derivatives
    !! are calculated without computing eigenvectors or performing matrix inversions, which
    !! stays well defined even for singular tensors (e.g., when \(\det(\Sigma) = 0\)).
    !!
    !! See the technical documentation for further algebraic details.
    !!
    !! @warning The eigenvector-free routines [[dEigenvalues_dTensor]] and
    !! [[d2Eigenvalues_dTensor2]] are superseded for new work, not withdrawn: finite-strain
    !! kinematics ([[muscle_kin_finite_base]]) still calls [[dEigenvalues_dTensor]], while
    !! [[d2Eigenvalues_dTensor2]] is now exercised only by the tests. Their repeated-root
    !! detection uses absolute thresholds on
    !! the characteristic-polynomial derivative (not scale invariant), so nearly repeated
    !! eigenvalues are misclassified as repeated, and the repair of a repeated root misses
    !! part of the second derivative. New code should use [[eigen_sym3]] with
    !! [[spectral_rotate]] / [[spectral_compose]].
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    implicit none
    private

    public :: dEigenvalues_dTensor
    public :: d2Eigenvalues_dTensor2
    public :: eigen_sym3
    public :: spectral_rotate
    public :: spectral_compose

    real(real64), parameter :: EPS_TOL = 1.0D-11
    real(real64), parameter :: EPS_TOL_SECOND = 1.0D-6
    !! Safe tolerance threshold to prevent division by zero at repeated roots (singularities)

    real(real64), parameter :: JACOBI_TOL = 1.0D-3*epsilon(1.0D0)
        !! Off-diagonal entries below JACOBI_TOL*||Sigma||_F are treated as zero (scale invariant)
    integer, parameter :: JACOBI_MAX_SWEEPS = 50
        !! Safety bound; cyclic Jacobi on 3x3 converges quadratically in a few sweeps
    real(real64), parameter :: JACOBI_HUGE_THETA = 1.0D150
        !! Above this |theta|, tan(phi) = 1/(2 theta) avoids overflow in theta**2

contains
    pure function dEigenvalues_dTensor(Sigma, eigenvalues) result(res)
        !! # First Derivative of Eigenvalues
        !!
        !! Computes the first analytical derivative of all three eigenvalues of a symmetric 
        !! second-order tensor \(\Sigma_{ij}\) with respect to itself:
        !!
        !! \[ H_{mn}^{(a)} = \frac{\partial \Sigma_a}{\partial \Sigma_{mn}} \]
        !!
        !! ## Mathematical Formulation
        !! The computation is based on the implicit differentiation of the characteristic 
        !! polynomial and the Cayley-Hamilton theorem. This approach yields a closed-form 
        !! expression free of eigenvectors:
        !!
        !! \[ H_{mn}^{(a)} = \frac{\Sigma^2_{mn} - (I_1 - \Sigma_a)\Sigma_{mn} + (\Sigma_a^2 - I_1\Sigma_a + I_2)\delta_{mn}}{3\Sigma_a^2 - 2I_1\Sigma_a + I_2} \]
        !!
        !! ## Singularity Handling (Repeated Roots)
        !! When eigenvalues coincide (e.g., uniaxial tension or hydrostatic states), the denominator 
        !! approaches zero, leading to an indeterminate form \(0/0\). This routine robustly handles 
        !! such states by exploiting the identity property \(\sum_{a=1}^3 \frac{\partial \Sigma_a}{\partial \pmb{\Sigma}} = \mathbf{I}\).
        !! 
        !! - **Two identical eigenvalues:** The remaining subspace is isotropically divided among the repeated roots.
        !! - **Three identical eigenvalues:** The identity tensor is isotropically divided by 3.
        !!
        !! @note This function is pure and returns an array of three symmetric second-order tensors.
        !!
        !! @warning Deprecated for new code (see the module note). An eigenvalue whose
        !! denominator \((\Sigma_a-\Sigma_b)(\Sigma_a-\Sigma_c)\) is below EPS_TOL = 1e-11 in
        !! absolute value is treated as exactly repeated, which returns averaged projectors
        !! instead of the true ones (for unit-sized tensors: a gap below ~1e-11 next to a
        !! distant eigenvalue, or below ~3e-6 when all three are clustered). If only one
        !! denominator falls below the threshold, that entry of `res` is left undefined.
        type(ten_3D2Osym), intent(in)  :: Sigma          !! Symmetric second-order tensor \(\Sigma_{ij}\)
        real(real64), intent(in)       :: eigenvalues(3) !! Evaluated eigenvalues \(\Sigma_1, \Sigma_2, \Sigma_3\)
        type(ten_3D2Osym)              :: res(3)         !! Resulting derivatives \(H_{ij}^{(a)}\) for \(a=1,2,3\)
        
        real(real64)      :: I1, I2, D_scal(3), lam
        type(ten_3D2Osym) :: F_tensor
        type(iden_2O)     :: I2O
        integer           :: a, num_sing
        integer           :: good_idx, sing_idx(3)

        ! 1. Compute Invariants
        ! First invariant: trace(\Sigma)
        I1 = Sigma%vals(1) + Sigma%vals(2) + Sigma%vals(3)
        ! Second invariant: 0.5*(I1^2 - tr(\Sigma^2))
        I2 = 0.5D0 * (I1**2 - (Sigma .ddot. Sigma))
        ! Second-order identity tensor (\delta_ij)

        num_sing = 0
        good_idx = 1 ! Default initialization

        ! 2. First Pass: Compute denominators and identify singular roots
        do a = 1, 3
            lam = eigenvalues(a)
            ! Characteristic polynomial derivative (scalar denominator)
            D_scal(a) = (3.0D0*lam - 2.0D0*I1)*lam + I2
            
            ! Check for distinct eigenvalues (well-conditioned denominator)
            if (abs(D_scal(a)) > EPS_TOL) then
                ! Numerator tensor using Cayley-Hamilton expression
                F_tensor = Sigma%square() - (I1 - lam)*Sigma + (lam*(lam - I1))*I2O + I2*I2O 
                res(a) = F_tensor / D_scal(a)
                good_idx = a
            else
                ! Mark as repeated/singular root for post-processing
                num_sing = num_sing + 1
                sing_idx(num_sing) = a
            end if
        end do

        ! 3. Second Pass: Repair singular roots using orthogonal projections
        if (num_sing == 2) then
            ! Two repeated eigenvalues (e.g., pure uniaxial tension/compression)
            ! The degenerate subspace is isotropically averaged from the remaining identity space.
            F_tensor = 0.5D0 * (I2O - res(good_idx))
            res(sing_idx(1)) = F_tensor
            res(sing_idx(2)) = F_tensor
            
        else if (num_sing == 3) then
            ! Three repeated eigenvalues (e.g., pure hydrostatic stress)
            ! All principal directions are isotropic.
            F_tensor = (1.0D0 / 3.0D0) * I2O
            res(1) = F_tensor
            res(2) = F_tensor
            res(3) = F_tensor
        end if

    end function dEigenvalues_dTensor



    pure function d2Eigenvalues_dTensor2(Sigma, eigenvalues, H_tensors) result(res)
        !! # Second Derivative (Hessian) of Eigenvalues
        !!
        !! Computes the second analytical derivative of all three eigenvalues of a symmetric 
        !! second-order tensor \(\Sigma_{ij}\) with respect to itself:
        !!
        !! \[ \mathcal{H}_{mnop}^{(a)} = \frac{\partial^2 \Sigma_a}{\partial \Sigma_{mn} \partial \Sigma_{op}} \]
        !!
        !! This function is pure and returns an array of three fully symmetric fourth-order tensors.
        !!
        !! @warning Deprecated for new code (see the module note). Same absolute-threshold
        !! classification as [[dEigenvalues_dTensor]] (EPS_TOL_SECOND = 1e-6); for repeated or
        !! nearly repeated eigenvalues the result is not the second derivative of the
        !! corresponding spectral function. Use the Daleckii-Krein form with [[eigen_sym3]].
        implicit none
        type(ten_3D2Osym), intent(in)   :: Sigma          !! Symmetric second-order tensor \(\Sigma_{ij}\)
        real(real64), intent(in)        :: eigenvalues(3) !! Evaluated eigenvalues \(\Sigma_1, \Sigma_2, \Sigma_3\)
        type(ten_3D2Osym), intent(in)   :: H_tensors(3)   !! Pre-computed first derivatives \(H_{ij}^{(a)}\)
        type(ten_3D4O3sym)              :: res(3)         !! Resulting Hessians for \(a=1,2,3\)
        
        type(ten_3D4O3sym) :: dF_dSigma, M_tensor, test
        type(ten_3D4O3sym) :: N_tensor
        type(ten_3D2Osym) :: H, K
        type(iden_2O) :: I2O
        type(iden_4O4T) :: I4O4T
        type(iden_4O3T) :: I4O3T
        real(real64) :: I1, I2, D_scal, lam
        integer :: a, num_sing, good_idx, sing_idx(3)

        ! First invariant: trace(\Sigma)
        I1 = Sigma%vals(1) + Sigma%vals(2) + Sigma%vals(3)
        ! Second invariant
        I2 = 0.5D0 * (I1**2 - (Sigma .ddot. Sigma))
       
        ! Fourth-order tensor M_ijkl (analytical derivative of \Sigma^2) 
        call M_tensor%init( &
            xxxx=2.0D0*Sigma%xx(), xxyy=0.0D0, xxzz=0.0D0, xxxy=Sigma%xy(), xxyz=0.0D0, xxxz=Sigma%xz(), &
            yyyy=2.0D0*Sigma%yy(), yyzz=0.0D0, yyxy=Sigma%xy(), yyyz=Sigma%yz(), yyxz=0.0D0, &
            zzzz=2.0D0*Sigma%zz(), zzxy=0.0D0, zzyz=Sigma%yz(), zzxz=Sigma%xz(), &
            xyxy=0.5D0*(Sigma%xx()+Sigma%yy()), xyyz=0.5D0*Sigma%xz(), xyxz=0.5D0*Sigma%yz(), &
            yzyz=0.5D0*(Sigma%yy()+Sigma%zz()), yzxz=0.5D0*Sigma%xy(), &
            xzxz=0.5D0*(Sigma%xx()+Sigma%zz()) &
        )

        num_sing = 0
        good_idx = 1

        do a = 1, 3
            lam = eigenvalues(a)
            D_scal = (3.0D0*lam - 2.0D0*I1)*lam + I2

            if (abs(D_scal) > EPS_TOL_SECOND) then
                H = H_tensors(a)                
                ! 2. Symmetrized fourth-order N_tensor (numerator) using Cayley-Hamilton-based expression
                N_tensor = M_tensor &
                           - (I1 - lam)*I4O4T &
                           + (I1 - lam)*I4O3T &
                           - 2.0D0 * (Sigma .tdotsym. I2O) &
                           + 2.0D0 * (Sigma .tdotsym. H) &
                           + 2.0D0 * (2.0D0*lam - I1) * (I2O .tdotsym. H) &
                           - 2.0D0 * (3.0D0*lam - I1) * (.tdotsym. H)
    
                res(a) = N_tensor / D_scal
                good_idx = a
            else 
                num_sing = num_sing + 1
                sing_idx(num_sing) = a
            end if
        end do

        ! Repair singular roots using orthogonal projections
        if (num_sing == 2) then
            ! Two repeated eigenvalues: Isotropic averaging of the degenerate subspace
            res(sing_idx(1)) = (-0.5D0) * res(good_idx)
            res(sing_idx(2)) = (-0.5D0) * res(good_idx)
        else if (num_sing == 3) then
            ! Three repeated eigenvalues: Isotropic distribution of the identity subspace
            res(1) = 0.0D0
            res(2) = 0.0D0
            res(3) = 0.0D0
        end if

    end function d2Eigenvalues_dTensor2


    pure subroutine eigen_sym3(Sigma, eigenvalues, V)
        !! # Eigen-decomposition of a symmetric second-order tensor (cyclic Jacobi)
        !!
        !! Computes \(\Sigma = V\,\mathrm{diag}(\Sigma_1,\Sigma_2,\Sigma_3)\,V^T\) with orthonormal
        !! eigenvectors in the columns of `V` and eigenvalues in descending order.
        !!
        !! Jacobi rotations are backward stable and scale invariant: the result does not
        !! depend on the stress units, and repeated or nearly repeated eigenvalues are
        !! handled without special cases (the eigenvectors of a cluster are an arbitrary
        !! orthonormal basis of its invariant subspace, which is all that spectral-function
        !! derivatives need). No division by a vanishing quantity occurs.
        !!
        !! Intended for spectral-function derivatives via [[spectral_rotate]] and
        !! [[spectral_compose]]. The existing eigenvector-free routines above are unchanged.
        type(ten_3D2Osym), intent(in) :: Sigma          !! Symmetric tensor (tensorial Voigt storage)
        real(real64), intent(out)     :: eigenvalues(3) !! Eigenvalues, descending
        real(real64), intent(out)     :: V(3,3)         !! Orthonormal eigenvectors (columns)

        real(real64) :: A(3,3), tol, theta, t, c, s, app, aqq, apq, arp, arq, tmp, col(3)
        integer :: sweep, k, p, q, r, i, j
        integer, parameter :: PP(3) = [1, 1, 2]
        integer, parameter :: QQ(3) = [2, 3, 3]
        logical :: rotated

        A(1,1) = Sigma%vals(1)
        A(2,2) = Sigma%vals(2)
        A(3,3) = Sigma%vals(3)
        A(1,2) = Sigma%vals(4)
        A(2,1) = A(1,2)
        A(2,3) = Sigma%vals(5)
        A(3,2) = A(2,3)
        A(1,3) = Sigma%vals(6)
        A(3,1) = A(1,3)

        V = 0.0D0
        V(1,1) = 1.0D0
        V(2,2) = 1.0D0
        V(3,3) = 1.0D0

        tol = JACOBI_TOL*sqrt(sum(A**2))

        do sweep = 1, JACOBI_MAX_SWEEPS
            rotated = .false.
            do k = 1, 3
                p = PP(k)
                q = QQ(k)
                r = 6 - p - q
                apq = A(p,q)
                if (abs(apq) <= tol) then
                    A(p,q) = 0.0D0
                    A(q,p) = 0.0D0
                    cycle
                end if
                rotated = .true.
                app = A(p,p)
                aqq = A(q,q)
                ! tan(phi) of the rotation that annihilates A(p,q); smaller root
                theta = (aqq - app)/(2.0D0*apq)
                if (abs(theta) > JACOBI_HUGE_THETA) then
                    t = 0.5D0/theta
                else
                    t = sign(1.0D0, theta)/(abs(theta) + sqrt(theta*theta + 1.0D0))
                end if
                c = 1.0D0/sqrt(t*t + 1.0D0)
                s = t*c
                A(p,p) = app - t*apq
                A(q,q) = aqq + t*apq
                A(p,q) = 0.0D0
                A(q,p) = 0.0D0
                arp = A(r,p)
                arq = A(r,q)
                A(r,p) = c*arp - s*arq
                A(p,r) = A(r,p)
                A(r,q) = s*arp + c*arq
                A(q,r) = A(r,q)
                do i = 1, 3
                    tmp = V(i,p)
                    V(i,p) = c*tmp - s*V(i,q)
                    V(i,q) = s*tmp + c*V(i,q)
                end do
            end do
            if (.not. rotated) exit
        end do

        eigenvalues = [A(1,1), A(2,2), A(3,3)]

        ! Descending order (columns of V follow their eigenvalues)
        do i = 1, 2
            j = i - 1 + maxloc(eigenvalues(i:3), dim=1)
            if (j /= i) then
                tmp = eigenvalues(i)
                eigenvalues(i) = eigenvalues(j)
                eigenvalues(j) = tmp
                col = V(:,i)
                V(:,i) = V(:,j)
                V(:,j) = col
            end if
        end do
    end subroutine eigen_sym3


    pure function spectral_rotate(V, E) result(z)
        !! Components of a symmetric tensor in the eigenbasis: \(z = V^T E V\).
        !! The diagonal gives the first-order eigenvalue changes along `E`
        !! (for a cluster, of the chosen basis); off-diagonals feed the
        !! Daleckii-Krein rotation terms.
        real(real64), intent(in)      :: V(3,3) !! Orthonormal eigenvectors (columns)
        type(ten_3D2Osym), intent(in) :: E      !! Symmetric direction tensor
        real(real64)                  :: z(3,3)

        real(real64) :: Em(3,3), W(3,3)
        integer :: i, j, m

        Em(1,1) = E%vals(1)
        Em(2,2) = E%vals(2)
        Em(3,3) = E%vals(3)
        Em(1,2) = E%vals(4)
        Em(2,1) = E%vals(4)
        Em(2,3) = E%vals(5)
        Em(3,2) = E%vals(5)
        Em(1,3) = E%vals(6)
        Em(3,1) = E%vals(6)

        ! W = E V
        do j = 1, 3
            do i = 1, 3
                W(i,j) = Em(i,1)*V(1,j) + Em(i,2)*V(2,j) + Em(i,3)*V(3,j)
            end do
        end do
        ! z = V^T W, symmetrized against round-off
        do j = 1, 3
            do i = 1, 3
                z(i,j) = 0.0D0
                do m = 1, 3
                    z(i,j) = z(i,j) + V(m,i)*W(m,j)
                end do
            end do
        end do
        do j = 2, 3
            do i = 1, j - 1
                z(i,j) = 0.5D0*(z(i,j) + z(j,i))
                z(j,i) = z(i,j)
            end do
        end do
    end function spectral_rotate


    pure function spectral_compose(V, M) result(res)
        !! Symmetric tensor \(V M V^T\) from its eigenbasis components `M` (symmetric).
        !! With `M = diag(g)` this is the gradient \(\sum_i g_i\, v_i \otimes v_i\) of a
        !! spectral function; with the Daleckii-Krein matrix it is the Hessian action.
        real(real64), intent(in) :: V(3,3) !! Orthonormal eigenvectors (columns)
        real(real64), intent(in) :: M(3,3) !! Symmetric components in the eigenbasis
        type(ten_3D2Osym)        :: res

        real(real64) :: T(3,3)
        integer :: a, b, i, j
        integer, parameter :: IA(6) = [1, 2, 3, 1, 2, 1]
        integer, parameter :: IB(6) = [1, 2, 3, 2, 3, 3]

        ! T = V M
        do j = 1, 3
            do i = 1, 3
                T(i,j) = V(i,1)*M(1,j) + V(i,2)*M(2,j) + V(i,3)*M(3,j)
            end do
        end do
        ! res_ab = (V M V^T)_ab for the six stored components
        do i = 1, 6
            a = IA(i)
            b = IB(i)
            res%vals(i) = 0.0D0
            do j = 1, 3
                res%vals(i) = res%vals(i) + T(a,j)*V(b,j)
            end do
        end do
    end function spectral_compose

end module muscle_math_spectral_derivs