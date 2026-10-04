! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_math_spectral_functions
    !! Module muscle_math_spectral_functions
    !! =====================================
    !! Gradient and Hessian, with respect to a symmetric tensor T, of a symmetric function of
    !! its eigenvalues, phi(T) = phi(lam_1, lam_2, lam_3). The caller supplies the scalar
    !! derivatives phi_a = dphi/dlam_a and phi_ab = d2phi/dlam_a dlam_b
    !! (de Souza Neto et al., 2008, Appendix A; Miehe, 1998):
    !! \[ \frac{\partial \phi}{\partial T} = \sum_a \phi_a E_a, \qquad
    !!    \frac{\partial^2 \phi}{\partial T^2} = \sum_{a,b} \phi_{ab} E_a \otimes E_b
    !!    + \sum_{a<b} 2\theta_{ab} N_{ab} \otimes N_{ab} \]
    !! with \(E_a = v_a \otimes v_a\), \(N_{ab} = \mathrm{sym}(v_a \otimes v_b)\) and
    !! \(\theta_{ab} = (\phi_a - \phi_b)/(\lambda_a - \lambda_b)\). For (nearly) repeated
    !! eigenvalues \(\theta_{ab}\) takes its limit \(\frac{1}{2}(\phi_{aa} + \phi_{bb}) - \phi_{ab}\).
    !! The eigenvalues come from `eigenvals`; one decomposition per call serves both derivatives.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    implicit none
    private

    public :: spectral_decomposition
    public :: spectral_gradient
    public :: spectral_hessian

    real(real64), parameter :: REPEATED_TOL = 1.0D-9
    !! Relative eigenvalue gap below which theta_ab takes its limit value. The limit assumes phi
    !! is twice differentiable between the two eigenvalues: 1e-9 keeps that band narrow where it
    !! is not (CPB06, a = 2, k /= 0: h'' jumps at lam = 0) and the quotient error below 1e-7.

    ! Voigt index pairs (i, j) of the components of ten_3D2Osym
    integer, parameter :: VI(6) = [1, 2, 3, 1, 2, 1]
    integer, parameter :: VJ(6) = [1, 2, 3, 2, 3, 3]

contains

    pure subroutine spectral_decomposition(T, lam, V)
        !! Eigenvalues and orthonormal eigenvectors of T. The eigenvector of the eigenvalue
        !! farthest from the middle one is the normalized cross product of two rows of
        !! T - lam I (Kopp, 2008); the other two come from one Jacobi rotation in the plane
        !! normal to it, which stays accurate when they (nearly) coincide. That eigenvalue is at
        !! least half the spread away from the other two, so the cross product needs no fallback.
        !! Called by the CPB06 derivatives.
        use muscle_math_operations, only : eigenvals
        type(ten_3D2Osym), intent(in) :: T   !! Symmetric tensor
        real(real64), intent(out) :: lam(3)  !! Eigenvalues, in descending order
        real(real64), intent(out) :: V(3,3)  !! Eigenvectors: column a belongs to lam(a)

        real(real64) :: rows(3,3), c(3,3), n(3), u(3), w(3), tu(3), tw(3)
        real(real64) :: p, q, r, theta, tg, cs, sn
        integer :: a, b, i

        lam = eigenvals(T)

        ! a: most separated eigenvalue; b, b+1: the remaining pair
        if (lam(1) - lam(2) >= lam(2) - lam(3)) then
            a = 1
            b = 2
        else
            a = 3
            b = 1
        end if

        rows(:,1) = [T%vals(1) - lam(a), T%vals(4), T%vals(6)]
        rows(:,2) = [T%vals(4), T%vals(2) - lam(a), T%vals(5)]
        rows(:,3) = [T%vals(6), T%vals(5), T%vals(3) - lam(a)]
        c(:,1) = cross(rows(:,1), rows(:,2))
        c(:,2) = cross(rows(:,2), rows(:,3))
        c(:,3) = cross(rows(:,3), rows(:,1))
        n = sum(c**2, dim=1)
        i = maxloc(n, dim=1)
        if (n(i) <= 0.0D0) then
            ! T is a multiple of the identity: any orthonormal basis is an eigenbasis
            V = 0.0D0
            V(1,1) = 1.0D0
            V(2,2) = 1.0D0
            V(3,3) = 1.0D0
            return
        end if
        V(:,a) = c(:,i)/sqrt(n(i))

        ! Orthonormal basis (u, w) of the plane normal to V(:,a)
        if (abs(V(1,a)) > abs(V(2,a))) then
            u = [-V(3,a), 0.0D0, V(1,a)]/sqrt(V(1,a)**2 + V(3,a)**2)
        else
            u = [0.0D0, V(3,a), -V(2,a)]/sqrt(V(2,a)**2 + V(3,a)**2)
        end if
        w = cross(V(:,a), u)

        ! Jacobi rotation that diagonalizes [[p, q], [q, r]], the restriction of T to (u, w)
        tu = [T%vals(1)*u(1) + T%vals(4)*u(2) + T%vals(6)*u(3), &
              T%vals(4)*u(1) + T%vals(2)*u(2) + T%vals(5)*u(3), &
              T%vals(6)*u(1) + T%vals(5)*u(2) + T%vals(3)*u(3)]
        tw = [T%vals(1)*w(1) + T%vals(4)*w(2) + T%vals(6)*w(3), &
              T%vals(4)*w(1) + T%vals(2)*w(2) + T%vals(5)*w(3), &
              T%vals(6)*w(1) + T%vals(5)*w(2) + T%vals(3)*w(3)]
        p = dot_product(u, tu)
        q = dot_product(u, tw)
        r = dot_product(w, tw)
        tg = 0.0D0
        if (abs(q) > 0.0D0) then
            theta = (r - p)/(2.0D0*q)
            tg = sign(1.0D0, theta)/(abs(theta) + sqrt(theta*theta + 1.0D0))
        end if
        cs = 1.0D0/sqrt(tg*tg + 1.0D0)
        sn = tg*cs
        ! Rotated eigenvalues: p - tg*q for cs*u - sn*w, and r + tg*q for sn*u + cs*w
        if (p - tg*q >= r + tg*q) then
            V(:,b) = cs*u - sn*w
            V(:,b+1) = sn*u + cs*w
        else
            V(:,b) = sn*u + cs*w
            V(:,b+1) = cs*u - sn*w
        end if
    end subroutine spectral_decomposition

    pure function spectral_gradient(V, dphi) result(res)
        !! dphi/dT = sum_a phi_a v_a (x) v_a (de Souza Neto et al., 2008, Appendix A).
        !! Called by the CPB06 gradient.
        real(real64), intent(in) :: V(3,3)   !! Eigenvectors from spectral_decomposition
        real(real64), intent(in) :: dphi(3)  !! phi_a = dphi/dlam_a
        type(ten_3D2Osym) :: res

        integer :: q

        do q = 1, 6
            res%vals(q) = sum(dphi*V(VI(q),:)*V(VJ(q),:))
        end do
    end function spectral_gradient

    pure function spectral_hessian(lam, V, dphi, d2phi) result(res)
        !! d2phi/dT2 = sum_ab phi_ab E_a (x) E_b + sum_(a<b) 2 theta_ab N_ab (x) N_ab
        !! (de Souza Neto et al., 2008, Appendix A; Miehe, 1998). Called by the CPB06 Hessian.
        real(real64), intent(in) :: lam(3)      !! Eigenvalues from spectral_decomposition
        real(real64), intent(in) :: V(3,3)      !! Eigenvectors from spectral_decomposition
        real(real64), intent(in) :: dphi(3)     !! phi_a = dphi/dlam_a
        real(real64), intent(in) :: d2phi(3,3)  !! phi_ab = d2phi/dlam_a dlam_b
        type(ten_3D4O3sym) :: res

        type(ten_3D2Osym) :: E(3), N
        real(real64) :: tol, theta
        integer :: a, b, q

        do a = 1, 3
            do q = 1, 6
                E(a)%vals(q) = V(VI(q),a)*V(VJ(q),a)
            end do
        end do
        tol = REPEATED_TOL*max(abs(lam(1)), abs(lam(3)))

        res = 0.0D0
        do a = 1, 3
            res = res + d2phi(a,a)*(.tdotsym. E(a))
            do b = a + 1, 3
                if (abs(lam(a) - lam(b)) > tol) then
                    theta = (dphi(a) - dphi(b))/(lam(a) - lam(b))
                else
                    theta = 0.5D0*(d2phi(a,a) + d2phi(b,b)) - d2phi(a,b)
                end if
                do q = 1, 6
                    N%vals(q) = 0.5D0*(V(VI(q),a)*V(VJ(q),b) + V(VI(q),b)*V(VJ(q),a))
                end do
                res = res + (2.0D0*d2phi(a,b))*(E(a) .tdotsym. E(b)) + (2.0D0*theta)*(.tdotsym. N)
            end do
        end do
    end function spectral_hessian

    pure function cross(x, y) result(res)
        ! Cross product x * y; used by spectral_decomposition
        real(real64), intent(in) :: x(3), y(3)
        real(real64) :: res(3)

        res = [x(2)*y(3) - x(3)*y(2), x(3)*y(1) - x(1)*y(3), x(1)*y(2) - x(2)*y(1)]
    end function cross

end module muscle_math_spectral_functions
