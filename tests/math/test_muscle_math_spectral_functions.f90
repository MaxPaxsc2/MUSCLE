program test_muscle_math_spectral_functions
    implicit none

    logical :: passed

    call test_spectral_trace_square_generic_1(passed)
    if (.not. passed) STOP 1

    call test_spectral_trace_square_repeated_2(passed)
    if (.not. passed) STOP 2

    call test_spectral_trace_cube_nearly_repeated_3(passed)
    if (.not. passed) STOP 3

    call test_spectral_coupled_repeated_4(passed)
    if (.not. passed) STOP 4

    print*, "Passed!", passed
    STOP 0
end program test_muscle_math_spectral_functions

subroutine rotated_tensor(lam, T)
    ! T = Q diag(lam) Q^T, with Q the rotation of 0.7 rad about the axis (1, 2, 3)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    real(real64), intent(in) :: lam(3)
    type(ten_3D2Osym), intent(out) :: T

    real(real64) :: k(3), Q(3,3), M(3,3), c, s
    integer :: i

    k = [1.0D0, 2.0D0, 3.0D0]/sqrt(14.0D0)
    c = cos(0.7D0)
    s = sin(0.7D0)
    Q = reshape([c + k(1)*k(1)*(1 - c), k(2)*k(1)*(1 - c) + k(3)*s, k(3)*k(1)*(1 - c) - k(2)*s, &
                 k(1)*k(2)*(1 - c) - k(3)*s, c + k(2)*k(2)*(1 - c), k(3)*k(2)*(1 - c) + k(1)*s, &
                 k(1)*k(3)*(1 - c) + k(2)*s, k(2)*k(3)*(1 - c) - k(1)*s, c + k(3)*k(3)*(1 - c)], [3,3])
    M = 0.0D0
    do i = 1, 3
        M(:,i) = lam(i)*Q(:,i)
    end do
    M = matmul(M, transpose(Q))
    call T%init(xx=M(1,1), yy=M(2,2), zz=M(3,3), xy=M(1,2), yz=M(2,3), xz=M(1,3))
end subroutine rotated_tensor

subroutine test_spectral_trace_square_generic_1(passed)
    ! phi = tr(T^2) = sum lam_a^2 in a generic rotated state: gradient 2T, Hessian 2 I4sym
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_math_spectral_functions
    implicit none
    logical, intent(out) :: passed

    real(real64), parameter :: TOL = 1.0D-12
    type(ten_3D2Osym) :: T, grad, grad_exp
    type(ten_3D4O3sym) :: hess, hess_exp
    real(real64) :: lam(3), V(3,3), d2phi(3,3)

    call rotated_tensor([1.3D0, 0.4D0, -0.9D0], T)
    call spectral_decomposition(T, lam, V)
    d2phi = 0.0D0
    d2phi(1,1) = 2.0D0
    d2phi(2,2) = 2.0D0
    d2phi(3,3) = 2.0D0
    grad = spectral_gradient(V, 2.0D0*lam)
    hess = spectral_hessian(lam, V, 2.0D0*lam, d2phi)

    grad_exp = 2.0D0*T
    hess_exp = 2.0D0*iden_4O4T()
    passed = grad%is_approx(grad_exp, tol=TOL) .and. hess%is_approx(hess_exp, tol=TOL)
    if (.not. passed) print*, "Spectral tr(T^2) generic failed:", grad%vals, hess%vals
end subroutine test_spectral_trace_square_generic_1

subroutine test_spectral_trace_square_repeated_2(passed)
    ! phi = tr(T^2) with an exactly repeated pair, diag(2,1,1), and triplet, diag(1,1,1):
    ! theta_ab takes its limit and the Hessian must still be exactly 2 I4sym
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_math_spectral_functions
    implicit none
    logical, intent(out) :: passed

    real(real64), parameter :: TOL = 1.0D-12
    type(ten_3D2Osym) :: T
    type(ten_3D4O3sym) :: hess, hess_exp
    real(real64) :: lam(3), V(3,3), d2phi(3,3)

    d2phi = 0.0D0
    d2phi(1,1) = 2.0D0
    d2phi(2,2) = 2.0D0
    d2phi(3,3) = 2.0D0
    hess_exp = 2.0D0*iden_4O4T()

    call T%init((/2D0, 1D0, 1D0, 0D0, 0D0, 0D0/))
    call spectral_decomposition(T, lam, V)
    hess = spectral_hessian(lam, V, 2.0D0*lam, d2phi)
    passed = hess%is_approx(hess_exp, tol=TOL)
    if (.not. passed) print*, "Spectral tr(T^2) diag(2,1,1) failed:", hess%vals
    if (.not. passed) return

    call T%init((/1D0, 1D0, 1D0, 0D0, 0D0, 0D0/))
    call spectral_decomposition(T, lam, V)
    hess = spectral_hessian(lam, V, 2.0D0*lam, d2phi)
    passed = hess%is_approx(hess_exp, tol=TOL)
    if (.not. passed) print*, "Spectral tr(T^2) diag(1,1,1) failed:", hess%vals
end subroutine test_spectral_trace_square_repeated_2

subroutine test_spectral_trace_cube_nearly_repeated_3(passed)
    ! phi = tr(T^3) = sum lam_a^3 in a rotated state with a nearly repeated pair (relative gap
    ! 1e-12, so theta_ab takes its limit, exact for a cubic): gradient 3 T^2 and Hessian
    ! H:D = 3 (T D + D T) = 3 ((T + D)^2 - T^2 - D^2).
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_math_spectral_functions
    implicit none
    logical, intent(out) :: passed

    real(real64), parameter :: TOL = 1.0D-11
    type(ten_3D2Osym) :: T, D, TD, grad, grad_exp, hd, hd_exp
    type(ten_3D4O3sym) :: hess
    real(real64) :: lam(3), V(3,3), d2phi(3,3)

    call rotated_tensor([1.0D0, 1.0D0 - 1.0D-12, -0.3D0], T)
    call spectral_decomposition(T, lam, V)
    d2phi = 0.0D0
    d2phi(1,1) = 6.0D0*lam(1)
    d2phi(2,2) = 6.0D0*lam(2)
    d2phi(3,3) = 6.0D0*lam(3)
    grad = spectral_gradient(V, 3.0D0*lam**2)
    hess = spectral_hessian(lam, V, 3.0D0*lam**2, d2phi)

    call D%init((/0.2D0, -0.5D0, 0.1D0, 0.3D0, -0.4D0, 0.6D0/))
    TD = T + D
    grad_exp = 3.0D0*T%square()
    hd = hess .ddot. D
    hd_exp = 3.0D0*(TD%square() - T%square() - D%square())
    passed = grad%is_approx(grad_exp, tol=TOL) .and. hd%is_approx(hd_exp, tol=TOL)
    if (.not. passed) print*, "Spectral tr(T^3) nearly repeated failed:", grad%vals, hd%vals
end subroutine test_spectral_trace_cube_nearly_repeated_3

subroutine test_spectral_coupled_repeated_4(passed)
    ! phi = tr(T) tr(T^2), not separable (phi_ab /= 0 for a /= b), in a rotated state with a
    ! repeated pair: gradient tr(T^2) I + 2 tr(T) T, Hessian 4 sym(T (x) I) + 2 tr(T) I4sym
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_math_spectral_functions
    implicit none
    logical, intent(out) :: passed

    real(real64), parameter :: TOL = 1.0D-12
    type(ten_3D2Osym) :: T, grad, grad_exp
    type(ten_3D4O3sym) :: hess, hess_exp
    real(real64) :: lam(3), V(3,3), d2phi(3,3), s1, s2
    integer :: a, b

    call rotated_tensor([0.8D0, -0.2D0, -0.2D0], T)
    call spectral_decomposition(T, lam, V)
    s1 = sum(lam)
    s2 = sum(lam**2)
    do a = 1, 3
        do b = 1, 3
            d2phi(a,b) = 2.0D0*(lam(a) + lam(b))
        end do
        d2phi(a,a) = d2phi(a,a) + 2.0D0*s1
    end do
    grad = spectral_gradient(V, s2 + 2.0D0*s1*lam)
    hess = spectral_hessian(lam, V, s2 + 2.0D0*s1*lam, d2phi)

    s1 = T%xx() + T%yy() + T%zz()
    grad_exp = (T .ddot. T)*iden_2O() + (2.0D0*s1)*T
    hess_exp = 4.0D0*(T .tdotsym. iden_2O()) + (2.0D0*s1)*iden_4O4T()
    passed = grad%is_approx(grad_exp, tol=TOL) .and. hess%is_approx(hess_exp, tol=TOL)
    if (.not. passed) print*, "Spectral tr(T) tr(T^2) repeated failed:", grad%vals, hess%vals
end subroutine test_spectral_coupled_repeated_4
