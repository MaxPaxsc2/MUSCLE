program test_muscle_yield_linear_transform
    implicit none

    logical :: passed

    call test_linear_transform_deviator_1(passed)
    if (.not. passed) STOP 1

    call test_linear_transform_pressure_2(passed)
    if (.not. passed) STOP 2

    call test_linear_transform_pull_back_3(passed)
    if (.not. passed) STOP 3

    call test_linear_transform_hessian_4(passed)
    if (.not. passed) STOP 4

    print*, "Passed!", passed
    STOP 0
end program test_muscle_yield_linear_transform

subroutine linear_transform_test_data(C, c_shear, stress)
    ! Non-symmetric normal block with a nonzero diagonal, three different shear factors and a
    ! stress with normal and shear components together
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    real(real64), intent(out) :: C(3,3), c_shear(3)
    type(ten_3D2Osym), intent(out) :: stress

    C(1,1) = 1.1D0
    C(1,2) = -0.4D0
    C(1,3) = 0.3D0
    C(2,1) = 0.2D0
    C(2,2) = 0.9D0
    C(2,3) = -0.6D0
    C(3,1) = -0.5D0
    C(3,2) = 0.7D0
    C(3,3) = 1.3D0
    c_shear(1) = 1.2D0
    c_shear(2) = 0.8D0
    c_shear(3) = 1.5D0
    call stress%init((/120D0, -40D0, 75D0, 30D0, -55D0, 20D0/))
end subroutine linear_transform_test_data

subroutine test_linear_transform_deviator_1(passed)
    ! S = L : stress with L = deviatoric_block(C) against C applied to .dev. stress, and with
    ! C = I (isotropic limit) against .dev. stress itself
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_linear_transform
    implicit none
    logical, intent(out) :: passed

    real(real64), parameter :: TOL = 1.0D-14
    real(real64) :: C(3,3), c_shear(3), one(3)
    type(ten_3D2Osym) :: stress, s, S_ref, S_new

    call linear_transform_test_data(C, c_shear, stress)
    s = .dev. stress
    call S_ref%init(xx=C(1,1)*s%xx() + C(1,2)*s%yy() + C(1,3)*s%zz(), &
                    yy=C(2,1)*s%xx() + C(2,2)*s%yy() + C(2,3)*s%zz(), &
                    zz=C(3,1)*s%xx() + C(3,2)*s%yy() + C(3,3)*s%zz(), &
                    xy=c_shear(1)*s%xy(), yz=c_shear(2)*s%yz(), xz=c_shear(3)*s%xz())
    S_new = linear_transform(deviatoric_block(C), c_shear, stress)
    passed = S_new%is_approx(S_ref, tol=TOL)
    if (.not. passed) print*, "Linear transform, general C failed:", S_new%vals, S_ref%vals
    if (.not. passed) return

    C = 0.0D0
    C(1,1) = 1.0D0
    C(2,2) = 1.0D0
    C(3,3) = 1.0D0
    one = 1.0D0
    S_new = linear_transform(deviatoric_block(C), one, stress)
    passed = S_new%is_approx(s, tol=TOL)
    if (.not. passed) print*, "Linear transform, C = I failed:", S_new%vals, s%vals
end subroutine test_linear_transform_deviator_1

subroutine test_linear_transform_pressure_2(passed)
    ! L : (p I + stress) against L : stress with a pressure p = 2**40, about 1e10 times the
    ! stress: L annihilates p I, so the pressure must not leave rounding of order 1e-16 p in S.
    ! The stress has integer components, so p I + stress is exact in real64
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_linear_transform
    implicit none
    logical, intent(out) :: passed

    real(real64), parameter :: TOL = 1.0D-14
    real(real64) :: C(3,3), c_shear(3), L(3,3), p
    type(ten_3D2Osym) :: stress, loaded, S_ref, S_new

    call linear_transform_test_data(C, c_shear, stress)
    L = deviatoric_block(C)
    p = 2.0D0**40
    call loaded%init(xx=p + stress%xx(), yy=p + stress%yy(), zz=p + stress%zz(), &
                     xy=stress%xy(), yz=stress%yz(), xz=stress%xz())
    S_ref = linear_transform(L, c_shear, stress)
    S_new = linear_transform(L, c_shear, loaded)
    passed = S_new%is_approx(S_ref, tol=TOL)
    if (.not. passed) print*, "Linear transform with a large pressure failed:", S_new%vals, &
                              S_ref%vals
end subroutine test_linear_transform_pressure_2

subroutine test_linear_transform_pull_back_3(passed)
    ! Adjoint identity (L^T : x) : y = x : (L : y) with x and y full tensors: every entry of L^T
    ! must be the transpose of the one used by linear_transform
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_linear_transform
    implicit none
    logical, intent(out) :: passed

    real(real64), parameter :: TOL = 1.0D-14
    real(real64) :: C(3,3), c_shear(3), L(3,3), lhs, rhs
    type(ten_3D2Osym) :: x, y

    call linear_transform_test_data(C, c_shear, x)
    call y%init((/-0.3D0, 0.8D0, 0.25D0, -0.6D0, 0.45D0, 0.9D0/))
    L = deviatoric_block(C)
    lhs = pull_back(L, c_shear, x) .ddot. y
    rhs = x .ddot. linear_transform(L, c_shear, y)
    passed = abs(lhs - rhs) < TOL*abs(rhs)
    if (.not. passed) print*, "Pull back adjoint failed:", lhs, rhs
end subroutine test_linear_transform_pull_back_3

subroutine test_linear_transform_hessian_4(passed)
    ! (L^T : H : L) : y = L^T : (H : (L : y)) for the six unit tensors y of the Voigt basis,
    ! with H having its 21 components different and nonzero: each y checks one column, so
    ! every normal, normal-shear and shear component of the result is checked
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_linear_transform
    implicit none
    logical, intent(out) :: passed

    real(real64), parameter :: TOL = 1.0D-14
    real(real64) :: C(3,3), c_shear(3), L(3,3)
    type(ten_3D2Osym) :: stress, y, hy, hy_ref
    type(ten_3D4O3sym) :: H
    integer :: i

    call linear_transform_test_data(C, c_shear, stress)
    L = deviatoric_block(C)
    call H%init((/2.1D0, 1.7D0, 2.6D0, 0.9D0, 1.3D0, 0.7D0, -0.4D0, 0.35D0, 0.12D0, -0.22D0, &
                  0.18D0, -0.6D0, -0.27D0, 0.41D0, 0.08D0, 0.33D0, -0.15D0, 0.5D0, -0.38D0, &
                  0.24D0, -0.11D0/))
    do i = 1, 6
        y%vals = 0.0D0
        y%vals(i) = 1.0D0
        hy = pull_back_hessian(L, c_shear, H) .ddot. y
        hy_ref = pull_back(L, c_shear, H .ddot. linear_transform(L, c_shear, y))
        passed = hy%is_approx(hy_ref, tol=TOL)
        if (.not. passed) print*, "Pull back of the Hessian failed, direction", i, &
                                  hy%vals, hy_ref%vals
        if (.not. passed) return
    end do
end subroutine test_linear_transform_hessian_4
