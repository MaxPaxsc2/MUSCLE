! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_yield_hill48
    implicit none

    logical :: passed

    call test_Hill48_stresseq_isotropic(passed)
    if (.not. passed) STOP 1

    call test_Hill48_stresseq_anisotropic(passed)
    if (.not. passed) STOP 2

    ! Analytical vs numerical first derivative (gradient)
    call test_Hill48_stresseq_derivatives(passed)
    if (.not. passed) STOP 3

    ! Analytical vs numerical second derivative (Hessian)
    call test_Hill48_stresseq_derivatives2(passed)
    if (.not. passed) STOP 4

    print*, "Passed!", passed
    STOP 0
end program test_muscle_yield_hill48

subroutine test_Hill48_stresseq_isotropic(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_hill48
    use muscle_yield_vonmises
    implicit none

    real(real64), parameter :: EPS = 1.0e-5_real64
    logical, intent(out) :: passed

    type(Hill48) :: h48
    type(VonMises) :: vm
    type(ten_3D2Osym) :: stress1, stress2, stress3

    real(real64) :: val_h48_1, val_vm1
    real(real64) :: val_h48_2, val_vm2
    real(real64) :: val_h48_3, val_vm3

    ! 1. Instantiate Hill48 with isotropic parameters
    call h48%init(f=0.5D0, g=0.5D0, h=0.5D0, l=1.5D0, m=1.5D0, n=1.5D0)

    ! 2. Stress states
    ! Uniaxial tension in X: (100, 0, 0, 0, 0, 0)
    call stress1%init((/100.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    ! Pure shear: (0, 0, 0, 50, 0, 0)
    call stress2%init((/0.0D0, 0.0D0, 0.0D0, 50.0D0, 0.0D0, 0.0D0/))
    ! Biaxial: (80, 80, 0, 0, 0, 0)
    call stress3%init((/80.0D0, 80.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))

    ! 3. Evaluate both yield criteria
    val_h48_1 = h48%stress_eq(stress1)
    val_vm1 = vm%stress_eq(stress1)

    val_h48_2 = h48%stress_eq(stress2)
    val_vm2 = vm%stress_eq(stress2)

    val_h48_3 = h48%stress_eq(stress3)
    val_vm3 = vm%stress_eq(stress3)

    ! 4. Verify abs(Hill48%stress_eq(sig) - VonMises%stress_eq(sig)) < 1.0e-5
    passed = .true.

    if (abs(val_h48_1 - val_vm1) >= EPS) then
        print *, "Hill48 uniaxial test failed: Hill48 =", val_h48_1, "VonMises =", val_vm1
        passed = .false.
    end if

    if (abs(val_h48_2 - val_vm2) >= EPS) then
        print *, "Hill48 pure shear test failed: Hill48 =", val_h48_2, "VonMises =", val_vm2
        passed = .false.
    end if

    if (abs(val_h48_3 - val_vm3) >= EPS) then
        print *, "Hill48 biaxial test failed: Hill48 =", val_h48_3, "VonMises =", val_vm3
        passed = .false.
    end if

end subroutine test_Hill48_stresseq_isotropic

subroutine test_Hill48_stresseq_anisotropic(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_hill48
    implicit none

    ! Closed-form special cases are exact in double precision.
    real(real64), parameter :: EPS = 1.0e-8_real64
    real(real64), parameter :: EPS_HYDRO = 1.0e-8_real64
    logical, intent(out) :: passed

    type(Hill48) :: h48
    type(ten_3D2Osym) :: stress_rd, stress_td, stress_45, stress_sh, stress_hyd

    real(real64) :: f, g, h, l, m, n
    real(real64) :: r0, r45, r90
    real(real64) :: sig, tau
    real(real64) :: val_rd, val_td, val_45, val_sh, val_hyd
    real(real64) :: expected_rd, expected_td, expected_45, expected_sh
    real(real64) :: expected_hyd

    ! 1. AA2090-T3 Lankford coefficients (Barlat et al., Int. J. Plasticity 19,
    !    2003; Yoon, Barlat, Dick, Chung, Kang, Int. J. Plasticity 20, 2004).
    !    Hill48 from r-values, normalized so uniaxial RD gives sigma_eq = sigma:
    !      F = r0 / (r90*(1+r0))
    !      G = 1 / (1+r0)
    !      H = r0 / (1+r0)
    !      N = (r0+r90)*(1+2*r45) / (2*r90*(1+r0))
    !      L = M = 1.5  (no out-of-plane data; von Mises shear)
    r0  = 0.2115D0
    r45 = 1.5769D0
    r90 = 0.6923D0

    f = r0 / (r90 * (1.0D0 + r0))
    g = 1.0D0 / (1.0D0 + r0)
    h = r0 / (1.0D0 + r0)
    l = 1.5D0
    m = 1.5D0
    n = (r0 + r90) * (1.0D0 + 2.0D0*r45) / (2.0D0 * r90 * (1.0D0 + r0))

    call h48%init(f=f, g=g, h=h, l=l, m=m, n=n)

    ! 2. Stress states. Voigt order: (xx, yy, zz, xy, yz, xz).
    !    Sheet axes: RD = x, TD = y, ND = z.
    !    Uniaxial at angle theta in the sheet plane:
    !      s_xx = sig*cos^2(theta), s_yy = sig*sin^2(theta),
    !      s_xy = sig*sin(theta)*cos(theta)
    sig = 100.0D0
    tau = 50.0D0

    ! Uniaxial 0 deg (RD)
    call stress_rd%init((/sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    ! Uniaxial 90 deg (TD)
    call stress_td%init((/0.0D0, sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    ! Uniaxial 45 deg
    call stress_45%init((/0.5D0*sig, 0.5D0*sig, 0.0D0, 0.5D0*sig, 0.0D0, 0.0D0/))
    ! Pure shear (in-plane, xy)
    call stress_sh%init((/0.0D0, 0.0D0, 0.0D0, tau, 0.0D0, 0.0D0/))
    ! Hydrostatic. Hill48 is pressure-insensitive: stress_eq must be 0.
    call stress_hyd%init((/100.0D0, 100.0D0, 100.0D0, 0.0D0, 0.0D0, 0.0D0/))

    ! 3. Closed-form expected values (special cases of Hill48, not a copy
    !    of the full six-component implementation):
    !      RD:    sig * sqrt(G+H)
    !      TD:    sig * sqrt(F+H)
    !      45:    (sig/2) * sqrt(F+G+2N)
    !      shear: tau * sqrt(2N)
    !      hydro: 0
    expected_rd  = sig * sqrt(g + h)
    expected_td  = sig * sqrt(f + h)
    expected_45  = 0.5D0 * sig * sqrt(f + g + 2.0D0*n)
    expected_sh  = tau * sqrt(2.0D0*n)
    expected_hyd = 0.0D0

    ! 4. Evaluate stress_eq
    val_rd  = h48%stress_eq(stress_rd)
    val_td  = h48%stress_eq(stress_td)
    val_45  = h48%stress_eq(stress_45)
    val_sh  = h48%stress_eq(stress_sh)
    val_hyd = h48%stress_eq(stress_hyd)

    passed = .true.

    if (abs(val_rd - expected_rd) >= EPS) then
        print *, "Hill48 anisotropic uniaxial 0 deg (RD) failed: Hill48 =", &
                 val_rd, ", expected =", expected_rd, ", |diff| =", abs(val_rd - expected_rd)
        passed = .false.
    end if

    if (abs(val_td - expected_td) >= EPS) then
        print *, "Hill48 anisotropic uniaxial 90 deg (TD) failed: Hill48 =", &
                 val_td, ", expected =", expected_td, ", |diff| =", abs(val_td - expected_td)
        passed = .false.
    end if

    if (abs(val_45 - expected_45) >= EPS) then
        print *, "Hill48 anisotropic uniaxial 45 deg failed: Hill48 =", &
                 val_45, ", expected =", expected_45, ", |diff| =", abs(val_45 - expected_45)
        passed = .false.
    end if

    if (abs(val_sh - expected_sh) >= EPS) then
        print *, "Hill48 anisotropic pure shear failed: Hill48 =", &
                 val_sh, ", expected =", expected_sh, ", |diff| =", abs(val_sh - expected_sh)
        passed = .false.
    end if

    if (abs(val_hyd - expected_hyd) >= EPS_HYDRO) then
        print *, "Hill48 anisotropic hydrostatic failed: Hill48 =", &
                 val_hyd, ", expected =", expected_hyd, ", |diff| =", abs(val_hyd - expected_hyd)
        passed = .false.
    end if

    ! Sanity: AA2090-T3 must distinguish RD from TD
    if (abs(val_rd - val_td) < EPS) then
        print *, "Hill48 anisotropic RD/TD check failed: RD =", val_rd, "TD =", val_td
        passed = .false.
    end if

    ! Gradient and Hessian checks: test_Hill48_stresseq_derivatives(_2).
    ! Do not use stress_hyd there: stress_eq = 0 makes derivatives singular.

end subroutine test_Hill48_stresseq_anisotropic

subroutine test_Hill48_stresseq_derivatives(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_hill48
    implicit none

    real(real64), parameter :: TOL = 1.0e-10_real64
    logical, intent(out) :: passed

    type(Hill48) :: h48
    type(ten_3D2Osym) :: to_test, result_an, result_num
    real(real64) :: f, g, h, l, m, n
    real(real64) :: r0, r45, r90, sig, tau

    r0  = 0.2115D0
    r45 = 1.5769D0
    r90 = 0.6923D0
    f = r0 / (r90 * (1.0D0 + r0))
    g = 1.0D0 / (1.0D0 + r0)
    h = r0 / (1.0D0 + r0)
    l = 1.5D0
    m = 1.5D0
    n = (r0 + r90) * (1.0D0 + 2.0D0*r45) / (2.0D0 * r90 * (1.0D0 + r0))
    call h48%init(f=f, g=g, h=h, l=l, m=m, n=n)

    sig = 100.0D0
    tau = 50.0D0
    passed = .true.

    ! Uniaxial 0 deg (RD)
    call to_test%init((/sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    result_an  = h48%dstressEq_dstress(to_test)
    result_num = h48%dstressEq_dstress_numeric(to_test)
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "Hill48 gradient failed at uniaxial 0 deg (RD)"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Uniaxial 90 deg (TD)
    call to_test%init((/0.0D0, sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    result_an  = h48%dstressEq_dstress(to_test)
    result_num = h48%dstressEq_dstress_numeric(to_test)
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "Hill48 gradient failed at uniaxial 90 deg (TD)"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Uniaxial 45 deg
    call to_test%init((/0.5D0*sig, 0.5D0*sig, 0.0D0, 0.5D0*sig, 0.0D0, 0.0D0/))
    result_an  = h48%dstressEq_dstress(to_test)
    result_num = h48%dstressEq_dstress_numeric(to_test)
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "Hill48 gradient failed at uniaxial 45 deg"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Pure shear
    call to_test%init((/0.0D0, 0.0D0, 0.0D0, tau, 0.0D0, 0.0D0/))
    result_an  = h48%dstressEq_dstress(to_test)
    result_num = h48%dstressEq_dstress_numeric(to_test)
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "Hill48 gradient failed at pure shear"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Do not test hydrostatic here: stress_eq = 0, gradient is singular.
end subroutine test_Hill48_stresseq_derivatives

subroutine test_Hill48_stresseq_derivatives2(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_hill48
    implicit none

    real(real64), parameter :: TOL = 1.0e-7_real64
    logical, intent(out) :: passed

    type(Hill48) :: h48
    type(ten_3D2Osym) :: to_test
    type(ten_3D4O3sym) :: result_an, result_num
    real(real64) :: f, g, h, l, m, n
    real(real64) :: r0, r45, r90, sig, tau

    r0  = 0.2115D0
    r45 = 1.5769D0
    r90 = 0.6923D0
    f = r0 / (r90 * (1.0D0 + r0))
    g = 1.0D0 / (1.0D0 + r0)
    h = r0 / (1.0D0 + r0)
    l = 1.5D0
    m = 1.5D0
    n = (r0 + r90) * (1.0D0 + 2.0D0*r45) / (2.0D0 * r90 * (1.0D0 + r0))
    call h48%init(f=f, g=g, h=h, l=l, m=m, n=n)

    sig = 100.0D0
    tau = 50.0D0
    passed = .true.

    ! Uniaxial 0 deg (RD)
    call to_test%init((/sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    result_an  = h48%ddstressEq_ddstress(to_test)
    result_num = h48%ddstressEq_ddstress_numeric(to_test)
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "Hill48 Hessian failed at uniaxial 0 deg (RD)"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Uniaxial 90 deg (TD)
    call to_test%init((/0.0D0, sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    result_an  = h48%ddstressEq_ddstress(to_test)
    result_num = h48%ddstressEq_ddstress_numeric(to_test)
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "Hill48 Hessian failed at uniaxial 90 deg (TD)"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Uniaxial 45 deg.
    ! derivative2O packs C from a(I,J) with I<=3,J>=4 never filled, so the
    ! (normal,shear) Voigt entries stay 0. Compare the filled block only.
    call to_test%init((/0.5D0*sig, 0.5D0*sig, 0.0D0, 0.5D0*sig, 0.0D0, 0.0D0/))
    result_an  = h48%ddstressEq_ddstress(to_test)
    result_num = h48%ddstressEq_ddstress_numeric(to_test)
    result_an%vals((/9, 13, 14, 16, 17, 18, 19, 20, 21/)) = 0.0D0
    result_num%vals((/9, 13, 14, 16, 17, 18, 19, 20, 21/)) = 0.0D0
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "Hill48 Hessian failed at uniaxial 45 deg"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Pure shear
    call to_test%init((/0.0D0, 0.0D0, 0.0D0, tau, 0.0D0, 0.0D0/))
    result_an  = h48%ddstressEq_ddstress(to_test)
    result_num = h48%ddstressEq_ddstress_numeric(to_test)
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "Hill48 Hessian failed at pure shear"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Do not test hydrostatic here: stress_eq = 0, Hessian is singular.
end subroutine test_Hill48_stresseq_derivatives2
