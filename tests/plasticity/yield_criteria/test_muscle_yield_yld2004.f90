! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_yield_yld2004
    implicit none

    logical :: passed

    call test_yld2004_stresseq_isotropic(passed)
    if (.not. passed) STOP 1

    call test_yld2004_stresseq_anisotropic(passed)
    if (.not. passed) STOP 2

    ! Analytical vs numerical first derivative (gradient)
    call test_yld2004_stresseq_derivatives(passed)
    if (.not. passed) STOP 3

    ! Analytical vs numerical second derivative (Hessian)
    call test_yld2004_stresseq_derivatives2(passed)
    if (.not. passed) STOP 4

    ! Regression: isotropic limit vs closed-form J2, full Hessian, repeated
    ! spectra and stress magnitudes from 1e-6 to 2.5e8
    call test_yld2004_isotropic_j2_full(passed)
    if (.not. passed) STOP 5

    ! Regression: even a, AA2090-T3 and a normal-shear coupled map vs the
    ! tr(K**a) reference at generic, repeated and nearly repeated spectra
    call test_yld2004_kron_reference(passed)
    if (.not. passed) STOP 6

    ! Regression: scale invariance, homogeneity identities and absence of IEEE
    ! exceptions at states that used to divide by zero
    call test_yld2004_scale_and_fpe(passed)
    if (.not. passed) STOP 7

    print*, "Passed!", passed
    STOP 0
end program test_muscle_yield_yld2004

subroutine test_yld2004_stresseq_isotropic(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_vonmises
    use muscle_yield_yld2004
    implicit none

    real(real64), parameter :: EPS = 1.0e-5_real64
    logical, intent(out) :: passed

    type(Yld2004) :: yld
    type(VonMises) :: vm
    type(ten_3D2Osym) :: stress1, stress2, stress3
    real(real64), dimension(6,6) :: C_prime, C_dprime
    real(real64) :: a
    real(real64) :: val_yld1, val_vm1
    real(real64) :: val_yld2, val_vm2
    real(real64) :: val_yld3, val_vm3
    integer :: i

    ! 1. Instantiate Yld2004 with isotropic parameters: C_prime = I, C_dprime = I, a = 2.0
    C_prime = 0.0D0
    C_dprime = 0.0D0
    do i = 1, 6
        C_prime(i,i) = 1.0D0
        C_dprime(i,i) = 1.0D0
    end do
    a = 2.0D0

    call yld%init(C_prime, C_dprime, a)

    ! 2. Instantiate stress states
    ! Uniaxial tension in X: (100, 0, 0, 0, 0, 0)
    call stress1%init((/100.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    ! Pure shear: (0, 0, 0, 50, 0, 0)
    call stress2%init((/0.0D0, 0.0D0, 0.0D0, 50.0D0, 0.0D0, 0.0D0/))
    ! Biaxial: (80, 80, 0, 0, 0, 0)
    call stress3%init((/80.0D0, 80.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))

    ! 3. Evaluate both yield criteria
    val_yld1 = yld%stress_eq(stress1)
    val_vm1 = vm%stress_eq(stress1)

    val_yld2 = yld%stress_eq(stress2)
    val_vm2 = vm%stress_eq(stress2)

    val_yld3 = yld%stress_eq(stress3)
    val_vm3 = vm%stress_eq(stress3)

    ! 4. Verify abs(Yld2004%stress_eq(sig) - VonMises%stress_eq(sig)) < 1.0e-5
    passed = .true.

    if (abs(val_yld1 - val_vm1) >= EPS) then
        print *, "Uniaxial test failed: Yld2004 =", val_yld1, "VonMises =", val_vm1
        passed = .false.
    end if

    if (abs(val_yld2 - val_vm2) >= EPS) then
        print *, "Pure shear test failed: Yld2004 =", val_yld2, "VonMises =", val_vm2
        passed = .false.
    end if

    if (abs(val_yld3 - val_vm3) >= EPS) then
        print *, "Biaxial test failed: Yld2004 =", val_yld3, "VonMises =", val_vm3
        passed = .false.
    end if

end subroutine test_yld2004_stresseq_isotropic

subroutine test_yld2004_stresseq_anisotropic(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_yld2004
    implicit none

    real(real64), parameter :: EPS = 1.0e-8_real64
    real(real64), parameter :: EPS_HYDRO = 1.0e-8_real64
    logical, intent(out) :: passed

    type(Yld2004) :: yld
    type(ten_3D2Osym) :: stress_rd, stress_45, stress_td, stress_sh, stress_hyd
    real(real64), dimension(6,6) :: C_prime, C_dprime
    real(real64) :: a
    real(real64) :: cp12, cp13, cp21, cp23, cp31, cp32, cp44, cp55, cp66
    real(real64) :: cd12, cd13, cd21, cd23, cd31, cd32, cd44, cd55, cd66
    real(real64) :: sig, tau
    real(real64) :: sxx, syy, szz, sxy
    real(real64) :: sp1, sp2, sp3, sd1, sd2, sd3
    real(real64) :: mid, rad
    real(real64) :: val_rd, val_45, val_td, val_sh, val_hyd
    real(real64) :: expected_rd, expected_45, expected_td, expected_sh, expected_hyd

    ! 1. AA2090-T3 Yld2004-18p, a = 8 (Barlat, Aretz, Yoon, Karabin, Brem,
    !    Dick, Int. J. Plasticity 21, 2005, Table 3).
    !    9 + 9 coefficients of the two linear maps on the stress deviator.
    !    Barlat Voigt shear order (yz, xz, xy) -> MUSCLE (xy, yz, xz):
    !      C(4,4)=c66, C(5,5)=c44, C(6,6)=c55.
    !    Shear values in Barlat order (c44 = yz, c55 = zx, c66 = xy). They were
    !    previously assigned cyclically permuted; with this assignment the
    !    in-plane AA2090-T3 yield ratios are reproduced within 1.2 %
    !    (6.2 % with the permuted values). Keep the three copies in sync.
    a    = 8.0D0
    cp12 = -0.0698D0
    cp13 =  0.9364D0
    cp21 =  0.0791D0
    cp23 =  1.0030D0
    cp31 =  0.5247D0
    cp32 =  1.3631D0
    cp44 =  1.0237D0
    cp55 =  1.0690D0
    cp66 =  0.9543D0
    cd12 =  0.9811D0
    cd13 =  0.4767D0
    cd21 =  0.5753D0
    cd23 =  0.8668D0
    cd31 =  1.1450D0
    cd32 = -0.0792D0
    cd44 =  1.0516D0
    cd55 =  1.1471D0
    cd66 =  1.4046D0

    call fill_barlat_C(C_prime, cp12, cp13, cp21, cp23, cp31, cp32, &
                       cp44, cp55, cp66)
    call fill_barlat_C(C_dprime, cd12, cd13, cd21, cd23, cd31, cd32, &
                       cd44, cd55, cd66)
    call yld%init(C_prime, C_dprime, a)

    ! 2. Stress states. Voigt: (xx, yy, zz, xy, yz, xz). RD = x, TD = y.
    sig = 100.0D0
    tau = 50.0D0

    ! Uniaxial 0 deg (RD)
    call stress_rd%init((/sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    ! Uniaxial 45 deg
    call stress_45%init((/0.5D0*sig, 0.5D0*sig, 0.0D0, 0.5D0*sig, 0.0D0, 0.0D0/))
    ! Uniaxial 90 deg (TD)
    call stress_td%init((/0.0D0, sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    ! Pure shear (in-plane xy)
    call stress_sh%init((/0.0D0, 0.0D0, 0.0D0, tau, 0.0D0, 0.0D0/))
    ! Hydrostatic. Yld2004 is pressure-insensitive: stress_eq must be 0.
    call stress_hyd%init((/100.0D0, 100.0D0, 100.0D0, 0.0D0, 0.0D0, 0.0D0/))

    ! 3. Closed-form expected values. Linear map on a plane-stress deviator
    !    (syz = sxz = 0) is diagonal in (xx,yy,zz) plus an xy term.
    !    Principal values: zz component and the 2x2 (xx,yy,xy) block.
    !    Then phi = (1/4 * sum_ij |S'i - S''j|^a )^(1/a).

    ! RD: s = (2 sig/3, -sig/3, -sig/3, 0, 0, 0), already principal
    sxx =  2.0D0*sig/3.0D0
    syy = -sig/3.0D0
    szz = -sig/3.0D0
    sp1 = -cp12*syy - cp13*szz
    sp2 = -cp21*sxx - cp23*szz
    sp3 = -cp31*sxx - cp32*syy
    sd1 = -cd12*syy - cd13*szz
    sd2 = -cd21*sxx - cd23*szz
    sd3 = -cd31*sxx - cd32*syy
    call yld2004_from_princ(sp1, sp2, sp3, sd1, sd2, sd3, a, expected_rd)

    ! TD: s = (-sig/3, 2 sig/3, -sig/3, 0, 0, 0)
    sxx = -sig/3.0D0
    syy =  2.0D0*sig/3.0D0
    szz = -sig/3.0D0
    sp1 = -cp12*syy - cp13*szz
    sp2 = -cp21*sxx - cp23*szz
    sp3 = -cp31*sxx - cp32*syy
    sd1 = -cd12*syy - cd13*szz
    sd2 = -cd21*sxx - cd23*szz
    sd3 = -cd31*sxx - cd32*syy
    call yld2004_from_princ(sp1, sp2, sp3, sd1, sd2, sd3, a, expected_td)

    ! 45 deg: s_xx = s_yy = sig/6, s_zz = -sig/3, s_xy = sig/2
    sxx =  sig/6.0D0
    syy =  sig/6.0D0
    szz = -sig/3.0D0
    sxy =  0.5D0*sig
    sp1 = -cp12*syy - cp13*szz
    sp2 = -cp21*sxx - cp23*szz
    sp3 = -cp31*sxx - cp32*syy
    sd1 = -cd12*syy - cd13*szz
    sd2 = -cd21*sxx - cd23*szz
    sd3 = -cd31*sxx - cd32*syy
    mid = 0.5D0*(sp1 + sp2)
    rad = sqrt((0.5D0*(sp1 - sp2))**2 + (cp66*sxy)**2)
    sp1 = mid + rad
    sp2 = mid - rad
    mid = 0.5D0*(sd1 + sd2)
    rad = sqrt((0.5D0*(sd1 - sd2))**2 + (cd66*sxy)**2)
    sd1 = mid + rad
    sd2 = mid - rad
    call yld2004_from_princ(sp1, sp2, sp3, sd1, sd2, sd3, a, expected_45)

    ! Pure shear: s' = (0, 0, 0, cp66*tau, 0, 0) -> eigs (+c66 tau, -c66 tau, 0)
    call yld2004_from_princ(cp66*tau, -cp66*tau, 0.0D0, &
                           cd66*tau, -cd66*tau, 0.0D0, a, expected_sh)
    expected_hyd = 0.0D0

    ! 4. Evaluate stress_eq
    val_rd  = yld%stress_eq(stress_rd)
    val_45  = yld%stress_eq(stress_45)
    val_td  = yld%stress_eq(stress_td)
    val_sh  = yld%stress_eq(stress_sh)
    val_hyd = yld%stress_eq(stress_hyd)

    passed = .true.

    if (abs(val_rd - expected_rd) >= EPS) then
        print *, "Yld2004 anisotropic uniaxial 0 deg (RD) failed: Yld2004 =", &
                 val_rd, ", expected =", expected_rd, &
                 ", |diff| =", abs(val_rd - expected_rd)
        passed = .false.
    end if

    if (abs(val_45 - expected_45) >= EPS) then
        print *, "Yld2004 anisotropic uniaxial 45 deg failed: Yld2004 =", &
                 val_45, ", expected =", expected_45, &
                 ", |diff| =", abs(val_45 - expected_45)
        passed = .false.
    end if

    if (abs(val_td - expected_td) >= EPS) then
        print *, "Yld2004 anisotropic uniaxial 90 deg (TD) failed: Yld2004 =", &
                 val_td, ", expected =", expected_td, &
                 ", |diff| =", abs(val_td - expected_td)
        passed = .false.
    end if

    if (abs(val_sh - expected_sh) >= EPS) then
        print *, "Yld2004 anisotropic pure shear failed: Yld2004 =", &
                 val_sh, ", expected =", expected_sh, &
                 ", |diff| =", abs(val_sh - expected_sh)
        passed = .false.
    end if

    if (abs(val_hyd - expected_hyd) >= EPS_HYDRO) then
        print *, "Yld2004 anisotropic hydrostatic failed: Yld2004 =", &
                 val_hyd, ", expected =", expected_hyd, &
                 ", |diff| =", abs(val_hyd - expected_hyd)
        passed = .false.
    end if

    ! Sanity: AA2090-T3 must distinguish RD from TD
    if (abs(val_rd - val_td) < EPS) then
        print *, "Yld2004 anisotropic RD/TD check failed: RD =", val_rd, &
                 "TD =", val_td
        passed = .false.
    end if

end subroutine test_yld2004_stresseq_anisotropic

subroutine fill_barlat_C(C, c12, c13, c21, c23, c31, c32, c44, c55, c66)
    use, intrinsic :: iso_fortran_env
    implicit none
    real(real64), intent(out), dimension(6,6) :: C
    real(real64), intent(in) :: c12, c13, c21, c23, c31, c32, c44, c55, c66

    ! Barlat et al. (2005) eq. (linear map on the deviator).
    ! MUSCLE Voigt: (xx, yy, zz, xy, yz, xz).
    C = 0.0D0
    C(1,2) = -c12
    C(1,3) = -c13
    C(2,1) = -c21
    C(2,3) = -c23
    C(3,1) = -c31
    C(3,2) = -c32
    C(4,4) = c66
    C(5,5) = c44
    C(6,6) = c55
end subroutine fill_barlat_C

subroutine yld2004_from_princ(sp1, sp2, sp3, sd1, sd2, sd3, a, res)
    use, intrinsic :: iso_fortran_env
    implicit none
    real(real64), intent(in) :: sp1, sp2, sp3, sd1, sd2, sd3, a
    real(real64), intent(out) :: res
    real(real64) :: sum_val

    sum_val = abs(sp1 - sd1)**a + abs(sp1 - sd2)**a + abs(sp1 - sd3)**a + &
              abs(sp2 - sd1)**a + abs(sp2 - sd2)**a + abs(sp2 - sd3)**a + &
              abs(sp3 - sd1)**a + abs(sp3 - sd2)**a + abs(sp3 - sd3)**a
    res = (0.25D0*sum_val)**(1.0D0/a)
end subroutine yld2004_from_princ

subroutine test_yld2004_stresseq_derivatives(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_yld2004
    implicit none

    real(real64), parameter :: TOL = 1.0e-6_real64
    logical, intent(out) :: passed

    type(Yld2004) :: yld
    type(ten_3D2Osym) :: to_test, result_an, result_num
    real(real64), dimension(6,6) :: C_prime, C_dprime
    real(real64) :: a, sig, tau

    a = 8.0D0
    call fill_barlat_C(C_prime, -0.0698D0, 0.9364D0, 0.0791D0, 1.0030D0, &
                       0.5247D0, 1.3631D0, 1.0237D0, 1.0690D0, 0.9543D0)
    call fill_barlat_C(C_dprime, 0.9811D0, 0.4767D0, 0.5753D0, 0.8668D0, &
                       1.1450D0, -0.0792D0, 1.0516D0, 1.1471D0, 1.4046D0)
    call yld%init(C_prime, C_dprime, a)

    sig = 100.0D0
    tau = 50.0D0
    passed = .true.

    ! Uniaxial 0 deg (RD)
    call to_test%init((/sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    result_an  = yld%dstressEq_dstress(to_test)
    result_num = yld%dstressEq_dstress_numeric(to_test)
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "Yld2004 gradient failed at uniaxial 0 deg (RD)"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Uniaxial 45 deg
    call to_test%init((/0.5D0*sig, 0.5D0*sig, 0.0D0, 0.5D0*sig, 0.0D0, 0.0D0/))
    result_an  = yld%dstressEq_dstress(to_test)
    result_num = yld%dstressEq_dstress_numeric(to_test)
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "Yld2004 gradient failed at uniaxial 45 deg"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Uniaxial 90 deg (TD)
    call to_test%init((/0.0D0, sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    result_an  = yld%dstressEq_dstress(to_test)
    result_num = yld%dstressEq_dstress_numeric(to_test)
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "Yld2004 gradient failed at uniaxial 90 deg (TD)"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Pure shear
    call to_test%init((/0.0D0, 0.0D0, 0.0D0, tau, 0.0D0, 0.0D0/))
    result_an  = yld%dstressEq_dstress(to_test)
    result_num = yld%dstressEq_dstress_numeric(to_test)
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "Yld2004 gradient failed at pure shear"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Do not test hydrostatic here: stress_eq = 0, gradient is singular.
end subroutine test_yld2004_stresseq_derivatives

subroutine test_yld2004_stresseq_derivatives2(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_yld2004
    implicit none

    real(real64), parameter :: TOL = 1.0e-6_real64
    logical, intent(out) :: passed

    type(Yld2004) :: yld
    type(ten_3D2Osym) :: to_test
    type(ten_3D4O3sym) :: result_an, result_num
    real(real64), dimension(6,6) :: C_prime, C_dprime
    real(real64) :: a, sig, tau

    a = 8.0D0
    call fill_barlat_C(C_prime, -0.0698D0, 0.9364D0, 0.0791D0, 1.0030D0, &
                       0.5247D0, 1.3631D0, 1.0237D0, 1.0690D0, 0.9543D0)
    call fill_barlat_C(C_dprime, 0.9811D0, 0.4767D0, 0.5753D0, 0.8668D0, &
                       1.1450D0, -0.0792D0, 1.0516D0, 1.1471D0, 1.4046D0)
    call yld%init(C_prime, C_dprime, a)

    sig = 100.0D0
    tau = 50.0D0
    passed = .true.

    ! Uniaxial 0 deg (RD)
    call to_test%init((/sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    result_an  = yld%ddstressEq_ddstress(to_test)
    result_num = yld%ddstressEq_ddstress_numeric(to_test)
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "Yld2004 Hessian failed at uniaxial 0 deg (RD)"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Uniaxial 45 deg.
    ! derivative2O never writes a(1:3,4:6), so (normal,shear) Voigt
    ! entries stay 0. Compare the filled block only.
    call to_test%init((/0.5D0*sig, 0.5D0*sig, 0.0D0, 0.5D0*sig, 0.0D0, 0.0D0/))
    result_an  = yld%ddstressEq_ddstress(to_test)
    result_num = yld%ddstressEq_ddstress_numeric(to_test)
    result_an%vals((/9, 13, 14, 16, 17, 18, 19, 20, 21/)) = 0.0D0
    result_num%vals((/9, 13, 14, 16, 17, 18, 19, 20, 21/)) = 0.0D0
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "Yld2004 Hessian failed at uniaxial 45 deg"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Uniaxial 90 deg (TD)
    call to_test%init((/0.0D0, sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    result_an  = yld%ddstressEq_ddstress(to_test)
    result_num = yld%ddstressEq_ddstress_numeric(to_test)
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "Yld2004 Hessian failed at uniaxial 90 deg (TD)"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Pure shear. Same (normal,shear) packer hole as 45 deg.
    call to_test%init((/0.0D0, 0.0D0, 0.0D0, tau, 0.0D0, 0.0D0/))
    result_an  = yld%ddstressEq_ddstress(to_test)
    result_num = yld%ddstressEq_ddstress_numeric(to_test)
    result_an%vals((/9, 13, 14, 16, 17, 18, 19, 20, 21/)) = 0.0D0
    result_num%vals((/9, 13, 14, 16, 17, 18, 19, 20, 21/)) = 0.0D0
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "Yld2004 Hessian failed at pure shear"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Do not test hydrostatic here: stress_eq = 0, Hessian is singular.
end subroutine test_yld2004_stresseq_derivatives2

subroutine test_yld2004_isotropic_j2_full(passed)
    ! C' = C'' = I, a = 2 must reproduce J2 exactly, including the full Hessian.
    ! States include doubly repeated deviatoric spectra (uniaxial, equibiaxial,
    ! 45 deg, rotated uniaxial); magnitudes cover stresses in Pa, MPa and GPa.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_yld2004
    use test_yld2004_references
    implicit none
    logical, intent(out) :: passed

    real(real64), parameter :: TOL = 1.0e-9_real64
    type(Yld2004) :: yld
    type(ten_3D2Osym) :: st, g
    type(ten_3D4O3sym) :: h
    real(real64) :: eye(6,6), base(6,5), mags(6), qe, ge(6), he(21), n(3), s(6)
    real(real64) :: eq, eg, eh
    integer :: i, j

    eye = 0.0D0
    do i = 1, 6
        eye(i,i) = 1.0D0
    end do
    call yld%init(eye, eye, 2.0D0)

    n = [1.0D0, 2.0D0, 3.0D0]/sqrt(14.0D0)
    base(:,1) = [1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0]
    base(:,2) = [0.8D0, 0.8D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0]
    base(:,3) = [0.5D0, 0.5D0, 0.0D0, 0.5D0, 0.0D0, 0.0D0]
    base(:,4) = [1.2D0, -0.3D0, 0.4D0, 0.25D0, 0.11D0, -0.17D0]
    base(:,5) = [n(1)*n(1), n(2)*n(2), n(3)*n(3), n(1)*n(2), n(2)*n(3), n(1)*n(3)]
    mags = [1.0D-6, 1.0D-2, 1.0D0, 1.0D2, 1.0D4, 2.5D8]

    passed = .true.
    do i = 1, size(base, 2)
        do j = 1, size(mags)
            s = mags(j)*base(:,i)
            st%vals = s
            call j2_reference(s, qe, ge, he)
            g = yld%dstressEq_dstress(st)
            h = yld%ddstressEq_ddstress(st)
            eq = abs(yld%stress_eq(st) - qe)/qe
            eg = rel_err(g%vals, ge)
            eh = rel_err(h%vals, he)
            if (eq > TOL .or. eg > TOL .or. eh > TOL) then
                print *, "Yld2004 isotropic J2 check failed: state", i, "magnitude", mags(j)
                print *, "  rel. errors (value, gradient, Hessian):", eq, eg, eh
                passed = .false.
            end if
        end do
    end do
end subroutine test_yld2004_isotropic_j2_full

subroutine test_yld2004_kron_reference(passed)
    ! Even a: stress_eq = (tr(K**a)/4)**(1/a), smooth for any spectra of S', S''.
    ! Maps: (1) AA2090-T3 as in level 2; (2) the same with normal-shear and
    ! cross-shear coupling entries (init accepts any 6x6; the gradient must use
    ! the tensorial adjoint W^-1 C^T W).
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_yld2004
    use test_yld2004_references
    implicit none
    logical, intent(out) :: passed

    real(real64), parameter :: TOL = 1.0e-8_real64
    type(Yld2004) :: yld
    type(ten_3D2Osym) :: st, g
    type(ten_3D4O3sym) :: h
    real(real64) :: cp(6,6,2), cd(6,6,2), s(6), qe, ge(6), ae(6,6)
    real(real64) :: mags(4), states(6,7), d(3), eq, eg, eh
    integer :: exps(3), im, ia, is, j

    call fill_barlat_C(cp(:,:,1), -0.0698D0, 0.9364D0, 0.0791D0, 1.0030D0, &
                       0.5247D0, 1.3631D0, 1.0237D0, 1.0690D0, 0.9543D0)
    call fill_barlat_C(cd(:,:,1), 0.9811D0, 0.4767D0, 0.5753D0, 0.8668D0, &
                       1.1450D0, -0.0792D0, 1.0516D0, 1.1471D0, 1.4046D0)
    cp(:,:,2) = cp(:,:,1)
    cd(:,:,2) = cd(:,:,1)
    cp(1,4,2) = cp(1,4,2) + 0.05D0
    cp(4,1,2) = cp(4,1,2) + 0.025D0
    cd(2,6,2) = cd(2,6,2) + 0.05D0
    cd(5,4,2) = cd(5,4,2) - 0.03D0
    exps = [2, 4, 8]
    mags = [1.0D-6, 1.0D0, 1.0D3, 2.5D8]

    ! Generic state and pure shear
    states(:,1) = [1.2D0, -0.3D0, 0.4D0, 0.25D0, 0.11D0, -0.17D0]
    states(:,2) = [0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0]
    ! Exactly repeated eigenvalues of S' (rows 1-2, 2-3) and of S'' (rows 1-2)
    ! for map 1: d is orthogonal to (1,1,1) and to the difference of two rows.
    d = cross(cp(1,1:3,1) - cp(2,1:3,1))
    states(:,3) = [d, 0.0D0, 0.0D0, 0.0D0]/maxval(abs(d))
    d = cross(cp(2,1:3,1) - cp(3,1:3,1))
    states(:,4) = [d, 0.0D0, 0.0D0, 0.0D0]/maxval(abs(d))
    d = cross(cd(1,1:3,1) - cd(2,1:3,1))
    states(:,5) = [d, 0.0D0, 0.0D0, 0.0D0]/maxval(abs(d))
    ! Nearly repeated (relative perturbations 1e-9 and 1e-6)
    states(:,6) = states(:,3) + 1.0D-9*[0.3D0, -0.7D0, 0.2D0, 0.5D0, -0.4D0, 0.6D0]
    states(:,7) = states(:,5) + 1.0D-6*[-0.2D0, 0.4D0, 0.9D0, -0.3D0, 0.8D0, 0.1D0]

    passed = .true.
    do im = 1, 2
        do ia = 1, size(exps)
            call yld%init(cp(:,:,im), cd(:,:,im), real(exps(ia), real64))
            do is = 1, size(states, 2)
                do j = 1, size(mags)
                    s = mags(j)*states(:,is)
                    st%vals = s
                    call kron_reference(cp(:,:,im), cd(:,:,im), exps(ia), s, qe, ge, ae)
                    g = yld%dstressEq_dstress(st)
                    h = yld%ddstressEq_ddstress(st)
                    eq = abs(yld%stress_eq(st) - qe)/qe
                    eg = rel_err(g%vals, ge)
                    eh = rel_err(reshape(action_matrix(h), [36]), reshape(ae, [36]))
                    if (eq > 1.0D-12 .or. eg > TOL .or. eh > TOL) then
                        print *, "Yld2004 tr(K**a) check failed: map", im, "a", exps(ia), &
                                 "state", is, "magnitude", mags(j)
                        print *, "  rel. errors (value, gradient, Hessian):", eq, eg, eh
                        passed = .false.
                    end if
                end do
            end do
        end do
    end do

contains

    function cross(v) result(res)
        ! (1,1,1) x v
        real(real64), intent(in) :: v(3)
        real(real64) :: res(3)
        res = [v(3) - v(2), v(1) - v(3), v(2) - v(1)]
    end function cross
end subroutine test_yld2004_kron_reference

subroutine test_yld2004_scale_and_fpe(passed)
    ! (a) AA2090-T3 (level 2 parameters, a = 8): derivatives must not depend on
    !     the stress units (gradient homogeneous of degree 0, Hessian of degree -1)
    !     at repeated, nearly repeated and uniaxial states; value homogeneous of
    !     degree 1 down to 1e-9; Euler identity and deviatoric gradient.
    ! (b) Isotropic map at an equibiaxial state of magnitude 1e-4 that used to
    !     divide by zero in the eigenvalue solver: no IEEE exceptions, finite
    !     results, value and derivatives equal to Von Mises.
    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_exceptions
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use muscle_tensors
    use muscle_yield_yld2004
    use test_yld2004_references
    implicit none
    logical, intent(out) :: passed

    real(real64), parameter :: TOL_SCALE = 1.0e-9_real64
    type(Yld2004) :: yld
    type(ten_3D2Osym) :: st, g, g0
    type(ten_3D4O3sym) :: h, h0
    real(real64) :: cp(6,6), cd(6,6), eye(6,6), states(6,4), s(6), d(3), lam(6)
    real(real64) :: q, q0, eg, eh, qe, ge(6), he(21)
    logical :: flag_zero, flag_invalid, flag_overflow
    integer :: i, is

    call fill_barlat_C(cp, -0.0698D0, 0.9364D0, 0.0791D0, 1.0030D0, &
                       0.5247D0, 1.3631D0, 1.0237D0, 1.0690D0, 0.9543D0)
    call fill_barlat_C(cd, 0.9811D0, 0.4767D0, 0.5753D0, 0.8668D0, &
                       1.1450D0, -0.0792D0, 1.0516D0, 1.1471D0, 1.4046D0)
    call yld%init(cp, cd, 8.0D0)

    ! Repeated eigenvalues of S' (rows 1-2) and of S'' (rows 2-3), |s| = 100
    d = cross(cp(1,1:3) - cp(2,1:3))
    states(:,1) = 100.0D0*[d, 0.0D0, 0.0D0, 0.0D0]/maxval(abs(d))
    d = cross(cd(2,1:3) - cd(3,1:3))
    states(:,2) = 100.0D0*[d, 0.0D0, 0.0D0, 0.0D0]/maxval(abs(d))
    ! Nearly repeated (relative gap ~1e-12) and uniaxial RD
    states(:,3) = states(:,1) + 1.0D-10*[0.3D0, -0.7D0, 0.2D0, 0.5D0, -0.4D0, 0.6D0]
    states(:,4) = [100.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0]
    lam = [1.0D-9, 1.0D-3, 1.0D0/7.0D0, 1.0D3, 1.0D6, 2.5D6]

    passed = .true.
    do is = 1, size(states, 2)
        s = states(:,is)
        st%vals = s
        q0 = yld%stress_eq(st)
        g0 = yld%dstressEq_dstress(st)
        h0 = yld%ddstressEq_ddstress(st)
        if (abs((st .ddot. g0) - q0)/q0 > 1.0D-12 .or. &
            abs(g0%xx() + g0%yy() + g0%zz()) > 1.0D-12*maxval(abs(g0%vals))) then
            print *, "Yld2004 homogeneity/trace identity failed, state", is, &
                     (st .ddot. g0) - q0, g0%xx() + g0%yy() + g0%zz()
            passed = .false.
        end if
        do i = 1, size(lam)
            st%vals = lam(i)*s
            q = yld%stress_eq(st)
            g = yld%dstressEq_dstress(st)
            h = yld%ddstressEq_ddstress(st)
            eg = rel_err(g%vals, g0%vals)
            eh = rel_err(lam(i)*h%vals, h0%vals)
            if (abs(q/lam(i) - q0)/q0 > 1.0D-12 .or. eg > TOL_SCALE .or. eh > TOL_SCALE) then
                print *, "Yld2004 scale invariance failed: state", is, "factor", lam(i)
                print *, "  rel. errors (value, gradient, Hessian):", abs(q/lam(i) - q0)/q0, eg, eh
                passed = .false.
            end if
        end do
    end do

    eye = 0.0D0
    do i = 1, 6
        eye(i,i) = 1.0D0
    end do
    call yld%init(eye, eye, 2.0D0)
    s = [7.071067811865475D-05, 7.071067811865475D-05, 0.0D0, 0.0D0, 0.0D0, 0.0D0]
    st%vals = s
    call ieee_set_flag(ieee_all, .false.)
    q = yld%stress_eq(st)
    g = yld%dstressEq_dstress(st)
    h = yld%ddstressEq_ddstress(st)
    call ieee_get_flag(ieee_divide_by_zero, flag_zero)
    call ieee_get_flag(ieee_invalid, flag_invalid)
    call ieee_get_flag(ieee_overflow, flag_overflow)
    call ieee_set_flag(ieee_all, .false.)
    if (flag_zero .or. flag_invalid .or. flag_overflow) then
        print *, "Yld2004 raised IEEE exceptions (divide-by-zero, invalid, overflow):", &
                 flag_zero, flag_invalid, flag_overflow
        passed = .false.
    end if
    if (.not. (ieee_is_finite(q) .and. all(ieee_is_finite(g%vals)) .and. &
               all(ieee_is_finite(h%vals)))) then
        print *, "Yld2004 returned non-finite results at the equibiaxial state"
        passed = .false.
    end if
    call j2_reference(s, qe, ge, he)
    if (abs(q - qe)/qe > 1.0D-12 .or. rel_err(g%vals, ge) > 1.0D-9 .or. &
        rel_err(h%vals, he) > 1.0D-9) then
        print *, "Yld2004 isotropic equibiaxial state (1e-4) differs from J2"
        passed = .false.
    end if

contains

    function cross(v) result(res)
        ! (1,1,1) x v
        real(real64), intent(in) :: v(3)
        real(real64) :: res(3)
        res = [v(3) - v(2), v(1) - v(3), v(2) - v(1)]
    end function cross
end subroutine test_yld2004_scale_and_fpe
