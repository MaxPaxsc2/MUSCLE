! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_yield_cpb06
    implicit none
    
    logical :: passed

    call test_CPB06_stresseq_isotropic(passed)
    if (.not. passed) STOP 1

    call test_CPB06_stresseq_anisotropic(passed)
    if (.not. passed) STOP 2

    ! Analytical vs numerical first derivative (gradient)
    call test_CPB06_stresseq_derivatives(passed)
    if (.not. passed) STOP 3

    ! Analytical vs numerical second derivative (Hessian)
    call test_CPB06_stresseq_derivatives2(passed)
    if (.not. passed) STOP 4

    ! Regression: isotropic limit vs closed-form J2, full Hessian, repeated
    ! spectra and stress magnitudes from 1e-6 to 2.5e8
    call test_CPB06_isotropic_j2_full(passed)
    if (.not. passed) STOP 5

    ! Regression: k = 0, even a, anisotropic C vs tr(Sigma**a) reference at
    ! generic, repeated and nearly repeated transformed spectra, several magnitudes
    call test_CPB06_tracepower_reference(passed)
    if (.not. passed) STOP 6

    ! Regression: published k /= 0 set (repeated transformed spectrum), scale
    ! invariance of the derivatives and absence of IEEE exceptions
    call test_CPB06_published_repeated(passed)
    if (.not. passed) STOP 7

    print*, "Passed!", passed
    STOP 0
end program test_muscle_yield_cpb06

subroutine test_CPB06_stresseq_isotropic(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    use muscle_yield_vonmises
    implicit none
    
    real(real64), parameter :: EPS = 1.0e-5_real64
    logical, intent(out) :: passed

    type(CPB06) :: cpb
    type(VonMises) :: vm
    type(ten_3D2Osym) :: stress1, stress2, stress3

    real(real64) :: val_cpb1, val_vm1
    real(real64) :: val_cpb2, val_vm2
    real(real64) :: val_cpb3, val_vm3

    ! 1. Instantiate CPB06 with isotropic parameters
    call cpb%init(c11=1.0D0, c12=0.0D0, c13=0.0D0, &
                  c21=0.0D0, c22=1.0D0, c23=0.0D0, &
                  c31=0.0D0, c32=0.0D0, c33=1.0D0, &
                  c44=1.0D0, c55=1.0D0, c66=1.0D0, &
                  k=0.0D0, a=2.0D0                 &
                  )

    ! 2. Stress states
    ! Uniaxial tension in X: (100, 0, 0, 0, 0, 0)
    call stress1%init((/100.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    ! Pure shear: (0, 0, 0, 50, 0, 0)
    call stress2%init((/0.0D0, 0.0D0, 0.0D0, 50.0D0, 0.0D0, 0.0D0/))
    ! Biaxial: (80, 80, 0, 0, 0, 0)
    call stress3%init((/80.0D0, 80.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))

    ! 3. Evaluate both yield criteria
    val_cpb1 = cpb%stress_eq(stress1)
    val_vm1 = vm%stress_eq(stress1)

    val_cpb2 = cpb%stress_eq(stress2)
    val_vm2 = vm%stress_eq(stress2)

    val_cpb3 = cpb%stress_eq(stress3)
    val_vm3 = vm%stress_eq(stress3)

    ! 4. Verify abs(CPB06%stress_eq(sig) - VonMises%stress_eq(sig)) < 1.0e-5
    passed = .true.

    if (abs(val_cpb1 - val_vm1) >= EPS) then
        print *, "CPB06 uniaxial test failed: CPB06 =", val_cpb1, "VonMises =", val_vm1
        passed = .false.
    end if

    if (abs(val_cpb2 - val_vm2) >= EPS) then
        print *, "CPB06 pure shear test failed: CPB06 =", val_cpb2, "VonMises =", val_vm2
        passed = .false.
    end if

    if (abs(val_cpb3 - val_vm3) >= EPS) then
        print *, "CPB06 biaxial test failed: CPB06 =", val_cpb3, "VonMises =", val_vm3
        passed = .false.
    end if

end subroutine test_CPB06_stresseq_isotropic

subroutine test_CPB06_stresseq_anisotropic(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    use muscle_math_operations, only: eigenvals
    implicit none

    real(real64), parameter :: EPS = 1.0e-8_real64
    real(real64), parameter :: EPS_HYDRO = 1.0e-8_real64
    real(real64), parameter :: EPS_EIG = 1.0e-8_real64
    logical, intent(out) :: passed

    type(CPB06) :: cpb
    type(ten_3D2Osym) :: stress_rd_t, stress_rd_c
    type(ten_3D2Osym) :: stress_td_t, stress_td_c
    type(ten_3D2Osym) :: stress_sh, stress_hyd

    real(real64) :: c11, c12, c13, c21, c22, c23, c31, c32, c33
    real(real64) :: c44, c55, c66, k, a
    real(real64) :: g1, g2, g3, d1, d2, d3, B
    real(real64) :: phi_t, phi_c
    real(real64) :: sig, tau
    real(real64) :: val_rd_t, val_rd_c, val_td_t, val_td_c, val_sh, val_hyd
    real(real64) :: expected_rd_t, expected_rd_c
    real(real64) :: expected_td_t, expected_td_c
    real(real64) :: expected_sh, expected_hyd
    real(real64), dimension(3) :: eigs
    integer :: n_zero

    ! 1. Ti-6Al-4V, CPB06 at initial yield (Tuninetti et al., Int. J. Plasticity
    !    67, 2015, Table 3, Wp = 1.857 J/cm3). Criterion of Cazacu, Plunkett,
    !    Barlat, Int. J. Plasticity 22, 2006. a = 2, k /= 0 (SD effect).
    !    C is symmetric; C44 = C55 = C66.
    c11 = 1.000D0
    c12 = 2.373D0
    c13 = 2.364D0
    c21 = 2.373D0
    c22 = 1.838D0
    c23 = 1.196D0
    c31 = 2.364D0
    c32 = 1.196D0
    c33 = 2.444D0
    c44 = 3.607D0
    c55 = 3.607D0
    c66 = 3.607D0
    k   = 0.136D0
    a   = 2.000D0

    call cpb%init(c11=c11, c12=c12, c13=c13, &
                  c21=c21, c22=c22, c23=c23, &
                  c31=c31, c32=c32, c33=c33, &
                  c44=c44, c55=c55, c66=c66, &
                  k=k, a=a                   &
                  )

    ! 2. Stress states. Voigt: (xx, yy, zz, xy, yz, xz). RD = x, TD = y.
    sig = 100.0D0
    tau = 50.0D0

    ! Uniaxial tension RD (also repeated Cauchy eigenvalues: sig, 0, 0)
    call stress_rd_t%init((/sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    ! Uniaxial compression RD
    call stress_rd_c%init((/-sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    ! Uniaxial tension TD
    call stress_td_t%init((/0.0D0, sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    ! Uniaxial compression TD
    call stress_td_c%init((/0.0D0, -sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    ! Pure shear (in-plane xy)
    call stress_sh%init((/0.0D0, 0.0D0, 0.0D0, tau, 0.0D0, 0.0D0/))
    ! Hydrostatic. CPB06 is pressure-insensitive: stress_eq must be 0.
    call stress_hyd%init((/100.0D0, 100.0D0, 100.0D0, 0.0D0, 0.0D0, 0.0D0/))

    ! 3. Closed-form special cases (not a copy of the six-component path).
    !    B normalizes so uniaxial RD tension of magnitude 1 gives stress_eq = 1:
    !      gamma_i = (2 C_i1 - C_i2 - C_i3)/3
    !      B = [sum (|gamma_i| - k gamma_i)^a ]^{-1/a}
    !    RD tension:     sig
    !    RD compression: B * sig * [sum (|gamma_i| + k gamma_i)^a]^{1/a}
    !    TD: same with delta_i = (-C_i1 + 2 C_i2 - C_i3)/3
    !    Pure shear:     B * [ (|c44 tau| - k c44 tau)^a
    !                        + (|c44 tau| + k c44 tau)^a ]^{1/a}
    !    Hydrostatic:    0
    g1 = (2.0D0*c11 - c12 - c13)/3.0D0
    g2 = (2.0D0*c21 - c22 - c23)/3.0D0
    g3 = (2.0D0*c31 - c32 - c33)/3.0D0
    phi_t = ((abs(g1) - k*g1)**a + (abs(g2) - k*g2)**a + (abs(g3) - k*g3)**a)**(1.0D0/a)
    B = 1.0D0/phi_t

    d1 = (-c11 + 2.0D0*c12 - c13)/3.0D0
    d2 = (-c21 + 2.0D0*c22 - c23)/3.0D0
    d3 = (-c31 + 2.0D0*c32 - c33)/3.0D0

    expected_rd_t = sig
    phi_c = ((abs(g1) + k*g1)**a + (abs(g2) + k*g2)**a + (abs(g3) + k*g3)**a)**(1.0D0/a)
    expected_rd_c = B*sig*phi_c

    expected_td_t = B*sig*((abs(d1) - k*d1)**a + (abs(d2) - k*d2)**a + &
                           (abs(d3) - k*d3)**a)**(1.0D0/a)
    expected_td_c = B*sig*((abs(d1) + k*d1)**a + (abs(d2) + k*d2)**a + &
                           (abs(d3) + k*d3)**a)**(1.0D0/a)

    expected_sh = B*((abs(c44*tau) - k*c44*tau)**a + &
                     (abs(c44*tau) + k*c44*tau)**a)**(1.0D0/a)
    expected_hyd = 0.0D0

    ! 4. Evaluate stress_eq
    val_rd_t = cpb%stress_eq(stress_rd_t)
    val_rd_c = cpb%stress_eq(stress_rd_c)
    val_td_t = cpb%stress_eq(stress_td_t)
    val_td_c = cpb%stress_eq(stress_td_c)
    val_sh   = cpb%stress_eq(stress_sh)
    val_hyd  = cpb%stress_eq(stress_hyd)

    passed = .true.

    if (abs(val_rd_t - expected_rd_t) >= EPS) then
        print *, "CPB06 anisotropic uniaxial RD tension failed: CPB06 =", &
                 val_rd_t, ", expected =", expected_rd_t, &
                 ", |diff| =", abs(val_rd_t - expected_rd_t)
        passed = .false.
    end if

    if (abs(val_rd_c - expected_rd_c) >= EPS) then
        print *, "CPB06 anisotropic uniaxial RD compression failed: CPB06 =", &
                 val_rd_c, ", expected =", expected_rd_c, &
                 ", |diff| =", abs(val_rd_c - expected_rd_c)
        passed = .false.
    end if

    if (abs(val_rd_t - val_rd_c) < EPS) then
        print *, "CPB06 SD check RD failed: tension =", val_rd_t, &
                 "compression =", val_rd_c
        passed = .false.
    end if

    if (abs(val_td_t - expected_td_t) >= EPS) then
        print *, "CPB06 anisotropic uniaxial TD tension failed: CPB06 =", &
                 val_td_t, ", expected =", expected_td_t, &
                 ", |diff| =", abs(val_td_t - expected_td_t)
        passed = .false.
    end if

    if (abs(val_td_c - expected_td_c) >= EPS) then
        print *, "CPB06 anisotropic uniaxial TD compression failed: CPB06 =", &
                 val_td_c, ", expected =", expected_td_c, &
                 ", |diff| =", abs(val_td_c - expected_td_c)
        passed = .false.
    end if

    if (abs(val_td_t - val_td_c) < EPS) then
        print *, "CPB06 SD check TD failed: tension =", val_td_t, &
                 "compression =", val_td_c
        passed = .false.
    end if

    if (abs(val_sh - expected_sh) >= EPS) then
        print *, "CPB06 anisotropic pure shear failed: CPB06 =", &
                 val_sh, ", expected =", expected_sh, &
                 ", |diff| =", abs(val_sh - expected_sh)
        passed = .false.
    end if

    if (abs(val_hyd - expected_hyd) >= EPS_HYDRO) then
        print *, "CPB06 anisotropic hydrostatic failed: CPB06 =", &
                 val_hyd, ", expected =", expected_hyd, &
                 ", |diff| =", abs(val_hyd - expected_hyd)
        passed = .false.
    end if

    ! Repeated Cauchy eigenvalues: uniaxial RD has principal stresses (sig, 0, 0).
    eigs = eigenvals(stress_rd_t)
    n_zero = 0
    if (abs(eigs(1)) < EPS_EIG) n_zero = n_zero + 1
    if (abs(eigs(2)) < EPS_EIG) n_zero = n_zero + 1
    if (abs(eigs(3)) < EPS_EIG) n_zero = n_zero + 1
    if (n_zero < 2) then
        print *, "CPB06 repeated-eigenvalue state failed: eigs =", eigs
        passed = .false.
    end if
    if (abs(val_rd_t - expected_rd_t) >= EPS) then
        print *, "CPB06 stress_eq on repeated-eigenvalue state failed: CPB06 =", &
                 val_rd_t, ", expected =", expected_rd_t
        passed = .false.
    end if

end subroutine test_CPB06_stresseq_anisotropic

subroutine test_CPB06_stresseq_derivatives(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    implicit none

    real(real64), parameter :: TOL = 1.0e-6_real64
    logical, intent(out) :: passed

    type(CPB06) :: cpb
    type(ten_3D2Osym) :: to_test, result_an, result_num
    real(real64) :: sig, tau

    call cpb%init(c11=1.000D0, c12=2.373D0, c13=2.364D0, &
                  c21=2.373D0, c22=1.838D0, c23=1.196D0, &
                  c31=2.364D0, c32=1.196D0, c33=2.444D0, &
                  c44=3.607D0, c55=3.607D0, c66=3.607D0, &
                  k=0.136D0, a=2.000D0                   &
                  )

    sig = 100.0D0
    tau = 50.0D0
    passed = .true.

    ! Uniaxial RD tension (repeated Cauchy eigenvalues)
    call to_test%init((/sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    result_an  = cpb%dstressEq_dstress(to_test)
    result_num = cpb%dstressEq_dstress_numeric(to_test)
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "CPB06 gradient failed at uniaxial RD tension"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Uniaxial RD compression
    call to_test%init((/-sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    result_an  = cpb%dstressEq_dstress(to_test)
    result_num = cpb%dstressEq_dstress_numeric(to_test)
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "CPB06 gradient failed at uniaxial RD compression"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Uniaxial TD tension
    call to_test%init((/0.0D0, sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    result_an  = cpb%dstressEq_dstress(to_test)
    result_num = cpb%dstressEq_dstress_numeric(to_test)
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "CPB06 gradient failed at uniaxial TD tension"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Pure shear
    call to_test%init((/0.0D0, 0.0D0, 0.0D0, tau, 0.0D0, 0.0D0/))
    result_an  = cpb%dstressEq_dstress(to_test)
    result_num = cpb%dstressEq_dstress_numeric(to_test)
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "CPB06 gradient failed at pure shear"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Do not test hydrostatic here: stress_eq = 0, gradient is singular.
end subroutine test_CPB06_stresseq_derivatives

subroutine test_CPB06_stresseq_derivatives2(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    implicit none

    real(real64), parameter :: TOL = 1.0e-6_real64
    logical, intent(out) :: passed

    type(CPB06) :: cpb
    type(ten_3D2Osym) :: to_test
    type(ten_3D4O3sym) :: result_an, result_num
    real(real64) :: sig, tau

    call cpb%init(c11=1.000D0, c12=2.373D0, c13=2.364D0, &
                  c21=2.373D0, c22=1.838D0, c23=1.196D0, &
                  c31=2.364D0, c32=1.196D0, c33=2.444D0, &
                  c44=3.607D0, c55=3.607D0, c66=3.607D0, &
                  k=0.136D0, a=2.000D0                   &
                  )

    sig = 100.0D0
    tau = 50.0D0
    passed = .true.

    ! Uniaxial RD tension (n_sing = 2 on Cauchy stress)
    call to_test%init((/sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    result_an  = cpb%ddstressEq_ddstress(to_test)
    result_num = cpb%ddstressEq_ddstress_numeric(to_test)
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "CPB06 Hessian failed at uniaxial RD tension"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Uniaxial RD compression
    call to_test%init((/-sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    result_an  = cpb%ddstressEq_ddstress(to_test)
    result_num = cpb%ddstressEq_ddstress_numeric(to_test)
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "CPB06 Hessian failed at uniaxial RD compression"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Uniaxial TD tension
    call to_test%init((/0.0D0, sig, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    result_an  = cpb%ddstressEq_ddstress(to_test)
    result_num = cpb%ddstressEq_ddstress_numeric(to_test)
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "CPB06 Hessian failed at uniaxial TD tension"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Pure shear.
    ! derivative2O never writes a(1:3,4:6), so (normal,shear) Voigt
    ! entries stay 0. Compare the filled block only.
    call to_test%init((/0.0D0, 0.0D0, 0.0D0, tau, 0.0D0, 0.0D0/))
    result_an  = cpb%ddstressEq_ddstress(to_test)
    result_num = cpb%ddstressEq_ddstress_numeric(to_test)
    result_an%vals((/9, 13, 14, 16, 17, 18, 19, 20, 21/)) = 0.0D0
    result_num%vals((/9, 13, 14, 16, 17, 18, 19, 20, 21/)) = 0.0D0
    if (.not. result_an%is_approx(result_num, tol=TOL)) then
        print *, "CPB06 Hessian failed at pure shear"
        print *, "  Analytical:", result_an%vals
        print *, "  Numerical :", result_num%vals
        print *, "  |diff|    :", abs(result_an%vals - result_num%vals)
        passed = .false.
    end if

    ! Do not test hydrostatic here: stress_eq = 0, Hessian is singular.
end subroutine test_CPB06_stresseq_derivatives2

subroutine test_CPB06_isotropic_j2_full(passed)
    ! C = I, k = 0, a = 2 must reproduce J2 exactly, including the full Hessian.
    ! States include doubly repeated deviatoric spectra (uniaxial, equibiaxial,
    ! 45 deg, rotated uniaxial); magnitudes cover stresses in Pa, MPa and GPa.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    use test_cpb06_references
    implicit none
    logical, intent(out) :: passed

    real(real64), parameter :: TOL = 1.0e-9_real64
    type(CPB06) :: cpb
    type(ten_3D2Osym) :: st, g
    type(ten_3D4O3sym) :: h
    real(real64) :: base(6,5), mags(6), qe, ge(6), he(21), n(3), s(6)
    real(real64) :: eq, eg, eh
    integer :: i, j

    call cpb%init(c11=1.0D0, c12=0.0D0, c13=0.0D0, &
                  c21=0.0D0, c22=1.0D0, c23=0.0D0, &
                  c31=0.0D0, c32=0.0D0, c33=1.0D0, &
                  c44=1.0D0, c55=1.0D0, c66=1.0D0, &
                  k=0.0D0, a=2.0D0)

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
            g = cpb%dstressEq_dstress(st)
            h = cpb%ddstressEq_ddstress(st)
            eq = abs(cpb%stress_eq(st) - qe)/qe
            eg = rel_err(g%vals, ge)
            eh = rel_err(h%vals, he)
            if (eq > TOL .or. eg > TOL .or. eh > TOL) then
                print *, "CPB06 isotropic J2 check failed: state", i, "magnitude", mags(j)
                print *, "  rel. errors (value, gradient, Hessian):", eq, eg, eh
                passed = .false.
            end if
        end do
    end do
end subroutine test_CPB06_isotropic_j2_full

subroutine test_CPB06_tracepower_reference(passed)
    ! k = 0 and even a: stress_eq = B*tr(Sigma**a)**(1/a), smooth for any spectrum.
    ! Anisotropic maps: Ti-6Al-4V (Tuninetti et al. 2013, Table 2, Wp = 1.857,
    ! with k set to 0 for this check) and a non-symmetric synthetic C.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    use test_cpb06_references
    implicit none
    logical, intent(out) :: passed

    real(real64), parameter :: TOL = 1.0e-8_real64
    type(CPB06) :: cpb
    type(ten_3D2Osym) :: st, g
    type(ten_3D4O3sym) :: h
    real(real64) :: c1(3,3,2), c2(3,2), s(6), d(3), qe, ge(6), ae(6,6)
    real(real64) :: mags(4), states(6,6), eq, eg, eh
    integer :: exps(3), im, ia, is, j

    c1(:,:,1) = reshape([1.000D0, -2.373D0, -2.364D0, &
                         -2.373D0, -1.838D0, 1.196D0, &
                         -2.364D0, 1.196D0, -2.444D0], [3,3])
    c2(:,1) = [3.607D0, 3.607D0, 3.607D0]
    c1(:,:,2) = reshape([1.10D0, 0.05D0, -0.10D0, &
                         -0.20D0, 0.90D0, 0.30D0, &
                         0.15D0, -0.25D0, 1.20D0], [3,3])
    c2(:,2) = [0.8D0, 1.3D0, 1.7D0]
    exps = [2, 4, 8]
    mags = [1.0D-6, 1.0D0, 1.0D3, 2.5D8]

    passed = .true.
    do im = 1, 2
        ! Generic state and pure shear
        states(:,1) = [1.2D0, -0.3D0, 0.4D0, 0.25D0, 0.11D0, -0.17D0]
        states(:,2) = [0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0]
        ! Exactly repeated transformed eigenvalues: d orthogonal to (1,1,1) and to
        ! the difference of two rows of C1, so Sigma_ii = Sigma_jj and no shear.
        d = cross(c1(1,:,im) - c1(2,:,im))
        states(:,3) = [d, 0.0D0, 0.0D0, 0.0D0]/maxval(abs(d))
        d = cross(c1(2,:,im) - c1(3,:,im))
        states(:,4) = [d, 0.0D0, 0.0D0, 0.0D0]/maxval(abs(d))
        ! Nearly repeated (relative perturbations 1e-9 and 1e-6)
        states(:,5) = states(:,3) + 1.0D-9*[0.3D0, -0.7D0, 0.2D0, 0.5D0, -0.4D0, 0.6D0]
        states(:,6) = states(:,4) + 1.0D-6*[-0.2D0, 0.4D0, 0.9D0, -0.3D0, 0.8D0, 0.1D0]
        do ia = 1, size(exps)
            call cpb%init(c11=c1(1,1,im), c12=c1(1,2,im), c13=c1(1,3,im), &
                          c21=c1(2,1,im), c22=c1(2,2,im), c23=c1(2,3,im), &
                          c31=c1(3,1,im), c32=c1(3,2,im), c33=c1(3,3,im), &
                          c44=c2(1,im), c55=c2(2,im), c66=c2(3,im), &
                          k=0.0D0, a=real(exps(ia), real64))
            do is = 1, size(states, 2)
                do j = 1, size(mags)
                    s = mags(j)*states(:,is)
                    st%vals = s
                    call tracepower_reference(c1(:,:,im), c2(:,im), exps(ia), s, qe, ge, ae)
                    g = cpb%dstressEq_dstress(st)
                    h = cpb%ddstressEq_ddstress(st)
                    eq = abs(cpb%stress_eq(st) - qe)/qe
                    eg = rel_err(g%vals, ge)
                    eh = rel_err(reshape(action_matrix(h), [36]), reshape(ae, [36]))
                    if (eq > 1.0D-12 .or. eg > TOL .or. eh > TOL) then
                        print *, "CPB06 tr(Sigma**a) check failed: map", im, "a", exps(ia), &
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
end subroutine test_CPB06_tracepower_reference

subroutine test_CPB06_published_repeated(passed)
    ! Ti-6Al-4V, Tuninetti et al. (2013 preprint, Table 2, row Wp = 48.66 J/cm3):
    ! k = -0.165, a = 2. The state below has a doubly repeated transformed spectrum
    ! at |sigma| ~ 1e3. Expected values: independent 80-digit reference
    ! (value-only central differences with two step sizes, agreement < 1e-16).
    ! The literals below are that reference truncated to real64. The high-precision
    ! computation is not part of this repository, so they cannot be regenerated from
    ! here if the Hessian packing or the rounding ever changes.
    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_exceptions
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use muscle_tensors
    use muscle_yield_cpb06
    use test_cpb06_references
    implicit none
    logical, intent(out) :: passed

    real(real64), parameter :: TOL_REF = 1.0e-8_real64
    real(real64), parameter :: TOL_SCALE = 1.0e-9_real64
    type(CPB06) :: cpb
    type(ten_3D2Osym) :: st, g, g0
    type(ten_3D4O3sym) :: h, h0
    real(real64) :: s(6), qe, ge(6), ae(6,6), lam(5), eg, eh, q
    logical :: flag_zero, flag_invalid, flag_overflow
    integer :: i

    call cpb%init(c11=1.000D0, c12=-2.428D0, c13=-2.920D0, &
                  c21=-2.428D0, c22=1.652D0, c23=-2.236D0, &
                  c31=-2.920D0, c32=-2.236D0, c33=1.003D0, &
                  c44=-3.996D0, c55=-3.996D0, c66=-3.996D0, &
                  k=-0.165D0, a=2.0D0)

    s = [-3.88263584394196585D+02, 8.16175634522563200D+02, -4.27912050128366673D+02, &
         0.0D0, 0.0D0, 0.0D0]
    qe = 1.19771508616519463D+03
    ge = [-5.21124512831125819D-01, 9.79333602114445645D-01, -4.58209089283319826D-01, &
          0.0D0, 0.0D0, 0.0D0]
    ae = 0.0D0
    ae(1:3,1:3) = reshape([ 4.12789900568876519D-04, -1.31554121385814856D-05, -3.99634488430295018D-04, &
                           -1.31554121385814856D-05,  4.19256547452909818D-07,  1.27361555911285750D-05, &
                           -3.99634488430295018D-04,  1.27361555911285750D-05,  3.86898332839166433D-04], &
                          [3,3])
    ae(4,4) = 1.41403169549418292D-03
    ae(5,5) = 1.41403169549418292D-03
    ae(6,6) = 8.29923238267711973D-04

    passed = .true.

    st%vals = s
    q = cpb%stress_eq(st)
    g0 = cpb%dstressEq_dstress(st)
    h0 = cpb%ddstressEq_ddstress(st)
    eg = rel_err(g0%vals, ge)
    eh = rel_err(reshape(action_matrix(h0), [36]), reshape(ae, [36]))
    if (abs(q - qe)/qe > 1.0D-12 .or. eg > TOL_REF .or. eh > TOL_REF) then
        print *, "CPB06 published repeated-spectrum check failed"
        print *, "  rel. errors (value, gradient, Hessian):", abs(q - qe)/qe, eg, eh
        passed = .false.
    end if

    ! Euler identity (homogeneity of degree 1) and pressure insensitivity
    if (abs((st .ddot. g0) - q)/q > 1.0D-12 .or. &
        abs(g0%xx() + g0%yy() + g0%zz()) > 1.0D-12*maxval(abs(g0%vals))) then
        print *, "CPB06 homogeneity/trace identity failed:", (st .ddot. g0) - q, &
                 g0%xx() + g0%yy() + g0%zz()
        passed = .false.
    end if

    ! Gradient is homogeneous of degree 0, Hessian of degree -1: results must not
    ! depend on the stress units.
    lam = [1.0D-9, 1.0D-3, 1.0D0/7.0D0, 1.0D3, 1.0D6]
    do i = 1, size(lam)
        st%vals = lam(i)*s
        g = cpb%dstressEq_dstress(st)
        h = cpb%ddstressEq_ddstress(st)
        eg = rel_err(g%vals, g0%vals)
        eh = rel_err(lam(i)*h%vals, h0%vals)
        if (eg > TOL_SCALE .or. eh > TOL_SCALE) then
            print *, "CPB06 scale invariance failed for factor", lam(i), ":", eg, eh
            passed = .false.
        end if
    end do

    ! Value-only evaluation that used to divide by zero in the eigenvalue solver
    ! (aborted with SIGFPE when traps were enabled).
    st%vals = [-81.61960844304869D0, 42.72804252855544D0, 38.89156591449324D0, &
               0.0D0, 0.0D0, 0.0D0]
    call ieee_set_flag(ieee_all, .false.)
    q = cpb%stress_eq(st)
    g = cpb%dstressEq_dstress(st)
    h = cpb%ddstressEq_ddstress(st)
    call ieee_get_flag(ieee_divide_by_zero, flag_zero)
    call ieee_get_flag(ieee_invalid, flag_invalid)
    call ieee_get_flag(ieee_overflow, flag_overflow)
    call ieee_set_flag(ieee_all, .false.)
    if (flag_zero .or. flag_invalid .or. flag_overflow) then
        print *, "CPB06 raised IEEE exceptions (divide-by-zero, invalid, overflow):", &
                 flag_zero, flag_invalid, flag_overflow
        passed = .false.
    end if
    if (.not. (ieee_is_finite(q) .and. all(ieee_is_finite(g%vals)) .and. &
               all(ieee_is_finite(h%vals)))) then
        print *, "CPB06 returned non-finite results at the repeated-spectrum state"
        passed = .false.
    end if
    if (abs(q - 1.1650051722899815D+02)/q > 1.0D-12) then
        print *, "CPB06 value at the repeated-spectrum state:", q, "expected", 1.1650051722899815D+02
        passed = .false.
    end if
end subroutine test_CPB06_published_repeated
