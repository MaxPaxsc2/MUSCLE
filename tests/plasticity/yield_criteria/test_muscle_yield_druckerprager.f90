! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_yield_dp
    implicit none
    logical :: passed

    call test_dp_tension_uniaxial(passed)
    if (.not. passed) stop 1

    call test_dp_compression_uniaxial(passed)
    if (.not. passed) stop 2

    call test_dp_shear_calibration(passed)
    if (.not. passed) stop 3

    call test_dp_pure_hydrostatic(passed)
    if (.not. passed) stop 4

    call test_dp_tension_uniaxial_comp_abaqus(passed)
    if (.not. passed) stop 5

    call test_dp_compression_uniaxial_comp_abaqus(passed)
    if (.not. passed) stop 6

    call test_dp_biaxial_comp_abaqus(passed)
    if (.not. passed) stop 7

    call test_dp_closed_form(passed)
    if (.not. passed) stop 8

    call test_dp_gradient_identities(passed)
    if (.not. passed) stop 9

    stop 0
end program test_muscle_yield_dp


subroutine test_dp_tension_uniaxial(passed)
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    implicit none
    logical, intent(out) :: passed

    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: stress
    real(real64) :: res, expected, beta, K_ratio
    real(real64), parameter :: EPS = 1.0D-7

    beta = 20.0D0
    K_ratio = 0.8D0   ! (0.778 <= K <= 1.0)
    expected = 150.0D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_TENSION)

    call stress%init(xx=expected, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    passed = abs(res - expected) < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Tension test"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", res - expected
    else
        print *, "PASS: Uniaxial Tension Test"
    end if
end subroutine test_dp_tension_uniaxial


subroutine test_dp_compression_uniaxial(passed)
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    implicit none
    logical, intent(out) :: passed

    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: stress
    real(real64) :: res, applied_stress, expected, beta, K_ratio
    real(real64), parameter :: EPS = 1.0D-7

    beta = 20.0D0
    K_ratio = 0.8D0
    applied_stress = -150.0D0
    expected = 150.0D0 

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_COMPRESSION)

    call stress%init(xx=applied_stress, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    passed = abs(res - expected) < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Compression test"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", res - expected
    else
        print *, "PASS: Uniaxial Compression Test"
    end if
end subroutine test_dp_compression_uniaxial


subroutine test_dp_shear_calibration(passed)
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    implicit none
    logical, intent(out) :: passed

    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: stress
    real(real64) :: res, expected, beta, K_ratio, tau, q_val
    real(real64), parameter :: EPS = 1.0D-7

    beta = 10.0D0
    K_ratio = 0.8D0
    tau = 80.0D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_SHEAR)

    call stress%init(xx=0.0D0, yy=0.0D0, zz=0.0D0, xy=tau, yz=0.0D0, xz=0.0D0)

    q_val = sqrt(3.0D0) * tau
    expected = 0.5D0 * q_val * (1.0D0 + 1.0D0 / K_ratio)

    res = dp%stress_eq(stress)

    passed = abs(res - expected) < EPS
    if (.not. passed) then
        print *, "FAIL: Pure Shear test"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", res - expected
    else
        print *, "PASS: Pure Shear Test"
    end if
end subroutine test_dp_shear_calibration


subroutine test_dp_pure_hydrostatic(passed)
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    implicit none
    logical, intent(out) :: passed

    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: stress
    real(real64) :: res, expected, beta, K_ratio, p_hydro, tanbeta, PI
    real(real64), parameter :: EPS = 1.0D-7

    PI = acos(-1.0D0)
    beta = 30.0D0
    K_ratio = 1.0D0
    p_hydro = 100.0D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_TENSION)

    call stress%init(xx=p_hydro, yy=p_hydro, zz=p_hydro, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    tanbeta = tan(beta * PI / 180.0D0)
    expected = (tanbeta * (1.0D0 + tanbeta / 3.0D0)**(-1)) * p_hydro

    res = dp%stress_eq(stress)

    passed = abs(res - expected) < EPS
    if (.not. passed) then
        print *, "FAIL: Pure Hydrostatic test"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", res - expected
    else
        print *, "PASS: Pure Hydrostatic Test"
    end if
end subroutine test_dp_pure_hydrostatic

subroutine test_dp_tension_uniaxial_comp_abaqus(passed)
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    use muscle_hard_bilinear
    implicit none
    logical, intent(out) :: passed
    logical :: ok

    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: stress, tangent
    type(Bilinear_hardening) :: bilinear
    real(real64) :: res, expected, expected_tan, beta, K_ratio
    real(real64), parameter :: EPS = 1.0D-3

    passed = .true.

    ! ---------------------------------------------------
    ! --------  Tensile type  ---------------------------
    ! ---------------------------------------------------
    bilinear = Bilinear_hardening(y0=300D0, k=1000D0)
    beta = 16.0D0
    K_ratio = 0.85D0   ! (0.778 <= K <= 1.0)
    expected = bilinear%stress(0.09344D0)
    expected_tan = -2.581929D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_TENSION)
    call stress%init(xx=393.4D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    ok = abs(res - expected)/expected < EPS
    if (.not. ok) then
        print *, "FAIL: Uniaxial Tension test (TENSIONDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Uniaxial Tension Test (TENSIONDEF)"
    end if
    passed = passed .and. ok
    
    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    ok = abs(res - expected_tan) < EPS
    if (.not. ok) then
        print *, "FAIL: Uniaxial Tangent test (TENSIONDEF)"
        print *, "  Expected:", expected_tan
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan)
    else
        print *, "PASS: Uniaxial Tangent Test (TENSIONDEF)"
    end if
    passed = passed .and. ok

    ! ---------------------------------------------------
    ! --------  Compression type  -----------------------
    ! ---------------------------------------------------
    bilinear = Bilinear_hardening(y0=300D0, k=1000D0)
    beta = 16.0D0
    K_ratio = 0.85D0   ! (0.778 <= K <= 1.0)
    expected = bilinear%stress(0.06688D0)
    expected_tan = -2.58221D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_COMPRESSION)
    call stress%init(xx=260.8D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    ok = abs(res - expected)/expected < EPS
    if (.not. ok) then
        print *, "FAIL: Uniaxial Tension test (COMPRESSIONDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Uniaxial Tension Test (COMPRESSIONDEF)"
    end if
    passed = passed .and. ok

    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    ok = abs(res - expected_tan) < EPS
    if (.not. ok) then
        print *, "FAIL: Uniaxial Tangent test (COMPRESSIONDEF)"
        print *, "  Expected:", expected_tan
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan)
    else
        print *, "PASS: Uniaxial Tangent Test (COMPRESSIONDEF)"
    end if
    passed = passed .and. ok

    ! ---------------------------------------------------
    ! --------  Shear type  -----------------------------
    ! ---------------------------------------------------
    bilinear = Bilinear_hardening(y0=300D0, k=1000D0)
    beta = 16.0D0
    K_ratio = 0.85D0   ! (0.778 <= K <= 1.0)
    expected = bilinear%stress(0.08032D0)
    expected_tan = -2.58221D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_SHEAR)
    call stress%init(xx=299.0D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    ok = abs(res - expected)/expected < EPS
    if (.not. ok) then
        print *, "FAIL: Uniaxial Tension test (SHEARDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Uniaxial Tension Test (SHEARDEF)"
    end if
    passed = passed .and. ok

    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    ok = abs(res - expected_tan) < EPS
    if (.not. ok) then
        print *, "FAIL: Uniaxial Tangent test (SHEARDEF)"
        print *, "  Expected:", expected_tan
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan)
    else
        print *, "PASS: Uniaxial Tangent Test (SHEARDEF)"
    end if
    passed = passed .and. ok
end subroutine test_dp_tension_uniaxial_comp_abaqus

subroutine test_dp_compression_uniaxial_comp_abaqus(passed)
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    use muscle_hard_bilinear
    implicit none
    logical, intent(out) :: passed
    logical :: ok

    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: stress, tangent
    type(Bilinear_hardening) :: bilinear
    real(real64) :: res, expected, expected_tan, beta, K_ratio
    real(real64), parameter :: EPS = 1.0D-3

    passed = .true.

    bilinear = Bilinear_hardening(y0=300D0, k=1000D0)
    beta = 16.0D0
    K_ratio = 0.85D0   ! (0.778 <= K <= 1.0)
    ! ---------------------------------------------------
    ! --------  Tensile type  ---------------------------
    ! ---------------------------------------------------
    expected = bilinear%stress(0.0112D0)
    expected_tan = -1.5180952D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_TENSION)
    call stress%init(xx=-437.7D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    ok = abs(res - expected)/expected < EPS
    if (.not. ok) then
        print *, "FAIL: Uniaxial Compression test (TENSIONDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Uniaxial Compression Test (TENSIONDEF)"
    end if
    passed = passed .and. ok
    
    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    ok = abs(res - expected_tan) < EPS
    if (.not. ok) then
        print *, "FAIL: Uniaxial Compression Tangent test (TENSIONDEF)"
        print *, "  Expected:", expected_tan
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan)
    else
        print *, "PASS: Uniaxial Compression Tangent Test (TENSIONDEF)"
    end if
    passed = passed .and. ok

    ! ---------------------------------------------------
    ! --------  Compression type  -----------------------
    ! ---------------------------------------------------
    expected = bilinear%stress(0.00858D0)
    expected_tan = -1.518584D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_COMPRESSION)
    call stress%init(xx=-308.58D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    ok = abs(res - expected)/expected < EPS
    if (.not. ok) then
        print *, "FAIL: Uniaxial Compression test (COMPRESSIONDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Uniaxial Compression Test (COMPRESSIONDEF)"
    end if
    passed = passed .and. ok

    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    ok = abs(res - expected_tan) < EPS
    if (.not. ok) then
        print *, "FAIL: Uniaxial Compression Tangent test (COMPRESSIONDEF)"
        print *, "  Expected:", expected_tan
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan)
    else
        print *, "PASS: Uniaxial Compression Tangent Test (COMPRESSIONDEF)"
    end if
    passed = passed .and. ok

    ! ---------------------------------------------------
    ! --------  Shear type  -----------------------------
    ! ---------------------------------------------------
    expected = bilinear%stress(0.01013D0)
    expected_tan = -1.5198555D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_SHEAR)
    call stress%init(xx=-342.9D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    ok = abs(res - expected)/expected < EPS
    if (.not. ok) then
        print *, "FAIL: Uniaxial Compression test (SHEARDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Uniaxial Compression Test (SHEARDEF)"
    end if
    passed = passed .and. ok

    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    ok = abs(res - expected_tan) < EPS*1.5 ! Falta de decimales
    if (.not. ok) then
        print *, "FAIL: Uniaxial Compression Tangent test (SHEARDEF)"
        print *, "  Expected:", expected_tan
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan)
    else
        print *, "PASS: Uniaxial Compression Tangent Test (SHEARDEF)"
    end if
    passed = passed .and. ok
end subroutine test_dp_compression_uniaxial_comp_abaqus


subroutine test_dp_biaxial_comp_abaqus(passed)
    ! Abaqus reference for a non-proportional (constrained) loading path.
    !
    ! The stress components and PEEQ below are the final state of the Abaqus run, so the
    ! equivalent stress is asserted against the hardening curve at that PEEQ.
    !
    ! The expected_tan1/expected_tan2 values are ratios of ACCUMULATED plastic strain
    ! components (PE11/PE22, PE11/PE33) printed by Abaqus with 3-4 digits, e.g.
    ! 0.0169/0.0303 = 0.55775577. Along a non-proportional path PE is the sum of the flow
    ! directions of all increments, so it is not the flow direction N at the final state.
    ! These ratios coincide with the N of this same criterion evaluated about 0.2-0.25 deg
    ! (Lode angle) before the final state, and differ from N at the final state by 1.3-2.6 %.
    ! They are printed for information only and are not asserted; the gradient at these
    ! states is checked against finite differences in test_dp_gradient_identities.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    use muscle_hard_bilinear
    implicit none
    logical, intent(out) :: passed
    logical :: ok

    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: stress, tangent
    type(Bilinear_hardening) :: bilinear
    real(real64) :: res, expected, expected_tan1, expected_tan2, beta, K_ratio
    real(real64), parameter :: EPS = 1.0D-4

    passed = .true.

    bilinear = Bilinear_hardening(y0=300D0, k=1000D0)
    beta = 16.0D0
    K_ratio = 0.85D0   ! (0.778 <= K <= 1.0)

    ! ---------------------------------------------------
    ! --------  Tensile type  ---------------------------
    ! ---------------------------------------------------
    expected = bilinear%stress(0.03497D0)
    expected_tan1 = -0.471306D0
    expected_tan2 = 0.55486D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_TENSION)
    call stress%init(xx=-994.63D0, yy=-1711.35D0, zz=-901.47D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    ok = abs(res - expected)/expected < EPS
    if (.not. ok) then
        print *, "FAIL: Biaxial test (TENSIONDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Biaxial Test (TENSIONDEF)"
    end if
    passed = passed .and. ok

    tangent = dp%dstressEq_dstress(stress)
    call print_pe_ratio("Biaxial Nxx/Nyy (TENSIONDEF)", tangent%xx()/tangent%yy(), expected_tan1)
    call print_pe_ratio("Biaxial Nxx/Nzz (TENSIONDEF)", tangent%xx()/tangent%zz(), expected_tan2)

    ! ---------------------------------------------------
    ! --------  Compression type  -----------------------
    ! ---------------------------------------------------
    expected = bilinear%stress(0.02531D0)
    expected_tan1 = -0.472822D0
    expected_tan2 = 0.55775577D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_COMPRESSION)
    call stress%init(xx=-1048.47D0, yy=-1647.97D0, zz=-970.69D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    ok = abs(res - expected)/expected < EPS
    if (.not. ok) then
        print *, "FAIL: Biaxial test (COMPRESSIONDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Biaxial Test (COMPRESSIONDEF)"
    end if
    passed = passed .and. ok

    tangent = dp%dstressEq_dstress(stress)
    call print_pe_ratio("Biaxial Nxx/Nyy (COMPRESSIONDEF)", tangent%xx()/tangent%yy(), expected_tan1)
    call print_pe_ratio("Biaxial Nxx/Nzz (COMPRESSIONDEF)", tangent%xx()/tangent%zz(), expected_tan2)

    ! ---------------------------------------------------
    ! --------  Shear type  -----------------------------
    ! ---------------------------------------------------
    expected = bilinear%stress(0.03031D0)
    expected_tan1 = -0.4724789D0
    expected_tan2 = 0.55682D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_SHEAR)
    call stress%init(xx=-1033.74D0, yy=-1665.32D0, zz=-951.75D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    ok = abs(res - expected)/expected < EPS
    if (.not. ok) then
        print *, "FAIL: Biaxial test (SHEARDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Biaxial Test (SHEARDEF)"
    end if
    passed = passed .and. ok

    tangent = dp%dstressEq_dstress(stress)
    call print_pe_ratio("Biaxial Nxx/Nyy (SHEARDEF)", tangent%xx()/tangent%yy(), expected_tan1)
    call print_pe_ratio("Biaxial Nxx/Nzz (SHEARDEF)", tangent%xx()/tangent%zz(), expected_tan2)

contains

    subroutine print_pe_ratio(label, got, abaqus_pe)
        character(len=*), intent(in) :: label
        real(real64), intent(in) :: got, abaqus_pe

        print *, "INFO: ", label, " (not asserted, Abaqus value is accumulated PE)"
        print *, "  Abaqus PE ratio:", abaqus_pe
        print *, "  Code N ratio   :", got
        print *, "  Rel. diff      :", abs(got - abaqus_pe)/abs(abaqus_pe)
    end subroutine print_pe_ratio

end subroutine test_dp_biaxial_comp_abaqus


subroutine test_dp_closed_form(passed)
    ! Closed-form values derived by hand from the formulation documented in
    ! muscle_yield_druckerprager (beta = 16 deg, K = 0.85, all hardening modes):
    !
    !   sigma_eq = f*[ q/2*(1 + 1/K - (1 - 1/K)*xi) + tan(beta)*p ],  xi = (r/q)**3
    !   N        = d(sigma_eq)/d(sigma)
    !
    ! xi = +1 for uniaxial/triaxial tension, -1 for uniaxial/triaxial compression and
    ! equibiaxial tension, 0 for pure shear, so every value reduces to an explicit
    ! expression in tb = tan(beta), K and the hardening-mode factor f.
    ! With c1 = (1 + 1/K)/2 and c2 = (1 - 1/K)/2:
    !
    !   uniaxial tension  (s,0,0):  sigma_eq = f*s*(1/K + tb/3)
    !                               N = f*(1/K + tb/3, -1/(2K) + tb/3, -1/(2K) + tb/3)
    !                               Nxx/Nyy = -2*(3 + K*tb)/(3 - 2*K*tb)
    !   uniaxial compr.  (-s,0,0):  sigma_eq = f*s*(1 - tb/3)
    !                               N = f*(tb/3 - 1, 1/2 + tb/3, 1/2 + tb/3)
    !                               Nxx/Nyy = -2*(3 - tb)/(3 + 2*tb)
    !   triaxial tension  (a,b,b):  sigma_eq = f*((a - b)/K + tb*(a + 2b)/3), N as uniaxial tension
    !   triaxial compr.   (b,a,a):  sigma_eq = f*((a - b) + tb*(2a + b)/3),   N as uniaxial compr.
    !   pure shear        xy = t:   sigma_eq = f*sqrt(3)*c1*t
    !                               N = f*(tb/3 - 3c2/2, tb/3 - 3c2/2, tb/3 + 3c2; xy = sqrt(3)*c1/2)
    !   equibiaxial       (s,s,0):  sigma_eq = f*s*(1 + 2*tb/3)
    !                               N = f*(1/2 + tb/3, 1/2 + tb/3, tb/3 - 1)
    !   plane shear      (t,-t,0):  sigma_eq = f*sqrt(3)*c1*t
    !                               N = f*(sqrt(3)*c1/2 + tb/3 - 3c2/2, -sqrt(3)*c1/2 + tb/3 - 3c2/2, tb/3 + 3c2)
    !
    ! N components are tensor components (xx, yy, zz, xy, yz, xz); those not listed are zero.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    implicit none
    logical, intent(out) :: passed

    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: stress, tangent
    real(real64) :: tb, c1, c2, f, fm(3), s, t, a, b
    real(real64) :: n_ten(6), n_comp(6), n_shear(6), n_equi(6), n_plane(6)
    integer :: m
    logical :: mode_passed
    real(real64), parameter :: EPS = 1.0D-12
    real(real64), parameter :: PI = acos(-1.0D0)
    real(real64), parameter :: BETA = 16.0D0, K_RATIO = 0.85D0
    integer, parameter :: MODES(3) = [DP_HARDENING_TENSION, DP_HARDENING_COMPRESSION, DP_HARDENING_SHEAR]
    character(len=14), parameter :: NAMES(3) = [character(len=14) :: "TENSIONDEF", "COMPRESSIONDEF", "SHEARDEF"]

    passed = .true.

    tb = tan(BETA * PI / 180.0D0)
    c1 = 0.5D0 * (1.0D0 + 1.0D0 / K_RATIO)
    c2 = 0.5D0 * (1.0D0 - 1.0D0 / K_RATIO)
    fm = [1.0D0 / (1.0D0 / K_RATIO + tb / 3.0D0), 1.0D0 / (1.0D0 - tb / 3.0D0), 1.0D0]

    s = 300.0D0
    t = 100.0D0
    a = 200.0D0
    b = 50.0D0

    do m = 1, 3
        mode_passed = .true.
        f = fm(m)
        call dp%init(beta_deg=BETA, K=K_RATIO, hardening_mode=MODES(m))

        n_ten   = f * [1.0D0/K_RATIO + tb/3.0D0, -0.5D0/K_RATIO + tb/3.0D0, -0.5D0/K_RATIO + tb/3.0D0, &
                       0.0D0, 0.0D0, 0.0D0]
        n_comp  = f * [tb/3.0D0 - 1.0D0, 0.5D0 + tb/3.0D0, 0.5D0 + tb/3.0D0, 0.0D0, 0.0D0, 0.0D0]
        n_shear = f * [tb/3.0D0 - 1.5D0*c2, tb/3.0D0 - 1.5D0*c2, tb/3.0D0 + 3.0D0*c2, &
                       0.5D0*sqrt(3.0D0)*c1, 0.0D0, 0.0D0]
        n_equi  = f * [0.5D0 + tb/3.0D0, 0.5D0 + tb/3.0D0, tb/3.0D0 - 1.0D0, 0.0D0, 0.0D0, 0.0D0]
        n_plane = f * [0.5D0*sqrt(3.0D0)*c1 + tb/3.0D0 - 1.5D0*c2, &
                       -0.5D0*sqrt(3.0D0)*c1 + tb/3.0D0 - 1.5D0*c2, &
                       tb/3.0D0 + 3.0D0*c2, 0.0D0, 0.0D0, 0.0D0]

        ! Uniaxial tension
        call stress%init(xx=s, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
        tangent = dp%dstressEq_dstress(stress)
        call check("uniaxial tension sigma_eq", dp%stress_eq(stress), f*s*(1.0D0/K_RATIO + tb/3.0D0))
        call check_n("uniaxial tension N", tangent, n_ten)
        call check("uniaxial tension Nxx/Nyy", tangent%xx()/tangent%yy(), &
                   -2.0D0*(3.0D0 + K_RATIO*tb)/(3.0D0 - 2.0D0*K_RATIO*tb))

        ! Uniaxial compression
        call stress%init(xx=-s, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
        tangent = dp%dstressEq_dstress(stress)
        call check("uniaxial compression sigma_eq", dp%stress_eq(stress), f*s*(1.0D0 - tb/3.0D0))
        call check_n("uniaxial compression N", tangent, n_comp)
        call check("uniaxial compression Nxx/Nyy", tangent%xx()/tangent%yy(), &
                   -2.0D0*(3.0D0 - tb)/(3.0D0 + 2.0D0*tb))

        ! Triaxial tension and compression (non-zero pressure)
        call stress%init(xx=a, yy=b, zz=b, xy=0.0D0, yz=0.0D0, xz=0.0D0)
        call check("triaxial tension sigma_eq", dp%stress_eq(stress), &
                   f*((a - b)/K_RATIO + tb*(a + 2.0D0*b)/3.0D0))
        call check_n("triaxial tension N", dp%dstressEq_dstress(stress), n_ten)

        call stress%init(xx=b, yy=a, zz=a, xy=0.0D0, yz=0.0D0, xz=0.0D0)
        call check("triaxial compression sigma_eq", dp%stress_eq(stress), &
                   f*((a - b) + tb*(2.0D0*a + b)/3.0D0))
        call check_n("triaxial compression N", dp%dstressEq_dstress(stress), n_comp)

        ! Pure shear
        call stress%init(xx=0.0D0, yy=0.0D0, zz=0.0D0, xy=t, yz=0.0D0, xz=0.0D0)
        call check("pure shear sigma_eq", dp%stress_eq(stress), f*sqrt(3.0D0)*c1*t)
        call check_n("pure shear N", dp%dstressEq_dstress(stress), n_shear)

        ! Equibiaxial tension
        call stress%init(xx=s, yy=s, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
        call check("equibiaxial sigma_eq", dp%stress_eq(stress), f*s*(1.0D0 + 2.0D0*tb/3.0D0))
        call check_n("equibiaxial N", dp%dstressEq_dstress(stress), n_equi)

        ! Plane shear in principal axes
        call stress%init(xx=t, yy=-t, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
        call check("plane shear sigma_eq", dp%stress_eq(stress), f*sqrt(3.0D0)*c1*t)
        call check_n("plane shear N", dp%dstressEq_dstress(stress), n_plane)

        if (mode_passed) print *, "PASS: Closed-form values (", trim(NAMES(m)), ")"
        passed = passed .and. mode_passed
    end do

contains

    subroutine check(label, got, want)
        character(len=*), intent(in) :: label
        real(real64), intent(in) :: got, want
        logical :: ok

        ok = abs(got - want) <= EPS * max(abs(want), 1.0D0)
        if (.not. ok) then
            print *, "FAIL: Closed-form ", label, " (", trim(NAMES(m)), ")"
            print *, "  Expected:", want
            print *, "  Got     :", got
            print *, "  Diff    :", got - want
        end if
        mode_passed = mode_passed .and. ok
    end subroutine check

    subroutine check_n(label, n, want)
        character(len=*), intent(in) :: label
        type(ten_3D2Osym), intent(in) :: n
        real(real64), intent(in) :: want(6)

        call check(label // " xx", n%xx(), want(1))
        call check(label // " yy", n%yy(), want(2))
        call check(label // " zz", n%zz(), want(3))
        call check(label // " xy", n%xy(), want(4))
        call check(label // " yz", n%yz(), want(5))
        call check(label // " xz", n%xz(), want(6))
    end subroutine check_n

end subroutine test_dp_closed_form


subroutine test_dp_gradient_identities(passed)
    ! Analytical gradient versus the finite-difference gradient of stress_eq, and the
    ! homogeneity (Euler) identity N:sigma = sigma_eq, at states with and without shear,
    ! including the three Abaqus biaxial final states.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    implicit none
    logical, intent(out) :: passed

    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: stress, analytic, numeric, diff
    real(real64) :: states(6, 6), res_fd, res_euler, seq
    integer :: m, i
    logical :: ok, mode_passed
    real(real64), parameter :: EPS_FD = 1.0D-8
    real(real64), parameter :: EPS_EULER = 1.0D-12
    integer, parameter :: MODES(3) = [DP_HARDENING_TENSION, DP_HARDENING_COMPRESSION, DP_HARDENING_SHEAR]
    character(len=14), parameter :: NAMES(3) = [character(len=14) :: "TENSIONDEF", "COMPRESSIONDEF", "SHEARDEF"]

    states(:, 1) = [-994.63D0, -1711.35D0, -901.47D0, 0.0D0, 0.0D0, 0.0D0]
    states(:, 2) = [-1048.47D0, -1647.97D0, -970.69D0, 0.0D0, 0.0D0, 0.0D0]
    states(:, 3) = [-1033.74D0, -1665.32D0, -951.75D0, 0.0D0, 0.0D0, 0.0D0]
    states(:, 4) = [120.0D0, -40.0D0, 75.0D0, 30.0D0, -22.0D0, 15.0D0]
    states(:, 5) = [-250.0D0, 60.0D0, -30.0D0, -45.0D0, 90.0D0, 70.0D0]
    states(:, 6) = [0.0D0, 0.0D0, 0.0D0, 80.0D0, 0.0D0, 0.0D0]

    passed = .true.

    do m = 1, 3
        mode_passed = .true.
        call dp%init(beta_deg=16.0D0, K=0.85D0, hardening_mode=MODES(m))
        do i = 1, size(states, 2)
            call stress%init(xx=states(1, i), yy=states(2, i), zz=states(3, i), &
                             xy=states(4, i), yz=states(5, i), xz=states(6, i))
            seq = dp%stress_eq(stress)
            analytic = dp%dstressEq_dstress(stress)
            numeric = dp%dstressEq_dstress_numeric(stress)
            diff = analytic - numeric

            res_fd = sqrt(diff .ddot. diff) / sqrt(analytic .ddot. analytic)
            res_euler = abs((analytic .ddot. stress) - seq) / abs(seq)

            ok = res_fd < EPS_FD .and. res_euler < EPS_EULER
            if (.not. ok) then
                print *, "FAIL: Gradient identities (", trim(NAMES(m)), ") state", i
                print *, "  |N - N_FD|/|N|      :", res_fd
                print *, "  |N:sigma - seq|/seq :", res_euler
            end if
            mode_passed = mode_passed .and. ok
        end do
        if (mode_passed) print *, "PASS: Gradient identities (", trim(NAMES(m)), ")"
        passed = passed .and. mode_passed
    end do

end subroutine test_dp_gradient_identities
