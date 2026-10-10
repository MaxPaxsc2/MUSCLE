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

    call test_DruckerPrager_hessian(passed)
    if (.not. passed) stop 8

    call test_dp_apex_slope(passed)
    if (.not. passed) stop 9

    call test_dp_apex_dev_gauge(passed)
    if (.not. passed) stop 10

    call test_dp_biaxial_comp_abaqus(passed)
    if (.not. passed) stop 7

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

    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: stress, tangent
    type(Bilinear_hardening) :: bilinear
    real(real64) :: res, expected, expected_tan, beta, K_ratio
    real(real64), parameter :: EPS = 1.0D-3

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

    passed = abs(res - expected)/expected < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Tension test (TENSIONDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Uniaxial Tension Test (TENSIONDEF)"
    end if
    
    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    passed = abs(res - expected_tan) < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Tangent test (TENSIONDEF)"
        print *, "  Expected:", expected_tan
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan)
    else
        print *, "PASS: Uniaxial Tangent Test (TENSIONDEF)"
    end if

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

    passed = abs(res - expected)/expected < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Tension test (COMPRESSIONDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Uniaxial Tension Test (COMPRESSIONDEF)"
    end if

    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    passed = abs(res - expected_tan) < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Tangent test (COMPRESSIONDEF)"
        print *, "  Expected:", expected_tan
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan)
    else
        print *, "PASS: Uniaxial Tangent Test (COMPRESSIONDEF)"
    end if

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

    passed = abs(res - expected)/expected < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Tension test (SHEARDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Uniaxial Tension Test (SHEARDEF)"
    end if

    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    passed = abs(res - expected_tan) < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Tangent test (SHEARDEF)"
        print *, "  Expected:", expected_tan
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan)
    else
        print *, "PASS: Uniaxial Tangent Test (SHEARDEF)"
    end if
end subroutine test_dp_tension_uniaxial_comp_abaqus

subroutine test_dp_compression_uniaxial_comp_abaqus(passed)
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    use muscle_hard_bilinear
    implicit none
    logical, intent(out) :: passed

    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: stress, tangent
    type(Bilinear_hardening) :: bilinear
    real(real64) :: res, expected, expected_tan, beta, K_ratio
    real(real64), parameter :: EPS = 1.0D-3

    bilinear = Bilinear_hardening(y0=300D0, k=1000D0)
    beta = 16.0D0
    K_ratio = 0.85D0   ! (0.778 <= K <= 1.0)
    ! ---------------------------------------------------
    ! --------  Tensile type  ---------------------------
    ! ---------------------------------------------------
    expected = bilinear%stress(0.0112D0)
    expected_tan = -1.5180952

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_TENSION)
    call stress%init(xx=-437.7D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    passed = abs(res - expected)/expected < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Compression test (TENSIONDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Uniaxial Compression Test (TENSIONDEF)"
    end if
    
    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    passed = abs(res - expected_tan) < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Compression Tangent test (TENSIONDEF)"
        print *, "  Expected:", expected_tan
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan)
    else
        print *, "PASS: Uniaxial Compression Tangent Test (TENSIONDEF)"
    end if

    ! ---------------------------------------------------
    ! --------  Compression type  -----------------------
    ! ---------------------------------------------------
    expected = bilinear%stress(0.00858D0)
    expected_tan = -1.518584D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_COMPRESSION)
    call stress%init(xx=-308.58D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    passed = abs(res - expected)/expected < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Compression test (COMPRESSIONDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Uniaxial Compression Test (COMPRESSIONDEF)"
    end if

    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    passed = abs(res - expected_tan) < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Compression Tangent test (COMPRESSIONDEF)"
        print *, "  Expected:", expected_tan
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan)
    else
        print *, "PASS: Uniaxial Compression Tangent Test (COMPRESSIONDEF)"
    end if

    ! ---------------------------------------------------
    ! --------  Shear type  -----------------------------
    ! ---------------------------------------------------
    expected = bilinear%stress(0.01013D0)
    expected_tan = -1.5198555

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_SHEAR)
    call stress%init(xx=-342.9D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    passed = abs(res - expected)/expected < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Compression test (SHEARDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Uniaxial Compression Test (SHEARDEF)"
    end if

    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    passed = abs(res - expected_tan) < EPS*1.5 ! Falta de decimales
    if (.not. passed) then
        print *, "FAIL: Uniaxial Compression Tangent test (SHEARDEF)"
        print *, "  Expected:", expected_tan
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan)
    else
        print *, "PASS: Uniaxial Compression Tangent Test (SHEARDEF)"
    end if
end subroutine test_dp_compression_uniaxial_comp_abaqus


subroutine test_dp_biaxial_comp_abaqus(passed)
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    use muscle_hard_bilinear
    implicit none
    logical, intent(out) :: passed

    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: stress, tangent
    type(Bilinear_hardening) :: bilinear
    real(real64) :: res, expected, expected_tan1, expected_tan2, beta, K_ratio
    real(real64), parameter :: EPS = 1.0D-4

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

    passed = abs(res - expected)/expected < EPS
    if (.not. passed) then
        print *, "FAIL: Biaxial test (TENSIONDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Biaxial Test (TENSIONDEF)"
    end if
    
    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    passed = abs(res - expected_tan1)/expected_tan1 < EPS
    if (.not. passed) then
        print *, "FAIL: Biaxial Tangent1 test (TENSIONDEF)"
        print *, "  Expected:", expected_tan1
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan1)
    else
        print *, "PASS: Biaxial Tangent1 Test (TENSIONDEF)"
    end if

    res = tangent%xx()/tangent%zz()
    passed = abs(res - expected_tan2)/expected_tan2 < EPS
    if (.not. passed) then
        print *, "FAIL: Biaxial Tangent2 test (TENSIONDEF)"
        print *, "  Expected:", expected_tan2
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan2)
    else
        print *, "PASS: Biaxial Tangent2 Test (TENSIONDEF)"
    end if

    ! ---------------------------------------------------
    ! --------  Compression type  -----------------------
    ! ---------------------------------------------------
    expected = bilinear%stress(0.02531D0)
    expected_tan1 = -0.472822D0
    expected_tan2 = 0.55775577D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_COMPRESSION)
    call stress%init(xx=-1048.47D0, yy=-1647.97D0, zz=-970.69D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    passed = abs(res - expected)/expected < EPS
    if (.not. passed) then
        print *, "FAIL: Biaxial test (COMPRESSIONDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Biaxial Test (COMPRESSIONDEF)"
    end if
    
    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    passed = abs(res - expected_tan1)/expected_tan1 < EPS
    if (.not. passed) then
        print *, "FAIL: Biaxial Tangent1 test (COMPRESSIONDEF)"
        print *, "  Expected:", expected_tan1
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan1)
    else
        print *, "PASS: Biaxial Tangent1 Test (COMPRESSIONDEF)"
    end if

    res = tangent%xx()/tangent%zz()
    passed = abs(res - expected_tan2)/expected_tan2 < EPS
    if (.not. passed) then
        print *, "FAIL: Biaxial Tangent2 test (COMPRESSIONDEF)"
        print *, "  Expected:", expected_tan2
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan2)
    else
        print *, "PASS: Biaxial Tangent2 Test (COMPRESSIONDEF)"
    end if

    ! ! ---------------------------------------------------
    ! ! --------  Shear type  -----------------------------
    ! ! ---------------------------------------------------
    expected = bilinear%stress(0.03031D0)
    expected_tan1 = -0.4724789D0
    expected_tan2 = 0.55682D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_SHEAR)
    call stress%init(xx=-1033.74D0, yy=-1665.32D0, zz=-951.75D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    passed = abs(res - expected)/expected < EPS
    if (.not. passed) then
        print *, "FAIL: Biaxial test (COMPRESSIONDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Biaxial Test (COMPRESSIONDEF)"
    end if
    
    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    passed = abs(res - expected_tan1)/expected_tan1 < EPS
    if (.not. passed) then
        print *, "FAIL: Biaxial Tangent1 test (COMPRESSIONDEF)"
        print *, "  Expected:", expected_tan1
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan1)
    else
        print *, "PASS: Biaxial Tangent1 Test (COMPRESSIONDEF)"
    end if

    res = tangent%xx()/tangent%zz()
    passed = abs(res - expected_tan2)/expected_tan2 < EPS
    if (.not. passed) then
        print *, "FAIL: Biaxial Tangent2 test (COMPRESSIONDEF)"
        print *, "  Expected:", expected_tan2
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan2)
    else
        print *, "PASS: Biaxial Tangent2 Test (COMPRESSIONDEF)"
    end if

end subroutine test_dp_biaxial_comp_abaqus

subroutine test_DruckerPrager_hessian(passed)
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    use muscle_yield_vonmises
    implicit none
    logical, intent(out) :: passed

    type(DruckerPrager) :: dp
    type(VonMises) :: vm
    type(ten_3D2Osym) :: stress
    type(ten_3D4O3sym) :: hess
    real(real64), parameter :: TOL_VM = 1.0D-12
    real(real64), parameter :: TOL_NUM = 1.0D-6

    ! K = 1 removes the third-invariant term and the shear calibration gives f = 1:
    ! the Hessian must equal the Von Mises one. Normals and shears coupled.
    call dp%init(beta_deg=16.0D0, K=1.0D0, hardening_mode=DP_HARDENING_SHEAR)
    call stress%init(xx=210.0D0, yy=-45.0D0, zz=80.0D0, xy=60.0D0, yz=-35.0D0, xz=25.0D0)
    hess = dp%ddstressEq_ddstress(stress)
    passed = hess%is_approx(vm%ddstressEq_ddstress(stress), tol=TOL_VM)
    if (.not. passed) then
        print *, "FAIL: Hessian K = 1 vs Von Mises"
        return
    end if
    print *, "PASS: Hessian K = 1 vs Von Mises"

    ! K = 0.78 (lower bound) makes the third-invariant term largest; tension-like Lode angle
    ! (r^3 > 0) with all six components nonzero.
    call dp%init(beta_deg=16.0D0, K=0.78D0, hardening_mode=DP_HARDENING_TENSION)
    call stress%init(xx=300.0D0, yy=40.0D0, zz=-20.0D0, xy=45.0D0, yz=-30.0D0, xz=55.0D0)
    hess = dp%ddstressEq_ddstress(stress)
    passed = hess%is_approx(dp%ddstressEq_ddstress_numeric(stress), tol=TOL_NUM)
    if (.not. passed) then
        print *, "FAIL: Hessian vs numeric (TENSIONDEF, K = 0.78)"
        return
    end if
    print *, "PASS: Hessian vs numeric (TENSIONDEF, K = 0.78)"

    ! Compression-like Lode angle (r^3 < 0) on top of a large hydrostatic part, as in the
    ! biaxial case above, with shears added.
    call dp%init(beta_deg=16.0D0, K=0.85D0, hardening_mode=DP_HARDENING_COMPRESSION)
    call stress%init(xx=-1048.47D0, yy=-1647.97D0, zz=-970.69D0, xy=-120.0D0, yz=85.0D0, xz=40.0D0)
    hess = dp%ddstressEq_ddstress(stress)
    passed = hess%is_approx(dp%ddstressEq_ddstress_numeric(stress), tol=TOL_NUM)
    if (.not. passed) then
        print *, "FAIL: Hessian vs numeric (COMPRESSIONDEF, K = 0.85)"
        return
    end if
    print *, "PASS: Hessian vs numeric (COMPRESSIONDEF, K = 0.85)"

    ! Shear calibration with K < 1 activates the third-invariant term missing from the K = 1 check.
    call dp%init(beta_deg=16.0D0, K=0.78D0, hardening_mode=DP_HARDENING_SHEAR)
    call stress%init(xx=300.0D0, yy=40.0D0, zz=-20.0D0, xy=45.0D0, yz=-30.0D0, xz=55.0D0)
    hess = dp%ddstressEq_ddstress(stress)
    passed = hess%is_approx(dp%ddstressEq_ddstress_numeric(stress), tol=TOL_NUM)
    if (.not. passed) then
        print *, "FAIL: Hessian vs numeric (SHEARDEF, K = 0.78)"
        return
    end if
    print *, "PASS: Hessian vs numeric (SHEARDEF, K = 0.78)"

    ! Pure hydrostatic stress (q = 0): the gradient is constant there, so the Hessian is
    ! exactly zero and must not be NaN (the return mapping reaches this state in tension).
    call stress%init(xx=500.0D0, yy=500.0D0, zz=500.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    hess = dp%ddstressEq_ddstress(stress)
    passed = (sum(abs(hess%vals)) < tiny(1.0D0))
    if (.not. passed) then
        print *, "FAIL: Hessian on the hydrostatic axis"
        return
    end if
    print *, "PASS: Hessian on the hydrostatic axis"
end subroutine test_DruckerPrager_hessian

subroutine test_dp_apex_slope(passed)
    ! The apex is the hydrostatic singularity of the cone: the slope must match the
    ! TENSION calibration of d = f*tan(beta) and the value of stress_eq at p = 1.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    implicit none
    logical, intent(out) :: passed

    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: stress
    real(real64) :: res, expected, tanbeta, PI
    real(real64), parameter :: EPS = 1.0D-12

    PI = acos(-1.0D0)
    tanbeta = tan(16.0D0 * PI / 180.0D0)
    expected = tanbeta / (1.0D0 / 0.85D0 + tanbeta / 3.0D0)

    call dp%init(beta_deg=16.0D0, K=0.85D0, hardening_mode=DP_HARDENING_TENSION)
    call stress%init(xx=1.0D0, yy=1.0D0, zz=1.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%apex_slope()
    passed = abs(res - expected) < EPS .and. abs(res - dp%stress_eq(stress)) < EPS
    if (.not. passed) then
        print *, "FAIL: Apex slope test"
        print *, "  Expected:", expected
        print *, "  Got     :", res
    else
        print *, "PASS: Apex slope test"
    end if
end subroutine test_dp_apex_slope

subroutine test_dp_apex_dev_gauge(passed)
    ! Gauge of the deviatoric subdifferential at the apex, max over s of (n : s)/t(s), for
    ! coupled directions n (normals and shears at once). References for K < 1 come from the
    ! maximization over the full deviatoric space with mpmath (40 digits).
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    use muscle_yield_vonmises
    use, intrinsic :: ieee_exceptions
    implicit none
    logical, intent(out) :: passed

    type(DruckerPrager) :: dp
    type(VonMises) :: vm
    type(ten_3D2Osym) :: n
    real(real64) :: res, expected
    logical :: invalid
    real(real64), parameter :: TOL = 1.0D-10

    ! K = 1, SHEAR (f = 1): circular cone, the gauge is (2/3) q(n)
    call dp%init(beta_deg=20.0D0, K=1.0D0, hardening_mode=DP_HARDENING_SHEAR)
    call n%init((/1.81536221119D-3, -1.11120307957D-3, -7.0415913162D-4, &
                  2.82442330595D-3, 1.21196306882D-3, 2.0174822516D-4/))
    expected = (2.0D0 / 3.0D0) * vm%stress_eq(n)
    res = dp%apex_dev_gauge(n)
    passed = abs(res - expected) < TOL * expected
    if (.not. passed) then
        print *, "FAIL: Apex gauge K = 1:", res, expected
        return
    end if

    ! K = 0.85, TENSION: Lode angle 0.6 of n, the maximizer is not at the Lode angle of n
    call dp%init(beta_deg=16.0D0, K=0.85D0, hardening_mode=DP_HARDENING_TENSION)
    call n%init((/2.02901221125D-4, 4.19317108241D-4, -6.22218329366D-4, &
                  6.27452450337D-4, 1.99507951108D-3, -1.06702722181D-4/))
    expected = 3.0469160045314333873D-3
    res = dp%apex_dev_gauge(n)
    passed = abs(res - expected) < TOL * expected
    if (.not. passed) then
        print *, "FAIL: Apex gauge K = 0.85:", res, expected
        return
    end if

    ! K = 0.778 (convexity limit), COMPRESSION: Lode angle 0.003, where an unguarded
    ! Newton step leaves [0, pi/3]
    call dp%init(beta_deg=30.0D0, K=0.778D0, hardening_mode=DP_HARDENING_COMPRESSION)
    call n%init((/1.1019546062D-2, -5.01214440894D-3, -6.00740165306D-3, &
                  8.19660642397D-3, -2.93972584491D-3, -7.04295564963D-3/))
    expected = 1.0683686244491726332D-2
    res = dp%apex_dev_gauge(n)
    passed = abs(res - expected) < TOL * expected
    if (.not. passed) then
        print *, "FAIL: Apex gauge K = 0.778:", res, expected
        return
    end if

    ! Zero direction (pure volumetric trial strain): no Lode angle, the gauge is zero and
    ! no invalid operation (0/0) is raised
    call n%init((/0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    call ieee_set_flag(ieee_invalid, .false.)
    res = dp%apex_dev_gauge(n)
    call ieee_get_flag(ieee_invalid, invalid)
    passed = abs(res) < TOL .and. .not. invalid
    if (.not. passed) then
        print *, "FAIL: Apex gauge of a zero direction:", res
    else
        print *, "PASS: Apex gauge test"
    end if
end subroutine test_dp_apex_dev_gauge
