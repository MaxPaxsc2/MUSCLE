! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_yield_cpb06ex2
    implicit none
    logical :: passed
    call test_CPB06ex2_reduces_to_CPB06(passed)
    if (.not. passed) STOP 1
    call test_CPB06ex2_isotropic_vonmises(passed)
    if (.not. passed) STOP 2
    call test_CPB06ex2_paper_tables(passed)
    if (.not. passed) STOP 3
    call test_CPB06ex2_tension_compression_asymmetry(passed)
    if (.not. passed) STOP 4
    call test_CPB06ex2_pressure_insensitive(passed)
    if (.not. passed) STOP 5
    call test_CPB06ex2_gradient_vs_numeric(passed)
    if (.not. passed) STOP 6
    call test_CPB06ex2_gradient_analytical(passed)
    if (.not. passed) STOP 7
    print*, "Passed!", passed
    STOP 0
end program test_muscle_yield_cpb06ex2

subroutine test_CPB06ex2_reduces_to_CPB06(passed)
    ! With Cp = C and kp = k (the primed map), CPB06ex2 is CPB06: both sums are equal and
    ! normalized the same way.
    ! General nonsymmetric normal block and three distinct shear factors expose map errors.
    ! Called by program test_muscle_yield_cpb06ex2 (STOP 1)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06ex2
    use muscle_yield_cpb06
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    ! Same eigenvalues and powers in both models: round-off only. Measured: 1.6e-16
    real(real64), parameter :: TOL = 1D-12
    type(CPB06ex2) :: yld
    type(CPB06) :: cpb
    type(ten_3D2Osym) :: s
    passed = .true.
    call cpb%init(1D0,-0.2D0,0.3D0,0.15D0,0.8D0,-0.4D0,-0.1D0,0.2D0,1.2D0, &
                  1.3D0,0.7D0,1.1D0,0.25D0,4D0)
    call yld%init(1D0,-0.2D0,0.3D0,0.15D0,0.8D0,-0.4D0,-0.1D0,0.2D0,1.2D0, &
                  1.3D0,0.7D0,1.1D0,0.25D0, &
                  1D0,-0.2D0,0.3D0,0.15D0,0.8D0,-0.4D0,-0.1D0,0.2D0,1.2D0, &
                  1.3D0,0.7D0,1.1D0,0.25D0,4D0)
    call s%init((/120.5D0,-35.25D0,60.75D0,42.5D0,-28D0,33.25D0/))
    call check_approx(passed, 'CPB06 reduction coupled', yld%stress_eq(s), cpb%stress_eq(s), &
                      tol=TOL)
    s = -s
    call check_approx(passed, 'CPB06 reduction compression', yld%stress_eq(s), cpb%stress_eq(s), &
                      tol=TOL)
end subroutine test_CPB06ex2_reduces_to_CPB06

subroutine test_CPB06ex2_isotropic_vonmises(passed)
    ! Isotropic a=2 checks all tensorial shear weights against the closed-form VonMises.
    ! Called by program test_muscle_yield_cpb06ex2 (STOP 2)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06ex2
    use muscle_yield_vonmises
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    ! Closed form against eigenvalues: round-off only. Measured: 1.7e-16
    real(real64), parameter :: TOL = 1D-12
    type(CPB06ex2) :: yld
    type(VonMises) :: vm
    type(ten_3D2Osym) :: s
    passed = .true.
    call yld%init(1D0,0D0,0D0,0D0,1D0,0D0,0D0,0D0,1D0,1D0,1D0,1D0,0D0, &
                  1D0,0D0,0D0,0D0,1D0,0D0,0D0,0D0,1D0,1D0,1D0,1D0,0D0,2D0)
    call s%init((/120.5D0,-35.25D0,60.75D0,42.5D0,-28D0,33.25D0/))
    call check_approx(passed, 'VonMises coupled', yld%stress_eq(s), vm%stress_eq(s), tol=TOL)
end subroutine test_CPB06ex2_isotropic_vonmises

subroutine CPB06ex2_test_materials(yld)
    ! Plunkett et al. (2008), Tables 2, 8, 1: steel, 6111-T4 and AZ31B, respectively.
    ! The shear order yz,xz,xy of the paper is permuted to library xy,yz,xz.
    ! AZ31B lacks out-of-plane shear data: unit factors there are synthetic test inputs.
    ! Called by test_CPB06ex2_paper_tables, test_CPB06ex2_tension_compression_asymmetry,
    ! test_CPB06ex2_pressure_insensitive and the gradient tests
    use muscle_yield_cpb06ex2
    implicit none
    type(CPB06ex2), intent(out) :: yld(3)
    call yld(1)%init(1.283D0,0.208D0,-0.421D0,0.208D0,0.535D0,-0.791D0,-0.421D0,-0.791D0,0.746D0, &
                     1.156D0,1.418D0,1.329D0,0.030D0, &
                     1.283D0,0.380D0,-0.001D0,0.380D0,2.145D0,0.792D0,-0.001D0,0.792D0,0.218D0, &
                     1.721D0,1.418D0,1.329D0,-0.027D0,5D0)
    call yld(2)%init(0.250D0,-1.399D0,-0.971D0,-1.399D0,-0.375D0,-2.254D0, &
                     -0.971D0,-2.254D0,-0.831D0, &
                     1.316D0,1.403D0,1.101D0,0D0, &
                     0.250D0,-0.840D0,-1.330D0,-0.840D0,0.780D0,-0.565D0, &
                     -1.330D0,-0.565D0,0.479D0, &
                     1.413D0,1.094D0,1.448D0,0D0,12D0)
    call yld(3)%init(0.771D0,-0.914D0,-1.312D0,-0.914D0,0.654D0,-1.209D0, &
                     -1.312D0,-1.209D0,-0.293D0, &
                     1.617D0,1D0,1D0,0.295D0, &
                     0.771D0,0.337D0,-0.764D0,0.337D0,0.646D0,-0.401D0,-0.764D0,-0.401D0,-0.022D0, &
                     0.499D0,1D0,1D0,0.900D0,4D0)
end subroutine CPB06ex2_test_materials

subroutine test_CPB06ex2_paper_tables(passed)
    ! Rounded steel parameters and out-of-plane 6111-T4 loadings detect the shear permutation.
    ! Called by program test_muscle_yield_cpb06ex2 (STOP 3)
    use muscle_tensors
    use muscle_yield_cpb06ex2
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(CPB06ex2) :: yld(3)
    type(ten_3D2Osym) :: s
    passed = .true.
    call CPB06ex2_test_materials(yld)
    ! Table 2 has three-digit parameters: 2e-3 allows their rounding, not solver error.
    ! Measured: 9.7e-4 or less (Y_C)
    call s%init((/1D0,0D0,0D0,0D0,0D0,0D0/))
    call check_approx(passed, 'steel Table 3 X_T', 1D0/yld(1)%stress_eq(s), 1D0, tol=2D-3)
    call s%init((/0D0,1D0,0D0,0D0,0D0,0D0/))
    call check_approx(passed, 'steel Table 3 Y_T', 1D0/yld(1)%stress_eq(s), 0.914D0, tol=2D-3)
    call s%init((/-1D0,0D0,0D0,0D0,0D0,0D0/))
    call check_approx(passed, 'steel Table 3 X_C', 1D0/yld(1)%stress_eq(s), 0.963D0, tol=2D-3)
    call s%init((/0D0,-1D0,0D0,0D0,0D0,0D0/))
    call check_approx(passed, 'steel Table 3 Y_C', 1D0/yld(1)%stress_eq(s), 0.960D0, tol=2D-3)
    call s%init((/0D0,0D0,-1D0,0D0,0D0,0D0/))
    call check_approx(passed, 'steel Table 3 Z_C', 1D0/yld(1)%stress_eq(s), 0.958D0, tol=2D-3)
    ! Table 9 has two decimals: half a last digit plus rounding of Table 8. Measured: 2.6e-4 or less
    call s%init((/0D0,0.5D0,0.5D0,0D0,0.5D0,0D0/))
    call check_approx(passed, '6111 Table 9 YZ 45', 1D0/yld(2)%stress_eq(s), 1.07D0, tol=6D-3)
    call s%init((/0.5D0,0D0,0.5D0,0D0,0D0,0.5D0/))
    call check_approx(passed, '6111 Table 9 ZX 45', 1D0/yld(2)%stress_eq(s), 1.04D0, tol=6D-3)
    call s%init((/0D0,0D0,0D0,0D0,1D0,0D0/))
    call check_approx(passed, '6111 Table 9 YZ shear', 1D0/yld(2)%stress_eq(s), 0.67D0, tol=6D-3)
    call s%init((/0D0,0D0,0D0,0D0,0D0,1D0/))
    call check_approx(passed, '6111 Table 9 ZX shear', 1D0/yld(2)%stress_eq(s), 0.65D0, tol=6D-3)
    ! Independent 90-digit trace-of-powers reference, all normal and shear components coupled.
    ! Measured: 1.7e-16
    call s%init((/120.5D0,-35.25D0,60.75D0,42.5D0,-28D0,33.25D0/))
    call check_approx(passed, '6111 independent scalar', yld(2)%stress_eq(s), &
                      170.33623760110963846252D0, tol=1D-12)
end subroutine test_CPB06ex2_paper_tables

subroutine test_CPB06ex2_tension_compression_asymmetry(passed)
    ! kp = 0.9 (the primed k) makes AZ31B compression much weaker than tension in RD (Table 1).
    ! The RD normalization of init is exact: f = 1 in unit tension along RD (measured: 0).
    ! Called by program test_muscle_yield_cpb06ex2 (STOP 4)
    use muscle_tensors
    use muscle_yield_cpb06ex2
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(CPB06ex2) :: yld(3)
    type(ten_3D2Osym) :: s
    passed = .true.
    call CPB06ex2_test_materials(yld)
    call s%init((/1D0,0D0,0D0,0D0,0D0,0D0/))
    call check_approx(passed, 'AZ31 RD normalization', yld(3)%stress_eq(s), 1D0, tol=1D-12)
    call check_true(passed, 'AZ31 tension compression asymmetry', &
                    yld(3)%stress_eq(-s) > 1.5D0*yld(3)%stress_eq(s))
end subroutine test_CPB06ex2_tension_compression_asymmetry

subroutine test_CPB06ex2_pressure_insensitive(passed)
    ! Both maps act on the deviator (Eq. 9): a hydrostatic stress gives zero, and adding a
    ! pressure to a coupled state leaves the equivalent stress unchanged. AZ31B, kp = 0.9.
    ! Called by program test_muscle_yield_cpb06ex2 (STOP 5)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06ex2
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    ! linear_transform subtracts the zz stress before the products, so the hydrostatic state
    ! maps to an exact zero and the shifted state differs only by the rounding of the
    ! differences. Measured: exact (0) in both checks
    real(real64), parameter :: TOL = 1D-12
    type(CPB06ex2) :: yld(3)
    type(ten_3D2Osym) :: s, p
    passed = .true.
    call CPB06ex2_test_materials(yld)
    call p%init((/250D0,250D0,250D0,0D0,0D0,0D0/))
    call check_approx(passed, 'AZ31 hydrostatic stress = 0', yld(3)%stress_eq(p), 0D0, tol=TOL)
    call s%init((/120.5D0,-35.25D0,60.75D0,42.5D0,-28D0,33.25D0/))
    call check_approx(passed, 'AZ31 coupled state plus pressure 250', yld(3)%stress_eq(s + p), &
                      yld(3)%stress_eq(s), tol=TOL)
end subroutine test_CPB06ex2_pressure_insensitive

subroutine test_CPB06ex2_gradient_vs_numeric(passed)
    ! Gradient against the finite-difference gradient of the base class, for the three
    ! materials in a state with normal and shear components together, and for 6111-T4 in a
    ! rotated state where two principal values of Sigma differ by about 1e-8 relatively
    ! Called by program test_muscle_yield_cpb06ex2 (STOP 6)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06ex2
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    ! Central differences of the base class: truncation and the rounding of eigenvals.
    ! Measured: 1.2e-10 or less
    real(real64), parameter :: TOL = 1D-8
    type(CPB06ex2) :: yld(3)
    type(ten_3D2Osym) :: s
    passed = .true.
    call CPB06ex2_test_materials(yld)
    call s%init((/120.5D0,-35.25D0,60.75D0,42.5D0,-28D0,33.25D0/))
    ! Steel: k and kp of opposite sign, a = 5
    call check_approx(passed, 'gradient vs numeric: steel', yld(1)%dstressEq_dstress(s), &
                      yld(1)%dstressEq_dstress_numeric(s), tol=TOL)
    ! 6111-T4: a = 12, the largest powers of the three
    call check_approx(passed, 'gradient vs numeric: 6111-T4', yld(2)%dstressEq_dstress(s), &
                      yld(2)%dstressEq_dstress_numeric(s), tol=TOL)
    ! AZ31B: kp = 0.9, the second transformation almost one-sided
    call check_approx(passed, 'gradient vs numeric: AZ31B', yld(3)%dstressEq_dstress(s), &
                      yld(3)%dstressEq_dstress_numeric(s), tol=TOL)
    call s%init((/197.95561360136386D0,-58.04950942270237D0,-139.90610417866154D0, &
                  -50.24571382486021D0,11.771438358965758D0,-51.20999705196427D0/))
    call check_approx(passed, 'gradient vs numeric: 6111-T4, close pair', &
                      yld(2)%dstressEq_dstress(s), yld(2)%dstressEq_dstress_numeric(s), tol=TOL)
end subroutine test_CPB06ex2_gradient_vs_numeric

subroutine test_CPB06ex2_gradient_analytical(passed)
    ! Gradient of 6111-T4 against Euler's identity for the homogeneous f (grad : stress = f)
    ! and against an independent 90-digit reference computed with traces of powers of
    ! Sigma and Sigmap (no eigenvalues), in a state with all components
    ! Called by program test_muscle_yield_cpb06ex2 (STOP 7)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06ex2
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    ! Exact references: round-off only. Measured: 2.8e-16 or less
    real(real64), parameter :: TOL = 1D-12
    type(CPB06ex2) :: yld(3)
    type(ten_3D2Osym) :: s, expected
    passed = .true.
    call CPB06ex2_test_materials(yld)
    call s%init((/120.5D0,-35.25D0,60.75D0,42.5D0,-28D0,33.25D0/))
    call check_approx(passed, 'gradient Euler identity: 6111-T4', &
                      yld(2)%dstressEq_dstress(s) .ddot. s, yld(2)%stress_eq(s), tol=TOL)
    call expected%init((/0.63491427307714846058430D0, -0.8014554330786033866660D0, &
                         0.16654116000145492608174D0, 0.29454772703368041095476D0, &
                         -0.2025437658446382337026D0, 0.28693805299310551585208D0/))
    call check_approx(passed, 'gradient vs reference: 6111-T4', yld(2)%dstressEq_dstress(s), &
                      expected, tol=TOL)
end subroutine test_CPB06ex2_gradient_analytical
