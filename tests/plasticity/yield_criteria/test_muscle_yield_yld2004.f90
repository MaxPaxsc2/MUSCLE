! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_yield_yld2004
    implicit none
    logical :: passed

    call test_Yld2004_stresseq_1(passed)
    if (.not. passed) STOP 1

    print*, "Passed!", passed
    STOP 0
end program test_muscle_yield_yld2004

subroutine Yld2004_test_materials(yld)
    ! yld(1): isotropic (all 18 coefficients 1, so S' = S'' = s), a = 2: von Mises.
    ! yld(2): AA2090-T3, a = 8 (Barlat et al., 2005, Table 2).
    use muscle_yield_yld2004
    implicit none
    type(Yld2004), intent(out) :: yld(2)

    call yld(1)%init(cp12=1D0, cp13=1D0, cp21=1D0, cp23=1D0, cp31=1D0, cp32=1D0, &
                     cp44=1D0, cp55=1D0, cp66=1D0, &
                     cpp12=1D0, cpp13=1D0, cpp21=1D0, cpp23=1D0, cpp31=1D0, cpp32=1D0, &
                     cpp44=1D0, cpp55=1D0, cpp66=1D0, a=2D0)
    call yld(2)%init(cp12=-0.069888D0, cp13=0.936408D0, cp21=0.079143D0, cp23=1.003060D0, &
                     cp31=0.524741D0, cp32=1.363180D0, cp44=1.023770D0, cp55=1.069060D0, cp66=0.954322D0, &
                     cpp12=0.981171D0, cpp13=0.476741D0, cpp21=0.575316D0, cpp23=0.866827D0, &
                     cpp31=1.145010D0, cpp32=-0.079294D0, cpp44=1.051660D0, cpp55=1.147100D0, &
                     cpp66=1.404620D0, a=8D0)
end subroutine Yld2004_test_materials

subroutine test_Yld2004_stresseq_1(passed)
    ! Equivalent stress against independent references. The AA2090-T3 values come from a
    ! Python script with mpmath at 90 digits that writes phi with traces of powers of S' and
    ! S'' (valid for even a), without eigenvalues.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_yld2004
    use muscle_yield_vonmises
    implicit none
    logical, intent(out) :: passed

    ! Measured errors are below 4e-16; 1e-12 leaves room for other compilers and still separates
    ! the shear assignments, which differ by more than 1e-3
    real(real64), parameter :: TOL = 1.0D-12
    type(Yld2004) :: yld(2)
    type(VonMises) :: vm
    type(ten_3D2Osym) :: s, s_aa(4)
    real(real64) :: expected(4), err
    integer :: i

    call Yld2004_test_materials(yld)

    ! Isotropic limit, general 3D state: 4 f**2 = sum_ij (s_i - s_j)**2, i.e. von Mises
    ! (three distinct eigenvalues with all shears present)
    call s%init((/180D0, -40D0, 25D0, 60D0, -30D0, 45D0/))
    err = abs(yld(1)%stress_eq(s) - vm%stress_eq(s))/vm%stress_eq(s)
    passed = err <= TOL
    if (.not. passed) print*, "Yld2004 isotropic vs von Mises, rel. error =", err
    if (.not. passed) return

    ! Uniaxial RD: S' and S'' diagonal, pairs lam'_i - lam''_j of both signs
    call s_aa(1)%init((/100D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    ! Uniaxial TD: the same with the other rows of C' and C''
    call s_aa(2)%init((/0D0, 100D0, 0D0, 0D0, 0D0, 0D0/))
    ! Uniaxial at 45 deg: the xy shear picks c66 (Barlat shear order yz, zx, xy)
    call s_aa(3)%init((/50D0, 50D0, 0D0, 50D0, 0D0, 0D0/))
    ! Rotated 3D state: yz and xz shears pick c44 and c55
    call s_aa(4)%init((/120.5D0, -35.25D0, 60.75D0, 42.5D0, -28.0D0, 33.25D0/))
    expected = (/99.932153929799531134D0, 110.32696825195229388D0, &
                 122.30436942533559492D0, 187.28212180970985109D0/)
    do i = 1, 4
        err = abs(yld(2)%stress_eq(s_aa(i)) - expected(i))/expected(i)
        passed = err <= TOL
        if (.not. passed) print*, "Yld2004 AA2090-T3 state", i, "rel. error =", err
        if (.not. passed) return
    end do

    ! Homogeneity of degree 1 and pressure insensitivity on the rotated state
    s = 2.5D0*s_aa(4)
    err = abs(yld(2)%stress_eq(s) - 2.5D0*expected(4))/expected(4)
    passed = err <= TOL
    if (.not. passed) print*, "Yld2004 homogeneity, rel. error =", err
    if (.not. passed) return
    call s%init(s_aa(4)%vals + (/300D0, 300D0, 300D0, 0D0, 0D0, 0D0/))
    err = abs(yld(2)%stress_eq(s) - expected(4))/expected(4)
    passed = err <= TOL
    if (.not. passed) print*, "Yld2004 pressure insensitivity, rel. error =", err
end subroutine test_Yld2004_stresseq_1
