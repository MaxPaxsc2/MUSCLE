! SPDX-License-Identifier: GPL-3.0-or-later
! Reproducer: all three isotropic limits must recover the derivatives of J2.
program test_isotropic_derivatives
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use muscle_tensors
    use muscle_yield_base, only: Base_yield_critera
    use muscle_yield_vonmises, only: VonMises
    use muscle_yield_hill48, only: Hill48
    use muscle_yield_cpb06, only: CPB06
    use muscle_yield_yld2004, only: Yld2004
    implicit none
    type(Hill48) :: hill
    type(CPB06) :: cpb
    type(Yld2004) :: yld
    type(VonMises) :: vm
    real(real64) :: eye(6,6)
    logical :: passed
    integer :: i
    passed=.true.
    call hill%init(.5_real64,.5_real64,.5_real64,1.5_real64,1.5_real64,1.5_real64)
    call cpb%init(1._real64,0._real64,0._real64,0._real64,1._real64,0._real64, &
                  0._real64,0._real64,1._real64,1._real64,1._real64,1._real64,0._real64,2._real64)
    eye=0._real64
    do i=1,6
        eye(i,i)=1._real64
    end do
    call yld%init(eye,eye,2._real64)
    call check(hill,'Hill48',passed)
    call check(cpb,'CPB06',passed)
    call check(yld,'Yld2004',passed)
    if (.not. passed) stop 1
    stop 0
contains
    subroutine check(criterion,name,all_passed)
        class(Base_yield_critera), intent(in) :: criterion
        character(len=*), intent(in) :: name
        logical, intent(inout) :: all_passed
        type(ten_3D2Osym) :: s,g,gr
        type(ten_3D4O3sym) :: h,hr
        real(real64) :: dq,dg,dh,x
        integer :: j
        dq=0._real64;dg=0._real64;dh=0._real64
        do j=1,4
            select case(j)
            case(1)
                s%vals=[100._real64,0._real64,0._real64,0._real64,0._real64,0._real64]
            case(2)
                s%vals=[80._real64,80._real64,0._real64,0._real64,0._real64,0._real64]
            case(3)
                s%vals=[50._real64,50._real64,0._real64,50._real64,0._real64,0._real64]
            case(4)
                s%vals=[120._real64,-30._real64,40._real64,25._real64,11._real64,-17._real64]
            end select
            x=criterion%stress_eq(s)
            g=criterion%dstressEq_dstress(s);gr=vm%dstressEq_dstress(s)
            h=criterion%ddstressEq_ddstress(s);hr=vm%ddstressEq_ddstress(s)
            if (.not. ieee_is_finite(x) .or. .not. all(ieee_is_finite(g%vals)) .or. &
                .not. all(ieee_is_finite(h%vals))) all_passed=.false.
            dq=max(dq,abs(x-vm%stress_eq(s)))
            dg=max(dg,maxval(abs(g%vals-gr%vals)))
            dh=max(dh,maxval(abs(h%vals-hr%vals)))
        end do
        print '(A,3ES16.7)',name//' max abs residuals (value, gradient, full Hessian): ',dq,dg,dh
        if (dq>1e-8_real64 .or. dg>1e-10_real64 .or. dh>1e-10_real64) all_passed=.false.
    end subroutine check
end program test_isotropic_derivatives
