program test_muscle_solver_cutting_plane
    implicit none
    logical :: passed

    call test_cutting_plane_vonmises(passed)
    if (.not. passed) STOP 1


    print*, "Passed!", passed
    STOP 0
end program test_muscle_solver_cutting_plane

subroutine test_cutting_plane_vonmises(passed)
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_solver_cutting_plane, only : cutting_plane
    implicit none

    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(VonMises) :: vm
    type(ten_3D2Osym) :: strain, strain_p, stress
    type(Swift_hardening) :: sw
    type(Elasticity_linear) :: elas
    real(real64) :: strain_pf
    logical :: error
    type(ten_3D2Osym) :: expected_stress
    real(real64), parameter :: TOL=1D-8
    ! Radial-return reference for uniaxial strain exx=0.8, E=1000, nu=0.3, Swift(100, 0.1, 1e-4)
    real(real64), parameter :: EXPECTED_STRAIN_PF=0.453258459032D0

    passed = .False.
    strain_pf = 0D0
    error = .False.
    

    call strain%init(xx=0.8D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)

    call cutting_plane(strain=strain, elasticity=elas, hardening=sw, yield=vm,    &
                       stress=stress, strain_pf=strain_pf, strain_p=strain_p, error=error)

    write(*,*) "stress"
    write(*,*) stress
    write(*,*) 

    write(*,*) "strain_pf:", strain_pf
    write(*,*) 

    write(*,*) "strain_p"
    write(*,*) strain_p
    write(*,*) 

    write(*,*) "Error:", error

    passed = .not. error
    if (.not. passed) return

    passed = abs(strain_pf - EXPECTED_STRAIN_PF) < TOL
    if (.not. passed) print*, "Effective plastic strain is not equal", strain_pf - EXPECTED_STRAIN_PF
    if (.not. passed) return

    call expected_stress%init(xx=728.2627238217D0, yy=635.8686380891D0, zz=635.8686380891D0, &
                              xy=0D0, yz=0D0, xz=0D0)
    passed = stress%is_approx(expected_stress, tol=TOL)
    if (.not. passed) print*, "Stress is not equal", stress - expected_stress
    if (.not. passed) return
end subroutine
