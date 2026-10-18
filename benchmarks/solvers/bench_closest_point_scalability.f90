! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program bench_closest_point_scalability
    !! Strong-scaling benchmark of the `Closest_point` return-mapping solver with OpenMP.
    !!
    !! Solves the same set of synthetic Gauss-point strain states sequentially and then with
    !! several thread counts, and reports for each thread count the wall time, the speedup
    !! \( S_p = t_1 / t_p \), the efficiency \( S_p / p \) and the Karp-Flatt metric (experimentally
    !! determined serial fraction)
    !! \[ e_p = \frac{1/S_p - 1/p}{1 - 1/p}, \]
    !! which stays constant if the loss of scaling is due to serial work and grows with p if it is
    !! due to parallel overhead (Karp & Flatt, Commun. ACM 33(5), 1990).
    !!
    !! It also reports how many points converged (solver `status`) and the mean number of
    !! iterations, and checks that every parallel run reproduces the sequential stresses
    !! (exit code 1 otherwise).
    !!
    !! Timing: one untimed sequential warm-up pass, then every configuration is measured
    !! `--repeat` times and the minimum wall time is reported (least affected by noise).
    !!
    !! This is not a unit test (it is not registered in CTest). Build it with
    !! `-DMUSCLE_BUILD_BENCHMARKS=ON` (preferably in a Release build) and run, e.g.:
    !!
    !! ```text
    !! bench_closest_point_scalability --points 20000 --threads 1,2,4,8 --yield vm --load tension
    !! ```
    use, intrinsic :: iso_fortran_env, only : real64, output_unit
    use omp_lib
    use muscle_tensors
    use muscle_plastic_history
    use muscle_solver_closest_point
    use muscle_elasticity_linear
    use muscle_hard_swift
    use muscle_yield_base, only : Base_yield_critera
    use muscle_yield_vonmises
    use muscle_yield_druckerprager
    implicit none

    integer, parameter :: MAX_CONFIGS = 64
        !! Maximum number of thread counts that can be requested.

    ! Options (command line)
    integer :: n_points = 200000
        !! Number of Gauss points.
    integer :: n_repeat = 3
        !! Repetitions per configuration (the minimum time is reported).
    integer :: thread_list(MAX_CONFIGS)
        !! Thread counts to run.
    integer :: n_configs = 0
        !! Number of entries in `thread_list`.
    character(len=16) :: yield_name = 'vm'
        !! Yield criterion: `vm` (von Mises) or `dp` (Drucker-Prager).
    character(len=16) :: load_name = 'tension'
        !! Load case: `tension` or `compression` (sign of the strain states).

    type(ten_3D2Osym), allocatable :: strains(:)
    type(Plastic_material_history), allocatable :: hist_seq(:), hist_par(:)
    integer, allocatable :: status_seq(:), iters_seq(:)

    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: hard
    type(VonMises) :: vm
    type(DruckerPrager) :: dp
    class(Base_yield_critera), allocatable :: yld
    type(Closest_point) :: solver, solver_loc

    integer :: i, k, r, nth, max_threads, status_loc, n_mismatch
    real(real64) :: t0, t_seq, t_par, t_run, speedup, efficiency, karp_flatt, sgn
    real(real64) :: val_xx, val_yy, val_xy
    logical :: all_ok

    max_threads = omp_get_max_threads()
    call parse_arguments()

    ! -------------------------------------------------------------------------
    ! 1. Material, solver and synthetic strain states
    ! -------------------------------------------------------------------------
    call elas%set_parameters(young=200000.0D0, poisson=0.3D0)
    hard = Swift_hardening(k=100.0D0, n=0.1D0, e0=1.0D0)
    select case (trim(yield_name))
    case ('vm')
        allocate(yld, source=vm)
    case ('dp')
        call dp%init(beta_deg=16D0, K=0.85D0, hardening_mode=DP_HARDENING_TENSION)
        allocate(yld, source=dp)
    end select
    call solver%init(elasticity=elas, hardening=hard, yield=yld)

    sgn = merge(1.0D0, -1.0D0, trim(load_name) == 'tension')
    allocate(strains(n_points), hist_seq(n_points), hist_par(n_points))
    allocate(status_seq(n_points), iters_seq(n_points))
    do i = 1, n_points
        val_xx = 0.02D0 + 0.01D0 * sin(dble(i) * 0.01D0)
        val_yy = -0.001D0 * cos(dble(i) * 0.01D0)
        val_xy = 0.002D0 * sin(dble(i) * 0.05D0)
        call strains(i)%init(xx=sgn*val_xx, yy=sgn*val_yy, zz=0.0D0, xy=sgn*val_xy, yz=0.0D0, xz=0.0D0)
        call hist_seq(i)%init()
    end do

    call print_header()

    ! -------------------------------------------------------------------------
    ! 2. Sequential baseline
    ! -------------------------------------------------------------------------
    call omp_set_num_threads(1)
    ! Untimed warm-up pass (first touch of the arrays, caches)
    do i = 1, n_points
        call solver%solve(strain=strains(i), history=hist_seq(i), status=status_seq(i), iters=iters_seq(i))
    end do
    t_seq = huge(1.0D0)
    do r = 1, n_repeat
        do i = 1, n_points
            call hist_seq(i)%init()
        end do
        t0 = omp_get_wtime()
        do i = 1, n_points
            call solver%solve(strain=strains(i), history=hist_seq(i), status=status_seq(i), iters=iters_seq(i))
        end do
        t_seq = min(t_seq, omp_get_wtime() - t0)
    end do
    call print_status_summary()
    write(*, '(a, f10.4, a, f10.2, a)') " Sequential baseline: ", t_seq, " s  (", &
        1.0D6 * t_seq / n_points, " us/point)"
    print '(a)', ''

    ! -------------------------------------------------------------------------
    ! 3. Strong scaling
    ! -------------------------------------------------------------------------
    print '(a)', " Threads |   Time (s) |  Speedup | Efficiency | Karp-Flatt e | Parallel == sequential"
    print '(a)', " --------+------------+----------+------------+--------------+-----------------------"
    all_ok = .true.
    do k = 1, n_configs
        nth = thread_list(k)
        if (nth > max_threads) then
            write(*, '(i8, a, i0, a)') nth, " | skipped (only ", max_threads, " threads available)"
            cycle
        end if

        call omp_set_num_threads(nth)
        t_par = huge(1.0D0)
        do r = 1, n_repeat
            do i = 1, n_points
                call hist_par(i)%init()
            end do
            t0 = omp_get_wtime()
!$OMP PARALLEL PRIVATE(i, status_loc, solver_loc) SHARED(strains, hist_par, elas, hard, yld)
            ! Each thread builds its own solver (allocatable polymorphic components)
            call solver_loc%init(elasticity=elas, hardening=hard, yield=yld)
!$OMP DO SCHEDULE(guided, 64)
            do i = 1, n_points
                call solver_loc%solve(strain=strains(i), history=hist_par(i), status=status_loc)
            end do
!$OMP END DO
!$OMP END PARALLEL
            t_run = omp_get_wtime() - t0
            t_par = min(t_par, t_run)
        end do

        n_mismatch = 0
        do i = 1, n_points
            if (.not. (hist_par(i)%state_np1%stress .approx. hist_seq(i)%state_np1%stress)) then
                n_mismatch = n_mismatch + 1
            end if
        end do
        all_ok = all_ok .and. (n_mismatch == 0)

        speedup = t_seq / t_par
        efficiency = 100.0D0 * speedup / dble(nth)
        if (nth > 1) then
            karp_flatt = (1.0D0/speedup - 1.0D0/dble(nth)) / (1.0D0 - 1.0D0/dble(nth))
            write(*, '(i8, " |", f11.4, " |", f8.2, "x |", f10.1, "% |", f13.3, " | ", a)') &
                nth, t_par, speedup, efficiency, karp_flatt, mismatch_text(n_mismatch)
        else
            write(*, '(i8, " |", f11.4, " |", f8.2, "x |", f10.1, "% |", a13, " | ", a)') &
                nth, t_par, speedup, efficiency, "-", mismatch_text(n_mismatch)
        end if
        flush(output_unit)
    end do
    print '(a)', ''

    if (.not. all_ok) then
        print '(a)', " ERROR: some parallel runs differ from the sequential baseline."
        stop 1
    end if
    print '(a)', " All parallel runs reproduce the sequential stresses."

contains

    subroutine parse_arguments()
        !! Reads the command-line options (see `print_usage`).
        character(len=256) :: arg, val
        integer :: n_args, j, ios

        n_args = command_argument_count()
        j = 1
        do while (j <= n_args)
            call get_command_argument(j, arg)
            select case (trim(arg))
            case ('-h', '--help')
                call print_usage()
                stop 0
            case ('--points', '--threads', '--yield', '--load', '--repeat')
                if (j == n_args) call usage_error("missing value for " // trim(arg))
                call get_command_argument(j + 1, val)
                select case (trim(arg))
                case ('--points')
                    read(val, *, iostat=ios) n_points
                    if (ios /= 0 .or. n_points < 1) call usage_error("invalid --points: " // trim(val))
                case ('--repeat')
                    read(val, *, iostat=ios) n_repeat
                    if (ios /= 0 .or. n_repeat < 1) call usage_error("invalid --repeat: " // trim(val))
                case ('--threads')
                    call parse_thread_list(val)
                case ('--yield')
                    if (trim(val) /= 'vm' .and. trim(val) /= 'dp') &
                        call usage_error("--yield must be vm or dp")
                    yield_name = val(1:len(yield_name))
                case ('--load')
                    if (trim(val) /= 'tension' .and. trim(val) /= 'compression') &
                        call usage_error("--load must be tension or compression")
                    load_name = val(1:len(load_name))
                end select
                j = j + 2
            case default
                call usage_error("unknown option: " // trim(arg))
            end select
        end do

        if (n_configs == 0) call default_thread_list()
    end subroutine parse_arguments

    subroutine parse_thread_list(text)
        !! Parses a comma-separated list of positive thread counts.
        character(len=*), intent(in) :: text
            !! List, e.g. "1,2,4,8".
        integer :: p0, p1, ios, nth_val

        n_configs = 0
        p0 = 1
        do while (p0 <= len_trim(text))
            p1 = index(text(p0:), ',')
            if (p1 == 0) then
                p1 = len_trim(text) + 1
            else
                p1 = p0 + p1 - 1
            end if
            read(text(p0:p1-1), *, iostat=ios) nth_val
            if (ios /= 0 .or. nth_val < 1) call usage_error("invalid --threads: " // trim(text))
            if (n_configs == MAX_CONFIGS) call usage_error("too many thread counts")
            n_configs = n_configs + 1
            thread_list(n_configs) = nth_val
            p0 = p1 + 1
        end do
    end subroutine parse_thread_list

    subroutine default_thread_list()
        !! 1, 2, 4, 8, ... up to the number of available threads, plus that number itself.
        integer :: nth_val

        nth_val = 1
        do while (nth_val <= max_threads)
            n_configs = n_configs + 1
            thread_list(n_configs) = nth_val
            nth_val = 2 * nth_val
        end do
        if (thread_list(n_configs) /= max_threads) then
            n_configs = n_configs + 1
            thread_list(n_configs) = max_threads
        end if
    end subroutine default_thread_list

    subroutine print_usage()
        !! Prints the command-line help.
        print '(a)', "Usage: bench_closest_point_scalability [options]"
        print '(a)', "  --points N          number of Gauss points (default 200000)"
        print '(a)', "  --repeat R          repetitions per configuration, minimum time reported (default 3)"
        print '(a)', "  --threads a,b,...   thread counts (default 1,2,4,... up to the available threads)"
        print '(a)', "  --yield vm|dp       von Mises (default) or Drucker-Prager"
        print '(a)', "  --load tension|compression   sign of the strain states (default tension)"
        print '(a)', "  -h, --help          show this help"
    end subroutine print_usage

    subroutine usage_error(msg)
        !! Prints an error and the help, and stops with exit code 2.
        character(len=*), intent(in) :: msg
            !! Error message.
        print '(a)', "ERROR: " // msg
        call print_usage()
        stop 2
    end subroutine usage_error

    subroutine print_header()
        !! Prints the configuration of the run.
        character(len=:), allocatable :: build_type

        build_type = MUSCLE_BUILD_TYPE
        print '(a)', " ========================================================================="
        print '(a)', "  MUSCLE benchmark: strong scaling of the closest-point solver (OpenMP)"
        print '(a)', " ========================================================================="
        write(*, '(a, i0)') "  Gauss points      : ", n_points
        write(*, '(a, i0, a)') "  Repetitions       : ", n_repeat, " (minimum time reported, after a warm-up pass)"
        print '(a)', "  Yield criterion   : " // merge('von Mises     ', 'Drucker-Prager', trim(yield_name) == 'vm')
        print '(a)', "  Load case         : " // trim(load_name)
        write(*, '(a, i0)') "  Available threads : ", max_threads
        print '(a)', "  Build type        : " // trim(build_type)
        if (trim(build_type) /= 'Release') then
            print '(a)', "  WARNING: not a Release build; timings are not representative."
        end if
        print '(a)', ''
    end subroutine print_header

    subroutine print_status_summary()
        !! Prints how many points were elastic, converged or not converged in the baseline.
        integer :: n_el, n_conv, n_nonconv

        n_el = count(status_seq == STATUS_ELASTIC_CASE)
        n_conv = count(status_seq == STATUS_CONVERGED)
        n_nonconv = n_points - n_el - n_conv
        write(*, '(a, i0, a, i0, a, i0)') " Solver status (sequential): elastic = ", n_el, &
            ", converged = ", n_conv, ", NOT converged = ", n_nonconv
        if (n_conv + n_nonconv > 0) then
            write(*, '(a, f8.2)') " Mean Newton iterations (plastic points): ", &
                dble(sum(iters_seq, mask=status_seq /= STATUS_ELASTIC_CASE)) / dble(n_conv + n_nonconv)
        end if
        if (n_nonconv > 0) then
            print '(a)', " WARNING: some points did not converge; timings include the maximum iterations."
        end if
    end subroutine print_status_summary

    function mismatch_text(n) result(txt)
        !! "yes" or the number of mismatching points.
        integer, intent(in) :: n
            !! Number of points that differ from the sequential baseline.
        character(len=:), allocatable :: txt
        character(len=32) :: buf

        if (n == 0) then
            txt = "yes"
        else
            write(buf, '("NO (", i0, " points differ)")') n
            txt = trim(buf)
        end if
    end function mismatch_text

end program bench_closest_point_scalability
