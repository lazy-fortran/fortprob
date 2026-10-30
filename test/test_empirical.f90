module empirical_tests
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
        ieee_positive_inf, ieee_negative_inf, ieee_is_nan, ieee_is_finite
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use fortnum_status, only: fortnum_status_t, status_ok
    use fortprob_empirical, only: empirical_calibration_t
    implicit none
contains
    subroutine collect(tests)
        type(unittest_type), allocatable, intent(out) :: tests(:)
        allocate(tests(4))
        tests(1) = new_unittest('ties-and-boundaries', ties)
        tests(2) = new_unittest('quantile-ranks', ranks)
        tests(3) = new_unittest('invalid-input', invalid)
        tests(4) = new_unittest('permutation-and-monotonicity', permutation)
    end subroutine

    subroutine ties(error)
        type(error_type), allocatable, intent(out) :: error
        type(empirical_calibration_t) :: cal
        type(fortnum_status_t) :: status
        real(dp) :: p, queries(7), expected(7)
        integer :: i
        call cal%fit([4.0_dp, 2.0_dp, 2.0_dp, -1.0_dp], status)
        call check(error, status_ok(status)); if (allocated(error)) return
        call check(error, cal%count() == 4); if (allocated(error)) return
        queries = [-2.0_dp, -1.0_dp, 0.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
        expected = [1.0_dp, 1.0_dp, 0.8_dp, 0.8_dp, 0.4_dp, 0.4_dp, 0.2_dp]
        do i = 1, size(queries)
            p = cal%tail_probability(queries(i), status)
            call check(error, status_ok(status)); if (allocated(error)) return
            call check(error, abs(p - expected(i)) < 1e-15_dp)
            if (allocated(error)) return
        end do
    end subroutine

    subroutine ranks(error)
        type(error_type), allocatable, intent(out) :: error
        type(empirical_calibration_t) :: cal
        type(fortnum_status_t) :: status
        real(dp) :: q
        call cal%fit([8.0_dp, 2.0_dp, 4.0_dp], status)
        q = cal%upper_quantile(0.5_dp, status)
        call check(error, status_ok(status)); if (allocated(error)) return
        call check(error, q == 4.0_dp); if (allocated(error)) return
        q = cal%upper_quantile(0.25_dp, status)
        call check(error, q == 8.0_dp); if (allocated(error)) return
        q = cal%upper_quantile(0.125_dp, status)
        call check(error, .not. ieee_is_finite(q) .and. q > 0.0_dp)
        if (allocated(error)) return
        call cal%fit([7.0_dp], status)
        q = cal%upper_quantile(0.05_dp, status)
        call check(error, status_ok(status)); if (allocated(error)) return
        call check(error, .not. ieee_is_finite(q) .and. q > 0.0_dp)
        if (allocated(error)) return
        q = cal%upper_quantile(0.5_dp, status)
        call check(error, q == 7.0_dp)
    end subroutine

    subroutine invalid(error)
        type(error_type), allocatable, intent(out) :: error
        type(empirical_calibration_t) :: cal
        type(fortnum_status_t) :: status
        real(dp) :: bad(3), alphas(5), q, pair(2), empty(0)
        integer :: i
        bad(1) = ieee_value(0.0_dp, ieee_quiet_nan)
        bad(2) = ieee_value(0.0_dp, ieee_positive_inf)
        bad(3) = ieee_value(0.0_dp, ieee_negative_inf)
        q = cal%tail_probability(1.0_dp, status)
        call check(error, .not. status_ok(status) .and. ieee_is_nan(q))
        if (allocated(error)) return
        q = cal%upper_quantile(0.5_dp, status)
        call check(error, .not. status_ok(status) .and. ieee_is_nan(q))
        if (allocated(error)) return
        call cal%fit(empty, status)
        call check(error, .not. status_ok(status)); if (allocated(error)) return
        do i = 1, size(bad)
            call cal%fit([1.0_dp], status)
            pair(1) = 1.0_dp
            pair(2) = bad(i)
            call cal%fit(pair, status)
            call check(error, .not. status_ok(status) .and. cal%count() == 0)
            if (allocated(error)) return
            call cal%fit([1.0_dp], status)
            q = cal%tail_probability(bad(i), status)
            call check(error, .not. status_ok(status) .and. ieee_is_nan(q))
            if (allocated(error)) return
        end do
        alphas(1:4) = [0.0_dp, 1.0_dp, -0.1_dp, 1.1_dp]
        alphas(5) = bad(1)
        do i = 1, size(alphas)
            q = cal%upper_quantile(alphas(i), status)
            call check(error, .not. status_ok(status) .and. ieee_is_nan(q))
            if (allocated(error)) return
        end do
        do i = 2, 3
            q = cal%upper_quantile(bad(i), status)
            call check(error, .not. status_ok(status) .and. ieee_is_nan(q))
            if (allocated(error)) return
        end do
        call cal%clear()
        call cal%clear()
        call check(error, cal%count() == 0)
    end subroutine

    subroutine permutation(error)
        type(error_type), allocatable, intent(out) :: error
        type(empirical_calibration_t) :: first, second
        type(fortnum_status_t) :: status
        real(dp) :: p, previous, other, q, previous_q
        integer :: i
        call first%fit([4.0_dp, -1.0_dp, 2.0_dp, 2.0_dp], status)
        call second%fit([2.0_dp, 2.0_dp, -1.0_dp, 4.0_dp], status)
        previous = 1.0_dp
        previous_q = huge(1.0_dp)
        do i = 1, 19
            p = first%tail_probability(real(i-5, dp)/2.0_dp, status)
            other = second%tail_probability(real(i-5, dp)/2.0_dp, status)
            call check(error, p == other .and. p <= previous)
            if (allocated(error)) return
            previous = p
            q = first%upper_quantile(real(i, dp)/20.0_dp, status)
            other = second%upper_quantile(real(i, dp)/20.0_dp, status)
            call check(error, q == other); if (allocated(error)) return
            if (ieee_is_finite(q)) then
                call check(error, q <= previous_q); if (allocated(error)) return
                previous_q = q
            end if
        end do
    end subroutine
end module

program test_empirical
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: testsuite_type, new_testsuite, run_testsuite
    use empirical_tests, only: collect
    implicit none
    type(testsuite_type) :: suite
    integer :: stat
    stat = 0
    suite = new_testsuite('empirical', collect)
    call run_testsuite(suite%collect, error_unit, stat)
    if (stat /= 0) error stop 1
end program
