module truncated_tests
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use fortnum_status, only: fortnum_status_t, status_ok
    use fortprob_truncated, only: normal_lower_truncated_moments
    implicit none
contains
    subroutine collect(tests)
        type(unittest_type), allocatable, intent(out) :: tests(:)
        allocate (tests(3))
        tests(1) = new_unittest('half-normal-and-affine', half_normal)
        tests(2) = new_unittest('quadrature-tail-oracle', quadrature)
        tests(3) = new_unittest('deterministic-and-invalid', invalid)
    end subroutine
    subroutine half_normal(error)
        type(error_type), allocatable, intent(out) :: error
        type(fortnum_status_t) :: status
        real(dp) :: m, v, l
        call normal_lower_truncated_moments(0.0_dp, 1.0_dp, 0.0_dp, m, v, l, status)
        call check(error, status_ok(status)); if (allocated(error)) return
        call check(error, abs(m - sqrt(2.0_dp/acos(-1.0_dp))) < 1.e-14_dp)
        if (allocated(error)) return
        call check(error, abs(v - (1.0_dp - 2.0_dp/acos(-1.0_dp))) < 1.e-14_dp)
        if (allocated(error)) return
        call check(error, abs(l + log(2.0_dp)) < 1.e-14_dp)
        if (allocated(error)) return
        call normal_lower_truncated_moments(3.0_dp, 4.0_dp, 3.0_dp, m, v, l, status)
        call check(error, abs(m - 3.0_dp - 2.0_dp*sqrt(2.0_dp/acos(-1.0_dp))) < 1.e-14_dp)
        if (allocated(error)) return
        call check(error, abs(v - 4.0_dp*(1.0_dp - 2.0_dp/acos(-1.0_dp))) < 1.e-14_dp)
        if (allocated(error)) return
        call normal_lower_truncated_moments(10.0_dp, 1.0_dp, 0.0_dp, m, v, l, status)
        call check(error, abs(m - 10.0_dp) < 1.e-14_dp .and. abs(v - 1.0_dp) < 1.e-14_dp)
    end subroutine
    subroutine quadrature(error)
        type(error_type), allocatable, intent(out) :: error
        type(fortnum_status_t) :: status
        real(dp), parameter :: cases(6) = [1.0_dp, 5.0_dp, 8.0_dp, 8.01_dp, 40.0_dp, 100.0_dp]
        real(dp) :: a, h, t, w, sums(3), m, v, l, expected_m, expected_v, expected_l
        integer :: j, i
        ! Independent Simpson integration of the excess density on [0,infinity).
        do j = 1, size(cases)
            a = cases(j)
            h = 40.0_dp/(a + 1.0_dp)/20000.0_dp
            sums = 0.0_dp
            do i = 0, 20000
                t = real(i, dp)*h
                w = 2.0_dp
                if (mod(i, 2) == 1) w = 4.0_dp
                if (i == 0 .or. i == 20000) w = 1.0_dp
                w = w*exp(-a*t - 0.5_dp*t*t)
                sums(1) = sums(1) + w
                sums(2) = sums(2) + w*t
                sums(3) = sums(3) + w*t*t
            end do
            sums = sums*h/3.0_dp
            expected_m = sums(2)/sums(1)
            expected_v = sums(3)/sums(1) - expected_m**2
            expected_l = -0.5_dp*a*a - 0.5_dp*log(2.0_dp*acos(-1.0_dp)) + log(sums(1))
            call normal_lower_truncated_moments(-a, 1.0_dp, 0.0_dp, m, v, l, status)
            call check(error, status_ok(status)); if (allocated(error)) return
            call check(error, abs(m/expected_m - 1.0_dp) < 1.e-9_dp)
            if (allocated(error)) return
            call check(error, abs(v/expected_v - 1.0_dp) < 1.e-9_dp)
            if (allocated(error)) return
            call check(error, abs(l - expected_l) < 2.e-11_dp)
            if (allocated(error)) return
            call normal_lower_truncated_moments(7.0_dp - 3.0_dp*a, 9.0_dp, &
                7.0_dp, m, v, l, status)
            call check(error, abs((m - 7.0_dp)/(3.0_dp*expected_m) - 1.0_dp) < 1.e-9_dp)
            if (allocated(error)) return
            call check(error, abs(v/(9.0_dp*expected_v) - 1.0_dp) < 1.e-9_dp)
            if (allocated(error)) return
        end do
    end subroutine
    subroutine invalid(error)
        type(error_type), allocatable, intent(out) :: error
        type(fortnum_status_t) :: status
        real(dp) :: m, v, l
        call normal_lower_truncated_moments(2.0_dp, 0.0_dp, 2.0_dp, m, v, l, status)
        call check(error, status_ok(status) .and. m == 2.0_dp .and. v == 0.0_dp .and. l == 0.0_dp)
        if (allocated(error)) return
        call normal_lower_truncated_moments(1.0_dp, 0.0_dp, 2.0_dp, m, v, l, status)
        call check(error, .not. status_ok(status) .and. ieee_is_nan(m))
        if (allocated(error)) return
        call normal_lower_truncated_moments(0.0_dp, -1.0_dp, 0.0_dp, m, v, l, status)
        call check(error, .not. status_ok(status)); if (allocated(error)) return
        call normal_lower_truncated_moments(ieee_value(0.0_dp, ieee_quiet_nan), &
            1.0_dp, 0.0_dp, m, v, l, status)
        call check(error, .not. status_ok(status))
    end subroutine
end module
program test_truncated
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: testsuite_type, new_testsuite, run_testsuite
    use truncated_tests, only: collect
    implicit none
    type(testsuite_type) :: suite
    integer :: stat
    stat = 0
    suite = new_testsuite('truncated', collect)
    call run_testsuite(suite%collect, error_unit, stat)
    if (stat /= 0) error stop 1
end program
