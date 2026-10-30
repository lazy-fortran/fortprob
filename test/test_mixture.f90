module mixture_tests
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan, &
        ieee_positive_inf
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use fortnum_status, only: fortnum_status_t, status_ok, FORTNUM_DOMAIN_ERROR
    use fortprob_mixture, only: probability_mean
    implicit none
contains
    subroutine collect(tests)
        type(unittest_type), allocatable, intent(out) :: tests(:)
        allocate(tests(4))
        tests(1) = new_unittest('exact-endpoints-and-single-component', endpoints)
        tests(2) = new_unittest('normalized-high-precision-reference', reference)
        tests(3) = new_unittest('finite-weight-overflow-and-scale-invariance', scaling)
        tests(4) = new_unittest('invalid-input-is-not-a-probability', invalid)
    end subroutine

    subroutine endpoints(error)
        type(error_type), allocatable, intent(out) :: error
        type(fortnum_status_t) :: s
        real(dp) :: probabilities(9), weights(9), value
        weights = 1.0_dp/9.0_dp
        probabilities = 1.0_dp
        call probability_mean(probabilities, weights, value, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call check(error, value == 1.0_dp); if (allocated(error)) return
        probabilities = 0.0_dp
        call probability_mean(probabilities, weights, value, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call check(error, value == 0.0_dp); if (allocated(error)) return
        call probability_mean([0.37_dp], [1e300_dp], value, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call check(error, value == 0.37_dp); if (allocated(error)) return
        call probability_mean([0.0_dp, 1.0_dp, 0.4_dp], &
            [0.0_dp, 0.0_dp, 1e-300_dp], value, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call check(error, value == 0.4_dp)
    end subroutine

    subroutine reference(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: hp = selected_real_kind(30, 300)
        type(fortnum_status_t) :: s
        real(dp) :: p(257), w(257), value, expected
        real(hp) :: numerator, denominator
        integer :: i, j
        do j = 1, 17
            do i = 1, size(p)
                p(i) = real(modulo(37*i + 19*j, 258), dp)/257.0_dp
                w(i) = 10.0_dp**(-300.0_dp + &
                    600.0_dp*real(modulo(23*i + 17*j, 257), dp)/256.0_dp)
            end do
            numerator = sum(real(p, hp)*real(w, hp))
            denominator = sum(real(w, hp))
            expected = real(numerator/denominator, dp)
            call probability_mean(p, w, value, s)
            call check(error, status_ok(s), s%msg); if (allocated(error)) return
            call check(error, value >= 0.0_dp .and. value <= 1.0_dp)
            if (allocated(error)) return
            call check(error, abs(value - expected) < 5e-15_dp)
            if (allocated(error)) return
            call probability_mean(p(size(p):1:-1), w(size(w):1:-1), value, s)
            call check(error, status_ok(s), s%msg); if (allocated(error)) return
            call check(error, value >= 0.0_dp .and. value <= 1.0_dp)
            if (allocated(error)) return
            call check(error, abs(value - expected) < 5e-15_dp)
            if (allocated(error)) return
        end do
    end subroutine

    subroutine scaling(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: hp = selected_real_kind(30, 300)
        type(fortnum_status_t) :: s
        real(dp) :: p(4), w(4), scaled(4), value, expected, factors(3)
        real(hp) :: numerator, denominator
        integer :: i
        p = [0.0_dp, 1.0_dp, 0.25_dp, 0.8_dp]
        w = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
        numerator = sum(real(p, hp)*real(w, hp))
        denominator = sum(real(w, hp))
        expected = real(numerator/denominator, dp)
        factors = [1e-300_dp, 1.0_dp, 1e300_dp]
        do i = 1, 3
            scaled = w*factors(i)
            call probability_mean(p, scaled, value, s)
            call check(error, status_ok(s), s%msg); if (allocated(error)) return
            call check(error, abs(value - expected) < 3e-16_dp)
            if (allocated(error)) return
        end do
        call probability_mean([0.0_dp, 1.0_dp], &
            [huge(1.0_dp), huge(1.0_dp)], value, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call check(error, value == 0.5_dp); if (allocated(error)) return
        call probability_mean([1.0_dp, 0.0_dp], [1e-300_dp, 1e300_dp], value, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call check(error, value == 0.0_dp)
    end subroutine

    subroutine invalid(error)
        type(error_type), allocatable, intent(out) :: error
        type(fortnum_status_t) :: s
        real(dp) :: p(2), w(2), value, nan, inf
        integer :: i
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        inf = ieee_value(0.0_dp, ieee_positive_inf)
        do i = 1, 9
            p = [0.2_dp, 0.7_dp]; w = [1.0_dp, 2.0_dp]
            select case(i)
            case(1)
                p(1) = -epsilon(1.0_dp)
            case(2)
                p(2) = nearest(1.0_dp, 1.0_dp)
            case(3)
                p(1) = nan
            case(4)
                p(2) = inf
            case(5)
                w(1) = -epsilon(1.0_dp)
            case(6)
                w(2) = nan
            case(7)
                w(1) = inf
            case(8)
                w = 0.0_dp
            case(9)
                p(1) = nan; w(1) = 0.0_dp
            end select
            value = 0.0_dp
            call probability_mean(p, w, value, s)
            call check(error, s%code == FORTNUM_DOMAIN_ERROR)
            if (allocated(error)) return
            call check(error, ieee_is_nan(value)); if (allocated(error)) return
        end do
        call probability_mean(p, w(:1), value, s)
        call check(error, s%code == FORTNUM_DOMAIN_ERROR)
        if (allocated(error)) return
        call check(error, ieee_is_nan(value)); if (allocated(error)) return
        call probability_mean(p(:0), w(:0), value, s)
        call check(error, s%code == FORTNUM_DOMAIN_ERROR)
        if (allocated(error)) return
        call check(error, ieee_is_nan(value))
    end subroutine
end module
program test_mixture
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: testsuite_type, new_testsuite, run_testsuite
    use mixture_tests, only: collect
    implicit none
    type(testsuite_type) :: suite
    integer :: stat
    stat = 0
    suite = new_testsuite('mixture', collect)
    call run_testsuite(suite%collect, error_unit, stat)
    if (stat /= 0) error stop 1
end program
