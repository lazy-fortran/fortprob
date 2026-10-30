module laplace_tests
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use fortnum_status, only: fortnum_status_t, status_ok
    use fortprob_laplace, only: laplace_at_mode
    use fortprob_gaussian, only: gaussian_logpdf
    implicit none
contains
    subroutine collect(tests)
        type(unittest_type), allocatable, intent(out) :: tests(:)
        allocate(tests(4))
        tests(1) = new_unittest('normalized-gaussian', gaussian)
        tests(2) = new_unittest('rescaled-regression', regression)
        tests(3) = new_unittest('nonlinear-quadrature', nonlinear)
        tests(4) = new_unittest('invalid-curvature', invalid)
    end subroutine
    subroutine gaussian(error)
        type(error_type), allocatable, intent(out) :: error
        type(fortnum_status_t) :: s
        real(dp), allocatable :: c(:, :)
        real(dp) :: h(2, 2), expected(2, 2), z, logjoint
        h = reshape([2.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], [2, 2])
        expected = reshape([1.0_dp, -0.5_dp, -0.5_dp, 2.0_dp], [2, 2])/1.75_dp
        logjoint = -log(2.0_dp*acos(-1.0_dp)) + 0.5_dp*log(1.75_dp)
        call laplace_at_mode(logjoint, h, c, z, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call check(error, abs(z) < 1e-14_dp); if (allocated(error)) return
        call check(error, maxval(abs(c - expected)) < 1e-14_dp)
    end subroutine
    subroutine regression(error)
        type(error_type), allocatable, intent(out) :: error
        type(fortnum_status_t) :: s
        real(dp), allocatable :: c(:, :)
        real(dp) :: x(3), y(3), noise(3, 3), marginal(3, 3), h(1, 1)
        real(dp) :: mode, lp, ll, reference, z, scale
        real(dp) :: obs_mean(3), prior_mode(1), prior_mean(1), prior_cov(1, 1)
        integer :: i, j, k
        x = [1.0_dp, -2.0_dp, 0.7_dp]; y = [0.2_dp, -0.6_dp, 1.1_dp]
        noise = 0.0_dp
        do i = 1, 3
            noise(i, i) = 0.4_dp
            do j = 1, 3
                marginal(i, j) = 1.7_dp*x(i)*x(j)
            end do
        end do
        marginal = marginal + noise
        obs_mean = 0.3_dp*x
        call gaussian_logpdf(y, obs_mean, marginal, reference, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        do k = 1, 3
            scale = 10.0_dp**(3*(k - 2))
            h(1, 1) = (1.0_dp/1.7_dp + sum(x*x)/0.4_dp)*scale**2
            mode = (0.3_dp/1.7_dp + sum(x*y)/0.4_dp)/h(1, 1)*scale
            prior_mode = mode; prior_mean = 0.3_dp/scale
            prior_cov(1, 1) = 1.7_dp/scale**2
            call gaussian_logpdf(prior_mode, prior_mean, prior_cov, lp, s)
            call check(error, status_ok(s), s%msg); if (allocated(error)) return
            obs_mean = x*scale*mode
            call gaussian_logpdf(y, obs_mean, noise, ll, s)
            call check(error, status_ok(s), s%msg); if (allocated(error)) return
            call laplace_at_mode(lp + ll, h, c, z, s)
            call check(error, status_ok(s), s%msg); if (allocated(error)) return
            call check(error, abs(z - reference) < 2e-13_dp)
            if (allocated(error)) return
        end do
    end subroutine
    recursive function integrate(a, b, fa, fm, fb, whole, tolerance) result(value)
        real(dp), intent(in) :: a, b, fa, fm, fb, whole, tolerance
        real(dp) :: value, m, fl, fr, left, right
        m = (a + b)/2.0_dp
        fl = integrand((a + m)/2.0_dp); fr = integrand((m + b)/2.0_dp)
        left = (m - a)*(fa + 4.0_dp*fl + fm)/6.0_dp
        right = (b - m)*(fm + 4.0_dp*fr + fb)/6.0_dp
        if (abs(left + right - whole) < 15.0_dp*tolerance) then
            value = left + right + (left + right - whole)/15.0_dp
        else
            value = integrate(a, m, fa, fl, fm, left, tolerance/2.0_dp) + &
                integrate(m, b, fm, fr, fb, right, tolerance/2.0_dp)
        end if
    end function
    function integrand(x) result(value)
        real(dp), intent(in) :: x
        real(dp) :: value
        ! N(x;0,0.25) N(1;x**2,1): full curvature 2, Gauss-Newton 4.
        value = exp(-0.5_dp*(4.0_dp*x*x + (1.0_dp - x*x)**2))/ &
            acos(-1.0_dp)
    end function
    subroutine nonlinear(error)
        type(error_type), allocatable, intent(out) :: error
        type(fortnum_status_t) :: s
        real(dp), allocatable :: c(:, :)
        real(dp) :: z, exact, a, b, fa, fm, fb, whole, relative_error
        a = -10.0_dp; b = 10.0_dp
        fa = integrand(a); fm = integrand(0.0_dp); fb = integrand(b)
        whole = (b - a)*(fa + 4.0_dp*fm + fb)/6.0_dp
        exact = integrate(a, b, fa, fm, fb, whole, 1e-12_dp)
        call laplace_at_mode(log(fm), reshape([2.0_dp], [1, 1]), c, z, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        relative_error = exp(z)/exact - 1.0_dp
        call check(error, relative_error > 0.18_dp .and. relative_error < 0.20_dp)
        if (allocated(error)) return
        call check(error, abs(c(1, 1) - 0.5_dp) < 1e-14_dp)
    end subroutine
    subroutine invalid(error)
        type(error_type), allocatable, intent(out) :: error
        type(fortnum_status_t) :: s
        real(dp), allocatable :: c(:, :)
        real(dp) :: z, h(2, 2), bad
        integer :: i
        bad = ieee_value(0.0_dp, ieee_quiet_nan)
        do i = 1, 5
            h = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
            select case(i)
            case(1)
                h(1, 2) = 0.2_dp
            case(2)
                h(1, 1) = bad
            case(3)
                h(1, 1) = 0.0_dp
            case(4)
                h(1, 1) = -1.0_dp
            end select
            if (i == 5) then
                call laplace_at_mode(bad, h, c, z, s)
            else
                call laplace_at_mode(0.0_dp, h, c, z, s)
            end if
            call check(error, .not. status_ok(s)); if (allocated(error)) return
            call check(error, .not. allocated(c)); if (allocated(error)) return
        end do
        call laplace_at_mode(0.0_dp, h(:, 1:1), c, z, s)
        call check(error, .not. status_ok(s))
    end subroutine
end module
program test_laplace
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: testsuite_type, new_testsuite, run_testsuite
    use laplace_tests, only: collect
    implicit none
    type(testsuite_type) :: suite
    integer :: stat
    stat = 0
    suite = new_testsuite('laplace', collect)
    call run_testsuite(suite%collect, error_unit, stat)
    if (stat /= 0) error stop 1
end program
