module gaussian_tests
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use fortnum_status, only: fortnum_status_t, status_ok
    use fortprob_gaussian, only: gaussian_condition, gaussian_condition_jvp, &
        gaussian_logpdf, linear_propagate, gaussian_intervals, validate_covariance, &
        gaussian_conditioner_t
    use fortprob_regression, only: bayesian_linear_fit, bayesian_linear_predict, &
        bayesian_linear_fit_diagonal, bayesian_linear_predict_diagonal
    implicit none
contains
    subroutine collect(tests)
        type(unittest_type), allocatable, intent(out) :: tests(:)
        allocate(tests(8))
        tests(1) = new_unittest('scalar-conditioning', scalar)
        tests(2) = new_unittest('correlated-regression', regression)
        tests(3) = new_unittest('conditioning-jvp', jvp)
        tests(4) = new_unittest('diagonal-regression', diagonal)
        tests(5) = new_unittest('density-propagation-intervals', helpers)
        tests(6) = new_unittest('invalid-and-singular', invalid)
        tests(7) = new_unittest('diagonal-prediction', diagonal_prediction)
        tests(8) = new_unittest('cached-conditioning-and-refit', cached)
    end subroutine
    subroutine cached(error)
        type(error_type), allocatable, intent(out) :: error
        type(gaussian_conditioner_t) :: conditioner
        type(fortnum_status_t) :: s
        real(dp) :: yy(2, 2), xy(1, 2), alpha(2), inverse(2, 2), solved(2)
        real(dp), allocatable :: mean(:), covariance(:, :)
        yy = reshape([2.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], [2, 2])
        xy = reshape([0.3_dp, -0.2_dp], [1, 2])
        inverse = reshape([1.0_dp, -0.5_dp, -0.5_dp, 2.0_dp], [2, 2])/1.75_dp
        alpha = matmul(inverse, [2.0_dp, -1.0_dp])
        call conditioner%fit([1.0_dp, 2.0_dp], yy, [3.0_dp, 1.0_dp], s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call conditioner%condition([0.7_dp], reshape([1.0_dp], [1, 1]), xy, &
            mean, covariance, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call check(error, abs(mean(1) - 0.7_dp - dot_product(xy(1, :), alpha)) < 1e-14_dp)
        if (allocated(error)) return
        solved = matmul(inverse, xy(1, :))
        call check(error, abs(covariance(1, 1) - 1.0_dp + &
            dot_product(xy(1, :), solved)) < 1e-14_dp)
        if (allocated(error)) return
        yy(1, 2) = 0.7_dp
        call conditioner%fit([1.0_dp, 2.0_dp], yy, [3.0_dp, 1.0_dp], s)
        call check(error, .not. status_ok(s)); if (allocated(error)) return
        call conditioner%condition([0.7_dp], reshape([1.0_dp], [1, 1]), xy, &
            mean, covariance, s)
        call check(error, .not. status_ok(s)); if (allocated(error)) return
        call check(error, .not. allocated(mean) .and. .not. allocated(covariance))
        if (allocated(error)) return
        call conditioner%fit([1.0_dp], reshape([2.0_dp], [1, 1]), [3.0_dp], s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call conditioner%condition([0.0_dp], reshape([1.0_dp], [1, 1]), &
            reshape([0.5_dp], [1, 1]), mean, covariance, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call check(error, abs(mean(1) - 0.5_dp) < 1e-14_dp)
        if (allocated(error)) return
        call check(error, abs(covariance(1, 1) - 0.875_dp) < 1e-14_dp)
    end subroutine cached

    subroutine scalar(error)
        type(error_type), allocatable, intent(out) :: error
        type(fortnum_status_t) :: s
        real(dp), allocatable :: m(:), c(:, :)
        call gaussian_condition([2.0_dp], reshape([4.0_dp], [1, 1]), [1.0_dp], &
            reshape([5.0_dp], [1, 1]), reshape([2.0_dp], [1, 1]), [3.0_dp], m, c, s)
        call check(error, status_ok(s)); if (allocated(error)) return
        call check(error, abs(m(1) - 2.8_dp) < 1e-14_dp); if (allocated(error)) return
        call check(error, abs(c(1, 1) - 3.2_dp) < 1e-14_dp)
    end subroutine
    subroutine regression(error)
        type(error_type), allocatable, intent(out) :: error
        type(fortnum_status_t) :: s
        real(dp) :: eye(2, 2), noise(2, 2), query(2, 2)
        real(dp), allocatable :: m(:), c(:, :), pm(:), pc(:, :), oc(:, :)
        eye = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
        noise = reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
        query = reshape([1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp], [2, 2])
        call bayesian_linear_fit(eye, [1.0_dp, 2.0_dp], noise, [0.0_dp, 0.0_dp], &
            eye, m, c, s)
        call check(error, status_ok(s)); if (allocated(error)) return
        call check(error, maxval(abs(m - [0.125_dp, 0.625_dp])) < 1e-14_dp)
        if (allocated(error)) return
        call check(error, maxval(abs(c - reshape([0.625_dp, 0.125_dp, &
            0.125_dp, 0.625_dp], [2, 2]))) < 1e-14_dp)
        if (allocated(error)) return
        call bayesian_linear_predict(query, m, c, pm, pc, s, eye, oc)
        call check(error, status_ok(s)); if (allocated(error)) return
        call check(error, maxval(abs(pm - [0.75_dp, -0.5_dp])) < 1e-14_dp)
        if (allocated(error)) return
        call check(error, maxval(abs(pc - reshape([1.5_dp, 0.0_dp, &
            0.0_dp, 1.0_dp], [2, 2]))) < 1e-14_dp)
        if (allocated(error)) return
        call check(error, maxval(abs(oc - pc - eye)) < 1e-14_dp)
    end subroutine
    subroutine diagonal_prediction(error)
        type(error_type), allocatable, intent(out) :: error
        type(fortnum_status_t) :: s
        real(dp) :: query(2, 2), cov(2, 2)
        real(dp), allocatable :: m(:), v(:), observation(:)
        query = reshape([1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp], [2, 2])
        cov = reshape([0.625_dp, 0.125_dp, 0.125_dp, 0.625_dp], [2, 2])
        call bayesian_linear_predict_diagonal(query, [0.125_dp, 0.625_dp], &
            cov, m, v, s, [1.0_dp, 2.0_dp], observation)
        call check(error, status_ok(s)); if (allocated(error)) return
        call check(error, maxval(abs(m-[0.75_dp, -0.5_dp])) < 1e-14_dp)
        if (allocated(error)) return
        call check(error, maxval(abs(v-[1.5_dp, 1.0_dp])) < 1e-14_dp)
        if (allocated(error)) return
        call check(error, maxval(abs(observation-[2.5_dp, 3.0_dp])) < 1e-14_dp)
        if (allocated(error)) return
        call bayesian_linear_predict_diagonal(query, [0.125_dp, 0.625_dp], &
            cov, m, v, s, observation_variance=observation)
        call check(error, .not. status_ok(s)); if (allocated(error)) return
        call check(error, .not. allocated(m) .and. .not. allocated(v))
    end subroutine
    subroutine diagonal(error)
        type(error_type), allocatable, intent(out) :: error
        type(fortnum_status_t) :: s
        real(dp) :: design(3, 1), noise(3), y(3)
        real(dp), allocatable :: m(:), c(:, :)
        design = 1.0_dp; noise = [1.0_dp, 2.0_dp, 4.0_dp]
        y = [0.0_dp, 2.0_dp, 4.0_dp]
        call bayesian_linear_fit_diagonal(design, y, noise, [2.0_dp], &
            reshape([4.0_dp], [1, 1]), m, c, s)
        call check(error, status_ok(s)); if (allocated(error)) return
        call check(error, abs(m(1)-1.25_dp) < 1e-14_dp)
        if (allocated(error)) return
        call check(error, abs(c(1, 1)-0.5_dp) < 1e-14_dp)
        if (allocated(error)) return
        noise(1) = 0.0_dp
        call bayesian_linear_fit_diagonal(design, y, noise, [2.0_dp], &
            reshape([4.0_dp], [1, 1]), m, c, s)
        call check(error, .not. status_ok(s))
    end subroutine
    subroutine jvp(error)
        type(error_type), allocatable, intent(out) :: error
        type(fortnum_status_t) :: s
        real(dp) :: x(2), xx(2, 2), y0(2), yy(2, 2), xy(2, 2), y(2)
        real(dp) :: dx(2), dxx(2, 2), dy0(2), dyy(2, 2), dxy(2, 2), dy(2)
        real(dp), parameter :: h = 2e-5_dp
        real(dp), allocatable :: dm(:), dc(:, :), mp(:), cp(:, :), mm(:), cm(:, :)
        x = [0.2_dp, -0.3_dp]; y0 = [0.5_dp, 0.8_dp]; y = [1.1_dp, -0.2_dp]
        xx = reshape([2.0_dp, 0.4_dp, 0.4_dp, 3.0_dp], [2, 2])
        yy = reshape([1.5_dp, 0.2_dp, 0.2_dp, 1.2_dp], [2, 2])
        xy = reshape([0.4_dp, 0.1_dp, 0.2_dp, 0.5_dp], [2, 2])
        dx = [0.1_dp, -0.1_dp]; dy0 = [0.2_dp, 0.3_dp]; dy = [-0.2_dp, 0.1_dp]
        dxx = reshape([0.2_dp, 0.1_dp, 0.1_dp, -0.2_dp], [2, 2])
        dyy = reshape([0.1_dp, 0.02_dp, 0.02_dp, -0.1_dp], [2, 2])
        dxy = reshape([0.03_dp, 0.01_dp, -0.02_dp, 0.04_dp], [2, 2])
        call gaussian_condition_jvp(x, xx, y0, yy, xy, y, dx, dxx, dy0, dyy, dxy, dy, &
            dm, dc, s)
        call check(error, status_ok(s)); if (allocated(error)) return
        call gaussian_condition(x + h*dx, xx + h*dxx, y0 + h*dy0, yy + h*dyy, &
            xy + h*dxy, y + h*dy, mp, cp, s)
        call check(error, status_ok(s)); if (allocated(error)) return
        call gaussian_condition(x - h*dx, xx - h*dxx, y0 - h*dy0, yy - h*dyy, &
            xy - h*dxy, y - h*dy, mm, cm, s)
        call check(error, status_ok(s)); if (allocated(error)) return
        call check(error, maxval(abs(dm - (mp - mm)/(2*h))) < 3e-10_dp)
        if (allocated(error)) return
        call check(error, maxval(abs(dc - (cp - cm)/(2*h))) < 3e-10_dp)
    end subroutine
    subroutine helpers(error)
        type(error_type), allocatable, intent(out) :: error
        type(fortnum_status_t) :: s
        real(dp) :: cov(2, 2), jac(1, 2), density
        real(dp), allocatable :: m(:), c(:, :), lo(:), hi(:)
        cov = reshape([4.0_dp, 1.0_dp, 1.0_dp, 9.0_dp], [2, 2])
        jac = reshape([1.0_dp, -2.0_dp], [1, 2])
        call linear_propagate([2.0_dp, 3.0_dp], cov, jac, m, c, s)
        call check(error, status_ok(s)); if (allocated(error)) return
        call check(error, abs(m(1) + 4.0_dp) < 1e-14_dp); if (allocated(error)) return
        call check(error, abs(c(1, 1) - 36.0_dp) < 1e-14_dp)
        if (allocated(error)) return
        call gaussian_intervals(m, [c(1, 1)], 2.0_dp, lo, hi, s)
        call check(error, status_ok(s)); if (allocated(error)) return
        call check(error, abs(lo(1) + 16.0_dp) < 1e-14_dp)
        if (allocated(error)) return
        call check(error, abs(hi(1) - 8.0_dp) < 1e-14_dp)
        if (allocated(error)) return
        call gaussian_intervals([7.0_dp], [0.0_dp], 2.0_dp, lo, hi, s)
        call check(error, status_ok(s)); if (allocated(error)) return
        call check(error, lo(1) == 7.0_dp .and. hi(1) == 7.0_dp)
        if (allocated(error)) return
        call gaussian_logpdf([1.0_dp, 2.0_dp], [0.0_dp, 0.0_dp], &
            reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2]), density, s)
        call check(error, status_ok(s)); if (allocated(error)) return
        ! det=3, residual^T inverse residual=2.
        call check(error, abs(density + log(2*acos(-1.0_dp)) &
            + 0.5_dp*log(3.0_dp) + 1.0_dp) < 1e-14_dp)
    end subroutine
    subroutine invalid(error)
        type(error_type), allocatable, intent(out) :: error
        type(fortnum_status_t) :: s
        real(dp) :: a(2, 2), nan
        real(dp), allocatable :: m(:), c(:, :)
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        a = reshape([1.0_dp, 0.1_dp, 0.2_dp, 1.0_dp], [2, 2])
        call validate_covariance(a, s)
        call check(error,.not. status_ok(s)); if (allocated(error)) return
        a = reshape([1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp], [2, 2])
        call validate_covariance(a, s)
        call check(error,.not. status_ok(s)); if (allocated(error)) return
        a = 0.0_dp; a(1, 1) = nan
        call validate_covariance(a, s)
        call check(error,.not. status_ok(s)); if (allocated(error)) return
        call gaussian_condition([0.0_dp], reshape([1.0_dp], [1, 1]), [0.0_dp], &
            reshape([1.0_dp], [1, 1]), reshape([2.0_dp], [1, 1]), [1.0_dp], m, c, s)
        call check(error,.not. status_ok(s)); if (allocated(error)) return
        call check(error,.not. allocated(m) .and. .not. allocated(c))
        if (allocated(error)) return
        call gaussian_condition([0.0_dp], reshape([0.0_dp], [1, 1]), [0.0_dp], &
            reshape([1.0_dp], [1, 1]), reshape([0.0_dp], [1, 1]), [1.0_dp], m, c, s)
        call check(error, status_ok(s)); if (allocated(error)) return
        call check(error, c(1, 1) == 0.0_dp); if (allocated(error)) return
        call gaussian_condition([0.0_dp], reshape([1.0_dp], [1, 1]), [0.0_dp], &
            reshape([0.0_dp], [1, 1]), reshape([0.0_dp], [1, 1]), [1.0_dp], m, c, s)
        call check(error,.not. status_ok(s)); if (allocated(error)) return
        call gaussian_condition([0.0_dp], reshape([1.0_dp], [1, 1]), [0.0_dp], &
            reshape([1.0_dp], [1, 1]), reshape([0.0_dp], [1, 1]), [nan], m, c, s)
        call check(error,.not. status_ok(s))
    end subroutine
end module
program test_gaussian
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: testsuite_type, new_testsuite, run_testsuite
    use gaussian_tests, only: collect
    implicit none
    type(testsuite_type) :: suite
    integer :: stat
    stat = 0
    suite = new_testsuite('fortprob', collect)
    call run_testsuite(suite%collect, error_unit, stat)
    if (stat /= 0) error stop 1
end program
