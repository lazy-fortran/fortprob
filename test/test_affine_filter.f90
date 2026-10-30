module affine_filter_tests
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use fortnum_status, only: fortnum_status_t, status_ok
    use fortprob_affine_filter, only: gaussian_interval_t, affine_gaussian_filter_t
    use fortprob_gaussian, only: gaussian_condition, gaussian_logpdf
    use fortprob_regression, only: bayesian_linear_fit_diagonal
    implicit none
contains
    subroutine collect(tests)
        type(unittest_type), allocatable, intent(out) :: tests(:)
        allocate(tests(5))
        tests(1) = new_unittest('dense-joint-whitening-and-prediction', dense)
        tests(2) = new_unittest('normalized-coefficient-evidence', evidence)
        tests(3) = new_unittest('invalid-atomic-and-reset', invalid)
        tests(4) = new_unittest('cancellation-stress', cancellation)
        tests(5) = new_unittest('symmetric-seven-state-recurrence', symmetric_recurrence)
    end subroutine

    subroutine dense(error)
        type(error_type), allocatable, intent(out) :: error
        type(gaussian_interval_t) :: interval
        type(affine_gaussian_filter_t) :: filter
        type(fortnum_status_t) :: s
        real(dp) :: state(2, 22), obs(5, 22), h(5, 2), y(5), f(2, 2), a(2)
        real(dp) :: l(3, 3), q(3, 3), initial(2, 2), noise(5, 5), prior(2, 2)
        real(dp) :: offset, heff(2), variance, input_y, wx(2, 2), wy(2), nv(2)
        real(dp) :: obs_work(22), state_work(2, 22), hprior(5, 2), coef_work(2)
        real(dp) :: scalar_mean(1), scalar_cov(1, 1), hfit(2, 2), beta_xy(2, 2)
        real(dp), allocatable :: zeros(:), fit_mean(:)
        real(dp) :: full_cov(5, 5), full_mean(5), z, reference, expected
        real(dp), allocatable :: m(:), c(:, :), cm(:), cc(:, :)
        real(dp), allocatable :: yy(:, :), xy(:, :), past_y(:), past_h(:, :)
        integer :: i, j, k, nfit, fit(2), start
        logical :: condition
        initial = reshape([1.0_dp, 0.2_dp, 0.2_dp, 0.53_dp], [2, 2])
        prior = reshape([1.3_dp, -0.2_dp, -0.2_dp, 0.7_dp], [2, 2])
        h(:, 1) = [1.0_dp, 0.3_dp, -0.4_dp, 1.1_dp, 0.6_dp]
        h(:, 2) = [0.2_dp, -0.7_dp, 0.8_dp, 0.1_dp, -0.4_dp]
        y = [0.8_dp, 99.0_dp, -0.3_dp, -77.0_dp, 66.0_dp]
        state = 0.0_dp; obs = 0.0_dp
        state(:, 1) = [1.0_dp, 0.2_dp]
        state(:, 2) = [0.0_dp, 0.7_dp]
        call filter%init(initial, 2, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        nfit = 0
        do i = 1, 5
            f = reshape([0.7_dp, 0.1_dp, -0.05_dp, 0.8_dp], [2, 2])
            f = f*(1.0_dp - 0.04_dp*real(i, dp))
            a(1) = 0.8_dp
            a(2) = -0.3_dp + 0.02_dp*real(i, dp)
            l = 0.0_dp
            l(:, 1) = [0.2_dp, 0.04_dp, 0.08_dp]
            l(:, 2) = [0.0_dp, 0.3_dp, -0.06_dp]
            l(3, 3) = 0.1_dp
            q = matmul(l, transpose(l))
            start = 3 + 4*(i - 1)
            obs_work = matmul(a, state)
            obs(i, :) = obs_work
            obs(i, start:start + 2) = l(3, :)
            obs(i, start + 3) = sqrt(0.17_dp)
            state_work = matmul(f, state)
            state = state_work
            state(:, start:start + 2) = l(:2, :)
            noise = matmul(obs, transpose(obs))
            condition = i == 1 .or. i == 3
            input_y = ieee_value(0.0_dp, ieee_quiet_nan)
            if (condition) input_y = y(i)
            call interval%init(f, a, q(:2, :2), q(3, 3), q(:2, 3), 0.17_dp, s)
            call check(error, status_ok(s), s%msg); if (allocated(error)) return
            call filter%step(interval, h(i, :), input_y, condition, offset, heff, variance, s)
            call check(error, status_ok(s), s%msg); if (allocated(error)) return
            if (nfit == 0) then
                call check(error, abs(offset) < 1e-14_dp)
                if (allocated(error)) return
                call check(error, maxval(abs(heff - h(i, :))) < 1e-14_dp)
                if (allocated(error)) return
                call check(error, abs(variance - noise(i, i)) < 1e-14_dp)
                if (allocated(error)) return
            else
                allocate(yy(nfit, nfit), xy(1, nfit), past_y(nfit), past_h(nfit, 2), zeros(nfit))
                zeros = 0.0_dp
                do j = 1, nfit
                    past_y(j) = y(fit(j)); past_h(j, :) = h(fit(j), :)
                    xy(1, j) = noise(i, fit(j))
                    do k = 1, nfit
                        yy(j, k) = noise(fit(j), fit(k))
                    end do
                end do
                scalar_mean = 0.0_dp; scalar_cov(1, 1) = noise(i, i)
                call gaussian_condition(scalar_mean, scalar_cov, &
                    zeros, yy, xy, past_y, m, c, s)
                call check(error, status_ok(s), s%msg); if (allocated(error)) return
                call check(error, abs(offset - m(1)) < 2e-13_dp)
                if (allocated(error)) return
                call check(error, abs(variance - c(1, 1)) < 2e-13_dp)
                if (allocated(error)) return
                do j = 1, 2
                    scalar_mean(1) = h(i, j)
                    call gaussian_condition(scalar_mean, scalar_cov, &
                        past_h(:, j), yy, xy, zeros, m, c, s)
                    call check(error, status_ok(s), s%msg); if (allocated(error)) return
                    call check(error, abs(heff(j) - m(1)) < 2e-13_dp)
                    if (allocated(error)) return
                end do
                deallocate(yy, xy, past_y, past_h, zeros)
            end if
            if (condition) then
                nfit = nfit + 1; fit(nfit) = i
                wx(nfit, :) = heff; wy(nfit) = y(i) - offset; nv(nfit) = variance
            end if
        end do
        call bayesian_linear_fit_diagonal(wx, wy, nv, [0.3_dp, -0.1_dp], prior, &
            cm, cc, s, z)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        hprior = matmul(h, prior)
        full_cov = matmul(hprior, transpose(h))
        full_cov = noise + full_cov
        full_mean = matmul(h, [0.3_dp, -0.1_dp])
        allocate(yy(2, 2), xy(1, 2), past_y(2))
        do j = 1, 2
            past_y(j) = y(fit(j)); xy(1, j) = full_cov(5, fit(j))
            do k = 1, 2
                yy(j, k) = full_cov(fit(j), fit(k))
            end do
        end do
        fit_mean = full_mean(fit)
        hfit = h(fit, :)
        beta_xy = matmul(prior, transpose(hfit))
        call gaussian_condition([0.3_dp, -0.1_dp], prior, fit_mean, yy, beta_xy, &
            past_y, m, c, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call check(error, maxval(abs(cm - m)) < 3e-13_dp)
        if (allocated(error)) return
        call check(error, maxval(abs(cc - c)) < 3e-13_dp)
        if (allocated(error)) return
        call gaussian_logpdf(past_y, fit_mean, yy, reference, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call check(error, abs(z - reference) < 3e-13_dp); if (allocated(error)) return
        scalar_mean(1) = full_mean(5); scalar_cov(1, 1) = full_cov(5, 5)
        call gaussian_condition(scalar_mean, scalar_cov, &
            fit_mean, yy, xy, past_y, m, c, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        expected = offset + dot_product(heff, cm)
        call check(error, abs(expected - m(1)) < 3e-13_dp); if (allocated(error)) return
        coef_work = matmul(cc, heff)
        expected = variance + dot_product(heff, coef_work)
        call check(error, abs(expected - c(1, 1)) < 3e-13_dp)
    end subroutine

    subroutine evidence(error)
        type(error_type), allocatable, intent(out) :: error
        type(fortnum_status_t) :: s
        real(dp) :: x(3, 2), scaled(3, 2), y(3), nv(3), prior(2, 2), pm(2)
        real(dp) :: cov(3, 3), z, reference, scale(2), prior_scaled(2, 2)
        real(dp) :: xprior(3, 2), obs_mean(3), scaled_mean(2)
        real(dp), allocatable :: m(:), c(:, :), old_m(:), old_c(:, :)
        integer :: i, j, k
        x(:, 1) = [1.0_dp, -2.0_dp, 0.7_dp]
        x(:, 2) = [-0.3_dp, 0.1_dp, 0.9_dp]
        y = [0.2_dp, -0.6_dp, 1.1_dp]; nv = [0.4_dp, 0.7_dp, 0.2_dp]
        pm = [0.3_dp, -0.4_dp]
        prior = reshape([1.7_dp, 0.3_dp, 0.3_dp, 0.8_dp], [2, 2])
        xprior = matmul(x, prior)
        cov = matmul(xprior, transpose(x))
        do i = 1, 3
            cov(i, i) = cov(i, i) + nv(i)
        end do
        obs_mean = matmul(x, pm)
        call gaussian_logpdf(y, obs_mean, cov, reference, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call bayesian_linear_fit_diagonal(x, y, nv, pm, prior, old_m, old_c, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        do k = 1, 3
            scale(1) = 10.0_dp**(3*(k - 2))
            scale(2) = 10.0_dp**(-2*(k - 2))
            do i = 1, 2
                scaled(:, i) = x(:, i)*scale(i)
                do j = 1, 2
                    prior_scaled(i, j) = prior(i, j)/(scale(i)*scale(j))
                end do
            end do
            scaled_mean = pm/scale
            call bayesian_linear_fit_diagonal(scaled, y, nv, scaled_mean, prior_scaled, &
                m, c, s, z)
            call check(error, status_ok(s), s%msg); if (allocated(error)) return
            call check(error, abs(z - reference) < 3e-13_dp)
            if (allocated(error)) return
            call check(error, maxval(abs(m*scale - old_m)) < 3e-13_dp)
            if (allocated(error)) return
            do i = 1, 2
                c(:, i) = c(:, i)*scale*scale(i)
            end do
            call check(error, maxval(abs(c - old_c)) < 3e-13_dp)
            if (allocated(error)) return
        end do
        x = 0.0_dp
        call bayesian_linear_fit_diagonal(x, y, nv, pm, prior, m, c, s, z)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call check(error, maxval(abs(m - pm)) < 3e-13_dp)
        if (allocated(error)) return
        call check(error, maxval(abs(c - prior)) < 3e-13_dp)
        if (allocated(error)) return
        cov = 0.0_dp
        do i = 1, 3
            cov(i, i) = nv(i)
        end do
        obs_mean = 0.0_dp
        call gaussian_logpdf(y, obs_mean, cov, reference, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call check(error, abs(z - reference) < 3e-13_dp)
    end subroutine

    subroutine symmetric_recurrence(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: hp = selected_real_kind(30, 300)
        type(gaussian_interval_t) :: interval
        type(affine_gaussian_filter_t) :: filter
        type(fortnum_status_t) :: s
        real(dp) :: f(7, 7), a(7), initial(7, 7), zero(7, 7), fd(7)
        real(dp) :: off, h(1), variance, expected, r
        real(hp) :: precision(6, 6), b(6), power(6), aug(6, 7), solution(6)
        real(hp) :: multiplier, pivot_row(7)
        integer :: i, j, k, row, pivot
        fd = [0.716531310573789271_dp, 0.846481724890614129_dp, &
            0.920044414629323293_dp, 0.972604477116348343_dp, &
            0.986207116743916257_dp, 0.993079612490316066_dp, &
            0.716531310573789271_dp]
        a = [0.850406068278632188_dp, 0.921109650656315226_dp, &
            0.959467024448121042_dp, 0.986238823811457888_dp, &
            0.993087594438032295_dp, 0.996535801394494047_dp, &
            0.850406068278632188_dp]
        f = 0.0_dp; initial = 0.0_dp; zero = 0.0_dp
        precision = 0.0_hp; power = 1.0_hp
        do j = 1, 7
            f(j, j) = fd(j)
        end do
        do j = 1, 6
            initial(j, j) = (400.0_dp/6.0_dp)*4.0_dp
            precision(j, j) = 1.0_hp/real(initial(j, j), hp)
        end do
        r = 0.0100000000000000019_dp
        call interval%init(f, a, zero, 0.0_dp, [0.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], r, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call filter%init(initial, 1, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        do k = 1, 151
            ! Condition the initial coordinates directly, independently of filtering.
            b = real(a(:6), hp)*power
            aug(:, :6) = precision; aug(:, 7) = b
            do j = 1, 6
                pivot = j - 1 + maxloc(abs(aug(j:6, j)), dim=1)
                pivot_row = aug(j, :); aug(j, :) = aug(pivot, :)
                aug(pivot, :) = pivot_row
                do row = j + 1, 6
                    multiplier = aug(row, j)/aug(j, j)
                    aug(row, j:7) = aug(row, j:7) - multiplier*aug(j, j:7)
                end do
            end do
            do j = 6, 1, -1
                solution(j) = (aug(j, 7) - sum(aug(j, j + 1:6) &
                    *solution(j + 1:6)))/aug(j, j)
            end do
            expected = real(real(r, hp) + sum(b*solution), dp)
            call filter%step(interval, [0.0_dp], 0.0_dp, .true., off, h, variance, s)
            call check(error, status_ok(s), s%msg); if (allocated(error)) return
            call check(error, abs(variance - expected) &
                < 1e-9_dp*max(r, abs(expected)))
            if (allocated(error)) return
            do j = 1, 6
                do i = 1, 6
                    precision(i, j) = precision(i, j) + b(i)*b(j)/real(r, hp)
                end do
            end do
            power = power*real(fd(:6), hp)
        end do
    end subroutine

    subroutine cancellation(error)
        type(error_type), allocatable, intent(out) :: error
        type(gaussian_interval_t) :: interval
        type(affine_gaussian_filter_t) :: filter
        type(fortnum_status_t) :: s
        real(dp) :: off, h(1), v, expected
        integer :: i
        call filter%init(reshape([1e12_dp], [1, 1]), 1, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call interval%init(reshape([1.0_dp], [1, 1]), [1.0_dp], &
            reshape([0.0_dp], [1, 1]), 0.0_dp, [0.0_dp], 1e-6_dp, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        do i = 1, 100
            call filter%step(interval, [1.0_dp], 0.0_dp, .true., off, h, v, s)
            call check(error, status_ok(s), s%msg); if (allocated(error)) return
            if (i > 1) then
                expected = 1e-6_dp + 1.0_dp/(1e-12_dp + real(i - 1, dp)/1e-6_dp)
                call check(error, abs(v - expected) < 1e-19_dp)
                if (allocated(error)) return
            end if
        end do
    end subroutine

    subroutine invalid(error)
        type(error_type), allocatable, intent(out) :: error
        type(gaussian_interval_t) :: interval
        type(affine_gaussian_filter_t) :: filter, saved
        type(fortnum_status_t) :: s
        real(dp) :: f(2, 2), q(2, 2), a(2), off, v, heff(1), off2, v2, h2(1), nan
        integer :: k
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        f = reshape([0.8_dp, 0.1_dp, 0.0_dp, 0.6_dp], [2, 2])
        q = reshape([0.2_dp, 0.0_dp, 0.0_dp, 0.1_dp], [2, 2]); a = [0.9_dp, 0.2_dp]
        call filter%init(q, 1, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call interval%init(f, a, q, 0.1_dp, [0.05_dp, 0.01_dp], 0.2_dp, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        saved = filter
        call filter%step(interval, [1.0_dp], nan, .true., off, heff, v, s)
        call check(error, .not. status_ok(s)); if (allocated(error)) return
        call filter%step(interval, [1.0_dp], 0.4_dp, .true., off, heff, v, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call saved%step(interval, [1.0_dp], 0.4_dp, .true., off2, h2, v2, s)
        call check(error, abs(off - off2) + abs(v - v2) + maxval(abs(heff - h2)) < 1e-14_dp)
        if (allocated(error)) return
        ! Reinitialization is the boundary/reset operation.
        call filter%init(q, 1, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call filter%step(interval, [1.0_dp], nan, .false., off, heff, v, s)
        call check(error, status_ok(s), s%msg); if (allocated(error)) return
        call check(error, abs(off) + abs(heff(1) - 1.0_dp) < 1e-14_dp)
        if (allocated(error)) return
        do k = 1, 6
            select case(k)
            case(1)
                q(1, 2) = 0.05_dp
            case(2)
                q(1, 2) = 0.0_dp; q(1, 1) = -1.0_dp
            case(3)
                q(1, 1) = nan
            case(4)
                q(1, 1) = 0.2_dp; f(1, 1) = nan
            case(5)
                f(1, 1) = 0.8_dp; a(1) = nan
            case(6)
                a(1) = 0.9_dp
            end select
            if (k == 6) then
                call interval%init(f, a, q, 0.1_dp, [1.0_dp, 0.0_dp], 0.2_dp, s)
            else
                call interval%init(f, a, q, 0.1_dp, [0.05_dp, 0.01_dp], 0.2_dp, s)
            end if
            call check(error, .not. status_ok(s)); if (allocated(error)) return
        end do
        call interval%init(f, a, q, 0.1_dp, [0.05_dp, 0.01_dp], 0.0_dp, s)
        call check(error, .not. status_ok(s)); if (allocated(error)) return
        call interval%init(f(:, :1), a, q, 0.1_dp, [0.05_dp, 0.01_dp], 0.2_dp, s)
        call check(error, .not. status_ok(s))
    end subroutine
end module
program test_affine_filter
    use, intrinsic :: iso_fortran_env, only: error_unit
    use testdrive, only: testsuite_type, new_testsuite, run_testsuite
    use affine_filter_tests, only: collect
    implicit none
    type(testsuite_type) :: suite
    integer :: stat
    stat = 0
    suite = new_testsuite('affine-filter', collect)
    call run_testsuite(suite%collect, error_unit, stat)
    if (stat /= 0) error stop 1
end program
