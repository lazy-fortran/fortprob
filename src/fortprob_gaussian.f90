module fortprob_gaussian
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, status_ok, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    use fortnum_cholesky, only: cholesky_factorization_t
    use fortnum_symmetric_eigen, only: symmetric_eigen
    implicit none
    private
    public :: gaussian_condition, gaussian_condition_jvp, gaussian_logpdf
    public :: validate_covariance, linear_propagate, gaussian_intervals
    interface
        pure subroutine fortprob_interval_leaf(mean, variance, z, lower, upper, ds)
            import dp
            real(dp), intent(in) :: mean, variance, z
            real(dp), intent(out) :: lower, upper, ds
        end subroutine
        pure subroutine fortprob_logpdf_leaf(n, logdet, quadratic, density)
            import dp
            real(dp), intent(in) :: n, logdet, quadratic
            real(dp), intent(out) :: density
        end subroutine
    end interface
contains
    subroutine validate_covariance(cov, status, positive_definite)
        real(dp), intent(in) :: cov(:, :)
        type(fortnum_status_t), intent(out) :: status
        logical, intent(in), optional :: positive_definite
        type(cholesky_factorization_t) :: factor
        real(dp), allocatable :: values(:), vectors(:, :)
        real(dp) :: tolerance
        integer :: n, info
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid covariance')
        n = size(cov, 1)
        if (n < 1 .or. size(cov, 2) /= n) return
        if (.not. all(ieee_is_finite(cov))) return
        tolerance = 128.0_dp*epsilon(1.0_dp)*real(n, dp)*maxval(abs(cov))
        if (maxval(abs(cov - transpose(cov))) > tolerance) return
        if (present(positive_definite)) then
            if (positive_definite) then
                call factor%factorize(cov, status)
                return
            end if
        end if
        allocate (values(n), vectors(n, n))
        call symmetric_eigen(cov, values, vectors, info)
        if (info /= 0) then
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                'covariance eigenvalue validation failed')
            return
        end if
        if (minval(values) < -tolerance) return
        call status_set(status, FORTNUM_OK, '')
    end subroutine

    subroutine gaussian_condition(mean_x, cov_xx, mean_y, cov_yy, cov_xy, y, &
            post_mean, post_cov, status)
        real(dp), intent(in) :: mean_x(:), cov_xx(:, :), mean_y(:), cov_yy(:, :)
        real(dp), intent(in) :: cov_xy(:, :), y(:)
        real(dp), allocatable, intent(out) :: post_mean(:), post_cov(:, :)
        type(fortnum_status_t), intent(out) :: status
        type(cholesky_factorization_t) :: factor
        real(dp), allocatable :: alpha(:), solved(:, :), candidate(:, :)
        integer :: nx, ny
        nx = size(mean_x); ny = size(mean_y)
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'conditioning shape/data')
        if (nx < 1 .or. ny < 1) return
        if (any(shape(cov_xx) /= [nx, nx])) return
        if (any(shape(cov_yy) /= [ny, ny])) return
        if (any(shape(cov_xy) /= [nx, ny]) .or. size(y) /= ny) return
        if (.not. all(ieee_is_finite(mean_x))) return
        if (.not. all(ieee_is_finite(mean_y))) return
        if (.not. all(ieee_is_finite(y))) return
        if (.not. all(ieee_is_finite(cov_xy))) return
        call validate_covariance(cov_xx, status)
        if (.not. status_ok(status)) return
        call validate_covariance(cov_yy, status, .true.)
        if (.not. status_ok(status)) return
        call factor%factorize(cov_yy, status)
        if (.not. status_ok(status)) return
        alpha = y - mean_y
        solved = transpose(cov_xy)
        call factor%solve(alpha, status)
        if (.not. status_ok(status)) return
        call factor%solve(solved, status)
        if (.not. status_ok(status)) return
        candidate = cov_xx - matmul(cov_xy, solved)
        call validate_covariance(candidate, status)
        if (.not. status_ok(status)) return
        post_mean = mean_x + matmul(cov_xy, alpha)
        post_cov = candidate
        if (.not. all(ieee_is_finite(post_mean))) then
            deallocate (post_mean, post_cov)
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'conditioning overflow')
        end if
    end subroutine

    subroutine gaussian_condition_jvp(mean_x, cov_xx, mean_y, cov_yy, cov_xy, y, &
            dmean_x, dcov_xx, dmean_y, dcov_yy, dcov_xy, dy, &
            dpost_mean, dpost_cov, status)
        real(dp), intent(in) :: mean_x(:), cov_xx(:, :), mean_y(:), cov_yy(:, :)
        real(dp), intent(in) :: cov_xy(:, :), y(:), dmean_x(:), dcov_xx(:, :)
        real(dp), intent(in) :: dmean_y(:), dcov_yy(:, :), dcov_xy(:, :), dy(:)
        real(dp), allocatable, intent(out) :: dpost_mean(:), dpost_cov(:, :)
        type(fortnum_status_t), intent(out) :: status
        type(cholesky_factorization_t) :: factor
        real(dp), allocatable :: m(:), c(:, :), alpha(:), beta(:, :), da(:), db(:, :)
        call gaussian_condition(mean_x, cov_xx, mean_y, cov_yy, cov_xy, y, m, c, status)
        if (.not. status_ok(status)) return
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'conditioning JVP shape/data')
        if (size(dmean_x) /= size(mean_x) .or. size(dy) /= size(y)) return
        if (size(dmean_y) /= size(mean_y)) return
        if (any(shape(dcov_xx) /= shape(cov_xx))) return
        if (any(shape(dcov_yy) /= shape(cov_yy))) return
        if (any(shape(dcov_xy) /= shape(cov_xy))) return
        if (.not. all(ieee_is_finite(dmean_x))) return
        if (.not. all(ieee_is_finite(dmean_y))) return
        if (.not. all(ieee_is_finite(dy))) return
        if (.not. all(ieee_is_finite(dcov_xx))) return
        if (.not. all(ieee_is_finite(dcov_yy))) return
        if (.not. all(ieee_is_finite(dcov_xy))) return
        call factor%factorize(cov_yy, status)
        if (.not. status_ok(status)) return
        alpha = y - mean_y; beta = transpose(cov_xy)
        call factor%solve(alpha, status)
        if (.not. status_ok(status)) return
        call factor%solve(beta, status)
        if (.not. status_ok(status)) return
        da = dy - dmean_y - matmul(dcov_yy, alpha)
        db = transpose(dcov_xy) - matmul(dcov_yy, beta)
        call factor%solve(da, status)
        if (.not. status_ok(status)) return
        call factor%solve(db, status)
        if (.not. status_ok(status)) return
        dpost_mean = dmean_x + matmul(dcov_xy, alpha) + matmul(cov_xy, da)
        dpost_cov = dcov_xx - matmul(dcov_xy, beta) - matmul(cov_xy, db)
    end subroutine

    subroutine gaussian_logpdf(value, mean, cov, log_density, status)
        real(dp), intent(in) :: value(:), mean(:), cov(:, :)
        real(dp), intent(out) :: log_density
        type(fortnum_status_t), intent(out) :: status
        type(cholesky_factorization_t) :: factor
        real(dp), allocatable :: residual(:), solved(:)
        real(dp) :: logdet
        log_density = 0.0_dp
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'Gaussian density shape/data')
        if (size(value) /= size(mean)) return
        if (any(shape(cov) /= [size(mean), size(mean)])) return
        if (.not. all(ieee_is_finite(value))) return
        if (.not. all(ieee_is_finite(mean))) return
        call validate_covariance(cov, status, .true.)
        if (.not. status_ok(status)) return
        call factor%factorize(cov, status)
        if (.not. status_ok(status)) return
        residual = value - mean; solved = residual
        call factor%solve(solved, status)
        if (.not. status_ok(status)) return
        call factor%log_determinant(logdet, status)
        if (.not. status_ok(status)) return
        call fortprob_logpdf_leaf(real(size(mean), dp), logdet, &
            dot_product(residual, solved), log_density)
        if (.not. ieee_is_finite(log_density)) &
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'Gaussian density overflow')
    end subroutine

    subroutine linear_propagate(mean, cov, jacobian, out_mean, out_cov, status)
        real(dp), intent(in) :: mean(:), cov(:, :), jacobian(:, :)
        real(dp), allocatable, intent(out) :: out_mean(:), out_cov(:, :)
        type(fortnum_status_t), intent(out) :: status
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'linear propagation shape/data')
        if (size(jacobian, 2) /= size(mean) .or. size(jacobian, 1) < 1) return
        if (any(shape(cov) /= [size(mean), size(mean)])) return
        if (.not. all(ieee_is_finite(mean))) return
        if (.not. all(ieee_is_finite(jacobian))) return
        call validate_covariance(cov, status)
        if (.not. status_ok(status)) return
        out_mean = matmul(jacobian, mean)
        out_cov = matmul(matmul(jacobian, cov), transpose(jacobian))
        if (.not. all(ieee_is_finite(out_mean)) .or. &
            .not. all(ieee_is_finite(out_cov))) then
            deallocate (out_mean, out_cov)
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'linear propagation overflow')
        end if
    end subroutine

    subroutine gaussian_intervals(mean, variance, z, lower, upper, status)
        real(dp), intent(in) :: mean(:), variance(:), z
        real(dp), allocatable, intent(out) :: lower(:), upper(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp) :: unused
        integer :: i
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'interval shape/data')
        if (size(mean) /= size(variance) .or. size(mean) < 1) return
        if (.not. all(ieee_is_finite(mean))) return
        if (.not. all(ieee_is_finite(variance))) return
        if (any(variance < 0.0_dp)) return
        if (.not. ieee_is_finite(z)) return
        if (z < 0.0_dp) return
        allocate (lower(size(mean)), upper(size(mean)))
        do i = 1, size(mean)
            if (variance(i) == 0.0_dp) then
                lower(i) = mean(i); upper(i) = mean(i)
            else
                call fortprob_interval_leaf(mean(i), variance(i), z, &
                    lower(i), upper(i), unused)
            end if
        end do
        call status_set(status, FORTNUM_OK, '')
    end subroutine
end module fortprob_gaussian
