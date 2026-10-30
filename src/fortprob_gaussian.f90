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
    type, public :: gaussian_conditioner_t
        private
        type(cholesky_factorization_t) :: factor
        real(dp), allocatable :: alpha(:)
        logical :: fitted = .false.
    contains
        procedure :: fit => conditioner_fit
        procedure :: condition => conditioner_condition
    end type gaussian_conditioner_t
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
    subroutine validate_covariance(cov, status, positive_definite, reference_scale)
        real(dp), intent(in) :: cov(:, :)
        type(fortnum_status_t), intent(out) :: status
        logical, intent(in), optional :: positive_definite
        real(dp), intent(in), optional :: reference_scale
        type(cholesky_factorization_t) :: factor
        real(dp), allocatable :: values(:), vectors(:, :)
        real(dp) :: tolerance
        integer :: n, info
        call covariance_structure(cov, status, tolerance, reference_scale)
        if (.not. status_ok(status)) return
        if (present(positive_definite)) then
            if (positive_definite) then
                call factor%factorize(cov, status)
                return
            end if
        end if
        n = size(cov, 1)
        allocate (values(n), vectors(n, n))
        call symmetric_eigen(cov, values, vectors, info)
        if (info /= 0) then
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                'covariance eigenvalue validation failed')
            return
        end if
        if (minval(values) < -tolerance) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid covariance')
            return
        end if
        call status_set(status, FORTNUM_OK, '')
    end subroutine

    subroutine covariance_structure(cov, status, tolerance, reference_scale)
        real(dp), intent(in) :: cov(:, :)
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(out) :: tolerance
        real(dp), intent(in), optional :: reference_scale
        integer :: n
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid covariance')
        n = size(cov, 1)
        if (n < 1 .or. size(cov, 2) /= n) return
        if (.not. all(ieee_is_finite(cov))) return
        tolerance = maxval(abs(cov))
        if (present(reference_scale)) then
            if (.not. ieee_is_finite(reference_scale)) return
            if (reference_scale < 0.0_dp) return
            tolerance = max(tolerance, reference_scale)
        end if
        tolerance = 128.0_dp*epsilon(1.0_dp)*real(n, dp)*tolerance
        if (maxval(abs(cov - transpose(cov))) > tolerance) return
        call status_set(status, FORTNUM_OK, '')
    end subroutine covariance_structure

    subroutine conditioner_fit(self, mean_y, cov_yy, y, status)
        class(gaussian_conditioner_t), intent(out) :: self
        real(dp), intent(in) :: mean_y(:), cov_yy(:, :), y(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp) :: tolerance
        integer :: n
        n = size(mean_y)
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'conditioning shape/data')
        if (n < 1 .or. size(y) /= n) return
        if (size(cov_yy, 1) /= n .or. size(cov_yy, 2) /= n) return
        if (.not. all(ieee_is_finite(mean_y))) return
        if (.not. all(ieee_is_finite(y))) return
        call covariance_structure(cov_yy, status, tolerance)
        if (.not. status_ok(status)) return
        call self%factor%factorize(cov_yy, status)
        if (.not. status_ok(status)) return
        self%alpha = y - mean_y
        call self%factor%solve(self%alpha, status)
        if (.not. status_ok(status)) return
        if (.not. all(ieee_is_finite(self%alpha))) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'conditioning overflow')
            return
        end if
        self%fitted = .true.
    end subroutine conditioner_fit

    subroutine conditioner_condition(self, mean_x, cov_xx, cov_xy, &
            post_mean, post_cov, status)
        class(gaussian_conditioner_t), intent(in) :: self
        real(dp), intent(in) :: mean_x(:), cov_xx(:, :), cov_xy(:, :)
        real(dp), allocatable, intent(out) :: post_mean(:), post_cov(:, :)
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: solved(:, :), candidate(:, :)
        integer :: nx
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'conditioning fit required')
        if (.not. self%fitted) return
        nx = size(mean_x)
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'conditioning shape/data')
        if (nx < 1) return
        if (size(cov_xx, 1) /= nx .or. size(cov_xx, 2) /= nx) return
        if (size(cov_xy, 1) /= nx .or. size(cov_xy, 2) /= size(self%alpha)) return
        if (.not. all(ieee_is_finite(mean_x))) return
        if (.not. all(ieee_is_finite(cov_xy))) return
        call validate_covariance(cov_xx, status)
        if (.not. status_ok(status)) return
        solved = transpose(cov_xy)
        call self%factor%solve(solved, status)
        if (.not. status_ok(status)) return
        candidate = matmul(cov_xy, solved)
        candidate = cov_xx - candidate
        call validate_covariance(candidate, status, &
            reference_scale=maxval(abs(cov_xx)))
        if (.not. status_ok(status)) return
        post_mean = matmul(cov_xy, self%alpha)
        post_mean = mean_x + post_mean
        post_cov = candidate
        if (.not. all(ieee_is_finite(post_mean))) then
            deallocate (post_mean, post_cov)
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'conditioning overflow')
        end if
    end subroutine conditioner_condition

    subroutine gaussian_condition(mean_x, cov_xx, mean_y, cov_yy, cov_xy, y, &
            post_mean, post_cov, status)
        real(dp), intent(in) :: mean_x(:), cov_xx(:, :), mean_y(:), cov_yy(:, :)
        real(dp), intent(in) :: cov_xy(:, :), y(:)
        real(dp), allocatable, intent(out) :: post_mean(:), post_cov(:, :)
        type(fortnum_status_t), intent(out) :: status
        type(gaussian_conditioner_t) :: conditioner
        call conditioner%fit(mean_y, cov_yy, y, status)
        if (.not. status_ok(status)) return
        call conditioner%condition(mean_x, cov_xx, cov_xy, post_mean, post_cov, status)
    end subroutine gaussian_condition

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
        real(dp), allocatable :: vector_work(:), matrix_work(:, :)
        call gaussian_condition(mean_x, cov_xx, mean_y, cov_yy, cov_xy, y, m, c, status)
        if (.not. status_ok(status)) return
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'conditioning JVP shape/data')
        if (size(dmean_x) /= size(mean_x) .or. size(dy) /= size(y)) return
        if (size(dmean_y) /= size(mean_y)) return
        if ((size(dcov_xx, 1) /= size(cov_xx, 1) .or. &
            size(dcov_xx, 2) /= size(cov_xx, 2))) return
        if ((size(dcov_yy, 1) /= size(cov_yy, 1) .or. &
            size(dcov_yy, 2) /= size(cov_yy, 2))) return
        if ((size(dcov_xy, 1) /= size(cov_xy, 1) .or. &
            size(dcov_xy, 2) /= size(cov_xy, 2))) return
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
        da = matmul(dcov_yy, alpha)
        da = dy - dmean_y - da
        db = matmul(dcov_yy, beta)
        db = transpose(dcov_xy) - db
        call factor%solve(da, status)
        if (.not. status_ok(status)) return
        call factor%solve(db, status)
        if (.not. status_ok(status)) return
        dpost_mean = matmul(dcov_xy, alpha)
        vector_work = matmul(cov_xy, da)
        dpost_mean = dmean_x + dpost_mean + vector_work
        dpost_cov = matmul(dcov_xy, beta)
        matrix_work = matmul(cov_xy, db)
        dpost_cov = dcov_xx - dpost_cov - matrix_work
        if (.not. all(ieee_is_finite(dpost_mean)) .or. &
            .not. all(ieee_is_finite(dpost_cov))) then
            deallocate(dpost_mean, dpost_cov)
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'conditioning JVP overflow')
        end if
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
        if ((size(cov, 1) /= size(mean) .or. size(cov, 2) /= size(mean))) return
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
        real(dp), allocatable :: work(:, :)
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'linear propagation shape/data')
        if (size(jacobian, 2) /= size(mean) .or. size(jacobian, 1) < 1) return
        if ((size(cov, 1) /= size(mean) .or. size(cov, 2) /= size(mean))) return
        if (.not. all(ieee_is_finite(mean))) return
        if (.not. all(ieee_is_finite(jacobian))) return
        call validate_covariance(cov, status)
        if (.not. status_ok(status)) return
        out_mean = matmul(jacobian, mean)
        work = matmul(jacobian, cov)
        out_cov = matmul(work, transpose(jacobian))
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
        if (.not. all(ieee_is_finite(lower)) .or. &
            .not. all(ieee_is_finite(upper))) then
            deallocate(lower, upper)
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'interval overflow')
            return
        end if
        call status_set(status, FORTNUM_OK, '')
    end subroutine
end module fortprob_gaussian
