module fortprob_regression
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, status_set, status_ok, &
        FORTNUM_DOMAIN_ERROR, FORTNUM_OK
    use fortnum_cholesky, only: cholesky_factorization_t
    use fortprob_gaussian, only: gaussian_condition, validate_covariance, &
        linear_propagate
    implicit none
    private
    public :: bayesian_linear_fit, bayesian_linear_predict, &
        bayesian_linear_fit_diagonal, bayesian_linear_predict_diagonal
    interface
        pure subroutine fortprob_logpdf_leaf(n, logdet, quadratic, density)
            import dp
            real(dp), intent(in) :: n, logdet, quadratic
            real(dp), intent(out) :: density
        end subroutine
    end interface
contains
    subroutine bayesian_linear_fit(design, y, noise_cov, prior_mean, prior_cov, &
            coef_mean, coef_cov, status)
        real(dp), intent(in) :: design(:, :), y(:), noise_cov(:, :)
        real(dp), intent(in) :: prior_mean(:), prior_cov(:, :)
        real(dp), allocatable, intent(out) :: coef_mean(:), coef_cov(:, :)
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: cross(:, :), obs_cov(:, :), obs_mean(:)
        integer :: n, p
        n = size(design, 1); p = size(design, 2)
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'regression shape/data')
        if (n < 1 .or. p < 1 .or. size(y) /= n) return
        if (size(prior_mean) /= p) return
        if ((size(prior_cov, 1) /= p .or. size(prior_cov, 2) /= p)) return
        if ((size(noise_cov, 1) /= n .or. size(noise_cov, 2) /= n)) return
        if (.not. all(ieee_is_finite(design))) return
        call validate_covariance(prior_cov, status, .true.)
        if (.not. status_ok(status)) return
        call validate_covariance(noise_cov, status, .true.)
        if (.not. status_ok(status)) return
        cross = matmul(prior_cov, transpose(design))
        obs_cov = matmul(design, cross)
        obs_cov = obs_cov + noise_cov
        obs_mean = matmul(design, prior_mean)
        call gaussian_condition(prior_mean, prior_cov, obs_mean, obs_cov, cross, y, &
            coef_mean, coef_cov, status)
    end subroutine

    subroutine bayesian_linear_fit_diagonal(design, y, noise_variance, &
            prior_mean, prior_cov, coef_mean, coef_cov, status, log_evidence)
        real(dp), intent(in) :: design(:, :), y(:), noise_variance(:)
        real(dp), intent(in) :: prior_mean(:), prior_cov(:, :)
        real(dp), allocatable, intent(out) :: coef_mean(:), coef_cov(:, :)
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(out), optional :: log_evidence
        type(cholesky_factorization_t) :: factor
        real(dp), allocatable :: precision(:, :), weighted(:, :), rhs(:)
        real(dp), allocatable :: identity(:, :), work(:, :), weighted_y(:), prior_rhs(:)
        real(dp), allocatable :: prior_precision(:, :), delta(:), residual(:), prior_work(:)
        real(dp) :: prior_logdet, precision_logdet, quadratic
        integer :: n, p, i
        if (present(log_evidence)) log_evidence = 0.0_dp
        n = size(design, 1); p = size(design, 2)
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'diagonal regression shape/data')
        if (n < 1 .or. p < 1 .or. size(y) /= n) return
        if (size(noise_variance) /= n .or. size(prior_mean) /= p) return
        if ((size(prior_cov, 1) /= p .or. size(prior_cov, 2) /= p)) return
        if (.not. all(ieee_is_finite(design))) return
        if (.not. all(ieee_is_finite(y))) return
        if (.not. all(ieee_is_finite(prior_mean))) return
        if (.not. all(ieee_is_finite(noise_variance))) return
        if (any(noise_variance <= 0.0_dp)) return
        call validate_covariance(prior_cov, status, .true.)
        if (.not. status_ok(status)) return
        allocate(identity(p, p), weighted(n, p))
        identity = 0.0_dp
        do i = 1, p
            identity(i, i) = 1.0_dp
        end do
        precision = identity
        call factor%factorize(prior_cov, status)
        if (.not. status_ok(status)) return
        if (present(log_evidence)) then
            call factor%log_determinant(prior_logdet, status)
            if (.not. status_ok(status)) return
        end if
        call factor%solve(precision, status)
        if (.not. status_ok(status)) return
        if (present(log_evidence)) prior_precision = precision
        ! Assemble coefficient-space normal equations without an n-by-n array.
        do i = 1, p
            weighted(:, i) = design(:, i)/noise_variance
        end do
        prior_rhs = matmul(precision, prior_mean)
        weighted_y = y/noise_variance
        rhs = matmul(transpose(design), weighted_y)
        rhs = rhs + prior_rhs
        work = matmul(transpose(design), weighted)
        precision = precision + work
        if (.not. all(ieee_is_finite(precision)) .or. &
            .not. all(ieee_is_finite(rhs))) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'diagonal regression overflow')
            return
        end if
        call factor%factorize(precision, status)
        if (.not. status_ok(status)) return
        if (present(log_evidence)) then
            call factor%log_determinant(precision_logdet, status)
            if (.not. status_ok(status)) return
        end if
        call factor%solve(rhs, status)
        if (.not. status_ok(status)) return
        call factor%solve(identity, status)
        if (.not. status_ok(status)) return
        if (.not. all(ieee_is_finite(rhs)) .or. &
            .not. all(ieee_is_finite(identity))) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'diagonal regression solve overflow')
            return
        end if
        if (present(log_evidence)) then
            residual = matmul(design, rhs)
            residual = y - residual
            delta = rhs - prior_mean
            prior_work = matmul(prior_precision, delta)
            quadratic = sum(residual*(residual/noise_variance)) &
                + dot_product(delta, prior_work)
            call fortprob_logpdf_leaf(real(n, dp), sum(log(noise_variance)) &
                + prior_logdet + precision_logdet, quadratic, log_evidence)
            if (.not. ieee_is_finite(log_evidence)) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, 'regression evidence overflow')
                return
            end if
        end if
        coef_mean = rhs
        coef_cov = identity
    end subroutine

    subroutine bayesian_linear_predict_diagonal(design, coef_mean, coef_cov, &
            mean, latent_variance, status, noise_variance, observation_variance)
        real(dp), intent(in) :: design(:, :), coef_mean(:), coef_cov(:, :)
        real(dp), allocatable, intent(out) :: mean(:), latent_variance(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(in), optional :: noise_variance(:)
        real(dp), allocatable, intent(out), optional :: observation_variance(:)
        real(dp), allocatable :: work(:, :)
        integer :: n, p, i
        n = size(design, 1); p = size(design, 2)
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'diagonal prediction shape/data')
        if (n < 1 .or. p < 1 .or. size(coef_mean) /= p) return
        if (size(coef_cov, 1) /= p .or. size(coef_cov, 2) /= p) return
        if (.not. all(ieee_is_finite(design))) return
        if (.not. all(ieee_is_finite(coef_mean))) return
        if (present(observation_variance)) then
            if (.not. present(noise_variance)) return
        end if
        if (present(noise_variance)) then
            if (size(noise_variance) /= n) return
            if (.not. all(ieee_is_finite(noise_variance))) return
            if (any(noise_variance < 0.0_dp)) return
        end if
        call validate_covariance(coef_cov, status)
        if (.not. status_ok(status)) return
        mean = matmul(design, coef_mean)
        work = matmul(design, coef_cov)
        allocate(latent_variance(n))
        !$omp parallel do default(none) shared(n, design, work, latent_variance) private(i)
        do i = 1, n
            latent_variance(i) = dot_product(design(i, :), work(i, :))
        end do
        !$omp end parallel do
        if (.not. all(ieee_is_finite(mean)) .or. &
            .not. all(ieee_is_finite(latent_variance))) then
            deallocate(mean, latent_variance)
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'diagonal prediction overflow')
            return
        end if
        if (present(observation_variance)) then
            observation_variance = latent_variance + noise_variance
            if (.not. all(ieee_is_finite(observation_variance))) then
                deallocate(mean, latent_variance, observation_variance)
                call status_set(status, FORTNUM_DOMAIN_ERROR, 'diagonal prediction overflow')
                return
            end if
        end if
        call status_set(status, FORTNUM_OK, '')
    end subroutine

    subroutine bayesian_linear_predict(design, coef_mean, coef_cov, mean, &
            latent_cov, status, noise_cov, observation_cov)
        real(dp), intent(in) :: design(:, :), coef_mean(:), coef_cov(:, :)
        real(dp), allocatable, intent(out) :: mean(:), latent_cov(:, :)
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(in), optional :: noise_cov(:, :)
        real(dp), allocatable, intent(out), optional :: observation_cov(:, :)
        if (present(observation_cov)) then
            if (.not. present(noise_cov)) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    'observation covariance requires explicit noise')
                return
            end if
        end if
        call linear_propagate(coef_mean, coef_cov, design, mean, latent_cov, status)
        if (.not. status_ok(status)) return
        if (present(noise_cov)) then
            if ((size(noise_cov, 1) /= size(latent_cov, 1) .or. &
                size(noise_cov, 2) /= size(latent_cov, 2))) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, 'prediction noise shape')
                deallocate (mean, latent_cov)
                return
            end if
            call validate_covariance(noise_cov, status)
            if (.not. status_ok(status)) then
                deallocate (mean, latent_cov)
                return
            end if
            if (present(observation_cov)) then
                observation_cov = latent_cov + noise_cov
                if (.not. all(ieee_is_finite(observation_cov))) then
                    deallocate(mean, latent_cov, observation_cov)
                    call status_set(status, FORTNUM_DOMAIN_ERROR, 'prediction overflow')
                    return
                end if
            end if
        end if
    end subroutine
end module fortprob_regression
