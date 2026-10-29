module fortprob_regression
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, status_set, status_ok, &
        FORTNUM_DOMAIN_ERROR
    use fortprob_gaussian, only: gaussian_condition, validate_covariance, &
        linear_propagate
    implicit none
    private
    public :: bayesian_linear_fit, bayesian_linear_predict
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
        if (any(shape(prior_cov) /= [p, p])) return
        if (any(shape(noise_cov) /= [n, n])) return
        if (.not. all(ieee_is_finite(design))) return
        call validate_covariance(prior_cov, status, .true.)
        if (.not. status_ok(status)) return
        call validate_covariance(noise_cov, status, .true.)
        if (.not. status_ok(status)) return
        cross = matmul(prior_cov, transpose(design))
        obs_cov = matmul(design, cross) + noise_cov
        obs_mean = matmul(design, prior_mean)
        call gaussian_condition(prior_mean, prior_cov, obs_mean, obs_cov, cross, y, &
            coef_mean, coef_cov, status)
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
            if (any(shape(noise_cov) /= shape(latent_cov))) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, 'prediction noise shape')
                deallocate (mean, latent_cov)
                return
            end if
            call validate_covariance(noise_cov, status)
            if (.not. status_ok(status)) then
                deallocate (mean, latent_cov)
                return
            end if
            if (present(observation_cov)) observation_cov = latent_cov + noise_cov
        end if
    end subroutine
end module fortprob_regression
