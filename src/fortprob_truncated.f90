module fortprob_truncated
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
        ieee_quiet_nan, ieee_negative_inf
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR
    implicit none
    private
    public :: normal_lower_truncated_moments
    interface
        pure subroutine fortprob_mills_step_leaf(a, denominator, derivative, k, nd, nv)
            import dp
            real(dp), intent(in) :: a, denominator, derivative, k
            real(dp), intent(out) :: nd, nv
        end subroutine
        pure subroutine fortprob_truncated_leaf(a, ratio, shift, variance_factor)
            import dp
            real(dp), intent(in) :: a, ratio
            real(dp), intent(out) :: shift, variance_factor
        end subroutine
    end interface
contains
    ! Moments and log P(X >= lower) for a scalar Gaussian, including a point mass.
    subroutine normal_lower_truncated_moments(mu, variance, lower, mean, &
            truncated_variance, log_mass, status)
        real(dp), intent(in) :: mu, variance, lower
        real(dp), intent(out) :: mean, truncated_variance, log_mass
        type(fortnum_status_t), intent(out) :: status
        real(dp) :: sigma, a, ratio, shift, factor, denominator, derivative, nd, nv
        real(dp), parameter :: log_sqrt_two_pi = 0.91893853320467274178_dp
        integer :: k
        mean = ieee_value(0.0_dp, ieee_quiet_nan)
        truncated_variance = mean
        log_mass = mean
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid truncated normal input')
        if (.not. ieee_is_finite(mu)) return
        if (.not. ieee_is_finite(variance)) return
        if (.not. ieee_is_finite(lower)) return
        if (variance < 0.0_dp) return
        if (variance == 0.0_dp) then
            if (mu < lower) then
                log_mass = ieee_value(0.0_dp, ieee_negative_inf)
                call status_set(status, FORTNUM_DOMAIN_ERROR, 'conditioning event has zero mass')
                return
            end if
            mean = mu
            truncated_variance = 0.0_dp
            log_mass = 0.0_dp
            call status_set(status, FORTNUM_OK, '')
            return
        end if
        sigma = sqrt(variance)
        a = (lower - mu)/sigma
        if (.not. ieee_is_finite(a)) return
        if (a > sqrt(huge(a))) return
        if (a > 8.0_dp) then
            ! Laplace's continued fraction also gives the derivative without
            ! subtracting the nearly equal terms in the usual variance formula.
            denominator = a
            derivative = 1.0_dp
            do k = 128, 2, -1
                call fortprob_mills_step_leaf(a, denominator, derivative, &
                    real(k, dp), nd, nv)
                denominator = nd
                derivative = nv
            end do
            shift = 1.0_dp/denominator
            ratio = a + shift
            factor = (derivative/denominator)/denominator
            mean = lower + sigma*shift
            log_mass = -0.5_dp*a*a - log_sqrt_two_pi - log(ratio)
        else
            log_mass = log(0.5_dp*erfc(a/sqrt(2.0_dp)))
            ratio = exp(-0.5_dp*a*a - log_sqrt_two_pi - log_mass)
            call fortprob_truncated_leaf(a, ratio, shift, factor)
            mean = mu + sigma*shift
        end if
        truncated_variance = variance*factor
        if (.not. ieee_is_finite(mean)) return
        if (.not. ieee_is_finite(truncated_variance)) return
        if (.not. ieee_is_finite(log_mass)) return
        call status_set(status, FORTNUM_OK, '')
    end subroutine
end module
