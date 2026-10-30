module fortprob_mixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR
    implicit none
    private
    public :: probability_mean
contains
    pure subroutine probability_mean(probabilities, weights, value, status)
        real(dp), intent(in) :: probabilities(:), weights(:)
        real(dp), intent(out) :: value
        type(fortnum_status_t), intent(out) :: status
        real(dp) :: weight_max, weight, mass, mean
        integer :: i, n
        value = ieee_value(0.0_dp, ieee_quiet_nan)
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'probability mean shape/data')
        n = size(probabilities)
        if (n < 1 .or. size(weights) /= n) return
        if (.not. all(ieee_is_finite(probabilities))) return
        if (any(probabilities < 0.0_dp) .or. any(probabilities > 1.0_dp)) return
        if (.not. all(ieee_is_finite(weights))) return
        if (any(weights < 0.0_dp)) return
        weight_max = maxval(weights)
        if (weight_max <= 0.0_dp) return
        mean = 0.0_dp
        mass = 0.0_dp
        do i = 1, n
            weight = weights(i)/weight_max
            if (weight == 0.0_dp) cycle
            mass = mass + weight
            ! Convex updates preserve the range despite rounded normalization.
            mean = mean + (probabilities(i) - mean)*(weight/mass)
        end do
        value = mean
        call status_set(status, FORTNUM_OK, '')
    end subroutine
end module fortprob_mixture
