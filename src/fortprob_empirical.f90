module fortprob_empirical
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
        ieee_quiet_nan, ieee_positive_inf
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_DOMAIN_ERROR, FORTNUM_OK
    use fortnum_multiroot, only: argsort
    implicit none
    private
    public :: empirical_calibration_t

    interface
        pure subroutine fortprob_empirical_tail_leaf(n, count_ge, probability)
            import dp
            real(dp), intent(in) :: n, count_ge
            real(dp), intent(out) :: probability
        end subroutine
        pure subroutine fortprob_empirical_rank_leaf(n, alpha, location)
            import dp
            real(dp), intent(in) :: n, alpha
            real(dp), intent(out) :: location
        end subroutine
    end interface

    ! Finite calibration scores sorted once; queries take logarithmic time.
    ! Ordering is discrete and has no derivative at rank boundaries.
    type :: empirical_calibration_t
        private
        real(dp), allocatable :: scores(:)
    contains
        procedure :: fit
        procedure :: tail_probability
        procedure :: upper_quantile
        procedure :: count => sample_count
        procedure :: clear
    end type
contains
    subroutine clear(self)
        class(empirical_calibration_t), intent(inout) :: self
        if (allocated(self%scores)) deallocate(self%scores)
    end subroutine

    integer function sample_count(self) result(n)
        class(empirical_calibration_t), intent(in) :: self
        n = 0
        if (allocated(self%scores)) n = size(self%scores)
    end function

    ! A rejected fit clears the previous calibration.
    subroutine fit(self, scores, status)
        class(empirical_calibration_t), intent(inout) :: self
        real(dp), intent(in) :: scores(:)
        type(fortnum_status_t), intent(out) :: status
        integer :: order(size(scores))
        call self%clear()
        call status_set(status, FORTNUM_DOMAIN_ERROR, &
            'finite nonempty calibration required')
        if (size(scores) == 0) return
        if (.not. all(ieee_is_finite(scores))) return
        call argsort(scores, order)
        self%scores = scores(order)
        call status_set(status, FORTNUM_OK, '')
    end subroutine

    ! Tie-inclusive tail: (1 + number of scores >= query)/(n + 1).
    function tail_probability(self, query, status) result(probability)
        class(empirical_calibration_t), intent(in) :: self
        real(dp), intent(in) :: query
        type(fortnum_status_t), intent(out) :: status
        real(dp) :: probability
        integer :: n
        integer(int64) :: lo, hi, middle
        probability = ieee_value(0.0_dp, ieee_quiet_nan)
        call status_set(status, FORTNUM_DOMAIN_ERROR, &
            'finite query and fitted calibration required')
        n = self%count()
        if (n == 0) return
        if (.not. ieee_is_finite(query)) return
        lo = 1_int64
        hi = int(n, int64)
        do while (lo <= hi)
            middle = lo + (hi - lo)/2_int64
            if (self%scores(middle) < query) then
                lo = middle + 1_int64
            else
                hi = middle - 1_int64
            end if
        end do
        call fortprob_empirical_tail_leaf(real(n, dp), &
            real(int(n, int64) - hi, dp), probability)
        call status_set(status, FORTNUM_OK, '')
    end function

    ! Rank ceil((n+1)*(1-alpha)); ranks above n return positive infinity.
    ! This empirical rule alone supplies no temporal coverage guarantee.
    function upper_quantile(self, alpha, status) result(threshold)
        class(empirical_calibration_t), intent(in) :: self
        real(dp), intent(in) :: alpha
        type(fortnum_status_t), intent(out) :: status
        real(dp) :: threshold, location
        integer :: rank, n
        threshold = ieee_value(0.0_dp, ieee_quiet_nan)
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'alpha must be finite and in (0,1)')
        n = self%count()
        if (n == 0) return
        if (.not. ieee_is_finite(alpha)) return
        if (alpha <= 0.0_dp .or. alpha >= 1.0_dp) return
        call fortprob_empirical_rank_leaf(real(n, dp), alpha, location)
        if (location > real(n, dp)) then
            threshold = ieee_value(0.0_dp, ieee_positive_inf)
        else
            rank = max(1, ceiling(location))
            threshold = self%scores(rank)
        end if
        call status_set(status, FORTNUM_OK, '')
    end function
end module fortprob_empirical
