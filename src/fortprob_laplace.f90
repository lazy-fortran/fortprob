module fortprob_laplace
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, status_set, status_ok, &
        FORTNUM_DOMAIN_ERROR, FORTNUM_OK
    use fortnum_cholesky, only: cholesky_factorization_t
    use fortprob_gaussian, only: validate_covariance
    implicit none
    private
    public :: laplace_at_mode
    interface
        pure subroutine fortprob_laplace_leaf(n, log_joint_mode, logdet, log_evidence)
            import dp
            real(dp), intent(in) :: n, log_joint_mode, logdet
            real(dp), intent(out) :: log_evidence
        end subroutine
    end interface
contains
    subroutine laplace_at_mode(log_joint_mode, hessian, covariance, &
            log_evidence, status)
        real(dp), intent(in) :: log_joint_mode, hessian(:, :)
        real(dp), allocatable, intent(out) :: covariance(:, :)
        real(dp), intent(out) :: log_evidence
        type(fortnum_status_t), intent(out) :: status
        type(cholesky_factorization_t) :: factor
        real(dp), allocatable :: inverse(:, :)
        real(dp) :: logdet
        integer :: n, i
        log_evidence = 0.0_dp
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'invalid Laplace log density')
        if (.not. ieee_is_finite(log_joint_mode)) return
        call validate_covariance(hessian, status, .true.)
        if (.not. status_ok(status)) return
        call factor%factorize(hessian, status)
        if (.not. status_ok(status)) return
        call factor%log_determinant(logdet, status)
        if (.not. status_ok(status)) return
        n = size(hessian, 1)
        allocate(inverse(n, n))
        inverse = 0.0_dp
        do i = 1, n
            inverse(i, i) = 1.0_dp
        end do
        call factor%solve(inverse, status)
        if (.not. status_ok(status)) return
        call fortprob_laplace_leaf(real(n, dp), log_joint_mode, logdet, log_evidence)
        if (.not. all(ieee_is_finite(inverse)) .or. &
            .not. ieee_is_finite(log_evidence)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'Laplace overflow')
            return
        end if
        call move_alloc(inverse, covariance)
        call status_set(status, FORTNUM_OK, '')
    end subroutine
end module
