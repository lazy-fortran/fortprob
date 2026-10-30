module fortprob_affine_filter
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, status_ok, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    use fortprob_gaussian, only: validate_covariance
    implicit none
    private
    type, public :: gaussian_interval_t
        private
        real(dp), allocatable :: f(:, :), a(:), qend(:, :), qcross(:)
        real(dp) :: qavg = 0.0_dp, measurement_variance = 0.0_dp
    contains
        procedure :: init => interval_init
    end type
    type, public :: affine_gaussian_filter_t
        private
        real(dp), allocatable :: c(:), d(:, :), p(:, :)
    contains
        procedure :: init => filter_init
        procedure :: step => filter_step
    end type
contains
    subroutine interval_init(self, f, a, qend, qavg, qcross, &
            measurement_variance, status)
        class(gaussian_interval_t), intent(inout) :: self
        real(dp), intent(in) :: f(:, :), a(:), qend(:, :), qavg, qcross(:)
        real(dp), intent(in) :: measurement_variance
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: joint(:, :)
        integer :: n
        n = size(a)
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'interval shape/data')
        if (n < 1) return
        if (size(f, 1) /= n .or. size(f, 2) /= n) return
        if (size(qend, 1) /= n .or. size(qend, 2) /= n) return
        if (size(qcross) /= n) return
        if (.not. all(ieee_is_finite(f))) return
        if (.not. all(ieee_is_finite(a))) return
        if (.not. ieee_is_finite(measurement_variance)) return
        if (measurement_variance <= 0.0_dp) return
        allocate(joint(n + 1, n + 1))
        joint(:n, :n) = qend
        joint(:n, n + 1) = qcross
        joint(n + 1, :n) = qcross
        joint(n + 1, n + 1) = qavg
        call validate_covariance(joint, status)
        if (.not. status_ok(status)) return
        self%f = f
        self%a = a
        self%qend = qend
        self%qavg = qavg
        self%qcross = qcross
        self%measurement_variance = measurement_variance
    end subroutine

    subroutine filter_init(self, p_initial, ncoeff, status)
        class(affine_gaussian_filter_t), intent(inout) :: self
        real(dp), intent(in) :: p_initial(:, :)
        integer, intent(in) :: ncoeff
        type(fortnum_status_t), intent(out) :: status
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'filter coefficient count')
        if (ncoeff < 1) return
        call validate_covariance(p_initial, status)
        if (.not. status_ok(status)) return
        self%p = p_initial
        if (allocated(self%c)) deallocate(self%c)
        if (allocated(self%d)) deallocate(self%d)
        allocate(self%c(size(p_initial, 1)), self%d(size(p_initial, 1), ncoeff))
        self%c = 0.0_dp
        self%d = 0.0_dp
    end subroutine

    subroutine filter_step(self, interval, design, y_minus_base, condition, &
            offset, effective_design, innovation_variance, status)
        class(affine_gaussian_filter_t), intent(inout) :: self
        type(gaussian_interval_t), intent(in) :: interval
        real(dp), intent(in) :: design(:), y_minus_base
        logical, intent(in) :: condition
        real(dp), intent(out) :: offset, effective_design(:), innovation_variance
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: pa(:), cross(:), c0(:), d0(:, :), p0(:, :)
        real(dp), allocatable :: work(:, :), gain(:), transition(:, :)
        real(dp) :: scale, tolerance, state_variance
        integer :: n, i, j
        offset = 0.0_dp
        effective_design = 0.0_dp
        innovation_variance = 0.0_dp
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'filter initialization required')
        if (.not. allocated(self%p)) return
        if (.not. allocated(interval%f)) return
        n = size(self%c)
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'filter step shape/data')
        if (size(interval%a) /= n) return
        if (size(design) /= size(self%d, 2)) return
        if (size(effective_design) /= size(design)) return
        if (.not. all(ieee_is_finite(design))) return
        if (condition) then
            if (.not. ieee_is_finite(y_minus_base)) return
        end if
        pa = matmul(self%p, interval%a)
        offset = dot_product(interval%a, self%c)
        effective_design = matmul(interval%a, self%d)
        effective_design = design - effective_design
        state_variance = dot_product(interval%a, pa)
        innovation_variance = state_variance + interval%qavg &
            + interval%measurement_variance
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'filter step overflow')
        if (.not. ieee_is_finite(offset)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'filter offset nonfinite')
            return
        end if
        if (.not. all(ieee_is_finite(effective_design))) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'filter effective design nonfinite')
            return
        end if
        if (.not. ieee_is_finite(innovation_variance)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'filter innovation variance nonfinite')
            return
        end if
        if (innovation_variance <= 0.0_dp) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'filter innovation variance nonpositive')
            return
        end if
        tolerance = 128.0_dp*epsilon(1.0_dp)*real(n, dp) &
            *maxval(abs(self%p))*sum(abs(interval%a))**2
        if (state_variance < -tolerance) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'filter projected state variance negative')
            return
        end if
        cross = matmul(interval%f, pa)
        cross = cross + interval%qcross
        c0 = matmul(interval%f, self%c)
        d0 = matmul(interval%f, self%d)
        work = matmul(interval%f, self%p)
        allocate(p0(n, n))
        ! Assemble each covariance entry once in its symmetric representation.
        do j = 1, n
            do i = 1, j
                p0(i, j) = dot_product(work(i, :), interval%f(j, :)) &
                    + interval%qend(i, j)
                p0(j, i) = p0(i, j)
            end do
        end do
        scale = maxval(abs(p0))
        if (condition) then
            c0 = c0 + cross*((y_minus_base - offset)/innovation_variance)
            do j = 1, size(design)
                d0(:, j) = d0(:, j) + cross*(effective_design(j)/innovation_variance)
            end do
            ! Joseph form avoids subtracting nearly equal state variances.
            gain = cross/innovation_variance
            transition = interval%f
            do j = 1, n
                transition(:, j) = transition(:, j) - gain*interval%a(j)
            end do
            work = matmul(transition, self%p)
            do j = 1, n
                do i = 1, j
                    p0(i, j) = dot_product(work(i, :), transition(j, :)) &
                        + interval%qend(i, j) &
                        - (gain(i)*interval%qcross(j) + interval%qcross(i)*gain(j)) &
                        + (gain(i)*gain(j))*(interval%qavg + interval%measurement_variance)
                    p0(j, i) = p0(i, j)
                end do
            end do
        end if
        if (.not. all(ieee_is_finite(c0))) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'filter propagated offset nonfinite')
            return
        end if
        if (.not. all(ieee_is_finite(d0))) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'filter propagated design nonfinite')
            return
        end if
        if (.not. all(ieee_is_finite(p0))) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'filter propagated covariance nonfinite')
            return
        end if
        ! Valid initial and joint process covariances preserve PSD analytically.
        tolerance = 128.0_dp*epsilon(1.0_dp)*real(n, dp)*scale
        if (maxval(abs(p0 - transpose(p0))) > tolerance) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'filter propagated covariance nonsymmetric')
            return
        end if
        do j = 1, n
            if (p0(j, j) < -tolerance) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, 'filter propagated covariance diagonal negative')
                return
            end if
        end do
        call move_alloc(c0, self%c)
        call move_alloc(d0, self%d)
        call move_alloc(p0, self%p)
        call status_set(status, FORTNUM_OK, '')
    end subroutine
end module fortprob_affine_filter
