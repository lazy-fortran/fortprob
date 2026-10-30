module fortprob_affine_filter
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, status_ok, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    use fortnum_psd_factor, only: psd_factor
    use fortnum_gram_factor, only: gram_factor_workspace_t
    implicit none
    private
    type, public :: gaussian_interval_t
        private
        real(dp), allocatable :: f(:, :), a(:), process_root(:, :), qcross(:)
        real(dp) :: qavg = 0.0_dp, measurement_variance = 0.0_dp
        logical :: endpoint_noise = .false.
    contains
        procedure :: init => interval_init
    end type
    type, public :: affine_gaussian_filter_t
        private
        real(dp), allocatable :: c(:), d(:, :), root(:, :)
        type(gram_factor_workspace_t) :: workspace
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
        real(dp), allocatable :: joint(:, :), candidate(:, :)
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
        allocate(joint(n + 1, n + 1), candidate(n + 1, n + 1))
        joint(:n, :n) = qend
        joint(:n, n + 1) = qcross
        joint(n + 1, :n) = qcross
        joint(n + 1, n + 1) = qavg
        call psd_factor(joint, candidate, status)
        if (.not. status_ok(status)) return
        self%f = f
        self%a = a
        call move_alloc(candidate, self%process_root)
        self%qavg = qavg
        self%qcross = qcross
        self%measurement_variance = measurement_variance
        self%endpoint_noise = any(qend /= 0.0_dp)
    end subroutine

    subroutine filter_init(self, p_initial, ncoeff, status)
        class(affine_gaussian_filter_t), intent(inout) :: self
        real(dp), intent(in) :: p_initial(:, :)
        integer, intent(in) :: ncoeff
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: candidate(:, :)
        integer :: n
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'filter coefficient count')
        if (ncoeff < 1) return
        n = size(p_initial, 1)
        allocate(candidate(n, n))
        call psd_factor(p_initial, candidate, status)
        if (.not. status_ok(status)) return
        call move_alloc(candidate, self%root)
        if (allocated(self%c)) deallocate(self%c)
        if (allocated(self%d)) deallocate(self%d)
        allocate(self%c(n), self%d(n, ncoeff))
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
        real(dp), allocatable :: averaged_root(:), cross(:), gain(:), c0(:), d0(:, :)
        real(dp), allocatable :: propagated(:, :), candidate(:, :), columns(:, :)
        real(dp) :: state_variance
        integer :: n, j
        offset = 0.0_dp
        effective_design = 0.0_dp
        innovation_variance = 0.0_dp
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'filter initialization required')
        if (.not. allocated(self%root)) return
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
        averaged_root = matmul(interval%a, self%root)
        offset = dot_product(interval%a, self%c)
        effective_design = matmul(interval%a, self%d)
        effective_design = design - effective_design
        ! A Gram representation preserves nonnegative variance through underflow.
        state_variance = dot_product(averaged_root, averaged_root)
        innovation_variance = state_variance + interval%qavg &
            + interval%measurement_variance
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
        c0 = matmul(interval%f, self%c)
        d0 = matmul(interval%f, self%d)
        propagated = matmul(interval%f, self%root)
        if (condition .or. interval%endpoint_noise) then
            allocate(columns(n, 2*n + 2), candidate(n, n))
            columns = 0.0_dp
            columns(:, :n) = propagated
            columns(:, n + 1:2*n + 1) = interval%process_root(:n, :)
            if (condition) then
                cross = matmul(propagated, averaged_root)
                cross = cross + interval%qcross
                gain = cross/innovation_variance
                c0 = c0 + gain*(y_minus_base - offset)
                do j = 1, size(design)
                    d0(:, j) = d0(:, j) + gain*effective_design(j)
                end do
                ! Joseph columns retain endpoint/average process correlation.
                do j = 1, n
                    columns(:, j) = columns(:, j) - gain*averaged_root(j)
                end do
                do j = 1, n + 1
                    columns(:, n + j) = columns(:, n + j) &
                        - gain*interval%process_root(n + 1, j)
                end do
                columns(:, 2*n + 2) = -gain*sqrt(interval%measurement_variance)
            end if
            call self%workspace%factor(columns, candidate, status)
            if (.not. status_ok(status)) return
        else
            call move_alloc(propagated, candidate)
        end if
        if (.not. all(ieee_is_finite(c0))) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'filter propagated offset nonfinite')
            return
        end if
        if (.not. all(ieee_is_finite(d0))) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'filter propagated design nonfinite')
            return
        end if
        if (.not. all(ieee_is_finite(candidate))) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'filter propagated covariance root nonfinite')
            return
        end if
        ! Finite covariance diagonals remain required even though storage is factored.
        do j = 1, n
            if (.not. ieee_is_finite(dot_product(candidate(j, :), candidate(j, :)))) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, 'filter propagated covariance overflow')
                return
            end if
        end do
        call move_alloc(c0, self%c)
        call move_alloc(d0, self%d)
        call move_alloc(candidate, self%root)
        call status_set(status, FORTNUM_OK, '')
    end subroutine
end module fortprob_affine_filter
