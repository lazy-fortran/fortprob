program generate
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
        operator(*), operator(/), operator(**), log
    use fortsym_diff, only: diff
    use fortsym_kernel_emit, only: kernel_emit_spec_t, emit_fortran_kernel_ir
    use fortsym_kernel_ir, only: kernel_ir_t, lower_kernel_ir
    use fortsym_string, only: chars, str, str_t
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: variance, mean, z, root(3), n, logdet, quadratic, density(1)
    type(expr_t) :: count_ge, alpha, empirical(1), log_joint_mode
    type(expr_t) :: a, denominator, derivative, k, ratio, leaves(2)
    type(kernel_emit_spec_t) :: spec
    type(kernel_ir_t) :: ir
    type(str_t) :: source
    character(:), allocatable :: message, leaf_source
    logical :: good
    integer :: unit
    call arena%init()
    variance = sym(arena, 'variance')
    mean = sym(arena, 'mean')
    z = sym(arena, 'z')
    root = [mean - z*variance**0.5_dp, mean + z*variance**0.5_dp, &
        diff(variance**0.5_dp, variance)]
    call lower_kernel_ir(root, ir, good, message)
    if (.not. good) error stop message
    spec%name = str('fortprob_interval_leaf')
    spec%args = [str('mean'), str('variance'), str('z')]
    spec%outputs = [str('lower'), str('upper'), str('dsigma_dvariance')]
    spec%pure_procedure = .true.
    spec%generator = str('gen/generate.f90')
    spec%generator_revision = str('fortsym 063008070086aaf2d78f6e9d9a8fddff44c74a46')
    spec%regenerate_command = str('gen/regenerate.sh')
    source = emit_fortran_kernel_ir(ir, spec, good, message)
    if (.not. good) error stop message
    open (newunit=unit, file='src/generated/fortprob_interval_leaf.f90', &
        status='replace')
    write (unit, '(a)') chars(source)
    close (unit)
    n = sym(arena, 'n')
    logdet = sym(arena, 'logdet')
    quadratic = sym(arena, 'quadratic')
    density = [-0.5_dp*(n*log(2.0_dp*acos(-1.0_dp)) + logdet + quadratic)]
    call lower_kernel_ir(density, ir, good, message)
    if (.not. good) error stop message
    spec%name = str('fortprob_logpdf_leaf')
    spec%args = [str('n'), str('logdet'), str('quadratic')]
    spec%outputs = [str('density')]
    source = emit_fortran_kernel_ir(ir, spec, good, message)
    if (.not. good) error stop message
    open (newunit=unit, file='src/generated/fortprob_logpdf_leaf.f90', &
        status='replace')
    write (unit, '(a)') chars(source)
    close (unit)
    log_joint_mode = sym(arena, 'log_joint_mode')
    density = [log_joint_mode + 0.5_dp*n*log(2.0_dp*acos(-1.0_dp)) &
        - 0.5_dp*logdet]
    call lower_kernel_ir(density, ir, good, message)
    if (.not. good) error stop message
    spec%policy%shape_fma = .false.
    spec%name = str('fortprob_laplace_leaf')
    spec%args = [str('n'), str('log_joint_mode'), str('logdet')]
    spec%outputs = [str('log_evidence')]
    source = emit_fortran_kernel_ir(ir, spec, good, message)
    if (.not. good) error stop message
    open (newunit=unit, file='src/generated/fortprob_laplace_leaf.f90', &
        status='replace')
    leaf_source = chars(source)
    do while (len(leaf_source) > 0)
        if (leaf_source(len(leaf_source):) /= new_line('a')) exit
        leaf_source = leaf_source(:len(leaf_source) - 1)
    end do
    write (unit, '(a)') leaf_source
    close (unit)
    spec%policy%shape_fma = .true.
    count_ge = sym(arena, 'count_ge')
    alpha = sym(arena, 'alpha')
    empirical = [(1.0_dp + count_ge)/(n + 1.0_dp)]
    call lower_kernel_ir(empirical, ir, good, message)
    if (.not. good) error stop message
    spec%name = str('fortprob_empirical_tail_leaf')
    spec%args = [str('n'), str('count_ge')]
    spec%outputs = [str('probability')]
    source = emit_fortran_kernel_ir(ir, spec, good, message)
    if (.not. good) error stop message
    open (newunit=unit, file='src/generated/fortprob_empirical_tail_leaf.f90', &
        status='replace')
    write (unit, '(a)') chars(source)
    close (unit)
    empirical = [(n + 1.0_dp)*(1.0_dp - alpha)]
    call lower_kernel_ir(empirical, ir, good, message)
    if (.not. good) error stop message
    spec%name = str('fortprob_empirical_rank_leaf')
    spec%args = [str('n'), str('alpha')]
    spec%outputs = [str('location')]
    source = emit_fortran_kernel_ir(ir, spec, good, message)
    if (.not. good) error stop message
    open (newunit=unit, file='src/generated/fortprob_empirical_rank_leaf.f90', &
        status='replace')
    write (unit, '(a)') chars(source)
    close (unit)
    a = sym(arena, 'a')
    denominator = sym(arena, 'denominator')
    derivative = sym(arena, 'derivative')
    k = sym(arena, 'k')
    leaves = [a + k/denominator, 1.0_dp - k*derivative/denominator**2]
    call lower_kernel_ir(leaves, ir, good, message)
    if (.not. good) error stop message
    spec%name = str('fortprob_mills_step_leaf')
    spec%args = [str('a'), str('denominator'), str('derivative'), str('k')]
    spec%outputs = [str('next_denominator'), str('next_derivative')]
    source = emit_fortran_kernel_ir(ir, spec, good, message)
    if (.not. good) error stop message
    open (newunit=unit, file='src/generated/fortprob_mills_step_leaf.f90', &
        status='replace')
    write (unit, '(a)') chars(source)
    close (unit)
    ratio = sym(arena, 'ratio')
    leaves = [ratio, 1.0_dp + a*ratio - ratio**2]
    call lower_kernel_ir(leaves, ir, good, message)
    if (.not. good) error stop message
    spec%name = str('fortprob_truncated_leaf')
    spec%args = [str('a'), str('ratio')]
    spec%outputs = [str('shift'), str('variance_factor')]
    source = emit_fortran_kernel_ir(ir, spec, good, message)
    if (.not. good) error stop message
    open (newunit=unit, file='src/generated/fortprob_truncated_leaf.f90', &
        status='replace')
    write (unit, '(a)') chars(source)
    close (unit)
end program generate
