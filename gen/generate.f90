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
    type(kernel_emit_spec_t) :: spec
    type(kernel_ir_t) :: ir
    type(str_t) :: source
    character(:), allocatable :: message
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
end program generate
