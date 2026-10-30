# fortprob

Initial implementation of [fortprob](https://github.com/lazy-fortran/fortprob),
providing dense Gaussian
conditioning and Bayesian linear regression with a proper Gaussian coefficient
prior. Models can use polynomial, spline, or other caller-supplied design
matrices. Every covariance retains its cross-correlations.

The runtime is native Fortran with real64 arrays. Numerical factorizations,
solves and symmetric eigenvalue checks come from
[fortnum](https://github.com/lazy-fortran/fortnum), pinned in `fpm.toml`.
There is no symbolic engine dependency at runtime.

```fortran
use fortprob_regression, only: bayesian_linear_fit, bayesian_linear_predict
use fortnum_status, only: fortnum_status_t, status_ok
! design(n_observations,n_coefficients), noise_cov(n_observations,n_observations)
call bayesian_linear_fit(design, y, noise_cov, prior_mean, prior_cov, &
    coefficient_mean, coefficient_cov, status)
if (.not. status_ok(status)) error stop status%msg
call bayesian_linear_predict(query_design, coefficient_mean, coefficient_cov, &
    mean, latent_cov, status, noise_cov=future_noise, &
    observation_cov=observation_cov)
```

For independent observation errors, use
`bayesian_linear_fit_diagonal(design,y,noise_variance,prior_mean,prior_cov,coef_mean,coef_cov,status)`.
It accepts a strictly positive variance per observation and solves in coefficient
space, storing only the design and coefficient-size matrices: `O(n*p+p*p)`
storage for `n` observations and `p` coefficients, without an `n` by `n` matrix.
The optional final argument `log_evidence` returns the normalized marginal log
density after integrating the coefficients against the supplied proper Gaussian
prior. Noise variances and prior parameters are fixed inputs; this does not
integrate their uncertainty or justify priors estimated from the same targets.
`bayesian_linear_predict_diagonal(design,coef_mean,coef_cov,mean,latent_variance,status,noise_variance,observation_variance)`
computes marginal prediction variances with `O(n*p+p*p)` storage for `n` query
rows and `p` coefficients. Optional
future `noise_variance` is a vector; an `observation_variance` result requires it.
It preserves the full coefficient covariance while returning only output
marginals; use the joint predictor when output cross-covariances are needed.

`prior_cov` and training `noise_cov` must be positive definite. Prediction
noise may be positive semidefinite and is independent of the latent prediction.
`latent_cov` describes the noiseless model. Requesting `observation_cov` requires
explicit `noise_cov`; the result adds that future noise to `latent_cov`. Outputs are allocatable and allocated on success. Failed
Gaussian operations leave array outputs unallocated; inspect `status%code` and
`status%msg` before accessing them.

`fortprob_gaussian` exports:

- `gaussian_condition(mean_x,cov_xx,mean_y,cov_yy,cov_xy,y,post_mean,post_cov,status)`.
  The cross-covariance has shape `(size(mean_x),size(mean_y))`. `cov_yy` must
  be positive definite; the joint distribution may be positive semidefinite.
- Use `gaussian_conditioner_t` for repeated queries with fixed observations.
  Call `fit(mean_y,cov_yy,y,status)` to store the covariance factor and solved
  observation residual, then call
  `condition(mean_x,cov_xx,cov_xy,post_mean,post_cov,status)` for each query.
  The object owns its fitted data, so changing the original arrays leaves the
  fit unchanged. A failed refit clears the previous fit. Covariance checks
  follow the same rules as `gaussian_condition`.
- `gaussian_condition_jvp` takes the same six inputs, followed by their six
  directional derivatives, then `dpost_mean,dpost_cov,status`. This analytical
  implicit derivative includes both covariance and observation perturbations.
  Tangent covariances should be symmetric; structures and dimensions stay fixed.
- `gaussian_logpdf(value,mean,cov,log_density,status)` requires positive
  definite covariance and returns the normalized joint log density.
- `linear_propagate(mean,cov,jacobian,out_mean,out_cov,status)` evaluates a
  linear map and its full covariance. For a nonlinear map, use its nominal
  output separately and pass its Jacobian for first-order covariance propagation.
- `gaussian_intervals(mean,variance,z,lower,upper,status)` returns marginal
  intervals using an explicit nonnegative normal multiplier `z`; for example
  `z=1.959963984540054` gives approximately 95% marginal Gaussian coverage.
  These are pointwise intervals, not simultaneous joint coverage regions.
- `validate_covariance(cov,status,positive_definite,reference_scale)` checks
  finite values, symmetry and definiteness. Its optional logical flag defaults
  to false. The optional nonnegative reference scale is for subtraction error.

`fortprob_truncated` exports
`normal_lower_truncated_moments(mu,variance,lower,mean,truncated_variance,log_mass,status)`.
It returns the exact scalar Gaussian conditional moments for `X >= lower`, and
`log_mass = log(P(X >= lower))`, with real64 scalar arguments and
`fortnum_status_t`. Far-tail evaluation avoids subtractive cancellation in the
variance and does not exponentiate the conditioning mass. Zero variance denotes
a point mass: conditioning succeeds when `mu >= lower`; otherwise it returns a
domain error, NaN moments and negative infinite log mass. Nonfinite inputs,
negative variance or unrepresentable standardized bounds/results return a domain
error. Check status before using outputs.

`fortprob_laplace` exports
`laplace_at_mode(log_joint_mode,hessian,covariance,log_evidence,status)`.
All numeric inputs and outputs use real64; `covariance` is allocatable and
`status` is `fortnum_status_t`. The caller supplies a stationary interior mode,
the log joint density there, and the full Hessian of the negative log joint
density in the same parameter coordinates. The result is the inverse Hessian
and the local Laplace approximation to the log integral. No optimizer or
stationarity check is provided. A normalized posterior gives a log integral
approximately zero; model evidence requires the joint density with every proper
prior and likelihood normalizer, including coordinate-transformation Jacobians.
The caller is responsible for these normalizers. Improper priors do not define
model evidence.

The Hessian must be finite, symmetric to the existing covariance roundoff
tolerance, and strictly positive definite. Failed calls leave `covariance`
unallocated; all outputs require a successful status. There is no jitter or
Gauss–Newton substitution. This is a local unimodal approximation, not an exact
integral or a treatment of boundaries or multiple modes. The nonlinear test uses
`x ~ N(0,0.25)` and `y=1 ~ N(x**2,1)`: its unique mode is zero, its full
curvature is 2 (Gauss–Newton gives 4), and Laplace overestimates independently
adaptively integrated evidence by about 19.25%. Exact Gaussian and independently
marginalized linear-regression tests also check normalization under parameter
rescaling.

`fortprob_affine_filter` provides conditional Gaussian filtering when the state
mean depends affinely on unknown regression coefficients. Use
`gaussian_interval_t%init(F,A,Qend,qavg,qcross,measurement_variance,status)` to
configure a transition. Here `F` is the state-transition matrix, `A` the
observation row, `Qend` the end-state process covariance, `qavg` the process
variance in the observation, and `qcross` their cross-covariance. The
matrix arguments `F`, `Qend` and `P_initial` have shape `(nstate,nstate)`;
`A` and `qcross` have shape `(nstate)`. The joint
process covariance must be positive semidefinite; measurement variance must be
strictly positive. This permits interval-averaged observations whose process
noise correlates with the state at the interval end.

Initialize `affine_gaussian_filter_t` with
`init(P_initial,ncoeff,status)`. Its conditional state has mean `c-D*beta`,
initially zero, and covariance `P_initial`. Call
`step(interval,design,y_minus_base,condition,offset,effective_design,innovation_variance,status)`
in observation order. The caller allocates `effective_design(ncoeff)`.
Before the optional observation update, the outputs describe the observation
mean `base+offset+effective_design*beta` and its conditional innovation variance.
If `condition` is false, the target is ignored and the state only propagates.
If true, the update conditions on that target while retaining its affine
coefficient dependence. Internally, the state covariance is stored as a factor
times its transpose. QR compression of the Joseph update's factor columns
preserves this representation without subtracting covariance entries. This
also avoids indefinite covariances caused by rounding extremely small entries
after long decay. Failed steps leave the filter unchanged.

For each conditioned row, whiten the residual as
`(y_minus_base-offset)/sqrt(innovation_variance)` and divide
every entry of `effective_design` by the same square root. Retain the sum of log innovation variances when
computing normalized evidence. Coefficient inference remains the caller's
responsibility. Marginal prediction variance additionally requires the
effective design times coefficient covariance times its transpose; it must not
add an independently estimated state variance again. The interface provides
conditional Gaussian moments, not nonlinear filtering or filter derivatives.

The affine filter initializes covariance factors with
[FortNum](https://github.com/lazy-fortran/fortnum)'s pivoted semidefinite
factorization and its symmetry and reconstruction checks. Successful
initialization establishes a semidefinite approximation within its documented
roundoff bound, not exact semidefiniteness of the input. Projected state
variance is a sum of squares; measurement variance must remain positive.
No diagonal jitter or eigenvalue clipping is introduced.

For the other covariance APIs, inputs are not repaired.
Symmetry and positive semidefiniteness checks allow roundoff of
`128*epsilon(real64)*dimension*max(abs(cov))`; conditional covariance uses
at least the prior covariance scale after subtraction. Returned covariances
are unchanged, so roundoff can leave a tiny negative variance, which the
interval helper rejects. Singular observation systems return an error;
rank-aware conditioning is not implemented.

Build and test with `fpm test`. The test-drive suite checks exact scalar and
correlated multivariate answers, negative cases, and directional derivatives
against independent central finite differences. All fixtures are synthetic.
On Homebrew systems, fo may need
`PKG_CONFIG_PATH=/opt/homebrew/opt/openblas/lib/pkgconfig` to select BLAS.

Scalar interval and density leaves are generated by
[FortSym](https://github.com/lazy-fortran/fortsym) through
[FortGen](https://github.com/lazy-fortran/fortgen). Regenerate using
`gen/regenerate.sh` with the pinned FortSym checkout and check
`git diff --exit-code -- src/generated`. Ordinary consumers build checked-in
sources without a computer algebra system.
