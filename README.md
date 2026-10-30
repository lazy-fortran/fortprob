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

There is no diagonal stabilization, eigenvalue clipping, or repair of inputs.
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
