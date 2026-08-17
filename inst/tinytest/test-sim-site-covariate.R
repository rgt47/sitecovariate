library(tinytest)

# Whitepaper issue 2.8: the simulation engine that generates every
# reported number in the manuscript previously had zero test
# coverage. These tests exercise `sim_site_covariate()` at a small,
# fixed-seed `n_sim` and check output shape, column types, and basic
# sanity (near-zero bias when the estimand is well identified).

script_path <- here::here("analysis", "scripts",
                           "sim_site_covariate.R")
if (!file.exists(script_path)) {
  exit_file("sim_site_covariate.R not found; skipping engine tests")
}
source(script_path)

RNGkind("L'Ecuyer-CMRG")
set.seed(20260816)
res <- sim_site_covariate(
  n_sim = 25, n_sites = 10,
  site_size = "equal", true_trt = 0.30,
  site_sd = 0.50, trt_site_sd = 0.0
)

expect_equal(
  nrow(res), 25 * 4,
  info = "one row per model (4) per replication"
)
expect_true(
  all(c("model", "estimate", "se", "pvalue", "converged",
        "status", "sim") %in% names(res)),
  info = "output carries model, estimate, se, pvalue, converged, status, sim"
)
expect_true(
  all(res$model %in% c("no_site", "site_fixed", "site_random",
                        "site_x_trt")),
  info = "model labels are restricted to the four analysis strategies"
)
expect_true(
  all(res$status %in% c("converged", "singular", "nonconverged",
                         "error")),
  info = "status is one of the four documented outcome categories (issue 2.6)"
)
expect_true(
  is.logical(res$converged),
  info = "converged is a logical fit-succeeded indicator"
)

fit_ok <- res[res$converged, ]
expect_true(
  all(is.finite(fit_ok$estimate)),
  info = "estimates are finite for every fit that did not error"
)
mean_bias <- mean(fit_ok$estimate - 0.30)
expect_true(
  abs(mean_bias) < 0.15,
  info = paste(
    "mean bias across 25 reps is within a generous tolerance of",
    "zero (unbiased estimator, small-n_sim Monte Carlo noise)"
  )
)

# Reproducibility: an identical seed and call must reproduce
# identical estimates (issue 2.2 remediation depends on this).
set.seed(20260816)
res2 <- sim_site_covariate(
  n_sim = 25, n_sites = 10,
  site_size = "equal", true_trt = 0.30,
  site_sd = 0.50, trt_site_sd = 0.0
)
expect_equal(
  res$estimate, res2$estimate,
  info = "identical seed reproduces identical treatment-effect estimates"
)

# Sparse-site scenario support added for whitepaper issue 2.5.
set.seed(1)
res_sparse <- sim_site_covariate(
  n_sim = 5, n_sites = 6,
  site_size = "sparse", true_trt = 0.30,
  site_sd = 0.50, trt_site_sd = 0.0
)
expect_equal(
  nrow(res_sparse), 5 * 4,
  info = "sparse site_size option runs and returns expected row count"
)
