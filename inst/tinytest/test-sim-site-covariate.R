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
  nrow(res), 25 * 5,
  info = "one row per model (5) per replication"
)
expect_true(
  all(c("model", "estimate", "se", "pvalue", "converged",
        "status", "sim") %in% names(res)),
  info = "output carries model, estimate, se, pvalue, converged, status, sim"
)
expect_true(
  all(res$model %in% c("no_site", "site_fixed", "site_random",
                        "site_slope", "site_x_trt")),
  info = "model labels are restricted to the five analysis strategies"
)
expect_true(
  all(res$cond_trt == 0.30),
  info = "conditional truth equals the marginal effect when trt_site_sd = 0"
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
  nrow(res_sparse), 5 * 5,
  info = "sparse site_size option runs and returns expected row count"
)

# Regression anchor: replicate 1 of the shipped null run
# (seed 20260816, L'Ecuyer-CMRG). Adding the site_slope model must not
# change the stream or the other four estimates.
RNGkind("L'Ecuyer-CMRG")
set.seed(20260816)
res_null <- suppressMessages(sim_site_covariate(
  n_sim = 1, n_sites = 20, site_size = "equal", true_trt = 0,
  site_sd = 0.50, trt_site_sd = 0
))
expect_equal(
  res_null$estimate[res_null$model == "site_fixed"], -0.02931353,
  tolerance = 1e-6,
  info = "null-run replicate 1 reproduces under the pinned stream"
)

# Allocation: alternation is exactly balanced for even sites;
# permuted blocks of 4 keep both arms in every site of size >= 4 and
# can leave a site imbalanced.
sizes <- c(4, 5, 6, 7, 8, 20)
alt <- allocate_arms(sizes, "alternate")
site_id <- rep(seq_along(sizes), times = sizes)
imb_alt <- tapply(alt, site_id, function(a) abs(sum(a) - sum(1 - a)))
expect_equal(
  as.numeric(imb_alt), c(0, 1, 0, 1, 0, 0),
  info = "alternation: imbalance 0 for even sizes, 1 for odd sizes"
)
set.seed(11)
imb_blk <- replicate(200, {
  b <- allocate_arms(sizes, "block", block_size = 4)
  expect_true(
    all(tapply(b, site_id, function(a) any(a == 0) && any(a == 1))),
    info = "block allocation leaves both arms in every site"
  )
  max(tapply(b, site_id, function(a) abs(sum(a) - sum(1 - a))))
})
expect_true(
  any(imb_blk >= 2),
  info = "block allocation can produce within-site imbalance >= 2"
)
expect_true(
  all(imb_blk <= 2),
  info = "block size 4 bounds within-site imbalance at 2"
)

# Conditional truth under treatment-by-site heterogeneity.
set.seed(3)
res_int <- suppressMessages(sim_site_covariate(
  n_sim = 2, n_sites = 6, site_size = "equal", true_trt = 0.30,
  site_sd = 0.50, trt_site_sd = 0.15
))
expect_true(
  all(res_int$cond_trt != 0.30),
  info = "conditional truth departs from 0.30 when trt_site_sd > 0"
)
