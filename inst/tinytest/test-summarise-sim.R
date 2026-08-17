library(tinytest)

# Whitepaper issue 2.8(b): verify summarise_sim() against
# hand-computed values on a toy input, independent of the simulation
# engine itself.

script_path <- here::here("analysis", "scripts",
                           "sim_site_covariate.R")
if (!file.exists(script_path)) {
  exit_file("sim_site_covariate.R not found; skipping summary tests")
}
source(script_path)

toy <- data.frame(
  model = rep("toy_model", 4),
  estimate = c(0.10, 0.20, 0.30, 0.40),
  se = c(0.05, 0.05, 0.05, 0.05),
  pvalue = c(0.20, 0.01, 0.03, 0.60),
  converged = c(TRUE, TRUE, TRUE, TRUE),
  status = c("converged", "converged", "converged", "converged")
)
true_trt <- 0.20

summ <- summarise_sim(toy, true_trt = true_trt)

expect_equal(
  summ$bias, mean(toy$estimate) - true_trt,
  info = "bias matches hand-computed mean(estimate) - true_trt"
)
expect_equal(
  summ$empirical_se, stats::sd(toy$estimate),
  info = "empirical_se matches hand-computed sd(estimate)"
)
expect_equal(
  summ$mse, mean((toy$estimate - true_trt)^2),
  info = "mse matches hand-computed mean squared error"
)
expect_equal(
  summ$power, mean(toy$pvalue < 0.05),
  info = "power matches hand-computed proportion with pvalue < 0.05 (0.5)"
)
expect_equal(summ$power, 0.5, info = "power is exactly 0.5 for this toy input")
expect_equal(
  summ$mean_model_se, mean(toy$se),
  info = "mean_model_se matches hand-computed mean(se)"
)
lo <- toy$estimate - 1.96 * toy$se
hi <- toy$estimate + 1.96 * toy$se
expect_equal(
  summ$coverage, mean(lo <= true_trt & hi >= true_trt),
  info = "coverage matches hand-computed CI-contains-truth proportion"
)
expect_equal(
  summ$convergence_rate, 1,
  info = "convergence_rate is 1 when every status is 'converged'"
)

# Distinct status categories are counted correctly (issue 2.6).
toy_status <- data.frame(
  model = rep("toy_model", 4),
  estimate = c(0.1, 0.2, 0.3, NA),
  se = c(0.05, 0.05, 0.05, NA),
  pvalue = c(0.2, 0.01, 0.03, NA),
  converged = c(TRUE, TRUE, TRUE, FALSE),
  status = c("converged", "singular", "nonconverged", "error")
)
summ_status <- summarise_sim(toy_status, true_trt = 0.20)
expect_equal(
  summ_status$convergence_rate, 0.25,
  info = "convergence_rate counts only status == 'converged' (1 of 4)"
)
expect_equal(
  summ_status$singular_rate, 0.25,
  info = "singular_rate counts status == 'singular' (1 of 4)"
)
expect_equal(
  summ_status$nonconverged_rate, 0.25,
  info = "nonconverged_rate counts status == 'nonconverged' (1 of 4)"
)
expect_equal(
  summ_status$error_rate, 0.25,
  info = "error_rate counts status == 'error' (1 of 4)"
)
expect_equal(
  summ_status$n_fit, 3,
  info = "n_fit counts fits that did not error (converged == TRUE), used as MCSE denominator"
)
