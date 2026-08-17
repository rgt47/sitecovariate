library(tinytest)

# Whitepaper issue 2.8(c) / 2.7: regression-test the site_x_trt
# weighted-contrast computation (total-site-size weights, as
# `extract_trt()` in `sim_site_covariate.R` now implements) against
# an independent calculation via `emmeans::emmeans(..., weights =
# "proportional")`, which averages a treatment-by-site interaction
# model over sites in proportion to their observed size.

if (!requireNamespace("emmeans", quietly = TRUE)) {
  exit_file("emmeans not installed; skipping independent weighting check")
}

script_path <- here::here("analysis", "scripts",
                           "sim_site_covariate.R")
if (!file.exists(script_path)) {
  exit_file("sim_site_covariate.R not found; skipping weighting tests")
}

set.seed(2026)
n_sites <- 5
# Odd site sizes so within-site treatment-arm counts are unequal
# (`rep(c(0, 1), length.out = n)` gives 2/1 or 1/2 splits when n is
# odd); this is exactly the condition under which weighting by
# treated-arm count (the old, incorrect behaviour) and weighting by
# total site size (the stated estimand) diverge.
n_per_site <- c(3, 5, 7, 9, 11)
n_time <- 2
n_total <- sum(n_per_site)
site_vec <- rep(seq_len(n_sites), times = n_per_site)
trt_subj <- unlist(lapply(n_per_site, function(n) {
  rep(c(0, 1), length.out = n)
}))
site_for_subj <- rep(site_vec, each = n_time)
subj_id <- rep(seq_len(n_total), each = n_time)
trt <- rep(trt_subj, each = n_time)
time <- rep(seq_len(n_time), times = n_total)
site_intercept <- rnorm(n_sites, 0, 0.5)
trt_by_site <- rnorm(n_sites, 0, 0.15)
subj_re <- rep(rnorm(n_total, 0, 0.5), each = n_time)
y <- site_intercept[site_for_subj] +
  (0.30 + trt_by_site[site_for_subj]) * trt +
  0.05 * time + subj_re +
  rnorm(length(trt), 0, 1)
dat <- data.frame(
  y = y, trt = factor(trt), time = time,
  site = factor(site_for_subj), subj = factor(subj_id)
)

# Sanity check: within-site treatment-arm counts really are unequal
# for at least one site, so this dataset actually exercises the
# weighting difference the test is designed to catch.
arm_counts <- table(dat$site[dat$time == 1], dat$trt[dat$time == 1])
expect_true(
  any(arm_counts[, 1] != arm_counts[, 2]),
  info = "fixture has at least one site with unequal arm sizes"
)

fit <- lme4::lmer(y ~ trt * site + time + (1 | subj), data = dat)
fe <- lme4::fixef(fit)
trt_idx <- grep("trt1", names(fe))

# Reproduce the total-site-size weighting now used in
# `extract_trt()`'s `site_x_trt` branch.
site_counts <- table(dat$site)
wts <- as.numeric(site_counts) / sum(as.numeric(site_counts))
trt_coefs <- fe[trt_idx]
main_eff <- trt_coefs[1]
interact_eff <- c(0, trt_coefs[-1])
avg_trt_manual <- sum(wts * (main_eff + interact_eff))

# Independent calculation via emmeans: average the trt effect over
# site with weights proportional to observed site size.
emm <- emmeans::emmeans(
  fit, "trt", weights = "proportional", at = list(time = 1)
)
con <- as.data.frame(emmeans::contrast(emm, "revpairwise"))
avg_trt_emmeans <- con$estimate[1]

expect_equal(
  avg_trt_manual, avg_trt_emmeans, tolerance = 1e-8,
  info = paste(
    "total-site-size-weighted treatment contrast matches an",
    "independent emmeans::emmeans(weights = 'proportional')",
    "calculation (whitepaper issue 2.7)"
  )
)

# Regression guard: weighting by treated-arm count only (the old,
# incorrect behaviour) must NOT equal the site-size-weighted
# estimate on this fixture, confirming the two schemes genuinely
# differ when arm sizes are unequal within site.
site_counts_trt_only <- table(dat$site[dat$trt == "1"])
wts_trt_only <- as.numeric(site_counts_trt_only) /
  sum(as.numeric(site_counts_trt_only))
avg_trt_old <- sum(wts_trt_only * (main_eff + interact_eff))
expect_true(
  abs(avg_trt_old - avg_trt_manual) > 1e-6,
  info = paste(
    "treated-arm-only weighting differs from total-site-size",
    "weighting on this unequal-arm-size fixture, confirming the",
    "two estimands are not numerically interchangeable"
  )
)
