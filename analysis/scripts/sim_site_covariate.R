sim_site_covariate <- function(
    n_sim = 1000,
    n_sites = 20,
    site_size = "equal",
    true_trt = 0.30,
    site_sd = 0.50,
    trt_site_sd = 0.0,
    n_time = 4,
    resid_sd = 1.0
) {

  # Morris, White, and Crowther (2019) §4.1: the RNG seed is set ONCE
  # by the caller; this function does NOT call `set.seed()`. Per-rep
  # RNG states are captured and attached as an attribute of the return
  # value for diagnostic reproducibility of any failing replication.
  rng_states <- vector("list", n_sim)

  # `n_per_site` is drawn once per scenario call and reused across all
  # `n_sim` replications (a fixed design factor, analogous to a fixed
  # total sample size), not resampled per replicate (disclosed in
  # Methods, "Simulation scenarios"; whitepaper issue 3.5).
  n_per_site <- if (site_size == "equal") {
    rep(20, n_sites)
  } else if (site_size == "sparse") {
    # Small, near-fixed sites per arm (4-8 patients per site total),
    # matching the sparse-site regime the Introduction cites via
    # @Kahan2015fewevents (median 3 events per centre-treatment arm)
    # as the condition under which fixed vs random effects are
    # expected to diverge (whitepaper issue 2.5).
    sample(4:8, n_sites, replace = TRUE)
  } else {
    sizes <- round(rexp(n_sites, rate = 1 / 20))
    pmax(sizes, 4)
  }
  n_total <- sum(n_per_site)

  results <- vector("list", n_sim)

  for (s in seq_len(n_sim)) {
    rng_states[[s]] <- .Random.seed

    site_intercept <- rnorm(n_sites, 0, site_sd)
    trt_by_site <- rnorm(n_sites, 0, trt_site_sd)

    site_vec <- rep(seq_len(n_sites), times = n_per_site)
    trt_subj <- unlist(lapply(n_per_site, function(n) {
      rep(c(0, 1), length.out = n)
    }))

    site_for_subj <- rep(site_vec, each = n_time)
    subj_id <- rep(seq_len(n_total), each = n_time)
    trt <- rep(trt_subj, each = n_time)
    time <- rep(seq_len(n_time), times = n_total)

    subj_re <- rep(rnorm(n_total, 0, 0.5), each = n_time)

    y <- site_intercept[site_for_subj] +
      (true_trt + trt_by_site[site_for_subj]) * trt +
      0.05 * time +
      subj_re +
      rnorm(length(site_for_subj), 0, resid_sd)

    dat <- data.frame(
      y = y,
      trt = factor(trt),
      time = time,
      site = factor(site_for_subj),
      subj = factor(subj_id)
    )

    # `lme4` signals optimizer failure-to-converge and singular-fit
    # conditions as *warnings*, not errors (see whitepaper issue 2.6);
    # a bare `tryCatch(error = ...)` therefore cannot detect them.
    # `fit_lmer_safely()` captures warnings with `withCallingHandlers()`
    # and additionally checks `lme4::isSingular()` and
    # `fit@optinfo$conv$opt`, classifying every fit into one of four
    # mutually exclusive outcome categories.
    fit_lmer_safely <- function(formula, data) {
      warn_msgs <- character(0)
      fit <- withCallingHandlers(
        tryCatch(
          lme4::lmer(formula, data = data),
          error = function(e) NULL
        ),
        warning = function(w) {
          warn_msgs <<- c(warn_msgs, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      if (is.null(fit)) {
        return(list(fit = NULL, status = "error"))
      }
      opt_failed <- isTRUE(fit@optinfo$conv$opt != 0) ||
        any(grepl("failed to converge", warn_msgs, fixed = TRUE))
      if (opt_failed) {
        return(list(fit = fit, status = "nonconverged"))
      }
      if (lme4::isSingular(fit, tol = 1e-4)) {
        return(list(fit = fit, status = "singular"))
      }
      list(fit = fit, status = "converged")
    }

    fit_none <- fit_lmer_safely(
      y ~ trt + time + (1 | subj), dat
    )
    fit_fixed <- fit_lmer_safely(
      y ~ trt + time + site + (1 | subj), dat
    )
    fit_random <- fit_lmer_safely(
      y ~ trt + time + (1 | site) + (1 | subj), dat
    )
    fit_interact <- fit_lmer_safely(
      y ~ trt * site + time + (1 | subj), dat
    )

    extract_trt <- function(fit_obj, label) {
      fit <- fit_obj$fit
      status <- fit_obj$status
      # `converged` retains its original meaning (fit produced a
      # usable estimate, i.e. did not error); `status` carries the
      # finer converged/singular/nonconverged/error classification
      # used by `summarise_sim()` for the convergence-rate metric.
      if (is.null(fit)) {
        return(data.frame(
          model = label, estimate = NA, se = NA,
          pvalue = NA, converged = FALSE, status = status
        ))
      }
      cf <- summary(fit)$coefficients
      fe <- lme4::fixef(fit)

      if (label == "site_x_trt") {
        trt_idx <- grep("trt1", names(fe))
        if (length(trt_idx) == 0) {
          return(data.frame(
            model = label, estimate = NA, se = NA,
            pvalue = NA, converged = TRUE, status = status
          ))
        }
        # Weighted by total site size (all subjects), matching the
        # "site-size-weighted average of site-specific effects"
        # estimand stated in the Methods section (report.Rmd, ADEMP
        # structure). Previously weighted by treated-arm counts only,
        # which is a distinct (and unstated) estimand whenever
        # within-site arm sizes are unequal (odd n per site).
        site_counts <- table(dat$site)
        wts <- as.numeric(site_counts) /
          sum(as.numeric(site_counts))
        trt_coefs <- fe[trt_idx]
        main_eff <- trt_coefs[1]
        interact_eff <- c(0, trt_coefs[-1])
        avg_trt <- sum(wts * (main_eff + interact_eff))
        vcov_mat <- as.matrix(
          vcov(fit)[trt_idx, trt_idx]
        )
        contrast <- wts[1] * c(1, rep(0, length(trt_idx) - 1))
        for (k in seq_along(trt_idx)[-1]) {
          cvec <- rep(0, length(trt_idx))
          cvec[1] <- wts[k]
          cvec[k] <- wts[k]
          contrast <- contrast + cvec
        }
        se_avg <- sqrt(
          as.numeric(t(contrast) %*% vcov_mat %*% contrast)
        )
        tval <- avg_trt / se_avg
        pval <- 2 * pt(
          abs(tval),
          df = nrow(dat) - length(fe),
          lower.tail = FALSE
        )
        return(data.frame(
          model = label, estimate = avg_trt,
          se = se_avg, pvalue = pval, converged = TRUE,
          status = status
        ))
      }

      idx <- grep("^trt1$", rownames(cf))
      if (length(idx) == 0) {
        return(data.frame(
          model = label, estimate = NA, se = NA,
          pvalue = NA, converged = TRUE, status = status
        ))
      }
      est <- cf[idx, "Estimate"]
      se_val <- cf[idx, "Std. Error"]
      tval <- cf[idx, "t value"]
      pval <- 2 * pt(abs(tval),
                      df = nrow(dat) - length(fe),
                      lower.tail = FALSE)
      data.frame(
        model = label, estimate = est, se = se_val,
        pvalue = pval, converged = TRUE, status = status
      )
    }

    row_none <- extract_trt(fit_none, "no_site")
    row_fixed <- extract_trt(fit_fixed, "site_fixed")
    row_random <- extract_trt(fit_random, "site_random")
    row_interact <- extract_trt(fit_interact, "site_x_trt")

    results[[s]] <- rbind(row_none, row_fixed,
                          row_random, row_interact)
  }

  out <- do.call(rbind, results)
  out$sim <- rep(seq_len(n_sim), each = 4)

  attr(out, "params") <- list(
    n_sim = n_sim, n_sites = n_sites,
    site_size = site_size, true_trt = true_trt,
    site_sd = site_sd, trt_site_sd = trt_site_sd,
    n_time = n_time, resid_sd = resid_sd
  )
  attr(out, "rng_states") <- rng_states

  out
}

summarise_sim <- function(sim_results, true_trt = 0.30) {
  # Morris et al. (2019) Table 6: report Monte Carlo SEs alongside
  # every performance estimate. Convergence is treated as a
  # first-class outcome per Morris §5.1 rather than silently
  # dropping non-converged replications.
  #
  # `n_fit` counts replications that produced a usable estimate (fit
  # did not error), whether or not the optimizer fully converged;
  # this is the denominator for bias/SE/power/coverage MCSEs, exactly
  # as before. `convergence_rate` (and its two companions,
  # `singular_rate` and `nonconverged_rate`) now come from the finer
  # `status` classification computed by `fit_lmer_safely()`
  # (converged / singular / nonconverged / error), which captures
  # `lme4` warnings and `isSingular()` rather than only R errors
  # (whitepaper issue 2.6).
  if (!"status" %in% names(sim_results)) {
    sim_results$status <- ifelse(
      sim_results$converged, "converged", "error"
    )
  }
  sim_results |>
    dplyr::group_by(model) |>
    dplyr::summarise(
      n_total = dplyr::n(),
      n_fit = sum(converged),
      n_converged = sum(status == "converged"),
      convergence_rate = n_converged / n_total,
      mcse_convergence = sqrt(
        convergence_rate * (1 - convergence_rate) / n_total
      ),
      singular_rate = sum(status == "singular") / n_total,
      nonconverged_rate = sum(status == "nonconverged") / n_total,
      error_rate = sum(status == "error") / n_total,
      bias = mean(estimate - true_trt, na.rm = TRUE),
      mcse_bias = stats::sd(estimate, na.rm = TRUE) /
        sqrt(n_fit),
      empirical_se = stats::sd(estimate, na.rm = TRUE),
      mcse_empirical_se = empirical_se /
        sqrt(2 * (n_fit - 1)),
      mse = mean((estimate - true_trt)^2, na.rm = TRUE),
      mcse_mse = sqrt(
        stats::var((estimate - true_trt)^2, na.rm = TRUE) /
          n_fit
      ),
      mean_model_se = mean(se, na.rm = TRUE),
      power = mean(pvalue < 0.05, na.rm = TRUE),
      mcse_power = sqrt(power * (1 - power) / n_fit),
      coverage = mean(
        (estimate - 1.96 * se <= true_trt) &
          (estimate + 1.96 * se >= true_trt),
        na.rm = TRUE
      ),
      mcse_coverage = sqrt(
        coverage * (1 - coverage) / n_fit
      ),
      .groups = "drop"
    )
}
