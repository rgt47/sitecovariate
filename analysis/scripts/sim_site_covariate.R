sim_site_covariate <- function(
    n_sim = 1000,
    n_sites = 20,
    site_size = "equal",
    true_trt = 0.30,
    site_sd = 0.50,
    trt_site_sd = 0.0,
    n_time = 4,
    resid_sd = 1.0,
    seed = 20260312
) {

  set.seed(seed)

  n_per_site <- if (site_size == "equal") {
    rep(20, n_sites)
  } else {
    sizes <- round(rexp(n_sites, rate = 1 / 20))
    pmax(sizes, 4)
  }
  n_total <- sum(n_per_site)

  results <- vector("list", n_sim)

  for (s in seq_len(n_sim)) {

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

    fit_none <- tryCatch(
      lme4::lmer(y ~ trt + time + (1 | subj), data = dat),
      error = function(e) NULL
    )
    fit_fixed <- tryCatch(
      lme4::lmer(y ~ trt + time + site + (1 | subj), data = dat),
      error = function(e) NULL
    )
    fit_random <- tryCatch(
      lme4::lmer(y ~ trt + time + (1 | site) + (1 | subj),
                 data = dat),
      error = function(e) NULL
    )
    fit_interact <- tryCatch(
      lme4::lmer(y ~ trt * site + time + (1 | subj), data = dat),
      error = function(e) NULL
    )

    extract_trt <- function(fit, label) {
      if (is.null(fit)) {
        return(data.frame(
          model = label, estimate = NA, se = NA,
          pvalue = NA, converged = FALSE
        ))
      }
      cf <- summary(fit)$coefficients
      fe <- lme4::fixef(fit)

      if (label == "site_x_trt") {
        trt_idx <- grep("trt1", names(fe))
        if (length(trt_idx) == 0) {
          return(data.frame(
            model = label, estimate = NA, se = NA,
            pvalue = NA, converged = TRUE
          ))
        }
        site_counts <- table(dat$site[dat$trt == "1"])
        wts <- as.numeric(site_counts) /
          sum(as.numeric(site_counts))
        trt_coefs <- fe[trt_idx]
        main_eff <- trt_coefs[1]
        interact_eff <- c(0, trt_coefs[-1])
        avg_trt <- sum(wts * (main_eff + interact_eff))
        vcov_mat <- as.matrix(
          vcov(fit)[trt_idx, trt_idx]
        )
        contrast <- c(1, rep(0, length(trt_idx) - 1))
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
          se = se_avg, pvalue = pval, converged = TRUE
        ))
      }

      idx <- grep("^trt1$", rownames(cf))
      if (length(idx) == 0) {
        return(data.frame(
          model = label, estimate = NA, se = NA,
          pvalue = NA, converged = TRUE
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
        pvalue = pval, converged = TRUE
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
    n_time = n_time, resid_sd = resid_sd, seed = seed
  )

  out
}

summarise_sim <- function(sim_results, true_trt = 0.30) {
  sim_results |>
    dplyr::filter(converged) |>
    dplyr::group_by(model) |>
    dplyr::summarise(
      bias = mean(estimate - true_trt, na.rm = TRUE),
      empirical_se = sd(estimate, na.rm = TRUE),
      mean_model_se = mean(se, na.rm = TRUE),
      power = mean(pvalue < 0.05, na.rm = TRUE),
      coverage = mean(
        (estimate - 1.96 * se <= true_trt) &
          (estimate + 1.96 * se >= true_trt),
        na.rm = TRUE
      ),
      convergence = mean(converged, na.rm = TRUE),
      n_converged = sum(converged),
      .groups = "drop"
    )
}
