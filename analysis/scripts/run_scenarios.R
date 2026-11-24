# Driver for every simulation scenario reported in
# analysis/report/report.Rmd. Each scenario runs in its own R process
# and writes one .rds file, which the manuscript reads.
#
# Usage, from the inner repository root:
#   Rscript analysis/scripts/run_scenarios.R <scenario_id> [n_sim] [out_dir]
#   Rscript analysis/scripts/run_scenarios.R all
#
# RNG: one master seed, `sim_seed`, under L'Ecuyer-CMRG. Scenario k
# uses the k-th independent substream, obtained by applying
# `parallel::nextRNGStream()` k times to the state after
# `set.seed(sim_seed)`. Scenarios can therefore run in parallel and
# in any order, and each result depends only on (sim_seed, k).

script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f) == 1 && basename(f) == "run_scenarios.R") {
    dirname(normalizePath(f))
  } else {
    "analysis/scripts"
  }
})
source(file.path(script_dir, "sim_site_covariate.R"))

sim_seed <- 20260312

scenarios <- list(
  eq_no = list(stream = 1, n_sim = 1000, site_size = "equal",
               true_trt = 0.30, trt_site_sd = 0.00),
  uneq_no = list(stream = 2, n_sim = 1000, site_size = "unequal",
                 true_trt = 0.30, trt_site_sd = 0.00),
  eq_int = list(stream = 3, n_sim = 1000, site_size = "equal",
                true_trt = 0.30, trt_site_sd = 0.15),
  uneq_int = list(stream = 4, n_sim = 1000, site_size = "unequal",
                  true_trt = 0.30, trt_site_sd = 0.15),
  null_a = list(stream = 5, n_sim = 1000, site_size = "equal",
                true_trt = 0.00, trt_site_sd = 0.00),
  null_b = list(stream = 6, n_sim = 1000, site_size = "equal",
                true_trt = 0.00, trt_site_sd = 0.00),
  sparse = list(stream = 7, n_sim = 1000, site_size = "sparse",
                true_trt = 0.30, trt_site_sd = 0.00),
  sparse_block = list(stream = 8, n_sim = 1000, site_size = "sparse",
                      true_trt = 0.30, trt_site_sd = 0.00,
                      allocation = "block", block_size = 4)
)

stream_seed <- function(seed, k) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  s <- .Random.seed
  for (i in seq_len(k)) s <- parallel::nextRNGStream(s)
  s
}

run_one <- function(id, n_sim = NULL, out_dir = NULL) {
  sc <- scenarios[[id]]
  stopifnot(!is.null(sc))
  if (is.null(n_sim)) n_sim <- sc$n_sim
  if (is.null(out_dir)) {
    out_dir <- file.path(script_dir, "..", "data", "derived_data")
  }
  assign(".Random.seed", stream_seed(sim_seed, sc$stream),
         envir = globalenv())
  t0 <- Sys.time()
  res <- suppressMessages(sim_site_covariate(
    n_sim = n_sim, n_sites = 20, site_size = sc$site_size,
    true_trt = sc$true_trt, site_sd = 0.50,
    trt_site_sd = sc$trt_site_sd,
    allocation = if (is.null(sc$allocation)) "alternate" else
      sc$allocation,
    block_size = if (is.null(sc$block_size)) 4 else sc$block_size
  ))
  attr(res, "scenario") <- id
  attr(res, "sim_seed") <- sim_seed
  attr(res, "stream") <- sc$stream
  attr(res, "rng_kind") <- RNGkind()
  stopifnot(length(unique(res$sim)) == n_sim)
  path <- file.path(out_dir, paste0("sim_", id, ".rds"))
  saveRDS(res, path)
  message(sprintf("%s: n_sim = %d, %.1f min -> %s", id, n_sim,
                  as.numeric(difftime(Sys.time(), t0, units = "mins")),
                  normalizePath(path)))
  invisible(res)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) stop("give a scenario id or 'all'")
  ids <- if (args[1] == "all") names(scenarios) else args[1]
  n_sim <- if (length(args) >= 2) as.integer(args[2]) else NULL
  out_dir <- if (length(args) >= 3) args[3] else NULL
  for (id in ids) run_one(id, n_sim = n_sim, out_dir = out_dir)
}
