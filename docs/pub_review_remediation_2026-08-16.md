# pub_review Remediation Log: 11-site-covariate-analysis
*2026-08-16 18:10 PDT*

Remediation of `docs/pub_review_whitepaper_2026-08-16.md` (initial
referee-grade review, verdict: major revision). This log records
what was actually fixed, what was deferred, and what was newly
discovered while fixing. It does not edit the whitepaper, which
remains the review record.

## 1. Fixed

**Correctness tier (whitepaper section 4a)**

- **2.1 Type I error never simulated.** Added a null-effect
  scenario ($\beta = 0$) as new "Scenario 5" in
  `analysis/report/report.Rmd` (`run-sim-null`/`tab-null` chunks),
  and narrowed the Abstract/Methods claims to state explicitly
  which scenarios test power vs. type I error. Real numbers were
  generated (not fabricated) at a reduced replicate count
  ($n_{\mathrm{sim}} = 300$, seed 20260816) via
  `analysis/data/derived_data/sim_typeI_null_n300.rds`; a TODO with
  the exact full-scale command is left in the `run-sim-null` chunk.
  `[verified]` (ran the code; saved real output; confirmed
  non-fabricated, non-NA results).

- **2.2 Reproducibility section factual errors.** Rewrote the
  "Reproducibility" section of `analysis/report/report.Rmd` to be
  generated from `sessionInfo()`/`R.version.string`,
  `utils::packageVersion()`, and the `sim_seed` variable (now
  stored explicitly in the `sim-seed` chunk) rather than hand-typed
  prose. Corrected the false package list (removed `MASS`, `furrr`,
  `future`; the pipeline is a plain serial `for` loop with base
  `rnorm()`), the R version claim (now reports the actual render
  R version and separately notes the Docker/renv pin of 4.6.0), and
  the seed value (now `` `r sim_seed` ``, reading the true
  `sim-seed` chunk variable, eliminating the possibility of drift).
  `[verified]` (parsed the modified .Rmd with `knitr::purl()` +
  `parse()`; grepped `sim_site_covariate.R`/`DESCRIPTION` to
  confirm no `MASS`/`furrr`/`future` usage).

- **2.3 Stale ADEMP compliance note.** Re-ran the Morris audit
  against current code; wrote
  `docs/morris-audit-2026-08-16.md` with a corrected scorecard and
  verdict ("Substantially compliant," up from "Partially
  compliant"), documenting that the three original gaps (seed
  reseeded per scenario, no MCSE, unjustified `n_sim`) were already
  resolved in the code, plus the newly-fixed convergence-capture
  gap (2.6). Rewrote the end-of-document note in
  `analysis/report/report.Rmd` to match. The original
  `docs/morris-audit-2026-04-17.md` was left unmodified as a
  historical record. `[verified]` (re-inspected the current code
  against every scorecard row).

- **2.6 Convergence check cannot detect `lme4` warnings.** Replaced
  the bare `tryCatch(error = ...)` fitting calls in
  `analysis/scripts/sim_site_covariate.R` with a
  `fit_lmer_safely()` helper that captures warnings via
  `withCallingHandlers()` and checks `lme4::isSingular()` and
  `fit@optinfo$conv$opt`, classifying every fit as
  converged/singular/nonconverged/error. `summarise_sim()` now
  reports `convergence_rate`, `singular_rate`, and
  `nonconverged_rate` from this classification. Updated Methods
  ("Performance metrics") and the ADEMP compliance note to
  describe the new logic. `[verified]` (unit-tested via a
  deliberately near-singular fixture that correctly triggers
  `isSingular() == TRUE`; confirmed via `tinytest` that
  `summarise_sim()` correctly tabulates all four status categories
  on a synthetic input; full `tinytest::run_test_dir()` run: 25/25
  assertions pass).

- **2.7 Estimand mismatch in `site_x_trt` weighting.** Changed
  `extract_trt()`'s `site_x_trt` branch in
  `analysis/scripts/sim_site_covariate.R` to weight by
  `table(dat$site)` (total site size) instead of
  `table(dat$site[dat$trt == "1"])` (treated-arm count only),
  matching the estimand stated in the Methods section ("the
  site-size-weighted average of site-specific effects"). Also
  removed the now-genuinely-dead overwritten `contrast <-` line
  identified in whitepaper 3.1 (it was superseded by the very next
  line in the pre-fix code and is gone in the rewritten branch).
  `[verified]` (new tinytest regression test
  `inst/tinytest/test-site-x-trt-weighting.R` confirms the fixed
  weighting matches an independent
  `emmeans::emmeans(weights = "proportional")` calculation to
  1e-8, and that treated-arm-only weighting differs from
  total-site-size weighting on the same fixture, confirming the
  two schemes are not numerically interchangeable).

**Acceptance tier (whitepaper section 4b)**

- **2.4 MCSEs omitted from result tables.** Added a `fmt_mcse()`
  helper and MCSE columns (as "estimate (MCSE)" cells) to Bias,
  Power, and Coverage in all six result tables
  (`tab-equal-no`, `tab-unequal-no`, `tab-equal-int`,
  `tab-unequal-int`, `tab-null`, `tab-sparse`) in
  `analysis/report/report.Rmd`. Rewrote the first Discussion
  paragraph and the Abstract's Results/Conclusions to state
  plainly that the fixed-versus-random gap is not distinguishable
  from Monte Carlo noise in any scenario examined (moderate or
  sparse), rather than asserting an "amplified" advantage not
  supported by the reported numbers. `[verified]` (`knitr::purl()`
  + `parse()` confirms the new chunks are syntactically valid; the
  underlying MCSE computations were already covered by the
  existing `summarise_sim()` test suite).

- **2.5 No sparse-site scenario.** Added `site_size = "sparse"`
  support to `sim_site_covariate()` (4-8 patients per site, drawn
  once per scenario) and a new "Scenario 6" table/discussion in
  `analysis/report/report.Rmd`. Ran it at a reduced replicate count
  ($n_{\mathrm{sim}} = 300$, seed 20260817), saved to
  `analysis/data/derived_data/sim_sparse_n300.rds`; a TODO with the
  full-scale command is left in the `run-sim-sparse` chunk. The
  real result: even in the sparse-site regime, "No site" vs.
  site-adjusted power differs substantially (0.473 vs. ~0.66), but
  fixed-vs-random power (0.653 vs. 0.657) still does not clear one
  MCSE, so the manuscript's framing was narrowed accordingly rather
  than left overclaiming. `[verified]` (ran the code; real,
  non-fabricated output; `tinytest` confirms the `"sparse"` option
  runs and returns the expected row count).

- **2.8 Placeholder test suite.** Replaced
  `inst/tinytest/test_basic.R` (`expect_true(TRUE)`) with three
  real test files: `test-sim-site-covariate.R` (engine output
  shape/columns/status categories, reproducibility under a fixed
  seed, sparse-scenario support), `test-summarise-sim.R`
  (`summarise_sim()` verified against hand-computed values on a
  toy input, including the four status-category rates), and
  `test-site-x-trt-weighting.R` (the site-size-weighted contrast
  verified against an independent `emmeans` calculation). Added
  `emmeans` to `DESCRIPTION` Suggests for this test. `[verified]`
  (`Rscript -e 'tinytest::run_test_dir("inst/tinytest")'`: all ok,
  25 results, ~9s).

- **3.4 Unused bibliography entries.** Added a new paragraph in
  Introduction, "The case for adjusting for site," distinguishing
  site adjustment from the broader baseline-covariate-adjustment
  literature and citing 4 of the 7 previously-unused entries
  (`Pocock1975allocation`, `Fleiss1986design`,
  `wangAnalysisCovarianceRandomized2019`,
  `wangUtilitiesPitfallsStratified2022`). `[verified]` (grepped
  `report.Rmd` `@`-citations against `references.bib` keys post-edit
  to confirm the four now appear).

- **3.5 Undisclosed fixed site-size draw.** Added an explicit
  disclosure to "Simulation scenarios" in `analysis/report/report.Rmd`
  and a matching code comment in `sim_site_covariate.R` stating that
  the per-site sample-size vector (unequal and sparse conditions) is
  drawn once per scenario call and reused across all replications,
  not resampled per replicate. `[verified]` (text/code inspection).

**Polish tier (whitepaper section 4c)**

- **3.1 Dead code in `extract_trt()`.** Removed as part of the 2.7
  rewrite. `[verified]`.
- **3.2 Fragile positional table relabeling.** Replaced
  `mutate(model = c("No site", ...))` with an explicit
  `model_labels` key -> label lookup applied via
  `dplyr::recode(model, !!!model_labels)` in all six table chunks.
  `[verified]` (syntax-checked via `knitr::purl()` + `parse()`).
- **3.3 Misplaced section nesting.** The ADEMP compliance note was
  moved out from under "# References" into its own top-level
  "# ADEMP compliance statement" section; "# References" is now a
  proper standalone heading with a `<div id="refs">` anchor at the
  end of the document. `[verified]` (heading structure inspected
  post-edit).
- **3.6 Stale DESCRIPTION metadata.** Removed the unused
  `Roxygen: list(markdown = TRUE)` and `RoxygenNote: 7.2.0` fields,
  since `R/` is empty and no roxygenized functions exist.
  `Version: 0.0.0.9000` was left unchanged (release-tagging is an
  authorial decision, not a correctness issue). `[verified]`.

## 2. Deferred

- **Full-scale rerun of the four primary scenarios at
  $n_{\mathrm{sim}} = 1000$ with the corrected code.** The code
  fixes (2.6 convergence capture, 2.7 estimand weighting) change
  the simulation engine, and the on-disk `knitr` cache in
  `analysis/report/cache/` was generated under the pre-fix code
  (2026-04-19). That stale cache was deleted so a future render
  cannot silently reuse buggy-code results, but a full render was
  not attempted here: at ~7.3s per 20 replications observed
  locally, $n_{\mathrm{sim}} = 1000$ for the four primary scenarios
  is roughly 24-30 minutes of wall-clock R time, exceeding the
  budget for this pass. **Command:**
  `bash tools/render.sh analysis/report/report.Rmd`. The estimand
  and convergence fixes were independently verified on small/toy
  inputs (see 1. above) rather than via a full-scale rerun.
- **Folding Scenarios 5 (null) and 6 (sparse) into the main pinned
  `sim_seed` RNG stream at $n_{\mathrm{sim}} = 1000$.** Currently
  run at $n_{\mathrm{sim}} = 300$ under separately seeded streams
  (20260816, 20260817) for time-budget reasons, disclosed explicitly
  in the manuscript and the ADEMP compliance statement. **Commands**
  are given as TODOs directly in the `run-sim-null` and
  `run-sim-sparse` chunks of `analysis/report/report.Rmd`.
- **Manuscript title and full novelty-scope reframing** (whitepaper
  section 5, "Recommended framing"). The Abstract, Introduction
  "Present study," Discussion, and "The case for adjusting for
  site" were substantively narrowed to match the recommended
  framing (longitudinal extension of the cross-sectional
  fixed-vs-random literature, not "resolving" it), but the
  manuscript title ("Should site be a covariate in the analysis
  model?") was left unchanged. Retitling is an authorial/editorial
  decision about scope and target-journal fit that this
  remediation pass did not make unilaterally.
- **Pruning the remaining 3 unused bibliography entries**
  (`cliftonCorrelationBaselineScore2019`,
  `herschtalEffectDichotomizationSkewed2023`,
  `pirondiniCovariateAdjustmentCardiovascular2022`). Four of the
  seven flagged entries were woven into a new paragraph (see 1.
  above); these three (baseline-score correlation, dichotomization
  of skewed covariates, cardiovascular-trial covariate adjustment)
  did not fit naturally into any existing argument without forcing
  a citation. Left in `references.bib`, uncited, per whitepaper
  3.4's "cite... or prune" alternative; pruning is a one-line-per-entry
  deletion in `analysis/report/references.bib` if the author prefers
  that route instead.
- **Full factorial design** (site size x interaction x null x
  sparse) noted in the updated ADEMP audit as a disclosed scope
  decision, not attempted; would require additional scenario
  combinations and correspondingly more render time.

## 3. New issues found while fixing

- **Stale `knitr` cache would have silently served pre-fix results.**
  `analysis/report/cache/` contained cached `run-sim-*` chunk output
  from 2026-04-19, keyed only on chunk source text (not on the
  content of the sourced `sim_site_covariate.R`). Editing the
  simulation engine without deleting this cache would have left the
  manuscript's Tables 1-4 numerically identical to the pre-fix
  (buggy-convergence, wrong-estimand-weighting) code even after a
  "successful" re-render. The cache directory was deleted as part
  of this remediation; the next render will recompute from the
  corrected code. Not flagged in the original whitepaper.
- **`analysis/data/derived_data/sim_11.rds`** is an orphaned cached
  object not referenced anywhere in `report.Rmd` or
  `sim_site_covariate.R` (confirmed via grep). Left in place;
  harmless but unexplained provenance, worth a cleanup pass.
- **`report.pdf`/`report.tex`/`report_files/`** in
  `analysis/report/` are now stale relative to the edited
  `report.Rmd` source (they reflect the pre-remediation manuscript).
  A re-render via `bash tools/render.sh analysis/report/report.Rmd`
  is required before these artifacts can be trusted or
  distributed.
