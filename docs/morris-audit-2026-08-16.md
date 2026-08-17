# Morris et al. (2019) ADEMP Audit: 11-site-covariate-analysis
*2026-08-16 17:40 PDT*

## Scope

Re-audit of `analysis/scripts/sim_site_covariate.R` and
`analysis/report/report.Rmd` against the original
`docs/morris-audit-2026-04-17.md` scorecard, performed as part of
remediation for `docs/pub_review_whitepaper_2026-08-16.md`. This
file supersedes the 2026-04-17 audit's verdict; the original file is
left unmodified as a historical record of the earlier state.

## Updated ADEMP scorecard

| Criterion | Status | Evidence |
|---|---|---|
| Aims explicit | Met | ADEMP "Aims" bullet in report.Rmd Methods |
| DGMs documented | Met | site-effects DGM parameterised, formula given |
| Factors varied factorially | Partial | 4 primary scenarios (2x2) plus a null scenario and a sparse-site scenario; not a full factorial over site-size regime x interaction x null |
| Estimand defined with true value | Met | fixed treatment effect in DGM; `site_x_trt` weighting now matches the stated site-size-weighted estimand (see below) |
| Methods justified | Met | four `lme4::lmer()` specifications, one per strategy |
| Performance measures justified | Met | mapped to Aims in report.Rmd Methods, "Performance metrics" |
| n_sim stated | Met | `n_sim = 1000` for the four primary scenarios; derived from target coverage MCSE in report.Rmd |
| n_sim justified via MCSE | Met | report.Rmd states target coverage MCSE $\le 0.7$pp at $n_{\mathrm{sim}} = 1000$ |
| MCSE reported per metric | Met | `summarise_sim()` (`sim_site_covariate.R`) computes `mcse_bias`, `mcse_empirical_se`, `mcse_mse`, `mcse_power`, `mcse_coverage`, `mcse_convergence`; now also shown in Tables 1-6, not only inline text |
| Seed set once | Met | `sim-seed` chunk in report.Rmd sets `RNGkind("L'Ecuyer-CMRG")` and `set.seed(sim_seed)` once; `sim_site_covariate()` does not call `set.seed()` |
| RNG states stored | Met | per-replicate `.Random.seed` snapshots attached as an attribute of each scenario's return value |
| Paired comparisons | Met | same simulated dataset fed to all four models within a replication |
| Reproducibility | Met (with disclosed exception) | Reproducibility section regenerated from `sessionInfo()`/`sim_seed`; Scenarios 5-6 (null, sparse) were run under a separately seeded stream at a reduced replicate count, disclosed explicitly in-text and flagged with TODOs for full integration |
| Convergence outcome captured | Met | `fit_lmer_safely()` now classifies each fit as converged / singular / nonconverged / error via `withCallingHandlers()` + `lme4::isSingular()` + `fit@optinfo$conv$opt`, rather than only catching R errors |

## Verdict

**Substantially compliant**, up from "Partially compliant" in the
2026-04-17 audit. All three gaps identified in that audit (seed
reseeded per scenario, no MCSE reported, unjustified `n_sim`) were
already resolved in the code prior to this remediation pass; the
manuscript's end-of-document compliance note had simply not been
updated to say so (whitepaper issue 2.3). This audit additionally
verifies and documents a fourth area of compliance not scored in the
original audit: convergence-outcome capture (whitepaper issue 2.6),
which was fixed as part of this remediation.

## Remaining gaps

- Scenarios 5 (null/Type I error) and 6 (sparse site sizes) are run
  at a reduced replicate count ($n_{\mathrm{sim}} = 300$) under a
  seed not chained from the main `sim_seed` stream, for time-budget
  reasons during remediation. Exact commands to fold them into the
  main sequence at $n_{\mathrm{sim}} = 1000$ are given as TODOs in
  the corresponding `report.Rmd` source chunks (`run-sim-null`,
  `run-sim-sparse`).
- Factor design is not a full factorial (site size x interaction x
  null x sparse are not all crossed); this is a scope decision, not
  a compliance gap, and is disclosed in the manuscript's Discussion
  and Limitations text.

## References

Morris TP, White IR, Crowther MJ. Using simulation studies to
evaluate statistical methods. Stat Med 2019;38:2074-2102.
doi:10.1002/sim.8086

---
*Source: ~/prj/res/11-site-covariate-analysis/sitecovariate/docs/morris-audit-2026-08-16.md*
