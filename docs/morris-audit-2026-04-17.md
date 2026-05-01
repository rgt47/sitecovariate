# Morris et al. (2019) ADEMP Audit: 11-site-covariate-analysis
*2026-04-17 09:02 PDT*

## Scope

Files audited:

- `analysis/scripts/sim_site_covariate.R`
- `analysis/report/report.Rmd`

Package `R/` directory is empty; simulation source lives entirely in
`analysis/scripts/`.

## ADEMP scorecard

| Criterion | Status | Evidence |
|---|---|---|
| Aims explicit | Partial | goals in prose; no ADEMP header |
| DGMs documented | Met | site-effects DGM parameterised |
| Factors varied factorially | Partial | 4 scenarios, not full factorial |
| Estimand defined with true value | Met | fixed treatment effect in DGM |
| Methods justified | Met | with / without site covariate methods |
| Performance measures justified | Partial | listed without mapping to aims |
| n_sim stated | Met | `n_sim = 1000` hardcoded at `sim_site_covariate.R:2` |
| n_sim justified via MCSE | Not met | no derivation |
| MCSE reported per metric | Not met | `summarise_sim()` at `sim_site_covariate.R:166-184` returns no MCSE |
| Seed set once | Not met | `set.seed()` at `sim_site_covariate.R:13` fires on every function call; `report.Rmd:306-418` calls the function 4 times (scenarios 1-4) — each re-seeds |
| RNG states stored | Not met | not stored |
| Paired comparisons | Met | same data fed to both methods within rep |
| Reproducibility | Partial | seed set; RNGkind not pinned |

## Overall verdict

**Partially compliant.**

## Gaps

- `set.seed()` at `sim_site_covariate.R:13` fires inside the worker;
  `report.Rmd:306-418` invokes the worker four times (once per
  scenario), so the RNG state is reset before each scenario. Morris
  §4.1: one seed at program start.
- No Monte Carlo SE on any metric (`summarise_sim()` L166-184).
- `n_sim = 1000` hardcoded without an MCSE target derivation.
- `RNGkind()` not pinned; no per-rep `.Random.seed` snapshot.

## Remediation plan

1. Move `set.seed()` out of the worker function. Add a single
   `set.seed(20260310)` at the top of the main `01_run_simulation.R`
   (create if missing) or at the first sim-using chunk of
   `report.Rmd`.
2. Generate per-replicate L'Ecuyer streams from the master seed up
   front and pass them into each scenario so that scenarios 1-4
   receive independent but reproducible streams.
3. Add `mcse_*` columns to `summarise_sim()` in
   `sim_site_covariate.R:166-184` per Morris Table 6.
4. Derive `n_sim` from target MCSE at script top. For coverage MCSE
   ≤ 1 pp at 0.95, need n_sim ≥ 475.
5. Create `R/performance.R` with exportable MCSE helpers
   (`mcse_bias`, `mcse_coverage`, `mcse_rejection`, `mcse_mse`).
6. Pin `RNGkind("L'Ecuyer-CMRG")` and store `.Random.seed` per rep.
7. Add ADEMP Methods section to `report.Rmd` and wire MCSE into the
   four result tables.

## References

Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate
statistical methods. Stat Med 2019;38:2074-2102. doi:10.1002/sim.8086

---
*Source: ~/prj/res/11-site-covariate-analysis/sitecovariate/docs/morris-audit-2026-04-17.md*
