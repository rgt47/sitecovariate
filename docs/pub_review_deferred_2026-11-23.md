# pub_review Deferred-Item Pass: 11-site-covariate-analysis
*2026-11-23 12:45 PST*

This pass reviews the Deferred section and the "New issues found
while fixing" section of `docs/pub_review_remediation_2026-08-16.md`
(which answers `docs/pub_review_whitepaper_2026-08-16.md`) against
the current working tree. Nothing was staged or committed.

Important context: the work that resolves most items is in the
working tree but is not committed. The last commit touching the
manuscript is `d743ab2` (2026-09-07 auto-backup). A later referee
cycle (`analysis/report/referee-report-2026-11-18.md`,
`analysis/report/referee-remediation-log-2026-11-19.md`) rewrote the
simulation design: a committed-in-tree driver
`analysis/scripts/run_scenarios.R` (untracked in git) produces eight
`sim_<id>.rds` files (also untracked), and the manuscript only reads
them. "RESOLVED-SINCE" below therefore cites files, not commits.

## 1. Inventory

| # | Source | Item | Class |
|---|---|---|---|
| D1 | Aug log, Deferred | Full-scale rerun of Scenarios 1-4 at n_sim = 1000 with the corrected engine | RESOLVED-SINCE |
| D2 | Aug log, Deferred | Fold null and sparse scenarios into the pinned seed at n_sim = 1000 | RESOLVED-SINCE |
| D3 | Aug log, Deferred | Manuscript title and novelty-scope reframing | AUTHOR-DECISION |
| D4 | Aug log, Deferred | Three uncited bibliography entries (Clifton, Herschtal, Pirondini) | AUTHOR-DECISION |
| D5 | Aug log, Deferred | Full factorial design | AUTHOR-DECISION |
| N1 | Aug log, New issues | Stale knitr cache could serve pre-fix results | RESOLVED-SINCE |
| N2 | Aug log, New issues | Orphaned `sim_11.rds` | AUTHOR-DECISION |
| N3 | Aug log, New issues | `report.pdf`/`report.tex` stale relative to source | DO-NOW |

Evidence for the RESOLVED-SINCE items:

- **D1.** `analysis/scripts/run_scenarios.R` defines `eq_no`,
  `uneq_no`, `eq_int`, `uneq_int` at `n_sim = 1000`. Structural check
  (I ran it): each of the four rds files has 1000 unique replicates,
  `sim_seed` attribute 20260312, streams 1-4, and no NA estimates.
  The files are dated 2026-11-19, after the engine fixes of
  2026-08-16. Per the November log they were produced with the
  current engine. [inspected; replicate counts verified]
- **D2.** The separate seeds 20260816 and 20260817 are gone. The
  null scenario now uses streams 5 and 6 (2 x 1000 = 2000
  replicates) and the sparse scenario stream 7 (1000), all under
  master seed 20260312 with L'Ecuyer-CMRG substreams. This replaces
  the "single pinned stream" proposed in August with per-scenario
  substreams, which is a sound design choice (each result depends
  only on seed and stream index). [verified by reading attributes of
  the rds files]
- **N1.** `analysis/report/cache/` does not exist. No chunk in
  `report.Rmd` sets `cache = TRUE`; the setup chunk retains an inert
  `cache.path = "cache/"` option. The manuscript runs no simulation.
  [inspected by grep]

## 2. Done now

- **N3. Re-render.** The previous render (2026-11-19 12:54) post-
  dated the Rmd edits, but `references.bib` was modified later the
  same day (19:31), including a corrected DOI for Kahan and Harhay
  (2015), `10.1016/j.jclinepi.2015.03.016`. I re-rendered through the
  wrapper:
  `bash tools/render.sh analysis/report/report.Rmd`.
  Output: `analysis/report/report.pdf`, `report.tex`, figures under
  `report_files/`, staged copy
  `analysis/report/share/report-2026-11-23-1244-d743ab2-wip.pdf`, and
  a new `share/MANIFEST.md` row. The corrected DOI now appears in
  `report.tex` (line 1234); no `??` unresolved references.
  [verified]
- **Tests.** `Rscript -e 'tinytest::run_test_dir("inst/tinytest")'`:
  all ok, 234 results. No code was changed in this pass. [verified]

## 3. Long runs

None of the August deferred items requires one. The November
referee log lists further deferred items that would; see section 5.

## 4. Author decisions

- **D3. Title.** The title is still "Should site be a covariate in
  the analysis model?". The Abstract now states that the analysis is
  equivalent to one of subject means and that the conclusions are
  conditional on near-balanced allocation. Options: keep the
  question title; or retitle to name the scope, for example
  "Fixed, random, and interaction site terms under near-balanced
  allocation: a simulation study". Recommendation: retitle. The
  question form promises a general answer that the conditional
  findings do not support.
- **D4. Uncited bibliography entries.**
  `cliftonCorrelationBaselineScore2019`,
  `herschtalEffectDichotomizationSkewed2023`, and
  `pirondiniCovariateAdjustmentCardiovascular2022` remain in
  `references.bib` and are still uncited. pandoc ignores them, so
  the PDF is unaffected. Options: prune, or leave. Recommendation:
  leave them unless the bib is shared with a journal submission; the
  entries appear to be Zotero-managed (`zotero-sync/`), so a manual
  prune may be reverted by the next sync.
- **D5. Full factorial.** The ADEMP statement discloses that
  Scenarios 5-7 vary one factor from a reference condition.
  Additional combinations (for example, heterogeneity x sparse, null
  x unequal) would each cost one run of up to about 10 minutes
  (November log timings). Recommendation: accept the disclosed
  one-factor design; the more valuable extra scenario is simple
  randomization (section 5), not a factorial fill-in.
- **N2. Orphaned result files.** `sim_11.rds` (April 2026, a list of
  four two-element lists, no attributes, unknown provenance) and the
  now-superseded `sim_typeI_null_n300.rds` and
  `sim_sparse_n300.rds` are read by nothing in `report.Rmd`,
  `analysis/scripts/`, `inst/`, or `Makefile` (grep). I did not
  delete them. Recommendation: delete all three, or move them to an
  `archive/` directory, before committing.

## 5. New issues found

- **Uncommitted provenance.** The driver `run_scenarios.R` and all
  eight result files that the manuscript reads are untracked, and
  the Rmd, engine, tests, and bib changes are uncommitted. The
  staged PDFs carry the version tag `d743ab2-wip`, which does not
  identify the code that produced them. The author should commit
  before any distribution.
- **Duplicate LaTeX label.** The render warns
  `Label 'tab:tab-status' multiply defined`. The Fit status table
  (chunk `tab-status`, `longtable = TRUE`) repeats its labelled
  caption on continuation pages. It is not cross-referenced, so the
  output is unaffected; a short caption or kableExtra's
  `repeat_header` handling would silence it. Polish only; not
  changed.
- **Related deferrals in the November referee log** (not named
  `pub_review_remediation_*`, so not inventoried above, but they are
  the live open items): simple randomization within site;
  Kenward-Roger or Satterthwaite degrees of freedom (requires a full
  rerun, eight runs of up to about 10 minutes each:
  `Rscript analysis/scripts/run_scenarios.R <id>` per scenario after
  the engine change); variance 0.15 versus SD 0.15 for the
  heterogeneity term; varying the number of sites; paired scenarios
  drawing different site-size vectors. The November simulations ran
  on host R 4.6.0, not in the renv container.
