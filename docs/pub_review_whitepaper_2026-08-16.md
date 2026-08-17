# Referee-Grade Review: "Should Site Be a Covariate in the Analysis
Model?"

*Review date: 2026-08-16 16:07 PDT*

Workspace: `11-site-covariate-analysis` (inner repository:
`sitecovariate`). This is the first `pub_review` audit of this
workspace; no prior `docs/pub_review_whitepaper_*.md` file exists.

## 1. Summary of the work under review

The repository contains a single manuscript,
`analysis/report/report.Rmd` ("Should site be a covariate in the
analysis model?"), rendered to `analysis/report/report.pdf` /
`report.tex`. The paper reports a Monte Carlo simulation study
comparing four linear mixed-model strategies for handling
"site" (centre) in the analysis of a longitudinal two-arm
multicentre trial: (1) omitting site, (2) site as a fixed effect,
(3) site as a random effect, and (4) a treatment-by-site
interaction model. Four scenarios cross equal versus unequal site
sizes with absence versus presence of a true treatment-by-site
interaction variance component, each run at `n_sim = 1000`
replications via `analysis/scripts/sim_site_covariate.R`. The
paper is framed as an applied simulation study supporting
regulatory guidance (ICH E9, FDA, EMA) that stratification
variables belong in the analysis model, and its headline
conclusion is that site should enter as a random effect, with an
interaction term added only on strong a priori grounds. The
manuscript follows the Morris, White and Crowther (2019) ADEMP
structure and includes an internal audit note
(`docs/morris-audit-2026-04-17.md`) that is reproduced in
abbreviated form at the end of the report itself.

## 2. Major issues

**2.1 The Abstract promises evidence on Type I error, but no
scenario simulates the null.** *Location: Abstract ("consequences
for power, type I error..."); Methods, all four `run-sim-*`
chunks (report.Rmd lines 357-387).* All four simulation calls use
`true_trt = 0.30`; none sets the treatment effect to zero. The
reported "power" columns in Tables 1-4 are therefore true-positive
rates, not Type I error rates, and no false-positive-rate evidence
exists anywhere in the pipeline (inspected: `sim_site_covariate.R`
has no null-effect call site). A referee would require either an
explicit null scenario (`true_trt = 0`) with the corresponding
Type I error table, or removal of the Type I error claim from the
Abstract and Introduction. **Remediation:** add at least one
null-effect scenario per site-size/interaction combination (or at
minimum one representative null scenario) and report empirical
Type I error with its Wald-based MCSE, or strike the Type I error
claim.

**2.2 The Reproducibility section contains multiple factual
errors, verified against the code and build files.** *Location:
report.Rmd, "Reproducibility" section (around lines 811-829).*
Three specific claims were checked and each is wrong: (a) it
states "Principal packages used by the simulation pipeline are
`lme4`..., `MASS` (multivariate normal data generation), `furrr`
and `future` (parallel processing)" — `grep` over
`sim_site_covariate.R` and `DESCRIPTION` returns zero matches for
`MASS`, `furrr`, or `future` (verified); the simulation is a plain
serial `for` loop with `rnorm()`, not `MASS::mvrnorm()`, and there
is no parallel backend anywhere in the repository. (b) it states
"This manuscript was rendered with R 4.4.x," but the project's own
`Dockerfile` pins `FROM rocker/verse:4.6.0` (verified by `grep`),
so the disclosed R version does not match the build environment
actually used. (c) the RNG paragraph states
`set.seed(20260411)`, but the actual seed set in the `sim-seed`
chunk (report.Rmd line 351) is `set.seed(20260312)` — the two
numbers differ (verified by direct read of both locations in the
same file). A reproducibility statement with three independently
wrong claims in one short section undermines the paper's own
"Data availability" and "Reproducibility" sections, which a
statistical journal treats as load-bearing for acceptance.
**Remediation:** regenerate this section programmatically from
`sessionInfo()`/`renv::lockfile` output and the actual seed
variable rather than hand-typing it, so it cannot drift from the
code again.

**2.3 The end-of-document "Morris et al. (2019) ADEMP Compliance"
note misrepresents the current state of the pipeline to the
reader.** *Location: report.Rmd, final section (lines 833-851),
citing `docs/morris-audit-2026-04-17.md`.* This note tells the
reader the simulation has three outstanding gaps: seed reseeded on
every scenario call, no Monte Carlo SE on any metric, and
`n_sim = 1000` hardcoded without an MCSE-based derivation. All
three are contradicted by the manuscript's own current Methods
section: the `sim-seed` chunk sets `RNGkind` and the seed exactly
once (report.Rmd lines 344-352, with an explicit code comment
stating the worker "no longer accepts a `seed` argument"); the
paragraph immediately below it derives `n_sim = 1000` from a
target coverage MCSE of 0.7 percentage points; and
`summarise_sim()` (`sim_site_covariate.R` lines 171-209) computes
`mcse_bias`, `mcse_empirical_se`, `mcse_mse`, `mcse_power`,
`mcse_coverage`, and `mcse_convergence` for every model. The
remediation described in the audit file has evidently been carried
out in the code, but the manuscript's own audit summary was never
updated to reflect it, so the paper currently tells a referee it
is "partially compliant" with gaps that no longer exist while
simultaneously exhibiting a different, real gap (2.4 below) that
the stale audit does not mention. **Remediation:** re-run the
Morris audit against the current code, replace the stale summary
with a corrected verdict, and either move the `## Morris et al.
(2019) ADEMP Compliance` subsection out from under the "# References"
heading (see 3.5) or retitle it as its own top-level section (for
example, "ADEMP compliance statement").

**2.4 Computed Monte Carlo standard errors are not shown in the
tables where they are most needed, and a direct check shows the
paper's central fixed-vs-random claim is not distinguishable from
Monte Carlo noise in this simulation.** *Location: report.Rmd,
Tables 1-4 (`tab-equal-no`, `tab-unequal-no`, `tab-equal-int`,
`tab-unequal-int`, lines 430-546); Discussion, first two
paragraphs (lines 700-719).* The Methods section promises
"[performance measures], each reported with its ... Monte Carlo
SE" and the headline paragraph does use `mcse_bias` inline, but
none of the four per-scenario result tables includes any `mcse_*`
column, even though `summarise_sim()` computes them. This matters
because the paper's central claim — that random effects
outperforms fixed effects, "amplified when site sizes are
unequal" — is not visibly supported once Monte Carlo error is
taken into account. I (Rscript, `verify_sim2.R`, `n_sim = 100`,
`site_size = "unequal"`, `trt_site_sd = 0.15`, i.e. the scenario
the Discussion singles out as the case where the advantage should
be "amplified") verified directly that the mean model-based SE for
`site_fixed` and `site_random` in that scenario differ by
approximately 3e-6 (0.083461 vs 0.083463), and the maximum
per-replicate absolute difference across 100 replications is
0.00035. The rendered manuscript's own Table 1 shows the same
pattern to three decimal places (Model SE 0.071 vs 0.071 for
`site_fixed` vs `site_random`; power 0.991 vs 0.991; coverage
0.953 vs 0.953) — in every one of the four tables, the fixed- and
random-effects rows agree to the displayed precision. With
`n_sim = 1000`, the MCSE of a power estimate near 0.99 is on the
order of 0.003; the paper's own numbers do not exhibit a
distinguishable fixed-vs-random advantage anywhere in the four
tables. The Discussion's claim of a precision advantage "amplified
when site sizes are unequal" is not supported by the reported
values; it is asserted from the cited literature rather than
demonstrated by this simulation. **Remediation:** either (a) add
MCSE columns to all four tables (or a companion table of paired
model differences with their MCSEs) and revise the Discussion to
report only differences that clear, say, two MCSEs, or (b)
redesign the scenarios to include the sparse-site regime discussed
in 2.5, where a real fixed/random gap is more likely to emerge.

**2.5 The scenario design does not test the regime the
Introduction identifies as where fixed vs random effects diverge.**
*Location: Methods, "Simulation scenarios" (lines 297-312);
Introduction, "Fixed versus random effects" (lines 169-197), citing
@Kahan2013continuous and @Kahan2015fewevents.* The Introduction
explicitly motivates the fixed-vs-random comparison by noting that
random effects show "notable advantages when the number of
patients per site was small (10 or fewer)" and that a review of
206 published trials found a "median number of events per
centre-treatment arm combination [of] just 3." The simulation,
however, uses `n_sites = 20` with a mean of 20 patients per site
(minimum enforced at 4 after rounding), which is well above the
sparse-site regime the paper cites as practically important. This
is consistent with the near-zero fixed/random SE gap documented in
2.4: an orthogonal, moderately sized design is exactly the
condition under which fixed and random effects are expected to
coincide (a classical ANCOVA orthogonality result — because
treatment is allocated at essentially 1:1 within every site, the
site design matrix is close to orthogonal to treatment, so
site-adjustment method has little effect on the treatment point
estimate; verified directly: in the *equal*-site-size, no-interaction
scenario, the treatment estimate is numerically identical to at
least six decimal places across all four models, `Rscript`
`verify_sim.R`, `n_sim = 20`). A referee would ask why the
simulation does not include a small-site-size scenario (e.g. a
factor for `n_per_site` at 4-8 patients, or a fifth scenario using
the empirical distribution reported by Kahan 2015) to actually
probe the regime motivating the study. **Remediation:** add a
sparse-site scenario (small, fixed sites per arm) crossed with the
existing factors, since this is the setting in which the paper's
own literature review claims the practical stakes are highest.

**2.6 Convergence is defined only as "did not raise an R error,"
which cannot detect the numerical failure modes `lme4` actually
produces.** *Location: `sim_site_covariate.R`, all four `tryCatch`
calls (lines 60-76); Methods, "Performance metrics" (lines 341-342,
"proportion of replications in which the model converged").* Each
model fit is wrapped only in `tryCatch(..., error = function(e)
NULL)`; there is no capture of `warning()`, no check of
`fit@optinfo$conv$opt`, and no `lme4::isSingular()` check. `lme4`
signals convergence problems (failure to converge, singular fit,
boundary estimates) as *warnings*, not errors, in the overwhelming
majority of cases, so this code structurally cannot register the
condition the Methods section says it measures. I verified this
directly: fitting all four models on 100 replications of the
scenario most likely to stress the optimizer (unequal site sizes,
non-zero interaction variance) raised zero warnings and zero
errors (`Rscript verify_sim2.R`), and the manuscript's own tables
report convergence rates of 1.000 (or 0.999) throughout. A
"convergence rate" metric that can only ever equal 1 (barring a
crash) is not informative and should not be presented as a
performance measure distinguishing the four models. **Remediation:**
wrap fits with `withCallingHandlers()` or check
`lme4::isSingular(fit)` and `fit@optinfo$conv$opt != 0` explicitly,
and record these as distinct outcome categories (converged /
singular / non-converged / error).

**2.7 The estimand implemented for the interaction model does not
match the estimand stated in the Methods section.** *Location:
Methods, "ADEMP structure" — Estimand ("the site-size-weighted
average of site-specific effects," line 269); `sim_site_covariate.R`,
`extract_trt()`, `site_x_trt` branch (lines 96-98).* The code
computes averaging weights as `table(dat$site[dat$trt == "1"])`,
i.e. the number of *treated-arm* observations at each site, not
the total site size. Because treatment is allocated at
approximately 1:1 within each site by construction, the two
weighting schemes are numerically close but not identical in the
unequal-site-size, odd-`n`-per-site scenarios (`rep(c(0, 1),
length.out = n)` gives unequal arm sizes whenever `n` is odd).
This is a genuine specification/implementation mismatch, not
merely an approximation, and should be resolved explicitly one way
or the other and the Methods text corrected to match the code (or
vice versa). **Remediation:** either reweight by total site `n`
(matching the stated estimand) or restate the estimand as the
treated-arm-size-weighted average and justify that choice.

**2.8 Test suite provides no coverage of the simulation engine.**
*Location: `inst/tinytest/test_basic.R`.* The entire test file is
`library(tinytest); expect_true(TRUE)` — a placeholder with zero
assertions about `sim_site_covariate()` or `summarise_sim()`. Given
that this package's Suggests field lists `tinytest` as the mandated
test framework and the manuscript's scientific conclusions rest
entirely on this function, the absence of any regression test
(e.g. a small fixed-seed run checked against a stored reference
output, a check that bias is near zero for `true_trt = 0`, or a
unit check of the `mcse_*` formulas in `summarise_sim()`) is a
material reproducibility gap. **Remediation:** add tests that (a)
run `sim_site_covariate()` at small `n_sim` with a fixed seed and
check output shape/column types, (b) verify `summarise_sim()`
against hand-computed values for a toy input, and (c) regression-test
the `site_x_trt` weighted-contrast computation (issue 2.7) against
an independent calculation (e.g. `emmeans`).

## 3. Minor issues

**3.1 Dead code in `extract_trt()`.** *Location:
`sim_site_covariate.R`, lines 106-107.* The line
`contrast <- c(1, rep(0, length(trt_idx) - 1))` is immediately
overwritten by the next line and never used; harmless
(verified: values are recomputed before any use) but should be
removed as vestigial cruft.

**3.2 Fragile positional relabeling of model rows in the display
tables.** *Location: report.Rmd, all four table chunks, e.g. lines
435-438 (`mutate(model = c("No site", "Site (fixed)", ...))`).*
The four-row rename relies on `dplyr::group_by(model) |>
summarise()` alphabetizing `no_site, site_fixed, site_random,
site_x_trt` into exactly the intended display order; this happens
to hold today but is not defended anywhere and would silently
mislabel rows if a model name or count ever changed. Recommend an
explicit `case_when()`/join on the `model` key rather than
positional assignment.

**3.3 Misplaced section nesting.** *Location: report.Rmd, lines
831-838.* The "Morris et al. (2019) ADEMP Compliance" subsection is
nested under the top-level "# References" heading, even though it
has nothing to do with the reference list; with
`number_sections: true` this yields an oddly numbered
sub-heading under "References" in the rendered PDF (verified by
reading the heading hierarchy in the source). See remediation
under 2.3.

**3.4 Unused, topically relevant bibliography entries.** *Location:
`analysis/report/references.bib`.* Seven of the twenty-four bib
entries are never cited in `report.Rmd` (verified by comparing
`grep -oE '@[A-Za-z0-9_]+'` against the bib file's `@`-keys):
`Pocock1975allocation`, `Fleiss1986design`,
`cliftonCorrelationBaselineScore2019`,
`herschtalEffectDichotomizationSkewed2023`,
`pirondiniCovariateAdjustmentCardiovascular2022`,
`wangAnalysisCovarianceRandomized2019`, and
`wangUtilitiesPitfallsStratified2022`. Several are canonical
covariate-adjustment references (Pocock 1975 on stratified
allocation; Fleiss 1986's textbook; Wang et al. 2019 on ANCOVA
precision) that bear directly on the paper's question and were
evidently already gathered but never integrated into the argument.
Either cite them where relevant (see section 5) or prune the bib
file to the works actually used.

**3.5 "n_per_site" is fixed once per scenario, not resampled per
replication, and this is not disclosed.** *Location:
`sim_site_covariate.R`, lines 18-23 (outside the `for (s in
seq_len(n_sim))` loop); report.Rmd, "Simulation scenarios" (lines
297-312).* The unequal site-size vector is drawn once and reused
for all 1000 replications of a scenario, so "unequal site sizes"
describes one fixed, arbitrary site-size pattern rather than a
distribution over designs. This is a defensible ADEMP choice
(treating site sizes as a fixed design factor, analogous to total
sample size) but the manuscript's prose ("exponentially distributed
with mean 20, minimum 4") reads as though sizes vary by
replicate. Clarify in the Methods text.

**3.6 Version metadata is stale relative to the pipeline described.**
*Location: `DESCRIPTION` (`RoxygenNote: 7.2.0`, `Version:
0.0.0.9000`).* Cosmetic but worth tidying before any release tag;
`RoxygenNote` is unrelated to an empty `R/` directory and should
either be removed or the package should gain the documented
functions it purports to version.

## 4. What remains to be done

**(a) Required for correctness**

- Resolve the Type I error gap: add a null-effect (`true_trt = 0`)
  scenario, or remove the Type I error claim from the Abstract
  (2.1).
- Correct the three factual errors in the Reproducibility section
  (package list, R version, seed value) (2.2).
- Update or remove the stale end-of-document ADEMP compliance note
  so it reflects the current code (2.3).
- Fix the estimand mismatch in `extract_trt()`'s `site_x_trt`
  weighting, or correct the stated estimand to match the code
  (2.7).
- Replace the `tryCatch(error = ...)`-only convergence check with
  one that also captures `lme4` convergence/singularity warnings,
  and re-report convergence rates (2.6).

**(b) Required for acceptance**

- Add Monte Carlo SEs to Tables 1-4 (or a paired-difference table)
  and revise the Discussion so that no cross-model claim is stated
  without a demonstrable effect relative to MCSE (2.4).
- Add a sparse-site scenario that actually tests the regime the
  Introduction cites as practically important, or reframe the
  paper's claims to match the moderate-site-size regime that was
  actually simulated (2.5).
- Add non-trivial tests of `sim_site_covariate()` and
  `summarise_sim()` (2.8).
- Reconcile citation of already-gathered ANCOVA/stratification
  literature with the paper's argument, or prune the bibliography
  (3.4).
- Disclose that site-size allocation is fixed once per scenario,
  not resampled per replicate (3.5).

**(c) Desirable polish**

- Remove dead code in `extract_trt()` (3.1).
- Replace positional model relabeling with an explicit join/key
  match (3.2).
- Move the ADEMP compliance note out from under "# References"
  (3.3).
- Tidy `DESCRIPTION` metadata (3.6).

## 5. Recommended framing

**Plausible framings.** (i) *New methodological guidance paper*:
present the simulation as new evidence resolving an open question
about fixed vs. random site effects in general. (ii) *Applied
simulation for the longitudinal-outcome case specifically*:
position the contribution narrowly as extending the existing
fixed-vs-random literature (dominated by cross-sectional binary
and continuous outcomes, per Kahan 2013, 2014, 2015) to a
longitudinal repeated-measures setting with a mixed-model
structure. (iii) *Pedagogical/practice-guidance note*: a shorter
piece aimed at trialists and statisticians that synthesizes
existing evidence (Kahan, Edgar, regulatory guidance) into
concrete recommendations, using a small illustrative simulation
rather than claiming new methodological findings.

**Recommendation: framing (ii), the longitudinal extension, but
scaled back from "new evidence resolving the question" to "a
targeted check of whether the established cross-sectional
fixed/random findings transfer to a longitudinal, repeated-measures
setting."** The existing literature already densely covers the
central claim the manuscript currently leads with — that random
effects perform at least as well as fixed effects, especially with
few patients per site (Kahan2013continuous, across 378 scenarios;
Kahan2014binary for binary outcomes; Edgar 2021 via re-analysis of
a real trial) — and the regulatory recommendation to include
stratification variables is likewise already settled and cited
correctly. What is not already covered, and what this simulation
is actually positioned to speak to, is whether that guidance
continues to hold once outcomes are longitudinal with a
subject-level random intercept layered under the site structure —
a genuine, comparatively narrow gap. However, as documented in
2.4-2.5, the simulation as designed (moderate, near-balanced site
sizes) does not actually probe the sparse-site conditions under
which the cross-sectional literature reports the fixed/random
distinction mattering, so under the current design the paper
cannot yet claim to have extended that literature's practical
conclusions — it has largely reproduced the "site should be
adjusted for at all" result (site-adjusted vs. unadjusted, a
much less novel comparison already well established by
Localio 2001 and Kahan 2013) rather than the "fixed vs. random"
result under conditions that matter.

**Implications of the recommended framing.**

- *Title*: reframe to signal the longitudinal-specific,
  narrower scope, e.g. "Fixed versus random site effects in
  longitudinal multicentre trials: does a sparse-site advantage
  persist?" rather than the current general "Should site be a
  covariate in the analysis model?", which oversells the scope
  given how well-trodden the "adjust vs. don't adjust" question
  already is.
- *Abstract/Introduction*: lead with the longitudinal gap
  specifically (repeated-measures correlation structure layered
  under site structure), not with the general regulatory-guidance
  motivation, which is already answered in the literature and adds
  little novelty; keep the regulatory context as brief background,
  not as the paper's main contribution.
- *Comparators*: add the sparse-site scenario (2.5) so the paper's
  central fixed-vs-random claim is actually tested where the
  literature says it matters; consider adding a Type I error
  scenario (2.1) since this is a routine referee expectation for
  any power/coverage simulation.
- *Target journal*: with the narrower, longitudinal-specific
  framing and the described revisions, this is realistic for a
  simulation-methods-focused outlet such as *Statistics in
  Medicine*, *Pharmaceutical Statistics*, or *Trials* (the latter
  already houses several of the paper's key citations, e.g. Edgar
  2021, Clifton & Clifton 2019); it is not yet ready, on rigor
  grounds (2.1-2.7), for a flagship venue such as *JASA* or
  *Biometrics*.
- *Material to de-emphasize or move to supplement*: the "no site"
  vs. "site-adjusted" comparison (the biggest, least novel effect
  in the paper's own tables — e.g. power 0.973 vs 0.991 in Scenario
  1) should be summarized briefly rather than treated as a
  headline finding; the regulatory-context subsection (currently a
  full subsection under Introduction) could be condensed to a
  paragraph. The currently unused ANCOVA/stratified-randomization
  references (3.4) should either be woven into a short paragraph
  distinguishing this paper's site-specific question from the
  broader baseline-covariate-adjustment literature, or dropped.

## 6. Assessment

**Verdict: major revision.**

The simulation infrastructure is sound in its broad strokes (ADEMP
structure is present, MCSE machinery exists in the code, the
`lme4` model specifications are individually correct as far as
they were inspected), and the topic is a legitimate, currently
active area of methodological interest. However, a referee would
not accept this manuscript in its present form for three
independent reasons, each individually sufficient to require major
revision: (i) the paper's central quantitative claim — a
meaningful fixed-vs-random precision advantage — is not
distinguishable from Monte Carlo noise in the reported tables, and
a direct check confirms the underlying point estimates and
standard errors are numerically almost identical across models in
these scenarios (2.4-2.5); (ii) the Reproducibility section
contains multiple verifiably false statements about the pipeline
that a referee would flag on inspection alone (2.2); and (iii) a
scientifically important condition promised in the Abstract (Type
I error) is never simulated (2.1). None of these requires new
methodology to fix — they require additional simulation factors,
corrected reporting, and tightened claims — but together they are
more than "minor revision" in scope.

## 7. Revision history

- 2026-08-16: Initial referee-grade review. No prior
  `pub_review_whitepaper_*.md` existed in this repository. Nine
  issues identified as major (Type I error scenario missing;
  Reproducibility section factually incorrect on package list, R
  version, and seed value; stale self-reported ADEMP audit;
  Monte Carlo SEs omitted from result tables while the central
  fixed-vs-random claim is not distinguishable from noise; scenario
  design does not probe the sparse-site regime the Introduction
  cites as important; convergence metric structurally cannot
  detect `lme4` warnings; estimand mismatch in the interaction
  model's averaging weights; and a placeholder test suite with no
  real coverage), six as minor (dead code; fragile positional
  table relabeling; misplaced section nesting; unused bibliography
  entries; undisclosed fixed site-size draw across replicates; and
  stale package metadata). Overall assessment: major revision.
