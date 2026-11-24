# Referee remediation log: report.Rmd
*2026-11-19 00:16 PST*

This log covers remediation of `referee-report-2026-11-18.md`, followed
by an abstract voice pass. It builds on the voice revision of
2026-11-18 and on `docs/pub_review_remediation_2026-08-16.md`. Nothing
was staged, committed, pushed, or rendered.

## Design change underlying most fixes

All seven scenarios are now produced by a committed driver,
`analysis/scripts/run_scenarios.R`. The driver runs one R process per
scenario and writes `analysis/data/derived_data/sim_<id>.rds`. The
manuscript only reads these files. It runs no simulation and has no
knitr cache. RNG: L'Ecuyer-CMRG with master seed 20260312. Scenario k
uses substream k (`parallel::nextRNGStream()` applied k times). This
replaces the single sequential stream, because sequential runs would
have exceeded the 10-minute limit per R run. Wall-clock times with all
eight runs in parallel were 5.5 to 9.9 minutes each.

Engine changes (`analysis/scripts/sim_site_covariate.R`):

- A fifth model, `site_slope`: `y ~ trt + time + (1 + trt | site) +
  (1 | subj)`.
- An `allocation` argument ('alternate', the default, or 'block'),
  implemented in a new `allocate_arms()` helper. 'alternate' draws no
  random numbers, so earlier streams are unchanged.
- A per-replicate `cond_trt` column holding the conditional truth
  (beta plus the size-weighted mean of gamma_i).
- A realized `n_per_site` attribute on the result.
- A guard that returns NA for the interaction model if a column is
  dropped for rank deficiency.
- `summarise_sim()` now also returns `se_ratio`, `coverage_cond`, and
  `mcse_coverage_cond`.

Regression anchor: replicate 1 of the old null run (seed 20260816)
still gives -0.02931353 for all four original models. The first three
replicates of both old `_n300` files reproduce exactly under the new
engine. The `allocate_arms()` refactor reproduces the first three
replicates of `sim_sparse_block.rds` exactly. [verified]

## Fixed

- **M1 (alternation, not randomization).** Methods now has a
  'Treatment allocation' subsection. It states that allocation is
  deterministic alternation within site, which gives exact or near
  balance, and calls this a limit of the design. The false
  'randomized 1:1' sentence is gone. Scenario 7 is new: sparse sites
  with permuted blocks of 4, drawn afresh each replicate. Fixed and
  random power were 0.626 and 0.633 (MCSE 0.015). The fixed-versus-
  random claim is now stated as conditional on near-balanced
  allocation, in the Discussion and the Abstract. Simple randomization
  was not run (see Deferred). [verified]
- **M2 (inert longitudinal layer).** I reproduced the equivalence on
  one simulated dataset. The site-fixed LMM and `lm(ybar ~ trt +
  site)` on subject means gave the same estimate (0.2553992) and SE
  (0.06917704). The equivalence is now stated in 'Present study', in
  the DGM subsection, at the start of the Discussion, and in the
  Abstract. The paper no longer claims a longitudinal contribution.
  [verified]
- **M3 ('interaction reduces power').** The claim is removed from the
  Abstract, the Discussion, and the Conclusions. The new results show
  the interaction model's power minus the fixed model's power ranging
  from -0.002 to 0.010. The text gives the two structural reasons
  (observation-level df; balance). It also states that the interaction
  test itself was not studied. [verified]
- **M4 (under-coverage; estimand).** The ADEMP estimand bullet now
  defines a marginal estimand (beta) and a conditional one (beta_c).
  Coverage is scored against both. The random-slope model was added
  as a fifth comparator. Results:
  - Marginal coverage for fixed, random, and interaction: 0.922,
    0.922, 0.920 in Scenario 3 and 0.898, 0.898, 0.895 in Scenario 4.
  - Conditional coverage (fixed): 0.942 and 0.950.
  - Random-slope marginal coverage: 0.949 and 0.939.

  Under-coverage is reported in the headline findings, in the
  Scenario 3-4 text, in the Discussion (fourth point), and in the
  Abstract. [verified]
- **M5 (mechanism).** 'SE ratio' is added to every table. The no-site
  SE ratio was 1.17 to 1.23 in the scenarios without heterogeneity,
  while its empirical SE was close to that of the adjusted models.
  The 'narrower error stratum' sentence is replaced, and the third
  Discussion point is rewritten around Kahan and Morris (2012).
  [verified]
- **M6 (type I error).** The null scenario was rerun at
  n_sim = 2,000 (two substreams of 1,000). Type I error was 0.049
  (fixed), 0.049 (random), and 0.051 (interaction), with MCSE 0.005.
  No-site was 0.018 and random-slope 0.036. The 'roughly at' wording
  is gone. [verified]
- **M7 (parameter text).**
  - (a) The text now says SD 0.15 (variance 0.0225), matching the
    code. I did not rerun at variance 0.15.
  - (b) The Scenario 5 sigma^2_s = 0.50 sentence is removed. Methods
    states sigma^2_s = 0.25.
  - (c) The time slope is renamed tau. The estimand is beta.
  - (d) Replicate counts are stated per scenario.
  - (e) The weighting of the interaction-model estimator is described
    in Methods.

  [applied, unverified; code values checked by inspection]
- **M8 (reproducibility).** There is now a committed driver, one rds
  per scenario, and no knitr cache dependence. The Reproducibility,
  Data availability, and ADEMP compliance sections are rewritten to
  match. The plural 'driver scripts' is corrected. [verified: the
  driver runs end to end; purl plus execution of every chunk and
  every inline R expression succeeds]
- **Minor 1.** A new 'Fit status' table covers every scenario and
  model. Methods states that singular and nonconverged fits are
  included in all measures. The random-slope model was singular in
  36 to 64 percent of fits. [verified]
- **Minor 2.** The averaged-power sentence is removed. [applied, unverified]
- **Minor 3.** The events-versus-patients wording for Kahan and Harhay
  (2015) is rephrased in Scenarios 6-7. [applied, unverified]
- **Minor 6.** The df and 1.96 inconsistency is disclosed, with
  Kenward and Roger cited. The code is unchanged (see Deferred).
  [applied, unverified]
- **Minor 7.** Process residue is removed from the prose and the Rmd
  setup chunk: whitepaper issue numbers, TODOs, the 'previously
  nested' sentence, and the audit history. [verified by grep]
- **Minor 8.** The Scenario 5 opening that argued with the Abstract is
  removed. [applied, unverified]
- **Minor 10.** 'randomised' and 'centre' are gone, except in the
  verbatim ICH E9 quotation and in bibliography keys and titles.
  [verified by grep]
- **Minor 11.** The Abstract has moved to the YAML `abstract:` field.
  [verified: yaml.load parses]
- **Minor 12, references.**
  - `Fleiss1986design` is now `@book`.
  - `ICHE9guideline` author is now 'International Conference on
    Harmonisation'.
  - `R-lme4` now cites Bates et al. 2015, JSS 67(1). DOI verified on
    Crossref; pages 1-48 not returned by Crossref.
- **Minor 13.** The 'no gap clears two MCSEs' sentence is replaced by
  scoped statements. [applied, unverified]
- **Minor 14.** The realized site sizes for Scenarios 2, 4, 6, and 7
  are printed in Methods through inline R. [verified]
- **Section 6, references.**
  - Added `Chu2011continuous` (DOI 10.1186/1471-2288-11-21) and
    `Pickering2007smallcentres` (10.1002/sim.3068). Both are cited only
    for what their titles state.
  - Added `Kenward1997small` (10.2307/2533558). Crossref confirms the
    first page, 983. The end page, 997, is not returned.
  - Added `ICHE9R1addendum`, a guideline with no DOI and no URL.
    [applied, unverified: no Crossref record exists to check]
  - Kahan and Morris (2012) is now cited for the mechanism.
  - `edgarIncludingRandomCentre2021` is removed from the
    continuous-outcome precision claim, the 'underpowered interaction'
    claim, and the fixed-versus-random advantage. It is now cited only
    for the distinction between outcome variation and treatment-effect
    variation, from its abstract.

  All DOIs were checked against Crossref on 2026-11-19. [verified]

## Deferred

- **Simple randomization within site.** Single-arm sites make the
  interaction model rank-deficient. The new guard returns NA, but the
  scenario was not designed or run. It is needed for a full test of
  M1 (Suggestion 1).
- **A truly longitudinal design** (treatment-by-time estimand with
  MMRM, site-by-time variation, dropout). This is an author decision
  on scope (M2, Suggestion 3).
- **Kenward-Roger or Satterthwaite df.** The code change would alter
  every p-value and need a full rerun (Minor 6, Suggestion 5).
- **Varying the number of sites** (Suggestion 4).
- **Interaction F-test or LRT for random slope** (M3 remedy). Not
  implemented.
- **Variance 0.15 for the heterogeneity term.** If the author meant
  variance rather than SD, rerun with `trt_site_sd = sqrt(0.15)` in
  `run_scenarios.R`.
- **Nested-loop figure across all scenarios** (Minor 9). The figures
  still cover Scenarios 1-4. They were updated to five models.
- **Unused bib entries** (Clifton, Herschtal, Pirondini). Left in
  place, since pandoc ignores them.
- **Render.** The PDF and tex files are stale. Run
  `bash tools/render.sh analysis/report/report.Rmd`.
- **Retiring the old files.** `sim_typeI_null_n300.rds` and
  `sim_sparse_n300.rds` are no longer read. Delete them at the
  author's discretion.

## Illustrative values for author confirmation

- Permuted-block size of 4 in Scenario 7 (labeled illustrative in
  Methods).
- 4 to 8 patients per site in the sparse scenarios (pre-existing; now
  labeled illustrative).
- n_sim = 2,000 for the null scenario, split into two substreams.

## New issues found while fixing

- Each scenario draws its own site sizes. Scenario 2 enrolled 318
  subjects and Scenario 4 enrolled 473, and Scenarios 6 and 7
  enrolled 126 and 118. Paired scenarios are therefore not
  one-factor contrasts. This is now disclosed in Methods and the
  Discussion. Deferred: pass one fixed size vector to paired
  scenarios.
- The Scenario 1-4 numbers no longer match the referee's rerun.
  The design moved to per-scenario substreams, and the referee's
  pinned sequential stream is no longer used. For the four original
  models, the results agree with the referee's within Monte Carlo
  error. Scenario 4 differs most, because of the new size draw.
- The simulations ran on host R 4.6.0 with lme4 2.0.1. renv was
  skipped by the project `.Rprofile`, so the runs did not use the
  container.

- The random-slope model is singular in 36 to 64 percent of fits. Its
  null type I error (0.036) is conservative.
- In Scenario 4 the no-site model's marginal coverage is 0.939, with
  an SE ratio of 0.976. The 'offset' explanation in the text is
  hedged ('It appears that'), and it is inferred rather than tested.
- Under 8-way parallel load, single runs took up to 9.9 minutes, close
  to the 10-minute limit. Serial runs would be faster per scenario.
- The Abstract's numbers (0.018; 'about 20 percent') are typed by hand
  because they sit in YAML. They match the rds files as of this run.

## Abstract voice pass

- Word count: 513 -> 328. The abstract moved to YAML `abstract: |`,
  keeping the Background/Methods/Results/Conclusions structure.
- Factual changes:
  1. Four models became five, with the random-slope model added (M4).
  2. 'Longitudinal' framing became an explicit statement that the
     analysis is equivalent to one of subject means (M2).
  3. The allocation scheme is now stated, including the permuted-
     block scenario (M1).
  4. 'Six scenarios' became seven, with the replicate counts changed
     from 1,000/300 to 1,000/2,000 (M6, M8).
  5. The claim that the interaction term 'reduces power' is deleted
     (M3).
  6. Type I error 'roughly at 0.05' became 0.018 (no site) against
     about 0.05 for the adjusted models (M6).
  7. Under-coverage of the marginal effect, and the random-slope
     remedy, are added (M4).
  8. The mechanism is now an inflated model-based SE (M5).
  9. Citations of Localio and of Edgar et al. for the
     fixed-versus-random advantage are dropped.
  10. The conclusions are restated as conditional on near-balanced
      allocation, with a recommendation to state the estimand first.
- No sentence exceeds 25 words. Subordinate clauses were reduced, and
  the two when-clauses in the Conclusions were split. A few short
  clauses remain ('so', 'if'). 'randomised' became 'randomized'.
- The no-site SE claim ('about 20 percent') is limited to the
  scenarios without heterogeneity. The type I error of about 0.05 is
  attributed only to the fixed, random-intercept, and interaction
  models. The random-slope model gave 0.036.

## Checks

- tinytest: 234 results, all pass.
- `knitr::purl` succeeds. Executing every chunk, with figures
  printed to a null device, succeeds against the new rds files. So
  does evaluating every inline R expression.
- Every `@key` resolves in `references.bib`. The one non-citation
  match, `@optinfo`, is a code slot name.
- Em dashes: 0 before and 0 after.
