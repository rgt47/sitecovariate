# Referee report: 'Should site be a covariate in the analysis model?'
*2026-11-18 16:21 PST*

Manuscript: `analysis/report/report.Rmd` (sitecovariate workspace).
Supplement: none. Target venue assumed: *Statistics in Medicine*.
Version reviewed: the committed text at HEAD (`d743ab2`), before the
voice revision of the same date. Quoted wording refers to that
version.

## 1. Summary of the manuscript

The manuscript reports a Monte Carlo simulation that compares four
linear mixed models for a two-arm multicenter trial with a continuous
outcome measured at four visits. The models differ only in how site
enters: omitted, fixed effect, random intercept, or a fixed
treatment-by-site interaction whose site-size-weighted average is
reported. Four primary scenarios cross equal versus unequal site sizes
(about 20 patients per site, 20 sites) with presence versus absence of
a random treatment-by-site effect, at 1,000 replicates. A null-effect
scenario and a sparse-site scenario (4-8 patients per site) are run at
300 replicates from separately seeded streams and loaded from `.rds`
files. The stated contributions are: (a) site adjustment of any kind
outperforms omission when site contributes outcome variability; (b) in
the longitudinal setting, fixed and random site effects are
indistinguishable, so the cross-sectional fixed-versus-random advantage
reported by Kahan and Morris (2013) is not reproduced; (c) an unneeded
interaction term reduces power; and (d) the regulatory recommendation
to include stratification factors is supported.

## 2. Overall assessment and recommendation

**Recommendation: major revision (bordering on reject for a top-tier
methods journal in its current form).**

The question is relevant and the reporting follows ADEMP with Monte
Carlo SEs, which is to the authors' credit. However, the central
'longitudinal' finding (b) is an artifact of the data-generating
mechanism rather than a property of longitudinal data. Treatment is
assigned by deterministic within-site alternation, which makes it
exactly (equal sizes) or nearly (odd sizes) orthogonal to site. Under
that design the fixed, random, and interaction estimators coincide
almost by construction. Moreover, with a treatment effect constant over
time, complete data, and a random subject intercept, the four-visit
mixed model is algebraically equivalent to a cross-sectional analysis
of subject means; the repeated-measures layer contributes nothing to
the comparison. Claim (c) is contradicted by the manuscript's own
tables, and claim (a) rests on a mechanism (an inflated model-based SE
under balanced allocation) that differs from the one the text gives.
In the interaction scenarios the site-adjusted models under-cover,
which the headline findings do not report. None of these defects is
fatal to a revised paper, but each changes what the paper may claim.

## 3. Significance and novelty

The comparison of fixed and random center effects for continuous
outcomes is well covered. Kahan and Morris (2013) ran 378 scenarios
and found random effects at least as good as fixed effects
everywhere, with gains when centers are small. Chu et al. (2011)
compared six methods and found that fixed effects lose precision with
many centers and chance within-center imbalance. Pickering and
Weatherall (2007) addressed small center sizes specifically. The
present paper's claimed novelty is the longitudinal extension. As
noted above, the specific longitudinal DGM used (no treatment-by-time
interaction, compound symmetry, no missing data) reduces exactly to
the cross-sectional problem with residual variance
$\sigma^2_u + \sigma^2_e/T$, so no new information about longitudinal
designs is produced. A genuinely longitudinal contribution would
require at least one of: a treatment-by-time estimand (for example, the
final-visit contrast in an MMRM), site-by-time effects, dropout, or
random allocation within site that permits chance imbalance. In its
current form the paper would fit an applied or tutorial venue (for
example *Contemporary Clinical Trials* or *Pharmaceutical Statistics*)
rather than *Statistics in Medicine*, and only after the corrections
below. The manuscript does not engage with the ICH E9(R1) estimand
framework, which a current referee would expect, given that the choice
of fixed versus random site and the handling of treatment-by-site
heterogeneity determine whether the target is conditional on the
observed sites or marginal over a population of sites.

## 4. Major comments

Evidence labels: **verified** (I executed code or re-derived the
result), **inspected** (I read the source and judged it inconsistent),
**suspected** (plausible but not confirmed).

**M1. Treatment is not randomized within site; the design removes the
mechanism under which fixed and random effects differ (verified).**
Methods, 'Simulation scenarios', states that 'subjects are randomized
1:1 to treatment and control'. The code
(`sim_site_covariate.R`, lines 46-48) assigns
`rep(c(0, 1), length.out = n)` within every site, which is
deterministic alternation. I computed the within-site arm imbalance:
it is exactly 0 in every site under equal sizes and at most 1 under the
sparse design. With exact within-site balance, site is orthogonal to
treatment, and the treatment estimates from the no-site, fixed, random,
and interaction models are identical to machine precision (confirmed
in `sim_typeI_null_n300.rds`, replicate 1, where all four estimates
equal -0.02931353). The Kahan and Morris (2013) and Chu et al. (2011)
advantage for random effects arises from chance within-center
imbalance and from centers with patients in only one arm, which this
DGM cannot produce. The Discussion's statement that the result is 'the
expected behavior under near-orthogonal allocation' is correct, but
the Abstract and Conclusions then present the null fixed-versus-random
result as a finding about longitudinal outcomes. Remedy: allocate by
simple or permuted-block randomization within site (with block sizes
that permit imbalance in small sites), correct the Methods text, and
re-run. Restate claim (b) as conditional on the allocation scheme.

**M2. The longitudinal layer is inert; the model is equivalent to a
cross-sectional analysis of subject means (verified).** The DGM in
Methods has a treatment effect constant over time, a subject random
intercept, a common linear time trend, and no missing data. For a
balanced complete design, generalized least squares under compound
symmetry for a between-subject contrast reduces to OLS on subject
means. I fitted `y ~ trt + time + site + (1|subj)` and
`lm(ybar ~ trt + site)` on subject means for one simulated equal-size
dataset. Both gave estimate 0.2892739 and SE 0.07247337, identical to
the printed precision. The Introduction and Abstract frame the
contribution as testing 'whether the same conclusions hold once
outcomes are longitudinal'. This design cannot answer that question.
Remedy: either reframe the paper as a cross-sectional study with
reduced residual variance, or introduce features that make the
longitudinal structure matter: a treatment-by-time estimand, site-by-
time random slopes, MAR or MNAR dropout, or an unstructured
within-subject covariance analyzed by MMRM.

**M3. The claim that the interaction term 'reduces power' is
contradicted by the manuscript's own output (verified).** The Abstract
(Results), the Discussion (second point), and the Conclusions state
that including the treatment-by-site interaction when none exists
'expends degrees of freedom without benefit and reduces power'. The
sparse-site run (`sim_sparse_n300.rds`) gives power 0.660 for
`site_x_trt` against 0.653 (fixed) and 0.657 (random). My re-run of
Scenarios 1-4 at 1,000 replicates on the pinned stream, under the
current code, gives interaction-model power of 0.991 against 0.991
for fixed (Scenario 1), 0.991 against 0.989 (Scenario 2), 0.983
against 0.982 (Scenario 3), and 0.921 against 0.919 (Scenario 4).
The interaction model is never below the additive models by more
than rounding. The reason is
structural. The p-value uses `df = nrow(dat) - length(fe)`, which is
counted in observations (about 1,600) rather than subjects, so 19
extra parameters cost almost nothing. Also, under exact balance the
site-size-weighted average of site-specific effects equals the
additive-model estimate. Remedy: delete the claim or replace it with
what the tables show. If the authors want to study the interaction
model, they should study the test of interaction itself (type I error
and power of the 19-df interaction F-test or a likelihood ratio test
for a random slope), which the paper never performs despite the
Abstract's reference to 'treatment-by-site interaction testing'.

**M4. Site-adjusted models under-cover when treatment-by-site
heterogeneity is present; the headline findings omit this, and the
estimand is inconsistently defined (verified for under-coverage;
inspected for the estimand).** My pinned-stream re-run of Scenario 3
at 1,000 replicates gives coverage 0.929 (MCSE 0.008) for the fixed
and random site models and 0.928 for the interaction model, against
0.966 for no-site. The model-based SE is 0.071 against an empirical
SE of 0.079. The coverage of the site-adjusted models is about 2.6
MCSEs below nominal. In Scenario 4 coverage is 0.916 (fixed, random)
and 0.915 (interaction), about 3.9 MCSEs below nominal, with model SE
0.084 against empirical SE 0.097. The no-site model attains 0.952
(Scenario 4) and 0.966 (Scenario 3, about 2.8 MCSEs above nominal).
The mechanism is that bias and
coverage are computed against the population effect $\beta = 0.30$,
while the site-adjusted models condition on the realized sites and
ignore the variance of the average $\gamma_i$. However, the ADEMP
'Estimand' bullet says the target is 'the site-size-weighted average of
site-specific effects', which is the conditional quantity. The paper
therefore evaluates a conditional estimator against a marginal target
while declaring a conditional estimand. The no-site model's SE barely
changes between Scenario 1 (0.0856) and Scenario 3 (0.0866), so it
does not model the heterogeneity. My reading (inferred, not tested)
is that its over-stated SE from ignoring a balanced factor (M5)
happens to offset the missing slope variance. Either way, the
no-site model has coverage at or above nominal in the interaction
scenarios while the adjusted models under-cover. This reverses the
'adjusting outperforms omitting' message for coverage. Remedy: define the estimand per ICH E9
(R1) (conditional on sites, or marginal over sites), evaluate each
model against the matching truth, and add a random treatment-by-site
slope model `(1 + trt | site)` as a fifth comparator, since that is
the model that targets the marginal effect.

**M5. The stated mechanism for the no-site model's poorer performance
is not the one operating (verified).** The Introduction and Discussion
(third point) attribute the benefit of adjustment to reduced residual
variance and more precise estimates. In the sparse run the empirical
SE is 0.126 for no-site against 0.123 to 0.124 for the adjusted
models, essentially identical. In the null run the empirical SEs are
identical (0.0635) for all four models. The difference lies in the
model-based SE: 0.152 (no-site) against 0.127, and 0.0856 against
0.0706. That is, the no-site model overstates its SE, giving type I
error 0.017 and conservative coverage (0.98). This is exactly the
Kahan and Morris (2012) result for analyses that ignore a balancing
stratification factor, which the manuscript cites in the Introduction
but does not connect to its results. The Scenario 5 sentence about the
site-adjusted models 'targeting a narrower error stratum' should be
replaced with this explanation. Remedy: rewrite the Discussion's
third point and the Scenario 5 interpretation; report the ratio of
mean model SE to empirical SE as a performance measure.

**M6. The Type I error claim in the Abstract overstates the evidence
(verified).** The Abstract says empirical type I error for the
site-adjusted models 'sits roughly at the nominal 0.05 level'. The
null run gives 0.030 (MCSE 0.0098) for all three site-adjusted models,
about two MCSEs below 0.05, and 0.017 for the no-site model. The
model-based SE (0.0706) exceeds the empirical SE (0.0635, MCSE about
0.0026) by about 11 percent even in the adjusted models. That is
suspicious and may indicate that the SE or degrees of freedom are
conservative, or simply Monte Carlo noise at 300 replicates. Remedy:
run the null scenario at $n_{\mathrm{sim}} \geq 2{,}000$ (type I error
MCSE about 0.005) within the pinned stream and report the result
without the qualifier 'roughly at'.

**M7. Parameter values in the text do not match the code (inspected).**
(a) Methods states $\sigma^2_{s \times t} = 0.15$, but the code passes
`trt_site_sd = 0.15` to `rnorm(..., sd)`, so the variance is 0.0225.
The heterogeneity studied is therefore about one-seventh of the
stated value. (b) Scenario 5 states $\sigma^2_s = 0.50$; the code uses
`site_sd = 0.50`, a variance of 0.25 (consistent with the Methods
value 0.25). (c) The symbol $\delta$ is the time slope in the model
equation but the overall treatment effect in the ADEMP 'Estimand'
bullet; the treatment effect is $\beta$ elsewhere. (d) 'Performance
metrics' says 'For each model across 1,000 replications', which is
false for Scenarios 5-6. (e) The 'Site x Trt' model in the Methods
list is described as an analysis model but its reported estimate is a
post-hoc weighted contrast; the weighting (by site size including
both arms) should be stated in Methods, not only in a code comment.
Remedy: correct each item; decide whether 0.15 was intended as SD or
variance and re-run if the latter.

**M8. Reproducibility: the rendered PDF does not correspond to the
source, and Scenarios 5-6 have no committed driver (verified).** The
committed `report.tex`/`report.pdf` predate the remediation: their
tables carry an 'N conv.' column, no MCSE, and no Tables 5-6. No
knitr cache is committed, so Scenarios 1-4 can be reproduced only by a
full re-run. (I did re-run them. For the no-site, fixed, and random
models, bias, power, and coverage agree with the stale render to
three decimals. The interaction row differs slightly in Scenarios 2
and 4 because its weighting changed after that render.) The two
`.rds` files were generated outside the
manuscript by code that is not in the repository. I was able to
reproduce replicate 1 of the null run with
`RNGkind("L'Ecuyer-CMRG"); set.seed(20260816)` and the shipped
function, which supports the stated provenance, but the driver should
be committed. The Data Availability statement refers to 'driver
scripts' (plural); only `sim_site_covariate.R` exists. Remedy: fold
both scenarios into the pinned stream at full size (the TODO comments
already give the commands), commit a driver, and re-render.

## 5. Minor comments

1. Abstract, Methods: the Abstract mentions 'converged/singular/
   nonconverged/error' classification, but no table reports singular
   or nonconverged rates (all four statuses are computed in
   `summarise_sim()` but only `convergence_rate` is displayed). My
   re-run shows a convergence rate of 0.989 to 0.994 for the
   random-site model with a singular rate of 0, so 0.6 to 1.1 percent
   of its fits are classed as nonconverged. These fits are still
   included in bias, power, and coverage (via `n_fit`). The text
   should say so, and the nonconverged rate should be tabulated.
2. Results, 'Headline findings': power averaged across scenarios with
   and without interaction is not an interpretable quantity. Report
   per-scenario values or drop the sentence.
3. Results, Scenario 6: '@Kahan2015fewevents: a median of 3 events per
   centre-treatment arm' concerns binary events, not patients. Using it
   to motivate 4-8 patients per site for a continuous outcome
   conflates events with sample size. Rephrase.
4. Introduction, 'Fixed versus random effects': the sentence 'Random
   effects performed as well as or better than fixed effects in all
   scenarios, with notable advantages when ... there was imbalance
   between treatment arms within sites' is consistent with Kahan and
   Morris (2013) as summarized in its abstract. It strengthens M1:
   the present DGM excludes that imbalance.
5. Introduction: 'Patients within a site tend to be more similar ...'
   and 'within-site correlation violates the independence assumption'
   are reasonable but uncited beyond Localio et al.; acceptable.
6. Methods, 'Performance metrics': coverage uses a normal 1.96
   quantile, while p-values use a $t$ distribution with observation-
   level df. Make them consistent. Observation-level df
   (`nrow(dat) - length(fe)`) are anti-conservative for a
   between-subject contrast; Kenward-Roger or Satterthwaite df
   (`lmerTest`, `pbkrtest`) are standard and matter for the sparse
   design with few subjects per site.
7. Results text retains process residue: 'whitepaper issue 2.1',
   'whitepaper issue 2.5', 'see TODO in the source chunk', 'This
   section was previously nested under the References heading', the
   audit history in the ADEMP compliance section, and code comments
   such as 'Whitepaper issue 2.4'. These are internal notes and should
   not appear in a submitted manuscript.
8. Scenario 5 opening quotes the Abstract ('consequences for power,
   type I error') to justify a section. A manuscript should not argue
   with its own abstract; state the aim of the scenario directly.
9. Figures 1-3 cover only Scenarios 1-4 and plot raw boxplots of bias
   with no MCSE. Consider a single nested-loop or lollipop plot with
   MCSE intervals across all six scenarios.
10. Spelling: 'randomised' (Abstract) and 'centre' (several places)
    are British; the CSL is *Statistics in Medicine*, which accepts
    either, but the paper should be consistent.
11. The 'Abstract' is a numbered top-level section in the PDF (and in
    the table of contents); use the YAML `abstract:` field.
12. `Fleiss1986design` is typed `@article` for a book; `R-lme4` should
    cite Bates et al. (2015, *Journal of Statistical Software* 67(1))
    rather than a manual with a hand-typed version; `ICHE9guideline`
    author was the 'International Conference on Harmonisation' in 1998.
    Three bib entries are unused (Clifton 2019, Herschtal 2023,
    Pirondini 2022).
13. Discussion first point: 'no cross-model gap in this simulation
    clears two MCSEs' is true for fixed versus random but false for
    no-site versus adjusted (power 0.473 versus 0.653 in the sparse
    run, a gap of more than four combined MCSEs). The sentence should
    be scoped.
14. The unequal-size draw (one draw per scenario) is disclosed, but
    with a single draw of 20 exponential sizes the results depend on
    that one configuration. Report the realized site sizes, or
    resample per replicate.

## 6. Missing and questionable references

| Claim location | Issue | Suggested source |
|---|---|---|
| Introduction, fixed vs random | Directly comparable simulation not cited; it identifies within-center imbalance as the mechanism | Comparing methods to estimate treatment effects on a continuous outcome in multicentre randomized controlled trials: a simulation study. *BMC Medical Research Methodology* 2011;11:21. doi:10.1186/1471-2288-11-21 (title, journal, and DOI from search results; author list recalled as Chu R, Thabane L, et al. and not verified) |
| Introduction, small centers | Small-center continuous-outcome work not cited | Pickering RM, Weatherall M. The analysis of continuous outcomes in multi-centre trials with small centre sizes. *Statistics in Medicine* 2007 (recalled from memory; volume, pages, and DOI not verified) |
| Methods, ADEMP estimand | No estimand framework cited | ICH E9(R1) Addendum on estimands and sensitivity analysis in clinical trials, 2019 |
| Methods, analysis models | Small-sample df not addressed | Kenward MG, Roger JH. Small sample inference for fixed effects from restricted maximum likelihood. *Biometrics* 1997;53(3):983-997 |
| Discussion, third point | Mechanism should cite the balancing-factor result | Kahan BC, Morris TP. Improper analysis of trials randomised using stratified blocks or minimisation. *Statistics in Medicine* 2012;31(4):328-340 (already in the bibliography) |
| Results, Scenario 6 | Events per center cited for patient counts | Kahan and Harhay (2015) is about binary events; cite a source on patients per center or rephrase |
| `@edgarIncludingRandomCentre2021` | Cited for continuous-outcome precision claims and 'underpowered interaction tests'; the paper is a re-analysis of CRASH-2 with a binary outcome | Check that each attached claim appears in Edgar et al.; restrict to what it reports |
| Methods, analysis models | `lme4` cited as manual | Bates D, Maechler M, Bolker B, Walker S. Fitting linear mixed-effects models using lme4. *Journal of Statistical Software* 2015;67(1):1-48 |

## 7. Suggestions for strengthening the paper

1. Replace alternation with permuted-block randomization within site
   and include a scenario with unequal allocation across sites, so that
   the fixed-versus-random comparison has a mechanism to detect.
2. Add a random treatment-by-site slope model, and evaluate each model
   against its own estimand (conditional and marginal).
3. Make the longitudinal structure matter: a treatment-by-time
   estimand at the final visit analyzed by MMRM, site-by-time
   variation, and MAR dropout.
4. Vary the number of sites (for example 5, 10, 40), since the
   fixed-versus-random distinction depends more on the number of
   sites than on their size.
5. Use Kenward-Roger df throughout, and report the ratio of mean
   model SE to empirical SE.
6. Present results across all six scenarios in one nested-loop plot
   with MCSE intervals.
7. Commit a single driver script and a cached-results file per
   scenario, and re-render so the PDF matches the source.

## 8. Scope of this review

**Verified (executed):**

- Summaries of `sim_typeI_null_n300.rds` and `sim_sparse_n300.rds`
  with the shipped `summarise_sim()` (numbers quoted in M3, M5, M6).
- Replicate 1 of the null run reproduces exactly under
  `RNGkind("L'Ecuyer-CMRG"); set.seed(20260816)` but not under
  Mersenne-Twister.
- Within-site arm imbalance under the shipped allocation code (0 for
  equal sizes; at most 1 for sparse).
- Equivalence of the four-visit LMM and the subject-mean OLS for one
  simulated dataset (M2).
- 150-replicate exploratory runs of Scenario 3 and of a sparse
  scenario with interaction under current code. These agreed in
  direction with the full re-run below and are not quoted elsewhere.
- Re-run of Scenarios 1-4 at 1,000 replicates each on the pinned
  stream (`RNGkind("L'Ecuyer-CMRG"); set.seed(20260312)`, scenarios
  in manuscript order, state carried between two R sessions by saving
  `.Random.seed`), using the current shipped code. Results, as
  power / coverage for no-site, fixed, random, interaction:
  Scenario 1, 0.973/0.982, 0.991/0.953, 0.991/0.953, 0.991/0.953;
  Scenario 2, 0.967/0.983, 0.989/0.943, 0.989/0.943, 0.991/0.941;
  Scenario 3, 0.954/0.966, 0.982/0.929, 0.982/0.929, 0.983/0.928;
  Scenario 4, 0.858/0.952, 0.919/0.916, 0.920/0.916, 0.921/0.915.
  These agree with the stale `report.tex` to three decimals except
  for the interaction model, whose weighting changed since that
  render. Absolute bias is at most 0.002 throughout.
  I assumed that no RNG draws occur between the `sim-seed` chunk and
  the first scenario in a knitr render; the agreement with the stale
  render supports that assumption.

**Inspected (read, not executed):** all prose, every code chunk,
`sim_site_covariate.R`, `references.bib`, the stale `report.tex`, and
the remediation log `docs/pub_review_remediation_2026-08-16.md`.

**Not checked:** a full knitr render of the manuscript (not
performed, per instructions); whether knitr caching preserves the sequential RNG stream across
cached chunks (suspected risk); full text of the regulatory documents
(the ICH E9 quotation and the FDA 2023 and EMA 2015 paraphrases were
not checked against the source); full text of Kahan and Morris (2013)
beyond its abstract, of Edgar et al. (2021), and of Kahan and Harhay
(2015); figures (not rendered).

**Literature searches run (WebSearch):** 'Kahan Morris 2013 analysis
of multicentre trials continuous outcomes centre effects fixed random
simulation results'; 'multicentre trial longitudinal repeated measures
centre effect fixed versus random mixed model simulation'; 'treatment
effect heterogeneity across centres estimand multicentre trial random
treatment-by-centre interaction coverage'. These confirmed the 378
scenarios and the main conclusion of Kahan and Morris (2013) and
located Chu et al. (2011) and Pickering and Weatherall (2007). I found
no prior simulation specific to longitudinal multicenter outcomes in
these searches, so the novelty gap is real in principle; the problem
is that the present DGM does not exploit it.
