# Voice and humanizer revision log: report.Rmd
*2026-11-18 16:21 PST*

Manuscript: `analysis/report/report.Rmd`. Supplement: none.
Pass 2 applied `~/.claude/commands/use_rgtvoice.md`; Pass 3 applied
the humanizer skill in File mode to the Pass 2 text. Both passes were
carried out together, passage by passage, with the Edit tool. The file
was not staged, committed, rendered, or pushed.

## Scope of edits

There were 34 Edit operations on prose, touching about 90 sentences
in the Introduction, Methods, Results prose, Discussion, Future
research, Reproducibility, and ADEMP compliance sections. The
Abstract, YAML header, all code chunks and chunk options, display and
inline math, inline R expressions, citation keys, table and figure
content, and the section structure are unchanged.

## Kinds of changes (Pass 2, voice)

1. **Sentence splitting.** Sentences over 25 words, or carrying
   subordinate clauses, were split into two or three sentences joined
   by a plain connective ('That is,', 'As a result,', 'By contrast,',
   'In that case,', 'therefore'). Examples: the four-strategy list in
   'Present study'; the parameter definitions after the model
   equation; the first and third Discussion points; the
   Reproducibility RNG paragraph.
2. **Common words.** 'demonstrated' became 'showed'; 'primary
   argument' became 'main argument'; 'require' became 'need' in
   places; 'consumes degrees of freedom proportional to' became 'uses
   degrees of freedom in proportion to'; 'facilitate' became
   'support'.
3. **Voice markers.** 'On the plus side ... On the minus side' for the
   fixed-effects trade-off; 'merit' for 'deserve'; 'would be useful'
   for 'would be valuable'; 'consider' as a signpost in the first
   Discussion point. The existing 'We note that', 'It appears that',
   'Unfortunately,' and 'We shall' were kept.
4. **US spelling.** 'centre' became 'center' in prose (Introduction,
   Agresti and Islam sentences, Kahan 2015 sentences, Scenario 6).
5. **Quotation marks.** Double quotes in prose became single quotes
   (global rule), including around 'No site', 'power', 'unequal site
   sizes', the Abstract quotations, the 'References' heading
   reference, and the audit verdict.

## Kinds of changes (Pass 3, humanizer)

1. **Not-X-but-Y.** 'The question is not merely academic' became 'The
   answer has consequences in practice.' 'It is not between the fixed-
   and random-effects models specifically, whose ...' became 'By
   contrast, the power estimates of the fixed- and random-effects
   models ... differ by less than one Monte Carlo SE.' 'rather than a
   novel finding of this study' became 'and it is not a novel finding
   of this study' (kept as a scope statement, since it carries
   information).
2. **Staged run-ups.** 'What the simulation does show consistently
   ... is the much larger ... effect' became a direct statement.
   'The practical relevance of this comparison is underscored by'
   became '@Kahan2015fewevents show the practical relevance of this
   comparison.'
3. **AI-typical words.** 'key question' became 'practical question';
   'underscored' removed; 'illuminate' became 'bear on'; 'valuable'
   became 'useful'.
4. **Decorative bold.** The inline '**not**' in the Discussion was
   removed.
5. **Dashes.** None present before or after.

## Passages deliberately left alone

- **Abstract (lines 77-146):** left byte-identical per instructions,
  including the British 'randomised'.
- **ICH E9 quotation:** 'randomisation' retained as verbatim quote.
- **Claims the referee report flags as wrong** were reworded only and
  never strengthened: 'subjects are randomized 1:1' (M1); 'This
  introduces an additional layer of within-subject correlation' (M2);
  'including an interaction term ... reduces power' (M3); the
  'narrower error stratum' explanation and 'inflated residual
  variance' (M5); 'in the vicinity of the nominal 0.05 level' (M6);
  $\sigma^2_s = 0.50$ in Scenario 5 and $\delta$ in the Estimand
  bullet (M7); 'For each model across 1,000 replications' (M7).
- **Process residue:** 'whitepaper issue 2.1/2.2/2.5/3.3', 'see TODO in
  the source chunk', and 'This section was previously nested under
  the References heading' would normally be cut under humanizer §25.
  They carry information about the manuscript's state, and the hard
  constraint forbids deleting information-bearing sentences, so they
  were kept and raised as minor comment 7 in the referee report.
- **Bold run-in labels:** the ADEMP bullets ('Aims.', 'Estimand.'),
  performance-metric labels ('Bias', 'Power'), Future-research item
  labels, 'Random-number generator.', 'Data and code availability.',
  'Verdict:', and the audit-gap labels were left. They act as
  structural run-in heads in a numbered or bulleted enumeration, and
  removing them would alter list structure.
- **LaTeX lists** (`itemize` for scenarios, `enumerate` for model
  formulas) were left untouched, as they are part of the scenario and
  model definitions.
- **Headings:** all are already in sentence case (after-colon
  capitals in 'Scenario N: ...' and 'Type I' are conventional), so
  none were changed.
- **Number words versus digits** ('one Monte Carlo SE', 'two MCSEs')
  and ranges ('1-4', '4-8') were kept as written.
- **No sentence merged.** No sentence was deleted as a literal repeat.

## Verification

Script: `verify.py` in the private scratch directory, comparing the
edited file against `git show HEAD:analysis/report/report.Rmd`.

| Check | Result |
|---|---|
| YAML header identical | PASS |
| Code-chunk lines (multiset) | PASS |
| Inline R expressions (whitespace-normalized multiset) | PASS |
| Numbers in prose (multiset, trailing `.,:` ignored) | PASS |
| `@` citation keys (multiset, trailing punctuation ignored) | PASS |
| Inline and display math spans (whitespace-normalized) | PASS |
| Heading list (case-insensitive, in order) | PASS |
| Abstract block byte-identical | PASS |
| Em or en dashes in file | 0 |

A Scenario 5 split ('In the 'No site' model, that variance is left
in the residual.') turned an implied point into a standalone claim
on a passage the referee report flags (M5). It was folded back into
the preceding sentence, and the checks were re-run with the same
result.

The 'no subordinate clauses' rule was applied to reduce subordinate
clauses, but it has not eliminated them. Several short when-, if-,
and which-clauses remain where splitting would have produced choppy
prose or separated a qualifier from its claim. Only the 25-word
ceiling was checked mechanically.

One intermediate run failed the number check, because 'Scenario 6
likewise shows' added a '6'. It was reworded to 'This scenario likewise
shows', and the check then passed.

The sentence-length scan (outside the abstract) flags three items.
Two are splitter artifacts: a sentence ending in a closing single
quote ('... type I error.' and '... an R error.') is joined to the
next sentence. Each real sentence is under 25 words. The third is
the 'Remaining gap' paragraph, at 27 words including its five-word
bold run-in label (22 without it). No sentence of revised prose
exceeds 25 words.
