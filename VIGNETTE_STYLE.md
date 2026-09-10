# Vignette prose style

How the prose in `vignettes/` is written, and the `desc:` blocks in
`_pkgdown.yml` with it. It covers register and rhetoric only — for what the
vignettes should contain and how they are grouped, see the "How these vignettes
are organised" section of `vignettes/biomasslmb.Rmd` and the `articles:`
section of `_pkgdown.yml`.

Not shipped in the built package: this file is listed in `.Rbuildignore`.

## The stance

The reader is a colleague who has to make a decision, not a student being
taught and not someone following a recipe. Three consequences run through
everything below:

- **Authority comes from measurement, not from consensus.** Where a claim
  matters, measure it on the packaged data or the facility's own experiments
  and report the number. Cite the literature for mechanism and for methods,
  not as a substitute for evidence you could produce.
- **The subject of an article is a decision, not a procedure.** A vignette that
  only shows what to type has not done its job. Say what the options cost.
- **State the limits of your own evidence before the reader finds them.**
  Sample sizes, confounds, and cases the numbers do not cover go in the text,
  unprompted.

## Rules

### 1. Do not narrate in the first person

Avoid `we`, `our`, `us`. The analysis is described, not performed together.
Where the sentence needs a person, it is the reader: `you`.

> Here, we will QC and filter the PSM level abundances from PD before
> summarising them to protein-level abundances.

becomes

> This vignette works through one typical experiment from end to end: a
> whole-proteome comparison between two conditions, three replicates each.

`we` is acceptable in the rare sentence about the facility as an institution
("the facility runs most of its whole-proteome work on TMT"), and in
`biomasslmb.Rmd` where the offer of help is personal.

### 2. No "Here, we will…" openers

A section opens on its subject, not on an announcement of itself. Delete the
announcement and the sentence usually already exists underneath it.

### 3. Say what it costs

Every option in these vignettes has a price. Naming it is the point of the
prose; without it the reader has a menu rather than a decision.

> Since the same amount of sample was labelled in each case, it's reasonable
> to use 'diff.median' normalisation.

is a decision with its cost left out. Add the condition under which it stops
being reasonable, and what happens then.

### 4. Show the failure mode, not just the correct path

For any step that looks the same whether or not it worked, show what a broken
result looks like and give the check that catches it. This is the whole
subject of the `gotcha_*` articles, but it belongs in the core workflow too.

> These are separate defences that do not always agree with one another, and
> when one of them silently fails the filtering still appears to work.

### 5. Bold the finding, then expand it

Where a paragraph reports a result, put the claim in a short bold sentence and
let the rest of the paragraph support it. This is what makes long analytical
articles navigable without callout boxes.

> **The y axis separates almost nothing.** Enrichment designs sit at a median
> of … on interquartile ranges that overlap across most of their length.

### 6. Inline the qualification; do not box it

Caveats belong in the sentence they qualify, usually behind an em-dash, not in
a callout or a footnote. Qualification that is worth reading is worth reading
in place.

> …on twelve TMT experiments and three LFQ-DDA ones, so read them as an order
> of magnitude rather than as precise figures.

This is what makes the sentences long. That is accepted: mean sentence length
around 24 words, with roughly one sentence in six over 35 words.

### 7. Audit your own numbers in the text

Where a figure or table carries an argument, name the confound and the sample
size in the prose next to it.

> The facility runs most of its whole-proteome work on TMT and most of its
> enrichment work label-free, so acquisition type and experimental design are
> close to collinear here, and every comparison between colours is partly a
> comparison between shapes.

### 8. Say what the article is, and what it is not

Open each vignette by placing it: what it covers, what it assumes has been
read, what it deliberately leaves to another article. Where an article is
reference material rather than a route to follow, say so.

> This article is a map of that path … Apart from the worked example
> immediately below, it contains little analysis itself.

### 9. Cross-link at the decision, not in a list

Link from the sentence where the reader would want the other article, with the
reason attached. Avoid "see also" collections at the foot of a section.

> An enrichment experiment — an IP, BioID or TurboID pulldown — needs different
> reasoning about normalisation and about missing values, and is covered in
> [enrichment designs](interactome_designs.html).

### 10. Hedge less and instruct less

Fewer `may`/`might`/`typically`, and fewer `should`/`must`/`need to`. Both are
replaced by the same thing: state the condition and the cost, and the reader
decides. If a recommendation really is unconditional, give it flatly and say
what it rests on.

### 11. Write as though this were the only version

No `corrected`, `now we`, `instead of X`, `updated to`, `previously`. The
editing history belongs in commit messages and PR descriptions. (This repeats
a global rule, and it is the one most often broken when migrating old prose.)

## Register check

Measured on prose only — YAML, code chunks and inline code stripped. Figures
are per 1000 words of prose, taken from the fifteen vignettes written in this
style. Treat them as a smell test on a finished draft, not as targets to write
against.

| Marker | Target | Old style, for contrast |
|---|---|---|
| `we` / `our` / `us` | ~2 | 19–36 |
| `you` / `your` | ~5 | 1–3 |
| hedges (*may, might, typically, usually*) | ~2 | ~4 |
| obligation (*should, must, need, important*) | ~2.5 | ~7 |
| cost/failure vocabulary | ~2.5 | ~0.9 |
| em-dashes | ~7 | 0–2 |
| bold spans | ~4 | 0–1 |
| "Here, we…" / "We will…" openers | ~0 | ~6 |

A draft that is far outside the first column in more than one or two rows is
usually in the old register rather than being unusual on purpose.

Two rows depend on what kind of article it is rather than on register, and
comparing against the aggregate will mislead. Second-person and bold density
are high in the reference and cautionary articles, which address the reader
directly and report findings, and low in the workflow articles, which walk
through one route: `you` runs at 9-13 per 1000 words in the first group and
0.7-2.9 in the second. Compare a draft against articles of its own kind.

## Files still in the old register

Four vignettes were carried over from `main` and have only partly converged.
Their sentence length and hedging already match; what lags is person, section
openers, and emphasis.

- `LFQ_DDA_Peptide_QC_Summarisation.Rmd`
- `LFQ_DIA_Precursor_QC_Summarisation.Rmd`
- `TMT_PSM_QC_Summarisation.Rmd`
- `summarisation_methods.Rmd`

`LFQ_DDA_Peptide_QC_Summarisation.Rmd` shows the seam clearly: new framing
paragraphs sit on top of inherited body text, and the person switches between
them. Rules 1, 2, 5 and 6 are what these four need; the analysis they contain
is sound and should not be rewritten to suit the prose.
