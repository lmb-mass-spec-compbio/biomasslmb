# Summarise the PTM groups formed at a given min_candidate_prob

Runs
[`add_ambiguous_ptm_group_rowdata()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_ambiguous_ptm_group_rowdata.md)
at one `min_candidate_prob` and reports the size and span of the pooled
groups it forms, alongside how often those groups exclude the true site.
Sweeping it over a range of thresholds is how `min_candidate_prob` is
chosen: raising it buys smaller, more interpretable groups and pays in
coverage.

The miss rate columns read directly off the probabilities. Because they
sum to one per modification, the mass a threshold leaves behind is the
chance that the group it forms does not contain the real residue, so
`pc_miss_over_5` is the percentage of pooled peptides whose group has
more than a 5\\

## Usage

``` r
summarise_ptm_groups(
  obj,
  candidates,
  min_candidate_prob,
  min_prob = 0.501,
  max_group_span = 50,
  ...
)
```

## Arguments

- obj:

  `SummarizedExperiment`. As passed to
  [`add_ambiguous_ptm_group_rowdata()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_ambiguous_ptm_group_rowdata.md)

- candidates:

  `data.frame` of candidate residues from one of the
  `parse_ptm_candidates_*()` functions

- min_candidate_prob:

  `numeric` Threshold to summarise

- min_prob:

  `numeric` Probability at or above which a residue counts as localised

- max_group_span:

  `numeric` Widest span in residues a group may cover

- ...:

  Further arguments for
  [`add_ambiguous_ptm_group_rowdata()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_ambiguous_ptm_group_rowdata.md)

## Value

One row `data.frame` of group size, span and coverage statistics
