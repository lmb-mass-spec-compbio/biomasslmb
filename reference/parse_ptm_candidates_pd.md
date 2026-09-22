# Parse Proteome Discoverer ptmRS site probabilities into candidate residues

Extracts every candidate residue for a modification from the ptmRS
output, rather than only the residues that are confidently localised.

This reads `ptmRS.Phospho.Site.Probabilities`, **not** the
`ptmRS.Best.Site.Probabilities` column that
[`parse_PTM_scores_pd()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/parse_PTM_scores_pd.md)
uses. The Best Site column reports only the winning isoform's sites, so
the alternatives a group is built from are absent from it. Pointing this
function at the Best Site column fails silently, producing groups that
contain only the sites ptmRS already preferred.

    Best Site Probabilities   : S3(Phospho): 39.61; S4(Phospho): 39.61
    Phospho Site Probabilities: S(1): 1.3; S(3): 39.6; S(4): 39.6; S(6): 6.5

ptmRS reports percentages; they are rescaled to a 0 to 1 probability so
that `min_prob` and `min_candidate_prob` in
[`add_ambiguous_ptm_group_rowdata()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_ambiguous_ptm_group_rowdata.md)
mean the same thing whichever search engine produced the data. Note that
[`parse_PTM_scores_pd()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/parse_PTM_scores_pd.md)
takes its `threshold` as a percentage instead.

The number of PTMs is recovered from the total, which is 100 per
modification. Rows reading `Inconclusive data` or `Too many isoforms`
match nothing and so yield no candidates.

## Usage

``` r
parse_ptm_candidates_pd(obj, prob_col = "ptmRS.Phospho.Site.Probabilities")
```

## Arguments

- obj:

  `SummarizedExperiment`. Proteomics dataset with PD PSM-level rowData

- prob_col:

  `character` Column holding the site probabilities

## Value

`data.frame` with one row per candidate residue and columns `row`,
`pep_pos`, `residue`, `prob` and `n_ptms`
