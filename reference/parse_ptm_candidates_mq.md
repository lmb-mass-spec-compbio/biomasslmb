# Parse MaxQuant PTM localisation probabilities into candidate residues

Extracts every candidate residue for a modification from the MaxQuant
probability column, rather than only the residues that are confidently
localised. MaxQuant writes the probabilities inline in the peptide
sequence, so a residue's position is the length of the text before its
value.

Probabilities are on a 0 to 1 scale, as MaxQuant reports them, so
`min_prob` and `min_candidate_prob` in
[`add_ambiguous_ptm_group_rowdata()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_ambiguous_ptm_group_rowdata.md)
mean the same thing whichever search engine produced the data. Note that
[`parse_PTM_scores_pd()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/parse_PTM_scores_pd.md)
takes its `threshold` as a percentage instead.

The output is the input to
[`add_ambiguous_ptm_group_rowdata()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_ambiguous_ptm_group_rowdata.md).
Use
[`parse_ptm_candidates_pd()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/parse_ptm_candidates_pd.md)
for Proteome Discoverer data.

## Usage

``` r
parse_ptm_candidates_mq(
  obj,
  prob_col = "Phospho..STY..Probabilities",
  n_ptm_col = "Phospho..STY.",
  sequence_col = "Sequence"
)
```

## Arguments

- obj:

  `SummarizedExperiment`. Proteomics dataset with MaxQuant evidence.txt
  rowData

- prob_col:

  `character` Column holding the probability string

- n_ptm_col:

  `character` Column holding the number of PTMs on the peptide

- sequence_col:

  `character` Column holding the peptide sequence

## Value

`data.frame` with one row per candidate residue and columns `row`,
`pep_pos`, `residue`, `prob` and `n_ptms`
