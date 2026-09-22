# Group ambiguously localised PTM sites by candidate residue overlap

Filtering on localisation probability discards a peptide unless every
PTM on it has a residue reaching `min_prob`, so a peptide with two
candidate serines tied at 0.5/0.5 is thrown away even though it
establishes that the region is modified and by how much. This function
keeps those peptides by pooling them into a quantifiable group.

Peptides whose sites all reach `min_prob` keep their own site as their
group, exactly as the localisation filter would report it, and take no
part in the graph. Peptides that do not resolve instead contribute every
candidate residue above `min_candidate_prob`. Nodes are candidate
residues in protein coordinates and edges join residues that are
candidates on the same peptide, so connected components chain a "S32 or
S33" peptide with a "S33 or S37" peptide into a single group over S32,
S33 and S37. The graph is built per protein and per PTM count, which
keeps a singly modified peptide from merging with a doubly modified one
over the same residues.

The localisation probabilities on a peptide sum to the number of PTMs on
it, since each modification has to sit somewhere. A candidate's
probability is therefore a share of one modification, and the summed
probability of a group is the chance that the group contains the true
site. That makes `min_candidate_prob` a coverage guarantee rather than
an arbitrary denoising threshold: the probability mass it leaves behind
is the chance the group excludes the real residue.
[`summarise_ptm_groups()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/summarise_ptm_groups.md)
reports that miss rate across a range of thresholds.

A peptide that occurs at more than one position in its protein has no
single set of protein coordinates, and is left unassigned, as is a
peptide with no candidates at all.

## Usage

``` r
add_ambiguous_ptm_group_rowdata(
  obj,
  candidates,
  master_protein_col = "Leading.razor.protein",
  start_col = "start",
  min_prob = 0.501,
  min_candidate_prob = 0.02,
  max_group_span = 50,
  verbose = TRUE
)
```

## Arguments

- obj:

  `SummarizedExperiment`. Proteomics dataset with peptide start
  positions from
  [`add_peptide_positions_from_cleavage()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_peptide_positions_from_cleavage.md)

- candidates:

  `data.frame` from
  [`parse_ptm_candidates_mq()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/parse_ptm_candidates_mq.md)
  or
  [`parse_ptm_candidates_pd()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/parse_ptm_candidates_pd.md)

- master_protein_col:

  `character` Column identifying the parent protein

- start_col:

  `character` Column holding the peptide start position

- min_prob:

  `numeric` Probability at or above which a residue counts as localised.
  On a 0 to 1 scale

- min_candidate_prob:

  `numeric` Probability below which a residue on an unresolved peptide
  is treated as ruled out rather than as a candidate. On a 0 to 1 scale

- max_group_span:

  `numeric` Widest span in residues a group may cover

- verbose:

  `logical` Report how many rows were assigned to a group

## Value

`SummarizedExperiment` with rowData columns `ptm_group_id`,
`ptm_group_members` (the protein positions, semi-colon separated),
`ptm_group_n_candidates` and `ptm_group_resolved`. All are `NA` for a
row that could not be assigned. `ptm_group_n_candidates` is 1 for a
peptide whose sites are all localised, since there is one possible
assignment, and otherwise the number of candidate residues pooled into
the group

## Limitations

A pooled group has no single residue, so it has no motif and cannot be
passed to
[`add_site_sequence()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_site_sequence.md)
or used for kinase and motif enrichment. Downstream labelling has to
carry `ptm_group_n_candidates` through to the site labels, so that a
pooled group is never rendered as if it were one residue. Pooling
intensity across residues can average away opposing changes at
neighbouring sites. And a pooled group's candidates may overlap a
residue that is also quantified as its own resolved site, so the same
modification event can contribute to two features; the consequences of
that for FDR and for enrichment analyses are not quantified.
