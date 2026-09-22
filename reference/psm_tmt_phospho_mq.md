# PSM-level MaxQuant output for a phospho-enriched TMT18plex fraction

MaxQuant `evidence.txt` PSM-level output for the phospho-enriched
fraction of a real TMT18plex experiment, subsetted to a random selection
of proteins carrying more than one phosphopeptide PSM, plus a sample of
contaminant, decoy and non-unique-master-protein PSMs for use as the
MaxQuant worked example in the `PTM site quantification` vignette.

Includes the `Phospho..STY..Probabilities` and `Phospho..STY.` columns
needed by
[`parse_ptm_candidates_mq`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/parse_ptm_candidates_mq.md),
which is the reason this dataset exists: MaxQuant writes localisation
probabilities inline in the peptide sequence, where Proteome Discoverer
writes them as a separate list of percentages, so
[`psm_tmt_phospho`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/psm_tmt_phospho.md)
cannot demonstrate both.

[`psm_tmt_factorial`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/psm_tmt_factorial.md)
is the total (non-enriched) fraction of the same labelled pool, acquired
from the same plex, so the two share
[`tmt_factorial_design`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/tmt_factorial_design.md)
and sample columns are named to match its rownames. The design is a 3 x
2 of cell line against treatment, three replicates each.

`inst/extdata/tmt_phospho_mq_proteome.fasta.gz` holds the sequences of
the retained master proteins, taken from the search database, so that
[`add_peptide_positions_from_cleavage`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/add_peptide_positions_from_cleavage.md)
can place peptides without querying UniProt.

## Usage

``` r
psm_tmt_phospho_mq
```

## Format

An object of class `data.frame` with 10170 rows and 68 columns.
