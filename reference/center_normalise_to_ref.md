# Center-median normalise the expression matrix in an MSnSet using medians from a reference dataset

Center-median normalisation is a simple normalisation method that is
appropriate for relative abundance proteomics such as isobaric tagging.
This can be achieved with `QFeatures::normalize(method='diff.median')`.
However, for some experimental designs, the normalisation should be
against the medians in another dataset. For example, for PTM studies,
one may wish to isobaric tag samples, pool, and then PTM-enriched, with
the enriched sample quantified in a separate run to the non-enriched
(total) sample. In this case, it may make more sense to center-median
normalise the PTM-enriched samples using the median from the total
samaples.

## Usage

``` r
center_normalise_to_ref(
  obj,
  medians,
  center_to_zero = FALSE,
  on_log_scale = FALSE
)
```

## Arguments

- obj:

  `SummarisedExperiment`

- medians:

  `vector, numeric`. Sample medians from reference dataset

- center_to_zero:

  `logical`. Centre the data range on zero. If FALSE, normalisation
  retains original data range.

- on_log_scale:

  `logical`. Input data is log-transformed

## Value

Returns a `SummarisedExperiment` with the assay data column
center-median normalised
