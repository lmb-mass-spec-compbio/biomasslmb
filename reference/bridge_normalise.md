# Put TMT plexes onto a common scale using bridge (pooled reference) channels

An experiment too large for a single TMT plex is split over several
plexes, each of which is a separate MS run. A feature is therefore
quantified from a different set of PSMs in each plex, so its abundances
are not directly comparable between plexes and merging them confounds
plex with biology.

A bridge channel is a pooled sample, identical in every plex, which
provides a per-feature reference point. For each feature, this function
takes the difference between its abundance in a plex's bridge channel(s)
and its abundance across the plexes' bridge channels, and removes that
difference from every sample in the plex.

## Usage

``` r
bridge_normalise(
  obj,
  plex_col,
  bridge_cols,
  fun = mean,
  reference = NULL,
  on_missing = c("na", "drop", "ignore"),
  min_bridge = 1,
  on_log_scale = TRUE,
  verbose = TRUE
)
```

## Arguments

- obj:

  `SummarizedExperiment` containing quantification for all plexes, e.g.
  the assay produced by
  [`QFeatures::joinAssays`](https://rformassspectrometry.github.io/QFeatures/reference/joinAssays.html).

- plex_col:

  `string`. Name of the `colData` column identifying the plex.

- bridge_cols:

  `logical` vector of length `ncol(obj)`, or `character` vector of
  column names, identifying the bridge channels.

- fun:

  `function`. Used to collapse several bridge channels within a plex to
  a single reference. Applied to the non-missing values only. Default is
  `mean`.

- reference:

  `string`. Plex to align the others to. Default (`NULL`) is to align
  every plex to the mean of the per-plex references, which preserves the
  overall abundance scale.

- on_missing:

  `string`. What to do with a feature that has no usable bridge value in
  a plex: `'na'` (default) sets that plex's values for the feature to
  `NA`, `'drop'` removes the feature, and `'ignore'` leaves that plex
  uncorrected.

- min_bridge:

  `numeric`. Minimum number of non-missing bridge values required in a
  plex for the reference to be used. Default is 1.

- on_log_scale:

  `logical`. Input data is log-transformed, so the correction is
  subtractive. If `FALSE`, the correction is a ratio.

- verbose:

  `logical`. Report the size of the correction and the number of
  features affected by missing bridge values. Default is TRUE.

## Value

Returns a `SummarizedExperiment` with the plexes on a common scale.

## Details

Unlike
[`center_normalise_to_ref()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/center_normalise_to_ref.md),
which removes a single offset per sample, the correction here is per
feature *and* per plex, because that is what a plex effect is: different
plexes quantify a feature from different PSMs. The function is agnostic
to the feature level, so it can be applied to a joined protein-level
assay or to a joined peptide- or site-level assay.

**Several bridge channels per plex.** These are collapsed to one
reference per plex with `fun`, ignoring missing values, before the
plexes are compared. Collapsing within a plex first matters when the
plexes carry unequal numbers of bridge channels, since pooling all
bridge channels together would give a plex with more of them more weight
in defining the common scale.

**Features without a usable bridge value in a plex.** These cannot be
placed on the common scale, and `on_missing` decides what happens to
them. Note that `on_missing='ignore'` does not assume that the plex has
no offset for the feature; because the reference is centred across
plexes, it assumes the plex's offset is the average one.

After correction the bridge channels agree by construction, so their
agreement is not evidence that the correction worked. The evidence is
what happened to the other samples, e.g. with
[`plot_pca()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/plot_pca.md)
coloured by plex.
