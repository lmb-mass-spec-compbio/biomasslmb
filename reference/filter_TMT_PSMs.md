# Filter a PSM-level summarizedExperiment to remove low quality PSMs

Filter PSMs from TMT quantification to remove the following:

1.  Missing values (`NA`) for all tags

2.  Interference/co-isolation above a set value (default=100, e.g no
    filtering)

3.  Signal:noise ratio below a set value (default=0, e.g no filtering)

4.  SPS-MM below a set value (default=0, e.g no filtering)

5.  Quan.Info is not empty (”)

6.  PSM.Ambiguity is not 'Selected' or 'Unambiguous'

## Usage

``` r
filter_TMT_PSMs(
  obj,
  inter_thresh = 100,
  sn_thresh = 0,
  spsmm_thresh = 0,
  master_protein_col = "Master.Protein.Accessions",
  inter_col = "Isolation.Interference.in.Percent",
  sn_col = "Average.Reporter.SN",
  spsmm_col = "SPS.Mass.Matches.in.Percent",
  from_PD = TRUE,
  verbose = TRUE
)
```

## Arguments

- obj:

  `summarizedExperiment`. Should contain PSMs-level TMT quantification

- inter_thresh:

  `numeric`. Maximum allowed interference/co-isolation

- sn_thresh:

  `numeric`. Minimum allowed average signal:

- spsmm_thresh:

  `numeric`. Minimum allowed value for Synchronous precursor selection
  mass matches (%) - SPS-MM, e.g the percentage of MS3 fragments that
  can be explicitly traced back to the precursor peptides

- master_protein_col:

  `string`. Name of column containing master proteins.

- inter_col:

  `string`. Name of column containing the interference value.

- sn_col:

  `string`. Name of column containing the signal:noise value.

- spsmm_col:

  `string`. Name of column containing the SPS-MM value.

- from_PD:

  `logical`. Input is from ProteomeDiscover. If set, will filter using
  Quan.Info and PSM.Ambiguity columns. Default is TRUE

- verbose:

  `logical`. Default is TRUE, use verbose output messages.

## Value

Returns an `summarizedExperiment` with the filtered PSMs.

## Examples

``` r
tmt_qf <- QFeatures::readQFeatures(assayData = psm_tmt_clock,
  colData = tmt_clock_design,
  quantCols = rownames(tmt_clock_design),
  name = "psms_raw")
#> Checking arguments.
#> Loading data as a 'SummarizedExperiment' object.
#> Formatting sample annotations (colData).
#> Formatting data as a 'QFeatures' object.

# a more accurate average signal:noise than PD reports
tmt_qf[["psms_raw"]] <- update_average_sn(tmt_qf[["psms_raw"]])

# drop PSMs with low signal:noise or high co-isolation
psms_filtered <- filter_TMT_PSMs(tmt_qf[["psms_raw"]],
                                 inter_thresh = 50, sn_thresh = 10)
#> Filtering PSMs...
#> 11281 features found from 790 master proteins => Initial PSMs
#> 10554 features found from 768 master proteins => PSMs with Quan.Info removed
#> 10554 features found from 768 master proteins => PSMs which are not selected or unambiguous removed
#> 6780 features found from 691 master proteins => Removing PSMs without quantification values
#> 6675 features found from 686 master proteins => Removing PSMs with high Co-isolation/interference
#> 5921 features found from 631 master proteins => Removing PSMs with low average S:N ratio
#> Not performing filtering by SPS-MM (`spsmm_thresh`=0)

nrow(psms_filtered)
#> [1] 5921
```
