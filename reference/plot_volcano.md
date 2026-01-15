# Volcano plot for differential abundance testing

Statistical testing results for 'omics' differential abundance are
frequently displayed in a 'volcano' plot, where x = log fold-change and
y = -log10(p-value). This function generates a volcano plot from a
`data.frame` containing the statistical testing results

## Usage

``` r
plot_volcano(
  obj,
  lf_col = "logFC",
  p_col = "adj.P.Val",
  sig_col = NULL,
  title = NULL,
  xlab = "log2 Fold Change"
)
```

## Arguments

- obj:

  `data.frame`. Statistical testing results

- lf_col:

  `string`. Column for log fold-change

- p_col:

  `string`. Column for p-value

- sig_col:

  `string`. Column for significance annotation (optional)

- title:

  `string`. Plot title

- xlab:

  `string`. Label for x axis

## Value

Returns a *ggplot* object.
