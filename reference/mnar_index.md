# Deprecated: use condition_miss_index()

`mnar_index()` has been renamed to
[`condition_miss_index()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/condition_miss_index.md)
to better reflect that the index summarises missingness predictable from
experimental condition, not MNAR in the classical sense. All arguments
are identical; this wrapper calls the new function and will be removed
in a future release.

## Usage

``` r
mnar_index(...)
```

## Arguments

- ...:

  arguments passed on to
  [`condition_miss_index`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/condition_miss_index.md)

## Details

Note:
[`condition_miss_index()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/condition_miss_index.md)
expects the `$summary` output from
[`condition_miss_score()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/condition_miss_score.md),
which uses updated column names (`condition_miss_score`,
`condition_miss_class`) rather than the old `mnar_score` / `mnar_class`
names.

## See also

[`condition_miss_index`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/condition_miss_index.md)
