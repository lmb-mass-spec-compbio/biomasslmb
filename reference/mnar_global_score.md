# Deprecated: use global_condition_miss_score()

`mnar_global_score()` has been renamed to
[`global_condition_miss_score()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/global_condition_miss_score.md)
to better reflect that the score measures missingness predictable from
experimental condition, not MNAR in the classical sense. All arguments
are identical; this wrapper calls the new function and will be removed
in a future release.

## Usage

``` r
mnar_global_score(...)
```

## Arguments

- ...:

  arguments passed on to
  [`global_condition_miss_score`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/global_condition_miss_score.md)

## See also

[`global_condition_miss_score`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/global_condition_miss_score.md)
