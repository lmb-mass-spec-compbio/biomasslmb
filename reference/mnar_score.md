# Deprecated: use condition_miss_score()

`mnar_score()` has been renamed to
[`condition_miss_score()`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/condition_miss_score.md)
to better reflect that the score measures missingness predictable from
experimental condition, not MNAR in the classical sense. All arguments
are identical; this wrapper calls the new function and will be removed
in a future release.

## Usage

``` r
mnar_score(...)
```

## See also

[`condition_miss_score`](https://lmb-mass-spec-compbio.github.io/biomasslmb/reference/condition_miss_score.md)
