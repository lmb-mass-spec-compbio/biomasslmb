# Predict method for tls_model

Predict method for tls_model

## Usage

``` r
# S3 method for class 'tls_model'
predict(object, newdata, se.fit = FALSE, ...)
```

## Arguments

- object:

  `tls_model` object as returned by `tls_fit`

- newdata:

  `data.frame` with an `x` column, or a `numeric` vector of x values

- se.fit:

  `logical` Not supported for TLS; an error is raised if TRUE

- ...:

  ignored, present for consistency with the `predict` generic

## Value

`numeric` vector of predicted y values.
