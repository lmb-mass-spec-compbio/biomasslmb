#' Total least squares regression (for use with geom_smooth)
#'
#' @param formula Formula of the form y ~ x.
#' @param data A data frame containing columns for x and y.
#' @param ... Ignored.
#' @export
tls <- function(formula, data, ...) {
  M <- model.frame(formula, data)
  x <- M[[2]]
  y <- M[[1]]

  if (length(unique(x)) < 2 || length(unique(y)) < 2) {
    slope <- NA_real_
    intercept <- NA_real_
  } else {
    pca <- prcomp(cbind(x, y))
    slope <- pca$rotation[2, 1] / pca$rotation[1, 1]
    intercept <- mean(y) - slope * mean(x)
  }

  # return a model-like object
  structure(
    list(intercept = intercept, slope = slope),
    class = "tls_model"
  )
}

#' Predict method for tls_model
#'
#' @export
predict.tls_model <- function(object, newdata, se.fit = FALSE, ...) {
  if (is.data.frame(newdata) && "x" %in% names(newdata)) {
    x <- newdata$x
  } else {
    x <- newdata
  }

  y <- object$intercept + object$slope * x
  if (se.fit) stop("Standard errors not supported for TLS.")
  y
}
