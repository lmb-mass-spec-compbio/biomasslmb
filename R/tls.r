#' Total least squares regression (for use with geom_smooth)
#'
#' @param formula Formula of the form y ~ x.
#' @param data A data frame containing columns for x and y.
#' @param ... Ignored.
#' @examples
#' # both axes are noisy estimates, so ordinary least squares would give a
#' # slope that depends on which variable is on the x axis
#' x <- rnorm(100)
#' df <- data.frame(x = x + rnorm(100, sd = 0.3),
#'                  y = x + rnorm(100, sd = 0.3))
#'
#' unlist(tls(y ~ x, df))
#'
#' library(ggplot2)
#' ggplot(df, aes(x, y)) +
#'   geom_point() +
#'   geom_smooth(method = tls, formula = y ~ x, se = FALSE) +
#'   theme_biomasslmb()
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
#' @param object `tls_model` object as returned by `tls_fit`
#' @param newdata `data.frame` with an `x` column, or a `numeric` vector of x values
#' @param se.fit `logical` Not supported for TLS; an error is raised if TRUE
#' @param ... ignored, present for consistency with the `predict` generic
#' @return `numeric` vector of predicted y values.
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
