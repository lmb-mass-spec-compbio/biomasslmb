test_that("tls recovers a known slope/intercept from noise-free data", {
  df <- data.frame(x = 1:10, y = 3 + 2 * (1:10))
  fit <- tls(y ~ x, df)

  expect_s3_class(fit, "tls_model")
  expect_equal(fit$slope, 2, tolerance = 1e-8)
  expect_equal(fit$intercept, 3, tolerance = 1e-8)
})

test_that("tls returns NA slope/intercept when x or y is constant", {
  df <- data.frame(x = rep(1, 5), y = 1:5)
  fit <- tls(y ~ x, df)

  expect_true(is.na(fit$slope))
  expect_true(is.na(fit$intercept))
})

test_that("predict.tls_model applies the fitted line to new data", {
  df <- data.frame(x = 1:10, y = 3 + 2 * (1:10))
  fit <- tls(y ~ x, df)

  expect_equal(predict(fit, data.frame(x = c(0, 5))), c(3, 13), tolerance = 1e-8)
  expect_equal(predict(fit, c(0, 5)), c(3, 13), tolerance = 1e-8)
})

test_that("predict.tls_model errors if se.fit is requested", {
  df <- data.frame(x = 1:10, y = 3 + 2 * (1:10))
  fit <- tls(y ~ x, df)

  expect_error(predict(fit, data.frame(x = 1), se.fit = TRUE), "not supported")
})
