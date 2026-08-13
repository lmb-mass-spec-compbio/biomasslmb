test_that("bh_with_alpha matches p.adjust(method='BH') and flags significance at alpha", {
  p <- c(0.001, 0.01, 0.02, 0.03, 0.5, 0.8)
  res <- biomasslmb:::bh_with_alpha(p, alpha = 0.05)

  expect_equal(res$adj, p.adjust(p, method = "BH"))
  expect_equal(res$sig, res$adj <= 0.05)
})

test_that("bb2014 validates inputs", {
  mtx_p <- matrix(c(0.01, 0.02, 0.5, 0.6), nrow = 2, dimnames = list(c("f1", "f2"), NULL))

  expect_error(bb2014(as.data.frame(mtx_p), p_screen = c(f1 = 0.01, f2 = 0.5)),
               "must be a matrix")
  expect_error(bb2014(mtx_p, p_screen = c(0.01, 0.5, 0.2)),
               "same length")
  expect_error(bb2014(mtx_p, p_screen = c(f1 = 0.01, f2 = 1.5)),
               "valid p-values")
  expect_error(bb2014(mtx_p, p_screen = c(0.01, 0.5)),
               "must be provided")
  expect_error(bb2014(mtx_p, p_screen = c(f2 = 0.01, f1 = 0.5)),
               "identical and in the same order")
})

test_that("bb2014 only adjusts p-values within selected families", {
  # f1 has a strong screening p-value (selected), f2 does not (not selected)
  mtx_p <- matrix(
    c(0.001, 0.002, 0.003, 0.6, 0.7, 0.8),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("f1", "f2"), c("h1", "h2", "h3"))
  )
  p_screen <- c(f1 = 0.001, f2 = 0.9)

  res <- bb2014(mtx_p, p_screen, alpha = 0.05)

  expect_equal(res$selected, c(f1 = 1))
  expect_true(all(is.na(res$adj_p["f2", ])))
  expect_false(any(is.na(res$adj_p["f1", ])))
  expect_equal(res$alpha_star, 0.05 * 1 / 2)
})

test_that("bb2014 selects no families and returns an all-NA matrix when nothing passes screening", {
  mtx_p <- matrix(
    c(0.01, 0.02, 0.03, 0.04),
    nrow = 2,
    dimnames = list(c("f1", "f2"), c("h1", "h2"))
  )
  p_screen <- c(f1 = 0.9, f2 = 0.95)

  res <- bb2014(mtx_p, p_screen, alpha = 0.05)

  expect_equal(unname(res$selected), integer(0))
  expect_true(all(is.na(res$adj_p)))
  expect_equal(res$alpha_star, 0)
})
