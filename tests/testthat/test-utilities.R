test_that("remove_dots collapses runs of duplicated full stops to a single one", {
  expect_equal(remove_dots(c("a..b", "c...d", "e.f")), c("a.b", "c.d", "e.f"))
})

test_that("remove_dots leaves single full stops and non-dotted strings unchanged", {
  expect_equal(remove_dots("no.dots.here"), "no.dots.here")
  expect_equal(remove_dots("nodots"), "nodots")
})

test_that("remove_x strips a single leading capital X", {
  expect_equal(remove_x(c("X1", "X2", "Y1")), c("1", "2", "Y1"))
})

test_that("remove_x is case-sensitive and only strips a leading X", {
  expect_equal(remove_x("x1"), "x1")
  expect_equal(remove_x("AX1"), "AX1")
})

test_that("message_parse reports feature and unique master protein counts", {
  x <- data.frame(protein = c("P1", "P1", "P2"))
  expect_message(
    message_parse(x, "protein", "test step"),
    "3 features found from 2 master proteins => test step"
  )
})
