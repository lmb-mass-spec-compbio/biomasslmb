test_that("get_cat_palette returns n hex codes from the 7-colour palette", {
  result <- get_cat_palette(5)
  expect_length(result, 5)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", result)))
})

test_that("get_cat_palette switches to the 12-colour palette above 7", {
  result <- get_cat_palette(10)
  expect_length(result, 10)
  expect_false(identical(result[1:7], get_cat_palette(7)))
})

test_that("get_cat_palette errors outside 1-12", {
  expect_error(get_cat_palette(0), "n must be")
  expect_error(get_cat_palette(13), "n must be")
})
