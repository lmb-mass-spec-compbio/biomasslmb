test_that("make_unique_all numbers all members of a duplicated group from 1", {
  expect_equal(make_unique_all(c("A", "A", "B"), always = FALSE), c("A.1", "A.2", "B"))
})

test_that("make_unique_all handles multiple duplicate groups", {
  expect_equal(make_unique_all(c("A", "A", "B", "B", "C"), always = FALSE),
               c("A.1", "A.2", "B.1", "B.2", "C"))
})

test_that("make_unique_all suffixes singletons when always = TRUE", {
  expect_equal(make_unique_all(c("A", "A", "B"), always = TRUE), c("A.1", "A.2", "B.1"))
})

test_that("make_unique_all respects a custom separator", {
  expect_equal(make_unique_all(c("X", "X", "Y"), sep = "_", always = FALSE), c("X_1", "X_2", "Y"))
})
