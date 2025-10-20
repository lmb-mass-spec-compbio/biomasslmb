#' Make Elements of a Character Vector Unique (Numbering From 1)
#'
#' A variant of [base::make.unique()] that appends numeric suffixes
#' to elements of a character vector. Unlike `make.unique()`, which
#' leaves the first occurrence unchanged, this function can number
#' all members of a duplicated group starting from 1. Optionally,
#' it can number singletons as well.
#'
#' @param x A character vector.
#' @param sep A character string to separate the original element
#'   from the suffix. Defaults to `"."`.
#' @param always Logical. If `TRUE`, all elements (including singletons)
#'   get a suffix (e.g. `"B"` becomes `"B.1"`). If `FALSE` (default),
#'   singletons remain unchanged.
#'
#' @return A character vector of the same length as `x`,
#'   with duplicates disambiguated by appending a numeric suffix.
#'
#' @examples
#' # Mimics make.unique() but numbers first duplicates as well
#' make_unique_all(c("A", "A", "B"))
#' # [1] "A.1" "A.2" "B"
#'
#' # Multiple duplicate groups
#' make_unique_all(c("A", "A", "B", "B", "C"))
#' # [1] "A.1" "A.2" "B.1" "B.2" "C"
#'
#' # Force suffixes on all elements
#' make_unique_all(c("A", "A", "B"), always = TRUE)
#' # [1] "A.1" "A.2" "B.1"
#'
#' # Custom separator
#' make_unique_all(c("X", "X", "Y"), sep = "_")
#' # [1] "X_1" "X_2" "Y"
#'
#' @seealso [base::make.unique()]
#'
#' @export
make_unique_all <- function(x, sep = ".", always = TRUE) {
  # Number each occurrence within its group
  counts <- ave(seq_along(x), x, FUN = seq_along)

  # Get group sizes
  dup_counts <- ave(seq_along(x), x, FUN = length)

  # Decide whether to suffix
  if (always) {
    paste(x, counts, sep = sep)
  } else {
    ifelse(dup_counts > 1, paste(x, counts, sep = sep), x)
  }
}
