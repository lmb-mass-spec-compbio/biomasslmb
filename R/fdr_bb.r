bh_with_alpha <- function(p, alpha) {
  k <- length(p)
  o <- order(p)
  p_sorted <- p[o]

  # Regular BH adjusted p-values
  adj_sorted <- p_sorted * k / seq_len(k)
  adj_sorted <- pmin(adj_sorted, 1)

  adj <- numeric(k)
  adj[o] <- adj_sorted

  list(adj = adj, sig = adj <= alpha)
}


#' Benjamini & Bogomolov (2014) selective-inference adjustment
#' using external screening p-values (e.g., F-test p-values)
#'
#' @param mtx_p     Matrix of raw p-values (families × hypotheses)
#' @param p_screen  Vector of screening p-values (one per family)
#' @param alpha     Target FDR level across selected families
#'
#' @return A list with:
#'   - adj_p: matrix of adjusted p-values (NA for unselected families)
#'   - selected: integer indices of selected families
#'   - alpha_star: the B&B reduced alpha
#' @export
bb2014 <- function(mtx_p, p_screen, alpha = 0.05) {

  # ---- Basic validation ----
  if (!is.matrix(mtx_p))
    stop("mtx_p must be a matrix of hypothesis-level p-values.")

  m <- nrow(mtx_p)

  if (length(p_screen) != m)
    stop("p_screen must have the same length as the number of families (nrow(mtx_p)).")

  if (!is.numeric(p_screen) || any(p_screen < 0 | p_screen > 1))
    stop("p_screen must be valid p-values.")

  # ---- Name requirement ----
  if (is.null(rownames(mtx_p)) || is.null(names(p_screen))) {
    stop("Both mtx_p rownames and names(p_screen) must be provided.")
  }

  # Names must match exactly and in order
  if (!identical(rownames(mtx_p), names(p_screen))) {
    stop("Row names of mtx_p and names(p_screen) must be identical and in the same order.")
  }

  # ---- Selection step ----
  sel_fam <- which(p.adjust(p_screen, method = "BH") <= alpha)
  R <- length(sel_fam)

  # B&B scaling
  alpha_star <- alpha * R / m

  # ---- Output matrix ----
  adj_mat <- matrix(
    NA,
    nrow = m,
    ncol = ncol(mtx_p),
    dimnames = dimnames(mtx_p)
  )

  # ---- Within-family BH with reduced alpha ----
  if (R > 0) {
    for (i in sel_fam) {
      p_i <- mtx_p[i, ]
      adj_i <- bh_with_alpha(p_i, alpha_star)$adj
      adj_mat[i, ] <- adj_i
    }
  }

  list(
    adj_p = adj_mat,
    selected = sel_fam,
    alpha_star = alpha_star
  )
}

