# Grouping ambiguously localised PTM sites.
#
# Everything search-engine specific lives in a parse_ptm_candidates_*()
# function, each returning the same long table of candidate residues:
#
#   row      index of the row in the SummarizedExperiment
#   pep_pos  position of the residue within the peptide
#   residue  the amino acid
#   prob     localisation probability, rescaled to 0 to 1
#   n_ptms   number of PTMs on the peptide
#
# add_ambiguous_ptm_group_rowdata() reads only that table and the peptide start
# positions, so supporting another search output means writing another parser.

#' Parse MaxQuant PTM localisation probabilities into candidate residues
#'
#' @description Extracts every candidate residue for a modification from the
#' MaxQuant probability column, rather than only the residues that are
#' confidently localised. MaxQuant writes the probabilities inline in the
#' peptide sequence, so a residue's position is the length of the text before
#' its value.
#'
#' Probabilities are on a 0 to 1 scale, as MaxQuant reports them, so `min_prob`
#' and `min_candidate_prob` in `add_ambiguous_ptm_group_rowdata()` mean the same
#' thing whichever search engine produced the data. Note that
#' `parse_PTM_scores_pd()` takes its `threshold` as a percentage instead.
#'
#' The output is the input to `add_ambiguous_ptm_group_rowdata()`. Use
#' `parse_ptm_candidates_pd()` for Proteome Discoverer data.
#'
#' @param obj `SummarizedExperiment`. Proteomics dataset with MaxQuant
#'   evidence.txt rowData
#' @param prob_col `character` Column holding the probability string
#' @param n_ptm_col `character` Column holding the number of PTMs on the peptide
#' @param sequence_col `character` Column holding the peptide sequence
#' @return `data.frame` with one row per candidate residue and columns `row`,
#'   `pep_pos`, `residue`, `prob` and `n_ptms`
#' @export
parse_ptm_candidates_mq <- function(obj,
                                    prob_col = 'Phospho..STY..Probabilities',
                                    n_ptm_col = 'Phospho..STY.',
                                    sequence_col = 'Sequence') {

  check_se(obj)

  prob_str <- rowData(obj)[[prob_col]]
  matches <- gregexpr('\\((0\\.\\d+|1|0)\\)', prob_str)
  values <- regmatches(prob_str, matches)
  between <- regmatches(prob_str, matches, invert = TRUE)

  pep_pos <- unlist(lapply(between, function(x) head(cumsum(nchar(x)), -1)))

  data.frame(row = rep(seq_along(values), lengths(values)),
             pep_pos = pep_pos,
             residue = substr(rep(rowData(obj)[[sequence_col]], lengths(values)),
                              pep_pos, pep_pos),
             prob = as.numeric(gsub('[()]', '', unlist(values))),
             n_ptms = rep(as.integer(rowData(obj)[[n_ptm_col]]), lengths(values)))
}


#' Parse Proteome Discoverer ptmRS site probabilities into candidate residues
#'
#' @description Extracts every candidate residue for a modification from the
#' ptmRS output, rather than only the residues that are confidently localised.
#'
#' This reads `ptmRS.Phospho.Site.Probabilities`, **not** the
#' `ptmRS.Best.Site.Probabilities` column that `parse_PTM_scores_pd()` uses. The
#' Best Site column reports only the winning isoform's sites, so the
#' alternatives a group is built from are absent from it. Pointing this function
#' at the Best Site column fails silently, producing groups that contain only
#' the sites ptmRS already preferred.
#'
#' ```
#' Best Site Probabilities   : S3(Phospho): 39.61; S4(Phospho): 39.61
#' Phospho Site Probabilities: S(1): 1.3; S(3): 39.6; S(4): 39.6; S(6): 6.5
#' ```
#'
#' ptmRS reports percentages; they are rescaled to a 0 to 1 probability so that
#' `min_prob` and `min_candidate_prob` in `add_ambiguous_ptm_group_rowdata()`
#' mean the same thing whichever search engine produced the data. Note that
#' `parse_PTM_scores_pd()` takes its `threshold` as a percentage instead.
#'
#' The number of PTMs is recovered from the total, which is 100 per
#' modification. Rows reading `Inconclusive data` or `Too many isoforms` match
#' nothing and so yield no candidates.
#'
#' @param obj `SummarizedExperiment`. Proteomics dataset with PD PSM-level
#'   rowData
#' @param prob_col `character` Column holding the site probabilities
#' @return `data.frame` with one row per candidate residue and columns `row`,
#'   `pep_pos`, `residue`, `prob` and `n_ptms`
#' @export
parse_ptm_candidates_pd <- function(obj,
                                    prob_col = 'ptmRS.Phospho.Site.Probabilities') {

  check_se(obj)

  prob_str <- rowData(obj)[[prob_col]]
  site <- '([A-Z])\\((\\d+)\\):\\s*([0-9.]+)'
  values <- regmatches(prob_str, gregexpr(site, prob_str))

  if (sum(lengths(values)) == 0) {
    return(data.frame(row = integer(), pep_pos = numeric(), residue = character(),
                      prob = numeric(), n_ptms = integer()))
  }

  fields <- do.call(rbind, regmatches(unlist(values), regexec(site, unlist(values))))

  data.frame(row = rep(seq_along(values), lengths(values)),
             pep_pos = as.numeric(fields[, 3]),
             residue = fields[, 2],
             prob = as.numeric(fields[, 4]) / 100) %>%
    group_by(.data$row) %>%
    mutate(n_ptms = as.integer(round(sum(.data$prob)))) %>%
    ungroup() %>%
    as.data.frame()
}


#' Label the connected components of an undirected graph
#'
#' Union-find over a list of edges. The package needs connected components and
#' nothing else from graph theory, and the components formed here are small, so
#' this avoids a dependency on igraph.
#'
#' @param n_nodes `numeric` number of nodes, indexed 1 to n_nodes.
#' @param edges two column `matrix` of node indices.
#' @return `integer` component of each node, numbered from 1.
#' @noRd
connected_components <- function(n_nodes, edges) {

  parent <- seq_len(n_nodes)

  root <- function(i) {
    while (parent[i] != i) i <- parent[i]
    i
  }

  for (k in seq_len(nrow(edges))) {
    a <- root(edges[k, 1])
    b <- root(edges[k, 2])
    if (a != b) parent[b] <- a
  }

  roots <- vapply(seq_len(n_nodes), root, integer(1))
  match(roots, unique(roots))
}


#' Split components that span too much of the protein
#'
#' A component wider than max_group_span is a chain of overlapping peptides
#' walking along the protein rather than one ambiguous site, so it is cut at its
#' largest internal gap, repeatedly until every piece fits. The span is a
#' distance in residues along the protein, not a count of sites.
#'
#' @param nodes `numeric` protein positions, one per node.
#' @param membership `integer` component of each node.
#' @param max_group_span `numeric` widest span in residues a component may cover.
#' @return `integer` updated component membership.
#' @noRd
split_wide_components <- function(nodes, membership, max_group_span) {
  repeat {
    spans <- tapply(nodes, membership, function(x) max(x) - min(x))
    wide <- as.integer(names(spans)[spans > max_group_span])
    if (length(wide) == 0) break
    for (component in wide) {
      positions <- sort(nodes[membership == component])
      cut_at <- positions[which.max(diff(positions))]
      membership[membership == component & nodes > cut_at] <- max(membership) + 1L
    }
  }
  membership
}


#' Group ambiguously localised PTM sites by candidate residue overlap
#'
#' @description Filtering on localisation probability discards a peptide unless
#' every PTM on it has a residue reaching `min_prob`, so a peptide with two
#' candidate serines tied at 0.5/0.5 is thrown away even though it establishes
#' that the region is modified and by how much. This function keeps those
#' peptides by pooling them into a quantifiable group.
#'
#' Peptides whose sites all reach `min_prob` keep their own site as their group,
#' exactly as the localisation filter would report it, and take no part in the
#' graph. Peptides that do not resolve instead contribute every candidate
#' residue above `min_candidate_prob`. Nodes are candidate residues in protein
#' coordinates and edges join residues that are candidates on the same peptide,
#' so connected components chain a "S32 or S33" peptide with a "S33 or S37"
#' peptide into a single group over S32, S33 and S37. The graph is built per
#' protein and per PTM count, which keeps a singly modified peptide from merging
#' with a doubly modified one over the same residues.
#'
#' The localisation probabilities on a peptide sum to the number of PTMs on it,
#' since each modification has to sit somewhere. A candidate's probability is
#' therefore a share of one modification, and the summed probability of a group
#' is the chance that the group contains the true site. That makes
#' `min_candidate_prob` a coverage guarantee rather than an arbitrary denoising
#' threshold: the probability mass it leaves behind is the chance the group
#' excludes the real residue. `summarise_ptm_groups()` reports that miss rate
#' across a range of thresholds.
#'
#' A peptide that occurs at more than one position in its protein has no single
#' set of protein coordinates, and is left unassigned, as is a peptide with no
#' candidates at all.
#'
#' @section Limitations:
#' A pooled group has no single residue, so it has no motif and cannot be passed
#' to `add_site_sequence()` or used for kinase and motif enrichment. Downstream
#' labelling has to carry `ptm_group_n_candidates` through to the site labels,
#' so that a pooled group is never rendered as if it were one residue. Pooling
#' intensity across residues can average away opposing changes at neighbouring
#' sites. And a pooled group's candidates may overlap a residue that is also
#' quantified as its own resolved site, so the same modification event can
#' contribute to two features; the consequences of that for FDR and for
#' enrichment analyses are not quantified.
#'
#' @param obj `SummarizedExperiment`. Proteomics dataset with peptide start
#'   positions from `add_peptide_positions_from_cleavage()`
#' @param candidates `data.frame` from `parse_ptm_candidates_mq()` or
#'   `parse_ptm_candidates_pd()`
#' @param master_protein_col `character` Column identifying the parent protein
#' @param start_col `character` Column holding the peptide start position
#' @param min_prob `numeric` Probability at or above which a residue counts as
#'   localised. On a 0 to 1 scale
#' @param min_candidate_prob `numeric` Probability below which a residue on an
#'   unresolved peptide is treated as ruled out rather than as a candidate. On a
#'   0 to 1 scale
#' @param max_group_span `numeric` Widest span in residues a group may cover
#' @param verbose `logical` Report how many rows were assigned to a group
#' @return `SummarizedExperiment` with rowData columns `ptm_group_id`,
#'   `ptm_group_members` (the protein positions, semi-colon separated),
#'   `ptm_group_n_candidates` and `ptm_group_resolved`. All are `NA` for a row
#'   that could not be assigned. `ptm_group_n_candidates` is 1 for a peptide
#'   whose sites are all localised, since there is one possible assignment, and
#'   otherwise the number of candidate residues pooled into the group
#' @export
add_ambiguous_ptm_group_rowdata <- function(obj,
                                            candidates,
                                            master_protein_col = 'Leading.razor.protein',
                                            start_col = 'start',
                                            min_prob = 0.501,
                                            min_candidate_prob = 0.02,
                                            max_group_span = 50,
                                            verbose = TRUE) {

  check_se(obj)

  rd <- data.frame(rowData(obj))

  if (!start_col %in% colnames(rd)) {
    stop(sprintf("column `%s` is not in the rowData. Run add_peptide_positions_from_cleavage() first", start_col))
  }

  # a peptide matching its protein in more than one place gets a semicolon
  # separated start, and so has no single set of protein coordinates
  start <- suppressWarnings(as.numeric(rd[[start_col]]))

  candidates <- candidates %>%
    mutate(pos = .data$pep_pos + start[.data$row] - 1,
           protein = rd[[master_protein_col]][.data$row],
           stratum = paste(.data$protein, .data$n_ptms, sep = '_n')) %>%
    filter(!is.na(.data$pos)) %>%
    group_by(.data$row) %>%
    mutate(resolved = sum(.data$prob >= min_prob) == .data$n_ptms[1]) %>%
    ungroup()

  resolved_groups <- candidates %>%
    filter(.data$resolved, .data$prob >= min_prob) %>%
    group_by(.data$row, .data$stratum) %>%
    summarise(group_members = paste(sort(unique(.data$pos)), collapse = ';'),
              group_n = 1L, .groups = 'drop')

  ambiguous <- candidates %>%
    filter(!.data$resolved, .data$prob >= min_candidate_prob)

  ambiguous_groups <- data.frame(row = numeric(), stratum = character(),
                                 group_members = character(), group_n = integer())

  if (nrow(ambiguous) > 0) {

    ambiguous_components <- lapply(split(ambiguous, ambiguous$stratum), function(stratum) {

      nodes <- sort(unique(stratum$pos))

      # every candidate on a peptide is joined to the peptide's first candidate,
      # which connects them all without enumerating every pair
      edges <- lapply(split(match(stratum$pos, nodes), stratum$row), function(node) {
        node <- unique(node)
        if (length(node) > 1) cbind(node[1], node[-1])
      })
      edges <- do.call(rbind, edges)
      if (is.null(edges)) edges <- matrix(integer(), ncol = 2)

      data.frame(pos = nodes,
                 component = split_wide_components(
                   nodes, connected_components(length(nodes), edges), max_group_span))
    }) %>%
      bind_rows(.id = 'stratum')

    ambiguous <- merge(ambiguous, ambiguous_components, by = c('stratum', 'pos'))

    component_members <- ambiguous %>%
      group_by(.data$stratum, .data$component) %>%
      summarise(group_members = paste(sort(unique(.data$pos)), collapse = ';'),
                group_n = n_distinct(.data$pos), .groups = 'drop')

    ambiguous_groups <- ambiguous %>%
      distinct(.data$row, .data$stratum, .data$component) %>%
      # a peptide's candidates are all in one component unless max_group_span
      # has split them, in which case it has no single group and is dropped
      group_by(.data$row) %>%
      filter(n() == 1) %>%
      ungroup() %>%
      merge(component_members, by = c('stratum', 'component')) %>%
      select('row', 'stratum', 'group_members', 'group_n')
  }

  groups <- bind_rows(resolved_groups, ambiguous_groups) %>%
    transmute(row = .data$row,
              ptm_group_id = sprintf('%s_%s', .data$stratum, .data$group_members),
              ptm_group_members = .data$group_members,
              ptm_group_n_candidates = .data$group_n)

  groups <- data.frame(row = seq_len(nrow(rd))) %>%
    left_join(groups, by = 'row') %>%
    mutate(ptm_group_resolved = .data$ptm_group_n_candidates == 1)

  if (verbose) {
    message(sprintf('%s of %s rows assigned to a group; %s single residue, %s pooled',
                    sum(!is.na(groups$ptm_group_id)), nrow(groups),
                    sum(groups$ptm_group_resolved, na.rm = TRUE),
                    sum(!groups$ptm_group_resolved, na.rm = TRUE)))
  }

  rowData(obj) <- cbind(rowData(obj),
                        groups[, c('ptm_group_id', 'ptm_group_members',
                                   'ptm_group_n_candidates', 'ptm_group_resolved')])
  obj
}


#' Summarise the PTM groups formed at a given min_candidate_prob
#'
#' @description Runs `add_ambiguous_ptm_group_rowdata()` at one
#' `min_candidate_prob` and reports the size and span of the pooled groups it
#' forms, alongside how often those groups exclude the true site. Sweeping it
#' over a range of thresholds is how `min_candidate_prob` is chosen: raising it
#' buys smaller, more interpretable groups and pays in coverage.
#'
#' The miss rate columns read directly off the probabilities. Because they sum
#' to one per modification, the mass a threshold leaves behind is the chance
#' that the group it forms does not contain the real residue, so
#' `pc_miss_over_5` is the percentage of pooled peptides whose group has more
#' than a 5\% chance of excluding it.
#'
#' @param obj `SummarizedExperiment`. As passed to
#'   `add_ambiguous_ptm_group_rowdata()`
#' @param candidates `data.frame` of candidate residues from one of the
#'   `parse_ptm_candidates_*()` functions
#' @param min_candidate_prob `numeric` Threshold to summarise
#' @param min_prob `numeric` Probability at or above which a residue counts as
#'   localised
#' @param max_group_span `numeric` Widest span in residues a group may cover
#' @param ... Further arguments for `add_ambiguous_ptm_group_rowdata()`
#' @return One row `data.frame` of group size, span and coverage statistics
#' @export
summarise_ptm_groups <- function(obj, candidates, min_candidate_prob,
                                 min_prob = 0.501, max_group_span = 50, ...) {

  obj <- suppressMessages(add_ambiguous_ptm_group_rowdata(
    obj, candidates, min_prob = min_prob,
    min_candidate_prob = min_candidate_prob,
    max_group_span = max_group_span, verbose = FALSE, ...))

  rd <- data.frame(rowData(obj))
  pooled_rows <- which(!rd$ptm_group_resolved %in% TRUE & !is.na(rd$ptm_group_id))

  pooled <- rd[pooled_rows, ] %>%
    distinct(.data$ptm_group_id, .data$ptm_group_members, .data$ptm_group_n_candidates)
  span <- sapply(strsplit(pooled$ptm_group_members, ';'),
                 function(pos) diff(range(as.numeric(pos))))

  retained <- candidates %>%
    filter(.data$row %in% pooled_rows) %>%
    group_by(.data$row, .data$n_ptms) %>%
    summarise(retained = sum(.data$prob[.data$prob >= min_candidate_prob]),
              .groups = 'drop') %>%
    mutate(fraction = .data$retained / .data$n_ptms)

  data.frame(min_candidate_prob = min_candidate_prob,
             pooled_groups = nrow(pooled),
             median_size = stats::median(pooled$ptm_group_n_candidates),
             p90_size = stats::quantile(pooled$ptm_group_n_candidates, 0.9, names = FALSE),
             max_size = max(pooled$ptm_group_n_candidates),
             pc_size_5_plus = round(100 * mean(pooled$ptm_group_n_candidates >= 5), 1),
             median_span = stats::median(span),
             p95_span = stats::quantile(span, 0.95, names = FALSE),
             pc_miss_over_1 = round(100 * mean(retained$fraction < 0.99), 2),
             pc_miss_over_5 = round(100 * mean(retained$fraction < 0.95), 2),
             pc_miss_over_10 = round(100 * mean(retained$fraction < 0.90), 2))
}
