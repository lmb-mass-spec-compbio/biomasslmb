test_that("get_tm_info parses start/end/length/type from TRANSMEM segments", {
  trans_info_string <- c("1..20; /note=Extracellular", "45..67; /note=Cytoplasmic")

  result <- get_tm_info(trans_info_string)

  expect_equal(result, c("1;45", "20;67", "19;22", "Extracellular;Cytoplasmic",
                          "Cytoplasmic;Extracellular"))
})

test_that("get_tm_info returns NA for empty input", {
  expect_equal(get_tm_info(character(0)), rep(NA, 5))
})

test_that("get_tm_info returns NA for beta-stranded segments", {
  result <- get_tm_info("1..20; /note=Beta stranded")
  expect_equal(result, rep(NA, 5))
})

test_that("add_tm_info adds per-protein transmembrane columns", {
  query_result <- data.frame(
    Transmembrane = c(
      "TRANSMEM 1..20; /note=Extracellular; TRANSMEM 45..67; /note=Cytoplasmic",
      "TRANSMEM 10..30; /note=Extracellular"
    ),
    stringsAsFactors = FALSE
  )

  result <- add_tm_info(query_result)

  expect_equal(result$n_tms, c(2, 1))
  expect_equal(result$tm_start, c("1;45", "10"))
  expect_equal(result$tm_end, c("20;67", "30"))
  expect_equal(result$tm_length, c("19;22", "20"))
  expect_equal(result$tm_types, c("Extracellular;Cytoplasmic", "Extracellular"))
  expect_equal(result$tm_types_set, c("Cytoplasmic;Extracellular", "Extracellular"))
})

test_that("get_topology parses topology, termini, and ranges", {
  # single-residue segment ("21", no "..") exercises the end = start + 1 branch
  topology_string <- c("1..20; /note=Extracellular", "21; /note=Helical", "22..50; /note=Cytoplasmic")

  result <- get_topology(topology_string)

  expect_equal(result, c(
    "Extracellular;Helical;Cytoplasmic",
    "Extracellular",
    "Cytoplasmic",
    "1;21;22",
    "20;22;50",
    "19;1;28"
  ))
})

test_that("get_topology returns NA for empty input", {
  expect_equal(get_topology(character(0)), rep(NA, 6))
})

test_that("add_topology_info adds per-protein topology columns", {
  result <- data.frame(
    Topological.domain = c(
      "TOPO_DOM 1..20; /note=Extracellular; TOPO_DOM 21..44; /note=Cytoplasmic; TOPO_DOM 45..70; /note=Extracellular",
      "TOPO_DOM 1..9; /note=Cytoplasmic; TOPO_DOM 10..30; /note=Extracellular"
    ),
    stringsAsFactors = FALSE
  )

  out <- add_topology_info(result)

  expect_equal(out$topology, c(
    "Extracellular;Cytoplasmic;Extracellular",
    "Cytoplasmic;Extracellular"
  ))
  expect_equal(out$n_term, c("Extracellular", "Cytoplasmic"))
  expect_equal(out$c_term, c("Extracellular", "Extracellular"))
  expect_equal(out$topology_start, c("1;21;45", "1;10"))
  expect_equal(out$topology_end, c("20;44;70", "9;30"))
  expect_equal(out$topology_length, c("19;23;25", "8;20"))
})
