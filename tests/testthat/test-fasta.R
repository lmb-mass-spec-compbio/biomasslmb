write_fasta <- function(headers, seqs) {
  path <- tempfile(fileext = ".fasta")
  lines <- character(0)
  for (i in seq_along(headers)) {
    lines <- c(lines, paste0(">", headers[[i]]), seqs[[i]])
  }
  writeLines(lines, path)
  path
}

test_that("get_crap_fasta_accessions extracts the accession between the first pair of pipes", {
  fasta_path <- write_fasta(
    c("sp|P00330|ADH1_YEAST Alcohol dehydrogenase 1", "sp|P62258|1433E_HUMAN 14-3-3 protein epsilon"),
    c("MSIPETQKGVIFYESHGKLEYK", "MDDREDLVYQAKLAEQAERYDE")
  )

  result <- get_crap_fasta_accessions(fasta_path)
  expect_equal(result, c("P00330", "P62258"))
})

test_that("get_maxquant_cont_accessions defaults to the bundled contaminants.fasta.gz", {
  result <- suppressMessages(get_maxquant_cont_accessions())

  expect_true(is.character(result))
  expect_true(length(result) > 0)
  expect_true(all(startsWith(result, "CON__")))
  expect_true("CON__P00761" %in% result)
})

test_that("get_maxquant_cont_accessions parses a custom fasta and prefixes CON__", {
  fasta_path <- write_fasta(
    c("P00330 SWISS-PROT:P00330|ADH1_YEAST Alcohol dehydrogenase 1"),
    c("MSIPETQKGVIFYESHGKLEYK")
  )

  result <- get_maxquant_cont_accessions(fasta_path)
  expect_equal(result, "CON__P00330")
})
