test_that("append_fasta adds sequences from file1 onto the end of file2", {
  file1 <- tempfile(fileext = ".fasta")
  file2 <- tempfile(fileext = ".fasta")

  Biostrings::writeXStringSet(
    Biostrings::AAStringSet(c(seq1 = "MSIPETQKGV")), filepath = file1
  )
  Biostrings::writeXStringSet(
    Biostrings::AAStringSet(c(seq2 = "MDDREDLVYQ")), filepath = file2
  )

  append_fasta(file1, file2)

  result <- Biostrings::readAAStringSet(file2)
  expect_equal(names(result), c("seq2", "seq1"))
  expect_equal(as.character(result[["seq1"]]), "MSIPETQKGV")
})

test_that("append_fasta adds cRAP numbers to file1 headers when is_crap = TRUE", {
  file1 <- tempfile(fileext = ".fasta")
  file2 <- tempfile(fileext = ".fasta")

  Biostrings::writeXStringSet(
    Biostrings::AAStringSet(c("|proteaseA|" = "MSIPETQKGV", "|proteaseB|" = "MDDREDLVYQ")),
    filepath = file1
  )
  Biostrings::writeXStringSet(
    Biostrings::AAStringSet(c(seq2 = "MDDREDLVYQ")), filepath = file2
  )

  append_fasta(file1, file2, is_crap = TRUE, crap_start = 5)

  result <- Biostrings::readAAStringSet(file2)
  expect_equal(names(result), c("seq2", "|cRAP005|proteaseA|", "|cRAP006|proteaseB|"))
})
