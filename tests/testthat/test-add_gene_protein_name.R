make_se <- function(descriptions) {
  row_data <- S4Vectors::DataFrame(Master.Protein.Descriptions = descriptions)
  assay_mat <- matrix(1, nrow = length(descriptions), ncol = 1, dimnames = list(NULL, "s1"))
  SummarizedExperiment::SummarizedExperiment(assays = list(counts = assay_mat), rowData = row_data)
}

test_that("add_gene_long_protein_name_pd extracts gene name and long protein name", {
  se <- make_se("Serum albumin OS=Bos taurus GN=ALB PE=1 SV=4")
  out <- add_gene_long_protein_name_pd(se)

  expect_equal(SummarizedExperiment::rowData(out)$Master.Protein.gene.name, "ALB")
  expect_equal(SummarizedExperiment::rowData(out)$Master.Protein.long.name, "Serum albumin")
})

test_that("add_gene_long_protein_name_pd returns an empty gene name when GN= is absent", {
  se <- make_se("Some protein OS=Homo sapiens PE=1 SV=1")
  out <- add_gene_long_protein_name_pd(se)

  expect_equal(SummarizedExperiment::rowData(out)$Master.Protein.gene.name, "")
  expect_equal(SummarizedExperiment::rowData(out)$Master.Protein.long.name, "Some protein")
})

test_that("add_gene_long_protein_name_pd handles multiple ';'-delimited descriptions", {
  se <- make_se("Protein A OS=Homo sapiens GN=A1 PE=1 SV=1;Protein B OS=Homo sapiens PE=1 SV=1")
  out <- add_gene_long_protein_name_pd(se)

  expect_equal(SummarizedExperiment::rowData(out)$Master.Protein.gene.name, "A1")
  expect_equal(SummarizedExperiment::rowData(out)$Master.Protein.long.name, "Protein A;Protein B")
})
