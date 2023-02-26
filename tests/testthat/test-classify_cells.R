context("test-classify_cells.R")

library(org.Hs.eg.db)
library(org.Mm.eg.db)
data(test_classifier)
data(test_cds)

set.seed(10)
new_cds <- garnett:::get_communities(test_cds)

test_that("get_communities works", {
  skip_on_travis()
  expect_identical(counts(new_cds), counts(test_cds))
  expect_equal(ncol(colData(new_cds)) - 1, ncol(colData(test_cds)))
  expect_identical(rowData(new_cds), rowData(test_cds))
  expect_equal(length(unique(colData(new_cds)$louv_cluster)), 5)
  expect_equal(sum(colData(new_cds)$louv_cluster == 2), 172)
  expect_equal(as.character(colData(new_cds)$louv_cluster[5]), "2")
})

# classify cells
set.seed(10)
new_cds <- classify_cells(test_cds, test_classifier,
                          db = org.Hs.eg.db,
                          rank_prob_ratio = 1.5,
                          cluster_extend = TRUE,
                          cds_gene_id_type = "SYMBOL")

test_that("classify_cells works", {
  expect_identical(counts(new_cds), counts(test_cds))
  expect_identical(rowData(new_cds), rowData(test_cds))
  expect_equal(sum(colData(new_cds)$cell_type == "B cells"), 210)
  expect_equal(sum(colData(new_cds)$cell_type == "CD4 T cells"), 76)
  expect_equal(sum(colData(new_cds)$cell_type == "CD8 T cells"), 47)
  expect_equal(sum(colData(new_cds)$cell_type == "T cells"), 103)
  expect_equal(sum(colData(new_cds)$cluster_ext_type == "B cells"), 401)
  expect_equal(sum(colData(new_cds)$cluster_ext_type == "CD4 T cells"), 200)
  expect_equal(sum(colData(new_cds)$cluster_ext_type == "T cells"), 199)
})

colData(test_cds)$garnett_cluster <- c(rep(1, times=200), rep(2, times=200),
                                     rep(3, times=200), rep(4, times=200))
set.seed(10)
new_cds <- classify_cells(test_cds, test_classifier,
                          db = org.Hs.eg.db,
                          rank_prob_ratio = 1.5,
                          cluster_extend = TRUE,
                          cds_gene_id_type = "SYMBOL")

test_that("classify_cells works with provided clusters", {
  expect_identical(counts(new_cds), counts(test_cds))
  expect_identical(rowData(new_cds), rowData(test_cds))
  expect_equal(sum(colData(new_cds)$cell_type == "B cells"), 210)
  expect_equal(sum(colData(new_cds)$cell_type == "CD4 T cells"), 76)
  expect_equal(sum(colData(new_cds)$cell_type == "CD8 T cells"), 47)
  expect_equal(sum(colData(new_cds)$cell_type == "T cells"), 103)
  expect_equal(sum(colData(new_cds)$cluster_ext_type == "B cells"), 400)
  expect_equal(sum(colData(new_cds)$cluster_ext_type == "CD4 T cells"), 0)
  expect_equal(sum(colData(new_cds)$cluster_ext_type == "T cells"), 400)
})

test_that("classify_cells error messages work", {
  expect_error(new_cds <- classify_cells(test_cds, test_classifier,
                                         db = org.Hs.eg.db,
                                         rank_prob_ratio = 1.5,
                                         cluster_extend = TRUE),
               paste("None of the model genes are in your CDS object. Did you",
                     "specify the correct cds_gene_id_type and the correct",
                     "db?"))
  expect_error(suppressWarnings(new_cds <- classify_cells(test_cds,
                                                  test_classifier,
                                                  db = org.Mm.eg.db,
                                                  rank_prob_ratio = 1.5,
                                                  cluster_extend = TRUE,
                                                  cds_gene_id_type = "SYMBOL")),
               paste("None of the model genes are in your CDS object. Did you",
                     "specify the correct cds_gene_id_type and the correct",
                     "db?"))
})

set.seed(10)
new_cds <- classify_cells(test_cds, test_classifier,
                          db = org.Hs.eg.db,
                          rank_prob_ratio = 1.5,
                          cluster_extend = TRUE,
                          cds_gene_id_type = "SYMBOL",
                          cluster_extend_max_frac_unknown = .4,
                          cluster_extend_max_frac_incorrect = .5)

set.seed(10)
new_cds2 <- classify_cells(test_cds, test_classifier,
                          db = org.Hs.eg.db,
                          rank_prob_ratio = 1.5,
                          cluster_extend = TRUE,
                          cds_gene_id_type = "SYMBOL",
                          cluster_extend_max_frac_unknown = .95,
                          cluster_extend_max_frac_incorrect = .5)

test_that("classify_cells extension parameters work", {
  expect_identical(counts(new_cds), counts(test_cds))
  expect_identical(rowData(new_cds), rowData(test_cds))
  expect_equal(sum(colData(new_cds)$cell_type == "B cells"), 210)
  expect_equal(sum(colData(new_cds)$cell_type == "CD4 T cells"), 76)
  expect_equal(sum(colData(new_cds)$cell_type == "CD8 T cells"), 47)
  expect_equal(sum(colData(new_cds)$cell_type == "T cells"), 103)
  expect_equal(sum(colData(new_cds)$cluster_ext_type == "B cells"), 210)
  expect_equal(sum(colData(new_cds)$cluster_ext_type == "CD4 T cells"), 76)
  expect_equal(sum(colData(new_cds)$cluster_ext_type == "T cells"), 103)

  expect_identical(counts(new_cds2), counts(test_cds))
  expect_identical(rowData(new_cds2), rowData(test_cds))
  expect_equal(sum(colData(new_cds2)$cell_type == "B cells"), 210)
  expect_equal(sum(colData(new_cds2)$cell_type == "CD4 T cells"), 76)
  expect_equal(sum(colData(new_cds2)$cell_type == "CD8 T cells"), 47)
  expect_equal(sum(colData(new_cds2)$cell_type == "T cells"), 103)
  expect_equal(sum(colData(new_cds2)$cluster_ext_type == "B cells"), 400)
  expect_equal(sum(colData(new_cds2)$cluster_ext_type == "CD4 T cells"), 400)
  expect_equal(sum(colData(new_cds2)$cluster_ext_type == "T cells"), 0)
})

