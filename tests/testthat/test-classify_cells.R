context("test-classify_cells.R")

library(org.Hs.eg.db)
library(org.Mm.eg.db)
data(test_classifier)
data(test_cds)

set.seed(10)
new_cds <- garnett:::get_communities(test_cds)

test_that("get_communities works", {
  skip_on_travis()
  expect_identical(exprs(new_cds), exprs(test_cds))
  expect_equal(ncol(pData(new_cds)) - 1, ncol(pData(test_cds)))
  expect_identical(fData(new_cds), fData(test_cds))
  expect_equal(length(unique(pData(new_cds)$louv_cluster)), 5)
  expect_equal(sum(pData(new_cds)$louv_cluster == 2), 172)
  expect_equal(as.character(pData(new_cds)$louv_cluster[5]), "2")
})

# classify cells
set.seed(10)
new_cds <- suppressWarnings(classify_cells(test_cds, test_classifier,
                          db = org.Hs.eg.db,
                          rank_prob_ratio = 1.5,
                          cluster_extend = TRUE,
                          cds_gene_id_type = "SYMBOL"))

test_that("classify_cells works", {
  expect_identical(exprs(new_cds), exprs(test_cds))
  expect_identical(fData(new_cds), fData(test_cds))
  expect_equal(sum(pData(new_cds)$cell_type == "B cells"), 160)
  expect_equal(sum(pData(new_cds)$cell_type == "CD4 T cells"), 82)
  expect_equal(sum(pData(new_cds)$cell_type == "CD8 T cells"), 48)
  expect_equal(sum(pData(new_cds)$cell_type == "T cells"), 142)
  expect_equal(sum(pData(new_cds)$cluster_ext_type == "B cells"), 401)
  expect_equal(sum(pData(new_cds)$cluster_ext_type == "CD4 T cells"), 0)
  expect_equal(sum(pData(new_cds)$cluster_ext_type == "T cells"), 399)
})

pData(test_cds)$garnett_cluster <- c(rep(1, times=200), rep(2, times=200),
                                     rep(3, times=200), rep(4, times=200))
set.seed(10)
new_cds <- suppressWarnings(classify_cells(test_cds, test_classifier,
                          db = org.Hs.eg.db,
                          rank_prob_ratio = 1.5,
                          cluster_extend = TRUE,
                          cds_gene_id_type = "SYMBOL"))

test_that("classify_cells works with provided clusters", {
  expect_identical(exprs(new_cds), exprs(test_cds))
  expect_identical(fData(new_cds), fData(test_cds))
  expect_equal(sum(pData(new_cds)$cell_type == "B cells"), 160)
  expect_equal(sum(pData(new_cds)$cell_type == "CD4 T cells"), 82)
  expect_equal(sum(pData(new_cds)$cell_type == "CD8 T cells"), 48)
  expect_equal(sum(pData(new_cds)$cell_type == "T cells"), 142)
  expect_equal(sum(pData(new_cds)$cluster_ext_type == "B cells"), 400)
  expect_equal(sum(pData(new_cds)$cluster_ext_type == "CD4 T cells"), 0)
  expect_equal(sum(pData(new_cds)$cluster_ext_type == "T cells"), 400)
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
  expect_identical(exprs(new_cds), exprs(test_cds))
  expect_identical(fData(new_cds), fData(test_cds))
  expect_equal(sum(pData(new_cds)$cell_type == "B cells"), 160)
  expect_equal(sum(pData(new_cds)$cell_type == "CD4 T cells"), 82)
  expect_equal(sum(pData(new_cds)$cell_type == "CD8 T cells"), 48)
  expect_equal(sum(pData(new_cds)$cell_type == "T cells"), 142)
  expect_equal(sum(pData(new_cds)$cluster_ext_type == "B cells"), 160)
  expect_equal(sum(pData(new_cds)$cluster_ext_type == "CD4 T cells"), 0)
  expect_equal(sum(pData(new_cds)$cluster_ext_type == "T cells"), 400)

  expect_identical(exprs(new_cds2), exprs(test_cds))
  expect_identical(fData(new_cds2), fData(test_cds))
  expect_equal(sum(pData(new_cds)$cell_type == "B cells"), 160)
  expect_equal(sum(pData(new_cds)$cell_type == "CD4 T cells"), 82)
  expect_equal(sum(pData(new_cds)$cell_type == "CD8 T cells"), 48)
  expect_equal(sum(pData(new_cds)$cell_type == "T cells"), 142)
  expect_equal(sum(pData(new_cds)$cluster_ext_type == "B cells"), 160)
  expect_equal(sum(pData(new_cds)$cluster_ext_type == "CD4 T cells"), 0)
  expect_equal(sum(pData(new_cds)$cluster_ext_type == "T cells"), 400)
})

