context("test-get_training_sample")

library(org.Mm.eg.db)
library(org.Hs.eg.db)
data("test_cds")

set.seed(260)
lams <- unique(c(100000, 50000,
                 seq(10000, 100, by=-200),
                 seq(100,10, by=-20),
                 seq(10, 1, by=-2),
                 seq(1, .1, by=-.2), .1, .05, .01, .005, .001, .0005, .0001))
test_that("right error message when no cells picked", {
  expect_error(test_classifier <- train_cell_classifier(cds = test_cds,
                                   marker_file = "../pbmc_get_training_bad.txt",
                                   db=org.Hs.eg.db,
                                   min_observations = 10,
                                   cds_gene_id_type = "SYMBOL",
                                   num_unknown = 50,
                                   marker_file_gene_id_type = "SYMBOL"),
               paste("Not enough training samples for any cell types at root",
                     "of cell type hierarchy!"))
  expect_message(test_classifier <- train_cell_classifier(cds = test_cds,
                                  marker_file = "../pbmc_get_training_bad2.txt",
                                  db=org.Hs.eg.db,
                                  min_observations = 10,
                                  cds_gene_id_type = "SYMBOL",
                                  num_unknown = 50,
                                  marker_file_gene_id_type = "SYMBOL"),
                 paste("Not enough training samples for children of T cells.",
                       "They will not be subclassified."))
})

set.seed(260)
test_classifier <- train_cell_classifier(cds = test_cds,
                                      marker_file = "../pbmc_get_training.txt",
                                      db=org.Hs.eg.db,
                                      min_observations = 10,
                                      cds_gene_id_type = "SYMBOL",
                                      num_unknown = 50,
                                      marker_file_gene_id_type = "SYMBOL")

test_that("get_training_sample follows rules", {
  expect_equal(
    igraph::V(test_classifier@classification_tree)[1]$model[[1]]$glmnet.fit$nobs,
    298)
})
