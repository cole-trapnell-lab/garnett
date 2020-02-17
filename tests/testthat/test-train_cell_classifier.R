context("test-train_cell_classifier.R")

library(org.Mm.eg.db)
file_str = paste0(readChar("../logic_test.txt",
                           file.info("../logic_test.txt")$size),"\n")
parse_list <- garnett:::parse_input(file_str)
name_order <- unlist(parse_list[["name_order"]])
rm("name_order", envir=parse_list)

test_that("error message for convert ids works",{
  expect_error(garnett:::make_name_map(parse_list,
                                       c("App", "Itga2", "Ncam1", "Lyve1",
                                         "Tek", "Ncam3", "Kdr", "Ramp2",
                                         "Flt1", "App2", "Itga1", "Itga2",
                                         "Itga3"),
                                       "SYMBOL",
                                       "ENSEMBL",
                                       org.Mm.eg.db))
                 })

# Check gene names and keywords
gene_table <- garnett:::make_name_map(parse_list,
                            c("App", "Itga2", "Ncam1", "Lyve1", "Tek", "Ncam3",
                              "Kdr", "Ramp2", "Flt1", "App2", "Itga1", "Itga2",
                              "Itga3"),
                            "SYMBOL",
                            "SYMBOL",
                            org.Mm.eg.db)

test_that("make_name_map works", {
  expect_is(gene_table, "data.frame")
  gene_table <- suppressWarnings(garnett:::make_name_map(parse_list,
                                  c("ENSMUSG00000022892", "ENSMUSG00000015533",
                                    "ENSMUSG00000039542", "ENSMUSG00000030787",
                                    "ENSMUSG00000006386", "ENSMUSG00000062960",
                                    "ENSMUSG00000001240", "ENSMUSG00000029648",
                                    "ENSMUSG00000042284", "ENSMUSG00000015533",
                                    "ENSMUSG00000001507"),
                                    "ENSEMBL",
                                    "SYMBOL",
                                    org.Mm.eg.db))
  expect_equal(gene_table$fgenes[1], 'ENSMUSG00000001240')
  expect_equal(gene_table$orig_fgenes[1], "Ramp2")
  expect_equal(sum(gene_table$in_cds), 11)
})


logic_list <- garnett:::assemble_logic(parse_list[["Endothelial"]], gene_table)

test_that("logic is assembled correctly", {
  logic <- unlist(logic_list)
  expect_equal(logic[1], "(x['App',] > 1) & (x['App',] < 500)")
  expect_equal(logic[2], "(x['App2',] > 1) & (x['App2',] < 500)")
  expect_equal(logic[3], "(x['Itga1',] > 2) & (x['Itga1',] < Inf)")
  expect_equal(logic[4], "(x['Itga2',] > 2) & (x['Itga2',] < Inf)")
  expect_equal(logic[5], "(x['Itga3',] > 2) & (x['Itga3',] < Inf)")
  expect_equal(logic[6], "(x['Ncam1',] > -1) & (x['Ncam1',] < 100)")
  expect_equal(logic[7], "(x['Ncam3',] > -1) & (x['Ncam3',] < 100)")
  expect_equal(logic[8], "assigns == 'Endothelial'")
  expect_equal(logic[9], "mouse.sex %in% c('F')")
  expect_equal(logic[10], "tissue %in% c('liver', 'kidney', 'brain')")
})
library(org.Hs.eg.db)
data(test_cds)

set.seed(260)
test_classifier <- train_cell_classifier(cds = test_cds,
                                         marker_file = "../ensembl_test.txt",
                                         db=org.Hs.eg.db,
                                         min_observations = 10,
                                         cds_gene_id_type = "SYMBOL",
                                         num_unknown = 50,
                                         marker_file_gene_id_type = "ENSEMBL")

test_that("training works with ensembl ids", {
  expect_is(test_classifier, "garnett_classifier")
  expect_equal(length(test_classifier@classification_tree), 10)
})

data(test_cds)
ens_cds <- garnett:::cds_to_other_id(test_cds, db = org.Hs.eg.db,
                                     new_gene_id_type = "ENSEMBL",
                                     input_file_gene_id_type = "SYMBOL")

set.seed(260)
test_classifier <- train_cell_classifier(cds = ens_cds,
                                         marker_file = "../ensembl_test.txt",
                                         db=org.Hs.eg.db,
                                         min_observations = 10,
                                         cds_gene_id_type = "ENSEMBL",
                                         num_unknown = 50,
                                         marker_file_gene_id_type = "ENSEMBL")

test_that("training works with both ids", {
  expect_is(test_classifier, "garnett_classifier")
  expect_equal(length(test_classifier@classification_tree), 10)
})


set.seed(260)
test_classifier <- train_cell_classifier(cds = test_cds,
                                         marker_file = "../pbmc_test.txt",
                                         db=org.Hs.eg.db,
                                         min_observations = 10,
                                         cds_gene_id_type = "SYMBOL",
                                         num_unknown = 50,
                                         marker_file_gene_id_type = "SYMBOL")
test_that("training works with symbol ids", {
  expect_is(test_classifier, "garnett_classifier")
  expect_equal(length(test_classifier@classification_tree), 10)
})

set.seed(260)
test_classifier2 <- train_cell_classifier(cds = test_cds,
                                         marker_file = "../pbmc_test_dos.txt",
                                         db=org.Hs.eg.db,
                                         min_observations = 10,
                                         cds_gene_id_type = "SYMBOL",
                                         num_unknown = 50,
                                         marker_file_gene_id_type = "SYMBOL")
test_that("training works with carriage return", {
  expect_is(test_classifier, "garnett_classifier")
  expect_equal(length(test_classifier@classification_tree), 10)
})

test_that("circular subtypes not allowed", {
  expect_error(
    test_classifier2 <- train_cell_classifier(cds = test_cds,
                                              marker_file = "../bad_parse2.txt",
                                              db=org.Hs.eg.db,
                                              min_observations = 10,
                                              cds_gene_id_type = "SYMBOL",
                                              num_unknown = 50,
                                              marker_file_gene_id_type = "SYMBOL"),
    "'test cell 2' cannot be a subtype of itself. Please modify marker file.")
})
pData(test_cds)$overall_type <- ifelse(pData(test_cds)$FACS_type == "B cells", "B cells", "T cells")

test_classifier <- train_cell_classifier(cds = test_cds,
                                        marker_file = "../marker_free_test.txt",
                                        db=org.Hs.eg.db,
                                        min_observations = 10,
                                        cds_gene_id_type = "SYMBOL",
                                        num_unknown = 50,
                                        marker_file_gene_id_type = "SYMBOL")
