context("test-utils.R")

library(org.Hs.eg.db)
library(org.Mm.eg.db)
data(test_cds)

ensembl <- garnett:::cds_to_other_id(test_cds, db=org.Hs.eg.db,
                                     "SYMBOL", "ENSEMBL")

test_that("cds_to_other_id works", {
  expect_is(ensembl, "CellDataSet")
  expect_equal(nrow(fData(ensembl)), 19149)
  expect_equal(row.names(fData(ensembl)[
    fData(ensembl)$gene_short_name == "FAM87B",]), "ENSG00000177757")
  expect_identical(exprs(test_cds)["FAM87B",],
                   exprs(ensembl)["ENSG00000177757",])
})

data(test_classifier)

out <- get_feature_genes(test_classifier, db=org.Hs.eg.db)
out2 <- get_feature_genes(test_classifier, convert_ids = F, db=org.Hs.eg.db)
out3 <- get_feature_genes(test_classifier, db=org.Hs.eg.db, node="T cells")

test_that("get_feature_genes works", {
  expect_error(get_feature_genes(test_classifier, convert_ids = TRUE),
               "If convert_ids = TRUE, db must be provided.")
  expect_is(out, "data.frame")
  expect_is(out2, "data.frame")
  expect_is(out3, "data.frame")
  expect_equal(nrow(out), 31)
  expect_equal(ncol(out3), 3)
  expect_equal(row.names(out2)[2], "ENSG00000196154")
})

test_that("get_classifier_references works", {
  out <- get_classifier_references(test_classifier, cell_type = NULL)
  expect_equal(length(out), 2)
  expect_equal(out[["T cells"]],
               "https://www.ncbi.nlm.nih.gov/pubmed/?term=1534551")
  out <- get_classifier_references(test_classifier, cell_type = "B cells")
  expect_equal(length(out), 2)
  expect_equal(out[[1]],
               "https://www.ncbi.nlm.nih.gov/pubmed/?term=23688120")
})


file_str = paste0(readChar("../logic_test.txt",
                           file.info("../logic_test.txt")$size),"\n")
parse_list <- garnett:::parse_input(file_str)
name_order <- unlist(parse_list[["name_order"]])
rm("name_order", envir=parse_list)

# Check gene names and keywords
gene_table <- garnett:::check_marker_conversion(parse_list,
                                                c("App", "Itga2", "Ncam1",
                                                  "Lyve1", "Tek", "Ncam3",
                                                  "Kdr", "Ramp2", "Flt1",
                                                  "App2", "Itga1", "Itga2",
                                                  "Itga3"),
                                                "ENSEMBL", "SYMBOL",
                                                org.Mm.eg.db)

test_that("check_marker_conversion works", {
  expect_equal(nrow(gene_table), 12)
  expect_equal(gene_table$fgenes[1], "ENSMUSG00000030787")
  expect_equal(sum(gene_table$in_cds), 0)
  gene_table <- garnett:::check_marker_conversion(parse_list,
                                                  c("App", "Itga2", "Ncam1",
                                                    "Lyve1", "Tek", "Ncam3",
                                                    "Kdr", "Ramp2", "Flt1",
                                                    "App2", "Itga1", "Itga2"),
                                                  "SYMBOL", "SYMBOL",
                                                  org.Mm.eg.db)
  expect_equal(sum(gene_table$in_cds), 11)
})


marker_check <- check_markers(test_cds,
                              "../pbmc_test.txt",
                              db=org.Hs.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")

marker_check_dos <- check_markers(test_cds,
                              "../pbmc_test_dos.txt",
                              db=org.Hs.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")

marker_check2 <- check_markers(test_cds, use_tf_idf = F,
                              "../pbmc_1mark.txt",
                              db=org.Hs.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")

test_that("check_markers works", {
  expect_identical(marker_check, marker_check_dos)
  expect_equal(sum(marker_check$marker_score), 413.8816, tol = 1e-4)
  expect_equal(nrow(marker_check), 18)
  expect_equal(sum(marker_check$summary != "Ok"), 3)
  sub <- subset(marker_check, parent != "root")
  expect_equal(sum(sub$nominates), 323)
  expect_equal(sum(sub$exclusion_dismisses), 233)
  expect_equal(sum(sub$inclusion_ambiguates), 63)
  expect_equal(sub$most_overlap[5], "CD4 T cells")
  expect_equal(sum(marker_check2$marker_score, na.rm=T), 388.5029, tol = 1e-4)
  expect_equal(marker_check2[18, "total_nominated"], 46)
})

test_that("plot_markers works", {
  skip_on_travis()
  vdiffr::expect_doppelganger("basic_marker_plot",
                              fig = plot_markers(marker_check))
  vdiffr::expect_doppelganger("basic_marker_plot change amb",
                              fig = plot_markers(marker_check,
                                                 amb_marker_cutoff = .7))
})

marker_check <- check_markers(test_cds,
                              "../pbmc_test.txt",
                              db="none",
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")

marker_check2 <- check_markers(test_cds, use_tf_idf = F,
                               "../pbmc_1mark.txt",
                               db="none",
                               cds_gene_id_type = "SYMBOL",
                               marker_file_gene_id_type = "SYMBOL")

test_that("check_markers works db='none'", {
  expect_equal(sum(marker_check$marker_score), 414.3958, tol = 1e-3)
  expect_equal(nrow(marker_check), 18)
  expect_equal(sum(marker_check$summary != "Ok"), 3)
  sub <- subset(marker_check, parent != "root")
  expect_equal(sum(sub$nominates), 323)
  expect_equal(sum(sub$exclusion_dismisses), 233)
  expect_equal(sum(sub$inclusion_ambiguates), 63)
  expect_equal(sub$most_overlap[5], "CD4 T cells")
  expect_equal(sum(marker_check2$marker_score, na.rm=T), 388.5029, tol = 1e-4)
  expect_equal(marker_check2[18, "total_nominated"], 46)
})

