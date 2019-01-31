library(garnett)
library(org.Hs.eg.db)

# load in the data
# NOTE: the 'system.file' file name is only necessary to read in 
# included package data
#
mat <- read.table(system.file("extdata", "exprs.txt", package = "garnett"), 
                  check.names = FALSE)
fdata <- read.table(system.file("extdata", "fdata.txt", package = "garnett"))
pdata <- read.table(system.file("extdata", "pdata.txt", package = "garnett"), 
                    sep="\t")

# create a new CDS object
pd <- new("AnnotatedDataFrame", data = pdata)
fd <- new("AnnotatedDataFrame", data = fdata)
pbmc_cds <- newCellDataSet(as(as.matrix(mat), "dgCMatrix"),
                             phenoData = pd,
                             featureData = fd)

# generate size factors for normalization later
pbmc_cds <- estimateSizeFactors(pbmc_cds)

marker_file_path <- system.file("extdata", "pbmc_bad_markers.txt", 
                                package = "garnett")
marker_check <- check_markers(pbmc_cds, marker_file_path, 
                              db=org.Hs.eg.db, 
                              cds_gene_id_type = "SYMBOL", 
                              marker_file_gene_id_type = "SYMBOL")

png("~/Documents/Software/garnett_pages/images/marker_check.png", res=300, width = 5.5, height=5, units="in")
plot_markers(marker_check, amb_label_size = 3) 
dev.off()

set.seed(260)

marker_file_path <- system.file("extdata", "pbmc_test.txt", 
                                package = "garnett")
pbmc_classifier <- train_cell_classifier(cds = pbmc_cds,
                                         marker_file = marker_file_path,
                                         db=org.Hs.eg.db,
                                         cds_gene_id_type = "SYMBOL",
                                         num_unknown = 50,
                                         marker_file_gene_id_type = "SYMBOL")

pbmc_cds <- classify_cells(pbmc_cds, pbmc_classifier,
                           db = org.Hs.eg.db,
                           classify_communities = TRUE,
                           cds_gene_id_type = "SYMBOL")

library(ggplot2)

png("~/Documents/Software/garnett_pages/images/cell_type_ex.png", res=300, width = 5, height=4, units="in")
qplot(tsne_1, tsne_2, color = cell_type, data = pData(pbmc_cds)) + theme_bw()
dev.off()

png("~/Documents/Software/garnett_pages/images/community_type_ex.png", res=300, width = 5.5, height=4.4, units="in")
qplot(tsne_1, tsne_2, color = community_type, data = pData(pbmc_cds)) + theme_bw()
dev.off()

png("~/Documents/Software/garnett_pages/images/FACS_type_ex.png", res=300, width = 5.5, height=4.4, units="in")
qplot(tsne_1, tsne_2, color = FACS_type, data = pData(pbmc_cds)) + theme_bw()
dev.off()
