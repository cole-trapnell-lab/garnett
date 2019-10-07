cds_to_other_id <- function(cds,
                            db,
                            input_file_gene_id_type,
                            new_gene_id_type,
                            verbose = FALSE) {
  matrix <- exprs(cds)
  fdata <- fData(cds)

  new_g <- convert_gene_ids(row.names(fdata),
                            db,
                            input_file_gene_id_type,
                            new_gene_id_type)
  lstart <- length(new_g)
  new_g <- new_g[!is.na(new_g)]
  new_g <- new_g[!duplicated(new_g)]
  lend <- length(new_g)
  if((lstart-lend)/lstart > .7) warning(paste("More than 70% of IDs were lost",
                                              "when converting to",
                                              new_gene_id_type, "IDs. Did you",
                                              "specify the correct gene ID",
                                              "types and the correct db?"))
  if(verbose) message(paste("After converting CDS to", new_gene_id_type,"IDs,",
                            lstart - lend, "IDs were lost"))

  matrix <- matrix[names(new_g),]
  fdata <- fdata[names(new_g),, drop=FALSE]
  row.names(matrix) <- new_g
  row.names(fdata) <- new_g

  pd = new("AnnotatedDataFrame", data = pData(cds))
  fd = new("AnnotatedDataFrame", data = fdata)
  cds = suppressWarnings(newCellDataSet(matrix,
                       phenoData=pd,
                       featureData=fd,
                       expressionFamily=cds@expressionFamily,
                       lowerDetectionLimit=cds@lowerDetectionLimit))

  return(cds)
}

convert_gene_ids <- function(gene_list,
                             db,
                             start_type,
                             end_type) {

  tryCatch({suppressMessages(AnnotationDbi::mapIds(db, keys = gene_list,
                                         column = end_type, start_type))},
           error = function(e) {
             msg <- paste0("Garnett cannot convert the gene IDs using the ",
                         "db and types provided. Please check that your db, ",
                         "cds_gene_id_type and marker_file_gene_id_type ",
                         "parameters are correct. Please note that the ", "
                         cds_gene_id_type refers to the type of the ",
                         "row.names of the feature (gene) table in your cds. ",
                         "Conversion error: ", e)
             stop(msg)
           })
}



#' Extract feature genes
#'
#' Extract the genes chosen as features in cell type classification from a
#' trained garnett_classifier
#'
#' @param classifier Trained garnett_classifier - output from
#'  \code{\link{train_cell_classifier}}.
#' @param node Character. The name of the parent node of the multinomial
#'  classifier you would like to view features for. If top level, use "root".
#' @param convert_ids Logical. Should classifier IDs be converted to SYMBOL?
#' @param db Bioconductor AnnotationDb-class package for converting gene IDs.
#'  For example, for humans use org.Hs.eg.db. See available packages at
#'  \href{http://bioconductor.org/packages/3.8/data/annotation/}{Bioconductor}.
#'  If \code{convert_ids = FALSE}, db can be \code{NULL}.
#'
#' @return A data.frame of coefficient values for each gene with non-zero
#'  coefficients in the classifier.
#' @export
#'
#' @examples
#' library(org.Hs.eg.db)
#' data(test_classifier)
#' featuresdf <- get_feature_genes(test_classifier, db=org.Hs.eg.db)
#' featuresdf2 <- get_feature_genes(test_classifier,
#'                                  convert_ids = FALSE,
#'                                  node = "T cells")
#'
get_feature_genes <- function(classifier,
                              node = "root",
                              convert_ids = FALSE,
                              db=NULL) {
  assertthat::assert_that(is(classifier, "garnett_classifier"))
  assertthat::assert_that(is.character(node))
  assertthat::assert_that(is.logical(convert_ids))
  if (convert_ids) {
    if (is.null(db)) stop("If convert_ids = TRUE, db must be provided.")
    if (is(db, "character") && db == "none")
      stop("Cannot convert IDs if db = 'none'.")
    assertthat::assert_that(is(db, "OrgDb"),
                            msg = paste0("db must be an 'AnnotationDb' object ",
                                         "see http://bioconductor.org/",
                                         "packages/3.8/data/annotation/ ",
                                         "for available"))
  }
  s = "lambda.min"
  cvfit <- igraph::V(classifier@classification_tree)[node]$model
  feature_genes = stats::coef(cvfit[[1]], s = s)

  all <- as.data.frame(as.matrix(do.call("cbind", feature_genes)))
  names(all) <- names(feature_genes)
  zeroes <- all != 0
  all <- all[rowSums(zeroes) > 0,]

  if (.hasSlot(classifier, "gene_id_type")) {
    classifier_gene_id_type <- classifier@gene_id_type
  } else {
    classifier_gene_id_type <- "ENSEMBL"
  }

  if (convert_ids) {
    convs <- convert_gene_ids(row.names(all)[2:length(row.names(all))],
                              db,
                              classifier_gene_id_type,
                              "SYMBOL")
    convs[is.na(convs)] <- names(convs[is.na(convs)])
    if(sum(duplicated(convs)) > 0) {
      convs[duplicated(convs)] <- paste0(convs[duplicated(convs)], "_2")
    }
    row.names(all)[2:length(row.names(all))] <- convs
  }
  all
}


#' Retrieve marker references from garnett_classifier
#'
#' @param classifier garnett_classifier created using train_cell_classifier.
#' @param cell_type Cell type name or \code{NULL}. References for which cell
#'  type should be printed? If \code{NULL}, all are printed.
#'
#' @return List of references included when garnett_classifier was trained.
#' @export
#'
#' @examples
#' data(test_classifier)
#' get_classifier_references(test_classifier)
#'
get_classifier_references <- function(classifier,
                                      cell_type = NULL) {
  assertthat::assert_that(is(classifier, "garnett_classifier"))
  if (!is.null(cell_type)) {
    assertthat::assert_that(cell_type %in% names(classifier@references))
  }
  if(is.null(cell_type)) {
    return(classifier@references)
  } else {
    return(classifier@references[[cell_type]])
  }
}


#' Check marker file
#'
#' Check the markers chosen for the marker file and generate a table of useful
#' statistics. The output of this function can be fed into
#' \code{\link{plot_markers}} to generate a diagnostic plot.
#'
#' @param cds Input CDS object.
#' @param marker_file A character path to the marker file to define cell types.
#'  See details and documentation for \code{\link{Parser}} by running
#'  \code{?Parser} for more information.
#' @param db Bioconductor AnnotationDb-class package for converting gene IDs.
#'  For example, for humans use org.Hs.eg.db. See available packages at
#'  \href{http://bioconductor.org/packages/3.8/data/annotation/}{Bioconductor}.
#'  If your organism does not have an AnnotationDb-class database available,
#'  you can specify "none", however then Garnett will not check/convert gene
#'  IDs, so your CDS and marker file must have the same gene ID type.
#' @param cds_gene_id_type The type of gene ID used in the CDS. Should be one
#'  of the values in \code{columns(db)}. Default is "ENSEMBL". Ignored if
#'  db = "none".
#' @param marker_file_gene_id_type The type of gene ID used in the marker file.
#'  Should be one of the values in \code{columns(db)}. Default is "SYMBOL".
#'  Ignored if db = "none".
#' @param propogate_markers Logical. Should markers from child nodes of a cell
#'  type be used in finding representatives of the parent type? Should
#'  generally be \code{TRUE}.
#' @param use_tf_idf Logical. Should TF-IDF matrix be calculated during
#'  estimation? If \code{TRUE}, estimates will be more accurate, but
#'  calculation is slower with very large datasets.
#' @param classifier_gene_id_type The type of gene ID that will be used in the
#'  classifier. If possible for your organism, this should be "ENSEMBL", which
#'  is the default. Ignored if db = "none".
#'
#' @return Data.frame of marker check results.
#'
#' @details This function checks the chosen cell type markers in the marker
#'  file provided to ensure they are good candidates for use in classification.
#'  The function works by estimating which cells will be chosen given each
#'  marker gene and returning some statistics for each marker. Note that this
#'  function does not take into account meta data information when calculating
#'  statistics.
#'
#'  \describe{
#'  The output data.frame has several columns:
#'  \item{marker_gene}{Gene name as provided in the marker file}
#'  \item{ENSEMBL}{The corresponding ensembl ID derived from db conversion}
#'  \item{parent}{The parent cell type in the cell type hierarchy - 'root' if
#'  top level}
#'  \item{cell_type}{The cell type the marker belongs to}
#'  \item{in_cds}{Whether the marker is present in the CDS}
#'  \item{nominates}{The number of cells the marker is estimated to nominate to
#'  the cell type}
#'  \item{total_nominated}{The total number of cells nominated by all the
#'  markers for that cell type}
#'  \item{exclusion_dismisses}{The number of cells no longer nominated to the
#'  cell type if this marker is excluded (i.e. not captured by other markers
#'  for the cell type)}
#'  \item{inclusion_ambiguates}{How many cells become ambiguous (i.e. are
#'  nominated to multiple cell types) if this marker is included}
#'  \item{most_overlap}{The cell type that most often shares this marker (i.e.
#'  is the other side of the ambiguity). If inclusion_ambiguates is 0,
#'  most_overlap is NA}
#'  \item{ambiguity}{inclusion_ambiguates/nominates - if high, consider
#'  excluding this marker}
#'  \item{marker_score}{(1/(ambiguity + .01)) * nominates/total_nominated - a
#'  general measure of the quality of a marker. Higher is better}
#'  \item{summary}{A summary column that identifies potential problems with the
#'  provided markers}
#'  }
#'
#' @export
#'
#' @examples
#' library(org.Hs.eg.db)
#' data(test_cds)
#'
#' # generate size factors for normalization later
#' test_cds <- estimateSizeFactors(test_cds)
#' marker_file_path <- system.file("extdata", "pbmc_bad_markers.txt",
#'                                 package = "garnett")
#' marker_check <- check_markers(test_cds, marker_file_path,
#'                               db=org.Hs.eg.db,
#'                               cds_gene_id_type = "SYMBOL",
#'                               marker_file_gene_id_type = "SYMBOL")
#'
check_markers <- function(cds,
                          marker_file,
                          db,
                          cds_gene_id_type = "SYMBOL",
                          marker_file_gene_id_type = "SYMBOL",
                          propogate_markers = TRUE,
                          use_tf_idf = TRUE,
                          classifier_gene_id_type = "ENSEMBL") {

  ##### Check inputs #####
  assertthat::assert_that(is(cds, "CellDataSet"))
  assertthat::assert_that(assertthat::has_name(pData(cds), "Size_Factor"),
                          msg = paste("Must run estimateSizeFactors() on cds",
                                      "before calling check_markers"))
  assertthat::assert_that(sum(is.na(pData(cds)$Size_Factor)) == 0,
                          msg = paste("Must run estimateSizeFactors() on cds",
                                      "before calling check_markers"))
  assertthat::assert_that(is.character(marker_file))
  assertthat::is.readable(marker_file)

  if (is(db, "character") && db == "none") {
    cds_gene_id_type <- 'custom'
    classifier_gene_id_type <- 'custom'
    marker_file_gene_id_type <- 'custom'
  } else {
    assertthat::assert_that(is(db, "OrgDb"),
                            msg = paste0("db must be an 'AnnotationDb' object ",
                                         "or 'none' see ",
                                         "http://bioconductor.org/packages/",
                                         "3.8/data/annotation/ for available"))
    assertthat::assert_that(is.character(cds_gene_id_type))
    assertthat::assert_that(is.character(marker_file_gene_id_type))
    assertthat::assert_that(cds_gene_id_type %in% AnnotationDbi::keytypes(db),
                            msg = paste("cds_gene_id_type must be one of",
                                        "keytypes(db)"))
    assertthat::assert_that(marker_file_gene_id_type %in%
                              AnnotationDbi::keytypes(db),
                            msg = paste("marker_file_gene_id_type must be one",
                                        "of keytypes(db)"))
  }


  assertthat::assert_that(is.logical(propogate_markers))
  assertthat::assert_that(is.logical(use_tf_idf))

  ##### Set internal parameters #####
  back_cutoff <- 0.25

  sf <- pData(cds)$Size_Factor

  ##### Normalize and rename CDS #####
  if (!is(exprs(cds), "dgCMatrix")) {
    pd <- new("AnnotatedDataFrame", data = pData(cds))
    fd <- new("AnnotatedDataFrame", data = fData(cds))
    cds <- suppressWarnings(newCellDataSet(as(exprs(cds), "dgCMatrix"),
                          phenoData = pd,
                          featureData = fd))
    pData(cds)$Size_Factor <- sf
  }


  if(cds_gene_id_type != classifier_gene_id_type)  {
    cds <- cds_to_other_id(cds, db=db, cds_gene_id_type,
                           classifier_gene_id_type)
    pData(cds)$Size_Factor <- sf
  }
  pData(cds)$num_genes_expressed <- Matrix::colSums(as(exprs(cds), "lgCMatrix"))
  cell_totals <-  Matrix::colSums(exprs(cds))


  pd <- new("AnnotatedDataFrame", data = pData(cds))
  fd <- new("AnnotatedDataFrame", data = fData(cds))
  orig_cds <- cds
  temp <- exprs(cds)
  temp@x <- temp@x / rep.int(pData(cds)$Size_Factor, diff(temp@p))
  cds <- suppressWarnings(newCellDataSet(temp,
                        phenoData = pd, featureData = fd))

  pData(cds)$Size_Factor <- sf

  ##### Parse Marker File #####
  file_str = paste0(readChar(marker_file, file.info(marker_file)$size),"\n")

  parse_list <- parse_input(file_str)
  orig_name_order <- unlist(parse_list[["name_order"]])
  rm("name_order", envir=parse_list)
  if(is.null(parse_list)) stop("Parse failed!")
  message(paste("There are", length(parse_list), "cell type definitions"))


  ranks <- lapply(orig_name_order, function(i) parse_list[[i]]@parenttype)
  names(ranks) <- orig_name_order
  if(length(unlist(unique(ranks[which(!ranks %in% names(ranks) & lengths(ranks) != 0L)])) != 0)) {
    stop(paste("Subtype", unlist(unique(ranks[which(!ranks %in% names(ranks) & lengths(ranks) != 0L)])), "is not defined in marker file."))
  }

  name_order <- names(ranks[lengths(ranks) == 0L])
  ranks <- ranks[!names(ranks) %in% name_order]
  while(length(ranks) != 0) {
    name_order <- c(name_order, names(ranks)[ranks %in% name_order])
    ranks <- ranks[!names(ranks) %in% name_order]
  }


  # Check gene names and keywords
  gene_table <- check_marker_conversion(parse_list,
                                        as.character(row.names(fData(cds))),
                                        classifier_gene_id_type,
                                        marker_file_gene_id_type,
                                        db)
  gene_table$nominates <- NA
  gene_table$exclusion_dismisses <- NA
  gene_table$inclusion_ambiguates <- NA
  gene_table$most_overlap <- NA
  gene_table$total_nominated <- NA

  ##### Make garnett_classifier #####
  classifier <- new_garnett_classifier()

  other_rules <- c()
  for(i in name_order) {
    logic_list <- assemble_logic(parse_list[[i]], gene_table)
    classifier <- add_cell_rule(parse_list[[i]], classifier, logic_list)
    other_rules_mult <- get_rule_multiplier(i, classifier, orig_cds)
    other_rules <- c(other_rules, other_rules_mult)
  }

  names(other_rules) <- name_order

  if(propogate_markers) {
    root <- propogate_func(curr_node = "root", parse_list, classifier)
    gene_table <- check_marker_conversion(parse_list,
                                          as.character(row.names(fData(cds))),
                                          classifier_gene_id_type,
                                          marker_file_gene_id_type,
                                          db)
    gene_table$nominates <- NA
    gene_table$exclusion_dismisses <- NA
    gene_table$inclusion_ambiguates <- NA
    gene_table$most_overlap <- NA
    gene_table$total_nominated <- NA
  }

  if(use_tf_idf) {
    dat <- tfidf(cds)
  } else{
    dat <- Matrix::t(exprs(cds))
  }

  ##### For each node #####
  for (v in igraph::V(classifier@classification_tree)){
    child_cell_types <-
      igraph::V(classifier@classification_tree)[suppressWarnings(outnei(v))]$name

    if(length(child_cell_types) == 0) next

    ##### For each child #####
    all_types <- list()
    for (i in child_cell_types) {
      agg <- aggregate_positive_markers(parse_list[[i]], dat,
                                        gene_table, back_cutoff, agg = F)
      if (!is.null(agg)){
        agg[agg > 0] <- 1
        agg <- sweep(agg,MARGIN=1,other_rules[[i]][,1],`*`)
        all_types <- c(all_types, agg)
        names(all_types)[length(all_types)] <- i
      }
    }

    other_cells <- lapply(all_types, function(x) {
      rs <- Matrix::rowSums(x)
      rs[rs > 0] <- 1
      rs
    })
    amb_cells <- Reduce(`+`, other_cells)
    total_amb <- sum(amb_cells[amb_cells > 1]-1)
    assigned <- lapply(other_cells, function(x) names(x[x!=0]))

    for(cell_type in names(all_types)) {
      for(cols in colnames(all_types[[cell_type]])) {
        temp_other <- other_cells
        r <- which(gene_table$fgenes == cols &
                     gene_table$cell_type == cell_type)
        gene_table$nominates[r] <- sum(all_types[[cell_type]][,cols])
        if(length(colnames(all_types[[cell_type]])) == 1) {
          gene_table$exclusion_dismisses[r] <- gene_table$nominates[r]
          gene_table$total_nominated[r] <- gene_table$nominates[r]
          temp_other <- temp_other[setdiff(names(temp_other), cell_type)]
          temp <- Reduce(`+`, temp_other)
          new_amb <- sum(temp[temp > 1]-1)
          gene_table$inclusion_ambiguates[r] <- total_amb - new_amb
          ambs <- names(which(amb_cells != temp))
          amb_count <- lapply(assigned, function(x) sum(ambs %in% x))
          amb_count[cell_type] <- 0
          gene_table$most_overlap[r] <- names(which.max(amb_count))
        } else{
          total_assign <- sum(temp_other[[cell_type]] > 0)
          gene_table$total_nominated[r] <- total_assign
          temp_other[[cell_type]] <- Matrix::rowSums(
            all_types[[cell_type]][,setdiff(colnames(all_types[[cell_type]]),
                                            cols), drop=FALSE])
          gene_table$exclusion_dismisses[r] <-
            total_assign - sum(temp_other[[cell_type]] > 0)
          temp_other[[cell_type]][temp_other[[cell_type]] > 0] <- 1
          temp <- Reduce(`+`, temp_other)
          new_amb <- sum(temp[temp > 1]-1)
          gene_table$inclusion_ambiguates[r] <- total_amb - new_amb
          ambs <- names(which(amb_cells != temp))
          amb_count <- lapply(assigned, function(x) sum(ambs %in% x))
          amb_count[cell_type] <- 0
          gene_table$most_overlap[r] <- names(which.max(amb_count))
        }
      }
    }
  }

  names(gene_table) <- c("gene_id", "parent", "cell_type", "marker_gene",
                         "in_cds", "nominates", "exclusion_dismisses",
                         "inclusion_ambiguates", "most_overlap",
                         "total_nominated")

  gene_table <- gene_table[,c("marker_gene", "gene_id", "parent", "cell_type",
                              "in_cds", "nominates", "total_nominated",
                              "exclusion_dismisses", "inclusion_ambiguates",
                              "most_overlap")]
  gene_table$most_overlap[gene_table$inclusion_ambiguates == 0] <- NA
  gene_table$ambiguity <-
    gene_table$inclusion_ambiguates/gene_table$nominates
  gene_table$marker_score <- (1/(gene_table$ambiguity + .01)) *
    gene_table$nominates/gene_table$total_nominated
  gene_table$summary <- NA

  gene_table$summary[is.na(gene_table$gene_id)] <- "Not in db"
  gene_table$summary[is.na(gene_table$summary) &
                       !gene_table$in_cds] <- "Not in CDS"
  gene_table$summary[is.na(gene_table$summary) &
                       gene_table$ambiguity > 0.25] <- "High ambiguity?"
  gene_table$summary[is.na(gene_table$summary) &
                       gene_table$nominates <
                       stats::quantile(gene_table$nominates,
                                       .1, na.rm=TRUE)] <- "Low nomination?"
  gene_table$summary[is.na(gene_table$summary)] <- "Ok"
  gene_table$summary[gene_table$summary == "Ok" & is.na(gene_table$marker_score)] <- "Ok - marker_score N/A"
  return(gene_table)
}

check_marker_conversion <- function(parse_list,
                                    possible_genes,
                                    cds_gene_id_type,
                                    marker_file_gene_id_type,
                                    db) {
  gene_start <- collect_gene_names(parse_list)
  gene_table <- data.frame(fgenes = gene_start[,1], parent = gene_start[,2],
                           cell_type = gene_start[,3])
  gene_table$parent <- as.character(gene_table$parent)
  gene_table$fgenes <- as.character(gene_table$fgenes)
  gene_table$cell_type <- as.character(gene_table$cell_type)
  gene_table$orig_fgenes <- gene_table$fgenes
  if(cds_gene_id_type != marker_file_gene_id_type) {
    gene_table$fgenes <- convert_gene_ids(gene_table$orig_fgenes,
                                          db,
                                          marker_file_gene_id_type,
                                          cds_gene_id_type)
    bad_convert <- sum(is.na(gene_table$fgenes))
  }

  gene_table$in_cds <- gene_table$fgenes %in% possible_genes
  gene_table$in_cds[is.na(gene_table$in_cds)] <- FALSE

  gene_table$fgenes <- as.character(gene_table$fgenes)
  gene_table
}

get_rule_multiplier <- function(i, classifier, orig_cds) {
  ##### Exclude possibles using other definitions #####
  cell_class_func <-
    igraph::V(classifier@classification_tree)[i]$classify_func[[1]]

  parent <- environment(cell_class_func)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent=parent)

  pData(orig_cds)$assigns <- igraph::V(classifier@classification_tree)[i]$name
  Biobase::multiassign(names(pData(orig_cds)), pData(orig_cds), envir=e1)
  environment(cell_class_func) <- e1

  type_res <- cell_class_func(exprs(orig_cds))
  if (length(type_res)!= ncol(orig_cds)){
    message(paste("Error: classification function for",
                  igraph::V(classifier@classification_tree)[i]$name,
                  "returned a malformed result."))
    stop()
  }

  type_res <- as(as(type_res,"sparseVector"), "sparseMatrix")
  row.names(type_res) <- row.names(pData(orig_cds))
  colnames(type_res) <- i
  type_res
}



#' Plot marker metrics
#'
#' This function takes as input the output of the \code{\link{check_markers}}
#' function and generates a plot to visualize the most important metrics.
#'
#' @param marker_check_df Marker check data.frame - output of check_markers.
#' @param amb_marker_cutoff Numeric. Cutoff at which to label ambiguous markers.
#'  Default 0.5.
#' @param label_size Numeric, size of the text labels for ambiguous markers and
#' unplotted markers.
#'
#' @return A ggplot object of the marker plot.
#' @export
#'
#' @examples
#' library(org.Hs.eg.db)
#'
#' marker_file_path <- system.file("extdata", "pbmc_test.txt",
#'                                 package = "garnett")
#' data(test_cds)
#' marker_check <- check_markers(test_cds,
#'                               marker_file_path,
#'                               db=org.Hs.eg.db,
#'                               cds_gene_id_type = "SYMBOL",
#'                               marker_file_gene_id_type = "SYMBOL")
#'
#' plot_markers(marker_check)
#'
#'
#'
plot_markers <- function(marker_check_df,
                         amb_marker_cutoff = .5,
                         label_size = 2) {
  assertthat::assert_that(is.data.frame(marker_check_df))
  assertthat::assert_that(sum(!c("marker_gene", "cell_type", "nominates",
                                 "total_nominated", "most_overlap", "ambiguity",
                                 "marker_score", "summary") %in%
                                names(marker_check_df)) == 0,
                          msg = paste("marker_check_df must be the output of",
                                      "the check_markers function. Input is",
                                      "missing key columns"))
  labeldf <- data.frame(cell_type = marker_check_df$cell_type,
                        cell_type_label = paste0(marker_check_df$cell_type,
                                                 ": ",
                                                 marker_check_df$total_nominated),
                        total_nominated = marker_check_df$total_nominated)
  labeldf <- labeldf[order(labeldf$total_nominated),]
  labeldf <- labeldf[!duplicated(labeldf$cell_type),]
  labeldf$cell_type <- as.factor(labeldf$cell_type)
  marker_check_df$cell_type <- as.factor(marker_check_df$cell_type)
  marker_check_df$cell_type <-
    factor(marker_check_df$cell_type,
           labels = as.character(
             labeldf[order(labeldf$cell_type),]$cell_type_label))
  marker_check_df$marker_gene <- as.factor(marker_check_df$marker_gene)
  marker_check_df$tempy <- as.factor(paste(marker_check_df$marker_gene,
                                           marker_check_df$cell_type, sep=">"))
  ggplot2::ggplot(marker_check_df,
                  ggplot2::aes(x=ambiguity,
                               y=forcats::fct_reorder2(tempy,
                                                       cell_type,
                                                       -marker_score),
                               fill=100 * nominates/total_nominated)) +
    ggplot2::geom_point(ggplot2::aes(size = marker_score), color = "black",
                        pch=21, data = marker_check_df, stroke = .1) +
    ggplot2::geom_text(ggplot2::aes(x = 0.4,
                                    label = ifelse(is.na(ambiguity),
                                                   as.character(summary), '')),
                       color = "firebrick4", size = label_size) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label=ifelse(ambiguity > amb_marker_cutoff,
                                as.character(paste0("High overlap with\n",
                                                    most_overlap)),'')),
      color = "black", point.padding = .01, size=label_size,
      min.segment.length = ggplot2::unit(0, 'lines'), segment.size = .2) +
    ggplot2::facet_grid(cell_type~., scales="free", space="free_y",
                        labeller=ggplot2::label_value) +
    ggplot2::scale_size(name = "Marker\nscore", range = c(1,3)) +
    viridis::scale_fill_viridis(position="top", name="% of\nassigned") +
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = 0.5,
                                                    barheight = 4)) +
    ggplot2::xlim(c(0,1)) +
    ggplot2::ylab("") +
    ggplot2::xlab("Ambiguity") +
    ggplot2::annotate("segment", x = -Inf, xend = Inf, y = Inf, yend =Inf,
                      color = "grey93") +
    ggplot2::annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,
                      color = "grey93") +
    ggplot2::theme_classic() + #base_size=6
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0),
                   strip.background = ggplot2::element_rect(color = "grey93",
                                                            fill="grey93"),
                   legend.key.height = ggplot2::unit(0.2, "cm")) +
    ggplot2::scale_y_discrete(labels = function(x) {
      stringr::str_split_fixed(x, ">", 2)[,1]
    })

}

