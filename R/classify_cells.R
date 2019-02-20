#' Classify cells from trained garnett_classifier
#'
#' This function uses a previously trained \code{\link{garnett_classifier}}
#' (trained using \code{\link{train_cell_classifier}}) to classify cell types
#' in a CDS object.
#'
#' @param cds Input CDS object.
#' @param classifier Trained garnett_classifier - output from
#'  \code{\link{train_cell_classifier}}.
#' @param db Bioconductor AnnotationDb-class package for converting gene IDs.
#'  For example, for humans use org.Hs.eg.db. See available packages at
#'  \href{http://bioconductor.org/packages/3.8/data/annotation/}{Bioconductor}.
#'  If your organism does not have an AnnotationDb-class database available,
#'  you can specify "none", however then Garnett will not check/convert gene
#'  IDs, so your CDS and marker file must have the same gene ID type.
#' @param cds_gene_id_type The type of gene ID used in the CDS. Should be one
#'  of the values in \code{columns(db)}. Default is "ENSEMBL". Ignored if
#'  db = "none".
#' @param rank_prob_ratio Numeric value greater than 1. This is the minimum
#'  odds ratio between the probability of the most likely cell type to the
#'  second most likely cell type to allow assignment. Default is 1.5. Higher
#'  values are more conservative.
#' @param cluster_extend Logical. When \code{TRUE}, the classifier
#'  provides a secondary cluster-extended classification, which assigns type
#'  for the entire cluster based on the assignments of the cluster members. If
#'  the pData table of the input CDS has a column called "garnett_cluster",
#'  this will be used for cluster-extended assignments. Otherwise, assignments
#'  are calculated using Louvain community detection in PCA space. This
#'  assignment is returned as a column in the output CDS pData table. For large
#'  datasets, if the "garnett_cluster" column is not provided and
#'  \code{cluster_extend = TRUE}, the function can be significantly slower the
#'  first time it is run. See details for more information.
#' @param verbose Logical. Should progress messages be printed.
#'
#' @details This function applies a previously trained multinomial glmnet
#'  classifier at each node of a previously defined garnett_classifier tree.
#'  The output is a CDS object with cell type classifications added to the
#'  pData table.
#'
#'  When \code{cluster_extend = TRUE}, louvain communities are calculated in
#'  PCA space. Any cluster where >90% of classified cells are of a single type,
#'  >5% of cells are classified, and a minimum of 5 cells are classified will
#'  be assigned that cluster-extended type. Both cluster-extended type and
#'  originally calculated cell type are reported.
#'
#' @return CDS object with classifications in the \code{pData} table.
#' @export
#'
#' @examples
#' library(org.Hs.eg.db)
#' data(test_classifier)
#' data(test_cds)
#'
#' # classify cells
#' test_cds <- classify_cells(test_cds, test_classifier,
#'                            db = org.Hs.eg.db,
#'                            rank_prob_ratio = 1.5,
#'                            cluster_extend = TRUE,
#'                            cds_gene_id_type = "SYMBOL")
#'
classify_cells <- function(cds,
                           classifier,
                           db,
                           cds_gene_id_type = "ENSEMBL",
                           rank_prob_ratio = 1.5,
                           cluster_extend = FALSE,
                           verbose = FALSE) {
  if(verbose) message("Starting classification")
  ##### Check inputs #####
  if(verbose) message("Checking inputs")
  assertthat::assert_that(is(cds, "CellDataSet"))
  assertthat::assert_that(assertthat::has_name(pData(cds), "Size_Factor"),
                          msg = paste("Must run estimateSizeFactors() on cds",
                                      "before calling classify_cells"))
  assertthat::assert_that(sum(is.na(pData(cds)$Size_Factor)) == 0,
                          msg = paste("Must run estimateSizeFactors() on cds",
                                      "before calling classify_cells"))
  assertthat::assert_that(is(classifier, "garnett_classifier"))
  if(is(db, "character") && db == "none") {
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
    assertthat::assert_that(cds_gene_id_type %in% AnnotationDbi::keytypes(db),
                            msg = paste("cds_gene_id_type must be one of",
                                        "keytypes(db)"))
  }

  assertthat::assert_that(is.numeric(rank_prob_ratio))
  assertthat::assert_that(rank_prob_ratio > 1,
                          msg = "rank_prob_ratio must be greater than 1")
  assertthat::assert_that(is.logical(cluster_extend))
  assertthat::assert_that(is.logical(verbose))

  ##### Set internal parameters #####
  s <- "lambda.min"

  ##### Normalize CDS #####
  orig_cds <- cds

  if(verbose) message("Normalizing CDS object\n")

  if (!is(exprs(cds), "dgCMatrix")) {
    sf <- pData(cds)$Size_Factor
    pd <- new("AnnotatedDataFrame", data = pData(cds))
    fd <- new("AnnotatedDataFrame", data = fData(cds))
    cds <- suppressWarnings(newCellDataSet(as(exprs(cds), "dgCMatrix"),
                          phenoData = pd,
                          featureData = fd))
    pData(cds)$Size_Factor <- sf
  }

  pData(cds)$num_genes_expressed <- Matrix::colSums(as(exprs(cds),
                                                       "lgCMatrix"))
  new_cell_totals <- Matrix::colSums(exprs(cds))
  sfs <- new_cell_totals/(classifier@cell_totals *
                            stats::median(pData(cds)$num_genes_expressed))
  sfs[is.na(sfs)] <- 1
  save_sf <- pData(cds)$Size_Factor
  pData(cds)$Size_Factor <- sfs
  pd <- new("AnnotatedDataFrame", data = pData(cds))
  fd <- new("AnnotatedDataFrame", data = fData(cds))
  temp <- exprs(cds)
  temp@x <- temp@x / rep.int(pData(cds)$Size_Factor, diff(temp@p))
  norm_cds <- suppressWarnings(newCellDataSet(temp,
                             phenoData = pd, featureData = fd))

  if (.hasSlot(classifier, "gene_id_type")) {
    classifier_gene_id_type <- classifier@gene_id_type
  } else {
    classifier_gene_id_type <- "ENSEMBL"
  }

  ### Convert to Classifier IDs ###
  if(cds_gene_id_type != classifier_gene_id_type) {
    if (verbose) message(paste("Converting CDS IDs to",
                               classifier_gene_id_type, "\n"))
    lstart <- nrow(fData(norm_cds))
    norm_cds <- cds_to_other_id(norm_cds,
                                db=db,
                                cds_gene_id_type,
                                classifier_gene_id_type,
                                verbose = FALSE)
    lend <- nrow(fData(norm_cds))
  }

  pData(norm_cds)$Size_Factor <- sfs
  cds <- orig_cds

  ##### Calculate cell communities #####
  if (cluster_extend) {
    if ("garnett_cluster" %in% names(pData(cds))) {
      pData(norm_cds)$louv_cluster <- pData(cds)$garnett_cluster
    } else {
      if(verbose) message(paste("No garnett_cluster column provided,",
                                "generating clusters for classification\n"))
      norm_cds <- get_communities(norm_cds)
      pData(cds)$garnett_cluster <- pData(norm_cds)$louv_cluster
    }
  }

  ##### Classify cells #####
  if(verbose) message("Predicting cell types\n")
  class_df <- run_classifier(classifier, norm_cds,
                             cluster_extend = cluster_extend,
                             s=s,
                             rank_prob_ratio = rank_prob_ratio)

  pData(cds)$cell_type <- NULL
  if("cluster_ext_type" %in% names(pData(cds)))
    pData(cds)$cluster_ext_type <- NULL

  pData(cds)$Size_Factor <- save_sf
  pData(cds) <- cbind(pData(cds), class_df)
  if(verbose) message("Complete!\n")
  cds
}

run_classifier <- function(classifier,
                           cds,
                           cluster_extend,
                           rank_prob_ratio,
                           s) {

  imputed_gate_res <- list()

  for (v in igraph::V(classifier@classification_tree)){

    child_cell_types <- igraph::V(classifier@classification_tree)[
      suppressWarnings(outnei(v))]$name

    if (length(child_cell_types) > 0){
      new_gate_res <- make_predictions(cds, classifier, v,
                                       rank_prob_ratio = rank_prob_ratio,
                                       s=s)
      imputed_gate_res <- append(imputed_gate_res, new_gate_res)
    }
  }

  assignments <- rep("Unknown", length(imputed_gate_res[[1]]))
  names(assignments) <- row.names(imputed_gate_res[[1]])

  tree_levels <- igraph::distances(classifier@classification_tree,
                                   to = "root")[,"root"]
  tree_depth <- max(tree_levels)
  level_table <- data.frame(cell = row.names(imputed_gate_res[[1]]),
                            level1 = "Unknown" ,
                            stringsAsFactors = FALSE)

  fill_in_assignments <- function(curr_assignments, classifier, v,
                                  imputed_gate_res, level_table){

    for (child in igraph::V(classifier@classification_tree)[
      suppressWarnings(outnei(v))]){
      curr_level <- paste0("level", tree_levels[child])
      if(!curr_level %in% names(level_table)) {
        level_table[[curr_level]] <- "Unknown"
      }
      all_parents <- igraph::V(classifier@classification_tree)[
        igraph::all_simple_paths(classifier@classification_tree, v, to = child,
                                 mode = "out")[[1]]]$name
      parents <- setdiff(all_parents, "root")
      if (length(intersect(parents, names(imputed_gate_res))) > 0 &
          sum(all_parents %in% union(names(imputed_gate_res), "root")) > 1) {
        type_res <- imputed_gate_res[parents]
        if (length(type_res) > 1 ){
          mat <- do.call(cbind,type_res)
          type_res <- apply(mat, 1, function(x) { prod(x) })
          cell_type <- igraph::V(classifier@classification_tree) [ child ]$name
        }else{
          cell_type <- names(type_res)[[1]]
          type_res <- type_res[[1]]
        }

        curr_assignments[Matrix::which(type_res == TRUE)] <- cell_type

        level_table[[curr_level]][Matrix::which(type_res == TRUE)] <- cell_type
        level_table <- fill_in_assignments(curr_assignments, classifier, child,
                                           imputed_gate_res, level_table)
      }
    }

    return (level_table)
  }

  level_table <- fill_in_assignments(assignments, classifier, v = 1,
                                     imputed_gate_res, level_table)

  cell_type <- level_table$level1
  names(cell_type) <- level_table$cell
  for(col in names(level_table)[2:ncol(level_table)]) {
    cell_type[level_table[[col]] != "Unknown"] <-
      level_table[[col]][level_table[[col]] != "Unknown"]
  }

  cell_type <- as.data.frame(cell_type)

  if (cluster_extend) {
    level_table$cluster <- pData(cds)$louv_cluster
    community_assign <- data.frame(cluster = unique(pData(cds)$louv_cluster),
                                   assign = "Unknown",
                                   stringsAsFactors = FALSE)
    for(col in names(level_table)[2:(ncol(level_table)-1)]) {
      for(clust in community_assign$cluster) {
        sub <- level_table[level_table$cluster == clust,]
        num_unk <- sum(sub[[col]] == "Unknown")
        freqs <- as.data.frame(table(sub[[col]]))
        if(nrow(freqs) == 1) next
        freqs <- freqs[freqs$Var1 != "Unknown",]
        putative_type <- as.character(freqs$Var1[which.max(freqs$Freq)])
        if (freqs$Freq[freqs$Var1 == putative_type]/sum(freqs$Freq) > .9 &
            num_unk/(sum(freqs$Freq) + num_unk) < .95 &
            freqs$Freq[freqs$Var1 == putative_type] > 5) {
          community_assign[["assign"]][
            community_assign[["cluster"]] == clust] <- putative_type
        }
      }
    }

    pData(cds)$louv_cluster <- plyr::mapvalues(x = pData(cds)$louv_cluster,
                                               from = community_assign$cluster,
                                               to = community_assign$assign)
    cell_type$cluster_ext_type <- as.character(pData(cds)$louv_cluster)
    cell_type$cell_type <- as.character(cell_type$cell_type)
    cell_type$cluster_ext_type[cell_type$cluster_ext_type == "Unknown"] <-
      cell_type$cell_type[cell_type$cluster_ext_type == "Unknown"]
  }

  return(cell_type)
}

make_predictions <- function(cds,
                             classifier,
                             curr_node,
                             rank_prob_ratio,
                             cores = 1,
                             s) {
  cvfit <- igraph::V(classifier@classification_tree)[curr_node]$model[[1]]

  predictions <- tryCatch({
    if(is.null(cvfit)) {
      child_cell_types <- igraph::V(classifier@classification_tree)[
        suppressWarnings(outnei(curr_node)) ]$name
      predictions <- matrix(FALSE, nrow=nrow(pData(cds)),
                            ncol=length(child_cell_types),
                            dimnames=list(row.names(pData(cds)),
                                          child_cell_types))
      predictions <- split(predictions, rep(1:ncol(predictions),
                                            each = nrow(predictions)))
      names(predictions) <- child_cell_types
      predictions
    } else {
      candidate_model_genes <- cvfit$glmnet.fit$beta[[1]]@Dimnames[[1]]
      good_genes <- intersect(row.names(exprs(cds)),
                              candidate_model_genes)
      if (length(good_genes) == 0) stop(paste("None of the model genes are in",
                                              "your CDS object. Did you",
                                              "specify the correct",
                                              "cds_gene_id_type and the",
                                              "correct db?"))
      x <- Matrix::t(exprs(cds[intersect(row.names(exprs(cds)),
                                         candidate_model_genes),])) #slow

      extra <- as(matrix(0, nrow = nrow(x),
                         ncol = length(setdiff(candidate_model_genes,
                                               colnames(x)))), "sparseMatrix")
      row.names(extra) <- row.names(x)
      colnames(extra) <- setdiff(candidate_model_genes, colnames(x))

      x <- cbind(x, extra)
      x <- x[,candidate_model_genes]

      # predict probabilities using fitted model
      nonz <- Matrix::rowSums(do.call(cbind,
                                      glmnet::coef.cv.glmnet(cvfit,
                                                             s="lambda.min")))
      nonz <- names(nonz[nonz != 0])
      nonz <- nonz[2:length(nonz)]
      if (sum(!nonz %in% row.names(exprs(cds))) > 0) {
        warning(paste("The following genes used in the classifier are not",
                      "present in the input CDS. Interpret with caution.",
                      nonz[!nonz %in% row.names(exprs(cds))]))
      }

      temp <- stats::predict(cvfit, #slow
                             newx = x,
                             s = s,
                             type = "response")
      prediction_probs <- as.matrix(as.data.frame(temp))

      # normalize probabilities by dividing by max
      prediction_probs <- prediction_probs/Biobase::rowMax(prediction_probs)

      # find the odds ratio of top prob over second best
      prediction_probs <- apply(prediction_probs, 1, function(x) {
        m <- names(which.max(x))
        s <- sort(x, decreasing = T)
        c(cell_type = m, odds_ratio = s[1]/s[2])
      })

      prediction_probs <- as.data.frame(t(prediction_probs))
      prediction_probs$cell_name <- row.names(prediction_probs)
      names(prediction_probs) <- c("cell_type", "odds_ratio", "cell_name")
      prediction_probs$odds_ratio <-
        as.numeric(as.character(prediction_probs$odds_ratio))

      # odds ratio has to be larger than rank_prob_ratio
      assignments <- prediction_probs[prediction_probs$odds_ratio >
                                        rank_prob_ratio,]

      # odds ratio also must be larger than expected by random guess
      # (1/number of cell types)
      random_guess_thresh <- 1.0 / length(cvfit$glmnet.fit$beta)
      assignments <- assignments[assignments$odds_ratio > random_guess_thresh,]

      not_assigned <- row.names(pData(cds))[ !row.names(pData(cds)) %in%
                                               assignments$cell_name]
      if(length(not_assigned) > 0) {
        assignments <- rbind(assignments,
                             data.frame(cell_name = not_assigned,
                                        cell_type = NA, odds_ratio = NA))
      }

      assignments$cell_type <- stringr::str_replace_all(assignments$cell_type,
                                                        "\\.1",
                                                        "")

      # reformat predictions
      predictions <- reshape2::dcast(assignments, cell_name ~ cell_type,
                                     value.var = "odds_ratio")
      row.names(predictions) <- predictions$cell_name

      if (ncol(predictions) > 2){
        predictions <- predictions[,setdiff(colnames(predictions), "NA")]
        predictions <- predictions[,-1, drop=FALSE]
        predictions <- predictions[rownames(pData(cds)),,drop=FALSE]
        predictions <- as.matrix(predictions)
        predictions[is.na(predictions)] <- FALSE
        predictions[predictions != 0] <- TRUE
        cell_type_names <- colnames(predictions)

        predictions <- split(predictions, rep(1:ncol(predictions),
                                              each = nrow(predictions)))
        names(predictions) <- cell_type_names

      } else {
        cell_type_names <- names(cvfit$glmnet.fit$beta)
        one_type <- names(predictions)[2]
        if (one_type == "NA") {
          names(predictions)[2] <- "Unknown"
          one_type <- "Unknown"
        }
        predictions <- matrix(FALSE, nrow=nrow(pData(cds)),
                              ncol=length(cell_type_names),
                              dimnames=list(row.names(pData(cds)),
                                            cell_type_names))
        predictions[,one_type] <- TRUE

        predictions <- split(predictions, rep(1:ncol(predictions),
                                              each = nrow(predictions)))
        names(predictions) <- cell_type_names
      }
      predictions
    }

  },
  #warning = function(w) print(w),
  error = function(e) {
    if (e$message == paste("None of the model genes are in your CDS object.",
                           "Did you specify the correct cds_gene_id_type and",
                           "the correct db?"))
      stop(e)
    print (e)
    cell_type_names <- names(cvfit$glmnet.fit$beta)
    predictions <- matrix(FALSE, nrow=nrow(pData(cds)),
                          ncol=length(cell_type_names),
                          dimnames=list(row.names(pData(cds)),
                                        cell_type_names))
    predictions <- split(predictions, rep(1:ncol(predictions),
                                          each = nrow(predictions)))
    names(predictions) <- cell_type_names
    predictions
  })

  for (i in 1:length(predictions)){
    p <- as(as(predictions[[i]], "sparseVector"), "sparseMatrix")
    row.names(p) <- row.names(pData(cds))
    predictions[[i]] <- p
  }

  return(predictions)
}


get_communities <- function(cds) {
  k <- 20

  fm_rowsums <- Matrix::rowSums(exprs(cds))
  FM <- exprs(cds)[is.finite(fm_rowsums) & fm_rowsums != 0, ]

  x <- Matrix::t(FM)
  n <- min(50, min(dim(FM)) - 1)
  args <- list(A = x, nv = n)

  x <- DelayedArray::DelayedArray(x)
  args$center <- round(DelayedMatrixStats::colMeans2(x), 10)
  args$scale <- sqrt(DelayedMatrixStats::colVars(x))

  s <- do.call(irlba::irlba, args = args)

  sdev <- s$d/sqrt(max(1, nrow(x) - 1))
  pcs <- sweep(s$u, 2, s$d, FUN = `*`)
  colnames(pcs) <- paste("PC", seq(1, n), sep = "")

  vars <- sdev^2
  imp <- rbind(`Standard deviation` = sdev,
               `Proportion of Variance` = round(vars, 5),
               `Cumulative Proportion` = round(cumsum(vars), 5))

  row.names(pcs) <- colnames(FM)
  cell_names <- colnames(FM)

  tmp <- RANN::nn2(pcs, pcs, k + 1, searchtype = "standard")
  neighborMatrix <- tmp[[1]][, -1]

  links <- monocle:::jaccard_coeff(neighborMatrix, FALSE)

  links <- links[links[, 1] > 0, ]
  relations <- as.data.frame(links)
  colnames(relations) <- c("from", "to", "weight")

  relations$from <- cell_names[relations$from]
  relations$to <- cell_names[relations$to]
  g <- igraph::graph.data.frame(relations, directed = FALSE)

  Q <- igraph::cluster_louvain(g)

  pData(cds)$louv_cluster <- factor(igraph::membership(Q))
  cds
}

collect_gene_names <- function(cellont_list) {
  genes <- lapply(cellont_list, function(x) {
    if(length(collect_genes(x)@gene_names) != 0) {
      pt <- ifelse(identical(x@parenttype, character(0)), "root", x@parenttype)
      data.frame(genes = collect_genes(x)@gene_names, parent = pt,
                 cell_type = x@name)
    }
  })
  all <- do.call("rbind",genes)
  return(all)
}

