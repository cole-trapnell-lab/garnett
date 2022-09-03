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
#'  the colData table of the input CDS has a column called "garnett_cluster",
#'  this will be used for cluster-extended assignments. Otherwise, assignments
#'  are calculated using Louvain community detection in PCA space. This
#'  assignment is returned as a column in the output CDS colData table. For large
#'  datasets, if the "garnett_cluster" column is not provided and
#'  \code{cluster_extend = TRUE}, the function can be significantly slower the
#'  first time it is run. See details for more information.
#' @param verbose Logical. Should progress messages be printed.
#' @param cluster_extend_max_frac_unknown Numeric between 0 and 1. The maximum
#'   fraction of a cluster allowed to be classified as 'Unknown' and still
#'   extend classifications to the cluster. Only used when
#'   \code{cluster_extend = TRUE}. Default is 0.95. See details.
#' @param cluster_extend_max_frac_incorrect Numeric between 0 and 1. The
#'   maximum fraction of classified cells in a cluster allowed to be
#'   incorrectly classified (i.e. assigned to a non-dominant type) and still
#'   extend classifications to the cluster. Fraction does not include 'Unknown'
#'   cells. Only used when \code{cluster_extend = TRUE}. Default is 0.1. See
#'   details.
#' @param return_type_levels Logical. When \code{TRUE}, the function additionally
#'   appends assignments from each hierarchical level in the classifier as columns
#'   in the pData table labeled \code{cell_type_li}, where "i" indicates the
#'   corresponding level index
#'
#' @details This function applies a previously trained multinomial glmnet
#'  classifier at each node of a previously defined garnett_classifier tree.
#'  The output is a CDS object with cell type classifications added to the
#'  colData table.
#'
#'  When \code{cluster_extend = TRUE}, louvain communities are calculated in
#'  PCA space. Any cluster where >\code{cluster_extend_max_frac_unknown},
#'  (default 90%) of classified cells are of a single type,
#'  >\code{1 - cluster_extend_max_frac_unknown} (default 5%) of cells are classified, and a minimum of 5 cells are classified will
#'  be assigned that cluster-extended type. Both cluster-extended type and
#'  originally calculated cell type are reported.
#'
#' @return CDS object with classifications in the \code{colData} table.
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
                           verbose = FALSE,
                           cluster_extend_max_frac_unknown = 0.95,
                           cluster_extend_max_frac_incorrect = 0.1,
                           return_type_levels = FALSE) {
  if(verbose) message("Starting classification")
  ##### Check inputs #####
  if(verbose) message("Checking inputs")
  assertthat::assert_that(is(cds, "cell_data_set"))
  assertthat::assert_that(assertthat::has_name(colData(cds), "Size_Factor"),
                          msg = paste("Must run estimate_size_factors() on cds",
                                      "before calling classify_cells"))
  assertthat::assert_that(sum(is.na(colData(cds)$Size_Factor)) == 0,
                          msg = paste("Must run estimate_size_factors() on cds",
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
    classifier_gene_id_type <- 'temp'
  }

  assertthat::assert_that(is.numeric(rank_prob_ratio))
  assertthat::assert_that(rank_prob_ratio > 1,
                          msg = "rank_prob_ratio must be greater than 1")
  assertthat::assert_that(is.logical(cluster_extend))
  assertthat::assert_that(is.logical(verbose))
  assertthat::assert_that(is.logical(return_type_levels))

  ##### Set internal parameters #####
  s <- "lambda.min"

  ##### Normalize CDS #####
  orig_cds <- cds

  if(verbose) message("Normalizing CDS object\n")

  if (!is(counts(cds), "dgCMatrix")) {
    sf <- colData(cds)$Size_Factor
    cds <- suppressWarnings(new_cell_data_set(as(counts(cds), "dgCMatrix"),
                                       cell_metadata = colData(cds),
                                       gene_metadata = rowData(cds)))
    colData(cds)$Size_Factor <- sf
  }

  colData(cds)$num_genes_expressed <- Matrix::colSums(as(counts(cds),
                                                       "lgCMatrix"))
  save_sf <- colData(cds)$Size_Factor
  new_cell_totals <- Matrix::colSums(counts(cds))

  excluded_cells <- NULL
  if(sum(new_cell_totals == 0) != 0) {
    warning(paste0(sum(new_cell_totals == 0), " cells in cds have no reads. These cells will be excluded from classification."))
    excluded_cells <- names(new_cell_totals == 0)
    cds <- cds[,new_cell_totals != 0]
    new_cell_totals <- new_cell_totals[new_cell_totals != 0]
  }
  sfs <- new_cell_totals/(classifier@cell_totals *
                            stats::median(colData(cds)$num_genes_expressed))
  sfs[is.na(sfs)] <- 1

  colData(cds)$Size_Factor <- sfs
  temp <- counts(cds)
  temp@x <- temp@x / rep.int(colData(cds)$Size_Factor, diff(temp@p))
  norm_cds <- suppressWarnings(new_cell_data_set(temp,
                             cell_metadata = colData(cds),
                             gene_metadata = rowData(cds)))

  if (classifier_gene_id_type != "custom") {
    if (methods::.hasSlot(classifier, "gene_id_type")) {
      classifier_gene_id_type <- classifier@gene_id_type
    } else {
      classifier_gene_id_type <- "ENSEMBL"
    }
  }

  ### Convert to Classifier IDs ###
  if(cds_gene_id_type != classifier_gene_id_type) {
    if (verbose) message(paste("Converting CDS IDs to",
                               classifier_gene_id_type, "\n"))
    lstart <- nrow(rowData(norm_cds))
    norm_cds <- cds_to_other_id(norm_cds,
                                db=db,
                                cds_gene_id_type,
                                classifier_gene_id_type,
                                verbose = FALSE)
    lend <- nrow(rowData(norm_cds))
  } else if (cds_gene_id_type == "ENSEMBL") {
    norm_cds <- fix_ensembl_id(norm_cds)
  }

  colData(norm_cds)$Size_Factor <- sfs
  cds <- orig_cds

  ##### Calculate cell communities #####
  if (cluster_extend) {
    if ("garnett_cluster" %in% names(colData(cds))) {
      colData(norm_cds)$louv_cluster <- colData(cds)$garnett_cluster
    } else {
      if(verbose) message(paste("No garnett_cluster column provided,",
                                "generating clusters for classification\n"))
      norm_cds <- get_communities(norm_cds)
      colData(cds)$garnett_cluster <- NA
      colData(cds)[row.names(colData(norm_cds)),]$garnett_cluster <- colData(norm_cds)$louv_cluster
    }
  }

  ##### Classify cells #####
  if(verbose) message("Predicting cell types\n")
  class_df <- run_classifier(classifier, norm_cds,
                             cluster_extend = cluster_extend,
                             s=s,
                             rank_prob_ratio = rank_prob_ratio,
                             cluster_extend_max_frac_unknown = cluster_extend_max_frac_unknown,
                             cluster_extend_max_frac_incorrect = cluster_extend_max_frac_incorrect,
                             return_type_levels = return_type_levels)
  if(!is.null(excluded_cells)) {
    ext <- matrix(ncol=ncol(class_df), nrow = length(excluded_cells),
                  dimnames = list(excluded_cells))
    colnames(ext) <- colnames(class_df)
    class_df <- rbind(class_df, ext)
    class_df <- class_df[row.names(colData(cds)),]
  }

  colData(cds)$cell_type <- NULL
  if("cluster_ext_type" %in% names(colData(cds)))
    colData(cds)$cluster_ext_type <- NULL

  colData(cds)$Size_Factor <- save_sf
  colData(cds) <- cbind(colData(cds), class_df)
  if(verbose) message("Complete!\n")
  cds
}

run_classifier <- function(classifier,
                           cds,
                           cluster_extend,
                           rank_prob_ratio,
                           s,
                           cluster_extend_max_frac_unknown,
                           cluster_extend_max_frac_incorrect,
                           return_type_levels) {

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

        new_assignment_mask <- type_res == 1
        if (length(parents) > 1) {
          new_assignment_mask <- new_assignment_mask & (curr_assignments == parents[[length(parents) - 1]])
        }
        curr_assignments[Matrix::which(new_assignment_mask)] <- cell_type
        level_table[[curr_level]][Matrix::which(new_assignment_mask)] <- cell_type

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
    level_table$cluster <- colData(cds)$louv_cluster
    community_assign <- data.frame(cluster = unique(colData(cds)$louv_cluster),
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
        if (freqs$Freq[freqs$Var1 == putative_type]/sum(freqs$Freq) >
            (1 - cluster_extend_max_frac_incorrect) &
            num_unk/(sum(freqs$Freq) + num_unk) < cluster_extend_max_frac_unknown &
            freqs$Freq[freqs$Var1 == putative_type] > 5) {
          community_assign[["assign"]][
            community_assign[["cluster"]] == clust] <- putative_type
        }
      }
    }

    colData(cds)$louv_cluster <- plyr::mapvalues(x = colData(cds)$louv_cluster,
                                               from = community_assign$cluster,
                                               to = community_assign$assign)
    cell_type$cluster_ext_type <- as.character(colData(cds)$louv_cluster)
    cell_type$cell_type <- as.character(cell_type$cell_type)
    cell_type$cluster_ext_type[cell_type$cluster_ext_type == "Unknown"] <-
      cell_type$cell_type[cell_type$cluster_ext_type == "Unknown"]
  }

  if (return_type_levels) {
    level_table <- level_table[, grep("level", colnames(level_table))]
    for(col in 2:ncol(level_table)) {
      unknown_mask <- (level_table[[col]] == "Unknown")
      level_table[[col]][unknown_mask] <- level_table[[col - 1]][unknown_mask]
    }

    cell_type[gsub("level", "cell_type_l", colnames(level_table))] <- level_table
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
      predictions <- matrix(FALSE, nrow=nrow(colData(cds)),
                            ncol=length(child_cell_types),
                            dimnames=list(row.names(colData(cds)),
                                          child_cell_types))
      predictions <- split(predictions, rep(1:ncol(predictions),
                                            each = nrow(predictions)))
      names(predictions) <- child_cell_types
      predictions
    } else {
      candidate_model_genes <- cvfit$glmnet.fit$beta[[1]]@Dimnames[[1]]
      good_genes <- intersect(row.names(counts(cds)),
                              candidate_model_genes)
      if (length(good_genes) == 0) stop(paste("None of the model genes are in",
                                              "your CDS object. Did you",
                                              "specify the correct",
                                              "cds_gene_id_type and the",
                                              "correct db?"))
      x <- Matrix::t(counts(cds[intersect(row.names(counts(cds)),
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
                                      glmnet::coef.glmnet(cvfit,
                                                             s="lambda.min")))
      nonz <- nonz[2:length(nonz)]
      nonz <- names(nonz[nonz != 0])

      if (sum(!nonz %in% row.names(counts(cds))) > 0) {
        warning(paste("The following genes used in the classifier are not",
                      "present in the input CDS. Interpret with caution.",
                      nonz[!nonz %in% row.names(counts(cds))]))
      }

      temp <- stats::predict(cvfit, #slow
                             newx = x,
                             s = s,
                             type = "response")
      temp[is.nan(temp)] <- 0
      prediction_probs <- as.matrix(as.data.frame(temp))


      # normalize probabilities by dividing by max
      prediction_probs <- prediction_probs/Biobase::rowMax(prediction_probs)

      prediction_probs[is.nan(prediction_probs)] <- 0

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

      not_assigned <- row.names(colData(cds))[ !row.names(colData(cds)) %in%
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
      predictions <- predictions[!is.na(predictions$cell_name),]
      row.names(predictions) <- predictions$cell_name

      if (ncol(predictions) > 2){
        predictions <- predictions[,setdiff(colnames(predictions), "NA")]
        predictions <- predictions[,-1, drop=FALSE]
        predictions <- predictions[rownames(colData(cds)),,drop=FALSE]
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
        predictions <- matrix(FALSE, nrow=nrow(colData(cds)),
                              ncol=length(cell_type_names),
                              dimnames=list(row.names(colData(cds)),
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
    predictions <- matrix(FALSE, nrow=nrow(colData(cds)),
                          ncol=length(cell_type_names),
                          dimnames=list(row.names(colData(cds)),
                                        cell_type_names))
    predictions <- split(predictions, rep(1:ncol(predictions),
                                          each = nrow(predictions)))
    names(predictions) <- cell_type_names
    predictions
  })

  for (i in 1:length(predictions)){
    p <- as(as(predictions[[i]], "sparseVector"), "sparseMatrix")
    row.names(p) <- row.names(colData(cds))
    predictions[[i]] <- p
  }

  return(predictions)
}


get_communities <- function(cds) {
  k <- 20

  fm_rowsums <- Matrix::rowSums(counts(cds))
  FM <- counts(cds)[is.finite(fm_rowsums) & fm_rowsums != 0, ]

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

  links <- monocle3:::jaccard_coeff(neighborMatrix, FALSE)

  links <- links[links[, 1] > 0, ]
  relations <- as.data.frame(links)
  colnames(relations) <- c("from", "to", "weight")

  relations$from <- cell_names[relations$from]
  relations$to <- cell_names[relations$to]
  g <- igraph::graph.data.frame(relations, directed = FALSE)

  Q <- igraph::cluster_louvain(g)

  colData(cds)$louv_cluster <- factor(igraph::membership(Q))
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

