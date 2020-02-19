get_training_sample <- function(cds,
                                orig_cds,
                                classifier,
                                tf_idf,
                                gene_table,
                                curr_node,
                                parse_list,
                                name_order,
                                max_training_samples,
                                num_unknown,
                                back_cutoff,
                                training_cutoff,
                                marker_scores,
                                return_initial_assign,
                                meta_only) {

  ##### Find type assignment from expressed/not expressed #####

  child_cell_types <- igraph::V(classifier@classification_tree)[
    suppressWarnings(outnei(curr_node)) ]$name
  parent <- igraph::V(classifier@classification_tree)[curr_node]$name
  if (length(child_cell_types) > 0) {
    if (length(intersect(child_cell_types, colnames(marker_scores))) == 0 &
        length(intersect(child_cell_types, meta_only)) == 0) {
      return(NULL)
    }
    if (length(intersect(child_cell_types, colnames(marker_scores))) > 0) {
      assigns <- assign_type(marker_scores[,intersect(child_cell_types,
                                                      colnames(marker_scores)),
                                           drop=FALSE],
                             training_cutoff, return_initial_assign)
      if(return_initial_assign) {
        return(assigns)
      }
      assigns <- as.data.frame(assigns)
      names(assigns) <- "assigns"
      colData(orig_cds)$assigns <- assigns[row.names(colData(orig_cds)),"assigns"]
      colData(orig_cds)$assigns <- as.character(colData(orig_cds)$assigns)
      colData(orig_cds)$assigns[is.na(colData(orig_cds)$assigns)] <- "None"
    }
  }



  ##### Add/Exclude possibles using other definitions #####
  child_rules <- list()
  for (child in child_cell_types) {
    cell_class_func <-
      igraph::V(classifier@classification_tree)[child]$classify_func[[1]]

    parent <- environment(cell_class_func)
    if (is.null(parent))
      parent <- emptyenv()
    e1 <- new.env(parent=parent)

    Biobase::multiassign(names(colData(orig_cds)),
                         as.data.frame(colData(orig_cds)), envir=e1)
    environment(cell_class_func) <- e1

    type_res <- cell_class_func(counts(orig_cds))
    if (length(type_res)!= ncol(orig_cds)){
      message(paste("Error: classification function for",
                    igraph::V(classifier@classification_tree)[child]$name,
                    "returned a malformed result."))
      stop()
    }

    type_res <- as(as(type_res,"sparseVector"), "sparseMatrix")
    row.names(type_res) <- row.names(colData(orig_cds))
    colnames(type_res) <- child
    child_rules[[ child ]] <- type_res
  }


  # Assign cell types by classifier rules and check for conflicts
  ctf_cell_type <- rep("Unknown", length(child_rules[[1]]))
  names(ctf_cell_type) <- row.names(child_rules[[1]])
  for (child in child_cell_types) {
    ctf_cell_type[Matrix::which(as.matrix(child_rules[child][[1]]))] <- child
  }

  # Find training sample cells and downsample if necessary
  good_cell_type <- ctf_cell_type[ctf_cell_type %in% child_cell_types]
  obs_counts <- table(good_cell_type)
  training_sample <- c()
  for(i in names(obs_counts)){
    num_obs_for_type_i <- min(max_training_samples, obs_counts[i])
    obs_for_type_i <- sample(which(good_cell_type == i), num_obs_for_type_i)
    training_sample <- append(training_sample, obs_for_type_i)
  }
  training_sample <- good_cell_type[training_sample]

  # Find outgroup samples and downsample if necessary
  outgroup_samples <- !ctf_cell_type %in% child_cell_types

  if (length(outgroup_samples) > 0){
    # if very few cells are available for outgroup, use whole cds to find outgroup
    if(sum(outgroup_samples) < (.1 * num_unknown)) {
      out_group_cds <- cds
      outgroup_samples <- rep("Unknown", nrow(pData(cds)))
      names(outgroup_samples) <- row.names(pData(cds))
    } else {
      outgroup_samples <- ctf_cell_type[outgroup_samples]

      out_group_cds <- cds[,names(outgroup_samples)]
    }

    out_group_cds <- out_group_cds[,sample(row.names(colData(out_group_cds)),
                                           min(nrow(colData(out_group_cds)),
                                               num_unknown * 10),
                                           replace = F)]

    out_group_cds <- get_communities(out_group_cds)

    per_clust <-
      floor(num_unknown/length(unique(colData(out_group_cds)$louv_cluster)))

    outg <- lapply(unique(colData(out_group_cds)$louv_cluster), function(x) {
      sub <- colData(out_group_cds)[colData(out_group_cds)$louv_cluster == x,]
      sample(row.names(sub), min(nrow(sub), per_clust))
    })

    outg <- unlist(outg)
    if(length(outg) < min(length(outgroup_samples), num_unknown)) {
      to_get <- min(length(outgroup_samples), num_unknown)
      outgroup_samples <- outgroup_samples[!names(outgroup_samples) %in% outg]
      outg <- c(outg, sample(names(outgroup_samples), (to_get - length(outg))))
    }
    outgroup <- rep("Unknown", length(outg))
    names(outgroup) <- outg
    training_sample <- training_sample[!names(training_sample) %in% names(outgroup)]
    training_sample <- append(training_sample, outgroup)
  }

  training_sample <- factor(training_sample)
  training_sample <- droplevels(training_sample)

  return(training_sample)
}



assign_type <- function(total_vals,
                        training_cutoff,
                        return_initial_assign) {
  cutoffs <- apply(total_vals, 2, stats::quantile, probs = training_cutoff)

  q <- as.data.frame(t(t(total_vals) > cutoffs))
  q$total <- rowSums(q)
  q$assign <- sapply(q$total, function(x) {
    if (x == 0) {
      out <- "None"
    } else if (x > 1) {
      out <- "Ambiguous"
    } else {
      out <- "Assign"
    }
    out
  })
  q$cell <- row.names(total_vals)
  if (return_initial_assign) {
    q$total <- NULL
    types <- colnames(total_vals)
    q$Assignment <- "None"
    for(type in types) {
      q$Assignment[q[,type]] <- type
    }
    q$Assignment[q$assign == "Ambiguous"] <- "Ambiguous"
    q$cell <- NULL
    q$assign <- NULL
    return(q)
  }

  out <- q[q$assign != "Assign",]$assign
  names(out) <- q[q$assign != "Assign",]$cell
  sub <- q[q$assign == "Assign",]
  sub$total <- NULL
  sub$assign <- NULL
  sub <- reshape2::melt(sub, id.vars = "cell")

  sub <- sub[sub$value,]
  out2 <- as.character(sub$variable)
  names(out2) <- sub$cell
  out <- c(out, out2)
  out <- out[out != "None"]
  out
}

aggregate_positive_markers <- function(cell_type,
                                       tf_idf,
                                       gene_table,
                                       back_cutoff,
                                       agg = TRUE) {
  gene_list <- cell_type@expressed
  if (length(gene_list) == 0) return("no_gene")
  bad_genes <- gene_table[!gene_table$in_cds,]$orig_fgenes
  gene_list <- gene_list[!gene_list %in% bad_genes]
  gene_list <- gene_table$fgenes[match(gene_list,gene_table$orig_fgenes)]

  rel <- intersect(unlist(gene_list), colnames(tf_idf))
  rel <- unique(rel)
  rel_genes <- tf_idf[,rel, drop=FALSE]
  if(length(rel) == 0) return(NULL)

  if(length(unlist(gene_list)) == 1) {
    background_cutoff <- back_cutoff * stats::quantile(rel_genes[,1],
                                                       prob = .95)
    rel_genes[rel_genes < background_cutoff] <- 0
    agg <- rel_genes
  } else {
    background_cutoff <- back_cutoff * apply(rel_genes, 2,
                                             stats::quantile,
                                             prob = .95)

    temp <- Matrix::t(rel_genes)
    temp[temp < background_cutoff] <- 0
    rel_genes <- Matrix::t(temp)

    if(agg) {
      agg <-  (Matrix::rowSums(rel_genes))
    } else {
      agg <- rel_genes
    }
  }
  agg
}

get_negative_markers <- function(cell_type,
                                 tf_idf,
                                 gene_table,
                                 back_cutoff) {
  not_list <- cell_type@not_expressed
  if (length(not_list) != 0) {
    bad_genes <- gene_table[!gene_table$in_cds,]$orig_fgenes
    not_list <- not_list[!not_list %in% bad_genes]
    not_list <- gene_table$fgenes[match(not_list,gene_table$orig_fgenes)]
    rel <- intersect(unlist(not_list), colnames(tf_idf))
    rel <- unique(rel)
    rel_genes <- tf_idf[,rel, drop=FALSE]
    if(length(rel) == 0) return("")

    if(length(unlist(not_list)) == 1) {
      background_cutoff <- back_cutoff * stats::quantile(rel_genes[,1],
                                                         prob = .95)
      bad_cells <- names(rel_genes[rel_genes[,1] > background_cutoff,])
    } else {
      background_cutoff <- back_cutoff * apply(rel_genes, 2,
                                               stats::quantile,
                                               prob = .95)

      temp <- Matrix::t(rel_genes)
      temp <- temp > background_cutoff
      temp <- Matrix::t(temp)
      bad_cells <- row.names(temp[Matrix::rowSums(temp) != 0,])
    }
  } else {
    bad_cells <- ""
  }
  return(bad_cells)
}

