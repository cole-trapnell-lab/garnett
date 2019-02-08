setOldClass(c("igraph"), prototype=structure(list(), class="igraph"))

#' The garnett_classifier class
#'
#' Classifies cells according to a hierarchy of types.
#'
#' Classifies cells according to a hierarchy of types via user-defined gating
#' functions.
#'
#'
#'@section Slots:
#'  \describe{
#'    \item{\code{classification_tree}:}{Object of class \code{"igraph"}}
#'    \item{\code{cell_totals}:}{Object of class \code{"numeric"}}
#'    \item{\code{gene_id_type}:}{Object of class \code{"character"}}
#'    \item{\code{references}:}{Object of class \code{"list"}}
#'  }
#'
#' @name garnett_classifier
#' @rdname garnett_classifier
#' @aliases garnett_classifier-class
#' @exportClass garnett_classifier
setClass( "garnett_classifier",
          slots = c(classification_tree="igraph",
                    cell_totals = "numeric",
                    gene_id_type = "character",
                    references = "list"))
