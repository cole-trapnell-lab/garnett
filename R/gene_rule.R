#' gene_rule class
#'
#' Class definition for gene_rule, which is a container to hold a gene name
#' along with upper and lower expression bounds.
#'
#' @slot gene_name character. The name of the gene
#' @slot lower numeric. Lower bound of expression - same units as CDS object.
#' @slot upper numeric. Upper bound of expression - same units as CDS object.
#'
#' @exportClass gene_rule
setClass("gene_rule", representation(gene_name = "character",
                                     lower = "numeric",
                                     upper = "numeric"))
