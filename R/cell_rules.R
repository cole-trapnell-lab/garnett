#' cell_rules class
#'
#' Representation of cell type derived from marker file.
#'
#' @slot name character. Name of the cell type.
#' @slot gene_names character. A list of all of the genes included in the
#'  definition.
#' @slot expressed character. A list of genes defined as "expressed:".
#' @slot not_expressed character. A list of genes defined as "not expressed:".
#' @slot gene_rules vector of GeneRules-class. A list of genes defined under
#'  specific rules using "expressed below:", "expressed above:", or
#'  "expressed between:".
#' @slot meta data.frame of meta data rules specified in marker file.
#' @slot parenttype character. The name of the parent type - specified by
#'  "subtype of:".
#' @slot references character. A list of references included in the definition.
#'
#' @name cell_rules
#' @rdname cell_rules
#' @aliases cell_rules-class
#' @exportClass cell_rules
setClass("cell_rules", representation(name = "character",
                                      gene_names = "character",
                                      expressed = "character",
                                      not_expressed = "character",
                                      gene_rules = "vector",
                                      meta = "data.frame",
                                      parenttype = "character",
                                      references = "character"))
