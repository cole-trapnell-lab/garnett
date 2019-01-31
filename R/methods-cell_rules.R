setGeneric("collect_genes", signature = "x", function(x)
  standardGeneric("collect_genes"))

setMethod("collect_genes", c("x" = "cell_rules"), function(x) {
  n1 <- as.character(c(x@expressed, x@not_expressed))
  n2 <- as.character(unlist(lapply(x@gene_rules, function(y) y@gene_name)))
  x@gene_names <- unique(c(n1, n2))
  x
})


