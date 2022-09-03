new_garnett_classifier <- function()
{
  garc <- new( "garnett_classifier",
               classification_tree = igraph::graph.empty())

  root_node_id <- "root"

  garc@classification_tree <- garc@classification_tree +
    igraph::vertex(root_node_id,
                   classify_func=list(function(x) {rep(TRUE, ncol(x))}),
                   model = NULL)

  return(garc)
}

add_cell_type <- function(classifier,
                          cell_type_name,
                          classify_func,
                          parent_cell_type_name="root")
{
  if (cell_type_name %in% igraph::V(classifier@classification_tree)$name){
    stop(paste("Error: cell type",cell_type_name, "already exists."))
  }

  classifier@classification_tree <- classifier@classification_tree +
    igraph::vertex(cell_type_name, classify_func=list(classify_func),
                   model=NULL)

  classifier@classification_tree <- classifier@classification_tree +
    igraph::edge(parent_cell_type_name, cell_type_name)
  return (classifier)
}
