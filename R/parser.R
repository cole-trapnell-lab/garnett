Lexer <- R6::R6Class(
  "Lexer",
  public = list(
    tokens = c("KEYWORD", "NAME", "NEWLINE", "NUM",
               "SIMP_KEY", "COM"),
    literals = c(">", ":", ","),
    t_COM = function(re='#.*',t) {}, # Comment character
    t_KEYWORD =
      c('expressed below|expressed above|expressed between|subtype of'),
    t_SIMP_KEY = c('celltype|expressed|not expressed|references'),
    t_NAME =
      '[a-zA-Z0-9_+/\\-\\.|=`~\\*&<^%?@!$;:]*[a-zA-Z][a-zA-Z0-9_+/\\-\\.|=`~\\*&<^%?@!$;]*',
    t_NUM = '([0-9]*\\.[0-9]+)|([0-9]+)',
    t_ignore = " \t",
    t_NEWLINE = "\n|\r\n",
    t_error = function(t) {
      cat(sprintf(" Marker file error. Illegal character '%s'\n", t$value[1]))
      t$lexer$skip(1)
      return(t)
    }
  )
)

#' Parsing the Garnett marker file
#'
#' Garnett uses a marker file to allow users to specify cell type definitions.
#' While the marker file is designed to be easy to construct and human-readable,
#' it is parsed by Garnett automatically, and so it needs to follow certain
#' formatting constraints.
#'
#' The following describes the constraints necessary in the input to the
#' \code{marker_file} argument of \code{\link{train_cell_classifier}} and
#' \code{\link{check_markers}}.
#'
#' @section Elements of a cell type description:
#' The basic structure of the Garnett marker file is a series of entries, each
#' describing elements of a cell type. After the cell name, each additional
#' line will be a descriptor, which begins with a keyword, followed by a colon
#' (':'). After the colon, a series of specifications can be added, separated
#' by commas (','). Descriptors may spill onto following lines so long as you
#' do not split a specification across multiple lines (i.e. if breaking up a
#' long descriptor across multiple lines, all but the last line should end with
#' a comma). Each new descriptor should begin on a new line. A generic cell
#' type entry looks like this:
#'
#' ```
#' > cell type name
#' descriptor: spec1, spec2,
#' spec3, spec4
#' descriptor2: spec1
#' ```
#'
#' The following are the potential descriptors:
#' \describe{
#'   \item{cell name}{\strong{Required} Each cell type must have a unique name,
#'   and the name should head the cell type description. To indicate a new cell
#'   type, use the \code{>} symbol, followed by the cell name, followed by a
#'   new line. For example, \code{> T cell}.}
#'   \item{expressed:}{\strong{Required} After the cell name, the minimal
#'   requirement for each cell type is the name of a single marker gene. The
#'   line in the marker file will begin with \code{expressed:}, followed by one
#'   or more gene names separated by commas. The last gene name of the
#'   descriptor is not followed by a comma. Gene IDs can be of any type
#'   (ENSEMBL, SYMBOL, etc.) that is present in the Bioconductor
#'   \code{\link[AnnotationDbi]{AnnotationDb-class}} package for your species.
#'   (See available packages on the
#'   \href{http://bioconductor.org/packages/3.8/data/annotation/}{Bioconductor website}).
#'   For example, for human, use \code{\link[org.Hs.eg.db]{org.Hs.eg.db}}. To
#'   see available gene ID types, you can run \code{columns(db)}. You will
#'   specify which gene ID type you used when calling
#'   \code{\link{train_cell_classifier}}.} If your species does not have an
#'   annotation dataset of type \code{\link[AnnotationDbi]{AnnotationDb-class}},
#'   you can set \code{db = 'none'}, however Garnett will then not convert gene
#'   ID types, so CDS and marker file gene ID types need to be the same.
#'   \item{not expressed:}{In addition to
#'   specifying genes that the cell type should express, you can also specify
#'   genes that your cell type should not express. Details on specifying genes
#'   are the same as for \code{expressed:}.}
#'   \item{subtype of:}{When present, this descriptor specifies that a cell
#'   type is a subtype of another cell type that is also described in the
#'   marker file. A biological example would be a CD4 T cell being a subtype of
#'   a T cell. This descriptor causes the cell type to be classified on a
#'   separate sub-level of the classification hierarchy, after the
#'   classification of its parent type is done (i.e. first T cells are
#'   discriminated from other cell types, then the T cells are subclassified
#'   using any cell types with the descriptor \code{subtype of: T cell}).
#'   \code{subtype of:} can only include a single specification, and the
#'   specification must be the exact name of another cell type specified in
#'   this marker file.}
#'   \item{references:}{This descriptor is not required, but is highly
#'   recommended. The specifications for this descriptor should be links/DOIs
#'   documenting how you chose your marker genes. While these specifications
#'   will not influence cell type classification, they will be packaged with
#'   the built classifier so that future users of the classifier can trace the
#'   origins of the markers/ }
#'   \item{*meta data:}{This wildcard descriptor allows you to specify any
#'   other property of a cell type that you wish to specify. The keyword will
#'   be the name of the column in your \code{colData} (meta data) table that you
#'   wish to specify, and the specifications will be a list of acceptable
#'   values for that meta data. An example use of this would be
#'   \code{tissue: liver, kidney}, which would specify that training cells for
#'   this cell type must have "liver" or "kidney" as their entry in the
#'   "tissue" column of the \code{colData} table.}
#'   \item{expressed below:}{While we recommend that you use \code{expressed:}
#'   and \code{not expressed:} to specify the cell type's marker genes, because
#'   these terms utilize the entirety of Garnett's built-in normalization and
#'   standardization, you can also specify expression using the following
#'   logical descriptors
#'   \code{expressed below:, expressed above:, expressed between:}.
#'   Note that no normalization occurs with these descriptors; they are used as
#'   logical gates only. To specify \code{expressed below:}, use the gene name,
#'   followed by a space, followed by a number. This will only allow training
#'   cells that have this gene expressed below the given value \strong{in the
#'   units of the expression matrix provided}. For example,
#'   \code{expressed below: MYOD1 7, MYH3 2}.}
#'   \item{expressed above:}{Similar to \code{expressed below:}, but will only
#'   allow training cells expressing the given gene above the value provided.}
#'   \item{expressed between:}{Similar to \code{expressed below:}, but provide
#'   two values separated by spaces. For example
#'   \code{expressed between: ACT5 2 5.5, ACT2 1 2.7}. This descriptor will
#'   only allow training cells expressing the given gene between the two values
#'   provided.}
#' }
#' @section Checking your marker file:
#'
#' Because only specific expressed markers are useful for Garnett
#' classification, we recommend that you always check your marker file for
#' ambiguity before proceeding with classification. To do this, we have
#' provided the functions \code{\link{check_markers}} and
#' \code{\link{plot_markers}} to facilitate marker checking. See that manual
#' pages for those functions for details.
#'
#' @seealso \code{\link{train_cell_classifier}}
#'
Parser <- R6::R6Class(
  "Parser",
  public = list(
    tokens = c("KEYWORD", "NAME", "NEWLINE",
               "SIMP_KEY", "NUM"),
    literals = c(">", ":", ","),

    # Parsing rules
    linenum = 1,
    p_file_1 = function(doc="file : NEWLINE file", p) {
      p$set(1, p$get(3))
    },
    p_file_2 = function(doc="file : celltype
                                | file celltype", p) {
      if(p$length() == 2) {
        cts = new.env(hash=TRUE)
        name_order <- list()
        cts[["name_order"]] <- name_order
        cell_ont <- p$get(2)
        cts[[cell_ont@name]] <- cell_ont
        cts[["name_order"]] <- c(cts[["name_order"]], cell_ont@name)
      } else {
        cts <- p$get(2)
        cell_ont <- p$get(3)
        if(cell_ont@name %in% names(cts)) {
          stop(paste0("Cell type '", cell_ont@name,
                      "' is defined a second time at or near line ",
                      self[["linenum"]], "."))
        }
        cts[[cell_ont@name]] <- cell_ont
        cts[["name_order"]] <- c(cts[["name_order"]], cell_ont@name)
      }
      p$set(1, cts)
    },
    p_header_start = function(doc="header_start : '>' NAME
                                                | '>' NUM
                                                | header_start NUM
                                                | header_start NAME", p) {
      if (is.character(p$get(2))) {
        ct <- new("cell_rules")
        ct@name <- p$get(3)
      } else {
        ct <- p$get(2)
        ct@name <- paste(ct@name, p$get(3))
      }
      p$set(1, ct)
    },

    p_header_1 = function(doc="header : header_start NEWLINE", p) {
      self[["linenum"]] <- self[["linenum"]] + 1
      p$set(1, p$get(2))
    },
    p_celltype_1 = function(doc="celltype : header", p) {
      p$set(1, p$get(2))
    },
    p_celltype_2 = function(doc="celltype : celltype rule", p) {
      ct <- p$get(2)
      rule <- p$get(3)

      if (rule[[1]] == "exp") {
        ct@expressed <- c(ct@expressed, rule[2:length(rule)])
      } else if (rule[[1]] == "nexp") {
        ct@not_expressed <- c(ct@not_expressed, rule[2:length(rule)])
      } else if (rule[[1]] == "ref") {
        ct@references <- c(ct@references, rule[2:length(rule)])
      } else if (rule[[1]]  == "sub") {
        ct@parenttype <- rule[[2]]
      } else if (rule[[1]] == "meta") {
        if (nrow(ct@meta) == 0) {
          ct@meta <- rule[[2]]
        } else {
          ct@meta <- rbind(ct@meta, rule[[2]])
        }
      } else if (rule[[1]] %in% c("exp_ab", "exp_bel", "exp_bet")) {
        ct@gene_rules <- c(ct@gene_rules, rule[2:length(rule)])
      }
      p$set(1, ct)
    },
    p_rule_head_1 = function(doc="rule_head : SIMP_KEY ':'", p) {
      p$set(1, p$get(2))
    },
    p_rule_head_comp_1 = function(doc="rule_head_comp : KEYWORD ':'", p) {
      p$set(1, p$get(2))
    },
    p_rule_head_2 = function(doc="rule_head : NAME ':'", p) {
      p$set(1, p$get(2))
    },
    p_rule_1 = function(doc="rule : rule_head expression NEWLINE", p) {
      self[["linenum"]] <- self[["linenum"]] + 1
      if(p$get(2) == "expressed") {
        p$set(1, c("exp", p$get(3)))
      } else if(p$get(2) == "not expressed") {
        p$set(1, c("nexp", p$get(3)))
      } else if(p$get(2) == "references") {
        p$set(1, c("ref", p$get(3)))
      } else {
        p$set(1, list("meta", data.frame(name = p$get(2),
                                         spec = p$get(3),
                                         stringsAsFactors = FALSE)))
      }
    },
    p_rule_2 = function(doc="rule : rule_head_comp comp_expression NEWLINE",
                        p) {
      self[["linenum"]] <- self[["linenum"]] + 1
      if(!is.list(p$get(3))) p$set(3, list(p$get(3)))
      if(p$get(2) == "expressed above") {
        grs <- lapply(p$get(3), function(x) {
          if (length(x) != 2) stop(paste0("Syntax error in marker file at or",
                                          " near line ", self[["linenum"]],
                                          ": expressed above needs one ",
                                          "value"))
          x <- new("gene_rule", gene_name = x[1],
                   lower = as.numeric(x[2]), upper = Inf)
        })
        p$set(1, c("exp_ab", grs))
      } else if(p$get(2) == "expressed below") {
        grs <- lapply(p$get(3), function(x) {
          if (length(x) != 2) stop(paste0("Syntax error in marker file at or",
                                          " near line ", self[["linenum"]],
                                          ": expressed below needs one ",
                                          "value"))
          x <- new("gene_rule", gene_name = x[1], upper = as.numeric(x[2]),
                   lower = -1)
        })
        p$set(1, c("exp_bel", grs))
      } else if(p$get(2) == "expressed between") {
        grs <- lapply(p$get(3), function(x) {
          if (length(x) != 3) stop(paste0("Syntax error in marker file at or",
                                          " near line ", self[["linenum"]],
                                          ": expressed between needs two",
                                          " values"))
          n1 <- as.numeric(x[2])
          n2 <- as.numeric(x[3])
          if (n1 > n2) stop(paste0("expressed between: ", x[1],
                                   " has expression values out of order ",
                                   "- swapping."))
          x <- new("gene_rule", gene_name = x[1], lower = min(n1, n2),
                   upper = max(n1, n2))
        })
        p$set(1, c("exp_bet", grs))
      } else if(p$get(2) == "subtype of") {
        p$set(1, list("sub", paste(unlist(p$get(3)), collapse = " ")))
      }
    },
    p_rule_3 = function(doc="rule : rule NEWLINE", p) {
      p$set(1, p$get(2))
      self[["linenum"]] <- self[["linenum"]] + 1
    },
    p_expression_1 = function(doc="expression : NAME
                                              | NUM", p) {
      p$set(1, p$get(2))
    },
    p_expression_2 = function(doc="expression : expression NAME
                                              | expression NUM", p) {
      if(length(p$get(2)) == 1) {
        p$set(1, paste(p$get(2), p$get(3)))
      } else {
        p$set(1, c(p$get(2)[1:(length(p$get(2)) - 1)],
                   paste(p$get(2)[length(p$get(2))], p$get(3))))
      }

    },
    p_expression_3 = function(doc="expression : expression ',' NAME
                                              | expression ',' NUM", p) {
      p$set(1, c(p$get(2), p$get(4)))
    },
    p_expression_4 = function(doc="expression : expression ',' NEWLINE NAME
                                              | expression ',' NEWLINE NUM",
                              p) {
      p$set(1, c(p$get(2), p$get(5)))
    },
    p_comp_expression_1 = function(doc="comp_expression : NAME
                                                        | NUM" , p) {
      p$set(1, p$get(2))
    },
    p_comp_expression_2 = function(doc="comp_expression : comp_expression NUM
                                                        | comp_expression NAME",
                                   p) {
      p$set(1, c(p$get(2), p$get(3)))
    },
    p_comp_expression_3 =
      function(doc="comp_expression : comp_expression ',' NAME NUM", p) {
        l <- p$get(2)
        if(!is.list(l)) l <- list(l)
        l[[length(l) + 1]] <- c(p$get(4), p$get(5))
        p$set(1, l)
      },
    p_comp_expression_4 =
      function(doc="comp_expression : comp_expression ',' NEWLINE NAME NUM",
               p) {
        l <- p$get(2)
        if(!is.list(l)) l <- list(l)
        l[[length(l) + 1]] <- c(p$get(5), p$get(6))
        p$set(1, l)
      },
    p_comp_expression_5 =
      function(doc="comp_expression : comp_expression ',' NAME NUM NUM", p) {
        l <- p$get(2)
        if(!is.list(l)) l <- list(l)
        l[[length(l) + 1]] <- c(p$get(4), p$get(5), p$get(6))
        p$set(1, l)
      },
    p_comp_expression_6 =
      function(doc="comp_expression : comp_expression ',' NEWLINE NAME NUM NUM",
               p) {
        l <- p$get(2)
        if(!is.list(l)) l <- list(l)
        l[[length(l) + 1]] <- c(p$get(5), p$get(6), p$get(7))
        p$set(1, l)
      },
    p_error = function(p) {
      if(is.null(p)) stop("Marker file error. Syntax error at EOF")
      else           stop(sprintf("Marker file error. Syntax error '%s' at or near line number %s\n",
                                  p$value, self[["linenum"]]))
    }
  )
)





