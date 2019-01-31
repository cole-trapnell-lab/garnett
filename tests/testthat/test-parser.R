context("test-parser.R")

test_parse <- function(input_file) {
  file_str = paste0(readChar(input_file, file.info(input_file)$size), "\n")
  garnett:::parse_input(file_str)
}

test_that("basic parsing works", {
  parse_list <- test_parse("../good_parse.txt")
  expect_equal(length(parse_list[["name_order"]]), 6)
  expect_is(parse_list[["test cell 1"]], "cell_rules")
  expect_equal(parse_list[["test cell 1"]]@name,"test cell 1")
  expect_equal(parse_list[["test cell 1"]]@expressed, "hannah")
  expect_equal(parse_list[["test cell 1"]]@not_expressed, "zach")
  expect_is(parse_list[["test cell 1"]]@gene_rules, "list")
  expect_is(parse_list[["test cell 1"]]@gene_rules[[1]], "gene_rule")
  expect_equal(parse_list[["cell 2 test"]]@parenttype, "test cell 1")
  expect_equal(parse_list[["cell 2 test"]]@gene_rules[[1]]@upper, 7)
  expect_equal(parse_list[["test cell 1"]]@references[[1]], "website1.htm")
  expect_equal(parse_list[["test cell 14"]]@references[[3]], "website2.htm")
})

test_that("error messages work", {
  expect_error(garnett:::parse_input("hannah"),
               "Syntax error 'hannah' at or near line number 1")
  expect_error(garnett:::parse_input(">hannah\nexpressed below: cole\n"),
               paste("Syntax error in marker file at or near line 3: expressed",
                     "below needs one value"))
  expect_error(garnett:::parse_input(">hannah\nexpressed above: cole 1 2\n"),
               paste("Syntax error in marker file at or near line 3: expressed",
                     "above needs one value"))
  expect_error(garnett:::parse_input(">hannah\nexpressed between: cole 6\n"),
               paste("Syntax error in marker file at or near line 3: expressed",
               "between needs two values"))
  expect_error(garnett:::parse_input(">hannah\nexpressed between: cole 6 3\n"),
               paste("expressed between: cole has expression values out of",
                     "order - swapping."))
  expect_error(test_parse("../bad_parse1.txt"),
               paste("Cell type 'test cell 1' is defined a second time at or",
                     "near line 20."))
  expect_error(garnett:::parse_input(paste(">hannah\nexpressed between: cole 3",
                                           "6\n han")),
               "Syntax error at EOF")
})

