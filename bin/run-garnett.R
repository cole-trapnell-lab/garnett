#!/usr/bin/env Rscript
predict_immune_cell_types <- function(input.path,marker.path, output.path){
  dat <- read.table(input.path)
  temp <- load(marker.path)
  markers <- get(temp)
  rm(temp)
  
  dat_cds <- newCellDataSet(as.matrix(dat))
  dat_cds <- estimateSizeFactors(dat_cds)
  
  dat_result <- classify_cells(dat_cds,markers,
                               db = org.Hs.eg.db,
                               cluster_extend = TRUE,
                               cds_gene_id_type = "SYMBOL")
  result <- data.frame("sample"=colnames(dat),
                       "prediction"=dat_result$cell_type,
                       "score"=1)
  write.csv(result,file=output.path,row.names = FALSE)
}

getArgs<-function(){
  require(optparse)
  option_list <- list(
    make_option(c("-i", "--input"), dest='input',default=NULL, help="Path to tab-delimited input matrix"),
    make_option(c("-o", "--output"), default="output.csv", dest='output',help = "Path to output file"),
    make_option(c('-m','--marker'),default=NULL, dest='marker',help="Path to markers file")#,
#    make_option(c('-t','--testmode'),default=FALSE,dest='testmode',action='store_true',help='Run in test mode')
  )
  
  args=parse_args(OptionParser(option_list = option_list))
  
  return(args)
}

main<-function(){
  args<-getArgs()
  input.path <- args$input
  output.path <- args$output
  marker.file <- args$marker
  
#  mode=ifelse(args$testmode,'mode','prod')
  
  require(garnett)
  require(org.Hs.eg.db)
  
  predict_immune_cell_types(input.path, marker.file, output.path)
  
}

main()