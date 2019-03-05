# loadCopyNumberCalls.CNVkit example.
cnvkit.file <-"testFiles/cl_seq.cns"
loadCopyNumberCalls.CNVkit(cnvkit.file)
# x <- read.table(file=ex.file,header = F, sep = "\t",quote = "\'")
# colnames(x) <- c("ID","chrom", "loc.start", "loc.end", "num.mark", "seg.mean")

loadCopyNumberCalls.CNVkit <- function(cnvkit.file, 
                                       chr.col = NULL, 
                                       start.col = NULL,
                                       end.col = NULL, 
                                       segment.value.col = "log2", 
                                       cn.col = NA,
                                       genome = NULL,
                                       verbose = TRUE){
  #Load file.seg. First, probe if file exist
  if(!file.exists(cnvkit.file)){
    stop(paste0(cnvkit.file, " does not exist or you are not in the correct directory."))
  }
  
  segs <- loadCopyNumberCalls(cnv.data = cnvkit.file, 
                              chr.col = chr.col, 
                              start.col = start.col,
                              end.col = end.col, 
                              cn.col = cn.col, 
                              segment.value.col = segment.value.col,
                              loh.col = NA, 
                              genome = genome,
                              verbose = verbose)
  
  
  return(segs)
}

# by this function its not posible to recognize hg19 as reference genome
x<-loadCopyNumberCalls.CNVkit(cnvkit.file, genome = "hg19")

