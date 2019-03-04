# loadCopyNumberCalls.CNVkit example.
cnvkit.file <-"testFiles/example.seg"
seg<-loadCopyNumberCalls.CNVkit(cnvkit.file)
# x <- read.table(file=ex.file,header = F, sep = "\t",quote = "\'")
# colnames(x) <- c("ID","chrom", "loc.start", "loc.end", "num.mark", "seg.mean")

loadCopyNumberCalls.CNVkit <- function(cnvkit.file, 
                                       chr.col=NULL, 
                                       start.col=NULL,
                                       end.col=NULL, 
                                       segment.value.col = NULL, # we can only have one value use match.arg()
                                       genome = NULL, verbose=TRUE){
  #Load file.seg. First, probe if file exist
  if(!file.exists(cnvkit.file)){
    stop(paste0(cnvkit.file, " does not exist or you are not in the correct directory."))
    
  }else{
    cnvkit.data <- read.table(file = cnvkit.file, sep="\t", skip = 1, stringsAsFactors = FALSE)
    colnames(cnvkit.data) <-  c("ID","chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
    cnvkit.data <- cnvkit.data[,c("chrom", "loc.start", "loc.end", "ID", "num.mark", "seg.mean")]
  }
  
  segs <- loadCopyNumberCalls(cnv.data = cnvkit.data, 
                                          chr.col=chr.col, 
                                          start.col=start.col,
                                          end.col=end.col, 
                                          cn.col=cn.col, 
                                          segment.value.col = segment.value.col, 
                                          genome = genome,
                                          verbose = verbose)
  return(seg)
}

x<-loadCopyNumberCalls.CNVkit(cnvkit.file)

