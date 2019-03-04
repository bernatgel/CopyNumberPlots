#loadCopyNumberCalls.DECoN

loadCopyNumberCalls.DECoN <- function(decon.file, 
                                      chr.col=NULL, 
                                      start.col= "Start",
                                      end.col= "End",
                                      cn.col="CNV.type",
                                      segment.value.col=NULL,
                                      genome=NULL, 
                                      verbose=TRUE){
  #load data
  #Load decon.file First, probe if file exist
  if(!file.exists(decon.file)){
    stop(paste0(decon.file, " does not exist or you are not in the correct directory."))
  }
    
  decon.data <- read.table(file = decon.file, sep = "\t",header = TRUE, stringsAsFactors = FALSE)
  #change the name of Start.b and End.b by "num.first.exon" and "num.last.exon"
  colnames(decon.data)[c(5,6)] <- c("Num.first.exon", "Num.last.exon")
  
  
  
  segs <- loadCopyNumberCalls(cnv.data = decon.data,
                              chr.col=chr.col, 
                              start.col=start.col,
                              end.col=end.col, 
                              cn.col=cn.col, 
                              segment.value.col = segment.value.col, 
                              genome = genome,
                              verbose = verbose)
  
  segs <- sort(segs)
  
  seqlevelsStyle(segs)<-"UCSC"
  
  return(segs)
}

#Example:

decon.file <- "testFiles/calls_all.txt"
decon.data <- read.table(file = decon.file, sep = "\t",header = TRUE, stringsAsFactors = FALSE)
decon.data <- head(decon.data, 10)
class(decon.data)

x <- loadCopyNumberCalls.DECoN(decon.file = decon.file,genome = "hg19")
  