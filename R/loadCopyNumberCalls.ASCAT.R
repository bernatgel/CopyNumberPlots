# library(ASCAT)
# 
# 
# install.packages("devtools")
# devtools::install_github(repo = "Crick-CancerGenomics/ascat", subdir = "ASCAT")
# 
# 
# library(ASCAT)
# ascat.bc = ascat.loadData("ASCAT_data/Tumor_LogR.txt","ASCAT_data/Tumor_BAF.txt","ASCAT_data/Germline_LogR.txt", "ASCAT_data/Germline_BAF.txt")
# ascat.plotRawData(ascat.bc)
# ascat.bc = ascat.aspcf(ascat.bc)
# ascat.plotSegmentedData(ascat.bc)
# ascat.output = ascat.runAscat(ascat.bc) 


loadCopyNumberCalls.ASCAT <- function(cnv.data) {
  #ASCAT output does not have a specific class. Test for its basic structure
  if(methods::is(cnv.data, "list") && c("nA", "nB", "segments") %in% names(cnv.data)) {
    #It looks like ASCAT output
    segments <- cnv.data$segments[,c("chr", "startpos", "endpos", "nMajor", "nMinor", "sample")]
    segments$cn <- segments$nMajor + segments$nMinor
    segments$loh <- ifelse(segments$nMinor==0, 1, 0)
    segments <- toGRanges(segments)
    #If there's more than one sample, split per sample to get a list of GRanges
    if(length(unique(segments$sample))>1) {
      segments <- split(x = segments, segments$sample)  
    }
    
    return(segments)
  }
}
