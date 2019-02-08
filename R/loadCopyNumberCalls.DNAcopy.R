# BiocManager::install("DNAcopy")
# 
# library(DNAcopy)
# 
# 
# data(coriell)
# 
# 
# CNA.object <- CNA(cbind(coriell$Coriell.05296),
#                   coriell$Chromosome,coriell$Position,
#                   data.type="logratio",sampleid="c05296")
# 
# 
# smoothed.CNA.object <- smooth.CNA(CNA.object)
# 
# segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
# 
# plot(segment.smoothed.CNA.object, plot.type="w")
# 
# names(segment.smoothed.CNA.object)
# 
# cnv.data <- segment.smoothed.CNA.object

loadCopyNumberCalls.DNAcopy <- function(cnv.data) {
  if(methods::is(cnv.data, "DNAcopy")) {
    segments <- cnv.data$output[,c("chrom", "loc.start", "loc.end", "seg.mean", "ID", "num.mark")]
    names(segments)[4] <- "segment.value"
    segments <- toGRanges(segments)
    return(segments)
  }
}

