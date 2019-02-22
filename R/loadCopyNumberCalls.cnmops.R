#CN.mops
#ejemplo:
# data(cn.mops, package = "cn.mops")
# resCNMOPSEX <- exomecn.mops(exomeCounts)
# resCNMOPSEX <- calcIntegerCopyNumbers(resCNMOPSEX)
# 
# resCNMOPS <- cn.mops(XRanges)
# resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)
# cn.mops.res <- resCNMOPSEX
# cn.mops.res <- "cn.mops.segmentation.csv"
# cnv.data <- "x"
# # rr <- tryCatch(toGRanges("1:1-10"), error = function(e){return(NULL)})
# # rr
# cnv.data <- "cn.mops.segmentation.csv"
# 
# 
# loadCopyNumberCalls.cnmops(cn.mops.res = cn.mops.res)  



loadCopyNumberCalls.cnmops<- function(cn.mops.res, 
                                      chr.col=NULL, 
                                      start.col=NULL,
                                      end.col=NULL, 
                                      cn.col="CN",
                                      segment.value.col = c("median", "mean"), # we can only have one value use match.arg()
                                      genome = NULL, verbose=TRUE){
  
  #if cn.mops.res is "CNVDetectionResult"
  if(methods::is(cn.mops.res, "CNVDetectionResult")){
    segs <- cn.mops::segmentation(cn.mops.res)
    segs$cn <- as.integer(gsub(segs$CN,pattern = "CN",replacement = ""))
    segs <- loadCopyNumberCalls(cnv.data = segs,
                                chr.col=chr.col, 
                                start.col=start.col, 
                                end.col=end.col, 
                                cn.col=cn.col, 
                                segment.value.col = segment.value.col,
                                genome = genome,
                                verbose = verbose)
  }
  
  #if cn.mops.res come from a file
  if(is.character(cn.mops.res)){
      segs <- loadCopyNumberCalls(cnv.data = cn.mops.res, 
                                  chr.col=chr.col, 
                                  start.col=start.col,
                                  end.col=end.col, 
                                  cn.col=cn.col, 
                                  segment.value.col = segment.value.col, 
                                  genome = genome,
                                  verbose = verbose)

      mcols(segs)$cn <- as.integer(gsub(segs$cn,pattern = "CN",replacement = ""))
      mcols(segs) <- mcols(segs)[,c("sampleName", "median","mean", "cn")]
  }
  
  #if there are more than one sample in cn.mops.res
  if(length(table(segs$sampleName))>1){
    return(split(segs, mcols(segs)$sampleName))
  }
  
  return(segs)
}  


