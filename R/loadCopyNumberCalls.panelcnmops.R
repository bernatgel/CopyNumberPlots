# loadCopyNumberCalls.panelcnmops 

# EXAMPLE: PANELCNMOPS
data(panelcn.mops, package = "panelcn.mops")
XandCB <- test
control
resultlist
elementMetadata(XandCB) <- cbind(elementMetadata(XandCB),
                                 elementMetadata(control))
result <- panelcn.mops::panelcn.mops(XandCB)
s <-integerCopyNumber(resultlist[[1]])

# result.count <- panelcn.mops::runPanelcnMops(XandCB,countWindows = countWindows)
# 
# sampleNames <- colnames(elementMetadata(test))
sampleNames <- colnames(elementMetadata(test))

resulttable <- createResultTable(resultlist = resultlist,
                                 XandCB = XandCB,
                                 countWindows = countWindows,
                                 sampleNames = sampleNames)
table(resulttable[[1]]$CN)
panelcn.mops.res <- resulttable[[1]]

loadCopyNumberCalls.panelcnmops <- function(panelcn.mops.res, 
                                            chr.col = NULL, 
                                            start.col = NULL,
                                            end.col = NULL, 
                                            cn.col = NA,
                                            segment.value.col = NULL, 
                                            genome = NULL, 
                                            verbose=TRUE){

 
  segs <- loadCopyNumberCalls(cnv.data = panelcn.mops.res, 
                              chr.col= chr.col,  
                              start.col = start.col,  
                              end.col = end.col,  
                              cn.col = cn.col,  
                              segment.value.col = segment.value.col, 
                              genome = genome, 
                              verbose = verbose) 
  
  GenomicRanges::mcols(segs)$cn <- as.integer(gsub(x = GenomicRanges::mcols(segs)$CN, pattern = "CN", replacement = ""))
  
  return(segs)
                                            
}