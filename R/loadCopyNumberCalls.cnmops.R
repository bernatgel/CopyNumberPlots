#' loadCopyNumberCalls.cnmops
#'
#' @description
#' Loads copy number calls from cn.mops and panelcn.mops in a tabular format
#'
#' @details
#' This function will load segments data from any "bed-like" data structure
#' in R or file. Internally it uses the \code{toGRanges} function from regioneR
#' package and can work with any format accepted by it, including R objects
#' and local or remote files.
#' If no column names are specified, it will use simple heuristics to try to
#' identify the relevant data columns.
#'
#'
#' @note
#' The returned GRanges will have the chromosome names following the UCSC style
#' irrespective of the original format.
#'
#' @usage loadCopyNumberCalls.cnmops(cn.mops.res, chr.col=NULL, start.col=NULL,end.col=NULL, cn.col="CN", segment.value.col = c("median", "mean"),genome = NULL, verbose=TRUE)
#'
#' @param cn.mops.res The name of the file with the data or the name of the variable with the data
#' @param chr.col The name or number of the column with chromosome information. If NULL, it is automatically identified. (default to NULL)
#' @param start.col The name or number of the column with start position information. If NULL, it is automatically identified. (default to NULL)
#' @param end.col The name or number of the column with end position information. If NULL, it is automatically identified. (default to NULL)
#' @param cn.col The name or number of the column with CN information. If NULL, it is automatically identified. (default to NULL)
#' @param segment.value.col The name or number of the column with SNP identifier information. If NULL, it is automatically identified. (default to NULL)
#' @param genome The name of the genome (default to NULL)
#' @param verbose Wether information messages should be generated. (default to TRUE)
#'
#' @return
#' A GRanges object with a range per copy number segment
#'
#' @examples
#'
#'
#'
#' @export loadCopyNumberCalls.cnmops
#'
#' @importFrom cn.mops segmentation


loadCopyNumberCalls.cnmops<- function(cn.mops.res, 
                                      chr.col = NULL, 
                                      start.col = NULL,
                                      end.col = NULL, 
                                      cn.col = "CN",
                                      segment.value.col = "median", 
                                      genome = NULL, 
                                      verbose=TRUE){
  
  #if cn.mops.res is "CNVDetectionResult"
  if(methods::is(cn.mops.res, "CNVDetectionResult")){
    segs <- cn.mops::segmentation(cn.mops.res)
    segs$cn <- as.integer(gsub(segs$CN,pattern = "CN",replacement = ""))
    segs <- loadCopyNumberCalls(cnv.data = segs,
                                chr.col = chr.col, 
                                start.col = start.col, 
                                end.col = end.col, 
                                cn.col = cn.col, 
                                segment.value.col = segment.value.col,
                                loh.col = NA,
                                genome = genome,
                                verbose = verbose)

  }else{
    segs <- loadCopyNumberCalls(cnv.data = cn.mops.res, 
                                chr.col = chr.col, 
                                start.col = start.col,
                                end.col = end.col, 
                                cn.col = cn.col, 
                                segment.value.col = segment.value.col, 
                                loh.col = NA,
                                genome = genome,
                                verbose = verbose)
    
    GenomicRanges::mcols(segs)$cn <- as.integer(gsub(segs$cn,pattern = "CN",replacement = ""))
    GenomicRanges::mcols(segs) <-  GenomicRanges::mcols(segs)[,-which(names(GenomicRanges::mcols(segs)) %in% c("width","strand"))]
  }
  
  #if there is more than one sample in cn.mops.res
  if(length(table(segs$sampleName)) > 1){
     segs <- split(segs, mcols(segs)$sampleName)
  }
  
  return(segs)
}  



#CN.mops
#ejemplo:

#PANELCNMOPS
data(panelcn.mops, package = "panelcn.mops")
XandCB <- test
control
resultlist
elementMetadata(XandCB) <- cbind(elementMetadata(XandCB),
                                 elementMetadata(control))
result <- panelcn.mops::panelcn.mops(XandCB)
s <-integerCopyNumber(resultlist[[1]])

result.count <- panelcn.mops::runPanelcnMops(XandCB,countWindows = countWindows)

sampleNames <- colnames(elementMetadata(test))
resulttable <- createResultTable(resultlist = resultlist, XandCB = XandCB,
                                 countWindows = countWindows,
                                 selectedGenes = selectedGenes,
                                 sampleNames = sampleNames)
sampleNames <- colnames(elementMetadata(test))

resulttable <- createResultTable(resultlist = resultlist,
                                 XandCB = XandCB,
                                 countWindows = countWindows,
                                 sampleNames = sampleNames)
table(resulttable[[1]]$CN)
cnv.data <- resulttable[[1]]
loadCopyNumberCalls.cnmops(cn.mops.res = cnv.data)  


# df <-cnv.data
# loadCopyNumberCalls(cnv.data = cnv.data)
# resulttable<-resulttable[[1]]
# loadCopyNumberCalls.cnmops(cn.mops.res = cnv.data)
# 
# cnv.data <- resulttable[[1]]
# class(resulttable)
# class(cnv.data)
# toGRanges(cnv.data)
# colnames(x = cnv.data)
################################

#CN.MOPS
data(cn.mops, package = "cn.mops")
resCNMOPSEX <- exomecn.mops(exomeCounts)
resCNMOPSEX <- calcIntegerCopyNumbers(resCNMOPSEX)
class(resCNMOPSEX)
resCNMOPS <- cn.mops(XRanges)
resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)
cn.mops.res <- resCNMOPS
cn.mops.res <- "cn.mops.segmentation.csv"
cnv.data <- "x"
# rr <- tryCatch(toGRanges("1:1-10"), error = function(e){return(NULL)})
# rr
resCNMOPS@individualCall
loadCopyNumberCalls.cnmops(cn.mops.res = cn.mops.res)  
cn.mops.res <- "cn.mops.segmentation.csv"
cnv.data <- cn.mops.res
loadCopyNumberCalls.cnmops(cn.mops.res = cn.mops.res)  
