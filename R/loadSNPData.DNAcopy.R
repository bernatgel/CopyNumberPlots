#' loadSNPDataDNAcopy
#'
#' @description
#' Loads SNP array data in a tabular format using data from DNAcopy rawdata or DNAcopy results
#'
#' @details
#' Given a DNAcopy rawdata or DNAcopy results object, 
#' the function load SNP array data in a tabular format.
#' It will try to identify the columns with the relevant information
#' (chr, position, BAF, LRR, etc...) or will use the column number or name
#' supplied by the user, if any. It will convert the tabular data into a
#' GRanges, with one range per SNP in the table.
#'
#' @usage loadSNPDataDNAcopy (snp.data, chr.col = "chrom", start.col = NULL, end.col = NULL, pos.col = "maploc", lrr.col = NULL, chr.transformation = "23:X,24:Y,25:MT", genome = NULL, na.rm = FALSE, verbose = TRUE)
#' 
#' @inheritParams loadSNPData
#' @param chr.transformation (character) The transformation of the chromosome 
#' names in a comma separated "key:value" format.(defaults to "23:X,24:Y,25:MT")
#' @param na.rm (logical) Whether to remove NA values. (defaults to FALSE)
#' 
#' @return
#' A list of GRanges with a range per SNP and one GRanges per sample.
#'
#' @examples
#' library(DNAcopy)
#' data(coriell)
#' #one sample
#' CNA.object <- CNA(cbind(coriell$Coriell.05296), coriell$Chromosome, coriell$Position, data.type="logratio")
#' CNA.object <-  smooth.CNA(CNA.object)
#' snps <- loadSNPDataDNAcopy(snp.data = CNA.object, na.rm = TRUE) 
#' 
#' # more than 1 sample
#' CNA.object <- CNA(genomdat = cbind(coriell$Coriell.05296, coriell$Coriell.13330), chrom = coriell$Chromosome, maploc = coriell$Position, data.type = "logratio", sampleid = c("c05296", "c13330"))
#' CNA.object <-  smooth.CNA(CNA.object)
#' snps <- loadSNPDataDNAcopy(snp.data = CNA.object, na.rm = TRUE)
#' 
#' #If the data come from DNAcopy results
#' DNAcopy.object <- segment(CNA.object, verbose=1)
#' snps <- loadSNPDataDNAcopy(snp.data = DNAcopy.object)
#'
#' @export loadSNPDataDNAcopy
#'

loadSNPDataDNAcopy <- function(snp.data,
                                chr.col = "chrom",
                                start.col = NULL,
                                end.col = NULL,
                                pos.col = "maploc",
                                lrr.col = NULL,
                                chr.transformation = "23:X,24:Y,25:MT",
                                genome = NULL,
                                na.rm = FALSE,
                                verbose = TRUE){
  
  if(!methods::is(snp.data, "DNAcopy") && !methods::is(snp.data, "CNA")){
    stop("snp.data must be either DNAcopy object or CNA")
  }
    
  #The results object from DNAcopy includes the input object. Extract it
  if(methods::is(snp.data, "DNAcopy")) snp.data <- data.frame(snp.data$data) 
  
  #If it's the input object from DNAcopy, extract the marker intensity values
  if(methods::is(snp.data, "CNA")) snp.data <- data.frame(snp.data) 
  
  #The LRR information is saved under sample name.
  samples <- colnames(snp.data)[!(colnames(snp.data) %in% c(chr.col, pos.col))]
  
  
  
  if(!is.null(chr.transformation)){
    chr.col <- names(snp.data)[getChrColumn(df = snp.data, col = chr.col, verbose = FALSE)]
    snp.data[,chr.col] <- transformChr(chr = snp.data[,chr.col], chr.transformation = chr.transformation)
  }

  if(length(samples)>=1){
    snps <- list()
    for(sample in samples){
      snp.data.samp <- snp.data[,c(chr.col, pos.col, sample)]
      snps[[sample]] <- loadSNPData(snp.data = snp.data.samp,
                                    chr.col = chr.col,
                                    start.col = start.col,
                                    end.col = end.col,
                                    pos.col = pos.col,
                                    baf.col = NULL, 
                                    lrr.col = sample, 
                                    id.col = NULL,
                                    genome = genome,
                                    verbose = verbose)
    }
  }
  #if the user wants to remove the NA values
  if(na.rm == TRUE) snps <- removeNAs(snp.data = snps, lrr.na = TRUE, baf.na = FALSE, id.na = FALSE, verbose = verbose)
  
  return(snps)
}
