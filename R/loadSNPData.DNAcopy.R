

loadSNPData.DNAcopy <- function(snp.data) {
  if(methods::is(snp.data, "DNAcopy")) { #The results object from DNAcopy includes the input object. Extract it
    snp.data <- snp.data$data
  }
  if(methods::is(snp.data, "CNA")) { #If it's the input object from DNAcopy, extract the marker intensity values
    samples <- names(snp.data)
    samples <- samples[!(samples %in% c("chrom", "maploc"))]

    sample.data <- setNames(data.frame(snp.data[[samples]], stringsAsFactors = FALSE), samples)

    snps <-  toGRanges(cbind(chr=snp.data$chrom, start=snp.data$maploc, end=snp.data$maploc, sample.data, stringsAsFactors=FALSE))



    return(snps)
  }

}
