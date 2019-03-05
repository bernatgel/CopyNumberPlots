# loadCopyNumberCalls.CNVkit example.
seg.file <-"testFiles/example.seg"
loadCopyNumberCalls.seg(seg.file)
# x <- read.table(file=ex.file,header = F, sep = "\t",quote = "\'")
# colnames(x) <- c("ID","chrom", "loc.start", "loc.end", "num.mark", "seg.mean")

loadCopyNumberCalls.seg <- function(seg.file, 
                                       chr.col = "chrom", 
                                       start.col = "loc.start",
                                       end.col = "loc.end", 
                                       segment.value.col = "seg.mean", 
                                       genome = NULL,
                                       # chr.transformation = "23:X,24:Y,25:MT",
                                       verbose = TRUE){
  #Load file.seg. First, probe if file exist
  if(!file.exists(seg.file)){
    stop(paste0(seg.file, " does not exist or you are not in the correct directory."))
    
  }else{
    seg.data <- read.table(file = seg.file, sep="\t", skip = 1, stringsAsFactors = FALSE)
    colnames(seg.file) <-  c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
    # we change the order of the columns to be transfromed in a G
    seg.data <- seg.file[,c("chrom", "loc.start", "loc.end", "ID", "num.mark", "seg.mean")] 
    
    ### añadir un parámetro que contenga por defecto 23:X,24:Y,25:MT 
    # chr.trans <- data.frame(limma::strsplit2(x=chr.transformation, split = ","))
    # chr.trans <- data.frame(limma::strsplit2(x=chr.trans, split = ":"))
    
    # una vez añadida la transformación lo que hago es convertir el seqlevelsStyle en UCSC.
  }
  
  segs <- loadCopyNumberCalls(cnv.data = seg.data, 
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

x<-loadCopyNumberCalls.CNVkit(cnvkit.file, genome = NULL)

