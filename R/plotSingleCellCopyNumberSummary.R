#' plotSingleCellCopyNumberSummary
#'
#' @description
#' Plot a summary of the copy number calls from a single-cell CNV data file. 
#' Plot a histogram-like plot with the number of cells with gains and losses. 
#' The file must be an HDF5 file with the same format as the
#' ones produced by 10X CellRanger software. 
#'
#' @details
#' This function will open the HDF5 file, extract the CNV values and counts,
#' for each genomic region, the number of cells with a gain or a loss. It 
#' then plots a histogram-like summary of these numbers.
#'
#' @note If the file is open by any other application the function will fail.
#' 
#' @usage plotSingleCellCopyNumberSummary(karyoplot, cnv.file, direction="in",  gain.color=NULL, normal.color=NULL, loss.color=NULL, r0=0, r1=1, ...)
#'
#' @param karyoplot A karyoplot object
#' @param cnv.file (an HDF5 file) The path to an HDF5 file containing the single-cell CNV data. 
#' @param direction The direction to which the coverage plot point, either "in" for inward or "out" for outward. (defaults to "in")
#' @param gain.color (color) The color assigned to gains (defaults to NULL)
#' @param normal.color (color) The color assigned to normal ploidy(defaults to NULL)
#' @param loss.color (colors) The color assigned to losses(defaults to NULL)
#' @param r0 (numeric) (karyoploteR parameter) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)(defaults to 0)
#' @param r1 (numeric) (karyoploteR parameter) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)(defaults to 1)
#' @param ... The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to karyoploteR functions.
#'
#'
#' @return
#' Invisibly returns the karyoplot object representing the plot. With it
#' it is possible to add other elements to the plot using standrad karyoploteR
#' functions. The returned object will have an additional "latest.plot" 
#' element with a list containing: the number of cells, the bin size used to 
#' partition the genome, the windows representing such partitioning with
#' the number of cells with gains and losses in that each window,
#' the number of bins per chromosome and a GRanges with the regions where 
#' no cell had any data, the no-call regions.
#' 
#' @examples
#'
#' kp <- plotKaryotype(plot.type=4, genome="hg38")
#' #NOT RUN - Using 10X example data from https://www.10xgenomics.com/resources/datasets/
#' #plotSingleCellCopyNumberSummary(kp, "breast_tissue_D_2k_cnv_data.h5")
#' 
#'
#' @export plotSingleCellCopyNumberSummary
#'
#' @importFrom GenomicRanges start end
#'
plotSingleCellCopyNumberSummary <- function(karyoplot, cnv.file,
                                            direction="in",  gain.color=NULL, 
                                            normal.color=NULL, loss.color=NULL, 
                                            r0=0, r1=1, ...) {

  if(missing(karyoplot)) stop("The parameter 'karyoplot' is required")
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  #data
  if(missing(cnv.file)) stop("The parameter 'cnv.file' is required")

  #Get the colors  
  if(any(is.null(gain.color), is.null(loss.color), is.null(normal.color))) {
    segment.colors <- getCopyNumberColors()
    if(is.null(gain.color)) gain.color <- segment.colors["3"]
    if(is.null(normal.color)) normal.color <- segment.colors["2"]
    if(is.null(loss.color)) loss.color <- segment.colors["1"]
  }
  
  #Open the file
  #TODO: Check it's a file and that it exists (can it be a connection?)
  cnv.data <- rhdf5::H5Fopen(cnv.file)
  
  #Check the genome size and the chromosomes match the one in the genome
    #1 - At least one of the chromosomes should be shared. If not, exit with a warning
      if(!any(cnv.data$constants$chroms %in% karyoplot$chromosomes)) {
        warning("None of the chromosomes in the single-cell CNV file is present 
in the plot. Nothing will be plot. 
Are you using the same genome used to create the data?
Chromosomes in the file: ", paste0(cnv.data$constants$chroms, collapse=", "), "
Chromosomes in the plot: ", paste0(karyoplot$chromosomes, collapse=", "))
        return(karyoplot)
      }
  
    #2 - Check the size of the shared chromosomes match
      shared.chrs <- which(cnv.data$constants$chroms %in% karyoplot$chromosomes)
      
      #Check if the chromosome sizes match those of the used genome reference
      if(!all(karyoplot$chromosome.lengths[cnv.data$constants$chroms[shared.chrs]] - (cnv.data$constants$num_bins_per_chrom[shared.chrs] * cnv.data$constants$bin_size) < cnv.data$constants$bin_size)) {
        stop("The size of the chromosomes in the data file is different than the size of the chromosomes in the genome. Are you using the same genome?")        
      }

  shared.chr.names <- cnv.data$constants$chroms[shared.chrs]
  
  #Create a tiling of the genome to plot the CN values
  windows <- GenomicRanges::tileGenome(karyoplot$chromosome.lengths[shared.chr.names],
                        tilewidth = cnv.data$constants$bin_size, cut.last.tile.in.chrom = TRUE)
  #and split the windows object by chr
  windows <- split(windows, GenomicRanges::seqnames(windows))
  
  #Get the total number of cells
  num.cells <- cnv.data$constants$num_cells
  
  #Prepare a list to hold the no-call (cn=-128) regions
  gains.and.losses <- list()
  
  #For each chromosome
  for(chr in shared.chr.names) {
    #Read the CN data of all cells
    chr.cnvs <- rhdf5::h5read(cnv.data, paste0("cnvs/", chr))[,seq_len(num.cells)]
    chr.cnvs[chr.cnvs == -128] <- NA
    chr.cnvs <- abs(chr.cnvs)
    
    chr.no.call.regions <- apply(chr.cnvs, 1, function(x) {all(is.na(x))})
    no.call.regions[[chr]] <- regioneR::joinRegions(windows[[chr]][chr.no.call.regions], min.dist = 1)
    
    valid.regions <- regioneR::subtractRegions(karyoplot$genome[chr], no.call.regions[[chr]])
    
    #compute the number of gains and losses
    w <- windows[[chr]]
    w$gains <- apply(chr.cnvs, 1, function(x) {length(which(x>2))})
    w$normal <- apply(chr.cnvs, 1, function(x) {length(which(x==2))})
    w$losses <- apply(chr.cnvs, 1, function(x) {length(which(x<2))})
    
    #To get kpArea to plot the real coverage (flat tops), we need to build a
    #GRanges with two elements per range, one at the start and one at the end
    starts <- w
    GenomicRanges::end(starts) <- GenomicRanges::start(starts)
    ends <- w
    GenomicRanges::start(ends) <- GenomicRanges::end(ends)
    w <- sort(c(starts, ends))
    
    gains.and.losses[[chr]] <- w
    
    if(direction=="in") {
      karyoploteR::kpRect(karyoplot, data=valid.regions, y0=0, y1=1, r0=r0, r1=r1, col=normal.color, ...)
      karyoploteR::kpArea(karyoplot, data=w, y=w$gains, col=gain.color,  r0=r1, r1=r0, ymax=num.cells, ...)
      karyoploteR::kpArea(karyoplot, data=w, y=w$losses, col=loss.color, r0=r0, r1=r1, ymax=num.cells, ...)
    } else {
      mid <- r0+(r1-r0)/2
      karyoploteR::kpArea(karyoplot, data=w, y=w$gains, col=gain.color,  r0=mid, r1=r1, ymax=num.cells, ...)
      karyoploteR::kpArea(karyoplot, data=w, y=w$losses, col=loss.color, r0=mid, r1=r0, ymax=num.cells, ...)
    }
    
  }
  
  #add the latest.plot list with: the computed hclust, the windows, the number of cells...
  karyoplot$latest.plot <- list(funct="plotSingleCellCopyNumberSummary", 
                                computed.values=list(
                                  num.cells=num.cells,
                                  bin.size=cnv.data$constants$bin_size,
                                  bins.per.chr=setNames(cnv.data$constants$num_bins_per_chrom, cnv.data$constants$chroms),
                                  gains.and.losses=gains.and.losses,
                                  no.call.regions=no.call.regions
                                  )
                                )

  rhdf5::H5Fclose(cnv.data)
  
    
  return(invisible(karyoplot))
}

