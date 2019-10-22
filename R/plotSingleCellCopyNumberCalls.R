#' plotSingleCellCopyNumberCalls
#'
#' @description
#' Plot the segments representing the copy number calls from a single-cell
#' CNV data file. The file must be an HDF5 file with the same format as the
#' ones produced by 10X CellRanger software. 
#'
#' @details
#' This function will open the HDF5 file, extract the CNV values for each 
#' cell and plot them as a different colored segments using 
#' \code{\link{plotCopyNumberCalls}}. 
#' By default cells will be reordered and plotted according the hierarchical 
#' clustering tree embedded in the file. If reorder is FALSE or the tree is
#' not present, cell will be plotted in the order they are stored in the file.
#'
#' @note If the file is open by any other application the function will fail.
#' 
#' @usage plotSingleCellCopyNumberCalls(karyoplot, cnv.file, reorder=TRUE, cn.colors=NULL, track.margin=0, r0=0, r1=1, ...) 
#'
#' @param karyoplot A karyoplote object
#' @param cnv.file (an HDF5 file) The path to an HDF5 file containing the single-cell CNV data. 
#' @param reorder (logical) If TRUE, the cells will be plotted inthe order specified by the hirearchical clustering tree embedded in the HDF5 file. If FALSE the cells will be plotted in the order they are stored in the file. (defaukts to TRUE)
#' @param cn.colors (colors) The colors assigned to gains and losses. \code{\link{getCopyNumberColors}} is used to determine them. If NULL, the default color scheme is used. (defaults to NULL)
#' @param track.margin (numeric) This is the margin between the cells CN, in portion of the total per cell space. (deafults to 0)
#' @param r0 (numeric) (karyoploteR parameter) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)(defaults to 0)
#' @param r1 (numeric) (karyoploteR parameter) r0 and r1 define the vertical range of the data panel to be used to draw this plot. They can be used to split the data panel in different vertical ranges (similar to tracks in a genome browser) to plot differents data. If NULL, they are set to the min and max of the data panel, it is, to use all the available space. (defaults to NULL)(defaults to 1)
#' @param ... The ellipsis operator can be used to specify any additional graphical parameters. Any additional parameter will be passed to the internal calls to karyoploteR functions.
#'
#'
#' @return
#' Invisibly returns the karyoplot object representing the plot. With it
#' it is possible to add other elements to the plot using standrad karyoploteR
#' functions. The returned pbject will have an additional "latest.plot" 
#' element with a list containing: the number of cells, the hierarchical 
#' clustering tree, the bin size used to partition the genome, the windows
#' representing such partitioning, the number of bins per chromosome and a
#' GRAnges with the regions where no cell had any data, the no-call regions.
#' 
#' 
#'
#' @examples
#'
#' kp <- plotKaryotype(plot.type=4, genome="hg38")
#' #NOT RUN - Using 10X example data from https://www.10xgenomics.com/resources/datasets/
#' #plotSingleCellCopyNumberCalls(kp, "breast_tissue_D_2k_cnv_data.h5")
#' 
#'
#' @export plotSingleCellCopyNumberCalls
#'
#' @importFrom rhdf5 H5Fopen h5read H5Fclose
#' @importFrom GenomicRanges tileGenome seqnames
#' @importFrom stats setNames
#'
plotSingleCellCopyNumberCalls <- function(karyoplot, cnv.file, 
                                          reorder=TRUE, cn.colors=NULL, 
                                          track.margin=0, r0=0, r1=1, ...) {

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
  
  #Read the ztree into an hclust object
  cells.tree <- readHDF5Ztree(cnv.data)
  
  #Prepare a list to hold the no-call (cn=-128) regions
  no.call.regions <- list()
  
  #For each chromosome
  for(chr in shared.chr.names) {
    #Read the CN data of all cells
    chr.cnvs <- rhdf5::h5read(cnv.data, paste0("cnvs/", chr))[,seq_len(num.cells)]
    chr.cnvs[chr.cnvs == -128] <- NA
    chr.cnvs <- abs(chr.cnvs)
    
    chr.no.call.regions <- apply(chr.cnvs, 1, function(x) {all(is.na(x))})
    no.call.regions[[chr]] <- regioneR::joinRegions(windows[[chr]][chr.no.call.regions], min.dist = 1)
    
    #And plot it for each cell individually in the order dictated by the hierarchical cluestering tree
    for(i in seq_len(num.cells)) {
      at <- autotrack(i, num.cells, r0=r0, r1=r1, margin = track.margin)
      plotCopyNumberCalls(karyoplot, windows[[chr]], 
                          cn.values=chr.cnvs[,cells.tree$order[i]], 
                          cn.colors = cn.colors,
                          loh.height = 0, labels=NA,
                          r0=at$r0, r1=at$r1, ...)
    }
  }
  
  #add the latest.plot list with: the computed hclust, the windows, the number of cells...
  karyoplot$latest.plot <- list(funct="plotSingleCellCopyNumberCalls", 
                                computed.values=list(
                                  cells.tree=cells.tree,
                                  num.cells=num.cells,
                                  bin.size=cnv.data$constants$bin_size,
                                  bins.per.chr=setNames(cnv.data$constants$num_bins_per_chrom, cnv.data$constants$chroms),
                                  windows=windows,
                                  no.call.regions=no.call.regions
                                  )
                                )
  rhdf5::H5Fclose(cnv.data)
  
  return(invisible(karyoplot))
}

