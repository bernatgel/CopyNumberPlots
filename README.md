[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)


# CopyNumberPlots: create copy-number plots using karyoploteR 

## Introduction

Data visualisation is a powerful tool used for data analysis and exploration in 
many fields. Genomics data analysis is one of these fields where good 
visualisation tools can be of great help. 
The aim of [CopyNumberPlots](http://bioconductor.org/packages/CopyNumberPlots/)
is to offer the user an 
easy way to create copy-number related plots using the infrastructure provided
by the R/Bioconductor package [karyoploteR](http://bioconductor.org/packages/karyoploteR/)

In addition to a set of specialized plotting functions for copy-number analysis
data and results (`plotBAF`, `plotCopyNumberCalls`, ...), 
CopyNumberPlots contains a number of data loading 
functions to help parsing and loading the results of widely used 
copy-number calling software such as [DNAcopy](	http://bioconductor.org/packages/DNAcopy/),
[DECoN](https://wellcomeopenresearch.org/articles/1-20/v1) or
[CNVkit](https://cnvkit.readthedocs.io/en/stable/).



Finally, since CopyNumberPlots extends the 
functionality of karyoploteR, it is possible 
to combine the plotting functions of both packages to get the perfect
figure for your data.

## Installation

CopyNumberPlots is a [Bioconductor](http://bioconductor.org) 
package and to install it we have to use 
[BiocManager](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html)

```{r getPackage, eval=FALSE}
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("CopyNumberPlots")
```

## More information

You can find more information on CopyNumberPlots in 
[its Bioconductor page](http://bioconductor.org/packages/CopyNumberPlots/) 
and specially in the 
[CopyNumberPlots vignette](https://bioconductor.org/packages/release/bioc/vignettes/CopyNumberPlots/inst/doc/CopyNumberPlots.html).

You can find more information on how to use karyoploteR in the
[karyoploteR tutorial](https://bernatgel.github.io/karyoploter_tutorial/).



