#
# library(mclust)
#
# pos <- sort(floor(runif(1000, 1, 10000)))
# baf.data <- toGRanges("chr1", pos, pos)
# baf.data$baf <- rnorm(1000, mean = 0.5, sd = 0.06)
# baf.data$baf[1:400] <- baf.data$baf[1:400] + c(0.18, -0.08)
#
# kp <- plotKaryotype(zoom=toGRanges("chr1", 1, 10000))
# plotBAF(kp, baf.data)
# kpAbline(kp, v=start(baf.data[400]))
#
#
# for(i in 1:10) {
#   st <- (i-1)*100
#   bd <- baf.data[st:(st+100)]
#
#   kk <- Mclust(bd$baf)
#
#   #For every identified cluster
#   for(clust.num in seq_len(kk$G)) {
#     #plot the points
#     clust.bd <- bd[kk$classification==clust.num]
#     plotBAF(kp, clust.bd, points.col = rainbow(5)[clust.num])
#     med <- median(clust.bd$baf)
#     kpSegments(kp, data=joinRegions(clust.bd, 1e6), y0=med, y1=med, col=darker(rainbow(5)[clust.num]))
#   }
#
#
#   # kk <- kmeans(bd$baf, 2)
#   #
#   # if(kk$betweenss/kk$totss > 0.75) {
#   #   #It's two clusters
#   #   meds <- c(median(bd[kk$cluster==1]$baf), median(bd[kk$cluster==2]$baf))
#   #   kpPoints(kp, data=bd, y=meds[1], col="red")
#   #   kpPoints(kp, data=bd, y=meds[2], col="blue")
#   # } else {
#   #   #only one cluster
#   #   meds <- median(bd$baf)
#   #   kpPoints(kp, data=bd, y=meds, col="green")
#   # }
#
# }
#
#
#
#
# kp <- plotKaryotype(zoom=toGRanges("chr1", 1, 10000))
# plotBAF(kp, bd[kk$cluster==1], points.col = "red")
# plotBAF(kp, bd[kk$cluster==2], points.col = "blue")
#
# kk$betweenss
#
#
#
#
# wss <- numeric(15)
# for (i in 2:15) wss[i] <- sum(kmeans(bd$baf, centers=i)$withinss)
# plot(1:15, wss, type="b", xlab="Number of Clusters",
#      ylab="Within groups sum of squares")
#
#
#
#
#
# shift <- function(n, x) {
#   if(n==0) return(x)
#   if(n>0) return(c(rep(0, n), x[1:(length(x)-n)]))
#   if(n<0) return(c(x[(abs(n)+1):length(x)], rep(0, abs(n))))
# }
#
# mergeContiguous <- function(peaks, valleys) {
#   merged.valleys <- numeric()
#   pp <- c(0, peaks, 100)
#   for(p in seq_len(length(pp)-1)) {
#     v <- valleys[valleys>pp[p] & valleys<pp[p+1]]
#     if(length(v)>0) {
#       merged.valleys <- c(merged.valleys, floor(median(v)))
#     }
#   }
#   return(merged.valleys)
# }
#
# filterValleys <- function(peaks, valleys, ss) {
#   filtered.valleys <- numeric()
#   for(nv in seq_len(length(valleys))) {
#     v <- valleys[nv]
#     left.p <- peaks[peaks<v][length(peaks[peaks<v])]
#     right.p <- peaks[peaks>v][1]
#
#     if((length(left.p)==0 || ss[left.p]>ss[v]) && (length(right.p)==0 || ss[right.p]>ss[v])) {
#       filtered.valleys <- c(filtered.valleys, v)
#     }
#   }
#   return(filtered.valleys)
# }
#
#
# filterPeaks <- function(peaks, valleys, ss) {
#   filtered.peaks <- numeric()
#   for(np in seq_len(length(peaks))) {
#     p <- peaks[np]
#     left.v <- valleys[valleys<p][length(valleys[valleys<p])]
#     right.v <- valleys[valleys>p][1]
#
#     if((length(left.v)==0 || ss[left.v]<ss[p]) && (is.na(right.v) || ss[right.v]<ss[p])) {
#       filtered.peaks <- c(filtered.peaks, p)
#     }
#   }
#   return(filtered.peaks)
# }
#
#
#
# partitionGenome <- function(chr, pos, genome="hg19") {
#   chr <- as.character(chr)
#   gg <- getGenome(genome)
#   gg <- filterChromosomes(gg, keep.chr = unique(chr))
#   gr <- GRanges()
#   for(cc in unique(chr)) {
#     pp <- pos[chr==cc]
#     pp <- sort(unique(pp))
#     if(1 %in% pp) pp <- pp[-1]
#     gr <- c(gr, toGRanges(data.frame(chr=cc, start=c(1, pp), end=c(pp-1, end(gg[as.character(seqnames(gg))==cc])))))
#   }
#   gr <- setGenomeToGRanges(gr, genome)
#   return(sort(gr))
# }
#
#
#
# bd <- seg.snps
# hh <- hist(bd$baf, breaks = 100)$counts
#
# ss <- (shift(-5, hh) + shift(-4, hh) + shift(-3, hh) + shift(-2, hh) + shift(-1, hh) + hh + shift(1, hh) + shift(2, hh) + shift(3, hh) + shift(4, hh) + shift(5, hh))/11
#
# peaks <- which(ss != 0 &
#          ss>=shift(-5, ss) & ss>=shift(-4, ss)  & ss>=shift(-3, ss)  & ss>=shift(-2, ss)  & ss>=shift(-1, ss) &
#          ss>=shift(1, ss) & ss>=shift(2, ss) & ss>=shift(3, ss) & ss>=shift(4, ss) & ss>=shift(5, ss))
#
# valleys <- which(ss<=shift(-5, ss) & ss<=shift(-4, ss)  & ss<=shift(-3, ss)  & ss<=shift(-2, ss)  & ss<=shift(-1, ss) &
#                  ss<=shift(1, ss) & ss<=shift(2, ss) & ss<=shift(3, ss) & ss<=shift(4, ss) & ss<=shift(5, ss))
#
# #merge.valleys
# valleys <- mergeValleys(peaks, valleys)
# peaks <- mergePeaks(peaks, valleys)
#
# #Valleys partition the BAF space into clusters
#
#
# plot(x=1:99, y=ss, type="l")
# points(x=peaks, y=ss[peaks], col="blue")
# points(x=valleys, y=ss[valleys], col="red")
#
