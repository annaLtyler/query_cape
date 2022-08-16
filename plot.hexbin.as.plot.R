#This plots a hexbin plot using a regular plotting function so 
#we can get multiple plots on one page. I'm still not sure how
#do to it when plotting an actual hexbin object.

plot.hexbin.as.plot <- function(x, y, xlab, ylab, main, min.cex = 1, max.cex = 3,
    n.bins = 10, legend.pos = "topright", round.legend = 1000){

    if(missing(main)){main = ""}
	if(missing(xlab)){xlab = deparse(substitute(x))}
	if(missing(ylab)){ylab = deparse(substitute(y))}

    hb <- hexbin(x, y, xlab = xlab, ylab = ylab)
    #plot(hb)
    xy <- hcell2xy(hb)
    pt.cex = scale.between.vals(hb@count, target.min = min.cex, target.max = max.cex)
    cols <- colors.from.values(hb@count, use.pheatmap.colors = TRUE)
    
    plot(xy, xlab = "Main Effect of Query Locus", 
        ylab = "Interaction Effects", pch = 16,
        main = "Query Effects vs. Interaction Effects", col = cols, 
        cex = pt.cex)
    
    bins <- sort(unique(hb@count))
    binned.bins <- round(segment.region(min(bins), max(bins), n.bins, "ends"))
    rounded.bins <- unique(round(binned.bins/round.legend) * round.legend)
    rounded.bins[which(rounded.bins == 0)] <- 1

    bin.cols <- colors.from.values(rounded.bins, use.pheatmap.colors = TRUE)
    bin.cex <- scale.between.vals(rounded.bins, target.min = min.cex, target.max = max.cex)
    legend(legend.pos, legend = rounded.bins, col = bin.cols, pch = 16)

}