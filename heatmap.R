file<-commandArgs(trailingOnly=TRUE)[1:2]
# R --vanilla <this --args comparison.out label.out

library(RColorBrewer)

marPad<-15
fileL=length(readLines(file[[1]]))
x<-read.table(file[[1]], sep="\t", nrows=(fileL-2), header=T, row.names=1)
pdf(file=paste(file[[2]], ".pdf", sep=""), width=8, height=8, bg="white", pointsize=12)
heatmap(as.matrix(x), Rowv=NA, Colv=NA, na.rm=T, revC=T, scale="none", margins=c(log2(nrow(x))+marPad,log2(ncol(x))+marPad), cexRow=1, cexCol=1, col=colorRampPalette(rev(brewer.pal(9, "Reds")))(32) )
legend(x="bottomright", legend=c("match", "uncertain", "mismatch"), fill=colorRampPalette(rev(brewer.pal(3, "Reds")))(3) )
dev.off()
