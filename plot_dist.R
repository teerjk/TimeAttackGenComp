file<-commandArgs(trailingOnly=TRUE)[1:2]

x<-read.table(file[[1]])
pdf(file=paste(file[[2]], "af.pdf", sep="."), width=8, height=8, bg="white", pointsize=12)
hist(x[[8]], breaks=50, main=file[[2]], xlab="Alt allele fraction")
dev.off()

pdf(file=paste(file[[2]], "af.chr.pdf", sep="."), width=8, height=8, bg="white", pointsize=12)
plot(x[[8]], col=x[[1]], pch=".", cex=2, main=file[[2]], ylab="Position", xlab="Alt allele fraction")
dev.off()

