#!/usr/bin/Rscript
#install.packages("seqinr")
#install.packages("data.table")
library("seqinr")
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
wd = args[1]
setwd(wd)

fasta.file = "454LargeContigs.fna"
fasta = read.fasta(fasta.file)
gc = sapply(fasta, function(x){GC(x)})
cov = fread("grep 'contig' 454ContigGraph.txt")
cov = cov[1:length(gc), ]

pdf("sim_portrait.pdf")
plot(cov$V4, gc, pch = 3, cex = 0.5, lwd = 1.5, xlab = "coverage", ylab = "GC",
     col = rgb(0.2,0.2,0.2,0.4), ylim=c(0,1), xlim = c(0,400))
plot(cov$V4, gc, xlab = "coverage", ylab = "GC", type="n")
text(cov$V4, gc, xlab = "coverage", ylab = "GC", labels = cov$V2,cex = 0.5,col = rgb(0.2,0.2,0.2,0.5))
dev.off()
