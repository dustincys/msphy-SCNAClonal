#!/usr/bin/env Rscript
library(latex2exp)
library(ggplot2)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
rTableFilePath <- args[1]
phiFigurePath <- args[2]

segmentsData <- read.table(rTableFilePath, head=T, sep = "\t")
segmentsData <- segmentsData[segmentsData$gc != 0,]
segmentsData$loga <- log(segmentsData$tReadNum + 1) - log(segmentsData$nReadNum + 1)
segmentsData <- segmentsData[segmentsData$gc != 0,]
segmentsData$length <- segmentsData$end - segmentsData$start

if('phiAnswer' %in% colnames(segmentsData)){
  segmentsData$Predict <- segmentsData$copyNumber == segmentsData$copyNumberAnswer
}

for(chr in unique(segmentsData$chrom)){
  tempSD <- segmentsData[segmentsData$chrom == chr,]
  tempSD <- tempSD[order(tempSD$start),]
  tempSD$mappedLength <- tempSD$length / min(tempSD$length)

  tempSD <- data.table(tempSD)
  tempSD <- tempSD[, mappedEnd:=cumsum(mappedLength),]
  tempSD$mappedStart <- tempSD$mappedEnd - tempSD$mappedLength

  g <- ggplot(data = tempSD)
  if('phiAnswer' %in% colnames(segmentsData)){
    g <- g + geom_segment(aes(x = mappedStart, y = phiAnswer, xend = mappedEnd,
                              yend = phiAnswer,color='Ground truth'),size=1,alpha=1)
    g <- g + geom_segment(aes(x = mappedStart, y = phi, xend = mappedEnd, yend = phi,
                              color='Predicted'),size=1,alpha=0.5)
  }else{
    g <- g + geom_segment(aes(x = mappedStart, y = phi, xend = mappedEnd, yend = phi),
                          size=1,alpha=0.5)
  }
  #g <- g + scale_x_continuous(breaks=c(50000,200000,350000),
  #                     labels=c("1","Baseline","2"))
  g <- g + theme_linedraw()
  g <- g + xlab("SCNA") + ylab(TeX('$\\phi$'))
  g <- g + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key = element_rect(fill = "white"),
    legend.title=element_text(size=10.5),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black"),
    axis.text=element_text(size=10.5),
    axis.title=element_text(size=10.5),
    axis.text.x = element_text(size=10.5),
    axis.text.y = element_text(size=10.5),
    plot.title = element_text(size=10.5))
  
  png(paste(phiFigurePath, "_chr", chr, ".png", sep = ""), pointsize = 10.5, family = "Times", width = 8, height= 2.5, units="in", res=600)
  print(g)
  dev.off()
}
