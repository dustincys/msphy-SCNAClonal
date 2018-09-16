#!/usr/bin/env Rscript
library(latex2exp)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

rTableFilePath <- args[1]
copyNumberFigurePath <- args[2]


segmentsData <- read.table(rTableFilePath, head=T, sep = "\t")
segmentsData <- segmentsData[segmentsData$gc != 0,]
segmentsData$loga <- log(segmentsData$tReadNum + 1) - log(segmentsData$nReadNum + 1)

if('copyNumberAnswer' %in% colnames(segmentsData)){
  segmentsData$Predict <- segmentsData$copyNumber == segmentsData$copyNumberAnswer

  g <- ggplot(data = segmentsData)
  g <- g + geom_text(aes( x = gc, y = loga, label = copyNumberAnswer, color=Predict), show.legend = T, alpha=0.8, size=5)
  g <- g + theme_linedraw()
  g <- g + xlab("GC content") + ylab(TeX('$\\log (\\hat{D}^T/D^N)$'))
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

  png(copyNumberFigurePath, pointsize = 10.5, family = "Times",
      width = 5, height= 2.5, units="in", res=600)
  print(g)
  dev.off()
}else{
  g <- ggplot(data = segmentsData)
  g <- g + geom_text(aes( x = gc, y = loga, label = copyNumber), alpha=0.8, size=5)
  g <- g + theme_linedraw()
  g <- g + xlab("GC content") + ylab(TeX('$\\log (\\hat{D}^T/D^N)$'))
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

  png(copyNumberFigurePath, pointsize = 10.5, family = "Times",
      width = 5, height= 2.5, units="in", res=600)
  print(g)
  dev.off()
}
