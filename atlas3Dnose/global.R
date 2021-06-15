library(shiny)
library(shinymanager)
library(reshape2)
library(plotly)
library(gtools)
library(ggplot2)
library(RColorBrewer)
library(locfit)
library(shinythemes)

#showNotification("loading data...")

load("3DMOE.Rdata")

plotXYcorrelation<-function(x, y, main, xlab, ylab, col, ylim=c(min(y, na.rm=T), max(y, na.rm=T)), xlim=c(min(x, na.rm=T), max(x, na.rm=T))){
  #x=values in x axis, y=values in y axis, xlab=x axis name, ylab=y axis name, col=vector of colors per sample
  plot(x, y, main=main, ylab=ylab, xlab=xlab, col=alpha(col, 0.3), cex.main=1, cex=0.8, pch=19, ylim=ylim, xlim=xlim)
  abline(lm(y~x))
  mtext(c(paste("R=", round(cor.test(x, y, method="spearman")$estimate, digits=2), "\t", "\t", "p=", cor.test(x, y, method="spearman")$p.value)), side=3, cex=0.6)
}

#plotXYcorrelation(x=as.numeric(DV["Acsm4",]), y=as.numeric(DV["Nqo1",]), col=1, ylab="Nqo1", 
#                  xlab="Acsm4", main="")
