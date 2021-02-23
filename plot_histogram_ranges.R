 plotter <- results_MSFragger_NX[, c('Peptide', 'Experimental.Mass')]
 plotter$ranges <- cut(as.numeric(plotter$Experimental.Mass), seq(0,max(plotter$Experimental.Mass),200))
plotter$labels <- cut(as.numeric(plotter$Experimental.Mass), seq(0,max(as.numeric(plotter$Experimental.Mass)),200), labels=c(1:346))
hist(plotter$labels)
dev.off()


subset["peptide","original.delta.mass"]
library(doby)



aa<-results_MSFragger_NX[(results_MSFragger_NX$Original.Delta.Mass < 50) && (results_MSFragger_NX$Original.Delta.Mass > -1*(50)),] ### rango

bb <- summaryBy(Original.Delta.Mass ~ Spectrum, data = aa, FUN = function(x) length(unique(x)))

bb[,1] <- as.numeric(paste(bb[,1]))
bb[,2] <- as.numeric(paste(bb[,2]))
plot(bb[,1], bb[,2], type = "p", pch = ".")
dev.off()
