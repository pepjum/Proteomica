##### filtrar JAMM
library(data.table)


#antiguo$V11<-(log10(antiguo$V7)-mean(log10(antiguo$V7)))/sd(log10(antiguo$V7))


#qnorm(0.01)
# -2.326348


#selected<-antiguo[which(antiguo$V11 > abs(qnorm(0.01))),]



##### JAMM_EGR1

JAMM_EGR1<-fread("~/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/EGR1/JAMM_RESULTS/peaks/filtered.peaks.narrowPeak")

JAMM_EGR1$V11<-(log10(JAMM_EGR1$V7)-mean(log10(JAMM_EGR1$V7)))/sd(log10(JAMM_EGR1$V7))
JAMM$V11<-(log10(JAMM$V7)-mean(log10(JAMM$V7)))/sd(log10(JAMM$V7))


selected_JAMM<-JAMM_EGR1[which(JAMM_EGR1$V11 > abs(qnorm(0.01))),]
selected_JAMM<-as.data.frame(selected_JAMM)

write.table(selected_JAMM, file="~/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/EGR1/JAMM_RESULTS/peaks/selected_JAMM_EGR1.narrowPeak", col.names=F, row.names=F, quote=F, sep="\t")


###JAMM_h3k4me3

JAMM_h3k4me3<-fread("~/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k4me3/JAMM_RESULTS/peaks/filtered.peaks.narrowPeak")

JAMM_h3k4me3$V11<-(log10(JAMM_h3k4me3$V7)-mean(log10(JAMM_h3k4me3$V7)))/sd(log10(JAMM_h3k4me3$V7))

selected_JAMM<-JAMM_h3k4me3[which(JAMM_h3k4me3$V11 > abs(qnorm(0.01))),]
selected_JAMM<-as.data.frame(selected_JAMM)

write.table(selected_JAMM, file="~/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k4me3/JAMM_RESULTS/peaks/selected_JAMM_h3k4me3.narrowPeak", col.names=F, row.names=F, quote=F, sep="\t")

#### JAMM JUND

JAMM_JUND<-fread("~/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/JUND/JAMM_RESULTS/peaks/filtered.peaks.narrowPeak")

JAMM_JUND$V11<-(log10(JAMM_JUND$V7)-mean(log10(JAMM_JUND$V7)))/sd(log10(JAMM_JUND$V7))

selected_JAMM<-JAMM_JUND[which(JAMM_JUND$V11 > abs(qnorm(0.01))),]
selected_JAMM<-as.data.frame(selected_JAMM)

write.table(selected_JAMM, file="~/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/JUND/JAMM_RESULTS/peaks/selected_JAMM_JUND.narrowPeak", col.names=F, row.names=F, quote=F, sep="\t")

#### JAMM SP1

JAMM_SP1<-fread("~/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/SP1/JAMM_RESULTS/peaks/filtered.peaks.narrowPeak")

JAMM_SP1$V11<-(log10(JAMM_SP1$V7)-mean(log10(JAMM_SP1$V7)))/sd(log10(JAMM_SP1$V7))

selected_JAMM<-JAMM_SP1[which(JAMM_SP1$V11 > abs(qnorm(0.01))),]
selected_JAMM<-as.data.frame(selected_JAMM)

write.table(selected_JAMM, file="~/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/SP1/JAMM_RESULTS/peaks/selected_JAMM_SP1.narrowPeak", col.names=F, row.names=F, quote=F, sep="\t")

### JAMM h3k36me3 v2

JAMM_h3k36me3<-fread("~/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k36me3/JAMM_RESULTS_2/peaks/filtered.peaks.narrowPeak")

JAMM_h3k36me3$V11<-(log10(JAMM_h3k36me3$V7)-mean(log10(JAMM_h3k36me3$V7)))/sd(log10(JAMM_h3k36me3$V7))

selected_JAMM<-JAMM_h3k36me3[which(JAMM_h3k36me3$V11 > abs(qnorm(0.01))),]
selected_JAMM<-as.data.frame(selected_JAMM)

write.table(selected_JAMM, file="~/data/pepe/20_SEÑALES_SIGNALS_BED_METHODS/h3k36me3/JAMM_RESULTS_2/peaks/selected_JAMM_h3k36me3.narrowPeak", col.names=F, row.names=F, quote=F, sep="\t")
