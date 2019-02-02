library(dplyr)
library(DESeq2)
library(ggplot2)
library(mygene)
library(biomaRt)
library(tidyr)
library(purrr)

ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")

tcrdata <- read.delim("HedrickTCR.txt")
colnames(tcrdata) <- c("TranscriptID","Chr","Start","End","Strand","Length","Copies","Annotation","Rep1_01_PCC","Rep1_100_K99A","Rep1_10_K99A","Rep1_10_PCC","Rep1_No","Rep2_01_PCC_MEKi","Rep2_01_PCC","Rep2_100_K99A_MEKi","Rep2_100_K99A","Rep2_10_K99A_MEKi","Rep2_10_K99A","Rep2_10_PCC_MEKi","Rep2_10_PCC","Rep2_No_MEKi","Rep2_No")
noMEK <- cbind(tcrdata[,1:13], tcrdata[,15],tcrdata[,17],tcrdata[,19],tcrdata[,21],tcrdata[,23])
colnames(noMEK)[14:18] <- c(colnames(tcrdata)[15],colnames(tcrdata)[17],colnames(tcrdata)[19],colnames(tcrdata)[21],colnames(tcrdata)[23])
noMEK$mean_01_PCC <- rowMeans(noMEK[c(9,14)])
noMEK$mean_100_K99A <- rowMeans(noMEK[c(10,15)])
noMEK$mean_10_K99A <- rowMeans(noMEK[c(11,16)])
noMEK$mean_10_PCC <- rowMeans(noMEK[c(12,17)])
noMEK$mean_No <- rowMeans(noMEK[c(13,18)])
filt_noMEK <- noMEK[noMEK$mean_10_K99A >0,]
filt_noMEK <- filt_noMEK[filt_noMEK$mean_No >0,]
filt_noMEK <- filt_noMEK[filt_noMEK$mean_10_PCC >0,]
filt_noMEK$Fc_K99AtoNO <- filt_noMEK[21]/filt_noMEK[23]
filt_noMEK$Fc_PCCtoK99A <- filt_noMEK[22]/filt_noMEK[21]

#No peptide FPKM over 10, less than 1 Fc K99A to No, less than 1 Fc PCC to K99A
filt_noMEK <- filt_noMEK[filt_noMEK$mean_No > 10,]
filt_noMEK <- filt_noMEK[filt_noMEK$Fc_K99AtoNO < 0.8,]
filt_noMEK <- filt_noMEK[filt_noMEK$Fc_PCCtoK99A < 0.8,]

#Identify ensembl gene id for the filtered dataset
filtnM_ensmbl <- queryMany(filt_noMEK[1],scopes="refseq", fields="ensembl.gene", species="mouse")

#Read in filtered data from in vivo TCR stimulation RNA-seq
JKTCR <- read.csv("global_results_subset.csv")
colnames(JKTCR)[1] <- "ensembl"

#Reformat filtered Glass dataset and overlap with in vivo dataset
a <- as.data.frame(filtnM_ensmbl)
a <- a[4]
a <- unlist2(a)
a <- as.data.frame(a)
colnames(a) <- "ensembl"
overlap <- semi_join(JKTCR, a, by="ensembl")