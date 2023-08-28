
##################################################################
#load the packages
library(edgeR)
library(DESeq2)
library(limma)
library(org.Hs.eg.db)
library(tximport)
library(csaw)
library(readr)
library(ggplot2)
library(ggrepel)
library(WGCNA)
library(ggpubr)
library(bestNormalize)
library(ComplexHeatmap)
library(reshape2)
library(tidyverse)
library(caret)
library(MASS)
library(scatterplot3d)
library(GenomicAlignments)
library(GenomicFeatures)
library(org.Pf.plasmo.db)
library(variancePartition)
library(elasticnet)
library("FactoMineR")
library("factoextra")
library(Seurat)
library(patchwork)
library(cowplot)
library(mixOmics)
library(readxl)
library("DiscoRhythm")
library(tximport)
library(csaw)
library(openxlsx)
library(matrixTests)
library(ISOpureR)
library(cqn)
library(NormExpression)
library(sva)
library(devtools)
library(QuantNorm)
library(scBatch)

#vst()
##make a tx2gene file for plasmodium
txdb<-makeTxDbFromGFF(file = "~/Documents/phd_after_upgrading_rnaseq/12_timepoints/mapped2pf3D7/PlasmoDB-55_Pfalciparum3D7.gff",
                      format = "gff",)
k <- keys(txdb, keytype = "TXNAME")
keytypes(txdb)
keys(txdb,keytype = "TXID")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
head(tx2gene)
dim(tx2gene)

###################################################################################################
annot<-read.xlsx("~/Documents/phd_after_upgrading_rnaseq/12_timepoints/mapped2pf3D7/Pf3D7_annotation_file.xlsx")
head(annot)
annot<-annot[!duplicated(annot$geneID),]
library(dplyr)
annot_filter<-annot%>%dplyr::filter(grepl(pattern = "PfEMP1|rifin|stevor|ribosomal RNA",geneName))
dim(annot_filter)
annot_filter$geneName
head(annot_filter)
annot_filter2<-read.xlsx("~/Documents/phd_after_upgrading_rnaseq/12_timepoints/mapped2pf3D7_version2/Pfalciparum_amino.acyltRNAs.xlsx")


######################################################################################################
annot_filter_bind<-rbind(annot_filter2,annot_filter)
annot_filter_bind<-annot_filter_bind[!duplicated(annot_filter_bind$geneID),]
head(annot_filter_bind)
dim(annot_filter_bind)

##########################################################
tx2gene<-tx2gene[!tx2gene$GENEID%in%annot_filter_bind$geneID,]
dim(tx2gene)

#set the path to all Kallisto files
files_all_data <- dir(path = "kallisto_output", 
                      pattern="*.h5", 
                      full.names =TRUE,recursive=TRUE)
files_all_data


merged_meta<-read.csv("~/Documents/phd_after_upgrading_rnaseq/12_timepoints/mapped2pf3D7_version2/all_phenodata.csv",
                      row.names = 2)
head(merged_meta)
Group <- factor(paste(merged_meta$condition,merged_meta$timepoint,sep="_"))
Group
merged_meta<-cbind(merged_meta,Group=Group)
merged_meta


head(files_all_data)
subjectid<-gsub(pattern = "/home/kmwikali/Documents/phd_after_upgrading_rnaseq/12_timepoints/mapped2pf3D7_version2/all_kallisto_files_version2/",
            replacement = "", x = files_all_data)
subjectid<-gsub(pattern = "/abundance.h5",
                replacement = "", x = subjectid)
subjectid


names(files_all_data)<-paste0(subjectid)
files_all_data



#make a tximport object
Retinopathy_allData_kallisto_tximport <- tximport(files_all_data, type = "kallisto",
                               txOut=FALSE,
                               tx2gene = tx2gene,
                               txIdCol = "TXNAME",
                               geneIdCol = "GENEID",
                               countsFromAbundance="lengthScaledTPM", 
                               ignoreAfterBar = TRUE)
save(Retinopathy_allData_kallisto_tximport, file="Retinopathy_allData_kallisto_tximport.Rdata")

#start the process of normalizing the counts by gene length and lib size
dim(kallisto_genelevel$abundance) #5317
head(kallisto_genelevel$abundance)
head(kallisto_genelevel$counts)



length_gc <- read.csv("~/Documents/phd_after_upgrading_rnaseq/12_timepoints/mapped2pf3D7/Pfalciparum_length_pGC_content.txt", sep="")
head(length_gc)
dim(length_gc)
length_gc$Name<-gsub("\\..*","", length_gc$Name)
head(length_gc)
names(length_gc)[1:3]<-c("IDs","length","pGC")
length_gc<-length_gc[!duplicated(length_gc$IDs),]
rownames(length_gc)<-length_gc$IDs
dim(length_gc)
head(length_gc)
length_gc$PropGC<-length_gc$pGC/100

library(cqn)
library(NormExpression)

inter_genes<-intersect(rownames(length_gc),rownames(kallisto_genelevel$counts))
inter_genes

###########################################################################################
#calculate tpm
tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}


TPMdata<-tpm3(counts = kallisto_genelevel$counts[inter_genes,],
              len = length_gc[inter_genes,]$length)
head(TPMdata)
dim(TPMdata)

###########################################################################################
nTPMdata<-(TPMdata+0.001)
head(nTPMdata)
nTPMdata<-normalizeVSN(nTPMdata)
truehist(nTPMdata)

############################################################################################
merged_meta<-merged_meta[order(merged_meta$timepoint),]
merged_meta<-merged_meta[order(merged_meta$Isolate),]
Merged_EV<-subset(merged_meta, condition=="EV")
Merged_EV
Merged_WP<-subset(merged_meta, condition=="WP")
Merged_WP

nTPM_EV<-as.data.frame((nTPMdata[,rownames(Merged_EV)]))
nTPM_EV_s<-t(scale(t(nTPM_EV)))
head(nTPM_EV_s)
#calculate signal to noise ratio
nTPM_EV_mvfft<-mvfft(as.matrix(t(nTPM_EV)))
Heatmap(t(scale(Re(nTPM_EV_mvfft[-1,]))), show_row_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,column_labels = 2:70)
nTPM_EV_mvfft[-c(1,7),]=0+0i
nTPM_EV_imvfft<-mvfft(nTPM_EV_mvfft, inverse = TRUE)/70
t(Re(nTPM_EV_imvfft))
dim(na.omit(t(nTPM_EV_imvfft)))
##########################################################################
#extract the fourier normalized data
nTPM_EV_e_imvfft<-na.omit(as.data.frame(t(Re(nTPM_EV_imvfft))))
head(nTPM_EV_e_imvfft)


nTPM_WP<-as.data.frame((nTPMdata[,rownames(Merged_WP)]))
nTPM_WP_s<-t(scale(t(nTPM_WP)))
head(nTPM_WP_s)
#calculate signal to noise ratio
nTPM_WP_mvfft<-mvfft(as.matrix(t(nTPM_WP)))
Heatmap(t(scale(Re(nTPM_WP_mvfft[-1,]))), show_row_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,column_labels = 2:70)
nTPM_WP_mvfft[-c(1,7),]=0+0i
nTPM_WP_imvfft<-mvfft(nTPM_WP_mvfft, inverse = TRUE)/70
t(Re(nTPM_WP_imvfft))
dim(na.omit(t(nTPM_WP_imvfft)))
##########################################################################
#extract the fourier normalized data
nTPM_WP_p_imvfft<-na.omit(as.data.frame(t(Re(nTPM_WP_imvfft))))
head(nTPM_WP_p_imvfft)


library(itsmr)
###########################################################################
#plot EBA175
plotc((as.data.frame(t(nTPM_WP_p_imvfft))$PF3D7_0731500),
      (as.data.frame(t(nTPM_EV_e_imvfft))$PF3D7_0731500))

#msp2
plotc((as.data.frame(t(nTPM_WP_p_imvfft))$PF3D7_0206800),
      (as.data.frame(t(nTPM_EV_e_imvfft))$PF3D7_0206800))

#kahrp
plotc((as.data.frame(t(nTPM_WP_p_imvfft))$PF3D7_1353200),
      (as.data.frame(t(nTPM_EV_e_imvfft))$PF3D7_1353200))

#msp1
plotc((as.data.frame(t(nTPM_WP_p_imvfft))$PF3D7_0930300),
      (as.data.frame(t(nTPM_EV_e_imvfft))$PF3D7_0930300))


#rh5
plotc(scale(as.data.frame(t(nTPM_WP_p_imvfft))$PF3D7_0424100),
      scale(as.data.frame(t(nTPM_EV_e_imvfft))$PF3D7_0424100))

plotc((as.data.frame(t(nTPM_WP_p_imvfft))$PF3D7_0304600),
      (as.data.frame(t(nTPM_EV_e_imvfft))$PF3D7_0304600))
cor((as.data.frame(t(nTPM_WP_p_imvfft))$PF3D7_0304600),
    (as.data.frame(t(nTPM_EV_e_imvfft))$PF3D7_0304600))

#############################################################################################
head(nTPM_WP_p_imvfft)
head(nTPM_EV_e_imvfft)

#correlate EVs and parasites
pf<-as.data.frame(t(nTPM_WP_p_imvfft))
pf[,1:4]

ev<-as.data.frame(t(nTPM_EV_e_imvfft[rownames(nTPM_WP_p_imvfft),]))
ev[,1:4]


rho <- pval <- NULL
for(i in 1:ncol(pf)){
  rho[i]  <- cor.test(pf[,i],ev[,i], method = "pearson")$estimate
  pval[i] <- cor.test(pf[,i],ev[,i], method = "pearson")$p.value
}
res1 <- cbind(rho,pval)
res1



rownames(res1)<-colnames(pf)

res1<-as.data.frame(res1)
head(res1)

res1$IDs<- rownames(res1)
head(res1)


#######################################################################################
####plothistogram to choose cut-off

ggplot(res1, aes(x=rho))+
  geom_histogram(color="darkblue", fill="lightblue")


cosinor_WP<-row_cosinor(nTPM_WP_p_imvfft,
                        t = Merged_WP$timepoint,
                        period = 48)
head(cosinor_WP)


cosinor_EV<-row_cosinor(nTPM_EV_e_imvfft,
                        t = Merged_EV$timepoint,
                        period = 48)
head(cosinor_EV)

plot(cosinor_EV$acrophase, cosinor_WP$acrophase)

##########################################################################################

phase1<-abs(cosinor_EV$acrophase - cosinor_WP$acrophase)
phase1
phase2<-48-abs(cosinor_EV$acrophase - cosinor_WP$acrophase)
phase2

phase_diff<-pmin(phase1,phase2)
table(phase_diff>17)

cosinor_EV$phase_diff<-phase_diff
cosinor_WP$phase_diff<-phase_diff
head(cosinor_WP)

cosinor_EV_sig<-subset(cosinor_EV, phase_diff>12)
dim(cosinor_EV_sig)
pc1<-pca(t(cbind(t(scale(t(nTPM_EV_e_imvfft))),t(scale(t(nTPM_WP_p_imvfft))))[rownames(cosinor_EV_sig),][,rownames(merged_meta)]),
         center = FALSE)
plotIndiv(pc1, ind.names = FALSE, group = merged_meta$timepoint)

comp<-pc1$variates$X[,1:2]
all(rownames(comp)==rownames(merged_meta))

merged_meta_comp<-cbind(comp,merged_meta)
head(merged_meta_comp)

design<-model.matrix(~(PC1+PC2), merged_meta_comp)
design

################################################################################################################
fit<-lmFit(cbind(t(scale(t(nTPM_EV))),t(scale(t(nTPM_WP))))
           [,rownames(design)],design = design)
resid<-residuals(object = fit,cbind(t(scale(t(nTPM_EV))),t(scale(t(nTPM_WP))))
                 [,rownames(design)],type=deviance)


library(RUVSeq)
###############################################################################
ntpmdata<-RUVr(x =cbind(t(scale(t(nTPM_EV))),t(scale(t(nTPM_WP))))
                        [,rownames(design)],isLog = TRUE,
                        cIdx = rownames(nTPM_WP_p_imvfft),k = 50,residuals = resid)
ntpmdata$normalizedCounts


pc1<-pca(t(cbind((ntpmdata$normalizedCounts))[rownames(cosinor_EV),][,rownames(merged_meta_comp)]),
         center = FALSE)
plotIndiv(pc1, ind.names = FALSE, legend = TRUE,
          group = merged_meta_comp$timepoint)

#

#correlate EVs and parasites
pf<-as.data.frame(t(ntpmdata$normalizedCounts[,rownames(Merged_WP)]))
pf[,1:4]

ev<-as.data.frame(t(ntpmdata$normalizedCounts[,rownames(Merged_EV)][rownames(ntpmdata$normalizedCounts[,rownames(Merged_WP)]),]))
ev[,1:4]


rho <- pval <- NULL
for(i in 1:ncol(pf)){
  rho[i]  <- cor.test(pf[,i],ev[,i], method = "pearson")$estimate
  pval[i] <- cor.test(pf[,i],ev[,i], method = "pearson")$p.value
}
res1 <- cbind(rho,pval)
res1



rownames(res1)<-colnames(pf)

res1<-as.data.frame(res1)
head(res1)

res1$IDs<- rownames(res1)
head(res1)



#############################################################################
#order the phenodata by
pheno_9106<-subset(merged_meta, Isolate=="KE01")
pheno_9106<-pheno_9106[order(pheno_9106$timepoint),]
pheno_9106
#select the whole parasite data
pheno_9106e<-subset(pheno_9106, condition=="EV")
pheno_9106e

nTPM_EV_9106<-as.data.frame((nTPMdata[,rownames(pheno_9106e)]))
nTPM_EV_9106_mvfft<-mvfft(as.matrix(t(nTPM_EV_9106)))
Heatmap(t(scale(Re(nTPM_EV_9106_mvfft))), show_row_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,column_labels = 1:11)
nTPM_EV_9106_mvfft[-c(1),]=0+0i
nTPM_EV_9106_imvfft<-mvfft(nTPM_EV_9106_mvfft, inverse = TRUE)/11
t(Re(nTPM_EV_9106_imvfft))
dim(na.omit(t(nTPM_EV_9106_imvfft)))

pheno_9106p<-subset(pheno_9106, condition=="WP")
pheno_9106p

############################################################################
nTPM_WP_9106<-as.data.frame((nTPMdata[,rownames(pheno_9106p)]))
##########################################################################
#extract the fourier normalized data
nTPM_EV_9106_e_imvfft<-na.omit(as.data.frame(t(Re(nTPM_EV_9106_imvfft))))
head(nTPM_EV_9106_e_imvfft)

#####################################################################################
#calculate signal to noise ratio
nTPM_WP_9106_mvfft<-mvfft(as.matrix(t(nTPM_WP_9106)))
Heatmap(t(scale(Re(nTPM_WP_9106_mvfft))), show_row_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,column_labels = 1:11)
nTPM_WP_9106_mvfft[-c(1),]=0+0i
nTPM_WP_9106_imvfft<-mvfft(nTPM_WP_9106_mvfft, inverse = TRUE)/11
t(Re(nTPM_WP_9106_imvfft))
dim(na.omit(t(nTPM_WP_9106_imvfft)))

nTPM_WP_9106_p_imvfft<-na.omit(as.data.frame(t(Re(nTPM_WP_9106_imvfft))))
head(nTPM_WP_9106_p_imvfft)


nTPM_mean9106<-cbind(nTPM_WP_9106_p_imvfft, nTPM_EV_9106_e_imvfft)
head(nTPM_mean9106)

###########################################################################
#order the phenodata by
pheno_9106HB<-subset(merged_meta, Isolate=="sKE01")
pheno_9106HB<-pheno_9106HB[order(pheno_9106HB$timepoint),]
pheno_9106HB
#select the whole parasite data
pheno_9106HBe<-subset(pheno_9106HB, condition=="EV")
pheno_9106HBe

nTPM_EV_9106HB<-as.data.frame((nTPMdata[,rownames(pheno_9106HBe)]))
nTPM_EV_9106HB_mvfft<-mvfft(as.matrix(t(nTPM_EV_9106HB)))
Heatmap(t(scale(Re(nTPM_EV_9106HB_mvfft))), show_row_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,column_labels = 1:12)
nTPM_EV_9106HB_mvfft[-c(1),]=0+0i
nTPM_EV_9106HB_imvfft<-mvfft(nTPM_EV_9106HB_mvfft, inverse = TRUE)/12
t(Re(nTPM_EV_9106HB_imvfft))
dim(na.omit(t(nTPM_EV_9106HB_imvfft)))

pheno_9106HBp<-subset(pheno_9106HB, condition=="WP")
pheno_9106HBp

############################################################################
nTPM_WP_9106HB<-as.data.frame((nTPMdata[,rownames(pheno_9106HBp)]))
##########################################################################
#extract the fourier normalized data
nTPM_EV_9106HB_e_imvfft<-na.omit(as.data.frame(t(Re(nTPM_EV_9106HB_imvfft))))
head(nTPM_EV_9106HB_e_imvfft)

#####################################################################################
#calculate signal to noise ratio
nTPM_WP_9106HB_mvfft<-mvfft(as.matrix(t(nTPM_WP_9106HB)))
Heatmap(t(scale(Re(nTPM_WP_9106HB_mvfft))), show_row_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,column_labels = 1:12)
nTPM_WP_9106HB_mvfft[-c(1),]=0+0i
nTPM_WP_9106HB_imvfft<-mvfft(nTPM_WP_9106HB_mvfft, inverse = TRUE)/11
t(Re(nTPM_WP_9106HB_imvfft))
dim(na.omit(t(nTPM_WP_9106HB_imvfft)))

nTPM_WP_9106HB_p_imvfft<-na.omit(as.data.frame(t(Re(nTPM_WP_9106HB_imvfft))))
head(nTPM_WP_9106HB_p_imvfft)


nTPM_mean9106HB<-cbind(nTPM_WP_9106HB_p_imvfft, nTPM_EV_9106HB_e_imvfft)
head(nTPM_mean9106HB)

###########################################################################
#order the phenodata by
pheno_9626<-subset(merged_meta, Isolate=="KE02")
pheno_9626<-pheno_9626[order(pheno_9626$timepoint),]
pheno_9626
#select the whole parasite data
pheno_9626e<-subset(pheno_9626, condition=="EV")
pheno_9626e

nTPM_EV_9626<-as.data.frame((nTPMdata[,rownames(pheno_9626e)]))
nTPM_EV_9626_mvfft<-mvfft(as.matrix(t(nTPM_EV_9626)))
Heatmap(t(scale(Re(nTPM_EV_9626_mvfft))), show_row_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,column_labels = 1:12)
nTPM_EV_9626_mvfft[-c(1),]=0+0i
nTPM_EV_9626_imvfft<-mvfft(nTPM_EV_9626_mvfft, inverse = TRUE)/12
t(Re(nTPM_EV_9626_imvfft))
dim(na.omit(t(nTPM_EV_9626_imvfft)))

pheno_9626p<-subset(pheno_9626, condition=="WP")
pheno_9626p

############################################################################
nTPM_WP_9626<-as.data.frame((nTPMdata[,rownames(pheno_9626p)]))
##########################################################################
#extract the fourier normalized data
nTPM_EV_9626_e_imvfft<-na.omit(as.data.frame(t(Re(nTPM_EV_9626_imvfft))))
head(nTPM_EV_9626_e_imvfft)

#####################################################################################
#calculate signal to noise ratio
nTPM_WP_9626_mvfft<-mvfft(as.matrix(t(nTPM_WP_9626)))
Heatmap(t(scale(Re(nTPM_WP_9626_mvfft))), show_row_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,column_labels = 1:12)
nTPM_WP_9626_mvfft[-c(1),]=0+0i
nTPM_WP_9626_imvfft<-mvfft(nTPM_WP_9626_mvfft, inverse = TRUE)/12
t(Re(nTPM_WP_9626_imvfft))
dim(na.omit(t(nTPM_WP_9626_imvfft)))

nTPM_WP_9626_p_imvfft<-na.omit(as.data.frame(t(Re(nTPM_WP_9626_imvfft))))
head(nTPM_WP_9626_p_imvfft)


nTPM_mean9626<-cbind(nTPM_WP_9626_p_imvfft, nTPM_EV_9626_e_imvfft)
head(nTPM_mean9626)

#########################################################################
########################################################################

#order the phenodata by
pheno_10668<-subset(merged_meta, Isolate=="KE04")
pheno_10668<-pheno_10668[order(pheno_10668$timepoint),]
pheno_10668
#select the whole parasite data
pheno_10668e<-subset(pheno_10668, condition=="EV")
pheno_10668e

nTPM_EV_10668<-as.data.frame((nTPMdata[,rownames(pheno_10668e)]))
nTPM_EV_10668_mvfft<-mvfft(as.matrix(t(nTPM_EV_10668)))
Heatmap(t(scale(Re(nTPM_EV_10668_mvfft))), show_row_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,column_labels = 1:12)
nTPM_EV_10668_mvfft[-c(1),]=0+0i
nTPM_EV_10668_imvfft<-mvfft(nTPM_EV_10668_mvfft, inverse = TRUE)/12
t(Re(nTPM_EV_10668_imvfft))
dim(na.omit(t(nTPM_EV_10668_imvfft)))

pheno_10668p<-subset(pheno_10668, condition=="WP")
pheno_10668p

############################################################################
nTPM_WP_10668<-as.data.frame((nTPMdata[,rownames(pheno_10668p)]))
##########################################################################
#extract the fourier normalized data
nTPM_EV_10668_e_imvfft<-na.omit(as.data.frame(t(Re(nTPM_EV_10668_imvfft))))
head(nTPM_EV_10668_e_imvfft)

#####################################################################################
#calculate signal to noise ratio
nTPM_WP_10668_mvfft<-mvfft(as.matrix(t(nTPM_WP_10668)))
Heatmap(t(scale(Re(nTPM_WP_10668_mvfft))), show_row_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,column_labels = 1:12)
nTPM_WP_10668_mvfft[-c(1),]=0+0i
nTPM_WP_10668_imvfft<-mvfft(nTPM_WP_10668_mvfft, inverse = TRUE)/12
t(Re(nTPM_WP_10668_imvfft))
dim(na.omit(t(nTPM_WP_10668_imvfft)))

nTPM_WP_10668_p_imvfft<-na.omit(as.data.frame(t(Re(nTPM_WP_10668_imvfft))))
head(nTPM_WP_10668_p_imvfft)


nTPM_mean10668<-cbind(nTPM_WP_10668_p_imvfft, nTPM_EV_10668_e_imvfft)
head(nTPM_mean10668)

##############################################################################
##############################################################################
#order the phenodata by
pheno_10936<-subset(merged_meta, Isolate=="KE06")
pheno_10936<-pheno_10936[order(pheno_10936$timepoint),]
pheno_10936
#select the whole parasite data
pheno_10936e<-subset(pheno_10936, condition=="EV")
pheno_10936e

nTPM_EV_10936<-as.data.frame((nTPMdata[,rownames(pheno_10936e)]))
nTPM_EV_10936_mvfft<-mvfft(as.matrix(t(nTPM_EV_10936)))
Heatmap(t(scale(Re(nTPM_EV_10936_mvfft))), show_row_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,column_labels = 1:12)
nTPM_EV_10936_mvfft[-c(1),]=0+0i
nTPM_EV_10936_imvfft<-mvfft(nTPM_EV_10936_mvfft, inverse = TRUE)/12
t(Re(nTPM_EV_10936_imvfft))
dim(na.omit(t(nTPM_EV_10936_imvfft)))

pheno_10936p<-subset(pheno_10936, condition=="WP")
pheno_10936p

############################################################################
nTPM_WP_10936<-as.data.frame((nTPMdata[,rownames(pheno_10936p)]))
##########################################################################
#extract the fourier normalized data
nTPM_EV_10936_e_imvfft<-na.omit(as.data.frame(t(Re(nTPM_EV_10936_imvfft))))
head(nTPM_EV_10936_e_imvfft)

#####################################################################################
#calculate signal to noise ratio
nTPM_WP_10936_mvfft<-mvfft(as.matrix(t(nTPM_WP_10936)))
Heatmap(t(scale(Re(nTPM_WP_10936_mvfft))), show_row_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,column_labels = 1:12)
nTPM_WP_10936_mvfft[-c(1),]=0+0i
nTPM_WP_10936_imvfft<-mvfft(nTPM_WP_10936_mvfft, inverse = TRUE)/12
t(Re(nTPM_WP_10936_imvfft))
dim(na.omit(t(nTPM_WP_10936_imvfft)))

nTPM_WP_10936_p_imvfft<-na.omit(as.data.frame(t(Re(nTPM_WP_10936_imvfft))))
head(nTPM_WP_10936_p_imvfft)


nTPM_mean10936<-cbind(nTPM_WP_10936_p_imvfft, nTPM_EV_10936_e_imvfft)
head(nTPM_mean10936)

#########################################################################
#########################################################################
#order the phenodata by
pheno_Dd2<-subset(merged_meta, Isolate=="Dd2")
pheno_Dd2<-pheno_Dd2[order(pheno_Dd2$timepoint),]
pheno_Dd2
#select the whole parasite data
pheno_Dd2e<-subset(pheno_Dd2, condition=="EV")
pheno_Dd2e

nTPM_EV_Dd2<-as.data.frame((nTPMdata[,rownames(pheno_Dd2e)]))
nTPM_EV_Dd2_mvfft<-mvfft(as.matrix(t(nTPM_EV_Dd2)))
Heatmap(t(scale(Re(nTPM_EV_Dd2_mvfft))), show_row_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,column_labels = 1:11)
nTPM_EV_Dd2_mvfft[-c(1),]=0+0i
nTPM_EV_Dd2_imvfft<-mvfft(nTPM_EV_Dd2_mvfft, inverse = TRUE)/11
t(Re(nTPM_EV_Dd2_imvfft))
dim(na.omit(t(nTPM_EV_Dd2_imvfft)))

pheno_Dd2p<-subset(pheno_Dd2, condition=="WP")
pheno_Dd2p

############################################################################
nTPM_WP_Dd2<-as.data.frame((nTPMdata[,rownames(pheno_Dd2p)]))
##########################################################################
#extract the fourier normalized data
nTPM_EV_Dd2_e_imvfft<-na.omit(as.data.frame(t(Re(nTPM_EV_Dd2_imvfft))))
head(nTPM_EV_Dd2_e_imvfft)

#####################################################################################
#calculate signal to noise ratio
nTPM_WP_Dd2_mvfft<-mvfft(as.matrix(t(nTPM_WP_Dd2)))
Heatmap(t(scale(Re(nTPM_WP_Dd2_mvfft))), show_row_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,column_labels = 1:11)
nTPM_WP_Dd2_mvfft[-c(1),]=0+0i
nTPM_WP_Dd2_imvfft<-mvfft(nTPM_WP_Dd2_mvfft, inverse = TRUE)/11
t(Re(nTPM_WP_Dd2_imvfft))
dim(na.omit(t(nTPM_WP_Dd2_imvfft)))

nTPM_WP_Dd2_p_imvfft<-na.omit(as.data.frame(t(Re(nTPM_WP_Dd2_imvfft))))
head(nTPM_WP_Dd2_p_imvfft)


nTPM_meanDd2<-cbind(nTPM_WP_Dd2_p_imvfft, nTPM_EV_Dd2_e_imvfft)
head(nTPM_meanDd2)

#############################################################################
#######################################################################################
####plothistogram to choose cut-off



frequncy1<-cbind(nTPM_mean9106,nTPM_mean9626, 
                 nTPM_mean9106HB,nTPM_mean10668,
                 nTPM_meanDd2,nTPM_mean10936)
head(frequncy1)

frequncy1<-frequncy1[,colnames(ntpmdata$normalizedCounts)]

all(rownames(frequncy1)==rownames(ntpmdata$normalizedCounts))
all(colnames(frequncy1)==colnames(ntpmdata$normalizedCounts))


nTPM_norm_final_data<-ntpmdata$normalizedCounts+frequncy1
head(nTPM_norm_final_data)





#plot EBA175
plotc((as.data.frame(t(nTPM_norm_final_data[,rownames(Merged_WP)]))$PF3D7_0731500),
      (as.data.frame(t(nTPM_norm_final_data[,rownames(Merged_EV)]))$PF3D7_0731500))


#plot EBA175
plotc((as.data.frame(t(nTPM_norm_final_data[,rownames(Merged_WP)]))$PF3D7_0930300),
      (as.data.frame(t(nTPM_norm_final_data[,rownames(Merged_EV)]))$PF3D7_0930300))



#plot EBA175
plotc((as.data.frame(t(nTPM_norm_final_data[,rownames(Merged_WP)]))$PF3D7_1353200),
      (as.data.frame(t(nTPM_norm_final_data[,rownames(Merged_EV)]))$PF3D7_1353200))

save(nTPM_norm_final_data,file = "~/Documents/phd_after_upgrading_rnaseq/12_timepoints/mapped2pf3D7_version2/Fourier_PCA_RUVr/nTPM_norm_final_data2.Rdata")



##################################################################################
#load("~/Documents/phd_after_upgrading_rnaseq/12_timepoints/mapped2pf3D7_version2/Fourier_PCA_RUVr/nTPM_norm_final_data.Rdata")


nTPMdata_meta<-merge(merged_meta, t(nTPMdata), by="row.names")
head(nTPMdata_meta[1:10])

library(dplyr)
median_nTPM<-nTPMdata_meta[,-c(1,2,4,6,7)] %>% 
  group_by(Isolate, condition) %>%
  summarise_all("mean")
head(median_nTPM)

median_nTPM_melt<-melt(median_nTPM,id.vars = c(1:2))
head(median_nTPM_melt)

ev_median<-subset(median_nTPM_melt, condition=="EV")
wp_median<-subset(median_nTPM_melt, condition=="WP")
head(wp_median)

merge_variable<-merge(ev_median, wp_median, by="variable")
head(merge_variable)
colnames(merge_variable)[1]<-"IDs"

cluster_df<-read.xlsx("Fourier_PCA_RUVr/Tables/Delta_mesor_withMaraLawnizack_withTranscription_DecayRates.xlsx")
head(cluster_df)
colnames(cluster_df)

merge_variable_cluster<-merge(merge_variable, cluster_df[,c(1,11), drop=FALSE], by="IDs")
head(merge_variable_cluster)


graphics.off()


jpeg(filename = "Fourier_PCA_RUVr/time_independent_analysis/Total_expression_by_cluster.jpeg",
     units = "in",width = 20,height = 4.5,res = 300)
ggplot(merge_variable_cluster, aes(value.x, value.y, color=Legend, shape=Legend))+
  geom_point()+
  facet_grid(cols = vars(Isolate.x))+
  theme_pubr(border = TRUE,legend = "right")+
  theme(text = element_text(size = 20))+
  ylab("Mean expression in WP")+
  xlab("Mean secretion in PfEVs")+
  labs(colour = "Cluster", shape="Cluster")+
  scale_x_continuous(breaks = c(3,6,9,12))+
  guides(colour = guide_legend(override.aes = list(size=4)))
dev.off()


ggplot(subset(merge_mesor,variable.x=="deltamesor_KE01"), aes(value.x, value.y))+
  geom_point()
facet_grid(cols = vars(variable.x))



merge_variable_cluster<-merge(merge_mesor, cluster_df[,c(1,11,12), drop=FALSE], by="IDs")
head(merge_variable_cluster)


merge_variable_cluster$Isolate.x<-factor(as.factor(merge_variable_cluster$Isolate.x),
                                         levels = c("KE01","sKE01","KE02" ,"KE04","KE06","Dd2"))
head(merge_variable_cluster)

############################################################################################

jpeg(filename = "Fourier_PCA_RUVr/time_independent_analysis/Total_expression_by_cluster.jpeg",
     units = "in",width = 6,height = 5,res = 600)
ggplot(subset(merge_variable_cluster,variable.x=="deltamesor_KE01"), aes(value.x, value.y, color=Legend, shape=Legend))+
  geom_point()+
  theme_pubr(border = TRUE,legend = "right")+
  theme(text = element_text(size = 20))+
  ylab("Mean expression in WP")+
  xlab("Mean secretion in PfEVs")+
  labs(colour = "Cluster", shape="Cluster")+
  scale_x_continuous(breaks = c(3,6,9,12))+
  guides(colour = guide_legend(override.aes = list(size=4)))
dev.off()

jpeg(filename = "Fourier_PCA_RUVr/time_independent_analysis/MA_plots.jpeg",
     units = "in",width = 10,height = 5,res = 600)
p1<-ggplot(subset(merge_variable_cluster,variable.x=="deltamesor_KE01"), aes(value.x, median_logFC, 
                                                                             color=Legend, shape=Legend))+
  geom_point()+
  theme_pubr(border = TRUE,legend = "right")+
  theme(text = element_text(size = 20))+
  xlab("Mean secretion in PfEVs")+
  ylab("")+
  labs(colour = "Cluster", shape="Cluster")+
  #scale_x_continuous(breaks = c(0,3,6,9,12))+
  guides(colour = guide_legend(override.aes = list(size=4)))
p2<-ggplot(subset(merge_variable_cluster,variable.x=="deltamesor_KE01"), aes(value.y, median_logFC, 
                                                                             color=Legend, shape=Legend))+
  geom_point()+
  theme_pubr(border = TRUE,legend = "right")+
  theme(text = element_text(size = 20))+
  xlab("Mean expression in WP")+
  ylab("Delta mesor")+
  labs(colour = "Cluster", shape="Cluster")+
  #scale_x_continuous(breaks = c(0,3,6,9,12))+
  guides(colour = guide_legend(override.aes = list(size=4)))+NoLegend()

p2+p1

dev.off()
































