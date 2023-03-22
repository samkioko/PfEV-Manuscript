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
library(tximport)
library(RUVSeq)

##################################################################
#load the phenodata
pheno_9106<-read.csv("~/Documents/phd_after_upgrading_rnaseq/12_timepoints/9106/9106_12_timepoint_metadata.csv",
                     row.names = 2, 
                     check.names = FALSE)

#load the normalized data
load("~/Documents/phd_after_upgrading_rnaseq/12_timepoints/mapped2pf3D7_version2/Fourier_PCA_RUVr/nTPM_norm_final_data2.Rdata")
pheno_9106<-pheno_9106[order(pheno_9106$timepoint),]
pheno_9106

################################################################

TPM_9106<-nTPM_norm_final_data[,rownames(pheno_9106)]

truehist(as.matrix(TPM_9106))


#select the whole parasite data
pheno_p<-subset(pheno_9106, condition=="WP")
pheno_p
tpm_9106_p<-TPM_9106[,rownames(pheno_p)]
head(tpm_9106_p)
tpm_9106_p_t<-as.data.frame(t(tpm_9106_p))
tpm_9106_p_t

##LOESS the whole data
vars <- colnames(tpm_9106_p_t)
vars
## covariate
id <- c(1:nrow(tpm_9106_p_t))
id


#id<-pheno_p$time_numeric
## define a loess filter function (fitting loess regression line)
loess.filter <- function (x, span) loess(formula = paste(x, "id", 
                                                         sep = "~"),
                                         data = tpm_9106_p_t,
                                         degree = 2,
                                         span = span)$fitted 
## apply filter column-by-column
new.dat_p <- as.data.frame(lapply(vars, loess.filter, span = 0.75),
                           col.names = colnames(tpm_9106_p_t))

log2_loess_p<-t(new.dat_p)
head(log2_loess_p)
colnames(log2_loess_p)<-colnames(tpm_9106_p)
head(tpm_9106_p)
head(log2_loess_p)

#fourier transform RNA-seq data
tpm_9106_p_t<-as.data.frame(t(log2_loess_p))
tpm_9106_p_t[1:10]

#5/12
#######################################################################################
#calculate signal to noise ratio
tpm_9106_p_t_mvfft<-mvfft(as.matrix(tpm_9106_p_t))
#Heatmap(t(scale(Re(tpm_9106_p_t_mvfft))), show_row_names = FALSE,cluster_columns = FALSE)
tpm_9106_p_t_mvfft[-c(1:3,11),]=0+0i
tpm_9106_p_t_imvfft<-mvfft(tpm_9106_p_t_mvfft, inverse = TRUE)/11
t(Re(tpm_9106_p_t_imvfft))

##########################################################################
#extract the fourier normalized data
tpm_9106_p_imvfft<-as.data.frame(t(Re(tpm_9106_p_t_imvfft)))
head(tpm_9106_p_imvfft)

############################################################################
library(matrixTests)
dim(tpm_9106_p_imvfft)
cosinor_p<-row_cosinor(x = (tpm_9106_p_imvfft),
                       t = c(4,8,12,16,20,24,28,32,36,40,44),
                       period = 44)
head(cosinor_p)
dim(cosinor_p)
table(cosinor_p$pvalue<0.05)
cosinor_p$acrophase<-(cosinor_p$acrophase/44)*48
head(cosinor_p)

############################################################
noise_p<-as.data.frame(t(as.data.frame(lapply(tpm_9106_p_t,function(x) sd(x)))))
head(noise_p)

signal_p<-as.data.frame(t(as.data.frame(lapply(as.data.frame(Re(tpm_9106_p_t_imvfft)),
                                               function(x) max(x)))))
head(signal_p)

merge_SNR_p<-merge(signal_p,noise_p,by="row.names")
head(merge_SNR_p)
names(merge_SNR_p)[1:3]<-c("IDs","signal","noise")
merge_SNR_p$SNR<-sqrt(merge_SNR_p$signal^2/merge_SNR_p$noise^2)
head(merge_SNR_p)

ggplot(merge_SNR_p, aes(x=SNR))+geom_histogram()

table(merge_SNR_p$SNR>5)


#############################################################
#extract genes with signal to noise ratio of greater than 10
ars_p_meta<-cosinor_p
ars_p_meta$IDs<-rownames(ars_p_meta)
head(ars_p_meta)

#merge with signal to noise ratio
ars_p_meta<-merge(ars_p_meta,merge_SNR_p, by="IDs")
ars_p_meta

ars_p_meta_sig<-subset(ars_p_meta, ars_p_meta$SNR>5)
dim(ars_p_meta_sig)

ars_p_meta_sig<-ars_p_meta_sig[order(ars_p_meta_sig$acrophase),]

lgr_sig_p<-(tpm_9106_p_imvfft)[ars_p_meta_sig$IDs,]
dim(lgr_sig_p)
head(lgr_sig_p)
###########################################################
library(circlize)
col_funp = colorRamp2(c(-2, 0, 2), c("green", "black", "red"))

#plot 9106 parasite heatmap
hab<-HeatmapAnnotation(time = anno_barplot(pheno_p$timepoint,
                                           gp = gpar(fontsize = 20,
                                                     fill=2),
                                           bar_width = 1),
                       show_annotation_name = TRUE,
                       annotation_name_side = "left",
                       annotation_name_rot = 90,
                       annotation_name_gp = gpar(fontsize=20),
                       height = unit(0.8,"cm"))

dim(lgr_sig_p)
#graphics.off()
jpeg(filename = "TU_cqn_loess_fft/9106/parasite_9106.jpeg",
     units = "in",width = 1.7,height = 10,res = 250)
htm_9106_p<-Heatmap(t(scale(t(lgr_sig_p[ars_p_meta_sig$IDs,]))),
                   clustering_method_rows = "average",
                   clustering_distance_rows = "canberra",
                   heatmap_legend_param = list(title = "mean z score", 
                                               direction="vertical"),
                   top_annotation = hab,
                   show_row_names = FALSE,
                   col = col_funp,cluster_rows = FALSE,
                   show_column_names = FALSE,
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   show_row_dend = FALSE, 
                   cluster_row_slices = FALSE,
                   row_gap = unit(0, "mm"),
                   row_names_gp = gpar(fontsize = 8),
                   row_names_max_width = unit(15,"cm"),
                   column_title = "KE01",
                   row_title="5233 genes",
                   column_title_gp = gpar(fontsize = 20),
                   row_title_gp = gpar(fontsize = 20),
                   row_title_rot = 90, use_raster = FALSE)
htm_9106_p

#+rowAnnotation(mark = anno)
htm_9106_pj<-draw(htm_9106_p, merge_legend = TRUE, 
                 heatmap_legend_side = "right", 
                 annotation_legend_side = "top",
                 show_heatmap_legend = FALSE)
dev.off()


######################################################
pheno_9106<-pheno_9106[order(pheno_9106$timepoint),]
pheno_9106
#select the whole parasite data
pheno_e<-subset(pheno_9106, condition=="EV")
pheno_e
tpm_9106_e<-RPKM_9106[,rownames(pheno_e)]
head(tpm_9106_e)
head(tpm_9106_e)
tpm_9106_e_t<-as.data.frame(t(tpm_9106_e))
tpm_9106_e_t

dim(tpm_9106_e_t)

##LOESS the whole data
vars <- colnames(tpm_9106_e_t)
vars
## covariate
id <- c(1:nrow(tpm_9106_e_t))
id
#id<-pheno_p$time_numeric
## define a loess filter function (fitting loess regression line)
loess.filter <- function (x, span) loess(formula = paste(x, "id", 
                                                         sep = "~"),
                                         data = tpm_9106_e_t,
                                         degree = 2,
                                         span = span)$fitted 
## apply filter column-by-column
new.dat_e_exvivo <- as.data.frame(lapply(vars, loess.filter, 
                                         span = 0.75),
                                  col.names = colnames(tpm_9106_e_t))

log2_loess_e<-t(new.dat_e_exvivo)
head(log2_loess_e)
colnames(log2_loess_e)<-colnames(tpm_9106_e)
head(log2_loess_e)
dim(log2_loess_e)

tpm_9106_e_t<-as.data.frame(t(log2_loess_e))
tpm_9106_e_t[1:10]


#######################################################################################
tpm_9106_e_t_mvfft<-mvfft(as.matrix(tpm_9106_e_t))
tpm_9106_e_t_mvfft[c(4:10),]=0+0i
tpm_9106_e_t_imvfft<-mvfft(tpm_9106_e_t_mvfft, inverse = TRUE)/11
t(Re(tpm_9106_e_t_imvfft))
############
#extract the fourier normalized data
tpm_9106_e_imvfft<-as.data.frame(t(Re(tpm_9106_e_t_imvfft)))
head(tpm_9106_e_imvfft)


#extract peak signal to noise ratio (PSNR)
##########################################################################
dim(tpm_9106_e_imvfft)

library(matrixTests)
cosinor_e<-row_cosinor(x = tpm_9106_e_imvfft,
                       t = c(4,8,12,16,20,24,28,32,36,40,44),
                       period = 44)
head(cosinor_e)
table(cosinor_e$pvalue<0.05)
cosinor_e$acrophase<-(cosinor_e$acrophase/44)*48
#extract peak signal to noise ratio (PSNR)
noise_e<-as.data.frame(t(as.data.frame(lapply(tpm_9106_e_t,function(x) sd(x)))))
head(noise_e)
signal_e<-as.data.frame(t(as.data.frame(lapply(as.data.frame((Re(tpm_9106_e_t_imvfft))),function(x) max(x)))))
head(signal_e)


merge_SNR_e<-merge(signal_e,noise_e,by="row.names")
head(merge_SNR_e)
names(merge_SNR_e)[1:3]<-c("IDs","signal","noise")
merge_SNR_e$SNR<-sqrt(merge_SNR_e$signal^2/merge_SNR_e$noise^2)
head(merge_SNR_e)

ggplot(merge_SNR_e, aes(x=SNR))+geom_histogram()

table(merge_SNR_e$SNR>5)#4785 genes meet the cut-off of signal to noise ratio


#extract genes with signal to noise ratio of greater than 10
ars_e_meta<-cosinor_e
ars_e_meta$IDs<-rownames(ars_e_meta)
head(ars_e_meta)

#merge with signal to noise ratio
ars_e_meta<-merge(ars_e_meta,merge_SNR_e, by="IDs")
ars_e_meta

ars_e_meta_sig<-subset(ars_e_meta, ars_e_meta$SNR>5)
dim(ars_e_meta_sig)

ars_e_meta_sig<-ars_e_meta_sig[order(ars_e_meta_sig$acrophase),]

lgr_sig_e<-tpm_9106_e_imvfft[ars_e_meta_sig$IDs,]
dim(lgr_sig_e)
head(lgr_sig_e)



library(circlize)
col_fune = colorRamp2(c(-2,0,2), c("dodgerblue", "black", "red"))


#plot 9106 parasite heatmap
hab<-HeatmapAnnotation(time = anno_barplot(pheno_e$timepoint,
                                           gp = gpar(fontsize = 20,
                                                     fill=2),
                                           bar_width = 1),
                       show_annotation_name = TRUE,
                       annotation_name_side = "left",
                       annotation_name_rot = 90,
                       annotation_name_gp = gpar(fontsize=20),
                       height = unit(0.8,"cm"))

#graphics.off()
jpeg(filename = "TU_cqn_loess_fft/9106/ev_9106.jpeg",
     units = "in",width = 1.7,height = 10,res = 250)
htm_9106_e<-Heatmap(t(scale(t(lgr_sig_e[ars_e_meta_sig$IDs,]))),
                   clustering_method_rows = "average",
                   clustering_distance_rows = "canberra",
                   heatmap_legend_param = list(title = "log2 exp", 
                                               direction="vertical"),
                   top_annotation = hab,
                   show_row_names = FALSE,
                   col = col_fune,
                   cluster_rows = FALSE,
                   show_column_names = FALSE,
                   show_column_dend = FALSE,
                   cluster_columns = FALSE,
                   show_row_dend = FALSE, 
                   cluster_row_slices = FALSE,
                   row_gap = unit(0, "mm"),
                   row_names_gp = gpar(fontsize = 8),
                   row_names_max_width = unit(15,"cm"),
                   column_title = "KE01",
                   row_title="5245 genes",
                   column_title_gp = gpar(fontsize = 20),
                   row_title_gp = gpar(fontsize = 20),
                   row_title_rot = 90, use_raster = FALSE)
htm_9106_e

#+rowAnnotation(mark = anno)
htm_9106_ej<-draw(htm_9106_e, merge_legend = TRUE, 
                 heatmap_legend_side = "right", 
                 annotation_legend_side = "top",
                 show_heatmap_legend = FALSE)
dev.off()

###################################################################
library(itsmr)
truehist(as.matrix(tpm_9106_p_imvfft))
truehist(as.matrix(tpm_9106_e_imvfft))
#plot EBA175
plotc((as.data.frame(t(tpm_9106_p_imvfft))$PF3D7_0731500),
      (as.data.frame(t(tpm_9106_e_imvfft))$PF3D7_0731500))

#msp2
plotc((as.data.frame(t(tpm_9106_p_imvfft))$PF3D7_0206800),
      (as.data.frame(t(tpm_9106_e_imvfft))$PF3D7_0206800))

#
plotc((as.data.frame(t(tpm_9106_p_imvfft))$PF3D7_1353200),
      (as.data.frame(t(tpm_9106_e_imvfft))$PF3D7_1353200))

#msp1
plotc((as.data.frame(t(tpm_9106_p_imvfft))$PF3D7_0930300),
      (as.data.frame(t(tpm_9106_e_imvfft))$PF3D7_0930300))

#rh5
plotc((as.data.frame(t(tpm_9106_p_imvfft))$PF3D7_0424100),
      (as.data.frame(t(tpm_9106_e_imvfft))$PF3D7_0424100))

#load annotation table
annot<-read.xlsx("~/Documents/phd_after_upgrading_rnaseq/12_timepoints/mapped2pf3D7/Pf3D7_annotation_file.xlsx")
head(annot)
colnames(annot)[1]<-"IDs"
annot<-annot[!duplicated(annot$IDs),]
annot

#find in phase genes <less than 12 hours
merge_data<-merge(ars_p_meta, ars_e_meta, by="IDs")
head(merge_data)
dim(merge_data)
names(merge_data)[1]<-"IDs"
merge_data_annot<-merge(annot,merge_data,by="IDs",all.y=TRUE)
head(merge_data_annot)
dim(merge_data_annot)

colnames(merge_data_annot)<-gsub(pattern = "\\.x",
                                 replacement = "_KE01_Parasite",
                                 colnames(merge_data_annot))

colnames(merge_data_annot)<-gsub(pattern = "\\.y",
                                 replacement = "_KE01_EVs",
                                 colnames(merge_data_annot))

head(merge_data_annot)
phase1<-abs(merge_data_annot$acrophase_KE01_Parasite - merge_data_annot$acrophase_KE01_EVs)
phase1
phase2<-48-abs(merge_data_annot$acrophase_KE01_Parasite - merge_data_annot$acrophase_KE01_EVs)
phase2

merge_data_annot$phase_diff<-pmin(phase1,phase2)
merge_data_annot$phase_diff
head(merge_data_annot)
dim(merge_data_annot)

#############################

#colnames(merge_data_annot)
merge_data_annot$pmin_SNR<-pmin(merge_data_annot$SNR_KE01_Parasite,
                                merge_data_annot$SNR_KE01_EVs)


#correlate EVs and parasites
pf<-as.data.frame(t(tpm_9106_p_imvfft))
pf[,1:4]

ev<-as.data.frame(t(tpm_9106_e_imvfft[rownames(tpm_9106_p_imvfft),]))
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

#merge with annotation
merge_data_annot_ortho_cor<-merge(res1, merge_data_annot, by="IDs")
head(merge_data_annot_ortho_cor)
dim(merge_data_annot_ortho_cor)

#########################plothistogram to choose cut-off
ggplot(merge_data_annot_ortho_cor, aes(x=rho))+
  geom_histogram(color="darkblue", fill="lightblue")

############################################################################
log2_loess_scale_p<-as.data.frame(t(scale(t(tpm_9106_p_imvfft))))
head(log2_loess_scale_p)
log2_loess_scale_e<-as.data.frame(t(scale(t(tpm_9106_e_imvfft))))
head(log2_loess_scale_e)

pf9106_cbind_scale<-cbind(log2_loess_scale_e,log2_loess_scale_p)
head(merge_data_annot_ortho_cor)

library(plotmm)
#find the cut-off of phase and out-phase using gausian mixture models
mixmdl <- mixtools::normalmixEM(merge_data_annot_ortho_cor$phase_diff, k = 2)
mixmdl
plot_cut_point(mixmdl, plot = TRUE, color = "amerika") 
plot_cut_point(mixmdl, plot = FALSE, color = "amerika") # cut-off is 12.2 hours

# produces plot
#visualize ggplot
plot_mm(mixmdl, 2) +
  ggplot2::labs(title = "gmm model of cut-off of phase-shift")


#write.xlsx(merge_data_annot_ortho_cor,
#file = "TU_cqn_loess_fft/9106/pf9106_merge_data_annot_ortho_cor.xlsx",
#overwrite = TRUE)

merge_data_annot_ortho_cor$trough_EV<-merge_data_annot_ortho_cor$acrophase_KE01_EVs-24
index<-merge_data_annot_ortho_cor$trough_EV<0
index
merge_data_annot_ortho_cor$trough_EV[index]<-(merge_data_annot_ortho_cor$trough_EV[index])+48


merge_phase<-subset(merge_data_annot_ortho_cor, 
                    phase_diff<4 &
                      merge_data_annot_ortho_cor$pmin_SNR>5)
dim(merge_phase)
table(merge_data_annot_ortho_cor$rho>0)


merge_phase$pmax<-pmax(merge_phase$acrophase_KE01_Parasite,
                       merge_phase$acrophase_KE01_EVs)

merge_phase<-merge_phase[order(merge_phase$pmax,
                               decreasing =FALSE),]
data_phase<-pf9106_cbind_scale[merge_phase$IDs,]
head(data_phase)
dim(data_phase)


library(circlize)
col_funev = colorRamp2(c(-2, 0, 2), c("blue", "black", "gold"))

habar<-HeatmapAnnotation(time = anno_barplot(pheno_9106$timepoint,
                                             gp = gpar(fontsize = 20,fill=2),
                                             bar_width = 1),
                         show_annotation_name = TRUE,
                         annotation_name_rot = 90,
                         annotation_name_side = "left",
                         height = unit(0.8,"cm"),
                         annotation_name_gp = gpar(fontsize=20))
habar
dim(data_phase)
head(data_phase)
graphics.off()
#jpeg(filename = "TU_cqn_loess_fft/9106/circabidian_inphase_9106_sig.jpeg",
# units = "in",width = 2.7,height = 7,res = 320)
htm_9106_phasic<-Heatmap(as.matrix(data_phase[,rownames(pheno_9106)]),
                        clustering_method_rows = "average",
                        clustering_distance_rows = "canberra",
                        heatmap_legend_param = list(title = "log2 exp", 
                                                    direction="vertical"),
                        top_annotation = habar,
                        show_row_names = FALSE,
                        col = col_funev,
                        cluster_rows = FALSE,
                        show_column_names = FALSE,
                        row_names_side = "right", 
                        show_column_dend = FALSE,
                        cluster_columns = FALSE,
                        row_title = "652 in phase genes",
                        show_row_dend = FALSE, 
                        column_split = factor(pheno_9106$condition,
                                              levels = c("WP","EV")), 
                        cluster_row_slices = FALSE,
                        cluster_column_slices = FALSE,
                        row_gap = unit(0, "mm"),use_raster = FALSE,
                        row_names_gp = gpar(fontsize = 12,col="black"),
                        row_names_max_width = unit(5,"cm"),
                        column_title_gp = gpar(fontsize = 20,col="black"),
                        row_title_gp = gpar(fontsize = 20,col="black"))

htm_9106_phasic
#graphics.off()
#+rowAnnotation(mark = anno)
htm_9106_phasicj<-draw(htm_9106_phasic, merge_legend = TRUE, 
                      heatmap_legend_side = "right", 
                      annotation_legend_side = "bottom",
                      show_heatmap_legend = FALSE,
                      column_title = "9106", 
                      column_title_gp = gpar(fontsize = 20, col="black"))

#dev.off()

merge_phase
#######################################################
#choose out of phase genes and plot
merge_ophase<-subset(merge_data_annot_ortho_cor,
                     merge_data_annot_ortho_cor$pmin_SNR>0&
                       phase_diff>4)
dim(merge_ophase)#4202

#data_outphase<-read.xlsx("TU_cqn_loess_fft/9106/pf9106_9106_cosinor_results_excelSheet.xlsx",
#sheet = 5)
merge_ophase$pmin<-pmin(merge_ophase$acrophase_KE01_Parasite, 
                        merge_ophase$trough_EV)

merge_ophase<-merge_ophase[order(merge_ophase$pmin,
                                 decreasing =FALSE),]
merge_ophase
data_ophase<-pf9106_cbind_scale[merge_ophase$IDs,]
head(data_ophase)
dim(data_ophase)

library(circlize)
col_funev = colorRamp2(c(-2, 0, 2), c("green","black", "purple"))
#
dim(data_ophase)
graphics.off()

jpeg(filename = "TU_cqn_loess_fft/9106/circabidian_outphase_9106_sigV3.jpeg",
     units = "in",width = 3.5,height = 11,res = 330)
htm_9106_phasic<-Heatmap(as.matrix(data_ophase[,rownames(pheno_9106)]),
                        clustering_method_rows = "average",
                        clustering_distance_rows = "canberra",
                        heatmap_legend_param = list(title = "log2 exp", 
                                                    direction="vertical"),
                        top_annotation = habar,
                        show_row_names = FALSE,
                        col = col_funev,
                        cluster_rows = FALSE,
                        show_column_names = FALSE,
                        row_names_side = "right", 
                        show_column_dend = FALSE,
                        cluster_columns = FALSE,
                        row_title = "4543 out of phase genes",
                        show_row_dend = FALSE, 
                        column_split = factor(pheno_9106$condition,
                                              levels = c("WP","EV")), 
                        cluster_row_slices = FALSE,
                        cluster_column_slices = FALSE,
                        row_gap = unit(0, "mm"),use_raster = FALSE,
                        row_names_gp = gpar(fontsize = 12,col="black"),
                        row_names_max_width = unit(5,"cm"),
                        column_title_gp = gpar(fontsize = 20,col="black"),
                        row_title_gp = gpar(fontsize = 20,col="black"))
htm_9106_phasic

#+rowAnnotation(mark = anno)
htm_9106_phasicj<-draw(htm_9106_phasic, merge_legend = TRUE, 
                      heatmap_legend_side = "right", 
                      annotation_legend_side = "bottom",
                      show_heatmap_legend = FALSE,
                      column_title = "KE01", 
                      column_title_gp = gpar(fontsize = 20, col="black"))

dev.off()
#?pmin
merge_phase
library(ENIGMA)

#plot RH1
plotc(scale(as.data.frame(t(tpm_9106_p_imvfft))$PF3D7_0424200),
      scale(as.data.frame(t(tpm_9106_e_imvfft))$PF3D7_0424200))

#plotEBA165
plotc((as.data.frame(t(tpm_9106_p_imvfft))$PF3D7_0424300),
      (as.data.frame(t(tpm_9106_e_imvfft))$PF3D7_0424300))


#plot EBA175
plotc(scale(as.data.frame(t(tpm_9106_p_imvfft))$PF3D7_0731500),
      scale(as.data.frame(t(tpm_9106_e_imvfft))$PF3D7_0731500))


#KAHRP
plotc(scale(as.data.frame(t(tpm_9106_p_imvfft))$PF3D7_0202000),
      scale(as.data.frame(t(tpm_9106_e_imvfft))$PF3D7_0202000))

#CSP
plotc(scale(as.data.frame(t(tpm_9106_p_imvfft))$PF3D7_0304600),
      scale(as.data.frame(t(tpm_9106_e_imvfft))$PF3D7_0304600))

#plot RESA
plotc((as.data.frame(t(tpm_9106_p_imvfft))$PF3D7_0102200),
      (as.data.frame(t(tpm_9106_e_imvfft))$PF3D7_0102200))

#ferrochelatase
plotc(scale(as.data.frame(t(tpm_9106_p_imvfft))$PF3D7_1364900),
      scale(as.data.frame(t(tpm_9106_e_imvfft))$PF3D7_1364900))


########################################################################
wb <- createWorkbook()
addWorksheet(wb, "cyclic_WP_all")
addWorksheet(wb, "cyclic_EVs_all")
addWorksheet(wb, "Cyclic_both_WP_and_EVs")
#addWorksheet(wb, "cyclic_inphase")
addWorksheet(wb, "cyclic_outphase")
addWorksheet(wb, "log2CPM_loess_WP_fft")
addWorksheet(wb, "log2CPM_loess_EV_fft")
addWorksheet(wb, "expressed_matrix_scaled")
addWorksheet(wb, "phenodata")

writeData(wb, "cyclic_WP_all", ars_p_meta, rowNames = TRUE)
writeData(wb, "cyclic_EVs_all", ars_e_meta, rowNames = TRUE)
writeData(wb, "Cyclic_both_WP_and_EVs", merge_data_annot_ortho_cor, rowNames = TRUE)
#writeData(wb, "cyclic_inphase", merge_phase, rowNames = TRUE)
writeData(wb, "cyclic_outphase", merge_ophase, rowNames = TRUE)
writeData(wb, "log2CPM_loess_WP_fft", tpm_9106_p_imvfft, rowNames = TRUE)
writeData(wb, "log2CPM_loess_EV_fft", tpm_9106_e_imvfft, rowNames = TRUE)
writeData(wb, "expressed_matrix_scaled", pf9106_cbind_scale, rowNames = TRUE)
writeData(wb, "phenodata", pheno_9106, rowNames = TRUE)

#export the data
saveWorkbook(wb, file = "TU_cqn_loess_fft/9106/pfKE01_KE01_cosinor_results_excelSheet.xlsx", 
             overwrite = TRUE)
###############################################################################
#run circacompare
#######################################################################
pf9106_cbind_scale<-as.data.frame(t(cbind(tpm_9106_p_imvfft,
                                           tpm_9106_e_imvfft[rownames(tpm_9106_p_imvfft),])))
pf9106_cbind_scale[1:30]

merge_with_phenodata<-merge(pheno_9106,pf9106_cbind_scale, by="row.names")
head(merge_with_phenodata[,1:10])
merge_with_phenodata$timepoint<-merge_with_phenodata$timepoint
head(merge_with_phenodata)[1:6]
dim(merge_with_phenodata)
library(circacompare)

#######run the circacompare algorithm
merge_with_phenodata$condition<-relevel(as.factor(merge_with_phenodata$condition),
                                        ref = "WP")

fit_un <- lapply(colnames(merge_with_phenodata[,c(5:5243)]), function(x) circacompare(x=merge_with_phenodata,period = 48,alpha_threshold = 1,
                                                                                      col_time="timepoint", col_group="condition", 
                                                                                      col_outcome= x,control = list(period_min=24,period_max=48))$summary)
names(fit_un)<-colnames(merge_with_phenodata[,c(5:5243)])
fit_df<-do.call(rbind.data.frame, fit_un)

fit_df$name<-rownames(fit_df)
head(fit_df)
fit_df$name<-gsub(pattern ="\\..*","" ,fit_df$name)

fit_df_wide<-pivot_wider(data = fit_df,names_from = "parameter",
                         id_cols = "name",values_from = "value")

head(fit_df_wide)
names(fit_df_wide)[1]<-"IDs"

library(metaseqR2)
#load annotation table
annot<-read.xlsx("../mapped2pf3D7/Pf3D7_annotation_file.xlsx")
head(annot)
colnames(annot)[1]<-"IDs"
annot<-annot[!duplicated(annot$IDs),]
annot

fit_df_wide$pvalue_meta_rhythmic<-apply(as.matrix(fit_df_wide[,2:3]),
                                        1,combineSimes)
fit_df_wide<-merge(annot,fit_df_wide, by="IDs")
dim(fit_df_wide)
head(fit_df_wide)
table(fit_df_wide$`P-value for mesor difference`<0.05)

#library(openxlsx)
#fit_df_wide<-read.xlsx("TU_cqn_loess_fft/9626/Pf9626_KE02_Circacompare.xlsx",
                       #colNames = TRUE,check.names = FALSE)

circa_WP<-subset(fit_df_wide,
                 fit_df_wide$`Presence of rhythmicity (p-value) for WP`<0.05)
circa_WP
######################################
circa_WP<-circa_WP[order(circa_WP$`WP peak time hours`),]
dim(circa_WP)

lgr_sig_p<-(tpm_9106_p_imvfft)[circa_WP$IDs,]
dim(lgr_sig_p)
head(lgr_sig_p)
###########################################################
library(circlize)
col_funp = colorRamp2(c(-2, 0, 2), c("dodgerblue", "black", "red"))

#plot 9106 parasite heatmap
hab<-HeatmapAnnotation(time = anno_barplot(pheno_p$timepoint,
                                           gp = gpar(fontsize = 35,
                                                     fill=2),
                                           bar_width = 1),
                       show_annotation_name = FALSE,
                       annotation_name_side = "left",
                       annotation_name_rot = 90,
                       annotation_name_gp = gpar(fontsize=35),
                       height = unit(1.5,"cm"))

library(circlize)
col_fun = colorRamp2(breaks = c(0,48), colors = c("grey99","grey0"))

column_ha = HeatmapAnnotation(time=pheno_p$timepoint,col = list(time=col_fun),
                              show_annotation_name = FALSE)


dim(lgr_sig_p)
#graphics.off()
jpeg(filename = "Fourier_PCA_RUVr2/9106/parasite_9106V2.jpeg",
     units = "in",width = 2,height = 8,res = 500)
htm_9106_p<-Heatmap(t(scale(t(lgr_sig_p[circa_WP$IDs,]))),
                     clustering_method_rows = "average",
                     clustering_distance_rows = "canberra",
                     heatmap_legend_param = list(title = "mean z score", 
                                                 direction="vertical"),
                     top_annotation = column_ha,
                     show_row_names = FALSE,
                     col = col_funp,cluster_rows = FALSE,
                     show_column_names = FALSE,
                     show_column_dend = FALSE,
                     cluster_columns = FALSE,
                     show_row_dend = FALSE, 
                     cluster_row_slices = FALSE,
                     row_gap = unit(0, "mm"),
                     row_names_gp = gpar(fontsize = 8),
                     row_names_max_width = unit(15,"cm"),
                     column_title = "KE01",
                     row_title="5098 genes",
                     column_title_gp = gpar(fontsize = 35),
                     row_title_gp = gpar(fontsize = 35),
                     row_title_rot = 90, use_raster = FALSE)
htm_9106_p

#+rowAnnotation(mark = anno)
htm_9106_pj<-draw(htm_9106_p, merge_legend = TRUE, 
                   heatmap_legend_side = "right", 
                   annotation_legend_side = "top",
                   show_heatmap_legend = FALSE)
dev.off()

##########################################################
###############################################################
##############################################################
circa_EV<-subset(fit_df_wide,
                 fit_df_wide$`Presence of rhythmicity (p-value) for EV`<0.05)

######################################
circa_EV<-circa_EV[order(circa_EV$`EV peak time hours`),]
dim(circa_EV)

lgr_sig_e<-(tpm_9106_e_imvfft)[circa_EV$IDs,]
dim(lgr_sig_e)
head(lgr_sig_e)
###########################################################
library(circlize)
col_fune = colorRamp2(c(-2, 0, 2), c("dodgerblue", "black", "gold"))

#plot 9106 parasite heatmap
hab<-HeatmapAnnotation(time = anno_barplot(pheno_e$timepoint,
                                           gp = gpar(fontsize = 20,
                                                     fill=2),
                                           bar_width = 1),
                       show_annotation_name = TRUE,
                       annotation_name_side = "left",
                       annotation_name_rot = 90,
                       annotation_name_gp = gpar(fontsize=20),
                       height = unit(0.8,"cm"))

dim(lgr_sig_e)
#graphics.off()
jpeg(filename = "Fourier_PCA_RUVr2/9106/ev_9106V2.jpeg",
     units = "in",width = 2,height = 8,res = 350)
htm_9106_e<-Heatmap(t(scale(t(lgr_sig_e[circa_EV$IDs,]))),
                     clustering_method_rows = "average",
                     clustering_distance_rows = "canberra",
                     heatmap_legend_param = list(title = "RPKM z score", 
                                                 direction="vertical"),
                     top_annotation = column_ha,
                     show_row_names = FALSE,
                     col = col_fune,cluster_rows = FALSE,
                     show_column_names = FALSE,
                     show_column_dend = FALSE,
                     cluster_columns = FALSE,
                     show_row_dend = FALSE, 
                     cluster_row_slices = FALSE,
                     row_gap = unit(0, "mm"),
                     row_names_gp = gpar(fontsize = 8),
                     row_names_max_width = unit(15,"cm"),
                     column_title = "KE01",
                     row_title="5048 genes",
                     column_title_gp = gpar(fontsize = 35),
                     row_title_gp = gpar(fontsize = 35),
                     row_title_rot = 90, use_raster = FALSE)
htm_9106_e

#+rowAnnotation(mark = anno)
htm_9106_ej<-draw(htm_9106_e, merge_legend = TRUE, 
                   heatmap_legend_side = "right", 
                   annotation_legend_side = "top",
                   show_heatmap_legend = FALSE)
dev.off()

##############################################################
fit_df_wide$trough_EV<-fit_df_wide$`EV peak time hours`-24
index<-fit_df_wide$trough_EV<0
index
fit_df_wide$trough_EV[index]<-(fit_df_wide$trough_EV[index])+48

#################################################
###############################################
#choose out of phase genes and plot
merge_ophase<-subset(fit_df_wide,
                     fit_df_wide$`Presence of rhythmicity (p-value) for EV`<0.05 &
                      fit_df_wide$`Presence of rhythmicity (p-value) for WP`<0.05 &
                       #abs(fit_df_wide$`Phase difference estimate`)>12 &
                       fit_df_wide$`P-value for difference in phase`<0.05
                     )
dim(merge_ophase)#
head(merge_ophase)

library(ggplot2)
library(viridis)

ggplot(fit_df_wide, aes(`EV peak time hours`, 
                         `WP peak time hours`))+
  geom_point(adjust=4)+
  scale_color_viridis(option="B")


merge_ophase$pmin<-pmin(merge_ophase$`WP peak time hours`, 
                        merge_ophase$trough_EV)

merge_ophase<-merge_ophase[order(merge_ophase$pmin,
                                 decreasing =FALSE),]
dim(merge_ophase)
log2_loess_scale_p<-as.data.frame(t(scale(t(tpm_9106_p_imvfft))))
head(log2_loess_scale_p)
log2_loess_scale_e<-as.data.frame(t(scale(t(tpm_9106_e_imvfft))))
head(log2_loess_scale_e)

pf9106_cbind_scale<-cbind(log2_loess_scale_e,log2_loess_scale_p)

data_ophase<-(pf9106_cbind_scale)[merge_ophase$IDs,]
pf9106_cbind_scale[1:10]
head(data_ophase)
dim(data_ophase)

library(circlize)
col_funev = colorRamp2(c(-2, 0, 2), c("green","black", "purple"))
#
dim(data_ophase)
#graphics.off()
habar<-HeatmapAnnotation(time = anno_barplot(pheno_9106$timepoint,
                                             gp = gpar(fontsize = 20,fill=2),
                                             bar_width = 1),
                         show_annotation_name = TRUE,
                         annotation_name_rot = 90,
                         annotation_name_side = "left",
                         height = unit(1.5,"cm"),
                         annotation_name_gp = gpar(fontsize=20))
habar

habar = HeatmapAnnotation(time=pheno_9106$timepoint,col = list(time=col_fun),
                          show_annotation_name = FALSE,
                          height = unit(2,"cm"))

#pheno_9106$condition<-gsub(pattern = "EV",replacement = "PfEVs",x = pheno_9106$condition)

#graphics.off()
jpeg(filename = "Fourier_PCA_RUVr2/9106/circabidian_outphase_9106_sigV3.jpeg",
     units = "in",width = 4,height = 8,res = 300)
htm_9106_phasic<-Heatmap(as.matrix(data_ophase[,rownames(pheno_9106)]),
                          clustering_method_rows = "average",
                          clustering_distance_rows = "canberra",
                          heatmap_legend_param = list(title = "log2 exp", 
                                                      direction="vertical"),
                          top_annotation = habar,
                          show_row_names = FALSE,
                          col = col_funev,
                          cluster_rows = FALSE,
                          show_column_names = FALSE,
                          row_names_side = "right", 
                          show_column_dend = FALSE,
                          cluster_columns = FALSE,
                          row_title = "4669 genes",
                          show_row_dend = FALSE, 
                          column_split = factor(pheno_9106$condition,
                                                levels = c("WP","PfEVs")), 
                          cluster_row_slices = FALSE,
                          cluster_column_slices = FALSE,
                          row_gap = unit(0, "mm"),use_raster = FALSE,
                          row_names_gp = gpar(fontsize = 12,col="black"),
                          row_names_max_width = unit(5,"cm"),
                          column_title_gp = gpar(fontsize = 35,col="black"),
                          row_title_gp = gpar(fontsize = 35,col="black"))
htm_9106_phasic

#+rowAnnotation(mark = anno)
htm_9106_phasicj<-draw(htm_9106_phasic, merge_legend = TRUE, 
                        heatmap_legend_side = "right", 
                        annotation_legend_side = "bottom",
                        show_heatmap_legend = FALSE,
                        column_title = "KE01", 
                        column_title_gp = gpar(fontsize = 35, col="black"))

dev.off()

peakdata<-merge(merge_SNR_p,merge_SNR_e, by="IDs")
head(peakdata)
colnames(peakdata)<-gsub(pattern = "\\.x",
                         replacement = "_9106_Parasite",
                         colnames(peakdata))

colnames(peakdata)<-gsub(pattern = "\\.y",
                         replacement = "_9106_EVs",
                         colnames(peakdata))
head(peakdata)

fit_df_wide_with_peakdata<-merge(peakdata, fit_df_wide, by="IDs")
head(fit_df_wide_with_peakdata)


l<-list(fit_df_wide_with_peakdata, tpm_9106_p_imvfft, tpm_9106_e_imvfft, pheno_9106)
l

write.xlsx(l,rowNames=TRUE, "Fourier_PCA_RUVr2/9106/Pf9106_KE01_Circacompare.xlsx")




