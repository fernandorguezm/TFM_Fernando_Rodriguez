# Author: Fernando Rodríguez Marín
# Contact: fernando.rodriguez.marin8@gmail.com
# Project: master's thesis

library(tidyverse)
library(cluster)
library(pheatmap)
library(eulerr)

##########################################################################
##########################################################################
############################ ARABIDOPSIS #################################
##########################################################################
##########################################################################

# Reading TPM and DEGs data

DEGs <- read.csv("DEGs_venny_mutants_at.csv", sep = ";")$All

mean_TPM_at <- read.csv2(file = "Atcol0/Mean_TPM_arabidopsis_Col0.csv", header = T)
rownames(mean_TPM_at) <- mean_TPM_at[,1]
mean_TPM_at <- mean_TPM_at[,2:5]
colnames(mean_TPM_at) <- c("16ºC & 10uE col0", "16ºC & 100uE col0", "24ºC & 10uE col0", "24ºC & 100uE col0")

mean_TPM_atpif <- read.csv2(file = "Atpif4/Mean_TPM_atpif4.csv", header = T)
rownames(mean_TPM_atpif) <- mean_TPM_atpif[,1]
mean_TPM_atpif <- mean_TPM_atpif[,2:5]
colnames(mean_TPM_atpif) <- c("16ºC & 10uE Atpif4", "16ºC & 100uE Atpif4", "24ºC & 10uE Atpif4", "24ºC & 100uE Atpif4")

mean_TPM_atphyb <- read.csv2(file = "Atphyb/Mean_TPM_atphyb.csv", header = T)
rownames(mean_TPM_atphyb) <- mean_TPM_atphyb[,1]
mean_TPM_atphyb <- mean_TPM_atphyb[,2:5]
colnames(mean_TPM_atphyb) <- c("16ºC & 10uE Atphyb", "16ºC & 100uE Atphyb", "24ºC & 10uE Atphyb", "24ºC & 100uE Atphyb")

mean_TPM_athy5 <- read.csv2(file = "Athy5/Mean_TPM_athy5.csv", header = T)
rownames(mean_TPM_athy5) <- mean_TPM_athy5[,1]
mean_TPM_athy5 <- mean_TPM_athy5[,2:5]
colnames(mean_TPM_athy5) <- c("16ºC & 10uE Athy5", "16ºC & 100uE Athy5", "24ºC & 10uE Athy5", "24ºC & 100uE Athy5")

TPMs <- cbind(mean_TPM_at, mean_TPM_atpif, mean_TPM_atphyb, mean_TPM_athy5)
TPMs <- TPMs[DEGs,]

trans_cts_at <- TPMs

hclust_matrix_at <- trans_cts_at
hclust_matrix_at <- mutate_all(hclust_matrix_at, function(x) as.numeric(as.character(x)))
##### get counts for candidate genes ####

hclust_matrix_at <- hclust_matrix_at %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

hclust_matrix_at <- subset(hclust_matrix_at, rownames(hclust_matrix_at) %in% DEGs)
nrow(hclust_matrix_at)

##########################################################################
##########################################################################
############################## HEATMAP ###################################
##########################################################################
##########################################################################

# Preparing data for annotation (data of samples)

lowtemp_lowlight <- c("16ºC", "10uE")
lowtemp_highlight <- c("16ºC", "100uE")
hightemp_lowlight <- c("24ºC", "10uE")
hightemp_highlight <- c("24ºC", "100uE")

pheno.datahm <- matrix(nrow = 4, ncol = 2)
pheno.datahm <- as.data.frame(pheno.datahm)
pheno.datahm[1,] <- lowtemp_lowlight
pheno.datahm[2,] <- lowtemp_highlight
pheno.datahm[3,] <- hightemp_lowlight
pheno.datahm[4,] <- hightemp_highlight

colnames(pheno.datahm) <- c("temperature", "light")
rownames(pheno.datahm) <- c("16ºC/10uE", "16ºC/100uE", "24ºC/10uE", "24ºC/100uE")
t(pheno.datahm)

DEGs <- read.csv("DEGs_venny_mutants_at.csv", sep = ";")

rowanno <- DEGs[,8:9]
rownames(rowanno) <- rowanno[,1]
rowanno <- rowanno[,-1,  drop=FALSE]

paletteLength <- 50
myColor <- colorRampPalette(c("navy", "white", "firebrick3"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(hclust_matrix_at), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(hclust_matrix_at)/paletteLength, max(hclust_matrix_at), length.out=floor(paletteLength/2)))
  
# Complete heatmap

png("Heatmap_venny_annotated.png", 500, 1500, pointsize = 50)
pheatmap(hclust_matrix_at,
         annotation_row = rowanno,
         show_rownames =F,
         cluster_cols = F,
         cluster_rows = F,
         fontsize = 7,
         cutree_cols = 3,
         color=myColor, 
         breaks=myBreaks,
         gaps_col=c(4,8,12))
dev.off()

# Heatmap of each part

png("Heatmap_venny_annotated_athy5.png", 500, length(na.omit(DEGs$Athy5))/8, pointsize = 50)
pheatmap(hclust_matrix_at[na.omit(DEGs$Athy5),],
         annotation_row = rowanno,
         show_rownames =F,
         cluster_cols = F,
         cluster_rows = T,
         fontsize = 7,
         cutree_cols = 3,
         color=myColor, 
         breaks=myBreaks,
         gaps_col=c(4,8,12))
dev.off()

png("Heatmap_venny_annotated_atphyb.png", 500, length(na.omit(DEGs$Atphyb))/8, pointsize = 50)
pheatmap(hclust_matrix_at[na.omit(DEGs$Atphyb),],
         annotation_row = rowanno,
         show_rownames =F,
         cluster_cols = F,
         cluster_rows = T,
         fontsize = 7,
         cutree_cols = 3,
         color=myColor, 
         breaks=myBreaks,
         gaps_col=c(4,8,12))
dev.off()

png("Heatmap_venny_annotated_atpif4.png", 500, length(na.omit(DEGs$Atpif4))/1.6, pointsize = 1)
pheatmap(hclust_matrix_at[na.omit(DEGs$Atpif4),],
         annotation_row = rowanno,
         show_rownames =F,
         show_colnames = T,
         cluster_cols = F,
         cluster_rows = T,
         fontsize = 7,
         cutree_cols = 3,
         color=myColor, 
         breaks=myBreaks,
         gaps_col=c(4,8,12))
dev.off()

png("Heatmap_venny_annotated_athy5_atphyb.png", 500, length(na.omit(DEGs$Athy5.Atphyb))/2, pointsize = 50)
pheatmap(hclust_matrix_at[na.omit(DEGs$Athy5.Atphyb),],
         annotation_row = rowanno,
         show_rownames =F,
         show_colnames = F,
         annotation_legend = F,
         legend = F,
         cluster_cols = F,
         cluster_rows = T,
         fontsize = 7,
         cutree_cols = 3,
         color=myColor, 
         breaks=myBreaks,
         gaps_col=c(4,8,12))
dev.off()
png("Heatmap_venny_annotated_athy5_atpif4.png", 500, length(na.omit(DEGs$Athy5.Atpif4))/2, pointsize = 50)
pheatmap(hclust_matrix_at[na.omit(DEGs$Athy5.Atpif4),],
         annotation_row = rowanno,
         show_rownames =F,
         show_colnames = F,
         annotation_legend = F,
         legend = F,
         cluster_cols = F,
         cluster_rows = T,
         fontsize = 7,
         cutree_cols = 3,
         color=myColor, 
         breaks=myBreaks,
         gaps_col=c(4,8,12))
dev.off()
png("Heatmap_venny_annotated_atphyb_atpif4.png", 500, length(na.omit(DEGs$Atphyb.Atpif4))/2, pointsize = 2)
pheatmap(hclust_matrix_at[na.omit(DEGs$Atphyb.Atpif4),],
         annotation_row = rowanno,
         show_rownames =F,
         show_colnames = F,
         annotation_legend = F,
         legend = F,
         cluster_cols = F,
         cluster_rows = T,
         border_color=NA,
         fontsize = 7,
         cutree_cols = 3,
         color=myColor,
         breaks=myBreaks,
         gaps_col=c(4,8,12))
dev.off()
png("Heatmap_venny_legend_at.png", 200, 200, pointsize = 50)
pheatmap(hclust_matrix_at[na.omit(DEGs$Athy5.Atphyb.Atpif4),],
         annotation_row = rowanno,
         show_rownames =F,
         show_colnames = F,
         annotation_legend = F,
         legend = T,
         cluster_cols = F,
         cluster_rows = T,
         fontsize = 7,
         cutree_cols = 3,
         color=myColor, 
         breaks=myBreaks,
         gaps_col=c(4,8,12))
dev.off()

##########################################################################
##########################################################################
############################ VENN DIAGRAM ################################
##########################################################################
##########################################################################

VennDiag <- euler(c("A" = 2570, "B" = 905, "C" = 151, "A&B" = 806, "B&C" = 118, 
        "A&C" = 160, "A&B&C" = 489))

png("Venn_Atmutants.png", 1000, 1000, pointsize=100)
plot(VennDiag, fontsize=1000, cex=1000, alpha=0.7,
     fill=c("grey","firebrick2", "navy"),
     #edges = list(col = c("black","black", "black"), lex = 2), 
     quantities = list(cex = 3),
     labels = c("", "", ""))
dev.off()

# ##########################################################################
# ##########################################################################
# ####################### PEARSON CORRELATION ##############################
# ##########################################################################
# ##########################################################################
# 
# 
# DEGs <- read.csv("DEGs_venny_mutants_at.csv", sep = ";")$All
# 
# mean_TPM_at <- read.csv2(file = "Atcol0/Mean_TPM_arabidopsis_Col0.csv", header = T)
# rownames(mean_TPM_at) <- mean_TPM_at[,1]
# mean_TPM_at <- mean_TPM_at[,2:5]
# colnames(mean_TPM_at) <- c("16ºC & 10uE col0", "16ºC & 100uE col0", "24ºC & 10uE col0", "24ºC & 100uE col0")
# 
# mean_TPM_atpif <- read.csv2(file = "Atpif4/Mean_TPM_atpif4.csv", header = T)
# rownames(mean_TPM_atpif) <- mean_TPM_atpif[,1]
# mean_TPM_atpif <- mean_TPM_atpif[,2:5]
# colnames(mean_TPM_atpif) <- c("16ºC & 10uE Atpif4", "16ºC & 100uE Atpif4", "24ºC & 10uE Atpif4", "24ºC & 100uE Atpif4")
# 
# mean_TPM_atphyb <- read.csv2(file = "Atphyb/Mean_TPM_atphyb.csv", header = T)
# rownames(mean_TPM_atphyb) <- mean_TPM_atphyb[,1]
# mean_TPM_atphyb <- mean_TPM_atphyb[,2:5]
# colnames(mean_TPM_atphyb) <- c("16ºC & 10uE Atphyb", "16ºC & 100uE Atphyb", "24ºC & 10uE Atphyb", "24ºC & 100uE Atphyb")
# 
# mean_TPM_athy5 <- read.csv2(file = "Athy5/Mean_TPM_athy5.csv", header = T)
# rownames(mean_TPM_athy5) <- mean_TPM_athy5[,1]
# mean_TPM_athy5 <- mean_TPM_athy5[,2:5]
# colnames(mean_TPM_athy5) <- c("16ºC & 10uE Athy5", "16ºC & 100uE Athy5", "24ºC & 10uE Athy5", "24ºC & 100uE Athy5")
# 
# TPMs <- cbind(mean_TPM_at, mean_TPM_atpif, mean_TPM_atphyb, mean_TPM_athy5)
# 
# library(rstatix)
# to<- as.tibble(TPMs)
# rownames(to) <- rownames(TPMs)
# gene.correlation <- cor_mat(to, method = "pearson")
# dim(gene.correlation)
# rstatix::cor_mark_significant(gene.correlation)
# gene.correlation[1:5, 1:5]
# gene.correlation <- as.data.frame(gene.correlation)
# rownames(gene.correlation) <- gene.correlation[,1]
# gene.correlation <- gene.correlation[,-1]
# gene.correlation <- mutate_all(gene.correlation, function(x) as.numeric(as.character(x)))
# gene.correlation[upper.tri(gene.correlation)] = NA
# 
# paletteLength <- 60
# myColor <- colorRampPalette(c("navy", "white", "firebrick3"))(paletteLength)
# # length(breaks) == length(paletteLength) + 1
# # use floor and ceiling to deal with even/odd length pallettelengths
# breaksList = seq(0.4, 1, by = 0.01)
# png("Correlation_Pearson_Matrix_2.png", 500, 500, pointsize = 50)
# pheatmap(gene.correlation, cluster_rows=F, cluster_cols=F, na_col="white", 
#          breaks = breaksList, color = myColor)
# dev.off()

##########################################################################
##########################################################################
############################# MARCHANTIA ########@@#######################
##########################################################################
##########################################################################


DEGs_mp <- read.csv("DEGs_venny_mutants_mp.csv", sep = ";")$All

mean_TPM_mp <- read.csv2(file = "Mptak1/Mean_TPM_marchantia_Tak1.csv", header = T)
rownames(mean_TPM_mp) <- mean_TPM_mp[,1]
mean_TPM_mp <- mean_TPM_mp[,2:5]
colnames(mean_TPM_mp) <- c("16ºC & 10uE tak1", "16ºC & 100uE tak1", "24ºC & 10uE tak1", "24ºC & 100uE tak1")

mean_TPM_mppif <- read.csv2(file = "Mppif/Mean_TPM_mppif.csv", header = T)
rownames(mean_TPM_mppif) <- mean_TPM_mppif[,1]
mean_TPM_mppif <- mean_TPM_mppif[,2:5]
colnames(mean_TPM_mppif) <- c("16ºC & 10uE Mppif", "16ºC & 100uE Mppif", "24ºC & 10uE Mppif", "24ºC & 100uE Mppif")

mean_TPM_mpphyb <- read.csv2(file = "Mpphy/Mean_TPM_mpphy.csv", header = T)
rownames(mean_TPM_mpphyb) <- mean_TPM_mpphyb[,1]
mean_TPM_mpphyb <- mean_TPM_mpphyb[,2:5]
colnames(mean_TPM_mpphyb) <- c("16ºC & 10uE Mpphy", "16ºC & 100uE Mpphy", "24ºC & 10uE Mpphy", "24ºC & 100uE Mpphy")

mean_TPM_mphy5 <- read.csv2(file = "Mphy5/Mean_TPM_mphy5.csv", header = T)
rownames(mean_TPM_mphy5) <- mean_TPM_mphy5[,1]
mean_TPM_mphy5 <- mean_TPM_mphy5[,2:5]
colnames(mean_TPM_mphy5) <- c("16ºC & 10uE Mphy5", "16ºC & 100uE Mphy5", "24ºC & 10uE Mphy5", "24ºC & 100uE Mphy5")

TPMs <- cbind(mean_TPM_mp, mean_TPM_mppif, mean_TPM_mpphyb, mean_TPM_mphy5)
TPMs <- TPMs[DEGs_mp,]

trans_cts_mp <- TPMs

hclust_matrix_mp <- trans_cts_mp
hclust_matrix_mp <- mutate_all(hclust_matrix_mp, function(x) as.numeric(as.character(x)))
##### get counts for candidate genes ####

hclust_matrix_mp <- hclust_matrix_mp %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

hclust_matrix_mp <- subset(hclust_matrix_mp, rownames(hclust_matrix_mp) %in% DEGs_mp)
nrow(hclust_matrix_mp)

##########################################################################
##########################################################################
############################## HEATMAP ###################################
##########################################################################
##########################################################################

# Preparing data for annotation (data of samples)

lowtemp_lowlight <- c("16ºC", "10uE")
lowtemp_highlight <- c("16ºC", "100uE")
hightemp_lowlight <- c("24ºC", "10uE")
hightemp_highlight <- c("24ºC", "100uE")

pheno.datahm <- matrix(nrow = 4, ncol = 2)
pheno.datahm <- as.data.frame(pheno.datahm)
pheno.datahm[1,] <- lowtemp_lowlight
pheno.datahm[2,] <- lowtemp_highlight
pheno.datahm[3,] <- hightemp_lowlight
pheno.datahm[4,] <- hightemp_highlight

colnames(pheno.datahm) <- c("temperature", "light")
rownames(pheno.datahm) <- c("16ºC/10uE", "16ºC/100uE", "24ºC/10uE", "24ºC/100uE")
t(pheno.datahm)

DEGs <- read.csv("DEGs_venny_mutants_mp.csv", sep = ";")

rowanno <- DEGs[,8:9]
rownames(rowanno) <- rowanno[,1]
rowanno <- rowanno[,-1,  drop=FALSE]

paletteLength <- 50
myColor <- colorRampPalette(c("navy", "white", "firebrick3"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(hclust_matrix_mp), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(hclust_matrix_mp)/paletteLength, max(hclust_matrix_mp), length.out=floor(paletteLength/2)))


png("Heatmap_venny_annotated_mp.png", 500, 1500, pointsize = 50)
pheatmap(hclust_matrix_mp,
         annotation_row = rowanno,
         show_rownames =F,
         cluster_cols = F,
         cluster_rows = F,
         fontsize = 7,
         cutree_cols = 3,
         color=myColor, 
         breaks=myBreaks,
         gaps_col=c(4,8,12))
dev.off()



png("Heatmap_venny_annotated_mphy5_mpphy.png", 500, length(na.omit(DEGs$Mphy5.Mpphy))/2, pointsize = 50)
pheatmap(hclust_matrix_mp[na.omit(DEGs$Mphy5.Mpphy),],
         annotation_row = rowanno,
         show_rownames =F,
         show_colnames = F,
         annotation_legend = F,
         legend = F,
         cluster_cols = F,
         cluster_rows = T,
         fontsize = 7,
         cutree_cols = 3,
         color=myColor, 
         breaks=myBreaks,
         gaps_col=c(4,8,12))
dev.off()
png("Heatmap_venny_annotated_mphy5_mppif.png", 500, length(na.omit(DEGs$Mphy5.Mppif))/2, pointsize = 50)
pheatmap(hclust_matrix_mp[na.omit(DEGs$Mphy5.Mppif),],
         annotation_row = rowanno,
         show_rownames =F,
         show_colnames = F,
         annotation_legend = F,
         legend = F,
         cluster_cols = F,
         cluster_rows = T,
         fontsize = 7,
         cutree_cols = 3,
         border_color = NA,
         color=myColor, 
         breaks=myBreaks,
         gaps_col=c(4,8,12))
dev.off()
png("Heatmap_venny_annotated_mpphy_mppif.png", 500, length(na.omit(DEGs$Mpphy.Mppif))/2, pointsize = 2)
pheatmap(hclust_matrix_mp[na.omit(DEGs$Mpphy.Mppif),],
         annotation_row = rowanno,
         show_rownames =F,
         show_colnames = F,
         annotation_legend = F,
         legend = F,
         cluster_cols = F,
         cluster_rows = T,
         border_color=NA,
         fontsize = 7,
         cutree_cols = 3,
         color=myColor,
         breaks=myBreaks,
         gaps_col=c(4,8,12))
dev.off()
png("Heatmap_venny_annotated_mphy5_mpphy_mppif.png", 500, length(na.omit(DEGs$Mphy5.Mpphy.Mppif))/2, pointsize = 50)
pheatmap(hclust_matrix_mp[na.omit(DEGs$Mphy5.Mpphy.Mppif),],
         annotation_row = rowanno,
         show_rownames =F,
         show_colnames = F,
         annotation_legend = F,
         legend = F,
         cluster_cols = F,
         cluster_rows = T,
         fontsize = 7,
         cutree_cols = 3,
         color=myColor, 
         breaks=myBreaks,
         gaps_col=c(4,8,12))
dev.off()

##########################################################################
##########################################################################
############################ VENN DIAGRAM ################################
##########################################################################
##########################################################################

VennDiag <- euler(c("A" = 123, "B" = 2010, "C" = 273, "A&B" = 187, "B&C" = 1077, 
                    "A&C" = 67, "A&B&C" = 746))

png("Venn_Mpmutants.png", 1000, 1000, pointsize=100)
plot(VennDiag, fontsize=1000, cex=1000, alpha=0.7,
     fill=c("grey","firebrick2", "navy"),
     #edges = list(col = c("black","black", "black"), lex = 2), 
     quantities = list(cex = 3),
     labels = c("", "", ""))
dev.off()


# ##########################################################################
# ##########################################################################
# ####################### PEARSON CORRELATION ##############################
# ##########################################################################
# ##########################################################################
# 
# 
# DEGs <- read.csv("DEGs_venny_mutants_mp.csv", sep = ";")$All
# 
# mean_TPM_mp <- read.csv2(file = "Mptak1/Mean_TPM_marchantia_Tak1.csv", header = T)
# rownames(mean_TPM_mp) <- mean_TPM_mp[,1]
# mean_TPM_mp <- mean_TPM_mp[,2:5]
# colnames(mean_TPM_mp) <- c("16ºC & 10uE tak1", "16ºC & 100uE tak1", "24ºC & 10uE tak1", "24ºC & 100uE tak1")
# 
# mean_TPM_mppif <- read.csv2(file = "Mppif/Mean_TPM_mppif.csv", header = T)
# rownames(mean_TPM_mppif) <- mean_TPM_mppif[,1]
# mean_TPM_mppif <- mean_TPM_mppif[,2:5]
# colnames(mean_TPM_mppif) <- c("16ºC & 10uE Mppif", "16ºC & 100uE Mppif", "24ºC & 10uE Mppif", "24ºC & 100uE Mppif")
# 
# mean_TPM_mpphyb <- read.csv2(file = "Mpphy/Mean_TPM_mpphy.csv", header = T)
# rownames(mean_TPM_mpphyb) <- mean_TPM_mpphyb[,1]
# mean_TPM_mpphyb <- mean_TPM_mpphyb[,2:5]
# colnames(mean_TPM_mpphyb) <- c("16ºC & 10uE Mpphy", "16ºC & 100uE Mpphy", "24ºC & 10uE Mpphy", "24ºC & 100uE Mpphy")
# 
# mean_TPM_mphy5 <- read.csv2(file = "Mphy5/Mean_TPM_mphy5.csv", header = T)
# rownames(mean_TPM_mphy5) <- mean_TPM_mphy5[,1]
# mean_TPM_mphy5 <- mean_TPM_mphy5[,2:5]
# colnames(mean_TPM_mphy5) <- c("16ºC & 10uE Mphy5", "16ºC & 100uE Mphy5", "24ºC & 10uE Mphy5", "24ºC & 100uE Mphy5")
# 
# TPMs <- cbind(mean_TPM_mp, mean_TPM_mppif, mean_TPM_mpphyb, mean_TPM_mphy5)
# 
# library(rstatix)
# to<- as.tibble(TPMs)
# rownames(to) <- rownames(TPMs)
# gene.correlation <- cor_mat(to, method = "pearson")
# dim(gene.correlation)
# rstatix::cor_mark_significant(gene.correlation)
# gene.correlation[1:5, 1:5]
# gene.correlation <- as.data.frame(gene.correlation)
# rownames(gene.correlation) <- gene.correlation[,1]
# gene.correlation <- gene.correlation[,-1]
# gene.correlation <- mutate_all(gene.correlation, function(x) as.numeric(as.character(x)))
# gene.correlation[upper.tri(gene.correlation)] = NA
# 
# paletteLength <- 30
# myColor <- colorRampPalette(c("navy", "white", "firebrick3"))(paletteLength)
# # length(breaks) == length(paletteLength) + 1
# # use floor and ceiling to deal with even/odd length pallettelengths
# breaksList = seq(0.7, 1, by = 0.01)
# 
# png("Correlation_Pearson_Matrix_mp_2.png", 500, 500, pointsize = 50)
# pheatmap(gene.correlation, cluster_rows=F, cluster_cols=F, na_col="white", 
#          breaks = breaksList, color = myColor)
# dev.off()
# 
# DEGs <- read.csv("DEGs_venny_mutants_mp.csv", sep = ";")
     
