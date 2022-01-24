# Author: Fernando Rodríguez Marín
# Contact: fernando.rodriguez.marin8@gmail.com
# Project: master's thesis

##########################################################################
##########################################################################
########################### LOADING PACKAGES #############################
##########################################################################
##########################################################################

library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(readr)

##########################################################################
##########################################################################
####################### IMPORT & PRE-PROCESSING ##########################
##########################################################################
##########################################################################

# Import data from stringtie
## Previous commands line:
## Transcript assembly
## stringtie -G $WD/annotation/annotation.gtf -o sample.gtf -l sample sample.bam
## Quantification
## stringtie -e -B -G $WD/annotation/annotation.gtf -o sample.gtf sample.bam
## The script prepDE.py give you an input for DESeq2 (gene_count_matrix.csv)
## You can download it from: http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual

##########################################################################
##########################################################################
######################## ANALYSIS WITH DESeq2 ############################
##########################################################################
##########################################################################

# Loading metadata and gene expression data of the experiment.
pheno.data <- read.csv(file="pheno_data.csv")
pheno.data
gene.expression <- read.csv(file = "gene_count_matrix.csv")
rownames(gene.expression) <- gene.expression[,1]
gene.expression <- gene.expression[,2:13]
head(gene.expression)

# Run the DESeq pipeline for each comparison

# 1st comparison (comparison of 24ºC vs 16ºC with 10uE fixed)
pheno.data1 <- pheno.data[c(1,2,3,7,8,9), c(1,2,4)]
gene.expression1 <- gene.expression[,c(1,2,3,7,8,9)]
dds_temperature_10uE <- DESeqDataSetFromMatrix(countData=gene.expression1, colData=pheno.data1, design = ~ temperature)
dds_temperature_10uE <- DESeq(dds_temperature_10uE)
res_temperature_10uE <- results(dds_temperature_10uE, contrast = c("temperature", "24C", "16C"))

# 2nd comparison (comparison of 24ºC vs 16ºC with 100uE fixed))
pheno.data2 <- pheno.data[c(4,5,6,10,11,12), c(1,2,4)]
gene.expression2 <- gene.expression[,c(4,5,6,10,11,12)]
dds_temperature_100uE <- DESeqDataSetFromMatrix(countData=gene.expression2, colData=pheno.data2, design = ~ temperature)
dds_temperature_100uE <- DESeq(dds_temperature_100uE)
res_temperature_100uE <- results(dds_temperature_100uE, contrast = c("temperature", "24C", "16C"))

# 3rd comparison (comparison of 100uE vs 10uE with 16ºC fixed))
pheno.data3 <- pheno.data[1:6, c(1,3,4)]
gene.expression3 <- gene.expression[,1:6]
dds_light_16ºC <- DESeqDataSetFromMatrix(countData=gene.expression3, colData=pheno.data3, design = ~light)
dds_light_16ºC <- DESeq(dds_light_16ºC)
res_light_16ºC <- results(dds_light_16ºC, contrast = c("light", "100uE", "10uE"))

# 4th comparison (comparison of 100uE vs 10uE with 24ºC fixed))
pheno.data4 <- pheno.data[7:12, c(1,3,4)]
gene.expression4 <- gene.expression[,7:12]
dds_light_24ºC <- DESeqDataSetFromMatrix(countData=gene.expression4, colData=pheno.data4, design = ~ light)
dds_light_24ºC <- DESeq(dds_light_24ºC)
res_light_24ºC <- results(dds_light_24ºC, contrast = c("light", "100uE", "10uE"))

# Comparison with all samples

dds_pca <- DESeqDataSetFromMatrix(countData=gene.expression, colData=pheno.data, design = ~ temperature + light + temperature:light)
dds_pca <- DESeq(dds_pca)

##########################################################################
##########################################################################
################################### TPM ##################################
##########################################################################
##########################################################################

# Import data of TPM from stringtie. Stringtie generates a table:
# gene_abun.tab per sample. Convert this table to a .csv and read it. 

# Select gene_ID and TPM column

gene_abun1 <- read.csv(file = "sample1/gene_abun1.csv", sep = "\t")[,c(1,9)]
gene_abun2 <- read.csv(file = "sample2/gene_abun2.csv", sep = "\t")[,c(1,9)]
gene_abun3 <- read.csv(file = "sample3/gene_abun3.csv", sep = "\t")[,c(1,9)]
gene_abun4 <- read.csv(file = "sample4/gene_abun4.csv", sep = "\t")[,c(1,9)]
gene_abun5 <- read.csv(file = "sample5/gene_abun5.csv", sep = "\t")[,c(1,9)]
gene_abun6 <- read.csv(file = "sample6/gene_abun6.csv", sep = "\t")[,c(1,9)]
gene_abun25 <- read.csv(file = "sample25/gene_abun25.csv", sep = "\t")[,c(1,9)]
gene_abun26 <- read.csv(file = "sample26/gene_abun26.csv", sep = "\t")[,c(1,9)]
gene_abun27 <- read.csv(file = "sample27/gene_abun27.csv", sep = "\t")[,c(1,9)]
gene_abun28 <- read.csv(file = "sample28/gene_abun28.csv", sep = "\t")[,c(1,9)]
gene_abun29 <- read.csv(file = "sample29/gene_abun29.csv", sep = "\t")[,c(1,9)]
gene_abun30 <- read.csv(file = "sample30/gene_abun30.csv", sep = "\t")[,c(1,9)]

# Merge all the genes (ignore Warnings message, change the colnames later)
TPM <- merge(gene_abun1, gene_abun2, by = "Gene.ID")
TPM <- merge(TPM, gene_abun3, by = "Gene.ID")
TPM <- merge(TPM, gene_abun4, by = "Gene.ID")
TPM <- merge(TPM, gene_abun5, by = "Gene.ID")
TPM <- merge(TPM, gene_abun6, by = "Gene.ID")
TPM <- merge(TPM, gene_abun25, by = "Gene.ID")
TPM <- merge(TPM, gene_abun26, by = "Gene.ID")
TPM <- merge(TPM, gene_abun27, by = "Gene.ID")
TPM <- merge(TPM, gene_abun28, by = "Gene.ID")
TPM <- merge(TPM, gene_abun29, by = "Gene.ID")
TPM <- merge(TPM, gene_abun30, by = "Gene.ID")

colnames(TPM) <- c("gene_ID", "sample1", "sample2", "sample3", "sample4", "sample5", "sample6", "sample25", "sample26", "sample27", "sample28", "sample29", "sample30")
rownames(TPM) <- sub("gene:*", "", TPM[,1])
TPM <- TPM[,2:13]
write.csv(TPM, file = "TPM.csv")

##########################################################################
##########################################################################
########################## MEAN TPM MATRIX ###############################
##########################################################################
##########################################################################

head(TPM)
lowtemp_lowlight <- (TPM$sample1 + TPM$sample2 + TPM$sample3)/3
lowtemp_highlight <- (TPM$sample3 + TPM$sample4 + TPM$sample6)/3
hightemp_lowlight <- (TPM$sample25 + TPM$sample26 + TPM$sample27)/3
hightemp_highlight <- (TPM$sample28 + TPM$sample29 + TPM$sample30)/3


mean_TPM <- matrix(c(rownames(TPM), lowtemp_lowlight, lowtemp_highlight, hightemp_lowlight, hightemp_highlight), ncol = 5)
colnames(mean_TPM) <- c("gene_id","16ºC & 10uE", "16ºC & 100uE", "24ºC & 10uE", "24ºC & 100uE")
head(mean_TPM)

rownames(mean_TPM) <- mean_TPM[,1]
mean_TPM <- mean_TPM[,2:5]
mean_TPM <- matrix(as.numeric(mean_TPM),    # Convert it to numeric matrix
                             ncol = ncol(mean_TPM))
colnames(mean_TPM) <- c("16ºC & 10uE", "16ºC & 100uE", "24ºC & 10uE", "24ºC & 100uE")
rownames(mean_TPM) <- rownames(TPM)
head(mean_TPM)

write.csv2(mean_TPM, file = "Mean_TPM.csv")

##########################################################################
##########################################################################
########################## REPRESENTING DDS ##############################
##########################################################################
##########################################################################

# Applying variance stabilizing transformation
vst <- vst(dds_pca)
head(assay(vst))
par(mfrow=c(1,1))
hist(assay(vst))

#################################################
## Select data for the 500 most variable genes ##
#################################################

# We estimate the variance for each row in the TPM matrix
var_genes <- apply(assay(vst), 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)
# Subset TPM matrix
highly_variable_TPM <- vst[select_var,]
dim(highly_variable_TPM)
head(highly_variable_TPM)

# Principal components analysis
## This function belongs to Stephen Turner, @genetics_blog: 
## https://gist.github.com/stephenturner/f60c1934405c127f09a6
                                            
vst_pca <- function (rld, intgroup ="hour", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)                        
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], cex=1.5, pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  # to label all the points
  # with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)),
  #                                 cex = 0.3))
  #legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #   rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
mycols <- c("#4575B4", "#74ADD1", "#D73027", "#FDAE61")

png("PCA.png", 1000, 1000, pointsize=30)
vst_pca(highly_variable_TPM, colors=mycols, intgroup=c("temperature", "light"),
        ylim=c(-40,40), xlim=c(-40,40))
dev.off()


