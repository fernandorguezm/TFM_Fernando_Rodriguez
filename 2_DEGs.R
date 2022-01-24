# Author: Fernando Rodríguez Marín
# Contact: fernando.rodriguez.marin8@gmail.com
# Project: master's thesis

library(qpcR)

##########################################################################
##########################################################################
##################### REPRESENTATION OF RESULTS ##########################
##########################################################################
##########################################################################

# MUTANTS

keyvals <- ifelse(
  res_temperature_16C10uEhy5$log2FoldChange < -1 & res_temperature_16C10uEhy5$padj < 0.01, 'navy',
  ifelse(res_temperature_16C10uEhy5$log2FoldChange > 1 & res_temperature_16C10uEhy5$padj < 0.01, 'red3',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'navy'] <- 'downregulated'
names(keyvals)[keyvals == 'black'] <- 'no regulated'
names(keyvals)[keyvals == 'red3'] <- 'upregulated'

png("volcanoplot_16C10uE.png", 800, 1000, pointsize=20)
EnhancedVolcano(res_temperature_16C10uEhy5,
                lab = NA,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Athy54 vs Atcol-0 16ºC 10uE',
                pCutoff = 0.01,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 0.0,
                colCustom = keyvals,
                colAlpha = 1, 
                gridlines.major = FALSE,
                gridlines.minor = FALSE,    
                legendPosition = 'right',
                legendLabSize = 0,
                legendIconSize = 0)
dev.off()


keyvals <- ifelse(
  res_temperature_16C100uEhy5$log2FoldChange < -1 & res_temperature_16C100uEhy5$padj < 0.01, 'navy',
  ifelse(res_temperature_16C100uEhy5$log2FoldChange > 1 & res_temperature_16C100uEhy5$padj < 0.01, 'red3',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'navy'] <- 'downregulated'
names(keyvals)[keyvals == 'black'] <- 'no regulated'
names(keyvals)[keyvals == 'red3'] <- 'upregulated'

png("volcanoplot_16C100uE.png", 800, 1000, pointsize=20)
EnhancedVolcano(res_temperature_16C100uEhy5,
                lab = NA,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Athy54 vs Atcol-0 16ºC 100uE',
                pCutoff = 0.01,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 0.0,
                colCustom = keyvals,
                colAlpha = 1, 
                gridlines.major = FALSE,
                gridlines.minor = FALSE,    
                legendPosition = 'right',
                legendLabSize = 0,
                legendIconSize = 0)
dev.off()



keyvals <- ifelse(
  res_temperature_24C10uEhy5$log2FoldChange < -1 & res_temperature_24C10uEhy5$padj < 0.01, 'navy',
  ifelse(res_temperature_24C10uEhy5$log2FoldChange > 1 & res_temperature_24C10uEhy5$padj < 0.01, 'red3',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'navy'] <- 'downregulated'
names(keyvals)[keyvals == 'black'] <- 'no regulated'
names(keyvals)[keyvals == 'red3'] <- 'upregulated'

png("volcanoplot_24C10uE.png", 800, 1000, pointsize=20)
EnhancedVolcano(res_temperature_24C10uEhy5,
                lab = NA,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Athy54 vs Atcol-0 24ºC 10uE',
                pCutoff = 0.01,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 0.0,
                colCustom = keyvals,
                colAlpha = 1, 
                gridlines.major = FALSE,
                gridlines.minor = FALSE,    
                legendPosition = 'right',
                legendLabSize = 0,
                legendIconSize = 0)
dev.off()


keyvals <- ifelse(
  res_temperature_24C100uEhy5$log2FoldChange < -1 & res_temperature_24C100uEhy5$padj < 0.01, 'navy',
  ifelse(res_temperature_24C100uEhy5$log2FoldChange > 1 & res_temperature_24C100uEhy5$padj < 0.01, 'red3',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'navy'] <- 'downregulated'
names(keyvals)[keyvals == 'black'] <- 'no regulated'
names(keyvals)[keyvals == 'red3'] <- 'upregulated'

png("volcanoplot_24C100uE.png", 800, 1000, pointsize=20)
EnhancedVolcano(res_temperature_24C100uEhy5,
                lab = NA,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Athy54 vs Atcol-0 24ºC 100uE',
                pCutoff = 0.01,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 0.0,
                colCustom = keyvals,
                colAlpha = 1, 
                gridlines.major = FALSE,
                gridlines.minor = FALSE,    
                legendPosition = 'right',
                legendLabSize = 0,
                legendIconSize = 0)
dev.off()


##########################################################################
##########################################################################
############################ LIST OF DEGs ################################
##########################################################################
##########################################################################

head(res_temperature_10uE)

#results 10uE fixed and 24ºC vs 16ºC

fold.change <- res_temperature_10uE$log2FoldChange
p.adj <- res_temperature_10uE$padj
length(fold.change)
genes.ids <- rownames(res_temperature_10uE)

fc<- log2(2)
up.genes_temperature_10uE <- na.omit(genes.ids[fold.change > fc & p.adj < 0.01])
down.genes_temperature_10uE <- na.omit(genes.ids[fold.change < -fc  & p.adj < 0.01])
DEG_temperature_10uE <- c(up.genes_temperature_10uE, down.genes_temperature_10uE)

length(up.genes_temperature_10uE)
length(down.genes_temperature_10uE)
length(DEG_temperature_10uE)

#results 100uE fixed and 24ºC vs 16ºC

fold.change <- res_temperature_100uE$log2FoldChange
p.adj <- res_temperature_100uE$padj
length(fold.change)
genes.ids <- rownames(res_temperature_100uE)

fc<- log2(2)
up.genes_temperature_100uE <- na.omit(genes.ids[fold.change > fc  & p.adj < 0.01])
down.genes_temperature_100uE <- na.omit(genes.ids[fold.change < -fc  & p.adj < 0.01])
DEG_temperature_100uE <- c(up.genes_temperature_100uE, down.genes_temperature_100uE)

length(up.genes_temperature_100uE)
length(down.genes_temperature_100uE)
length(DEG_temperature_100uE)

#results 16ºC fixed and hour 100uE vs 10uE

fold.change <- res_light_16ºC$log2FoldChange
p.adj <- res_light_16ºC$padj
length(fold.change)
genes.ids <- rownames(res_light_16ºC)

fc<- log2(2)
up.genes_light_16ºC <- na.omit(genes.ids[fold.change > fc  & p.adj < 0.01])
down.genes_light_16ºC <- na.omit(genes.ids[fold.change < -fc  & p.adj < 0.01])
DEG_light_16ºC <- c(up.genes_light_16ºC, down.genes_light_16ºC)

length(up.genes_light_16ºC)
length(down.genes_light_16ºC)
length(DEG_light_16ºC)

res_light_16ºC[DEG_light_16ºC,]

#results 24ºC fixed and hour 100uE vs 10uE

fold.change <- res_light_24ºC$log2FoldChange
p.adj <- res_light_24ºC$padj
length(fold.change)
genes.ids <- rownames(res_light_24ºC)

fc<- log2(2)
up.genes_light_24ºC <- na.omit(genes.ids[fold.change > fc  & p.adj < 0.01])
down.genes_light_24ºC <- na.omit(genes.ids[fold.change < -fc  & p.adj < 0.01])
DEG_light_24ºC <- c(up.genes_light_24ºC, down.genes_light_24ºC)

length(up.genes_light_24ºC)
length(down.genes_light_24ºC)
length(DEG_light_24ºC)

# temperature DEGs (combining the 2 fixed hour)

up.genes_temperature <- unique(c(up.genes_temperature_10uE,up.genes_temperature_100uE))
down.genes_temperature <- unique(c(down.genes_temperature_10uE,down.genes_temperature_100uE))
DEGs_temperature <- unique(c(DEG_temperature_10uE,DEG_temperature_100uE))

# light DEGs (combining the 2 fixed NaOH)

up.genes_light <- unique(c(up.genes_light_16ºC,up.genes_light_24ºC))
down.genes_light <- unique(c(down.genes_light_16ºC,down.genes_light_24ºC))
DEGs_light <- unique(c(DEG_light_16ºC,DEG_light_24ºC))

# All DEGs (combining all)

up.genes_all <- unique(c(up.genes_temperature,up.genes_light))
down.genes_all <- unique(c(down.genes_temperature,down.genes_light))
DEGs_all <- unique(c(DEGs_temperature,DEGs_light))

# Making table with all DEGs in all comparisons

table_DEGs <- qpcR:::cbind.na(up.genes_temperature_10uE, down.genes_temperature_10uE, DEG_temperature_10uE,
                                 up.genes_temperature_100uE, down.genes_temperature_100uE, DEG_temperature_100uE,
                                 up.genes_light_16ºC, down.genes_light_16ºC, DEG_light_16ºC,
                                 up.genes_light_24ºC, down.genes_light_24ºC, DEG_light_24ºC,
                                 up.genes_temperature, down.genes_temperature, DEGs_temperature,
                                 up.genes_light, down.genes_light, DEGs_light,
                                 up.genes_all, down.genes_all, DEGs_all)

table_DEGs

write.csv2(table_DEGs, file = "DEGs_Arabidopsis.csv")

