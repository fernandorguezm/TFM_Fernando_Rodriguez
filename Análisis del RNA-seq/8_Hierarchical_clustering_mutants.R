library(tidyverse)
library(cluster)
library(pheatmap)

##########################################################################
##########################################################################
################################# Atpif4 #################################
##########################################################################
##########################################################################

mean_TPM_atpif4 <- read.csv2(file = "Atpif4/Mean_TPM_atpif4.csv", header = T)
rownames(mean_TPM_atpif4) <- mean_TPM_atpif4[,1]
mean_TPM_atpif4 <- mean_TPM_atpif4[,2:5]
colnames(mean_TPM_atpif4) <- c("16ºC & 10uE", "16ºC & 100uE", "24ºC & 10uE", "24ºC & 100uE")

# Hierarchical clustering based on the Hugo Tavares' code
# https://tavareshugo.github.io/data-carpentry-rnaseq/04b_rnaseq_clustering.html
#trans_cts_atpif4 <- mean_expression_atpif4

trans_cts_atpif4 <- mean_TPM_atpif4

hclust_matrix_atpif4 <- trans_cts_atpif4

##### get counts for candidate genes ####

hclust_matrix_atpif4 <- hclust_matrix_atpif4 %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

DEGs_col0 <- read.csv2(file = "Atcol0/DEGs_Arabidopsis.csv")
DEGs_col0 <- DEGs_col0$DEGs_all_at

# Visualise gene expression trends per cluster

trans_cts_atpif4 <- as.data.frame(mean_TPM_atpif4)

# Adding Gene_id column

trans_cts_atpif4 <- add_column(trans_cts_atpif4, rownames(trans_cts_atpif4), .before = 1)

#Changing names of conditions to condition* to make easier the next steps.
colnames(trans_cts_atpif4)
colnames(trans_cts_atpif4) <- c("gene_id", "condition1", "condition2", "condition3", "condition4")
colnames(trans_cts_atpif4)

lowtemp_lowlight <- c("condition1", "16C", "10uE")
lowtemp_highlight <- c("condition2", "16C", "100uE")
hightemp_lowlight <- c("condition3", "24C", "10uE")
hightemp_highlight <- c("condition4", "24C", "100uE")

sample_info <- matrix(nrow = 4, ncol = 3)
sample_info <- as.data.frame(sample_info)
sample_info[1,] <- lowtemp_lowlight
sample_info[2,] <- lowtemp_highlight
sample_info[3,] <- hightemp_lowlight
sample_info[4,] <- hightemp_highlight

colnames(sample_info) <- c("condition", "temperature", "light")
rownames(sample_info) <- c("condition1", "condition2", "condition3", "condition4")

##########################################################################
##########################################################################
############################### Z-SCORE ##################################
##########################################################################
##########################################################################

trans_cts_mean_atpif4 <- trans_cts_atpif4 %>% 
  # convert to long format. Grouping all the samples columns in one
  pivot_longer(cols = condition1:condition4, names_to = "condition", values_to = "cts")  %>% 
  # join with sample info table. Adding the rest of colums (NaOH, hour, replicate...)
  full_join(sample_info, by = ("condition")) %>% 
  # filter to retain only genes of interest
  filter(gene_id %in% DEGs_col0) %>% ##CHANGE THIS IF Nº CLUSTER DOESN'T FIT
  # for each gene
  group_by(gene_id) %>%
  # scale the cts column
  mutate(cts_scaled = (cts - mean(cts))/sd(cts)) %>% 
  # for each gene, strain and minute
  group_by(gene_id, light, temperature) %>%
  # calculate the mean (scaled) cts
  summarise(mean_cts_scaled = mean(cts_scaled),
            nrep = n()) %>% 
  ungroup()

head(trans_cts_mean_atpif4)


indx <- c(2,3,4,5,6,7,8,10,11,19,21)
clusters <- read.csv("Atcol0/clusters_genes_at.csv", sep = ";")

i <- 1
for (i in 1:length(indx)) {
  
trans_cts_mean_atpif4_cluster <- trans_cts_mean_atpif4 %>%
filter(gene_id %in% clusters[,indx[i]+1])

lowtemperature_atpif4 <- trans_cts_mean_atpif4_cluster[trans_cts_mean_atpif4_cluster$temperature=="16C",]
highttemperature_atpif4 <- trans_cts_mean_atpif4_cluster[trans_cts_mean_atpif4_cluster$temperature=="24C",]

# Doing representation with an average of the expression

figure <- ggplot(NULL, aes(x=factor(light, level = c("10uE", "100uE")), y=mean_cts_scaled, group = gene_id)) + 
  geom_line(data=highttemperature_atpif4, color='#FDAE61', alpha = 0.2) + 
  geom_line(data=lowtemperature_atpif4, color='#4575B4', alpha = 0.2)+
  #geom_point(data = trans_cts_cluster, aes(x = hour, y = mean_cts_scaled, group = interaction(cluster,NaOH)), colour = "black", cex=1, alpha =0.1) +
  #stat_summary(data = highttemperature_atpif4,aes(group=1), fun=mean, geom="line", colour="#FDAE61", size = 2) + 
  #stat_summary(data = lowtemperature_atpif4,aes(group=1), fun=mean, geom="line", colour="#4575B4", size = 2)  +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=20),
                     axis.title=element_text(size=25),
                     strip.background = element_blank(),
                     strip.text = element_blank()) +
  labs(x = "Intensidad lumínica", y = "Z-score") +
  scale_colour_manual(values=c("orangered1", "dodgerblue1"), 
                      name  ="Temperature", 
                      labels=c("100mM", "0mM"))
ggsave(
  paste("Atcol0/cluster_images_Atpif4/cluster", indx[i], ".png", sep = ""),
  figure,
  width = 5,
  height = 5,
  dpi = 500
)
}

##########################################################################
##########################################################################
############################## Atphyb ####################################
##########################################################################
##########################################################################

mean_TPM_atphyb <- read.csv2(file = "Atphyb/Mean_TPM_atphyb.csv", header = T)
rownames(mean_TPM_atphyb) <- mean_TPM_atphyb[,1]
mean_TPM_atphyb <- mean_TPM_atphyb[,2:5]
colnames(mean_TPM_atphyb) <- c("16ºC & 10uE", "16ºC & 100uE", "24ºC & 10uE", "24ºC & 100uE")

# Hierarchical clustering based on the Hugo Tavares' code
# https://tavareshugo.github.io/data-carpentry-rnaseq/04b_rnaseq_clustering.html
#trans_cts_atphyb <- mean_expression_atphyb

trans_cts_atphyb <- mean_TPM_atphyb

hclust_matrix_atphyb <- trans_cts_atphyb

##### get counts for candidate genes ####

hclust_matrix_atphyb <- hclust_matrix_atphyb %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

DEGs_col0 <- read.csv2(file = "Atcol0/DEGs_Arabidopsis.csv")
DEGs_col0 <- DEGs_col0$DEGs_all_at

# Visualise gene expression trends per cluster

trans_cts_atphyb <- as.data.frame(mean_TPM_atphyb)

# Adding Gene_id column

trans_cts_atphyb <- add_column(trans_cts_atphyb, rownames(trans_cts_atphyb), .before = 1)

#Changing names of conditions to condition* to make easier the next steps.
colnames(trans_cts_atphyb)
colnames(trans_cts_atphyb) <- c("gene_id", "condition1", "condition2", "condition3", "condition4")
colnames(trans_cts_atphyb)

lowtemp_lowlight <- c("condition1", "16C", "10uE")
lowtemp_highlight <- c("condition2", "16C", "100uE")
hightemp_lowlight <- c("condition3", "24C", "10uE")
hightemp_highlight <- c("condition4", "24C", "100uE")

sample_info <- matrix(nrow = 4, ncol = 3)
sample_info <- as.data.frame(sample_info)
sample_info[1,] <- lowtemp_lowlight
sample_info[2,] <- lowtemp_highlight
sample_info[3,] <- hightemp_lowlight
sample_info[4,] <- hightemp_highlight

colnames(sample_info) <- c("condition", "temperature", "light")
rownames(sample_info) <- c("condition1", "condition2", "condition3", "condition4")

##########################################################################
##########################################################################
############################### Z-SCORE ##################################
##########################################################################
##########################################################################

trans_cts_mean_atphyb <- trans_cts_atphyb %>% 
  # convert to long format. Grouping all the samples columns in one
  pivot_longer(cols = condition1:condition4, names_to = "condition", values_to = "cts")  %>% 
  # join with sample info table. Adding the rest of colums (NaOH, hour, replicate...)
  full_join(sample_info, by = ("condition")) %>% 
  # filter to retain only genes of interest
  filter(gene_id %in% DEGs_col0) %>% ##CHANGE THIS IF Nº CLUSTER DOESN'T FIT
  # for each gene
  group_by(gene_id) %>%
  # scale the cts column
  mutate(cts_scaled = (cts - mean(cts))/sd(cts)) %>% 
  # for each gene, strain and minute
  group_by(gene_id, light, temperature) %>%
  # calculate the mean (scaled) cts
  summarise(mean_cts_scaled = mean(cts_scaled),
            nrep = n()) %>% 
  ungroup()

head(trans_cts_mean_atphyb)

indx <- c(2,3,4,5,6,7,8,10,11,19,21)
clusters <- read.csv("Atcol0/clusters_genes_at.csv", sep = ";")

i <- 1
for (i in 1:length(indx)) {
  
trans_cts_mean_atphyb_cluster <- trans_cts_mean_atphyb %>%
    filter(gene_id %in% clusters[,indx[i]+1])
  

lowtemperature_atphyb <- trans_cts_mean_atphyb_cluster[trans_cts_mean_atphyb_cluster$temperature=="16C",]
highttemperature_atphyb <- trans_cts_mean_atphyb_cluster[trans_cts_mean_atphyb_cluster$temperature=="24C",]

# Doing representation with an average of the expression

figure <-ggplot(NULL, aes(x=factor(light, level = c("10uE", "100uE")), y=mean_cts_scaled, group = gene_id)) + 
  geom_line(data=highttemperature_atphyb, color='#FDAE61', alpha = 0.2) + 
  geom_line(data=lowtemperature_atphyb, color='#4575B4', alpha = 0.2)+
  #geom_point(data = trans_cts_cluster, aes(x = hour, y = mean_cts_scaled, group = interaction(cluster,NaOH)), colour = "black", cex=1, alpha =0.1) +
  #stat_summary(data = highttemperature_atphyb,aes(group=1), fun=mean, geom="line", colour="#FDAE61", size = 2) + 
  #stat_summary(data = lowtemperature_atphyb,aes(group=1), fun=mean, geom="line", colour="#4575B4", size = 2)  +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=25),
                     axis.title=element_text(size=30),
                     strip.background = element_blank(),
                     strip.text = element_blank()) +
  labs(x = "light", y = "Z-score") +
  scale_colour_manual(values=c("orangered1", "dodgerblue1"), 
                      name  ="Temperature", 
                      labels=c("100mM", "0mM"))

ggsave(
  paste("Atcol0/cluster_images_Atphyb/cluster", indx[i], ".png", sep = ""),
  figure,
  width = 5,
  height = 5,
  dpi = 500
)
}

##########################################################################
##########################################################################
############################## Athy5 ####################################
##########################################################################
##########################################################################

mean_TPM_athy5 <- read.csv2(file = "Athy5/Mean_TPM_athy5.csv", header = T)
rownames(mean_TPM_athy5) <- mean_TPM_athy5[,1]
mean_TPM_athy5 <- mean_TPM_athy5[,2:5]
colnames(mean_TPM_athy5) <- c("16ºC & 10uE", "16ºC & 100uE", "24ºC & 10uE", "24ºC & 100uE")

# Hierarchical clustering based on the Hugo Tavares' code
# https://tavareshugo.github.io/data-carpentry-rnaseq/04b_rnaseq_clustering.html
#trans_cts_athy5 <- mean_expression_athy5

trans_cts_athy5 <- mean_TPM_athy5

hclust_matrix_athy5 <- trans_cts_athy5

##### get counts for candidate genes ####

hclust_matrix_athy5 <- hclust_matrix_athy5 %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

DEGs_col0 <- read.csv2(file = "Atcol0/DEGs_Arabidopsis.csv")
DEGs_col0 <- DEGs_col0$DEGs_all_at

# Visualise gene expression trends per cluster

trans_cts_athy5 <- as.data.frame(mean_TPM_athy5)

# Adding Gene_id column

trans_cts_athy5 <- add_column(trans_cts_athy5, rownames(trans_cts_athy5), .before = 1)

#Changing names of conditions to condition* to make easier the next steps.
colnames(trans_cts_athy5)
colnames(trans_cts_athy5) <- c("gene_id", "condition1", "condition2", "condition3", "condition4")
colnames(trans_cts_athy5)

lowtemp_lowlight <- c("condition1", "16C", "10uE")
lowtemp_highlight <- c("condition2", "16C", "100uE")
hightemp_lowlight <- c("condition3", "24C", "10uE")
hightemp_highlight <- c("condition4", "24C", "100uE")

sample_info <- matrix(nrow = 4, ncol = 3)
sample_info <- as.data.frame(sample_info)
sample_info[1,] <- lowtemp_lowlight
sample_info[2,] <- lowtemp_highlight
sample_info[3,] <- hightemp_lowlight
sample_info[4,] <- hightemp_highlight

colnames(sample_info) <- c("condition", "temperature", "light")
rownames(sample_info) <- c("condition1", "condition2", "condition3", "condition4")

##########################################################################
##########################################################################
############################### Z-SCORE ##################################
##########################################################################
##########################################################################

trans_cts_mean_athy5 <- trans_cts_athy5 %>% 
  # convert to long format. Grouping all the samples columns in one
  pivot_longer(cols = condition1:condition4, names_to = "condition", values_to = "cts")  %>% 
  # join with sample info table. Adding the rest of colums (NaOH, hour, replicate...)
  full_join(sample_info, by = ("condition")) %>% 
  # filter to retain only genes of interest
  filter(gene_id %in% DEGs_col0) %>% ##CHANGE THIS IF Nº CLUSTER DOESN'T FIT
  # for each gene
  group_by(gene_id) %>%
  # scale the cts column
  mutate(cts_scaled = (cts - mean(cts))/sd(cts)) %>% 
  # for each gene, strain and minute
  group_by(gene_id, light, temperature) %>%
  # calculate the mean (scaled) cts
  summarise(mean_cts_scaled = mean(cts_scaled),
            nrep = n()) %>% 
  ungroup()

head(trans_cts_mean_athy5)

indx <- c(2,3,4,5,6,7,8,10,11,19,21)
clusters <- read.csv("Atcol0/clusters_genes_at.csv", sep = ";")

i <- 1
for (i in 1:length(indx)) {
  
  trans_cts_mean_athy5_cluster <- trans_cts_mean_athy5 %>%
    filter(gene_id %in% clusters[,indx[i]+1])
  
  
  lowtemperature_athy5 <- trans_cts_mean_athy5_cluster[trans_cts_mean_athy5_cluster$temperature=="16C",]
  highttemperature_athy5 <- trans_cts_mean_athy5_cluster[trans_cts_mean_athy5_cluster$temperature=="24C",]
  
  # Doing representation with an average of the expression
  
  figure <-ggplot(NULL, aes(x=factor(light, level = c("10uE", "100uE")), y=mean_cts_scaled, group = gene_id)) + 
    geom_line(data=highttemperature_athy5, color='#FDAE61', alpha = 0.2) + 
    geom_line(data=lowtemperature_athy5, color='#4575B4', alpha = 0.2)+
    #geom_point(data = trans_cts_cluster, aes(x = hour, y = mean_cts_scaled, group = interaction(cluster,NaOH)), colour = "black", cex=1, alpha =0.1) +
    #stat_summary(data = highttemperature_athy5,aes(group=1), fun=mean, geom="line", colour="#FDAE61", size = 2) + 
    #stat_summary(data = lowtemperature_athy5,aes(group=1), fun=mean, geom="line", colour="#4575B4", size = 2)  +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                       axis.text=element_text(size=25),
                       axis.title=element_text(size=30),
                       strip.background = element_blank(),
                       strip.text = element_blank()) +
    labs(x = "light", y = "Z-score") +
    scale_colour_manual(values=c("orangered1", "dodgerblue1"), 
                        name  ="Temperature", 
                        labels=c("100mM", "0mM"))
  
  ggsave(
    paste("Atcol0/cluster_images_Athy5/cluster", indx[i], ".png", sep = ""),
    figure,
    width = 5,
    height = 5,
    dpi = 500
  )
}
