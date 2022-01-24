# Author: Fernando Rodríguez Marín
# Contact: fernando.rodriguez.marin8@gmail.com
# Project: master's thesis

library(tidyverse)
library(cluster)
library(pheatmap)
library(qpcR)

##########################################################################
##########################################################################
###################### HIERARCHICAL CLUSTERING ###########################
##########################################################################
##########################################################################

mean_TPM <- read.csv2(file = "Mean_TPM.csv", header = T)
rownames(mean_TPM) <- mean_TPM[,1]
mean_TPM <- mean_TPM[,2:5]
colnames(mean_TPM) <- c("16ºC & 10uE", "16ºC & 100uE", "24ºC & 10uE", "24ºC & 100uE")

# Hierarchical clustering based on the Hugo Tavares' code
# https://tavareshugo.github.io/data-carpentry-rnaseq/04b_rnaseq_clustering.html
#trans_cts <- mean_expression

trans_cts <- mean_TPM

hclust_matrix <- trans_cts

##### get counts for candidate genes ####

hclust_matrix <- hclust_matrix %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

DEGs_table <- read.csv2(file = "DEGs.csv")
DEGs <- DEGs_table$DEGs_all
hclust_matrix <- subset(hclust_matrix, rownames(hclust_matrix) %in% DEGs)
nrow(hclust_matrix)
gene_dist <- dist(hclust_matrix)

# Perform hierarchical clustering using hclust()
gene_hclust <- hclust(gene_dist)

# The default `plot()` function can be used to produce a simple dendrogram
plot(gene_hclust, labels = FALSE)
abline(h = 1.75, col = "green", lwd = 2) # add horizontal line to illustrate cutting dendrogram

# Now we have if th silhouette of the clustering is not too much affected

i <- 2
for (i in 2:30) {
  variable <- cutree(gene_hclust,k=i)
  nam = paste("hclust.",i, sep = "")
  nam2 = paste("sil",i, sep = "")
  assign(nam, variable)
  variable2 <- silhouette(variable, gene_dist)
  assign(nam2, variable2)
  
}

hclust.sil.values <- c(summary(sil2)[["avg.width"]],summary(sil3)[["avg.width"]],summary(sil4)[["avg.width"]],summary(sil5)[["avg.width"]],summary(sil6)[["avg.width"]],summary(sil7)[["avg.width"]],summary(sil8)[["avg.width"]],summary(sil9)[["avg.width"]],summary(sil10)[["avg.width"]],summary(sil11)[["avg.width"]],summary(sil12)[["avg.width"]],summary(sil13)[["avg.width"]],summary(sil14)[["avg.width"]],summary(sil15)[["avg.width"]],summary(sil16)[["avg.width"]],summary(sil17)[["avg.width"]],summary(sil18)[["avg.width"]],summary(sil19)[["avg.width"]],summary(sil20)[["avg.width"]],summary(sil21)[["avg.width"]],summary(sil22)[["avg.width"]],summary(sil23)[["avg.width"]],summary(sil24)[["avg.width"]],summary(sil25)[["avg.width"]],summary(sil26)[["avg.width"]],summary(sil27)[["avg.width"]],summary(sil28)[["avg.width"]],summary(sil29)[["avg.width"]],summary(sil30)[["avg.width"]])
plot(2:30,hclust.sil.values,type="o",col="blue",pch=0,ylim=c(0,1),xlab="Number of clusters",ylab="Silhouette",lwd=3)
max(hclust.sil.values)

# Continue with K = clustering number 

gene_cluster <- cutree(gene_hclust, k = 30) %>% 
  # turn the named vector into a tibble
  enframe() %>%
  # rename some of the columns
  rename(gene_id = name, cluster = value)

head(gene_cluster)

# Visualise gene expression trends per cluster

trans_cts <- as.data.frame(mean_TPM)

# Adding Gene_id column

trans_cts <- add_column(trans_cts, rownames(trans_cts), .before = 1)

#Changing names of conditions to condition* to make easier the next steps.
colnames(trans_cts)
colnames(trans_cts) <- c("gene_id", "condition1", "condition2", "condition3", "condition4")
colnames(trans_cts)

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

trans_cts_mean <- trans_cts %>% 
  # convert to long format. Grouping all the samples columns in one
  pivot_longer(cols = condition1:condition4, names_to = "condition", values_to = "cts")  %>% 
  # join with sample info table. Adding the rest of colums (NaOH, hour, replicate...)
  full_join(sample_info, by = ("condition")) %>% 
  # filter to retain only genes of interest
  filter(gene_id %in% DEGs) %>% ##CHANGE THIS IF Nº CLUSTER DOESN'T FIT
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

trans_cts_cluster <- trans_cts_mean %>% 
  inner_join(gene_cluster, by = "gene_id") 

head(trans_cts_cluster)

##########################################################################
##########################################################################
####################### CLUSTER REPRESENTATION ###########################
##########################################################################
##########################################################################

# Doing representation with an average of the expression

lowtemperature <- trans_cts_cluster[trans_cts_cluster$temperature=="16C",]
highttemperature <- trans_cts_cluster[trans_cts_cluster$temperature=="24C",]

png("clustering_30.png", 2000, 1000, pointsize=30)
ggplot(NULL, aes(x=factor(light, level = c("10uE", "100uE")), y=mean_cts_scaled, group = gene_id)) + 
  geom_line(data=highttemperature, color='orangered1', alpha = 0.02) + 
  geom_line(data=lowtemperature, color='dodgerblue1', alpha = 0.02)+
  #geom_point(data = trans_cts_cluster, aes(x = hour, y = mean_cts_scaled, group = interaction(cluster,NaOH)), colour = "black", cex=1, alpha =0.1) +
  facet_grid(cols = vars(cluster)) + 
  stat_summary(data = highttemperature,aes(group=1), fun=mean, geom="line", colour="orangered1", size = 1) + 
  stat_summary(data = lowtemperature,aes(group=1), fun=mean, geom="line", colour="dodgerblue1", size = 1)  +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap( ~ cluster, ncol = 10) +
  labs(x = "light", y = "Z-score") +
  scale_colour_manual(values=c("orangered1", "dodgerblue1"), 
                      name  ="Temperature", 
                      labels=c("100mM", "0mM"))
dev.off()

# List of DEGs with their cluster

cluster_annotation <- unique(trans_cts_cluster[,c("gene_id", "cluster")])
write.csv(cluster_annotation, file = "cluster_annotation_GRN.csv")

# List of DEGs for each cluster

clusters <- data.frame()
variable1 <- cluster_annotation[cluster_annotation$cluster==1,]$gene_id
variable2 <- cluster_annotation[cluster_annotation$cluster==2,]$gene_id
clusters <- qpcR:::cbind.na(variable1, variable2)
i <- 3
for (i in 3:30) {
  variable <- cluster_annotation[cluster_annotation$cluster==i,]$gene_id
  clusters <- qpcR:::cbind.na(clusters, variable)
}
write.csv(clusters, file = "clusters_genes.csv")

# Calculating the percentage of DEGs that represents the 10 more populated clusters

indx <- tail(names(sort(table(gene_cluster$cluster))),10)
indx
suma <- sum(gene_cluster$cluster==indx[1])
i <- 2
for (i in 2:length(indx)) {
  suma <- suma + sum(gene_cluster$cluster==indx[i])
  porcentaje <- suma/nrow(gene_cluster)
}

# Representation of each cluster one by one

png("cluster_images/cluster1.png", 400, 400, pointsize=30)
ggplot(NULL, aes(x=factor(light, level = c("10uE", "100uE")), y=mean_cts_scaled, group = gene_id)) + 
  geom_line(data=highttemperature[highttemperature$cluster==1,], color='#FDAE61', alpha = 0.05) + 
  geom_line(data=lowtemperature[lowtemperature$cluster==1,], color='#4575B4', alpha = 0.05)+
  #geom_point(data = trans_cts_cluster, aes(x = hour, y = mean_cts_scaled, group = interaction(cluster,NaOH)), colour = "black", cex=1, alpha =0.1) +
  facet_grid(cols = vars(cluster)) + 
  stat_summary(data = highttemperature[highttemperature$cluster==1,],aes(group=1), fun=mean, geom="line", colour="#FDAE61", size = 2) + 
  stat_summary(data = lowtemperature[lowtemperature$cluster==1,],aes(group=1), fun=mean, geom="line", colour="#4575B4", size = 2)  +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=25),
                     axis.title=element_text(size=30),
                     strip.background = element_blank(),
                     strip.text = element_blank()) +
  facet_wrap( ~ cluster, ncol = 10) +
  labs(x = "light", y = "Z-score") +
  scale_colour_manual(values=c("orangered1", "dodgerblue1"), 
                      name  ="Temperature", 
                      labels=c("100mM", "0mM"))
dev.off()
