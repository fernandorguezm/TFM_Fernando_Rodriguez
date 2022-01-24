# Author: Fran Romero Campero
# Edited by: Fernando Rodríguez Marín
# Contact: fernando.rodriguez.marin8@gmail.com
# Project: master's thesis

library(igraph)
library(cluster)
library(ggpubr)
library(rstatix)
library(tidyverse)

##########################################################################
##########################################################################
##################### CO-EXPRESSION NETWORKS #############################
##########################################################################
##########################################################################

DEGs <- read.csv(file = "DEGs.csv", sep = ";")
DEGs_all <- na.omit(DEGs$DEGs_all)

mean_TPM <- read.csv("Mean_TPM.csv", sep = ";")
rownames(mean_TPM) <- mean_TPM[,1]
mean_TPM <- mean_TPM[,2:5]

# Correlation analysis between expression levels. Calculating correlation matrixs between DEGs.

gene.correlation <- cor(t(mean_TPM[,1:2][DEGs_all,]))
dim(gene.correlation)
gene.correlation[1:5, 1:5]

# Choosing a threshold that produces a good aproach to a free scale network.
# All genes have (in average) more than 10 neighbors 

thresholds <- seq(from=0.70,to=0.99,by=0.01)
mean.connectivities <- vector(length=length(thresholds))
scale.free.R2 <- vector(length=length(thresholds))

## Go through all the correlation thresholds
for(i in 1:length(thresholds))
{
  print(thresholds[i])
  ## Construct the network corresponding to the current threshold being checked
  ## Adjacency matrix
  current.adjacency <- (gene.correlation > thresholds[i] & gene.correlation < 1) 
  ## Network
  threshold.network <- graph.adjacency(current.adjacency, mode="undirected")
  
  ## Compute node degrees
  node.degrees <- degree(threshold.network)
  
  ## Keep track of the mean connectivity or mean node degree of the current network
  mean.connectivities[i] <- mean(node.degrees)
  
  ## Check scale free property
  h <- hist(node.degrees)
  ## Compute degree frequencies
  degree.frequencies <- table(node.degrees)
  ## Determine linear regression for logarithmic transformed degree frequencies
  lm.r <- lm(log10(as.numeric(names(degree.frequencies[-1]))) ~ log10(degree.frequencies[-1]))
  ## Extract R squared as a measure of adjustment to the scale free property
  s.lm <- summary(lm.r)
  scale.free.R2[i] <- s.lm[["adj.r.squared"]]
}

plot(thresholds,mean.connectivities,type="o",col="red",lwd=3,xlab="Correlation Threshold",ylab="Mean connectivity")
plot(thresholds,scale.free.R2,type="o",col="blue",lwd=3,xlim=c(0.70,0.99),xlab="Correlation Threshold",ylab="Scale Free Model Fit (R²)")

names(mean.connectivities) <- thresholds
names(scale.free.R2) <- thresholds

mean.connectivities
scale.free.R2

# Choose the threshold with a bigger Scale free model fit (R2)

adjacency <- (gene.correlation > 0.99) & (gene.correlation < 1)
gene.coexpression.network <- graph.adjacency(adjacency, mode="undirected")
write.graph(gene.coexpression.network,file="_GCN.gml",format="gml")

######################################################################
######################################################################
######################################################################

# To see how many nodes have your network: 
gene.coexpression.network

## To read .gml file
gene.coexpression.network <- read.graph(file="_gene_coexpression_network.gml",format="gml")

# Anotattion of DEGs regulated by light, temperature or both

light_temperature_both_data <- data.frame(DEGs_all)

i <- 1
for (i in 1:length(DEGs_all)) {
  
if (DEGs_all[i] %in% DEGs$DEGs_temperature) {
  light_temperature_both_data$data[light_temperature_both_data[,1] == DEGs_all[i]] <- "temperature"
} 
if (DEGs_all[i] %in% DEGs$DEGs_light) {
  light_temperature_both_data$data[light_temperature_both_data[,1] == DEGs_all[i]] <- "light"
}
if (DEGs_all[i] %in% DEGs$DEGs_temperature & DEGs_all[i] %in% DEGs$DEGs_light) {
    light_temperature_both_data$data[light_temperature_both_data[,1] == DEGs_all[i]] <- "both"
} 
}

write.table(light_temperature_both_data, file = "info_network.txt")

write.csv(mean_TPM[DEGs_all,], file = "mean_expression_DEGs.csv")
