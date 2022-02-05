# Author: Fernando Rodríguez Marín
# Contact: fernando.rodriguez.marin8@gmail.com
# Project: master's thesis

library(stringr)
library(pheatmap)
library(tidyr)
library(MetBrewer)

##########################################################################
##########################################################################
############################## ORTHOGROUPS ###############################
##########################################################################
##########################################################################

# Reading file with Orthogroups from orthofinder

orthogroups <- read.csv(file = "Orthogroups.csv", sep = ";", header = T)
head(orthogroups)
# Adding column where we can add T or F if a gen from our DEGs is in 
# an orthogroup
orthogroups[,6] <- NA
orthogroups[,7] <- NA
orthogroups[,8] <- NA
orthogroups[,9] <- NA
orthogroups[,10] <- NA
orthogroups[,11] <- NA
orthogroups[,12] <- NA
orthogroups[,13] <- NA
orthogroups[,14] <- NA
orthogroups[,15] <- NA
orthogroups[,16] <- NA
orthogroups[,17] <- NA

colnames(orthogroups) <- c("Orthogroups", "Arabidopsis", "Marchantia",
                           "Mesotaenium", "Ostreococcus",
                           "temperature_Arabidopsis", "light_Arabidopsis", "both_Arabidopsis",
                           "temperature_Marchantia", "light_Marchantia", "both_Marchantia",
                           "temperature_Mesotaenium", "light_Mesotaenium", "both_Mesotaenium",
                           "temperature_Ostreococcus", "light_Ostreococcus", "both_Ostreococcus")
                          # Change it according your table

head(orthogroups)
orthogroups[1,2]

DEG_arabidopsis <- read.csv(file = "download_atcol0/DEGs_Arabidopsis.csv", sep = ";")
DEG_marchantia <- read.csv(file = "download_mptak1/DEGs_Marchantia_Tak1.csv", sep = ";")
DEG_mesotaenium <- read.csv(file = "download_Mesotaenium_v3/DEGs_Mesotaenium.csv", sep = ";")
DEG_ostreococcus <- read.csv(file = "download_Ostreococcus_v2/DEGs_Ostreococcus.csv", sep = ";")


degat_temperature <- na.omit(DEG_arabidopsis$DEGs_temperature_at)
degat_light <- na.omit(DEG_arabidopsis$DEGs_light_at)

degmp_temperature <- na.omit(DEG_marchantia$DEGs_temperature_mp)
degmp_light <- na.omit(DEG_marchantia$DEGs_light_mp)

degme_temperature <- na.omit(DEG_mesotaenium$DEGs_temperature_me)
degme_light <- na.omit(DEG_mesotaenium$DEGs_light_me)

degot_temperature <- na.omit(DEG_ostreococcus$DEGs_temperature_ot)
degot_light <- na.omit(DEG_ostreococcus$DEGs_light_ot)

##########################################################################
##########################################################################
############################## ARABIDOPSIS ###############################
##########################################################################
##########################################################################

# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(degat_light)) {
    if (str_detect(orthogroups$Arabidopsis[j], degat_light[k])) {
      orthogroups$light_Arabidopsis[j] <- str_detect(orthogroups$Arabidopsis[j], degat_light[k])
      break
    }
    else {
      orthogroups$light_Arabidopsis[j] <- FALSE
    }
  }
  k <- 1
}

# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(degat_temperature)) {
    if (str_detect(orthogroups$Arabidopsis[j], degat_temperature[k])) {
      orthogroups$temperature_Arabidopsis[j] <- str_detect(orthogroups$Arabidopsis[j], degat_temperature[k])
      break
    }
    else {
      orthogroups$temperature_Arabidopsis[j] <- FALSE
    }
  }
  k <- 1
}

# Adding TRUE or FALSE to the column: both (temperature & light)

j <- 1

for (j in 1:nrow(orthogroups)) {
    if (orthogroups$temperature_Arabidopsis[j] & 
        orthogroups$light_Arabidopsis[j]) {
      orthogroups$both_Arabidopsis[j] <- TRUE
    }
    else {
      orthogroups$both_Arabidopsis[j] <- FALSE
    }
}

##########################################################################
##########################################################################
############################## MARCHANTIA ################################
##########################################################################
##########################################################################

# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(degmp_light)) {
    if (str_detect(orthogroups$Marchantia[j], degmp_light[k])) {
      orthogroups$light_Marchantia[j] <- str_detect(orthogroups$Marchantia[j], degmp_light[k])
      break
    }
    else {
      orthogroups$light_Marchantia[j] <- FALSE
    }
  }
  k <- 1
}

# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(degmp_temperature)) {
    if (str_detect(orthogroups$Marchantia[j], degmp_temperature[k])) {
      orthogroups$temperature_Marchantia[j] <- str_detect(orthogroups$Marchantia[j], degmp_temperature[k])
      break
    }
    else {
      orthogroups$temperature_Marchantia[j] <- FALSE
    }
  }
  k <- 1
}

# Adding TRUE or FALSE to the column: both (temperature & light)

j <- 1

for (j in 1:nrow(orthogroups)) {
  if (orthogroups$temperature_Marchantia[j] & 
      orthogroups$light_Marchantia[j]) {
    orthogroups$both_Marchantia[j] <- TRUE
  }
  else {
    orthogroups$both_Marchantia[j] <- FALSE
  }
}

##########################################################################
##########################################################################
############################## MESOTAENIUM ###############################
##########################################################################
##########################################################################

# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(degme_light)) {
    if (str_detect(orthogroups$Mesotaenium[j], degme_light[k])) {
      orthogroups$light_Mesotaenium[j] <- str_detect(orthogroups$Mesotaenium[j], degme_light[k])
      break
    }
    else {
      orthogroups$light_Mesotaenium[j] <- FALSE
    }
  }
  k <- 1
}

# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(degme_temperature)) {
    if (str_detect(orthogroups$Mesotaenium[j], degme_temperature[k])) {
      orthogroups$temperature_Mesotaenium[j] <- str_detect(orthogroups$Mesotaenium[j], degme_temperature[k])
      break
    }
    else {
      orthogroups$temperature_Mesotaenium[j] <- FALSE
    }
  }
  k <- 1
}

# Adding TRUE or FALSE to the column: both (temperature & light)

j <- 1

for (j in 1:nrow(orthogroups)) {
  if (orthogroups$temperature_Mesotaenium[j] & 
      orthogroups$light_Mesotaenium[j]) {
    orthogroups$both_Mesotaenium[j] <- TRUE
  }
  else {
    orthogroups$both_Mesotaenium[j] <- FALSE
  }
}

##########################################################################
##########################################################################
############################ OSTREOCOCCUS ################################
##########################################################################
##########################################################################

# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(degot_light)) {
    if (str_detect(orthogroups$Ostreococcus[j], degot_light[k])) {
      orthogroups$light_Ostreococcus[j] <- str_detect(orthogroups$Ostreococcus[j], degot_light[k])
      break
    }
    else {
      orthogroups$light_Ostreococcus[j] <- FALSE
    }
  }
  k <- 1
}

# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(degot_temperature)) {
    if (str_detect(orthogroups$Ostreococcus[j], degot_temperature[k])) {
      orthogroups$temperature_Ostreococcus[j] <- str_detect(orthogroups$Ostreococcus[j], degot_temperature[k])
      break
    }
    else {
      orthogroups$temperature_Ostreococcus[j] <- FALSE
    }
  }
  k <- 1
}

# Adding TRUE or FALSE to the column: both (temperature & light)

j <- 1

for (j in 1:nrow(orthogroups)) {
  if (orthogroups$temperature_Ostreococcus[j] & 
      orthogroups$light_Ostreococcus[j]) {
    orthogroups$both_Ostreococcus[j] <- TRUE
  }
  else {
    orthogroups$both_Ostreococcus[j] <- FALSE
  }
}

write.csv2(orthogroups, file = "Orthogroups_species.csv")

# Deleting orthogroups with only FALSES
keep <- orthogroups[,6] | orthogroups[,7]| orthogroups[,8]| orthogroups[,9]| orthogroups[,10]| orthogroups[,11]| orthogroups[,12]| orthogroups[,13]| orthogroups[,14]| orthogroups[,15]| orthogroups[,16]| orthogroups[,17]
orthogroups <- orthogroups[keep,]

write.csv2(orthogroups, file = "Orthogroups_filtered.csv")

##########################################################################
##########################################################################
################################ HEATMAP #################################
##########################################################################
##########################################################################

# Calculating how many organism responding to temperature, light or both

orthogroups <- read.csv(file = "Orthogroups_filtered.csv", sep = ";", header = T)

# Mantain only booleans columns
orthogroups <- orthogroups[,c(2,7,8,9,10,11,12,13,14,15,16,17,18)]
head(orthogroups)
heatmap_orthogroups <- data.frame()
i <- 1

for (i in 1:length(orthogroups$Orthogroups)) {
  vector_temperature <- c(orthogroups[i,2], orthogroups[i,5], orthogroups[i,8], orthogroups[i,11])
  if (sum(vector_temperature, na.rm = TRUE)>=1) {
    heatmap_orthogroups[i,1] <- orthogroups[,1][i]
    heatmap_orthogroups[i,2] <- sum(vector_temperature, na.rm = TRUE) 
    i <- i + 1
  }
  else { i <- i + 1}
}

i <- 1

for (i in 1:length(orthogroups$Orthogroups)) {
  vector_light <- c(orthogroups[i,3], orthogroups[i,6], orthogroups[i,9], orthogroups[i,12])
  if (sum(vector_light, na.rm = TRUE)>=1) {
    heatmap_orthogroups[i,1] <- orthogroups[,1][i]
    heatmap_orthogroups[i,3] <- sum(vector_light, na.rm = TRUE) 
    i <- i + 1
  }
  else { i <- i + 1}
}

i <- 1

for (i in 1:length(orthogroups$Orthogroups)) {
  vector_both <- c(orthogroups[i,4], orthogroups[i,7], orthogroups[i,10], orthogroups[i,13])
  if (sum(vector_both, na.rm = TRUE)>=1) {
    heatmap_orthogroups[i,1] <- orthogroups[,1][i]
    heatmap_orthogroups[i,4] <- sum(vector_both, na.rm = TRUE) 
    i <- i + 1
  }
  else { i <- i + 1}
}

head(heatmap_orthogroups)

heatmap_orthogroups <- heatmap_orthogroups %>% drop_na(V1)

heatmap_orthogroups[is.na(heatmap_orthogroups)] <- 0

rownames(heatmap_orthogroups) <- heatmap_orthogroups[,1]
heatmap_orthogroups <- heatmap_orthogroups[,2:4]
colnames(heatmap_orthogroups) <- c("temperature", "light", "temperature & light")

heatmap_orthogroups <- heatmap_orthogroups[1:40,]
paletteLength <- 4
myColor <- colorRampPalette(c("white", "darkgreen"))(paletteLength)

pheatmap(heatmap_orthogroups, color = myColor, show_rownames =F, cluster_rows = T, cluster_cols = F)

# Selecting the orthogroups where there are DEGs for the 2 conditions in 
# at least 3 species

i <- 1
both_2_GO <- c()
for (i in 1:nrow(heatmap_orthogroups)) {
  if (((heatmap_orthogroups[i,3] >= 3))) {
    both_2_GO <- append(both_2_GO, rownames(heatmap_orthogroups)[i])
  }
}

heatmap2 <- heatmap_orthogroups[both_2_GO,]

pheatmap(heatmap2, color = c("palegreen", "palegreen4"), show_rownames =T, cluster_rows = T, cluster_cols = F)

##########################################################################
##########################################################################
################################ HEATMAP 2 ###############################
##########################################################################
##########################################################################

# Calculating how many organism responding to temperature, light or both

orthogroups <- read.csv(file = "Orthogroups_filtered.csv", sep = ";", header = T)

i <- 1
heatmap2_orthologs <- data.frame(orthogroups$Orthogroups)
for (i in 1:nrow(orthogroups)) {
if (orthogroups$temperature_Arabidopsis[i]) {
  heatmap2_orthologs[i,2] <- c("Arabidopsis")
}
}
i <- 1
for (i in 1:nrow(orthogroups)) {
  if (orthogroups$temperature_Marchantia[i]) {
    heatmap2_orthologs[i,2] <- paste(heatmap2_orthologs[i,2], "Marchantia")
  }
}
i <- 1
for (i in 1:nrow(orthogroups)) {
  if (orthogroups$temperature_Mesotaenium[i]) {
    heatmap2_orthologs[i,2] <- paste(heatmap2_orthologs[i,2], "Mesotaenium")
  }
}
i <- 1
for (i in 1:nrow(orthogroups)) {
  if (orthogroups$temperature_Ostreococcus[i]) {
    heatmap2_orthologs[i,2] <- paste(heatmap2_orthologs[i,2], "Ostreococcus")
  }
}

# Light

i <- 1
for (i in 1:nrow(orthogroups)) {
  if (orthogroups$light_Arabidopsis[i]) {
    heatmap2_orthologs[i,3] <- c("Arabidopsis")
  }
}
i <- 1
for (i in 1:nrow(orthogroups)) {
  if (orthogroups$light_Marchantia[i]) {
    heatmap2_orthologs[i,3] <- paste(heatmap2_orthologs[i,3], "Marchantia")
  }
}
i <- 1
for (i in 1:nrow(orthogroups)) {
  if (orthogroups$light_Mesotaenium[i]) {
    heatmap2_orthologs[i,3] <- paste(heatmap2_orthologs[i,3], "Mesotaenium")
  }
}
i <- 1
for (i in 1:nrow(orthogroups)) {
  if (orthogroups$light_Ostreococcus[i]) {
    heatmap2_orthologs[i,3] <- paste(heatmap2_orthologs[i,3], "Ostreococcus")
  }
}

# Both

i <- 1
for (i in 1:nrow(orthogroups)) {
  if (orthogroups$both_Arabidopsis[i]) {
    heatmap2_orthologs[i,4] <- c("Arabidopsis")
  }
}
i <- 1
for (i in 1:nrow(orthogroups)) {
  if (orthogroups$both_Marchantia[i]) {
    heatmap2_orthologs[i,4] <- paste(heatmap2_orthologs[i,4], "Marchantia")
  }
}
i <- 1
for (i in 1:nrow(orthogroups)) {
  if (orthogroups$both_Mesotaenium[i]) {
    heatmap2_orthologs[i,4] <- paste(heatmap2_orthologs[i,4], "Mesotaenium")
  }
}
i <- 1
for (i in 1:nrow(orthogroups)) {
  if (orthogroups$both_Ostreococcus[i]) {
    heatmap2_orthologs[i,4] <- paste(heatmap2_orthologs[i,4], "Ostreococcus")
  }
}

heatmap2_orthologs[] <- lapply(heatmap2_orthologs, gsub, pattern="NA", replacement='')
heatmap2_orthologs[is.na(heatmap2_orthologs)] <- c("None")

colnames(heatmap2_orthologs) <- c("Orthogroups", "temperature", "light", "both")

write.csv(heatmap2_orthologs, file = "Orthogroups_with_each_specie.csv")

# Reading file

heatmap2_orthologs <- read.csv("Orthogroups_with_each_specie.csv")
rownames(heatmap2_orthologs) <- heatmap2_orthologs[,2]

i <- 1
both_2_GO <- c()
for (i in 1:nrow(heatmap_orthogroups)) {
  if (((heatmap_orthogroups[i,1] >= 2 &heatmap_orthogroups[i,2] >= 2&heatmap_orthogroups[i,3] >= 2))) {
    both_2_GO <- append(both_2_GO, rownames(heatmap_orthogroups)[i])
  }
}

# Asign number to each combination

heatmap_filtered <- heatmap2_orthologs[both_2_GO,][,3:5]
head(heatmap_filtered)
heatmap_filtered[heatmap_filtered=="Arabidopsis Marchantia"] <- 6
heatmap_filtered[heatmap_filtered=="Arabidopsis Mesotaenium"] <- 5
heatmap_filtered[heatmap_filtered=="Arabidopsis Ostreococcus"] <- 3
heatmap_filtered[heatmap_filtered==" Marchantia Mesotaenium"] <- 4
heatmap_filtered[heatmap_filtered==" Marchantia Ostreococcus"] <- 2
heatmap_filtered[heatmap_filtered==" Mesotaenium Ostreococcus"] <- 1
heatmap_filtered[heatmap_filtered=="Arabidopsis Marchantia Mesotaenium"] <- 10
heatmap_filtered[heatmap_filtered=="Arabidopsis Marchantia Ostreococcus"] <- 9
heatmap_filtered[heatmap_filtered=="Arabidopsis Mesotaenium Ostreococcus"] <- 8
heatmap_filtered[heatmap_filtered==" Marchantia Mesotaenium Ostreococcus"] <- 7
heatmap_filtered[heatmap_filtered=="Arabidopsis Marchantia Mesotaenium Ostreococcus"] <- 11

heatmap_filtered_v2 <- data.frame(apply(heatmap_filtered, 2, function(x) as.numeric(as.character(x))))
rownames(heatmap_filtered_v2) <- rownames(heatmap_filtered)

#met.brewer("Robert",n=15,type="continuous")
color <- c("#041D2C", "#0D3E61", "#296FA2", "#053C29", "#2F6856", "#487952", "#556219", "#AB8929", "#8F6527", "#6D481B", "#472C0B")
color <- c("#A20000", "#FF0000", "#FF3300","#FF6600",
           "#FFC000", "#92D050", "#00B050",
           "#9DC3E6", "#0070C0", "#002060", "#7030A0")
color <- c("#FF4D4D", "#7C9DD6", "#9BC67F","#9B6FBD",
           "#FFAA4C", "#FFFF4D", "#D78C59",
           "#FF94B8", "#CBD8EE", "#4D4D4D", "#B9FFE8")


pheatmap(heatmap_filtered_v2, show_rownames =F, cluster_rows = T, cluster_cols = F,
         color = color, width = 100)

png("Heatmap_orthogroups_coloured_rainbow2.png", 1000, 4000, pointsize=30)
pheatmap(heatmap_filtered_v2, show_rownames =F, cluster_rows = T, cluster_cols = F,
         color = color)
dev.off()

write.csv(heatmap_filtered_v2, file = "Orthogroups_table_figure.csv")


##########################################################################
##########################################################################
################################ DOTPLOT #################################
##########################################################################
##########################################################################


tabla <- read.csv(file = "Dotplot_orthogroups.csv", sep = ";", header = T)

library("reshape2")
library(ggplot2)
#convert data frame from a "wide" format to a "long" format
pcm = melt(tabla, id = c("X"))
pcm$X <- factor(pcm$X,levels=unique(pcm$X))


color <- c("#FF4D4D", "#7C9DD6", "#9BC67F","#9B6FBD",
           "#FFAA4C", "#FFFF4D", "#D78C59",
           "#FF94B8", "#CBD8EE", "#4D4D4D", "#B9FFE8")

tiff("Plot2.tif", res = 50)

plot <- ggplot(pcm, aes(x = X, y = variable)) +
  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) +
  scale_size_continuous(limits = c(1, 105), range = c(1,20), breaks = c(1,15,30,45,60,75,90,105)) +
  labs( x= "", y = "", size =  "Abundancia", fill = "")  +
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, angle = 90), 
        axis.text.y = element_text(colour = "black", size = 7, face="italic"), 
        legend.text = element_text(size = 7, face ="bold", colour ="black" ), 
        legend.title = element_text(size = 10, face = "bold"), 
        panel.background = element_blank(), 
        legend.position = "left",
        rect = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = color)  +
  scale_y_discrete(limits = rev(levels(pcm$variable))) 

ggsave(
  "ggtest.png",
  plot,
  width = 4,
  height = 7,
  dpi = 1000
)

##########################################################################
##########################################################################
############################## GO ORTHOGROUPS ############################
##########################################################################
##########################################################################

orthogroups <- read.csv(file = "Orthogroups.csv", sep = ";", header = T)

DEGs_cluster <- read.csv("Mptak1/clusters_genes.csv")
DEGs_cluster <- na.omit(DEGs_cluster$variable28)


j<-1
k<-1
l<-1

DEGs_orthogroups <- c()
for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEGs_cluster)) {
    if (str_detect(orthogroups$Marchantia_polymorpha[j], DEGs_cluster[k])) {
      DEGs_orthogroups[l] <- orthogroups$Arabidopsis_thaliana[j]
      l <- l+1
      break
    }
    else {
      k <- k+1
    }
  }
  k <- 1
}

DEGs_orthogroups <- data.frame(DEGs_orthogroups)
nrow(DEGs_orthogroups)
library(dplyr)
DEGs_orthogroups <- DEGs_orthogroups %>% 
  mutate(DEGs_orthogroups = strsplit(as.character(DEGs_orthogroups), ",")) %>%
  unnest(DEGs_orthogroups)
DEGs_orthogroups$DEGs_orthogroups <- substr(DEGs_orthogroups$DEGs_orthogroups, 0, 12)

write.csv(DEGs_orthogroups, "cluster28_Marchantia_to_Arabidopsis.csv")
DEGs_cluster <- read.csv("cluster28_Marchantia_to_Arabidopsis.csv")
