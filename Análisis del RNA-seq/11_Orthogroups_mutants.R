# Author: Fernando Rodríguez Marín
# Contact: fernando.rodriguez.marin8@gmail.com
# Project: master's thesis

library(stringr)
library(pheatmap)
library(tidyr)
library(MetBrewer)
library(VennDiagram)

##########################################################################
##########################################################################
############################## ORTHOGROUPS ###############################
##########################################################################
##########################################################################

# Reading file with Orthogroups from orthofinder

orthogroups <- read.csv(file = "Orthogroups_at_mp_common.csv", sep = ";", header = T)[1:3]
head(orthogroups)
# Adding column where we can add T or F if a gen from our DEGs is in 
# an orthogroup
orthogroups[,4] <- NA
orthogroups[,5] <- NA
orthogroups[,6] <- NA
orthogroups[,7] <- NA
orthogroups[,8] <- NA
orthogroups[,9] <- NA

colnames(orthogroups) <- c("Orthogroups", 
                           "Arabidopsis", "Marchantia",
                           "Athy5", "Mphy5",
                           "Atphyb", "Mpphy",
                           "Atpif4", "Mppif")
                          # Change it according your table

head(orthogroups)
orthogroups[1,2]

DEG_athy5 <- read.csv(file = "Athy5_vs_Atcol0/DEGs_Arabidopsis.csv", sep = ";")$DEGs_all_at
DEG_mphy5 <- read.csv(file = "Mphy5_vs_Mptak1/DEGs_Marchantia.csv", sep = ";")$DEGs_all_mp
DEG_atphyb <- read.csv(file = "Atphyb_vsAtcol0/DEGs_Arabidopsis.csv", sep = ";")$DEGs_all_at
DEG_mpphy <- read.csv(file = "Mpphy_vs_Mptak1/DEGs_Marchantia.csv", sep = ";")$DEGs_all_mp
DEG_atpif4 <- read.csv(file = "Atpif4_vs_Atcol0/DEGs_Arabidopsis.csv", sep = ";")$DEGs_all_at
DEG_mppif <- read.csv(file = "Mppif_vs_Mptak1/DEGs_Marchantia.csv", sep = ";")$DEGs_all_mp

# Athy5
# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEG_athy5)) {
    if (str_detect(orthogroups$Arabidopsis[j], DEG_athy5[k])) {
      orthogroups$Athy5[j] <- str_detect(orthogroups$Arabidopsis[j], DEG_athy5[k])
      break
    }
    else {
      orthogroups$Athy5[j] <- FALSE
    }
  }
  k <- 1
}

# Atphyb
# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEG_atphyb)) {
    if (str_detect(orthogroups$Arabidopsis[j], DEG_atphyb[k])) {
      orthogroups$Atphyb[j] <- str_detect(orthogroups$Arabidopsis[j], DEG_atphyb[k])
      break
    }
    else {
      orthogroups$Atphyb[j] <- FALSE
    }
  }
  k <- 1
}

# Atpif4
# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEG_atpif4)) {
    if (str_detect(orthogroups$Arabidopsis[j], DEG_atpif4[k])) {
      orthogroups$Atpif4[j] <- str_detect(orthogroups$Arabidopsis[j], DEG_atpif4[k])
      break
    }
    else {
      orthogroups$Atpif4[j] <- FALSE
    }
  }
  k <- 1
}


# Mphy5
# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEG_mphy5)) {
    if (str_detect(orthogroups$Marchantia[j], DEG_mphy5[k])) {
      orthogroups$Mphy5[j] <- str_detect(orthogroups$Marchantia[j], DEG_mphy5[k])
      break
    }
    else {
      orthogroups$Mphy5[j] <- FALSE
    }
  }
  k <- 1
}
# Mpphy
# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEG_mpphy)) {
    if (str_detect(orthogroups$Marchantia[j], DEG_mpphy[k])) {
      orthogroups$Mpphy[j] <- str_detect(orthogroups$Marchantia[j], DEG_mpphy[k])
      break
    }
    else {
      orthogroups$Mpphy[j] <- FALSE
    }
  }
  k <- 1
}

# Mppif
# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEG_mppif)) {
    if (str_detect(orthogroups$Marchantia[j], DEG_mppif[k])) {
      orthogroups$Mppif[j] <- str_detect(orthogroups$Marchantia[j], DEG_mppif[k])
      break
    }
    else {
      orthogroups$Mppif[j] <- FALSE
    }
  }
  k <- 1
}

write.csv(orthogroups, file = "Orthogroups_mutants.csv")

orthogroups <- read.csv("Orthogroups_mutants.csv")
head(orthogroups)

orthogroups[,11] <- NA
orthogroups[,12] <- NA
orthogroups[,13] <- NA


# Mantain only booleans columns
orthogroups_mutants <- orthogroups[,c(2,5,6,7,8,9,10,11,12,13)]
head(orthogroups_mutants)


# HY5

i <- 1
for (i in 1:length(orthogroups_mutants$Orthogroups)) {
  vector_hy5 <- c(orthogroups_mutants[i,2], orthogroups_mutants[i,3])
  if (sum(vector_hy5, na.rm = TRUE)==2) {
    orthogroups_mutants[i,8] <- TRUE
    i <- i + 1
  }
  else { 
    orthogroups_mutants[i,8] <- FALSE
    i <- i + 1}
}

# PHY

i <- 1
for (i in 1:length(orthogroups_mutants$Orthogroups)) {
  vector_hy5 <- c(orthogroups_mutants[i,4], orthogroups_mutants[i,5])
  if (sum(vector_hy5, na.rm = TRUE)==2) {
    orthogroups_mutants[i,9] <- TRUE
    i <- i + 1
  }
  else { 
    orthogroups_mutants[i,9] <- FALSE
    i <- i + 1}
}

# PIF

i <- 1
for (i in 1:length(orthogroups_mutants$Orthogroups)) {
  vector_hy5 <- c(orthogroups_mutants[i,6], orthogroups_mutants[i,7])
  if (sum(vector_hy5, na.rm = TRUE)==2) {
    orthogroups_mutants[i,10] <- TRUE
    i <- i + 1
  }
  else { 
    orthogroups_mutants[i,10] <- FALSE
    i <- i + 1}
}

write.csv(orthogroups_mutants, file = "Orthogroups_mutants_comparations.csv")

library(VennDiagram)

orthogroups_mutants <- read.csv("Orthogroups_mutants_comparations.csv")
head(orthogroups_mutants)

sum(orthogroups_mutants$Athy5, na.rm = TRUE)
sum(orthogroups_mutants$Mphy5, na.rm = TRUE)
sum(orthogroups_mutants$HY5, na.rm = TRUE)

png("Venn_HY5_AtMp.png", 1000, 1000, pointsize=50)
draw.pairwise.venn(area1 = 1352, area2 = 246, cross.area = 163, 
                   category = c("", ""),
                  col=c("#A5C68F","#BFE498"),
                  fill = c("#A5C68F","#BFE498"),
                  alpha = 0.4,
                  lwd =5)
dev.off()

VennDiag <- euler(c("A" = 1189, "B" = 83, "A&B" = 163))

png("Venn_HY5_AtMp.png", 1000, 1000, pointsize=100)
plot(VennDiag, fontsize=1000, cex=1000, alpha=0.9,
     fill=c("#A5C68F","#BFE498"),
     edges = list(col = c("black","black", "black"), lex = 2), 
     quantities = list(cex = 4),
     labels = c("", "", ""))
dev.off()


sum(orthogroups_mutants$Atphyb, na.rm = TRUE)
sum(orthogroups_mutants$Mpphy, na.rm = TRUE)
sum(orthogroups_mutants$PHY, na.rm = TRUE)

png("Venn_PHY_AtMp.png", 1000, 1000, pointsize=50)
draw.pairwise.venn(area1 = 771, area2 = 948, cross.area = 315, 
                   category = c("", ""),
                   col=c("#A5C68F","#BFE498"),
                   fill = c("#A5C68F","#BFE498"),
                   alpha = 0.4,
                   lwd =5)
dev.off()

VennDiag <- euler(c("A" = 456, "B" = 633, "A&B" = 315))

png("Venn_PHY_AtMp.png", 1000, 1000, pointsize=100)
plot(VennDiag, fontsize=1000, cex=1000, alpha=0.9,
     fill=c("#A5C68F","#BFE498"),
     edges = list(col = c("black","black", "black"), lex = 2), 
     quantities = list(cex = 4),
     labels = c("", "", ""))
dev.off()


sum(orthogroups_mutants$Atpif4, na.rm = TRUE)
sum(orthogroups_mutants$Mppif, na.rm = TRUE)
sum(orthogroups_mutants$PIF, na.rm = TRUE)

png("Venn_PIF_AtMp.png", 1000, 1000, pointsize=50)
draw.pairwise.venn(area1 = 357, area2 = 427, cross.area = 109, 
                   category = c("", ""),
                   col=c("#A5C68F","#BFE498"),
                   fill = c("#A5C68F","#BFE498"),
                   alpha = 0.4,
                   lwd =5)
dev.off()

VennDiag <- euler(c("A" = 248, "B" = 318, "A&B" = 109))

png("Venn_PIF_AtMp.png", 1000, 1000, pointsize=100)
plot(VennDiag, fontsize=1000, cex=1000, alpha=0.9,
     fill=c("#A5C68F","#BFE498"),
     edges = list(col = c("black","black", "black"), lex = 2), 
     quantities = list(cex = 4),
     labels = c("", "", ""))
dev.off()

##########################################################################
##########################################################################
############################## COMBINATIONS ##############################
##########################################################################
##########################################################################

orthogroups <- read.csv("Orthogroups_mutants.csv")
head(orthogroups)

orthogroups[,11] <- NA
orthogroups[,12] <- NA
orthogroups[,13] <- NA
orthogroups[,14] <- NA
orthogroups[,15] <- NA
orthogroups[,16] <- NA
orthogroups[,17] <- NA
orthogroups[,18] <- NA
orthogroups[,19] <- NA
orthogroups[,20] <- NA
orthogroups[,21] <- NA
orthogroups[,22] <- NA
orthogroups[,23] <- NA
orthogroups[,24] <- NA

colnames(orthogroups) <- c("X","Orthogroups", 
                           "Arabidopsis", "Marchantia",
                           "Athy5", "Mphy5",
                           "Atphyb", "Mpphy",
                           "Atpif4", "Mppif",
                           "Athy5_only", "Mphy5_only",
                           "Atphyb_only", "Mpphy_only",
                           "Atpif4_only", "Mppif_only",
                           "Athy5_Atphy", "Mphy5_Mpphy",
                           "Athy5_Atpif4", "Mphy5_Mppif",
                           "Atphyb_Atpif4", "Mpphy_Mppif",
                           "Athy5_Atphyb_Atpif4", "Mphy5_Mpphy_Mppif")

DEG_athy5_only <- na.omit(read.csv(file = "DEGs_venny_mutants_at.csv", sep = ";")$Athy5)
DEG_mphy5_only <- na.omit(read.csv(file = "DEGs_venny_mutants_mp.csv", sep = ";")$Mphy5)
DEG_atphy_only <- na.omit(read.csv(file = "DEGs_venny_mutants_at.csv", sep = ";")$Atphyb)
DEG_mpphy_only <- na.omit(read.csv(file = "DEGs_venny_mutants_mp.csv", sep = ";")$Mpphy)
DEG_atpif_only <- na.omit(read.csv(file = "DEGs_venny_mutants_at.csv", sep = ";")$Atpif4)
DEG_mppif_only <- na.omit(read.csv(file = "DEGs_venny_mutants_mp.csv", sep = ";")$Mppif)

DEG_athy5_atphy <- na.omit(read.csv(file = "DEGs_venny_mutants_at.csv", sep = ";")$Athy5_Atphyb)
DEG_mphy5_mpphy <- na.omit(read.csv(file = "DEGs_venny_mutants_mp.csv", sep = ";")$Mphy5_Mpphy)
DEG_athy5_atpif <- na.omit(read.csv(file = "DEGs_venny_mutants_at.csv", sep = ";")$Athy5_Atpif4)
DEG_mphy5_mppif <- na.omit(read.csv(file = "DEGs_venny_mutants_mp.csv", sep = ";")$Mphy5_Mppif)
DEG_atphy_atpif <- na.omit(read.csv(file = "DEGs_venny_mutants_at.csv", sep = ";")$Atphyb_Atpif4)
DEG_mpphy_mppif <- na.omit(read.csv(file = "DEGs_venny_mutants_mp.csv", sep = ";")$Mpphy_Mppif)

DEG_athy5_atphy_atpif <- na.omit(read.csv(file = "DEGs_venny_mutants_at.csv", sep = ";")$Athy5_Atphyb_Atpif4)
DEG_mphy5_mpphy_mppif <- na.omit(read.csv(file = "DEGs_venny_mutants_mp.csv", sep = ";")$Mphy5_Mpphy_Mppif)


# Athy5 only
# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEG_athy5_only)) {
    if (str_detect(orthogroups$Arabidopsis[j], DEG_athy5_only[k])) {
      orthogroups$Athy5_only[j] <- str_detect(orthogroups$Arabidopsis[j], DEG_athy5_only[k])
      break
    }
    else {
      orthogroups$Athy5_only[j] <- FALSE
    }
  }
  k <- 1
}

# Mphy5 only
# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEG_mphy5_only)) {
    if (str_detect(orthogroups$Marchantia[j], DEG_mphy5_only[k])) {
      orthogroups$Mphy5_only[j] <- str_detect(orthogroups$Marchantia[j], DEG_mphy5_only[k])
      break
    }
    else {
      orthogroups$Mphy5_only[j] <- FALSE
    }
  }
  k <- 1
}

# Atphyb only
# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEG_atphy_only)) {
    if (str_detect(orthogroups$Arabidopsis[j], DEG_atphy_only[k])) {
      orthogroups$Atphyb_only[j] <- str_detect(orthogroups$Arabidopsis[j], DEG_atphy_only[k])
      break
    }
    else {
      orthogroups$Atphyb_only[j] <- FALSE
    }
  }
  k <- 1
}

# Mpphy only
# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEG_mpphy_only)) {
    if (str_detect(orthogroups$Marchantia[j], DEG_mpphy_only[k])) {
      orthogroups$Mpphy_only[j] <- str_detect(orthogroups$Marchantia[j], DEG_mpphy_only[k])
      break
    }
    else {
      orthogroups$Mpphy_only[j] <- FALSE
    }
  }
  k <- 1
}

# Atpif4 only
# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEG_atpif_only)) {
    if (str_detect(orthogroups$Arabidopsis[j], DEG_atpif_only[k])) {
      orthogroups$Atpif4_only[j] <- str_detect(orthogroups$Arabidopsis[j], DEG_atpif_only[k])
      break
    }
    else {
      orthogroups$Atpif4_only[j] <- FALSE
    }
  }
  k <- 1
}

# Mppif only
# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEG_mppif_only)) {
    if (str_detect(orthogroups$Marchantia[j], DEG_mppif_only[k])) {
      orthogroups$Mppif_only[j] <- str_detect(orthogroups$Marchantia[j], DEG_mppif_only[k])
      break
    }
    else {
      orthogroups$Mppif_only[j] <- FALSE
    }
  }
  k <- 1
}

###################### 
## COMBINATION OF 2 ##
######################

# Athy5_Atphyb
# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEG_athy5_atphy)) {
    if (str_detect(orthogroups$Arabidopsis[j], DEG_athy5_atphy[k])) {
      orthogroups$Athy5_Atphy[j] <- str_detect(orthogroups$Arabidopsis[j], DEG_athy5_atphy[k])
      break
    }
    else {
      orthogroups$Athy5_Atphy[j] <- FALSE
    }
  }
  k <- 1
}

# Mphy5_Mpphy
# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEG_mphy5_mpphy)) {
    if (str_detect(orthogroups$Marchantia[j], DEG_mphy5_mpphy[k])) {
      orthogroups$Mphy5_Mpphy[j] <- str_detect(orthogroups$Marchantia[j], DEG_mphy5_mpphy[k])
      break
    }
    else {
      orthogroups$Mphy5_Mpphy[j] <- FALSE
    }
  }
  k <- 1
}

# Athy5_Atpif4
# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEG_athy5_atpif)) {
    if (str_detect(orthogroups$Arabidopsis[j], DEG_athy5_atpif[k])) {
      orthogroups$Athy5_Atpif4[j] <- str_detect(orthogroups$Arabidopsis[j], DEG_athy5_atpif[k])
      break
    }
    else {
      orthogroups$Athy5_Atpif4[j] <- FALSE
    }
  }
  k <- 1
}

# Mphy5_Mppif
# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEG_mphy5_mppif)) {
    if (str_detect(orthogroups$Marchantia[j], DEG_mphy5_mppif[k])) {
      orthogroups$Mphy5_Mppif[j] <- str_detect(orthogroups$Marchantia[j], DEG_mphy5_mppif[k])
      break
    }
    else {
      orthogroups$Mphy5_Mppif[j] <- FALSE
    }
  }
  k <- 1
}

# Atphyb_Atpif4
# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEG_atphy_atpif)) {
    if (str_detect(orthogroups$Arabidopsis[j], DEG_atphy_atpif[k])) {
      orthogroups$Atphyb_Atpif4[j] <- str_detect(orthogroups$Arabidopsis[j], DEG_atphy_atpif[k])
      break
    }
    else {
      orthogroups$Atphyb_Atpif4[j] <- FALSE
    }
  }
  k <- 1
}

# Mpphy_Mpphy
# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEG_mpphy_mppif)) {
    if (str_detect(orthogroups$Marchantia[j], DEG_mpphy_mppif[k])) {
      orthogroups$Mpphy_Mppif[j] <- str_detect(orthogroups$Marchantia[j], DEG_mpphy_mppif[k])
      break
    }
    else {
      orthogroups$Mpphy_Mppif[j] <- FALSE
    }
  }
  k <- 1
}


###################### 
## COMBINATION OF 3 ##
######################

# Athy5_Atphyb_Atpif4
# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEG_athy5_atphy_atpif)) {
    if (str_detect(orthogroups$Arabidopsis[j], DEG_athy5_atphy_atpif[k])) {
      orthogroups$Athy5_Atphyb_Atpif4[j] <- str_detect(orthogroups$Arabidopsis[j], DEG_athy5_atphy_atpif[k])
      break
    }
    else {
      orthogroups$Athy5_Atphyb_Atpif4[j] <- FALSE
    }
  }
  k <- 1
}

# Mphy5_Mpphy_Mppif
# Loop to add TRUE if a gen of the DEGs list is in a specific orthogroup

j<-1
k<-1

for (j in 1:nrow(orthogroups)) {
  for (k in 1:length(DEG_mphy5_mpphy_mppif)) {
    if (str_detect(orthogroups$Marchantia[j], DEG_mphy5_mpphy_mppif[k])) {
      orthogroups$Mphy5_Mpphy_Mppif[j] <- str_detect(orthogroups$Marchantia[j], DEG_mphy5_mpphy_mppif[k])
      break
    }
    else {
      orthogroups$Mphy5_Mpphy_Mppif[j] <- FALSE
    }
  }
  k <- 1
}


write.csv(orthogroups, file = "Orthogroups_mutants_venny.csv")
orthogroups_mutants <- read.csv(file = "Orthogroups_mutants_venny.csv")

# HY5-PHY

i <- 1
for (i in 1:length(orthogroups_mutants$Orthogroups)) {
  vector_hy5_phy <- c(orthogroups_mutants[i,18], orthogroups_mutants[i,19])
  if (sum(vector_hy5_phy, na.rm = TRUE)==2) {
    orthogroups_mutants[i,26] <- TRUE
    i <- i + 1
  }
  else { 
    orthogroups_mutants[i,26] <- FALSE
    i <- i + 1}
}

# HY5-PIF

i <- 1
for (i in 1:length(orthogroups_mutants$Orthogroups)) {
  vector_hy5_pif <- c(orthogroups_mutants[i,20], orthogroups_mutants[i,21])
  if (sum(vector_hy5_pif, na.rm = TRUE)==2) {
    orthogroups_mutants[i,27] <- TRUE
    i <- i + 1
  }
  else { 
    orthogroups_mutants[i,27] <- FALSE
    i <- i + 1}
}

# PHY-PIF

i <- 1
for (i in 1:length(orthogroups_mutants$Orthogroups)) {
  vector_phy_pif <- c(orthogroups_mutants[i,22], orthogroups_mutants[i,23])
  if (sum(vector_phy_pif, na.rm = TRUE)==2) {
    orthogroups_mutants[i,28] <- TRUE
    i <- i + 1
  }
  else { 
    orthogroups_mutants[i,28] <- FALSE
    i <- i + 1}
}

# HY5-PHY-PIF

i <- 1
for (i in 1:length(orthogroups_mutants$Orthogroups)) {
  vector_hy5_phy_pif <- c(orthogroups_mutants[i,24], orthogroups_mutants[i,25])
  if (sum(vector_hy5_phy_pif, na.rm = TRUE)==2) {
    orthogroups_mutants[i,29] <- TRUE
    i <- i + 1
  }
  else { 
    orthogroups_mutants[i,29] <- FALSE
    i <- i + 1}
}


write.csv(orthogroups_mutants, file = "Orthogroups_mutants_venny.csv")
orthogroups_mutants <- read.csv(file = "Orthogroups_mutants_venny.csv", sep = ";")

sum(orthogroups_mutants$Athy5_Atphy, na.rm = TRUE)
sum(orthogroups_mutants$Mphy5_Mpphy, na.rm = TRUE)
sum(orthogroups_mutants$hy5_phy, na.rm = TRUE)

png("Venn_HY5_PHY.png", 1000, 1000, pointsize=50)
draw.pairwise.venn(area1 = 303, area2 = 65, cross.area = 23, 
                   category = c("", ""),
                   col=c("#A5C68F","#BFE498"),
                   fill = c("#A5C68F","#BFE498"),
                   alpha = 0.4,
                   lwd =5)
dev.off()

VennDiag <- euler(c("A" = 280, "B" = 42, "A&B" = 23))

png("Venn_HY5_PHY.png", 1000, 1000, pointsize=100)
plot(VennDiag, fontsize=1000, cex=1000, alpha=0.9,
     fill=c("#A5C68F","#BFE498"),
     edges = list(col = c("black","black", "black"), lex = 2), 
     quantities = list(cex = 4),
     labels = c("", "", ""))
dev.off()


sum(orthogroups_mutants$Athy5_Atpif4, na.rm = TRUE)
sum(orthogroups_mutants$Mphy5_Mppif, na.rm = TRUE)
sum(orthogroups_mutants$hy5_pif, na.rm = TRUE)

png("Venn_HY5_PIF.png", 1000, 1000, pointsize=50)
draw.pairwise.venn(area1 = 79, area2 = 24, cross.area = 6, 
                   category = c("", ""),
                   col=c("#A5C68F","#BFE498"),
                   fill = c("#A5C68F","#BFE498"),
                   alpha = 0.4,
                   lwd =5)
dev.off()

VennDiag <- euler(c("A" = 73, "B" = 18, "A&B" = 6))

png("Venn_HY5_PIF.png", 1000, 1000, pointsize=100)
plot(VennDiag, fontsize=1000, cex=1000, alpha=0.9,
     fill=c("#A5C68F","#BFE498"),
     edges = list(col = c("black","black", "black"), lex = 2), 
     quantities = list(cex = 4),
     labels = c("", "", ""))
dev.off()

sum(orthogroups_mutants$Atphyb_Atpif4, na.rm = TRUE)
sum(orthogroups_mutants$Mpphy_Mppif, na.rm = TRUE)
sum(orthogroups_mutants$phy_pif, na.rm = TRUE)

png("Venn_PHY_PIF.png", 1000, 1000, pointsize=50)
draw.pairwise.venn(area1 = 61, area2 = 259, cross.area = 25, 
                   category = c("", ""),
                   col=c("#A5C68F","#BFE498"),
                   fill = c("#A5C68F","#BFE498"),
                   alpha = 0.4,
                   lwd =5)
dev.off()

VennDiag <- euler(c("A" = 36, "B" = 234, "A&B" = 25))

png("Venn_PHY_PIF.png", 1000, 1000, pointsize=100)
plot(VennDiag, fontsize=1000, cex=1000, alpha=0.9,
     fill=c("#A5C68F","#BFE498"),
     edges = list(col = c("black","black", "black"), lex = 2), 
     quantities = list(cex = 4),
     labels = c("", "", ""))
dev.off()


sum(orthogroups_mutants$Athy5_Atphyb_Atpif4, na.rm = TRUE)
sum(orthogroups_mutants$Mphy5_Mpphy_Mppif, na.rm = TRUE)
sum(orthogroups_mutants$hy5_phy_pif, na.rm = TRUE)

png("Venn_HY5_PHY_PIF.png", 1000, 1000, pointsize=50)
draw.pairwise.venn(area1 = 167, area2 = 132, cross.area = 27, 
                   category = c("", ""),
                   col=c("#A5C68F","#BFE498"),
                   fill = c("#A5C68F","#BFE498"),
                   alpha = 0.4,
                   lwd =5)
dev.off()

VennDiag <- euler(c("A" = 140, "B" = 105, "A&B" = 27))

png("Venn_HY5_PHY_PIF.png", 1000, 1000, pointsize=100)
plot(VennDiag, fontsize=1000, cex=1000, alpha=0.9,
     fill=c("#A5C68F","#BFE498"),
     edges = list(col = c("black","black", "black"), lex = 2), 
     quantities = list(cex = 4),
     labels = c("", "", ""))
dev.off()