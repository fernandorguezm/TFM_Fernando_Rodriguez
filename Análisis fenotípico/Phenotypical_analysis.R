# Author: Fernando Rodríguez Marín
# Contact: fernando.rodriguez.marin8@gmail.com
# Project: master's thesis

library(ggplot2)
library(RColorBrewer)

#####################################################################
#####################################################################
################ Arabidopsis thaliana mutants #######################
#####################################################################
#####################################################################

display.brewer.pal(n = 11, name = 'RdYlBu')
brewer.pal(n = 11, name = 'RdYlBu')

par(mfrow=c(4,1))

# Col-0
data <- read.csv(file = "Archivos/At_Hypocotyl_length.csv", sep = ";")
data <- data[1:358,]
tail(data)
data$Light.intensity <- factor(data$Light.intensity, levels=c("1", "10", "100"))

# Option 3

col0 <- ggplot(data, aes(x=Light.intensity, y=Hypocotyl.length..cm., group=interaction(Light.intensity, Temperature))) +
  stat_boxplot(geom="errorbar", color = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_boxplot(alpha=1, fill=c("#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D"), 
               color=c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(width=0.75),
               binwidth = 0.015, alpha = 1, 
               aes(fill=Temperature)) +
  
  ggtitle("Col-0") +
  scale_fill_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  scale_color_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  xlab("Light intensity (uE)")+ 
  ylab("Hypocotyl length (cm)")+
  theme_classic()  +
  ylim(c(0,1.5))

col0 <- col0 + theme(legend.position="none",
               text = element_text(size=30)
               )

ggsave(
  "Col0.png",
  col0,
  width = 7,
  height = 7,
  dpi = 500
)

# Atpif4

data <- read.csv(file = "Archivos/At_Hypocotyl_length.csv", sep = ";")
data <- data[359:715,]
tail(data)
data$Light.intensity <- factor(data$Light.intensity, levels=c("1", "10", "100"))

# Option 3

atpif4 <- ggplot(data, aes(x=Light.intensity, y=Hypocotyl.length..cm., group=interaction(Light.intensity, Temperature))) +
  stat_boxplot(geom="errorbar", color = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_boxplot(alpha=1, fill=c("#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D"), 
               color=c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(width=0.75),
               binwidth = 0.015, alpha = 1, 
               aes(fill=Temperature))+
  ggtitle("Atpif4-2") +
  scale_fill_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  scale_color_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  xlab("Light intensity (uE)")+ 
  ylab("Hypocotyl length (cm)")+
  theme_classic() +
  ylim(c(0,1.5)) 

atpif4 <- atpif4 + theme(legend.position="none",
                     text = element_text(size=30)
)

ggsave(
  "AtPIF4.png",
  atpif4,
  width = 7,
  height = 7,
  dpi = 500
)

# Atphyb

data <- read.csv(file = "Archivos/At_Hypocotyl_length.csv", sep = ";")
data <- data[716:1066,]
tail(data)
data$Light.intensity <- factor(data$Light.intensity, levels=c("1", "10", "100"))

# Option 3

atphyb <- ggplot(data, aes(x=Light.intensity, y=Hypocotyl.length..cm., group=interaction(Light.intensity, Temperature))) +
  stat_boxplot(geom="errorbar", color = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_boxplot(alpha=1, fill=c("#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D"), 
               color=c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(width=0.75),
               binwidth = 0.015, alpha = 1, 
               aes(fill=Temperature))+
  ggtitle("Atphyb-9") +
  scale_fill_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  scale_color_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  xlab("Light intensity (uE)")+ 
  ylab("Hypocotyl length (cm)")+
  theme_classic() +
  ylim(c(0,1.5)) 

atphyb <- atphyb + theme(legend.position="none",
                         text = element_text(size=30)
)
ggsave(
  "AtPHYB.png",
  atphyb,
  width = 7,
  height = 7,
  dpi = 500
)

# Athy5

data <- read.csv(file = "Archivos/At_Hypocotyl_length.csv", sep = ";")
data <- data[1037:1396,]
tail(data)
data$Light.intensity <- factor(data$Light.intensity, levels=c("1", "10", "100"))

# Option 3

athy5 <- ggplot(data, aes(x=Light.intensity, y=Hypocotyl.length..cm., group=interaction(Light.intensity, Temperature))) +
  stat_boxplot(geom="errorbar", color = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_boxplot(alpha=1, fill=c("#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D"), 
               color=c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(width=0.75),
               binwidth = 0.015, alpha = 1, 
               aes(fill=Temperature))+
  ggtitle("Athy5-2") +
  scale_fill_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  scale_color_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  xlab("Light intensity (uE)")+ 
  ylab("Hypocotyl length (cm)")+
  theme_classic() +
  ylim(c(0,1.5)) 

athy5 <- athy5 + theme(legend.position="none",
                         text = element_text(size=30)
)
ggsave(
  "AtHY5.png",
  athy5,
  width = 7,
  height = 7,
  dpi = 500
)

#####################################################################
#####################################################################
################### Marchantia polymorpha  mutants ##################
#####################################################################
#####################################################################

display.brewer.pal(n = 11, name = 'RdYlBu')
brewer.pal(n = 11, name = 'RdYlBu')

# Tak-1
data <- read.csv(file = "Archivos/Marchantia_thallus_surface_area.csv", sep = ";")
data <- data[1:167,]
tail(data)
data$Light.intensity <- factor(data$Light.intensity, levels=c("1", "10", "100"))

tak1 <- ggplot(data, aes(x=Light.intensity, y=Thallus.surface.area..mm2., group=interaction(Light.intensity, Temperature))) +
  stat_boxplot(geom="errorbar", color = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_boxplot(alpha=1, fill=c("#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D"), 
               color=c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(width=0.75),
               binwidth = 0.8, alpha = 1, 
               aes(fill=Temperature))+
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("TAK-1") +
  scale_fill_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  scale_color_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  xlab("Light intensity (uE)")+ 
  ylab("Thallus surface area (mm2)")+
  theme_classic()+
  ylim(c(0,75))

tak1 <- tak1 + theme(legend.position="none",
                       text = element_text(size=30)
)


ggsave(
  "tak1.png",
  tak1,
  width = 7,
  height = 7,
  dpi = 500
)

# Mppif-1
data <- read.csv(file = "Archivos/Marchantia_thallus_surface_area.csv", sep = ";")
data <- data[168:332,]
tail(data)
data$Light.intensity <- factor(data$Light.intensity, levels=c("1", "10", "100"))

mppif <- ggplot(data, aes(x=Light.intensity, y=Thallus.surface.area..mm2., group=interaction(Light.intensity, Temperature))) +
  stat_boxplot(geom="errorbar", color = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_boxplot(alpha=1, fill=c("#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D"), 
               color=c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(width=0.75),
               binwidth = 0.8, alpha = 1, 
               aes(fill=Temperature))+
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Mppif-1") +
  scale_fill_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  scale_color_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  xlab("Light intensity (uE)")+ 
  ylab("Thallus surface area (mm2)")+
  theme_classic()+
  ylim(c(0,75))

mppif <- mppif + theme(legend.position="none",
                     text = element_text(size=30)
)
ggsave(
  "mppif.png",
  mppif,
  width = 7,
  height = 7,
  dpi = 500
)
# Mpphy-1
data <- read.csv(file = "Archivos/Marchantia_thallus_surface_area.csv", sep = ";")
data <- data[333:498,]
tail(data)
data$Light.intensity <- factor(data$Light.intensity, levels=c("1", "10", "100"))

mpphy <- ggplot(data, aes(x=Light.intensity, y=Thallus.surface.area..mm2., group=interaction(Light.intensity, Temperature))) +
  stat_boxplot(geom="errorbar", color = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_boxplot(alpha=1, fill=c("#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D"), 
               color=c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(width=0.75),
               binwidth = 0.8, alpha = 1, 
               aes(fill=Temperature))+
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Mpphy-1") +
  scale_fill_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  scale_color_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  xlab("Light intensity (uE)")+ 
  ylab("Thallus surface area (mm2)")+
  theme_classic() +
  ylim(c(0,75))

mpphy <- mpphy + theme(legend.position="none",
                       text = element_text(size=30)
)
ggsave(
  "mpphy.png",
  mpphy,
  width = 7,
  height = 7,
  dpi = 500
)

# Mphy5-1
data <- read.csv(file = "Archivos/Marchantia_thallus_surface_area.csv", sep = ";")
data <- data[499:653,]
tail(data)
data$Light.intensity <- factor(data$Light.intensity, levels=c("1", "10", "100"))

mphy5 <- ggplot(data, aes(x=Light.intensity, y=Thallus.surface.area..mm2., group=interaction(Light.intensity, Temperature))) +
  stat_boxplot(geom="errorbar", color = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_boxplot(alpha=1, fill=c("#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D"), 
               color=c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(width=0.75),
               binwidth = 0.8, alpha = 1, 
               aes(fill=Temperature))+
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Mphy5-1") +
  scale_fill_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  scale_color_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  xlab("Light intensity (uE)")+ 
  ylab("Thallus surface area (mm2)")+
  theme_classic() +
  ylim(c(0,75))

mphy5 <- mphy5 + theme(legend.position="none",
                       text = element_text(size=30)
)
ggsave(
  "mphy5.png",
  mphy5,
  width = 7,
  height = 7,
  dpi = 500
)

#####################################################################
#####################################################################
############################ Wild-types #############################
#####################################################################
#####################################################################

# Arabidopsis thaliana

data <- read.csv(file = "Archivos/At_Hypocotyl_length.csv", sep = ";")
data <- data[1:358,]
tail(data)
data$Light.intensity <- factor(data$Light.intensity, levels=c("1", "10", "100"))

# Option 3

col0_without_ylim <- ggplot(data, aes(x=Light.intensity, y=Hypocotyl.length..cm., group=interaction(Light.intensity, Temperature))) +
  stat_boxplot(geom="errorbar", color = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_boxplot(alpha=1, fill=c("#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D"), 
               color=c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(width=0.75),
               binwidth = 0.008, alpha = 1, 
               aes(fill=Temperature)) +
  
  ggtitle("Col-0") +
  scale_fill_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  scale_color_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  xlab("Light intensity (uE)")+ 
  ylab("Hypocotyl length (cm)")+
  theme_classic() 

col0_without_ylim <- col0_without_ylim + theme(legend.position="none",
                     text = element_text(size=30)
)

ggsave(
  "Col0_without_ylim.png",
  col0_without_ylim,
  width = 6,
  height = 7,
  dpi = 500
)

# Marchantia polymorpha
data <- read.csv(file = "Archivos/Marchantia_thallus_surface_area.csv", sep = ";")
data <- data[1:167,]
tail(data)
data$Light.intensity <- factor(data$Light.intensity, levels=c("1", "10", "100"))

tak1_whitout_ylim <- ggplot(data, aes(x=Light.intensity, y=Thallus.surface.area..mm2., group=interaction(Light.intensity, Temperature))) +
  stat_boxplot(geom="errorbar", color = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_boxplot(alpha=1, fill=c("#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D"), 
               color=c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(width=0.75),
               binwidth = 0.9, alpha = 1, 
               aes(fill=Temperature))+
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("TAK-1") +
  scale_fill_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  scale_color_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  xlab("Light intensity (uE)")+ 
  ylab("Thallus surface area (mm2)")+
  theme_classic()+
  ylim(c(0,75))

tak1_whitout_ylim <- tak1_whitout_ylim + theme(legend.position="none",
                     text = element_text(size=30)
)


ggsave(
  "tak1_whitout_ylim.png",
  tak1_whitout_ylim,
  width = 6,
  height = 7,
  dpi = 500
)

# Mesotaenium endlicherianum

data <- read.csv(file = "Archivos/Mesotaenium_cell_number.csv", sep = ";")
data$Light.intensity <- factor(data$Light.intensity, levels=c("1", "10", "100"))

me <- ggplot(data, aes(x=Light.intensity, y=Cell.number, group=interaction(Light.intensity, Temperature))) +
  stat_boxplot(geom="errorbar", color = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_boxplot(alpha=1, fill=c("#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D"), 
               color=c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(width=0.75),
               binwidth = 3.5, alpha = 1, 
               aes(fill=Temperature))+
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Mesotaenium endlicherianum") +
  scale_fill_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  scale_color_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  xlab("Light intensity (uE)")+ 
  ylab("Cell number")+
  theme_classic() 

me <- me + theme(legend.position="none",
                 text = element_text(size=30)
)
ggsave(
  "Mesotaenium.png",
  me,
  width = 6,
  height = 7,
  dpi = 500
)

# Ostreococcus tauri

data <- read.csv(file = "Archivos/Ot_Number_Cells.csv", sep = ";")
data$Light <- factor(data$Light, levels=c("10", "100"))

ot <- ggplot(data, aes(x=Light, y=Cells, group=interaction(Light, Temperature))) +
  stat_boxplot(geom="errorbar", color = c("#4575B4", "#74ADD1", "#FDAE61","#4575B4", "#74ADD1", "#FDAE61")) +
  geom_boxplot(alpha=1, fill=c("#AEC3E0", "#D4ECF4", "#FED5AC","#AEC3E0", "#D4ECF4", "#FED5AC"), 
               color=c("#4575B4", "#74ADD1", "#FDAE61","#4575B4", "#74ADD1", "#FDAE61")) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(width=0.75),
               binwidth = 75, alpha = 1, 
               aes(fill=Temperature))+
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Ostreococcus tauri") +
  scale_fill_manual(values = c("#4575B4", "#74ADD1", "#FDAE61"))+
  scale_color_manual(values = c("#4575B4", "#74ADD1", "#FDAE61"))+
  xlab("Light intensity (uE)")+ 
  ylab("Cell number")+
  theme_classic() 

ot <- ot + theme(legend.position="none",
                 text = element_text(size=30)
)
ggsave(
  "Ostreococcus.png",
  ot,
  width = 6,
  height = 7,
  dpi = 500
)

#####################################################################
#####################################################################
####################### Statistical Analysis ########################
#####################################################################
#####################################################################

library(multcompView)

# Arabidopsis thaliana 

data <- read.csv(file = "Archivos/Arabidopsis_statistical.csv", sep = ";")
data <- data[1:1396,]
tail(data)

model <- aov(data = data, Hypocotyl.length..cm.~Sample)
hsd <- TukeyHSD(model)

letters_at <- multcompLetters(extract_p(hsd$Sample))
letters_at

# Arabidopsis thaliana 

data <- read.csv(file = "Archivos/Arabidopsis_statistical.csv", sep = ";")
data <- data[1:358,]
tail(data)

model <- aov(data = data, Hypocotyl.length..cm.~Sample)
hsd <- TukeyHSD(model)

letters_at <- multcompLetters(extract_p(hsd$Sample))
letters_at

# Marchantia polymorpha

data <- read.csv(file = "Archivos/Marchantia_statistical.csv", sep = ";")
data <- data[1:653,]
tail(data)

model <- aov(data = data, Thallus.surface.area..mm2.~Sample)
hsd <- TukeyHSD(model)

letters_mp <- multcompLetters(extract_p(hsd$Sample))

# Marchantia polymorpha WT

data <- read.csv(file = "Archivos/Marchantia_statistical.csv", sep = ";")
data <- data[1:167,]
tail(data)

model <- aov(data = data, Thallus.surface.area..mm2.~Sample)
hsd <- TukeyHSD(model)

letters_mp <- multcompLetters(extract_p(hsd$Sample))

# Mesotaenium endlicherianum

data <- read.csv(file = "Archivos/Mesotaenium_statistical.csv", sep = ";")

model <- aov(data = data,  Cell.number~Sample)
hsd <- TukeyHSD(model)

letters_me <- multcompLetters(extract_p(hsd$Sample))

# Ostreococcus tauri

data <- read.csv(file = "Archivos/Ot_statistical.csv", sep = ";")

model <- aov(data = data, Cells~Sample)
hsd <- TukeyHSD(model)

letters_ot <- multcompLetters(extract_p(hsd$Sample))


#####################################################################
#####################################################################
########################### Complementary ###########################
#####################################################################
#####################################################################

# Mesotaenium endlicherianum

data <- read.csv(file = "Archivos/Mesotaenium_cell_length.csv", sep = ";")
data$Light.intensity <- factor(data$Light.intensity, levels=c("1", "10", "100"))

me <- ggplot(data, aes(x=Light.intensity, y=Cell.length, group=interaction(Light.intensity, Temperature))) +
  stat_boxplot(geom="errorbar", color = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_boxplot(alpha=1, fill=c("#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D","#AEC3E0", "#D4ECF4", "#FED5AC", "#EB918D"), 
               color=c("#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027","#4575B4", "#74ADD1", "#FDAE61", "#D73027")) +
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge(width=0.75),
               binwidth = 0.5, alpha = 1, 
               aes(fill=Temperature, color=Temperature))+
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Mesotaenium endlicherianum") +
  scale_fill_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  scale_color_manual(values = c("#4575B4", "#74ADD1", "#FDAE61", "#D73027"))+
  xlab("Light intensity (uE)")+ 
  ylab("Cell length (um)")+
  theme_classic() 

me <- me + theme(legend.position="none",
                 text = element_text(size=30)
)
ggsave(
  "Mesotaenium_cell_length.png",
  me,
  width = 7,
  height = 7,
  dpi = 500
)

# Tests

data <- read.csv(file = "Archivos/Mesotaenium_cell_length_dens.csv", sep = ";")
tail(data)
data$Light.intensity <- factor(data$Light.intensity, levels=c("1", "10", "100"))


ggplot(data, aes(x=Light.intensity, y=Cell.length, group=interaction(Light.intensity, Temperature))) + 
  geom_violin()

data <- data[3565:6353,]

ggplot(data, aes(x=Cell.length, group=interaction(Light.intensity, Temperature))) + 
  geom_density(alpha=0.6)
