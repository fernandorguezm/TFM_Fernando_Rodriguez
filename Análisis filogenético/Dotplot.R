# Author: Fernando Rodríguez Marín
# Contact: fernando.rodriguez.marin8@gmail.com
# Project: master's thesis

tabla <- read.csv(file = "Tabla_dotplot.csv", sep = ";", header = T)

library("reshape2")
library(ggplot2)
#convert data frame from a "wide" format to a "long" format
pcm = melt(tabla, id = c("ï.."))
pcm$ï.. <- factor(pcm$ï..,levels=unique(pcm$ï..))


colours = c( "#7F7F7F",  "#A20000", "#FF0000", "#FF3300","#FF6600",
             "#F4B183", "#FFC000", "#92D050", "#00B050",
             "#9DC3E6", "#0070C0", "#002060", "#7030A0")

tiff("Plot2.tif", res = 50)

plot <- ggplot(pcm, aes(x = ï.., y = variable)) +
  geom_point(aes(size = value, fill = ï..), alpha = 0.75, shape = 21) +
  scale_size_continuous(limits = c(1, 50), range = c(1,18), breaks = c(1,4,7,10,13,16,19)) +
  labs( x= "", y = "", size =  "Abundancia", fill = "")  +
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, angle = 45, vjust = 1, hjust = 1), 
        axis.text.y = element_text(colour = "black", size = 7, face="italic"), 
        legend.text = element_text(size = 7, face ="bold", colour ="black" ), 
        legend.title = element_text(size = 10, face = "bold"), 
        panel.background = element_blank(), 
        legend.position = "right",
        rect = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = colours)  +
  scale_y_discrete(limits = rev(levels(pcm$variable))) 

ggsave(
  "ggtest_es.png",
  plot,
  width = 7,
  height = 7,
  dpi = 1000
)
