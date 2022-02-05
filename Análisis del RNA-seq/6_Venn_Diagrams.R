# Author: Fernando Rodríguez Marín
# Contact: fernando.rodriguez.marin8@gmail.com
# Project: master's thesis

library(VennDiagram)

VennDiag <- euler(c("A" = 1405, "B" = 1807, "A&B" = 1380))

png("Venn_lighttemp_At.png", 1000, 1000, pointsize=100)
plot(VennDiag, fontsize=1000, cex=1000, alpha=0.7,
     fill=c("#7BB8DD","#E33228"),
     edges = list(col = c("black","black", "black"), lex = 2), 
     quantities = list(cex = 4),
     labels = c("", "", ""))
dev.off()

VennDiag <- euler(c("A" = 1089, "B" = 748, "A&B" = 1217))

png("Venn_lighttemp_Mp.png", 1000, 1000, pointsize=100)
plot(VennDiag, fontsize=1000, cex=1000, alpha=0.7,
     fill=c("#7BB8DD","#E33228"),
     edges = list(col = c("black","black", "black"), lex = 2), 
     quantities = list(cex = 4),
     labels = c("", "", ""))
dev.off()


VennDiag <- euler(c("A" = 205, "B" = 1489, "A&B" = 1264))

png("Venn_lighttemp_Me.png", 1000, 1000, pointsize=100)
plot(VennDiag, fontsize=1000, cex=1000, alpha=0.7,
     fill=c("#7BB8DD","#E33228"),
     edges = list(col = c("black","black", "black"), lex = 2), 
     quantities = list(cex = 4),
     labels = c("", "", ""))
dev.off()

VennDiag <- euler(c("A" = 1370, "B" = 324, "A&B" = 1049))

png("Venn_lighttemp_Ot.png", 1000, 1000, pointsize=100)
plot(VennDiag, fontsize=1000, cex=1000, alpha=0.7,
     fill=c("#7BB8DD","#E33228"),
     edges = list(col = c("black","black", "black"), lex = 2), 
     quantities = list(cex = 4),
     labels = c("", "", ""))
dev.off()
