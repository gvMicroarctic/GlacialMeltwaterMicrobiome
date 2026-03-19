##########
###rDNA###
##########

library(ggplot2)
library(svglite)

setwd("~/Desktop/R_icevirome/v3")

rdna <- read.table("./rDNA.txt", header = T, check.names = FALSE)
rdna

p_rdna <- ggplot(rdna, aes(x = sample, y = perc_rDNA, fill = vb)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        plot.title = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black")) +
  labs(x = "Sample", y = "Reads mapped to rRNA genes (%)") +
  scale_fill_manual(values = c("b" = "darkgreen", "v" = "darkorange")) + 
  geom_hline(yintercept = 0.483, color = "black", linetype = "dashed", linewidth = 0.7) +
  geom_hline(yintercept = 0.271, color = "black", linetype = "dotted", linewidth = 0.7)
p_rdna

#save plot
svglite(filename="rDNA_barplot.svg", width=5.5, height=5.5, fix_text_size=FALSE)
p_rdna
dev.off()
