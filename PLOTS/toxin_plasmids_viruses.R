###################
###Toxins, args###
##################

setwd("~/Desktop/R_icevirome/v3")

library(readr)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(svglite)

#Toxin

#import data
toxin_plasmid <- read_tsv("./PathoFact_plasmids_p/PathoFact_toxins_unique_counts.tsv")
toxin_virus <- read_tsv("./PathoFact_viruses/PathoFact_toxins_unique_counts.tsv")

#how many
nrow(toxin_plasmid)
#[1] 60
nrow(toxin_virus)
#[1] 48

#how many with > 1
nrow(toxin_plasmid %>% filter(Count > 1))
#[1] 34
nrow(toxin_virus %>% filter(Count > 1))
#[1] 28

#how many with > 2
nrow(toxin_plasmid %>% filter(Count > 2))
#[1] 22
nrow(toxin_virus %>% filter(Count > 2))
#[1] 22

#how many with >= 5
nrow(toxin_plasmid %>% filter(Count >= 5))
#[1] 16
nrow(toxin_virus %>% filter(Count >= 5))
#[1] 9

#Decided to go for 5 as threshold
toxin_plasmid5 <- data.frame(toxin_plasmid %>% filter(Count >= 5))
toxin_virus5 <- data.frame(toxin_virus %>% filter(Count >= 5))
toxin_plasmid5
toxin_virus5

#barplot
svglite("barplot_toxin.svg", width = 4, height = 8, fix_text_size=FALSE)

par(mfrow = c(2,1))   # 2 plots, 1 column

toxin_plasmid5 <- toxin_plasmid5[order(toxin_plasmid5$Count, decreasing = TRUE), ]
barplot(
  toxin_plasmid5$Count,
  names.arg = toxin_plasmid5$Description,
  las = 2,
  cex.names = 0.7,
  col = "black",
  ylab = "Count",
  main = "Plasmid toxins"
)

toxin_virus5 <- toxin_virus5[order(toxin_virus5$Count, decreasing = TRUE), ]
barplot(
  toxin_virus5$Count,
  names.arg = toxin_virus5$Description,
  las = 2,
  cex.names = 0.7,
  col = "black",
  ylab = "Count",
  main = "Virus toxins"
)

dev.off()

#AMR category

#import data
amr_cat_plasmid <- read_tsv("./PathoFact_plasmids_p/args_AMR_category_counts.tsv")[-1,]
amr_cat_virus <- read_tsv("./PathoFact_viruses/args_AMR_category_counts.tsv")[-1,]

#how many
nrow(amr_cat_plasmid)
#[1] 7
nrow(amr_cat_virus)
#[1] 14

#how many with > 1
nrow(amr_cat_plasmid %>% filter(Count > 1))
#[1] 5
nrow(amr_cat_virus %>% filter(Count > 1))
#[1] 8

#how many with > 2
nrow(amr_cat_plasmid %>% filter(Count > 2))
#[1] 3
nrow(amr_cat_virus %>% filter(Count > 2))
#[1] 7

#how many with >= 5
nrow(amr_cat_plasmid %>% filter(Count >= 5))
#[1] 0
nrow(amr_cat_virus %>% filter(Count >= 5))
#[1] 5

#Decided to go for 1 as threshold
amr_cat_plasmid1 <- data.frame(amr_cat_plasmid %>% filter(Count > 1))
amr_cat_virus1 <- data.frame(amr_cat_virus %>% filter(Count > 1))
amr_cat_plasmid1
amr_cat_virus1

#barplot
svglite("barplot_amr_cat.svg", width = 4, height = 8, fix_text_size=FALSE)

par(mfrow = c(2,1))   # 2 plots, 1 column

amr_cat_plasmid1 <- amr_cat_plasmid1[order(amr_cat_plasmid1$Count, decreasing = TRUE), ]
barplot(
  amr_cat_plasmid1$Count,
  names.arg = amr_cat_plasmid1$AMR_category,
  las = 2,
  cex.names = 0.7,
  col = "black",
  ylab = "Count",
  main = "Plasmid ARGs"
)

amr_cat_virus1 <- amr_cat_virus1[order(amr_cat_virus1$Count, decreasing = TRUE), ]
barplot(
  amr_cat_virus1$Count,
  names.arg = amr_cat_virus1$AMR_category,
  las = 2,
  cex.names = 0.7,
  col = "black",
  ylab = "Count",
  main = "Virus ARGs"
)
dev.off()

#AMR resistance mechanism

#import data
amr_res_plasmid <- read_tsv("./PathoFact_plasmids_p/args_Resistance_mechanism_counts.tsv")[-1,]
amr_res_virus <- read_tsv("./PathoFact_viruses/args_Resistance_mechanism_counts.tsv")[-1,]

#how many
nrow(amr_res_plasmid)
#[1] 5
nrow(amr_res_virus)
#[1] 4

#how many with > 1
nrow(amr_res_plasmid %>% filter(Count > 1))
#[1] 4
nrow(amr_res_virus %>% filter(Count > 1))
#[1] 2

#how many with > 2
nrow(amr_res_plasmid %>% filter(Count > 2))
#[1] 2
nrow(amr_res_virus %>% filter(Count > 2))
#[1] 4

#how many with >= 5
nrow(amr_res_plasmid %>% filter(Count >= 5))
#[1] 1
nrow(amr_res_virus %>% filter(Count >= 5))
#[1] 3

#Decided to go for 1 as threshold
amr_res_plasmid1 <- data.frame(amr_res_plasmid %>% filter(Count > 1))
amr_res_virus1 <- data.frame(amr_res_virus %>% filter(Count > 1))
amr_res_plasmid1
amr_res_virus1

#barplot
svglite("barplot_amr_res.svg", width = 4, height = 8, fix_text_size=FALSE)

par(mfrow = c(2,1))   # 2 plots, 1 column

amr_res_plasmid1 <- amr_res_plasmid1[order(amr_res_plasmid1$Count, decreasing = TRUE), ]
barplot(
  amr_res_plasmid1$Count,
  names.arg = amr_res_plasmid1$Resistance_mechanism,
  las = 2,
  cex.names = 0.7,
  col = "black",
  ylab = "Count",
  main = "Plasmid resistence mechanisms"
)

amr_res_virus1 <- amr_res_virus1[order(amr_res_virus1$Count, decreasing = TRUE), ]
barplot(
  amr_res_virus1$Count,
  names.arg = amr_res_virus1$Resistance_mechanism,
  las = 2,
  cex.names = 0.7,
  col = "black",
  ylab = "Count",
  main = "Virus resistence mechanism"
)

dev.off()


