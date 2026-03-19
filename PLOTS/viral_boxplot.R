###########################################################
####Analyses on entire dataset (only amplified samples)####
###########################################################

library(factoextra)
library(svglite)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(vegan)
library(gridExtra)
library(Maaslin2)
library(tibble)

setwd("~/Desktop/R_icevirome/v3")

#############
#count table#
#############

#get count information
votu_count <- read.table("vOTU_abundance_10.txt", sep = "\t", header = T, row.names = 1)

#get only columns with tpm counts
votu_tpm0 <- votu_count[, grepl("TPM", colnames(votu_count))]
#format sample names (columns)
colnames(votu_tpm0) <- gsub("^a\\.|_si\\.TPM$", "", colnames(votu_tpm0))
dim(votu_tpm0)
#[1] 7474     23
votu_tpm <- votu_tpm0[,-(ncol(votu_tpm0))]
dim(votu_tpm)
#[1] 7474     22
votu_tpm_name <- cbind(votu_tpm, contig = row.names(votu_tpm)) #I need it for below steps

#remove big data.frame
rm(votu_count)

#add contig information to count_entire_tpm
tpm_contig <- cbind(votu_tpm, contig = row.names(votu_tpm_name))
head(tpm_contig)

############
##Taxonomy##
############

#get data from file
taxonomy <- read.table("votu_taxonomy_clean.txt", row.names = 1)
colnames(taxonomy) <- c("domain", "phylum", "class", "order", "family")
dim(taxonomy)
#[1] 7474    5

#format taxonomic names
taxonomy$domain <-gsub("d__","",as.character(taxonomy$domain))
taxonomy$phylum <-gsub("p__","",as.character(taxonomy$phylum))
taxonomy$class <-gsub("c__","",as.character(taxonomy$class))
taxonomy$order <-gsub("o__","",as.character(taxonomy$order))
taxonomy$family <-gsub("f__","",as.character(taxonomy$family))

#get taxonomy correspondences to each family
path <- as.data.frame(taxonomy[!duplicated(taxonomy[,5]),])

#check how many vOTUs taxonomically assigned for each taxonomic level
nrow(taxonomy %>% filter(phylum != "Unclassified"))/nrow(taxonomy)*100
#[1] 35.34921
nrow(taxonomy %>% filter(class != "Unclassified"))/nrow(taxonomy)*100
#[1] 34.72036
nrow(taxonomy %>% filter(family != "Unclassified"))/nrow(taxonomy)*100
#[1] 13.27268

#add contig information to taxonomy
taxonomy_c <- cbind(taxonomy, contig = row.names(taxonomy))
taxonomy_contig <- merge(tpm_contig, taxonomy_c, by = "contig", all = FALSE)

#####################
##Taxonomy - phylum##
#####################

#N.B. violin plots do not work for this dataser because not enouth data to compute kernel density extimation; so boxplots

#sum counts
phylum_contig_sum <- taxonomy_contig %>%
  group_by(phylum) %>%
  dplyr::summarize(sum_BR1 = sum(as.numeric(BR1)),sum_BR2 = sum(as.numeric(BR2)),sum_BR3 = sum(as.numeric(BR3)),
                   sum_BR4 = sum(as.numeric(BR4)),sum_BR5 = sum(as.numeric(BR5)),
                   
                   sum_BM1a = sum(as.numeric(BM1a)),sum_BM1b = sum(as.numeric(BM1b)),sum_BM2 = sum(as.numeric(BM2)),
                   sum_BM3 = sum(as.numeric(BM3)),sum_BM4 = sum(as.numeric(BM4)),sum_BM5 = sum(as.numeric(BM5)),
                   
                   sum_VR1 = sum(as.numeric(VR1)),sum_VR2 = sum(as.numeric(VR2)),sum_VR3 = sum(as.numeric(VR3)),
                   sum_VR4 = sum(as.numeric(VR4)),sum_VR5 = sum(as.numeric(VR5)),
                   
                   sum_VM1a = sum(as.numeric(VM1a)),sum_VM1b = sum(as.numeric(VM1b)),sum_VM3a = sum(as.numeric(VM3a)),
                   sum_VM3b = sum(as.numeric(VM3b)),sum_VM4 = sum(as.numeric(VM4)),sum_VM5 = sum(as.numeric(VM5)))

phylum_contig_sum_b <- taxonomy_contig %>%
  group_by(phylum) %>%
  dplyr::summarize(sum_BR1 = sum(as.numeric(BR1)),sum_BR2 = sum(as.numeric(BR2)),sum_BR3 = sum(as.numeric(BR3)),
                   sum_BR4 = sum(as.numeric(BR4)),sum_BR5 = sum(as.numeric(BR5)),
                   
                   sum_BM1a = sum(as.numeric(BM1a)),sum_BM1b = sum(as.numeric(BM1b)),sum_BM2 = sum(as.numeric(BM2)),
                   sum_BM3 = sum(as.numeric(BM3)),sum_BM4 = sum(as.numeric(BM4)),sum_BM5 = sum(as.numeric(BM5)))

phylum_contig_sum_v <- taxonomy_contig %>%
  group_by(phylum) %>%
  dplyr::summarize(sum_VR1 = sum(as.numeric(VR1)),sum_VR2 = sum(as.numeric(VR2)),sum_VR3 = sum(as.numeric(VR3)),
                   sum_VR4 = sum(as.numeric(VR4)),sum_VR5 = sum(as.numeric(VR5)),
                   
                   sum_VM1a = sum(as.numeric(VM1a)),sum_VM1b = sum(as.numeric(VM1b)),sum_VM3a = sum(as.numeric(VM3a)),
                   sum_VM3b = sum(as.numeric(VM3b)),sum_VM4 = sum(as.numeric(VM4)),sum_VM5 = sum(as.numeric(VM5)))

#reshape the data into long format
phylum_contig_long <- phylum_contig_sum %>%
  pivot_longer(cols = starts_with("sum"), names_to = "sample", values_to = "value")

#add I vs II filters, and M vs R
phylum_contig_long_dt <- phylum_contig_long  %>% 
  mutate(dataset = case_when(grepl("BM", sample) ~ "B", grepl("BR", sample) ~ "B",
                             grepl("VM", sample) ~ "V", grepl("VR", sample) ~ "V", TRUE ~ NA_character_))

#order + get top 15 phyla (excluding Unclassified)
top <- phylum_contig_long_dt %>%
  filter(phylum != "Unclassified") %>%
  group_by(phylum) %>%
  summarise(total = sum(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 15) %>%
  pull(phylum)

#filter + factor order
phylum_contig_long_dt_top <- phylum_contig_long_dt %>%
  filter(phylum %in% top) %>%                         # remove non-top phyla
  mutate(
    phylum  = factor(phylum, levels = rev(top)),      # ordered taxa
    dataset = factor(dataset, levels = c("B", "V"))   # dataset order
  )

#plot
p_boxplot_phylum_dataset <- phylum_contig_long_dt_top %>%
  ggplot(aes(x = phylum, y = value, colour = dataset)) +
  geom_boxplot(outlier.shape = NA, width = 0.8) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 2,
    fill = "white",
    position = position_dodge(width = 0.8)
  ) +
  coord_flip() +
  labs(x = "Phylum", y = "TPM") +
  theme_bw() +
  theme(
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  scale_colour_manual(values = c("B" = "darkgreen", "V" = "darkorange"))

p_boxplot_phylum_dataset

#####################
##Taxonomy - class##
#####################

#N.B. violin plots do not work for this dataser because not enouth data to compute kernel density extimation; so boxplots

#sum counts
class_contig_sum <- taxonomy_contig %>%
  group_by(class) %>%
  dplyr::summarize(sum_BR1 = sum(as.numeric(BR1)),sum_BR2 = sum(as.numeric(BR2)),sum_BR3 = sum(as.numeric(BR3)),
                   sum_BR4 = sum(as.numeric(BR4)),sum_BR5 = sum(as.numeric(BR5)),
                   
                   sum_BM1a = sum(as.numeric(BM1a)),sum_BM1b = sum(as.numeric(BM1b)),sum_BM2 = sum(as.numeric(BM2)),
                   sum_BM3 = sum(as.numeric(BM3)),sum_BM4 = sum(as.numeric(BM4)),sum_BM5 = sum(as.numeric(BM5)),
                   
                   sum_VR1 = sum(as.numeric(VR1)),sum_VR2 = sum(as.numeric(VR2)),sum_VR3 = sum(as.numeric(VR3)),
                   sum_VR4 = sum(as.numeric(VR4)),sum_VR5 = sum(as.numeric(VR5)),
                   
                   sum_VM1a = sum(as.numeric(VM1a)),sum_VM1b = sum(as.numeric(VM1b)),sum_VM3a = sum(as.numeric(VM3a)),
                   sum_VM3b = sum(as.numeric(VM3b)),sum_VM4 = sum(as.numeric(VM4)),sum_VM5 = sum(as.numeric(VM5)))

class_contig_sum_b <- taxonomy_contig %>%
  group_by(class) %>%
  dplyr::summarize(sum_BR1 = sum(as.numeric(BR1)),sum_BR2 = sum(as.numeric(BR2)),sum_BR3 = sum(as.numeric(BR3)),
                   sum_BR4 = sum(as.numeric(BR4)),sum_BR5 = sum(as.numeric(BR5)),
                   
                   sum_BM1a = sum(as.numeric(BM1a)),sum_BM1b = sum(as.numeric(BM1b)),sum_BM2 = sum(as.numeric(BM2)),
                   sum_BM3 = sum(as.numeric(BM3)),sum_BM4 = sum(as.numeric(BM4)),sum_BM5 = sum(as.numeric(BM5)))

class_contig_sum_v <- taxonomy_contig %>%
  group_by(class) %>%
  dplyr::summarize(sum_VR1 = sum(as.numeric(VR1)),sum_VR2 = sum(as.numeric(VR2)),sum_VR3 = sum(as.numeric(VR3)),
                   sum_VR4 = sum(as.numeric(VR4)),sum_VR5 = sum(as.numeric(VR5)),
                   
                   sum_VM1a = sum(as.numeric(VM1a)),sum_VM1b = sum(as.numeric(VM1b)),sum_VM3a = sum(as.numeric(VM3a)),
                   sum_VM3b = sum(as.numeric(VM3b)),sum_VM4 = sum(as.numeric(VM4)),sum_VM5 = sum(as.numeric(VM5)))

#reshape the data into long format
class_contig_long <- class_contig_sum %>%
  pivot_longer(cols = starts_with("sum"), names_to = "sample", values_to = "value")

#add I vs II filters, and M vs R
class_contig_long_dt <- class_contig_long  %>% 
  mutate(dataset = case_when(grepl("BM", sample) ~ "B", grepl("BR", sample) ~ "B",
                             grepl("VM", sample) ~ "V", grepl("VR", sample) ~ "V", TRUE ~ NA_character_))

#order + get top 15 phyla (excluding Unclassified)
top <- class_contig_long_dt %>%
  filter(class != "Unclassified") %>%
  group_by(class) %>%
  summarise(total = sum(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 15) %>%
  pull(class)

#filter + factor order
class_contig_long_dt_top <- class_contig_long_dt %>%
  filter(class %in% top) %>%                         # remove non-top phyla
  mutate(
    class  = factor(class, levels = rev(top)),      # ordered taxa
    dataset = factor(dataset, levels = c("B", "V"))   # dataset order
  )

#plot
p_boxplot_class_dataset <- class_contig_long_dt_top %>%
  ggplot(aes(x = class, y = value, colour = dataset)) +
  geom_boxplot(outlier.shape = NA, width = 0.8) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 2,
    fill = "white",
    position = position_dodge(width = 0.8)
  ) +
  coord_flip() +
  labs(x = "Class", y = "TPM") +
  theme_bw() +
  theme(
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  scale_colour_manual(values = c("B" = "darkgreen", "V" = "darkorange"))

p_boxplot_class_dataset

#####################
##Taxonomy - family##
#####################

#N.B. violin plots do not work for this dataser because not enouth data to compute kernel density extimation; so boxplots

#sum counts
family_contig_sum <- taxonomy_contig %>%
  group_by(family) %>%
  dplyr::summarize(sum_BR1 = sum(as.numeric(BR1)),sum_BR2 = sum(as.numeric(BR2)),sum_BR3 = sum(as.numeric(BR3)),
                   sum_BR4 = sum(as.numeric(BR4)),sum_BR5 = sum(as.numeric(BR5)),
                   
                   sum_BM1a = sum(as.numeric(BM1a)),sum_BM1b = sum(as.numeric(BM1b)),sum_BM2 = sum(as.numeric(BM2)),
                   sum_BM3 = sum(as.numeric(BM3)),sum_BM4 = sum(as.numeric(BM4)),sum_BM5 = sum(as.numeric(BM5)),
                   
                   sum_VR1 = sum(as.numeric(VR1)),sum_VR2 = sum(as.numeric(VR2)),sum_VR3 = sum(as.numeric(VR3)),
                   sum_VR4 = sum(as.numeric(VR4)),sum_VR5 = sum(as.numeric(VR5)),
                   
                   sum_VM1a = sum(as.numeric(VM1a)),sum_VM1b = sum(as.numeric(VM1b)),sum_VM3a = sum(as.numeric(VM3a)),
                   sum_VM3b = sum(as.numeric(VM3b)),sum_VM4 = sum(as.numeric(VM4)),sum_VM5 = sum(as.numeric(VM5)))

family_contig_sum_b <- taxonomy_contig %>%
  group_by(family) %>%
  dplyr::summarize(sum_BR1 = sum(as.numeric(BR1)),sum_BR2 = sum(as.numeric(BR2)),sum_BR3 = sum(as.numeric(BR3)),
                   sum_BR4 = sum(as.numeric(BR4)),sum_BR5 = sum(as.numeric(BR5)),
                   
                   sum_BM1a = sum(as.numeric(BM1a)),sum_BM1b = sum(as.numeric(BM1b)),sum_BM2 = sum(as.numeric(BM2)),
                   sum_BM3 = sum(as.numeric(BM3)),sum_BM4 = sum(as.numeric(BM4)),sum_BM5 = sum(as.numeric(BM5)))

family_contig_sum_v <- taxonomy_contig %>%
  group_by(family) %>%
  dplyr::summarize(sum_VR1 = sum(as.numeric(VR1)),sum_VR2 = sum(as.numeric(VR2)),sum_VR3 = sum(as.numeric(VR3)),
                   sum_VR4 = sum(as.numeric(VR4)),sum_VR5 = sum(as.numeric(VR5)),
                   
                   sum_VM1a = sum(as.numeric(VM1a)),sum_VM1b = sum(as.numeric(VM1b)),sum_VM3a = sum(as.numeric(VM3a)),
                   sum_VM3b = sum(as.numeric(VM3b)),sum_VM4 = sum(as.numeric(VM4)),sum_VM5 = sum(as.numeric(VM5)))

#reshape the data into long format
family_contig_long <- family_contig_sum %>%
  pivot_longer(cols = starts_with("sum"), names_to = "sample", values_to = "value")

#add I vs II filters, and M vs R
family_contig_long_dt <- family_contig_long  %>% 
  mutate(dataset = case_when(grepl("BM", sample) ~ "B", grepl("BR", sample) ~ "B",
                             grepl("VM", sample) ~ "V", grepl("VR", sample) ~ "V", TRUE ~ NA_character_))

#order + get top 15 phyla (excluding Unclassified)
top <- family_contig_long_dt %>%
  filter(family != "Unclassified") %>%
  group_by(family) %>%
  summarise(total = sum(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 15) %>%
  pull(family)

#filter + factor order
family_contig_long_dt_top <- family_contig_long_dt %>%
  filter(family %in% top) %>%                         # remove non-top phyla
  mutate(
    family  = factor(family, levels = rev(top)),      # ordered taxa
    dataset = factor(dataset, levels = c("B", "V"))   # dataset order
  )

#plot
p_boxplot_family_dataset <- family_contig_long_dt_top %>%
  ggplot(aes(x = family, y = value, colour = dataset)) +
  geom_boxplot(outlier.shape = NA, width = 0.8) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 2,
    fill = "white",
    position = position_dodge(width = 0.8)
  ) +
  coord_flip() +
  labs(x = "Family", y = "TPM") +
  theme_bw() +
  theme(
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  scale_colour_manual(values = c("B" = "darkgreen", "V" = "darkorange"))

p_boxplot_family_dataset

#plot
svglite(filename="viral_barplots_SI.svg", width=6, height=22, fix_text_size=FALSE)
grid.arrange(p_boxplot_phylum_dataset, p_boxplot_class_dataset, 
             p_boxplot_family_dataset, nrow = 3)
dev.off()

p_boxplot_phylum_dataset <- p_boxplot_phylum_dataset +
  theme(legend.position = "none",axis.title.y = element_blank(),
        axis.text.y  = element_blank())

p_boxplot_class_dataset <- p_boxplot_class_dataset +
  theme(legend.position = "none",axis.title.y = element_blank(),
        axis.text.y  = element_blank())

p_boxplot_family_dataset <-p_boxplot_family_dataset +
  theme(legend.position = "none",axis.title.y = element_blank(),
        axis.text.y  = element_blank())

svglite(filename="viral_barplots_SI_no_legend.svg", width=4, height=22, fix_text_size=FALSE)
grid.arrange(p_boxplot_phylum_dataset, p_boxplot_class_dataset, 
             p_boxplot_family_dataset, nrow = 3)
dev.off()

############
##Maaslin2##
############

#At the end I am using Maaslin2 because I can use normalised and transformed data with this!!!! Whereas I couldn't not with others.
#Results also make sense!
#ANCOMBC is the most powerful and the one I wanted to use at the beginning but not enough taxa to run it (< 50) so I had to find an anternative.. and better this way so then I can use TPMs!

#create metadata file
metadata <- data.frame(vb = c("b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v"),
                       glacier = c("r", "r", "r", "r", "r", "m", "m", "m", "m", "m", "m", "r", "r", "r", "r", "r","m", "m", "m", "m", "m", "m"))
row.names(metadata) <- colnames(class_contig_sum)[-1]

###################
#Maaslin2 - phylum#
###################

dt_phylum <- as.data.frame(phylum_contig_sum[,-1])
row.names(dt_phylum) <- phylum_contig_sum$phylum

maaslin_phy <- Maaslin2(
  input_data = t(sqrt(dt_phylum)), #transpose: samples = rows
  input_metadata = metadata,
  output = "maaslin_phy_v", #output folder
  fixed_effects = c("vb", "glacier"), #test both variables
  normalization = "NONE", #none because already TPM
  transform = "NONE", #none because already square root
  analysis_method = "LM", #linear model
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  max_significance = 0.05)

#significant taxa
res_phy <- read.table("maaslin_phy_v/significant_results.tsv", header = TRUE, sep = "\t")
res_phy

#intersect with taxa in plots
intersect(res_phy[res_phy$metadata == "vb",]$feature, unique(phylum_contig_long_dt_top$phylum))
#[1] "Pisuviricota"       "Artverviricota"     "Nucleocytoviricota" "Uroviricota" 

##################
#Maaslin2 - class#
##################

dt_class <- as.data.frame(class_contig_sum[,-1])
row.names(dt_class) <- class_contig_sum$class

maaslin_cla <- Maaslin2(
  input_data = t(sqrt(dt_class)), #transpose: samples = rows
  input_metadata = metadata,
  output = "maaslin_cla_v", #output folder
  fixed_effects = c("vb", "glacier"), #test both variables
  normalization = "NONE", #none because already TPM
  transform = "NONE", #none because already square root
  analysis_method = "LM", #linear model
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  max_significance = 0.05)

#significant taxa
res_cla <- read.table("maaslin_cla_v/significant_results.tsv", header = TRUE, sep = "\t")
res_cla

#intersect with taxa in plots
intersect(res_cla[res_cla$metadata == "vb",]$feature, unique(class_contig_long_dt_top$class))
#[1] "Pisoniviricetes"  "Stelpaviricetes"  "Caudoviricetes"   "Megaviricetes"   
#[5] "Polintoviricetes" "Revtraviricetes"  "Maveriviricetes"  "Tectiliviricetes"

###################
#Maaslin2 - family#
###################

dt_family <- as.data.frame(family_contig_sum[,-1])
row.names(dt_family) <- family_contig_sum$family

maaslin_fam <- Maaslin2(
  input_data = t(sqrt(dt_family)), #transpose: samples = rows
  input_metadata = metadata,
  output = "maaslin_fam_v", #output folder
  fixed_effects = c("vb", "glacier"), #test both variables
  normalization = "NONE", #none because already TPM
  transform = "NONE", #none because already square root
  analysis_method = "LM", #linear model
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  max_significance = 0.05)

#significant taxa
res_fam <- read.table("maaslin_fam_v/significant_results.tsv", header = TRUE, sep = "\t")
res_fam

#intersect with taxa in plots
intersect(res_fam[res_fam$metadata == "vb",]$feature, unique(family_contig_long_dt_top$family))
#[1] "Phycodnaviridae" "Paulinoviridae"  "Mimiviridae"     "Straboviridae"   "Caulimoviridae" 