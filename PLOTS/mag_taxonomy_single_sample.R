##################
###MAG taxonomy###
##################

library(vegan)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(grid)  #needed for arrow plotting
library(Maaslin2)
library(tibble)
library(svglite)
library(Hmisc)
library(igraph)
library(reshape2)
library(compositions)
library(gtools)
library(ComplexHeatmap)

setwd("~/Desktop/R_icevirome/v3")

######################
#import MAG specifics#
######################

#import MAG specifics
specifics <- read.table("specifics_corrected_drep_new.txt", sep = "\t", row.names = 1, header = T)
dim(specifics)
#[1] 112 17

######################
#import MAG abundance#
######################

#read file
mag_tpm0 <- read.table("MAG_abundance_10.txt", sep = "\t", header = T, row.names = 1)[-1,]
dim(mag_tpm0)
#[1] 112 88
mag_tpm0
#I am using MAG_abundance_10.txt so then cleaner results

#get only columns with "TPM"
mag_tpm1 <- mag_tpm0[, grepl("TPM", colnames(mag_tpm0))]
mag_tpm1

#format sample names (columns)
colnames(mag_tpm1) <- gsub("^a\\.|_si\\.TPM$", "", colnames(mag_tpm1))

#format MAG names (rows)
rownames(mag_tpm1) <- sub("\\.(strict|orig|permissive)$", "", rownames(mag_tpm1))
mag_tpm1

#check that same MAGs
setdiff(row.names(mag_tpm1), row.names(specifics))
#character(0)

#get only medium and high quality MAGs
mag_tpm <- mag_tpm1[row.names(specifics),]
dim(mag_tpm)
#[1] 112  22
row.names(mag_tpm) == row.names(specifics)

#correct MAG names
row.names(mag_tpm) <- specifics$Mag

#split by dataset
mag_tpm_I <- mag_tpm[,c(1:11)]
mag_tpm_II <- mag_tpm[,c(12:22)]

#split by glacier
mag_tpm_m <- mag_tpm[,c(1:6, 12:17)]
mag_tpm_r <- mag_tpm[,c(7:11, 18:22)]

colnames(mag_tpm_I)
colnames(mag_tpm_II)
colnames(mag_tpm_m)
colnames(mag_tpm_r)

###########
#permanova#
###########

metadata <- data.frame(vb = c("b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b",
                              "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v"),
                       glacier = c("m", "m", "m", "m", "m", "m",
                                   "r", "r", "r", "r", "r",
                                   "m", "m", "m", "m", "m", "m",
                                   "r", "r", "r", "r", "r"))
row.names(metadata) <- colnames(mag_tpm)[1:22]

dist_matrix <- vegdist(t(sqrt(mag_tpm[1:22])), method = "bray")  # Bray-Curtis distance

permanova_result <- adonis2(dist_matrix ~ vb * glacier, data = metadata, permutations = 999, by = "terms")
permanova_result
#           Df SumOfSqs      R2       F Pr(>F)    
#vb          1  0.51967 0.23695 11.5997  0.001 ***
#glacier     1  0.69117 0.31515 15.4277  0.001 ***
#vb:glacier  1  0.17589 0.08020  3.9261  0.007 ** 
#Residual   18  0.80641 0.36770                   
#Total      21  2.19315 1.00000                                 

#############
#bubble plot#
#############

#get column with mags
df <- data.frame(mag_tpm) %>% rownames_to_column(var = "mag")

#get column with mag + class
row.names(data.frame(mag_tpm)) == specifics$Mag
df$mag_class <- paste(specifics$class, specifics$family, specifics$Mag)

# Change sample names
colnames(df) <- c("mag",
                  "FI-M-1a", "FI-M-1b", "FI-M-2", "FI-M-3", "FI-M-4", "FI-M-5",
                  "FI-R-1", "FI-R-2", "FI-R-3", "FI-R-4", "FI-R-5",
                  "FII-M-1a", "FII-M-1b", "FII-M-3a", "FII-M-3b", "FII-M-4", "FII-M-5",
                  "FII-R-1", "FII-R-2", "FII-R-3", "FII-R-4", "FII-R-5",
                  "merged","mag_class")

# Convert to long format and filter zeros
df_long <- df %>%
  pivot_longer(cols = -c(mag, mag_class, merged), names_to = "sample", values_to = "value") %>%
  filter(value > 0)

# Order MAGs (row names) alphanumerically within each phylum
specifics_ordered <- specifics %>%
  arrange(phylum, class, family) %>%                     # first sort by higher taxonomic ranks
  group_by(phylum, class, family) %>%                   # group by taxonomic hierarchy
  arrange(.by_group = TRUE) %>%                         # make sure grouping order is preserved
  slice(mixedorder(Mag)) %>%                             # natural alphanumeric sort of Mag column within each group
  ungroup()

# Add group column based on glacier
df_long <- df_long %>%
  mutate(group = ifelse(grepl("M", sample), "M", "R"))

# Re-factor samples and MAGs in sorted order
df_long$sample <- factor(df_long$sample, levels = c(
  "FI-M-1a", "FI-M-1b", "FI-M-2", "FI-M-3", "FI-M-4", "FI-M-5",
  "FII-M-1a", "FII-M-1b", "FII-M-3a", "FII-M-3b", "FII-M-4", "FII-M-5",
  "FI-R-1", "FI-R-2", "FI-R-3", "FI-R-4", "FI-R-5",
  "FII-R-1", "FII-R-2", "FII-R-3", "FII-R-4", "FII-R-5"))
df_long$mag <- factor(df_long$mag, levels =rev(specifics_ordered$Mag))

#plot
label_map <- setNames(df_long$mag_class, df_long$mag)
p_bubble <- ggplot(df_long, aes(x = sample, y = mag, size = value, color = group)) +
  geom_point(alpha = 0.6) +  # bubbles with transparency
  scale_size_continuous(range = c(0, 20), limits = c(1, 650000)) +
  scale_color_manual(values = c("M" = "darkblue", "R" = "darkred")) +  # custom colors
  scale_y_discrete(labels = label_map) +  # use mag_class for y-axis labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # rotates x-axis labels
  labs(x = "Sample", y = "MAG", size = "TPM counts", color = "Group")
p_bubble

#save
svglite(filename="mag_bubble_plot.svg", width=10.5, height=15, fix_text_size=FALSE)
p_bubble
dev.off()

###########################
#Maaslin2 - entire dataset#
###########################

maaslin <- Maaslin2(
  input_data = t(sqrt(mag_tpm[1:22])), #transpose: samples = rows
  input_metadata = metadata,
  output = "maaslin_mag_tpm", #output folder
  fixed_effects = c("vb", "glacier"), #test both variables
  normalization = "NONE", #none because already TPM
  transform = "NONE", #none because already square root
  analysis_method = "LM", #linear model
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  max_significance = 0.05)

#significant taxa
res <- read.table("maaslin_mag_tpm/significant_results.tsv", header = TRUE, sep = "\t")
res

#- -> more in i
#+ -> more in ii
nrow(res[res$metadata == "vb",])
#[1] 84
mag_i <- res %>%
  filter(coef < 0) %>%
  filter(metadata == "vb") %>%
  mutate(feature = gsub("^X", "", feature)) %>%
  pull(feature)
length(mag_i)
#[1] 40
mag_ii <- res %>%
  filter(coef > 0) %>%
  filter(metadata == "vb") %>%
  mutate(feature = gsub("^X", "", feature)) %>%
  pull(feature)
length(mag_ii)
#[1] 25

#- -> more in m
#+ -> more in r
nrow(res[res$metadata == "glacier",])
#[1] 71
mag_m <- res %>%
  filter(coef < 0) %>%
  filter(metadata == "glacier") %>%
  mutate(feature = gsub("^X", "", feature)) %>%
  pull(feature)
length(mag_m)
#[1] 17

mag_r <- res %>%
  filter(coef > 0) %>%
  filter(metadata == "glacier") %>%
  mutate(feature = gsub("^X", "", feature)) %>%
  pull(feature)
length(mag_r)
#[1] 54

######################
#Maaslin2 - dataset I#
######################

metadata_I <- data.frame(glacier = c("m", "m", "m", "m", "m", "m",
                                   "r", "r", "r", "r", "r"))
row.names(metadata_I) <- colnames(mag_tpm_I)

maaslin_I <- Maaslin2(
  input_data = t(sqrt(mag_tpm_I)), #transpose: samples = rows
  input_metadata = metadata_I,
  output = "maaslin_mag_tpm_I", #output folder
  fixed_effects = c("glacier"), #test both variables
  normalization = "NONE", #none because already TPM
  transform = "NONE", #none because already square root
  analysis_method = "LM", #linear model
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  max_significance = 0.05)

#significant taxa
res <- read.table("maaslin_mag_tpm_I/significant_results.tsv", header = TRUE, sep = "\t")
res

#- -> more in m
#+ -> more in r
nrow(res[res$metadata == "glacier",])
#[1] 68
mag_m_I <- res %>%
  filter(coef < 0) %>%
  filter(metadata == "glacier") %>%
  mutate(feature = gsub("^X", "", feature)) %>%
  pull(feature)
length(mag_m_I)
#[1] 10

mag_r_I <- res %>%
  filter(coef > 0) %>%
  filter(metadata == "glacier") %>%
  mutate(feature = gsub("^X", "", feature)) %>%
  pull(feature)
length(mag_r_I)
#[1] 58

#######################
#Maaslin2 - dataset II#
#######################

metadata_II <- data.frame(glacier = c("m", "m", "m", "m", "m", "m",
                                     "r", "r", "r", "r", "r"))
row.names(metadata_II) <- colnames(mag_tpm_II)

maaslin_II <- Maaslin2(
  input_data = t(sqrt(mag_tpm_II)), #transpose: samples = rows
  input_metadata = metadata_II,
  output = "maaslin_mag_tpm_II", #output folder
  fixed_effects = c("glacier"), #test both variables
  normalization = "NONE", #none because already TPM
  transform = "NONE", #none because already square root
  analysis_method = "LM", #linear model
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  max_significance = 0.05)

#significant taxa
res <- read.table("maaslin_mag_tpm_II/significant_results.tsv", header = TRUE, sep = "\t")
res

#- -> more in m
#+ -> more in r
nrow(res[res$metadata == "glacier",])
#[1] 31
mag_m_II <- res %>%
  filter(coef < 0) %>%
  filter(metadata == "glacier") %>%
  mutate(feature = gsub("^X", "", feature)) %>%
  pull(feature)
length(mag_m_II)
#[1] 8

mag_r_II <- res %>%
  filter(coef > 0) %>%
  filter(metadata == "glacier") %>%
  mutate(feature = gsub("^X", "", feature)) %>%
  pull(feature)
length(mag_r_II)
#[1] 23

#######################################
#heatmap with metadata for bubble plot#
#######################################

#all mags
mag_all <- specifics_ordered$Mag
mag_all

#differentially abundant in morterasch
mag_m_I
mag_m_II

#differentially abundant in rhone
mag_r_I
mag_r_II

#create the list of groups
mag_groups <- list(m_I = mag_m_I, m_II = mag_m_II,
                   r_I = mag_r_I, r_II = mag_r_II)

#create a presence/absence data frame
mag_matrix <- sapply(mag_groups, function(vec) as.integer(mag_all %in% vec))

#convert to data frame and add row names
mag_df <- as.data.frame(mag_matrix)
rownames(mag_df) <- mag_all
mag_df

#check if contrasts
c1 <- mag_df[,c("m_I", "r_II")]
c2 <- mag_df[,c("r_I", "m_II")]
c1[rowSums(c1) > 1,]
c2[rowSums(c2) > 1,]
#<0 rows> (or 0-length row.names)

#format diffentially abundant data.frame
da <- data.frame(cbind(m = rowSums(mag_df[,c(1,2)]), r = (rowSums(mag_df[,c(3,4)]))))
da
da$r[da$r == 1] <- 3
da$r[da$r == 2] <- 4
da
da_sum <- data.frame(cbind(da, sum = rowSums(da)))
da_sum

#get data.frame for phylum
phylum <- specifics_ordered %>%
  filter(Mag %in% row.names(da_sum)) %>%
  select(Mag, phylum)

#final dataset
final_hm <- cbind(da = da_sum$sum, phylum = phylum)
row.names(final_hm) <- row.names(da_sum)
final_hm

#heatmap
colors <- c("white", "#ADD8E6", "darkblue", "brown1", "darkred")
ht1 = Heatmap(as.matrix(final_hm[,1]),
        cluster_rows = FALSE,
        col = colors, # Pass your declared vector of 5 colors here
        show_heatmap_legend = TRUE) # Set to TRUE to verify the legend, then FALSE if desired

colors <- c("#cab2d6", "orange", "cornflowerblue", "blue", "aquamarine1", "#666666", "#b2df8a",
            "yellow", "#1b9e77", "darkgoldenrod4", "red", "#e7298a", "mistyrose", "brown")
ht2 = Heatmap(as.matrix(final_hm[,3]),
              cluster_rows = FALSE,
              col = colors, # Pass your declared vector of 5 colors here
              show_heatmap_legend = TRUE) # Set to TRUE to verify the legend, then FALSE if desired
ht1 + ht2

#save
svglite(filename="heatmap_da_phylum.svg", width=3, height=18, fix_text_size=FALSE)
ht1 + ht2
dev.off()

