######################
####Viral taxonomy####
######################

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(gridExtra)
#library(Maaslin2)

setwd("~/Desktop/R_icevirome/v4")

##############
#votu count 0#
##############

#Get number of reads to vOTUs

#get count information
votu_count_0 <- read.table("vOTU_abundance_0.txt", sep = "\t", header = T, row.names = 1)

#get only columns with tpm counts
votu_tpm0_0 <- votu_count_0[, grepl(".Read.Count", colnames(votu_count_0))]
#format sample names (columns)
colnames(votu_tpm0_0) <- gsub("^a\\.|_si\\.Read.Count$", "", colnames(votu_tpm0_0))
dim(votu_tpm0_0)
#[1] 7474   23
votu_tpm_0 <- as.data.frame(cbind(contig = row.names(votu_tpm0_0), votu_tpm0_0))[,-1]
dim(votu_tpm_0)
#[1] 7474   23

#number of reads mapped to each sample
colSums(votu_tpm_0)
#BM1a     BM1b      BM2      BM3      BM4      BM5      BR1      BR2      BR3      BR4      BR5     VM1a     VM1b 
#1579150  1356760  1478216  2036580  2050750  1638652  2103026  3457774  2150207   604282  2119921  1333601  2915691 
#VM3a     VM3b      VM4      VM5      VR1      VR2      VR3      VR4      VR5   merged 
#3109799  4963609  3109044  1859132  1601955  2703011  3403853  1260002  1754455 48589470

############
#votu count#
############

#get count information
votu_count <- read.table("vOTU_abundance_10.txt", sep = "\t", header = T, row.names = 1)

#get only columns with tpm counts
votu_tpm0 <- votu_count[, grepl("TPM", colnames(votu_count))]
#format sample names (columns)
colnames(votu_tpm0) <- gsub("^a\\.|_si\\.TPM$", "", colnames(votu_tpm0))
dim(votu_tpm0)
#[1] 7474   23
votu_tpm <- as.data.frame(cbind(contig = row.names(votu_tpm0), votu_tpm0))
dim(votu_tpm)
#[1] 7474   24

#add mean and and sd
votu_tpm_m_sd <- votu_tpm[,-24] %>%
  rowwise() %>%
  mutate(
    mean = mean(c_across(-contig), na.rm = TRUE),
    sd = sd(c_across(-contig), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  as.data.frame()

#save tpm counts
write.csv(votu_tpm_m_sd, "votu_tpm_m_sd.csv", row.names = FALSE)

##########
#Taxonomy#
##########

#get data from file
taxonomy0 <- read.table("votu_taxonomy_clean.txt", row.names = 1)
colnames(taxonomy0) <- c("domain", "phylum", "class", "order", "family")
dim(taxonomy0)
#[1] 7474    5
taxonomy <- taxonomy0[row.names(votu_tpm),]
dim(taxonomy)
#[1] 7474    5

table(taxonomy$phylum)

#format taxonomic names
taxonomy$domain <-gsub("d__","",as.character(taxonomy$domain))
taxonomy$phylum <-gsub("p__","",as.character(taxonomy$phylum))
taxonomy$class <-gsub("c__","",as.character(taxonomy$class))
taxonomy$order <-gsub("o__","",as.character(taxonomy$order))
taxonomy$family <-gsub("f__","",as.character(taxonomy$family))

#format
taxonomy_v <- cbind(contig = row.names(taxonomy), taxonomy)

#add contig information to taxonomy
taxonomy_contig <- merge(votu_tpm, taxonomy_v, by = "contig", all = FALSE)

########
#Tables#
########

rel_abundance_table <- function(df, column, filename) {
  counts <- table(df[[column]])
  rel_abund <- prop.table(counts)
  
  out <- data.frame(
    taxon = names(counts),
    count = as.integer(counts),
    relative_abundance = as.numeric(rel_abund)*100
  )
  
  write.csv(out, filename, row.names = FALSE)
  return(out)
}

phylum_df <- rel_abundance_table(taxonomy, "phylum", "phylum_vout_RA.csv")
class_df  <- rel_abundance_table(taxonomy, "class",  "class_vout_RA.csv")
family_df <- rel_abundance_table(taxonomy, "family", "family_vout_RA.csv")

#################
#alluvional plot#
#################

#get only phlyum
phylum_contig_sum <- taxonomy_contig %>%
  group_by(phylum) %>%
  filter(phylum != "Unclassified") %>%
  dplyr::summarize(sum_merged = sum(as.numeric(merged))) #%>%
  #filter(sum_merged > 1000)
phylum_contig_sum$phylum

phylum_to_family_sums <- taxonomy_contig %>%
  filter(phylum %in% phylum_contig_sum$phylum) %>%
  group_by(phylum, class, family) %>%
  summarise(sum_merged = sum(as.numeric(merged)), .groups = "drop")
phylum_to_family_sums

#plot
p_viruses <- ggplot(phylum_to_family_sums,
                       aes(axis1 = phylum, axis2 = class, axis3 = family, y = sum_merged)) +
  scale_x_discrete(limits = c("Phylum", "Class", "Family"),
                   expand = c(.05, .05)) +
  geom_alluvium(aes(fill = phylum), width = 1/2, knot.pos = 0.1, alpha = .8) +   # thinner flows
  geom_stratum(width = 1/2, fill = "grey95", colour = "grey50", alpha = 0.1) + # wider boxes
  geom_text(stat = "stratum", aes(label = ifelse(
    after_stat(count) > 5000,
    as.character(after_stat(stratum)), "")),
    size = 3.5, color = "black") +
  scale_fill_manual(values = c("Artverviricota" = "lightskyblue1", "Hofneiviricota" = "lightpink1", "Nucleocytoviricota" = "#8da0cb",
                               "Preplasmiviricota" = "darkseagreen2","Uroviricota" = "brown1")) +
  labs(y = "TPM", x = NULL, fill = "Phylum") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.title = element_text(color = "black"),      # axis titles
        axis.text = element_text(color = "black"),       # tick labels
        legend.title = element_text(color = "black"),    # legend title
        legend.text = element_text(color = "black"),     # legend labels
        plot.title = element_text(color = "black"))
p_viruses
#width must be the same, and the higher the wider the box with labels

#save plot
svglite(filename="alluvional_plot_viruses.svg", width=6.5, height=12, fix_text_size=FALSE)
p_viruses
dev.off()

explore_class <- taxonomy_contig %>%
  group_by(class) %>%
  filter(class != "Unclassified") %>%
  dplyr::summarize(sum_merged = sum(as.numeric(merged))) %>%
  filter(sum_merged > 1000)
explore_class$class

explore_family <- taxonomy_contig %>%
  group_by(family) %>%
  filter(family != "Unclassified") %>%
  dplyr::summarize(sum_merged = sum(as.numeric(merged))) %>%
  filter(sum_merged > 1000)
explore_family$family

#########
#barplot#
#########

#Format samples
colnames(taxonomy_contig) <- c("contig",
                               "FI-M-1a", "FI-M-1b", "FI-M-2", "FI-M-3", "FI-M-4", "FI-M-5",
                               "FI-R-1", "FI-R-2", "FI-R-3", "FI-R-4", "FI-R-5",
                               "FII-M-1a", "FII-M-1b", "FII-M-3a", "FII-M-3b", "FII-M-4", "FII-M-5",
                               "FII-R-1", "FII-R-2", "FII-R-3", "FII-R-4", "FII-R-5",
                               "merged","domain","phylum","class", "order", "family")

#Reshape to long format (all sample columns)
df_long <- taxonomy_contig %>%
  pivot_longer(
    cols = starts_with("FI"), # all TPM/sample columns
    names_to = "Sample",
    values_to = "TPM"
  )

#Aggregate by phylum per sample
df_agg <- df_long %>%
  group_by(Sample, phylum) %>%
  summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

#Check which phyla
unique(df_agg$phylum)
#[1] "Artverviricota"     "Hofneiviricota"     "Negarnaviricota"    "Nucleocytoviricota" "Peploviricota"      "Pisuviricota"       "Preplasmiviricota" 
#[8] "Produgelaviricota"  "Unclassified"       "Uroviricota" 

#Order taxa
unique(df_agg$phylum)
df_agg$phylum <- factor(df_agg$phylum, levels  = c("Unclassified", "Uroviricota", "Produgelaviricota", "Preplasmiviricota", "Pisuviricota", 
                                                   "Peploviricota", "Nucleocytoviricota", "Negarnaviricota","Hofneiviricota", "Artverviricota"))

#Stacked barplot by phylum with your colors
p_barplot_phylum <- ggplot(df_agg, aes(x = Sample, y = TPM, fill = phylum)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black")
  ) +
  labs(x = "Sample", y = "TPM", fill = "Phylum") +
  scale_fill_manual(values = c(
    "Artverviricota" = "lightskyblue1",
    "Hofneiviricota" = "lightpink1",
    "Nucleocytoviricota" = "#8da0cb",
    "Preplasmiviricota" = "darkseagreen2",
    "Uroviricota" = "brown1",
    "Peploviricota" = "#e6ab02",
    "Unclassified" = "grey70",
    "Negarnaviricota" = "#17becf",
    "Pisuviricota" = "darkgreen",
    "Produgelaviricota" = "#ff7f0e"
  ))
p_barplot_phylum

#Aggregate by phylum per sample
df_agg <- df_long %>%
  group_by(Sample, class) %>%
  summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

#Check which class
unique(df_agg$class)
#[1] "Ainoaviricetes"   "Caudoviricetes"      "Faserviricetes"   "Herviviricetes"   "Maveriviricetes"  "Megaviricetes"    "Naldaviricetes"  
#[9]   "Pokkesviricetes"  "Polintoviricetes" "Revtraviricetes"    "Tectiliviricetes" "Unclassified"

#Order taxa
df_agg$class <- factor(df_agg$class, levels  = rev(unique(df_agg$class)))

class_colors <- c(
  "Alsuviricetes"     = "#1f77b4",  # blue #
  "Caudoviricetes"    = "#ff7f0e",  # orange #
  "Ellioviricetes"    = "yellow",
  "Faserviricetes"    = "#2ca02c",  # green #
  "Herviviricetes"    = "salmon",  # red #
  "Maveriviricetes"   = "#9467bd",  # purple #
  "Megaviricetes"     = "#8c564b",  # brown #
  "Naldaviricetes"    = "#e377c2",  # pink #
  "Pisoniviricetes"   = "green",
  "Pokkesviricetes"   = "#7f7f7f",  # grey #
  "Polintoviricetes"  = "#bcbd22",  # olive #
  "Revtraviricetes"   = "#17becf",  # cyan #
  "Stelpaviricetes"   = "red",
  "Tectiliviricetes"  = "#aec7e8",  # light blue
  "Unclassified"      = "grey70"    # light grey
)

#Stacked barplot by phylum with your colors
p_barplot_class <- ggplot(df_agg, aes(x = Sample, y = TPM, fill = class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = class_colors) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black")
  ) +
  labs(x = "Sample", y = "TPM", fill = "Class")
p_barplot_class

svglite(filename="barplot_viruses.svg", width=7, height=12, fix_text_size=FALSE)
grid.arrange(p_barplot_phylum, p_barplot_class, ncol = 1)
dev.off()