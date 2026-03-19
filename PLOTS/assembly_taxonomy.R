#####################################
####Viral and eukaryotic taxonomy####
#####################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(gridExtra)

setwd("~/Desktop/R_icevirome/v3")

################
#Assembly count#
################

#get count information
assembly_count <- read.table("assembly_abundance_10.txt", sep = "\t", header = T, row.names = 1)

#get only columns with tpm counts
assembly_tpm0 <- assembly_count[, grepl("TPM", colnames(assembly_count))]
#format sample names (columns)
colnames(assembly_tpm0) <- gsub("^a\\.|_si\\.TPM$", "", colnames(assembly_tpm0))
dim(assembly_tpm0)
#[1] 780852     23
assembly_tpm <- as.data.frame(cbind(contig = row.names(assembly_tpm0), assembly_tpm0))
dim(assembly_tpm)
#[1] 780852     24

##########
#Taxonomy#
##########

taxonomy <- read.csv("cat_taxonomy_only_official_formatted_trimmed.csv")
dim(taxonomy)
#[1] 780852      8

#number of unclassified contigs
table(taxonomy$domain)

#add contig information to taxonomy
taxonomy_contig <- merge(assembly_tpm, taxonomy, by = "contig", all = FALSE)

table(taxonomy_contig$domain)

###############################
###Prokaryotes vs eukaryotes###
###############################

#get only domain
domain_contig_sum <- taxonomy_contig %>%
  group_by(domain) %>%
  summarise(across(starts_with(c("V", "B")), ~ sum(as.numeric(.)), .names = "sum_{.col}"))
dim(domain_contig_sum)
write.csv(domain_contig_sum, "domain_contig_sum.csv")

##################################
###Eukaryotes - alluvional plot###
##################################

#get only phlum
phylum_contig_sum_e <- taxonomy_contig %>% filter(domain %in% c("Eukaryota")) %>%
  group_by(phylum) %>%
  filter(phylum != "Unclassified") %>%
  dplyr::summarize(sum_merged = sum(as.numeric(merged))) #%>%
  #filter(sum_merged > 100)
phylum_contig_sum_e$phylum

phylum_to_family_sums <- taxonomy_contig %>%
  filter(phylum %in% phylum_contig_sum_e$phylum) %>%
  group_by(phylum, class, family) %>%
  summarise(sum_merged = sum(as.numeric(merged)), .groups = "drop")
phylum_to_family_sums

#plot
p_eukaryotes <- ggplot(phylum_to_family_sums,
  aes(axis1 = phylum, axis2 = class, axis3 = family, y = sum_merged)) +
  scale_x_discrete(limits = c("Phylum", "Class", "Family"),
    expand = c(.05, .05)) +
  geom_alluvium(aes(fill = phylum), width = 1/3, alpha = .8, knot.pos = 0.1) +   # thinner flows
  geom_stratum(width = 1/3, fill = "grey95", colour = "grey50", alpha = 0.1) + # wider boxes
  geom_text(stat = "stratum", aes(label = ifelse(
      after_stat(count) > 100,
      as.character(after_stat(stratum)), "")),
    size = 3.5, color = "black") +
  scale_fill_manual(values = c("Arthropoda" = "grey","Ascomycota" = "mistyrose",
    "Basidiomycota" = "orange", "Chlorophyta" = "mediumspringgreen",
    "Chytridiomycota" = "burlywood3", "Ciliophora" = "#ffd92f",
    "Streptophyta" = "#a6d854", "Tardigrada" = "moccasin")) +
  labs(y = "TPM", x = NULL, fill = "Phylum") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
p_eukaryotes

#save plot
svglite(filename="alluvional_plot_eukaryotes.svg", width=6.5, height=12, fix_text_size=FALSE)
p_eukaryotes
dev.off()

##########################
###Eukaryotes - barplot###
##########################

# Prepare phylum summary
phylum_contig_sum_e <- taxonomy_contig %>%
  filter(domain == "Eukaryota") %>%
  filter(phylum != "Unclassified") %>%
  group_by(phylum) %>%
  summarise(sum_merged = sum(as.numeric(merged)), .groups = "drop") %>%
  filter(sum_merged > 100) %>%   # optional threshold
  arrange(desc(sum_merged))       # sort descending for abundance

# Convert phylum to factor for plotting (preserves order)
phylum_contig_sum_e <- phylum_contig_sum_e %>%
  mutate(phylum = factor(phylum, levels = rev(phylum)))

# Horizontal bar/segment plot with colored squares
p_barplot_phylum <- ggplot(phylum_contig_sum_e, aes(x = phylum, y = sum_merged)) +
  geom_segment(aes(x = phylum, xend = phylum, y = 0, yend = sum_merged), color = "black") +
  geom_point(aes(fill = phylum), size = 6, shape = 22, color = "black") +  # square with border
  coord_flip() +
  scale_fill_manual(values = c(
    "Ascomycota" = "red",
    "Basidiomycota" = "orange",
    "Chytridiomycota" = "brown",
    
    "Chlorophyta" = "mediumspringgreen",
    "Streptophyta" = "#a6d854",
    
    "Ciliophora" = "#ffd92f",
    
    "Arthropoda" = "grey",
    
    "Tardigrada" = "moccasin"
  )) +
  labs(y = "TPM", x = NULL, fill = "Phylum") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.title = element_text(color = "black"),      # axis titles
        axis.text = element_text(color = "black"),       # tick labels
        legend.title = element_text(color = "black"),    # legend title
        legend.text = element_text(color = "black"),     # legend labels
        plot.title = element_text(color = "black"))

# Display plot
p_barplot_phylum

# Summarize merged counts per class
class_contig_sum_e <- taxonomy_contig %>%
  filter(domain == "Eukaryota") %>%
  filter(class != "Unclassified") %>%
  group_by(class) %>%
  summarise(sum_merged = sum(as.numeric(merged)), .groups = "drop") %>%
  filter(sum_merged > 100) %>%          # optional threshold
  arrange(desc(sum_merged))             # sort descending for abundance

# Convert class to factor for plotting (preserves order)
class_contig_sum_e <- class_contig_sum_e %>%
  mutate(class = factor(class, levels = rev(class)))

unique(class_contig_sum_e$class)

# Horizontal bar/segment plot with colored squares
p_barplot_class <- ggplot(class_contig_sum_e, aes(x = class, y = sum_merged)) +
  geom_segment(aes(x = class, xend = class, y = 0, yend = sum_merged), color = "black") +
  geom_point(aes(fill = class), size = 6, shape = 22, color = "black") +  # square with border
  coord_flip() +
  # Optional: custom colors for top classes, adjust as needed
  scale_fill_manual(values = c(
    "Dothideomycetes" = "red",
    "Lecanoromycetes" = "red",
    "Leotiomycetes" = "red",
    "Taphrinomycetes" = "red",
    
    "Microbotryomycetes" = "orange",
    
    "Chytridiomycetes" = "brown",
    
    "Trebouxiophyceae" = "mediumspringgreen", #
    "Oligohymenophorea" = "#ffd92f", #
    "Collembola" = "grey", #
    "Eutardigrada" = "moccasin" #
  )) +
  labs(y = "TPM", x = NULL, fill = "Class") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    plot.title = element_text(color = "black")
  )

# Display plot
p_barplot_class

svglite(filename="barplot_eukaryotes.svg", width=8, height=4, fix_text_size=FALSE)
grid.arrange(p_barplot_phylum, p_barplot_class, ncol = 2)
dev.off()

#################
###Prokaryotes###
#################

#get only phlum
phylum_contig_sum_b <- taxonomy_contig %>% filter(domain %in% c("Bacteria")) %>%
  group_by(phylum) %>%
  filter(phylum != "Unclassified") %>%
  dplyr::summarize(sum_merged = sum(as.numeric(merged))) #%>%
#filter(sum_merged > 100)
phylum_contig_sum_b$phylum

phylum_to_family_sums_b <- taxonomy_contig %>%
  filter(phylum %in% phylum_contig_sum_b$phylum) %>%
  group_by(phylum, class, family) %>%
  summarise(sum_merged = sum(as.numeric(merged)), .groups = "drop")
phylum_to_family_sums_b

#########
#barplot#
#########

#Format samples
colnames(taxonomy_contig) <- c("contig",
                               "FI-M-1a", "FI-M-1b", "FI-M-2", "FI-M-3", "FI-M-4", "FI-M-5",
                               "FI-R-1", "FI-R-2", "FI-R-3", "FI-R-4", "FI-R-5",
                               "FII-M-1a", "FII-M-1b", "FII-M-3a", "FII-M-3b", "FII-M-4", "FII-M-5",
                               "FII-R-1", "FII-R-2", "FII-R-3", "FII-R-4", "FII-R-5",
                               "merged","domain","phylum","class", "order", "family", "genus", "species")

#Reshape to long format (all sample columns)
df_long <- taxonomy_contig %>% filter(domain %in% c("Bacteria")) %>%
  pivot_longer(
    cols = starts_with(c("FI", "merged")), # all TPM/sample columns
    names_to = "Sample",
    values_to = "TPM"
  )

#Aggregate by phylum per sample
df_agg <- df_long %>%
  group_by(Sample, phylum) %>%
  summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

#get first 20 phyla
top20 <- df_agg %>%
  filter(Sample == "merged") %>%
  arrange(desc(TPM)) %>%
  filter(phylum != "Unclassified") %>%
  slice_head(n = 20)
top20$phylum
df_agg_top20 <- df_agg %>% filter(phylum %in% top20$phylum)

#correct nomenclature (consistent with)
df_agg_top20 <- df_agg_top20 %>%
  mutate(phylum = ifelse(phylum == "Candidatus Saccharibacteria", "Patescibacteriota", phylum))

#Order taxa
taxa_vector <- c("Acidobacteriota", "Actinomycetota", "Armatimonadota", "Bacteroidota", "Bdellovibrionota", 
                 "Chlamydiota", "Chloroflexota", "Patescibacteriota", "Cyanobacteriota", "Deinococcota", 
                 "Gemmatimonadota", "Planctomycetota", "Pseudomonadota", "Vulcanimicrobiota", 
                 "Candidatus Paceibacterota", "Myxococcota", "Nitrospirota", "Thermodesulfobacteriota", "Verrucomicrobiota")
df_agg_top20$phylum <- factor(df_agg_top20$phylum, levels  = taxa_vector)

#########
#barplot#
#########

#Reshape to long format (all sample columns)
df_long <- taxonomy_contig %>% filter(domain %in% c("Bacteria")) %>%
  pivot_longer(
    cols = starts_with(c("FI", "merged")), # all TPM/sample columns
    names_to = "Sample",
    values_to = "TPM"
  )

#Aggregate by phylum per sample
df_agg <- df_long %>%
  group_by(Sample, phylum) %>%
  summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

#get first 20 phyla
top20 <- df_agg %>%
  filter(Sample == "merged") %>%
  arrange(desc(TPM)) %>%
  filter(phylum != "Unclassified") %>%
  slice_head(n = 20)
top20$phylum
df_agg_top20 <- df_agg %>% filter(phylum %in% top20$phylum)

#correct nomenclature (consistent with)
df_agg_top20 <- df_agg_top20 %>%
  mutate(phylum = ifelse(phylum == "Candidatus Saccharibacteria", "Patescibacteriota", phylum))

#order taxa
taxa_vector <- c("Acidobacteriota", "Actinomycetota", "Armatimonadota", "Bacteroidota", "Bdellovibrionota", 
                 "Chlamydiota", "Chloroflexota", "Patescibacteriota", "Cyanobacteriota", "Deinococcota", 
                 "Gemmatimonadota", "Planctomycetota", "Pseudomonadota", "Vulcanimicrobiota", "Bacillota",
                 "Candidatus Paceibacterota", "Myxococcota", "Nitrospirota", "Thermodesulfobacteriota", "Verrucomicrobiota")
df_agg_top20$phylum <- factor(df_agg_top20$phylum, levels  = taxa_vector)
unique(df_agg_top20$phylum)

#remove "merged"
df_agg_top20_trim <- df_agg_top20 %>% filter(!Sample %in% "merged")

#Colors
phylum_colors <- c(
  "Acidobacteriota" = "#cab2d6",
  "Actinomycetota" = "orange",
  "Armatimonadota" = "#1f78b4",
  "Bacteroidota" = "blue",
  "Bdellovibrionota" = "aquamarine1",
  "Chlamydiota" = "#666666",
  "Chloroflexota" = "#b2df8a",
  "Patescibacteriota" = "yellow",
  "Cyanobacteriota" = "#1b9e77",
  "Deinococcota" = "darkgoldenrod4",
  "Gemmatimonadota" = "red",
  "Planctomycetota" = "#e7298a",
  "Pseudomonadota" = "mistyrose",
  "Vulcanimicrobiota" = "brown",
  
  "Bacillota" = "grey", #different color this as well because interesting from iphop
  
  "Candidatus Paceibacterota" = "black",
  "Myxococcota" = "black",
  "Nitrospirota" = "black",
  "Thermodesulfobacteriota" = "black",
  "Verrucomicrobiota" = "black"
)

#Stacked barplot by phylum with your colors
p_barplot_phylum <- ggplot(df_agg_top20_trim, aes(x = Sample, y = TPM, fill = phylum)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black")
  ) + scale_fill_manual(values = phylum_colors) +
  labs(x = "Sample", y = "TPM", fill = "Phylum")
p_barplot_phylum

###########
###class###
###########

#Aggregate by phylum per sample
df_agg <- df_long %>%
  group_by(Sample, class) %>%
  summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

#get first 20 classes
top20 <- df_agg %>%
  filter(Sample == "merged") %>%
  arrange(desc(TPM)) %>%
  filter(class != "Unclassified") %>%
  slice_head(n = 20)
top20$class
df_agg_top20 <- df_agg %>% filter(class %in% top20$class)

#Order taxa
taxa_vector <- c("Terriglobia", "Actinomycetes", "Armatimonadia", "Chthonomonadia", "Bacteroidia", "Chitinophagia", 
                 "Cytophagia", "Saprospiria", "Sphingobacteriia", "Bdellovibrionia", "Oligoflexia", 
                 "Anaerolineae", "Cyanophyceae", "Deinococci", "Alphaproteobacteria", "Betaproteobacteria", 
                 "Gammaproteobacteria", "Desulfobacteria", "Flavobacteriia", "Ktedonobacteria")
df_agg_top20$class <- factor(df_agg_top20$class, levels  = taxa_vector)
unique(df_agg_top20$class)

#remove "merged"
df_agg_top20_trim <- df_agg_top20 %>% filter(!Sample %in% "merged")

#Colors
class_colors <- c(
  
  ##Acido
  "Terriglobia"        = "#cab2d6",
  
  ##Actinomycetes
  "Actinomycetes"      = "#ff7f00",  # orange
  
  ##Arm
  "Armatimonadia"      = "#1f78b4",  # turquoise
  "Chthonomonadia"     = "#a6cee3",  # mint
  
  ##Bactero
  "Bacteroidia"        = "blue",  # gray
  "Chitinophagia"      = "#6a3d9a",  # purple
  "Cytophagia"         = "darkblue",  # red
  "Saprospiria"        = "violet",  # light pink
  "Sphingobacteriia"   = "darkviolet",  # pink
  
  ##Bdello
  "Bdellovibrionia"    = "cyan",  # peach
  "Oligoflexia"        = "cyan4",  # light green
  
  ##Chloro
  "Anaerolineae"       = "darkgreen",  # yellow
  
  ##Cyanoba
  "Cyanophyceae"       = "#33a02c",  # blue
  
  ##Deino
  "Deinococci"         = "#b15928",  # brown
  
  ##Pseudo
  "Alphaproteobacteria"= "pink",  # light blue
  "Betaproteobacteria" = "#fb8072",  # green
  "Gammaproteobacteria"= "darkred",  # lavender
  
  ##Others
  "Desulfobacteria"    = "black",
  "Flavobacteriia"     = "black",
  "Ktedonobacteria"    = "black"
)

#Stacked barplot by phylum with your colors
p_barplot_class <- ggplot(df_agg_top20_trim, aes(x = Sample, y = TPM, fill = class)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black")
  ) + scale_fill_manual(values = class_colors) +
  labs(x = "Sample", y = "TPM", fill = "Class")
p_barplot_class

############
###family###
############

#Aggregate by phylum per sample
df_agg <- df_long %>%
  group_by(Sample, family) %>%
  summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

#get first 20 phyla
top20 <- df_agg %>%
  filter(Sample == "merged") %>%
  arrange(desc(TPM)) %>%
  filter(family != "Unclassified") %>%
  slice_head(n = 20)
top20$family
df_agg_top20 <- df_agg %>% filter(family %in% top20$family)

#Order taxa
taxa_vector <- c("Acidobacteriaceae", "Microbacteriaceae", "Capsulimonadaceae", "Chthonomonadaceae", "Chitinophagaceae", 
                 "Cytophagaceae", "Hymenobacteraceae", "Sphingobacteriaceae", "Chamaesiphonaceae", "Leptolyngbyaceae", 
                 "Deinococcaceae", "Acetobacteraceae", "Comamonadaceae", "Methylophilaceae", "Oxalobacteraceae", 
                 "Sphingomonadaceae", "Burkholderiaceae", "Candidatus Magnetomoraceae", "Flectobacillaceae", "Pseudanabaenaceae")

df_agg_top20$family <- factor(df_agg_top20$family, levels  = taxa_vector)
unique(df_agg_top20$family)

#remove "merged"
df_agg_top20_trim <- df_agg_top20 %>% filter(!Sample %in% "merged")

#Colors
family_colors <- c(
  
  # Acidobacteriota
  "Acidobacteriaceae"          = "#cab2d6",  # gray
  
  # Actinobacteria
  "Microbacteriaceae"          = "#ff7f00",  # orange
  
  # Armatimonadota
  "Capsulimonadaceae"          = "#1f78b4",  # teal
  "Chthonomonadaceae"          = "#a6cee3",  # turquoise
  
  # Bacteroidota
  "Chitinophagaceae"           = "blue",  # light green
  "Cytophagaceae"              = "darkblue",  # red
  "Hymenobacteraceae"          = "violet",  # green
  "Sphingobacteriaceae"        = "darkviolet",  # pink
  
  # Cyanobacteria
  "Chamaesiphonaceae"          = "#b2df8a",  # blue
  "Leptolyngbyaceae"           = "darkgreen",  # light blue
  
  # Deinococci
  "Deinococcaceae"             = "#b15928",  # brown
  
  # Proteobacteria
  "Acetobacteraceae"           = "#fb8072",  # coral
  "Comamonadaceae"             = "pink",  # lavender
  "Methylophilaceae"           = "#e31a1c",  # yellow
  "Oxalobacteraceae"           = "#fdbf6f",  # peach
  "Sphingomonadaceae"          = "darkred",  # purple
  
  # Others
  "Burkholderiaceae"           = "black",  # violet
  "Candidatus Magnetomoraceae" = "black",  # mint
  "Flectobacillaceae"          = "black",  # light pink
  "Pseudanabaenaceae"          = "black"   # light yellow
)

#Stacked barplot by phylum with your colors
p_barplot_family <- ggplot(df_agg_top20_trim, aes(x = Sample, y = TPM, fill = family)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black")
  ) + scale_fill_manual(values = family_colors) +
  labs(x = "Sample", y = "TPM", fill = "Family")
p_barplot_family

#save
svglite(filename="barplot_prokaryotic_assembly.svg", width=8, height=16, fix_text_size=FALSE)
grid.arrange(p_barplot_phylum, p_barplot_class, p_barplot_family, ncol = 1)
dev.off()

