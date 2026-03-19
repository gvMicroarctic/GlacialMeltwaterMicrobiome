#############################
####MAG phylogenetic tree####
#############################

setwd("~/Desktop/R_icevirome/v3")

library(ape)
library(ggtree)
library(tidytree)
library(dplyr)
library(svglite)
library(ggplot2)

#info on mags
specifics <- read.table("specifics_corrected_drep_new.txt", sep = "\t", row.names = 1, header = TRUE)

#read phylogenetic tree from the file
tree <- read.tree("accession_aa.tre")

#add node names
tree$node.label <- paste0("Node_", seq_len(tree$Nnode))

#write updated Newick
write.tree(tree, file = "accession_aa_with_nodes.tree")

#read new tree
tree <- read.tree("accession_aa_with_nodes.tree")

# Convert tree to tibble for easy node access
tree_df <- as_tibble(tree)

# Extract node IDs by label (internal nodes)
node_id_1  <- tree_df$node[tree_df$label == "Node_1787"]  # Bacteroidota
node_id_2  <- tree_df$node[tree_df$label == "Node_4384"]  # Patescibacteriota
node_id_3  <- tree_df$node[tree_df$label == "Node_5056"]    # Actinomycetota; but finish by had to colour in inkscape (Node_20 but then the entire tree was getting colored)
node_id_4  <- tree_df$node[tree_df$label == "Node_5127"]    # Actinomycetota; but finish by had to colour in inkscape (Node_20 but then the entire tree was getting colored)
node_id_5  <- tree_df$node[tree_df$label == "Node_466"]   # Pseudomonadota
node_id_6  <- tree_df$node[tree_df$label == "Node_283"]   # Acidobacteriota
node_id_7  <- tree_df$node[tree_df$label == "Node_4997"]  # Vulcanimicrobiota
node_id_8  <- tree_df$node[tree_df$label == "Node_3985"]  # Deinococcota
node_id_9  <- tree_df$node[tree_df$label == "Node_3231"]  # Cyanobacteriota
node_id_10  <- tree_df$node[tree_df$label == "Node_2769"]  # Chlamydiota
node_id_11 <- tree_df$node[tree_df$label == "Node_33"]    # Bdellovibrionota
node_id_12 <- tree_df$node[tree_df$label == "Node_260"] #Bdellovibrionota
node_id_13 <- tree_df$node[tree_df$label == "Node_2234"]  # Gemmatimonadota
node_id_14 <- tree_df$node[tree_df$label == "Node_4948"]  # Armatimonadota
node_id_15 <- tree_df$node[tree_df$label == "Node_3994"]  # Chloroflexota
node_id_16 <- tree_df$node[tree_df$label == "Node_2367"]  # Planctomycetota

# Combine into a single vector
nodes_to_group <- c(
  node_id_1, node_id_2, node_id_3, node_id_4, node_id_5,
  node_id_6, node_id_7, node_id_8, node_id_9, node_id_10,
  node_id_11, node_id_12, node_id_13, node_id_14, node_id_15, node_id_16
)

# Group clades by node
tree_grouped <- groupClade(tree, .node = nodes_to_group)

# Identify MAGs (bins)
bins <- tree_df$label[grepl("bin_", tree_df$label)]

# Plot tree
p <- ggtree(tree_grouped, layout = "equal_angle", aes(color = group)) +
  
  # Color by clade
  scale_color_manual(
    values=c("0" = "grey93", 
             "1" = "blue", "2" = "yellow", "3" = "orange", "4" = "orange", "5" = "mistyrose",
             "6" = "#cab2d6", "7" = "brown", "8" = "darkgoldenrod4", "9" = "#1b9e77", "10" = "#666666", 
             "11" = "aquamarine1", "12" = "aquamarine1", "13" = "red", "14" = "cornflowerblue", "15" = "#b2df8a", "16" = "#e7298a"), guide="none"
  ) +
  
  # Add points for MAGs
  geom_point(
    data = . %>% filter(grepl("bin_", label) & isTip),
    aes(x = x, y = y),
    color = "black",
    size = 2.5,
    shape = 8
  ) +
  
  # Add tree scale
  geom_treescale(x = -0.8, y = -0.4, width = 0.1, color = "black", linesize = 0.25, offset = 0.1) +
  
  # Remove tip labels
  theme(legend.position="none")
p

#save plot
svglite(filename="mag_phylogenetic_tree.svg", width=10, height=9, fix_text_size=FALSE)
p
dev.off()

#########################################
####Add TPM abundances to prokaryotes####
#########################################

#get count information
MAG_count <- read.table("MAG_abundance_10.txt", sep = "\t", header = T, row.names = 1)[-1,]

#get only columns with tpm counts
MAG_tpm0 <- MAG_count[, grepl("TPM", colnames(MAG_count))]
#format sample names (columns)
colnames(MAG_tpm0) <- gsub("^a\\.|_si\\.TPM$", "", colnames(MAG_tpm0))
dim(MAG_tpm0)
#[1] 780852     23
MAG_tpm <- as.data.frame(cbind(File = row.names(MAG_tpm0), MAG_tpm0))
dim(MAG_tpm)
#[1] 112     24

##########
#Taxonomy#
##########

specifics <- read.table("specifics_corrected_drep_new.txt", sep = "\t", header = TRUE)
dim(specifics)
#[1] 112  18

#add contig information to taxonomy
taxonomy_contig <- merge(MAG_tpm, specifics, by = "File", all = FALSE)

##########
###Plot###
##########

# Summarize merged counts per class
class_contig_sum_e <- taxonomy_contig %>%
  filter(class != "Unclassified") %>%
  group_by(class, phylum) %>%
  summarise(sum_merged = sum(as.numeric(merged)), .groups = "drop") %>%
  arrange(sum_merged)   #order by abundance

# Convert class to factor for plotting (preserves order)
class_contig_sum_e$class <- factor(
  class_contig_sum_e$class,
  levels = class_contig_sum_e$class
)

unique(class_contig_sum_e$class)

# Horizontal bar/segment plot with colored squares
p_barplot_class <- ggplot(
  class_contig_sum_e,
  aes(x = class, y = sum_merged)
) +
  geom_segment(
    aes(x = class, xend = class, y = 0, yend = sum_merged),
    color = "black"
  ) +
  geom_point(
    aes(fill = phylum),      # <<< THIS IS THE KEY CHANGE
    size = 6,
    shape = 21,
    color = "black"
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Acidobacteriota"      = "#cab2d6",
      "Actinomycetota"       = "orange",
      "Armatimonadota"       = "#1f78b4",
      "Bacteroidota"         = "blue",
      "Bdellovibrionota"     = "aquamarine1",
      "Chlamydiota"          = "#666666",
      "Chloroflexota"        = "#b2df8a",
      "Patescibacteriota"      = "yellow",
      "Cyanobacteriota"      = "#1b9e77",
      "Deinococcota"         = "darkgoldenrod4",
      "Gemmatimonadota"      = "red",
      "Planctomycetota"      = "#e7298a",
      "Pseudomonadota"       = "mistyrose",
      "Vulcanimicrobiota"    = "brown"
    )
  ) +
  labs(y = "TPM", x = NULL, fill = "Phylum") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

# Display plot
p_barplot_class

svglite(filename="barplot_prokaryotes.svg", width=4, height=8, fix_text_size=FALSE)
p_barplot_class
dev.off()
