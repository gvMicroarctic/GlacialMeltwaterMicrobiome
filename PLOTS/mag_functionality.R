#########################
####MAG functionality####
#########################

setwd("~/Desktop/R_icevirome/v3")

#/storage/varliero/icevirome_prj/MAG/kofamscan/phylophlan/db.tre_formatted.treefile
#/storage/varliero/icevirome_prj/MAG/kofamscan/keggdecoder.txt

library(ape)
library(dendextend)
library(dplyr)
library(svglite)
library(circlize)
library(ComplexHeatmap)
library(ape)

##############
#Get the tree#
##############

# Load and filter specifics
specifics <- read.table("specifics_corrected_drep_new.txt", sep = "\t", row.names = 1, header = TRUE)
#specifics_combined <- specifics %>% mutate(combined = paste(Mag, class, family, sep = "_"))
specifics_combined <- specifics %>% mutate(combined = paste(Mag, family, sep = ", "))

# Filter near-complete MAGs
nc <- specifics_combined %>%
  filter(Completeness >= 70, Contamination <= 10)
dim(nc)
#[1] 87 18

# Generate enough distinct colors for all phyla
phyla <- sort(unique(nc$phylum))
colors <- c(
  "Bacteroidota" = "blue",
  "Patescibacteriota" = "yellow",
  "Actinomycetota" = "orange",
  "Pseudomonadota" = "mistyrose",
  "Acidobacteriota" = "#cab2d6",
  "Vulcanimicrobiota" = "brown",
  "Deinococcota" = "darkgoldenrod4",
  "Cyanobacteriota" = "#1b9e77",
  "Chlamydiota" = "#666666",
  "Bdellovibrionota" = "aquamarine1",
  "Armatimonadota" = "cornflowerblue",
  "Chloroflexota" = "#b2df8a",
  "Planctomycetota" = "#e7298a"
)

# Tree 2: From KEGG abundance
abundance <- as.matrix(read.table("keggdecoder.txt", sep = "\t", row.names = 1, header = TRUE, check.names = FALSE))
row.names(abundance) <- gsub("\\.", "_", row.names(abundance))
abundance_nc <- abundance[row.names(nc), ]
dist_mat <- dist(abundance_nc, method = "euclidean")
hc_tree <- hclust(dist_mat, method = "average")
tree <- as.dendrogram(hc_tree)

# Get dots colors at tree leaves
phylum_lookup <- setNames(nc$phylum, row.names(nc))
line_colors <- colors[phylum_lookup[labels(tree)]]
tip_colors <- colors[phylum_lookup[labels(tree)]]

# Replace labels with combined names
tree_labels <- labels(tree)
name <- specifics_combined$combined[match(tree_labels, row.names(specifics_combined))]
labels(tree) <- name

# Add dots on leaves
tree <- set(tree, "leaves_pch", 19)
tree <- set(tree, "leaves_cex", 1.5)
tree <- set(tree, "leaves_col", tip_colors)

# Plot
plot(tree, horiz = TRUE)
svglite(filename="tree.svg", width=9, height=13, fix_text_size=FALSE)
plot(tree, horiz = TRUE)
dev.off()

#########################
#Get functional profiles#
#########################

# Tree: From KEGG abundance
abundance <- as.matrix(read.table("keggdecoder.txt", sep = "\t", row.names = 1, header = TRUE, check.names = FALSE))
row.names(abundance) <- gsub("\\.", "_", row.names(abundance))
abundance_nc <- abundance[row.names(nc), ]

keep <- read.table("keggdecoder_to_keep.txt", sep = "\t", header = T) ###start from here - only patterns and remove

abundance_nc_keep0 <- abundance_nc[,keep$Function]
dim(abundance_nc_keep0)
abundance_nc_keep1 <- abundance_nc_keep0[,colSums(abundance_nc_keep0) > 0]
dim(abundance_nc_keep1)

#order by tree labels
#abundance_nc_keep <- abundance_nc_keep1
abundance_nc_keep <- abundance_nc_keep1[rev(tree_labels),]

#format data for plot
t1 <- keep %>% filter(category %in% "core metabolism")
t2 <- keep %>% filter(category %in% "o2 cytochromes")
t3 <- keep %>% filter(category %in% "carbon fixation")

t4 <- keep %>% filter(category %in% "carbon degradation")
t5 <- keep %>% filter(category %in% "nitrogen cycle")
t6 <- keep %>% filter(category %in% "sufur cycle")

t7 <- keep %>% filter(category %in% "hydrogen cycle")
t8 <- keep %>% filter(category %in% "vitamins and transporters")
t9 <- keep %>% filter(category %in% "photosynthesis")

t10 <- keep %>% filter(category %in% "Fermentation")
t11 <- keep %>% filter(category %in% "competence")
t12 <- keep %>% filter(category %in% "Amino_acid")

#plot

# ht1 = Heatmap(as.matrix(abundance_nc_keep[,intersect(t1$Function, colnames(abundance_nc_keep))]), name = "core metabolism", 
#               cluster_rows = FALSE, col = colorRamp2(c(0, 1), c("white", "red")))
# ht10 = Heatmap(as.matrix(abundance_nc_keep[,intersect(t10$Function, colnames(abundance_nc_keep))]), name = "Fermentation", 
#                cluster_rows = FALSE, col = colorRamp2(c(0, 1), c("white", "salmon")))
# ht11 = Heatmap(as.matrix(abundance_nc_keep[,intersect(t11$Function, colnames(abundance_nc_keep))]), name = "competence", 
#                cluster_rows = FALSE, col = colorRamp2(c(0, 1), c("white", "brown4")))
# ht8 = Heatmap(as.matrix(abundance_nc_keep[,intersect(t8$Function, colnames(abundance_nc_keep))]), name = "vitamins and transporters", 
#               cluster_rows = FALSE, col = colorRamp2(c(0, 1), c("white", "brown4")))
# ht2 = Heatmap(as.matrix(abundance_nc_keep[,intersect(t2$Function, colnames(abundance_nc_keep))]), name = "o2 cytochromes", 
#               cluster_rows = FALSE, col = colorRamp2(c(0, 1), c("white", "green")))

# ht3 = Heatmap(as.matrix(abundance_nc_keep[,intersect(t3$Function, colnames(abundance_nc_keep))]), name = "carbon fixation", 
#               col = colorRamp2(c(0, 1), c("white", "blueviolet")), show_heatmap_legend = FALSE)
# ht4 = Heatmap(as.matrix(abundance_nc_keep[,intersect(t4$Function, colnames(abundance_nc_keep))]), name = "carbon degradation", 
#               col = colorRamp2(c(0, 1), c("white", "darkblue")), show_heatmap_legend = FALSE)
# ht9 = Heatmap(as.matrix(abundance_nc_keep[,intersect(t9$Function, colnames(abundance_nc_keep))]), name = "photosynthesis", 
#               col = colorRamp2(c(0, 1), c("white", "darkgreen")), show_heatmap_legend = FALSE)
# 
# ht5 = Heatmap(as.matrix(abundance_nc_keep[,intersect(t5$Function, colnames(abundance_nc_keep))]), name = "nitrogen cycle", 
#               col = colorRamp2(c(0, 1), c("white", "orange")), show_heatmap_legend = FALSE)
# ht6 = Heatmap(as.matrix(abundance_nc_keep[,intersect(t6$Function, colnames(abundance_nc_keep))]), name = "sufur cycle", 
#               col = colorRamp2(c(0, 1), c("white", "darkred")), show_heatmap_legend = FALSE)
# ht7 = Heatmap(as.matrix(abundance_nc_keep[,intersect(t7$Function, colnames(abundance_nc_keep))]), name = "hydrogen cycle", 
#               col = colorRamp2(c(0, 1), c("white", "darkgoldenrod4")), show_heatmap_legend = FALSE)
# 
# ht12 = Heatmap(as.matrix(abundance_nc_keep[,intersect(t12$Function, colnames(abundance_nc_keep))]), name = "Amino_acid", 
#               col = colorRamp2(c(0, 1), c("white", "black")), show_heatmap_legend = FALSE)

ht3 = Heatmap(as.matrix(abundance_nc_keep[,intersect(t3$Function, colnames(abundance_nc_keep))]), name = "carbon fixation", 
              cluster_rows = FALSE, col = colorRamp2(c(0, 1), c("white", "black")), show_heatmap_legend = FALSE)
ht4 = Heatmap(as.matrix(abundance_nc_keep[,intersect(t4$Function, colnames(abundance_nc_keep))]), name = "carbon degradation", 
              cluster_rows = FALSE, col = colorRamp2(c(0, 1), c("white", "black")), show_heatmap_legend = FALSE)
ht9 = Heatmap(as.matrix(abundance_nc_keep[,intersect(t9$Function, colnames(abundance_nc_keep))]), name = "photosynthesis", 
              cluster_rows = FALSE, col = colorRamp2(c(0, 1), c("white", "black")), show_heatmap_legend = FALSE)

ht5 = Heatmap(as.matrix(abundance_nc_keep[,intersect(t5$Function, colnames(abundance_nc_keep))]), name = "nitrogen cycle", 
              cluster_rows = FALSE, col = colorRamp2(c(0, 1), c("white", "black")), show_heatmap_legend = FALSE)
ht6 = Heatmap(as.matrix(abundance_nc_keep[,intersect(t6$Function, colnames(abundance_nc_keep))]), name = "sufur cycle", 
              cluster_rows = FALSE, col = colorRamp2(c(0, 1), c("white", "black")), show_heatmap_legend = FALSE)
ht7 = Heatmap(as.matrix(abundance_nc_keep[,intersect(t7$Function, colnames(abundance_nc_keep))]), name = "hydrogen cycle", 
              cluster_rows = FALSE, col = colorRamp2(c(0, 1), c("white", "black")), show_heatmap_legend = FALSE)

ht12 = Heatmap(as.matrix(abundance_nc_keep[,intersect(t12$Function, colnames(abundance_nc_keep))]), name = "Amino_acid", 
               cluster_rows = FALSE, col = colorRamp2(c(0, 1), c("white", "black")), show_heatmap_legend = FALSE)

ht3 + ht4 + ht9 + ht5 + ht6 + ht7 + ht12

svglite(filename="mag_functions.svg", width=9, height=13, fix_text_size=FALSE)
ht3 + ht4 + ht9 + ht5 + ht6 + ht7 + ht12
dev.off()

#plot legend

Heatmap(as.matrix(abundance_nc_keep[,intersect(t6$Function, colnames(abundance_nc_keep))]), name = "sufur cycle", 
              cluster_rows = FALSE, col = colorRamp2(c(0, 1), c("white", "black")))
#png : 1500x700



