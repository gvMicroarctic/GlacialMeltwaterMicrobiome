######################################
#votu-eukaryotic matches via proteins#
######################################

library(igraph)
library(dplyr)
library(readr)
library(svglite)

setwd("~/Desktop/R_icevirome/v3")

# -------------------------------
# Read data
# -------------------------------
data <- read.table("votu_euk_counts.tsv",
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

unique(data$contig_phylum)
unique(data$votu_phylum)

# Count number associated to votu-taxa and contig-taxa

data_votu <- data[!duplicated(data$votu), ]
table(data_votu$votu_phylum)

data_contig <- data[!duplicated(data$contig), ]
table(data_contig$contig_phylum)

# -------------------------------
# Create nodes
# -------------------------------
contig_nodes <- data %>%
  select(name = contig, phylum = contig_phylum) %>% distinct() %>%
  mutate(type = "contig")

votu_nodes <- data %>%
  select(name = votu, phylum = votu_phylum) %>% distinct() %>%
  mutate(type = "votu")

nodes <- bind_rows(contig_nodes, votu_nodes)

# -------------------------------
# Create edges with weights
# -------------------------------
edges <- data %>%
  select(contig, votu, count)

# -------------------------------
# Build igraph object
# -------------------------------
g <- graph_from_data_frame(d = edges[,1:2], vertices = nodes, directed = FALSE)

# Assign edge weights from count
E(g)$weight <- edges$count

# Optionally scale edge widths for plotting
# For example, scale counts to width range 1-5
E(g)$width <- scales::rescale(E(g)$weight, to = c(0.5, 15))

# -------------------------------
# Define colors for each phylum manually
# -------------------------------
phylum_colors <- c(
  "Ascomycota" = "red",
  "Basidiomycota" = "orange",
  "Chytridiomycota" = "brown",
  "Ciliophora" = "#ffd92f",
  "Tardigrada" = "moccasin",
  
  "Nucleocytoviricota" = "#8da0cb",
  "Uroviricota" = "brown1",
  "Artverviricota" = "lightskyblue1",
  
  "Unclassified" = "white"
)

V(g)$color <- phylum_colors[V(g)$phylum]

# -------------------------------
# Set shapes and sizes
# -------------------------------
V(g)$shape <- ifelse(V(g)$type == "contig", "circle", "square")
V(g)$size <- ifelse(V(g)$type == "contig", 7, 5)

# -------------------------------
# Plot network
# -------------------------------
plot(g,
     vertex.label = NA,
     edge.width = E(g)$width,   # use scaled counts
     layout = layout_with_fr)
legend("topleft",
       legend = names(phylum_colors),
       col = phylum_colors,
       pch = 19,
       pt.cex = 1.5,
       bty = "n",
       title = "Phylum")

svglite(filename="eukaryotes2votus_protein.svg", width=9, height=13, fix_text_size=FALSE)
plot(g,
     vertex.label = NA,
     edge.width = E(g)$width,   # use scaled counts
     layout = layout_with_fr)
legend("topleft",
       legend = names(phylum_colors),
       col = phylum_colors,
       pch = 19,
       pt.cex = 1.5,
       bty = "n",
       title = "Phylum")
dev.off()


svglite(filename="eukaryotes2votus_protein_labels.svg", width=9, height=13, fix_text_size=FALSE)
plot(g,
     edge.width = E(g)$width,   # use scaled counts
     layout = layout_with_fr)
legend("topleft",
       legend = names(phylum_colors),
       col = phylum_colors,
       pch = 19,
       pt.cex = 1.5,
       bty = "n",
       title = "Phylum")
dev.off()


###################################
#votu-eukaryotic matches via tRNAs#
###################################

# -------------------------------
# Read data
# -------------------------------
data <- read.table("votu_euk_trna_counts.tsv",
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

unique(data$contig_phylum)
unique(data$votu_phylum)

# Count number associated to votu-taxa and contig-taxa

data_votu <- data[!duplicated(data$votu), ]
table(data_votu$votu_phylum)

data_contig <- data[!duplicated(data$contig), ]
table(data_contig$contig_phylum)

# -------------------------------
# Create nodes
# -------------------------------
contig_nodes <- data %>%
  select(name = contig, phylum = contig_phylum) %>% distinct() %>%
  mutate(type = "contig")

votu_nodes <- data %>%
  select(name = votu, phylum = votu_phylum) %>% distinct() %>%
  mutate(type = "votu")

nodes <- bind_rows(contig_nodes, votu_nodes)

# -------------------------------
# Create edges with weights
# -------------------------------
edges <- data %>%
  select(contig, votu, count)

# -------------------------------
# Build igraph object
# -------------------------------
g <- graph_from_data_frame(d = edges[,1:2], vertices = nodes, directed = FALSE)

# Assign edge weights from count
E(g)$weight <- edges$count

# Optionally scale edge widths for plotting
# For example, scale counts to width range 1-5
E(g)$width <- scales::rescale(E(g)$weight, to = c(0.5, 15))

# -------------------------------
# Define colors for each phylum manually
# -------------------------------
phylum_colors <- c(
  "Chytridiomycota" = "brown", #
  "Ciliophora" = "#ffd92f", #
  "Chlorophyta" = "mediumspringgreen",
  "Streptophyta" = "#a6d854",
  
  "Cryptomycota" = "darkorange", #change to yellow for next time
  
  "Nucleocytoviricota" = "#8da0cb",
  "Uroviricota" = "brown1",
  
  "Unclassified" = "white"
)

V(g)$color <- phylum_colors[V(g)$phylum]

# -------------------------------
# Set shapes and sizes
# -------------------------------
V(g)$shape <- ifelse(V(g)$type == "contig", "circle", "square")
V(g)$size <- ifelse(V(g)$type == "contig", 7, 5)

# -------------------------------
# Plot network
# -------------------------------
plot(g,
     vertex.label = NA,
     edge.width = E(g)$width,   # use scaled counts
     layout = layout_with_fr)
legend("topleft",
       legend = names(phylum_colors),
       col = phylum_colors,
       pch = 19,
       pt.cex = 1.5,
       bty = "n",
       title = "Phylum")

svglite(filename="eukaryotes2votus_trna.svg", width=9, height=13, fix_text_size=FALSE)
plot(g,
     vertex.label = NA,
     edge.width = E(g)$width,   # use scaled counts
     layout = layout_with_fr)
legend("topleft",
       legend = names(phylum_colors),
       col = phylum_colors,
       pch = 19,
       pt.cex = 1.5,
       bty = "n",
       title = "Phylum")
dev.off()
