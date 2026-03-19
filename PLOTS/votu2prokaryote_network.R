######################################################
###Rhone glacier ice core analyses - virus-bacteria###
######################################################

library(tidyr)
library(dplyr)
library(ggplot2)
library(svglite)
library(igraph)

setwd("~/Desktop/R_icevirome/v4")

#viral taxonomy
taxonomy_v <- read.table(text = gsub(";", "\t", readLines("votu_taxonomy_clean.txt")))
row.names(taxonomy_v) <- taxonomy_v[,1]
colnames(taxonomy_v) <- c("contig", "domain_v", "phylum_v", "class_v", "order_v", "family_v")
dim(taxonomy_v)
#[1] 7474    6

#they need to be curated _curated
#I removed "_X" lettering; and redundant matches if then I had them 
#I changed changed from "Gammaproteobacteria" to "Betaproteobacteria" for order Burkholderiales
#I also removed Order and Family because it was a mess to sort out

#virus - microbe association (curated file; updated taxonomy)
#I decided not to use genome file because iphop for genus level, genome not good (I think)
host_prediction <- read.table("Host_prediction_to_genus_curated.txt", header = T, check.names = FALSE, sep = "\t")
dim(host_prediction)
#[1] 388    7

#number of vOTUs associated to microorganisms
length(unique(host_prediction$contig))
#[1] 344

#merge information for viral and microbial taxonomy
taxonomy_host_prediction <- merge(taxonomy_v, host_prediction)
dim(taxonomy_host_prediction)
#[1] 388   12

##################################################
#number of potentially infected genera per family#
##################################################

#how many viral contigs associated to each phylum
host_prediction_count_phylum <- as.data.frame(table(host_prediction$phylum))
host_prediction_count_phylum
host_prediction_count_phylum[host_prediction_count_phylum$Freq >= 10,]
#Var1 Freq
#1   Acidobacteriota   14
#2    Actinomycetota   35
#5         Bacillota   43
#6      Bacteroidota   73
#7  Bdellovibrionota   11
#12  Cyanobacteriota   40
#19   Pseudomonadota  130

#save table for SI
write.table(host_prediction_count_phylum, "host_prediction_count_phylum_SI.txt")

##################
#Network - phylum#
##################

#data.frame for microbial phyla
nt_phylum0 <- taxonomy_host_prediction[,c("contig", "phylum_v", "phylum")]
dim(nt_phylum0)
#[1] 388   3

# keep only unique lines
nt_phylum <- nt_phylum0 %>%
  add_count(contig, phylum, name = "n_occurrences") %>%   # count occurrences
  distinct(contig, phylum, .keep_all = TRUE)
dim(nt_phylum)
#[1] 360   4

#remove phyla where only one/two vOTU because too chaotic
trim <- as.vector(host_prediction_count_phylum[host_prediction_count_phylum$Freq >= 3,] %>% select(Var1))
nt_phylum_trim <- nt_phylum %>% filter(phylum %in% trim$Var1)

#create network
network <- graph_from_data_frame(nt_phylum_trim[,c("contig", "phylum")], directed = FALSE)
str(network)

#get microbial and viral taxon information
mapped_phylum0 <- nt_phylum$phylum_v[
  match(V(network)$name, nt_phylum$contig)
]
mapped_phylum <- ifelse(is.na(mapped_phylum0), V(network)$name, mapped_phylum0)

#add microbial and viral taxa metadata to vertices
V(network)$phylum <- mapped_phylum
unique(V(network)$phylum)

#colors for taxa with 5 or more viruses
phylum_colors <- c(
  "Actinomycetota" = "orange", #in mags
  "Pseudomonadota" = "mistyrose", #in mags
  "Bacteroidota" = "blue", #in mags
  "Cyanobacteriota" = "#1b9e77", #in mags
  "Bdellovibrionota" = "aquamarine1", #in mags
  
  "Acidobacteriota" = "#cab2d6", #in mags
  "Patescibacteriota" = "yellow", #in mags
  "Chloroflexota" = "#b2df8a", #in mags
  "Bacillota" = "lightcyan", 
  "Methanobacteriota" = "darkred",
  
  "Asgardarchaeota" = "#e7298a",
  "Campylobacterota" = "grey",
  "Spirochaetota" = "purple",
  
  "p__Uroviricota" = "brown1",
  "p__Hofneiviricota" = "lightpink1",
  "p__Nucleocytoviricota" = "#8da0cb",
  "p__Artverviricota" = "lightskyblue1",
  
  "p__Unclassified" = "white"
  
  #low abundance (< 2)
  #"Vulcanimicrobiota" = "brown", #in mags
  #"Planctomycetota" = "#e7298a", #in mags
  #"Chlamydiota" = "#666666", #in mags
  #"Deinococcota" = "darkgoldenrod4", #in mags
  #"Caldisericota" = "black",
  #"Armatimonadota" = "cornflowerblue", #in mags
  #"Gemmatimonadota" = "red", #in mags
  #"Halobacteriota" = "gold",
  #"Thermoproteota" = "violet",
)

#shapes
phylum_shapes <- c(
  # microbial phyla → now **square**
  "Actinomycetota" = "square",
  "Pseudomonadota" = "square",
  "Bacteroidota" = "square",
  "Cyanobacteriota" = "square",
  "Planctomycetota" = "square",
  "Bacillota" = "square",
  "Methanobacteriota" = "square",
  "Bdellovibrionota" = "square",
  "Acidobacteriota" = "square",
  "Patescibacteriota" = "square",
  "Gemmatimonadota" = "square",
  "Armatimonadota" = "square",
  "Vulcanimicrobiota" = "square",
  "Chloroflexota" = "square",
  "Asgardarchaeota" = "square",
  "Caldisericota" = "square",
  "Thermoproteota" = "square",
  "Campylobacterota" = "square",
  "Spirochaetota" = "square",
  "Chlamydiota" = "square",
  "Deinococcota" = "square",
  "Halobacteriota" = "square",
  
  
  "Methanobacteriota" = "square",
  "Thermoplasmatota" = "square",
  "Verrucomicrobiota" = "square",
  "Myxococcota" = "square",
  "Hadarchaeota" = "square",
  
  # viral phyla → now **circle**
  "p__Uroviricota"        = "circle",
  "p__Hofneiviricota"     = "circle",
  "p__Nucleocytoviricota" = "circle",
  "p__Artverviricota"     = "circle",
  "p__Unclassified"       = "circle"
)

#sizes
s_m <- 7
s_v <- 5
phylum_sizes <- c(
  # microbial phyla → now **3**
  "Actinomycetota" = s_m,
  "Pseudomonadota" = s_m,
  "Bacteroidota" = s_m,
  "Cyanobacteriota" = s_m,
  "Planctomycetota" = s_m,
  "Bacillota" = s_m,
  "Methanobacteriota" = s_m,
  "Bdellovibrionota" = s_m,
  "Acidobacteriota" = s_m,
  "Patescibacteriota" = s_m,
  "Gemmatimonadota" = s_m,
  "Armatimonadota" = s_m,
  "Vulcanimicrobiota" = s_m,
  "Chloroflexota" = s_m,
  "Asgardarchaeota" = s_m,
  "Caldisericota" = s_m,
  "Thermoproteota" = s_m,
  "Campylobacterota" = s_m,
  "Spirochaetota" = s_m,
  "Chlamydiota" = s_m,
  "Deinococcota" = s_m,
  "Halobacteriota" = s_m,
  
  "Methanobacteriota" = s_m,
  "Thermoplasmatota" = s_m,
  "Verrucomicrobiota" = s_m,
  "Myxococcota" = s_m,
  "Hadarchaeota" = s_m,
  
  # viral phyla → now **20**
  "p__Uroviricota"        = s_v,
  "p__Hofneiviricota"     = s_v,
  "p__Nucleocytoviricota" = s_v,
  "p__Artverviricota"     = s_v,
  "p__Unclassified"       = s_v
)

#assign colors, shapes and sizes to network
V(network)$color <- phylum_colors[V(network)$phylum]
V(network)$shape <- phylum_shapes[V(network)$phylum]
V(network)$size <- phylum_sizes[V(network)$phylum]

# plot the network
plot(network,
     vertex.label = NA)
legend("topleft",
       legend = names(phylum_colors),
       col    = phylum_colors,
       pch    = 19,
       pt.cex = 1.5,
       bty    = "n",
       title  = "Phylum")

svglite(filename="phylum2votus.svg", width=9, height=13, fix_text_size=FALSE)
plot(network,
     vertex.label = NA)
legend("topleft",
       legend = names(phylum_colors),
       col    = phylum_colors,
       pch    = 19,
       pt.cex = 1.5,
       bty    = "n",
       title  = "Phylum")
dev.off()

##########################
#Network - genus - entire#
##########################

#data.frame for microbial phyla
nt_genus <- taxonomy_host_prediction %>%
  select(contig, phylum_v, phylum, genus, method_code) %>%
  filter(!is.na(genus))
dim(nt_genus)
#[1] 374   5

#remove phyla where only one/two vOTU because too chaotic
trim <- as.vector(host_prediction_count_phylum[host_prediction_count_phylum$Freq >= 3,] %>% select(Var1))
nt_genus_trim <- nt_genus %>% filter(phylum %in% trim$Var1)

#unique links
nt <- nt_genus_trim[,c("contig", "genus","method_code")]
dim(nt)
#[1] 362   3
nt_unique <- nt %>%
  group_by(contig, genus) %>%
  filter(n() == 1) %>%
  ungroup()
dim(nt_unique)
#[1] 352   3

#create network
network <- graph_from_data_frame(nt_unique, directed = FALSE)
str(network)

# names of nodes
node_names <- V(network)$name

# map viral phylum using contig names ----
mapped_phylum_v <- nt_genus$phylum_v[ match(node_names, nt_genus$contig) ]

# map microbial phylum using genus names ----
mapped_phylum_microbe <- nt_genus$phylum[ match(node_names, nt_genus$genus) ]

# combine: use viral phylum if contig node, else microbial phylum ----
mapped_phylum <- ifelse(
  !is.na(mapped_phylum_v),
  mapped_phylum_v,
  mapped_phylum_microbe
)

#add microbial and viral taxa metadata to vertices
V(network)$phylum <- mapped_phylum
unique(V(network)$phylum)

#colors for taxa with 5 or more viruses
phylum_colors <- c(
  "Actinomycetota" = "orange", #in mags
  "Pseudomonadota" = "mistyrose", #in mags
  "Bacteroidota" = "blue", #in mags
  "Cyanobacteriota" = "#1b9e77", #in mags
  "Bdellovibrionota" = "aquamarine1", #in mags
  
  "Acidobacteriota" = "#cab2d6", #in mags
  "Patescibacteriota" = "yellow", #in mags
  "Chloroflexota" = "#b2df8a", #in mags
  "Bacillota" = "lightcyan", 
  "Methanobacteriota" = "darkred",
  
  "Asgardarchaeota" = "#e7298a",
  "Campylobacterota" = "grey",
  "Spirochaetota" = "purple",
  
  "p__Uroviricota" = "brown1",
  "p__Hofneiviricota" = "lightpink2",
  "p__Nucleocytoviricota" = "#8da0cb",
  "p__Artverviricota" = "lightskyblue1",
  
  "p__Unclassified" = "white"
)

#shapes
phylum_shapes <- c(
  # microbial phyla → now **square**
  "Actinomycetota" = "square",
  "Pseudomonadota" = "square",
  "Bacteroidota" = "square",
  "Cyanobacteriota" = "square",
  "Planctomycetota" = "square",
  "Bacillota" = "square",
  "Methanobacteriota" = "square",
  "Bdellovibrionota" = "square",
  "Acidobacteriota" = "square",
  "Patescibacteriota" = "square",
  "Gemmatimonadota" = "square",
  "Armatimonadota" = "square",
  "Vulcanimicrobiota" = "square",
  "Chloroflexota" = "square",
  "Asgardarchaeota" = "square",
  "Caldisericota" = "square",
  "Thermoproteota" = "square",
  "Campylobacterota" = "square",
  "Spirochaetota" = "square",
  "Chlamydiota" = "square",
  "Deinococcota" = "square",
  "Halobacteriota" = "square",
  
  
  "Methanobacteriota" = "square",
  "Thermoplasmatota" = "square",
  "Verrucomicrobiota" = "square",
  "Myxococcota" = "square",
  "Hadarchaeota" = "square",
  
  # viral phyla → now **circle**
  "p__Uroviricota"        = "circle",
  "p__Hofneiviricota"     = "circle",
  "p__Nucleocytoviricota" = "circle",
  "p__Artverviricota"     = "circle",
  "p__Unclassified"       = "circle"
)

#shapes
s_m <- 5
s_v <- 4
phylum_sizes <- c(
  # microbial phyla → now **3**
  "Actinomycetota" = s_m,
  "Pseudomonadota" = s_m,
  "Bacteroidota" = s_m,
  "Cyanobacteriota" = s_m,
  "Planctomycetota" = s_m,
  "Bacillota" = s_m,
  "Methanobacteriota" = s_m,
  "Bdellovibrionota" = s_m,
  "Acidobacteriota" = s_m,
  "Patescibacteriota" = s_m,
  "Gemmatimonadota" = s_m,
  "Armatimonadota" = s_m,
  "Vulcanimicrobiota" = s_m,
  "Chloroflexota" = s_m,
  "Asgardarchaeota" = s_m,
  "Caldisericota" = s_m,
  "Thermoproteota" = s_m,
  "Campylobacterota" = s_m,
  "Spirochaetota" = s_m,
  "Chlamydiota" = s_m,
  "Deinococcota" = s_m,
  "Halobacteriota" = s_m,
  
  "Methanobacteriota" = s_m,
  "Thermoplasmatota" = s_m,
  "Verrucomicrobiota" = s_m,
  "Myxococcota" = s_m,
  "Hadarchaeota" = s_m,
  
  # viral phyla → now **20**
  "p__Uroviricota"        = s_v,
  "p__Hofneiviricota"     = s_v,
  "p__Nucleocytoviricota" = s_v,
  "p__Artverviricota"     = s_v,
  "p__Unclassified"       = s_v
)

V(network)$color <- phylum_colors[V(network)$phylum]

# assign shapes 
V(network)$shape <- phylum_shapes[V(network)$phylum]

# assign sizes
V(network)$size <- phylum_sizes[V(network)$phylum] * 0.7

method_colors <- c(
  "1" = "black",    # iPHoP-RF or multiple
  "2" = "#E69F00",  # CRISPR
  "3" = "deeppink",   # blast
  "4" = "blue"  # RaFAH
)

# Get edge list
edges <- as_data_frame(network, what="edges")

# Map method_code using contig (from column 'contig' in nt_genus_trim)
edges$method_code <- nt_unique$method_code[ match(edges$from, nt_unique$contig) ]

# Assign edge colors based on method_code
E(network)$color <- method_colors[as.character(edges$method_code)]

E(network)$width <- 1.7

# plot the network
plot(network,
     vertex.label = NA)
#I have the one above
# legend("topleft",
#        legend = names(phylum_colors),
#        col    = phylum_colors,
#        pch    = 19,
#        pt.cex = 1.5,
#        bty    = "n",
#        title  = "Phylum")

svglite(filename="genus2votus.svg", width=12, height=16, fix_text_size=FALSE)
plot(network,
     vertex.label = NA)
dev.off()

svglite(filename="genus2votus_names.svg", width=12, height=16, fix_text_size=FALSE)
plot(network)
dev.off()

#genera to display
#Var1 Freq
#1   Acidobacteriota   14
#2    Actinomycetota   35
#5         Bacillota   43
#6      Bacteroidota   73
#7  Bdellovibrionota   11
#12  Cyanobacteriota   40
#19   Pseudomonadota  130

###################
#vOTUs exploration#
###################

#Eukaryotic vOTUs
cont_euk1 <- nt_phylum %>% filter(phylum_v == "p__Nucleocytoviricota") %>% select(contig)
cont_euk1
cont_euk2 <- nt_phylum %>% filter(phylum_v == "p__Artverviricota") %>% select(contig)
cont_euk2
#saved to eukaryotic_votus.txt

#Bacillota linked vOTUs
cont_bacillota <- nt_phylum %>% filter(phylum == "Bacillota") %>% select(contig)
cont_bacillota