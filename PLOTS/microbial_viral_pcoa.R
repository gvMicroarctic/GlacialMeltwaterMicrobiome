############
####PcoA####
############

library(svglite)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(gridExtra)
library(Maaslin2)

setwd("~/Desktop/R_icevirome/v3")

##################
##microorganisms##
##################

#get count information
assembly_count <- read.table("assembly_abundance_10.txt", sep = "\t", header = T, row.names = 1)

#get only columns with tpm counts
assembly_tpm0 <- assembly_count[, grepl("TPM", colnames(assembly_count))]
#format sample names (columns)
colnames(assembly_tpm0) <- gsub("^a\\.|_si\\.TPM$", "", colnames(assembly_tpm0))
dim(assembly_tpm0)
#[1] 780851     23
assembly_tpm <- assembly_tpm0[,-(ncol(assembly_tpm0))]
dim(assembly_tpm)
#[1] 780851     22
assembly_tpm_name <- cbind(assembly_tpm, contig = row.names(assembly_tpm)) #I need it for below steps

#remove big data.frame
rm(assembly_count)

#add contig information to count_entire_tpm
tpm_contig <- cbind(assembly_tpm, contig = row.names(assembly_tpm_name))
head(tpm_contig)

######
#PCoA#
######

#transpose and transform data for Bray-Curtis
dt <- t(sqrt(assembly_tpm))  #square root transformation

#calculate Bray-Curtis distance
bc_dist <- vegdist(dt, method = "bray")

#run PCoA
pcoa_coords <- cmdscale(bc_dist, k = 2, eig = TRUE)

#extract coordinates and variance
pcoa_scores <- as.data.frame(pcoa_coords$points)
colnames(pcoa_scores) <- c("PCoA1", "PCoA2")
explained_var <- pcoa_coords$eig / sum(pcoa_coords$eig)

#samples
samples <- c("BM1a", "BM1b", "BM2", "BM3", "BM4", "BM5", "BR1", "BR2", "BR3", "BR4", 
             "BR5", "VM1a", "VM1b", "VM3a", "VM3b", "VM4", "VM5", "VR1", "VR2", 
             "VR3", "VR4", "VR5")

#create factors for plotting
glacier <- ifelse(grepl("BM|VM", samples), "M",
                  ifelse(grepl("BR|VR", samples), "R", "N"))
dataset <- ifelse(grepl("BM|BR", samples), "B",
                  ifelse(grepl("VM|VR", samples), "V", "N"))
pcoa_scores$Sample <- samples
pcoa_scores$Glacier <- factor(glacier)
pcoa_scores$Dataset <- factor(dataset)

#plot PCoA
p_pcoa_m <- ggplot(pcoa_scores, aes(x = PCoA1, y = PCoA2, fill = Dataset, shape = Glacier)) +
  geom_point(size = 5) +
  # geom_text(aes(label = Sample), vjust = -1, size = 3) +
  scale_shape_manual(values = c("M" = 21, "R" = 24)) +
  scale_fill_manual(values = c("B" = "darkgreen", "V" = "darkorange")) +
  xlab(paste0("PCoA1 (", round(explained_var[1] * 100, 1), "%)")) +
  ylab(paste0("PCoA2 (", round(explained_var[2] * 100, 1), "%)")) +
  theme_bw() +
  theme(legend.position = "right",
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        plot.title = element_text(color = "black"))
p_pcoa_m

###########
#permanova#
###########

metadata <- data.frame(vb = c("b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b",
                              "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v"),
                       glacier = c("m", "m", "m", "m", "m", "m",
                                   "r", "r", "r", "r", "r",
                                   "m", "m", "m", "m", "m", "m",
                                   "r", "r", "r", "r", "r"))

dist_matrix <- vegdist(t(sqrt(assembly_tpm)), method = "bray")  #Bray-Curtis distance
#clr and euclidean worst results.. so keeping Hellinger with BC; also BC so then constistent with PCOA

permanova_result <- adonis2(dist_matrix ~ vb * glacier, data = metadata, permutations = 999, by = "terms")
permanova_result
# Df SumOfSqs      R2       F Pr(>F)    
# vb          1   1.1165 0.25324 10.4933  0.001 ***
# glacier     1   1.1044 0.25050 10.3795  0.001 ***
# vb:glacier  1   0.2727 0.06186  2.5633  0.020 *  
# Residual   18   1.9152 0.43440                   
# Total      21   4.4089 1.00000 

#########
#viruses#
#########

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

######
#PCoA#
######

#transpose and transform data for Bray-Curtis
dt <- t(sqrt(votu_tpm))  #square root transformation

#calculate Bray-Curtis distance
bc_dist <- vegdist(dt, method = "bray")

#run PCoA
pcoa_coords <- cmdscale(bc_dist, k = 2, eig = TRUE)

#extract coordinates and variance
pcoa_scores <- as.data.frame(pcoa_coords$points)
colnames(pcoa_scores) <- c("PCoA1", "PCoA2")
explained_var <- pcoa_coords$eig / sum(pcoa_coords$eig)

#samples
samples <- c("BM1a", "BM1b", "BM2", "BM3", "BM4", "BM5", "BR1", "BR2", "BR3", "BR4", 
             "BR5", "VM1a", "VM1b", "VM3a", "VM3b", "VM4", "VM5", "VR1", "VR2", 
             "VR3", "VR4", "VR5")

#create factors for plotting
glacier <- ifelse(grepl("BM|VM", samples), "M",
                  ifelse(grepl("BR|VR", samples), "R", "N"))
dataset <- ifelse(grepl("BM|BR", samples), "B",
                  ifelse(grepl("VM|VR", samples), "V", "N"))
pcoa_scores$Sample <- samples
pcoa_scores$Glacier <- factor(glacier)
pcoa_scores$Dataset <- factor(dataset)

#plot PCoA
p_pcoa_v <- ggplot(pcoa_scores, aes(x = PCoA1, y = PCoA2, fill = Dataset, shape = Glacier)) +
  geom_point(size = 5) +
  # geom_text(aes(label = Sample), vjust = -1, size = 3) +
  scale_shape_manual(values = c("M" = 21, "R" = 24)) +
  scale_fill_manual(values = c("B" = "darkgreen", "V" = "darkorange")) +
  xlab(paste0("PCoA1 (", round(explained_var[1] * 100, 1), "%)")) +
  ylab(paste0("PCoA2 (", round(explained_var[2] * 100, 1), "%)")) +
  theme_bw() +
  theme(legend.position = "right",
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        plot.title = element_text(color = "black"))
p_pcoa_v

###########
#permanova#
###########

metadata <- data.frame(vb = c("b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b",
                              "v", "v", "v", "v", "v", "v", "v", "v", "v", "v", "v"),
                       glacier = c("m", "m", "m", "m", "m", "m",
                                   "r", "r", "r", "r", "r",
                                   "m", "m", "m", "m", "m", "m",
                                   "r", "r", "r", "r", "r"))

dist_matrix <- vegdist(t(sqrt(votu_tpm)), method = "bray")  #Bray-Curtis distance
#clr and euclidean worst results.. so keeping Hellinger with BC; also BC so then constistent with PCOA

permanova_result <- adonis2(dist_matrix ~ vb * glacier, data = metadata, permutations = 999, by = "terms")
permanova_result
# Df SumOfSqs      R2       F Pr(>F)    
# vb          1   1.4516 0.28561 12.9228  0.001 ***
# glacier     1   1.2577 0.24745 11.1960  0.001 ***
# vb:glacier  1   0.3513 0.06912  3.1272  0.012 *  
# Residual   18   2.0220 0.39783                   
# Total      21   5.0825 1.00000                   

######
#Plot#
######

svglite(filename="microbial_viral_pcoa.svg", width=13, height=5.5, fix_text_size=FALSE)
grid.arrange(p_pcoa_m, p_pcoa_v, nrow = 1)
dev.off()
