#################
###Trasposases###
#################

setwd("~/Desktop/R_icevirome/v3")

library(data.table)
library(dplyr)
library(svglite)

#N.B.
#No_ANNOTATION: proteins in flanks but no annotation in eggnog
#No_EGGNOG: in eggnog but no value
#No_PROTEIN: "-" in flank file

##############
#COG category#
##############

#import tables
cog_category = read.table("transposase_clusters_flanks_counts.cog_category.txt", row.names = 1, h = T, check.names = FALSE)
dim(cog_category)
#[1] 106  22
cog_category_trimmed <- cog_category[row.names(cog_category) != "No_PROTEIN", ]
dim(cog_category_trimmed)
#[1] 105  22

#calculate relative abundance
cog_category_rel <- sweep(cog_category_trimmed, 2, colSums(cog_category_trimmed), "/")

#Check which ones vary

# Initialize storage for the final table, keeping only rho and p_value
results <- data.frame(
  COG_category = character(),
  Spearman_rho = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Define distance vector
distance <- c(1,2,3,4,5,6,7,8,9,10,11,1,2,3,4,5,6,7,8,9,10,11)

# Loop through each row of the 'cog_category_rel' data frame
for (i in 1:nrow(cog_category_rel)) {
  # Extract values for the current COG category
  cat_values <- as.numeric(cog_category_rel[i, ])
  
  # Skip if the number of values doesn't match the distance vector length
  if(length(cat_values) != length(distance)) next
  
  # Spearman correlation for direction (rho)
  spearman_test <- cor.test(cat_values, distance, method = "spearman", exact = FALSE)
  rho <- spearman_test$estimate
  
  # Linear model to get p-value
  model <- lm(cat_values ~ distance)
  model_sum <- summary(model)
  
  # Calculate p-value from F-statistic
  # F-stat is model_sum$fstatistic[1], DF1 is model_sum$fstatistic[2], DF2 is model_sum$fstatistic[3]
  p_val <- pf(model_sum$fstatistic[1], model_sum$fstatistic[2], model_sum$fstatistic[3], lower.tail = FALSE) 
  
  # Only include significant lm results (p-value < 0.05)
  if (!is.na(p_val) && p_val < 0.05) {
    # Bind only the required columns to the results data frame
    results <- rbind(results, data.frame(
      COG_category = rownames(cog_category_rel)[i],
      Spearman_rho = rho,
      p_value = p_val,
      stringsAsFactors = FALSE
    ))
  }
}

# Print final table (contains only significant results)
print(results)

results[results$Spearman_rho > 0.5,]
# COG_category Spearman_rho      p_value
# rho              C    0.8949350 2.539618e-06
# rho6             E    0.5256327 4.278714e-03
# rho8            EQ    0.6704963 2.037710e-03
# rho9            EU    0.6220470 2.145235e-04
# rho10            F    0.5777428 6.852805e-03
# rho11            G    0.5165701 7.056686e-03
# rho13           GM    0.6003994 1.479805e-03
# rho16            H    0.8088399 1.556751e-05
# rho18           HP    0.7101911 3.325335e-04
# rho19            I    0.6479782 4.481527e-04
# rho21            J    0.5641489 1.580377e-03
# rho23           JM    0.6865150 2.077550e-03
# rho30           MU    0.6139934 5.456384e-03
# rho31            N    0.8224339 1.426239e-06
# rho33           NU    0.8246996 5.916646e-07
# rho36           OP    0.5353202 1.021899e-02
# rho37            P    0.5346953 2.444154e-03
# rho38            S    0.5414923 5.450248e-03

results[results$Spearman_rho < -0.5,]
# COG_category Spearman_rho      p_value
# rho1            CG   -0.5638703 0.0277647620
# rho5            DT   -0.5098911 0.0425761018
# rho20           IU   -0.5151194 0.0241263162
# rho25            L   -0.8926693 0.0038379801
# rho27           LO   -0.6877099 0.0004865645
# rho28           LU   -0.5230326 0.0268120121
# rho32         NPTU   -0.7096265 0.0046558630

#####
#COG#
#####

cog <- fread("transposase_clusters_flanks_counts.cog.txt",h = T)
cog <- as.data.frame(cog)
rownames(cog) <- cog[[1]]    # first column
cog[[1]] <- NULL
dim(cog)
#[1] 13709    24

cog_trimmed <- cog[row.names(cog) != "No_PROTEIN", ]
dim(cog_trimmed)
#[1] 13708    24

#calculate relative abundance
cog_rel <- sweep(cog_trimmed[,-c(1,2)], 2, colSums(cog_trimmed[,-c(1,2)]), "/")

#descriptions
description <- cbind(COG = row.names(cog), description=cog$Description, category=cog$COG_CAT)

# Initialize storage for the final table, keeping only rho and p_value
results <- data.frame(
  COG = character(),
  Spearman_rho = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Define distance vector
distance <- c(1,2,3,4,5,6,7,8,9,10,11,1,2,3,4,5,6,7,8,9,10,11)

# Loop through each row of the 'cog_category_rel' data frame
for (i in 1:nrow(cog_rel)) {
  # Extract values for the current COG category
  cat_values <- as.numeric(cog_rel[i, ])
  
  # Skip if the number of values doesn't match the distance vector length
  if(length(cat_values) != length(distance)) next
  
  # Spearman correlation for direction (rho)
  spearman_test <- cor.test(cat_values, distance, method = "spearman", exact = FALSE)
  rho <- spearman_test$estimate
  
  # Linear model to get p-value
  model <- lm(cat_values ~ distance)
  model_sum <- summary(model)
  
  # Calculate p-value from F-statistic
  # F-stat is model_sum$fstatistic[1], DF1 is model_sum$fstatistic[2], DF2 is model_sum$fstatistic[3]
  p_val <- pf(model_sum$fstatistic[1], model_sum$fstatistic[2], model_sum$fstatistic[3], lower.tail = FALSE) 
  
  # Only include significant lm results (p-value < 0.05)
  if (!is.na(p_val) && p_val < 0.05) {
    # Bind only the required columns to the results data frame
    results <- rbind(results, data.frame(
      COG = rownames(cog_rel)[i],
      Spearman_rho = rho,
      p_value = p_val,
      stringsAsFactors = FALSE
    ))
  }
}

# Print final table (contains only significant results)
print(results)

#only rho < -0.5
results_neg <- results[results$Spearman_rho < -0.5,]
dim(results_neg)
#[1] 447   3

#add description to cog
results_neg_description <- merge(results_neg, description, by = "COG", all.x = TRUE)
results_neg_description 

#checck how many cogs to each categors
data.frame(table(results_neg_description$category))

# Print final table (contains only significant results)
write.table(results_neg_description,"results_neg_description.txt")

#######################
##COG category - plot##
#######################

total_cog_cat <- read.table("total_COG_categories.txt", h = T, row.names = 1)
total_cog <- read.table("total_COGs.txt", h = T, row.names = 1)

total_number <- 325440

#create relative abundances  
total_cog_cat$rel_count <- (total_cog_cat$Count / total_number)

#reorder data
cog_category_rel_ord <- cog_category_rel[,c(11,10,9,8,7,6,5,4,3,2,1,
                                            12,13,14,15,16,17,18,19,20,21,22)]
#X-axis labels
x_labels <- colnames(cog_category_rel_ord)

# Vector of all categories to plot (ordered by type)
categories <- c("J","K","L", 
                "D","M", "N", "O", "T", "U", "V", 
                "C","E","F", "G", "H","I","P","Q", 
                "S")

#not present: "R", "Y", "Z", "W" 
#not abundant in general: "A", "B"

#format x axis
colnames(cog_category_rel_ord) <- c("-10-11", "-9-10", "-8-9", "-7-8", "-6-7", "-5-6", "-4-5", "-3-4", "-2-3", "-1-2", "-0-1", 
                                    "+0-1", "+1-2", "+2-3", "+3-4", "+4-5", "+5-6", "+6-7", "+7-8", "+8-9", "+9-10", "+10-11")


# Set up plotting area: one column, number of rows = number of categories
svglite(filename="transposon_clusters_eggnog.svg", width=8, height=11, fix_text_size=FALSE)
par(mfrow = c(7, 3), mar = c(2, 4, 2, 1)) # bottom, left, top, right margins
# Loop over each category
for(cat in categories) {
  plot(
    as.numeric(cog_category_rel_ord[cat, ]),
    type = "b",
    xaxt = "n",
    xlab = "",
    ylab = "Relative abundance",
    main = cat,
    ylim = c(0, 0.2)
  )
  axis(1, at = 1:length(colnames(cog_category_rel_ord)), labels = colnames(cog_category_rel_ord), las = 2, cex.axis = 0.8)
  abline(h = total_cog_cat[cat, 2], col = "red", lty = 2, lwd = 2)
}
dev.off()

#A	RNA processing and modification	Information Storage and Processing		
#B	Chromatin structure and dynamics	Information Storage and Processing		
#J	Translation, ribosomal structure and biogenesis	Information Storage and Processing		
#K	Transcription	Information Storage and Processing		
#L	Replication, recombination and repair	Information Storage and Processing		

#D	Cell cycle control, cell division, chromosome partitioning	Cellular Processes and Signaling		
#M	Cell wall/membrane/envelope biogenesis	Cellular Processes and Signaling		
#N	Cell motility	Cellular Processes and Signaling		
#O	Post-translational modification, protein turnover, chaperones	Cellular Processes and Signaling		
#T	Signal transduction mechanisms	Cellular Processes and Signaling		
#U	Intracellular trafficking, secretion, and vesicular transport	Cellular Processes and Signaling		
#V	Defense mechanisms	Cellular Processes and Signaling		
#W	Extracellular structures	Cellular Processes and Signaling		
#Y	Nuclear structure	Cellular Processes and Signaling		
#Z	Cytoskeleton	Cellular Processes and Signaling		

#C	Energy production and conversion	Metabolism		
#E	Amino acid transport and metabolism	Metabolism		
#F	Nucleotide transport and metabolism	Metabolism		
#G	Carbohydrate transport and metabolism	Metabolism		
#H	Coenzyme transport and metabolism	Metabolism		
#I	Lipid transport and metabolism	Metabolism		
#P	Inorganic ion transport and metabolism	Metabolism		
#Q	Secondary metabolites biosynthesis, transport and catabolism	Metabolism		

#R	General function prediction only	Poorly Characterized		
#S	Function unknown	Poorly Characterized		

##############
##COG - plot##
##############

# COGs you want to highlight
highlight_list <- list(
  M = c("4NJR9", "COG0438", "1HRW6",
        "1G1P3", 
        "47KJK"),
  O = c("47PAY", 
        "47RNU",
        "1INQ7"),
  T = c("1INN4",
        "1J0UC", "1GPYK", "1NU8D", "47XEZ",
        "47XFN", "47MJ1"),
  V = c("2JP11", 
        "1KNTE", 
        "1IPQ1", "1IQYT"),
  C = c("2WJWS", 
        "47JTF",
        "1INQK"),
  E = c("4ABRN", 
        "47JF7", "1IRQU",
        "4BDGG"),
  G = c("47NGW", "47Q36", "2FNA4", "2KER6", "47KZ5", "47U0P", "1WJPP", "1X4BZ", "47NZE",
        "47JBJ", "47Y9D",
        "47KAR"),
  P = c("1IQNY", "1IPNC",
        "1IPAN", "1IP00", "1IWGM", "47JGS", "2KKQN", "1IPBX", "47JTW",
        "4C44Z", "2K10Q")
)

# Corresponding colors for highlighted COGs
highlight_colors_list <- list(
  M = c("red", "red", "red", "blue", "yellow"),
  O = c("orange", "purple", "green"),
  T = c("brown",
        "pink", "pink", "pink", "pink",
        "darkgreen", "darkgreen"),
  V = c("cyan", "magenta", "darkblue", "darkblue"),
  C = c("mediumblue", "seagreen2", "palevioletred1"),
  E = c("darkred", 
        "lightgoldenrod1","lightgoldenrod1",
        "skyblue"),
  G = c("salmon", "salmon", "salmon", "salmon", "salmon", "salmon", "salmon", "salmon", "salmon",
        "aquamarine1","aquamarine1", "mediumorchid4"),
  P = c("firebrick", "firebrick",
        "pink1", "pink1", "pink1", "pink1", "pink1", "pink1", "pink1",
        "steelblue", "steelblue")
)

# Reorder data
cog_rel_ord <- cog_rel[, c(11,10,9,8,7,6,5,4,3,2,1,12,13,14,15,16,17,18,19,20,21,22)]
x_labels <- c("-10-11", "-9-10", "-8-9", "-7-8", "-6-7", "-5-6", "-4-5", "-3-4", "-2-3", "-1-2", "-0-1", 
              "+0-1", "+1-2", "+2-3", "+3-4", "+4-5", "+5-6", "+6-7", "+7-8", "+8-9", "+9-10", "+10-11")
ymax <- 0.0025

# Categories to loop through
categories <- c("M", "O", "T", "V", "C","E","G","P")

# Multi-panel layout
svglite(filename="transposon_clusters_cogs.svg", width=8, height=11, fix_text_size=FALSE)
par(mfrow = c(4, 2), mar = c(4, 4, 3, 1))

for(cat in categories) {
  
  # Extract COGs in this category
  cogs <- results_neg_description %>%
    filter(category == cat) %>% select(COG)
  
  if(nrow(cogs) == 0) next
  
  # Highlighted COGs and their colors
  highlight_cogs <- highlight_list[[cat]]
  highlight_colors <- highlight_colors_list[[cat]]
  
  # Empty plot
  plot(
    1:ncol(cog_rel_ord), type = "n",
    xaxt = "n", xlab = "", ylab = "RA",
    ylim = c(0, ymax), main = cat
  )
  axis(1, at = 1:length(x_labels), labels = x_labels, las = 2, cex.axis = 0.7)
  
  # Plot all COGs in black first
  for(i in 1:nrow(cogs)) {
    cog_name <- cogs$COG[i]
    if(!(cog_name %in% highlight_cogs)) {
      lines(as.numeric(cog_rel_ord[cog_name, ]), type = "l", col = "black", lwd = 1)
      #optional: add horizontal reference line for each COG
      abline(h = total_cog[cog_name, 2], col = "black", lty = 2)
    }
  }
  
  # Plot highlighted COGs in their specified colors
  for(i in seq_along(highlight_cogs)) {
    cog_name <- highlight_cogs[i]
    if(!(cog_name %in% cogs$COG)) next
    lines(as.numeric(cog_rel_ord[cog_name, ]), type = "l", col = highlight_colors[i], lwd = 3)
    #optional: add horizontal reference line for each COG
    abline(h = total_cog[cog_name, 2], col = highlight_colors[i], lty = 2)
  }
  
  # Legend for highlighted COGs
  if(length(highlight_cogs) > 0) {
    legend(
      "topright", legend = highlight_cogs,
      col = highlight_colors, lty = 1, lwd = 3, cex = 0.7
    )
  }
}
dev.off()


