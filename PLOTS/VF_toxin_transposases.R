#######################################
###Trasposases - virulence and toxin###
#######################################

setwd("~/Desktop/R_icevirome/v3")

library(svglite)
library(ggplot2)
library(gridExtra)
library(svglite)
library(tidyr)
library(data.table)

###########
#virulence#
###########

#import data
virulence <- read.table("virulence.txt", header = T)
dim(virulence)
#40 4
virulence$interval

#reorder
virulence_ord <- virulence[c(20,19,18,17,16,15,14,13,12,11,
                             10,9,8,7,6,5,4,3,2,1,21:40),]
virulence_ord$interval

virulence_ord$interval <- factor(virulence_ord$interval, levels = virulence_ord$interval)

p_virulence <- ggplot(virulence_ord, aes(x = interval, y = virulence_perc)) +
  geom_line(aes(group = 1)) +
  geom_point() +
  geom_hline(yintercept = 0.209418019, linetype = "dashed") +   # ← horizontal line
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(
    x = "Interval",
    y = "Virulence RA"
  )

#######
#toxin#
#######

#import data
toxin <- read.table("toxin.txt", header = T)
dim(toxin)
#40 4
toxin$interval

#reorder
toxin_ord <- toxin[c(20,19,18,17,16,15,14,13,12,11,
                             10,9,8,7,6,5,4,3,2,1,21:40),]
toxin_ord$interval

toxin_ord$interval <- factor(toxin_ord$interval, levels = toxin_ord$interval)

p_toxin <- ggplot(toxin_ord, aes(x = interval, y = toxin_perc)) +
  geom_line(aes(group = 1)) +
  geom_point() +
  geom_hline(yintercept = 0.020264872, linetype = "dashed") +   # ← horizontal line
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(
    x = "Interval",
    y = "Toxin RA"
  )

#save
svglite(filename="transposase_virulence_toxin.svg", width=6, height=7, fix_text_size=FALSE)
grid.arrange(p_virulence, p_toxin, ncol = 1)
dev.off()

################################################################################
# Analysis of Single Toxins: Abundance, Significance, and Enrichment
################################################################################

# 1. DATA IMPORT & NORMALIZATION (Assuming previous steps remain same)
tox_n_total <- fread("toxin_pathoFact_total_abundance_by_Description.tsv", h = TRUE)
tox_n <- fread("toxin_counts_by_Description.tsv", h = TRUE)
tox_n <- as.data.frame(tox_n)
rownames(tox_n) <- tox_n[[1]]
tox_n[[1]] <- NULL

# Match your 20 interval columns (11:30)
tox_n_trimmed <- tox_n[, c(11:30)] 

total <- fread("flank_interval_counts_total.tsv", h = TRUE)
total_trimmed <- total %>% filter(interval %in% colnames(tox_n_trimmed))

totals_vec <- setNames(total_trimmed$protein_count, total_trimmed$interval)
totals_vec <- totals_vec[colnames(tox_n_trimmed)]

tox_n_rel <- sweep(as.matrix(tox_n_trimmed), 2, totals_vec, FUN = "/")
tox_n_rel_trimmed <- as.data.table(tox_n_rel, keep.rownames = "NAME")

# 2. STATISTICAL TESTING
distance <- c(10,9,8,7,6,5,4,3,2,1,1,2,3,4,5,6,7,8,9,10)
results <- data.frame(NAME = character(), Rho = numeric(), P = numeric(), stringsAsFactors = FALSE)

for (i in 1:nrow(tox_n_rel_trimmed)) {
  cat_values <- as.numeric(tox_n_rel_trimmed[i, -1])
  if (length(cat_values) != length(distance) || all(is.na(cat_values))) next
  
  spearman_test <- cor.test(cat_values, distance, method = "spearman", exact = FALSE)
  rho <- as.numeric(spearman_test$estimate)
  
  model <- lm(cat_values ~ distance)
  model_sum <- summary(model)
  
  if (!is.null(model_sum$fstatistic)) {
    p_val <- pf(model_sum$fstatistic[1], model_sum$fstatistic[2], model_sum$fstatistic[3], lower.tail = FALSE)
    results <- rbind(results, data.frame(NAME = tox_n_rel_trimmed$NAME[i], Rho = rho, P = p_val))
  }
}

# explore significant

results[results$P < 0.05,]


# 3. FILTERING & JOINING
sig_names_corr <- results$NAME[results$P < 0.05]
sig_names_corr
sig_names_abundance <- tox_n_total$Description[1:17]
sig_names_all <- unique(c(sig_names_corr, sig_names_abundance))

plot_long <- tox_n_rel_trimmed %>% 
  filter(NAME %in% sig_names_all) %>%
  pivot_longer(cols = -NAME, names_to = "interval", values_to = "rel_abundance") %>%
  left_join(data.frame(interval = colnames(tox_n_rel_trimmed)[-1], distance = distance), by = "interval") %>%
  left_join(tox_n_total %>% select(Description, baseline_rel = rel_abundance) %>% rename(NAME = Description), by = "NAME") %>%
  mutate(direction = case_when(grepl("^-", interval) ~ "Negative", grepl("^\\+", interval) ~ "Positive", TRUE ~ "Other"))

# Apply enrichment filter (mean > baseline)
plot_long <- plot_long %>%
  group_by(NAME, baseline_rel) %>%
  mutate(mean_abundance = mean(rel_abundance, na.rm = TRUE)) %>%
  filter(mean_abundance > baseline_rel) %>%
  ungroup()

# 4. PREPARE TEXT LABELS FOR RHO & P
label_data <- plot_long %>%
  group_by(NAME) %>%
  summarize(
    y_pos = max(rel_abundance) * 0.95,
    baseline_rel = unique(baseline_rel)
  ) %>%
  left_join(results, by = "NAME") %>%
  mutate(label = paste0("rho = ", round(Rho, 2), "\np = ", format.pval(P, digits = 2)))

# 5. VISUALIZATION
line_data <- plot_long %>% filter(NAME %in% sig_names_corr)

p <- ggplot(plot_long, aes(x = distance, y = rel_abundance)) +
  # Baseline (Black dashed line) - Genome Average
  geom_hline(aes(yintercept = baseline_rel), 
             linetype = "dashed", color = "black", linewidth = 0.8, alpha = 1) +
  
  # NEW: Mean Abundance (Solid Green line) - Average across intervals
  geom_hline(aes(yintercept = mean_abundance), 
             linetype = "dashed", color = "red", linewidth = 0.8, alpha = 1) +
  
  # Raw data points (Red/Blue)
  geom_point(aes(color = direction), size = 2, alpha = 0.8) +
  
  # Conditional Trend Line (Significant only)
  geom_smooth(data = line_data, method = "lm", se = FALSE, color = "darkblue", linewidth = 1) +
  
  # Rho and P-value Text
  geom_text(data = label_data, aes(x = 10, y = y_pos, label = label), 
            hjust = 1, vjust = 1, size = 3, color = "black", inherit.aes = FALSE) +
  
  # Formatting
  scale_color_manual(values = c("Negative" = "red", "Positive" = "blue", "Other" = "gray")) +
  facet_wrap(~ NAME, scales = "free_y") +
  scale_x_continuous(breaks = seq(1, 10, 1)) +
  theme_bw() +
  labs(
    x = "Distance from Gene",
    y = "Relative Abundance",
    color = "Interval Sign"
  ) +
  theme(
    # General text colors
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    
    # FIX: Use 'panel.spacing.y' instead of 'panel.spacing.t'
    panel.spacing.y = unit(1.5, "lines"), 
    
    # Tall grey boxes: increase 't' and 'b' margins to expand the box height
    strip.text = element_text(
      color = "black", 
      margin = margin(t = 12, b = 12) #to increase margins in grey plot title boxes
    ),
    
    # Ensure the background and border of the box are defined
    strip.background = element_rect(fill = "grey90", color = "black")
  )

p

#save
svglite(filename="transposon_toxin_single.svg", width=10, height=11, fix_text_size=FALSE)
p
dev.off()
