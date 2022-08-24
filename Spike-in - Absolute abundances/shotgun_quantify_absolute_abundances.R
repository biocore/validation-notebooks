####################################################################################################
# Quantification of absolute abundances for shotgun metagenomics sequence data
# Adapted from Livia Zaramela
# Justin P. Shaffer
# justinparkshaffer@gmail.com
####################################################################################################


####################################################################################################
# Part 1: Generate linear models for absolute abundance calculation (synDNA-metagenomic)
####################################################################################################

# Set up working environment
install.packages("data.table") # For data wrangling
install.packages("reshape2") # For data wrangling
install.packages("ggplot2") # For making ggplots
install.packages("plyr") # For data wrangling
install.packages("quantreg") # For generating linear models
install.packages("ggpmisc") # For annotating ggplots

library(data.table)
library(reshape2)
library(ggplot2)
library(plyr)
library(quantreg)
library(ggpmisc)


# Set working directory
setwd("/path/to/inputs/")


# Import data
table_synthetic_hits <- read.delim("table_plasmid_sequence_hits.txt", h=T) # This is a feature-table with per-sample counts for each synthetic sequence. Feature IDs are plasmid IDs (e.g., 'p126')
metadata_pools <- read.delim("metadata_samples_plasmid_sequences.txt", h=T) # This is a sample metadata file that denotes which plasmid pool went into each sample (see the example file)
dilutions <- read.delim("metadata_dilutions_pool.tsv") # This is metadata file that describes how each plasmid sequence was diluted for each pool (see the example file)


# Parsing the files
## Calculate the total number of reads aligned to the plasmid sequences for each sample
table_synthetic_hit_totals <- data.frame(apply(table_synthetic_hits[,-1], 2, sum))
table_synthetic_hit_totals <- data.frame(sample_name_r = rownames(table_synthetic_hit_totals), table_synthetic_hit_totals = table_synthetic_hit_totals$apply.table_synthetic_hits....1...2..sum.)


## Convert feature-table to long-form and merge with total number of hits to plasmid sequences
table_synthetic_hits_long_with_totals_beta <- melt(table_synthetic_hits)
colnames(table_synthetic_hits_long_with_totals_beta) <- c("plasmid_id", "sample_name_r", "count")
table_synthetic_hits_long_with_totals <- merge(table_synthetic_hits_long_with_totals_beta, table_synthetic_hit_totals, by="sample_name_r")


# Merge updated feature-table with sample-pool information
table_synthetic_hits_long_with_totals_and_pools <- merge(table_synthetic_hits_long_with_totals, metadata_pools, by = "sample_name_r") 


# Aggregate counts from different sequencing lanes. If you have data from multiple sequencing lanes, uncomment the next line and comment out the line after the next.
#table_synthetic_hits_long_with_totals_and_pools_all_lanes <- aggregate(cbind(count, table_synthetic_hit_totals) ~ sample_name_r + plasmid_id + pool, table_synthetic_hits_long_with_totals_and_pools, sum)
table_synthetic_hits_long_with_totals_and_pools_all_lanes <- table_synthetic_hits_long_with_totals_and_pools


# Calculate counts per million 
table_synthetic_hits_long_with_totals_and_pools_all_lanes$counts_per_million <- (table_synthetic_hits_long_with_totals_and_pools_all_lanes$count/table_synthetic_hits_long_with_totals_and_pools_all_lanes$read_count_total)*1000000


# Calculate percentages
table_synthetic_hits_long_with_totals_and_pools_all_lanes$percentage <- (table_synthetic_hits_long_with_totals_and_pools_all_lanes$count/table_synthetic_hits_long_with_totals_and_pools_all_lanes$read_count_total)/100


# Calculate log10 of counts per million
table_synthetic_hits_long_with_totals_and_pools_all_lanes$counts_per_million_log10 <- log10(table_synthetic_hits_long_with_totals_and_pools_all_lanes$counts_per_million)


# Merge table with plasmid dilution information
table_synthetic_hits_long_with_totals_and_pools_all_lanes_with_dilutions <- merge(table_synthetic_hits_long_with_totals_and_pools_all_lanes, dilutions, by = c("plasmid_id", "pool"))


# Rename final table for tidyness
final_table_1 <- table_synthetic_hits_long_with_totals_and_pools_all_lanes_with_dilutions


# Adjust the dilution factor. The pooling uses a dilution of 1:20.
dilution_log <- final_table_1$dilution_id
dilution_log[which(dilution_log == 1)] <- 1/20
dilution_log[which(dilution_log == 2)] <- 0.1/20
dilution_log[which(dilution_log == 3)] <- 0.01/20
dilution_log[which(dilution_log == 4)] <- 0.001/20
dilution_log[which(dilution_log == 5)] <- 0.0001/20

final_table_1$dilution_log <- dilution_log


# Plot log10(counts per million) by the dilutions
ggplot(final_table_1, aes(counts_per_million_log10, log10((dilution_log)))) +
  geom_point(size = 1) +
  facet_wrap(~pool + sample_name, 
             ncol = 8) +
  geom_smooth(method = "lm", 
              formula = y ~ x, 
              se = FALSE, 
              size = 0.5) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = ..rr.label..), 
               parse = TRUE, 
               label.x.npc = "left", 
               size = 3) +
  stat_fit_glance(method = 'lm', 
                  geom = 'text', 
                  aes(label = paste0('p = ', format(..p.value.., 3))), 
                  label.x = 1.5, 
                  label.y = -1.8, 
                  size = 3) + 
  labs(y = expression(paste("dilutions ", 10^"n")), 
       x = "log10(counts per million)") +  
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_bw() + 
  guides(fill = guide_legend(ncol = 1))  +
  theme(panel.border = element_rect(fill = NA, 
                                    size = 1),
        strip.text.x = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_blank())


# Generate a linear model for each sample
v1 <- unique(final_table_1$pool)
v3 <- unique(final_table_1$sample_name_pool)
v1v2v3 <- expand.grid(v1,v3)


# Extract unique pool-sample combinations
combinations <- v1v2v3[!duplicated(v1v2v3),]
sel <- c(intersect(grep("pool1", combinations$Var1), grep("pool1", combinations$Var2)),
         intersect(grep("pool2", combinations$Var1), grep("pool2", combinations$Var2)))

combinations <- combinations[sel,]
colnames(combinations) <- c("pool", "sample_name_pool")


# Create a function for pulling variables from linear models for each sample
coefall <- NULL

for(i in 1:dim(combinations)[1]){
  criteria <- as.vector(unlist(combinations[i,]))
  temp <- subset(final_table_1, pool %in% criteria & sample_name_pool %in% criteria)
  temp <- temp[!(temp$count == 0), ]
  lmtemp <- lm(temp$dilution_id ~ temp$counts_per_million_log10)
  coef <- c(criteria, lmtemp$coefficients[1], lmtemp$coefficients[2])
  coefall <- rbind(coefall, coef)
}

colnames(coefall) <- c("pool", "sample_name_pool", "a_intercept", "b_slope")
rownames(coefall) <- NULL
coefall <- as.data.frame(coefall)


# Export table to file
write.table(coefall, "results_synDNA_linear_model_variables.txt", 
            quote = F, 
            row.names = F, 
            sep="\t")


####################################################################################################
# Part 2: Perform analysis
####################################################################################################

# Set up environment
library(tidyverse)
library(plyr)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(broom)

library(Hmisc)
library(corrplot)
library(ggpmisc)
library(scales)
library(FactoMineR)
library(factoextra)


# Set working directory
setwd("~/Mycelium/R/2022_absolute_abundances_synDNA-metagenomic/")


# Load synDNA metadata_pools, linear model variables, flow data, and feature-table
linear_models <- read.delim("results_synDNA_linear_model_variables.txt", h=T)
metadata_pools <- read.delim("rna_integrity_metaG_synDNA_metadata.txt", h=T)
table_community_beta <- fread("table_community_hits.txt", head=T)
table_community <- as.data.frame(table_community_beta)


# Clean sample names
linear_models$sample_name <- gsub("_[^_]+$", "", linear_models$sample_name_pool)


# Read in feature metadata including genome sizes per feature
metadata_features <- fread("feature_metadata.tsv", sep = "\t", h=T)
metadata_features <- as.data.frame(metadata_features)


# Calculate raw counts for each feature
feature_counts_raw <- NULL

for(i in 2:dim(table_community)[2]){
  df <- data.frame(OTUID = table_community$OTUID, count_raw = table_community[,i])
  SampleIDs <- colnames(table_community)
  
  metaspecies <- data.frame(OTUID = metadata_features$OTUID,
                            Species = metadata_features$species, 
                            Length = metadata_features$total_length)
  
  metaspecies <- merge(df, metaspecies, by = "OTUID")
  metaspecies <- metaspecies[which(metaspecies$count_raw != 0),]
  dfaggreg <- aggregate(metaspecies$count_raw, by=list(metaspecies$OTUID), FUN=sum)
  colnames(dfaggreg) <- c("OGU", "count_raw")
  dfout <- data.frame(SampleID = SampleIDs[i], AlignedReads = sum(dfaggreg$count_raw), dfaggreg)
  feature_counts_raw <- rbind(feature_counts_raw, dfout)
}


# Estimate genome length (average) and coverage
genome_length_average <- NULL
total_reads <- NULL
for(i in 1:dim(feature_counts_raw)[1]){
  genome_length_average <- c(genome_length_average, mean(metadata_features$total_length[which(metadata_features$OTUID == feature_counts_raw$OGU[i])]))
  total_reads <- c(total_reads, metadata_pools$read_count_total[which(metadata_pools$sample_name == feature_counts_raw$SampleID[i])])
}

feature_counts_raw$genome_length <- genome_length_average
feature_counts_raw$total_reads <- total_reads
feature_counts_raw$coverage <- (feature_counts_raw$count_raw*150)/feature_counts_raw$genome_length


# Remove lowly abundant species (coverage filter)
# NOTE: If you would like to use a coverage filter uncomment the next line and set your threshold there, and comment out the line after the next one.
#feature_counts_filtered <- feature_counts_raw[which(feature_counts_raw$coverage >= 1),]
feature_counts_filtered <- feature_counts_raw


# Apply linear models
new_and_old_counts <- NULL

for(i in 1:dim(linear_models)[1]){
  df <- feature_counts_filtered[feature_counts_filtered$SampleID %in% linear_models$sample_name[i],]
  df$count_reads_normalized_rowsums <- ((df$count_raw)/sum(df$count_raw)) * 100 
  df$count_reads_normalized_rowsums_log10 <- log10((df$count_raw/sum(df$count_raw)) * 1000000) # Because we used log10 values to build the model
  convert_reads <- linear_models$a_intercept[i] + (linear_models$b_slope[i]*df$count_reads_normalized_rowsums_log10)
  convert_reads <- 10^(-convert_reads) ## counts were in log scale.
  genome_size <- df$genome_length
  df$count_cells <- (convert_reads * (6.022 * (10^23))) / (genome_size * (1 * (10^9)) * 650)
  df$count_cells_normalized_rowsums <- (df$count_cells/sum(df$count_cells)) * 100
  
  df <- df[order(df$count_cells_normalized_rowsums, decreasing = TRUE), ]
  
  new_and_old_counts <- rbind(new_and_old_counts, df)
}

new_and_old_counts_with_metadata <- merge(metadata_pools, new_and_old_counts, by.x = "sample_name", by.y = "SampleID")


# Create feature-table with cell counts
feature_table_cell_counts <- pivot_wider(new_and_old_counts_with_metadata,
                                         id_cols = sample_name,
                                         names_from = OGU,
                                         values_from = count_cells,
                                         values_fill = 0)

write_tsv(feature_table_cell_counts,
          file = "feature_table_cell_counts.txt")


# Create feature-table with predicted relative abundances
feature_table_cell_counts_relative_abundance <- pivot_wider(new_and_old_counts_with_metadata,
                                                            id_cols = sample_name,
                                                            names_from = OGU,
                                                            values_from = count_cells_normalized_rowsums,
                                                            values_fill = 0)

write_tsv(feature_table_cell_counts_relative_abundance,
          file = "feature_table_cell_counts_normalized.txt")

## NOTE: If all you need are the estimates of cell counts, you can stop here.


# Select pool. Use whichever pool showed the best relationships in correlation analysis above
new_and_old_counts_with_metadata <- new_and_old_counts_with_metadata[grep("pool2", new_and_old_counts_with_metadata$pool), ]


# Correlation analysis: synthetic-based predicted vs. relative abundance
analysis <- new_and_old_counts_with_metadata %>%
  group_by(sample_name) %>%
  nest() %>%
  mutate(model = map(data, ~lm(count_cells_normalized_rowsums ~ count_reads_normalized_rowsums, data = .)),
         cor = map(data, ~tidy(cor.test(.x$count_cells_normalized_rowsums, .x$count_reads_normalized_rowsums, method = "pearson"), 3)))

stats <- analysis %>%
  unnest(cor)


# Plotting correlation graphs
ggplot(new_and_old_counts_with_metadata, aes(count_cells_normalized_rowsums, count_reads_normalized_rowsums)) +
  geom_point(shape = 21, 
             size = 1, 
             fill = "blue") +
  facet_wrap(~sample_name, 
             ncol = 8) +
  geom_smooth(method = "lm", 
              color = "black", 
              se = TRUE, 
              formula = y ~ x, 
              size = 0.5) +
  geom_text(data = stats, 
            aes(label = sprintf("r = %s", round(estimate, 3)), 
                x = 7, 
                y = 3), 
            size = 2.5, 
            hjust = 0) +
  geom_text(data = stats, 
            aes(label = sprintf("p = %s", formatC(p.value, format = "e", digits = 2)),  
                x = 7, 
                y = 1.5), 
            size = 2.5, 
            hjust = 0) +
  labs(x = "relative abundance (%)\nbased on cell count estimates", 
       y = "relative abundance (%)\nbased on read counts") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_bw() + 
  guides(fill = guide_legend(ncol  =1)) +
  theme(axis.title = element_text(size = 16))

