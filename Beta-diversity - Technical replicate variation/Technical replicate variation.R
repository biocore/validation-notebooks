#########################################################################################
# 2022
# Knight Lab - Validation notebooks
# Beta-diversity - Technical replicate distance analysis
# Justin Shaffer
# justinparkshaffer@gmail.com
#########################################################################################


# Set working directory
#########################################################################################
getwd()
setwd("/path/to/inputs/")


# Install and load libraries needed for analysis
#########################################################################################

install.packages("tidyverse")
install.packages("cowplot")
install.packages("ggpubr")
install.packages("svglite")

library(tidyverse)
library(cowplot)
library(ggpubr)
library(svglite)


# Beta-diversity - Technical replicate analysis
#########################################################################################

# Read in data
tech_rep_jaccard <- read_tsv("data_tech_reps_jaccard.txt")
tech_rep_rpca <- read_tsv("data_tech_reps_rpca.txt")
tech_rep_unifrac <- read_tsv("data_tech_reps_unifrac.txt")
tech_rep_wunifrac <- read_tsv("data_tech_reps_weighted_unifrac.txt")


# Re-order levels of the grouping variable (optional)
## Note: Change the levels to those for your variable
tech_rep_jaccard$sample1_extraction_kit_round <- factor(tech_rep_jaccard$sample1_extraction_kit_round, levels = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX Microbiome","NucleoMag Food", "Zymo MagBead"))
tech_rep_rpca$sample1_extraction_kit_round <- factor(tech_rep_rpca$sample1_extraction_kit_round, levels = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX Microbiome","NucleoMag Food", "Zymo MagBead"))
tech_rep_unifrac$sample1_extraction_kit_round <- factor(tech_rep_unifrac$sample1_extraction_kit_round, levels = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX Microbiome","NucleoMag Food", "Zymo MagBead"))
tech_rep_wunifrac$sample1_extraction_kit_round <- factor(tech_rep_wunifrac$sample1_extraction_kit_round, levels = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX Microbiome","NucleoMag Food", "Zymo MagBead"))


# Re-name levels of the grouping variable (optional)
## Note: Change the levels to those for your variable
levels(tech_rep_jaccard$sample1_extraction_kit_round) = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX","NucleoMag", "Zymo")
levels(tech_rep_rpca$sample1_extraction_kit_round) = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX","NucleoMag", "Zymo")
levels(tech_rep_unifrac$sample1_extraction_kit_round) = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX","NucleoMag", "Zymo")
levels(tech_rep_wunifrac$sample1_extraction_kit_round) = c("PowerSoil r1", "PowerSoil Pro", "Norgen", "PowerSoil r2", "MagMAX","NucleoMag", "Zymo")


# Re-order levels of the sample type variable (optional)
## Note: Change the levels to those for your variable
tech_rep_jaccard$sample_type_2 <- factor(tech_rep_jaccard$sample_type_2, levels = c("keyboard", "floor", "kimchi", "yogurt", "seawater", "freshwater", "bare soil", "mouse jejunum tissue", "mouse feces", "cat feces", "human feces", "human female urine", "human male urine", "human foot", "human armpit", "human saliva", "human milk"))
tech_rep_rpca$sample_type_2 <- factor(tech_rep_rpca$sample_type_2, levels = c("keyboard", "floor", "kimchi", "yogurt", "seawater", "freshwater", "bare soil", "mouse jejunum tissue", "mouse feces", "cat feces", "human feces", "human female urine", "human male urine", "human foot", "human armpit", "human saliva", "human milk"))
tech_rep_unifrac$sample_type_2 <- factor(tech_rep_unifrac$sample_type_2, levels = c("keyboard", "floor", "kimchi", "yogurt", "seawater", "freshwater", "bare soil", "mouse jejunum tissue", "mouse feces", "cat feces", "human feces", "human female urine", "human male urine", "human foot", "human armpit", "human saliva", "human milk"))
tech_rep_wunifrac$sample_type_2 <- factor(tech_rep_wunifrac$sample_type_2, levels = c("keyboard", "floor", "kimchi", "yogurt", "seawater", "freshwater", "bare soil", "mouse jejunum tissue", "mouse feces", "cat feces", "human feces", "human female urine", "human male urine", "human foot", "human armpit", "human saliva", "human milk"))


# Re-name levels for the sample type variable (optional)
## Note: Change the levels to those for your variable
levels(tech_rep_jaccard$sample_type_2) <- c("Surface, keyboard", "Surface, floor tile", "Food, kimchi", "Food, yogurt", "Water, saline", "Water, non-saline", "Soil, bare", "Mouse tissue, jejunum ", "Mouse feces", "Cat feces", "Human feces", "Human urine, female", "Human urine, male", "Human skin, foot", "Human skin, armpit", "Human saliva", "Human milk")
levels(tech_rep_rpca$sample_type_2) <- c("Surface, keyboard", "Surface, floor tile", "Food, kimchi", "Food, yogurt", "Water, saline", "Water, non-saline", "Soil, bare", "Mouse tissue, jejunum ", "Mouse feces", "Cat feces", "Human feces", "Human urine, female", "Human urine, male", "Human skin, foot", "Human skin, armpit", "Human saliva", "Human milk")
levels(tech_rep_unifrac$sample_type_2) <- c("Surface, keyboard", "Surface, floor tile", "Food, kimchi", "Food, yogurt", "Water, saline", "Water, non-saline", "Soil, bare", "Mouse tissue, jejunum ", "Mouse feces", "Cat feces", "Human feces", "Human urine, female", "Human urine, male", "Human skin, foot", "Human skin, armpit", "Human saliva", "Human milk")
levels(tech_rep_wunifrac$sample_type_2) <- c("Surface, keyboard", "Surface, floor tile", "Food, kimchi", "Food, yogurt", "Water, saline", "Water, non-saline", "Soil, bare", "Mouse tissue, jejunum ", "Mouse feces", "Cat feces", "Human feces", "Human urine, female", "Human urine, male", "Human skin, foot", "Human skin, armpit", "Human saliva", "Human milk")


# Create a plot for each distance metric
boxplot_jaccard <- ggplot(tech_rep_jaccard, aes(x = sample1_extraction_kit_round, y = value)) +
  facet_wrap(~sample_type_2) +
  geom_boxplot(outlier.alpha = 0, aes(fill = sample1_extraction_kit_round)) +
  geom_jitter(width = 0.2, size = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  ylab("Within sample\nJaccard distance") +
  geom_vline(xintercept = 3.5,
             color = "black",
             size = 0.25) +
  ylim(c(-0.01,5.5)) +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.background =element_rect(fill="white"))

boxplot_rpca <- ggplot(tech_rep_rpca, aes(x = sample1_extraction_kit_round, y = value)) +
  facet_wrap(~sample_type_2,
             ncol = 5) +
  geom_boxplot(outlier.alpha = 0, aes(fill = sample1_extraction_kit_round)) +
  geom_jitter(width = 0.2, size = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  ylab("Within sample\nRobust Aitchison distance") +
  geom_vline(xintercept = 3.5,
             color = "black",
             size = 0.25) +
  ylim(c(-0.01,5.5)) +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.background =element_rect(fill="white"))

boxplot_unifrac <- ggplot(tech_rep_unifrac, aes(x = sample1_extraction_kit_round, y = value)) +
  facet_wrap(~sample_type_2) +
  geom_boxplot(outlier.alpha = 0, aes(fill = sample1_extraction_kit_round)) +
  geom_jitter(width = 0.2, size = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  ylab("Within sample\nunweighted UniFrac distance") +
  geom_vline(xintercept = 3.5,
             color = "black",
             size = 0.25) +
  ylim(c(-0.01,5.5)) +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.background =element_rect(fill="white"))

boxplot_wunifrac <- ggplot(tech_rep_wunifrac, aes(x = sample1_extraction_kit_round, y = value)) +
  facet_wrap(~sample_type_2) +
  geom_boxplot(outlier.alpha = 0, aes(fill = sample1_extraction_kit_round)) +
  geom_jitter(width = 0.2, size = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  ylab("Within sample\nweighted UniFrac distance") +
  geom_vline(xintercept = 3.5,
             color = "black",
             size = 0.25) +
  ylim(c(-0.01,5.5)) +
  scale_fill_manual(values = c("blue", "steelblue1", "chocolate4", "blue", "green4", "white", "yellow")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        text = element_text(size = 14),
        strip.background =element_rect(fill="white"))


# Create multi-paneled plot showing all distance metrics
#########################################################################################

# Create plot
figure_tech_rep <- plot_grid(boxplot_jaccard,
                                 boxplot_rpca,
                                 boxplot_unifrac,
                                 boxplot_wunifrac,
                                 labels = c("A", "B", "C", "D"),
                                 label_size = 24,
                                 label_fontfamily = 'sans',
                                 ncol = 2,
                                 hjust = -0.25)


# View plot
figure_tech_rep


# Save plot as TIFF
save_plot('figure_tech_rep_20x15.tiff',
          figure_tech_rep,
          base_width = 20,
          base_height = 15)


# Save plot as SVG
save_plot('figure_tech_rep_20x15.svg',
          figure_tech_rep_16S,
          base_width = 20,
          base_height = 15)


