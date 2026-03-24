# 16S Metagenomic Sequence Analysis - Bee Gut Bacteria
# Project: Anuja's Lysate Caging Experiment
# Analyst: Yva
# Date: December 2025
# Description: Analysis of hoeny bee gut microbiome response to ASx lysate
#              tetracycline, and Terra-Pro treatments

# packages
# uncomment & run these lines only once:
#install.packages("tidyverse")
#install.packages("vegan")
#install.packages("openxlsx")


#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("phyloseq")
#BiocManager::install("microbiome")
#BiocManager::install("DESeq2")
#BiocManager::install("ANCOMBC")


library(phyloseq) # v1.54.0 - microbiome analysis
library(vegan)  # used for diversity calculations?
library(readxl) 
library(DESeq2) # v1.50.2 - differential abundance
library(ANCOMBC) # v2.12.0 - composition analysis
library(microbiome) # additional microbiome tools (part of ancombc)
library(openxlsx)
library(tidyverse) # DESeq2 library will mask dplyr from tidyverse
#                   don't mess with the library order!


# directories
project_dir <- "C:/Users/elin1473/Documents/r_projects/anuja_lacto"
data_dir <- file.path(project_dir, "data")
metadata_path <- file.path(project_dir, "data/sample_metadata.xlsx")
output_dir <- file.path(project_dir, "output")
figures_dir <- file.path(project_dir, "figures")
notes_dir <- file.path(project_dir, "notes")


# read taxonomy data files
csv_files <- list.files(data_dir,
                        pattern = "\\.csv$",
                        recursive = TRUE,
                        full.names = TRUE)

# exclude the excluded_samples folder
csv_files <- csv_files[!str_detect(csv_files, "excluded_samples")]

# this section is just for standardizing the sample id from anuja's metadata sheet
# with the sample id from plasmidsaurus - the sample look-up table is for showing the
# original names from both files and the changes implemented
# the metadata_clean.csv merges anuja's metadata with the new sample name which has:
# Day Treatment Replicate number (1-3)

# clean metadata starts here
metadata <- read_excel(metadata_path, skip = 3) # column headers are buried on row 4

# clean metadata and create standardized sample names
metadata_clean <- metadata %>%
  # extract day & cage/rep numbers from tube name
  mutate(Day = str_extract(`Tube Name`, "Day-\\d+") %>%
           str_remove("Day-"),
         Cage = str_extract(`Tube Name`, "(Cage|Rep)-(\\d+)") %>%
           str_extract("\\d+"),
         Treatment = Treatments) %>%
  # fix apiary-marked samples (tubes 39-41)
  mutate(Day = case_when(
    `SAMPLE Tube ID` %in% c("39", "40", "41") ~ "15",
    TRUE ~ Day),
    Treatment = case_when(
      `SAMPLE Tube ID` %in% c("39", "40", "41") ~ "Apiary-Marked",
      TRUE ~ Treatment
    )
  ) %>%
  
  group_by(Day, Treatment) %>%
  mutate(Rep = row_number()) %>%
  ungroup() %>% # leave it for the double mutate??
# create standardized sample names - check the spaces & adjust
  mutate(NewSampleName = paste("Day ", Day, " ",
                               str_replace_all(Treatment, " ", " "),
                               "Rep ", Rep))

# read & merge sequencing data
raw_data <- map_dfr(csv_files, function(file) {
  sample_name <- basename(dirname(file))
  tube_id <- str_extract(sample_name, "\\d+$")
  read_csv(file, show_col_types = FALSE) %>%
    mutate(OriginalSampleName = sample_name, TubeID = tube_id)
})

raw_data_named <- raw_data %>%
  left_join(metadata_clean, by = c("TubeID" = "SAMPLE Tube ID")) %>%
  mutate(NewSampleName = str_remove(NewSampleName, " rerun"))

# create lookup table & save
sample_name_lookup <- raw_data_named %>%
  select(OriginalSampleName, TubeID, NewSampleName, Day, Treatment, Rep) %>%
  distinct() %>%
  arrange(as.numeric(TubeID)) # is this correct??

write_csv(sample_name_lookup, file.path(output_dir, "sample_name_lookup.csv"))
write_csv(metadata_clean, file.path(output_dir, "metadata_clean.csv")) # save for supplemental xlxs
# clean metadata finished


# sequencing stats - each file has ont_stats.tsv which has details for the number of reads etc.
# use this to generate the read stats for the dataset

# find ont stats files
ont_stats_files <- list.files(data_dir, 
                              pattern = "_ont-stats.tsv",
                              recursive = TRUE, 
                              full.names = TRUE)
# exclude the excluded file
ont_stats_files <- ont_stats_files[!str_detect(ont_stats_files, "excluded_samples")]

ont_stats_all <- map_dfr(ont_stats_files, function(file) {
  sample_name <- basename(dirname(file))
  tube_id <- str_extract(sample_name, "\\d+$")
  
  read_tsv(file, col_names = c("metric", "value"),
           show_col_types = FALSE) %>%
    pivot_wider(names_from = metric, values_from = value) %>%
    mutate(OriginalSampleName = sample_name, TubeID = tube_id)
})

# merge with metadata for NewSampleName association etc.
ont_stats_named <- ont_stats_all %>%
  left_join(metadata_clean %>%
              select(`SAMPLE Tube ID`, NewSampleName, Day, Treatment, Rep),
            by = c("TubeID" = "SAMPLE Tube ID")) %>%
  mutate(NewSampleName = str_remove(NewSampleName, " rerun"))

write_csv(ont_stats_named, file.path(output_dir, 
                                     "sequencing_stats_all_samples.csv")) # supplemental xlxs

# summary statistics - all samples
sequencing_summary <- ont_stats_named %>%
  summarize(
    Total_Samples = n(),
    Mean_Reads = mean(total_num_reads, na.rm = TRUE),
    SD_Reads = sd(total_num_reads, na.rm = TRUE),
    Min_Reads = min(total_num_reads, na.rm = TRUE),
    Max_Reads = max(total_num_reads, na.rm = TRUE),
    
    Mean_Bases = mean(total_num_bases, na.rm = TRUE),
    
    Mean_Read_Length = mean(mean_read_length, na.rm = TRUE),
    
    Mean_Quality = mean(mean_quality, na.rm = TRUE),
    Mean_N50 = mean(read_n50, na.rm = TRUE)
  )

# summary by treatment - all samples
treatment_summary <- ont_stats_named %>%
  group_by(Treatment) %>%
  summarize(
    N_Samples = n(),
    Mean_Reads = mean(total_num_reads, na.rm = TRUE),
    SD_Reads = sd(total_num_reads, na.rm = TRUE),
    Mean_Quality = mean(mean_quality, na.rm = TRUE)
  ) %>%
  arrange(desc(Mean_Reads))

write_csv(treatment_summary, file.path(output_dir, "sequencing_stats_by_treatment.csv"))
# end of sequencing stats chonk

# let's start some analysis!
# phyloseq used to import, store, and analyze taxonomic & abundance data

# create phyloseq objects
otu_table_wide <- raw_data_named %>%
  select(`Taxonomy ID`, NewSampleName, `Relative Abundance`) %>%
  pivot_wider(names_from = NewSampleName,
              values_from = `Relative Abundance`,
              values_fill = 0) %>%
  column_to_rownames("Taxonomy ID")

# create taxonomy table
tax_table_df <- raw_data_named %>%
  select(`Taxonomy ID`, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  distinct()
# fix column consistency issue
names(tax_table_df) <- tolower(names(tax_table_df))

# taxonomy corrections - reclass lactos Firm5
# motta & moran 2025 paper has bombilactobacillus for mellis & mellifer

# lactobacillus are dominant in 40/44 samples
# however, changing taxa names for mellis & mellifer to bombilactobacillus
# changed the dominance to 38/44
# to show lactobacillus dominance in 40/44 samples, comment or skip 
# this section - or ignore and go with 38/44
tax_table_df <- tax_table_df %>%
  mutate(
    genus = case_when(
      species == "Lactobacillus mellis" ~ "Bombilactobacillus",
      species == "Lactobacillus mellifer" ~ "Bombilactobacillus",
      TRUE ~ genus
    ),
    
    species = case_when(
      species == "Lactobacillus mellis" ~ "Bombilactobacillus mellis",
      species == "Lactobacillus mellifer" ~ "Bombilactobacillus mellifer",
      TRUE ~ species
    )
  )

# convert to rownames
tax_table_df <- tax_table_df %>% column_to_rownames("taxonomy id")

# create sample metadata
sample_metadata <- raw_data_named %>%
  select(NewSampleName, Day, Treatment, Rep) %>%
  distinct() %>%
  mutate(Day = as.numeric(Day)) %>%
  column_to_rownames("NewSampleName")

# convert to matrices for phyloseq
otu_mat <- as.matrix(otu_table_wide)
tax_mat <- as.matrix(tax_table_df)

# create phyloseq object with relative abundance
ps <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(sample_metadata)
)
# check it - remove this when no longer needed
print(ps)

## ANALYSIS SECTION ##
# start with taxonomic composition plots (genus & species lvl)
# genus taxa bar graphs
ps_genus <- tax_glom(ps, taxrank = "genus", NArm = FALSE)
top_genera <- names(sort(taxa_sums(ps_genus), decreasing = TRUE) [1:20])

# transform to relative abundance & prepare for plots
ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x/sum(x))

ps_genus_melt <- psmelt(ps_genus_rel) %>%
  mutate(
    Day_numeric = as.numeric(Day),
    Day_for_sorting = if_else(is.na(Day_numeric), 999, Day_numeric),
    Genus_display = ifelse(OTU %in% top_genera, genus, "Other"),
    Treatment_display = case_when(
      Treatment == "For baseline data" ~ "NEWs",
      Treatment == "Control" ~ "Control",
      Treatment == "Tetracycline-treated" ~ "Tetracycline",
      Treatment == "Terra-Pro-treated" ~ "Terra-Pro",
      Treatment == "ASx-lysate-treated" ~ "ASx lysate",
      Treatment == "Apiary-Marked" ~ "Hive Reared",
      Treatment == "Microbiome in the comb samples" ~ "Comb",
      TRUE ~ Treatment
    ),
    Treatment_order = factor(Treatment_display,
                             levels = c("NEWs",
                                        "Comb",
                                        "Control",
                                        "Tetracycline",
                                        "Terra-Pro",
                                        "ASx lysate",
                                        "Hive Reared")),
    Sample_label = paste("Day", Day, "", Treatment_display, "Rep-", Rep)
  )

# fiiiigure tiiiime!
# genus bar plot
p_genus <- ps_genus_melt %>%
  ggplot(aes(x = reorder(Sample_label, Day_for_sorting),
             y = Abundance, 
             fill = Genus_display)) +
  geom_bar(stat= "identity", width = 0.8) +
  facet_grid(~ Treatment_order, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(
    axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(face = "italic"),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    x = "Sample",
    y = "Relative Abundance",
    title = "Taxonomic Composition by Genus",
    fill = "Genus"
  ) +
  scale_fill_viridis_d(option = "turbo")

ggsave(file.path(figures_dir, "taxonomy_by_genus.pdf"), 
       plot = p_genus, width = 18, height = 8)

ggsave(file.path(figures_dir, "taxonomy_by_genus.png"), 
       plot = p_genus, width = 18, height = 8)

# reevaluate for species
# top 20 species
top_species <- names(sort(taxa_sums(ps), decreasing = TRUE)[1:20])
ps_species_rel <- transform_sample_counts(ps, function(x) x/sum(x))

ps_species_melt <- psmelt(ps_species_rel) %>%
  mutate(
    Day_numeric = as.numeric(Day),
    Day_for_sorting = if_else(is.na(Day_numeric), 999, Day_numeric),
    Species_display = ifelse(OTU %in% top_species, species, "Other"),
    Treatment_display = case_when(
      Treatment == "For baseline data" ~ "NEWs",
      Treatment == "Control" ~ "Control",
      Treatment == "Tetracycline-treated" ~ "Tetracycline",
      Treatment == "Terra-Pro-treated" ~ "Terra-Pro",
      Treatment == "ASx-lysate-treated" ~ "ASx lysate",
      Treatment == "Apiary-Marked" ~ "Hive Reared",
      Treatment == "Microbiome in the comb samples" ~ "Comb",
      TRUE ~ Treatment
    ),
    Day = if_else(Treatment_display == "Comb", "17", as.character(Day)), # comb day
    Treatment_order = factor(Treatment_display,
                             levels = c("NEWs",
                                        "Comb",
                                        "Control",
                                        "Tetracycline",
                                        "Terra-Pro",
                                        "ASx lysate",
                                        "Hive Reared")),
    Sample_label = paste("Day", Day, "", Treatment_display, "Rep ", Rep)
  )

# species bar plot
p_species <- ps_species_melt %>%
  ggplot(aes(x = reorder(Sample_label, Day_for_sorting),
             y = Abundance, 
             fill = Species_display)) +
  geom_bar(stat= "identity", width = 0.8) +
  facet_grid(~ Treatment_order, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(
    axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(face = "italic"),
    legend.direction = "vertical", # make single coumn
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    panel.spacing = unit(1, "lines")
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(
    x = "Sample",
    y = "Relative Abundance",
    title = "Taxonomic Composition by Species",
    fill = "Species"
  ) +
  scale_fill_viridis_d(option = "turbo")

ggsave(file.path(figures_dir, "taxonomy_by_species.pdf"), 
       plot = p_species, width = 18, height = 8)

ggsave(file.path(figures_dir, "taxonomy_by_species.png"), 
       plot = p_species, width = 18, height = 8)
# end of taxonomy relative abundance figures

# this section looks at the composition of taxa/each sample
# taxonomic composition statistics
# dominant genus per sample - which genus is in the most samples and by how much?
# visually, based on bar plots, lactobacillus is dominant
dominant_genus_per_sample <- ps_genus_melt %>%
  group_by(Sample) %>%
  slice_max(Abundance, n = 1) %>%
  ungroup() %>%
  select(Sample, Sample_label, Treatment, Day, genus, Abundance) %>%
  arrange(Treatment, Day)

# lacto dominance - output to screen
lacto_dominated <- sum(dominant_genus_per_sample$genus == "Lactobacillus")
cat("Samples dominated by Lactobacillus:", lacto_dominated, "out of", 
    nrow(dominant_genus_per_sample), 
    "(", round(lacto_dominated/nrow(dominant_genus_per_sample)*100, 1), "%)\n\n")
print(table(dominant_genus_per_sample$genus))

# samples NOT dominated by lactobacillus - second runner up
non_lacto_dominated <- dominant_genus_per_sample %>%
  filter(genus != "Lactobacillus")
if (nrow(non_lacto_dominated) > 0) {
  cat("\nSamples NOT dominated by Lacto\n")
  print(non_lacto_dominated)
  
  # is lactobacillus in these samples? could be duds??
  lacto_in_non_dominated <- ps_genus_melt %>%
    filter(Sample %in% non_lacto_dominated$Sample, 
           genus == "Lactobacillus") %>%
    select(Sample, Sample_label, genus, Abundance) %>%
    arrange(Sample)
  print(lacto_in_non_dominated)
  # lacto rank in non-dominated samples - is it a runner up? 
  lacto_rank_in_non_dominated <- ps_genus_melt %>%
    filter(Sample %in% non_lacto_dominated$Sample) %>%
    group_by(Sample, Sample_label) %>%
    arrange(desc(Abundance)) %>%
    mutate(Rank = row_number()) %>%
    filter(genus == "Lactobacillus") %>% # can I add bombilactos here?
    ungroup() %>%
    select(Sample_label, genus, Abundance, Rank)
  
  cat("\n=== LACTOBACILLUS RANK IN NON-DOMINATED SAMPLES ===\n")
  print(lacto_rank_in_non_dominated)
}

# how diverse is the lactobacillus in the samples? 
# looking specifically for firm4 & firm5 species 
# to include bombilactobacillus, comment out taxa naming section 
# and re-run from beginning

# lacto species diversity - print to screen
lacto_diversity <- ps_species_melt %>%
  filter(genus == "Lactobacillus", Abundance > 0) %>%
  group_by(Sample, Sample_label, Treatment) %>%
  summarise(
    N_Lacto_Species = n_distinct(species),
    Total_Lacto_Abundance = sum(Abundance),
    .groups = "drop"
  ) %>%
  arrange(desc(N_Lacto_Species))

cat("\n=== LACTOBACILLUS SPECIES DIVERSITY ===\n")
cat("Mean Lactobacillus species per sample:", 
    round(mean(lacto_diversity$N_Lacto_Species), 1), 
    "±", round(sd(lacto_diversity$N_Lacto_Species), 1), "\n\n")

cat("Top 10 samples with most Lactobacillus species:\n")
print(head(lacto_diversity, 10))


write_csv(dominant_genus_per_sample, 
          file.path(output_dir, "dominant_genus_per_sample.csv"))

write_csv(lacto_diversity, 
          file.path(output_dir, "lactobacillus_diversity_per_sample.csv"))
# end of taxa analysis

# beta diversity - pcoa plots
cat("\n Starting Beta Diversity Analysis...\n")

# treatment colors consistent w/anuja's figures
sample_colors <- c(
  "NEWs" = "#984ea3",
  "Control" = "#377eb8",
  "Tetracycline" = "#4daf4a",
  "Terra-Pro" = "#ff7f00",
  "ASx lysate" = "#e41a1c",
  "Hive Reared" = "#ffff33",
  "Comb" = "#a65628"
)

# calculate bray-curtis distance
ps_bray <- phyloseq::distance(ps, method = "bray")

# pcoa
ps_pcoa <- ordinate(ps, method = "PCoA", distance = ps_bray)
# prep data for plotting
pcoa_df <- plot_ordination(ps, ps_pcoa, justDF = TRUE) %>%
  rownames_to_column("Sample") %>%
  mutate(
    Treatment_display = case_when(
      Treatment == "For baseline data" ~ "NEWs",
      Treatment == "Control" ~ "Control",
      Treatment == "Tetracycline-treated" ~ "Tetracycline",
      Treatment == "Terra-Pro-treated" ~ "Terra-Pro",
      Treatment == "ASx-lysate-treated" ~ "ASx lysate",
      Treatment == "Apiary-Marked" ~ "Hive Reared",
      Treatment == "Microbiome in the comb samples" ~ "Comb",
      TRUE ~ Treatment
    ),
    Day = if_else(Treatment_display == "Comb", "17", as.character(Day)), # comb day
    Treatment_order = factor(Treatment_display,
                             levels = c("NEWs",
                                        "Comb",
                                        "Control",
                                        "Tetracycline",
                                        "Terra-Pro",
                                        "ASx lysate",
                                        "Hive Reared")),
    rep_factor = factor(Rep)
  )
# pcoa of all samples
p_pcoa <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2,
                              color = Treatment_order,
                              shape = rep_factor)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = Treatment_order),
               type = "norm", 
               level = 0.95) +
  scale_color_manual(values = sample_colors) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(
    title = "PCoA of Bee Gut Bacterial Communities",
    subtitle = "Bray-Curtis dissimilarity",
    x = paste0("PC1 (", round(ps_pcoa$values$Relative_eig[1]*100, 1), "%)"),
    y = paste0("PC1 (", round(ps_pcoa$values$Relative_eig[2]*100, 1), "%)"),
    color = "Treatment",
    shape = "Replicate"
  )

ggsave(file.path(figures_dir, "pcoa_bray_curtis.pdf"),
       plot = p_pcoa, width = 10, height = 7)


ggsave(file.path(figures_dir, "pcoa_bray_curtis.png"),
       plot = p_pcoa, width = 10, height = 7)

# pcoa of all samples - with labels
p_pcoa_label <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2,
                              color = Treatment_order,
                              shape = rep_factor)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(aes(label = paste0("D", Day, "-R", Rep)),
            vjust = -1, size = 2.5, show.legend = FALSE) +
  stat_ellipse(aes(group = Treatment_order),
               type = "norm", 
               level = 0.95) +
  scale_color_manual(values = sample_colors) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(
    title = "PCoA of Bee Gut Bacterial Communities - Labeled",
    subtitle = "Bray-Curtis dissimilarity",
    x = paste0("PC1 (", round(ps_pcoa$values$Relative_eig[1]*100, 1), "%)"),
    y = paste0("PC1 (", round(ps_pcoa$values$Relative_eig[2]*100, 1), "%)"),
    color = "Treatment",
    shape = "Replicate"
  )

ggsave(file.path(figures_dir, "pcoa_bray_curtis_label.pdf"),
       plot = p_pcoa_label, width = 10, height = 7)


ggsave(file.path(figures_dir, "pcoa_bray_curtis_label.png"),
       plot = p_pcoa_label, width = 10, height = 7)

# end of pcoa

## DIFFERENTIAL ABUNDANCE STARTS HERE ##
# DESeq2 & ANCOM-BC2 require count data, not relative abundance
# read taxonomy-detailed.tsv files with estimated read counts

tsv_files <- list.files(data_dir, 
                        pattern = "taxonomy-detailed\\.tsv",
                        recursive = TRUE,
                        full.names = TRUE)

# exclude the excluded_samples folder
tsv_files <- tsv_files[!str_detect(tsv_files, "excluded_samples")]

# read & merge all count data
raw_counts <- map_dfr(tsv_files, function(file) {
  sample_name <- basename(dirname(file))
  tube_id <- str_extract(sample_name, "\\d+$")
  
  read_tsv(file, show_col_types = FALSE) %>%
    mutate(OriginalSampleName = sample_name,
           TubeID = tube_id)
})
# merge with metadata
raw_counts_named <- raw_counts %>%
  left_join(metadata_clean, by = c("TubeID" = "SAMPLE Tube ID")) %>%
  mutate(NewSampleName = str_remove(NewSampleName, " rerun"))

# otu table with counts
otu_table_counts <- raw_counts_named %>%
  select(tax_id, NewSampleName, `estimated counts`) %>%
  pivot_wider(names_from = NewSampleName,
              values_from = `estimated counts`,
              values_fill = 0) %>%
  column_to_rownames("tax_id")

# create taxonomy table
tax_table_counts_df <- raw_counts_named %>%
  select(tax_id, superkingdom, phylum, class, order, family, genus, species) %>%
  distinct()

# taxonomy corrections - reclass lactos Firm5
tax_table_counts_df <- tax_table_counts_df %>%
  mutate(
    genus = case_when(
      species == "Lactobacillus mellis" ~ "Bombilactobacillus",
      species == "Lactobacillus mellifer" ~ "Bombilactobacillus",
      TRUE ~ genus
    ),
    
    species = case_when(
      species == "Lactobacillus mellis" ~ "Bombilactobacillus mellis",
      species == "Lactobacillus mellifer" ~ "Bombilactobacillus mellifer",
      TRUE ~ species
    )
  ) %>%
  column_to_rownames("tax_id")

# create sample metadata
sample_metadata_counts <- raw_counts_named %>%
  select(NewSampleName, Day, Treatment, Rep) %>%
  distinct() %>%
  mutate(Day = as.numeric(Day)) %>%
  column_to_rownames("NewSampleName")

# convert to matrices
otu_mat_counts <- as.matrix(otu_table_counts)
tax_mat_counts <- as.matrix(tax_table_counts_df)

# phyloseq object with count data
ps_counts <- phyloseq(
  otu_table(otu_table_counts, taxa_are_rows = TRUE),
  tax_table(tax_mat_counts),
  sample_data(sample_metadata_counts)
)

# did this effing work?
cat("Samples:", nsamples(ps_counts), "\n")
cat("Taxa:", ntaxa(ps_counts), "\n")

# deseq2 starts here
# treatments and controls only - exclude others
ps_counts_treatments <- subset_samples(ps_counts,
                                       !Treatment %in% c("For baseline data",
                                                         "Apiary-Marked",
                                                         "Microbiome in the comb samples"))

# set control as the reference level
sample_data(ps_counts_treatments)$Treatment <- factor(
  sample_data(ps_counts_treatments)$Treatment,
  levels = c("Control", "Tetracycline-treated", "Terra-Pro-treated", "ASx-lysate-treated")
)

# run deseq2
deseq_treatments <- phyloseq_to_deseq2(ps_counts_treatments, ~ Treatment)
deseq_treatments <- estimateSizeFactors(deseq_treatments, type = "poscounts")
deseq_treatments <- DESeq(deseq_treatments, fitType = "parametric")

# extract results - control vs treatment for comparison
res_control_tet <- results(deseq_treatments,
                           contrast = c("Treatment", "Control",
                                        "Tetracycline-treated"), alpha = 0.5)

res_control_terra <- results(deseq_treatments,
                             contrast = c("Treatment", "Control",
                                          "Terra-Pro-treated"), alpha = 0.5)

res_control_asx <- results(deseq_treatments,
                           contrast = c("Treatment", "Control",
                                        "ASx-lysate-treated"), alpha = 0.5)

# extract significant taxa
extract_sig_taxa <- function(res_obj, comparison_name) {
  as.data.frame(res_obj) %>%
    rownames_to_column("tax_id") %>%
    left_join(tax_table_counts_df %>% 
                rownames_to_column("tax_id"), by = "tax_id") %>%
    filter(padj < 0.05) %>%
    mutate(Comparison = comparison_name,
           direction = ifelse(log2FoldChange > 0, 
                              "Enriched_in_Control", "Enriched_in_Treatment")) %>%
    select(Comparison, species, genus, family, log2FoldChange, padj, 
           baseMean, direction) %>%
    arrange(padj)
}

# extract significant taxa for all comparisons 
sig_tet <- extract_sig_taxa(res_control_tet, "Control_vs_Tetracycline")
sig_terra <- extract_sig_taxa(res_control_terra, "Control_vs_Terra-Pro")
sig_asx <- extract_sig_taxa(res_control_asx, "Control_vs_ASx_lysate")
# combine
all_sig_taxa_treatments <- rbind(sig_tet, sig_terra, sig_asx)

write_csv(all_sig_taxa_treatments, file.path(output_dir, 
                                             "deseq2_significant_taxa_treatments.csv"))
## end of deseq2

# ancom-bc2 starts here
# run ancombc2 (control is already set as reference)
ancom_treatments <- ancombc2(
  data = ps_counts_treatments,
  fix_formula = "Treatment",
  p_adj_method = "holm",
  prv_cut = 0.10,
  lib_cut = 1000,
  group = "Treatment",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  global = TRUE
)

# extract results
ancom_res <- ancom_treatments$res

sig_ancom_tet <- ancom_res %>%
  filter(`diff_TreatmentTetracycline-treated` == TRUE) %>%
  select(taxon, `lfc_TreatmentTetracycline-treated`, `q_TreatmentTetracycline-treated`) %>%
  rename(lfc = `lfc_TreatmentTetracycline-treated`,
         q_val = `q_TreatmentTetracycline-treated`) %>%
  mutate(Comparison = "Tetracycline vs Control") %>%
  left_join(tax_table_df %>% rownames_to_column("taxon"), by = "taxon")

sig_ancom_terra <- ancom_res %>%
  filter(`diff_TreatmentTerra-Pro-treated` == TRUE) %>%
  select(taxon, `lfc_TreatmentTerra-Pro-treated`, `q_TreatmentTerra-Pro-treated`) %>%
  rename(lfc = `lfc_TreatmentTerra-Pro-treated`, 
         q_val = `q_TreatmentTerra-Pro-treated`) %>%
  mutate(Comparison = "Terra-Pro vs Control") %>%
  left_join(tax_table_df %>% rownames_to_column("taxon"), by = "taxon")

sig_ancom_asx <- ancom_res %>%
  filter(`diff_TreatmentASx-lysate-treated` == TRUE) %>%
  select(taxon, `lfc_TreatmentASx-lysate-treated`, `q_TreatmentASx-lysate-treated`) %>%
  rename(lfc = `lfc_TreatmentASx-lysate-treated`,
         q_val = `q_TreatmentASx-lysate-treated`) %>%
  mutate(Comparison = "ASx lysate vs Control") %>%
  left_join(tax_table_df %>% rownames_to_column("taxon"), by = "taxon")

# combine results
all_sig_ancom <- rbind(
  sig_ancom_tet %>% select(Comparison, taxon, species, genus, lfc, q_val),
  sig_ancom_terra %>% select(Comparison, taxon, species, genus, lfc, q_val),
  sig_ancom_asx %>% select(Comparison, taxon, species, genus, lfc, q_val)
)

write_csv(all_sig_ancom, 
          file.path(output_dir, "ancombc2_significant_taxa.csv"))


# top candidates (lowest q-values)
ancom_borderline <- ancom_res %>%
  left_join(tax_table_counts_df %>% rownames_to_column("taxon"),
            by = "taxon", suffix = c("", ".y")) %>%  select(taxon, species, genus,
         `lfc_TreatmentTetracycline-treated`, `q_TreatmentTetracycline-treated`,
         `lfc_TreatmentTerra-Pro-treated`, `q_TreatmentTerra-Pro-treated`,
         `lfc_TreatmentASx-lysate-treated`, `q_TreatmentASx-lysate-treated`) %>%
  select(-ends_with(".y")) %>%
  filter(`q_TreatmentTetracycline-treated` < 0.1 | 
           `q_TreatmentTerra-Pro-treated` < 0.1 | 
           `q_TreatmentASx-lysate-treated` < 0.1)

if (nrow(ancom_borderline) > 0) {
  cat("Taxa with q < 0.1 (borderline):", nrow(ancom_borderline), "\n")
  print(ancom_borderline)
  write_csv(ancom_borderline, file.path(output_dir, "ancombc2_borderline_taxa.csv"))
} 
# end of ancombc2


# compare deseq2 & ancombc2

# prep deseq2 results
deseq2_for_compare <- all_sig_taxa_treatments %>%
  mutate(Method = "DESeq2") %>%
  rename(lfc = log2FoldChange, q_val = padj) %>%
  select(Method, Comparison, species, genus, lfc, q_val)

# prep ancombc2 (including borderline)
ancombc2_for_compare <- rbind(
  ancom_res %>% #tetracycline
    left_join(tax_table_counts_df %>% rownames_to_column("taxon"), by = "taxon") %>%
    select(taxon, species, genus,
           `lfc_TreatmentTetracycline-treated`, `q_TreatmentTetracycline-treated`) %>%
    rename(lfc = `lfc_TreatmentTetracycline-treated`, 
           q_val = `q_TreatmentTetracycline-treated`) %>%
    mutate(Comparison = "Control_vs_Tetracycline", Method = "ANCOM-BC2") %>%
    filter(q_val < 0.1),
  
  ancom_res %>% # terra-pro
    left_join(tax_table_counts_df %>% rownames_to_column("taxon"), by = "taxon") %>%
    select(taxon, species, genus,
           `lfc_TreatmentTerra-Pro-treated`, `q_TreatmentTerra-Pro-treated`) %>%
    rename(lfc = `lfc_TreatmentTerra-Pro-treated`, 
           q_val = `q_TreatmentTerra-Pro-treated`) %>%
    mutate(Comparison = "Control_vs_Terra-Pro", Method = "ANCOM-BC2") %>%
    filter(q_val < 0.1),
  
  ancom_res %>% # asx lysate
    left_join(tax_table_counts_df %>% rownames_to_column("taxon"), by = "taxon") %>%
    select(taxon, species, genus,
           `lfc_TreatmentASx-lysate-treated`, `q_TreatmentASx-lysate-treated`) %>%
    rename(lfc = `lfc_TreatmentASx-lysate-treated`, 
           q_val = `q_TreatmentASx-lysate-treated`) %>%
    mutate(Comparison = "Control_vs_ASx_lysate", Method = "ANCOM-BC2") %>%
    filter(q_val < 0.1)
) %>%
  left_join(tax_table_counts_df %>% rownames_to_column("taxon"), 
            by = "taxon", suffix = c("", ".y")) %>%
  select(Method, Comparison, species, genus, lfc, q_val)
  
comparison_table <- rbind(
  deseq2_for_compare, ancombc2_for_compare) %>% 
  arrange(Comparison, species, Method)

write_csv(comparison_table, file.path(output_dir, 
                                      "deseq2_vs_ancombc2_comparison.csv"))


# summary table
summary_comparison <- data.frame(
  Comparison = c("Control_vs_Tetracycline", "Control_vs_Terra-Pro", "Control_vs_ASx_lysate"),
  DESeq2_Significant = c(
    sum(all_sig_taxa_treatments$Comparison == "Control_vs_Tetracycline"),
    sum(all_sig_taxa_treatments$Comparison == "Control_vs_Terra-Pro"),
    sum(all_sig_taxa_treatments$Comparison == "Control_vs_ASx_lysate")
  ),
  ANCOMBC2_Significant = c(
    sum(sig_ancom_tet$Comparison == "Tetracycline vs Control"),
    sum(sig_ancom_terra$Comparison == "Terra-Pro vs Control"),
    sum(sig_ancom_asx$Comparison == "ASx lysate vs Control")
  ),
  ANCOMBC2_Borderline = c(
    sum(ancom_borderline$`q_TreatmentTetracycline-treated` < 0.1, na.rm = TRUE),
    sum(ancom_borderline$`q_TreatmentTerra-Pro-treated` < 0.1, na.rm = TRUE),
    sum(ancom_borderline$`q_TreatmentASx-lysate-treated` < 0.1, na.rm = TRUE)
  )
) # tired of typing

write_csv(summary_comparison, 
          file.path(output_dir, "methods_comparison_summary.csv"))

# end of deseq2 & ancombc2 comparison

# abundance rank analysis - are sig taxa bee gut specific?
# top 20 most abundant species
top_20_species <- names(sort(taxa_sums(ps_counts), decreasing = TRUE)[1:20])
top_20_names <- tax_table_counts_df[top_20_species, "species"]

# are deseq2 significant taxa in top 20? - print to screen
deseq_sig_taxa <- all_sig_taxa_treatments %>%
  mutate(In_Top_20 = species %in% top_20_names) %>%
  select(Comparison, species, genus, In_Top_20, log2FoldChange, padj) %>%
  arrange(Comparison, desc(In_Top_20))

print(table(deseq_sig_taxa$In_Top_20, deseq_sig_taxa$Comparison))

# taxa in top20
deseq_in_top20 <- deseq_sig_taxa %>% filter(In_Top_20 == TRUE)

if (nrow(deseq_in_top20) > 0) {
  cat("\nDESeq2 significant taxa IN top 20:\n")
  for (i in 1:nrow(deseq_in_top20)) {
    cat("  ", deseq_in_top20$species[i], " (", deseq_in_top20$Comparison[i], ")\n", sep="")
  }
} else {
  cat("\nNone of the DESeq2 significant taxa are in the top 20 most abundant.\n")
}

# taxa not in top 20
deseq_not_top20 <- deseq_sig_taxa %>% filter(In_Top_20 == FALSE)

cat("\nDESeq2 significant taxa NOT in top 20:", nrow(deseq_not_top20), "\n")
if (nrow(deseq_not_top20) > 0) {
  for (i in 1:nrow(deseq_not_top20)) {
    cat("  ", deseq_not_top20$species[i], " (", deseq_not_top20$Comparison[i], ")\n", sep="")
  }
}

# ancombc2 borderline taxon
if (nrow(ancom_borderline) > 0) {
  ancom_borderline_species <- unique(ancom_borderline$species)
  ancom_in_top20 <- ancom_borderline_species %in% top_20_names
  for (i in 1:length(ancom_borderline_species)) {
    cat(ancom_borderline_species[i], ": ", 
        ifelse(ancom_in_top20[i], "In top 20", "Not in top 20"), "\n", sep="")
  }
}

write_csv(deseq_sig_taxa, 
          file.path(output_dir, "significant_taxa_abundance_rank.csv"))



# community similarity - based on taxa in each sample, how similar/dissimilar
# are the samples 
# use jaccard index (0 = no overlap 1 = identical sets)
# use treatment_taxa_lists_all from before 
if (!exists("treatment_taxa_lists_all")) {
  treatment_taxa_lists_all <- list()
  
  all_treatments <- c("For baseline data", 
                      "Control", 
                      "Tetracycline-treated", 
                      "Terra-Pro-treated", 
                      "ASx-lysate-treated", 
                      "Apiary-Marked", 
                      "Microbiome in the comb samples")
  
  for(treat in all_treatments) {
    samples_in_treatment <- data.frame(sample_data(ps_counts)) %>%
      rownames_to_column("SampleID") %>%
      filter(Treatment == treat) %>%
      pull(SampleID)
    
    taxa_present <- otu_table(ps_counts)[, samples_in_treatment] %>%
      as.matrix() %>%
      rowSums() > 0
    
    treatment_taxa_lists_all[[treat]] <- names(taxa_present[taxa_present])
  }
  
  # clean up names
  names(treatment_taxa_lists_all) <- c("NEWs", 
                                       "Control", 
                                       "Tetracycline", 
                                       "Terra-Pro", 
                                       "ASx lysate", 
                                       "Hive Reared", 
                                       "Comb")
}

# jaccard distances between sample types
pa_matrix <- matrix(0, 
                    nrow = length(unique(unlist(treatment_taxa_lists_all))),
                    ncol = length(treatment_taxa_lists_all))
colnames(pa_matrix) <- names(treatment_taxa_lists_all)
rownames(pa_matrix) <- unique(unlist(treatment_taxa_lists_all))

for(i in 1:length(treatment_taxa_lists_all)) {
  pa_matrix[treatment_taxa_lists_all[[i]], i] <- 1
}

jaccard_dist <- vegdist(t(pa_matrix), method = "jaccard")
jaccard_matrix <- as.matrix(jaccard_dist)

write.csv(jaccard_matrix, 
          file.path(output_dir, "jaccard_distances_sample_types.csv"), row.names = TRUE)

# summary statistics
core_all <- Reduce(intersect, treatment_taxa_lists_all)

summary_stats <- data.frame(
  Sample_Type = names(treatment_taxa_lists_all),
  Total_Taxa = sapply(treatment_taxa_lists_all, length),
  Unique_Taxa = sapply(names(treatment_taxa_lists_all), function(x) {
    length(setdiff(treatment_taxa_lists_all[[x]], 
                   unlist(treatment_taxa_lists_all[names(treatment_taxa_lists_all) != x])))
  }),
  Core_Taxa_All = length(core_all)
)

# caged core taxa count
caged_names <- c("Control", "Tetracycline", "Terra-Pro", "ASx lysate")
core_caged <- Reduce(intersect, treatment_taxa_lists_all[caged_names])

summary_stats$Core_Taxa_Caged <- ifelse(
  summary_stats$Sample_Type %in% caged_names,
  length(core_caged), NA)

write.csv(summary_stats, 
          file.path(output_dir, "taxa_summary_statistics.csv"),
          row.names = FALSE)


# comb as bacterial source & hive vs control comparison starts here
# want to know how efficiently the comb seeded the control & treatment samples
# want to know how the hive reared compares to comb and controls (no treatments)

# comb seeding effectiveness - how many shared taxa with comb and treatments/control
comb_taxa <- treatment_taxa_lists_all$Comb

for(treat in c("Control", "Tetracycline", "Terra-Pro", "ASx lysate")) {
  shared <- intersect(comb_taxa, treatment_taxa_lists_all[[treat]])
  percent_shared <- (length(shared) / length(treatment_taxa_lists_all[[treat]])) * 100
  
  cat(treat, ": ", length(shared), " taxa shared with Comb (", 
      round(percent_shared, 1), "% of microbiome)\n", sep="")
}

# hive vs comb vs control comparison - print to screen
hive_taxa <- treatment_taxa_lists_all$`Hive Reared`
control_taxa <- treatment_taxa_lists_all$Control

# shared between all three
shared_all_three <- Reduce(intersect, list(hive_taxa, comb_taxa, control_taxa))
cat("Shared by Hive, Comb, AND Control:", length(shared_all_three), "taxa\n")

# two-way intersections (excluding the third)
shared_hive_comb_only <- setdiff(intersect(hive_taxa, comb_taxa), control_taxa)
shared_hive_control_only <- setdiff(intersect(hive_taxa, control_taxa), comb_taxa)
shared_comb_control_only <- setdiff(intersect(comb_taxa, control_taxa), hive_taxa)

cat("Hive + Comb only (not Control):", length(shared_hive_comb_only), "taxa\n")
cat("Hive + Control only (not Comb):", length(shared_hive_control_only), "taxa\n")
cat("Comb + Control only (not Hive):", length(shared_comb_control_only), "taxa\n")

# unique to each
unique_hive <- setdiff(hive_taxa, c(comb_taxa, control_taxa))
unique_comb <- setdiff(comb_taxa, c(hive_taxa, control_taxa))
unique_control <- setdiff(control_taxa, c(hive_taxa, comb_taxa))

cat("\nUnique to Hive only:", length(unique_hive), "taxa\n")
cat("Unique to Comb only:", length(unique_comb), "taxa\n")
cat("Unique to Control only:", length(unique_control), "taxa\n")

# taxa in hive but not from comb
hive_not_from_comb <- setdiff(hive_taxa, comb_taxa)

if (length(hive_not_from_comb) > 0) {
  cat("\nTaxa in Hive-Reared but NOT in Comb:", length(hive_not_from_comb), "\n")
  hive_not_comb_species <- tax_table_counts_df[hive_not_from_comb, c("genus", "species")]
  print(hive_not_comb_species)
}

# overlap percentages
hive_from_comb_pct <- (length(intersect(hive_taxa, comb_taxa)) / length(hive_taxa)) * 100
control_from_comb_pct <- (length(intersect(control_taxa, comb_taxa)) / length(control_taxa)) * 100

cat("\n=== MICROBIOME SOURCE COMPARISON ===\n")
cat("Hive-Reared:", round(hive_from_comb_pct, 1), "% overlap with Comb\n")
cat("Control (caged):", round(control_from_comb_pct, 1), "% overlap with Comb\n")
cat("Difference:", round(abs(hive_from_comb_pct - control_from_comb_pct), 1), "percentage points\n")

# save summary
venn_summary <- data.frame(
  Category = c("All three", "Hive + Comb only", "Hive + Control only", 
               "Comb + Control only", "Hive only", "Comb only", "Control only"),
  Count = c(length(shared_all_three), length(shared_hive_comb_only), 
            length(shared_hive_control_only), length(shared_comb_control_only),
            length(unique_hive), length(unique_comb), length(unique_control))
)

write.csv(venn_summary, 
          file.path(output_dir, "hive_vs_cage_taxa_overlap.csv"),
          row.names = FALSE)

### differential abundance comb vs hive & controls
# subset to comb, hive, control 
ps_counts_comb_comparison <- subset_samples(ps_counts, 
                                            Treatment %in% c("Microbiome in the comb samples",
                                                             "Apiary-Marked",
                                                             "Control"))

# comb as reference
sample_data(ps_counts_comb_comparison)$Treatment <- factor(
  sample_data(ps_counts_comb_comparison)$Treatment,
  levels = c("Microbiome in the comb samples", "Apiary-Marked", "Control")
)

# deseq2
deseq_comb <- phyloseq_to_deseq2(ps_counts_comb_comparison, ~ Treatment)
deseq_comb <- estimateSizeFactors(deseq_comb, type = "poscounts")
deseq_comb <- DESeq(deseq_comb, fitType = "parametric")

# extract results
res_comb_hive <- results(deseq_comb,
                         contrast = c("Treatment", 
                                      "Microbiome in the comb samples", 
                                      "Apiary-Marked"),
                         alpha = 0.05)

res_comb_control <- results(deseq_comb,
                            contrast = c("Treatment", 
                                         "Microbiome in the comb samples", 
                                         "Control"),
                            alpha = 0.05)

# significant taxa
sig_comb_hive <- extract_sig_taxa(res_comb_hive, "Comb_vs_Hive_Reared")
sig_comb_control <- extract_sig_taxa(res_comb_control, "Comb_vs_Control")

# combine
all_sig_comb_comparison <- rbind(sig_comb_hive, sig_comb_control)

write_csv(all_sig_comb_comparison, 
          file.path(output_dir, "deseq2_comb_comparison.csv"))

cat("Comb vs Hive-Reared:", nrow(sig_comb_hive), "significant taxa\n")
cat("Comb vs Control:", nrow(sig_comb_control), "significant taxa\n")

# ancombc2 
ancom_comb <- ancombc2(
  data = ps_counts_comb_comparison,
  fix_formula = "Treatment",
  p_adj_method = "holm",
  prv_cut = 0.10,
  lib_cut = 1000,
  group = "Treatment",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  global = TRUE
)

ancom_comb_res <- ancom_comb$res

# ancombc2 significant taxa
sig_ancom_comb_hive <- ancom_comb_res %>%
  left_join(tax_table_counts_df %>% rownames_to_column("taxon"), by = "taxon") %>%
  filter(`diff_TreatmentApiary-Marked` == TRUE) %>%
  select(taxon, species, genus, `lfc_TreatmentApiary-Marked`, `q_TreatmentApiary-Marked`) %>%
  rename(lfc = `lfc_TreatmentApiary-Marked`, q_val = `q_TreatmentApiary-Marked`) %>%
  mutate(Comparison = "Comb vs Hive-Reared")

sig_ancom_comb_control <- ancom_comb_res %>%
  left_join(tax_table_counts_df %>% rownames_to_column("taxon"), by = "taxon") %>%
  filter(`diff_TreatmentControl` == TRUE) %>%
  select(taxon, species, genus, `lfc_TreatmentControl`, `q_TreatmentControl`) %>%
  rename(lfc = `lfc_TreatmentControl`, q_val = `q_TreatmentControl`) %>%
  mutate(Comparison = "Comb vs Control")

# combine
all_sig_ancom_comb <- rbind(sig_ancom_comb_hive, sig_ancom_comb_control)

write_csv(all_sig_ancom_comb, 
          file.path(output_dir, "ancombc2_comb_comparison.csv"))
# end of comb vs hive vs control differential


# are the bee gut specific bacteria in all the samples??
# gut microbiome genus
core_bee_bacteria <- c(
  "Snodgrassella", 
  "Gilliamella",
  "Lactobacillus",
  "Bombilactobacillus",
  "Bifidobacterium",
  "Frischella",
  "Bartonella"
)

# presence/absence
core_bacteria_presence <- ps_genus_melt %>%
  filter(genus %in% core_bee_bacteria) %>%
  group_by(Treatment_display, genus) %>%
  summarize(
    N_Samples_Present = sum(Abundance > 0),
    Total_Samples = n_distinct(Sample),
    Percent_Present = (N_Samples_Present / Total_Samples) * 100,
    Mean_Abundance = mean(Abundance[Abundance > 0], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Treatment_display, genus)

print(core_bacteria_presence)

write_csv(core_bacteria_presence, file.path(output_dir, "core_bacteria_presence_by_treatment.csv"))

#fml forgot the overlap percentages csv
comb_taxa <- treatment_taxa_lists_all$Comb
hive_taxa <- treatment_taxa_lists_all$`Hive Reared`
control_taxa <- treatment_taxa_lists_all$Control
tet_taxa <- treatment_taxa_lists_all$Tetracycline
terra_taxa <- treatment_taxa_lists_all$`Terra-Pro`
asx_taxa <- treatment_taxa_lists_all$`ASx lysate`

comb_seeding_table <- data.frame(
  Sample_Type = c("Control", "Tetracycline", "Terra-Pro", "ASx lysate", "Hive Reared"),
  Taxa_Shared_With_Comb = c(
    length(intersect(control_taxa, comb_taxa)),
    length(intersect(tet_taxa, comb_taxa)),
    length(intersect(terra_taxa, comb_taxa)),
    length(intersect(asx_taxa, comb_taxa)),
    length(intersect(hive_taxa, comb_taxa))
  ),
  Total_Taxa = c(
    length(control_taxa),
    length(tet_taxa),
    length(terra_taxa),
    length(asx_taxa),
    length(hive_taxa)
  )
) %>%
  mutate(Percent_From_Comb = round((Taxa_Shared_With_Comb / Total_Taxa) * 100, 1))

write.csv(comb_seeding_table, 
          file.path(output_dir, "comb_seeding_effectiveness.csv"),
          row.names = FALSE)
# end bee gut specific genus


# summaries section 
data_summary_text <- paste0(
  "16S Microbiome Analysis - Data Import & Sequencing Summary\n",
  "Date: ", Sys.Date(), "\n\n",
  
  "Dataset Overview\n\n",
  "Total Samples: ", length(unique(raw_data_named$NewSampleName)), "\n",
  "Total Unique Taxa: ", length(unique(raw_data_named$`Taxonomy ID`)), "\n",
  "Total Observations: ", nrow(raw_data), "\n\n",
  
  "Sequencing Statistics\n\n",
  "Mean Reads per Sample: ", round(sequencing_summary$Mean_Reads, 0),
  " ± ", round(sequencing_summary$SD_Reads, 0), "\n",
  "Range: ", round(sequencing_summary$Min_Reads, 0), " - ",
  round(sequencing_summary$Max_Reads, 0), " reads\n",
  "Mean Read Length: ", round(sequencing_summary$Mean_Read_Length, 1), " bp\n",
  "Mean Quality Score: ", round(sequencing_summary$Mean_Quality, 2), "\n",
  "Mean N50: ", round(sequencing_summary$Mean_N50, 0), " bp\n\n",
  
  "Experimental Design\n\n",
  "Treatments:\n",
  paste(" - ", unique(metadata_clean$Treatment), collapse = "\n"), "\n\n",
  "Time Points (days): ", paste(sort(unique(as.numeric(metadata_clean$Day))), 
                                collapse = ", "), "\n\n",
  
  "Samples per Treatment\n"
)

# add sample counts per treatment & day
treatment_counts <- metadata_clean %>%
  dplyr::count(Treatment, Day) %>%
  arrange(Day, Treatment)

for (i in 1:nrow(treatment_counts)) {
  data_summary_text <- paste0(data_summary_text, 
                              " Day ", treatment_counts$Day[i], " - ",
                              treatment_counts$Treatment[i], ": ",
                              treatment_counts$n[i], " replicates\n")
}

# add sequencing stats by treatment
data_summary_text <- paste0(data_summary_text, 
                            "\n\n Sequencing Stats by Treatment\n\n")

for(i in 1:nrow(treatment_summary)) {
  data_summary_text <- paste0(data_summary_text,
                              "  ", treatment_summary$Treatment[i], ":\n",
                              "    N = ", treatment_summary$N_Samples[i], "\n",
                              "    Mean Reads = ", round(treatment_summary$Mean_Reads[i], 0), 
                              " ± ", round(treatment_summary$SD_Reads[i], 0), "\n",
                              "    Mean Quality = ", round(treatment_summary$Mean_Quality[i], 2), "\n\n"
  )
}

# output files
data_summary_text <- paste0(data_summary_text,
                            "Output Files\n",
                            "Sample Information:\n",
                            " - output/sample_name_lookup.csv (sample ID crosswalk)\n",
                            " - output/metadata_clean.csv (cleaned metadata)\n\n",
                            "Sequencing Statistics:\n",
                            " - output/sequencing_stats_all_samples.csv (all samples)\n",
                            " - output/sequencing_stats_by_treatment.csv (by treatment)\n\n",
                            "Summary:\n",
                            " - notes/data_import_summary.txt (ths file)"
)

writeLines(data_summary_text, file.path(notes_dir, "data_import_summary.txt"))

analysis_summary <- paste0(
  "16S Microbiome Analysis - Complete Summary \n",
  "Date: ", Sys.Date(), "\n",

  "DATASET OVERVIEW\n",
  "----------------\n",
  "Total samples: ", nsamples(ps_counts), "\n",
  "Total unique taxa: ", ntaxa(ps_counts), "\n",
  "Mean reads per sample: ", round(mean(ont_stats_named$total_num_reads), 0), 
  " ± ", round(sd(ont_stats_named$total_num_reads), 0), "\n",
  "Mean read length: ", round(mean(ont_stats_named$mean_read_length), 1), " bp\n",
  "Mean quality score: ", round(mean(ont_stats_named$mean_quality), 2), "\n\n",
  
  "Sample Types\n",
  "------------\n",
  "NEWs: Newly emerged bees (sterile baseline)\n",
  "Control: Caged, no treatment\n",
  "Tetracycline: Caged, antibiotic treatment\n",
  "Terra-Pro: Caged, probiotic treatment\n",
  "ASx lysate: Caged, phage lysate treatment\n",
  "Hive Reared: Marked and returned to hive\n",
  "Comb: Comb material provided in cages\n\n",
  
  "Figures Generated\n",
  
  "Taxonomic Composition\n",
  "---------------------\n",
  "1. taxonomy_by_genus.pdf/.png\n",
  "   - Top 20 genera, relative abundance by treatment and day\n\n",
  
  "2. taxonomy_by_species.pdf/.png\n",
  "   - Top 20 species, relative abundance by treatment and day\n",
  "   - Includes Bombilactobacillus reclassification\n\n",
  
  "Beta Diversity (PCoA)\n",
  "---------------------\n",
  "3. pcoa_bray_curtis.pdf/.png\n",
  "   - All samples, Bray-Curtis dissimilarity\n",
  "   - PC1: ", round(ps_pcoa$values$Relative_eig[1]*100, 1), "% variance\n",
  "   - PC2: ", round(ps_pcoa$values$Relative_eig[2]*100, 1), "% variance\n",
  "   - Total captured: ", round(sum(ps_pcoa$values$Relative_eig[1:2])*100, 1), "%\n\n",
  
  "4. pcoa_bray_curtis_labeled.pdf/.png\n",
  "   - Same as #3 with sample labels (Day-Rep)\n\n",
  
  "Statistical Analysis\n",
  
  "Taxonomic Composition\n",
  
  "Lactobacillus dominance: ", lacto_dominated, "/", nrow(dominant_genus_per_sample), 
  " samples (", round(lacto_dominated/nrow(dominant_genus_per_sample)*100, 1), "%)\n",
  "Mean Lactobacillus species per sample: ", round(mean(lacto_diversity$N_Lacto_Species), 1), 
  " ± ", round(sd(lacto_diversity$N_Lacto_Species), 1), "\n",
  "Gilliamella: Second most abundant genus\n\n",
  
  "DESeq2 Results (vs Control)\n",
  "Tetracycline vs Control: ", nrow(sig_tet), " significant taxa\n",
  "Terra-Pro vs Control: ", nrow(sig_terra), " significant taxa\n",
  "ASx lysate vs Control: ", nrow(sig_asx), " significant taxa\n",
  "Significant taxa abundance: All are rare (not in top 20)\n\n",
  
  "ANCOM-BC2 Results (vs Control)\n",
  "More conservative compositional method\n",
  "Tetracycline vs Control: ", nrow(sig_ancom_tet), " significant taxa\n",
  "Terra-Pro vs Control: ", nrow(sig_ancom_terra), " significant taxa\n",
  "ASx lysate vs Control: ", nrow(sig_ancom_asx), " significant taxa\n",
  "Borderline (q < 0.1): ", nrow(ancom_borderline), " taxa\n\n",
  
  "Community Similarity (Jaccard Distance)\n",
  "Most similar:\n",
  "  - Terra-Pro & Tetracycline: 0.269 (both antimicrobials)\n",
  "  - Control & ASx lysate: 0.372\n\n",
  "Most different:\n",
  "  - NEWs: 0.65-0.72 from all others (sterile baseline)\n",
  "  - Hive Reared: 0.43-0.69 from caged (environment matters)\n\n",
  
  "Comb as Bacterial Source\n",
  "Caged bees acquired from comb:\n",
  "  - Control: 58.1%\n",
  "  - Tetracycline: 61.9%\n",
  "  - Terra-Pro: 62.5%\n",
  "  - ASx lysate: 51.3%\n\n",
  
  "Hive-reared vs Caged:\n",
  "  - Hive-Reared: 91.7% overlap with comb\n",
  "  - Control: 58.1% overlap with comb\n",
  "  - Difference: 33.6 percentage points\n",
  "  - Despite high overlap, abundances differ significantly\n\n",
  
  "Core Bee Gut Bacteria\n",
  "Analyzed genera: Lactobacillus, Bombilactobacillus, Gilliamella,\n",
  "                 Snodgrassella, Bifidobacterium, Frischella, Bartonella\n",
  "See: core_bacteria_presence_by_treatment.csv for details\n\n",
  
 
  "KEY FINDINGS\n",
  
  "1. Antibiotic Controls Worked:\n",
  "   - Tetracycline: ", nrow(sig_tet), " significantly different taxa\n",
  "   - Terra-Pro: ", nrow(sig_terra), " significantly different taxa\n",
  "   - Cluster together in PCoA (similar antimicrobial effect)\n\n",
  
  "2. ASx Lysate Had Comparable Effect:\n",
  "   - ASx lysate: ", nrow(sig_asx), " significantly different taxon (DESeq2)\n",
  "   - ", nrow(sig_ancom_asx), " significant taxa (ANCOM-BC2)\n",
  "   - Microbiome composition similar to untreated controls\n\n",
  
  "3. Caging Affects Microbiome:\n",
  "   - Hive-reared bees cluster separately from caged controls\n",
  "   - Environment influences community structure beyond colonization\n\n",
  
  "4. Comb Successfully Seeded Microbiome:\n",
  "   - 50-62% of caged bee gut bacteria came from comb\n",
  "   - Validates experimental design for controlled establishment\n\n",
  
  "5. Hive Environment Matters:\n",
  "   - Hive-reared: 91.7% taxa overlap with comb\n",
  "   - Caged: 58.1% taxa overlap with comb\n",
  "   - Same taxa, different abundances = environment shapes development\n\n",
  

  "Output Files\n",
  
  "Figures (figures/)\n",
  "------------------\n",
  "All figures saved as both PDF and PNG\n",
  "- Taxonomic composition (genus and species)\n",
  "- PCoA plots (labeled and unlabeled)\n\n",
  
  "Data Tables (output/)\n",
  "Sample Information:\n",
  "  - sample_name_lookup.csv\n",
  "  - metadata_clean.csv\n\n",
  
  "Sequencing Statistics:\n",
  "  - sequencing_stats_all_samples.csv\n",
  "  - sequencing_stats_by_treatment.csv\n\n",
  
  "Taxonomic Composition:\n",
  "  - dominant_genus_per_sample.csv\n",
  "  - lactobacillus_diversity_per_sample.csv\n",
  "  - core_bacteria_presence_by_treatment.csv\n",
  "  - core_bacteria_hive_vs_cage.csv\n\n",
  
  "Differential Abundance:\n",
  "  - deseq2_significant_taxa.csv\n",
  "  - ancombc2_significant_taxa.csv\n",
  "  - ancombc2_borderline_taxa.csv\n",
  "  - deseq2_vs_ancombc2_comparison.csv\n",
  "  - methods_comparison_summary.csv\n",
  "  - significant_taxa_abundance_rank.csv\n",
  "  - deseq2_comb_comparison.csv\n",
  "  - ancombc2_comb_comparison.csv\n\n",
  
  "Community Analysis:\n",
  "  - jaccard_distances_sample_types.csv\n",
  "  - taxa_summary_statistics.csv\n",
  "  - hive_vs_cage_taxa_overlap.csv\n\n",
  
  "Summary Notes (notes/)\n",
  "  - data_import_summary.txt\n",
  "  - analysis_complete_summary.txt (this file)\n\n",
  "R Objects\n",
  "  - phyloseq_object.rds (relative abundance)\n\n"
)

# save
writeLines(analysis_summary, file.path(notes_dir, "analysis_complete_summary.txt"))


# consolidate all of the csv files for supplemental
wb <- createWorkbook()

# readme sheet with table descriptions
table_descriptions <- tibble::tibble(
  Sheet_Name = c(
    "sample_lookup", "metadata_clean", "sequencing_stats_all", "sequencing_stats_treat",
    "dominant_genus", "lacto_diversity", "core_bacteria_presence", "core_bacteria_hive_cage",
    "deseq2_sig_taxa", "ancombc2_sig_taxa", "ancombc2_borderline", 
    "deseq2_v_ancombc2", "methods_comparison", "sig_taxa_abundance",
    "deseq2_comb_comp", "ancombc2_comb_comp",
    "jaccard_dist", "taxa_summary_stat", "hive_vs_cage", "comb_seeding_effectiveness"
  ),
  Description = c(
    "Sample name lookup table matching original sample names to cleaned names with metadata",
    "Cleaned metadata with standardized sample names, treatment groups, days, and replicates",
    "Complete sequencing statistics for all samples including read counts, length, and quality",
    "Summary of sequencing statistics grouped by treatment type",
    "Dominant genus per sample with abundance percentages",
    "Lactobacillus species diversity per sample",
    "Presence and abundance of core bee gut bacteria by treatment",
    "Core bacteria comparison between hive-reared, caged, and comb samples",
    "Significantly different taxa identified by DESeq2 (vs Control, padj < 0.05)",
    "Significantly different taxa identified by ANCOM-BC2 (vs Control, q < 0.05)",
    "ANCOM-BC2 borderline significant taxa (q < 0.1)",
    "Side-by-side comparison of DESeq2 vs ANCOM-BC2 results",
    "Summary comparing number of significant taxa found by each statistical method",
    "DESeq2 significant taxa with abundance rank (top 20 vs rare)",
    "DESeq2 differential abundance: Comb vs Hive-Reared and Comb vs Control",
    "ANCOM-BC2 differential abundance: Comb vs Hive-Reared and Comb vs Control",
    "Jaccard distance matrix showing community similarity (0=identical, 1=completely different)",
    "Summary statistics: total taxa, unique taxa, and core microbiome counts by sample type",
    "Venn-style summary of taxa overlap between Hive-Reared, Comb, and Control",
    "Table showing Comb as bacterial source analysis"
  )
)

# add readme sheet
addWorksheet(wb, "README", gridLines = FALSE)
writeData(wb, "README", "SUPPLEMENTARY DATA - TABLE DESCRIPTIONS", startRow = 1)
titleStyle <- createStyle(fontSize = 14, textDecoration = "bold")
addStyle(wb, "README", titleStyle, rows = 1, cols = 1)

writeData(wb, "README", table_descriptions, startRow = 3)
headerStyle <- createStyle(textDecoration = "bold", border = "bottom", fontSize = 11)
addStyle(wb, "README", headerStyle, rows = 3, cols = 1:2, gridExpand = TRUE)
setColWidths(wb, "README", cols = 1, widths = 25)
setColWidths(wb, "README", cols = 2, widths = 90)


# csv files to combine
csv_files_to_combine <- list(
  "sample_lookup" = "sample_name_lookup.csv",
  "metadata_clean" = "metadata_clean.csv",
  "sequencing_stats_all" = "sequencing_stats_all_samples.csv",
  "sequencing_stats_treat" = "sequencing_stats_by_treatment.csv",
  "dominant_genus" = "dominant_genus_per_sample.csv",
  "lacto_diversity" = "lactobacillus_diversity_per_sample.csv",
  "core_bacteria_presence" = "core_bacteria_presence_by_treatment.csv",
  "core_bacteria_hive_cage" = "hive_vs_cage_taxa_overlap.csv",
  "deseq2_sig_taxa" = "deseq2_significant_taxa_treatments.csv",
  "ancombc2_sig_taxa" = "ancombc2_significant_taxa.csv",
  "ancombc2_borderline" = "ancombc2_borderline_taxa.csv",
  "deseq2_v_ancombc2" = "deseq2_vs_ancombc2_comparison.csv",
  "methods_comparison" = "methods_comparison_summary.csv",
  "sig_taxa_abundance" = "significant_taxa_abundance_rank.csv",
  "deseq2_comb_comp" = "deseq2_comb_comparison.csv",
  "ancombc2_comb_comp" = "ancombc2_comb_comparison.csv",
  "jaccard_dist" = "jaccard_distances_sample_types.csv",
  "taxa_summary_stat" = "taxa_summary_statistics.csv",
  "hive_vs_cage" = "hive_vs_cage_taxa_overlap.csv",
  "comb_seeding" = "comb_seeding_effectiveness.csv"
)

# add csv as sheet
for (sheet_name in names(csv_files_to_combine)) {
  file_path <- file.path(output_dir, csv_files_to_combine[[sheet_name]])
  
  if (file.exists(file_path)) {
    data <- read.csv(file_path)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, data)
    
    headerStyle <- createStyle(textDecoration = "bold", border = "bottom")
    addStyle(wb, sheet_name, headerStyle, rows = 1, cols = 1:ncol(data), gridExpand = TRUE)
    setColWidths(wb, sheet_name, cols = 1:ncol(data), widths = "auto")
    
    cat("  Added sheet:", sheet_name, "\n")
  } else { 
    cat("  Achtung!  File not found:", csv_files_to_combine[[sheet_name]], "\n")
  }
}

# save
saveWorkbook(wb, 
             file.path(output_dir, "Supplementary_Data_All_Tables.xlsx"), 
             overwrite = TRUE)


# lol, need to make the jaccard into a table for the figures
jaccard_matrix <- read.csv(file.path(output_dir, "jaccard_distances_sample_types.csv"), 
                           row.names = 1)
# ugh, make it fancy
rownames(jaccard_matrix) <- c("NEWs", "Control", "Tetracycline", 
                              "Terra-Pro", "ASx lysate", "Hive Reared", "Comb")
colnames(jaccard_matrix) <- c("NEWs", "Control", "Tetracycline", 
                              "Terra-Pro", "ASx lysate", "Hive Reared", "Comb")

# round to 3
jaccard_matrix_rounded <- round(jaccard_matrix, 3)

# display - make it a lower triangle with dashes on upper
jaccard_display <- jaccard_matrix_rounded
jaccard_display[upper.tri(jaccard_display)] <- NA

write.csv(jaccard_display, 
          file.path(figures_dir, "jaccard_distances_formatted.csv"),
          row.names = TRUE, na = "—")

# check jaccard distances once more for ASx lysate
jaccard_matrix <- read.csv(file.path(output_dir, "jaccard_distances_sample_types.csv"), 
                           row.names = 1)

# extract asx lysate row
asx_distances <- jaccard_matrix["ASx lysate", ]
print(asx_distances)

# compare asx to treatments vs asx to control
cat("\nASx lysate distances:\n")
cat("  vs Control:", asx_distances$Control, "\n")
cat("  vs Tetracycline:", asx_distances$Tetracycline, "\n")
cat("  vs Terra-Pro:", asx_distances$Terra.Pro, "\n")
cat("  vs Hive Reared:", asx_distances$Hive.Reared, "\n")
cat("  vs Comb:", asx_distances$Comb, "\n")

