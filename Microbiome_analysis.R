# Install BiocManager if not installed; BiocManager manages Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install curatedMetagenomicData (Bioconductor package)
if (!requireNamespace("curatedMetagenomicData", quietly = TRUE)) {
  BiocManager::install("curatedMetagenomicData", ask = FALSE)
}
# Install other required packages from CRAN
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("vegan", quietly = TRUE)) {
  install.packages("vegan")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

# Load the libraries
library(curatedMetagenomicData)  # Curated microbiome datasets
library(dplyr)                  # Data manipulation (filter, mutate, etc.)
library(ggplot2)                # Plotting library for graphs
library(vegan)                  # Ecological analysis (diversity measures)
library(pheatmap)               # Heatmap generation

# Set a random seed for reproducibility of any random processes
set.seed(123)


# Retrieve the sample metadata data frame provided by curatedMetagenomicData
metadata <- curatedMetagenomicData::sampleMetadata

# Filter to keep only stool samples (gut microbiome) and samples with known age
metadata_filtered <- metadata %>%
  filter(body_site == "stool") %>%    # select only stool samples
  filter(!is.na(age))                # exclude samples without age information

# Create an 'AgeGroup' variable with specified age brackets using case_when
metadata_filtered <- metadata_filtered %>%
  mutate(AgeGroup = case_when(
    age <= 10              ~ "0-10",    # ages 0-10
    age > 10  & age <= 20  ~ "11-20",   # ages 11-20
    age > 20  & age <= 30  ~ "21-30",   # ages 21-30
    age > 30  & age <= 40  ~ "31-40",   # ages 31-40
    age > 40  & age <= 50  ~ "41-50",   # ages 41-50
    age > 50  & age <= 60  ~ "51-60",   # ages 51-60
    age >= 61             ~ "61+",      # ages 61 and above
    TRUE                   ~ NA_character_
  ))
# Remove any columns that are entirely NA to clean up the metadata
metadata_filtered <- metadata_filtered %>%
  select(where(~ !all(is.na(.x))))

# Use returnSamples() to fetch the taxonomic abundance data for these samples
# We want the 'relative_abundance' data type; rownames = "short" uses short taxonomy labels.
se_data <- returnSamples(metadata_filtered, dataType = "relative_abundance", rownames = "short")

# 'se_data' is a (Tree)SummarizedExperiment object with assays and metadata
# Extract the abundance matrix (rows: species, columns: samples) from the assay
otu_matrix <- assay(se_data, "relative_abundance")  # matrix of relative abundances

# Convert the sample metadata (colData of se_data) into a data frame
sample_metadata <- as.data.frame(colData(se_data))

# Ensure 'AgeGroup' is a factor with a specified order of levels for consistent plotting
age_levels <- c("0-10","11-20","21-30","31-40","41-50","51-60","61+")
sample_metadata$AgeGroup <- factor(sample_metadata$AgeGroup, levels = age_levels)


# Calculate richness (number of species) and Shannon diversity index for each sample.
# Richness is the count of species with non-zero abundance in a sample.
richness <- colSums(otu_matrix > 0)  # sum species with abundance >0 in each sample 

# Shannon diversity (use vegan::diversity, which computes the Shannon index by default)
shannon <- apply(otu_matrix, 2, function(x) {
  diversity(x, index = "shannon")
})
# Add these diversity measures to the sample metadata
sample_metadata$Richness <- richness
sample_metadata$Shannon <- shannon


# Count the number of distinct species present in each age group.
# A taxon is considered present in an age group if it has non-zero abundance in any sample of that group.
species_count_by_age <- sapply(levels(sample_metadata$AgeGroup), function(g) {
  # Identify which samples belong to this age group
  samples_in_group <- colnames(otu_matrix)[sample_metadata$AgeGroup == g]
  # Count species (rows) that have any non-zero entries in those samples
  sum(rowSums(otu_matrix[, samples_in_group, drop = FALSE] > 0) > 0)
})
# Put the results in a data frame for plotting
species_count_df <- data.frame(AgeGroup = levels(sample_metadata$AgeGroup),
                            speciesCount = as.numeric(species_count_by_age))


# Perform PCA on the taxonomic abundance data (transpose to have samples as rows).
# Log-transform the abundances (with a small pseudo-count) to stabilize variance.
otu_log <- log(otu_matrix + 1e-6)  # add a small constant to avoid log(0)
pca_res <- prcomp(t(otu_log), center = TRUE, scale. = FALSE)

# Extract PCA scores (principal component coordinates) for the first two components
pca_scores <- data.frame(PC1 = pca_res$x[, 1],
                         PC2 = pca_res$x[, 2],
                         AgeGroup = sample_metadata$AgeGroup)

# 1. Bar plot: Distinct Species count by age group
ggplot(species_count_df, aes(x = AgeGroup, y = speciesCount, fill = AgeGroup)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Number of Distinct Species by Age Group",
       x = "Age Group", y = "species Count") +
  theme_minimal() +
  theme(legend.position = "none")

# 2. Boxplots of alpha diversity (Richness and Shannon) by age group
# Richness by Age Group
ggplot(sample_metadata, aes(x = AgeGroup, y = Richness, fill = AgeGroup)) +
  geom_boxplot() +
  labs(title = "Richness (Species Count) by Age Group",
       x = "Age Group", y = "Richness") +
  theme_minimal() +
  theme(legend.position = "none")

# Shannon diversity by Age Group
ggplot(sample_metadata, aes(x = AgeGroup, y = Shannon, fill = AgeGroup)) +
  geom_boxplot() +
  labs(title = "Shannon Diversity by Age Group",
       x = "Age Group", y = "Shannon Index") +
  theme_minimal() +
  theme(legend.position = "none")

# 3. PCA scatter plot (PC1 vs PC2) colored by age group
ggplot(pca_scores, aes(x = PC1, y = PC2, color = AgeGroup)) +
  geom_point(alpha = 0.7) +
  labs(title = "PCA of Microbiome Samples by Age Group",
       x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()

# 4. Heatmap of species abundance (top 20 species across all samples)
# Identify top 20 species by mean abundance
species_mean_abundance <- rowMeans(otu_matrix)
top20_species <- names(sort(species_mean_abundance, decreasing = TRUE))[1:20]
otu_top20 <- log10(otu_matrix[top20_species, , drop = FALSE] + 1e-5)


# Prepare annotation for age group (columns) in heatmap
annotation_col <- data.frame(AgeGroup = sample_metadata$AgeGroup)
rownames(annotation_col) <- colnames(otu_top20)

# Draw the heatmap (rows = species, columns = samples), scaling rows to emphasize differences
pheatmap(otu_top20,
         scale = "row",               # normalize each taxon (row) to zero mean and unit variance
         annotation_col = annotation_col,
         show_colnames = FALSE,
         main = "Heatmap of Top 20 species Abundance")


# Extract top 20 species by average abundance
species_mean_abundance <- rowMeans(otu_matrix)
top20_species <- names(sort(species_mean_abundance, decreasing = TRUE))[1:20]
otu_top20 <- otu_matrix[top20_species, , drop = FALSE]

# Normalize each row (taxon) to range [-1, 1]
normalize_to_minus1_1 <- function(x) {
  rng <- max(x) - min(x)
  if (rng == 0) return(rep(0, length(x)))  # avoid division by zero
  return(2 * (x - min(x)) / rng - 1)
}

otu_top20_scaled <- t(apply(otu_top20, 1, normalize_to_minus1_1))

# Prepare annotation for age group
annotation_col <- data.frame(AgeGroup = sample_metadata$AgeGroup)
rownames(annotation_col) <- colnames(otu_top20_scaled)

# Create the heatmap with normalized data
pheatmap(otu_top20_scaled,
         scale = "none",               # Already scaled manually
         annotation_col = annotation_col,
         show_colnames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Top 20 species (Normalized to [-1, 1])")

alpha_val <- 0.01
#  one-way ANOVA
anova_result <- aov(Shannon ~ AgeGroup, data = sample_metadata)
anova_summary <- summary(anova_result)

# Print the summary
print(anova_summary)

# Check if the p-value is significant
p_value <- anova_summary[[1]][["Pr(>F)"]][1]
if (p_value < alpha_val) {
  cat("The result is significant (p < 0.01).\n")
} else {
  cat("The result is not significant (p >= 0.01).\n")
}
