# Analysis of Gut Microbiome Diversity by Age Group

## Description

This project analyzes how gut microbiome diversity changes across different age groups using data from the `curatedMetagenomicData` package available on Bioconductor. Only stool (fecal) samples were included in the analysis.

## Objectives

* To analyze how gut microbiome diversity varies with age
* To compare age groups using alpha diversity metrics (Richness and Shannon Index)
* To visualize sample distribution by age group using PCA
* To illustrate the distribution of the top 20 most abundant microbial species by age group with a heatmap

## Technologies Used

* **Programming Language:** R
* **Development Environment:** RStudio (2025.05.0+496)
* **R Version:** 4.5.0 (2025-04-11)
* **Libraries:**

  * `curatedMetagenomicData` - Microbiome dataset
  * `dplyr` - Data manipulation
  * `ggplot2` - Plotting
  * `vegan` - Alpha diversity calculation
  * `pheatmap` - Heatmap generation

## Outputs

* Bar plot of species count (Richness) by age group
* Box plots of alpha diversity (Richness and Shannon Index)
* PCA scatter plot by age group
* Heatmap of the top 20 most abundant species
* ANOVA test results for statistical significance

## Execution Instructions

1. Make sure R or RStudio is installed.
2. Install required R packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("curatedMetagenomicData")
install.packages(c("dplyr", "ggplot2", "vegan", "pheatmap"))
```

3. Run the code in `analysis.R` to perform the analysis.

## Notes

* Only stool samples with available age information were used.
* Bray-Curtis analysis was not included in this version.

## License

This project was developed for educational purposes and is released as open source.

---


