# Gerekli paketleri yüklüyoruz
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Gerekli Bioconductor paketlerini yüklüyoruz
BiocManager::install(c("scater", "mia"))
BiocManager::install("curatedMetagenomicData")
BiocManager::install("phyloseq", force = TRUE)

# Kütüphaneleri yüklüyoruz
library(dplyr)
library(DT)
library(curatedMetagenomicData)
library(phyloseq)

# Metadata setini yüklüyoruz
data("sampleMetadata")

# Yaş bilgisi olan ve dışkı örneği olanları filtreliyoruz
sample_data <- sampleMetadata %>%
  filter(!is.na(age)) %>%
  filter(body_site == "stool")

# Kullanılan çalışmalar
studies <- unique(sample_data$study_name)
print(studies)

# Sistemdeki tüm mevcut dataset isimlerini alıyoruz (dryrun TRUE sadece isim getirir)
available_datasets <- curatedMetagenomicData(".*", dryrun = TRUE)

# Sadece gene_families veya pathway_abundance içeren datasetleri filtreliyoruz
marker_datasets <- grep("gene_families|pathway_abundance", available_datasets, value = TRUE)

# marker_abundance dataset isimlerinden study adlarını ayıklıyoruz
marker_studies <- sub("\\.(gene_families|pathway_abundance)$", "", marker_datasets)

# Bizim çalışmamızdaki study isimlerinden geçerli olanları buluyoruz
valid_studies <- studies[studies %in% marker_studies]

# Bu geçerli study isimleri için dataset adları oluşturuyoruz
valid_datasets <- paste0(valid_studies, ".gene_families")  # veya ".pathway_abundance"

# Geçerli datasetleri yüklüyoruz ve phyloseq nesneleri oluşturuyoruz
physeq_list <- lapply(valid_datasets, function(dataset) {
  data <- curatedMetagenomicData(dataset, counts = TRUE, rownames = "long")
  
  # Assay verisini kontrol et
  if (!is.null(assay(data))) {
    assay_data <- assay(data)
    
    if (nrow(assay_data) > 0) {
      # Gene families veya pathway abundance kullanarak phyloseq nesnesi oluştur
      # Burada doğrudan assay_data'yı kullanıyoruz
      sample_data_info <- sample_data(data)
      return(phyloseq(otu_table(assay_data, taxa_are_rows = TRUE), sample_data_info))
    } else {
      warning(paste("Data is empty for dataset:", dataset))
      return(NULL)
    }
  } else {
    warning(paste("No assay data for dataset:", dataset))
    return(NULL)
  }
})

# NULL girişlerini listeden çıkar
physeq_list <- physeq_list[!sapply(physeq_list, is.null)]

# Geçerli phyloseq nesnelerinin sayısını kontrol et
cat("Number of valid phyloseq objects:", length(physeq_list), "\n")

# Eğer geçerli nesneler varsa birleştir
if (length(physeq_list) > 0) {
  physeq <- do.call(merge_phyloseq, physeq_list)
} else {
  stop("No valid phyloseq objects to merge.")
}

# Oluşan physeq nesnesinin yapısını kontrol et
str(physeq)