# Gerekli paketleri yüklüyoruz
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Gerekli Bioconductor paketlerini yüklüyoruz
BiocManager::install(c("scater", "mia"))
BiocManager::install("curatedMetagenomicData")
BiocManager::install("phyloseq", force = TRUE)

#  kütüphaneler
library(dplyr)
library(DT)
library(curatedMetagenomicData)
library(phyloseq)

# metadata setini yüklüyoruz
data("sampleMetadata")


sample_data <- sampleMetadata %>%
  filter(!is.na(age)) %>%              # Yaş bilgisi eksik olanlar çıkarılır
  filter(body_site == "stool")         # Sadece dışkı örnekleri alınır

# İlk 10 satırı göster
head(sample_data, 10)

# Yaş bilgilerine göre yaş grubu oluşturuyoruz
sample_data <- sample_data %>%
  mutate(
    age_group = case_when(
      age <= 10 ~ "0-10",
      age > 10 & age <= 20 ~ "11-20",
      age > 20 & age <= 30 ~ "21-30",
      age > 30 & age <= 40 ~ "31-40",
      age > 40 & age <= 50 ~ "41-50",
      age > 50 & age <= 60 ~ "51-60",
      age >= 61 ~ "61+"
    )
  )

# Yaş grubu dağılımı
table(sample_data$age_group)

# Kullanılan çalışmalar
studies <- unique(sample_data$study_name)
print(studies)

# Sistemdeki tüm mevcut dataset isimlerini alıyoruz (dryrun TRUE sadece isim getirir)
available_datasets <- curatedMetagenomicData(".*", dryrun = TRUE)

# Sadece marker_abundance içeren datasetleri filtreliyoruz
marker_datasets <- grep("marker_abundance", available_datasets, value = TRUE)

# marker_abundance dataset isimlerinden study adlarını ayıklıyoruz
marker_studies <- sub("\\.marker_abundance$", "", marker_datasets)

# Bizim çalışmamızdaki study isimlerinden geçerli olanları buluyoruz
valid_studies <- studies[studies %in% marker_studies]

# Bu geçerli study isimleri için dataset adları oluşturuyoruz
valid_datasets <- paste0(valid_studies, ".marker_abundance")

# Geçerli datasetleri yüklüyoruz
physeq_list <- curatedMetagenomicData(
  valid_datasets,
  counts = TRUE,
  rownames = "long"
)
str(physeq_list)
physeq_list<-as.list(physeq_list)
# Listeyi birleştirip tek bir phyloseq nesnesi oluşturuyoruz
physeq <- mergeData(physeq_list)
