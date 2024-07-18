# 필요한 패키지 설치 및 로드
library(DESeq2)
library(BiocGenerics)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(SummarizedExperiment)
library(Biobase)
library(tximport)
library(ggplot2)
library(cluster)

# 데이터 준비

folder_path <- "/home/aeri/bulk" 
files <- list.files(folder_path, pattern = "quant.genes.sf", full.names = TRUE)
sample_names <- gsub("_quant.genes.sf", "", basename(files))
sample_names
sample_info <- data.frame(sampleName = sample_names, condition = gsub("_[0-9]+$", "", sample_names))
sample_info
rownames(sample_info) <- sample_info$sampleName
txi <- tximport(files, type = "salmon", txOut = TRUE)

# DESeqDataSet 생성 및 필터링
dds <- DESeqDataSetFromTximport(txi, colData = sample_info, design = ~ condition)

# DESeq2 정규화 수행
dds <- DESeq(dds)

# 정규화된 카운트 데이터 추출 (사이즈 팩터 정규화)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts

epsilon <- .Machine$double.eps  # 아주 작은 값 (머신의 최소 양수값)
cv <- apply(normalized_counts, 1, function(x) {
  m <- mean(x)
  return(sd(x) / (m + epsilon))
})


# 무작위 퍼뮤테이션을 통한 배경 분포 생성
set.seed(123)
num_permutations <- 10 # 퍼뮤테이션 수를 설정
perm_cv <- replicate(num_permutations, {
  perm_counts <- normalized_counts
  for (i in 1:nrow(perm_counts)) {
    perm_counts[i, ] <- sample(perm_counts[i, ])
  }
  apply(perm_counts, 1, function(x) {
    m <- mean(x)
    return(sd(x) / (m + epsilon))
  })
})

# 배경 분포를 벡터 형태로 변환
perm_cv_flatten <- as.vector(perm_cv)

# 데이터 프레임 생성
cv_data <- data.frame(
  cv = perm_cv_flatten,
  type = rep("Permutation", length(perm_cv_flatten))
)

# 히스토그램 시각화
ggplot(cv_data, aes(x = cv, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  labs(title = "Distribution of Coefficient of Variation (Permutation)", x = "Coefficient of Variation", y = "Frequency") +
  theme_minimal()

















# 히스토그램 시각화
ggplot(cv_data, aes(x = cv, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  labs(title = "Distribution of Coefficient of Variation (Permutation)", x = "Coefficient of Variation", y = "Frequency") +
  theme_minimal()





# K-means 클러스터링 수행
fviz_nbclust(high_var_genes, kmeans, method = "wss")
set.seed(123)
k <- 5
kmeans_result <- kmeans(high_var_genes, centers = k)

# 클러스터 할당 결과 추가
cluster_assignment <- kmeans_result$cluster

# 클러스터링 결과 시각화
fviz_cluster(kmeans_result, data = high_var_genes, geom = "point", stand = FALSE)


