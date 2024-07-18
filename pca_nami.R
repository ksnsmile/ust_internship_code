
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
library(BiocManager)
library(limma)
library(ggplot2)
library(ggrepel) 

# 데이터 준비
folder_path <- "C:/Users/sn714/OneDrive/바탕 화면/bulk" 
files <- list.files(folder_path, pattern = "quant.genes.sf", full.names = TRUE)
sample_names <- gsub("_quant.genes.sf", "", basename(files))
sample_names
sample_info <- data.frame(sampleName = sample_names, condition = gsub("_[0-9]+$", "", sample_names))
sample_info
rownames(sample_info) <- sample_info$sampleName
txi <- tximport(files, type = "salmon", txOut = TRUE)


# DESeqDataSet 생성 및 필터링
dds <- DESeqDataSetFromTximport(txi, colData = sample_info, design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# VST 변환
vsd <- vst(dds, blind = FALSE)

#배치효과 제거
vsd_matrix <- assay(vsd)
batch <- sample_info$batch
design <- model.matrix(~ condition, data=sample_info)

# batch 제거 함수 사용
vsd_batch_corrected <- removeBatchEffect(vsd_matrix,batch=batch, design=design)

# vsd 객체에 배치 효과 제거된 데이터 다시 할당
assay(vsd) <- vsd_batch_corrected

# PCA 플롯 생성
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# 그룹별 중심 좌표 계산
centers <- aggregate(cbind(PC1, PC2) ~ condition, data = pca_data, FUN = mean)

# 그룹별 원을 그리기 위한 함수
create_circle <- function(center, radius, npoints = 100) {
  angles <- seq(0, 2 * pi, length.out = npoints)
  data.frame(
    PC1 = center[1] + radius * cos(angles),
    PC2 = center[2] + radius * sin(angles)
  )
}

# 그룹별 원 데이터 생성
radius <- 10  # 원하는 반지름 값
circles <- do.call(rbind, lapply(1:nrow(centers), function(i) {
  cbind(create_circle(centers[i, 2:3], radius), condition = centers[i, 1])
}))



ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  geom_polygon(data = circles, aes(x = PC1, y = PC2, group = condition), fill = NA, color = "black", linetype = "dotted") +
  geom_text(data = centers, aes(x = PC1, y = PC2, label = condition), size = 5, vjust = -1) +
  labs(title = "Bulk RNA-seq",
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "lightgrey", color = NA)
  ) +
  xlim(-50, 75) +  # x축 범위 설정 (예시)
  ylim(-40,40) 
















