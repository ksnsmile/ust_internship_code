

library(tximport)
library(biomaRt)
library(dplyr)

# Salmon 결과 파일 경로 지정
files <- list.files(path = "C:/Users/sn714/OneDrive/바탕 화면/volcano", pattern = "quant.genes.sf", full.names = TRUE)

# tximport를 사용하여 데이터 읽기
txi <- tximport(files, type = "salmon", txOut = TRUE)
txi
# 데이터 확인
head(txi$counts)


# biomaRt를 사용하여 Ensembl에서 유전자 정보를 가져오기
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Ensembl Gene ID를 일반적인 유전자 이름으로 변환
genes <- rownames(txi$counts)
gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = genes,
                   mart = ensembl)

# 유전자 이름 추가
txi$counts <- as.data.frame(txi$counts)
txi$counts$ensembl_gene_id <- rownames(txi$counts)
#colnames(txi$counts) <-c("ensembl_gene_id", paste0("V", 1:14))
colnames(gene_info)
txi_merged <- merge(x = txi$counts, y = gene_info, by = "ensembl_gene_id", all.x = TRUE)
head(txi)

str(txi$counts)
str(gene_info)

head(txi$counts)
head(gene_info)

colnames(txi$counts)
colnames(gene_info)
# 결과 확인
head(txi_merged)
# 필요한 열만 선택
txi_merged_selected <- txi_merged[, c("external_gene_name", "ER_high_1", "ER_high_2", "ER_high_3", "ER_low_1", "ER_low_2", "ER_low_3 ")]

# CSV 파일로 저장
write.csv(txi_merged_selected, file = "txi_merged_selected2.csv", row.names = FALSE)

# 저장된 파일 확인
# read.csv("txi_merged_selected.csv")  # 파일을 읽어서 확인할 수 있음
read.csv("txi_merged_selected2.csv")


