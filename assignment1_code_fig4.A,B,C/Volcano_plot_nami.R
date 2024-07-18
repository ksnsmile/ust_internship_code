# DEG 분석을 수행하면 두 개의 결과를 반드시 얻는다. (Fold-change와 P-value)요소이다. 
# Fold change란, 어떤 유전자에 대하여 실험군에서의 평균발현량이 대조군에서의 평균발현량의 몇 배인지를 나타내고,
# P-value는 두 군의 평균발현량 차이가 통계적으로 유의미한 값인지를 알려준다.


#패키지 준비
library(DESeq2)
library(BiocGenerics)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(SummarizedExperiment)
library(Biobase)
library(tximport)

# 데이터 준비
folder_path <- "C:/Users/sn714/OneDrive/바탕 화면/volcano" 
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

# DESeq 함수 실행
dds <- DESeq(dds)

# 결과 추출
res <- results(dds, pAdjustMethod = "BH")

# 결과 확인
head(res)

#volcano_plot 그리기

cut_lfc <- 1
cut_pvalue <- 0.05

par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(res)

# Adjusted P values

y_range <- c(0,300)
x_range <- c(-12.5,12.5)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="ERhigh vs. ERlow MPP3\nDEGs", col='grey', cex=1.0, xlab=bquote(~Log[2]~FC), ylab=bquote(~-log[10]~adjusted~p-value),xlim=x_range,ylim=y_range,xaxt='n',yaxt='n'))
axis(1, at=c(-10,0,10), labels=c("-10","0","10"))
axis(2, at=c(0, 100, 200, 300), labels=c("0","100","200", "300"))



with(subset(topT, padj<cut_pvalue & log2FoldChange>cut_lfc), points(log2FoldChange, -log10(padj), pch=20, col='red', cex=1.5))
with(subset(topT, padj<cut_pvalue & log2FoldChange<(-cut_lfc)), points(log2FoldChange, -log10(padj), pch=20, col='red', cex=1.5))

## Add lines for FC and P-value cut-off

abline(v=-cut_lfc, col='grey', lty=3, lwd=2.0)
abline(v=cut_lfc, col='grey', lty=3, lwd=2.0)
