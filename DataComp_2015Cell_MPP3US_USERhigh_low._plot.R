library(Seurat)
library(dplyr)
library(Matrix)
library(cowplot)
setwd("/Users/hyojung/workspace/HJ_backupData/buttelab_s01/proj_cytokine/GSE72857_LK_2015Cell")
rawCount=read.csv(file="GSE72857_umitab.txt",sep="\t", header=T, row.names = 1)
#rawCount --> MARS-seq data (colum=wellID, row--> geneName=MGI id;GeneSymbol; Human ontology or Entrez_id;GeneSymbol)
rawGenes=data.frame(rowN= row.names(rawCount))

#strsplit(row.names(rawCount), split=";")
rawGenes$GeneID=unlist((lapply(row.names(rawCount), function(x){unlist(unlist(strsplit(x, split=";")))[1]})))
rawGenes$GeneSymbOld=unlist((lapply(row.names(rawCount), function(x){unlist(unlist(strsplit(x, split=";")))[2]})))

cellTypeDf=read.table(file="GSE72857_experimental_design.txt", sep="\t", header=T, skip=19)
unique(cellTypeDf$Batch_desc)
#Experimental design table
#Each row represents a single cell.
#Field:  Description:
#  Well_ID well identifier
#Seq_batch_ID    Illumina suquencing batch ID
#Amp_batch_ID    amplification batch ID (of 192 single cells)
#Mouse_ID        subject mouse ID
#Plate_ID        sorting batch ID (384 wells plates)
#Batch_description       "Information of the experiment, sorting scheme and genetic background of mouse"                                                                  
#Pool_barcode    Pool barcode (may overlap between different sorting batches)
#Cell_barcode    Single cell barcode
#RMT_sequence    random molecular tag sequence
#Number_of_cells Number of cells sorted to well (0/1)
#CD34_measurement        CD34 protein levels (if index-sorted)
#FcgR3_measurement       FcgR3 protein levels (if index-sorted)


#filter out the given set
nrow(cellTypeDf) #10368
cellTypeDf=subset(cellTypeDf, cellTypeDf$Batch_desc=='Unsorted myeloid')
nrow(cellTypeDf) #3072
cellTypeDf=subset(cellTypeDf, cellTypeDf$Number_of_cells==1)
nrow(cellTypeDf) #3040

nrow(rawGenes) #27297
rawCount=rawCount[, cellTypeDf$Well_ID]
#filter out rowSum
nrow(rawCount) #27297
rawCount=rawCount[rowSums(rawCount)>0,]
nrow(rawCount) #14260 -->52.24%

nrow(rawGenes) #27297
rawGenes=rawGenes[rawGenes$rowN%in%intersect(rawGenes$rowN, row.names(rawCount)),]
nrow(rawGenes) #14260
head(rawGenes)

library(biomaRt)
library('org.Mm.eg.db')

ensembl=useMart("ENSEMBL_MART_ENSEMBL")
ensembl_mm10_mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
coordGeneNames=getBM(attributes =c("ensembl_gene_id","mgi_symbol", "entrezgene_id","entrezgene_accession"),  mart = ensembl_mm10_mart)

length(intersect(coordGeneNames$mgi_symbol, rawGenes$GeneID)) #11538 --> 80%
length(intersect(coordGeneNames$entrezgene_id, rawGenes$GeneID)) #1
length(intersect(coordGeneNames$entrezgene_accession, rawGenes$GeneID)) #11180 -->78%
length(intersect(intersect(coordGeneNames$entrezgene_accession, rawGenes$GeneID), intersect(coordGeneNames$mgi_symbol, rawGenes$GeneID))) #11171
11538+11180-11171 #11547 80.9%

#load myData
load("/Users/hyojung/workspace/HJ_backupData/buttelab_s01/proj_cytokine/scRNA_3rdMPP3ERhighlow_STUS/MPP3US:MPP3ST:MPP3USERL:MPP3STERL:MPP3USERH:MPP3STERH_SeuratObjList.RData") #my_sampleCtrl.list
my_sampleCtrl.list
selectedSet=c("MPP3US", "MPP3USERL","MPP3USERH")
my_sampleCtrl.list=my_sampleCtrl.list[selectedSet]
my_sampleCtrl.list$MPP3US@assays$RNA[1:10,1:10] #rowNames are geneSymbol --> Let's convert rawCount rowNames as geneSymbol
length(intersect(coordGeneNames$mgi_symbol, row.names(my_sampleCtrl.list$MPP3US@assays$RNA))) #14865 of 15347 it's 96%...
length(intersect(coordGeneNames$entrezgene_accession,row.names(my_sampleCtrl.list$MPP3US@assays$RNA)))
length(intersect(intersect(coordGeneNames$entrezgene_accession,row.names(my_sampleCtrl.list$MPP3US@assays$RNA)), intersect(coordGeneNames$mgi_symbol, row.names(my_sampleCtrl.list$MPP3US@assays$RNA)))) #13174
length(setdiff(intersect(coordGeneNames$mgi_symbol, row.names(my_sampleCtrl.list$MPP3US@assays$RNA)), intersect(coordGeneNames$entrezgene_accession,row.names(my_sampleCtrl.list$MPP3US@assays$RNA))))


rawCount2=rawCount
rownames(rawCount2)=rawGenes[match(row.names(rawCount2), rawGenes$rowN),]$GeneID #update by GeneID
length(intersect(rownames(rawCount2),row.names(my_sampleCtrl.list$MPP3US@assays$RNA) )) #10678 -->74% of rawCount2
length(intersect(rownames(rawCount2),row.names(my_sampleCtrl.list$MPP3USERL@assays$RNA) )) #10780 --> 75%
length(intersect(rownames(rawCount2),row.names(my_sampleCtrl.list$MPP3USERH@assays$RNA) )) #10692 -->74%
coordGeneNames[match(rawGenes$GeneID, coordGeneNames$entrezgene_accession),]
subset(coordGeneNames, coordGeneNames$mgi_symbol!=coordGeneNames$entrezgene_accession)
#rownames(rawCount2[rownames(rawCount2)%in%rawGenes$rowN,])=
#  coordGeneNames[match(rawGenes[rownames(rawCount2)%in%rawGenes$rowN,]$GeneID,coordGeneNames$entrezgene_accession),]$entrezgene_accession

PublicD <- CreateSeuratObject(counts = rawCount2, min.cells = 3, min.features = 100, project = "Cell2015Ref")
PublicD@meta.data #orig.ident, nCount_RNA, nFeature_RNA, rowname ==> Wall ID (Cell ID)
#for mouse (^mt-), for human (^MT-)
PublicD[["percent.mt"]] <- PercentageFeatureSet(PublicD, pattern = "^mt-")
VlnPlot(PublicD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # mitochonrial gene expression null??

#test data
table(PublicD@meta.data$orig.ident)
PublicD.norm <- NormalizeData(PublicD , normalization.method = "LogNormalize", scale.factor = 10000)
# Identify the 10 most highly variable genes
PublicD.norm <- FindVariableFeatures(PublicD.norm, selection.method = "vst", nfeatures = 2000) 
top10PublicD <- head(VariableFeatures(PublicD.norm), 10)
top10PublicD #[1] "Cd74"    "H2-Aa"   "Prss34"  "Prg2"    "Pf4"     "Ccl5"    "Epx"     "Atpase6" "Hba-a2"  "H2-Eb1" 


#add Meta data --> index sort value of FACS?
dim(PublicD) #[1] 12485  2930 (row--> genes,col--> cell count)
head(cellTypeDf)

cellTypeDf=cellTypeDf[rownames(PublicD@meta.data)%in%cellTypeDf$Well_ID,]
dim(cellTypeDf) #2930, 13
tmp_d<-gsub(' ', '', cellTypeDf$Batch_desc)
names(tmp_d)<-cellTypeDf$Well_ID

tmp_d1=cellTypeDf$CD34_measurement
names(tmp_d1)=cellTypeDf$Well_ID

tmp_d2=cellTypeDf$FcgR3_measurement
names(tmp_d2)=cellTypeDf$Well_ID

PublicD=AddMetaData(object=PublicD, metadata=tmp_d, col.name="group")
PublicD=AddMetaData(object=PublicD, metadata=tmp_d1, col.name="CD34_measurement")
PublicD=AddMetaData(object=PublicD, metadata=tmp_d2, col.name="FcgR3_measurement")



#Add PublicD to my Set & Normalize it again...
my_sampleCtrl.list$PublicD=PublicD

#addMeta
#a=rep("10x", dim(my_sampleCtrl.list[[1]])[2]+dim(my_sampleCtrl.list[[2]])[2]+dim(my_sampleCtrl.list[[3]])[2])
#a=append(a, rep("MARSseq", dim(my_sampleCtrl.list[[4]])[2]))
#my_sampleCtrl.list=AddMetaData(my_sampleCtrl.list, metadata=a, col.name="Method")

#a2=rep("Cultured", dim(my_sampleCtrl.list[[1]])[2]+dim(my_sampleCtrl.list[[2]])[2]+dim(my_sampleCtrl.list[[3]])[2])
#a2=append(a2, rep("Fresh", dim(my_sampleCtrl.list[[4]])[2]))
#my_sampleCtrl.list=AddMetaData(my_sampleCtrl.list, metadata=a2, col.name="Cond")

#make normalize & point default assay as 'integrated'
for (i in 1:length(my_sampleCtrl.list)) {
  my_sampleCtrl.list[[i]] <- NormalizeData(my_sampleCtrl.list[[i]], verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 10000)
  my_sampleCtrl.list[[i]] <- FindVariableFeatures(my_sampleCtrl.list[[i]], selection.method = "vst", 
                                                  nfeatures = 2000, verbose = FALSE)
}

for(i in 1:3){
  print(i)
  a=rep("10x", dim(my_sampleCtrl.list[[i]])[2])
  my_sampleCtrl.list[[i]]=AddMetaData(my_sampleCtrl.list[[i]], metadata=a, col.name="Method")
  a2=rep("Cultured", dim(my_sampleCtrl.list[[i]])[2])
  my_sampleCtrl.list[[i]]=AddMetaData(my_sampleCtrl.list[[i]], metadata=a2, col.name="Cond")
  #my_sampleCtrl.list=AddMetaData
}

i=4
a=rep("MARSseq", dim(my_sampleCtrl.list[[i]])[2])
my_sampleCtrl.list[[i]]=AddMetaData(my_sampleCtrl.list[[i]], metadata=a, col.name="Method")
a2=rep("Fresh", dim(my_sampleCtrl.list[[i]])[2])
my_sampleCtrl.list[[i]]=AddMetaData(my_sampleCtrl.list[[i]], metadata=a2, col.name="Cond")


my_sampleCtrl.anchors <- FindIntegrationAnchors(object.list = my_sampleCtrl.list, dims = 1:30)
my_sampleCtrl.integrated <- IntegrateData(anchorset = my_sampleCtrl.anchors, dims = 1:30)
DefaultAssay(my_sampleCtrl.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
#my_sampleCtrl.integrated <- ScaleData(my_sampleCtrl.integrated, verbose = FALSE,  vars.to.regress = c('Method', 'group','Cond'))
my_sampleCtrl.integrated <- ScaleData(my_sampleCtrl.integrated, verbose = FALSE,  vars.to.regress = c('Method', 'Cond'))
my_sampleCtrl.integrated <- RunPCA(my_sampleCtrl.integrated, npcs = 30, verbose = FALSE)
# t-SNE, UMAP  and Clustering
my_sampleCtrl.integrated<- RunUMAP(my_sampleCtrl.integrated, reduction = "pca", dims = 1:20)
#Constructs a Shared Nearest Neighbor (SNN) Graph for a given dataset.
my_sampleCtrl.integrated <- FindNeighbors(my_sampleCtrl.integrated, reduction = "pca", dims = 1:20) 
my_sampleCtrl.integrated <- FindClusters(my_sampleCtrl.integrated, resolution = 0.5)


# Visualization
p1 <- DimPlot(my_sampleCtrl.integrated, reduction = "umap", group.by = "group")
p1_1 <- DimPlot(my_sampleCtrl.integrated, reduction = "umap", group.by = "group", split.by = "group")
p2 <- DimPlot(my_sampleCtrl.integrated, reduction = "umap", label = TRUE)
p3 <- DimPlot(my_sampleCtrl.integrated, reduction = "umap", group.by = "Method")
p4 <- DimPlot(my_sampleCtrl.integrated, reduction = "umap", group.by = "Cond")
plot_grid(p1, p1_1, p2, p3, p4, ncol=2)

my_sampleCtrl.integrated
mylist=unique(my_sampleCtrl.integrated@meta.data$group) #"MPP3US"          "MPP3USERL"       "MPP3USERH"       "Unsortedmyeloid"

p1 <- DimPlot(my_sampleCtrl.integrated, reduction = "umap", group.by = "group")
p1_1 <- DimPlot(my_sampleCtrl.integrated, reduction = "umap", group.by = "group", split.by = "group")
p2 <- DimPlot(my_sampleCtrl.integrated, reduction = "umap", label = TRUE)
p3 <- DimPlot(my_sampleCtrl.integrated, reduction = "umap", group.by = "Method")
p4 <- DimPlot(my_sampleCtrl.integrated, reduction = "umap", group.by = "Cond")
plot_grid(p1, p1_1, p2, p3, p4, ncol=2)
my_tags=paste(mylist, collapse = ":")
save(file=paste(my_tags,"_integratedUMAP.RData", sep=""), my_sampleCtrl.integrated)

my_sampleCtrl.integratedCopy=my_sampleCtrl.integrated
#specific feature plot for pulicD data
my_sampleCtrl.integrated@meta.data[is.na(my_sampleCtrl.integrated@meta.data$CD34_measurement),]$CD34_measurement=0
my_sampleCtrl.integrated@meta.data[(my_sampleCtrl.integrated@meta.data$CD34_measurement)<0,]$CD34_measurement=0
my_sampleCtrl.integrated@meta.data[is.na(my_sampleCtrl.integrated@meta.data$FcgR3_measurement),]$FcgR3_measurement=0
my_sampleCtrl.integrated@meta.data[(my_sampleCtrl.integrated@meta.data$FcgR3_measurement)<0,]$FcgR3_measurement=0

p1 <- DimPlot(my_sampleCtrl.integrated, reduction = "umap", group.by = "group")
p1_1 <- DimPlot(my_sampleCtrl.integrated, reduction = "umap", group.by = "group", split.by = "group")
p2 <- DimPlot(my_sampleCtrl.integrated, reduction = "umap", label = TRUE)
p3 <- DimPlot(my_sampleCtrl.integrated, reduction = "umap", group.by = "Method")
p4 <- DimPlot(my_sampleCtrl.integrated, reduction = "umap", group.by = "Cond")
plot_grid(p1, p1_1, p2, p3, p4, ncol=2)

library(gridExtra)

library(ggplot2)
DimplotList=list()
#for DimPlot version
for(i in seq(length(mylist))){
  myCells=subset(x=my_sampleCtrl.integrated, subset = group%in%mylist[i])
  p1=DimPlot(my_sampleCtrl.integrated, label=FALSE, group.by = "group", cells.highlight = list(rownames(myCells@meta.data)), cols.highlight = c("darkblue"), sizes.highlight = 0.5, pt.size=0.5)+ ggplot2::ggtitle(mylist[i])
  DimplotList[[i]]=p1
}
grid.arrange(DimplotList[[1]],DimplotList[[2]],DimplotList[[3]], DimplotList[[4]],ncol = 2)
ggsave(paste(my_tags,"_ByGroup.pdf", sep=""),  grid.arrange(DimplotList[[1]],DimplotList[[2]],DimplotList[[3]], DimplotList[[4]],ncol = 2), device='pdf', scale=1, width=15, height=12, units="in")



FeatplotList=list()

ii=4
myCells=subset(x=my_sampleCtrl.integrated, subset = group%in%mylist[ii])
p2=FeaturePlot(my_sampleCtrl.integrated, cells=(rownames(myCells@meta.data)), features="CD34_measurement",cols=c("light gray", "blue", "red"),  min.cutoff = 'q10', max.cutoff = 'q90')
FeatplotList[[1]]=p2
p2=FeaturePlot(my_sampleCtrl.integrated, features="CD34_measurement",cols=c("light gray", "blue", "red"),   min.cutoff = 'q10', max.cutoff = 'q90')
FeatplotList[[2]]=p2
p2=FeaturePlot(my_sampleCtrl.integrated, cells=(rownames(myCells@meta.data)), features="FcgR3_measurement",cols=c("light gray", "blue", "red"),   min.cutoff = 'q10', max.cutoff = 'q90')
FeatplotList[[3]]=p2
p2=FeaturePlot(my_sampleCtrl.integrated, features="FcgR3_measurement",cols=c("light gray", "blue", "red"),  min.cutoff = 'q10', max.cutoff = 'q90')
FeatplotList[[4]]=p2
grid.arrange(FeatplotList[[1]],FeatplotList[[2]],FeatplotList[[3]], FeatplotList[[4]], ncol = 2)
ggsave(paste(mylist[ii],"_FACSFeaturPlot2.pdf", sep=""),  grid.arrange(FeatplotList[[1]],FeatplotList[[2]],FeatplotList[[3]], FeatplotList[[4]], ncol = 2), device='pdf', scale=1, width=17, height=10, units="in")



my_sampleCtrl.integrated@meta.data$CD34_measurement=log(my_sampleCtrl.integrated@meta.data$CD34_measurement+1, base=10)
my_sampleCtrl.integrated@meta.data$FcgR3_measurement=log(my_sampleCtrl.integrated@meta.data$FcgR3_measurement+1, base=10)
myCells=subset(x=my_sampleCtrl.integrated, subset = group%in%mylist[ii])
p2=FeaturePlot(my_sampleCtrl.integrated, cells=(rownames(myCells@meta.data)), features="CD34_measurement",cols=c("light gray", "blue", "red"))
FeatplotList[[5]]=p2
p2=FeaturePlot(my_sampleCtrl.integrated, features="CD34_measurement",cols=c("light gray", "blue", "red"))
FeatplotList[[6]]=p2
p2=FeaturePlot(my_sampleCtrl.integrated, cells=(rownames(myCells@meta.data)), features="FcgR3_measurement",cols=c("light gray", "blue", "red"))
FeatplotList[[7]]=p2
p2=FeaturePlot(my_sampleCtrl.integrated, features="FcgR3_measurement",cols=c("light gray", "blue", "red"))
FeatplotList[[8]]=p2
grid.arrange(FeatplotList[[5]],FeatplotList[[6]],FeatplotList[[7]], FeatplotList[[8]], ncol = 2)
ggsave(paste(mylist[ii],"_FACSFeaturPlot_log10Scale.pdf", sep=""),  grid.arrange(FeatplotList[[5]],FeatplotList[[6]],FeatplotList[[7]], FeatplotList[[8]], ncol = 2), device='pdf', scale=1, width=15, height=12, units="in")
