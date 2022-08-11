#######################################################################################
##############         Preprocessing snRNA from Lau et al,                  ###########
#######################################################################################
library(DoubletFinder)
setwd("/Projects/deng/Alzheimer/Lau/data/Sample/") #samples were downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157827
fileList <- dir()
list.Object <- lapply(fileList, FUN = function(x){
  data <- Read10X(data.dir = paste0(x,"/",sep=""))
  project <- CreateSeuratObject(counts = data, project = x, min.cells = 3, min.features = 200)
  project[["percent.mt"]] <- PercentageFeatureSet(project, pattern = "^MT-")
  project <- subset(project, subset = nFeature_RNA > 200  & percent.mt < 20) #As the high energy consumption and the global distribution of mitochodrial in neurons, we set the percent.mt as 20%. 
  project <- NormalizeData(project, normalization.method = "LogNormalize", scale.factor = 10000)
  project <- FindVariableFeatures(project, selection.method = "vst", nfeatures = 2000)
  project <- ScaleData(project, verbose = FALSE)
  project <- RunPCA(project, features = VariableFeatures(object = project))
  project<- RunUMAP(project, reduction = "pca", dims = 1:20)
  nExp <- round(ncol(project) * 0.05)  # expect 5% doublets
  project <- doubletFinder_v3(project, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:20)
  DF.name = colnames(project@meta.data)[grepl("DF.classification", colnames(project@meta.data))]
  project = project[, project@meta.data[, DF.name] == "Singlet"]
}
)

LauAD<- merge(list.Object[[1]], y = c(list.Object[-1], project = "LauADALL")
VlnPlot(LauAD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0, ncol = 3)
LauAD=DietSeurat(LauAD, counts = TRUE, data = TRUE,scale.data = FALSE, features = NULL,assays = NULL,dimreducs = NULL,graphs = NULL)
#LauAD <- subset(LauAD, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mt < 20)
LauAD <- NormalizeData(LauAD, normalization.method = "LogNormalize", scale.factor = 10000)
LauAD <- FindVariableFeatures(LauAD, selection.method = "vst", nfeatures = 2000)
LauAD <- ScaleData(LauAD, verbose = FALSE)
LauAD <- RunPCA(LauAD, features = VariableFeatures(object = LauAD))
#ElbowPlot(LauAD,ndims = 50) #confirmed the number of PCs
# t-SNE and Clustering
LauAD<- RunUMAP(LauAD, reduction = "pca", dims = 1:20)
LauAD<- FindNeighbors(LauAD, reduction = "pca", dims = 1:20)
LauAD<- FindClusters(LauAD, resolution = 0.5)
#DimPlot(LauAD,label=TRUE)+NoLegend()
#FeaturePlot(LauAD,c("NRGN","CAMK2A","GAD1","MBP","GFAP","VCAN","CD74","FLT1","PDGFRA"))
#FeaturePlot(LauAD,features=c("CLDN5","FLT1","ABCB1","EBF1"))
saveRDS(LauAD8Author, file="../mergeWithSampleBatchExisted.Rds")

#--------------divide the whole dataset into control and alzheimer group---------------------
setwd("/Projects/deng/Alzheimer/syn18485175/cellCycle/LauAllData")
LauAD=readRDS("/Projects/deng/Alzheimer/Lau/data/mergeWithSampleExisted.Rds")
SampleInfo=read.table("/Projects/deng/Alzheimer/Lau/sampleInfo.txt",header=T,sep="\t") 
SampleInfo=SampleInfo[,c('GSMNumber',"CONDITION","AGE","SEX")]
metaInfo=LauAD@meta.data
meta=merge(metaInfo,SampleInfo,by.x='orig.ident',by.y='GSMNumber')
all(metaInfo$nCount_RNA==meta$nCount_RNA)
LauAD$CONDITION=meta$CONDITION
LauAD$AGE=meta$AGE
LauAD$SEX=meta$SEX

pdf("LauADDimPlot.pdf")
DimPlot(LauAD, label=TRUE)&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("LauADMarkerFeaturePlot.pdf") #Mapping cells to known cell types using privous validated markers
FeaturePlot(LauAD,c("NRGN","GAD1","AQP4","MBP","CD74","VCAN","FLT1","SLC6A1"))&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
LauAD.markers <- FindAllMarkers(LauAD, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #Double checked the cell types by their marker genes
LauAD.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
pdf("LauADMarkerHeatmapPlot.pdf") 
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
dev.off()

#--------Isolated the control samples from Lau et,al.----------
LauAllCtrl=subset(LauAD,CONDITION %in% "Control")
LauAllCtrl<- RunTSNE(LauAllCtrl, reduction = "pca", dims = 1:30)
LauAllCtrl<- FindNeighbors(LauAllCtrl, reduction = "pca", dims = 1:30)
LauAllCtrl<- FindClusters(LauAllCtrl, resolution = 0.5)
pdf("Control/LauAllCtrlMarkerFeaturePlot.pdf")
FeaturePlot(LauAllCtrl,reduction="tsne",features=c("NRGN","GAD1","AQP4","MBP","CD74","VCAN","FLT1","SLC6A1"))&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("Control/LauAllCtrlDimPlot.pdf")
DimPlot(LauAllCtrl, reduction="tsne",label=TRUE)&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
names=c("Oli","Ast","Ex","Ex","Oli","Opc","In","Ex","Mic","In","Ex","Ex","In","In","Ex","Ex","Ex","Ex","Oli","End","In","Ast","Ex","Ex","Oli","Opc","Ex")
new.cluster.ids <- names
names(new.cluster.ids) <- levels(LauAllCtrl)
LauAllCtrl <- RenameIdents(LauAllCtrl, new.cluster.ids)
LauAllCtrl$cellType=Idents(LauAllCtrl)
pdf("tmp/LauAllCtrlcellType.pdf")
TSNEPlot(LauAllCtrl, label=TRUE,group.by="cellType")&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("tmp/LauAllCtrlcellType8UMAP.pdf")
DimPlot(LauAllCtrl, label=TRUE)&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
saveRDS(LauAllCtrl,file="/Projects/deng/Alzheimer/syn18485175/cellCycle/LauAllData/LauAllCtrl.rds")

#--------Isolated the AD samples from Lau et,al.----------
LauAllAD=subset(LauAD,CONDITION %in% "Alzheimer")
LauAllAD<- RunTSNE(LauAllAD, reduction = "pca", dims = 1:30)
LauAllAD<- FindNeighbors(LauAllAD, reduction = "pca", dims = 1:30)
LauAllAD<- FindClusters(LauAllAD, resolution = 0.5)
pdf("tmp/LauAllADMarkerFeaturePlot.pdf")
FeaturePlot(LauAllAD,reduction="tsne",features=c("NRGN","GAD1","AQP4","MBP","CD74","VCAN","FLT1","SLC6A1"))&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("tmp/LauAllADDimPlot.pdf")
DimPlot(LauAllAD, reduction="tsne",label=TRUE)&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
names=c("Ex","Oli","Ast","Opc","Oli","In","Ex","Mic","Ex","In","Oli","Ex","Ex","In","Ast","End","In","Ex","Ex","Ex","Ex","In","Ex","End","Oli","Ex","Opc")
new.cluster.ids <- names
names(new.cluster.ids) <- levels(LauAllAD)
LauAllAD <- RenameIdents(LauAllAD, new.cluster.ids)
LauAllAD$cellType=Idents(LauAllAD)
pdf("tmp/LauAllADcellType.pdf")
TSNEPlot(LauAllAD, label=TRUE,group.by="cellType")&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("tmp/LauAllADcellType8UMAP.pdf")
DimPlot(LauAllAD, label=TRUE)&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
saveRDS(LauAllAD,file="/Projects/deng/Alzheimer/syn18485175/cellCycle/LauAllData/LauAllAD.rds")

#--------Extract the Ex from all the Lau et,al. dataset and add the cell cycle score-----------------
Ex=subset(LauAD,idents=c(2,3,7,8,11,14,16,18,19,22))
Ex <- Ex[!grepl("^MT-", rownames(Ex)), ]
Ex<- RunTSNE(Ex, reduction = "pca", dims = 1:30)
Ex<- FindNeighbors(Ex, reduction = "pca", dims = 1:30)
Ex<- FindClusters(Ex, resolution = 0.1)
saveRDS(ExCycle,file="/Projects/deng/Alzheimer/syn18485175/cellCycle/LauAllData/LauEx.rds")




#######################################################################################
##############         Preprocessing snRNA from Mathy et al.    #######################
#######################################################################################
setwd("/Projects/deng/Alzheimer/syn18485175")
data <- Read10X(data.dir = "data/") #data were downloaded from https://www.synapse.org/#!Synapse:syn18485175
AD <- CreateSeuratObject(counts = data, project = "AD")
AD <- NormalizeData(AD, normalization.method = "LogNormalize", scale.factor = 10000)
AD <- FindVariableFeatures(AD, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(AD)
AD <- ScaleData(AD, features = all.genes)
AD <- RunPCA(AD, features = VariableFeatures(object = AD))
DimHeatmap(AD , dims = 20:45, cells = 500, balanced = TRUE)
ElbowPlot(AD,ndims = 50)
# t-SNE and Clustering
AD<- RunUMAP(AD, reduction = "pca", dims = 1:30)
AD<- FindNeighbors(AD, reduction = "pca", dims = 1:30)
AD<- FindClusters(AD, resolution = 0.5)
saveRDS(AD,file="/Projects/deng/Alzheimer/syn18485175/AD.rds")

library(Seurat)
library(ggplot2)
library(clustree)
setwd("/Projects/deng/Alzheimer/syn18485175/cellCycle/MathyAllData")
AD=readRDS("/Projects/deng/Alzheimer/syn18485175/AD.rds")
Info=read.table("/Projects/deng/Alzheimer/syn18485175/Phenotype4Cell.txt",header=TRUE,row.names=1,sep="\t")
AD$Gender=Info$Gender
AD$Age=Info$age_death
AD=subset(AD,idents=c(0:21))
Ex=subset(AD,OldCellType %in% c("Ex"))
pdf("ExElbowPlot.pdf")
ElbowPlot(Ex,ndims = 50)
dev.off()
Ex=subset(ExAll,Statues=="Control")
Ex<- RunTSNE(Ex, reduction = "pca", dims = 1:30)
Ex<- FindNeighbors(Ex, reduction = "pca", dims = 1:30)
Ex<- FindClusters(Ex, resolution = 0.1) #label 4

pdf("labelR0.1Ctrl.pdf")
DimPlot(Ex, label=TRUE,reduction="tsne")&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
Ex.markers <- FindAllMarkers(Ex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Ex.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(Ex, features = top5$gene) + NoLegend()
saveRDS(Ex,file="/Projects/deng/Alzheimer/syn18485175/cellCycle/MathyAllData/AllEx8Mathy.rds")


#######################################################################################
##############         Preprocessing snRNA from Yang et al.     #######################
#######################################################################################
setwd("/Projects/deng/Aging/Ex/Yang_GSE159812")
library("GEOquery")
gse <- GEOquery::getGEO("GSE159812", GSEMatrix = TRUE)
phenotype=gse[[1]]@phenoData@data[,c("age:ch1","Sex:ch1","batch:ch1","biogroup:ch1","disease group:ch1","source_name_ch1")]
colnames(phenotype)=c("Age","Sex","BAtch","BioGroup","Group","Region")
write.table(phenotype,file="sampleInfo.txt",sep="\t",quote=F)

setwd("/Projects/deng/Aging/Ex/Yang_GSE159812/rawData")
filenames <- dir(full.names = F,pattern="GSM")
list.Object <- lapply(filenames, FUN = function(x) {
    data <- Read10X(data.dir = x)
    project <- CreateSeuratObject(counts = data, project = x)
})
setwd("/Projects/deng/Aging/Ex/Yang_GSE159812")
COVID19 <- merge(list.Object[[1]], y = list.Object[-1], project = "COVID19_Yang")
min(COVID19$nFeature_RNA) #43
table(COVID19$orig.ident)
COVID19[["percent.mt"]] <- PercentageFeatureSet(COVID19, pattern = "^MT-") #all with 0
pdf("QC.pdf",width=30)
VlnPlot(COVID19, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
COVID19 <- subset(COVID19, subset = nFeature_RNA > 200 & nFeature_RNA < 6000) #
COVID19 <- NormalizeData(COVID19, normalization.method = "LogNormalize", scale.factor = 10000)
COVID19 <- FindVariableFeatures(COVID19, selection.method = "vst", nfeatures = 2000)
COVID19 <- ScaleData(COVID19)
COVID19 <- RunPCA(COVID19, features = VariableFeatures(object = COVID19))
SampleInfo=read.table("sampleInfo.txt",header=T,row.names=1,sep="\t")
SampleInfo$Sample=rownames(SampleInfo)
MetaInfo=COVID19@meta.data
MetaInfoResult=merge(MetaInfo,SampleInfo,by.x="orig.ident",by.y="Sample")
COVID19$Age=MetaInfoResult$Age
COVID19$Group=MetaInfoResult$BioGroup
COVID19$Gender=MetaInfoResult$Sex
COVID19$Region=MetaInfoResult$Region

COVID19Ctrl=subset(COVID19,Group %in% "Control"& Region %in% "Medial prefrontal cortex") #extracted the control samples 
COVID19Ctrl<- RunTSNE(COVID19Ctrl, reduction = "pca", dims = 1:30)
COVID19Ctrl<- FindNeighbors(COVID19Ctrl, reduction = "pca", dims = 1:30)
COVID19Ctrl<- FindClusters(COVID19Ctrl, resolution = 0.5)
pdf("COVID19CtrlLabel.pdf",width=9)
DimPlot(COVID19Ctrl, label=TRUE,reduction="tsne")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("CellTypeCtrlMarker.pdf",width=10,height=9)
FeaturePlot(COVID19Ctrl,features=c("NRGN","GAD1","AQP4","MBP","CD74","VCAN","FLT1","AMBP"),reduction="tsne")&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
new.cluster.ids <- c("Oli","Oli","Oli","Ast","Ex","Ex","Ex","In","In","Opc","Ex","Mic","Ex","In","Ex","Ast","Ex","In","Ex","Ex","Ex","Ex","End")
names(new.cluster.ids) <- levels(COVID19Ctrl)
COVID19Ctrl <- RenameIdents(COVID19Ctrl, new.cluster.ids)
COVID19Ctrl$cellType=Idents(COVID19Ctrl)
pdf("COVID19CtrlCellType.pdf",width=10,height=9)
DimPlot(COVID19Ctrl, label=TRUE,reduction="tsne")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

saveRDS(COVID19,file="/Projects/deng/Aging/Ex/Yang_GSE159812/YangCOVID19.rds")
Ex=subset(COVID19,idents=c("2","15","17","18","19","20","22"))

pdf("COVID19ExLabelInitial.pdf",width=9)
DimPlot(Ex, label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

Ex <- NormalizeData(Ex, normalization.method = "LogNormalize", scale.factor = 10000)
Ex <- FindVariableFeatures(Ex, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Ex)
Ex <- ScaleData(Ex, features = all.genes)
Ex <- RunPCA(Ex, features = VariableFeatures(object = Ex))
pdf("ExElbowPlot.pdf")
ElbowPlot(Ex,ndims = 50)
dev.off()
# t-SNE and Clustering
Ex<- RunTSNE(Ex, reduction = "pca", dims = 1:50)
Ex<- FindNeighbors(Ex, reduction = "pca", dims = 1:50)
Ex<- FindClusters(Ex, resolution = 0.5)

SampleInfo=read.table("sampleInfo.txt",header=T,row.names=1)
SampleInfo$Sample=rownames(SampleInfo)
MetaInfo=Ex@meta.data
MetaInfoResult=merge(MetaInfo,SampleInfo,by.x="orig.ident",by.y="Sample")
Ex$Age=MetaInfoResult$Age
Ex$Group=MetaInfoResult$Group
Ex$Gender=MetaInfoResult$Gender
saveRDS(Ex,file="/Projects/deng/Aging/Ex/Yang_GSE159812/YangEx.rds")


#######################################################################################
##############         Preprocessing snRNA from Nagy et al.     #######################
#######################################################################################
#GSE144136 Major depressive disorder (MDD)  34 sample post-mortem dosolateral prefrontal cortex   All subjects were male. 

library(Seurat)
setwd("/Projects/deng/Aging/Ex/Nagy_GSE144136")
MDD.data <- Read10X(data.dir = "rawData/",gene.column = 1)
MDD <- CreateSeuratObject(counts = MDD.data, project = "MDD8Nagy")
summary(MDD$nFeature_RNA) #43
table(MDD$orig.ident)
MDD[["percent.mt"]] <- PercentageFeatureSet(MDD, pattern = "^MT-") #all with 0
pdf("QC.pdf",width=30)
VlnPlot(MDD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
MDD <- NormalizeData(MDD, normalization.method = "LogNormalize", scale.factor = 10000)
MDD <- FindVariableFeatures(MDD, selection.method = "vst", nfeatures = 2000)
MDD <- ScaleData(MDD)
MDD <- RunPCA(MDD, features = VariableFeatures(object = MDD))
pdf("NagyElbowPlot.pdf")
ElbowPlot(MDD,ndims = 50)
dev.off()
# t-SNE and Clustering
MDD<- RunTSNE(MDD, reduction = "pca", dims = 1:50)
MDD<- FindNeighbors(MDD, reduction = "pca", dims = 1:50)
MDD<- FindClusters(MDD, resolution = 0.5)
# Visualization
pdf("MDDLabelTNSE.pdf",width=9)
DimPlot(MDD, reduction="tsne",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("CellTypeMarkerTNSE.pdf",width=10,height=9)
FeaturePlot(MDD,reduction="tsne",features=c("NRGN","GAD1","AQP4","MBP","CD74","VCAN","FLT1","AMBP"))&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("MDDcellType8AuthorTNSE.pdf",width=13)
DimPlot(MDD, reduction="tsne",group.by="orig.ident")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

cellInfo=read.table("cellInfo.txt",header=T)
MDD$Group=cellInfo$Group
MDD$sampleID=cellInfo$sampleID
saveRDS(MDD,file="NagyMDD.rds")

MDDCtrl=subset(MDD,Group %in% "Control")
NagyAllCtrl <- RunTSNE(NagyAllCtrl, dims = 1:30)
pdf("MDDCtrlCellTypeMarkerTNSE.pdf",width=10,height=9)
FeaturePlot(NagyAllCtrl,reduction="tsne",features=c("NRGN","GAD1","AQP4","MBP","CD74","VCAN","FLT1","AMBP"))&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("MDDCtrlLabelTNSE.pdf",width=9)
DimPlot(NagyAllCtrl, reduction="tsne",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

cellType=read.table("meta.data.cellType.txt",header=T)
NagyAllCtrl$cellType=cellType$cellType
pdf("MDDCtrlcellType8AuthorTNSE.pdf")
DimPlot(NagyAllCtrl, reduction="tsne",group.by="cellType")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

saveRDS(NagyAllCtrl,file="/Projects/deng/Aging/Ex/Nagy_GSE144136/MDDCtrl.rds")

Ex=subset(MDD,orig.ident %in% "Ex")
pdf("ExLabelInitial.pdf",width=9)
DimPlot(Ex, reduction="tsne",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

Ex <- NormalizeData(Ex, normalization.method = "LogNormalize", scale.factor = 10000)
Ex <- FindVariableFeatures(Ex, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Ex)
Ex <- ScaleData(Ex, features = all.genes)
Ex <- RunPCA(Ex, features = VariableFeatures(object = Ex))
pdf("ExElbowPlot.pdf")
ElbowPlot(Ex,ndims = 50)
dev.off()
# t-SNE and Clustering
Ex<- RunTSNE(Ex, reduction = "pca", dims = 1:50)
Ex<- FindNeighbors(Ex, reduction = "pca", dims = 1:50)
Ex<- FindClusters(Ex, resolution = 0.2)
pdf("ExLabel0.3.pdf",width=9)
DimPlot(ExTmp, reduction = "tsne", label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("ExnCount_RNA.pdf",width=8)
FeaturePlot(Ex,c("nCount_RNA"),reduction="tsne")&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("ExGroup.pdf",width=8)
DimPlot(Ex, reduction = "tsne", group.by="Group",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("ExSample.pdf",width=8)
DimPlot(Ex, reduction = "tsne", group.by="sampleID")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

saveRDS(Ex,file="/Projects/deng/Aging/Ex/Nagy_GSE144136/NagyEx.rds")