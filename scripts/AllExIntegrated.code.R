library(Seurat)
library(ggplot2)
library(clustree)

#######################################################################################
##############      Merge snRNA profiles from all samples       #######################
#######################################################################################

MathyEx=readRDS("/Projects/deng/Alzheimer/syn18485175/cellCycle/MathyAllData/AllEx8Mathy.rds") #prefrontal cortex, For AD, age>75 years old  BA10
MathyExCtrl=subset(MathyEx,Statues %in% "Control") #24 sample, 16952 cells
MathyExCtrl$Gender=MathyExCtrl$Gender
MathyExCtrl$Age=MathyExCtrl$Age
MathyExCtrl$Source="Mathy"
MathyExCtrl$orig.ident=MathyExCtrl$ProjectID
MathyExCtrl$Statues="Control"

MathyExAD=subset(MathyEx,Statues %in% "Alzheimer") #24 sample, 16952 cells
MathyExAD$Gender=MathyExAD$Gender
MathyExAD$Age=MathyExAD$Age
MathyExAD$Source="Mathy"
MathyExAD$orig.ident=MathyExAD$ProjectID
MathyExAD$Statues="Alzheimer"

YangEx=readRDS("/Projects/deng/Aging/Ex/Yang_GSE159812/YangIntegratedEx.rds") #Medial prefrontal cortex, For COVID, age>60 
YangExCtrl=subset(YangEx,Group %in% "Control") #8 samples, 8360 cells
tmp=ifelse(YangExCtrl$Gender == "F", "Female", "Male")
YangExCtrl$Gender=tmp
YangExCtrl$Age=YangExCtrl$Age
YangExCtrl$Source="Yang"
YangExCtrl$orig.ident=YangExCtrl$orig.ident
YangExCtrl$Statues="Control"

NagyEx=readRDS("/Projects/deng/Aging/Ex/Nagy_GSE144136/NagyEx.rds") #post-mortem dorsolateral prefrontal cortex (BA9), For MDD, age= 38.4(4)
NagyExCtrl=subset(NagyEx,Group %in% "Control") #17 samples, 18555 cells
NagyExCtrl$Gender="Male" #all male
NagyExCtrl$Age=38
NagyExCtrl$Source="Nagy"
NagyExCtrl$orig.ident=paste0("B",NagyExCtrl$sampleID,sep="")
NagyExCtrl$Statues="Control"

LauEx=readRDS("/Projects/deng/Alzheimer/syn18485175/cellCycle/LauAllData/LauEx.rds")#post-mortem dorsolateral prefrontal cortex (BA9), For AD, age= 60-95
LauExCtrl=subset(LauEx,CONDITION %in% "Control")#12 samples, 30215 cells
LauExCtrl$Gender=LauExCtrl$SEX
LauExCtrl$Age=LauExCtrl$AGE
LauExCtrl$Source="Lau"
LauExCtrl$orig.ident=LauExCtrl$orig.ident
LauExCtrl$Statues="Control"

LauEx=readRDS("/Projects/deng/Alzheimer/syn18485175/cellCycle/LauAllData/LauEx.rds")#post-mortem dorsolateral prefrontal cortex (BA9), For AD, age= 60-95
LauExAD=subset(LauEx,CONDITION %in% "Alzheimer")#
LauExAD$Gender=LauExAD$SEX
LauExAD$Age=LauExAD$AGE
LauExAD$Source="Lau"
LauExAD$orig.ident=LauExAD$orig.ident
LauExAD$Statues="Alzheimer"

setwd("/Projects/deng/Aging/Ex/AllEx")
ExAll=merge(MathyExCtrl,y=c(MathyExAD,LauExCtrl,LauExAD,YangExCtrl,NagyExCtrl))
table(paste0(ExAll$Source,"_",ExAll$Statues,sep=""))

saveRDS(ExAll,"ExAllMerge.rds")



#######################################################################################
##############   Integrated different datasets into a whole datasets      #############
#######################################################################################

ExAll=readRDS("ExAllMerge.rds")
table(paste0(ExAll$Source,"_",ExAll$Statues,sep=""))
ExAll[["percent.mt"]] <- PercentageFeatureSet(ExAll, pattern = "^MT-")
selected_mito <- WhichCells(ExAll, expression = percent.mt < 25)
ExAllTmp <- subset(ExAll, cells = selected_mito)
ExAllTmp <- ExAllTmp[!grepl("^MT-", rownames(ExAllTmp)), ]
length(selected_mito)
memory.limit(size=400000)
ExAll=ExAllTmp
ExAll.list <- SplitObject(ExAll, split.by = "Source")
ExAll.list <- lapply(X = ExAll.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = ExAll.list)
ExAll.list <- lapply(X = ExAll.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = ExAll.list, reference = c(1), reduction = "rpca",dims = 1:50)
ExAll.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
ExAll.integrated <- ScaleData(ExAll.integrated, verbose = FALSE)
ExAll.integrated <- RunPCA(ExAll.integrated, verbose = FALSE)
pdf("ExAll.integratedElbowPlot.pdf")
ElbowPlot(ExAll.integrated,ndims = 50)
dev.off()
# t-SNE and Clustering
ExAll.integrated <- RunTSNE(ExAll.integrated, dims = 1:50,check_duplicates = FALSE)
ExAll.integrated<- FindNeighbors(ExAll.integrated, reduction = "pca", dims = 1:30)


for(i in colnames(ExAll.integrated@meta.data)) {
  ExAll.integrated[[i]] <- NULL
}
ExAll.integrated$orig.ident=ExAll.integrated$orig.ident
ExAll.integrated$nCount_RNA=ExAll.integrated$nCount_RNA
ExAll.integrated$nFeature_RNA=ExAll.integrated$nFeature_RNA
ExAll.integrated$seurat_clusters=ExAll.integrated$seurat_clusters
ExAll.integrated$Gender=ExAll.integrated$Gender
ExAll.integrated$Age=ExAll.integrated$Age
ExAll.integrated$Source=ExAll.integrated$Source
ExAll.integrated$Statues=ExAll.integrated$Statues

#----confirmation the resolution------------------
ExAll.integratedTmp <- FindClusters(ExAll.integrated, resolution = c(seq(0.1,1,0.1)))
clus.tree.out <- clustree(ExAll.integratedTmp)
pdf(file = "ExAll.integrated.Resolution.Tree.pdf", width = 12, height = 10)
print(clus.tree.out)
dev.off()

ExAll.integrated<- FindClusters(ExAll.integrated, resolution = 0.5)
saveRDS(ExAll.integrated,"/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")


#######################################################################################
##########    Initial analysis of the Ex contain both AD and control      #############
#######################################################################################

ExAll.integrated=subset(ExAll.integrated,idents=c(0:15)) #remove cluster 16 by their smaller cell number

#-----Figure S8C---------------
tiff("ExAll.integratedLabel.tiff",width=300,height=270)
DimPlot(ExAll.integrated, raster=FALSE,reduction = "tsne")&NoLegend()&theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()


#-----Figure 3A---------------
#-----relationship between subcluster---------------
Ex=ExAll.integrated
DefaultAssay(Ex)="RNA"
ExExpr <- AverageExpression(Ex)[["RNA"]]
hclust=hclust(as.dist(1-cor(ExExpr)),method="ward.D2")
pdf("clusterRelationship.pdf",width=10,height=3)
plot(hclust,hang=-1)
dev.off()
#-----expression of layer marker---------------
order=hclust$labels[hclust$order]
layerMarker=c("GLRA3","TLE1","TLE3","MDGA1","LUX2","UNC5D","GPR6","MEF2C","DTX4","CUX1","CUX2","KITLG","SATB2","RASGRF2","PVRL3","PRSS12","RORB","NEFH","CNTN6","FOXO1","OPN3","LIX1","SYT9","S100A10","LDB2","CRIM1","KCNK2","SULF2","PCP4","HTR2C","FEZF2","BCL11B","CRYM","OTX1","SOX5","TOX","ETV1","RPRM","RXFP1","FOXP2","SEMA3E","NR4A3","LXN","PPP1R1B","SYT6","OPRK1","NR4A2","SYNPR","TLE","NTNG2","ADRA2A")
t=DotPlot(Ex,features=layerMarker)
t=ggplot(t$data, aes(factor(features.plot,levels=unique(features.plot)),factor(id,levels=order),size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_color_gradient2(low="white",mid="white",high = "red")+
labs(color="Expression",size="Percent",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0),face="bold"),axis.text.x = element_text(size=rel(1.0),face="bold",angle=90,vjust=1,hjust=1),axis.text.y = element_text(size=rel(1.0),face="bold")) 
pdf("layerMarker.pdf",width=12,height=5)
print(t)
dev.off()


