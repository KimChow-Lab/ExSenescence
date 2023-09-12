library(Seurat)
library(ggplot2)
library(clustree)
library(RColorBrewer)
library(dplyr)
library(scCustomize)

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
saveRDS(ExAll.integrated,"ExAll.integratedTSNE.rds")


#######################################################################################
##########    Phenotype visualization      #############
#######################################################################################
setwd("/Projects/deng/Aging/Ex/AllEx")
ExAll.integrated=readRDS("ExAll.integratedTSNE.rds")
ExAll.integrated=subset(ExAll.integrated,idents=c(0:15)) #remove cluster 16 by their smaller cell number

colorPalette=colorRampPalette(brewer.pal(8, "Set3"))(length(unique(ExAll.integrated$seurat_clusters))) 
colorPalette[4]="Orange"
colorPalette[6]="Firebrick3"
colorPalette[11]="LightCyan"

pdf("ExAll.Cluster.pdf",width=4,height=3.5)
DimPlot(ExAll.integrated, raster=TRUE,reduction = "tsne",label=T,label.size = 6,cols =colorPalette)&NoLegend()&theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

pdf("ExAll.Source.pdf",width=12,height=3.5)
DimPlot(ExAll.integrated,split.by="Source", raster=TRUE,reduction = "tsne",label=T,label.size = 6,cols =colorPalette)&NoLegend()&theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

pdf("ExAll.Status.pdf",width=7,height=3.5)
DimPlot(ExAll.integrated,split.by="Statues", raster=TRUE,reduction = "tsne",cols =colorPalette)&NoLegend()&theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

ExAll.integrated$SourceStatus=paste0(ExAll.integrated$Statues,"_",ExAll.integrated$Source,sep="")
ExAll.integrated$SourceStatus=factor(ExAll.integrated$SourceStatus,levels=c("Control_Mathy","Control_Lau","Control_Nagy","Control_Yang","Alzheimer_Mathy","Alzheimer_Lau"))
pdf("ExAll.SourceStatus.pdf",width=6,height=8)
DimPlot(ExAll.integrated,split.by="SourceStatus", raster=TRUE,reduction = "tsne",cols =colorPalette,ncol=2)&NoLegend()&theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("ExAll.SourceStatus_group.pdf",width=6,height=8)
DimPlot(ExAll.integrated,split.by="SourceStatus", raster=TRUE,reduction = "tsne",group.by="SourceStatus",ncol=2,cols=c("Goldenrod","DarkOrange","Yellow","OliveDrab","PaleGoldenrod","PeachPuff"))&NoLegend()&NoAxes()
dev.off()
table(ExAll.integrated$SourceStatus)
#Control_Mathy     Control_Lau    Control_Nagy    Control_Yang Alzheimer_Mathy
#  16936           29937           18390            8084           17863
#Alzheimer_Lau
#  32002

ExAll.integrated$GenderStatus=paste0(ExAll.integrated$Statues,"_",ExAll.integrated$Gender,sep="")
ExAll.integrated$GenderStatus=factor(ExAll.integrated$GenderStatus,levels=c("Control_Male","Control_Female","Alzheimer_Male","Alzheimer_Female"))
pdf("ExAll.GenderStatus.pdf",width=6,height=5)
DimPlot(ExAll.integrated,split.by="GenderStatus", raster=TRUE,reduction = "tsne",cols =colorPalette,ncol=2)&theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
#Control_Male   Control_Female   Alzheimer_Male Alzheimer_Female
#51923            21424            32778            17087

ExAll.integrated$geneCount=ExAll.integrated$nCount_RNA/1000
pdf("ExAll.integrated.nCount_RNA.pdf",width=12,height=4)
VlnPlot(ExAll.integrated,features=c("geneCount"),ncol=1,pt.size=0,cols=colorPalette)&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

ExAll.integrated$geneNumber=ExAll.integrated$nFeature_RNA/1000
pdf("ExAll.integrated.nFeature_RNA.pdf",width=12,height=4)
VlnPlot(ExAll.integrated,features=c("geneNumber"),ncol=1,pt.size=0,cols=colorPalette)&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

metaData=ExAll.integrated@meta.data
geneCount <- metaData %>%
          group_by(seurat_clusters) %>%
          summarize(mean=mean(nCount_RNA))
geneNumber <- metaData %>%
          group_by(seurat_clusters) %>%
          summarize(mean=mean(nFeature_RNA))
write.table(data.frame(geneCount=geneCount,geneNumber=geneNumber),file="geneCountAndNumber_mean.txt",sep="\t",quote=F)


phenotype=unique(ExAll.integrated@meta.data[,c("orig.ident","Source","Age","Gender","Statues")])
rownames(phenotype)=phenotype$orig.ident
geneCount <- metaData %>% group_by(orig.ident) %>% summarize(mean=mean(nCount_RNA))
phenotype=phenotype[geneCount$orig.ident,]
all(phenotype$orig.ident==geneCount$orig.ident)

geneNumber <- metaData %>% group_by(orig.ident) %>% summarize(mean=mean(nFeature_RNA))
all(phenotype$orig.ident==geneNumber$orig.ident)

cellNumber=data.frame(table(ExAll.integrated$orig.ident))
rownames(cellNumber)=cellNumber[,1]
cellNumber=cellNumber[geneCount$orig.ident,]
all(phenotype$orig.ident==cellNumber[,1])

phenotype$AveGeneNumber=geneNumber$mean
phenotype$AveGeneCount=geneCount$mean
phenotype$nucleiNumber=cellNumber[,2]
write.table(phenotype,file="PhenotypeForAllSample.txt",sep="\t",quote=F)


#### cell number from different source in  each cluster####
#Figure S8D
setwd("/Projects/deng/Aging/Ex/AllEx")
ExAll.integrated=readRDS("ExAll.integratedTSNE.rds")
ExAll.integrated=subset(ExAll.integrated,idents=c(0:15)) #remove cluster 16 by their smaller cell number

t=as.matrix(table(ExAll.integrated$Source))
tmp=paste0(ExAll.integrated$seurat_clusters,"_",ExAll.integrated$Statues,"-",ExAll.integrated$Source,sep="")
tmp=as.matrix(table(tmp))
TmpInfo=as.matrix(do.call(rbind, strsplit(as.character(rownames(tmp)),'_')))
result=data.frame("Cluster"=TmpInfo[,1],"SourceStatues"=TmpInfo[,2],"Count"=tmp[,1])
ClusterList=factor(result$Cluster,levels=rev(c(0:15)))
result$SourceStatues=factor(result$SourceStatues,levels=rev(c("Control-Mathy","Control-Lau","Control-Nagy","Control-Yang","Alzheimer-Mathy","Alzheimer-Lau")))
colorPalette=rev(c("Goldenrod","DarkOrange","Yellow","OliveDrab","PaleGoldenrod","PeachPuff"))
g=ggplot(result, aes(ClusterList, Count, fill=SourceStatues)) +
  geom_bar(stat="identity",position="fill") +
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=colorPalette)+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 0))+
  theme(axis.title.x = element_blank())
pdf("ExAll.integrated.SourceStatuesDistribution.pdf",height=8,width=5)
print(g)
dev.off()
result=result[order(result$Cluster),]
write.table(result,file="ExAll.integrated.cellNumberRatioInEachSource.txt",sep="\t",quote=F)


#######################################################################################
##########   Identify the cell cycle re-entry neurons    ###########################
#######################################################################################
setwd("/Projects/deng/Aging/Ex/AllEx")
ExAll.integrated=readRDS("ExAll.integratedTSNE.rds")
ExAll.integrated=subset(ExAll.integrated,idents=c(0:15)) #remove cluster 16 by their smaller cell number
DefaultAssay(ExAll.integrated)="RNA"

selected_f_10xGenomic <- rownames(ExAll.integrated)[Matrix::rowSums(ExAll.integrated) > ncol(ExAll.integrated)*0.001]
length(selected_f_10xGenomic)
#25834
############ cell cycle phase score calculation for integrated cells  ##############
s.genes <- intersect(cc.genes$s.genes,selected_f)
g2m.genes <- intersect(cc.genes$g2m.genes,selected_f)

ExAll.integrated <- CellCycleScoring(ExAll.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
colorPalette=colorRampPalette(brewer.pal(8, "Set3"))(length(unique(ExAll.integrated$seurat_clusters))) 
colorPalette[4]="Orange"
colorPalette[6]="Firebrick3"
colorPalette[11]="LightCyan"
pdf("IntegratedADCtrlCells_Phase_PvalueDis/ExAll.integrated.cellCycleScoreyBySeurat.pdf",width=6,height=3)
VlnPlot(ExAll.integrated,c("S.Score","G2M.Score"),group.by="seurat_clusters",ncol=1,pt.size=0,cols =colorPalette)&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", linewidth=1.5, linetype="solid"))
dev.off()


DefaultAssay(ExAll.integrated)="RNA"
CCGene=read.table("/Projects/deng/Aging/Ex/cellCycleGene/CCGeneCombine.txt",header=T,sep="\t") #500
genes=intersect(selected_f_10xGenomic,CCGene$NAME)
length(genes) #341
write.table(CCGene[CCGene$NAME%in%genes,],file="341CellCycleGeneIn10XGenomic.txt",sep="\t",quote=F)
CCGene=CCGene[CCGene$NAME%in%genes,]
CCGenelist=split(CCGene$NAME,CCGene$PHASE)
ExAll.integrated <- AddModuleScore(
  object = ExAll.integrated,
  features = CCGenelist,
  name = names(CCGenelist),
  ctrl = 100,
)
colorPalette=colorRampPalette(brewer.pal(8, "Set3"))(length(unique(ExAll.integrated$seurat_clusters))) 
colorPalette[4]="Orange"
colorPalette[6]="Firebrick3"
colorPalette[11]="LightCyan"
pdf("IntegratedADCtrlCells_Phase_PvalueDis/ExAll.integrated.cellCycleScore.pdf",width=6)
VlnPlot(ExAll.integrated,c("G1.S1","S.phase5","G22","G2.M3","M.G14"),group.by="seurat_clusters",ncol=1,pt.size=0,cols =colorPalette)&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", linewidth=1.5, linetype="solid"))
dev.off()


cluster=c(0:15)
PhaseList=c("G1.S1","S.phase5","G22","G2.M3","M.G14")
CycleScore=ExAll.integrated@meta.data[,c("seurat_clusters","G1.S1","S.phase5","G22","G2.M3","M.G14")]

for(p in PhaseList){
 PvalueBetCluster=matrix(data = NA, nrow = length(cluster), ncol = length(cluster), dimnames = list(cluster, cluster))
 for(i in 1:length(cluster)){
  for( j in 1:length(cluster)){
    C1_score=CycleScore[CycleScore$seurat_clusters==cluster[i],p]
    C2_score=CycleScore[CycleScore$seurat_clusters==cluster[j],p]
    delta=mean(C1_score)-mean(C2_score) #the row represent the current cluster
    pval=t.test(C1_score,C2_score)$p.value
    PvalueBetCluster[i,j]=delta*(-log10(pval))/abs(delta)
    if(pval==0){
        PvalueBetCluster[i,j]=500;
    }
    if(i==j||delta<0){
      PvalueBetCluster[i,j]=0;
    }
  }
 }
 PvalueBetCluster=data.frame(PvalueBetCluster,check.names=F)
 PvalueBetCluster$Cluster=as.character(rownames(PvalueBetCluster))
 PvalueBetCluster.df=reshape2::melt(PvalueBetCluster,id=(length(cluster)+1))
 PvalueBetCluster.df$Cluster=factor(PvalueBetCluster.df$Cluster,levels=cluster)

 t=ggplot(PvalueBetCluster.df, aes(x=Cluster, y=value, color=Cluster))+
 geom_boxplot(outlier.shape = NA,outlier.colour=NA)+labs(title="",x="", y = "")+
 geom_jitter(size=1)+
 scale_color_manual(values=colorPalette)+
 theme_classic()+
 theme(legend.position="none")+
 theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
 pdf(paste0("IntegratedADCtrlCells_Phase_PvalueDis/IntegratedADCtrlCells_",p,"_PvalueDis.pdf",sep=""),height=2,width=6)
 print(t)
 dev.off()
 
 colnames(PvalueBetCluster.df)=c("CurrentCluster","OtherClusters","-log10(pvalue).adjust")
 write.table(PvalueBetCluster.df,file=paste0("IntegratedADCtrlCells_Phase_PvalueDis/IntegratedADCtrlCells_",p,"_PvalueDis.txt",sep=""),sep="\t",quote=F,row.names=F)
}




############ InferCNV to define the copy number variation  ##############
tmux new -s InferCNV
tmux attach-session -t InferCNV
setwd("/data2/deng/Aging/Ex/AllEx/InferCNV")
ExAll.integrated=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
DefaultAssay(ExAll.integrated)="RNA"
ExAll.Downsampling=ExAll.integrated[, sample(colnames(ExAll.integrated), size = 10000, replace=F)]
saveRDS(ExAll.Downsampling,file="ExAll.Downsampling.rds")
ExAll.Downsampling <- ExAll.Downsampling[!grepl("^MT-", rownames(ExAll.Downsampling)), ] #remove the genes coded by mitocondrial
counts_matrix = as.matrix(ExAll.Downsampling@assays$RNA@counts)
pData=ExAll.Downsampling@meta.data$seurat_clusters
names(pData)=rownames(ExAll.Downsampling@meta.data)
pData=data.frame(pData)
annotation=read.table("/Projects/deng/refGene/gencode.v38.annotation.removeDupliSymbol.geneType.txt",header=TRUE,row.names=1)
anno=data.frame("Chr"=annotation$Chr,"Start"=annotation$Start,"End"=annotation$End)
rownames(anno)=annotation$Symbol
genes=intersect(rownames(counts_matrix),rownames(anno))
length(genes)#31566
anno=anno[genes,]
anno$Chr=factor(anno$Chr,levels=paste0("chr",c(1:22,"X","Y"),sep=""))
anno=anno[order(anno$Chr,anno$Star,anno$End),]
genes=intersect(rownames(anno),rownames(counts_matrix))
length(genes)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix[genes,],
                                      annotations_file=pData,
                                      delim="\t",
                                      gene_order_file=anno[genes,],
                                      ref_group_names=NULL)
infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics   0.1= 8037 genes and 5645 cells
                               denoise=TRUE,
                               cluster_by_groups = TRUE,
                               HMM=TRUE,
                               output_format = "png",
                               out_dir="/data2/deng/Aging/Ex/AllEx/InferCNV/ExAllDownsampling"
                               )


####perform the infercnv analysis with dataset dependent pattern ############
library(infercnv)
ExAll.integrated=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
ExAll.integrated=subset(ExAll.integrated,seurat_clusters%in%c(0:15))
ExAll.integrated$SourceStatus=paste0(ExAll.integrated$Source,"_",ExAll.integrated$Statues,sep="")
table(ExAll.integrated$SourceStatus)
DefaultAssay(ExAll.integrated)="RNA"
for(s in unique(ExAll.integrated$SourceStatus)){
  Source.ctrl=subset(ExAll.integrated,SourceStatus==s)
  DefaultAssay(Source.ctrl)="RNA"
  Source.ctrl <- Source.ctrl[!grepl("^MT-", rownames(Source.ctrl)), ] #remove the genes coded by mitocondrial
  if(ncol(Source.ctrl)>10000){
     Source.ctrl=Source.ctrl[, sample(colnames(Source.ctrl), size = 10000, replace=F)]
  }
  cellCounts=table(Source.ctrl$seurat_clusters)
  Source.ctrl <- subset(Source.ctrl,seurat_clusters%in%names(cellCounts[cellCounts>10]))#removed the subtype with cells smaller than 10 cells, this subcluster (cells=1) could lead to the 
  counts_matrix = as.matrix(Source.ctrl@assays$RNA@counts)
  pData=Source.ctrl@meta.data$seurat_clusters
  names(pData)=rownames(Source.ctrl@meta.data)
  pData=data.frame(pData)

  annotation=read.table("/Projects/deng/refGene/gencode.v38.annotation.removeDupliSymbol.geneType.txt",header=TRUE,row.names=1)
  anno=data.frame("Chr"=annotation$Chr,"Start"=annotation$Start,"End"=annotation$End)
  rownames(anno)=annotation$Symbol
  genes=intersect(rownames(counts_matrix),rownames(anno))
  print(paste0(s,":",length(genes),sep=""))#31566
  anno=anno[genes,]
  anno$Chr=factor(anno$Chr,levels=paste0("chr",c(1:22,"X","Y"),sep=""))
  anno=anno[order(anno$Chr,anno$Star,anno$End),]
  genes=intersect(rownames(anno),rownames(counts_matrix))

  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix[genes,rownames(pData)],
                                      annotations_file=pData,
                                      delim="\t",
                                      gene_order_file=anno[genes,],
                                      ref_group_names=NULL)
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics   0.1= 8037 genes and 5645 cells
                               denoise=TRUE,
                               cluster_by_groups = TRUE,
                               HMM=TRUE,
                               output_format = "png",
                               out_dir=paste0("/data2/deng/Aging/Ex/AllEx/InferCNV/",s,sep="")
                               )
}



#####################################################################################################
##########    Exploring the features of senescence neurons    ######################################
#####################################################################################################
setwd("/data2/deng/Aging/Ex/AllEx/C5Feature")
ExAll.integrated=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
ExAll.integrated=subset(ExAll.integrated,seurat_clusters%in%c(0:15))

#-----layer information for each subcluster---------------
DefaultAssay(ExAll.integrated)="RNA"
ExExpr <- AverageExpression(ExAll.integrated)[["RNA"]]
hclust=hclust(as.dist(1-cor(ExExpr)),method="ward.D2")
pdf("clusterRelationship.pdf",width=10)
plot(hclust,hang=-1)
dev.off()
order=hclust$labels[hclust$order]

layerMarker=c("GLRA3","TLE1","MDGA1","CARTPT","CUX1","CUX2","KITLG","RASGRF2","PVRL3","PRSS12","DTX4","GPR6","UNC5D","MEF2C","RORB","SATB2","DKK3","LDB2","SYT9","OPN3","NEFH","CNTN6","FOXO1","LIX1","S100A10","CRIM1","KCNK2","SULF2","PCP4","HTR2C","FEZF2","TOX","ETV1","RPRM","RXFP1","FOXP2","CRYM","OTX1","SOX5","SEMA3E","NR4A3","LXN","PPP1R1B","SYT6","OPRK1","NR4A2","SYNPR","NTNG2","ADRA2A")
t=DotPlot(ExAll.integrated,features=layerMarker)
order=rev(c(0,2,3,8,7,9,6,13,1,4,10,15,11,14,12,5))
t=ggplot(t$data, aes(factor(features.plot,levels=unique(features.plot)),factor(id,levels=order),size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_color_gradient2(low="white",mid="white",high = "red")+
labs(color="Expression",size="Percent",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0),face="bold"),axis.text.x = element_text(size=rel(1.0),face="bold",angle=90,vjust=1,hjust=1),axis.text.y = element_text(size=rel(1.0),face="bold")) 
pdf("layerMarkerExpression_test.pdf",width=12,height=5)
print(t)
dev.off()


setwd("/data2/deng/Aging/Ex/AllEx/C5Feature")
ExAll.integrated=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
ExAll.integrated$rpcaCluster=ExAll.integrated$seurat_clusters
DefaultAssay(ExAll.integrated)="integrated"
ExAll.integrated <- RunUMAP(ExAll.integrated, dims = 1:50)
ExAll.integrated[["Layer"]]=NA
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(0,2,8),]="L2/3"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(1,6,7,9,13),]="L4/5"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(4,10),]="L4/6"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(15),]="L5"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(11,14),]="L5/6"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(12),]="L6"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(3),]="ES"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(5),]="LS"
ExAll.integrated$Layer=factor(ExAll.integrated$Layer,levels=c("ES","LS","L2/3","L4/5","L4/6","L5","L5/6","L6"))
ExAll.integrated=subset(ExAll.integrated,seurat_clusters%in%c(0:15))

LayerPalette=c("Orange","Firebrick3","#DBE419","#74CF56","#21D66A","#3DBC75","#238B8E","#2E718E")
pdf("/data2/deng/Aging/Ex/AllEx/C5Feature/ExAll.Layer.tsne.pdf",width=4,height=3.5)
DimPlot(ExAll.integrated, raster=TRUE,reduction = "tsne",group.by="Layer",label=T,label.size = 6,cols =LayerPalette)&NoLegend()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("/data2/deng/Aging/Ex/AllEx/C5Feature/ExAll.Layer.umap.pdf",width=4,height=3.5)
DimPlot(ExAll.integrated, raster=TRUE,reduction = "umap",group.by="Layer",label=T,label.size = 6,cols =LayerPalette)&NoLegend()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

DefaultAssay(ExAll.integrated)="RNA"
SelectedLayerMarker=c("GLRA3","CUX2","RORB","RXFP1","PCP4","SEMA3E","NTNG2")
pdf("/data2/deng/Aging/Ex/AllEx/C5Feature/SelectedLayerMarker.VlnPlot.pdf",height=6,width=4)
scCustomize::Stacked_VlnPlot(seurat_object = ExAll.integrated, group.by="Layer",features = SelectedLayerMarker, colors_use=LayerPalette,x_lab_rotate = TRUE, plot_spacing = 0.05)&
theme(axis.title.x=element_blank())
dev.off()

#https://www.nature.com/articles/s41593-020-00787-0
DefaultAssay(ExAll.integrated)="RNA"
pdf("/data2/deng/Aging/Ex/AllEx/C5Feature/ExAll.LayerMarker.pdf",width=20,height=3)
FeaturePlot(ExAll.integrated, raster=TRUE,reduction = "tsne",features=c("GLRA3","CUX2","RORB","PCP4","SEMA3E","RXFP1","NTNG2"),ncol=7)&NoLegend()&theme(plot.title=element_text(face="italic"),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

pdf("/data2/deng/Aging/Ex/AllEx/C5Feature/ExAll.LayerMarker.umap.pdf",width=20,height=3)
FeaturePlot(ExAll.integrated, raster=TRUE,reduction = "umap",features=c("GLRA3","CUX2","RORB","PCP4","SEMA3E","RXFP1","NTNG2"),ncol=7)&NoLegend()&theme(plot.title=element_text(face="italic"),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

#https://www.nature.com/articles/nrm3629
# Exrepssion of cell cycle assocaited genes
DefaultAssay(ExAll.integrated)="RNA"
cellCycleGene=c("CCND1","CCND2","CCND3","CDK4","CDK6","CCNE1","CDK2","CCNA1","CCNA2","CCNB1","CDK1","CDKN1A","CDKN1B","CDKN2A","CDKN2B","CDKN2C","ATR","ATM","CHEK1","CHEK2")
t=DotPlot(ExAll.integrated,features=cellCycleGene,group.by="Layer")
result=t$data
g=ggplot(result, aes(id, factor(features.plot,levels=rev(cellCycleGene)), size= pct.exp,color=avg.exp.scaled)) +geom_point(alpha = 0.9)+
scale_colour_viridis_c()+
labs(color="avg.exp.scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(angle=90,vjust=1,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("CellCycleGene.pdf",width=4,height=5)
print(g)
dev.off()
write.table(result,file="CellCycleGeneExprAcrossSubCluster.txt",sep="\t",quote=F,row.names=F)

gene=c("ATR","ATM","CHEK1","CHEK2")
t=DotPlot(ExAll.integrated,features=gene,group.by="Layer")
result=t$data
g=ggplot(result, aes(id, factor(features.plot,levels=rev(gene)), size= pct.exp,color=avg.exp.scaled)) +geom_point(alpha = 0.9)+
scale_colour_viridis_c()+
labs(color="avg.exp.scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(angle=90,vjust=1,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("ATM.pdf",width=3.5,height=2)
print(g)
dev.off()

gene=c("CCND1","CCND2","CCND3","CDK4","CDK6","CCNE1","CDK2","CCNA1","CCNA2","CCNB1","CDK1")
t=DotPlot(ExAll.integrated,features=gene,group.by="Layer")
result=t$data
g=ggplot(result, aes(id, factor(features.plot,levels=rev(gene)), size= pct.exp,color=avg.exp.scaled)) +geom_point(alpha = 0.9)+
scale_colour_viridis_c()+
labs(color="avg.exp.scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(angle=90,vjust=1,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("cellCycle.pdf",width=3.5,height=3)
print(g)
dev.off()

gene=c("CDKN1A","CDKN1B","CDKN2A","CDKN2B","CDKN2C")
t=DotPlot(ExAll.integrated,features=gene,group.by="Layer")
result=t$data
g=ggplot(result, aes(id, factor(features.plot,levels=rev(gene)), size= pct.exp,color=avg.exp.scaled)) +geom_point(alpha = 0.9)+
scale_colour_viridis_c()+
labs(color="avg.exp.scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(angle=90,vjust=1,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("CDKN.pdf",width=3.5,height=2)
print(g)
dev.off()

#https://www.nature.com/articles/nrm3629 #23877564
gene=c("E2F1","E2F2","E2F3","E2F4","E2F5","E2F6","E2F7","E2F8","RB1","RBL1","RBL2")
t=DotPlot(ExAll.integrated,features=gene,group.by="Layer")
result=t$data
g=ggplot(result, aes(id, factor(features.plot,levels=rev(gene)), size= pct.exp,color=avg.exp.scaled)) +geom_point(alpha = 0.9)+
scale_colour_viridis_c()+
labs(color="avg.exp.scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(angle=90,vjust=1,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("cellCycleTFs.pdf",width=3.5,height=3)
print(g)
dev.off()


#-----expression of cell cycle pathways (KEGG) associated genes : Heatmap---------------
library(msigdbr)
setwd("/data2/deng/Aging/Ex/AllEx/C5Feature")

Idents(ExAll.integrated)=ExAll.integrated$Layer
ExAllxExpr <- AverageExpression(ExAll.integrated)[["RNA"]]

m_df<- msigdbr(species = "Homo sapiens", category = "C2",subcategory="CP:KEGG")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name) 
cellCycle=intersect(fgsea_sets$KEGG_CELL_CYCLE,rownames(ExAllxExpr))
exprTmp=ExAllxExpr[cellCycle,]
exprTmp=exprTmp[rowMeans(exprTmp)>0,]
pdf("cellCycle8KEGG/cellCycleHeatmapTmp.pdf",height=5,wid=4)
pheatmap(exprTmp,clustering_method="ward.D2",scale="row",border_color=NA,show_rownames=F,clustering_distance_cols="correlation",clustering_distance_rows="correlation",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()
t=pheatmap(exprTmp,clustering_method="ward.D2",scale="row",border_color=NA,show_rownames=F,clustering_distance_cols="correlation",clustering_distance_rows="correlation",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()
tree=cutree(t$tree_row,k=5)
anno=data.frame(tree)
anno$Group=paste0("G",anno$tree,sep="")
anno$tree=NULL
ann_colors = list(
    Group = c(G1 = "RoyalBlue", G2 = "Firebrick3", G3 = "CornflowerBlue",G4="Orange",G5="Tomato")
)
pdf("cellCycle8KEGG/cellCycleHeatmapGroup.pdf",height=5,wid=4)
pheatmap(exprTmp,clustering_method="ward.D2",scale="row",annotation_colors = ann_colors,annotation_row=anno,border_color=NA,show_rownames=F,clustering_distance_cols="correlation",clustering_distance_rows="correlation",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

tmp=t$tree_row
order=tmp$labels[tmp$order]
anno$Symbo=rownames(anno)
write.table(anno[order,],file="cellCycle8KEGG/cellCycleGeneGroup.txt",col.names=T,row.names=F,quote=F,sep="\t")

pdf("cellCycle8KEGG/cellCycleHeatmapGroup4Name.pdf",height=20,wid=6)
pheatmap(exprTmp,clustering_method="ward.D2",scale="row",annotation_row=anno,border_color=NA,show_rownames=T,clustering_distance_cols="correlation",clustering_distance_rows="correlation",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

exprTmp=exprTmp[rownames(anno),]
all(rownames(exprTmp)==rownames(anno))
exprTmp=data.frame(exprTmp,check.names=F)
exprTmp$Cluster=anno$Group
write.table(exprTmp,file="cellCycleFromKEGGExpr.txt",sep="\t",quote=F)



#-----senescence associated signatures: multiple sets---------------
library(dplyr)
library(tibble)
library(fgsea)

setwd("/Projects/deng/Aging/Ex/AllControlEx/Senescence")
ExAll.integrated=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
ExAll.integrated=subset(ExAll.integrated,seurat_clusters%in%c(0:15))
DefaultAssay(ExAll.integrated)="RNA"

ExSmall <- ExAll.integrated[, sample(colnames(ExAll.integrated), size = 10000, replace=F)]
t=DotPlot(ExSmall,features=rownames(ExSmall))
ExExpr8DotPlot=na.omit(t$data)

CSFromSegura=read.table("SenescenceAssociatedGeneList/55markersBySegura.txt",header=T)
CSFromCellAge=read.table("SenescenceAssociatedGeneList/GeneListFromCellAge.txt",header=T)
CSFromCellAge$Group="CSFromCellAge"
CSFromGabrielCasella=read.table("SenescenceAssociatedGeneList/ConservedGene8GabrielCasella.txt",header=T)
CSFromBinZhang=read.table("SenescenceAssociatedGeneList/Senescence8BinZhang.txt",header=T)
CSFromBinZhang$Group="CSFromBinZhang"
CSFromKEGG=read.table("SenescenceAssociatedGeneList/hsa04218_gene.txt",header=T)
CSFromReactome=read.table("SenescenceAssociatedGeneList/R-HSA-2559582_gene.txt",header=T)
SASPFromSenMayo=read.table("SenescenceAssociatedGeneList/SenMayoList.txt",header=T)

CombineCSList=rbind(CSFromSegura,CSFromCellAge,CSFromGabrielCasella,CSFromBinZhang,CSFromKEGG,CSFromReactome,SASPFromSenMayo)
CombineCSList=CombineCSList[CombineCSList$GeneSymbol%in%ExExpr8DotPlot$features.plot,]
CombineCSList=split(CombineCSList$GeneSymbol,CombineCSList$Group)
cluster="5"
clusterCell<- ExExpr8DotPlot %>% dplyr::filter(id == cluster) %>% arrange(desc(avg.exp.scaled)) %>% dplyr::select(features.plot, avg.exp.scaled)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(CombineCSList, stats = ranks,eps=0)
ranks=na.omit(ranks)
fgseaRes=fgseaRes[order(fgseaRes$NES,decreasing=T),]


pdf("/data2/deng/Aging/Ex/AllEx/C5Feature/CSGeneList8fgsea.pdf",width=10,height=5)
plotGseaTable(CombineCSList[fgseaRes$pathway], ranks, fgseaRes, gseaParam=3)
dev.off()

pdf("/data2/deng/Aging/Ex/AllEx/C5Feature/Reactome_SASP_R-HSA-2559582.pdf",height=4)
plotEnrichment(CombineCSList[["R-HSA-2559582"]],ranks)
dev.off()
    

p<-ggplot(data=fgseaRes, aes(x=factor(pathway,levels=rev(pathway)), y=-log10(pval),fill=-log10(pval))) +
  scale_fill_gradient(low = "MistyRose", high = "red")+ theme_bw()+
  geom_bar(stat="identity",width=0.6)+coord_flip()
pdf("/data2/deng/Aging/Ex/AllEx/C5Feature/CSGeneList8fgseaPvalue.pdf",width=5,height=5)
print(p)
dev.off()

fwrite(fgseaRes, file="/data2/deng/Aging/Ex/AllEx/C5Feature/CSGeneList8fgsea.txt", sep="\t", sep2=c("", " ", ""))


#-----senescence associated signatures: core signatures---------------
ExAll.integrated[["Layer"]]=NA
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(0,2,8),]="L2/3"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(1,6,7,9,13),]="L4/5"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(4,10),]="L4/6"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(15),]="L5"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(11,14),]="L5/6"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(12),]="L6"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(3),]="ES"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(5),]="LS"
ExAll.integrated$Layer=factor(ExAll.integrated$Layer,levels=c("ES","LS","L2/3","L4/5","L4/6","L5","L5/6","L6"))

coreCellularGene=read.table("/Projects/deng/Aging/Ex/AllControlEx/Senescence/SenescenceAssociatedGeneList/55markersBySegura.txt",header=T)
DefaultAssay(ExAll.integrated)="RNA"
ExExpr <- AverageExpression(ExAll.integrated)[["RNA"]]
hclust=hclust(as.dist(1-cor(ExExpr)),method="ward.D2")
order=hclust$labels[hclust$order]
t=DotPlot(ExAll.integrated,features=coreCellularGene[,1],group.by="Layer")
result=merge(t$data,coreCellularGene,by.x="features.plot",by.y="GeneSymbol")
LS=result[result$id=="LS",]
LS=LS[order(LS$Group,LS$pct.exp),]
geneOrder=LS$features.plot
g=ggplot(result, aes(factor(features.plot,levels=rev(geneOrder)),factor(id,levels=c("ES","LS","L2/3","L4/5","L4/6","L5","L5/6","L6")), size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_color_gradient2(midpoint=0, low="Blue", mid="white",high="Red")+
labs(color="Average",size="Percent",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),face="italic",angle=90,vjust=0.5,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("/data2/deng/Aging/Ex/AllEx/C5Feature/SignatureCellularGene.pdf",width=10,height=3)
print(g)
dev.off()



# Exrepssion of the SASP assocaited genes
SASPmarker=read.table("/Projects/deng/Aging/Ex/AllEx/Senescence/SASPmarker.txt",header=T,sep="\t")
t=DotPlot(ExAll.integrated,features=SASPmarker$GeneSymobl,group.by="Layer")
result=merge(t$data,SASPmarker,by.x="features.plot",by.y="GeneSymobl")
LS=result[result$id=="LS",]
LS=LS[order(LS$pct.exp),]
geneOrder=LS$features.plot
#result=result[result$id %in% c(2:6),]
result$Group=ifelse(result$Group%in%c("Group5_1","Group5_2","Group5_3"),"Group5",result$Group)
g=ggplot(result, aes(factor(features.plot,levels=rev(geneOrder)), id, size= pct.exp,color=avg.exp.scaled)) +
geom_point(alpha = 0.9)+
scale_colour_viridis_c()+
labs(color="avg.exp.scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+
facet_grid(~Group, scale="free_x",space = "free")+
theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(face="italic",color="black",angle=90,vjust=0.5,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("SASPmarkerGene.pdf",width=22,height=4)
print(g)
dev.off()
write.table(result,file="SASPmarkerGeneExprAcrossLayers.txt",row.names=F,sep="\t",quote=F)

#######################################################################################
##############         Expression of DNA damage associated genes          #############
#######################################################################################
#Figure S5
setwd("/data2/deng/Aging/Ex/AllEx/C5Feature/")
ExTmpExpr <- AverageExpression(ExAll.integrated,group.by="Layer")[["RNA"]]
DNADamage=read.table("/Projects/deng/Aging/Ex/AllEx/DNArepair/DNADamage8Kim.txt",header=T,row.names=1,sep="\t")
DNADamage$GeneSymobl=rownames(DNADamage)


tmpCluster=DNADamage[DNADamage$Group %in% c("Mismatch excision repair (MMR)"),]
tmpClusterExpr=ExTmpExpr[intersect(rownames(ExTmpExpr),tmpCluster$GeneSymobl),]
pdf("DDRHeatmap/Mismatch excision repair (MMR).pdf",width=6,height=4)
pheatmap(t(tmpClusterExpr),clustering_method="ward.D2",scale="column",color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

DNADamageExpr=ExTmpExpr[intersect(rownames(ExTmpExpr),rownames(DNADamage)),]
DNADamageExpr=data.frame(DNADamageExpr,check.names=F)
DNADamageExpr$GeneSymobl=rownames(DNADamageExpr)
DNADamageExprAddInfo=merge(DNADamage,DNADamageExpr,by="GeneSymobl")
write.table(DNADamageExprAddInfo,file="DDRHeatmap/Mismatch excision repair (MMR).txt",row.names=F,sep="\t",quote=F)



#######marker genes for C3 and C5#################################
setwd("/data2/deng/Aging/Ex/AllEx/C5Feature")
ExAll.integrated=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
ExAll.integrated=subset(ExAll.integrated,seurat_clusters%in%c(0:15))
DefaultAssay(ExAll.integrated)="RNA"
C5Marker=FindMarkers(ExAll.integrated,test.use ="MAST",latent.vars=c("Source","Gender","Age","Statues"),ident.1=5)
C5Marker=C5Marker[C5Marker$p_val<0.01,]
write.table(C5Marker,file="C5Marker_MAST.txt",sep="\t",quote=F)

C5MarkerTmp=C5Marker[abs(C5Marker$avg_log2FC)>1,]
C5MarkerTmp$Pattern=ifelse(C5MarkerTmp$avg_log2FC>0,"UpInLS","DnInLS")
write.table(C5MarkerTmp,file="C5MarkerFC1_MAST.txt",sep="\t",quote=F)


C3Marker=FindMarkers(ExAll.integrated,ident.1=3,ident.2=c(0,2,8),,test.use ="MAST",latent.vars=c("Source","Gender","Age","Statues"))
C3Marker=C3Marker[C3Marker$p_val<0.01,]
dim(C3Marker)
C3Marker$Pattern=ifelse(C3Marker$avg_log2FC>0,"UpInES","DnInES")
write.table(C3Marker,file="C3MarkerSuperiorLayer_MAST.txt",sep="\t",quote=F)


C3C5Marker=FindMarkers(ExAll.integrated,ident.1=c(3,5),,test.use ="MAST",latent.vars=c("Source","Gender","Age","Statues"))
C3C5Marker=C3C5Marker[C3C5Marker$p_val<0.01,]
C3C5Marker$Pattern=ifelse(C3C5Marker$avg_log2FC>0,"UpInSene","DnInSene")
dim(C3C5Marker)
write.table(C3C5Marker,file="/data2/deng/Aging/Ex/AllEx/C5Feature/C3C5Marker_MAST.txt",sep="\t",quote=F)



#######################################################################################
##############         Effect of TFs on C5 markers:           #############
#######################################################################################
See transcriptionFactorAnalysis_modified.code.R


#######################################################################################
##############         AD associated genes in C5 (CNV)           ###########
#######################################################################################
See InferCNV.code.R



#######################################################################################
##############         Expression of the AD associated genes in C5           ##########
#######################################################################################
library(dplyr)
library(tibble)
library(fgsea)

setwd("/data2/deng/Aging/Ex/AllEx/InferCNV/ADRiskPlot/")

ExSmall <- ExAll.integrated[, sample(colnames(ExAll.integrated), size = 10000, replace=F)]
t=DotPlot(ExSmall,features=rownames(ExSmall))
ExExpr8DotPlot=na.omit(t$data)

CNVGainLoss=read.table("/data2/deng/Aging/Ex/AllEx/InferCNV/ADRiskPlot/ExprFrame4GraphAddBandSum.txt",header=T)
CNVGainLoss=unique(CNVGainLoss[,c("Band","CNVSum")])
CNVGainLoss$Pattern=ifelse(CNVGainLoss$CNVSum>0.1,"Gain",ifelse(CNVGainLoss$CNVSum<=-0.1,"Loss","Normal"))

ADriskGene=read.table("/data2/deng/Alzheimer/GWAS/ADLociBandInfor.txt",header=T)
ADRiskGenePattern=merge(ADriskGene,CNVGainLoss,by.x="ChLocation",by.y="Band")
ADRiskGenePattern=ADRiskGenePattern[ADRiskGenePattern$Pattern%in%c("Gain","Loss"),]

ADriskGeneList=split(ADRiskGenePattern$Symbol,ADRiskGenePattern$Pattern)

cluster="5"
clusterCell<- ExExpr8DotPlot %>% dplyr::filter(id == cluster) %>% arrange(desc(avg.exp.scaled)) %>% dplyr::select(features.plot, avg.exp.scaled)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(ADriskGeneList, stats = ranks,eps=0)
ranks=na.omit(ranks)
fgseaRes=fgseaRes[order(fgseaRes$NES,decreasing=T),]
pdf("/data2/deng/Aging/Ex/AllEx/InferCNV/ADRiskPlot/ADRiskGene_Loss.pdf",height=4)
plotEnrichment(ADriskGeneList[["Loss"]],ranks)
dev.off()


ADGene=c("BIN1","HLA-A","HLA-B","HLA-C","HLA-DQB1","HLA-DRB1","NYAP1","CLU","BCKDK","DOC2A","KAT8","IL34","KLF16","APOE") #,"FOXF1","APOE","ADAMTS4","KLF16","MADD","ABI3","SHARPIN","PTK2B",
t=DotPlot(ExAll.integrated,features=ADGene,group.by="Layer")
result=t$data
g=ggplot(result, aes(features.plot, id, size= pct.exp,color=avg.exp.scaled)) +
geom_point(alpha = 0.9)+
scale_colour_viridis_c()+
labs(color="avg.exp.scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+
theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(face="italic",color="black",angle=90,vjust=0.5,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("ADGeneExprAcrossLayers.pdf",width=6,height=3)
print(g)
dev.off()
write.table(result,file="ADGeneExprAcrossLayers.txt",row.names=F,sep="\t",quote=F)

tiff("ADGeneExpr.tiff",width=2000,height=2000)
FeaturePlot(ExAll.integrated,features=ADGene)&NoLegend()&NoAxes()
dev.off()


LayerPalette=c("Orange","Firebrick3","#DBE419","#74CF56","#21D66A","#3DBC75","#238B8E","#2E718E")
pdf("ADGeneExprVlnPlots.pdf",width=20,height=16)
VlnPlot(ExAll.integrated,features=ADGene,group.by="Layer",pt.size=0.01,cols=LayerPalette,raster =TRUE)&NoLegend()&NoAxes()&theme(title = element_text(size=rel(1.0),face="italic"))
dev.off()




#################################################################################################
##############        cell number variation for ES and LS across phenotypes           ##########
#################################################################################################
See cellNumberVariation.code.R



#######################################################################################
############################### Deg between AD and Normal in cell cycling ##############
#######################################################################################
setwd("/Projects/deng/Aging/Ex/AllEx/DEGBetADNorInCycing")
ExAll.integrated=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
DefaultAssay(ExAll.integrated)="RNA"
ExAll.integrated=subset(ExAll.integrated,idents=c(0:15))
ExAllTmp=subset(ExAll.integrated,Age>60) #remove Nagy et al,
Idents(ExAllTmp)="Statues"
table(ExAllTmp$Source)
#Source=subset(ExAllTmp,Source=="Yang") #Lau Mathy 
for(i in c(0:15)){
  SubCluster=subset(ExAllTmp,seurat_clusters==i)
  ADDeg=FindMarkers(SubCluster,ident.1="Alzheimer",ident.2="Control",test.use ="MAST",latent.vars=c("Source","Gender","Age"))
  ADDeg=ADDeg[ADDeg$p_val<0.01,]
  write.table(ADDeg,file=paste0("/Projects/deng/Aging/Ex/AllEx/DEGBetADNorInCycing/ADDeg4Old/AllSource_Mast/C",i,".txt",sep=""),sep="\t",quote=F)
}


ExAllTmp=subset(ExAll.integrated,Age>60) #remove Nagy et al,
Idents(ExAllTmp)="Statues"
table(ExAllTmp$Statues,ExAllTmp$Layer)
#Source=subset(ExAllTmp,Source=="Yang") #Lau Mathy 
for(l in unique(ExAllTmp$Layer)){
  SubCluster=subset(ExAllTmp,Layer==l)
  ADDeg=FindMarkers(SubCluster,ident.1="Alzheimer",ident.2="Control",test.use ="MAST",latent.vars=c("Source","Gender","Age"))
  ADDeg=ADDeg[ADDeg$p_val<0.01,]
  lname=stringr::str_replace(l,"/","_")
  write.table(ADDeg,file=paste0("/Projects/deng/Aging/Ex/AllEx/DEGBetADNorInCycing/ADDeg4Old/AllSource_Layer_Mast/",lname,".txt",sep=""),sep="\t",quote=F)
}

#------------------------count the DEGs number in each cluster between AD and Ctrl----------------------------

ADDeg4Old=read.table("D:/Aging/Ex/AllEx/DEGBetADNorInCycing/ADDeg4Old/ADDegByMastInAllSource_Layer.txt",header=T)
ADDeg4OldCount=data.frame(table(ADDeg4Old$Cluster,ADDeg4Old$Pattern))
colnames(ADDeg4OldCount)=c("Cluster","Pattern","Count")
#ADDeg4OldCount=ADDeg4OldCount[ADDeg4OldCount$Cluster%in%c("C0","C3","C5"),]
ClusterOrder=factor(ADDeg4OldCount$Cluster,levels=c("ES","LS","L2_3","L4_5","L4_6","L5","L5_6","L6"))
g=ggplot(ADDeg4OldCount, aes(ClusterOrder, Count, fill=factor(Pattern,levels=c("Up","Down")))) +
  geom_bar(stat="identity",position=position_dodge()) +
  #scale_y_continuous(expand=c(0,0))+
  theme_bw()+
  scale_fill_manual(values=c("Violet","CornflowerBlue"))+
  guides(fill = guide_legend(title = "Pattern", title.position = "top"),col = guide_legend(nrow = 1))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
pdf("D:/Aging/Ex/AllEx/DEGBetADNorInCycing/DEGbetADInOld_Layer.pdf",height=4,width=6)
print(g)
dev.off()

write.table(ADDeg4OldCount,file="D:/Aging/Ex/AllEx/DEGBetADNorInCycing/DEGNumberbetADInOld_Layer.txt",sep="\t",quote=F)


UpInES=ADDeg4Old[ADDeg4Old$Cluster=="ES"&ADDeg4Old$Pattern=="Up","Gene"]
DnInES=ADDeg4Old[ADDeg4Old$Cluster=="ES"&ADDeg4Old$Pattern=="Down","Gene"]
write.table(UpInES,file="D:/Aging/Ex/AllEx/DEGBetADNorInCycing/ADDegByMast/UpDEGbetADAndCtrlInESByMast.txt",row.names=F,col.names=F,quote=F)
write.table(DnInES,file="D:/Aging/Ex/AllEx/DEGBetADNorInCycing/ADDegByMast/DownDEGbetADAndCtrlInESByMast.txt",row.names=F,col.names=F,quote=F)

UpInLS=ADDeg4Old[ADDeg4Old$Cluster=="LS"&ADDeg4Old$Pattern=="Up","Gene"]
DnInLS=ADDeg4Old[ADDeg4Old$Cluster=="LS"&ADDeg4Old$Pattern=="Down","Gene"]
write.table(UpInLS,file="D:/Aging/Ex/AllEx/DEGBetADNorInCycing/ADDegByMast/UpDEGbetADAndCtrlInLSByMast.txt",row.names=F,col.names=F,quote=F)
write.table(DnInLS,file="D:/Aging/Ex/AllEx/DEGBetADNorInCycing/ADDegByMast/DownDEGbetADAndCtrlInLSByMast.txt",row.names=F,col.names=F,quote=F)

ADDEGList=list(UpInADES=UpInES,DnInADES=DnInES,UpInADLS=UpInLS,DnInADLS=DnInLS)


LSMarker=read.table("D:/Aging/Ex/Manuscript/Plos/2.Revised/2.Graph/1.RawGraph/C5Feature/C5MarkerFC1_MAST.txt",header=T)
LSMarkerList=split(rownames(LSMarker),LSMarker$Pattern)
names(LSMarkerList)=c("InhibitedMarkersInLS","ActivatedMarkersInLS") #"DnInLS" "UpInLS"
AllList=c(ADDEGList,LSMarkerList)
pdf("D:/Aging/Ex/Manuscript/Plos/2.Revised/2.Graph/1.RawGraph/ADDEG/ADESAndLS.pdf",width=7,height=5)
upset(fromList(AllList),nsets=6,
       queries = list(
         list(query = intersects, params = list("UpInADES", "ActivatedMarkersInLS"), color = "Red", active = T),
         list(query = intersects, params = list("DnInADES", "InhibitedMarkersInLS"), color = "Blue", active = T)
       )
     )
dev.off()

1-phyper(length(intersect(AllList$UpInADES,AllList$ActivatedMarkersInLS))-1,length(AllList$UpInADES),25834-length(AllList$UpInADES),length(AllList$ActivatedMarkersInLS))
1-phyper(length(intersect(AllList$DnInADES,AllList$InhibitedMarkersInLS))-1,length(AllList$DnInADES),25834-length(AllList$DnInADES),length(AllList$InhibitedMarkersInLS))



ExAll.integrated=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
DefaultAssay(ExAll.integrated)="RNA"
ExAll.integrated=subset(ExAll.integrated,idents=c(0:15))
ExAllTmp=subset(ExAll.integrated,Age>60) 
ExAllTmp$Group=paste0(ExAllTmp$Layer,"_",ExAllTmp$Statues,sep="")
pdf("ADGeneExpr.pdf",width=10,height=4)
DotPlot(subset(ExAllTmp,Layer%in%c("ES","LS")),features=c("AGAP3","ELK1","PDZD4","PSD","GUK1","SYN1","YPEL3","MAPK8IP1","FTL","RASD1","C9orf16","CHCHD10","PEA15","NEUROD2"),group.by="Group")&RotatedAxis()
dev.off()


pdf("ADGeneExpr.pdf",width=4,height=25)
FeaturePlot(ExAllTmp,raster=TRUE,features=c("AGAP3","ELK1","PDZD4","PSD","GUK1","SYN1","YPEL3","MAPK8IP1","FTL","RASD1","C9orf16","CHCHD10","PEA15","NEUROD2"),split.by="Statues")&NoLegend()&NoAxes()
dev.off()


#-----expression of immune associated genes---------------
setwd("/data2/deng/Aging/Ex/AllEx/C5Feature")
ExTmp=subset(ExAll.integrated,idents=c(0:6)) #only show the top 7 clusters
table(ExTmp$seurat_clusters) #validate the custers
DefaultAssay(ExTmp)="RNA"
ExTmp$Group=paste0(ExTmp$Layer,"_",ExTmp$Statues,sep="")
ImmuneAssociateGene=read.table("ImmuneAssociateGene.txt",header=T,sep="\t")
t=DotPlot(ExTmp,features=ImmuneAssociateGene$Symbol)
result=merge(t$data,ImmuneAssociateGene,by.x="features.plot",by.y="Symbol")
C0=result[result$id=="0_Control",] #order by their expression in cluster 0
C0=C0[order(C0$Group,C0$pct.exp),]
geneOrder=C0$features.plot
g=ggplot(result, aes(factor(features.plot,levels=unique(rev(geneOrder))), id, size= pct.exp,color=avg.exp.scaled)) +geom_point(alpha = 0.9)+
scale_colour_viridis_c()+
labs(color="avg.exp.scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+
facet_grid(~Group, scale="free_x",space = "free")+
theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(angle=90,face = "italic",vjust=1,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("ImmuneAssociateGeneSplitVersion.pdf",width=15,height=4)
print(g)
dev.off()


#figure 4B
#-----expression of immune associated genes---------------
immuneGene=c("BIN1","COX7C","HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPA2","HLA-DPA3","HLA-DPB1","HLA-DPB2","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DQB2","HLA-DQB3","HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DRB6","HLA-DRB9","HLA-E","HLA-F","HLA-G","HLA-H","HLA-J","HLA-K","HLA-L","HLA-N","HLA-P","HLA-S","HLA-T","HLA-U","HLA-V","HLA-W","HLA-Z","NYAP1","ZCWPW1","ECHDC3","RABEP1","IL34","SCIMP","KLF16","ABCA7","APOE","SIGLEC11","KAT8")
ExAllTmp=subset(ExAll.integrated,idents=c(0,3,5))
Idents(ExAllTmp)=factor(ExAllTmp$seurat_clusters,levels=c(0,3,5))
ExAllTmp$Statues=factor(ExAllTmp$Statues,levels=c("Control","Alzheimer"))
tiff("immune.Gene_VlnPlot.tiff",width=1800,height=1800)
VlnPlot(ExAllTmp,features=immuneGene,pt.size=0.01,split.by="Statues")&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))+theme(legend.position="right")
dev.off()
counts <- ExAllTmp@assays$RNA@counts
counts=data.frame(counts,check.names=F)
immuneGeneExpr=counts[intersect(immuneGene,rownames(ExAllTmp)),]
all(colnames(ExAllTmp)==colnames(immuneGeneExpr))
immuneGeneExprInfo=cbind(ExAllTmp@meta.data[,c("seurat_clusters","Gender")],t(immuneGeneExpr))
write.table(immuneGeneExprInfo,file="immune.Gene.expr.txt",sep="\t",quote=F)
