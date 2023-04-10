library(Seurat)
library(ggplot2)
library(clustree)

#Preprocessing for snRNA from Mathys et,al can be found in snRNAPreprocessing.code.R 
#cell type identification
MathyAllCtrl=readRDS("/Projects/deng/Aging/Ex/MathyEx/MathyAllCtrl.rds")
Idents(MathyAllCtrl)=factor(MathyAllCtrl$OldCellType,levels=c("Ex","Oli","In","Ast","Opc","Mic","End","Per"))

#---Figure 1A---
tiff("MathyAllCtrl.CellType.tiff",width=280,height=200)
DimPlot(MathyAllCtrl, reduction = "tsne",label=TRUE)&NoLegend()&theme(plot.title=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

#---Figure 1B---
tiff("MathyAllCtrl.CellTypeMarker.tiff",width=1000,height=150)
FeaturePlot(MathyAllCtrl, reduction = "tsne",col=c("lightgrey",viridis(10)),features=c("NRGN","MBP","GAD1","AQP4","VCAN","CSF1R","FLT1"),ncol=7)&NoLegend()&NoAxes()&theme(plot.title=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())
dev.off()

####################################################################################
############ Differentiated cell type isolation and visualiation ############
####################################################################################

MathyCtrlDifferentiatedCells=subset(MathyAllCtrl,OldCellType %in% c("Ex","In","Ast","Oli"))
MathyCtrlDifferentiatedCells <- RunTSNE(MathyCtrlDifferentiatedCells, dims = 1:30)
MathyCtrlDifferentiatedCells<- FindNeighbors(MathyCtrlDifferentiatedCells, reduction = "pca", dims = 1:30)
MathyCtrlDifferentiatedCells<- FindClusters(MathyCtrlDifferentiatedCells, resolution = 0.5)

MathyTmp=MathyCtrlDifferentiatedCells
#---Figure 1C---
Idents(MathyTmp)=factor(MathyTmp$OldCellType,levels=c("Ex","Oli","In","Ast"))
tiff("MathyCtrlDifferentiatedCells.cellType.tiff",width=230,height=200)
DimPlot(MathyTmp, reduction = "tsne")&NoLegend()&theme(plot.title=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

Idents(MathyCtrlDifferentiatedCells)=factor(MathyCtrlDifferentiatedCells$seurat_clusters,levels=c(0:19))
tiff("MathyCtrlDifferentiatedCells.cellClusters.tiff",width=230,height=200)
DimPlot(MathyCtrlDifferentiatedCells, reduction = "tsne")&NoLegend()&theme(plot.title=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

####################################################################################
############ cell cycle phase score calculation for Differentiated cells  ##########
####################################################################################
MathyCtrlDifferentiatedCells=readRDS("/Projects/deng/Aging/Ex/MathyEx/MathyCtrlDifferentiatedCells.rds")

#cell phase score for each cell cluster from differentiated cell type
CCGene=read.table("/Projects/deng/Aging/Ex/cellCycleGene/CCGeneCombine.txt",header=T,sep="\t")

DefaultAssay(MathyCtrlDifferentiatedCells)="RNA"
genes=intersect(rownames(MathyCtrlDifferentiatedCells),CCGene$NAME)
length(genes)
CCGene=CCGene[CCGene$NAME%in%genes,]
CCGenelist=split(CCGene$NAME,CCGene$PHASE)
MathyCtrlDifferentiatedCells <- AddModuleScore(
  object = MathyCtrlDifferentiatedCells,
  features = CCGenelist,
  name = names(CCGenelist),
  ctrl = 100,
)
cluster=c("2","3","4","5","6","8","10","13","14","16","17","18","0","1","7","11","12","15","19","9")
#---Figure 1E---
MathyCtrlDifferentiatedCells$seurat_clusters=factor(MathyCtrlDifferentiatedCells$seurat_clusters,levels=cluster)
pdf("MathyCtrlDifferentiatedCells.cellCycleScore.pdf",width=6)
VlnPlot(MathyCtrlDifferentiatedCells,c("G1.S1","S.phase5","G22","G2.M3","M.G14"),ncol=1,pt.size=0.1,group.by="seurat_clusters")&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

CycleScore=MathyCtrlDifferentiatedCells@meta.data[,c("seurat_clusters","G1.S1","S.phase5","G22","G2.M3","M.G14")]
cluster=c(2:6,8,10,13,14,16,17,18,0,1,7,11,12,15,19,9)
PvalueBetCluster=matrix(data = NA, nrow = length(cluster), ncol = length(cluster), dimnames = list(cluster, cluster))
for(i in 1:length(cluster)){
  for( j in 1:length(cluster)){
    C1_score=CycleScore[CycleScore$seurat_clusters==cluster[i],"M.G14"]
    C2_score=CycleScore[CycleScore$seurat_clusters==cluster[j],"M.G14"]
    delta=mean(C1_score)-mean(C2_score) #the row represent the current cluster
    pval=wilcox.test(C1_score,C2_score)$p.value
    PvalueBetCluster[i,j]=delta*(-log10(pval))/abs(delta)
    if(i==j||delta<0){
      PvalueBetCluster[i,j]=0;
    }
  }
}
PvalueBetCluster=data.frame(PvalueBetCluster,check.names=F)
PvalueBetCluster$Cluster=as.character(rownames(PvalueBetCluster))
PvalueBetCluster.df=reshape2::melt(PvalueBetCluster,id=21)
PvalueBetCluster.df$Cluster=factor(PvalueBetCluster.df$Cluster,levels=cluster)
t=ggplot(PvalueBetCluster.df, aes(x=Cluster, y=value, color=Cluster))+
geom_boxplot(outlier.shape = NA,outlier.colour=NA)+labs(title="",x="", y = "")+
geom_jitter(size=1)+
theme_classic()+
theme(legend.position="none")+
theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
pdf("/Projects/deng/Aging/Ex/MathyEx/MG1_PvalueDis.pdf",height=3,width=6)
print(t)
dev.off()

####################################################################################
############   Mitosis cell type isolation and visualiation   ############
####################################################################################
setwd("/Projects/deng/Aging/Ex/MathyEx/")
MathyAllCtrl=readRDS("MathyAllCtrl.rds")
MathyCtrlMitosisCells=subset(MathyAllCtrl,OldCellType %in% c("Mic","Opc","Per","End"))
MathyCtrlMitosisCells <- RunTSNE(MathyCtrlMitosisCells, dims = 1:30)
MathyCtrlMitosisCells<- FindNeighbors(MathyCtrlMitosisCells, reduction = "pca", dims = 1:30)
MathyCtrlMitosisCells<- FindClusters(MathyCtrlMitosisCells, resolution = 0.5)

saveRDS(MathyCtrlMitosisCells,file="MathyCtrlMitosisCells.rds")
#---Figure S1A---
tiff("MathyCtrlMitosisCells.cellType.tiff",width=230,height=200)
DimPlot(MathyCtrlMitosisCells, group.by="OldCellType",reduction = "tsne")&NoLegend()&theme(plot.title=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

#---Figure S1B---
tiff("MathyCtrlMitosisCells.cellClusters.tiff",width=230,height=200)
DimPlot(MathyCtrlMitosisCells,group.by="seurat_clusters", reduction = "tsne")&NoLegend()&theme(plot.title=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()


####################################################################################
############ cell cycle phase score calculation for Differentiated cells  ##########
####################################################################################
setwd("/Projects/deng/Aging/Ex/MathyEx/")
MathyCtrlMitosisCells=readRDS("MathyCtrlMitosisCells.rds")
cluster=c("0","2","4","1","3","6","5")
#---Figure S1C---
DefaultAssay(MathyCtrlMitosisCells)="RNA"
genes=intersect(rownames(MathyCtrlMitosisCells),CCGene$NAME)
length(genes)
CCGene=CCGene[CCGene$NAME%in%genes,]
CCGenelist=split(CCGene$NAME,CCGene$PHASE)
MathyCtrlMitosisCells <- AddModuleScore(
  object = MathyCtrlMitosisCells,
  features = CCGenelist,
  name = names(CCGenelist),
  ctrl = 100,
)
Idents(MathyCtrlMitosisCells)=factor(MathyCtrlMitosisCells$seurat_clusters,levels=cluster)
pdf("MathyMitosisCells.cellCycleScore.pdf",width=4)
VlnPlot(MathyCtrlMitosisCells,c("G1.S1","S.phase5","G22","G2.M3","M.G14"),ncol=1,pt.size=0)&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

CycleScore=MathyCtrlMitosisCells@meta.data[,c("seurat_clusters","G1.S1","S.phase5","G22","G2.M3","M.G14")]
PvalueBetCluster=matrix(data = NA, nrow = length(cluster), ncol = length(cluster), dimnames = list(cluster, cluster))
for(i in 1:length(cluster)){
  for( j in 1:length(cluster)){
    C1_score=CycleScore[CycleScore$seurat_clusters==cluster[i],"G22"]
    C2_score=CycleScore[CycleScore$seurat_clusters==cluster[j],"G22"]
    delta=mean(C1_score)-mean(C2_score) #the row represent the current cluster
    pval=wilcox.test(C1_score,C2_score)$p.value
    PvalueBetCluster[i,j]=delta*(-log10(pval))/abs(delta)
    if(i==j||delta<0){
      PvalueBetCluster[i,j]=0;
    }
  }
}
PvalueBetCluster=data.frame(PvalueBetCluster,check.names=F)
PvalueBetCluster$Cluster=as.character(rownames(PvalueBetCluster))
PvalueBetCluster.df=reshape2::melt(PvalueBetCluster,id=8)
PvalueBetCluster.df$Cluster=factor(PvalueBetCluster.df$Cluster,levels=cluster)
t=ggplot(PvalueBetCluster.df, aes(x=Cluster, y=value, color=Cluster))+
geom_boxplot(outlier.shape = NA,outlier.colour=NA)+labs(title="",x="", y = "")+
geom_jitter(size=1)+
theme_classic()+
theme(legend.position="none")+
theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
pdf("MathyCtrlMitosisCells_Phase_PvalueDis/MathyCtrlMitosisCells_G2_PvalueDis.pdf",height=2,width=4)
print(t)
dev.off()

