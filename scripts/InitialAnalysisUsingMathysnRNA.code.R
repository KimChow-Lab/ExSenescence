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
Idents(MathyCtrlDifferentiatedCells)=factor(Idents(MathyCtrlDifferentiatedCells),levels=cluster)
#---Figure 1E---
pdf("MathyCtrlDifferentiatedCells.cellCycleScore.pdf",width=6)
VlnPlot(MathyCtrlDifferentiatedCells,c("G1.S1","S.phase5","G22","G2.M3","M.G14"),ncol=1,pt.size=0)&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()











####################################################################################
############   Mitosis cell type isolation and visualiation   ############
####################################################################################

MathyAllCtrl=readRDS("/Projects/deng/Aging/Ex/MathyEx/MathyAllCtrl.rds")
MathyCtrlMitosisCells=subset(MathyAllCtrl,OldCellType %in% c("Mic","Opc","Per","End"))
MathyCtrlMitosisCells <- RunTSNE(MathyCtrlMitosisCells, dims = 1:30)
MathyCtrlMitosisCells<- FindNeighbors(MathyCtrlMitosisCells, reduction = "pca", dims = 1:30)
MathyCtrlMitosisCells<- FindClusters(MathyCtrlMitosisCells, resolution = 0.5)

#---Figure S1A---
Idents(MathyTmp)=factor(MathyTmp$OldCellType,levels=c("Ex","Oli","In","Ast"))
tiff("MathyCtrlDifferentiatedCells.cellType.tiff",width=230,height=200)
DimPlot(MathyTmp, reduction = "tsne")&NoLegend()&theme(plot.title=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

#---Figure S1B---
Idents(MathyCtrlDifferentiatedCells)=factor(MathyCtrlDifferentiatedCells$seurat_clusters,levels=c(0:19))
tiff("MathyCtrlDifferentiatedCells.cellClusters.tiff",width=230,height=200)
DimPlot(MathyCtrlDifferentiatedCells, reduction = "tsne")&NoLegend()&theme(plot.title=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()


####################################################################################
############ cell cycle phase score calculation for Differentiated cells  ##########
####################################################################################

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
cluster=c("0","2","4","1","3","6","5")
Idents(MathyCtrlMitosisCells)=factor(Idents(MathyCtrlMitosisCells),levels=cluster)
pdf("MathyMitosisCells.cellCycleScore.pdf",width=4)
VlnPlot(MathyCtrlMitosisCells,c("G1.S1","S.phase5","G22","G2.M3","M.G14"),ncol=1,pt.size=0)&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

