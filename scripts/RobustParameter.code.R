library(Seurat)
library(harmony) #0.1.1
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(ggalluvial)
library(clustree)
packageVersion('Seurat') #4.2.0

setwd("/data2/deng/Aging/Ex/AllEx")
ExAll.integrated=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
ExAll.integrated$rpcaCluster=ExAll.integrated$seurat_clusters
ExAll.integrated=subset(ExAll.integrated,idents=c(0:15))

####### validate by different integration methods (hormony)#############################

DefaultAssay(ExAll.integrated)="RNA"
ExAll.integrated <- NormalizeData(ExAll.integrated)
ExAll.integrated <- FindVariableFeatures(ExAll.integrated)
ExAll.integrated <- ScaleData(ExAll.integrated)
ExAll.integrated <- RunPCA(ExAll.integrated)
ExAll.hormony.integrated <- RunHarmony(ExAll.integrated, group.by.vars=c("Source"))
ExAll.hormony.integrated <- RunUMAP(ExAll.hormony.integrated, reduction = "harmony", dims = 1:30)
ExAll.hormony.integrated <- FindNeighbors(ExAll.hormony.integrated, reduction = "harmony", dims = 1:30)
ExAll.hormony.integrated <- FindClusters(ExAll.hormony.integrated, resolution = 0.5) 

hormonyColorPalette=colorRampPalette(brewer.pal(8, "Set3"))(length(unique(ExAll.hormony.integrated$seurat_clusters))) 
hormonyColorPalette[4]="Orange"
hormonyColorPalette[6]="Firebrick3"
pdf("hormony/ExAll.hormony.integratedLabel.pdf",width=4,height=3.5)
DimPlot(ExAll.hormony.integrated, raster=TRUE,reduction = "tsne",label=T,label.size = 6,cols =hormonyColorPalette)&NoLegend()&theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", linewidth=1.5, linetype="solid"))
dev.off()

rpcaColorPalette=colorRampPalette(brewer.pal(8, "Set3"))(length(unique(ExAll.hormony.integrated$rpcaCluster))) 
rpcaColorPalette[4]="Orange"
rpcaColorPalette[6]="Firebrick3"
rpcaColorPalette[11]="LightCyan"
pdf("hormony/ExAll.hormony.integrated.rpcaCluster.pdf",width=4,height=3.5)
DimPlot(ExAll.hormony.integrated, raster=TRUE,reduction = "tsne",group.by="rpcaCluster",label=T,label.size = 6,cols =rpcaColorPalette)&NoLegend()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", linewidth=1.5, linetype="solid"))
dev.off()

saveRDS(ExAll.hormony.integrated,file="hormony/ExAll.hormony.integrated.rds")

ExAll.hormony.integrated=readRDS("/data2/deng/Aging/Ex/AllEx/hormony/ExAll.hormony.integrated.rds")
Link=data.frame(table(ExAll.hormony.integrated$seurat_clusters,ExAll.hormony.integrated$rpcaCluster))
colnames(Link)=c("HormonyCluster","rPCACluster","SharedCells")
Link=Link[Link$rPCACluster%in%c(0:15),]
#https://r-charts.com/flow/ggalluvial/#google_vignette
p<-ggplot(data = Link,aes(axis1 = HormonyCluster, axis2 = rPCACluster, y = SharedCells)) +
  geom_alluvium(aes(fill = rPCACluster)) +
  geom_stratum() +
  scale_fill_manual(values = rpcaColorPalette)+
  geom_text(stat = "stratum",aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Survey", "Response"),
                   expand = c(0.15, 0.05)) +
  theme_void()
pdf("hormony/LinkBetweenHormonyAndrPCA.pdf",width=6,height=5)
print(p)
dev.off()


####### check the influence of UMAP #############################
ExAll.integrated=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
ExAll.integrated$rpcaCluster=ExAll.integrated$seurat_clusters
DefaultAssay(ExAll.integrated)="integrated"
ExAll.integrated <- RunUMAP(ExAll.integrated, dims = 1:50)
ExAll.integrated<- FindNeighbors(ExAll.integrated, reduction = "pca", dims = 1:30)
ExAll.integrated <- FindClusters(ExAll.integrated, resolution = 0.5) 

rpcaColorPalette=colorRampPalette(brewer.pal(8, "Set3"))(length(unique(ExAll.integrated$rpcaCluster))) 
rpcaColorPalette[4]="Orange"
rpcaColorPalette[6]="Firebrick3"
rpcaColorPalette[11]="LightCyan"
pdf("umap/ExAll.umap.integrated.tsneCluster.pdf",width=4,height=3.5)
DimPlot(ExAll.integrated, raster=TRUE,reduction = "umap",group.by="rpcaCluster",label=T,label.size = 6,cols =rpcaColorPalette)&NoLegend()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", linewidth=1.5, linetype="solid"))
dev.off()

umapColorPalette=colorRampPalette(brewer.pal(8, "Set3"))(length(unique(ExAll.integrated$seurat_clusters))) 
umapColorPalette[3]="Orange"
umapColorPalette[8]="Firebrick3"
pdf("umap/ExAll.umap.integrated.umapCluster.pdf",width=4,height=3.5)
DimPlot(ExAll.integrated, raster=TRUE,reduction = "umap",label=T,label.size = 6,cols =umapColorPalette)&NoLegend()&theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()


ExAll.integrated=subset(ExAll.integrated,rpcaCluster%in%c(0:15))
ExAll.integrated[["Layer"]]=NA
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(0,2,8),]="L2/3"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(1,6,7,9,13),]="L5/6"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(4,10),]="L4/5"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(11,14,15),]="L5"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(12),]="L6"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(3),]="ES"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(5),]="LS"
ExAll.integrated$Layer=factor(ExAll.integrated$Layer,levels=c("ES","LS","L2/3","L4/5","L5","L5/6","L6"))

LayerPalette=c("Orange","Firebrick3","#DBE419","#74CF56","#238B8E","#2E718E")
pdf("/data2/deng/Aging/Ex/AllEx/C5Feature/ExAll.Layer.pdf",width=4,height=3.5)
DimPlot(ExAll.integrated, raster=TRUE,reduction = "tsne",group.by="Layer",label=T,label.size = 6,cols =LayerPalette)&NoLegend()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

pdf("/data2/deng/Aging/Ex/AllEx/C5Feature/ExAll.Layer.umap.pdf",width=4,height=3.5)
DimPlot(ExAll.integrated, raster=TRUE,reduction = "umap",group.by="Layer",label=T,label.size = 6,cols =LayerPalette)&NoLegend()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

####### check the influence of resolution #############################

ExAll.integrated=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
ExAll.integratedTmp <- FindClusters(ExAll.integrated, resolution = c(seq(0,0.9,0.1),seq(1,2,0.2)))
clus.tree.out <- clustree(ExAll.integratedTmp,edge_width = 0.5)
pdf(file = "Resoution/ExAll.integrated.Resolution.Tree.pdf", width = 15, height = 15)
print(clus.tree.out)
dev.off()

ExAll.integratedTmp=FindClusters(ExAll.integrated, resolution = 1)
colorPalette=colorRampPalette(brewer.pal(8, "Set3"))(length(unique(ExAll.integratedTmp$seurat_clusters))) 
pdf("Resoution/ExAll.integrated.R1.pdf",width=4,height=3.5)
DimPlot(ExAll.integratedTmp, raster=TRUE,reduction = "tsne",label=T,label.size = 6,cols =colorPalette)&NoLegend()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", linewidth=1.5, linetype="solid"))
dev.off()

ExAll.integratedTmp=FindClusters(ExAll.integrated, resolution = 2)
colorPalette=colorRampPalette(brewer.pal(8, "Set3"))(length(unique(ExAll.integratedTmp$seurat_clusters))) 
pdf("Resoution/ExAll.integrated.R2.pdf",width=4,height=3.5)
DimPlot(ExAll.integratedTmp, raster=TRUE,reduction = "tsne",label=T,label.size = 6,cols =colorPalette)&NoLegend()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", linewidth=1.5, linetype="solid"))
dev.off()



