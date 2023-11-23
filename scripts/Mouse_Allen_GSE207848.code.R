library(Seurat)
library(ggplot2)
library(homologene)
library(RColorBrewer)
library(scCustomize)

setwd("/Projects/deng/Aging/Ex/Mouse_Allen_GSE207848")

GSM6321077.count=Read10X_h5("rawData/GSM6321077_PFC_90wk_1_matrix.h5")
GSM6321077 <- CreateSeuratObject(counts = GSM6321077.count, project = "GSM6321077")
GSM6321078.count=Read10X_h5("rawData/GSM6321078_PFC_90wk_2_matrix.h5")
GSM6321078 <- CreateSeuratObject(counts = GSM6321078.count, project = "GSM6321078")
GSM6321079.count=Read10X_h5("rawData/GSM6321079_PFC_90wk_3_matrix.h5")
GSM6321079 <- CreateSeuratObject(counts = GSM6321079.count, project = "GSM6321079")
GSM6321080.count=Read10X_h5("rawData/GSM6321080_PFC_90wk_4_matrix.h5")
GSM6321080 <- CreateSeuratObject(counts = GSM6321080.count, project = "GSM6321080")

AgingMouse=merge(x=GSM6321077,y = c(GSM6321078,GSM6321079,GSM6321080), project = "AgingMouse_PFC")
AgingMouse[["percent.mt"]] <- PercentageFeatureSet(AgingMouse, pattern = "^mt-")
summary(AgingMouse$percent.mt)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.00000 0.00000 0.00000 0.02051 0.01099 8.71022
summary(AgingMouse$nFeature_RNA)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#394    2000    3146    3507    5010   13789
pdf("AgingMouseQC.pdf",width=30)
VlnPlot(AgingMouse, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
#remove the cells with AgingMousetremely high nFeature_RNA
AgingMouse=subset(AgingMouse, subset = nFeature_RNA > 200  & nFeature_RNA < 10000)

AgingMouse <- NormalizeData(AgingMouse, normalization.method = "LogNormalize", scale.factor = 10000)
AgingMouse <- FindVariableFeatures(AgingMouse, selection.method = "vst", nfeatures = 2000)
AgingMouse <- ScaleData(AgingMouse)
AgingMouse <- RunPCA(AgingMouse, features = VariableFeatures(object = AgingMouse))
AgingMouse<- RunUMAP(AgingMouse, reduction = "pca", dims = 1:50)
AgingMouse<- RunTSNE(AgingMouse, reduction = "pca", dims = 1:50)
AgingMouse<- FindNeighbors(AgingMouse, reduction = "pca", dims = 1:50)
AgingMouse<- FindClusters(AgingMouse, resolution = 0.5)

pdf("AgingMouseLabel.pdf",width=9)
DimPlot(AgingMouse,reduction ="umap",label=TRUE)&NoLegend()&theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

pdf("AgingMouseSample.pdf",width=10)
DimPlot(AgingMouse, reduction ="umap",group.by="orig.ident")&theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.5, linetype="solid"))
dev.off()

pdf("NeuronMarker.pdf",width=5,height=7)
FeaturePlot(AgingMouse, features=c("Camk2a","Syt1","Slc17a7","Nrgn","Gad1","Gad2"),raster=TRUE,ncol=2)&NoLegend()&theme(plot.title=element_text(size=10,face ="italic"),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("AgingMouseCellTypeMarker.pdf",width=10,height=13)
FeaturePlot(AgingMouse,reduction ="umap",raster=TRUE,features=c("Camk2a","Syt1","Slc17a7","Nrgn","Gad1","Gad2","Aqp4","Slc1a2","Mbp","Mobp","Plp1","Csf1r","Ctss","Vcan","Flt1","Vtn","Cldn5","F13a1","Cd3e"))&NoLegend()&theme(plot.title=element_text(size=10,face ="italic"),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()


new.cluster.ids <- c("Oli","Ex","Oli","Ex","Ast","Ex","Mic","Ex","Ex","In","Opc","In","Ex","In","In","Ex","Ex","In","In","Ex","In","In","Ex","Mac","Ex","Vascular","Vascular","Ast","In","Vascular","Ex","Ast")
names(new.cluster.ids) <- levels(AgingMouse)
AgingMouse <- RenameIdents(AgingMouse, new.cluster.ids)
AgingMouse$cellType=Idents(AgingMouse)

colorPalette=colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(AgingMouse$cellType))) 
pdf("AgingMouseCellType_UMAP.pdf",width=4,height=3.5)
DimPlot(AgingMouse, reduction="umap",group.by="cellType",label=T,label.size = 6,raster=TRUE,cols =colorPalette)&NoLegend()&NoAxes()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.0, linetype="solid"))
dev.off()
pdf("AgingMouseCellType_TSNE.pdf",width=4,height=3.5)
DimPlot(AgingMouse, reduction="tsne",group.by="cellType",label=T,label.size = 6,raster=TRUE,cols =colorPalette)&NoLegend()&NoAxes()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.0, linetype="solid"))
dev.off()
pdf("AgingMouseCellType_TSNE_nolabel.pdf",width=4,height=3.5)
DimPlot(AgingMouse, reduction="tsne",group.by="cellType",label=F,raster=TRUE,cols =colorPalette)&NoLegend()&NoAxes()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.0, linetype="solid"))
dev.off()

saveRDS(AgingMouse,file="/data2/deng/Aging/Ex/Mouse_Allen_GSE207848/Mouse_Allen_GSE207848.rds")
#32285 features across 48052 samples within 1 assay

setwd("/Projects/deng/Aging/Ex/Mouse_Allen_GSE207848")
AgingMouse=readRDS("/data2/deng/Aging/Ex/Mouse_Allen_GSE207848/Mouse_Allen_GSE207848.rds")

Ex=subset(AgingMouse,seurat_clusters%in%c(1,3,5,7,8,12,15,16,19,22,24,30))
pdf("Ex_AgingMouse.pdf",width=9)
DimPlot(Ex,reduction ="umap",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.5, linetype="solid"))
dev.off()
Ex <- NormalizeData(Ex, normalization.method = "LogNormalize", scale.factor = 10000)
Ex <- FindVariableFeatures(Ex, selection.method = "vst", nfeatures = 2000)
Ex <- ScaleData(Ex)
Ex <- RunPCA(Ex, features = VariableFeatures(object = Ex))
Ex<- RunTSNE(Ex, reduction = "pca", dims = 1:30)
Ex<- RunUMAP(Ex, reduction = "pca", dims = 1:30)
Ex<- FindNeighbors(Ex, reduction = "pca", dims = 1:30)
Ex<- FindClusters(Ex, resolution = 0.3) 



pdf("Ex_cyclingMarker.pdf",width=5,height=5)
FeaturePlot(Ex,reduction = "tsne",raster=TRUE,features=c("Nrgn","Hsp90aa1","Cdkn1a","Cdkn2a"),ncol=2)&NoLegend()&theme(plot.title=element_text(size=10,face ="italic"),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("Ex_cyclingMarkerVln.pdf",width=7,height=10)
VlnPlot(Ex,features=c("Nrgn","Hsp90aa1","Cdkn1a","Cdkn2a"),ncol=1)&NoLegend()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

colorPalette=colorRampPalette(brewer.pal(8, "Set3"))(length(unique(Ex$seurat_clusters))) 
colorPalette[3]="Firebrick3"
pdf("Ex_Label_TSNE.pdf",width=4,height=3.5)
DimPlot(Ex,label=T,reduction="tsne",label.size = 6,raster=TRUE,cols =colorPalette)&NoLegend()&NoAxes()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.0, linetype="solid"))
dev.off()
pdf("Ex_Label_TSNE_NoLabel.pdf",width=4,height=3.5)
DimPlot(Ex,reduction="tsne",raster=TRUE,cols =colorPalette)&NoLegend()&NoAxes()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.0, linetype="solid"))
dev.off()
pdf("Ex_Label_UMAP.pdf",width=4,height=3.5)
DimPlot(Ex,label=T,reduction="umap",label.size = 6,raster=TRUE,cols =colorPalette)&NoLegend()&NoAxes()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.0, linetype="solid"))
dev.off()
pdf("Ex_nFeatureRNA.pdf",width=4.5,height=3.5)
FeaturePlot(Ex,raster=TRUE,,reduction="tsne",features="nFeature_RNA")&NoAxes()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.0, linetype="solid"))
dev.off()
pdf("Ex_nFeatureCount.pdf",width=4.5,height=3.5)
FeaturePlot(Ex,raster=TRUE,,reduction="tsne",features="nCount_RNA")&NoAxes()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.0, linetype="solid"))
dev.off()


layerMarker=c("GLRA3","TLE1","TLE3","MDGA1","CARTPT","CUX1","CUX2","KITLG","RASGRF2","PVRL3","PRSS12","LUX2","DTX4","GPR6","UNC5D","MEF2C","RORB","SATB2","DKK3","LDB2","SYT9","OPN3","NEFH","CNTN6","FOXO1","LIX1","S100A10","CRIM1","KCNK2","SULF2","PCP4","HTR2C","FEZF2","TOX","ETV1","RPRM","RXFP1","FOXP2","CRYM","OTX1","SOX5","SEMA3E","NR4A3","LXN","PPP1R1B","SYT6","OPRK1","NR4A2","SYNPR","TLE","NTNG2","ADRA2A")
library(homologene)
layerMarkerForMouse=human2mouse(layerMarker)

Ex=subset(Ex,seurat_clusters%in%c(0:15)) #removed C16 and C17 with small cells
ExExpr <- AverageExpression(Ex)[["RNA"]]
hclust=hclust(as.dist(1-cor(ExExpr)),method="ward.D2")
pdf("clusterRelationship.pdf",width=10,height=5)
plot(hclust,hang=-1)
dev.off()

order=hclust$labels[hclust$order]
t=DotPlot(Ex,features=layerMarkerForMouse$mouseGene)
t=ggplot(t$data, aes(factor(features.plot,levels=unique(features.plot)),factor(id,levels=order),size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_color_gradient2(low="blue",mid="white",high = "red")+
labs(color="Expression",size="Percent",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),face="italic",angle=90,vjust=1,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("layerMarker.pdf",width=12,height=5)
print(t)
dev.off()


pdf("LayerMarkerFeture.pdf",width=8,height=8)
FeaturePlot(Ex,raster=TRUE,,reduction="umap",features=c("Glra3","Cux2","Rorb","Pcp4","Crim1","Rxfp1","Ppp1r1b","Ntng2"))&NoLegend()&theme(plot.title=element_text(size=10,face ="italic"),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()


Ex[["Layer"]]=NA
Ex[["Layer"]][Ex$seurat_clusters%in%c(1,4,6,9),]="L2/3"
Ex[["Layer"]][Ex$seurat_clusters%in%c(3,5,7,14,15),]="L4/5"
Ex[["Layer"]][Ex$seurat_clusters%in%c(0,8,10),]="L4/6"
Ex[["Layer"]][Ex$seurat_clusters%in%c(11,12),]="L5/6"
Ex[["Layer"]][Ex$seurat_clusters%in%c(13),]="L6"
Ex[["Layer"]][Ex$seurat_clusters%in%c(2),]="C2"
Ex$Layer=factor(Ex$Layer,levels=c("C2","L2/3","L4/5","L4/6","L5/6","L6"))

LayerPalette=c("Firebrick3","#DBE419","#74CF56","#21D66A","#238B8E","#2E718E")
pdf("ExAll.Layer.tsne.pdf",width=4,height=3.5)
DimPlot(Ex, raster=TRUE,reduction = "tsne",group.by="Layer",label=T,label.size = 6,cols =LayerPalette)&NoLegend()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
dev.off()
pdf("ExAll.Layer.umap.pdf",width=4,height=3.5)
DimPlot(Ex, raster=TRUE,reduction = "umap",group.by="Layer",label=T,label.size = 6,cols =LayerPalette)&NoLegend()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
dev.off()



saveRDS(Ex,file="/data2/deng/Aging/Ex/Mouse_Allen_GSE207848/Mouse_Allen_GSE207848_Ex.rds")


setwd("/Projects/deng/Aging/Ex/Mouse_Allen_GSE207848")
Ex=readRDS("/data2/deng/Aging/Ex/Mouse_Allen_GSE207848/Mouse_Allen_GSE207848_Ex.rds")
#32285 features across 18623 samples within 1 assay

colorPalette=colorRampPalette(brewer.pal(8, "Set3"))(length(unique(Ex$seurat_clusters))) 
colorPalette[3]="Firebrick3"

CCGene=read.table("/Projects/deng/Aging/Ex/cellCycleGene/CCGeneCombine.txt",header=T,sep="\t")
MouseID=human2mouse(CCGene$NAME)
CCGeneCombine=merge(CCGene,MouseID,by.x="NAME",by.y="humanGene")
genes=intersect(rownames(Ex),CCGeneCombine$mouseGene)
length(genes) #347
CCGeneCombine=CCGeneCombine[CCGeneCombine$mouseGene%in%genes,]
CCGenelist=split(CCGeneCombine$mouseGene,CCGeneCombine$PHASE)
Idents(Ex)=Ex$seurat_clusters
Ex <- AddModuleScore(
  object = Ex,
  features = CCGenelist,
  name = names(CCGenelist),
  ctrl = 100,
)
pdf("Cells_Phase_PvalueDis/Ex_cellCycleScore.pdf",width=6)
VlnPlot(Ex,c("G1.S1","S.phase5","G22","G2.M3","M.G14"),ncol=1,pt.size=0,cols=colorPalette)&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()



cluster=c(0:15)
CycleScore=Ex@meta.data[,c("seurat_clusters","G1.S1","S.phase5","G22","G2.M3","M.G14")]
PhaseList=c("G1.S1","S.phase5","G22","G2.M3","M.G14")

for(p in PhaseList){
 PvalueBetCluster=matrix(data = NA, nrow = length(cluster), ncol = length(cluster), dimnames = list(cluster, cluster))
 for(i in 1:length(cluster)){
  for( j in 1:length(cluster)){
    C1_score=CycleScore[CycleScore$seurat_clusters==cluster[i],p]
    C2_score=CycleScore[CycleScore$seurat_clusters==cluster[j],p]
    #delta=mean(C1_score)-mean(C2_score) #the row represent the current cluster
    delta=quantile(C1_score,0.9)[[1]]-quantile(C2_score,0.9)[[1]]
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

 pdf(paste0("Cells_Phase_PvalueDis/Ex_",p,"_PvalueDis.pdf",sep=""),height=1.5,width=6)
 print(t)
 dev.off()

}




table(Ex$seurat_clusters)
Marker=FindMarkers(Ex,ident.1=2,test.use ="MAST")
Marker=Marker[Marker$p_val<0.01,]
Marker$pattern=ifelse(Marker$avg_log2FC>0,"UpInLS_Mouse","DnInLS_Mouse")
table(Marker$pattern)
#DnInLS_Mouse UpInLS_Mouse
#534         1076
write.table(Marker,file="C2Marker.txt",quote=F,sep="\t")
pdf("C2Marker.pdf",width=10,height=13)
FeaturePlot(Ex,reduction ="tsne",raster=TRUE,features=c(rownames(Marker)[1:20]))&NoLegend()&theme(plot.title=element_text(size=10,face ="italic"),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()


library(homologene) #1.4.68.19.3.27
C2MarkerFromMouse=read.table("/Projects/deng/Aging/Ex/Mouse_Allen_GSE207848/C2Marker.txt",header=T)
C2MarkerFromMouse$Symbol=rownames(C2MarkerFromMouse)
LSMarkerFromHuman=read.table("/data2/deng/Aging/Ex/AllEx/C5Feature/C5MarkerFC1_MAST.txt",header=T)
LSMarkerFromHuman$Symbol=rownames(LSMarkerFromHuman)
#767

#https://www.genenames.org/tools/hcop/
OrthologyGene=read.csv("/Projects/deng/Public/HGNC/human_mouse_hcop_fifteen_column.txt",header=T,sep="\t")
MouseGeneFromHumanLS=OrthologyGene[OrthologyGene$human_symbol%in%LSMarkerFromHuman$Symbol,c("human_symbol","mouse_symbol")]
MouseGeneFromHumanLS=LSMarkerFromHumanOrthology[!duplicated(LSMarkerFromHumanOrthology$human_symbol),]
MouseGeneFromHumanLS=MouseGeneFromHumanLS[!duplicated(MouseGeneFromHumanLS$mouse_symbol),] #Slc25a4 and Uqcrh
colnames(MouseGeneFromHumanLS)=c("humanGene","mouseGene")
dim(MouseGeneFromHumanLS)
#717
setdiff(LSMarkerFromHuman$Symbol,MouseGeneFromHumanLS$humanGene)

MouseGeneFromHumanLS.df=merge(MouseGeneFromHumanLS,LSMarkerFromHuman,by.x="humanGene",by.y="Symbol")
MouseGeneFromHumanLS.df$Symbol=MouseGeneFromHumanLS.df$mouseGene

MouseGeneFromHumanLS.df$pattern=ifelse(MouseGeneFromHumanLS.df$avg_log2FC>0,"UpInLS_Human","DnInLS_Human")

MouseGeneFromHumanLS.df=MouseGeneFromHumanLS.df[,c("p_val","avg_log2FC","pct.1","pct.2","p_val_adj","pattern","Symbol")]
rownames(MouseGeneFromHumanLS.df)=MouseGeneFromHumanLS.df$Symbol

MarkerGeneCombine=rbind(C2MarkerFromMouse,MouseGeneFromHumanLS.df)
#MarkerGeneCombine=MarkerGeneCombine[abs(MarkerGeneCombine$avg_log2FC)>0.5,]
table(MarkerGeneCombine$pattern)
#DnInLS_Human DnInLS_Mouse UpInLS_Human UpInLS_Mouse
#         345          534          372         1076
MarkerGeneCombine.list=split(MarkerGeneCombine$Symbol,MarkerGeneCombine$pattern)
library(UpSetR)
pdf("OverlapMarker.pdf",width=6,height=4)
upset(fromList(MarkerGeneCombine.list),)
dev.off()

SharedGene=c(
            intersect(MarkerGeneCombine.list$UpInLS_Human,MarkerGeneCombine.list$UpInLS_Mouse),
            intersect(MarkerGeneCombine.list$DnInLS_Human,MarkerGeneCombine.list$DnInLS_Mouse)
            )

ShareGeneInHuman=MouseGeneFromHumanLS.df[SharedGene,]
colnames(ShareGeneInHuman)=paste0(colnames(ShareGeneInHuman),"_hsa")
ShareGeneInMouse=C2MarkerFromMouse[SharedGene,]
all(rownames(ShareGeneInHuman)==rownames(ShareGeneInMouse))
colnames(ShareGeneInMouse)=paste0(colnames(ShareGeneInMouse),"_mmu")

ShareGeneInfor=data.frame(cbind(ShareGeneInHuman,ShareGeneInMouse),check.names=F)
ShareGeneInfor$FCTotal=abs(ShareGeneInfor[,2]+ShareGeneInfor[,9])
ShareGeneInfor=ShareGeneInfor[order(ShareGeneInfor$FCTotal,decreasing=T),]

ShareGeneInfor$mouseGene=rownames(ShareGeneInfor)
ShareGeneInfor.gene=merge(MouseGeneFromHumanLS,ShareGeneInfor,by.x="mouseGene",check.names=F)

ShareGeneInfor.gene=ShareGeneInfor.gene[order(ShareGeneInfor.gene$FCTotal,decreasing=T),]
write.table(ShareGeneInfor.gene,file="SharedMarkerGeneInfo.txt",sep="\t",quote=F,row.names=F)


pdf("SharedMarkerForLSInMouse.pdf",width=10,height=13)
FeaturePlot(Ex,reduction ="tsne",raster=TRUE,features=c(ShareGeneInfor.gene$mouseGene[1:20]))&NoLegend()&theme(plot.title=element_text(size=10,face ="italic"),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()


ShareGeneInfor.gene=read.table("SharedMarkerGeneInfo.txt",header=T,sep="\t")

ExAll.integrated=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
ExAll.integrated=subset(ExAll.integrated,seurat_clusters%in%c(0:15))
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

DefaultAssay(ExAll.integrated)="RNA"
pdf("SharedMarkerForLSInHuman.pdf",width=10,height=13)
FeaturePlot(ExAll.integrated,reduction ="tsne",raster=TRUE,features=c(ShareGeneInfor.gene$humanGene[1:20]))&NoLegend()&theme(plot.title=element_text(size=10,face ="italic"),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

ShareGeneInfor.gene.show=ShareGeneInfor.gene[c(1:100),]
ShareGeneInfor.gene.show=ShareGeneInfor.gene.show[order(ShareGeneInfor.gene.show$pattern_hsa,ShareGeneInfor.gene.show$FCTotal,decreasing=TRUE),]

geneOrder=ShareGeneInfor.gene.show$humanGene
t=DotPlot(ExAll.integrated,features=ShareGeneInfor.gene.show$humanGene,group.by="Layer")
result=t$data
g=ggplot(result, aes(factor(features.plot,levels=geneOrder), id, size= pct.exp,color=avg.exp.scaled)) +
geom_point(alpha = 0.9)+
scale_colour_viridis_c()+
labs(color="avg.exp.scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+
theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(face="italic",color="black",angle=90,vjust=0.5,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("SharedMarkerForLSInHuman_dotplot.pdf",width=15,height=4)
print(g)
dev.off()

geneOrder=ShareGeneInfor.gene.show$mouseGene
t=DotPlot(Ex,features=ShareGeneInfor.gene.show$mouseGene,group.by="Layer")
result=t$data
g=ggplot(result, aes(factor(features.plot,levels=geneOrder), id, size= pct.exp,color=avg.exp.scaled)) +
geom_point(alpha = 0.9)+
scale_colour_viridis_c()+
labs(color="avg.exp.scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+
theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(face="italic",color="black",angle=90,vjust=0.5,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("SharedMarkerForLSInMouse_dotplot.pdf",width=15,height=4)
print(g)
dev.off()



pdf("NeuronMarkerByHuman.pdf",width=20,height=3)
FeaturePlot(ExAll.integrated, raster=TRUE,reduction = "tsne",features=c("NRGN","RBFOX3","NEFM","NEFH","SYP","MAP2","DLG4"),ncol=7)&NoLegend()&theme(plot.title=element_text(face="italic"),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.0, linetype="solid"))
dev.off()
#NeuN: RBFOX3
#Recombinant Anti-160 kD Neurofilament Medium antibody: NEFM
#200kDa neurofilament heavy: NEFH
#synaptophysin: SYP
#PSD95: DLG4
pdf("NeuronMarkerByAgingHuman.pdf",width=20,height=3)
FeaturePlot(subset(ExAll.integrated,Age>60), raster=TRUE,reduction = "tsne",features=c("NRGN","RBFOX3","NEFM","NEFH","SYP","MAP2","DLG4"),ncol=7)&NoLegend()&theme(plot.title=element_text(face="italic"),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.0, linetype="solid"))
dev.off()

pdf("NeuronMarkerByMouse.pdf",width=20,height=3)
FeaturePlot(Ex,reduction = "tsne",raster=TRUE,features=c("Nrgn","Rbfox3","Nefm","Nefh","Syp","Map2","Dlg4"),ncol=7)&NoLegend()&theme(plot.title=element_text(face="italic"),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.0, linetype="solid"))
dev.off()




pdf("LDHAMarkerByMouse.pdf",width=6,height=3)
FeaturePlot(Ex,reduction = "tsne",raster=TRUE,features=c("Ldha","Ldhb"))&NoAxes()&NoLegend()&theme(plot.title=element_text(face="italic"),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.0, linetype="solid"))
dev.off()

pdf("LDHAMarkerByByHuman.pdf",width=6,height=3)
FeaturePlot(ExAll.integrated, raster=TRUE,reduction = "tsne",features=c("LDHA","LDHB"))&NoAxes()&NoLegend()&theme(plot.title=element_text(face="italic"),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.0, linetype="solid"))
dev.off()


ExSmall <- Ex[, sample(colnames(Ex), size = 10000, replace=F)]
ExSmall <- ExSmall[!grepl("^mt-", rownames(ExSmall)), ] #remove mitocohndrial genes
table(ExSmall$seurat_clusters)
library("infercnv")
counts_matrix = as.matrix(ExSmall@assays$RNA@counts)
pData=ExSmall@meta.data$seurat_clusters
names(pData)=rownames(ExSmall@meta.data)
pData=data.frame(pData)
annotation=read.table("/Projects/deng/refGene/gencode.vM27.FullName.removeDupliSymbol.geneType.txt",header=TRUE,row.names=1)
anno=data.frame("Chr"=annotation$Chr,"Start"=annotation$Start,"End"=annotation$End)
rownames(anno)=annotation$Symbol
genes=intersect(rownames(counts_matrix),rownames(anno))
length(genes)#31624
anno=anno[genes,]
anno$Chr=factor(anno$Chr,levels=paste0("chr",c(1:22,"X","Y"),sep=""))
anno=anno[order(anno$Chr,anno$Star,anno$End),]

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix[genes,],
                                      annotations_file=pData,
                                      delim="\t",
                                      gene_order_file=anno[genes,],
                                      ref_group_names=NULL)

infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics   
                               denoise=TRUE,
                               cluster_by_groups = TRUE,
                               HMM=TRUE,
                               out_dir="/data2/deng/Aging/Ex/Mouse_Allen_GSE207848/InferCNV"
                               )
