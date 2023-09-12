#######################################################################################
##############         Preprocessing snRNA from Mathy et al.    #######################
#######################################################################################
library(CellChat)
library(ggplot2)
setwd("/Projects/deng/Aging/Ex/AllEx/cellChat")

MathysAll=readRDS("/Projects/deng/Alzheimer/syn18485175/AD.rds")
Glia=subset(MathysAll,OldCellType %in% c("Mic","Ast","Oli"))
Glia=subset(Glia,seurat_clusters%in%c(0,2,6,13))
Glia$Type=Glia$OldCellType

ExAll.integrated=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
table(ExAll.integrated$Source)
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

Mathys=subset(ExAll.integrated,Source=="Mathy")
Mathys=subset(Mathys,idents=c(0:15))

Mathys$Type=Mathys$Layer
Combined <- merge(Glia, y = Mathys, add.cell.ids = c("Glia", "Mathys"), project = "Combined")
Combined <- NormalizeData(Combined, normalization.method = "LogNormalize", scale.factor = 10000)
Combined <- FindVariableFeatures(Combined, selection.method = "vst", nfeatures = 2000)
Combined <- ScaleData(Combined)
Combined <- RunPCA(Combined, features = VariableFeatures(object = Combined),npcs = 100)
Combined<- RunTSNE(Combined, reduction = "pca", dims = 1:20)
Combined<- RunUMAP(Combined, reduction = "pca", dims = 1:20)
Combined<- FindNeighbors(Combined, reduction = "pca", dims = 1:20)

table(Combined$Type)
Idents(Combined)=factor(Combined$Type,levels=c("ES","LS","L2/3","L4/5","L4/6","L5","L5/6","L6","Oli","Ast","Mic"))
LayerPalette=c("Orange","Firebrick3","#DBE419","#74CF56","#21D66A","#3DBC75","#238B8E","#2E718E")
colorPalette=c(LayerPalette,"#B385FF","#EF67EB","#FF63B6")
pdf("MathyData.Type.pdf",width=6.5,height=5)
DimPlot(Combined, reduction="umap",label=T,label.size = 6,cols =colorPalette,raster=TRUE)&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", linewidth=1.5, linetype="solid"))
dev.off()
pdf("MathyData_tsne.Type.pdf",width=6.5,height=5)
DimPlot(Combined, reduction="tsne",label=T,label.size = 6,cols =colorPalette,raster=TRUE)&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", linewidth=1.5, linetype="solid"))
dev.off()


selected_f_10xGenomic <- rownames(Combined)[Matrix::rowSums(Combined) > 3]
length(selected_f_10xGenomic)
Combined=subset(Combined,features=selected_f_10xGenomic)
saveRDS(Combined,file="Combined.Mathys.seurat.rds")

exprMat  <-  as.matrix(Combined@assays$RNA@data)
meta = Combined@meta.data
table(meta$Type)

cellchat <- createCellChat(object = Combined, group.by = "Type", assay = "RNA")
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
options(stringsAsFactors = FALSE)
#manually download the interaction from cellchatDB
interaction_input <- read.table(file = '/Projects/deng/Alzheimer/syn18485175/cellCycle/cellchatDB/Interaction.txt', row.names = 1,header=T,sep="\t") #add the celltalkDB into database
MyCellChatDB <- list()
MyCellChatDB$interaction <- interaction_input
MyCellChatDB$complex <- CellChatDB$complex
MyCellChatDB$cofactor <- CellChatDB$cofactor
MyCellChatDB$geneInfo <- CellChatDB$geneInfo
table(MyCellChatDB$interaction$annotation)
CellChatDB.use <- MyCellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#cellchat <- projectData(cellchat, PPI.human) #(optional)
cellchat <- computeCommunProb(cellchat) #very slow, waiting.....
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = "cellchatExCAndGlia_Layer.rds")



cellchat=readRDS("cellchatExCAndGlia_Layer.rds") 

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

LayerPalette=c("Orange","#DBE419","#74CF56","#21D66A","#3DBC75","#238B8E","#2E718E","Firebrick3")
colorPalette=c("#EF67EB",LayerPalette,"#FF63B6","#B385FF")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",width = 6,cluster.cols=TRUE,color.use=colorPalette)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",width = 6,cluster.cols=TRUE,color.use=colorPalette)
pdf("netAnalysis_signalingRole_heatmap_Layers.pdf")
print(ht1+ht2)
dev.off()

levels(cellchat@idents)
cellchat@idents <- factor(cellchat@idents,levels=c("ES","LS","L2/3","L4/5","L4/6","L5","L5/6","L6","Oli","Ast","Mic"))
#Figure S10
pdf("GliaAsSourcesBubble_Layers.pdf",height=15,width=6)
netVisual_bubble(cellchat, targets.use = c("ES","LS","L2/3","L4/5","L4/6","L5","L5/6","L6"),  sources.use = c("Oli","Ast","Mic"), remove.isolate = FALSE)
dev.off()
pdf("GliaAsTargetsBubble_Layers.pdf",height=15,width=6)
netVisual_bubble(cellchat, sources.use = c("ES","LS","L2/3","L4/5","L4/6","L5","L5/6","L6"),  targets.use = c("Oli","Ast","Mic"), remove.isolate = FALSE)
dev.off()



############Expression for specific gene list###################

setwd("/Projects/deng/Aging/Ex/AllEx/cellChat")
DefaultAssay(ExAll.integrated)="RNA"
pdf("ReceptorLost.pdf",width=10,height=10)
FeaturePlot(ExAll.integrated,features=c("LRP1","LRP4","LRP6","LRP8","CD46","CD81","GRM7"),raster=TRUE)&NoLegend()&NoAxes()&theme(plot.title=element_text(face="italic"),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

pdf("LigandForMic.pdf",width=5,height=5)
FeaturePlot(ExAll.integrated,features=c("CD47","CD200","CX3CL1"),raster=TRUE)&NoLegend()&NoAxes()&theme(plot.title=element_text(face="italic"),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("LigandForMic_DotPlot.pdf",width=5,height=5)
DotPlot(ExAll.integrated,features=c("CD47","CD200","CX3CL1"))
dev.off()

pdf("ReceptorTest.pdf",width=10,height=10)
FeaturePlot(ExAll.integrated,features=c("CXCL14","RELN","CNR1","CHRNA7","INPP4B"),raster=TRUE)&NoLegend()&NoAxes()&theme(plot.title=element_text(face="italic"),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()


Combined=readRDS("Combined.Mathys.seurat.rds")
Combined$Layer=factor(Combined$Type,levels=c("ES","LS","L2/3","L4/5","L4/6","L5","L5/6","L6","Oli","Ast","Mic"))


t=DotPlot(ExAll.integrated,features=c("CD47","CD200","CX3CL1"),group.by="Layer")
result=t$data
result$id=factor(result$id,levels=rev(c("ES","LS","L2/3","L4/5","L4/6","L5","L5/6","L6","Oli","Ast","Mic")))
g=ggplot(result, aes(features.plot, id, size= pct.exp,color=avg.exp.scaled)) +
geom_point(alpha = 0.9)+
scale_color_gradient(low = "Blue",high = "Orange")+
labs(color="avg.exp.scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+
theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(face="italic",color="black",angle=90,vjust=0.5,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("LigandForMic_DotPlot.pdf",width=3,height=3)
print(g)
dev.off()
