#######################################################################################
##############         Preprocessing snRNA from Mathy et al.    #######################
#######################################################################################

MathysAll=readRDS("/Projects/deng/Alzheimer/syn18485175/MathyAll.rds")
Glia=subset(AD,OldCellType %in% c("Mic","Ast","Oli"))
Glia=subset(Glia,seurat_clusters%in%c(0,2,6,13))
Glia$Type=Glia$OldCellType

ExAll.integrated=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
table(ExAll.integrated$Source)
Mathys=subset(ExAll.integrated,Source=="Mathy")
Mathys=subset(Mathys,idents=c(0:15))
Mathys$Type=paste0("C",Mathys$seurat_clusters)
Combined <- merge(Glia, y = Mathys, add.cell.ids = c("Glia", "Mathys"), project = "Combined")
Combined <- NormalizeData(Combined, normalization.method = "LogNormalize", scale.factor = 10000)
Combined <- FindVariableFeatures(Combined, selection.method = "vst", nfeatures = 2000)
Combined <- ScaleData(Combined)
Combined <- RunPCA(Combined, features = VariableFeatures(object = Combined),npcs = 100)
pdf("CombinedElbowPlot.pdf")
ElbowPlot(Combined,ndims = 100)
dev.off()
Combined<- RunTSNE(Combined, reduction = "pca", dims = 1:20)
Combined<- FindNeighbors(Combined, reduction = "pca", dims = 1:20)

Idents(Combined)=factor(Combined$Type,levels=c(paste("C",c(0:15),sep=""),"Oli","Ast","Mic"))
pdf("MathyData.Type.pdf",width=6,height=5)
DimPlot(Combined, reduction="tsne",label=T,label.size = 6)&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

library(CellChat)
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
saveRDS(cellchat, file = "cellchatExCAndGlia.rds")

library(CellChat)
cellchat=readRDS("cellchatExCAndGlia.rds") 
levels(cellchat@idents)
cellchat@idents <- factor(cellchat@idents,levels=c(paste0("C",c(3,5,0:2,4,6:15)),"Oli","Ast","Mic"))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

#Figure 4C
pdf("GliaAsSourcesBubble.pdf",height=10,width=6)
netVisual_bubble(cellchat, targets.use = c("C0","C1","C2","C3","C4","C5","C6"),  sources.use = c("Oli","Ast","Mic"), remove.isolate = FALSE)
dev.off()
pdf("GliaAsTargetsBubble.pdf",height=12,width=6)
netVisual_bubble(cellchat, sources.use = c("C0","C1","C2","C3","C4","C5","C6"),  targets.use = c("Oli","Ast","Mic"), remove.isolate = FALSE)
dev.off()

#Figure S10
pdf("GliaAsSourcesBubble.pdf",height=10,width=6)
netVisual_bubble(cellchat, targets.use = c("C0","C1","C2","C3","C4","C5","C6"),  sources.use = c("Oli","Ast","Mic"), remove.isolate = FALSE)
dev.off()
pdf("GliaAsTargetsBubble.pdf",height=12,width=6)
netVisual_bubble(cellchat, sources.use = c("C0","C1","C2","C3","C4","C5","C6"),  targets.use = c("Oli","Ast","Mic"), remove.isolate = FALSE)
dev.off()
