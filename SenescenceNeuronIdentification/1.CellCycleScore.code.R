library(Seurat)
library(RColorBrewer)
setwd("/boot3/bixm/08_HNLongevity/pipline")
#using your own pipine to generated this neuron dataset
ExIntegrated=readRDS("Data/ExIntegrated.Downsampling.rds")
#set the default assay as RNA if you using the integrated method
DefaultAssay(ExIntegrated)="RNA"
#obtain the cell cycle assocaited genes (you can download it from sources/CCGeneCombine.txt)
CCGene=read.table("sources/CCGeneCombine.txt",header=T,sep="\t") #500
genes=intersect(CCGene$NAME,rownames(ExIntegrated))
length(genes) #364

#set the clusteras the defult category
Idents(ExIntegrated)=ExIntegrated$seurat_clusters

CCGene=CCGene[CCGene$NAME%in%genes,]
CCGenelist=split(CCGene$NAME,CCGene$PHASE)
ExAll.integrated <- AddModuleScore(
  object = ExIntegrated,
  features = CCGenelist,
  name = names(CCGenelist),
  ctrl = 100,
)

#set your own color palette
colorPalette=colorRampPalette(brewer.pal(8, "Set3"))(length(unique(ExAll.integrated$seurat_clusters))) 
colorPalette[4]="Orange"
colorPalette[6]="Firebrick3"
colorPalette[11]="LightCyan"

pdf("ResultGraph/cellCycleScore.pdf",width=6)
VlnPlot(ExAll.integrated,c("G1.S1","S.phase5","G22","G2.M3","M.G14"),group.by="seurat_clusters",ncol=1,pt.size=0,cols =colorPalette)&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", linewidth=1.5, linetype="solid"))
dev.off()



