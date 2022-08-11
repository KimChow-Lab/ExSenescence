library(monocle) #local monocle2
packageVersion("monocle")
library(Seurat)

####################################################################################
############                 monocle analysis                ##############
####################################################################################

ExCtrl=readRDS("/Projects/deng/Aging/Ex/AllControlEx/ExCtrl.integrated.rds")
MathyCtrl=subset(ExCtrl,Source %in% "Mathy") # Replaced by Nagy, Yang, and Lau when needed

DefaultAssay(MathyCtrl)="RNA"
data  <-  as.matrix(MathyCtrl@assays$RNA@counts)
celldata <- as.data.frame(MathyCtrl@meta.data)
genedata <- as.data.frame(x = row.names(MathyCtrl), row.names = row.names(MathyCtrl))
colnames(genedata) <- "gene_short_name"
pd <- new("AnnotatedDataFrame", data = celldata)
fd <- new("AnnotatedDataFrame", data = genedata)
MathyCtrl.cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
MathyCtrl.cds <- estimateSizeFactors(MathyCtrl.cds)
MathyCtrl.cds <- estimateDispersions(MathyCtrl.cds)
MathyCtrl.cds <- detectGenes(MathyCtrl.cds, min_expr = 0.1)
head(fData(MathyCtrl.cds))
expressed_genes <- row.names(subset(fData(MathyCtrl.cds),num_cells_expressed >= 10))
#diff_test_res <- differentialGeneTest(MathyCtrl.cds[expressed_genes,],fullModelFormulaStr = "~seurat_clusters")
#ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
MathyCtrl.cds <- reduceDimension(MathyCtrl.cds, max_components = 2,method = 'DDRTree')
MathyCtrl.cds <- orderCells(MathyCtrl.cds)
saveRDS(MathyCtrl.cds,file="MathyCtrl.cds.rds")

#----Figure 3B, S4----
pData(MathyCtrl.cds)$cellType=ifelse(pData(MathyCtrl.cds)$seurat_clusters %in% c(3,4),paste0("C",pData(MathyCtrl.cds)$seurat_clusters,sep=""),"Others")
pdf("MathyCtrlTrajectoryCellType.pdf",height=6)
plot_cell_trajectory(MathyCtrl.cds, color_by = "cellType",cell_size = 1.5)+
scale_color_manual(values=c("#A3A500","#6BB100","Gainsboro"))+
theme(plot.title=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()


tiff("MathyCtrlTrajectory.tiff",width=350,height=250)
plot_cell_trajectory(MathyCtrl.cds, color_by = "State",cell_size = 1.5)&NoLegend()&theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()


pdf("MathyCtrlPseudotime.pdf",height=6)
plot_cell_trajectory(MathyCtrl.cds, color_by = "Pseudotime",cell_size = 1)+scale_color_viridis_c()
dev.off()




changed_genes <- row.names(subset(fData(MathyCtrl.cds),gene_short_name %in% c("HSP90AA1","ALDOA","YWHAG","YWHAH","YWHAB")))
valid_cells <- row.names(subset(pData(MathyCtrl.cds),seurat_clusters %in% c(0,3,4)))
MathyCtrlSubset.cds <- MathyCtrl.cds[,valid_cells]
pdf("plot_genes_branched_pseudotime4DEG.pdf",height=9)
plot_genes_branched_pseudotime(MathyCtrlSubset.cds[changed_genes,],
                       color_by = "seurat_clusters",
                       ncol = 1)+scale_color_manual(values=c("#F8766D","#6BB100","#A3A500"))
dev.off()
