library(SCENIC) 
library(Seurat)
library(pheatmap)
library(viridis)

####################################################################################
############         activated transcription factors analysis         ##############
####################################################################################

ExAll=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
ExAllTmp=subset(ExAll,idents=c(0:15))
ExSmall <- ExAllTmp[, sample(colnames(ExAllTmp), size = 10000, replace=F)] # Randomly selected 10,000 cells to identify the activated TFs
table(ExSmall$seurat_clusters)
pdf("ExSmallLabel.pdf",width=9) #Check the coverage of small set
DimPlot(ExSmall,label=TRUE,reduction = "tsne")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

hg38_dbs <- list('500bp'= 'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather', 
                 '10kb' = 'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather')
scenicOptions <- initializeScenic(org='hgnc',
                                datasetTitle='Rectal', 
                                dbDir="databases",
                                dbs=hg38_dbs,
                                nCores=30) 
exprMat  <-  as.matrix(ExSmall@assays$RNA@data)

cellInfo <-  ExSmall@meta.data
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,minCountsPerGene=3*.01*ncol(exprMat),minSamples=ncol(exprMat)*.01)
length(genesKept)
exprMat_filtered <- exprMat[genesKept, ]
runGenie3(exprMat_filtered, scenicOptions)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
scenicOptions <- initializeScenic(org="hgnc",dbDir="databases" , dbs=hg38_dbs, datasetTitle='ExIntegrated', nCores=20)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat) #(Non-detected genes are shuffled at the end of the ranking. Keep it in mind when choosing the threshold for calculating the AUC)
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)

cellInfo <-  ExSmall@meta.data
cellInfo <- data.frame(seuratCluster=ExSmall$seurat_clusters)
scenicOptions <- initializeScenic(org="hgnc",dbDir="databases" , dbs=hg38_dbs, datasetTitle='ExIntegrated', nCores=1) #reset the nCores with 1
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster), function(cells) rowMeans(getAUC(regulonAUC)[,cells]))


regulonActivity_byCellType=as.matrix(regulonActivity_byCellType)
regulonActivity_byCellType=regulonActivity_byCellType[Biobase::rowMax(regulonActivity_byCellType)>0.01,] #only focused on the TFs with auc socre greater than 0.01 in at least one sub cluster
pdf("regulonActivity_byCellType4Graph.pdf",height=4,width=4)
pheatmap(regulonActivity_byCellType,scale="row",clustering_method="ward.D2",color = colorRampPalette(c( "white","white","white","pink","red"))(50),show_rownames = T)
dev.off()

