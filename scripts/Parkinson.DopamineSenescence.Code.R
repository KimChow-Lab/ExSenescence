library(Seurat)
library(ggplot2)

setwd("/Projects/deng/Aging/DANeuron/Kamath_GSE178265")
GSE178265.data <- Read10X(data.dir = "rawData/")

MetaData=read.table("rawData/METADATA_PD.tsv",header=T,row.names=1) #https://singlecell.broadinstitute.org/single_cell/study/SCP1768/single-cell-genomic-profiling-of-human-dopamine-neurons-identifies-a-population-that-selectively-degenerates-in-parkinsons-disease-single-nuclei-data#study-download
MetaData=MetaData[,c("species","donor_id","disease__ontology_label","organ__ontology_label","sex","Donor_Age","Donor_PMI","Status","Cause_of_Death","FACS_Classification")]
sampleInfo=unique(MetaData[,c("donor_id","disease__ontology_label","Cause_of_Death","organ__ontology_label","sex","Donor_Age","Status")])
write.table(sampleInfo,file="sampleInfo.txt",quote=F,sep="\t",row.names=F)

cells=intersect(colnames(GSE178265.data),rownames(MetaData))
MetaData=MetaData[cells,]
HumanDANeuron <- CreateSeuratObject(counts = GSE178265.data, project = "HumanDANeuron")
all(rownames(MetaData)==colnames(GSE178265.data)) #TRUE

HumanDANeuron$Sex=MetaData$sex
HumanDANeuron$SampleId=MetaData$donor_id
HumanDANeuron$Age=MetaData$Donor_Age
HumanDANeuron$Status=MetaData$Status
HumanDANeuron$CauseDeath=MetaData$Cause_of_Death

NR4A2Pos=subset(HumanDANeuron,FACSGroup=="Positive") #only focused on the neurons
#obtain the DA neurons
NR4A2Pos <- NormalizeData(NR4A2Pos, normalization.method = "LogNormalize", scale.factor = 10000)
NR4A2Pos <- FindVariableFeatures(NR4A2Pos, selection.method = "vst", nfeatures = 2000)
NR4A2Pos <- ScaleData(NR4A2Pos)
NR4A2Pos <- RunPCA(NR4A2Pos, features = VariableFeatures(object = NR4A2Pos))
NR4A2Pos<- RunUMAP(NR4A2Pos, reduction = "pca", dims = 1:50)
NR4A2Pos<- FindNeighbors(NR4A2Pos, reduction = "pca", dims = 1:50)
NR4A2Pos<- FindClusters(NR4A2Pos, resolution = 0.5)
pdf("QC/NR4A2PosCluster.pdf",width=12,height=9)
DimPlot(NR4A2Pos,label=TRUE,label.size = 7)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
tiff("QC/NR4A2PosMarkerGene.tiff",width=1000,height=800)
FeaturePlot(NR4A2Pos,features=c("NRGN","SLC18A2","TH","RBFOX3","GAD1","AQP4","CX3CR1","VCAN","OLIG1","CLDN5"),raster=FALSE)&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
DANeuron=subset(NR4A2Pos,idents=c(0,8,24)) #resolution = 0.5
DANeuron <- RunPCA(DANeuron, features = VariableFeatures(object = DANeuron))
DANeuron <- RunTSNE(DANeuron, dims = 1:30)
DANeuron<- FindNeighbors(DANeuron, reduction = "pca", dims = 1:30)
DANeuron<- FindClusters(DANeuron, resolution = 0.5)
pdf("DANeuron/DANeuronCluster.pdf",width=11,height=9)
DimPlot(DANeuron,label=TRUE,reduction="tsne")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
#check the batch effect by sample
pdf("DANeuron/DANeuronSample.pdf",width=11,height=9)
DimPlot(DANeuron,group.by="SampleId",reduction="tsne")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
#================intergrated control and disease samples===================
DANeuron$Group=ifelse(DANeuron$Status%in%("Ctrl"),"Control","LBD_PD")
DANeuron.list <- SplitObject(DANeuron, split.by = "Group")
DANeuron.list <- lapply(X = DANeuron.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = DANeuron.list)
DANeuron.list <- lapply(X = DANeuron.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = DANeuron.list, dims = 1:30)
DANeuron.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
DANeuron.integrated <- ScaleData(DANeuron.integrated, verbose = FALSE)
DANeuron.integrated <- RunPCA(DANeuron.integrated, verbose = FALSE)
DANeuron.integrated <- RunTSNE(DANeuron.integrated, dims = 1:30)
DANeuron.integrated<- FindNeighbors(DANeuron.integrated, reduction = "pca", dims = 1:30)
DANeuron.integrated<- FindClusters(DANeuron.integrated, resolution = 0.2)

pdf("DANeuronIntegrated/DANeuronIntegratedCluster.pdf",width=11,height=9)
DimPlot(DANeuron.integrated,label=TRUE,label.size=7,reduction="tsne")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("DANeuronIntegrated/DANeuronIntegratedSample.pdf",width=11,height=9)
DimPlot(DANeuron.integrated,group.by="SampleId",reduction="tsne")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("DANeuronIntegrated/DANeuronIntegratedStatues.pdf",width=18,height=6)
DimPlot(DANeuron.integrated,split.by="Status",reduction="tsne")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("DANeuronIntegrated/DANeuronIntegratedGroup.pdf",width=14,height=6)
DimPlot(DANeuron.integrated,split.by="Group",reduction="tsne")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()


CCGene=read.table("CCGeneCombine.txt",header=T,sep="\t") #711
DefaultAssay(DANeuron.integrated)="RNA"
genes=intersect(rownames(DANeuron.integrated),CCGene$NAME)
length(genes) #364
CCGene=CCGene[CCGene$NAME%in%genes,]
CCGenelist=split(CCGene$NAME,CCGene$PHASE)
DANeuron.integrated <- AddModuleScore(
  object = DANeuron.integrated,
  features = CCGenelist,
  name = names(CCGenelist),
  ctrl = 100,
)
pdf("DANeuronIntegrated/DANeuronIntegrated.cellCycleScore.pdf",width=6)
VlnPlot(DANeuron.integrated,c("G1.S1","S.phase5","G22","G2.M3","M.G14"),ncol=1,pt.size=0)&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

cluster=c(0:8)
CycleScore=DANeuron.integrated@meta.data[,c("seurat_clusters","G1.S1","S.phase5","G22","G2.M3","M.G14")]
PvalueBetCluster=matrix(data = NA, nrow = length(cluster), ncol = length(cluster), dimnames = list(cluster, cluster))
for(i in 1:length(cluster)){
  for( j in 1:length(cluster)){
    C1_score=CycleScore[CycleScore$seurat_clusters==cluster[i],"G2.M3"]
    C2_score=CycleScore[CycleScore$seurat_clusters==cluster[j],"G2.M3"]
    delta=mean(C1_score)-mean(C2_score) #the row represent the current cluster
    pval=wilcox.test(C1_score,C2_score)$p.value
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
theme_classic()+
theme(legend.position="none")+
theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
pdf("DANeuronIntegrated/IntegratedDANeuron_Phase_PvalueDis/DANeuron_G2M_PvalueDisR0.2.pdf",height=2,width=4)
print(t)
dev.off()
colnames(PvalueBetCluster.df)=c("CurrentCluster","OtherClusters","-log10(pvalue).adjust")
write.table(PvalueBetCluster.df,file="IntegratedADCtrlCells_Phase_PvalueDis/IntegratedADCtrlCells_M.G14_PvalueDis.txt",sep="\t",quote=F,row.names=F)




coreCellularGene=read.table("55markersBySegura.txt",header=T)
DefaultAssay(DANeuron.integrated)="RNA"
DANeuronExpr <- AverageExpression(DANeuron.integrated)[["RNA"]]
hclust=hclust(as.dist(1-cor(DANeuronExpr)),method="ward.D2")
pdf("DANeuronIntegrated/hclust.pdf",height=5)
plot(hclust,hang=-1)
dev.off()
order=hclust$labels[hclust$order]
t=DotPlot(DANeuron.integrated,features=coreCellularGene[,1])
result=merge(t$data,coreCellularGene,by.x="features.plot",by.y="GeneNames")
C0=result[result$id==0,]
C0=C0[order(C0$Group,C0$pct.exp),]
geneOrder=C0$features.plot
g=ggplot(result, aes(factor(features.plot,levels=rev(unique(geneOrder))),factor(id,levels=order), size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red")+
labs(color="Average",size="Percent",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0),face="bold"),axis.text.x = element_text(size=rel(1.0),face="bold",angle=90),axis.text.y = element_text(size=rel(1.0),face="bold")) 
pdf("DANeuronIntegrated/SignatureCellularGene.pdf",width=10,height=5)
print(g)
dev.off()
result$features.plot=factor(result$features.plot,levels=geneOrder)
result$id=factor(result$id,levels=order)
result=result[order(result$features.plot,result$id),]
write.table(result,file="DANeuronIntegrated/SignatureCellularGeneExpr.txt",sep="\t",quote=F)

library("infercnv")
EPC=DANeuron.integrated
counts_matrix = as.matrix(EPC@assays$RNA@counts)
pData=EPC@meta.data$seurat_clusters
names(pData)=rownames(EPC@meta.data)
pData=data.frame(pData)
annotation=read.table("gencode.v38.annotation.geneType.txt",header=TRUE,row.names=1) #gene annotation
anno=data.frame("Chr"=annotation$Chr,"Start"=annotation$Start,"End"=annotation$End)
rownames(anno)=annotation$Symbol
genes=intersect(rownames(counts_matrix),rownames(anno))
anno=anno[genes,]
anno$Chr=factor(anno$Chr,levels=paste0("chr",c(1:22,"X","Y"),sep=""))
anno=anno[order(anno$Chr,anno$Star,anno$End),]
genes=intersect(rownames(anno),rownames(counts_matrix))
genes=genes[-which(genes%in%grep("MT-", genes,value=T))]
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
                               #out_dir="/Projects/deng/Aging/DANeuron/Kamath_GSE178265/DANeuronIntegrated/InferCNV"
                               out_dir="/data2/deng/Aging/Ex/DANeuron/InferCNVR02"
                               )

#------------------immune genes between PD and Ctrl----------------------------
DANeuron=DANeuron.integrated 
table(DANeuron$seurat_clusters) #validate the custers
DefaultAssay(DANeuron)="RNA"
DANeuron$Group=paste0(DANeuron$seurat_clusters,"_",DANeuron$Group,sep="")
Idents(DANeuron)=factor(DANeuron$Group,levels=c(paste0(c(0:8),"_Control",sep=""),paste0(c(0:8),"_LBD_PD",sep="")))
ImmuneAssociateGene=read.table("ImmuneAssociateGene.txt",header=T,sep="\t")
t=DotPlot(DANeuron,features=ImmuneAssociateGene$Symbol)
result=merge(t$data,ImmuneAssociateGene,by.x="features.plot",by.y="Symbol")
C0=result[result$id=="0_Control",] #order by their expression in cluster 0
C0=C0[order(C0$Group,C0$pct.exp),]
geneOrder=C0$features.plot
g=ggplot(result, aes(factor(features.plot,levels=unique(rev(geneOrder))), id, size= pct.exp,color=avg.exp.scaled)) +geom_point(alpha = 0.9)+
scale_colour_viridis_c()+
labs(color="avg.exp.scaled",size="pct.exp",x="",y="",title="")+
facet_grid(~Group, scale="free_x",space = "free")+
theme_bw()+
theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(angle=90,face = "italic",vjust=1,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("DANeuronIntegrated/Immune/ImmuneAssociateGeneInPDSplitVersion.pdf",width=15,height=5)
print(g)
dev.off()
result$features.plot=factor(result$features.plot,levels=rev(geneOrder))
result=result[order(result$features.plot),]
write.table(result,file="DANeuronIntegrated/Immune/ImmuneAssociateGene.Expr.txt",sep="\t",row.names=F,quote=F)


Idents(DANeuron)=DANeuron$seurat_clusters
G2M=c("CETN2","DCTN1","DYNC1H1","DYNC1I2","HSP90AA1","HSP90AB1","PAFAH1B1","PPP2CA","PPP2R1A","PRKACA","PSMA3","PSMA4","PSMA7","PSMB1","PSMB2","PSMB3","PSMB4","PSMB5","PSMB6","PSMB7","PSMC1","PSMC3","PSMC4","PSMC5","PSMC6","PSMD1","PSMD2","PSMD3","PSMD4","PSMD7","PSMD8","PSMD12","PSMD13","PSME1","RBBP4","RPS27A","SKP1","TUBA4A","TUBB2A","UBA52","UBB","UBC","YWHAE","YWHAG","TUBA1A","SSNA1","DYNLL1","PSMD6","RBX1","ACTR1A","OPTN","PSME3","PSMD14","TUBA1B","TUBB4A","TUBB4B","DCTN2","TUBGCP2","DCTN3","PPME1","MZT2B","TUBB","TUBB2B","MZT2A")
t=DotPlot(DANeuron,features=G2M)
result=t$data
C0=result[result$id==2,]
C0=C0[order(C0$pct.exp),]
geneOrder=C0$features.plot
g=ggplot(result, aes(factor(features.plot,levels=rev(unique(geneOrder))),id, size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red")+
labs(color="avg.exp.scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90,face = "italic",vjust=1,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("DANeuronIntegrated/Marker/G2MSignature.pdf",width=10,height=4)
print(g)
dev.off()
result$features.plot=factor(result$features.plot,levels=rev(geneOrder))
result=result[order(result$features.plot),]
write.table(result,file="DANeuronIntegrated/Marker/G2MSignature.Expr.txt",sep="\t",row.names=F,quote=F)


#---------------monocyte---------------------------------------------
library(monocle) #local monocle2
packageVersion("monocle")
library(Seurat)
setwd("/Projects/deng/Aging/DANeuron/Kamath_GSE178265")
DANeuron.integrated=readRDS("DANeuronIntegrated/DANeuron.integrated.rds")

table(DANeuron.integrated$seurat_clusters)
DANeuron.integrated.small <- DANeuron.integrated[, sample(colnames(DANeuron.integrated), size = 5000, replace=F)]

DefaultAssay(DANeuron.integrated.small)="RNA"
data  <-  as.matrix(DANeuron.integrated.small@assays$RNA@counts)
celldata <- as.data.frame(DANeuron.integrated.small@meta.data)
genedata <- as.data.frame(x = row.names(DANeuron.integrated.small), row.names = row.names(DANeuron.integrated.small))
colnames(genedata) <- "gene_short_name"
pd <- new("AnnotatedDataFrame", data = celldata)
fd <- new("AnnotatedDataFrame", data = genedata)
DANeuron.cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
DANeuron.cds <- estimateSizeFactors(DANeuron.cds)
DANeuron.cds <- estimateDispersions(DANeuron.cds)
DANeuron.cds <- detectGenes(DANeuron.cds, min_expr = 0.1)
head(fData(DANeuron.cds))
gc()
expressed_genes <- row.names(subset(fData(DANeuron.cds),num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(DANeuron.cds[expressed_genes,],fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
length(ordering_genes)
DANeuron.cds <- setOrderingFilter(DANeuron.cds, ordering_genes)
DANeuron.cds <- reduceDimension(DANeuron.cds, max_components = 2,method = 'DDRTree')
DANeuron.cds <- orderCells(DANeuron.cds,reverse=TRUE)

pData(DANeuron.cds)$cellType=ifelse(pData(DANeuron.cds)$seurat_clusters %in% c(2,5),paste0("C",pData(DANeuron.cds)$seurat_clusters,sep=""),"Others")
pData(DANeuron.cds)$cellType=factor(pData(DANeuron.cds)$cellType,levels=c("C2","C5","Others"))
pdf("DANeuronIntegrated/Monocle/DANeuronTrajectoryCellType.pdf",height=6,width=6)
plot_cell_trajectory(DANeuron.cds, color_by = "cellType",cell_size = 1.5)+
scale_color_manual(values=c("#93AA00","#00B9E3","Gainsboro"))+
theme(plot.title=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
#"#DB72FB", C7
pdf("DANeuronIntegrated/Monocle/DANeuronTrajectoryCellTypeXYTitle.pdf",height=6,width=6)
plot_cell_trajectory(DANeuron.cds, color_by = "cellType",cell_size = 1.5)+
scale_color_manual(values=c("#93AA00","#00B9E3","Gainsboro"))+
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

pdf("DANeuronIntegrated/Monocle/DANeuronTrajectory.pdf",height=6,width=6)
plot_cell_trajectory(DANeuron.cds, color_by = "seurat_clusters")+
theme(plot.title=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("DANeuronIntegrated/Monocle/DANeuronState.pdf",height=6,width=6)
plot_cell_trajectory(DANeuron.cds, color_by = "State",cell_size = 1)+
theme(plot.title=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
library(viridis)
pdf("DANeuronIntegrated/Monocle/DANeuronPseudotime.pdf",height=6,width=6)
plot_cell_trajectory(DANeuron.cds, color_by = "Pseudotime",cell_size = 1)+scale_color_viridis_c()+
theme(plot.title=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("DANeuronIntegrated/Monocle/DANeuronPseudotimeXYTitle.pdf",height=6,width=6)
plot_cell_trajectory(DANeuron.cds, color_by = "Pseudotime",cell_size = 1)+scale_color_viridis_c()+
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
saveRDS(DANeuron.cds,file="DANeuronIntegrated/Monocle/DANeuronIntegrated.cds.rds")
