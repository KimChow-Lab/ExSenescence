library(Seurat)
library(ggplot2)
library(clustree)

#######################################################################################
##############      Merge snRNA profiles from all samples       #######################
#######################################################################################

MathyEx=readRDS("/Projects/deng/Alzheimer/syn18485175/cellCycle/MathyAllData/AllEx8Mathy.rds") #prefrontal cortex, For AD, age>75 years old  BA10
MathyExCtrl=subset(MathyEx,Statues %in% "Control") #24 sample, 16952 cells
MathyExCtrl$Gender=MathyExCtrl$Gender
MathyExCtrl$Age=MathyExCtrl$Age
MathyExCtrl$Source="Mathy"
MathyExCtrl$orig.ident=MathyExCtrl$ProjectID
MathyExCtrl$Statues="Control"

MathyExAD=subset(MathyEx,Statues %in% "Alzheimer") #24 sample, 16952 cells
MathyExAD$Gender=MathyExAD$Gender
MathyExAD$Age=MathyExAD$Age
MathyExAD$Source="Mathy"
MathyExAD$orig.ident=MathyExAD$ProjectID
MathyExAD$Statues="Alzheimer"

YangEx=readRDS("/Projects/deng/Aging/Ex/Yang_GSE159812/YangIntegratedEx.rds") #Medial prefrontal cortex, For COVID, age>60 
YangExCtrl=subset(YangEx,Group %in% "Control") #8 samples, 8360 cells
tmp=ifelse(YangExCtrl$Gender == "F", "Female", "Male")
YangExCtrl$Gender=tmp
YangExCtrl$Age=YangExCtrl$Age
YangExCtrl$Source="Yang"
YangExCtrl$orig.ident=YangExCtrl$orig.ident
YangExCtrl$Statues="Control"

NagyEx=readRDS("/Projects/deng/Aging/Ex/Nagy_GSE144136/NagyEx.rds") #post-mortem dorsolateral prefrontal cortex (BA9), For MDD, age= 38.4(4)
NagyExCtrl=subset(NagyEx,Group %in% "Control") #17 samples, 18555 cells
NagyExCtrl$Gender="Male" #all male
NagyExCtrl$Age=38
NagyExCtrl$Source="Nagy"
NagyExCtrl$orig.ident=paste0("B",NagyExCtrl$sampleID,sep="")
NagyExCtrl$Statues="Control"

LauEx=readRDS("/Projects/deng/Alzheimer/syn18485175/cellCycle/LauAllData/LauEx.rds")#post-mortem dorsolateral prefrontal cortex (BA9), For AD, age= 60-95
LauExCtrl=subset(LauEx,CONDITION %in% "Control")#12 samples, 30215 cells
LauExCtrl$Gender=LauExCtrl$SEX
LauExCtrl$Age=LauExCtrl$AGE
LauExCtrl$Source="Lau"
LauExCtrl$orig.ident=LauExCtrl$orig.ident
LauExCtrl$Statues="Control"

LauEx=readRDS("/Projects/deng/Alzheimer/syn18485175/cellCycle/LauAllData/LauEx.rds")#post-mortem dorsolateral prefrontal cortex (BA9), For AD, age= 60-95
LauExAD=subset(LauEx,CONDITION %in% "Alzheimer")#
LauExAD$Gender=LauExAD$SEX
LauExAD$Age=LauExAD$AGE
LauExAD$Source="Lau"
LauExAD$orig.ident=LauExAD$orig.ident
LauExAD$Statues="Alzheimer"

setwd("/Projects/deng/Aging/Ex/AllEx")
ExAll=merge(MathyExCtrl,y=c(MathyExAD,LauExCtrl,LauExAD,YangExCtrl,NagyExCtrl))
table(paste0(ExAll$Source,"_",ExAll$Statues,sep=""))

saveRDS(ExAll,"ExAllMerge.rds")



#######################################################################################
##############   Integrated different datasets into a whole datasets      #############
#######################################################################################

ExAll=readRDS("ExAllMerge.rds")
table(paste0(ExAll$Source,"_",ExAll$Statues,sep=""))
ExAll[["percent.mt"]] <- PercentageFeatureSet(ExAll, pattern = "^MT-")
selected_mito <- WhichCells(ExAll, expression = percent.mt < 25)
ExAllTmp <- subset(ExAll, cells = selected_mito)
ExAllTmp <- ExAllTmp[!grepl("^MT-", rownames(ExAllTmp)), ]
length(selected_mito)
memory.limit(size=400000)
ExAll=ExAllTmp
ExAll.list <- SplitObject(ExAll, split.by = "Source")
ExAll.list <- lapply(X = ExAll.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = ExAll.list)
ExAll.list <- lapply(X = ExAll.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = ExAll.list, reference = c(1), reduction = "rpca",dims = 1:50)
ExAll.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
ExAll.integrated <- ScaleData(ExAll.integrated, verbose = FALSE)
ExAll.integrated <- RunPCA(ExAll.integrated, verbose = FALSE)
pdf("ExAll.integratedElbowPlot.pdf")
ElbowPlot(ExAll.integrated,ndims = 50)
dev.off()
# t-SNE and Clustering
ExAll.integrated <- RunTSNE(ExAll.integrated, dims = 1:50,check_duplicates = FALSE)
ExAll.integrated<- FindNeighbors(ExAll.integrated, reduction = "pca", dims = 1:30)


for(i in colnames(ExAll.integrated@meta.data)) {
  ExAll.integrated[[i]] <- NULL
}
ExAll.integrated$orig.ident=ExAll.integrated$orig.ident
ExAll.integrated$nCount_RNA=ExAll.integrated$nCount_RNA
ExAll.integrated$nFeature_RNA=ExAll.integrated$nFeature_RNA
ExAll.integrated$seurat_clusters=ExAll.integrated$seurat_clusters
ExAll.integrated$Gender=ExAll.integrated$Gender
ExAll.integrated$Age=ExAll.integrated$Age
ExAll.integrated$Source=ExAll.integrated$Source
ExAll.integrated$Statues=ExAll.integrated$Statues

#----confirmation the resolution------------------
ExAll.integratedTmp <- FindClusters(ExAll.integrated, resolution = c(seq(0.1,1,0.1)))
clus.tree.out <- clustree(ExAll.integratedTmp)
pdf(file = "ExAll.integrated.Resolution.Tree.pdf", width = 12, height = 10)
print(clus.tree.out)
dev.off()

ExAll.integrated<- FindClusters(ExAll.integrated, resolution = 0.5)
saveRDS(ExAll.integrated,"ExAll.integratedTSNE.rds")


#######################################################################################
##########    Initial analysis of the Ex contain both AD and control      #############
#######################################################################################
setwd("/Projects/deng/Aging/Ex/AllEx")
ExAll.integrated=readRDS("ExAll.integratedTSNE.rds")
ExAll.integrated=subset(ExAll.integrated,idents=c(0:15)) #remove cluster 16 by their smaller cell number

#-----Figure S8C---------------
tiff("ExAll.integratedLabel.tiff",width=300,height=270)
DimPlot(ExAll.integrated, raster=FALSE,reduction = "tsne")&NoLegend()&theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

#### cell number from different source in  each cluster####
#Figure S8D
t=as.matrix(table(Ex$Source))
tmp=paste0(Ex$seurat_clusters,"_",Ex$Statues,"-",Ex$Source,sep="")
tmp=as.matrix(table(tmp))
TmpInfo=as.matrix(do.call(rbind, strsplit(as.character(rownames(tmp)),'_')))
result=data.frame("Cluster"=TmpInfo[,1],"SourceStatues"=TmpInfo[,2],"Count"=tmp[,1])
ClusterList=factor(result$Cluster,levels=c(0:15))
g=ggplot(result, aes(ClusterList, Count, fill=SourceStatues)) +
  geom_bar(stat="identity",position="fill") +
  scale_y_continuous(expand=c(0,0))+
  scale_fill_brewer(palette="Set3")+
  theme(axis.text.x = element_text(angle = 0))+
  theme(axis.title.x = element_blank())
pdf("SourceStatuesDistribution.pdf",height=2.5,width=8)
print(g)
dev.off()
result=result[order(result$Cluster),]
write.table(result,file="cellNumberRatioInEachSource.txt",sep="\t",quote=F)



#-----Figure 4L---------------
#-----expression of immune associated genes---------------
ExTmp=subset(ExAll.integrated,idents=c(0:6)) #only show the top 7 clusters
table(ExTmp$seurat_clusters) #validate the custers
DefaultAssay(ExTmp)="RNA"
ExTmp$Group=paste0(ExTmp$seurat_clusters,"_",ExTmp$Statues,sep="")
Idents(ExTmp)=factor(ExTmp$Group,levels=c(paste0(c(0:6),"_Control",sep=""),paste0(c(0:6),"_Alzheimer",sep="")))
ImmuneAssociateGene=read.table("Immune/ImmuneAssociateGene.txt",header=T,sep="\t")
t=DotPlot(ExTmp,features=ImmuneAssociateGene$Symbol)
result=merge(t$data,ImmuneAssociateGene,by.x="features.plot",by.y="Symbol")
C0=result[result$id=="0_Control",] #order by their expression in cluster 0
C0=C0[order(C0$Group,C0$pct.exp),]
geneOrder=C0$features.plot
g=ggplot(result, aes(factor(features.plot,levels=unique(rev(geneOrder))), id, size= pct.exp,color=avg.exp.scaled)) +geom_point(alpha = 0.9)+
scale_colour_viridis_c()+
labs(color="avg.exp.scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+
facet_grid(~Group, scale="free_x",space = "free")+
theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(angle=90,face = "italic",vjust=1,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("Immune/ImmuneAssociateGeneSplitVersion.pdf",width=15,height=4)
print(g)
dev.off()

#figure 4B
#-----expression of immune associated genes---------------
immuneGene=c("BIN1","COX7C","HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPA2","HLA-DPA3","HLA-DPB1","HLA-DPB2","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DQB2","HLA-DQB3","HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DRB6","HLA-DRB9","HLA-E","HLA-F","HLA-G","HLA-H","HLA-J","HLA-K","HLA-L","HLA-N","HLA-P","HLA-S","HLA-T","HLA-U","HLA-V","HLA-W","HLA-Z","NYAP1","ZCWPW1","ECHDC3","RABEP1","IL34","SCIMP","KLF16","ABCA7","APOE","SIGLEC11","KAT8")
ExAllTmp=subset(ExAll.integrated,idents=c(0,3,5))
Idents(ExAllTmp)=factor(ExAllTmp$seurat_clusters,levels=c(0,3,5))
ExAllTmp$Statues=factor(ExAllTmp$Statues,levels=c("Control","Alzheimer"))
tiff("immune.Gene_VlnPlot.tiff",width=1800,height=1800)
VlnPlot(ExAllTmp,features=immuneGene,pt.size=0.01,split.by="Statues")&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))+theme(legend.position="right")
dev.off()
counts <- ExAllTmp@assays$RNA@counts
counts=data.frame(counts,check.names=F)
immuneGeneExpr=counts[intersect(immuneGene,rownames(ExAllTmp)),]
all(colnames(ExAllTmp)==colnames(immuneGeneExpr))
immuneGeneExprInfo=cbind(ExAllTmp@meta.data[,c("seurat_clusters","Gender")],t(immuneGeneExpr))
write.table(immuneGeneExprInfo,file="immune.Gene.expr.txt",sep="\t",quote=F)

####################################################################################
############ cell cycle phase score calculation for integrated cells  ##############
####################################################################################
setwd("/Projects/deng/Aging/Ex/AllEx")
ExAll.integrated=readRDS("ExAll.integratedTSNE.rds")
ExAll.integrated=subset(ExAll.integrated,idents=c(0:15)) #remove cluster 16 by their smaller cell number
CCGene=read.table("/Projects/deng/Aging/Ex/cellCycleGene/CCGeneCombine.txt",header=T,sep="\t") #711
DefaultAssay(ExAll.integrated)="RNA"
genes=intersect(rownames(ExAll.integrated),CCGene$NAME)
length(genes) #364
CCGene=CCGene[CCGene$NAME%in%genes,]
CCGenelist=split(CCGene$NAME,CCGene$PHASE)
ExAll.integrated <- AddModuleScore(
  object = ExAll.integrated,
  features = CCGenelist,
  name = names(CCGenelist),
  ctrl = 100,
)
pdf("IntegratedADCtrlCells_Phase_PvalueDis/ExAll.integrated.cellCycleScore.pdf",width=6)
VlnPlot(ExAll.integrated,c("G1.S1","S.phase5","G22","G2.M3","M.G14"),ncol=1,pt.size=0)&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

cluster=c(0:15)
CycleScore=ExAll.integrated@meta.data[,c("seurat_clusters","G1.S1","S.phase5","G22","G2.M3","M.G14")]
PvalueBetCluster=matrix(data = NA, nrow = length(cluster), ncol = length(cluster), dimnames = list(cluster, cluster))
for(i in 1:length(cluster)){
  for( j in 1:length(cluster)){
    C1_score=CycleScore[CycleScore$seurat_clusters==cluster[i],"G2.M3"]
    C2_score=CycleScore[CycleScore$seurat_clusters==cluster[j],"G2.M3"]
    delta=mean(C1_score)-mean(C2_score) #the row represent the current cluster
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
theme_classic()+
theme(legend.position="none")+
theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
pdf("IntegratedADCtrlCells_Phase_PvalueDis/IntegratedADCtrlCells_G2M_PvalueDis.pdf",height=2,width=6)
print(t)
dev.off()
colnames(PvalueBetCluster.df)=c("CurrentCluster","OtherClusters","-log10(pvalue).adjust")
write.table(PvalueBetCluster.df,file="IntegratedADCtrlCells_Phase_PvalueDis/IntegratedADCtrlCells_M.G14_PvalueDis.txt",sep="\t",quote=F,row.names=F)




#######################################################################################
#################################   infercnv      ######################
#######################################################################################
setwd("/Projects/deng/Aging/Ex/AllEx/InferCNV")
library("infercnv")
ExAll=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
table(ExAll$seurat_clusters)
ExSmall <- ExAll[, sample(colnames(ExAll), size = 10000, replace=F)]
table(ExSmall$seurat_clusters)
saveRDS(ExSmall,"/data2/deng/Aging/Ex/AllEx/InferCNV/ExSmallSetCells.rds")

ExSmall=readRDS("/data2/deng/Aging/Ex/AllEx/InferCNV/ExSmallSetCells.rds")
ExSmall=subset(ExSmall,idents=c(0:15))
pdf("/data2/deng/Aging/Ex/AllEx/ExAllSmall.integratedLabel.pdf",width=9)
DimPlot(ExSmall, reduction = "tsne", label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
counts_matrix = as.matrix(ExSmall@assays$RNA@counts)
pData=ExSmall@meta.data$seurat_clusters
names(pData)=rownames(ExSmall@meta.data)
pData=data.frame(pData)
annotation=read.table("/Projects/deng/refGene/gencode.v38.annotation.removeDupliSymbol.geneType.txt",header=TRUE,row.names=1)
anno=data.frame("Chr"=annotation$Chr,"Start"=annotation$Start,"End"=annotation$End)
rownames(anno)=annotation$Symbol
rownames(anno)=annotation$Symbol
genes=intersect(rownames(counts_matrix),rownames(anno))
anno=anno[genes,]
anno$Chr=factor(anno$Chr,levels=paste0("chr",c(1:22,"X","Y"),sep=""))
anno=anno[order(anno$Chr,anno$Star,anno$End),]
genes=intersect(rownames(anno),rownames(counts_matrix))
length(genes)
gc()
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
                               out_dir="/data2/deng/Aging/Ex/AllEx/InferCNV/InferCNV"
                               )



#######################################################################################
############################### Deg between AD and Normal in cell cycling ##############
#######################################################################################
setwd("/Projects/deng/Aging/Ex/AllEx/DEGBetADNorInCycing")
ExAll.integrated=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
DefaultAssay(ExAll.integrated)="RNA"
ExAll.integrated=subset(ExAll.integrated,idents=c(0:15))
ExAllTmp=subset(ExAll.integrated,Age>60)
Idents(ExAllTmp)="Statues"
table(ExAllTmp$Source)
#Source=subset(ExAllTmp,Source=="Yang") #Lau Mathy 
for(i in c(0:15)){
  SubCluster=subset(ExAllTmp,seurat_clusters==i)
  ADDeg=FindMarkers(SubCluster,ident.1="Alzheimer",ident.2="Control",pct.min=0.1)
  ADDeg=ADDeg[ADDeg$p_val<0.01,]
  write.table(ADDeg,file=paste0("/Projects/deng/Aging/Ex/AllEx/DEGBetADNorInCycing/ADDeg4Old/AllSource/C",i,".txt",sep=""),sep="\t",quote=F)
}

#------------------------count the DEGs number in each cluster between AD and Ctrl----------------------------
ADDeg4Old=read.table("D:/Aging/Ex/AllEx/DEGBetADNorInCycing/ADDeg4Old/ADDegInAllSource.txt",header=T)
ADDeg4Old=ADDeg4Old[ADDeg4Old$p_val_adj<0.01,]
ADDeg4OldCount=data.frame(table(ADDeg4Old$Cluster,ADDeg4Old$Pattern))
colnames(ADDeg4OldCount)=c("Cluster","Pattern","Count")
#ADDeg4OldCount=ADDeg4OldCount[ADDeg4OldCount$Cluster%in%c("C0","C3","C5"),]
ClusterOrder=factor(ADDeg4OldCount$Cluster,levels=paste0("C",c(0:15),sep=""))
g=ggplot(ADDeg4OldCount, aes(ClusterOrder, Count, fill=factor(Pattern,levels=c("Up","Down")))) +
  geom_bar(stat="identity",position=position_dodge()) +
  scale_y_continuous(expand=c(0,0))+theme_bw()+
  scale_fill_manual(values=c("Violet","CornflowerBlue"))+
  guides(fill = guide_legend(title = "Pattern", title.position = "top"),col = guide_legend(nrow = 1))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
pdf("D:/Aging/Ex/AllEx/DEGBetADNorInCycing/DEGbetADInOld.pdf",height=4,width=8)
print(g)
dev.off()

write.table(ADDeg4OldCount,file="D:/Aging/Ex/AllEx/DEGBetADNorInCycing/DEGNumberbetADInOld.txt",sep="\t",quote=F)

pdf("D:/Aging/Ex/AllEx/DEGBetADNorInCycing/C0C3C5DEGbetADInOld.pdf",height=3,width=5)
print(g)
dev.off()
