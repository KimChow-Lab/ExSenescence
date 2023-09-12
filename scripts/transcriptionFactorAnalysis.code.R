library(Seurat)
library(SeuratDisk)
library(SeuratData)

ExAll.integrated=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
ExAll.integrated=subset(ExAll.integrated,idents=c(0:15))

ExAll.integrated$SourceStatus=paste0(ExAll.integrated$Source,"_",ExAll.integrated$Statues,sep="")
table(ExAll.integrated$SourceStatus)
DefaultAssay(ExAll.integrated)="RNA"
unique(ExAll.integrated$SourceStatus)

selected_f_10XGenomice <- rownames(ExAll.integrated)[Matrix::rowSums(ExAll.integrated) > ncol(ExAll.integrated)*0.001] #ncol(ExAll.integrated)*0.001=123
length(selected_f_10XGenomice)

#Extract the protein coding genes
UniqueGene=read.table("/Projects/deng/Public/Ensemble/gencode.v38.annotation.removeDupliSymbol.geneType.txt",header=T)
Protein_coding_Gene=UniqueGene[UniqueGene$Type%in%"protein_coding",] #19920     7
GeneList=intersect(Protein_coding_Gene$Symbol,rownames(ExAll.integrated)) #18509
length(GeneList) #18509

GeneRetain=intersect(GeneList,selected_f_10XGenomice) #16141 remained
ExAll.integrated <- subset(ExAll.integrated, features = GeneRetain)

#remove the mitochodrial, ribosome and noncoding gene
ExAll.integrated <- ExAll.integrated[!grepl("^MT-", rownames(ExAll.integrated)), ] #remove mitochdrial genes
ExAll.integrated <- ExAll.integrated[!grepl("^HB[^(P)]", rownames(ExAll.integrated)), ] #remove homoglobin genes
dim(ExAll.integrated) #  16133 123212

#https://www.nature.com/articles/nmeth.4463: SCENIC is robust to down-sampling of cells and sparse expression matrices.
setwd("/data2/deng/Aging/Ex/AllEx/SCENIC")
for(s in unique(ExAll.integrated$SourceStatus)){ #
  setwd(paste0("/data2/deng/Aging/Ex/AllEx/SCENIC/",s,"/"))
  Source.seurat.Downsample=subset(ExAll.integrated,SourceStatus==s)
  DefaultAssay(Source.seurat.Downsample)="RNA"
  if(ncol(Source.seurat.Downsample)>10000){
     Source.seurat.Downsample=Source.seurat.Downsample[, sample(colnames(Source.seurat.Downsample), size = 10000, replace=F)]
  }
    print(paste0(s,"  Raw:",ncol(Source.seurat.Downsample),"  DownSampling:",ncol(Source.seurat.Downsample),sep=""))
    Source.seurat.Downsample.loom <- as.loom(Source.seurat.Downsample, filename = paste0(s,".downsampling.loom",sep=""), verbose = FALSE)
    Source.seurat.Downsample.loom$close_all()	
}
#"Mathy  Raw:34799  DownSampling:10000"
#"Lau  Raw:61939  DownSampling:10000"
#"Yang  Raw:8084  DownSampling:8084"
#"Nagy  Raw:18390  DownSampling:10000"
ExAll.integrated.Downsample=ExAll.integrated[, sample(colnames(ExAll.integrated), size = 30000, replace=F)]
ExAll.integrated.loom <- as.loom(ExAll.integrated.Downsample, filename = "ExAll.integrated.downsampling.loom", verbose = FALSE)
ExAll.integrated.loom$close_all()

#using the pyscenic
#raw D:\Aging\Glia\scRNAseq\SCENIC.code.R
#create a tmux env
tmux new -s scenic
tmux attach-session -t scenic
conda activate pyscenic 
pyscenic -h

cd /data2/deng/Aging/Ex/AllEx/SCENICSourceStatus

source=ExAll.integrated

for source in Lau_Control Nagy_Control Yang_Control Mathy_Alzheimer Lau_Alzheimer
do
echo ${source}

pyscenic grn \
${source}.downsampling.loom \
--num_workers 5 \
--output ${source}.downsampling.tsv \
--method grnboost2 \
/Projects/deng/Public/SCENIC/allTFs_hg38.txt

pyscenic ctx \
${source}.downsampling.tsv \
/Projects/deng/Public/SCENIC/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather //Projects/deng/Public/SCENIC/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
--annotations_fname /Projects/deng/Aging/Ex/AllEx/SCENIC/databases/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname ${source}.downsampling.loom \
--transpose \
--mask_dropouts \
--mode "dask_multiprocessing" \
--output ${source}.reg.csv \
--num_workers 5

pyscenic aucell \
${source}.downsampling.loom \
${source}.reg.csv \
--output ${source}.downsampling.final.loom \
--num_workers 5

done




###### Define the cell type specific TFs ####################################
rm(list = ls())  
setwd("/data2/deng/Aging/Ex/AllEx/SCENICSourceStatus")
library(SCENIC)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(Seurat)
packageVersion("SCENIC")  
library(SCopeLoomR)
ExAll=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
DefaultAssay(ExAll)="RNA"
ExAll=subset(ExAll,idents=c(0:15))
ExAll[["Layer"]]=NA
ExAll[["Layer"]][ExAll$seurat_clusters%in%c(0,2,8),]="L2/3"
ExAll[["Layer"]][ExAll$seurat_clusters%in%c(1,6,7,9,13),]="L4/5"
ExAll[["Layer"]][ExAll$seurat_clusters%in%c(4,10),]="L4/6"
ExAll[["Layer"]][ExAll$seurat_clusters%in%c(15),]="L5"
ExAll[["Layer"]][ExAll$seurat_clusters%in%c(11,14),]="L5/6"
ExAll[["Layer"]][ExAll$seurat_clusters%in%c(12),]="L6"
ExAll[["Layer"]][ExAll$seurat_clusters%in%c(3),]="ES"
ExAll[["Layer"]][ExAll$seurat_clusters%in%c(5),]="LS"
ExAll$Layer=factor(ExAll$Layer,levels=c("ES","LS","L2/3","L4/5","L4/6","L5","L5/6","L6"))
ExAll$SourceStatus=paste0(ExAll$Source,"_",ExAll$Statues,sep="")



s="ExAll.integrated"
scenicLoomPath=paste0(s,".downsampling.final.loom",sep="")
loom <- open_loom(scenicLoomPath)
# Read information from loom file: 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC") 
dim(regulonAUC)

subSource=subset(ExAll,cells = colnames(regulonAUC)) #for ExAll.integrated
selected_f_10xGenomic <- rownames(subSource)[Matrix::rowSums(subSource) > ncol(subSource)*0.1] #expressed in at lest 10% cells
print(paste0(s,":Expressed Genes ",length(selected_f_10xGenomic)))
#7665
cellInfo=data.frame(clusters=subSource$seurat_clusters,Layers=subSource$Layer,Gender=subSource$Gender,Age=subSource$Age,Status=subSource$Statues)
clusterNumber=table(cellInfo$Layers)
cellInfo=cellInfo[cellInfo$Layers%in%names(clusterNumber[clusterNumber>10]),]
#https://www.jianshu.com/p/7ab2d6c8f764
print(paste0(s,":TFNumber ",dim(regulonAUC)[1]))
#507 10000
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
cellId=intersect(rownames(cellInfo),colnames(regulonAUC))
regulonAUC <- regulonAUC[,cellId]
cellInfo=cellInfo[cellId,]
all(rownames(cellInfo)==colnames(regulonAUC))
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$Layers), function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType <- regulonActivity_byCellType[,unique(cellInfo$Layers)]
dim(regulonActivity_byCellType)
regulonActivity_byCellType=regulonActivity_byCellType[Biobase::rowMax(regulonActivity_byCellType)>0.01,] #only focused on the TFs with auc socre greater than 0.01 in at least one sub cluster
regulonActivity_byCellType.df=data.frame(regulonActivity_byCellType)
regulonActivity_byCellType.df$TFName=stringr::str_sub(rownames(regulonActivity_byCellType.df), end=-4)
regulonActivity_byCellType.df=regulonActivity_byCellType.df[regulonActivity_byCellType.df$TF%in%selected_f_10xGenomic,]
regulonActivity_byCellType=data.frame(regulonActivity_byCellType,check.names=F)
regulonActivity_byCellType=regulonActivity_byCellType[rownames(regulonActivity_byCellType.df),]
print(paste0(s,":TFNumber_AfterFilter: ",dim(regulonActivity_byCellType)[1]))
dim(regulonActivity_byCellType)
pdf(paste0("Heatmap/",s,"_FullName.pdf",sep=""),height=20)
t=pheatmap(regulonActivity_byCellType,scale="row",clustering_method="ward.D2")
dev.off()

TFGroup=data.frame(Cluster=cutree(t$tree_row,k=5))
TFGroup$Cluster=paste0("C",TFGroup$Cluster,sep="")
TFGroup$Pattern=ifelse(TFGroup$Cluster%in%c("C1","C2","C4"),"Activated","Inhibited")

Cluster_colors = colorRampPalette(brewer.pal(8, "Set2"))(length(unique(TFGroup$Cluster))) 
names(Cluster_colors)=unique(TFGroup$Cluster)

ann_colors=list(Cluster=Cluster_colors,Pattern=c(Activated="Orange",Inhibited="RoyalBlue"))

pdf(paste0("Heatmap/",s,"_Show.pdf",sep=""),height=6,width=5)
pheatmap(regulonActivity_byCellType,scale="row",clustering_method="ward.D2",show_rownames=F,annotation_row=TFGroup,color=viridis(10),annotation_colors = ann_colors)
dev.off()
pdf(paste0("Heatmap/",s,"_ShowGenes.pdf",sep=""),height=15,width=6)
pheatmap(regulonActivity_byCellType,scale="row",clustering_method="ward.D2",fontsize = 6,show_rownames=T,annotation_row=TFGroup,color=viridis(10),annotation_colors = ann_colors)
dev.off()

all(rownames(regulonActivity_byCellType)==rownames(TFGroup))
regulonActivity_byCellType$Group=TFGroup$Cluster
regulonActivity_byCellType$Pattern=TFGroup$Pattern
tmp=t$tree_row
TFOrder=tmp$labels[tmp$order]
regulonActivity_byCellType=regulonActivity_byCellType[TFOrder,]
write.table(regulonActivity_byCellType,file=paste0("Heatmap/",s,".txt",sep=""),sep="\t",quote=F)


###### TF list for each dataset #############
Integrated=read.table("Heatmap/ExAll.integrated.txt",header=T)
Activate_TFs=rownames(Integrated[Integrated$Pattern%in%c("Activated"),])
Inhibite_TFs=rownames(Integrated[Integrated$Pattern%in%c("Inhibited"),])

s="Mathy_Control"
unique(ExAll$SourceStatus)
for(s in unique(ExAll$SourceStatus)){
scenicLoomPath=paste0(s,".downsampling.final.loom",sep="")
loom <- open_loom(scenicLoomPath)
# Read information from loom file: 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC") 
dim(regulonAUC)

ExAll$SourceStatus=paste0(ExAll$Source,"_",ExAll$Statues,sep="")
subSource=subset(ExAll,SourceStatus==s)

selected_f_10xGenomic <- rownames(subSource)[Matrix::rowSums(subSource) > ncol(subSource)*0.1] #expressed in at lest 10% cells
print(paste0(s,":Expressed Genes ",length(selected_f_10xGenomic)))
#7665

cellInfo=data.frame(clusters=subSource$seurat_clusters,Layers=subSource$Layer,Gender=subSource$Gender,Age=subSource$Age,Status=subSource$Statues)
clusterNumber=table(cellInfo$Layers)
cellInfo=cellInfo[cellInfo$Layers%in%names(clusterNumber[clusterNumber>10]),]

#https://www.jianshu.com/p/7ab2d6c8f764
print(paste0(s,":TFNumber ",dim(regulonAUC)[1]))
#507 10000
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
cellId=intersect(rownames(cellInfo),colnames(regulonAUC))
regulonAUC <- regulonAUC[,cellId]
cellInfo=cellInfo[cellId,]
all(rownames(cellInfo)==colnames(regulonAUC))
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$Layers), function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType <- regulonActivity_byCellType[,unique(cellInfo$Layers)]
dim(regulonActivity_byCellType)
regulonActivity_byCellType=regulonActivity_byCellType[Biobase::rowMax(regulonActivity_byCellType)>0.01,] #only focused on the TFs with auc socre greater than 0.01 in at least one sub cluster
regulonActivity_byCellType.df=data.frame(regulonActivity_byCellType)
regulonActivity_byCellType.df$TFName=stringr::str_sub(rownames(regulonActivity_byCellType.df), end=-4)
regulonActivity_byCellType.df=regulonActivity_byCellType.df[regulonActivity_byCellType.df$TF%in%selected_f_10xGenomic,]
regulonActivity_byCellType=data.frame(regulonActivity_byCellType,check.names=F)
regulonActivity_byCellType=regulonActivity_byCellType[rownames(regulonActivity_byCellType.df),]

print(paste0(s,":TFNumber_AfterFilter: ",dim(regulonActivity_byCellType)[1]))
dim(regulonActivity_byCellType)
pdf(paste0("Heatmap/",s,"_FullName.pdf",sep=""),height=20)
t=pheatmap(regulonActivity_byCellType,scale="row",clustering_method="ward.D2")
dev.off()

TFGroup=data.frame(Cluster=cutree(t$tree_row,k=5))
TFGroup$Cluster=paste0("C",TFGroup$Cluster,sep="")
TFGroup$Pattern=ifelse(rownames(TFGroup)%in%Activate_TFs,"Activate_TFs",ifelse(rownames(TFGroup)%in%Inhibite_TFs,"Inhibite_TFs","Novel"))

Cluster_colors = colorRampPalette(brewer.pal(8, "Set2"))(length(unique(TFGroup$Cluster))) 
names(Cluster_colors)=unique(TFGroup$Cluster)
ann_colors=list(Cluster=Cluster_colors,Pattern=c(Activate_TFs="Orange",Inhibite_TFs="RoyalBlue",Novel="LightGrey"))

pdf(paste0("Heatmap/",s,"_Show.pdf",sep=""),height=6,width=5)
pheatmap(regulonActivity_byCellType,scale="row",clustering_method="ward.D2",show_rownames=F,annotation_row=TFGroup,color=viridis(10),annotation_colors = ann_colors)
dev.off()

all(rownames(regulonActivity_byCellType)==rownames(TFGroup))
regulonActivity_byCellType$Group=TFGroup$Cluster
regulonActivity_byCellType$Pattern=TFGroup$Pattern
tmp=t$tree_row
TFOrder=tmp$labels[tmp$order]
regulonActivity_byCellType=regulonActivity_byCellType[TFOrder,]
write.table(regulonActivity_byCellType,file=paste0("Heatmap/",s,".txt",sep=""),sep="\t",quote=F)

}

[1] "Mathy_Control:Expressed Genes 8107"
[1] "Mathy_Control:TFNumber 507"
[1] "Mathy_Control:TFNumber_AfterFilter: 182"
[1] "Mathy_Alzheimer:Expressed Genes 7338"
[1] "Mathy_Alzheimer:TFNumber 500"
[1] "Mathy_Alzheimer:TFNumber_AfterFilter: 167"
[1] "Lau_Control:Expressed Genes 7563"
[1] "Lau_Control:TFNumber 517"
[1] "Lau_Control:TFNumber_AfterFilter: 185"
[1] "Lau_Alzheimer:Expressed Genes 7483"
[1] "Lau_Alzheimer:TFNumber 508"
[1] "Lau_Alzheimer:TFNumber_AfterFilter: 180"
[1] "Yang_Control:Expressed Genes 7216"
[1] "Yang_Control:TFNumber 496"
[1] "Yang_Control:TFNumber_AfterFilter: 165"
[1] "Nagy_Control:Expressed Genes 7665"
[1] "Nagy_Control:TFNumber 518"
[1] "Nagy_Control:TFNumber_AfterFilter: 168"


###Deine the conserved TFs 
MathysCtrl=read.table("Heatmap/Mathy_Control.txt",header=T)
Activate_MathysCtrl=rownames(MathysCtrl[MathysCtrl$Group%in%c("C1","C4"),])
Inhibite_MathysCtrl=rownames(MathysCtrl[MathysCtrl$Group%in%c("C2","LS","C5"),])

LauCtrl=read.table("Heatmap/Lau_Control.txt",header=T)
Activate_LauCtrl=rownames(LauCtrl[LauCtrl$Group%in%c("C1","LS","C4"),])
Inhibite_LauCtrl=rownames(LauCtrl[LauCtrl$Group%in%c("C2","C5"),])

NagyCtrl=read.table("Heatmap/Nagy_Control.txt",header=T)
Activate_NagyCtrl=rownames(NagyCtrl[NagyCtrl$Group%in%c("C1","C4"),])
Inhibite_NagyCtrl=rownames(NagyCtrl[NagyCtrl$Group%in%c("C2","LS","C5"),])

YangCtrl=read.table("Heatmap/Yang_Control.txt",header=T)
Activate_YangCtrl=rownames(YangCtrl[YangCtrl$Group%in%c("C2","C4","C5"),])
Inhibite_YangCtrl=rownames(YangCtrl[YangCtrl$Group%in%c("C1","LS"),])

MathyAD=read.table("Heatmap/Mathy_Alzheimer.txt",header=T)
Activate_MathyAD=rownames(MathyAD[MathyAD$Group%in%c("C1","C2","C4"),])
Inhibite_MathyAD=rownames(MathyAD[MathyAD$Group%in%c("LS","C5"),])

LauAD=read.table("Heatmap/Lau_Alzheimer.txt",header=T)
Activate_LauAD=rownames(LauAD[LauAD$Group%in%c("C2","C4","C5"),])
Inhibite_LauAD=rownames(LauAD[LauAD$Group%in%c("C1","LS"),])



ActivateList=c(Activate_MathysCtrl,Activate_LauCtrl,Activate_NagyCtrl,Activate_YangCtrl,Activate_MathyAD,Activate_LauAD)
InhibiteList=c(Inhibite_MathysCtrl,Inhibite_LauCtrl,Inhibite_NagyCtrl,Inhibite_YangCtrl,Inhibite_MathyAD,Inhibite_LauAD)
ActivateTF.df=data.frame(table(ActivateList))
colnames(ActivateTF.df)=c("Symbol","ActivatedSupport")

InhibiteTF.df=data.frame(table(InhibiteList))
colnames(InhibiteTF.df)=c("Symbol","InhibitedSupport")

#barplot of the conserved TFS
IntegratedInfo=read.table("Heatmap/ExAll.integrated.txt",header=T)
IntegratedInfo=IntegratedInfo[,c("Group","Pattern")]
IntegratedInfo$Symbol=rownames(IntegratedInfo)
IntegratedInfo.df=merge(IntegratedInfo,ActivateTF.df,by="Symbol",all.x = TRUE,sort=FALSE)
IntegratedInfo.df=merge(IntegratedInfo.df,InhibiteTF.df,by="Symbol",all.x = TRUE,sort=FALSE)
rownames(IntegratedInfo.df)=IntegratedInfo.df$Symbol
IntegratedInfo.df=IntegratedInfo.df[rownames(IntegratedInfo),]
all(rownames(IntegratedInfo.df)==rownames(IntegratedInfo))
#TRUE
IntegratedInfo.ldf=reshape2::melt(IntegratedInfo.df,id=c(1:3))
IntegratedInfo.ldf$Symbol=factor(IntegratedInfo.ldf$Symbol,levels=rev(rownames(IntegratedInfo)))
g=ggplot(IntegratedInfo.ldf, aes(Symbol, value, fill=factor(Pattern,levels=c("Activated","Inhibited")))) +
  geom_bar(stat="identity") +
  scale_y_continuous(expand=c(0.1,0.1))+theme_bw()+
  scale_fill_manual(values=c("Orange","RoyalBlue"))+
  facet_grid(.~variable)+
  coord_flip()+
  guides(fill = guide_legend(title = "Pattern", title.position = "top"),col = guide_legend(nrow = 1))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
pdf("Heatmap/SuppotDatasets.pdf",height=20,width=7)
print(g)
dev.off()

#output the conserved TFS
IntegratedInfo.df[is.na(IntegratedInfo.df)] <- 0
IntegratedInfo.df.activated=IntegratedInfo.df[IntegratedInfo.df$ActivatedSupport>=4&IntegratedInfo.df$InhibitedSupport<=2&IntegratedInfo.df$Pattern=="Activated",]
IntegratedInfo.df.inhibited=IntegratedInfo.df[IntegratedInfo.df$InhibitedSupport>=4&IntegratedInfo.df$ActivatedSupport<=2&IntegratedInfo.df$Pattern=="Inhibited",]
IntegratedInfo.df.activated=IntegratedInfo.df.activated[!is.na(IntegratedInfo.df.activated$Symbol),]
IntegratedInfo.df.inhibited=IntegratedInfo.df.inhibited[!is.na(IntegratedInfo.df.inhibited$Symbol),]
dim(IntegratedInfo.df.activated)
#69  5
dim(IntegratedInfo.df.inhibited)
#23  5
intersect(rownames(IntegratedInfo.df.activated),rownames(IntegratedInfo.df.inhibited))
write.table(rbind(IntegratedInfo.df.activated,IntegratedInfo.df.inhibited),file="Heatmap/conservedTFs.txt",sep="\t",quote=F)



#targets of the conserved TFS
s="ExAll.integrated"
scenicLoomPath=paste0(s,".downsampling.final.loom",sep="")
loom <- open_loom(scenicLoomPath)
# Read information from loom file: 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
#myTFlist=c(rownames(IntegratedInfo.df.activated),rownames(IntegratedInfo.df.inhibited))
myTFlist=IntegratedInfo.df$Symbol
myTFTarget.df=data.frame(regulons_incidMat[myTFlist,],check.names=F)
myTFTarget.df$Symbol=rownames(myTFTarget.df)
myTFTarget.pattern=merge(unique(IntegratedInfo.df[,c("Symbol","Pattern")]),myTFTarget.df,by.x="Symbol")

myTFTarget.pattern.ldf=reshape2::melt(myTFTarget.pattern,id=c(1:2))
colnames(myTFTarget.pattern.ldf)=c("TFName","TFStatus","TFTarget","Target")
myTFTarget.pattern.ldf=myTFTarget.pattern.ldf[myTFTarget.pattern.ldf$Target==1,]
myTFTarget.pattern.ldf$Symbol=stringr::str_sub(myTFTarget.pattern.ldf$TFName,end=-4)

myTFTarget.pattern.ldf$TFName=factor(myTFTarget.pattern.ldf$TFName,levels=rownames(IntegratedInfo.df))
myTFTarget.pattern.ldf=myTFTarget.pattern.ldf[order(myTFTarget.pattern.ldf$TFName),]

write.table(myTFTarget.pattern.ldf,file="AllTFTargetGene.txt",row.names=F,sep="\t",quote=F)


#=========== TFs effect on DEG ===============================
LSMarker=read.table("/data2/deng/Aging/Ex/AllEx/C5Feature/C5MarkerFC1_MAST.txt",header=T)
LSMarkerList=split(rownames(LSMarker),LSMarker$Pattern)
names(LSMarkerList)
#"DnInLS" "UpInLS"

C3C5Marker=read.table("/data2/deng/Aging/Ex/AllEx/C5Feature/C3C5Marker_MAST.txt",header=T)
C3C5MarkerList=split(rownames(C3C5Marker),C3C5Marker$Pattern)
names(C3C5MarkerList)
#"DnInSene" "UpInSene"

MarkerList=LSMarkerList
names(MarkerList)=c("Down","Up")

TFTarget=read.table("/data2/deng/Aging/Ex/AllEx/SCENICSourceStatus/AllTFTargetGene.txt",header=T)
TFlist=unique(TFTarget$Symbol)
length(TFlist)
#162 TFs



TotalGeneNumber=dim(regulons_incidMat)[2]
df=matrix(data = NA, nrow = length(TFlist), ncol = 6, byrow = FALSE, dimnames = NULL)
for(i in 1:length(TFlist)){
    print(TFlist[i])
    SingleTFTargetList=TFTarget[TFTarget$Symbol%in%TFlist[i],"TFTarget"]
    UpTarget=intersect(SingleTFTargetList,MarkerList$Up)
    DownTarget=intersect(SingleTFTargetList,MarkerList$Dn)
    df[i,1]=length(SingleTFTargetList)
    df[i,2]=TotalGeneNumber
    df[i,3]=length(MarkerList$Up)
    df[i,4]=length(UpTarget)
    df[i,5]=length(MarkerList$Dn)
    df[i,6]=length(DownTarget)
}
df=data.frame(df)
colnames(df)=c("TFTargetN","TotalN","Up","TargetUp","Dn","TargetDn")
rownames(df)=TFlist

#1 - phyper(overlap - 1, sampleb, totala - sampleb, samplec)
df$UpP=1-phyper(df$TargetUp-1, df$TFTargetN, df$TotalN-df$TFTargetN, df$Up, lower.tail = TRUE)
df$DnP=1-phyper(df$TargetDn-1, df$TFTargetN, df$TotalN-df$TFTargetN, df$Dn, lower.tail = TRUE)

df$Symbol=rownames(df)
df.status=merge(unique(TFTarget[,c("Symbol","TFStatus")]),df,by="Symbol")
df.status$Symbol=factor(df.status$Symbol,levels=unique(TFTarget$Symbol))
df.status=df.status[order(df.status$Symbol),]
write.table(df.status,file="C5LSMarkerTarget8TF_Pvalue.txt",sep="\t",quote=F,row.names=F)


p<-ggplot(data=df.status, aes(x=factor(Symbol,levels=rev(Symbol)), y=-log10(UpP+1e-16),fill=-log10(UpP+1e-16))) +
  scale_fill_gradient(low = "Wheat", high = "DarkOrange")+
  theme_bw()+
  geom_bar(stat="identity",width=0.6)+
  coord_flip()
pdf("UpGeneTarget8TFPvalue_C5LS.pdf",height=15,width=4)
print(p)
dev.off()


p<-ggplot(data=df.status, aes(x=factor(Symbol,levels=rev(Symbol)), y=-log10(DnP+1e-16),fill=-log10(DnP+1e-16))) +
  scale_fill_gradient(low = "LightSteelBlue", high = "RoyalBlue")+ 
  theme_bw()+
  geom_bar(stat="identity",width=0.6)+
  coord_flip()
pdf("DnGeneTarget8TFPvalue_C5LS.pdf",height=15,width=4)
print(p)
dev.off()










#=========== conserved TFs across sub source ===============================
ActivateList=list(MathysCtrl=Activate_MathysCtrl,
             LauCtrl=Activate_LauCtrl,
             NagyCtrl=Activate_NagyCtrl,
             YangCtrl=Activate_YangCtrl,
             MathysAD=Activate_MathyAD,
             LauAD=Activate_LauAD
    )

library(UpSetR)
listOrder=c("MathysCtrl","LauCtrl","NagyCtrl","YangCtrl","MathysAD","LauAD")
pdf("Heatmap/IntersectForActivateList.pdf",height=4,width=7)
upset(fromList(ActivateList),sets =listOrder,keep.order = TRUE,nsets =6,nintersects=NA,
    queries = list(
    list(
      query = intersects,
      params = list(listOrder),
      active = T,
      color = "red"
     )
    ))
dev.off()

InhibiteList=list(MathysCtrl=Inhibite_MathysCtrl,
             LauCtrl=Inhibite_LauCtrl,
             NagyCtrl=Inhibite_NagyCtrl,
             YangCtrl=Inhibite_YangCtrl,
             MathysAD=Inhibite_MathyAD,
             LauAD=Inhibite_LauAD
    )
pdf("Heatmap/IntersectForInhibiteList.pdf",height=4,width=7)
upset(fromList(InhibiteList),nsets =6,nintersects=NA,sets =listOrder,keep.order = TRUE,
     queries = list(
     list(
      query = intersects,
      params = list(listOrder),
      active = T,
      color = "blue"
     )
    ))
dev.off()


ActivateOverlap=Reduce(intersect,list(ActivateList$MathysCtrl,ActivateList$LauCtrl,ActivateList$NagyCtrl,ActivateList$YangCtrl,ActivateList$MathysAD,ActivateList$LauAD))
toString(ActivateOverlap)
"CEBPB(+), SIN3A(+), FOXO3(+), KLF13(+), KLF6(+), POLR2A(+), NR2C2(+), MAZ(+), BCLAF1(+), RBBP5(+), BHLHE40(+), HDAC2(+), TAF7(+), YY1(+), ATF4(+), KLF9(+), EGR3(+), JUN(+)"

InhibiteOverlap=Reduce(intersect,list(InhibiteList$MathysCtrl,InhibiteList$LauCtrl,InhibiteList$NagyCtrl,InhibiteList$YangCtrl,InhibiteList$MathysAD,InhibiteList$LauAD))
"NR3C1(+), CUX1(+), FOXP2(+), ZEB1(+)"

