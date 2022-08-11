library("Seurat")
library("infercnv")
library("dplyr")
library("ggplot2")
library("RColorBrewer")
library("ggrepel")

####################################################################################
############       copy number variation for single cells       ####################
####################################################################################

# the output figures were showen in Figure 1F, 1G, 1H, 1I, 2A, S1D, S1E, S1F

MathyAllCtrlDifCell=readRDS("/Projects/deng/Aging/Ex/MathyEx/MathyCtrlDifferentiatedCells.rds") #The file will be changed to other dataset to obtain their copy number variation 
MathyAllCtrlEx=subset(MathyAllCtrlDifCell,seurat_clusters %in% c(2,3,4,5,6,8,10,13,14,16,17,18))

counts_matrix = as.matrix(MathyAllCtrlEx@assays$RNA@counts)
pData=MathyAllCtrlEx@meta.data$seurat_clusters
names(pData)=rownames(MathyAllCtrlEx@meta.data)
pData=data.frame(pData)
annotation=read.table("/Projects/deng/refGene/gencode.v38.annotation.Symbol.location.txt",header=TRUE,row.names=1)
anno=data.frame("Chr"=annotation$Chr,"Start"=annotation$Start,"End"=annotation$End)
rownames(anno)=annotation$Symbol
genes=intersect(rownames(counts_matrix),rownames(anno))
anno=anno[genes,]
anno$Chr=factor(anno$Chr,levels=paste0("chr",c(1:22,"X","Y"),sep=""))
anno=anno[order(anno$Chr,anno$Star,anno$End),]
genes=intersect(rownames(anno),rownames(counts_matrix))
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix[genes,],
                                      annotations_file=pData,
                                      delim="\t",
                                      gene_order_file=anno[genes,],
                                      ref_group_names=NULL)
  
infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics   0.1= 8037 genes and 5645 cells
                               denoise=TRUE,
                               cluster_by_groups = TRUE,
                               HMM=FALSE,
                               output_format = "png",
                               out_dir="/Projects/deng/Aging/Ex/MathyEx/InferCNV4Ex"
                               )



####################################################################################
############     Ribon plot for each cluster across  choromosome      ##############
####################################################################################

# the output figures were showen in Figure 1F, 2A, 2B, 2C, S1E, S2G

# copy number estimated by inferCNV
ExData=read.table("/Projects/deng/Aging/Ex/MathyEx/InferCNV4Ex/infercnv.observations.txt",check.names=F,row.names=1,header=T)

# Defination the sub cluster of the cell type to generate their global CNV plot, here use cluster 13 from Ex as example
MathyAllCtrlDifCell=readRDS("/Projects/deng/Aging/Ex/MathyEx/MathyCtrlDifferentiatedCells.rds")
MathyAllCtrlEx=subset(MathyAllCtrlDifCell,seurat_clusters %in% c(2,3,4,5,6,8,10,11,13,14,16,17,18)) #
C13=subset(MathyAllCtrlEx,seurat_clusters==13)

data=ExData
AllCellMean=apply(data,1,mean)
C13Cells=intersect(colnames(C13),colnames(data))
C13CellsExpr=data[,C13Cells]
C13CellMean=apply(C13CellsExpr,1,mean)
head(C13CellMean)
ExprFrame=cbind("AllCellMean"=AllCellMean,"C13CellMean"=C13CellMean,"CNV"=(C13CellMean-1)) # we used the average CNV scores across all cells minus 1 as the CNV score for each gene
rownames(ExprFrame)=rownames(data)
annotation=read.table("/Projects/deng/refGene/gencode.v38.annotation.removeDupliSymbol.geneType.txt",header=TRUE,row.names=1)
anno=data.frame("Chr"=annotation$Chr,"Start"=annotation$Start,"End"=annotation$End)
rownames(anno)=annotation$Symbol
anno=anno[anno$Chr%in%paste0("chr",c(1:22),sep=""),]
genes=intersect(rownames(data),rownames(anno))
anno=anno[genes,]
anno$Chr=factor(anno$Chr,levels=paste0("chr",c(1:22),sep=""))
anno=anno[order(anno$Chr,anno$Star,anno$End),]
geneOrder=intersect(rownames(anno),rownames(data))
ExprFrame4Graph=data.frame(ExprFrame[geneOrder,])
ExprFrame4Graph$index=c(1:dim(ExprFrame4Graph)[1])
all(rownames(anno)==rownames(ExprFrame4Graph))
ExprFrame4Graph=cbind(anno,ExprFrame4Graph)
ChrIndex <- ExprFrame4Graph %>% group_by(Chr) %>% filter (! duplicated(Chr))
t=ggplot(ExprFrame4Graph, aes(x=index,y = CNV)) + 
  geom_ribbon(aes(ymin=pmin(CNV,0), ymax=0), fill="blue", col="blue", alpha=0.5) +
  geom_ribbon(aes(ymin=0, ymax=pmax(CNV,0)), fill="red", col="red",alpha=0.5) +
  geom_vline(
    aes(xintercept = as.numeric(index)), 
    data = ChrIndex,
    colour = "grey50", alpha = 0.5
  )+
  geom_line(aes(y=0))+
  theme_bw()+ylim(-0.1,0.1)+
  theme(axis.title.x = element_blank(),axis.text.x = element_blank())+
  theme(panel.spacing.x = unit(0,"line"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

tiff("ExCluster13CNV.tiff",width=1000,height=200)
print(t)
dev.off()



####################################################################################
############        Arm-level aneuploidy estimated     ##############
####################################################################################


ExData=read.table("/Projects/deng/Aging/Ex/AllControlEx/InferCNV/infercnv.observations.txt",check.names=F,row.names=1,header=T)
ExCtrl=readRDS("/Projects/deng/Aging/Ex/AllControlEx/ExCtrl.integrated.rds")
SubType=subset(ExCtrl,seurat_clusters==0)

data=ExData

#obtain the reference value
AllCellMean=apply(data,1,mean) 
summary(AllCellMean)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.9972  1.0015  1.0033  1.0035  1.0057  1.0088
quantile(AllCellMean,0.95) #C0 1.007714
quantile(AllCellMean,0.05) #C0 0.9994199

SubType=subset(ExCtrl,seurat_clusters==3)
SubTypeCells=intersect(colnames(SubType),colnames(data))
SubTypeCellsExpr=data[,SubTypeCells]
SubTypeCellMean=apply(SubTypeCellsExpr,1,mean)
ExprFrame=cbind("SubTypeCellMean"=SubTypeCellMean,"CNV"=(SubTypeCellMean-1))
summary(ExprFrame[,2])


GainGene=Ratio[Ratio$GainRatio>0.25,]#SubTypeCellsExpr=SubTypeCellsExpr-1

cellNumber=dim(SubTypeCellsExpr)[2]
geneNumber=dim(SubTypeCellsExpr)[1]
Ratio=matrix(data = NA, nrow = geneNumber, ncol = 2, dimnames = list(rownames(ExprFrame), c("Gain","Loss")))
for(i in 1:dim(ExprFrame)[1]){
  Ratio[i,1]=sum(SubTypeCellsExpr[i,]>1.007714)
  Ratio[i,2]=sum(SubTypeCellsExpr[i,]<0.9994199)
  if(!i%%100){
    print(i)
  }
}
tmp=Ratio
Ratio=data.frame(Ratio)
Ratio$GainRatio=Ratio$Gain/cellNumber
Ratio$LossRatio=Ratio$Loss/cellNumber

LossGene=Ratio[Ratio$LossRatio>0.25,]
write.table(GainGene,file="CNV8Arm/GainGeneC3.txt",sep="\t",quote=F)
write.table(LossGene,file="CNV8Arm/LossGeneC3.txt",sep="\t",quote=F)

head(SubTypeCellMean)
ExprFrame=cbind("GainRatio"=Ratio$GainRatio,"LossRatio"=Ratio$LossRatio)
rownames(ExprFrame)=rownames(Ratio)

AllGeneChrLocation=read.table("/Projects/deng/Public/HGNC/GeneChrLocation.txt",header=T,row.names=1,sep="\t")
AllGeneChrLocation=AllGeneChrLocation[,c("Symbol","ChLocation")]
AllGeneChrLocation$BandTotal=sub('(\\d+[p|q])(\\d+.*)',"\\1", AllGeneChrLocation$ChLocation)
BandTotal=data.frame(table(AllGeneChrLocation$BandTotal))
colnames(BandTotal)=c("Band","TotalBand")
genes=intersect(genes,AllGeneChrLocation$Symbol) 
AllGeneChrLocation=AllGeneChrLocation[AllGeneChrLocation$Symbol%in%genes,]
rownames(AllGeneChrLocation)=AllGeneChrLocation$Symbol
gene=intersect(rownames(ExprFrame),rownames(AllGeneChrLocation))
ExprFrame4Graph=ExprFrame[gene,]
AllGeneChrLocation=AllGeneChrLocation[gene,]
all(rownames(ExprFrame4Graph)==rownames(AllGeneChrLocation))
ExprFrame4Graph=data.frame(ExprFrame4Graph)
ExprFrame4Graph$Band=AllGeneChrLocation$BandTotal

ExprFrame4Graph=ExprFrame4Graph[!ExprFrame4Graph$Band=="q33.3",]
library(reshape2)
ExprFrame4GraphTmp=melt(ExprFrame4Graph,id=c(3))
colnames(ExprFrame4GraphTmp)=c("Band","Group","Ratio")


t=ggplot(ExprFrame4GraphTmp, aes(Group,Ratio,color=Group)) +
  geom_point(aes(fill=Group),size=0.3, show.legend=TRUE, alpha=0.5, position = position_jitterdodge(jitter.width =0.75,dodge.width = 0.75))+
  scale_shape_manual(values=c(21, 22, 24))+scale_size(range = c(2, 4))+
  scale_color_manual(values=c("red","Blue"))+
  scale_fill_manual(values=c("red","Blue"))+
  ylim(0,1)+
  theme_bw()+facet_wrap(factor(Band,levels=unique(ExprFrame4GraphTmp$Band))~.,ncol=6)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
pdf("CNVsplit8Arm4Cluster4.pdf",width=10)
print(t)
dev.off()



####################################################################################
############        Integrated the AD loci into infercnv barplot      ##############
####################################################################################
ExData=read.table("/Projects/deng/Aging/Ex/AllEx/InferCNV/infercnv.observations.txt",check.names=F,row.names=1,header=T)
ExAll=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
table(ExAll$Statues)
SubType=subset(ExAll,seurat_clusters==3)
SubType
#SubType=subset(ExAll,seurat_clusters==6&Source%in%"Nagy") #for each dataset
data=ExData
#data=OliData
SubTypeCells=intersect(colnames(SubType),colnames(data))
SubTypeCellsExpr=data[,SubTypeCells]
SubTypeCellMean=apply(SubTypeCellsExpr,1,mean)
head(SubTypeCellMean)
ExprFrame=cbind("SubTypeCellMean"=SubTypeCellMean,"CNV"=(SubTypeCellMean-1))
quantile(ExprFrame[,2],0.95) #C0 0.009026389
rownames(ExprFrame)=rownames(data)
annotation=read.table("/Projects/deng/refGene/gencode.v38.annotation.removeDupliSymbol.geneType.txt",header=TRUE,row.names=1)
anno=data.frame("Chr"=annotation$Chr,"Start"=annotation$Start,"End"=annotation$End)
rownames(anno)=annotation$Symbol
anno=anno[anno$Chr%in%paste0("chr",c(1:22),sep=""),]
genes=intersect(rownames(data),rownames(anno))
AllGeneChrLocation=read.table("/Projects/deng/Public/HGNC/GeneChrLocation.txt",header=T,row.names=1,sep="\t")
AllGeneChrLocation=AllGeneChrLocation[,c("Symbol","ChLocation")]
#AllGeneChrLocation$BandTotal=sub('(\\d+[p|q]\\d+)(.*)',"\\1", AllGeneChrLocation$ChLocation)
AllGeneChrLocation$BandTotal=AllGeneChrLocation$ChLocation
BandTotal=data.frame(table(AllGeneChrLocation$BandTotal))
colnames(BandTotal)=c("Band","TotalBand")

genes=intersect(genes,GeneChrLocation$Symbol) 
length(genes)
anno=anno[genes,]
anno$Chr=factor(anno$Chr,levels=paste0("chr",c(1:22),sep=""))
anno=anno[order(anno$Chr,anno$Star,anno$End),]
geneOrder=intersect(rownames(anno),rownames(data))
ExprFrame4Graph=data.frame(ExprFrame[geneOrder,])
ExprFrame4Graph$index=c(1:dim(ExprFrame4Graph)[1]) #build the index for each gene
all(rownames(anno)==rownames(ExprFrame4Graph))
ExprFrame4Graph=cbind(anno,ExprFrame4Graph)
GeneChrLocation=AllGeneChrLocation[AllGeneChrLocation$Symbol%in%genes,]
rownames(GeneChrLocation)=GeneChrLocation$Symbol
GeneChrLocation=GeneChrLocation[genes,]
all(rownames(GeneChrLocation)==rownames(ExprFrame4Graph)) #TRUE
#ExprFrame4Graph$Band=sub('(\\d+[p|q]\\d+)(.*)',"\\1", GeneChrLocation$ChLocation)
ExprFrame4Graph$Band=GeneChrLocation$ChLocation
ExprFrame4Graph$CNVSig=ifelse(abs(ExprFrame4Graph$CNV)>0.01,"Sig","NonSig")
BandTargetInSingleCell=data.frame(table(ExprFrame4Graph$Band))
colnames(BandTargetInSingleCell)=c("Band","TotalBandInSC")
BandSum=ExprFrame4Graph %>% group_by(Band) %>% summarise(CNVSum = sum(CNV))

Bands=intersect(BandTargetInSingleCell$Band,BandTotal$Band)
BandTargetInSingleCell=BandTargetInSingleCell[BandTargetInSingleCell$Band%in%Bands,]
BandTotal=BandTotal[BandTotal$Band%in%Bands,]
BandInfo=merge(BandTotal,BandTargetInSingleCell,by="Band")

SigBand=ExprFrame4Graph[ExprFrame4Graph$CNVSig=="Sig",]
BandSig=data.frame(table(SigBand$Band))
colnames(BandSig)=c("Band","SigBand")
BandInfoTotal=merge(BandInfo,BandSig,by="Band",all.x = TRUE)
BandInfoTotal[is.na(BandInfoTotal)]=0
BandInfoTotal$Ratio=BandInfoTotal$SigBand/BandInfoTotal$TotalBand

ExprFrame4Graph$Symbol=rownames(ExprFrame4Graph)
ExprFrame4GraphBandMerge=merge(ExprFrame4Graph,BandInfoTotal,by="Band")
ExprFrame4GraphBandMerge=ExprFrame4GraphBandMerge[order(ExprFrame4GraphBandMerge$index),]
rownames(ExprFrame4GraphBandMerge)=ExprFrame4GraphBandMerge$Symbol

ChrBandIndex <- ExprFrame4Graph %>% group_by(Band) %>% filter (! duplicated(Band)) #extract the start location for each chr
AllResult=merge(ChrBandIndex,BandInfoTotal,by="Band")
AllResult=merge(AllResult,BandSum,by="Band")
AllResult=AllResult[order(AllResult$index),]
BandMidIndex=array()
for(i in 2:(dim(AllResult)[1]+1)){
    BandMidIndex[i-1]=AllResult[i,"index"]-AllResult[(i-1),"index"]
}
AllResult$barWidth=BandMidIndex


#Add the GWAS loci information on the graph
GWAS=read.table("/Projects/deng/Aging/Ex/AllControlEx/GWASLoci.txt")
ADLoci=intersect(GWAS[,1],AllResult$Band)
AllResult$ADLoci=ifelse(AllResult$Band%in%ADLoci,"Yes","No")
AllResult$colorLabel=ifelse(AllResult$Band%in%ADLoci,"ADLoci",ifelse(AllResult$Chr%in%c(paste0("chr",seq(1,21,2))),"Single","Double"))

t=ggplot(AllResult, aes(index, abs(CNVSum),fill=colorLabel,color=colorLabel,width=barWidth))+geom_bar(stat="identity")+
scale_fill_manual(values = c("DarkOrange","DarkGray","LightGrey"))+
scale_color_manual(values = c("DarkOrange","DarkGray","LightGrey"))+
theme_bw()+ylim(0,6)+
geom_label_repel(  data = AllResult[(AllResult$Band%in%ADLoci&abs(AllResult$CNVSum)>0.5),],
                   aes(x = index, y = abs(CNVSum), label = Band),
                   size = 5,
                   colour="black",
                   force_pull   = 0, 
                   nudge_x = 0.5,
                   box.padding = 0.5,
                   nudge_y = 0.5,
                   min.segment.length = 0, # draw all lines no matter how short
                   segment.size = 0.2,
                   segment.colour="black",
                   segment.curvature = -0.1,
                   segment.ncp = 3,
                   segment.angle = 90,
                   label.size=NA, #no border/box
                   fill = "WhiteSmoke")+  #no background
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=rel(1.0),colour = "black"))
tiff("ChrBandADLociInC3.tiff",width=1000,height=200)
print(t)
dev.off()




