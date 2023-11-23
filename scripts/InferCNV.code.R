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
annotation=read.table("/Projects/deng/refGene/gencode.v38.annotation.removeDupliSymbol.geneType.txt",header=TRUE,row.names=1)
anno=data.frame("Chr"=annotation$Chr,"Start"=annotation$Start,"End"=annotation$End)
rownames(anno)=annotation$Symbol
genes=intersect(rownames(counts_matrix),rownames(anno))
length(genes)#17046
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
                               HMM=TRUE,
                               output_format = "png",
                               out_dir="/data2/deng/Aging/Ex/MathyEx/InferCNVEx"
                               )


##########################################################################################################
############     pvalue between cell cycle gene re-expressed cluster and normal ones      ##############
##########################################################################################################

MathyEx.infercnv.observations=read.table("/Projects/deng/Aging/Ex/MathyEx/InferCNV4Ex/infercnv.observations.txt",check.names=F,row.names=1,header=T)
MathyAllCtrlDifCell=readRDS("/Projects/deng/Aging/Ex/MathyEx/MathyCtrlDifferentiatedCells.rds") #The file will be changed to other dataset to obtain their copy number variation 
C2=subset(MathyAllCtrlDifCell,seurat_clusters==2)
C4=subset(MathyAllCtrlDifCell,seurat_clusters==4)
C5=subset(MathyAllCtrlDifCell,seurat_clusters==5)
C13=subset(MathyAllCtrlDifCell,seurat_clusters==13)

data=MathyEx.infercnv.observations

C2Cells=intersect(colnames(C2),colnames(data))
C2CellsExpr=data[,C2Cells]
C2CellMean=apply(C2CellsExpr,1,mean)
C4Cells=intersect(colnames(C4),colnames(data))
C4CellsExpr=data[,C4Cells]
C4CellMean=apply(C4CellsExpr,1,mean)
C5Cells=intersect(colnames(C5),colnames(data))
C5CellsExpr=data[,C5Cells]
C5CellMean=apply(C5CellsExpr,1,mean)
C13Cells=intersect(colnames(C13),colnames(data))
C13CellsExpr=data[,C13Cells]
C13CellMean=apply(C13CellsExpr,1,mean)

t.test(abs(C13CellMean-1),abs(C2CellMean-1),alternative = c("greater")) #p-value < 2.2e-16
t.test(abs(C13CellMean-1),abs(C5CellMean-1),alternative = c("greater")) #p-value < 2.2e-16
t.test(abs(C5CellMean-1),abs(C2CellMean-1),alternative = c("greater")) #p-value < 2.2e-16
t.test(abs(C4CellMean-1),abs(C2CellMean-1),alternative = c("greater")) #p-value = 1
####################################################################################
############     Ribon plot for each cluster across  choromosome      ##############
####################################################################################

# the output figures were showen in Figure 2E
ExAll.integrated=readRDS("ExAll.integratedTSNE.rds")
ExAll.integrated=subset(ExAll.integrated,idents=c(0:15))
Group="ExAllDownsampling"
setwd("/data2/deng/Aging/Ex/AllEx/InferCNV/ExAllDownsampling")
infercnv.observations=read.table("infercnv.observations.txt",check.names=F,row.names=1,header=T)

for(targetCluster in c(0,3,5)){
subType=subset(ExAll.integrated,seurat_clusters==targetCluster)
data=infercnv.observations
AllCellMean=apply(data,1,mean)
subTypeCells=intersect(colnames(subType),colnames(data))
subTypeCellsExpr=data[,subTypeCells]
subTypeCellMean=apply(subTypeCellsExpr,1,mean)
head(subTypeCellMean)
ExprFrame=cbind("AllCellMean"=AllCellMean,"subTypeCellMean"=subTypeCellMean,"CNV"=(subTypeCellMean-1)) # we used the average CNV scores across all cells minus 1 as the CNV score for each gene
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

  pdf(paste0("/data2/deng/Aging/Ex/AllEx/InferCNV/RibonPlot/RibonPlot_ExCtrl_Cluster",targetCluster,"_",Group,".pdf"),width=8,height=2)
  print(t)
  dev.off()
  
  #output meta data
  write.table(ExprFrame4Graph,file=paste0("/data2/deng/Aging/Ex/AllEx/InferCNV/RibonPlot/RibonPlot_ExCtrl_Cluster",targetCluster,"_",Group,".txt",sep=""),quote=F,sep="\t")
 }




ExAll.integrated$Group=paste0(ExAll.integrated$Source,"_",ExAll.integrated$Statues,sep="")

for(Group in unique(ExAll.integrated$Group)){
setwd(paste0("/data2/deng/Aging/Ex/AllEx/InferCNV/",Group,sep=""))
infercnv.observations=read.table("infercnv.observations.txt",check.names=F,row.names=1,header=T)
for(targetCluster in c(0,3,5)){
subType=subset(ExAll.integrated,seurat_clusters==targetCluster)
data=infercnv.observations
AllCellMean=apply(data,1,mean)
subTypeCells=intersect(colnames(subType),colnames(data))
subTypeCellsExpr=data[,subTypeCells]
subTypeCellMean=apply(subTypeCellsExpr,1,mean)
head(subTypeCellMean)
ExprFrame=cbind("AllCellMean"=AllCellMean,"subTypeCellMean"=subTypeCellMean,"CNV"=(subTypeCellMean-1)) # we used the average CNV scores across all cells minus 1 as the CNV score for each gene
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

  pdf(paste0("/data2/deng/Aging/Ex/AllEx/InferCNV/RibonPlot/RibonPlot_ExCtrl_Cluster",targetCluster,"_",Group,".pdf"),width=8,height=2)
  print(t)
  dev.off()

  #output meta data
  write.table(ExprFrame4Graph,file=paste0("/data2/deng/Aging/Ex/AllEx/InferCNV/RibonPlot/RibonPlot_ExCtrl_Cluster",targetCluster,"_",Group,".txt",sep=""),quote=F,sep="\t")
 }
}



####################################################################################
############        Arm-level aneuploidy estimated     ##############
####################################################################################
#Figure 2F

setwd("/data2/deng/Aging/Ex/AllEx/InferCNV/ExAllDownsampling")
ExAll=readRDS("/data2/deng/Aging/Ex/AllEx/InferCNV/ExAll.Downsampling.rds")
#infercnv.20_HMM_predHMMi6.hmm_mode-samples.Pnorm_0.5.repr_intensities.observations.txt
ExAllData=read.table("/data2/deng/Aging/Ex/AllEx/InferCNV/ExAllDownsampling/infercnv.20_HMM_predHMMi6.hmm_mode-samples.Pnorm_0.5.repr_intensities.observations.txt",check.names=F,row.names=1,header=T)
data=ExAllData

C0=subset(ExAll,seurat_clusters==0) #49835 genes across 17195 cells within 2 assays
C3=subset(ExAll,seurat_clusters==3) #49835 genes across 7147 cells within 2 assays
C5=subset(ExAll,seurat_clusters==5) #49835 genes across 6612 genes within 2 assays

C0Cells=intersect(colnames(C0),colnames(data))
C0CellsIntensities=data[,C0Cells]
C0CellsIntensMedian=apply(C0CellsIntensities,1,mean) #median score for the single genes across C0

C3Cells=intersect(colnames(C3),colnames(data))
C3CellsIntensities=data[,C3Cells]
C3CellsIntensMedian=apply(C3CellsIntensities,1,mean) #median score for the single genes across C3

C5Cells=intersect(colnames(C5),colnames(data))
C5CellsIntensities=data[,C5Cells]
C5CellsIntensMedian=apply(C5CellsIntensities,1,mean) #median score for the single genes across C5

all(names(C0CellsIntensMedian)==names(C5CellsIntensMedian))
all(names(C0CellsIntensMedian)==names(C3CellsIntensMedian))
CNVIntensFrame=data.frame(cbind(C0CNV=C0CellsIntensMedian,C3CNV=C3CellsIntensMedian,C5CNV=C5CellsIntensMedian))

#obtain the location for each genes from HGNC database
AllGeneChrLocation=read.table("/Projects/deng/Public/HGNC/GeneChrLocation.txt",header=T,row.names=1,sep="\t")
AllGeneChrLocation=AllGeneChrLocation[,c("Symbol","ChLocation")]
AllGeneChrLocation$BandTotal=sub('(\\d+[p|q])(\\d+.*)',"\\1", AllGeneChrLocation$ChLocation)
rownames(AllGeneChrLocation)=AllGeneChrLocation$Symbol

gene=intersect(rownames(CNVIntensFrame),rownames(AllGeneChrLocation))
length(gene)
#7981
CNVIntensFrame4Graph=CNVIntensFrame[gene,]
AllGeneChrLocation=AllGeneChrLocation[gene,]
all(rownames(CNVIntensFrame4Graph)==rownames(AllGeneChrLocation))
CNVIntensFrame4Graph$Band=AllGeneChrLocation$BandTotal
CNVIntensFrame4Graph=CNVIntensFrame4Graph[!CNVIntensFrame4Graph$Band=="q33.3",]
table(CNVIntensFrame4Graph$Band)
#10p 10q 11p 11q 12p 12q 13q 14q 15q 16p 16q 17p 17q 18p 18q 19p 19q  1p  1q 20p
#80 264 152 277 112 350 176 289 285 176 150 111 309  37 102 222 219 437 339  74
#20q 21q 22q  2p  2q  3p  3q  4p  4q  5p  5q  6p  6q  7p  7q  8p  8q  9p  9q
#144  73 187 250 374 252 271 104 241  74 361 181 236 144 267 110 206  94 250

#output meta data
write.table(CNVIntensFrame4Graph,file="CNVIntensFrame4Graph_Genedata.txt",sep="\t",quote=F)

GainInC5=CNVIntensFrame4Graph[CNVIntensFrame4Graph$C5CNV>1,]
dim(GainInC5)#1798    4
write.table(GainInC5,file="GainInC5.txt",sep="\t",quote=F)
LossInC5=CNVIntensFrame4Graph[CNVIntensFrame4Graph$C5CNV<1,]
write.table(LossInC5,file="LossInC5.txt",sep="\t",quote=F)
dim(LossInC5) #911   4


CNVIntensFrame4GraphTmp=reshape2::melt(CNVIntensFrame4Graph,id=c(4))
colnames(CNVIntensFrame4GraphTmp)=c("Band","Group","MedianIntensities")
CNVIntensFrame4GraphTmp$Pattern=ifelse(CNVIntensFrame4GraphTmp$MedianIntensities>1,"Gain",ifelse(CNVIntensFrame4GraphTmp$MedianIntensities<1,"Loss","Normal"))
CNVIntensFrame4GraphTmp$Pattern=factor(CNVIntensFrame4GraphTmp$Pattern,levels=c("Gain","Normal","Loss"))
CNVIntensFrame4GraphTmp$Band=factor(CNVIntensFrame4GraphTmp$Band,levels=c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q","11p","11q","12p","12q","13q","14q","15q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q"))
CNVIntensFrame4GraphTmp$Group=factor(CNVIntensFrame4GraphTmp$Group,levels=c("C5CNV","C3CNV","C0CNV"))

write.table(CNVIntensFrame4GraphTmp,file="CNVIntensFrame4Graph_data.txt",row.names=F,sep="\t",quote=F)


t=ggplot(CNVIntensFrame4GraphTmp, aes(Group,MedianIntensities,color=Pattern)) +
  geom_point(aes(fill=Pattern),size=0.3, show.legend=TRUE, alpha=0.5, position = position_jitterdodge(jitter.width =0.75,dodge.width = 0.75))+
  scale_color_manual(values=c("Red","Lightgrey","Blue"))+
  scale_fill_manual(values=c("Red","Lightgrey","Blue"))+
  #ylim(0,1)+
  theme_bw()+facet_wrap(factor(Band,levels=unique(CNVIntensFrame4GraphTmp$Band))~.,ncol=6)+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),axis.title.x = element_blank(),axis.title.y = element_blank())
pdf("CNVHMMIntensByBand.pdf",width=8,height=9)
print(t)
dev.off()


t=ggplot(CNVIntensFrame4GraphTmp, aes(Group,MedianIntensities,color=Group)) +
  geom_point(aes(fill=Group),size=0.3, show.legend=TRUE, alpha=0.5, position = position_jitterdodge(jitter.width =0.85,dodge.width = 0.85))+
  scale_color_manual(values=c("#8DD3C7","#FFA500","#CD2626"))+
  scale_fill_manual(values=c("#8DD3C7","#FFA500","#CD2626"))+
  #ylim(0,1)+
  theme_bw()+facet_wrap(factor(Band,levels=unique(CNVIntensFrame4GraphTmp$Band))~.,ncol=8)+
  theme(
    axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.spacing = unit(0,'lines'))
pdf("CNVHMMIntensByBand.Cluster.pdf",width=9,height=6)
print(t)
dev.off()


t=ggplot(CNVIntensFrame4GraphTmp, aes(Band,MedianIntensities,color=Pattern)) +
  geom_point(aes(fill=Pattern),size=0.3, show.legend=TRUE, alpha=0.5, position = position_jitterdodge(jitter.width =0.85,dodge.width = 0.85))+
  scale_color_manual(values=c("Red","Lightgrey","Blue"))+
  scale_fill_manual(values=c("Red","Lightgrey","Blue"))+
  #ylim(0,1)+
  theme_bw()+facet_grid(Group~.,)+
  theme(
    axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.spacing = unit(0.1,'lines'))
pdf("CNVHMMIntensByBand.Pattern.pdf",width=9,height=3)
print(t)
dev.off()





####################################################################################
############        Integrated the AD loci into infercnv barplot      ##############
####################################################################################
#ExData=read.table("/Projects/deng/Aging/Ex/AllEx/InferCNV/infercnv.observations.txt",check.names=F,row.names=1,header=T)
ExAll.integrated=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
DefaultAssay(ExAll.integrated)="RNA"
ExAll.integrated$Group=paste0(ExAll.integrated$Source,"_",ExAll.integrated$Statues,sep="")

for(g in unique(ExAll.integrated$Group)){
  setwd(paste0("/data2/deng/Aging/Ex/AllEx/InferCNV/",g,sep=""))
  ExData=read.table("infercnv.observations.txt",check.names=F,row.names=1,header=T)

  for(targetCluster in c(0,3,5)){
     SubType=subset(ExAll.integrated,Group==g&seurat_clusters==targetCluster)
     #SubType=subset(ExAll.integrated,seurat_clusters==targetCluster) #only for the integrated dataset
     data=ExData

     SubTypeCells=intersect(colnames(SubType),colnames(data))
     SubTypeCellsExpr=data[,SubTypeCells]
     SubTypeCellMean=apply(SubTypeCellsExpr,1,mean)
     #obtain the SNV score, CNV score with 1 reprent no gain or loss
     ExprFrame=cbind("SubTypeCellMean"=SubTypeCellMean,"CNV"=(SubTypeCellMean-1))
     quantile(ExprFrame[,2],0.95) #C0 0.009026389
     rownames(ExprFrame)=rownames(data)

     #this gene anotation file was used to order the genes based on their chromosome location from chr1 to chr22
     annotation=read.table("/Projects/deng/refGene/gencode.v38.annotation.removeDupliSymbol.geneType.txt",header=TRUE,row.names=1)
     anno=data.frame("Chr"=annotation$Chr,"Start"=annotation$Start,"End"=annotation$End)
     rownames(anno)=annotation$Symbol
     anno=anno[anno$Chr%in%paste0("chr",c(1:22),sep=""),] #only consider the genes from chr1 to chr22, remove genes from chrM
     genes=intersect(rownames(data),rownames(anno))

     #to get the band information for each genes, this information were downloaded from HGNC database
     GeneChrLocation=read.table("/Projects/deng/Public/HGNC/GeneChrLocation.txt",header=T,row.names=1,sep="\t")
     GeneChrLocation=GeneChrLocation[,c("Symbol","ChLocation")]
     #GeneChrLocation$BandTotal=sub('(\\d+[p|q]\\d+)(.*)',"\\1", GeneChrLocation$ChLocation)
     GeneChrLocation$BandTotal=GeneChrLocation$ChLocation
     BandTotal=data.frame(table(GeneChrLocation$BandTotal)) #get the total number genes for each Band
     colnames(BandTotal)=c("Band","TotalBand")

     #order CNV matrix based on their gene order on the chromosome
     genes=intersect(genes,GeneChrLocation$Symbol) 
     length(genes)
     anno=anno[genes,]
     anno$Chr=factor(anno$Chr,levels=paste0("chr",c(1:22),sep=""))
     anno=anno[order(anno$Chr,anno$Star,anno$End),]
     geneOrder=intersect(rownames(anno),rownames(data))
     ExprFrame4Graph=data.frame(ExprFrame[geneOrder,])
     
     #Add the gene chromosome information to CNV score matrix
     ExprFrame4Graph$index=c(1:dim(ExprFrame4Graph)[1]) #build the index for each gene
     all(rownames(anno)==rownames(ExprFrame4Graph))
     ExprFrame4Graph=cbind(anno,ExprFrame4Graph)

     #Add the gene band information to CNV score matrix
     rownames(GeneChrLocation)=GeneChrLocation$Symbol
     GeneChrLocation=GeneChrLocation[genes,]
     all(rownames(GeneChrLocation)==rownames(ExprFrame4Graph)) #TRUE
     ExprFrame4Graph$Band=GeneChrLocation$ChLocation

     #got the CNV score for each band (as the sum score within each Band)
     BandSum=ExprFrame4Graph %>% group_by(Band) %>% summarise(CNVSum = sum(CNV))
     
     #ExprFrame4Graph$Symbol=rownames(ExprFrame4Graph)
     ExprFrame4GraphAddBandSum=merge(ExprFrame4Graph,BandSum,by="Band") #this si the total CNV score matrix (keep the gene information)
     
     ExprFrame4GraphAddBandSum.uni=ExprFrame4GraphAddBandSum[!duplicated(ExprFrame4GraphAddBandSum$Band),] #randomly seleced one band (not this step removed the gene information)
     ExprFrame4GraphAddBandSum.uni=ExprFrame4GraphAddBandSum.uni[order(ExprFrame4GraphAddBandSum.uni$index),c("Band","index","Chr","CNVSum")] #remove the gene dependent columns
     #obtain the width of the band based on their expressed gene number
     BandMidIndex=array()
     for(i in 2:(dim(ExprFrame4GraphAddBandSum.uni)[1]+1)){
         BandMidIndex[i-1]=ExprFrame4GraphAddBandSum.uni[i,"index"]-ExprFrame4GraphAddBandSum.uni[(i-1),"index"]
     }
     ExprFrame4GraphAddBandSum.uni$barWidth=BandMidIndex
     
     #Add the GWAS loci information on the graph, D:/Alzheimer/GWAS
     GWAS=read.table("/data2/deng/Alzheimer/GWAS/ADLociBandInfor.txt",header=T)
     ADLoci=intersect(GWAS$ChLocation,ExprFrame4GraphAddBandSum.uni$Band)
     ExprFrame4GraphAddBandSum.uni$ADLoci=ifelse(ExprFrame4GraphAddBandSum.uni$Band%in%ADLoci,"Yes","No")
     ExprFrame4GraphAddBandSum.uni$colorLabel=ifelse(ExprFrame4GraphAddBandSum.uni$Band%in%ADLoci,"ADLoci",ifelse(ExprFrame4GraphAddBandSum.uni$Chr%in%c(paste0("chr",seq(1,21,2))),"Single","Double"))

     #ExprFrame4GraphAddBandSum$ADGene=ifelse(ExprFrame4GraphAddBandSum$Symbol%in%GWAS$Symbol,"ADRisk","NormalGene")
     #ExprFrame4GraphAddBandSum=ExprFrame4GraphAddBandSum[order(ExprFrame4GraphAddBandSum$index),]
     #write.table(ExprFrame4GraphAddBandSum,file="/data2/deng/Aging/Ex/AllEx/InferCNV/ADRiskPlot/ExprFrame4GraphAddBandSum.txt",sep="\t",quote=F,row.names=F)
     
     t=ggplot(ExprFrame4GraphAddBandSum.uni, aes(index, CNVSum,fill=colorLabel,color=colorLabel,width=barWidth))+
       geom_bar(stat="identity")+
       scale_fill_manual(values = c("DarkOrange","DarkGray","LightGrey"))+
       scale_color_manual(values = c("DarkOrange","DarkGray","LightGrey"))+
       theme_bw()+ 
       ylim(-3,6)+
       geom_label_repel(  data = ExprFrame4GraphAddBandSum.uni[(ExprFrame4GraphAddBandSum.uni$Band%in%ADLoci&abs(ExprFrame4GraphAddBandSum.uni$CNVSum)>0.5),],
                   aes(x = index, y = CNVSum, label = Band),
                   size = 2,
                   colour="black",
                   max.overlaps=20,
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
                   fill = NA
                   )+  #no background
                   theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=rel(1.0),colour = "black"))
     pdf(paste0("/data2/deng/Aging/Ex/AllEx/InferCNV/ADRiskPlot/ADRiskPlotFor_Cluster",targetCluster,"_",g,".pdf"),width=10,height=1.5)
     print(t)
     dev.off()

     
     #output meta data
     write.table(ExprFrame4GraphAddBandSum.uni,file=paste0("/data2/deng/Aging/Ex/AllEx/InferCNV/ADRiskPlot/ADRiskPlotFor_Cluster",targetCluster,"_",g,".txt",sep=""),quote=F,sep="\t")
  }
}





ExData=read.table("/Projects/deng/Aging/Ex/AllEx/InferCNV/infercnv.observations.txt",check.names=F,row.names=1,header=T)
ExAll.integrated=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
DefaultAssay(ExAll.integrated)="RNA"
ExAll.integrated$Group=paste0(ExAll.integrated$Source,"_",ExAll.integrated$Statues,sep="")

g="ExAllDownsampling"
setwd(paste0("/data2/deng/Aging/Ex/AllEx/InferCNV/",g,sep=""))

for(targetCluster in c(0,3,5)){
     SubType=subset(ExAll.integrated,seurat_clusters==targetCluster)
     #SubType=subset(ExAll.integrated,seurat_clusters==targetCluster) #only for the integrated dataset
     data=ExData

     SubTypeCells=intersect(colnames(SubType),colnames(data))
     SubTypeCellsExpr=data[,SubTypeCells]
     SubTypeCellMean=apply(SubTypeCellsExpr,1,mean)
     #obtain the SNV score, CNV score with 1 reprent no gain or loss
     ExprFrame=cbind("SubTypeCellMean"=SubTypeCellMean,"CNV"=(SubTypeCellMean-1))
     quantile(ExprFrame[,2],0.95) #C0 0.009026389
     rownames(ExprFrame)=rownames(data)

     #this gene anotation file was used to order the genes based on their chromosome location from chr1 to chr22
     annotation=read.table("/Projects/deng/refGene/gencode.v38.annotation.removeDupliSymbol.geneType.txt",header=TRUE,row.names=1)
     anno=data.frame("Chr"=annotation$Chr,"Start"=annotation$Start,"End"=annotation$End)
     rownames(anno)=annotation$Symbol
     anno=anno[anno$Chr%in%paste0("chr",c(1:22),sep=""),] #only consider the genes from chr1 to chr22, remove genes from chrM
     genes=intersect(rownames(data),rownames(anno))

     #to get the band information for each genes, this information were downloaded from HGNC database
     GeneChrLocation=read.table("/Projects/deng/Public/HGNC/GeneChrLocation.txt",header=T,row.names=1,sep="\t")
     GeneChrLocation=GeneChrLocation[,c("Symbol","ChLocation")]
     #GeneChrLocation$BandTotal=sub('(\\d+[p|q]\\d+)(.*)',"\\1", GeneChrLocation$ChLocation)
     GeneChrLocation$BandTotal=GeneChrLocation$ChLocation
     BandTotal=data.frame(table(GeneChrLocation$BandTotal)) #get the total number genes for each Band
     colnames(BandTotal)=c("Band","TotalBand")

     #order CNV matrix based on their gene order on the chromosome
     genes=intersect(genes,GeneChrLocation$Symbol) 
     length(genes)
     anno=anno[genes,]
     anno$Chr=factor(anno$Chr,levels=paste0("chr",c(1:22),sep=""))
     anno=anno[order(anno$Chr,anno$Star,anno$End),]
     geneOrder=intersect(rownames(anno),rownames(data))
     ExprFrame4Graph=data.frame(ExprFrame[geneOrder,])
     
     #Add the gene chromosome information to CNV score matrix
     ExprFrame4Graph$index=c(1:dim(ExprFrame4Graph)[1]) #build the index for each gene
     all(rownames(anno)==rownames(ExprFrame4Graph))
     ExprFrame4Graph=cbind(anno,ExprFrame4Graph)

     #Add the gene band information to CNV score matrix
     rownames(GeneChrLocation)=GeneChrLocation$Symbol
     GeneChrLocation=GeneChrLocation[genes,]
     all(rownames(GeneChrLocation)==rownames(ExprFrame4Graph)) #TRUE
     ExprFrame4Graph$Band=GeneChrLocation$ChLocation

     #got the CNV score for each band (as the sum score within each Band)
     BandSum=ExprFrame4Graph %>% group_by(Band) %>% summarise(CNVSum = sum(CNV))
     
     #ExprFrame4Graph$Symbol=rownames(ExprFrame4Graph)
     ExprFrame4GraphAddBandSum=merge(ExprFrame4Graph,BandSum,by="Band") #this si the total CNV score matrix (keep the gene information)
     
     ExprFrame4GraphAddBandSum.uni=ExprFrame4GraphAddBandSum[!duplicated(ExprFrame4GraphAddBandSum$Band),] #randomly seleced one band (not this step removed the gene information)
     ExprFrame4GraphAddBandSum.uni=ExprFrame4GraphAddBandSum.uni[order(ExprFrame4GraphAddBandSum.uni$index),c("Band","index","Chr","CNVSum")] #remove the gene dependent columns
     #obtain the width of the band based on their expressed gene number
     BandMidIndex=array()
     for(i in 2:(dim(ExprFrame4GraphAddBandSum.uni)[1]+1)){
         BandMidIndex[i-1]=ExprFrame4GraphAddBandSum.uni[i,"index"]-ExprFrame4GraphAddBandSum.uni[(i-1),"index"]
     }
     ExprFrame4GraphAddBandSum.uni$barWidth=BandMidIndex
     
     #Add the GWAS loci information on the graph, D:/Alzheimer/GWAS
     #GWAS=read.table("/data2/deng/Alzheimer/GWAS/ADLociBandInfor.txt",header=T)
     #ADLoci=intersect(GWAS$ChLocation,ExprFrame4GraphAddBandSum.uni$Band)
     #ExprFrame4GraphAddBandSum.uni$ADLoci=ifelse(ExprFrame4GraphAddBandSum.uni$Band%in%ADLoci,"Yes","No")
     #ExprFrame4GraphAddBandSum.uni$colorLabel=ifelse(ExprFrame4GraphAddBandSum.uni$Band%in%ADLoci,"ADLoci",ifelse(ExprFrame4GraphAddBandSum.uni$Chr%in%c(paste0("chr",seq(1,21,2))),"Single","Double"))

     targetBand=c("1p36.33","1q21.1","11q14.2","14q24.3","15q11.2","15q13.1","15q13.3","16p13.2","16p13.3","17q21.1","17q21.31-q21.32","19q13.41","19q13.41","1q32.2","21q21.3","21q22.2","2p23.3","2q14.3","2q33.3","2q33.3","3p11.2","3p12.3","3q11.2","3q22.1","4p16.2","4q13.1","6p21.32","7p21.2","7q22.1","9p24.2","9p24.1","9p24.3")
     #Add the GWAS loci information on the graph, D:/Alzheimer/GWAS
     #GWAS=read.table("/data2/deng/Alzheimer/GWAS/CNVRegionBandInfor.txt",header=T)
     CNVRegion=intersect(targetBand,ExprFrame4GraphAddBandSum.uni$Band)
     ExprFrame4GraphAddBandSum.uni$CNVRegion=ifelse(ExprFrame4GraphAddBandSum.uni$Band%in%CNVRegion,"Yes","No")
     ExprFrame4GraphAddBandSum.uni$colorLabel=ifelse(ExprFrame4GraphAddBandSum.uni$Band%in%CNVRegion,"CNVRegion",ifelse(ExprFrame4GraphAddBandSum.uni$Chr%in%c(paste0("chr",seq(1,21,2))),"Single","Double"))


     t=ggplot(ExprFrame4GraphAddBandSum.uni, aes(index, CNVSum,fill=colorLabel,width=barWidth))+
       geom_bar(stat="identity",alpha=0.8)+
       scale_fill_manual(values = c("DarkOrange","DarkGray","LightGrey"))+
       theme_bw()+ 
       ylim(-3,6)+
       ggrepel::geom_label_repel(  data = ExprFrame4GraphAddBandSum.uni[(ExprFrame4GraphAddBandSum.uni$Band%in%CNVRegion&abs(ExprFrame4GraphAddBandSum.uni$CNVSum)>0.2),],
                   aes(x = index, y = CNVSum, label = Band),
                   size = 2,
                   colour="black",
                   max.overlaps=20,
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
                   fill = NA
                   )+  #no background
                   theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=rel(1.0),colour = "black"))
     pdf(paste0("/data2/deng/Aging/Ex/AllEx/InferCNV/ADRiskPlot/CNVRegionPlotFor_Cluster",targetCluster,"_",g,".pdf"),width=10,height=1.5)
     print(t)
     dev.off()
  }




for(g in unique(ExAll.integrated$Group)){
  setwd(paste0("/data2/deng/Aging/Ex/AllEx/InferCNV/",g,sep=""))
  ExData=read.table("infercnv.observations.txt",check.names=F,row.names=1,header=T)
  for(targetCluster in c(0,3,5)){
     SubType=subset(ExAll.integrated,seurat_clusters==targetCluster)
     #SubType=subset(ExAll.integrated,seurat_clusters==targetCluster) #only for the integrated dataset
     data=ExData

     SubTypeCells=intersect(colnames(SubType),colnames(data))
     SubTypeCellsExpr=data[,SubTypeCells]
     SubTypeCellMean=apply(SubTypeCellsExpr,1,mean)
     #obtain the SNV score, CNV score with 1 reprent no gain or loss
     ExprFrame=cbind("SubTypeCellMean"=SubTypeCellMean,"CNV"=(SubTypeCellMean-1))
     quantile(ExprFrame[,2],0.95) #C0 0.009026389
     rownames(ExprFrame)=rownames(data)

     #this gene anotation file was used to order the genes based on their chromosome location from chr1 to chr22
     annotation=read.table("/Projects/deng/refGene/gencode.v38.annotation.removeDupliSymbol.geneType.txt",header=TRUE,row.names=1)
     anno=data.frame("Chr"=annotation$Chr,"Start"=annotation$Start,"End"=annotation$End)
     rownames(anno)=annotation$Symbol
     anno=anno[anno$Chr%in%paste0("chr",c(1:22),sep=""),] #only consider the genes from chr1 to chr22, remove genes from chrM
     genes=intersect(rownames(data),rownames(anno))

     #to get the band information for each genes, this information were downloaded from HGNC database
     GeneChrLocation=read.table("/Projects/deng/Public/HGNC/GeneChrLocation.txt",header=T,row.names=1,sep="\t")
     GeneChrLocation=GeneChrLocation[,c("Symbol","ChLocation")]
     #GeneChrLocation$BandTotal=sub('(\\d+[p|q]\\d+)(.*)',"\\1", GeneChrLocation$ChLocation)
     GeneChrLocation$BandTotal=GeneChrLocation$ChLocation
     BandTotal=data.frame(table(GeneChrLocation$BandTotal)) #get the total number genes for each Band
     colnames(BandTotal)=c("Band","TotalBand")

     #order CNV matrix based on their gene order on the chromosome
     genes=intersect(genes,GeneChrLocation$Symbol) 
     length(genes)
     anno=anno[genes,]
     anno$Chr=factor(anno$Chr,levels=paste0("chr",c(1:22),sep=""))
     anno=anno[order(anno$Chr,anno$Star,anno$End),]
     geneOrder=intersect(rownames(anno),rownames(data))
     ExprFrame4Graph=data.frame(ExprFrame[geneOrder,])
     
     #Add the gene chromosome information to CNV score matrix
     ExprFrame4Graph$index=c(1:dim(ExprFrame4Graph)[1]) #build the index for each gene
     all(rownames(anno)==rownames(ExprFrame4Graph))
     ExprFrame4Graph=cbind(anno,ExprFrame4Graph)

     #Add the gene band information to CNV score matrix
     rownames(GeneChrLocation)=GeneChrLocation$Symbol
     GeneChrLocation=GeneChrLocation[genes,]
     all(rownames(GeneChrLocation)==rownames(ExprFrame4Graph)) #TRUE
     ExprFrame4Graph$Band=GeneChrLocation$ChLocation

     #got the CNV score for each band (as the sum score within each Band)
     BandSum=ExprFrame4Graph %>% group_by(Band) %>% summarise(CNVSum = sum(CNV))
     
     #ExprFrame4Graph$Symbol=rownames(ExprFrame4Graph)
     ExprFrame4GraphAddBandSum=merge(ExprFrame4Graph,BandSum,by="Band") #this si the total CNV score matrix (keep the gene information)
     
     ExprFrame4GraphAddBandSum.uni=ExprFrame4GraphAddBandSum[!duplicated(ExprFrame4GraphAddBandSum$Band),] #randomly seleced one band (not this step removed the gene information)
     ExprFrame4GraphAddBandSum.uni=ExprFrame4GraphAddBandSum.uni[order(ExprFrame4GraphAddBandSum.uni$index),c("Band","index","Chr","CNVSum")] #remove the gene dependent columns
     #obtain the width of the band based on their expressed gene number
     BandMidIndex=array()
     for(i in 2:(dim(ExprFrame4GraphAddBandSum.uni)[1]+1)){
         BandMidIndex[i-1]=ExprFrame4GraphAddBandSum.uni[i,"index"]-ExprFrame4GraphAddBandSum.uni[(i-1),"index"]
     }
     ExprFrame4GraphAddBandSum.uni$barWidth=BandMidIndex
     
     #Add the GWAS loci information on the graph, D:/Alzheimer/GWAS
     #GWAS=read.table("/data2/deng/Alzheimer/GWAS/ADLociBandInfor.txt",header=T)
     #ADLoci=intersect(GWAS$ChLocation,ExprFrame4GraphAddBandSum.uni$Band)
     #ExprFrame4GraphAddBandSum.uni$ADLoci=ifelse(ExprFrame4GraphAddBandSum.uni$Band%in%ADLoci,"Yes","No")
     #ExprFrame4GraphAddBandSum.uni$colorLabel=ifelse(ExprFrame4GraphAddBandSum.uni$Band%in%ADLoci,"ADLoci",ifelse(ExprFrame4GraphAddBandSum.uni$Chr%in%c(paste0("chr",seq(1,21,2))),"Single","Double"))

     targetBand=c("1p36.33","1q21.1","11q14.2","14q24.3","15q11.2","15q13.1","15q13.3","16p13.2","16p13.3","17q21.1","17q21.31-q21.32","19q13.41","19q13.41","1q32.2","21q21.3","21q22.2","2p23.3","2q14.3","2q33.3","2q33.3","3p11.2","3p12.3","3q11.2","3q22.1","4p16.2","4q13.1","6p21.32","7p21.2","7q22.1","9p24.2","9p24.1","9p24.3")
     #Add the GWAS loci information on the graph, D:/Alzheimer/GWAS
     #GWAS=read.table("/data2/deng/Alzheimer/GWAS/CNVRegionBandInfor.txt",header=T)
     CNVRegion=intersect(targetBand,ExprFrame4GraphAddBandSum.uni$Band)
     ExprFrame4GraphAddBandSum.uni$CNVRegion=ifelse(ExprFrame4GraphAddBandSum.uni$Band%in%CNVRegion,"Yes","No")
     ExprFrame4GraphAddBandSum.uni$colorLabel=ifelse(ExprFrame4GraphAddBandSum.uni$Band%in%CNVRegion,"CNVRegion",ifelse(ExprFrame4GraphAddBandSum.uni$Chr%in%c(paste0("chr",seq(1,21,2))),"Single","Double"))


     t=ggplot(ExprFrame4GraphAddBandSum.uni, aes(index, CNVSum,fill=colorLabel,width=barWidth))+
       geom_bar(stat="identity",alpha=0.8)+
       scale_fill_manual(values = c("DarkOrange","DarkGray","LightGrey"))+
       theme_bw()+ 
       ylim(-3,6)+
       ggrepel::geom_label_repel(  data = ExprFrame4GraphAddBandSum.uni[(ExprFrame4GraphAddBandSum.uni$Band%in%CNVRegion&abs(ExprFrame4GraphAddBandSum.uni$CNVSum)>0.2),],
                   aes(x = index, y = CNVSum, label = Band),
                   size = 2,
                   colour="black",
                   max.overlaps=20,
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
                   fill = NA
                   )+  #no background
                   theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=rel(1.0),colour = "black"))
     pdf(paste0("/data2/deng/Aging/Ex/AllEx/InferCNV/ADRiskPlot/CNVRegionPlotFor_Cluster",targetCluster,"_",g,".pdf"),width=10,height=1.5)
     print(t)
     dev.off()
  }
}


####https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5115612/





     SubType=subset(ExAll.integrated,Group==g&seurat_clusters==targetCluster)
     #SubType=subset(ExAll.integrated,seurat_clusters==targetCluster) #only for the integrated dataset
     data=ExData

     SubTypeCells=intersect(colnames(SubType),colnames(data))
     SubTypeCellsExpr=data[,SubTypeCells]
     SubTypeCellMean=apply(SubTypeCellsExpr,1,mean)
     #obtain the SNV score, CNV score with 1 reprent no gain or loss
     ExprFrame=cbind("SubTypeCellMean"=SubTypeCellMean,"CNV"=(SubTypeCellMean-1))
     quantile(ExprFrame[,2],0.95) #C0 0.009026389
     rownames(ExprFrame)=rownames(data)

     #this gene anotation file was used to order the genes based on their chromosome location from chr1 to chr22
     annotation=read.table("/Projects/deng/refGene/gencode.v38.annotation.removeDupliSymbol.geneType.txt",header=TRUE,row.names=1)
     anno=data.frame("Chr"=annotation$Chr,"Start"=annotation$Start,"End"=annotation$End)
     rownames(anno)=annotation$Symbol
     anno=anno[anno$Chr%in%paste0("chr",c(1:22),sep=""),] #only consider the genes from chr1 to chr22, remove genes from chrM
     genes=intersect(rownames(data),rownames(anno))

     #to get the band information for each genes, this information were downloaded from HGNC database
     GeneChrLocation=read.table("/Projects/deng/Public/HGNC/GeneChrLocation.txt",header=T,row.names=1,sep="\t")
     GeneChrLocation=GeneChrLocation[,c("Symbol","ChLocation")]
     #GeneChrLocation$BandTotal=sub('(\\d+[p|q]\\d+)(.*)',"\\1", GeneChrLocation$ChLocation)
     GeneChrLocation$BandTotal=GeneChrLocation$ChLocation
     BandTotal=data.frame(table(GeneChrLocation$BandTotal)) #get the total number genes for each Band
     colnames(BandTotal)=c("Band","TotalBand")

     #order CNV matrix based on their gene order on the chromosome
     genes=intersect(genes,GeneChrLocation$Symbol) 
     length(genes)
     anno=anno[genes,]
     anno$Chr=factor(anno$Chr,levels=paste0("chr",c(1:22),sep=""))
     anno=anno[order(anno$Chr,anno$Star,anno$End),]
     geneOrder=intersect(rownames(anno),rownames(data))
     ExprFrame4Graph=data.frame(ExprFrame[geneOrder,])
     
     #Add the gene chromosome information to CNV score matrix
     ExprFrame4Graph$index=c(1:dim(ExprFrame4Graph)[1]) #build the index for each gene
     all(rownames(anno)==rownames(ExprFrame4Graph))
     ExprFrame4Graph=cbind(anno,ExprFrame4Graph)

     #Add the gene band information to CNV score matrix
     rownames(GeneChrLocation)=GeneChrLocation$Symbol
     GeneChrLocation=GeneChrLocation[genes,]
     all(rownames(GeneChrLocation)==rownames(ExprFrame4Graph)) #TRUE
     ExprFrame4Graph$Band=GeneChrLocation$ChLocation

     #got the CNV score for each band (as the sum score within each Band)
     BandSum=ExprFrame4Graph %>% group_by(Band) %>% summarise(CNVSum = sum(CNV))
     
     #ExprFrame4Graph$Symbol=rownames(ExprFrame4Graph)
     ExprFrame4GraphAddBandSum=merge(ExprFrame4Graph,BandSum,by="Band") #this si the total CNV score matrix (keep the gene information)
     
     ExprFrame4GraphAddBandSum.uni=ExprFrame4GraphAddBandSum[!duplicated(ExprFrame4GraphAddBandSum$Band),] #randomly seleced one band (not this step removed the gene information)
     ExprFrame4GraphAddBandSum.uni=ExprFrame4GraphAddBandSum.uni[order(ExprFrame4GraphAddBandSum.uni$index),c("Band","index","Chr","CNVSum")] #remove the gene dependent columns
     #obtain the width of the band based on their expressed gene number
     BandMidIndex=array()
     for(i in 2:(dim(ExprFrame4GraphAddBandSum.uni)[1]+1)){
         BandMidIndex[i-1]=ExprFrame4GraphAddBandSum.uni[i,"index"]-ExprFrame4GraphAddBandSum.uni[(i-1),"index"]
     }
     ExprFrame4GraphAddBandSum.uni$barWidth=BandMidIndex
     
     targetBand=c("1p36.33","1q21.1","11q14.2","14q24.3","15q11.2","15q13.1","15q13.3","16p13.2","16p13.3","17q21.1","17q21.31-q21.32","19q13.41","19q13.41","1q32.2","21q21.3","21q22.2","2p23.3","2q14.3","2q33.3","2q33.3","3p11.2","3p12.3","3q11.2","3q22.1","4p16.2","4q13.1","6p21.32","7p21.2","7q22.1","9p24.2","9p24.1","9p24.3")

     #Add the GWAS loci information on the graph, D:/Alzheimer/GWAS
     #GWAS=read.table("/data2/deng/Alzheimer/GWAS/CNVRegionBandInfor.txt",header=T)
     CNVRegion=intersect(targetBand,ExprFrame4GraphAddBandSum.uni$Band)
     ExprFrame4GraphAddBandSum.uni$CNVRegion=ifelse(ExprFrame4GraphAddBandSum.uni$Band%in%CNVRegion,"Yes","No")
     ExprFrame4GraphAddBandSum.uni$colorLabel=ifelse(ExprFrame4GraphAddBandSum.uni$Band%in%CNVRegion,"CNVRegion",ifelse(ExprFrame4GraphAddBandSum.uni$Chr%in%c(paste0("chr",seq(1,21,2))),"Single","Double"))


     #ExprFrame4GraphAddBandSum$ADGene=ifelse(ExprFrame4GraphAddBandSum$Symbol%in%GWAS$Symbol,"ADRisk","NormalGene")
     #ExprFrame4GraphAddBandSum=ExprFrame4GraphAddBandSum[order(ExprFrame4GraphAddBandSum$index),]
     #write.table(ExprFrame4GraphAddBandSum,file="/data2/deng/Aging/Ex/AllEx/InferCNV/ADRiskPlot/ExprFrame4GraphAddBandSum.txt",sep="\t",quote=F,row.names=F)
     
     t=ggplot(ExprFrame4GraphAddBandSum.uni, aes(index, CNVSum,fill=colorLabel,width=barWidth))+
       geom_bar(stat="identity",alpha=0.8)+
       scale_fill_manual(values = c("DarkOrange","DarkGray","LightGrey"))+
       theme_bw()+ 
       ylim(-3,6)+
       ggrepel::geom_label_repel(  data = ExprFrame4GraphAddBandSum.uni[(ExprFrame4GraphAddBandSum.uni$Band%in%CNVRegion&abs(ExprFrame4GraphAddBandSum.uni$CNVSum)>0.2),],
                   aes(x = index, y = CNVSum, label = Band),
                   size = 2,
                   colour="black",
                   max.overlaps=20,
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
                   fill = NA
                   )+  #no background
                   theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=rel(1.0),colour = "black"))
     pdf(paste0("/data2/deng/Aging/Ex/AllEx/InferCNV/ADRiskPlot/CNVRegionPlotFor_Cluster",targetCluster,"_",g,".pdf"),width=10,height=1.5)
     print(t)
     dev.off()
