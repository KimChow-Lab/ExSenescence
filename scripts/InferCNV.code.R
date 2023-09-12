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



####################################################################################
############     Ribon plot for each cluster across  choromosome      ##############
####################################################################################

# the output figures were showen in Figure 1F, 2A, 2B, 2C, S1E, S2G
ExAll.integrated=readRDS("ExAll.integratedTSNE.rds")
ExAll.integrated=subset(ExAll.integrated,idents=c(0:15))

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
}
}



####################################################################################
############        Arm-level aneuploidy estimated     ##############
####################################################################################

setwd("/Projects/deng/Aging/Ex/AllControlEx/InferCNV")
ExData=read.table("/Projects/deng/Aging/Ex/AllControlEx/InferCNV/infercnv.observations.txt",check.names=F,row.names=1,header=T)
ExCtrl=readRDS("/Projects/deng/Aging/Ex/AllControlEx/ExCtrl.integrated.rds")
data=ExData
#obtain the reference value
AllCellMean=apply(data,1,mean) #check the background CNA infor
summary(AllCellMean)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.9972  1.0015  1.0033  1.0035  1.0057  1.0088
quantile(AllCellMean,0.95) #C0 1.007714
quantile(AllCellMean,0.05) #C0 0.9994199

SubType=subset(ExCtrl,seurat_clusters==0)
SubTypeCells=intersect(colnames(SubType),colnames(data))
SubTypeCellsExpr=data[,SubTypeCells]
SubTypeCellMean=apply(SubTypeCellsExpr,1,mean)
ExprFrame=cbind("SubTypeCellMean"=SubTypeCellMean,"CNV"=(SubTypeCellMean-1))
summary(ExprFrame[,2])

cellNumber=dim(SubTypeCellsExpr)[2] #7094
geneNumber=dim(SubTypeCellsExpr)[1] #8109
Ratio=matrix(data = NA, nrow = geneNumber, ncol = 2, dimnames = list(rownames(ExprFrame), c("Gain","Loss")))
for(i in 1:dim(ExprFrame)[1]){
  Ratio[i,1]=sum(SubTypeCellsExpr[i,]>1.007714)  #the top 95%
  Ratio[i,2]=sum(SubTypeCellsExpr[i,]<0.9994199) #the bottom 5%
  if(!i%%100){
    print(i)
  }
}
tmp=Ratio
Ratio=data.frame(Ratio)
Ratio$GainRatio=Ratio$Gain/cellNumber
Ratio$LossRatio=Ratio$Loss/cellNumber

GainGene=Ratio[Ratio$GainRatio>0.25,]#SubTypeCellsExpr=SubTypeCellsExpr-1
LossGene=Ratio[Ratio$LossRatio>0.25,]
write.table(GainGene,file="CNV8Arm/GainGeneC0.txt",sep="\t",quote=F)
write.table(LossGene,file="CNV8Arm/LossGeneC0.txt",sep="\t",quote=F)


head(SubTypeCellMean)
ExprFrame=cbind("GainRatio"=Ratio$GainRatio,"LossRatio"=Ratio$LossRatio)
rownames(ExprFrame)=rownames(Ratio)
GeneChrLocation=read.table("/Projects/deng/Public/HGNC/GeneChrLocation.txt",header=T,row.names=1,sep="\t")
GeneChrLocation=GeneChrLocation[,c("Symbol","ChLocation")]
GeneChrLocation$BandTotal=sub('(\\d+[p|q])(\\d+.*)',"\\1", GeneChrLocation$ChLocation)
BandTotal=data.frame(table(GeneChrLocation$BandTotal))
colnames(BandTotal)=c("Band","TotalBand")
rownames(GeneChrLocation)=GeneChrLocation$Symbol
gene=intersect(rownames(ExprFrame),rownames(GeneChrLocation))
ExprFrame4Graph=ExprFrame[gene,]
GeneChrLocation=GeneChrLocation[gene,]
all(rownames(ExprFrame4Graph)==rownames(GeneChrLocation))
ExprFrame4Graph=data.frame(ExprFrame4Graph)
ExprFrame4Graph$Band=GeneChrLocation$BandTotal

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
pdf("CNVsplit8Arm4Cluster0.pdf",width=10)
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

  }
}
