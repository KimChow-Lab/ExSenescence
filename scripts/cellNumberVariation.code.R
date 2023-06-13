library(Seurat)
library(ggplot2)
library(ggpubr)
setwd("/Projects/deng/Aging/Ex/AllEx/cellNumberVariation")
ExAll.integrated=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")
ExAll.integrated=subset(ExAll.integrated,idents=c(0:15))

C3C5=subset(ExAll.integrated,idents=c(3,5))
C3C5Number4EachSample=data.frame(table(C3C5$orig.ident))
SampleTotalNumber=data.frame(table(ExAll.integrated$orig.ident))
C3C5Ratio=merge(C3C5Number4EachSample,SampleTotalNumber,by="Var1")
colnames(C3C5Ratio)=c("Sample","C3C5Number","Total")
C3C5Ratio$Ratio=C3C5Ratio$C3C5Number/C3C5Ratio$Total
SampleAge=data.frame("Sample"=ExAll.integrated$orig.ident,"Age"=ExAll.integrated$Age,"Statues"=ExAll.integrated$Statues)
SampleAge=unique(SampleAge)
C3C5RatioAge=merge(C3C5Ratio,SampleAge,by="Sample")
C3C5RatioAge=C3C5RatioAge[C3C5RatioAge$Age>69,]
ClusterInfo=C3C5RatioAge
dim(ClusterInfo[ClusterInfo$Ratio>0.4,]) #10 samples 1 C5, 9 C3
ClusterInfo=ClusterInfo[ClusterInfo$Ratio<0.4,]
ClusterInfo$Statues=factor(ClusterInfo$Statues,levels=c("Control","Alzheimer"))
t=ggplot(ClusterInfo, aes(x=as.numeric(Age), y=as.numeric(Ratio), color=Statues,fill=Statues)) +
  geom_point() + 
  geom_smooth(method="lm",aes(fill = Statues),alpha = 0.2)+stat_cor(data=ClusterInfo, method = "pearson")+
  scale_color_manual(values=c("#28A9A1", "#E71F19"))+
  scale_fill_manual(values=c("#28A9A1", "#E71F19"))+theme_bw()
pdf("C3C5CombineAgeRatioInStatues.pdf",height=4,width=6)
print(t)
dev.off()
write.table(ClusterInfo,file="C3C5CombineAgeRatioInStatues.txt",sep="\t",quote=F,row.names=F)

ClusterInfo$AgeBin=ifelse(ClusterInfo$Age<70,"60-69",ifelse(ClusterInfo$Age<80,"70-79",ifelse(ClusterInfo$Age<90,"80-89","90-99")))
t=ggplot(ClusterInfo, aes(x=AgeBin, y=Ratio, colour=Statues)) +
  geom_boxplot(outlier.color = "white") + ylim(0,0.4)+
  geom_point(aes(fill=Statues),size=3, show.legend=TRUE, alpha=0.5, position = position_jitterdodge(dodge.width = 0.55))+
  scale_color_manual(values=c("RoyalBlue", "DarkViolet"))+
  scale_fill_manual(values=c("RoyalBlue", "DarkViolet"))+
  #stat_compare_means(comparisons = my_comparisons,method = "t.test")+
  theme_bw()+theme(axis.title.x = element_blank(),axis.title.y = element_blank())
pdf("C3C5CombineAgeBinRatioInStatues.pdf",height=4,width=5)
print(t)
dev.off()


tmp=paste0("C",ExAll.integrated$seurat_clusters,"_",ExAll.integrated$Gender,"_",ExAll.integrated$Age,"_",ExAll.integrated$orig.ident,"_",ExAll.integrated$Source,"_",ExAll.integrated$Statues,sep="")
tmp=as.matrix(table(tmp))
TmpInfo=as.matrix(do.call(rbind, strsplit(as.character(rownames(tmp)),'_')))
result=data.frame(SampleID=TmpInfo[,4],Cluster=TmpInfo[,1],"Gender"=TmpInfo[,2],"Age"=TmpInfo[,3],"Source"=TmpInfo[,5],"Statues"=TmpInfo[,6],"Number"=tmp[,1])
result=result[order(result$Statues,result$Gender,result$Age,result$SampleID),]
sampleList=factor(result$SampleID,levels=unique(result$SampleID))
result=result[result$Age>69,]
Cluster_order=factor(result$Cluster,levels=paste0("C",c(5,3,c(0:2),4,c(6:15))))
sampleNumber=data.frame(table(ExAll.integrated$orig.ident))
colnames(sampleNumber)=c("SampleID","SampleNumber")
resultTmp=merge(result,sampleNumber,by="SampleID")
resultTmp$Ratio=resultTmp$Number/resultTmp$SampleNumber
resultTmp=resultTmp[order(resultTmp$Statues,resultTmp$Gender,resultTmp$Age,resultTmp$SampleID),]


C3C5Info=resultTmp[resultTmp$Cluster %in% c("C3","C5"),]
ClusterInfo=C3C5Info
dim(ClusterInfo[ClusterInfo$Ratio>0.4,]) #10 samples 1 C5, 9 C3
ClusterInfo=ClusterInfo[ClusterInfo$Ratio<0.4,]
ClusterInfo$Statues=factor(ClusterInfo$Statues,levels=c("Control","Alzheimer"))

t=ggplot(ClusterInfo, aes(x=as.numeric(Age), y=as.numeric(Ratio), color=Statues,fill=Statues)) +
  geom_point() + 
  facet_grid(~Cluster,scales="free_x")+
  geom_smooth(method="lm",aes(fill = Statues),alpha = 0.2)+stat_cor(data=ClusterInfo, method = "pearson")+
  scale_color_manual(values=c("#28A9A1", "#E71F19"))+
  scale_fill_manual(values=c("#28A9A1", "#E71F19"))+theme_bw()
pdf("C3C5AgeRatioInStatues.pdf",height=4,width=8)
print(t)
dev.off()

ClusterInfo$AgeBin=ifelse(ClusterInfo$Age<70,"60-69",ifelse(ClusterInfo$Age<80,"70-79",ifelse(ClusterInfo$Age<90,"80-89","90-99")))
t=ggplot(ClusterInfo, aes(x=AgeBin, y=Ratio, colour=Statues)) +
  geom_boxplot(outlier.color = "white") + 
  geom_point(aes(fill=Statues),size=3, show.legend=TRUE, alpha=0.5, position = position_jitterdodge(dodge.width = 0.75))+
  scale_color_manual(values=c("RoyalBlue", "DarkViolet"))+
  scale_fill_manual(values=c("RoyalBlue", "DarkViolet"))+
  facet_wrap(~Cluster,scales="free_y")+
  theme_bw()+theme(axis.title.x = element_blank(),axis.title.y = element_blank())
pdf("C3C5SplitAgeBinRatioInStatues.pdf",height=4,width=8)
print(t)
dev.off()

wilcox.test(ClusterInfo[ClusterInfo$Cluster=="C3"&ClusterInfo$Statues=="Alzheimer"&ClusterInfo$AgeBin=="70-79",]$Ratio,ClusterInfo[ClusterInfo$Cluster=="C3"&ClusterInfo$Statues=="Control"&ClusterInfo$AgeBin=="70-79",]$Ratio)
p-value = 0.04798
wilcox.test(ClusterInfo[ClusterInfo$Cluster=="C3"&ClusterInfo$Statues=="Alzheimer"&ClusterInfo$AgeBin=="80-89",]$Ratio,ClusterInfo[ClusterInfo$Cluster=="C3"&ClusterInfo$Statues=="Control"&ClusterInfo$AgeBin=="80-89",]$Ratio)
p-value = 0.2181
wilcox.test(ClusterInfo[ClusterInfo$Cluster=="C3"&ClusterInfo$Statues=="Alzheimer"&ClusterInfo$AgeBin=="90-99",]$Ratio,ClusterInfo[ClusterInfo$Cluster=="C3"&ClusterInfo$Statues=="Control"&ClusterInfo$AgeBin=="90-99",]$Ratio)
p-value = 0.1259


wilcox.test(ClusterInfo[ClusterInfo$Cluster=="C3"&ClusterInfo$Statues=="Alzheimer"&ClusterInfo$AgeBin=="70-79",]$Ratio,ClusterInfo[ClusterInfo$Cluster=="C3"&ClusterInfo$Statues=="Alzheimer"&ClusterInfo$AgeBin=="90-99",]$Ratio)


C3C5Info=resultTmp[resultTmp$Cluster %in% c("C0","C3","C5"),]
ClusterInfo=C3C5Info
ClusterInfo=ClusterInfo[ClusterInfo$Ratio<0.4,]
ClusterInfo=ClusterInfo[ClusterInfo$Age>60,]
ClusterInfo$Group=paste0(ClusterInfo$Gender,"_",ClusterInfo$Statues,sep="")
StatuesList=factor(ClusterInfo$Statues,levels=c("Control","Alzheimer"))
t=ggplot(ClusterInfo, aes(x=StatuesList, y=Ratio, colour=StatuesList)) +
  geom_boxplot(outlier.color = "white") + ylim(0,0.4)+
  geom_point(aes(fill=StatuesList),size=3, show.legend=TRUE, alpha=0.5, position = position_jitterdodge(dodge.width = 0.25))+
  scale_color_manual(values=c("RoyalBlue", "DarkViolet"))+
  scale_fill_manual(values=c("RoyalBlue", "DarkViolet"))+
  #stat_compare_means(comparisons = my_comparisons,method = "t.test")+
  theme_bw()+facet_wrap(.~Cluster,scales="free_x")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
pdf("C0C3C5RatioBetStatuesPvalue.pdf",height=4,width=8)
print(t)
dev.off()
write.table(C3C5Info,file="C0C3C5RatioBetStatues.txt",sep="\t",quote=F,row.names=F)


ClusterInfo=resultTmp
ClusterInfo=ClusterInfo[ClusterInfo$Ratio<0.4,]
ClusterInfo=ClusterInfo[ClusterInfo$Age>60,]
ClusterInfo$Group=paste0(ClusterInfo$Gender,"_",ClusterInfo$Statues,sep="")
StatuesList=factor(ClusterInfo$Statues,levels=c("Control","Alzheimer"))
ClusterInfo$Cluster=factor(ClusterInfo$Cluster,levels=paste0("C",c(0:15),sep=""))
t=ggplot(ClusterInfo, aes(x=StatuesList, y=Ratio, colour=StatuesList)) +
  geom_boxplot(outlier.color = "white") + ylim(0,0.4)+
  geom_point(aes(fill=StatuesList),size=3, show.legend=TRUE, alpha=0.5, position = position_jitterdodge(dodge.width = 0.25))+
  scale_color_manual(values=c("RoyalBlue", "DarkViolet"))+
  scale_fill_manual(values=c("RoyalBlue", "DarkViolet"))+
  #stat_compare_means(comparisons = list(c("Alzheimer","Control")),method = "t.test")+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  theme_bw()+facet_wrap(.~Cluster,scales="free_y",ncol=8)+
  #theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank())
pdf("CellRatioBetStatuesPvalue.pdf",height=4,width=15)
print(t)
dev.off()
write.table(ClusterInfo,file="AllClustersRatioBetStatues.txt",sep="\t",quote=F,row.names=F)



metaTmp=ExAll.integrated@meta.data
PhenoType=read.table("/Projects/deng/Aging/Ex/AllEx/SampleInformation/PhenotypeInformation.txt",header=T)
#BraakStage,Cogdx,Ceradsc,Plaq_n,NFT,Amyloid,Gpath,Tangles,Cogn_global_lv
result=merge(metaTmp, PhenoType,by="orig.ident")
tmp=paste0("C",result$seurat_clusters,"_",result$Gender,"_",result$orig.ident,"_",result$Statues,"_",result$Age,"_",result$Ceradsc,sep="")
tmp=as.matrix(table(tmp))
TmpInfo=as.matrix(do.call(rbind, strsplit(as.character(rownames(tmp)),'_')))
resultTmp=data.frame(orig.ident=TmpInfo[,3],Cluster=TmpInfo[,1],"Gender"=TmpInfo[,2],"Statues"=TmpInfo[,4],"PhenoType"=TmpInfo[,6],"Age"=TmpInfo[,5],"Number"=tmp[,1])
sampleNumber=data.frame(table(result$orig.ident))
colnames(sampleNumber)=c("orig.ident","SampleNumber")
resultTmp=merge(resultTmp,sampleNumber,by="orig.ident")
resultTmp$Ratio=resultTmp$Number/resultTmp$SampleNumber

#Cogdx
C3C5Info=resultTmp[resultTmp$Cluster %in% c("C0","C3","C5"),]
ClusterInfo=C3C5Info
ClusterInfo=ClusterInfo[ClusterInfo$Ratio<0.4,]
ClusterInfo=ClusterInfo[ClusterInfo$Age>60,]
ClusterInfo=ClusterInfo[ClusterInfo$PhenoType %in% names(which(table(ClusterInfo$PhenoType)>5)),]
t=ggplot(ClusterInfo, aes(x=PhenoType, y=Ratio, colour=PhenoType)) +
  geom_boxplot(outlier.color = "white") + 
  geom_point(aes(fill=PhenoType),size=3, show.legend=TRUE, alpha=0.5, position = position_jitterdodge(dodge.width = 0.25))+
  theme_bw()+facet_wrap(~Cluster,scales="free_y",ncol=8)+
  scale_color_manual(values=c("#28A9A1","#C9A77C","#F4A016","#E71F19"))+
  scale_fill_manual(values=c("#28A9A1","#C9A77C","#F4A016","#E71F19"))+
  #stat_compare_means(method = "wilcox.test",comparisons = list(c(1, 4)),label.y.npc=0.8,label = "p.format")+
  #scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank())
pdf("C0C3C5CogdxRatioIn.pdf",height=4,width=8)
print(t)
dev.off()
write.table(ClusterInfo,file="C0C3C5CogdxRatio.txt",sep="\t",quote=F,row.names=F)


ClusterOrder=paste0("C",c(0:15))
ClusterInfo=resultTmp
dim(ClusterInfo[ClusterInfo$Ratio>0.4,]) #10 samples 1 C5, 9 C3
ClusterInfo=ClusterInfo[ClusterInfo$Ratio<0.4,]
ClusterInfo=ClusterInfo[ClusterInfo$PhenoType %in% names(which(table(ClusterInfo$PhenoType)>16)),]
ClusterInfo$Cluster=factor(ClusterInfo$Cluster,levels=ClusterOrder)
ClusterInfo=ClusterInfo[order(ClusterInfo$Cluster),]
write.table(ClusterInfo,file="AllClusterCogdxRatio.txt",sep="\t",quote=F,row.names=F)

#Cogdx
t=ggplot(ClusterInfo, aes(x=PhenoType, y=Ratio, colour=PhenoType)) +
  geom_boxplot(outlier.color = "white") + 
  geom_point(aes(fill=PhenoType),size=3, show.legend=TRUE, alpha=0.5, position = position_jitterdodge(dodge.width = 0.25))+
  theme_bw()+facet_wrap(~factor(Cluster,levels=ClusterOrder),scales="free_y",ncol=8)+
  scale_color_manual(values=c("#28A9A1","#C9A77C","#F4A016","#E71F19"))+
  scale_fill_manual(values=c("#28A9A1","#C9A77C","#F4A016","#E71F19"))+
  #stat_compare_means(method = "wilcox.test",comparisons = list(c(1, 4)),label.y.npc=0.8,label = "p.format")+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  #theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank())
pdf("CogdxRatioAllYAxis.pdf",width=15,height=4)
print(t)
dev.off()

#ceradx
ClusterInfo=resultTmp[resultTmp$Cluster %in% c("C0","C3","C5"),]
ClusterInfo=resultTmp
t=ggplot(ClusterInfo, aes(x=factor(PhenoType,levels=c(4,3,2,1)), y=Ratio, colour=PhenoType)) +
  geom_boxplot(outlier.color = "white") + 
  geom_point(aes(fill=PhenoType),size=3, show.legend=TRUE, alpha=0.5, position = position_jitterdodge(dodge.width = 0.25))+
  theme_bw()+facet_wrap(~factor(Cluster,levels=ClusterOrder),scales="free_y",ncol=8)+
  scale_color_manual(values=rev(c("#28A9A1","#C9A77C","#F4A016","#E71F19")))+
  scale_fill_manual(values=rev(c("#28A9A1","#C9A77C","#F4A016","#E71F19")))+
  stat_compare_means(method = "wilcox.test",comparisons = list(c(1, 4)),label.y.npc=0.8,label = "p.format")+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  #theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank())
pdf("C0C3C5CeradxRatio.pdf",height=4,width=8)
print(t)
dev.off()
write.table(ClusterInfo,file="C0C3C5CeradRatio.txt",sep="\t",quote=F,row.names=F)


#braaksc
t=ggplot(ClusterInfo, aes(x=PhenoType, y=Ratio, colour=PhenoType)) +
  geom_boxplot(outlier.color = "white") + 
  geom_point(aes(fill=PhenoType),size=3, show.legend=TRUE, alpha=0.5, position = position_jitterdodge(dodge.width = 0.25))+
  theme_bw()+facet_wrap(~factor(Cluster,levels=ClusterOrder),scales="free_y",ncol=8)+
  scale_color_manual(values=c("#046586","#28A9A1","#C9A77C","#F4A016","#F6BBC6","#E71F19"))+
  scale_fill_manual(values=c("#046586","#28A9A1","#C9A77C","#F4A016","#F6BBC6","#E71F19"))+
  #stat_compare_means(method = "wilcox.test",comparisons = list(c(1, 6)),label.y.npc=0.8,label = "p.format")+
  #scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  #theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank())
pdf("BraakscRatioYAxis.pdf",width=15,height=4)
print(t)
dev.off()



#apoe genotype
metaTmp=ExAll.integrated@meta.data
APOE=read.table("/Projects/deng/Aging/Ex/AllEx/SampleInformation/APOEInformation.txt",header=T)
result=merge(metaTmp, APOE,by="orig.ident")
tmp=paste0("C",result$seurat_clusters,"_",result$Gender,"_",result$orig.ident,"_",result$Statues,"_",result$Age,"_",result$APOE,sep="")
tmp=as.matrix(table(tmp))
TmpInfo=as.matrix(do.call(rbind, strsplit(as.character(rownames(tmp)),'_')))
resultTmp=data.frame(orig.ident=TmpInfo[,3],Cluster=TmpInfo[,1],"Gender"=TmpInfo[,2],"Statues"=TmpInfo[,4],"APOE"=TmpInfo[,6],"Age"=TmpInfo[,5],"Number"=tmp[,1])

sampleNumber=data.frame(table(result$orig.ident))
colnames(sampleNumber)=c("orig.ident","SampleNumber")
resultTmp=merge(resultTmp,sampleNumber,by="orig.ident")
resultTmp$Ratio=resultTmp$Number/resultTmp$SampleNumber
resultTmp=resultTmp[resultTmp$Ratio<0.4,]
ClusterOrder=paste0("C",c(0:15))
ClusterInfo=resultTmp
dim(ClusterInfo[ClusterInfo$Ratio>0.4,]) #10 samples 1 C5, 9 C3
ClusterInfo=ClusterInfo[ClusterInfo$APOE %in% c("23","33","34","44"),]
ClusterOrder=paste0("C",c(0:15))
#APOE
t=ggplot(ClusterInfo, aes(x=APOE, y=Ratio, colour=APOE)) +
  geom_boxplot(outlier.color = "white") + 
  geom_point(aes(fill=APOE),size=3, show.legend=TRUE, alpha=0.5, position = position_jitterdodge(dodge.width = 0.25))+
  theme_bw()+facet_wrap(~factor(Cluster,levels=ClusterOrder),scales="free_y",ncol=8)+
  scale_color_manual(values=c("#28A9A1","#C9A77C","#F4A016","#E71F19"))+
  scale_fill_manual(values=c("#28A9A1","#C9A77C","#F4A016","#E71F19"))+
  #stat_compare_means(method = "t.test",comparisons = list(c("23", "44")),label.y.npc=0.8,label = "p.format")+
  #scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank())
pdf("APOEGenotypeAllNoLabel.pdf",width=15,height=4)
print(t)
dev.off()
ClusterInfo$Cluster=factor(ClusterInfo$Cluster,levels=ClusterOrder)
ClusterInfo=ClusterInfo[order(ClusterInfo$Cluster),]
write.table(ClusterInfo,file="AllClusterAPOEGenotypeRatio.txt",sep="\t",quote=F,row.names=F)



ExAll.integrated=readRDS("/Projects/deng/Aging/Ex/AllEx/ExAll.integratedTSNE.rds")

ExTmp <- FindClusters(ExAll.integrated, resolution = c(seq(0.1,1,0.1)))
library(clustree)
clus.tree.out <- clustree(ExTmp)
pdf(file = "ExAll.integrated.Resolution.Tree.pdf", width = 12, height = 10)
print(clus.tree.out)
dev.off()


ExAllHighR=FindClusters(ExAll.integrated, resolution = 1)
tiff("ExAllHighRLabel.tiff",width=300,height=270)
DimPlot(ExAllHighR, raster=FALSE,reduction = "tsne",label=TRUE,label.size = 6)&NoLegend()&theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
tiff("ExAllHighRNoLabel.tiff",width=300,height=270)
DimPlot(ExAllHighR, raster=FALSE,reduction = "tsne")&NoLegend()&theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

tmp=paste0("C",ExAllHighR$seurat_clusters,"_",ExAllHighR$Gender,"_",ExAllHighR$Age,"_",ExAllHighR$orig.ident,"_",ExAllHighR$Source,"_",ExAllHighR$Statues,sep="")
tmp=as.matrix(table(tmp))
TmpInfo=as.matrix(do.call(rbind, strsplit(as.character(rownames(tmp)),'_')))
result=data.frame(SampleID=TmpInfo[,4],Cluster=TmpInfo[,1],"Gender"=TmpInfo[,2],"Age"=TmpInfo[,3],"Source"=TmpInfo[,5],"Statues"=TmpInfo[,6],"Number"=tmp[,1])
result=result[order(result$Statues,result$Gender,result$Age,result$SampleID),]

sampleNumber=data.frame(table(ExAllHighR$orig.ident))
colnames(sampleNumber)=c("SampleID","SampleNumber")
resultTmp=merge(result,sampleNumber,by="SampleID")
resultTmp$Ratio=resultTmp$Number/resultTmp$SampleNumber
resultTmp=resultTmp[order(resultTmp$Statues,resultTmp$Gender,resultTmp$Age,resultTmp$SampleID),]
C1Info=resultTmp[resultTmp$Cluster %in% c("C1","C10"),]
ClusterInfo=C1Info
ClusterInfo=ClusterInfo[ClusterInfo$Ratio<0.4,]
ClusterInfo=ClusterInfo[ClusterInfo$Age>60,]
#ClusterInfo=ClusterInfo[ClusterInfo$Cluster %in% "C3",]
ClusterInfo$Group=ClusterInfo$Statues
GroupList=factor(ClusterInfo$Group,levels=c("Control","Alzheimer"))
t=ggplot(ClusterInfo, aes(x=GroupList, y=Ratio, colour=GroupList)) +
  geom_boxplot(outlier.color = "white") + 
  geom_point(aes(fill=GroupList),size=3, show.legend=TRUE, alpha=0.5, position = position_jitterdodge(dodge.width = 0.25))+
  scale_color_manual(values=c("RoyalBlue", "DarkViolet"))+
  scale_fill_manual(values=c("RoyalBlue", "DarkViolet"))+
  stat_compare_means(comparisons = list(c("Control","Alzheimer")),method = "t.test")+
  theme_bw()+facet_wrap(~Cluster,scales="free_x")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
pdf("C3SubTypeRatioBetGroupPvalue.pdf",height=4,width=6)
print(t)
dev.off()

