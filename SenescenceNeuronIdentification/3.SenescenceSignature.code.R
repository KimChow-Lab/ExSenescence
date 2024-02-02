library(Seurat)
library(fgsea)
library(dplyr)
library(tibble)

setwd("/boot3/bixm/08_HNLongevity/pipline")
#using your own pipine to generated this neuron dataset
ExIntegrated=readRDS("Data/ExIntegrated.Downsampling.rds")
#set the default assay as RNA if you using the integrated method
DefaultAssay(ExIntegrated)="RNA"

t=DotPlot(ExIntegrated,features=rownames(ExIntegrated))
ExExpr8DotPlot=na.omit(t$data)

CSFromSegura=read.table("sources/55markersBySegura.txt",header=T)
CSFromCellAge=read.table("sources/GeneListFromCellAge.txt",header=T)
CSFromCellAge$Group="CSFromCellAge"
CSFromGabrielCasella=read.table("sources/ConservedGene8GabrielCasella.txt",header=T)
CSFromBinZhang=read.table("sources/Senescence8BinZhang.txt",header=T)
CSFromBinZhang$Group="CSFromBinZhang"
CSFromKEGG=read.table("sources/hsa04218_gene.txt",header=T)
CSFromReactome=read.table("sources/R-HSA-2559582_gene.txt",header=T)
SASPFromSenMayo=read.table("sources/SenMayoList.txt",header=T)

CombineCSList=rbind(CSFromSegura,CSFromCellAge,CSFromGabrielCasella,CSFromBinZhang,CSFromKEGG,CSFromReactome,SASPFromSenMayo)
CombineCSList=CombineCSList[CombineCSList$GeneSymbol%in%ExExpr8DotPlot$features.plot,]
CombineCSList=split(CombineCSList$GeneSymbol,CombineCSList$Group)

#define the target cluster 
cluster="5"
clusterCell<- ExExpr8DotPlot %>% dplyr::filter(id == cluster) %>% arrange(desc(avg.exp.scaled)) %>% dplyr::select(features.plot, avg.exp.scaled)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(CombineCSList, stats = ranks,eps=0)
ranks=na.omit(ranks)
fgseaRes=fgseaRes[order(fgseaRes$NES,decreasing=T),]

#generate multiple results as a table graph
pdf("ResultGraph/CSGeneList8fgsea.pdf",width=10,height=5)
plotGseaTable(CombineCSList[fgseaRes$pathway], ranks, fgseaRes, gseaParam=3)
dev.off()

#output the p value for each pathway
p<-ggplot(data=fgseaRes, aes(x=factor(pathway,levels=rev(pathway)), y=-log10(pval),fill=-log10(pval))) +
  scale_fill_gradient(low = "MistyRose", high = "red")+ theme_bw()+
  geom_bar(stat="identity",width=0.6)+coord_flip()
pdf("ResultGraph/CSGeneList8fgseaPvalue.pdf",width=5,height=5)
print(p)
dev.off()


#output single pathway
pdf("ResultGraph/Reactome_SASP_R-HSA-2559582.pdf",height=4)
plotEnrichment(CombineCSList[["R-HSA-2559582"]],ranks)
dev.off()

