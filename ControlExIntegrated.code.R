library(Seurat)
library(ggplot2)
library(clustree)
library(pheatmap)
library(presto)
library(dplyr)
library(msigdbr)
library(tibble)
library(fgsea)

#######################################################################################
##############    Merge snRNA profiles from control samples     #######################
#######################################################################################

MathyEx=readRDS("/Projects/deng/Alzheimer/syn18485175/cellCycle/MathyAllData/AllEx8Mathy.rds") #prefrontal cortex, For AD, age>75 years old  BA10
MathyExCtrl=subset(MathyEx,Statues %in% "Control") #24 sample, 16952 cells
MathyExCtrl$Gender=MathyExCtrl$Gender
MathyExCtrl$Age=MathyExCtrl$Age
MathyExCtrl$Source="Mathy"
MathyExCtrl$orig.ident=MathyExCtrl$ProjectID

YangEx=readRDS("/Projects/deng/Aging/Ex/Yang_GSE159812/YangIntegratedEx.rds") #Medial prefrontal cortex, For COVID, age>60 
YangExCtrl=subset(YangEx,Group %in% "Control") #8 samples, 8360 cells
tmp=ifelse(YangExCtrl$Gender == "F", "Female", "Male")
YangExCtrl$Gender=tmp
YangExCtrl$Age=YangExCtrl$Age
YangExCtrl$Source="Yang"
YangExCtrl$orig.ident=YangExCtrl$orig.ident

NagyEx=readRDS("/Projects/deng/Aging/Ex/Nagy_GSE144136/NagyEx.rds") #post-mortem dorsolateral prefrontal cortex (BA9), For MDD, age= 38.4(4)
NagyExCtrl=subset(NagyEx,Group %in% "Control") #17 samples, 18555 cells
NagyExCtrl$Gender="Male" #all male
NagyExCtrl$Age=38
NagyExCtrl$Source="Nagy"
NagyExCtrl$orig.ident=paste0("B",NagyExCtrl$sampleID,sep="")


LauEx=readRDS("/Projects/deng/Alzheimer/syn18485175/cellCycle/LauAllData/LauEx.rds")#post-mortem dorsolateral prefrontal cortex (BA9), For AD, age= 60-95
LauExCtrl=subset(LauEx,CONDITION %in% "Control")#12 samples, 30215 cells
LauExCtrl$Gender=LauExCtrl$SEX
LauExCtrl$Age=LauExCtrl$AGE
LauExCtrl$Source="Lau"
LauExCtrl$orig.ident=LauExCtrl$orig.ident


ExCtrl=merge(MathyExCtrl,y=c(YangExCtrl,NagyExCtrl,LauExCtrl))
ExCtrl=ExCtrl
for(i in colnames(ExCtrl@meta.data)) {
  ExCtrl[[i]] <- NULL
}
ExCtrl$orig.ident=ExCtrl$orig.ident
ExCtrl$nCount_RNA=ExCtrl$nCount_RNA
ExCtrl$nFeature_RNA=ExCtrl$nFeature_RNA
ExCtrl$Gender=ExCtrl$Gender
ExCtrl$Age=ExCtrl$Age
ExCtrl$Source=ExCtrl$Source

setwd("/Projects/deng/Aging/Ex/AllControlEx")
ExCtrl=ExCtrl
ExCtrl <- NormalizeData(ExCtrl, normalization.method = "LogNormalize", scale.factor = 10000)
ExCtrl <- FindVariableFeatures(ExCtrl, selection.method = "vst", nfeatures = 2000)
ExCtrl <- ScaleData(ExCtrl)
ExCtrl <- RunPCA(ExCtrl, features = VariableFeatures(object = ExCtrl),npcs = 100)
pdf("ExCtrlElbowPlot.pdf",width=5,height=5)
ElbowPlot(ExCtrl,ndims = 100)
dev.off()
# t-SNE and Clustering
ExCtrl<- RunUMAP(ExCtrl, reduction = "pca", dims = 1:50)
ExCtrl<- FindNeighbors(ExCtrl, reduction = "pca", dims = 1:50)
ExCtrl<- FindClusters(ExCtrl, resolution = 0.5)
pdf("ExCtrlLabel.pdf",width=9)
DimPlot(ExCtrl, label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("ExCtrlSource.pdf",width=8) 
DimPlot(ExCtrl, group.by="Source")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid")) #extremly source determined
dev.off()

saveRDS(ExCtrl,"ExCtrlMerge.rds")

#######################################################################################
##############   Integrated different datasets into a whole datasets      #############
#######################################################################################

setwd("/Projects/deng/Aging/Ex/AllControlEx")
library("Seurat")
ExCtrl=readRDS("ExCtrlMerge.rds")
ExCtrl[["percent.mt"]] <- PercentageFeatureSet(ExCtrl, pattern = "^MT-")
selected_mito <- WhichCells(ExCtrl, expression = percent.mt < 25)
ExCtrl <- subset(ExCtrl, cells = selected_mito)
ExCtrl <- ExCtrl[!grepl("^MT-", rownames(ExCtrl)), ]
memory.limit(size=400000)
ExCtrl=ExCtrl
ExCtrl.list <- SplitObject(ExCtrl, split.by = "Source")
ExCtrl.list <- lapply(X = ExCtrl.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = ExCtrl.list)
ExCtrl.list <- lapply(X = ExCtrl.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = ExCtrl.list, reference = c(1), reduction = "rpca",dims = 1:50)
ExCtrl.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
ExCtrl.integrated <- ScaleData(ExCtrl.integrated, verbose = FALSE)
ExCtrl.integrated <- RunPCA(ExCtrl.integrated, verbose = FALSE)
pdf("ExCtrl.integratedElbowPlot.pdf")
ElbowPlot(ExCtrl.integrated,ndims = 50)
dev.off()
# t-SNE and Clustering
ExCtrl.integrated <- RunTSNE(ExCtrl.integrated, dims = 1:50)
ExCtrl.integrated<- FindNeighbors(ExCtrl.integrated, reduction = "pca", dims = 1:30)

#----------- resolution selected based on cluster tree distribution -------------
setwd("/Projects/deng/Aging/Ex/AllControlEx")
ExCtrl.integrated=readRDS("/Projects/deng/Aging/Ex/AllControlEx/ExCtrl.integrated.rds")
ExTmp <- FindClusters(ExCtrl.integrated, resolution = c(seq(0.1,2.0,0.1)))
clus.tree.out <- clustree(ExTmp)
pdf(file = "ExCtrl.integrated.Resolution.Tree.pdf", width = 15, height = 20)
print(clus.tree.out)
dev.off()

ExCtrl.integrated<- FindClusters(ExCtrl.integrated, resolution = 0.2)
saveRDS(ExCtrl.integrated,file="/Projects/deng/Aging/Ex/AllControlEx/ExCtrl.integrated.rds")

tiff("ExCtrl.integratedLabelTrue.tiff",width=350,height=270)
DimPlot(ExCtrl.integrated, raster=FALSE,reduction = "tsne",label=T,label.size=6)&theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()


#######################################################################################
##############   phenotype visualization      #############
#######################################################################################

#Figure S2B
tiff("ExCtrl.integratedLabel.tiff",width=350,height=270)
DimPlot(ExCtrl.integrated, raster=FALSE,reduction = "tsne")&theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
#Figure S2C
tiff("ExCtrl.integratedSource.tiff",width=500,height=540)
DimPlot(ExCtrl.integrated, raster=FALSE,reduction = "tsne",split.by="Source",ncol=2)&NoLegend()&theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

#Figure S2D
t=as.matrix(table(ExCtrl.integrated$Source))
tmp=paste0(ExCtrl.integrated$seurat_clusters,"_",ExCtrl.integrated$Source,sep="")
tmp=as.matrix(table(tmp))
TmpInfo=as.matrix(do.call(rbind, strsplit(as.character(rownames(tmp)),'_')))
result=data.frame("Cluster"=TmpInfo[,1],"Source"=TmpInfo[,2],"Count"=tmp[,1])
ClusterList=factor(result$Cluster,levels=c(0:14))
g=ggplot(result, aes(ClusterList, Count, fill=Source)) +
  geom_bar(stat="identity",position="fill") +
  scale_y_continuous(expand=c(0,0))+
  scale_fill_brewer(palette="Set3")+
  theme(axis.text.x = element_text(angle = 0))+
  theme(axis.title.x = element_blank())
pdf("SourceDistribution.pdf",height=2.5,width=8)
print(g)
dev.off()


#Figure S2E
tiff("ExCtrl.integratedGender.tiff",width=500,height=270)
DimPlot(ExCtrl.integrated, raster=FALSE,reduction = "tsne",split.by="Gender")&NoLegend()&theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()


#Figure S2F
CCGene=read.table("/Projects/deng/Aging/Ex/cellCycleGene/CCGeneCombine.txt",header=T,sep="\t") #711
DefaultAssay(ExCtrl.integrated)="RNA"
genes=intersect(rownames(ExCtrl.integrated),CCGene$NAME)
length(genes) #364
CCGene=CCGene[CCGene$NAME%in%genes,]
CCGenelist=split(CCGene$NAME,CCGene$PHASE)
ExCtrl.integrated <- AddModuleScore(
  object = ExCtrl.integrated,
  features = CCGenelist,
  name = names(CCGenelist),
  ctrl = 100,
)
pdf("ExCtrl.integrated.cellCycleScore.pdf",width=6)
VlnPlot(ExCtrl.integrated,c("G1.S1","S.phase5","G22","G2.M3","M.G14"),ncol=1,pt.size=0)&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()



#######################################################################################
##############         Expression of DNA damage associated genes          #############
#######################################################################################
#Figure S5
ExTmpExpr <- AverageExpression(ExCtrl.integrated)[["RNA"]]
DNADamage=read.table("/Projects/deng/Aging/Ex/AllEx/DNArepair/DNArepairGeneList.txt",header=T,row.names=1,sep="\t")
DNADamage$GeneSymobl=rownames(DNADamage)

tmpCluster=DNADamage[DNADamage$Group %in% c("Nucleotide excision repair"),]
tmpClusterExpr=ExTmpExpr[intersect(rownames(ExTmpExpr),tmpCluster$GeneSymobl),]
pdf("Pheatmap/NER.pdf",width=6,height=4)
pheatmap(t(tmpClusterExpr),clustering_method="ward.D2",scale="column")
dev.off()


#######################################################################################
##############      Expression of cellular senescence associated genes    #############
#######################################################################################
#Figure 2D
coreCellularGene=read.table("55markersBySegura.txt",header=T)
DefaultAssay(ExCtrl)="RNA"
ExExpr <- AverageExpression(ExCtrl)[["RNA"]]
hclust=hclust(as.dist(1-cor(ExExpr)),method="ward.D2")
order=hclust$labels[hclust$order]
t=DotPlot(ExCtrl,features=coreCellularGene[,1])
result=merge(t$data,coreCellularGene,by.x="features.plot",by.y="GeneNames")
C3=result[result$id==0,]
C3=C3[order(C3$Group,C3$pct.exp),]
geneOrder=C3$features.plot
g=ggplot(result, aes(factor(features.plot,levels=rev(unique(geneOrder))),factor(id,levels=order), size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red")+
labs(color="Average",size="Percent",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0),face="bold"),axis.text.x = element_text(size=rel(1.0),face="bold",angle=90),axis.text.y = element_text(size=rel(1.0),face="bold")) 
pdf("SignatureCellularGene.pdf",width=10,height=4)
print(g)
dev.off()


#Figure 2E
ExSmall <- ExCtrl.integrated[, sample(colnames(ExCtrl.integrated), size = 10000, replace=F)]
DefaultAssay(ExSmall)="RNA"
t=DotPlot(ExSmall,features=rownames(ExSmall))
ExExpr8DotPlot=t$data

# Three list of cellular senescence gene list from BinZhang et,al. Reactome database, and KEGG
ConservedSenescenceGene=read.table("/Projects/deng/Aging/Ex/AllEx/Senescence/Senescence8BinZhang.txt")
genes=intersect(rownames(ExSmall),ConservedSenescenceGene[,1])
fgsea_sets=list("ConservedSeneGene"=genes)

m_df<- msigdbr(species = "Homo sapiens", category = "C2",subcategory="CP:REACTOME")  #R-HSA-2559583
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_sets=list("ReacomeSeneGene"=fgsea_sets$REACTOME_CELLULAR_SENESCENCE)

geneList=read.table("/Projects/deng/Aging/Ex/AllEx/Senescence/hsa04218_gene.txt",sep="\t") #hsa04218
fgsea_sets=list("KEGGSeneGene"=geneList[,2])


cluster="3"
clusterCell<- ExExpr8DotPlot %>% dplyr::filter(id == cluster) %>% arrange(desc(avg.exp.scaled)) %>% dplyr::select(features.plot, avg.exp.scaled)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0)
ranks=na.omit(ranks)

#Figure 2E
#GSEA plot for the cellular senescence assocaited genes
pdf("ConservedSeneGene.pdf",height=4)
plotEnrichment(fgsea_sets[["ConservedSeneGene"]],ranks)
dev.off()
pdf("SenescenceFromReactome.pdf",height=4)
plotEnrichment(fgsea_sets[["ReacomeSeneGene"]],ranks)
dev.off()
pdf("SenescenceFromKEGG.pdf",height=4)
plotEnrichment(fgsea_sets[["KEGGSeneGene"]],ranks)
dev.off()


# Output the  enrichment result of three cellular senescence gene lists from BinZhang et,al. Reactome database, and KEGG
fwrite(fgseaRes, file="ConservedSeneGene.txt", sep="\t", sep2=c("", " ", ""))
fwrite(fgseaRes, file="ReactomeSeneGene.txt", sep="\t", sep2=c("", " ", ""))
fwrite(fgseaRes, file="KEGGSeneGene.txt", sep="\t", sep2=c("", " ", ""))


# Figure S6A-C
# Exrepssion of the leadingEdge genes for the the enrichment result
GSEAResult=read.table("ConservedSeneGene.txt",header=T,sep="\t")
GSEAResult=read.table("ReactomeSeneGene.txt",header=T,sep="\t")
GSEAResult=read.table("KEGGSeneGene.txt",header=T,sep="\t")
leadingEdge=unlist(strsplit(GSEAResult$leadingEdge," "))
t=DotPlot(ExCtrl,features=leadingEdge)
result=t$data
C3=result[result$id==3,]
C3=C3[order(C3$pct.exp,decreasing=T),]
geneOrder=C3$features.plot
g=ggplot(result, aes(factor(features.plot,levels=unique(geneOrder)),factor(id,levels=order), size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red")+
labs(color="Average",size="Percent",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0),face="bold"),axis.text.x = element_text(size=rel(1.0),face="bold",angle=90),axis.text.y = element_text(size=rel(1.0),face="bold")) 
pdf("LeadingEdgeGeneInKEGGSeneGene.pdf",width=9,height=4)
print(g)
dev.off()


# Figure S6D
# Exrepssion of the SASP assocaited genes
SASPmarker=read.table("/Projects/deng/Aging/Ex/AllEx/Senescence/SASPmarker.txt",header=T,sep="\t")
t=DotPlot(ExCtrl,features=SASPmarker$GeneSymobl)
result=merge(t$data,SASPmarker,by.x="features.plot",by.y="GeneSymobl")
C5=result[result$id==0,]
C5=C5[order(C5$Group,C5$pct.exp),]
geneOrder=C5$features.plot
#result=result[result$id %in% c(2:6),]
g=ggplot(result, aes(factor(features.plot,levels=unique(rev(geneOrder))), factor(id,levels=order), size= pct.exp,color=avg.exp.scaled)) +geom_point(alpha = 0.9)+
scale_colour_viridis_c()+
labs(color="Average",size="Percent",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0),face="bold"),axis.text.x = element_text(face="bold",angle=90),axis.text.y = element_text(size=rel(1.0),face="bold")) 
pdf("SASPmarkerGene.pdf",width=15,height=4)
print(g)
dev.off()


