library(Seurat)
library(RColorBrewer)
setwd("/boot3/bixm/08_HNLongevity/pipline")
#using your own pipine to generated this neuron dataset
ExIntegrated=readRDS("ExIntegrated.Downsampling.rds")
#set the default assay as RNA if you using the integrated method
DefaultAssay(ExIntegrated)="RNA"
ExIntegrated <- ExIntegrated[!grepl("^MT-", rownames(ExIntegrated)), ] #remove the genes coded by mitocondrial
#obtain the count matrix (as infercnv was robust to the cell numbers, and then you can down sampling the cells to obtain a quicker result)
counts_matrix = as.matrix(ExIntegrated@assays$RNA@counts)
#set the cluster as phenotype
pData=ExIntegrated@meta.data$seurat_clusters
names(pData)=rownames(ExIntegrated@meta.data)
pData=data.frame(pData)
#prepair this file 
annotation=read.table("gencode.v38.annotation.geneType.txt",header=TRUE,row.names=1)
anno=data.frame("Chr"=annotation$Chr,"Start"=annotation$Start,"End"=annotation$End)
rownames(anno)=annotation$Symbol
genes=intersect(rownames(counts_matrix),rownames(anno))
length(genes)#31566
anno=anno[genes,]
anno$Chr=factor(anno$Chr,levels=paste0("chr",c(1:22,"X","Y"),sep=""))
anno=anno[order(anno$Chr,anno$Star,anno$End),]
genes=intersect(rownames(anno),rownames(counts_matrix))
length(genes)
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
                               out_dir="/data2/deng/Aging/Ex/AllEx/InferCNV/ExAllDownsampling"
                               )
