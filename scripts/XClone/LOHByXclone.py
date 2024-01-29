conda create -n xclone python=3.7
conda activate xclone
pip install git+https://github.com/single-cell-genetics/XClone
pip install requests
pip install importlib-metadata

### validate the xclone on TNBC
import xclone
import anndata as an
import pandas as pd
import numpy as np
import scipy
print("scipy", scipy.__version__)
#scipy 1.7.3

xclone.pp.efficiency_preview()
#multiprocessing cpu total count in your device 48

dataset_name = "TNBC1_scRNA"
## output results dir
outdir = "/data2/deng/Aging/Ex/Lau_GSE157827/xclone/tutorials/"
RDR_adata = xclone.data.tnbc1_rdr()
BAF_adata = xclone.data.tnbc1_baf()

xconfig = xclone.XCloneConfig(dataset_name = dataset_name, module = "RDR")
xconfig.set_figure_params(xclone= True, fontsize = 18)
xconfig.outdir = outdir
xconfig.cell_anno_key = "cluster.pred"
xconfig.ref_celltype = "N"
xconfig.marker_group_anno_key = "cluster.pred"
xconfig.xclone_plot= True
xconfig.plot_cell_anno_key = "cluster"
xconfig.display()
RDR_Xdata = xclone.model.run_RDR(RDR_adata,config_file = xconfig)

RDR_Xdata

import anndata
filename = outdir+"data/RDR_adata_KNN_HMM_post.h5ad"
adata = anndata.read_h5ad(filename)
#Access the posterior_mtx layer:
cnv_prob = adata.layers["posterior_mtx"]
copy_states = cnv_prob.argmax(axis=-1)

xconfig = xclone.XCloneConfig(dataset_name = dataset_name, module = "BAF")
xconfig.set_figure_params(xclone= True, fontsize = 18)
xconfig.outdir = outdir
xconfig.cell_anno_key = "cluster.pred"
xconfig.ref_celltype = "N"
xconfig.xclone_plot= True
xconfig.plot_cell_anno_key = "cluster"
xconfig.display()

BAF_merge_Xdata = xclone.model.run_BAF(BAF_adata,config_file = xconfig)



### using XClone detect LOH of LS neurons

####prepair the data used by XClone using xcltk
pip install -U xcltk
library(Seurat)
setwd("/Projects/deng/Aging/Ex/AllEx")
ExAll.integrated=readRDS("ExAll.integratedTSNE.rds")
ExAll.integrated=subset(ExAll.integrated,idents=c(0:15)) #remove cluster 16 by their smaller cell number
ExAll.integrated[["Layer"]]=NA
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(0,2,8),]="L2/3"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(1,6,7,9,13),]="L4/5"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(4,10),]="L4/6"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(15),]="L5"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(11,14),]="L5/6"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(12),]="L6"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(3),]="ES"
ExAll.integrated[["Layer"]][ExAll.integrated$seurat_clusters%in%c(5),]="LS"
ExAll.integrated$Layer=factor(ExAll.integrated$Layer,levels=c("ES","LS","L2/3","L4/5","L4/6","L5","L5/6","L6"))

### output the normal neurons for target samples
for (sample in c("GSM4775561","GSM4775562","GSM4775578","GSM4775572","GSM4775580","GSM4775581")){

  targetSample=subset(ExAll.integrated,orig.ident%in%sample)
  cellAnno=targetSample@meta.data[,c("seurat_clusters","Statues","Layer")]
  cellAnno$cell=as.matrix(do.call(rbind, strsplit(as.character(rownames(cellAnno)),'_')))[,2]
  cellAnno$clusterPred=ifelse(cellAnno$seurat_clusters=="5","LS","Others")

  targetPath=paste0("/data2/deng/Aging/Ex/Lau_GSE157827/xcltk/",sample,sep="")
    if(!dir.exists(targetPath)){
	 dir.create(targetPath)
  }
  write.table(cellAnno,file=paste0("/data2/deng/Aging/Ex/Lau_GSE157827/xcltk/",sample,"/ExClusterByIntegrated.tsv",sep=""),sep="\t",quote=F,row.names=F)

  normalCells.seurat=subset(targetSample,seurat_clusters%in%c(0:2,4,6:15))
  normalCells=as.matrix(do.call(rbind, strsplit(as.character(colnames(normalCells.seurat)),'_')))[,2]
  write.table(normalCells,file=paste0("/data2/deng/Aging/Ex/Lau_GSE157827/xcltk/",sample,"/NormalExIDByIntegrated.tsv",sep=""),quote=F,row.names=F,col.names=F)	
}

#################################################################################################
#################################################################################################
#################################################################################################
###need cellsnp-lite (which depends on several external libraries such as htslib. )
conda create -n cellsnp
conda activate cellsnp
conda install -c bioconda cellsnp-lite
conda install python=3.7

#export PATH=$PATH:/Softwares/deng/samtools/samtools-1.9/bin;

wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x /Softwares/deng/liftOver/liftOver
wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2
./configure --prefix=/Softwares/deng/bcftools/bcftools-1.19
make
make install

wget https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2
cd samtools-1.19
./configure --prefix=/Softwares/deng/samtools/samtools-1.19
pip install -U xcltk

#################################################################################################
conda activate cellsnp
export PATH=$PATH:/Softwares/deng/liftOver;
export PATH=$PATH:/Softwares/deng/bcftools/bcftools-1.19
export PATH=$PATH:/Softwares/deng/samtools/samtools-1.19

sample="GSM4775572"
mkdir /data2/deng/Aging/Ex/Lau_GSE157827/xcltk/${sample}/RDRResult
cd /data2/deng/Aging/Ex/Lau_GSE157827/xcltk/${sample}/RDRResult
#################################################################################################
#RDR part: count sparse matrix obtain
xcltk basefc               \
    -s  /data2/deng/Aging/Ex/Lau_GSE157827/rawData/${sample}/${sample}/outs/possorted_genome_bam.bam          \
    -b  /data2/deng/Aging/Ex/Lau_GSE157827/rawData/${sample}/${sample}/outs/filtered_feature_bc_matrix/barcodes.tsv      \
    -r  /Softwares/deng/xcltk/xcltk-master/preprocess/data/annotate_genes_hg38_update_20230126.txt       \
    -T  tsv                 \
    -O  /data2/deng/Aging/Ex/Lau_GSE157827/xcltk/${sample}/RDRResult           \
    -p  10                  \
    --cellTAG  CB           \
    --UMItag   UB           \
    --minLEN   30           \
    --minMAPQ  20           \
    --maxFLAG  4096


#################################################################################################
#BAF part
cd /data2/deng/Aging/Ex/Lau_GSE157827/xcltk/${sample}
#https://github.com/hxj5/xcltk/blob/master/preprocess/README.md

#Pre-Phasing
#-b only for normal cells
/Softwares/deng/xcltk/xcltk-master/preprocess/baf_pre_phase.sh            \
    -N  ${sample}         \
    -s  /data2/deng/Aging/Ex/Lau_GSE157827/rawData/${sample}/${sample}/outs/possorted_genome_bam.bam            \
    -b  /data2/deng/Aging/Ex/Lau_GSE157827/xcltk/${sample}/NormalExIDByIntegrated.tsv       \
    -F  /Softwares/deng/xcltk/xcltk-master/preprocess/data/human_g1k_v37.fasta   \
    -O  /data2/deng/Aging/Ex/Lau_GSE157827/xcltk/${sample}/PrePhasing             \
    -g  38                    \
    -C  CB                    \
    -u  UB                    \
    -p  10
# Filter cells. Only keep cells with valid barcodes

PrePhasing/result/GSM4775572.hg38.raw.het.qc.vcf.gz
#Keep heterozygous SNPs only, done
#Filter SNPs with low GQ score. 
#function asort never defined and then the following step need to be manually modified

cd /data2/deng/Aging/Ex/Lau_GSE157827/xcltk/${sample}/PrePhasing/result
#Convert hg38 to hg19
python /Softwares/deng/xcltk/xcltk-master/preprocess/liftOver_vcf.py -c /Softwares/deng/xcltk/xcltk-master/preprocess/data/hg38ToHg19.over.chain.gz -i ${sample}.hg38.raw.het.qc.vcf.gz -o /data2/deng/Aging/Ex/Lau_GSE157827/xcltk/${sample}/PrePhasing/result/${sample}.hg38.raw.het.qc.gq.hg19.vcf.gz -P liftOver
#remove reference allele mismatch
xcltk fixref -i ${sample}.hg38.raw.het.qc.gq.hg19.vcf.gz -r /Softwares/deng/xcltk/xcltk-master/preprocess/data/human_g1k_v37.fasta -o ${sample}.hg38.raw.het.qc.gq.hg19.ref.vcf.gz
#sort for index
bcftools sort -Oz ${sample}.hg38.raw.het.qc.gq.hg19.ref.vcf.gz -o ${sample}.hg38.raw.het.qc.gq.hg19.sort.vcf.gz
#check the reference allele mismatch (may be used by imputation.sanger
bcftools norm --check-ref e -f /Softwares/deng/xcltk/xcltk-master/preprocess/data/human_g1k_v37.fasta ${sample}.hg38.raw.het.qc.gq.hg19.sort.vcf.gz -Ou -o /dev/null


#https://imputation.sanger.ac.uk/
ll /data2/deng/Aging/Ex/Lau_GSE157827/xcltk/${sample}/PrePhasing/result
#${sample}.hg38.raw.het.qc.gq.hg19.sort.vcf.gz was downloaded for sanger phasing
#Haplotype Reference Consortium (r1.1)
#phase with EAGLE2, no imputation

ll /data2/deng/Aging/Ex/Lau_GSE157827/xcltk/${sample}/SangerPhasing/
cd /data2/deng/Aging/Ex/Lau_GSE157827/rawData/${sample}/${sample}/outs/filtered_feature_bc_matrix
gzip -d barcodes.tsv.gz
#Post-Phasing
#GSM4775572.hg38.raw.het.qc.gq.hg19.SangerPhasing.vcf.gz  was generated by sanger web server. not the -b should contains the whole cells (consistent to the RDR part)
cd /data2/deng/Aging/Ex/Lau_GSE157827/xcltk/${sample}/PostPhasing
/Softwares/deng/xcltk/xcltk-master/preprocess/baf_post_phase.sh        \
    -N  ${sample}      \
    -s  /data2/deng/Aging/Ex/Lau_GSE157827/rawData/${sample}/${sample}/outs/possorted_genome_bam.bam         \
    -b  /data2/deng/Aging/Ex/Lau_GSE157827/rawData/${sample}/${sample}/outs/filtered_feature_bc_matrix/barcodes.tsv     \
    -v  /data2/deng/Aging/Ex/Lau_GSE157827/xcltk/${sample}/SangerPhasing/${sample}.hg38.raw.het.qc.gq.hg19.SangerPhasing.vcf.gz      \
    -f  /Softwares/deng/xcltk/xcltk-master/preprocess/data/human_g1k_v37.fasta       \
    -O  /data2/deng/Aging/Ex/Lau_GSE157827/xcltk/${sample}/PostPhasing       \
    -g  38                 \
    -C  CB                 \
    -u  UB                 \
    -p  10



#output_region_file = /data2/deng/Aging/Ex/Lau_GSE157827/xcltk/GSM4775572/PostPhasing/result/baf/xcltk.region.tsv
#output_sample_file = /data2/deng/Aging/Ex/Lau_GSE157827/xcltk/GSM4775572/PostPhasing/result/baf/xcltk.samples.tsv
#output_ad_file = /data2/deng/Aging/Ex/Lau_GSE157827/xcltk/GSM4775572/PostPhasing/result/baf/xcltk.AD.mtx
#output_dp_file = /data2/deng/Aging/Ex/Lau_GSE157827/xcltk/GSM4775572/PostPhasing/result/baf/xcltk.DP.mtx
#output_oth_file = /data2/deng/Aging/Ex/Lau_GSE157827/xcltk/GSM4775572/PostPhasing/result/baf/xcltk.OTH.mtx

#################################################################################################
#################################################################################################
#################################################################################################

#xclone preprocessing
#build anndata from the files: 1, need to install xclone from GitHub repository (for development version), hg38_genes files were provided in this version.  
conda activate xclone
pip install git+https://github.com/single-cell-genetics/XClone
pip install importlib-metadata

python
import xclone

hg38_genes = xclone.pp.load_anno(genome_mode = "hg38_genes")
hg38_blocks = xclone.pp.load_anno(genome_mode = "hg38_blocks")
anno_file="/data2/deng/Aging/Ex/Lau_GSE157827/xcltk/GSM4775572/ExClusterByIntegrated.tsv"

data_dir = "/data2/deng/Aging/Ex/Lau_GSE157827/xcltk/GSM4775572/RDRResult/"
RDR_file = data_dir + "matrix.mtx"
mtx_barcodes_file = data_dir + "barcodes.tsv" # cell barcodes
regions_anno_file = data_dir + "features.tsv" # feature annnotation
RDR_adata=xclone.pp.xclonedata(RDR_file,'RDR',mtx_barcodes_file,genome_mode = "hg38_genes")
RDR_adata = xclone.pp.extra_anno(RDR_adata, anno_file, barcodes_key = "cell",cell_anno_key = ["seurat_clusters", "Statues","clusterPred","Layer"], sep = "\t")

data_dir = "/data2/deng/Aging/Ex/Lau_GSE157827/xcltk/GSM4775572/PostPhasing/baf/"
AD_file = data_dir + "xcltk.AD.mtx"
DP_file = data_dir + "xcltk.DP.mtx"
mtx_barcodes_file = data_dir + "xcltk.samples.tsv" # cell barcodes
# use default gene annotation
BAF_adata = xclone.pp.xclonedata([AD_file, DP_file], 'BAF',mtx_barcodes_file,genome_mode = "hg38_genes")
BAF_adata = xclone.pp.extra_anno(BAF_adata, anno_file, barcodes_key = "cell", cell_anno_key = ["seurat_clusters", "Statues","clusterPred","Layer"], sep = "\t")


import anndata as an
import pandas as pd
import numpy as np
import scipy
print("scipy", scipy.__version__)
#scipy 1.7.3

xclone.pp.efficiency_preview()
#multiprocessing cpu total count in your device 48

dataset_name = "GSM4775572"
## output results dir
outdir = "/data2/deng/Aging/Ex/Lau_GSE157827/xclone/GSM4775572/"

xconfig = xclone.XCloneConfig(dataset_name = dataset_name, module = "RDR")
xconfig.set_figure_params(xclone= True, fontsize = 15)
xconfig.outdir = outdir
xconfig.cell_anno_key = "clusterPred"
xconfig.ref_celltype = "Others"
xconfig.marker_group_anno_key = "clusterPred"
xconfig.xclone_plot= True
xconfig.plot_cell_anno_key = "Layer"
xconfig.display()
RDR_adata.var
RDR_adata.obs
RDR_Xdata = xclone.model.run_RDR(RDR_adata,config_file = xconfig)
#Keep valid cells: Filter out 4408 cells / 6376 total cells, remain 1968 valid cells with annotation
#RDR CNV states chrs guiding(copy loss, copy neutral, copy gain): ['Yp', '2p', '19p']
#RDR CNV states ratio guiding(copy loss, copy neutral, copy gain): [0.54241007 1.23754538 1.60363657]



import anndata
filename = outdir+"data/RDR_adata_KNN_HMM_post.h5ad"
adata = anndata.read_h5ad(filename)
#Access the posterior_mtx layer:
cnv_prob = adata.layers["posterior_mtx"]
copy_states = cnv_prob.argmax(axis=-1)

xconfig = xclone.XCloneConfig(dataset_name = dataset_name, module = "BAF")
xconfig.set_figure_params(xclone= True, fontsize = 15)
xconfig.outdir = outdir
xconfig.cell_anno_key = "clusterPred"
xconfig.ref_celltype = "Others"
xconfig.xclone_plot= True
xconfig.plot_cell_anno_key = "Layer"
xconfig.display()

BAF_merge_Xdata = xclone.model.run_BAF(BAF_adata,config_file = xconfig)
#Keep valid cells: Filter out 8996 cells / 14842 total cells, remain 5846 valid cells with annotation

BAF_merge_Xdata.var
RDR_Xdata.var

flag = ~(RDR_Xdata.var["chr"] == "Y")
RDR_Xdata = RDR_Xdata[:, flag]

flag = ~(BAF_merge_Xdata.var["chr"] == "Y")
BAF_merge_Xdata = BAF_merge_Xdata[:, flag]


xconfig = xclone.XCloneConfig(dataset_name = dataset_name, module = "Combine")
xconfig.set_figure_params(xclone= True, fontsize = 15)
xconfig.outdir = outdir
xconfig.cell_anno_key = "clusterPred"
xconfig.ref_celltype = "Others"
xconfig.xclone_plot= True
xconfig.plot_cell_anno_key = "Layer"
xconfig.BAF_denoise = True

xconfig.display()

xclone.model.run_combine(RDR_Xdata,
                         BAF_merge_Xdata,
                         verbose = True,
                         run_verbose = True,
                         config_file = xconfig)

