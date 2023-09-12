R (v4.2.1)
#https://www.nature.com/articles/s41593-020-00787-0
#https://research.libd.org/spatialLIBD/articles/spatialLIBD.html
library("spatialLIBD") 
packageVersion("spatialLIBD")
#‘1.13.4’

## Connect to ExperimentHub
###############this code buck was used to delete the old cache (As I downloaded them in R4.0.5)#########################
#https://www.bioconductor.org/packages/devel/bioc/vignettes/BiocFileCache/inst/doc/BiocFileCache.html#default-caching-location-update
#7.1Option 1: Moving Files
#########################################################################################################################
setwd("D:/Aging/Glia/spatialLIBD")

## Connect to ExperimentHub
ehub <- ExperimentHub::ExperimentHub()
#C:\Users\dengw\AppData\Local/R/cache/R/BiocFileCache

spe <- fetch_data(type = "spe", eh = ehub)
# download the full real data 
#C:\Users\dengw\AppData\Local\Temp\Rtmp6hmoyY/BiocFileCache/902cda43110_Human_DLPFC_Visium_processedData_sce_scran_spatialLIBD.Rdata%3Fdl%3D1

sce_layer <- fetch_data(type = "sce_layer", eh = ehub)
#C:\Users\dengw\AppData\Local\Temp\Rtmp6hmoyY/BiocFileCache/902c1de5178f_Human_DLPFC_Visium_processedData_sce_scran_sce_layer_spatialLIBD.Rdata%3Fdl%3D1

modeling_results <- fetch_data("modeling_results", eh = ehub)


## View our LIBD layers for one sample
pdf("151507_Layers.pdf",width=6,height=7)
vis_clus(
    spe = spe,
    clustervar = "layer_guess_reordered",
    sampleid = "151507",
    colors = libd_layer_colors,
    ... = " LIBD Layers"
)
dev.off()


#"GLRA3","CUX2","RORB","PCP4","SEMA3E","AQP4","NTNG2"

pdf("AQP4_New.pdf",width=4,height=5)
vis_gene(
    spe = spe,
    sampleid = "151507",
    viridis=FALSE,
    geneid = grep("^AQP4; ",rowData(spe)$gene_search,value=T),
    spatial = TRUE
)
dev.off()

vis_gene(
    spe = spe,
    sampleid = "151507",
    viridis=TRUE,
    geneid = grep("^AQP4; ",rowData(spe)$gene_search,value=T),
    spatial = TRUE
)

sig_genes <-
     sig_genes_extract_all(
            n = nrow(sce_layer),
            modeling_results = modeling_results,
            sce_layer = sce_layer
        )

pdf("AQP4_Boxplot.pdf",width=7,height=6)
layer_boxplot(
    i = which(sig_genes$gene == "AQP4")[1],
    sig_genes = sig_genes,
    sce_layer = sce_layer,
    col_low_box = "palegreen3",
    col_low_point = "springgreen2",
    col_high_box = "palegreen3",
    col_high_point = "springgreen2",
    short_title=1
)
dev.off()

