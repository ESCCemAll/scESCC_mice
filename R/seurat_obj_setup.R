####################################################################  
#
# Seurat object setup
#
####################################################################  


library(Seurat)


###### Create Seurat object ######
setwd(workpath)
samplename  # chose the samplename basing on celltype. For example, Epi for  epithelial cells.
seuratobj.data <- readRDS("rawdata.rds")  # load the counts data
minimun_cells = ncol(seuratobj.data)/1000
seuratobj <- CreateSeuratObject(seuratobj.data, min.cells = minimun_cells, 
    min.genes = 500, project = "10X_seuratobj")


###### QC, Normalization and scaling ######
mito.genes <- grep(pattern = "^mt-", x = rownames(x = seuratobj@data), 
    value = TRUE)
percent.mito <- Matrix::colSums(seuratobj@raw.data[mito.genes, ])/Matrix::colSums(seuratobj@raw.data)
seuratobj <- AddMetaData(object = seuratobj, metadata = percent.mito, 
    col.name = "percent.mito")
seuratobj <- FilterCells(object = seuratobj, subset.names = c("nGene", 
    "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(Inf, 
    0.2))
seuratobj <- NormalizeData(object = seuratobj, normalization.method = "LogNormalize", 
    scale.factor = 10000)
seuratobj <- FindVariableGenes(object = seuratobj, mean.function = ExpMean, 
    dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, 
    y.cutoff = 0.5)
seuratobj <- ScaleData(object = seuratobj, vars.to.regress = c("nUMI", 
    "percent.mito"))
seuratobj <- RunPCA(object = seuratobj, pc.genes = seuratobj@var.genes, 
    do.print = TRUE, pcs.print = 1:5, genes.print = 5)


###### PCA ######
seuratobj <- ProjectPCA(object = seuratobj, do.print = FALSE)
JackStrawPlot(object = seuratobj, PCs = 1:12)
PCElbowPlot(object = seuratobj)


###### Clustering ######
pcs  # select the number of PCs where there is a clear elbow in the output of PCElbowPlot()
res  # set resolution value
seuratobj <- FindClusters(object = seuratobj, reduction.type = "pca", 
    dims.use = 1:pcs, resolution = res, print.output = 0, save.SNN = TRUE)
seuratobj <- RunTSNE(object = seuratobj, dims.use = 1:pcs, do.fast = TRUE)


###### Visualization ######
TSNEPlot(seuratobj, do.label = T)