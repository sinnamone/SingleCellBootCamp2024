library(Seurat)


#to use Seurat by new v5 assay

options(Seurat.object.assay.version = 'v5')


#####load and create Seurat object from 10x cellranger count output directory

  sample_1 <- Read10X("Patient_01/raw_feature_bc_matrix/")
sample_2 <- Read10X("Patient_02/raw_feature_bc_matrix/")
#sample_3 <- Read10X("corso_sc/Patient_03/raw_feature_bc_matrix/")

#if you have a metadata table by cell barcode as row and feature as column can provide by meta.data option
sample_1 <- CreateSeuratObject(counts = sample_1)
sample_2 <- CreateSeuratObject(counts = sample_2)
#sample_3 <- CreateSeuratObject(counts = sample_3)


#####define Sample variable for each sample
sample_1$Sample <- "Patient_1"
sample_2$Sample <- "Patient_2"

#sample_3$Sample <- "Patient 3"


#####in droplet based single-cell analysis we obtain many empty droplet. 
####To remove them we build a curve that shows the log-count of umi associated with each ranked barcode 

seuratList<-list(sample_1,sample_2)

seuratList<- lapply(seuratList,function(x){CalculateBarcodeInflections(x,threshold.high = 10000,threshold.low = 1000)})

pdf("curve_barcode.pdf")
BarcodeInflectionsPlot(seuratList[[1]]) ##dotted line is the inflection point computed in sample_1. 
dev.off()
###We should zoom to appreciate the inversion of the curve

seuratList<- lapply(seuratList,SubsetByBarcodeInflections)


##compute Mt and ribo percent

seuratList<- lapply(seuratList, function(x){PercentageFeatureSet(x, pattern = "^RP[SL]", col.name = "percent.ribo")})
seuratList<- lapply(seuratList,function(x){PercentageFeatureSet(x, pattern = "^MT-", col.name = "percent.mito")})

#####now let's display some of the most used metrics for qc
pdf("ViolinPlotRibo.pdf")
VlnPlot(seuratList[[1]], features = c("percent.ribo","percent.mito","nCount_RNA","nFeature_RNA"))
dev.off()

pdf("ViolinPlotRibo_ptz2.pdf")
VlnPlot(seuratList[[2]], features = c("percent.ribo","percent.mito","nCount_RNA","nFeature_RNA"))
dev.off()
#VlnPlot(seuratList[[3]], features = c("percent.ribo","percent.mito","nCount_RNA","nFeature_RNA"))


####the tresholds for cutting low-quality cells can be chosen arbitrarily or by trying to identify which cells are outliers for each 
####Feature. A widely used method is to remove all cells that are less than three times the standard deviation of the 
####cell distribution for a feature. Furthermore you can set high treshold to filtered out captured droplets.

###generally it is advisable to consider cells with a percentage of mitochondrial genes lower than 14 and at least 200 genes as good.
###The values we calculate here can be very stringent indeed and exclude rare cell populations that have a different distribution 
###due to their biology. Therefore you need to adapt this filter to each dataset

seuratList<- lapply(seuratList,function(x){subset(x, subset = percent.mito <= 14 & nFeature_RNA>= 200)})


#####now we can visualize the cleaned dataset
#####now let's display some of the most used metrics for qc
pdf("ViolinPlot_postQC_ptz1.pdf")
VlnPlot(seuratList[[1]], features = c("percent.ribo","percent.mito","nCount_RNA","nFeature_RNA"))
dev.off()
pdf("ViolinPlot_postQC_ptz2.pdf")
VlnPlot(seuratList[[2]], features = c("percent.ribo","percent.mito","nCount_RNA","nFeature_RNA"))
dev.off()

#VlnPlot(seuratList[[3]], features = c("percent.ribo","percent.mito","nCount_RNA","nFeature_RNA"))


###remove low expressed genes which could affect the normalization step. We filtered out all genes which are expressed
### less than 1% of all cells

### Version 5
for(i in 1:length(seuratList)){
  counts<-LayerData(seuratList[[i]],layer = "counts")
  cell_percent<- rowMeans(counts>0)*100
  expressed_genes<-names(cell_percent[cell_percent>1])
  seuratList[[i]]<-seuratList[[i]][expressed_genes,]
}

### Version 4
for(i in 1:length(seuratList)){
  counts<-GetAssayData(seuratList[[i]],slot = "counts")
  cell_percent<- rowMeans(counts>0)*100
  expressed_genes<-names(cell_percent[cell_percent>1])
  seuratList[[i]]<-seuratList[[i]][expressed_genes,]
}




#####merge all filtered sample into one Seurat object  

seuratObj <- merge(seuratList[[1]], c(seuratList[[2]]), add.cell.ids=c("sample_1","sample_2"))

##we can remove the 3 objects no longer needed from the environment

remove(sample_1,sample_2,seuratList)
gc()

#####let's look closely at some characteristics of the created object

##you can view the available assay of the object by

seuratObj@assays

#or by function Assays

Assays(seuratObj)

##for the moment only the "RNA" assay is present. You can view the assay in use by

seuratObj@active.assay

###each array is divided in layers which contained the raw count matrix of the samples
Layers(seuratObj)


##you can view the metadata by

seuratObj@meta.data[1:5,]

##you can view the number of total cell by
ncol(seuratObj)

##you can view the number of total genes by
nrow(seuratObj)

###we can check the number of cell per sample by
table(seuratObj$Sample)

###now we can proceed with the normalization of the counts. 
###Seurat by default normalizes by dividing all counts by the same factor, 
###10000 but this can be changed by indicating it inside the function, and running the natural logarithm. 
###Other normalizations are available such as SCT, relative count or centered log ratio.

seuratObj<-NormalizeData(seuratObj)

##now we can see that the data layers containing the normalized counts has been added
Layers(seuratObj)

seuratObj@assays$RNA@data[1:10,1:5]

##to save an object
saveRDS(seuratObj,"QC_seurat_object.rds")
#to reload the object

seuratObj<- readRDS("seurat_object.rds")

##to save the entire envinroment
save.image("seurat_analysis.Rda")

#to reload the env
load("seurat_analysis.Rda")

