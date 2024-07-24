setwd("/home/pkf/gendata2/bioinformatica/rcarriero/corso")

#Upload libraries

library(Seurat)
library(dplyr)
library(clustree)
library(ggraph)

#Upload Seurat object

tcells<- readRDS('Treg.rds')

#Find variable features

tcells <- FindVariableFeatures(tcells, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(tcells), 10)

# plot variable features with and without labels
pdf("plot_variable_features.pdf", width = 15, height = 10)
plot1 <- VariableFeaturePlot(tcells)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

#Scale data

all.genes <- rownames(tcells)
tcells <- ScaleData(tcells, features = all.genes)

#Run PCA

tcells <- RunPCA(tcells, features = VariableFeatures(object = tcells))

# Examine and visualize PCA results a few different ways
print(tcells[["pca"]], dims = 1:5, nfeatures = 5)


pdf("vizdim_pca_plot.pdf")
VizDimLoadings(tcells, dims = 1:2, reduction = "pca")
DimPlot(tcells, reduction = "pca")
DimHeatmap(tcells, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(tcells, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()


tcells <- JackStraw(tcells, num.replicate = 100)
tcells <- ScoreJackStraw(tcells, dims = 1:20)

pdf("JackStrawPlot.pdf")
JackStrawPlot(tcells, dims = 1:20)
ElbowPlot(tcells)
dev.off()

#Select significant components
tcells <- FindNeighbors(tcells, dims = 1:14)

cells_clust_tree <- FindClusters(tcells, resolution = seq(from=0, to=2, by=0.1), print.output = 0, save.SNN = T)

pdf("CL7_tregs_Resolution_tree_up_to_2.0_pca_1_17.pdf",  width=17, height=17)
clustree(cells_clust_tree)
dev.off()

tcells <- FindClusters(tcells, resolution = 0.4)
tcells <- RunUMAP(tcells, dims = 1:14)

pdf("umap.pdf")
DimPlot(tcells, reduction = "umap", label = T)
dev.off()

tcells.markers <- FindAllMarkers(tcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tcells.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)


tcells.markers %>% group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

pdf("heatmap.pdf")
DoHeatmap(tcells, features = top10$gene) + NoLegend()
dev.off()

write.table(tcells.markers, "tcells.markers.txt", sep="\t")

#Plot of marker genes

pdf("VlnPlot_FOXP3.pdf")
VlnPlot(tcells, features = 'FOXP3')
dev.off()

pdf("FeaturePlot_CXCL13.pdf")
FeaturePlot(tcells, features = 'CXCL13', label = T)
dev.off()


pdf("DotPlot_IL1R2.pdf")
DotPlot(tcells, features = 'IL1R2')
dev.off()

#Responder vs Non responder --> all cells
Idents(tcells)

levels(tcells@meta.data$`tcells@active.ident`)[1] <-"Non_responder"
levels(tcells@meta.data$`tcells@active.ident`)[2] <-"Responder"
table(tcells@meta.data$`tcells@active.ident`)

Idents(tcells)<- tcells@meta.data$`tcells@active.ident`
R_vs_NR<- FindMarkers(tcells, ident.1 = 'Responder', ident.2 = 'Non_responder')
head(R_vs_NR)
write.table(R_vs_NR, file='R_vs_NR.txt', sep="\t")

#Responder vs Non responder--> cluster 0
Idents(tcells)<- tcells$RNA_snn_res.0.4
tcells_c0<- subset(tcells, idents = '0')
Idents(tcells_c0)<-tcells_c0@meta.data$`tcells@active.ident`
Idents(tcells_c0)
c0_R_vs_NR<- FindMarkers(tcells_c0, ident.1 = 'Responder', ident.2 = 'Non_responder')
head(c0_R_vs_NR)
write.table(c0_R_vs_NR, file='c0_R_vs_NR.txt', sep="\t")
