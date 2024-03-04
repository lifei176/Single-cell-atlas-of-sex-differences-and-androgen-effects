################# 1. Prepare required softwares
#################
#################
#################
#################
#################
library(Seurat)
library(ggsci)
library(viridis)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(BiocParallel)
library(harmony) 
library(DoubletFinder)
library(DropletUtils) 
library(clustree)
library(future)
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', "#b20000",'#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')
mypal <- union(pal_npg("nrc", alpha = 0.7)(10),my36colors)
mypal
################# 2. Load expression matrix data and creat seurat object
#################
#################
#################
#################
#################

################# Adrenal_FD_1
counts_Adrenal_FD_1<- Read10X_h5("Adrenal_FD_1.h5")
seurat_Adrenal_FD_1 <- CreateSeuratObject(counts = counts_Adrenal_FD_1,project = "Adrenal_FD_1")
seurat_Adrenal_FD_1[["percent.mt"]] <- PercentageFeatureSet(seurat_Adrenal_FD_1, pattern = "^mt-")
################# Adrenal_FD_2
counts_Adrenal_FD_2<- Read10X_h5("Adrenal_FD_2.h5")
seurat_Adrenal_FD_2 <- CreateSeuratObject(counts = counts_Adrenal_FD_2,project = "Adrenal_FD_2")
seurat_Adrenal_FD_2[["percent.mt"]] <- PercentageFeatureSet(seurat_Adrenal_FD_2, pattern = "^mt-")
################# Adrenal_FD_3
counts_Adrenal_FD_3<- Read10X_h5("Adrenal_FD_3.h5")
seurat_Adrenal_FD_3 <- CreateSeuratObject(counts = counts_Adrenal_FD_3,project = "Adrenal_FD_3")
seurat_Adrenal_FD_3[["percent.mt"]] <- PercentageFeatureSet(seurat_Adrenal_FD_3, pattern = "^mt-")
################# Adrenal_FS_1
counts_Adrenal_FS_1<- Read10X_h5("Adrenal_FS_1.h5")
seurat_Adrenal_FS_1 <- CreateSeuratObject(counts = counts_Adrenal_FS_1,project = "Adrenal_FS_1")
seurat_Adrenal_FS_1[["percent.mt"]] <- PercentageFeatureSet(seurat_Adrenal_FS_1, pattern = "^mt-")
################# Adrenal_FS_2
counts_Adrenal_FS_2<- Read10X_h5("Adrenal_FS_2.h5")
seurat_Adrenal_FS_2 <- CreateSeuratObject(counts = counts_Adrenal_FS_2,project = "Adrenal_FS_2")
seurat_Adrenal_FS_2[["percent.mt"]] <- PercentageFeatureSet(seurat_Adrenal_FS_2, pattern = "^mt-")
################# Adrenal_FS_3
counts_Adrenal_FS_3<- Read10X_h5("Adrenal_FS_3.h5")
seurat_Adrenal_FS_3 <- CreateSeuratObject(counts = counts_Adrenal_FS_3,project = "Adrenal_FS_3")
seurat_Adrenal_FS_3[["percent.mt"]] <- PercentageFeatureSet(seurat_Adrenal_FS_3, pattern = "^mt-")
################# Adrenal_MC_1
counts_Adrenal_MC_1<- Read10X_h5("Adrenal_MC_1.h5")
seurat_Adrenal_MC_1 <- CreateSeuratObject(counts = counts_Adrenal_MC_1,project = "Adrenal_MC_1")
seurat_Adrenal_MC_1[["percent.mt"]] <- PercentageFeatureSet(seurat_Adrenal_MC_1, pattern = "^mt-")
################# Adrenal_MC_2
counts_Adrenal_MC_2<- Read10X_h5("Adrenal_MC_2.h5")
seurat_Adrenal_MC_2 <- CreateSeuratObject(counts = counts_Adrenal_MC_2,project = "Adrenal_MC_2")
seurat_Adrenal_MC_2[["percent.mt"]] <- PercentageFeatureSet(seurat_Adrenal_MC_2, pattern = "^mt-")
################# Adrenal_MC_3
counts_Adrenal_MC_3<- Read10X_h5("Adrenal_MC_3.h5")
seurat_Adrenal_MC_3 <- CreateSeuratObject(counts = counts_Adrenal_MC_3,project = "Adrenal_MC_3")
seurat_Adrenal_MC_3[["percent.mt"]] <- PercentageFeatureSet(seurat_Adrenal_MC_3, pattern = "^mt-")
################# Adrenal_MS_1
counts_Adrenal_MS_1<- Read10X_h5("Adrenal_MS_1.h5")
seurat_Adrenal_MS_1 <- CreateSeuratObject(counts = counts_Adrenal_MS_1,project = "Adrenal_MS_1")
seurat_Adrenal_MS_1[["percent.mt"]] <- PercentageFeatureSet(seurat_Adrenal_MS_1, pattern = "^mt-")
################# Adrenal_MS_2
counts_Adrenal_MS_2<- Read10X_h5("Adrenal_MS_2.h5")
seurat_Adrenal_MS_2 <- CreateSeuratObject(counts = counts_Adrenal_MS_2,project = "Adrenal_MS_2")
seurat_Adrenal_MS_2[["percent.mt"]] <- PercentageFeatureSet(seurat_Adrenal_MS_2, pattern = "^mt-")
################# Adrenal_MS_3
counts_Adrenal_MS_3<- Read10X_h5("Adrenal_MS_3.h5")
seurat_Adrenal_MS_3 <- CreateSeuratObject(counts = counts_Adrenal_MS_3,project = "Adrenal_MS_3")
seurat_Adrenal_MS_3[["percent.mt"]] <- PercentageFeatureSet(seurat_Adrenal_MS_3, pattern = "^mt-")

################# Creat a combined seurat based on multiple samples
Adrenal <-merge(seurat_Adrenal_FD_1, y = c(seurat_Adrenal_FD_2,seurat_Adrenal_FD_3,seurat_Adrenal_FS_1,seurat_Adrenal_FS_2,seurat_Adrenal_FS_3,seurat_Adrenal_MC_1,seurat_Adrenal_MC_2,seurat_Adrenal_MC_3,seurat_Adrenal_MS_1,seurat_Adrenal_MS_2,seurat_Adrenal_MS_3), add.cell.ids = c("Adrenal_FD_1","Adrenal_FD_2","Adrenal_FD_3","Adrenal_FS_1","Adrenal_FS_2","Adrenal_FS_3","Adrenal_MC_1","Adrenal_MC_2","Adrenal_MC_3","Adrenal_MS_1","Adrenal_MS_2","Adrenal_MS_3"),project = "Adrenal")

################# 2. Evaluate the potential batch effect and identify the potential low-quality scRNA-seq libraries by using non-integration processing pipeline
#################
#################
#################
#################
#################

################# 2.1 Perform non-integration processing method
Adrenal<- NormalizeData(Adrenal, normalization.method = "LogNormalize", scale.factor = 10000)
Adrenal <- FindVariableFeatures(Adrenal, selection.method = "vst", nfeatures = 2000)
Adrenal <- ScaleData(Adrenal)
Adrenal <- RunPCA(Adrenal, features = VariableFeatures(object =Adrenal))
Adrenal <- RunUMAP(Adrenal, dims = 1:50)
Adrenal <- FindNeighbors(Adrenal, dims = 1:50)
Adrenal <- FindClusters(Adrenal, resolution = 0.10,verbose = TRUE)

################# 2.2 Visualize clustering results
pdf(".../Adrenal/Raw_UMAP_Adrenal.pdf",width=8,height=6)
DimPlot(Adrenal, reduction = "umap",cols = my36colors,raster=FALSE,label = TRUE)
dev.off()

pdf(".../Adrenal/Raw_UMAP_Adrenal_split.pdf",width=14,height=10)
DimPlot(Adrenal, reduction = "umap",cols = my36colors,raster=FALSE,split.by="orig.ident",ncol=4,label = TRUE)
dev.off()

pdf(".../Adrenal/Raw_Featureplot_QC.pdf",width=12,height=5)
FeaturePlot(Adrenal, features = c("nFeature_RNA", "nCount_RNA"),raster=FALSE,min.cutoff =0,reduction = "umap",cols= c("gray", "red"))
dev.off()

pdf(".../Adrenal/Raw_Vln_QC.pdf",width=20,height=6)
VlnPlot(Adrenal, features = c("nFeature_RNA", "nCount_RNA","percent.mt"),group.by="orig.ident",pt.size=0.01)
dev.off()

################# 2.3 Perform analysis of cell compositions
cell.prop.sample<-as.data.frame(prop.table(table(Idents(Adrenal),Adrenal$orig.ident)))
colnames(cell.prop.sample)<-c("cluster","sample","proportion")
pdf(file=".../Adrenal/Raw_cell_proportion_sample.pdf",width=10,height=4)
ggplot(cell.prop.sample,aes(sample,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),axis.text.x = element_text(angle = 45))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=my36colors)
dev.off()

################# 3. Identify and remove empty droplets
#################
#################
#################
#################
#################
set.seed(100)
################# Adrenal_FD_1
# (1). Identify total RNA count for each cell
n_Adrenal_FD_1 <- seurat_Adrenal_FD_1$nCount_RNA 
# (2). Define the second least RNA count as the lower bound 
l_Adrenal_FD_1 <- sort(n_Adrenal_FD_1)[2] 
# (3). Detect the empty droplets
e_Adrenal_FD_1 <- emptyDrops(seurat_Adrenal_FD_1@assays$RNA@counts, lower = l_Adrenal_FD_1,BPPARAM=MulticoreParam())
seurat_Adrenal_FD_1_noempty <- seurat_Adrenal_FD_1[,which(e_Adrenal_FD_1$FDR <= 0.01)]
################# Adrenal_FD_2
# (1). Identify total RNA count for each cell
n_Adrenal_FD_2 <- seurat_Adrenal_FD_2$nCount_RNA
# (2). Define the second least RNA count as the lower bound 
l_Adrenal_FD_2 <- sort(n_Adrenal_FD_2)[2]
# (3). Detect the empty droplets
e_Adrenal_FD_2 <- emptyDrops(seurat_Adrenal_FD_2@assays$RNA@counts, lower = l_Adrenal_FD_2,BPPARAM=MulticoreParam())
seurat_Adrenal_FD_2_noempty <- seurat_Adrenal_FD_2[,which(e_Adrenal_FD_2$FDR <= 0.01)]
################# Adrenal_FD_3
# (1). Identify total RNA count for each cell
n_Adrenal_FD_3 <- seurat_Adrenal_FD_3$nCount_RNA
# (2). Define the second least RNA count as the lower bound 
l_Adrenal_FD_3 <- sort(n_Adrenal_FD_3)[2]
# (3). Detect the empty droplets
e_Adrenal_FD_3 <- emptyDrops(seurat_Adrenal_FD_3@assays$RNA@counts, lower = l_Adrenal_FD_3,BPPARAM=MulticoreParam()) 
seurat_Adrenal_FD_3_noempty <- seurat_Adrenal_FD_3[,which(e_Adrenal_FD_3$FDR <= 0.01)]
#################Adrenal_FS_1
# (1). Identify total RNA count for each cell
n_Adrenal_FS_1 <- seurat_Adrenal_FS_1$nCount_RNA
# (2). Define the second least RNA count as the lower bound 
l_Adrenal_FS_1 <- sort(n_Adrenal_FS_1)[2]
# (3). Detect the empty droplets
e_Adrenal_FS_1 <- emptyDrops(seurat_Adrenal_FS_1@assays$RNA@counts, lower = l_Adrenal_FS_1,BPPARAM=MulticoreParam()) 
seurat_Adrenal_FS_1_noempty <- seurat_Adrenal_FS_1[,which(e_Adrenal_FS_1$FDR <= 0.01)]
#################Adrenal_FS_2
# (1). Identify total RNA count for each cell
n_Adrenal_FS_2 <- seurat_Adrenal_FS_2$nCount_RNA
# (2). Define the second least RNA count as the lower bound 
l_Adrenal_FS_2 <- sort(n_Adrenal_FS_2)[2]
# (3). Detect the empty droplets
e_Adrenal_FS_2 <- emptyDrops(seurat_Adrenal_FS_2@assays$RNA@counts, lower = l_Adrenal_FS_2,BPPARAM=MulticoreParam()) 
seurat_Adrenal_FS_2_noempty <- seurat_Adrenal_FS_2[,which(e_Adrenal_FS_2$FDR <= 0.01)]
#################Adrenal_FS_3
# (1). Identify total RNA count for each cell
n_Adrenal_FS_3 <- seurat_Adrenal_FS_3$nCount_RNA
# (2). Define the second least RNA count as the lower bound 
l_Adrenal_FS_3 <- sort(n_Adrenal_FS_3)[2]
# (3). Detect the empty droplets
e_Adrenal_FS_3 <- emptyDrops(seurat_Adrenal_FS_3@assays$RNA@counts, lower = l_Adrenal_FS_3,BPPARAM=MulticoreParam()) 
seurat_Adrenal_FS_3_noempty <- seurat_Adrenal_FS_3[,which(e_Adrenal_FS_3$FDR <= 0.01)]
#################Adrenal_MC_1
# (1). Identify total RNA count for each cell
n_Adrenal_MC_1 <- seurat_Adrenal_MC_1$nCount_RNA
# (2). Define the second least RNA count as the lower bound 
l_Adrenal_MC_1 <- sort(n_Adrenal_MC_1)[2]
# (3). Detect the empty droplets
e_Adrenal_MC_1 <- emptyDrops(seurat_Adrenal_MC_1@assays$RNA@counts, lower = l_Adrenal_MC_1,BPPARAM=MulticoreParam())
seurat_Adrenal_MC_1_noempty <- seurat_Adrenal_MC_1[,which(e_Adrenal_MC_1$FDR <= 0.01)]
#################Adrenal_MC_2
# (1). Identify total RNA count for each cell
n_Adrenal_MC_2 <- seurat_Adrenal_MC_2$nCount_RNA
# (2). Define the second least RNA count as the lower bound 
l_Adrenal_MC_2 <- sort(n_Adrenal_MC_2)[2]
# (3). Detect the empty droplets
e_Adrenal_MC_2 <- emptyDrops(seurat_Adrenal_MC_2@assays$RNA@counts, lower = l_Adrenal_MC_2,BPPARAM=MulticoreParam())
seurat_Adrenal_MC_2_noempty <- seurat_Adrenal_MC_2[,which(e_Adrenal_MC_2$FDR <= 0.01)]
#################Adrenal_MC_3
# (1). Identify total RNA count for each cell
n_Adrenal_MC_3 <- seurat_Adrenal_MC_3$nCount_RNA
# (2). Define the second least RNA count as the lower bound 
l_Adrenal_MC_3 <- sort(n_Adrenal_MC_3)[2]
# (3). Detect the empty droplets
e_Adrenal_MC_3 <- emptyDrops(seurat_Adrenal_MC_3@assays$RNA@counts, lower = l_Adrenal_MC_3,BPPARAM=MulticoreParam())
seurat_Adrenal_MC_3_noempty <- seurat_Adrenal_MC_3[,which(e_Adrenal_MC_3$FDR <= 0.01)]
#################Adrenal_MS_1
# (1). Identify total RNA count for each cell
n_Adrenal_MS_1 <- seurat_Adrenal_MS_1$nCount_RNA
# (2). Define the second least RNA count as the lower bound 
l_Adrenal_MS_1 <- sort(n_Adrenal_MS_1)[2]
# (3). Detect the empty droplets
e_Adrenal_MS_1 <- emptyDrops(seurat_Adrenal_MS_1@assays$RNA@counts, lower = l_Adrenal_MS_1,BPPARAM=MulticoreParam())
seurat_Adrenal_MS_1_noempty <- seurat_Adrenal_MS_1[,which(e_Adrenal_MS_1$FDR <= 0.01)]
#################Adrenal_MS_2
# (1). Identify total RNA count for each cell
n_Adrenal_MS_2 <- seurat_Adrenal_MS_2$nCount_RNA
# (2). Define the second least RNA count as the lower bound 
l_Adrenal_MS_2 <- sort(n_Adrenal_MS_2)[2]
# (3). Detect the empty droplets
e_Adrenal_MS_2 <- emptyDrops(seurat_Adrenal_MS_2@assays$RNA@counts, lower = l_Adrenal_MS_2,BPPARAM=MulticoreParam())
seurat_Adrenal_MS_2_noempty <- seurat_Adrenal_MS_2[,which(e_Adrenal_MS_2$FDR <= 0.01)]
#################Adrenal_MS_3
# (1). Identify total RNA count for each cell
n_Adrenal_MS_3 <- seurat_Adrenal_MS_3$nCount_RNA
# (2). Define the second least RNA count as the lower bound 
l_Adrenal_MS_3 <- sort(n_Adrenal_MS_3)[2]
# (3). Detect the empty droplets
e_Adrenal_MS_3 <- emptyDrops(seurat_Adrenal_MS_3@assays$RNA@counts, lower = l_Adrenal_MS_3,BPPARAM=MulticoreParam())
seurat_Adrenal_MS_3_noempty <- seurat_Adrenal_MS_3[,which(e_Adrenal_MS_3$FDR <= 0.01)]

#################4. Identify and remove low-RNA content cells
#################
#################
#################
#################
#################

seurat_Adrenal_FD_1_high<-subset(seurat_Adrenal_FD_1_noempty,subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 15)
seurat_Adrenal_FD_2_high<-subset(seurat_Adrenal_FD_2_noempty,subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 15)
seurat_Adrenal_FD_3_high<-subset(seurat_Adrenal_FD_3_noempty,subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 15)
seurat_Adrenal_FS_1_high<-subset(seurat_Adrenal_FS_1_noempty,subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 15)
seurat_Adrenal_FS_2_high<-subset(seurat_Adrenal_FS_2_noempty,subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 15)
seurat_Adrenal_FS_3_high<-subset(seurat_Adrenal_FS_3_noempty,subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 15)
seurat_Adrenal_MC_1_high<-subset(seurat_Adrenal_MC_1_noempty,subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 15)
seurat_Adrenal_MC_2_high<-subset(seurat_Adrenal_MC_2_noempty,subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 15)
seurat_Adrenal_MC_3_high<-subset(seurat_Adrenal_MC_3_noempty,subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 15)
seurat_Adrenal_MS_1_high<-subset(seurat_Adrenal_MS_1_noempty,subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 15)
seurat_Adrenal_MS_2_high<-subset(seurat_Adrenal_MS_2_noempty,subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 15)
seurat_Adrenal_MS_3_high<-subset(seurat_Adrenal_MS_3_noempty,subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 15)

#################5. Identify and remove doublets
#################
#################
#################
#################
#################

#################Adrenal_FD_1
# (1). Pre-process Seurat object (standard) 
seurat_Adrenal_FD_1_high <- NormalizeData(seurat_Adrenal_FD_1_high)
seurat_Adrenal_FD_1_high <- FindVariableFeatures(seurat_Adrenal_FD_1_high, selection.method = "vst", nfeatures = 2000)
seurat_Adrenal_FD_1_high <- ScaleData(seurat_Adrenal_FD_1_high)
seurat_Adrenal_FD_1_high <- RunPCA(seurat_Adrenal_FD_1_high)
seurat_Adrenal_FD_1_high <- RunUMAP(seurat_Adrenal_FD_1_high, dims = 1:10)
# (2). pK Identification (no ground-truth)
sweep.res.list_Adrenal_FD_1 <- paramSweep_v3(seurat_Adrenal_FD_1_high, PCs = 1:10, sct = FALSE)
sweep.stats_Adrenal_FD_1 <- summarizeSweep(sweep.res.list_Adrenal_FD_1, GT = FALSE)
pdf(file=".../Adrenal/doublet/pK_Adrenal_FD_1.pdf",width=6,height=8)
bcmvn_Adrenal_FD_1 <- find.pK(sweep.stats_Adrenal_FD_1)
dev.off()
mpK_Adrenal_FD_1<-as.numeric(as.vector(bcmvn_Adrenal_FD_1$pK[which.max(bcmvn_Adrenal_FD_1$BCmetric)]))
# (3). Run DoubletFinder
DoubletRate_Adrenal_FD_1 = ncol(seurat_Adrenal_FD_1_high)*8*1e-6 # ~1000 cells, and a multiplet rate of ~0.8% according to 10XGenomics
nExp_poi_Adrenal_FD_1 <- round(DoubletRate_Adrenal_FD_1*ncol(seurat_Adrenal_FD_1_high)) 
seurat_Adrenal_FD_1_high <- doubletFinder_v3(seurat_Adrenal_FD_1_high, PCs = 1:10, pN = 0.25, pK = mpK_Adrenal_FD_1, nExp = nExp_poi_Adrenal_FD_1, reuse.pANN = FALSE, sct = F)
# (4). Remove doublets
table(seurat_Adrenal_FD_1_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_FD_1,"_",nExp_poi_Adrenal_FD_1)])
seurat_Adrenal_FD_1_high$doubletfinder<-seurat_Adrenal_FD_1_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_FD_1,"_",nExp_poi_Adrenal_FD_1)]
seurat_Adrenal_FD_1_singlet<-subset(seurat_Adrenal_FD_1_high,doubletfinder=="Singlet")
#################Adrenal_FD_2
# (1). Pre-process Seurat object (standard) 
seurat_Adrenal_FD_2_high <- NormalizeData(seurat_Adrenal_FD_2_high)
seurat_Adrenal_FD_2_high <- FindVariableFeatures(seurat_Adrenal_FD_2_high, selection.method = "vst", nfeatures = 2000)
seurat_Adrenal_FD_2_high <- ScaleData(seurat_Adrenal_FD_2_high)
seurat_Adrenal_FD_2_high <- RunPCA(seurat_Adrenal_FD_2_high)
seurat_Adrenal_FD_2_high <- RunUMAP(seurat_Adrenal_FD_2_high, dims = 1:10)
# (2). pK Identification (no ground-truth)
sweep.res.list_Adrenal_FD_2 <- paramSweep_v3(seurat_Adrenal_FD_2_high, PCs = 1:10, sct = FALSE)
sweep.stats_Adrenal_FD_2 <- summarizeSweep(sweep.res.list_Adrenal_FD_2, GT = FALSE)
pdf(file=".../Adrenal/doublet/pK_Adrenal_FD_2.pdf",width=6,height=8)
bcmvn_Adrenal_FD_2 <- find.pK(sweep.stats_Adrenal_FD_2)
dev.off()
mpK_Adrenal_FD_2<-as.numeric(as.vector(bcmvn_Adrenal_FD_2$pK[which.max(bcmvn_Adrenal_FD_2$BCmetric)]))
# (3). Run DoubletFinder
DoubletRate_Adrenal_FD_2 = ncol(seurat_Adrenal_FD_2_high)*8*1e-6 # ~1000 cells, and a multiplet rate of ~0.8% according to 10XGenomics
nExp_poi_Adrenal_FD_2 <- round(DoubletRate_Adrenal_FD_2*ncol(seurat_Adrenal_FD_2_high)) 
seurat_Adrenal_FD_2_high <- doubletFinder_v3(seurat_Adrenal_FD_2_high, PCs = 1:10, pN = 0.25, pK = mpK_Adrenal_FD_2, nExp = nExp_poi_Adrenal_FD_2, reuse.pANN = FALSE, sct = F)
# (4). Remove doublets
table(seurat_Adrenal_FD_2_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_FD_2,"_",nExp_poi_Adrenal_FD_2)])
seurat_Adrenal_FD_2_high$doubletfinder<-seurat_Adrenal_FD_2_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_FD_2,"_",nExp_poi_Adrenal_FD_2)]
seurat_Adrenal_FD_2_singlet<-subset(seurat_Adrenal_FD_2_high,doubletfinder=="Singlet")
#################Adrenal_FD_3
# (1). Pre-process Seurat object (standard) 
seurat_Adrenal_FD_3_high <- NormalizeData(seurat_Adrenal_FD_3_high)
seurat_Adrenal_FD_3_high <- FindVariableFeatures(seurat_Adrenal_FD_3_high, selection.method = "vst", nfeatures = 2000)
seurat_Adrenal_FD_3_high <- ScaleData(seurat_Adrenal_FD_3_high)
seurat_Adrenal_FD_3_high <- RunPCA(seurat_Adrenal_FD_3_high)
seurat_Adrenal_FD_3_high <- RunUMAP(seurat_Adrenal_FD_3_high, dims = 1:10)
# (2). pK Identification (no ground-truth)
sweep.res.list_Adrenal_FD_3 <- paramSweep_v3(seurat_Adrenal_FD_3_high, PCs = 1:10, sct = FALSE)
sweep.stats_Adrenal_FD_3 <- summarizeSweep(sweep.res.list_Adrenal_FD_3, GT = FALSE)
pdf(file=".../Adrenal/doublet/pK_Adrenal_FD_3.pdf",width=6,height=8)
bcmvn_Adrenal_FD_3 <- find.pK(sweep.stats_Adrenal_FD_3)
dev.off()
mpK_Adrenal_FD_3<-as.numeric(as.vector(bcmvn_Adrenal_FD_3$pK[which.max(bcmvn_Adrenal_FD_3$BCmetric)]))
# (3). Run DoubletFinder
DoubletRate_Adrenal_FD_3 = ncol(seurat_Adrenal_FD_3_high)*8*1e-6 # ~1000 cells, and a multiplet rate of ~0.8% according to 10XGenomics
nExp_poi_Adrenal_FD_3 <- round(DoubletRate_Adrenal_FD_3*ncol(seurat_Adrenal_FD_3_high)) 
seurat_Adrenal_FD_3_high <- doubletFinder_v3(seurat_Adrenal_FD_3_high, PCs = 1:10, pN = 0.25, pK = mpK_Adrenal_FD_3, nExp = nExp_poi_Adrenal_FD_3, reuse.pANN = FALSE, sct = F)
# (4). Remove doublets
table(seurat_Adrenal_FD_3_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_FD_3,"_",nExp_poi_Adrenal_FD_3)])
seurat_Adrenal_FD_3_high$doubletfinder<-seurat_Adrenal_FD_3_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_FD_3,"_",nExp_poi_Adrenal_FD_3)]
seurat_Adrenal_FD_3_singlet<-subset(seurat_Adrenal_FD_3_high,doubletfinder=="Singlet")
#################Adrenal_FS_1
# (1). Pre-process Seurat object (standard) 
seurat_Adrenal_FS_1_high <- NormalizeData(seurat_Adrenal_FS_1_high)
seurat_Adrenal_FS_1_high <- FindVariableFeatures(seurat_Adrenal_FS_1_high, selection.method = "vst", nfeatures = 2000)
seurat_Adrenal_FS_1_high <- ScaleData(seurat_Adrenal_FS_1_high)
seurat_Adrenal_FS_1_high <- RunPCA(seurat_Adrenal_FS_1_high)
seurat_Adrenal_FS_1_high <- RunUMAP(seurat_Adrenal_FS_1_high, dims = 1:10)
# (2). pK Identification (no ground-truth)
sweep.res.list_Adrenal_FS_1 <- paramSweep_v3(seurat_Adrenal_FS_1_high, PCs = 1:10, sct = FALSE)
sweep.stats_Adrenal_FS_1 <- summarizeSweep(sweep.res.list_Adrenal_FS_1, GT = FALSE)
pdf(file=".../Adrenal/doublet/pK_Adrenal_FS_1.pdf",width=6,height=8)
bcmvn_Adrenal_FS_1 <- find.pK(sweep.stats_Adrenal_FS_1)
dev.off()
mpK_Adrenal_FS_1<-as.numeric(as.vector(bcmvn_Adrenal_FS_1$pK[which.max(bcmvn_Adrenal_FS_1$BCmetric)]))
# (3). Run DoubletFinder
DoubletRate_Adrenal_FS_1 = ncol(seurat_Adrenal_FS_1_high)*8*1e-6 # ~1000 cells, and a multiplet rate of ~0.8% according to 10XGenomics
nExp_poi_Adrenal_FS_1 <- round(DoubletRate_Adrenal_FS_1*ncol(seurat_Adrenal_FS_1_high)) 
seurat_Adrenal_FS_1_high <- doubletFinder_v3(seurat_Adrenal_FS_1_high, PCs = 1:10, pN = 0.25, pK = mpK_Adrenal_FS_1, nExp = nExp_poi_Adrenal_FS_1, reuse.pANN = FALSE, sct = F)
# (4). Remove doublets
table(seurat_Adrenal_FS_1_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_FS_1,"_",nExp_poi_Adrenal_FS_1)])
seurat_Adrenal_FS_1_high$doubletfinder<-seurat_Adrenal_FS_1_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_FS_1,"_",nExp_poi_Adrenal_FS_1)]
seurat_Adrenal_FS_1_singlet<-subset(seurat_Adrenal_FS_1_high,doubletfinder=="Singlet")
#################Adrenal_FS_2
# (1). Pre-process Seurat object (standard) 
seurat_Adrenal_FS_2_high <- NormalizeData(seurat_Adrenal_FS_2_high)
seurat_Adrenal_FS_2_high <- FindVariableFeatures(seurat_Adrenal_FS_2_high, selection.method = "vst", nfeatures = 2000)
seurat_Adrenal_FS_2_high <- ScaleData(seurat_Adrenal_FS_2_high)
seurat_Adrenal_FS_2_high <- RunPCA(seurat_Adrenal_FS_2_high)
seurat_Adrenal_FS_2_high <- RunUMAP(seurat_Adrenal_FS_2_high, dims = 1:10)
# (2). pK Identification (no ground-truth)
sweep.res.list_Adrenal_FS_2 <- paramSweep_v3(seurat_Adrenal_FS_2_high, PCs = 1:10, sct = FALSE)
sweep.stats_Adrenal_FS_2 <- summarizeSweep(sweep.res.list_Adrenal_FS_2, GT = FALSE)
pdf(file=".../Adrenal/doublet/pK_Adrenal_FS_2.pdf",width=6,height=8)
bcmvn_Adrenal_FS_2 <- find.pK(sweep.stats_Adrenal_FS_2)
dev.off()
mpK_Adrenal_FS_2<-as.numeric(as.vector(bcmvn_Adrenal_FS_2$pK[which.max(bcmvn_Adrenal_FS_2$BCmetric)]))
# (3). Run DoubletFinder
DoubletRate_Adrenal_FS_2 = ncol(seurat_Adrenal_FS_2_high)*8*1e-6 # ~1000 cells, and a multiplet rate of ~0.8% according to 10XGenomics
nExp_poi_Adrenal_FS_2 <- round(DoubletRate_Adrenal_FS_2*ncol(seurat_Adrenal_FS_2_high)) 
seurat_Adrenal_FS_2_high <- doubletFinder_v3(seurat_Adrenal_FS_2_high, PCs = 1:10, pN = 0.25, pK = mpK_Adrenal_FS_2, nExp = nExp_poi_Adrenal_FS_2, reuse.pANN = FALSE, sct = F)
# (4). Remove doublets
table(seurat_Adrenal_FS_2_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_FS_2,"_",nExp_poi_Adrenal_FS_2)])
seurat_Adrenal_FS_2_high$doubletfinder<-seurat_Adrenal_FS_2_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_FS_2,"_",nExp_poi_Adrenal_FS_2)]
seurat_Adrenal_FS_2_singlet<-subset(seurat_Adrenal_FS_2_high,doubletfinder=="Singlet")
#################Adrenal_FS_3
# (1). Pre-process Seurat object (standard) 
seurat_Adrenal_FS_3_high <- NormalizeData(seurat_Adrenal_FS_3_high)
seurat_Adrenal_FS_3_high <- FindVariableFeatures(seurat_Adrenal_FS_3_high, selection.method = "vst", nfeatures = 2000)
seurat_Adrenal_FS_3_high <- ScaleData(seurat_Adrenal_FS_3_high)
seurat_Adrenal_FS_3_high <- RunPCA(seurat_Adrenal_FS_3_high)
seurat_Adrenal_FS_3_high <- RunUMAP(seurat_Adrenal_FS_3_high, dims = 1:10)
# (2). pK Identification (no ground-truth)
sweep.res.list_Adrenal_FS_3 <- paramSweep_v3(seurat_Adrenal_FS_3_high, PCs = 1:10, sct = FALSE)
sweep.stats_Adrenal_FS_3 <- summarizeSweep(sweep.res.list_Adrenal_FS_3, GT = FALSE)
pdf(file=".../Adrenal/doublet/pK_Adrenal_FS_3.pdf",width=6,height=8)
bcmvn_Adrenal_FS_3 <- find.pK(sweep.stats_Adrenal_FS_3)
dev.off()
mpK_Adrenal_FS_3<-as.numeric(as.vector(bcmvn_Adrenal_FS_3$pK[which.max(bcmvn_Adrenal_FS_3$BCmetric)]))
# (3). Run DoubletFinder
DoubletRate_Adrenal_FS_3 = ncol(seurat_Adrenal_FS_3_high)*8*1e-6 # ~1000 cells, and a multiplet rate of ~0.8% according to 10XGenomics
nExp_poi_Adrenal_FS_3 <- round(DoubletRate_Adrenal_FS_3*ncol(seurat_Adrenal_FS_3_high)) 
seurat_Adrenal_FS_3_high <- doubletFinder_v3(seurat_Adrenal_FS_3_high, PCs = 1:10, pN = 0.25, pK = mpK_Adrenal_FS_3, nExp = nExp_poi_Adrenal_FS_3, reuse.pANN = FALSE, sct = F)
# (4). Remove doublets
table(seurat_Adrenal_FS_3_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_FS_3,"_",nExp_poi_Adrenal_FS_3)])
seurat_Adrenal_FS_3_high$doubletfinder<-seurat_Adrenal_FS_3_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_FS_3,"_",nExp_poi_Adrenal_FS_3)]
seurat_Adrenal_FS_3_singlet<-subset(seurat_Adrenal_FS_3_high,doubletfinder=="Singlet")
#################Adrenal_MC_1
# (1). Pre-process Seurat object (standard) 
seurat_Adrenal_MC_1_high <- NormalizeData(seurat_Adrenal_MC_1_high)
seurat_Adrenal_MC_1_high <- FindVariableFeatures(seurat_Adrenal_MC_1_high, selection.method = "vst", nfeatures = 2000)
seurat_Adrenal_MC_1_high <- ScaleData(seurat_Adrenal_MC_1_high)
seurat_Adrenal_MC_1_high <- RunPCA(seurat_Adrenal_MC_1_high)
seurat_Adrenal_MC_1_high <- RunUMAP(seurat_Adrenal_MC_1_high, dims = 1:10)
# (2). pK Identification (no ground-truth)
sweep.res.list_Adrenal_MC_1 <- paramSweep_v3(seurat_Adrenal_MC_1_high, PCs = 1:10, sct = FALSE)
sweep.stats_Adrenal_MC_1 <- summarizeSweep(sweep.res.list_Adrenal_MC_1, GT = FALSE)
pdf(file=".../Adrenal/doublet/pK_Adrenal_MC_1.pdf",width=6,height=8)
bcmvn_Adrenal_MC_1 <- find.pK(sweep.stats_Adrenal_MC_1)
dev.off()
mpK_Adrenal_MC_1<-as.numeric(as.vector(bcmvn_Adrenal_MC_1$pK[which.max(bcmvn_Adrenal_MC_1$BCmetric)]))
# (3). Run DoubletFinder
DoubletRate_Adrenal_MC_1 = ncol(seurat_Adrenal_MC_1_high)*8*1e-6 # ~1000 cells, and a multiplet rate of ~0.8% according to 10XGenomics
nExp_poi_Adrenal_MC_1 <- round(DoubletRate_Adrenal_MC_1*ncol(seurat_Adrenal_MC_1_high)) 
seurat_Adrenal_MC_1_high <- doubletFinder_v3(seurat_Adrenal_MC_1_high, PCs = 1:10, pN = 0.25, pK = mpK_Adrenal_MC_1, nExp = nExp_poi_Adrenal_MC_1, reuse.pANN = FALSE, sct = F)
# (4). Remove doublets
table(seurat_Adrenal_MC_1_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_MC_1,"_",nExp_poi_Adrenal_MC_1)])
seurat_Adrenal_MC_1_high$doubletfinder<-seurat_Adrenal_MC_1_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_MC_1,"_",nExp_poi_Adrenal_MC_1)]
seurat_Adrenal_MC_1_singlet<-subset(seurat_Adrenal_MC_1_high,doubletfinder=="Singlet")
#################Adrenal_MC_2
# (1). Pre-process Seurat object (standard) 
seurat_Adrenal_MC_2_high <- NormalizeData(seurat_Adrenal_MC_2_high)
seurat_Adrenal_MC_2_high <- FindVariableFeatures(seurat_Adrenal_MC_2_high, selection.method = "vst", nfeatures = 2000)
seurat_Adrenal_MC_2_high <- ScaleData(seurat_Adrenal_MC_2_high)
seurat_Adrenal_MC_2_high <- RunPCA(seurat_Adrenal_MC_2_high)
seurat_Adrenal_MC_2_high <- RunUMAP(seurat_Adrenal_MC_2_high, dims = 1:10)
# (2). pK Identification (no ground-truth)
sweep.res.list_Adrenal_MC_2 <- paramSweep_v3(seurat_Adrenal_MC_2_high, PCs = 1:10, sct = FALSE)
sweep.stats_Adrenal_MC_2 <- summarizeSweep(sweep.res.list_Adrenal_MC_2, GT = FALSE)
pdf(file=".../Adrenal/doublet/pK_Adrenal_MC_2.pdf",width=6,height=8)
bcmvn_Adrenal_MC_2 <- find.pK(sweep.stats_Adrenal_MC_2)
dev.off()
mpK_Adrenal_MC_2<-as.numeric(as.vector(bcmvn_Adrenal_MC_2$pK[which.max(bcmvn_Adrenal_MC_2$BCmetric)]))
# (3). Run DoubletFinder
DoubletRate_Adrenal_MC_2 = ncol(seurat_Adrenal_MC_2_high)*8*1e-6 # ~1000 cells, and a multiplet rate of ~0.8% according to 10XGenomics
nExp_poi_Adrenal_MC_2 <- round(DoubletRate_Adrenal_MC_2*ncol(seurat_Adrenal_MC_2_high)) 
seurat_Adrenal_MC_2_high <- doubletFinder_v3(seurat_Adrenal_MC_2_high, PCs = 1:10, pN = 0.25, pK = mpK_Adrenal_MC_2, nExp = nExp_poi_Adrenal_MC_2, reuse.pANN = FALSE, sct = F)
# (4). Remove doublets
table(seurat_Adrenal_MC_2_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_MC_2,"_",nExp_poi_Adrenal_MC_2)])
seurat_Adrenal_MC_2_high$doubletfinder<-seurat_Adrenal_MC_2_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_MC_2,"_",nExp_poi_Adrenal_MC_2)]
seurat_Adrenal_MC_2_singlet<-subset(seurat_Adrenal_MC_2_high,doubletfinder=="Singlet")
#################Adrenal_MC_3
# (1). Pre-process Seurat object (standard) 
seurat_Adrenal_MC_3_high <- NormalizeData(seurat_Adrenal_MC_3_high)
seurat_Adrenal_MC_3_high <- FindVariableFeatures(seurat_Adrenal_MC_3_high, selection.method = "vst", nfeatures = 2000)
seurat_Adrenal_MC_3_high <- ScaleData(seurat_Adrenal_MC_3_high)
seurat_Adrenal_MC_3_high <- RunPCA(seurat_Adrenal_MC_3_high)
seurat_Adrenal_MC_3_high <- RunUMAP(seurat_Adrenal_MC_3_high, dims = 1:10)
# (2). pK Identification (no ground-truth)
sweep.res.list_Adrenal_MC_3 <- paramSweep_v3(seurat_Adrenal_MC_3_high, PCs = 1:10, sct = FALSE)
sweep.stats_Adrenal_MC_3 <- summarizeSweep(sweep.res.list_Adrenal_MC_3, GT = FALSE)
pdf(file=".../Adrenal/doublet/pK_Adrenal_MC_3.pdf",width=6,height=8)
bcmvn_Adrenal_MC_3 <- find.pK(sweep.stats_Adrenal_MC_3)
dev.off()
mpK_Adrenal_MC_3<-as.numeric(as.vector(bcmvn_Adrenal_MC_3$pK[which.max(bcmvn_Adrenal_MC_3$BCmetric)]))
# (3). Run DoubletFinder
DoubletRate_Adrenal_MC_3 = ncol(seurat_Adrenal_MC_3_high)*8*1e-6 # ~1000 cells, and a multiplet rate of ~0.8% according to 10XGenomics
nExp_poi_Adrenal_MC_3 <- round(DoubletRate_Adrenal_MC_3*ncol(seurat_Adrenal_MC_3_high)) 
seurat_Adrenal_MC_3_high <- doubletFinder_v3(seurat_Adrenal_MC_3_high, PCs = 1:10, pN = 0.25, pK = mpK_Adrenal_MC_3, nExp = nExp_poi_Adrenal_MC_3, reuse.pANN = FALSE, sct = F)
# (4). Remove doublets
table(seurat_Adrenal_MC_3_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_MC_3,"_",nExp_poi_Adrenal_MC_3)])
seurat_Adrenal_MC_3_high$doubletfinder<-seurat_Adrenal_MC_3_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_MC_3,"_",nExp_poi_Adrenal_MC_3)]
seurat_Adrenal_MC_3_singlet<-subset(seurat_Adrenal_MC_3_high,doubletfinder=="Singlet")
#################Adrenal_MS_1
# (1). Pre-process Seurat object (standard) 
seurat_Adrenal_MS_1_high <- NormalizeData(seurat_Adrenal_MS_1_high)
seurat_Adrenal_MS_1_high <- FindVariableFeatures(seurat_Adrenal_MS_1_high, selection.method = "vst", nfeatures = 2000)
seurat_Adrenal_MS_1_high <- ScaleData(seurat_Adrenal_MS_1_high)
seurat_Adrenal_MS_1_high <- RunPCA(seurat_Adrenal_MS_1_high)
seurat_Adrenal_MS_1_high <- RunUMAP(seurat_Adrenal_MS_1_high, dims = 1:10)
# (2). pK Identification (no ground-truth)
sweep.res.list_Adrenal_MS_1 <- paramSweep_v3(seurat_Adrenal_MS_1_high, PCs = 1:10, sct = FALSE)
sweep.stats_Adrenal_MS_1 <- summarizeSweep(sweep.res.list_Adrenal_MS_1, GT = FALSE)
pdf(file=".../Adrenal/doublet/pK_Adrenal_MS_1.pdf",width=6,height=8)
bcmvn_Adrenal_MS_1 <- find.pK(sweep.stats_Adrenal_MS_1)
dev.off()
mpK_Adrenal_MS_1<-as.numeric(as.vector(bcmvn_Adrenal_MS_1$pK[which.max(bcmvn_Adrenal_MS_1$BCmetric)]))
# (3). Run DoubletFinder
DoubletRate_Adrenal_MS_1 = ncol(seurat_Adrenal_MS_1_high)*8*1e-6 # ~1000 cells, and a multiplet rate of ~0.8% according to 10XGenomics
nExp_poi_Adrenal_MS_1 <- round(DoubletRate_Adrenal_MS_1*ncol(seurat_Adrenal_MS_1_high)) 
seurat_Adrenal_MS_1_high <- doubletFinder_v3(seurat_Adrenal_MS_1_high, PCs = 1:10, pN = 0.25, pK = mpK_Adrenal_MS_1, nExp = nExp_poi_Adrenal_MS_1, reuse.pANN = FALSE, sct = F)
# (4). Remove doublets
table(seurat_Adrenal_MS_1_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_MS_1,"_",nExp_poi_Adrenal_MS_1)])
seurat_Adrenal_MS_1_high$doubletfinder<-seurat_Adrenal_MS_1_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_MS_1,"_",nExp_poi_Adrenal_MS_1)]
seurat_Adrenal_MS_1_singlet<-subset(seurat_Adrenal_MS_1_high,doubletfinder=="Singlet")
#################Adrenal_MS_2
# (1). Pre-process Seurat object (standard) 
seurat_Adrenal_MS_2_high <- NormalizeData(seurat_Adrenal_MS_2_high)
seurat_Adrenal_MS_2_high <- FindVariableFeatures(seurat_Adrenal_MS_2_high, selection.method = "vst", nfeatures = 2000)
seurat_Adrenal_MS_2_high <- ScaleData(seurat_Adrenal_MS_2_high)
seurat_Adrenal_MS_2_high <- RunPCA(seurat_Adrenal_MS_2_high)
seurat_Adrenal_MS_2_high <- RunUMAP(seurat_Adrenal_MS_2_high, dims = 1:10)
# (2). pK Identification (no ground-truth)
sweep.res.list_Adrenal_MS_2 <- paramSweep_v3(seurat_Adrenal_MS_2_high, PCs = 1:10, sct = FALSE)
sweep.stats_Adrenal_MS_2 <- summarizeSweep(sweep.res.list_Adrenal_MS_2, GT = FALSE)
pdf(file=".../Adrenal/doublet/pK_Adrenal_MS_2.pdf",width=6,height=8)
bcmvn_Adrenal_MS_2 <- find.pK(sweep.stats_Adrenal_MS_2)
dev.off()
mpK_Adrenal_MS_2<-as.numeric(as.vector(bcmvn_Adrenal_MS_2$pK[which.max(bcmvn_Adrenal_MS_2$BCmetric)]))
# (3). Run DoubletFinder
DoubletRate_Adrenal_MS_2 = ncol(seurat_Adrenal_MS_2_high)*8*1e-6 # ~1000 cells, and a multiplet rate of ~0.8% according to 10XGenomics
nExp_poi_Adrenal_MS_2 <- round(DoubletRate_Adrenal_MS_2*ncol(seurat_Adrenal_MS_2_high)) 
seurat_Adrenal_MS_2_high <- doubletFinder_v3(seurat_Adrenal_MS_2_high, PCs = 1:10, pN = 0.25, pK = mpK_Adrenal_MS_2, nExp = nExp_poi_Adrenal_MS_2, reuse.pANN = FALSE, sct = F)
# (4). Remove doublets
table(seurat_Adrenal_MS_2_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_MS_2,"_",nExp_poi_Adrenal_MS_2)])
seurat_Adrenal_MS_2_high$doubletfinder<-seurat_Adrenal_MS_2_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_MS_2,"_",nExp_poi_Adrenal_MS_2)]
seurat_Adrenal_MS_2_singlet<-subset(seurat_Adrenal_MS_2_high,doubletfinder=="Singlet")
#################Adrenal_MS_3
# (1). Pre-process Seurat object (standard) 
seurat_Adrenal_MS_3_high <- NormalizeData(seurat_Adrenal_MS_3_high)
seurat_Adrenal_MS_3_high <- FindVariableFeatures(seurat_Adrenal_MS_3_high, selection.method = "vst", nfeatures = 2000)
seurat_Adrenal_MS_3_high <- ScaleData(seurat_Adrenal_MS_3_high)
seurat_Adrenal_MS_3_high <- RunPCA(seurat_Adrenal_MS_3_high)
seurat_Adrenal_MS_3_high <- RunUMAP(seurat_Adrenal_MS_3_high, dims = 1:10)
# (2). pK Identification (no ground-truth)
sweep.res.list_Adrenal_MS_3 <- paramSweep_v3(seurat_Adrenal_MS_3_high, PCs = 1:10, sct = FALSE)
sweep.stats_Adrenal_MS_3 <- summarizeSweep(sweep.res.list_Adrenal_MS_3, GT = FALSE)
pdf(file=".../Adrenal/doublet/pK_Adrenal_MS_3.pdf",width=6,height=8)
bcmvn_Adrenal_MS_3 <- find.pK(sweep.stats_Adrenal_MS_3)
dev.off()
mpK_Adrenal_MS_3<-as.numeric(as.vector(bcmvn_Adrenal_MS_3$pK[which.max(bcmvn_Adrenal_MS_3$BCmetric)]))
# (3). Run DoubletFinder
DoubletRate_Adrenal_MS_3 = ncol(seurat_Adrenal_MS_3_high)*8*1e-6 # ~1000 cells, and a multiplet rate of ~0.8% according to 10XGenomics
nExp_poi_Adrenal_MS_3 <- round(DoubletRate_Adrenal_MS_3*ncol(seurat_Adrenal_MS_3_high)) 
seurat_Adrenal_MS_3_high <- doubletFinder_v3(seurat_Adrenal_MS_3_high, PCs = 1:10, pN = 0.25, pK = mpK_Adrenal_MS_3, nExp = nExp_poi_Adrenal_MS_3, reuse.pANN = FALSE, sct = F)
# (4). Remove doublets
table(seurat_Adrenal_MS_3_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_MS_3,"_",nExp_poi_Adrenal_MS_3)])
seurat_Adrenal_MS_3_high$doubletfinder<-seurat_Adrenal_MS_3_high@meta.data[,paste0("DF.classifications_0.25_",mpK_Adrenal_MS_3,"_",nExp_poi_Adrenal_MS_3)]
seurat_Adrenal_MS_3_singlet<-subset(seurat_Adrenal_MS_3_high,doubletfinder=="Singlet")

#################5. Integrate data
#################
#################
#################
#################
#################

#################5.1 Merge multiple objects into a combined object
Adrenal.combine <-merge(seurat_Adrenal_FD_1_singlet, y = c(seurat_Adrenal_FD_2_singlet,seurat_Adrenal_FD_3_singlet,seurat_Adrenal_FS_1_singlet,seurat_Adrenal_FS_2_singlet,seurat_Adrenal_FS_3_singlet,seurat_Adrenal_MC_1_singlet,seurat_Adrenal_MC_2_singlet,seurat_Adrenal_MC_3_singlet,seurat_Adrenal_MS_1_singlet,seurat_Adrenal_MS_2_singlet,seurat_Adrenal_MS_3_singlet), add.cell.ids = c("Adrenal_FD_1","Adrenal_FD_2","Adrenal_FD_3","Adrenal_FS_1","Adrenal_FS_2","Adrenal_FS_3","Adrenal_MC_1","Adrenal_MC_2","Adrenal_MC_3","Adrenal_MS_1","Adrenal_MS_2","Adrenal_MS_3"),project = "Adrenal.combine")
table(Adrenal.combine$orig.ident)
#################5.2 Integrate data
Adrenal.combine.list <- SplitObject(Adrenal.combine, split.by = "orig.ident")
Adrenal.combine.list <- lapply(X = Adrenal.combine.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = Adrenal.combine.list)
Adrenal.combine.list <- lapply(X = Adrenal.combine.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
rm(Adrenal.combine)
gc()
anchors.Adrenal <- FindIntegrationAnchors(object.list = Adrenal.combine.list, reduction = "rpca", 
                                  dims = 1:50)
Adrenal.combine.integrated <- IntegrateData(anchorset = anchors.Adrenal, dims = 1:50)
rm(Adrenal.combine.list)
gc()
DefaultAssay(Adrenal.combine.integrated) <- "integrated"
Adrenal.combine.integrated <- ScaleData(Adrenal.combine.integrated, verbose = FALSE)
Adrenal.combine.integrated <- RunPCA(Adrenal.combine.integrated, verbose = FALSE)
Adrenal.combine.integrated <- RunUMAP(Adrenal.combine.integrated, dims = 1:50)
Adrenal.combine.integrated <- FindNeighbors(Adrenal.combine.integrated, dims = 1:50)
seq<-seq(0.05,1.0,by=0.05)
for (res in seq){
  Adrenal.combine.integrated <- FindClusters(Adrenal.combine.integrated, resolution = res,verbose = TRUE)
}
#################5.3 Determine resolution by clustress
pdf(".../Adrenal/Final_clustree_Adrenal.pdf",width=8,height=16)
clustree(Adrenal.combine.integrated,prefix="integrated_snn_res.")
dev.off()
Idents(Adrenal.combine.integrated)<-Adrenal.combine.integrated$integrated_snn_res.0.2
Idents(Adrenal.combine.integrated)<-factor(Idents(Adrenal.combine.integrated),levels=c(0:(length(unique(Adrenal.combine.integrated$integrated_snn_res.0.2))-1)))
#################5.4 Performe primary visualization 
pdf(".../Adrenal/Final_UMAP_Adrenal.pdf",width=8,height=6)
DimPlot(Adrenal.combine.integrated, reduction = "umap",cols = my36colors,raster=FALSE,label = TRUE)
dev.off()

pdf(".../Adrenal/Final_UMAP_Adrenal_split.pdf",width=14,height=10)
DimPlot(Adrenal.combine.integrated, reduction = "umap",cols = my36colors,raster=FALSE,split.by="orig.ident",ncol=4,label = TRUE)
dev.off()

pdf(".../Adrenal/Final_Featureplot_QC.pdf",width=12,height=5)
FeaturePlot(Adrenal.combine.integrated, features = c("nFeature_RNA", "nCount_RNA"),raster=FALSE,min.cutoff =0,reduction = "umap",cols= c("gray", "red"))
dev.off()

pdf(".../Adrenal/Final_Vln_QC.pdf",width=20,height=6)
VlnPlot(Adrenal.combine.integrated, features = c("nFeature_RNA", "nCount_RNA","percent.mt"),group.by="orig.ident",pt.size=0.01)
dev.off()

cell.prop.sample<-as.data.frame(prop.table(table(Idents(Adrenal.combine.integrated),Adrenal.combine.integrated$orig.ident)))
colnames(cell.prop.sample)<-c("cluster","sample","proportion")
pdf(file=".../Adrenal/Final_cell_proportion_sample.pdf",width=10,height=4)
ggplot(cell.prop.sample,aes(sample,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),axis.text.x = element_text(angle = 45))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=my36colors)
dev.off()
#################5.4 Find marker genes of each cell cluster
plan()
plan(multisession, workers = 4)
plan()
DefaultAssay(Adrenal.combine.integrated)<-"RNA"
Adrenal.combine.integrated.markers <- FindAllMarkers(Adrenal.combine.integrated,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(Adrenal.combine.integrated.markers,".../Adrenal/Adrenal.combine.integrated.markers.csv")
plan("sequential")#change to sequential
options(future.globals.maxSize = 1000 * 1024^3)

saveRDS(Adrenal.combine.integrated,".../Adrenal/Adrenal.combine.integrated.rds")

