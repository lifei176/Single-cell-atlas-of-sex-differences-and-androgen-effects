################# 1. Prepare required softwares and data
#################
#################
#################
#################
#################
library(Seurat)
library(dplyr)
library(BiocParallel)
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', "#b20000",'#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')

Adipose<-readRDS("Adipose.rds")
Adrenal<-readRDS("Adrenal.rds")
Bonemarrow<-readRDS("Bonemarrow.rds")
Brain<-readRDS("Brain.rds")
Colon<-readRDS("Colon.rds")
Heart<-readRDS("Heart.rds")
Intestine<-readRDS("Intestine.rds")
Kidney<-readRDS("Kidney.rds")
Lacrimal<-readRDS("Lacrimal.rds")
Liver<-readRDS("Liver.rds")
Lung<-readRDS("Lung.rds")
Pancreas<-readRDS("Pancreas.rds")
Salivary<-readRDS("Salivary.rds")
Skeletalmuscle<-readRDS("Skeletalmuscle.rds")
Spleen<-readRDS("Spleen.rds")
Stomach<-readRDS("Stomach.rds")
Thymus<-readRDS("Thymus.rds")

################### 2. Define DEGs among differenct conditions across multiple cell types in each tissue
###################
###################
###################
###################
###################

################### 2.1. Adipose
###################
###################
###################
###################
###################
###################(1). Construct object for DEGs definition across multiple cell types and multiple conditions
DefaultAssay(Adipose) <- "RNA"
Adipose$Cell_name_condition <- paste(Adipose$Cell_name, Adipose$Condition, sep = "_")
Idents(Adipose) <- "Cell_name_condition"
table(Adipose$Cell_name_condition)
length(table(Adipose$Cell_name_condition))#107
length(table(Adipose$Cell_name))#28
###################(2). Select cell types for downstream DEG analysis
cell.number.condition.Adipose<-as.data.frame.matrix(table(Adipose$Cell_name,Adipose$Condition))
celltypes_filter_Adipose<-row.names(cell.number.condition.Adipose[which(cell.number.condition.Adipose$FD>2&cell.number.condition.Adipose$FS>2&cell.number.condition.Adipose$MS>2&cell.number.condition.Adipose$MC>2),])

VS <- data.frame(compare=c('MSvsFS','FDvsFS','MCvsMS'),
                 ident1=c('_MS','_FD','_MC'),
                 ident2=c('_FS','_FS','_MS'))
###################(3). Define DEGs among differenct conditions across multiple cell types
library(future)
plan()
plan(multisession, workers = 4)
plan()
DEG_Adipose <- data.frame()
for (j in (1:dim(VS)[1])){
  eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
  for (i in (1:length(celltypes_filter_Adipose))){
    assign(paste0("deg_markers_",celltypes_filter_Adipose[i],"_",format(VS[j,1]),"_Adipose",sep=""),FindMarkers(Adipose, ident.1 = paste0(celltypes_filter_Adipose[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_filter_Adipose[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
    eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_filter_Adipose[i]),"_",format(VS[j,1]),"_Adipose)[1]>0){deg_markers_",format(celltypes_filter_Adipose[i]),"_",format(VS[j,1]),"_Adipose <- subset(deg_markers_",format(celltypes_filter_Adipose[i]),"_",format(VS[j,1]),"_Adipose,p_val_adj<0.05)};
                             if(dim(deg_markers_",format(celltypes_filter_Adipose[i]),"_",format(VS[j,1]),"_Adipose)[1]>0){deg_markers_",format(celltypes_filter_Adipose[i]),"_",format(VS[j,1]),"_Adipose$Gene_symbol <- row.names(deg_markers_",format(celltypes_filter_Adipose[i]),"_",format(VS[j,1]),"_Adipose);
                             deg_markers_",format(celltypes_filter_Adipose[i]),"_",format(VS[j,1]),"_Adipose$Cell_type <- '",format(celltypes_filter_Adipose[i]),"';
                             deg_markers_",format(celltypes_filter_Adipose[i]),"_",format(VS[j,1]),"_Adipose$Change <- ifelse(deg_markers_",format(celltypes_filter_Adipose[i]),"_",format(VS[j,1]),"_Adipose$avg_log2FC > 0 ,'Up','Down');
                             deg_markers_",format(celltypes_filter_Adipose[i]),"_",format(VS[j,1]),"_Adipose$Tissue <- 'Adipose';
                             deg_markers_",format(celltypes_filter_Adipose[i]),"_",format(VS[j,1]),"_Adipose$Comparison <- '",format(VS[j,1]),"';",
                             format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_filter_Adipose[i]),"_",format(VS[j,1]),"_Adipose)};
                             rm(deg_markers_",format(celltypes_filter_Adipose[i]),"_",format(VS[j,1]),"_Adipose)")))
  }
  gc()
  eval(parse(text = paste0("DEG_Adipose <- bind_rows(DEG_Adipose,",format(VS[j,1]),");
                           rm(",format(VS[j,1]),")")))
}

write.csv(DEG_Adipose,'DEG_Adipose.csv')

{

################### 2.2. Adrenal
###################
###################
###################
###################
###################
###################(1). Construct object for DEGs definition across multiple cell types and multiple conditions
DefaultAssay(Adrenal) <- "RNA"
Adrenal$Cell_name_condition <- paste(Adrenal$Cell_name, Adrenal$Condition, sep = "_")
Idents(Adrenal) <- "Cell_name_condition"
table(Adrenal$Cell_name_condition)
length(table(Adrenal$Cell_name_condition))#107
length(table(Adrenal$Cell_name))#28
#################(2). Select cell types for downstream DEG analysis
cell.number.condition.Adrenal<-as.data.frame.matrix(table(Adrenal$Cell_name,Adrenal$Condition))
celltypes_filter_Adrenal<-row.names(cell.number.condition.Adrenal[which(cell.number.condition.Adrenal$FD>2&cell.number.condition.Adrenal$FS>2&cell.number.condition.Adrenal$MS>2&cell.number.condition.Adrenal$MC>2),])

VS <- data.frame(compare=c('MSvsFS','FDvsFS','MCvsMS'),
                 ident1=c('_MS','_FD','_MC'),
                 ident2=c('_FS','_FS','_MS'))
#################(3). Define DEGs among differenct conditions across multiple cell types
library(future)
plan()
plan(multisession, workers = 4)
plan()
DEG_Adrenal <- data.frame()
for (j in (1:dim(VS)[1])){
  eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
  for (i in (1:length(celltypes_filter_Adrenal))){
    assign(paste0("deg_markers_",celltypes_filter_Adrenal[i],"_",format(VS[j,1]),"_Adrenal",sep=""),FindMarkers(Adrenal, ident.1 = paste0(celltypes_filter_Adrenal[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_filter_Adrenal[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
    eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_filter_Adrenal[i]),"_",format(VS[j,1]),"_Adrenal)[1]>0){deg_markers_",format(celltypes_filter_Adrenal[i]),"_",format(VS[j,1]),"_Adrenal <- subset(deg_markers_",format(celltypes_filter_Adrenal[i]),"_",format(VS[j,1]),"_Adrenal,p_val_adj<0.05)};
                             if(dim(deg_markers_",format(celltypes_filter_Adrenal[i]),"_",format(VS[j,1]),"_Adrenal)[1]>0){deg_markers_",format(celltypes_filter_Adrenal[i]),"_",format(VS[j,1]),"_Adrenal$Gene_symbol <- row.names(deg_markers_",format(celltypes_filter_Adrenal[i]),"_",format(VS[j,1]),"_Adrenal);
                             deg_markers_",format(celltypes_filter_Adrenal[i]),"_",format(VS[j,1]),"_Adrenal$Cell_type <- '",format(celltypes_filter_Adrenal[i]),"';
                             deg_markers_",format(celltypes_filter_Adrenal[i]),"_",format(VS[j,1]),"_Adrenal$Change <- ifelse(deg_markers_",format(celltypes_filter_Adrenal[i]),"_",format(VS[j,1]),"_Adrenal$avg_log2FC > 0 ,'Up','Down');
                             deg_markers_",format(celltypes_filter_Adrenal[i]),"_",format(VS[j,1]),"_Adrenal$Tissue <- 'Adrenal';
                             deg_markers_",format(celltypes_filter_Adrenal[i]),"_",format(VS[j,1]),"_Adrenal$Comparison <- '",format(VS[j,1]),"';",
                             format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_filter_Adrenal[i]),"_",format(VS[j,1]),"_Adrenal)};
                             rm(deg_markers_",format(celltypes_filter_Adrenal[i]),"_",format(VS[j,1]),"_Adrenal)")))
                             }
  gc()
  eval(parse(text = paste0("DEG_Adrenal <- bind_rows(DEG_Adrenal,",format(VS[j,1]),");
                           rm(",format(VS[j,1]),")")))
  }

write.csv(DEG_Adrenal,'DEG_Adrenal.csv')

################### 2.3. Brain
###################
###################
###################
###################
###################
###################(1). Construct object for DEGs definition across multiple cell types and multiple conditions
DefaultAssay(Brain) <- "RNA"
Brain$Cell_name_condition <- paste(Brain$Cell_name, Brain$Condition, sep = "_")
Idents(Brain) <- "Cell_name_condition"
table(Brain$Cell_name_condition)
length(table(Brain$Cell_name_condition))#107
length(table(Brain$Cell_name))#28
#################(2). Select cell types for downstream DEG analysis
cell.number.condition.Brain<-as.data.frame.matrix(table(Brain$Cell_name,Brain$Condition))
celltypes_filter_Brain<-row.names(cell.number.condition.Brain[which(cell.number.condition.Brain$FD>2&cell.number.condition.Brain$FS>2&cell.number.condition.Brain$MS>2&cell.number.condition.Brain$MC>2),])

VS <- data.frame(compare=c('MSvsFS','FDvsFS','MCvsMS'),
                 ident1=c('_MS','_FD','_MC'),
                 ident2=c('_FS','_FS','_MS'))
#################(3). Define DEGs among differenct conditions across multiple cell types
library(future)
plan()
plan(multisession, workers = 4)
plan()
DEG_Brain <- data.frame()
for (j in (1:dim(VS)[1])){
  eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
  for (i in (1:length(celltypes_filter_Brain))){
    assign(paste0("deg_markers_",celltypes_filter_Brain[i],"_",format(VS[j,1]),"_Brain",sep=""),FindMarkers(Brain, ident.1 = paste0(celltypes_filter_Brain[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_filter_Brain[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
    eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_filter_Brain[i]),"_",format(VS[j,1]),"_Brain)[1]>0){deg_markers_",format(celltypes_filter_Brain[i]),"_",format(VS[j,1]),"_Brain <- subset(deg_markers_",format(celltypes_filter_Brain[i]),"_",format(VS[j,1]),"_Brain,p_val_adj<0.05)};
                             if(dim(deg_markers_",format(celltypes_filter_Brain[i]),"_",format(VS[j,1]),"_Brain)[1]>0){deg_markers_",format(celltypes_filter_Brain[i]),"_",format(VS[j,1]),"_Brain$Gene_symbol <- row.names(deg_markers_",format(celltypes_filter_Brain[i]),"_",format(VS[j,1]),"_Brain);
                             deg_markers_",format(celltypes_filter_Brain[i]),"_",format(VS[j,1]),"_Brain$Cell_type <- '",format(celltypes_filter_Brain[i]),"';
                             deg_markers_",format(celltypes_filter_Brain[i]),"_",format(VS[j,1]),"_Brain$Change <- ifelse(deg_markers_",format(celltypes_filter_Brain[i]),"_",format(VS[j,1]),"_Brain$avg_log2FC > 0 ,'Up','Down');
                             deg_markers_",format(celltypes_filter_Brain[i]),"_",format(VS[j,1]),"_Brain$Tissue <- 'Brain';
                             deg_markers_",format(celltypes_filter_Brain[i]),"_",format(VS[j,1]),"_Brain$Comparison <- '",format(VS[j,1]),"';",
                             format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_filter_Brain[i]),"_",format(VS[j,1]),"_Brain)};
                             rm(deg_markers_",format(celltypes_filter_Brain[i]),"_",format(VS[j,1]),"_Brain)")))
                             }
  gc()
  eval(parse(text = paste0("DEG_Brain <- bind_rows(DEG_Brain,",format(VS[j,1]),");
                           rm(",format(VS[j,1]),")")))
  }

write.csv(DEG_Brain,'DEG_Brain.csv')

################### 2.4. Bonemarrow
###################
###################
###################
###################
###################
###################(1). Construct object for DEGs definition across multiple cell types and multiple conditions
DefaultAssay(Bonemarrow) <- "RNA"
Bonemarrow$Cell_name_condition <- paste(Bonemarrow$Cell_name, Bonemarrow$Condition, sep = "_")
Idents(Bonemarrow) <- "Cell_name_condition"
table(Bonemarrow$Cell_name_condition)
length(table(Bonemarrow$Cell_name_condition))#107
length(table(Bonemarrow$Cell_name))#28
#################(2). Select cell types for downstream DEG analysis
cell.number.condition.Bonemarrow<-as.data.frame.matrix(table(Bonemarrow$Cell_name,Bonemarrow$Condition))
celltypes_filter_Bonemarrow<-row.names(cell.number.condition.Bonemarrow[which(cell.number.condition.Bonemarrow$FD>2&cell.number.condition.Bonemarrow$FS>2&cell.number.condition.Bonemarrow$MS>2&cell.number.condition.Bonemarrow$MC>2),])

VS <- data.frame(compare=c('MSvsFS','FDvsFS','MCvsMS'),
                 ident1=c('_MS','_FD','_MC'),
                 ident2=c('_FS','_FS','_MS'))
#################(3). Define DEGs among differenct conditions across multiple cell types
library(future)
plan()
plan(multisession, workers = 4)
plan()
DEG_Bonemarrow <- data.frame()
for (j in (1:dim(VS)[1])){
  eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
  for (i in (1:length(celltypes_filter_Bonemarrow))){
    assign(paste0("deg_markers_",celltypes_filter_Bonemarrow[i],"_",format(VS[j,1]),"_Bonemarrow",sep=""),FindMarkers(Bonemarrow, ident.1 = paste0(celltypes_filter_Bonemarrow[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_filter_Bonemarrow[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
    eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_filter_Bonemarrow[i]),"_",format(VS[j,1]),"_Bonemarrow)[1]>0){deg_markers_",format(celltypes_filter_Bonemarrow[i]),"_",format(VS[j,1]),"_Bonemarrow <- subset(deg_markers_",format(celltypes_filter_Bonemarrow[i]),"_",format(VS[j,1]),"_Bonemarrow,p_val_adj<0.05)};
                             if(dim(deg_markers_",format(celltypes_filter_Bonemarrow[i]),"_",format(VS[j,1]),"_Bonemarrow)[1]>0){deg_markers_",format(celltypes_filter_Bonemarrow[i]),"_",format(VS[j,1]),"_Bonemarrow$Gene_symbol <- row.names(deg_markers_",format(celltypes_filter_Bonemarrow[i]),"_",format(VS[j,1]),"_Bonemarrow);
                             deg_markers_",format(celltypes_filter_Bonemarrow[i]),"_",format(VS[j,1]),"_Bonemarrow$Cell_type <- '",format(celltypes_filter_Bonemarrow[i]),"';
                             deg_markers_",format(celltypes_filter_Bonemarrow[i]),"_",format(VS[j,1]),"_Bonemarrow$Change <- ifelse(deg_markers_",format(celltypes_filter_Bonemarrow[i]),"_",format(VS[j,1]),"_Bonemarrow$avg_log2FC > 0 ,'Up','Down');
                             deg_markers_",format(celltypes_filter_Bonemarrow[i]),"_",format(VS[j,1]),"_Bonemarrow$Tissue <- 'Bonemarrow';
                             deg_markers_",format(celltypes_filter_Bonemarrow[i]),"_",format(VS[j,1]),"_Bonemarrow$Comparison <- '",format(VS[j,1]),"';",
                             format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_filter_Bonemarrow[i]),"_",format(VS[j,1]),"_Bonemarrow)};
                             rm(deg_markers_",format(celltypes_filter_Bonemarrow[i]),"_",format(VS[j,1]),"_Bonemarrow)")))
                             }
  gc()
  eval(parse(text = paste0("DEG_Bonemarrow <- bind_rows(DEG_Bonemarrow,",format(VS[j,1]),");
                           rm(",format(VS[j,1]),")")))
  }

write.csv(DEG_Bonemarrow,'DEG_Bonemarrow.csv')

################### 2.5. Colon
###################
###################
###################
###################
###################
###################(1). Construct object for DEGs definition across multiple cell types and multiple conditions
DefaultAssay(Colon) <- "RNA"
Colon$Cell_name_condition <- paste(Colon$Cell_name, Colon$Condition, sep = "_")
Idents(Colon) <- "Cell_name_condition"
table(Colon$Cell_name_condition)
length(table(Colon$Cell_name_condition))#107
length(table(Colon$Cell_name))#28
#################(2). Select cell types for downstream DEG analysis
cell.number.condition.Colon<-as.data.frame.matrix(table(Colon$Cell_name,Colon$Condition))
celltypes_filter_Colon<-row.names(cell.number.condition.Colon[which(cell.number.condition.Colon$FD>2&cell.number.condition.Colon$FS>2&cell.number.condition.Colon$MS>2&cell.number.condition.Colon$MC>2),])

VS <- data.frame(compare=c('MSvsFS','FDvsFS','MCvsMS'),
                 ident1=c('_MS','_FD','_MC'),
                 ident2=c('_FS','_FS','_MS'))
#################(3). Define DEGs among differenct conditions across multiple cell types
library(future)
plan()
plan(multisession, workers = 4)
plan()
DEG_Colon <- data.frame()
for (j in (1:dim(VS)[1])){
  eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
  for (i in (1:length(celltypes_filter_Colon))){
    assign(paste0("deg_markers_",celltypes_filter_Colon[i],"_",format(VS[j,1]),"_Colon",sep=""),FindMarkers(Colon, ident.1 = paste0(celltypes_filter_Colon[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_filter_Colon[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
    eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_filter_Colon[i]),"_",format(VS[j,1]),"_Colon)[1]>0){deg_markers_",format(celltypes_filter_Colon[i]),"_",format(VS[j,1]),"_Colon <- subset(deg_markers_",format(celltypes_filter_Colon[i]),"_",format(VS[j,1]),"_Colon,p_val_adj<0.05)};
                             if(dim(deg_markers_",format(celltypes_filter_Colon[i]),"_",format(VS[j,1]),"_Colon)[1]>0){deg_markers_",format(celltypes_filter_Colon[i]),"_",format(VS[j,1]),"_Colon$Gene_symbol <- row.names(deg_markers_",format(celltypes_filter_Colon[i]),"_",format(VS[j,1]),"_Colon);
                             deg_markers_",format(celltypes_filter_Colon[i]),"_",format(VS[j,1]),"_Colon$Cell_type <- '",format(celltypes_filter_Colon[i]),"';
                             deg_markers_",format(celltypes_filter_Colon[i]),"_",format(VS[j,1]),"_Colon$Change <- ifelse(deg_markers_",format(celltypes_filter_Colon[i]),"_",format(VS[j,1]),"_Colon$avg_log2FC > 0 ,'Up','Down');
                             deg_markers_",format(celltypes_filter_Colon[i]),"_",format(VS[j,1]),"_Colon$Tissue <- 'Colon';
                             deg_markers_",format(celltypes_filter_Colon[i]),"_",format(VS[j,1]),"_Colon$Comparison <- '",format(VS[j,1]),"';",
                             format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_filter_Colon[i]),"_",format(VS[j,1]),"_Colon)};
                             rm(deg_markers_",format(celltypes_filter_Colon[i]),"_",format(VS[j,1]),"_Colon)")))
  }
  gc()
  eval(parse(text = paste0("DEG_Colon <- bind_rows(DEG_Colon,",format(VS[j,1]),");
                           rm(",format(VS[j,1]),")")))
}

write.csv(DEG_Colon,'DEG_Colon.csv')

################### 2.6. Heart
###################
###################
###################
###################
###################
###################(1). Construct object for DEGs definition across multiple cell types and multiple conditions
DefaultAssay(Heart) <- "RNA"
Heart$Cell_name_condition <- paste(Heart$Cell_name, Heart$Condition, sep = "_")
Idents(Heart) <- "Cell_name_condition"
table(Heart$Cell_name_condition)
length(table(Heart$Cell_name_condition))#107
length(table(Heart$Cell_name))#28
#################(2). Select cell types for downstream DEG analysis
cell.number.condition.Heart<-as.data.frame.matrix(table(Heart$Cell_name,Heart$Condition))
celltypes_filter_Heart<-row.names(cell.number.condition.Heart[which(cell.number.condition.Heart$FD>2&cell.number.condition.Heart$FS>2&cell.number.condition.Heart$MS>2&cell.number.condition.Heart$MC>2),])

VS <- data.frame(compare=c('MSvsFS','FDvsFS','MCvsMS'),
                 ident1=c('_MS','_FD','_MC'),
                 ident2=c('_FS','_FS','_MS'))
#################(3). Define DEGs among differenct conditions across multiple cell types
library(future)
plan()
plan(multisession, workers = 4)
plan()
DEG_Heart <- data.frame()
for (j in (1:dim(VS)[1])){
  eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
  for (i in (1:length(celltypes_filter_Heart))){
    assign(paste0("deg_markers_",celltypes_filter_Heart[i],"_",format(VS[j,1]),"_Heart",sep=""),FindMarkers(Heart, ident.1 = paste0(celltypes_filter_Heart[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_filter_Heart[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
    eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_filter_Heart[i]),"_",format(VS[j,1]),"_Heart)[1]>0){deg_markers_",format(celltypes_filter_Heart[i]),"_",format(VS[j,1]),"_Heart <- subset(deg_markers_",format(celltypes_filter_Heart[i]),"_",format(VS[j,1]),"_Heart,p_val_adj<0.05)};
                             if(dim(deg_markers_",format(celltypes_filter_Heart[i]),"_",format(VS[j,1]),"_Heart)[1]>0){deg_markers_",format(celltypes_filter_Heart[i]),"_",format(VS[j,1]),"_Heart$Gene_symbol <- row.names(deg_markers_",format(celltypes_filter_Heart[i]),"_",format(VS[j,1]),"_Heart);
                             deg_markers_",format(celltypes_filter_Heart[i]),"_",format(VS[j,1]),"_Heart$Cell_type <- '",format(celltypes_filter_Heart[i]),"';
                             deg_markers_",format(celltypes_filter_Heart[i]),"_",format(VS[j,1]),"_Heart$Change <- ifelse(deg_markers_",format(celltypes_filter_Heart[i]),"_",format(VS[j,1]),"_Heart$avg_log2FC > 0 ,'Up','Down');
                             deg_markers_",format(celltypes_filter_Heart[i]),"_",format(VS[j,1]),"_Heart$Tissue <- 'Heart';
                             deg_markers_",format(celltypes_filter_Heart[i]),"_",format(VS[j,1]),"_Heart$Comparison <- '",format(VS[j,1]),"';",
                             format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_filter_Heart[i]),"_",format(VS[j,1]),"_Heart)};
                             rm(deg_markers_",format(celltypes_filter_Heart[i]),"_",format(VS[j,1]),"_Heart)")))
                             }
  gc()
  eval(parse(text = paste0("DEG_Heart <- bind_rows(DEG_Heart,",format(VS[j,1]),");
                           rm(",format(VS[j,1]),")")))
  }

write.csv(DEG_Heart,'DEG_Heart.csv')



################### 2.7. Intestine
###################
###################
###################
###################
###################
###################(1). Construct object for DEGs definition across multiple cell types and multiple conditions
DefaultAssay(Intestine) <- "RNA"
Intestine$Cell_name_condition <- paste(Intestine$Cell_name, Intestine$Condition, sep = "_")
Idents(Intestine) <- "Cell_name_condition"
table(Intestine$Cell_name_condition)
length(table(Intestine$Cell_name_condition))#107
length(table(Intestine$Cell_name))#28
#################(2). Select cell types for downstream DEG analysis
cell.number.condition.Intestine<-as.data.frame.matrix(table(Intestine$Cell_name,Intestine$Condition))
celltypes_filter_Intestine<-row.names(cell.number.condition.Intestine[which(cell.number.condition.Intestine$FD>2&cell.number.condition.Intestine$FS>2&cell.number.condition.Intestine$MS>2&cell.number.condition.Intestine$MC>2),])

VS <- data.frame(compare=c('MSvsFS','FDvsFS','MCvsMS'),
                 ident1=c('_MS','_FD','_MC'),
                 ident2=c('_FS','_FS','_MS'))
#################(3). Define DEGs among differenct conditions across multiple cell types
library(future)
plan()
plan(multisession, workers = 4)
plan()
DEG_Intestine <- data.frame()
for (j in (1:dim(VS)[1])){
  eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
  for (i in (1:length(celltypes_filter_Intestine))){
    assign(paste0("deg_markers_",celltypes_filter_Intestine[i],"_",format(VS[j,1]),"_Intestine",sep=""),FindMarkers(Intestine, ident.1 = paste0(celltypes_filter_Intestine[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_filter_Intestine[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
    eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_filter_Intestine[i]),"_",format(VS[j,1]),"_Intestine)[1]>0){deg_markers_",format(celltypes_filter_Intestine[i]),"_",format(VS[j,1]),"_Intestine <- subset(deg_markers_",format(celltypes_filter_Intestine[i]),"_",format(VS[j,1]),"_Intestine,p_val_adj<0.05)};
                             if(dim(deg_markers_",format(celltypes_filter_Intestine[i]),"_",format(VS[j,1]),"_Intestine)[1]>0){deg_markers_",format(celltypes_filter_Intestine[i]),"_",format(VS[j,1]),"_Intestine$Gene_symbol <- row.names(deg_markers_",format(celltypes_filter_Intestine[i]),"_",format(VS[j,1]),"_Intestine);
                             deg_markers_",format(celltypes_filter_Intestine[i]),"_",format(VS[j,1]),"_Intestine$Cell_type <- '",format(celltypes_filter_Intestine[i]),"';
                             deg_markers_",format(celltypes_filter_Intestine[i]),"_",format(VS[j,1]),"_Intestine$Change <- ifelse(deg_markers_",format(celltypes_filter_Intestine[i]),"_",format(VS[j,1]),"_Intestine$avg_log2FC > 0 ,'Up','Down');
                             deg_markers_",format(celltypes_filter_Intestine[i]),"_",format(VS[j,1]),"_Intestine$Tissue <- 'Intestine';
                             deg_markers_",format(celltypes_filter_Intestine[i]),"_",format(VS[j,1]),"_Intestine$Comparison <- '",format(VS[j,1]),"';",
                             format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_filter_Intestine[i]),"_",format(VS[j,1]),"_Intestine)};
                             rm(deg_markers_",format(celltypes_filter_Intestine[i]),"_",format(VS[j,1]),"_Intestine)")))
                             }
  gc()
  eval(parse(text = paste0("DEG_Intestine <- bind_rows(DEG_Intestine,",format(VS[j,1]),");
                           rm(",format(VS[j,1]),")")))
  }

write.csv(DEG_Intestine,'DEG_Intestine.csv')

################### 2.8. Kidney
###################
###################
###################
###################
###################
###################(1). Construct object for DEGs definition across multiple cell types and multiple conditions
DefaultAssay(Kidney) <- "RNA"
Kidney$Cell_name_condition <- paste(Kidney$Cell_name, Kidney$Condition, sep = "_")
Idents(Kidney) <- "Cell_name_condition"
table(Kidney$Cell_name_condition)
length(table(Kidney$Cell_name_condition))#107
length(table(Kidney$Cell_name))#28
#################(2). Select cell types for downstream DEG analysis
cell.number.condition.Kidney<-as.data.frame.matrix(table(Kidney$Cell_name,Kidney$Condition))
celltypes_filter_Kidney<-row.names(cell.number.condition.Kidney[which(cell.number.condition.Kidney$FD>2&cell.number.condition.Kidney$FS>2&cell.number.condition.Kidney$MS>2&cell.number.condition.Kidney$MC>2),])

VS <- data.frame(compare=c('MSvsFS','FDvsFS','MCvsMS'),
                 ident1=c('_MS','_FD','_MC'),
                 ident2=c('_FS','_FS','_MS'))
#################(3). Define DEGs among differenct conditions across multiple cell types
library(future)
plan()
plan(multisession, workers = 4)
plan()
DEG_Kidney <- data.frame()
for (j in (1:dim(VS)[1])){
  eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
  for (i in (1:length(celltypes_filter_Kidney))){
    assign(paste0("deg_markers_",celltypes_filter_Kidney[i],"_",format(VS[j,1]),"_Kidney",sep=""),FindMarkers(Kidney, ident.1 = paste0(celltypes_filter_Kidney[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_filter_Kidney[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
    eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_filter_Kidney[i]),"_",format(VS[j,1]),"_Kidney)[1]>0){deg_markers_",format(celltypes_filter_Kidney[i]),"_",format(VS[j,1]),"_Kidney <- subset(deg_markers_",format(celltypes_filter_Kidney[i]),"_",format(VS[j,1]),"_Kidney,p_val_adj<0.05)};
                             if(dim(deg_markers_",format(celltypes_filter_Kidney[i]),"_",format(VS[j,1]),"_Kidney)[1]>0){deg_markers_",format(celltypes_filter_Kidney[i]),"_",format(VS[j,1]),"_Kidney$Gene_symbol <- row.names(deg_markers_",format(celltypes_filter_Kidney[i]),"_",format(VS[j,1]),"_Kidney);
                             deg_markers_",format(celltypes_filter_Kidney[i]),"_",format(VS[j,1]),"_Kidney$Cell_type <- '",format(celltypes_filter_Kidney[i]),"';
                             deg_markers_",format(celltypes_filter_Kidney[i]),"_",format(VS[j,1]),"_Kidney$Change <- ifelse(deg_markers_",format(celltypes_filter_Kidney[i]),"_",format(VS[j,1]),"_Kidney$avg_log2FC > 0 ,'Up','Down');
                             deg_markers_",format(celltypes_filter_Kidney[i]),"_",format(VS[j,1]),"_Kidney$Tissue <- 'Kidney';
                             deg_markers_",format(celltypes_filter_Kidney[i]),"_",format(VS[j,1]),"_Kidney$Comparison <- '",format(VS[j,1]),"';",
                             format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_filter_Kidney[i]),"_",format(VS[j,1]),"_Kidney)};
                             rm(deg_markers_",format(celltypes_filter_Kidney[i]),"_",format(VS[j,1]),"_Kidney)")))
                             }
  gc()
  eval(parse(text = paste0("DEG_Kidney <- bind_rows(DEG_Kidney,",format(VS[j,1]),");
                           rm(",format(VS[j,1]),")")))
  }

write.csv(DEG_Kidney,'DEG_Kidney.csv')

################### 2.9. Lacrimal
###################
###################
###################
###################
###################
###################(1). Construct object for DEGs definition across multiple cell types and multiple conditions
DefaultAssay(Lacrimal) <- "RNA"
Lacrimal$Cell_name_condition <- paste(Lacrimal$Cell_name, Lacrimal$Condition, sep = "_")
Idents(Lacrimal) <- "Cell_name_condition"
table(Lacrimal$Cell_name_condition)
length(table(Lacrimal$Cell_name_condition))#107
length(table(Lacrimal$Cell_name))#28
#################(2). Select cell types for downstream DEG analysis
cell.number.condition.Lacrimal<-as.data.frame.matrix(table(Lacrimal$Cell_name,Lacrimal$Condition))
celltypes_filter_Lacrimal<-row.names(cell.number.condition.Lacrimal[which(cell.number.condition.Lacrimal$FD>2&cell.number.condition.Lacrimal$FS>2&cell.number.condition.Lacrimal$MS>2&cell.number.condition.Lacrimal$MC>2),])

VS <- data.frame(compare=c('MSvsFS','FDvsFS','MCvsMS'),
                 ident1=c('_MS','_FD','_MC'),
                 ident2=c('_FS','_FS','_MS'))
#################(3). Define DEGs among differenct conditions across multiple cell types
library(future)
plan()
plan(multisession, workers = 4)
plan()
DEG_Lacrimal <- data.frame()
for (j in (1:dim(VS)[1])){
  eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
  for (i in (1:length(celltypes_filter_Lacrimal))){
    assign(paste0("deg_markers_",celltypes_filter_Lacrimal[i],"_",format(VS[j,1]),"_Lacrimal",sep=""),FindMarkers(Lacrimal, ident.1 = paste0(celltypes_filter_Lacrimal[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_filter_Lacrimal[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
    eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_filter_Lacrimal[i]),"_",format(VS[j,1]),"_Lacrimal)[1]>0){deg_markers_",format(celltypes_filter_Lacrimal[i]),"_",format(VS[j,1]),"_Lacrimal <- subset(deg_markers_",format(celltypes_filter_Lacrimal[i]),"_",format(VS[j,1]),"_Lacrimal,p_val_adj<0.05)};
                             if(dim(deg_markers_",format(celltypes_filter_Lacrimal[i]),"_",format(VS[j,1]),"_Lacrimal)[1]>0){deg_markers_",format(celltypes_filter_Lacrimal[i]),"_",format(VS[j,1]),"_Lacrimal$Gene_symbol <- row.names(deg_markers_",format(celltypes_filter_Lacrimal[i]),"_",format(VS[j,1]),"_Lacrimal);
                             deg_markers_",format(celltypes_filter_Lacrimal[i]),"_",format(VS[j,1]),"_Lacrimal$Cell_type <- '",format(celltypes_filter_Lacrimal[i]),"';
                             deg_markers_",format(celltypes_filter_Lacrimal[i]),"_",format(VS[j,1]),"_Lacrimal$Change <- ifelse(deg_markers_",format(celltypes_filter_Lacrimal[i]),"_",format(VS[j,1]),"_Lacrimal$avg_log2FC > 0 ,'Up','Down');
                             deg_markers_",format(celltypes_filter_Lacrimal[i]),"_",format(VS[j,1]),"_Lacrimal$Tissue <- 'Lacrimal';
                             deg_markers_",format(celltypes_filter_Lacrimal[i]),"_",format(VS[j,1]),"_Lacrimal$Comparison <- '",format(VS[j,1]),"';",
                             format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_filter_Lacrimal[i]),"_",format(VS[j,1]),"_Lacrimal)};
                             rm(deg_markers_",format(celltypes_filter_Lacrimal[i]),"_",format(VS[j,1]),"_Lacrimal)")))
                             }
  gc()
  eval(parse(text = paste0("DEG_Lacrimal <- bind_rows(DEG_Lacrimal,",format(VS[j,1]),");
                           rm(",format(VS[j,1]),")")))
  }

write.csv(DEG_Lacrimal,'DEG_Lacrimal.csv')

################### 2.10. Liver
###################
###################
###################
###################
###################
###################(1). Construct object for DEGs definition across multiple cell types and multiple conditions
DefaultAssay(Liver) <- "RNA"
Liver$Cell_name_condition <- paste(Liver$Cell_name, Liver$Condition, sep = "_")
Idents(Liver) <- "Cell_name_condition"
table(Liver$Cell_name_condition)
length(table(Liver$Cell_name_condition))#107
length(table(Liver$Cell_name))#28
#################(2). Select cell types for downstream DEG analysis
cell.number.condition.Liver<-as.data.frame.matrix(table(Liver$Cell_name,Liver$Condition))
celltypes_filter_Liver<-row.names(cell.number.condition.Liver[which(cell.number.condition.Liver$FD>2&cell.number.condition.Liver$FS>2&cell.number.condition.Liver$MS>2&cell.number.condition.Liver$MC>2),])

VS <- data.frame(compare=c('MSvsFS','FDvsFS','MCvsMS'),
                 ident1=c('_MS','_FD','_MC'),
                 ident2=c('_FS','_FS','_MS'))
#################(3). Define DEGs among differenct conditions across multiple cell types
library(future)
plan()
plan(multisession, workers = 4)
plan()
DEG_Liver <- data.frame()
for (j in (1:dim(VS)[1])){
  eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
  for (i in (1:length(celltypes_filter_Liver))){
    assign(paste0("deg_markers_",celltypes_filter_Liver[i],"_",format(VS[j,1]),"_Liver",sep=""),FindMarkers(Liver, ident.1 = paste0(celltypes_filter_Liver[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_filter_Liver[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
    eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_filter_Liver[i]),"_",format(VS[j,1]),"_Liver)[1]>0){deg_markers_",format(celltypes_filter_Liver[i]),"_",format(VS[j,1]),"_Liver <- subset(deg_markers_",format(celltypes_filter_Liver[i]),"_",format(VS[j,1]),"_Liver,p_val_adj<0.05)};
                             if(dim(deg_markers_",format(celltypes_filter_Liver[i]),"_",format(VS[j,1]),"_Liver)[1]>0){deg_markers_",format(celltypes_filter_Liver[i]),"_",format(VS[j,1]),"_Liver$Gene_symbol <- row.names(deg_markers_",format(celltypes_filter_Liver[i]),"_",format(VS[j,1]),"_Liver);
                             deg_markers_",format(celltypes_filter_Liver[i]),"_",format(VS[j,1]),"_Liver$Cell_type <- '",format(celltypes_filter_Liver[i]),"';
                             deg_markers_",format(celltypes_filter_Liver[i]),"_",format(VS[j,1]),"_Liver$Change <- ifelse(deg_markers_",format(celltypes_filter_Liver[i]),"_",format(VS[j,1]),"_Liver$avg_log2FC > 0 ,'Up','Down');
                             deg_markers_",format(celltypes_filter_Liver[i]),"_",format(VS[j,1]),"_Liver$Tissue <- 'Liver';
                             deg_markers_",format(celltypes_filter_Liver[i]),"_",format(VS[j,1]),"_Liver$Comparison <- '",format(VS[j,1]),"';",
                             format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_filter_Liver[i]),"_",format(VS[j,1]),"_Liver)};
                             rm(deg_markers_",format(celltypes_filter_Liver[i]),"_",format(VS[j,1]),"_Liver)")))
                             }
  gc()
  eval(parse(text = paste0("DEG_Liver <- bind_rows(DEG_Liver,",format(VS[j,1]),");
                           rm(",format(VS[j,1]),")")))
  }

write.csv(DEG_Liver,'DEG_Liver.csv')


################### 2.11. Lung
###################
###################
###################
###################
###################
###################(1). Construct object for DEGs definition across multiple cell types and multiple conditions
DefaultAssay(Lung) <- "RNA"
Lung$Cell_name_condition <- paste(Lung$Cell_name, Lung$Condition, sep = "_")
Idents(Lung) <- "Cell_name_condition"
table(Lung$Cell_name_condition)
length(table(Lung$Cell_name_condition))#107
length(table(Lung$Cell_name))#28
#################(2). Select cell types for downstream DEG analysis
cell.number.condition.Lung<-as.data.frame.matrix(table(Lung$Cell_name,Lung$Condition))
celltypes_filter_Lung<-row.names(cell.number.condition.Lung[which(cell.number.condition.Lung$FD>2&cell.number.condition.Lung$FS>2&cell.number.condition.Lung$MS>2&cell.number.condition.Lung$MC>2),])

VS <- data.frame(compare=c('MSvsFS','FDvsFS','MCvsMS'),
                 ident1=c('_MS','_FD','_MC'),
                 ident2=c('_FS','_FS','_MS'))
#################(3). Define DEGs among differenct conditions across multiple cell types
library(future)
plan()
plan(multisession, workers = 4)
plan()
DEG_Lung <- data.frame()
for (j in (1:dim(VS)[1])){
  eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
  for (i in (1:length(celltypes_filter_Lung))){
    assign(paste0("deg_markers_",celltypes_filter_Lung[i],"_",format(VS[j,1]),"_Lung",sep=""),FindMarkers(Lung, ident.1 = paste0(celltypes_filter_Lung[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_filter_Lung[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
    eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_filter_Lung[i]),"_",format(VS[j,1]),"_Lung)[1]>0){deg_markers_",format(celltypes_filter_Lung[i]),"_",format(VS[j,1]),"_Lung <- subset(deg_markers_",format(celltypes_filter_Lung[i]),"_",format(VS[j,1]),"_Lung,p_val_adj<0.05)};
                             if(dim(deg_markers_",format(celltypes_filter_Lung[i]),"_",format(VS[j,1]),"_Lung)[1]>0){deg_markers_",format(celltypes_filter_Lung[i]),"_",format(VS[j,1]),"_Lung$Gene_symbol <- row.names(deg_markers_",format(celltypes_filter_Lung[i]),"_",format(VS[j,1]),"_Lung);
                             deg_markers_",format(celltypes_filter_Lung[i]),"_",format(VS[j,1]),"_Lung$Cell_type <- '",format(celltypes_filter_Lung[i]),"';
                             deg_markers_",format(celltypes_filter_Lung[i]),"_",format(VS[j,1]),"_Lung$Change <- ifelse(deg_markers_",format(celltypes_filter_Lung[i]),"_",format(VS[j,1]),"_Lung$avg_log2FC > 0 ,'Up','Down');
                             deg_markers_",format(celltypes_filter_Lung[i]),"_",format(VS[j,1]),"_Lung$Tissue <- 'Lung';
                             deg_markers_",format(celltypes_filter_Lung[i]),"_",format(VS[j,1]),"_Lung$Comparison <- '",format(VS[j,1]),"';",
                             format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_filter_Lung[i]),"_",format(VS[j,1]),"_Lung)};
                             rm(deg_markers_",format(celltypes_filter_Lung[i]),"_",format(VS[j,1]),"_Lung)")))
                             }
  gc()
  eval(parse(text = paste0("DEG_Lung <- bind_rows(DEG_Lung,",format(VS[j,1]),");
                           rm(",format(VS[j,1]),")")))
  }

write.csv(DEG_Lung,'DEG_Lung.csv')

################### 2.12. Pancreas
###################
###################
###################
###################
###################
###################(1). Construct object for DEGs definition across multiple cell types and multiple conditions
DefaultAssay(Pancreas) <- "RNA"
Pancreas$Cell_name_condition <- paste(Pancreas$Cell_name, Pancreas$Condition, sep = "_")
Idents(Pancreas) <- "Cell_name_condition"
table(Pancreas$Cell_name_condition)
length(table(Pancreas$Cell_name_condition))#107
length(table(Pancreas$Cell_name))#28
#################(2). Select cell types for downstream DEG analysis
cell.number.condition.Pancreas<-as.data.frame.matrix(table(Pancreas$Cell_name,Pancreas$Condition))
celltypes_filter_Pancreas<-row.names(cell.number.condition.Pancreas[which(cell.number.condition.Pancreas$FD>2&cell.number.condition.Pancreas$FS>2&cell.number.condition.Pancreas$MS>2&cell.number.condition.Pancreas$MC>2),])

VS <- data.frame(compare=c('MSvsFS','FDvsFS','MCvsMS'),
                 ident1=c('_MS','_FD','_MC'),
                 ident2=c('_FS','_FS','_MS'))
#################(3). Define DEGs among differenct conditions across multiple cell types
library(future)
plan()
plan(multisession, workers = 4)
plan()
DEG_Pancreas <- data.frame()
for (j in (1:dim(VS)[1])){
  eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
  for (i in (1:length(celltypes_filter_Pancreas))){
    assign(paste0("deg_markers_",celltypes_filter_Pancreas[i],"_",format(VS[j,1]),"_Pancreas",sep=""),FindMarkers(Pancreas, ident.1 = paste0(celltypes_filter_Pancreas[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_filter_Pancreas[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
    eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_filter_Pancreas[i]),"_",format(VS[j,1]),"_Pancreas)[1]>0){deg_markers_",format(celltypes_filter_Pancreas[i]),"_",format(VS[j,1]),"_Pancreas <- subset(deg_markers_",format(celltypes_filter_Pancreas[i]),"_",format(VS[j,1]),"_Pancreas,p_val_adj<0.05)};
                             if(dim(deg_markers_",format(celltypes_filter_Pancreas[i]),"_",format(VS[j,1]),"_Pancreas)[1]>0){deg_markers_",format(celltypes_filter_Pancreas[i]),"_",format(VS[j,1]),"_Pancreas$Gene_symbol <- row.names(deg_markers_",format(celltypes_filter_Pancreas[i]),"_",format(VS[j,1]),"_Pancreas);
                             deg_markers_",format(celltypes_filter_Pancreas[i]),"_",format(VS[j,1]),"_Pancreas$Cell_type <- '",format(celltypes_filter_Pancreas[i]),"';
                             deg_markers_",format(celltypes_filter_Pancreas[i]),"_",format(VS[j,1]),"_Pancreas$Change <- ifelse(deg_markers_",format(celltypes_filter_Pancreas[i]),"_",format(VS[j,1]),"_Pancreas$avg_log2FC > 0 ,'Up','Down');
                             deg_markers_",format(celltypes_filter_Pancreas[i]),"_",format(VS[j,1]),"_Pancreas$Tissue <- 'Pancreas';
                             deg_markers_",format(celltypes_filter_Pancreas[i]),"_",format(VS[j,1]),"_Pancreas$Comparison <- '",format(VS[j,1]),"';",
                             format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_filter_Pancreas[i]),"_",format(VS[j,1]),"_Pancreas)};
                             rm(deg_markers_",format(celltypes_filter_Pancreas[i]),"_",format(VS[j,1]),"_Pancreas)")))
  }
  gc()
  eval(parse(text = paste0("DEG_Pancreas <- bind_rows(DEG_Pancreas,",format(VS[j,1]),");
                           rm(",format(VS[j,1]),")")))
}

write.csv(DEG_Pancreas,'DEG_Pancreas.csv')

################### 2.13. Salivary
###################
###################
###################
###################
###################
###################(1). Construct object for DEGs definition across multiple cell types and multiple conditions
DefaultAssay(Salivary) <- "RNA"
Salivary$Cell_name_condition <- paste(Salivary$Cell_name, Salivary$Condition, sep = "_")
Idents(Salivary) <- "Cell_name_condition"
table(Salivary$Cell_name_condition)
length(table(Salivary$Cell_name_condition))#107
length(table(Salivary$Cell_name))#28
#################(2). Select cell types for downstream DEG analysis
cell.number.condition.Salivary<-as.data.frame.matrix(table(Salivary$Cell_name,Salivary$Condition))
celltypes_filter_Salivary<-row.names(cell.number.condition.Salivary[which(cell.number.condition.Salivary$FD>2&cell.number.condition.Salivary$FS>2&cell.number.condition.Salivary$MS>2&cell.number.condition.Salivary$MC>2),])

VS <- data.frame(compare=c('MSvsFS','FDvsFS','MCvsMS'),
                 ident1=c('_MS','_FD','_MC'),
                 ident2=c('_FS','_FS','_MS'))
#################(3). Define DEGs among differenct conditions across multiple cell types
library(future)
plan()
plan(multisession, workers = 4)
plan()
DEG_Salivary <- data.frame()
for (j in (1:dim(VS)[1])){
  eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
  for (i in (1:length(celltypes_filter_Salivary))){
    assign(paste0("deg_markers_",celltypes_filter_Salivary[i],"_",format(VS[j,1]),"_Salivary",sep=""),FindMarkers(Salivary, ident.1 = paste0(celltypes_filter_Salivary[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_filter_Salivary[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
    eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_filter_Salivary[i]),"_",format(VS[j,1]),"_Salivary)[1]>0){deg_markers_",format(celltypes_filter_Salivary[i]),"_",format(VS[j,1]),"_Salivary <- subset(deg_markers_",format(celltypes_filter_Salivary[i]),"_",format(VS[j,1]),"_Salivary,p_val_adj<0.05)};
                             if(dim(deg_markers_",format(celltypes_filter_Salivary[i]),"_",format(VS[j,1]),"_Salivary)[1]>0){deg_markers_",format(celltypes_filter_Salivary[i]),"_",format(VS[j,1]),"_Salivary$Gene_symbol <- row.names(deg_markers_",format(celltypes_filter_Salivary[i]),"_",format(VS[j,1]),"_Salivary);
                             deg_markers_",format(celltypes_filter_Salivary[i]),"_",format(VS[j,1]),"_Salivary$Cell_type <- '",format(celltypes_filter_Salivary[i]),"';
                             deg_markers_",format(celltypes_filter_Salivary[i]),"_",format(VS[j,1]),"_Salivary$Change <- ifelse(deg_markers_",format(celltypes_filter_Salivary[i]),"_",format(VS[j,1]),"_Salivary$avg_log2FC > 0 ,'Up','Down');
                             deg_markers_",format(celltypes_filter_Salivary[i]),"_",format(VS[j,1]),"_Salivary$Tissue <- 'Salivary';
                             deg_markers_",format(celltypes_filter_Salivary[i]),"_",format(VS[j,1]),"_Salivary$Comparison <- '",format(VS[j,1]),"';",
                             format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_filter_Salivary[i]),"_",format(VS[j,1]),"_Salivary)};
                             rm(deg_markers_",format(celltypes_filter_Salivary[i]),"_",format(VS[j,1]),"_Salivary)")))
                             }
  gc()
  eval(parse(text = paste0("DEG_Salivary <- bind_rows(DEG_Salivary,",format(VS[j,1]),");
                           rm(",format(VS[j,1]),")")))
  }

write.csv(DEG_Salivary,'DEG_Salivary.csv')


################### 2.14. Skeletalmuscle
###################
###################
###################
###################
###################
###################(1). Construct object for DEGs definition across multiple cell types and multiple conditions
DefaultAssay(Skeletalmuscle) <- "RNA"
Skeletalmuscle$Cell_name_condition <- paste(Skeletalmuscle$Cell_name, Skeletalmuscle$Condition, sep = "_")
Idents(Skeletalmuscle) <- "Cell_name_condition"
table(Skeletalmuscle$Cell_name_condition)
length(table(Skeletalmuscle$Cell_name_condition))#107
length(table(Skeletalmuscle$Cell_name))#28
#################(2). Select cell types for downstream DEG analysis
cell.number.condition.Skeletalmuscle<-as.data.frame.matrix(table(Skeletalmuscle$Cell_name,Skeletalmuscle$Condition))
celltypes_filter_Skeletalmuscle<-row.names(cell.number.condition.Skeletalmuscle[which(cell.number.condition.Skeletalmuscle$FD>2&cell.number.condition.Skeletalmuscle$FS>2&cell.number.condition.Skeletalmuscle$MS>2&cell.number.condition.Skeletalmuscle$MC>2),])

VS <- data.frame(compare=c('MSvsFS','FDvsFS','MCvsMS'),
                 ident1=c('_MS','_FD','_MC'),
                 ident2=c('_FS','_FS','_MS'))
#################(3). Define DEGs among differenct conditions across multiple cell types
library(future)
plan()
plan(multisession, workers = 4)
plan()
DEG_Skeletalmuscle <- data.frame()
for (j in (1:dim(VS)[1])){
  eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
  for (i in (1:length(celltypes_filter_Skeletalmuscle))){
    assign(paste0("deg_markers_",celltypes_filter_Skeletalmuscle[i],"_",format(VS[j,1]),"_Skeletalmuscle",sep=""),FindMarkers(Skeletalmuscle, ident.1 = paste0(celltypes_filter_Skeletalmuscle[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_filter_Skeletalmuscle[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
    eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_filter_Skeletalmuscle[i]),"_",format(VS[j,1]),"_Skeletalmuscle)[1]>0){deg_markers_",format(celltypes_filter_Skeletalmuscle[i]),"_",format(VS[j,1]),"_Skeletalmuscle <- subset(deg_markers_",format(celltypes_filter_Skeletalmuscle[i]),"_",format(VS[j,1]),"_Skeletalmuscle,p_val_adj<0.05)};
                             if(dim(deg_markers_",format(celltypes_filter_Skeletalmuscle[i]),"_",format(VS[j,1]),"_Skeletalmuscle)[1]>0){deg_markers_",format(celltypes_filter_Skeletalmuscle[i]),"_",format(VS[j,1]),"_Skeletalmuscle$Gene_symbol <- row.names(deg_markers_",format(celltypes_filter_Skeletalmuscle[i]),"_",format(VS[j,1]),"_Skeletalmuscle);
                             deg_markers_",format(celltypes_filter_Skeletalmuscle[i]),"_",format(VS[j,1]),"_Skeletalmuscle$Cell_type <- '",format(celltypes_filter_Skeletalmuscle[i]),"';
                             deg_markers_",format(celltypes_filter_Skeletalmuscle[i]),"_",format(VS[j,1]),"_Skeletalmuscle$Change <- ifelse(deg_markers_",format(celltypes_filter_Skeletalmuscle[i]),"_",format(VS[j,1]),"_Skeletalmuscle$avg_log2FC > 0 ,'Up','Down');
                             deg_markers_",format(celltypes_filter_Skeletalmuscle[i]),"_",format(VS[j,1]),"_Skeletalmuscle$Tissue <- 'Skeletalmuscle';
                             deg_markers_",format(celltypes_filter_Skeletalmuscle[i]),"_",format(VS[j,1]),"_Skeletalmuscle$Comparison <- '",format(VS[j,1]),"';",
                             format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_filter_Skeletalmuscle[i]),"_",format(VS[j,1]),"_Skeletalmuscle)};
                             rm(deg_markers_",format(celltypes_filter_Skeletalmuscle[i]),"_",format(VS[j,1]),"_Skeletalmuscle)")))
                             }
  gc()
  eval(parse(text = paste0("DEG_Skeletalmuscle <- bind_rows(DEG_Skeletalmuscle,",format(VS[j,1]),");
                           rm(",format(VS[j,1]),")")))
  }

write.csv(DEG_Skeletalmuscle,'DEG_Skeletalmuscle.csv')
                             
################### 2.15. Spleen
###################
###################
###################
###################
###################
###################(1). Construct object for DEGs definition across multiple cell types and multiple conditions
DefaultAssay(Spleen) <- "RNA"
Spleen$Cell_name_condition <- paste(Spleen$Cell_name, Spleen$Condition, sep = "_")
Idents(Spleen) <- "Cell_name_condition"
table(Spleen$Cell_name_condition)
length(table(Spleen$Cell_name_condition))#107
length(table(Spleen$Cell_name))#28
#################(2). Select cell types for downstream DEG analysis
cell.number.condition.Spleen<-as.data.frame.matrix(table(Spleen$Cell_name,Spleen$Condition))
celltypes_filter_Spleen<-row.names(cell.number.condition.Spleen[which(cell.number.condition.Spleen$FD>2&cell.number.condition.Spleen$FS>2&cell.number.condition.Spleen$MS>2&cell.number.condition.Spleen$MC>2),])

VS <- data.frame(compare=c('MSvsFS','FDvsFS','MCvsMS'),
                 ident1=c('_MS','_FD','_MC'),
                 ident2=c('_FS','_FS','_MS'))
#################(3). Define DEGs among differenct conditions across multiple cell types
library(future)
plan()
plan(multisession, workers = 4)
plan()
DEG_Spleen <- data.frame()
for (j in (1:dim(VS)[1])){
  eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
  for (i in (1:length(celltypes_filter_Spleen))){
    assign(paste0("deg_markers_",celltypes_filter_Spleen[i],"_",format(VS[j,1]),"_Spleen",sep=""),FindMarkers(Spleen, ident.1 = paste0(celltypes_filter_Spleen[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_filter_Spleen[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
    eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_filter_Spleen[i]),"_",format(VS[j,1]),"_Spleen)[1]>0){deg_markers_",format(celltypes_filter_Spleen[i]),"_",format(VS[j,1]),"_Spleen <- subset(deg_markers_",format(celltypes_filter_Spleen[i]),"_",format(VS[j,1]),"_Spleen,p_val_adj<0.05)};
                             if(dim(deg_markers_",format(celltypes_filter_Spleen[i]),"_",format(VS[j,1]),"_Spleen)[1]>0){deg_markers_",format(celltypes_filter_Spleen[i]),"_",format(VS[j,1]),"_Spleen$Gene_symbol <- row.names(deg_markers_",format(celltypes_filter_Spleen[i]),"_",format(VS[j,1]),"_Spleen);
                             deg_markers_",format(celltypes_filter_Spleen[i]),"_",format(VS[j,1]),"_Spleen$Cell_type <- '",format(celltypes_filter_Spleen[i]),"';
                             deg_markers_",format(celltypes_filter_Spleen[i]),"_",format(VS[j,1]),"_Spleen$Change <- ifelse(deg_markers_",format(celltypes_filter_Spleen[i]),"_",format(VS[j,1]),"_Spleen$avg_log2FC > 0 ,'Up','Down');
                             deg_markers_",format(celltypes_filter_Spleen[i]),"_",format(VS[j,1]),"_Spleen$Tissue <- 'Spleen';
                             deg_markers_",format(celltypes_filter_Spleen[i]),"_",format(VS[j,1]),"_Spleen$Comparison <- '",format(VS[j,1]),"';",
                             format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_filter_Spleen[i]),"_",format(VS[j,1]),"_Spleen)};
                             rm(deg_markers_",format(celltypes_filter_Spleen[i]),"_",format(VS[j,1]),"_Spleen)")))
                             }
  gc()
  eval(parse(text = paste0("DEG_Spleen <- bind_rows(DEG_Spleen,",format(VS[j,1]),");
                           rm(",format(VS[j,1]),")")))
  }

write.csv(DEG_Spleen,'DEG_Spleen.csv')

################### 2.16. Stomach
###################
###################
###################
###################
###################
###################(1). Construct object for DEGs definition across multiple cell types and multiple conditions
DefaultAssay(Stomach) <- "RNA"
Stomach$Cell_name_condition <- paste(Stomach$Cell_name, Stomach$Condition, sep = "_")
Idents(Stomach) <- "Cell_name_condition"
table(Stomach$Cell_name_condition)
length(table(Stomach$Cell_name_condition))#107
length(table(Stomach$Cell_name))#28
#################(2). Select cell types for downstream DEG analysis
cell.number.condition.Stomach<-as.data.frame.matrix(table(Stomach$Cell_name,Stomach$Condition))
celltypes_filter_Stomach<-row.names(cell.number.condition.Stomach[which(cell.number.condition.Stomach$FD>2&cell.number.condition.Stomach$FS>2&cell.number.condition.Stomach$MS>2&cell.number.condition.Stomach$MC>2),])

VS <- data.frame(compare=c('MSvsFS','FDvsFS','MCvsMS'),
                 ident1=c('_MS','_FD','_MC'),
                 ident2=c('_FS','_FS','_MS'))
#################(3). Define DEGs among differenct conditions across multiple cell types
library(future)
plan()
plan(multisession, workers = 4)
plan()
DEG_Stomach <- data.frame()
for (j in (1:dim(VS)[1])){
  eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
  for (i in (1:length(celltypes_filter_Stomach))){
    assign(paste0("deg_markers_",celltypes_filter_Stomach[i],"_",format(VS[j,1]),"_Stomach",sep=""),FindMarkers(Stomach, ident.1 = paste0(celltypes_filter_Stomach[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_filter_Stomach[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
    eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_filter_Stomach[i]),"_",format(VS[j,1]),"_Stomach)[1]>0){deg_markers_",format(celltypes_filter_Stomach[i]),"_",format(VS[j,1]),"_Stomach <- subset(deg_markers_",format(celltypes_filter_Stomach[i]),"_",format(VS[j,1]),"_Stomach,p_val_adj<0.05)};
                             if(dim(deg_markers_",format(celltypes_filter_Stomach[i]),"_",format(VS[j,1]),"_Stomach)[1]>0){deg_markers_",format(celltypes_filter_Stomach[i]),"_",format(VS[j,1]),"_Stomach$Gene_symbol <- row.names(deg_markers_",format(celltypes_filter_Stomach[i]),"_",format(VS[j,1]),"_Stomach);
                             deg_markers_",format(celltypes_filter_Stomach[i]),"_",format(VS[j,1]),"_Stomach$Cell_type <- '",format(celltypes_filter_Stomach[i]),"';
                             deg_markers_",format(celltypes_filter_Stomach[i]),"_",format(VS[j,1]),"_Stomach$Change <- ifelse(deg_markers_",format(celltypes_filter_Stomach[i]),"_",format(VS[j,1]),"_Stomach$avg_log2FC > 0 ,'Up','Down');
                             deg_markers_",format(celltypes_filter_Stomach[i]),"_",format(VS[j,1]),"_Stomach$Tissue <- 'Stomach';
                             deg_markers_",format(celltypes_filter_Stomach[i]),"_",format(VS[j,1]),"_Stomach$Comparison <- '",format(VS[j,1]),"';",
                             format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_filter_Stomach[i]),"_",format(VS[j,1]),"_Stomach)};
                             rm(deg_markers_",format(celltypes_filter_Stomach[i]),"_",format(VS[j,1]),"_Stomach)")))
  }
  gc()
  eval(parse(text = paste0("DEG_Stomach <- bind_rows(DEG_Stomach,",format(VS[j,1]),");
                           rm(",format(VS[j,1]),")")))
}

write.csv(DEG_Stomach,'DEG_Stomach.csv')

################### 2.17. Thymus
###################
###################
###################
###################
###################
###################(1). Construct object for DEGs definition across multiple cell types and multiple conditions
DefaultAssay(Thymus) <- "RNA"
Thymus$Cell_name_condition <- paste(Thymus$Cell_name, Thymus$Condition, sep = "_")
Idents(Thymus) <- "Cell_name_condition"
table(Thymus$Cell_name_condition)
length(table(Thymus$Cell_name_condition))#107
length(table(Thymus$Cell_name))#28
#################(2). Select cell types for downstream DEG analysis
cell.number.condition.Thymus<-as.data.frame.matrix(table(Thymus$Cell_name,Thymus$Condition))
celltypes_filter_Thymus<-row.names(cell.number.condition.Thymus[which(cell.number.condition.Thymus$FD>2&cell.number.condition.Thymus$FS>2&cell.number.condition.Thymus$MS>2&cell.number.condition.Thymus$MC>2),])

VS <- data.frame(compare=c('MSvsFS','FDvsFS','MCvsMS'),
                 ident1=c('_MS','_FD','_MC'),
                 ident2=c('_FS','_FS','_MS'))
#################(3). Define DEGs among differenct conditions across multiple cell types
library(future)
plan()
plan(multisession, workers = 4)
plan()
DEG_Thymus <- data.frame()
for (j in (1:dim(VS)[1])){
  eval(parse(text = paste0(format(VS[j,1])," <- data.frame()")))
  for (i in (1:length(celltypes_filter_Thymus))){
    assign(paste0("deg_markers_",celltypes_filter_Thymus[i],"_",format(VS[j,1]),"_Thymus",sep=""),FindMarkers(Thymus, ident.1 = paste0(celltypes_filter_Thymus[i],format(VS[j,2]),"",sep=""), ident.2 = paste0(celltypes_filter_Thymus[i],format(VS[j,3]),"",sep=""), logfc.threshold = 0.5, verbose = FALSE))
    eval(parse(text = paste0("if(dim(deg_markers_",format(celltypes_filter_Thymus[i]),"_",format(VS[j,1]),"_Thymus)[1]>0){deg_markers_",format(celltypes_filter_Thymus[i]),"_",format(VS[j,1]),"_Thymus <- subset(deg_markers_",format(celltypes_filter_Thymus[i]),"_",format(VS[j,1]),"_Thymus,p_val_adj<0.05)};
                             if(dim(deg_markers_",format(celltypes_filter_Thymus[i]),"_",format(VS[j,1]),"_Thymus)[1]>0){deg_markers_",format(celltypes_filter_Thymus[i]),"_",format(VS[j,1]),"_Thymus$Gene_symbol <- row.names(deg_markers_",format(celltypes_filter_Thymus[i]),"_",format(VS[j,1]),"_Thymus);
                             deg_markers_",format(celltypes_filter_Thymus[i]),"_",format(VS[j,1]),"_Thymus$Cell_type <- '",format(celltypes_filter_Thymus[i]),"';
                             deg_markers_",format(celltypes_filter_Thymus[i]),"_",format(VS[j,1]),"_Thymus$Change <- ifelse(deg_markers_",format(celltypes_filter_Thymus[i]),"_",format(VS[j,1]),"_Thymus$avg_log2FC > 0 ,'Up','Down');
                             deg_markers_",format(celltypes_filter_Thymus[i]),"_",format(VS[j,1]),"_Thymus$Tissue <- 'Thymus';
                             deg_markers_",format(celltypes_filter_Thymus[i]),"_",format(VS[j,1]),"_Thymus$Comparison <- '",format(VS[j,1]),"';",
                             format(VS[j,1])," <- bind_rows(",format(VS[j,1]),",deg_markers_",format(celltypes_filter_Thymus[i]),"_",format(VS[j,1]),"_Thymus)};
                             rm(deg_markers_",format(celltypes_filter_Thymus[i]),"_",format(VS[j,1]),"_Thymus)")))
                             }
  gc()
  eval(parse(text = paste0("DEG_Thymus <- bind_rows(DEG_Thymus,",format(VS[j,1]),");
                           rm(",format(VS[j,1]),")")))
  }

write.csv(DEG_Thymus,'DEG_Thymus.csv')
}

################### 3. Constrcut a DEG list
###################
###################
###################
###################
###################
################### 3.1. Combine all the DEGs from each tissue
DEG_Adipose<-read.csv("DEG_Adipose.csv")
DEG_Adrenal<-read.csv("DEG_Adrenal.csv")
DEG_Bonemarrow<-read.csv("DEG_Bonemarrow.csv")
DEG_Brain<-read.csv("DEG_Brain.csv")
DEG_Colon<-read.csv("DEG_Colon.csv")
DEG_Heart<-read.csv("DEG_Heart.csv")
DEG_Intestine<-read.csv("DEG_Intestine.csv")
DEG_Kidney<-read.csv("DEG_Kidney.csv")
DEG_Lacrimal<-read.csv("DEG_Lacrimal.csv")
DEG_Liver<-read.csv("DEG_Liver.csv")
DEG_Lung<-read.csv("DEG_Lung.csv")
DEG_Pancreas<-read.csv("DEG_Pancreas.csv")
DEG_Salivary<-read.csv("DEG_Salivary.csv")
DEG_Skeletalmuscle<-read.csv("DEG_Skeletalmuscle.csv")
DEG_Spleen<-read.csv("DEG_Spleen.csv")
DEG_Stomach<-read.csv("DEG_Stomach.csv")
DEG_Thymus<-read.csv("DEG_Thymus.csv")

DEG_all<-rbind(DEG_Adipose,DEG_Adrenal,DEG_Bonemarrow,DEG_Brain,
               DEG_Colon,DEG_Heart,DEG_Intestine,DEG_Kidney,
               DEG_Lacrimal,DEG_Liver,DEG_Lung,DEG_Pancreas,
               DEG_Salivary,DEG_Skeletalmuscle,DEG_Spleen,DEG_Stomach,DEG_Thymus)

table(DEG_all$Tissue)
table(DEG_all$Cell_type)

dategene<-read.table("genename.txt",header = T)
dategene1<-intersect(DEG_all$Gene_symbol,dategene$gene_symbol_old)
for (i in 1:length(dategene1)){
 DEG_all$Gene_symbol<-gsub(dategene1[i],dategene[which(dategene$gene_symbol_old==dategene1[i]),]$gene_symbol_new,DEG_all$Gene_symbol)
}

write.csv(DEG_all,"DEG_all.csv")

################### 3.2. Remove cell types
DEG_filter<-DEG_all[which(DEG_all$Cell_type!="Erythroblast"&DEG_all$Cell_type!="Erythrocyte"&DEG_all$Cell_type!="Proliferating"&
                            DEG_all$Cell_type!="B_proliferating"&
                            DEG_all$Cell_type!="Neutrophil_proliferating"),]

table(DEG_filter$Tissue)
write.csv(DEG_filter,"DEG_filter.csv")
