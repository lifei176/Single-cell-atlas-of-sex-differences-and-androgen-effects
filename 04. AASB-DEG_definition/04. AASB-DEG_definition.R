################### Identify AASB-DEGs
###################
###################
###################
###################
###################
DEG<-read.csv("/.../DEG_filter.csv")
tissue<-names(table(DEG$Tissue))
data.frame.AASB<-data.frame()
AASB_count<-data.frame()
AASB_tissue<-data.frame()
for (i in 1:length(tissue)){
  DEG_1<-DEG[which(DEG$Tissue==tissue[i]),]
  DEG_1_MSvsFS<-droplevels(DEG_1[DEG_1$Comparison=="MSvsFS",])
  DEG_1_MSvsFS_up<-droplevels(DEG_1_MSvsFS[DEG_1_MSvsFS$Change=="Up",])
  DEG_1_MSvsFS_down<-droplevels(DEG_1_MSvsFS[DEG_1_MSvsFS$Change=="Down",])
  
  DEG_1_MCvsMS<-droplevels(DEG_1[DEG_1$Comparison=="MCvsMS",])
  DEG_1_MCvsMS_up<-droplevels(DEG_1_MCvsMS[DEG_1_MCvsMS$Change=="Up",])
  DEG_1_MCvsMS_down<-droplevels(DEG_1_MCvsMS[DEG_1_MCvsMS$Change=="Down",])
  
  DEG_1_FDvsFS<-droplevels(DEG_1[DEG_1$Comparison=="FDvsFS",])
  DEG_1_FDvsFS_up<-droplevels(DEG_1_FDvsFS[DEG_1_FDvsFS$Change=="Up",])
  DEG_1_FDvsFS_down<-droplevels(DEG_1_FDvsFS[DEG_1_FDvsFS$Change=="Down",])
  
  ################### Define the cell types with up-regulated DEGs or down-regulated DEGs
  ###################
  ###################
  ###################
  ###################
  ###################
  cell_up<-intersect(names(table(DEG_1_MSvsFS_up$Cell_type)),intersect(names(table(DEG_1_MCvsMS_down$Cell_type)),
                                                                       names(table(DEG_1_FDvsFS_up$Cell_type))))
  
  cell_down<-intersect(names(table(DEG_1_MSvsFS_down$Cell_type)),intersect(names(table(DEG_1_MCvsMS_up$Cell_type)),
                                                                           names(table(DEG_1_FDvsFS_down$Cell_type))))
  
  ################### Identify the positive AASB-DEGs
  ###################
  ###################
  ###################
  ###################
  ###################
  data.frame.up<-data.frame()
  for (j in 1:length(cell_up)){
    cell_up_MSvsFS<-as.character(droplevels(DEG_1_MSvsFS_up[which(DEG_1_MSvsFS_up$Cell_type==cell_up[j]),])$Gene_symbol)
    cell_down_MCvsMS<-as.character(droplevels(DEG_1_MCvsMS_down[which(DEG_1_MCvsMS_down$Cell_type==cell_up[j]),])$Gene_symbol)
    cell_up_FDvsFS<-as.character(droplevels(DEG_1_FDvsFS_up[which(DEG_1_FDvsFS_up$Cell_type==cell_up[j]),])$Gene_symbol)
    positive_AASB_DEG<-intersect(cell_up_MSvsFS,intersect(cell_down_MCvsMS,cell_up_FDvsFS))
    data.frame.up1<-data.frame(Tissue=rep(tissue[i],length(positive_AASB_DEG)),
                               Cell_type=rep(cell_up[j],length(positive_AASB_DEG)),
                               Change=rep("Positive",length(positive_AASB_DEG)),
                               Gene=positive_AASB_DEG)
    data.frame.up<-rbind(data.frame.up,data.frame.up1)
  }
  
  ################### Identify the negative AASB-DEGs
  ###################
  ###################
  ###################
  ###################
  ###################
  data.frame.down<-data.frame()
  for (h in 1:length(cell_down)){
    cell_down_MSvsFS<-as.character(droplevels(DEG_1_MSvsFS_down[which(DEG_1_MSvsFS_down$Cell_type==cell_down[h]),])$Gene_symbol)
    cell_up_MCvsMS<-as.character(droplevels(DEG_1_MCvsMS_up[which(DEG_1_MCvsMS_up$Cell_type==cell_down[h]),])$Gene_symbol)
    cell_down_FDvsFS<-as.character(droplevels(DEG_1_FDvsFS_down[which(DEG_1_FDvsFS_down$Cell_type==cell_down[h]),])$Gene_symbol)
    negative_AASB_DEG<-intersect(cell_down_MSvsFS,intersect(cell_up_MCvsMS,cell_down_FDvsFS))
    data.frame.down1<-data.frame(Tissue=rep(tissue[i],length(negative_AASB_DEG)),
                                 Cell_type=rep(cell_down[h],length(negative_AASB_DEG)),
                                 Change=rep("Negative",length(negative_AASB_DEG)),
                                 Gene=negative_AASB_DEG)
    data.frame.down<-rbind(data.frame.down,data.frame.down1)
  }
  data.frame.AASB1<-rbind(data.frame.up,data.frame.down)
  data.frame.AASB<-rbind(data.frame.AASB,data.frame.AASB1)
  
  AASB_count1<-data.frame(Tissue=tissue[i],Positive=length(unique(data.frame.up$Gene)),Negative=length(unique(data.frame.down$Gene)))
  AASB_count<-rbind(AASB_count,AASB_count1)
  
  AASB_tissue_positive<-data.frame(Tissue=rep(tissue[i],length(unique(data.frame.up$Gene))),
                                   Change= rep("Positive",length(unique(data.frame.up$Gene))),          
                                   Gene= unique(data.frame.up$Gene))
  AASB_tissue_negative<-data.frame(Tissue=rep(tissue[i],length(unique(data.frame.down$Gene))),
                                   Change= rep("Negative",length(unique(data.frame.down$Gene))),          
                                   Gene=unique(data.frame.down$Gene))                               
  AASB_tissue1<-rbind(AASB_tissue_positive,AASB_tissue_negative)                                
  
  AASB_tissue<-rbind(AASB_tissue,AASB_tissue1)#unique gene
}
AASB_DEG<-data.frame.AASB
write.csv(AASB_DEG,"/.../AASB_DEG.csv")
write.csv(AASB_count,"/.../AASB_count.csv")
write.csv(AASB_tissue,"/.../AASB_tissue.csv")
