library(Seurat)
library(cowplot)
library(tidyverse)
library(fgsea)
library(msigdbr)
library(clustree)
library(xtable)
library(flextable)
library(ggprism)
library(ggthemes)
library(ggpubr)
library(RColorBrewer)




control.data<-Read10X('Control_MG_10X_Matrix/')
treatment.data<-Read10X('Treatment_MG_10X_Matrix/')

control<-CreateSeuratObject(counts = control.data, project='mammary gland', min.cells = 3, min.features = 200)
treatment<-CreateSeuratObject(counts = treatment.data, project='mammary gland', min.cells = 3, min.features = 200)


#28360*7930
control[["percent.mt"]] <- PercentageFeatureSet(control, pattern = "^mt-")
control$sample='control'
#VlnPlot(control,features = c("nCount_RNA", "nFeature_RNA", 'percent.mt'))
control <- subset(control, subset = nFeature_RNA > 300 & nFeature_RNA < 7500 
                  & percent.mt < 5)
#28360*7706
control<-NormalizeData(control, normalization.method = "LogNormalize", scale.factor = 10000)
control<-FindVariableFeatures(control, selection.method = "vst", nfeatures = 2000)


#28199*8624
treatment[["percent.mt"]] <- PercentageFeatureSet(treatment, pattern = "^mt")
treatment$sample='treatment'
VlnPlot(treatment,features = c("nCount_RNA", "nFeature_RNA", 'percent.mt'))
treatment <- subset(treatment, subset = nFeature_RNA > 300 & nFeature_RNA < 7500 
                    & percent.mt < 5)
#28199*8516
treatment<-NormalizeData(treatment, normalization.method = "LogNormalize", scale.factor = 10000)
treatment<-FindVariableFeatures(treatment, selection.method = "vst", nfeatures = 2000)


features <- SelectIntegrationFeatures(object.list = list(control,treatment))
anchors <- FindIntegrationAnchors(object.list = list(control,treatment), anchor.features = features)
mg.combined<-IntegrateData(anchorset = anchors)
#RNA 29715*16222
#integrated 2000*16222
DefaultAssay(mg.combined) <- "integrated"


Idents(mg.combined) <- ""
VlnPlot(mg.combined,features = c("nCount_RNA", "nFeature_RNA", 'percent.mt'),
        group.by = 'sample')


mg.combined <- ScaleData(mg.combined)
mg.combined <- RunPCA(mg.combined)
ElbowPlot(mg.combined)


mg.combined <- RunUMAP(mg.combined, reduction = "pca", dims = 1:30)
mg.combined <- FindNeighbors(mg.combined, reduction = "pca", dims = 1:30)
mg.combined <- FindClusters(mg.combined, resolution = c(0.1,0.2,0.3,0.5,0.7,0.8,0.9))



clustree(mg.combined,prefix='integrated_snn_res.')
Idents(mg.combined) <- "integrated_snn_res.0.3"


p1 <- DimPlot(mg.combined, reduction = "umap",label=T,
              cols =  brewer.pal(11,'Paired'))
p2 <- DimPlot(mg.combined, reduction = "umap",group.by ='sample',label=T)
plot_grid(p1, p2)
p3 <- DimPlot(mg.combined, reduction = "umap", split.by = 'sample',
              cols =  brewer.pal(11,'Paired'))


DefaultAssay(mg.combined) <- "RNA"
#luminal 4,11,12
VlnPlot(mg.combined, features = c('Krt8','Krt18'))
VlnPlot(mg.combined, features = c("Pgr"))
VlnPlot(mg.combined, features = c("Wfdc18", "Elf5"))
#basal 6
VlnPlot(mg.combined, features = c('Krt5','Krt14','Lgr5'))
#Adig Mrap Abhd15 Adipoq 0,2,8
#Endothelial Sox17 9
#Fibroblast Dcn,Col1a1 1,5
#Smooth muscle cell Irag1 13
#Macrophage (B) Ms4a7 14
#Macrophage (A) Cd163 3
#T cell Cd3d 10
#Dendritic cell Flt3  7

#marker.0 <- FindMarkers(mg.combined, ident.1 = 0, min.pct = 0.25)
#marker.1 <- FindMarkers(mg.combined, ident.1 = 1, min.pct = 0.25)
#marker.7 <- FindMarkers(mg.combined, ident.1 = 0,ident.2 = 1, min.pct = 0.25)


mg.combined <- RenameIdents(mg.combined, 
                            '0'='Adipocyte','1'= 'Fibroblast','2'='Adipocyte',
                            '3'='Macrophage A','4'='Luminal AV','5'= 'Fibroblast',
                            '6'='Basal', '7'='Dendritic cell','8'='Adipocyte',
                            '9'='Endothelial','10'='T cell','11'='Luminal HS',
                            '12'='Luminal AV','13'='Smooth muscle','14'='Macrophage B')


DotPlot(mg.combined,features = 
          rev(c('Adig','Adipoq','Abhd15','Dcn','Col1a1','Cd163','Mrc1','Cd209f',
                'Elf5','Csn3','Wfdc18','Krt5','Krt14','Lgr5','Flt3','Traf1','Cd209a',
                'Sox17','Flt4','Cd3d','Cd3e','Cd3g','Pgr','Prlr','Esr1',
                'Irag1','Flna','Mmp12','Mmp13','Spic')),
        cols = c('lightgrey', 'red'))+
  theme(axis.text.x = element_text(angle = -45, hjust = -0.1))




ctrl_cluster <- table(Idents(mg.combined)[mg.combined$sample == "control"])
trt_cluster <- table(Idents(mg.combined)[mg.combined$sample == "treatment"])


ctrl_prop <- as.numeric(ctrl_cluster)/ sum(as.numeric(ctrl_cluster)) * 100
trt_prop <- as.numeric(trt_cluster)/ sum(as.numeric(trt_cluster)) * 100


cluster_df<- rbind((ctrl_cluster), (trt_cluster))
rownames(cluster_df) <- c('Control','Treatment')
df<-as.data.frame(cluster_df)
df[1,] <- paste0(as.numeric(ctrl_cluster),' (', round(ctrl_prop,2),'%)')
df[2,] <- paste0(as.numeric(trt_cluster),' (', round(trt_prop,2),'%)')


as_flextable(xtable(df))
doc = read_docx()
doc = body_add_flextable(doc,as_flextable(xtable(cluster_df)))
print(doc,"tt.docx")


chisq.test(cluster_df,correct = F)


cluster_df.2<- cbind((ctrl_cluster), (trt_cluster))
colnames(cluster_df.2) <- c('Control','Treatment')
cluster_df.2 <- as.data.frame(cluster_df.2)
cluster_df.2$Cluster<-rownames(cluster_df.2)
rownames(cluster_df.2)<- NULL




p1 <- ggplot(cluster_df.2) + 
  geom_col(
    aes(x = Cluster,y = Control,fill = Control),
    position = "dodge2",show.legend = TRUE,alpha = 0.9)+
  geom_hline(
    aes(yintercept = y), data.frame(y = c(0:6) * 500),color = "lightgrey") +
  scale_y_continuous(limits = c(-1000, 3500)) +
  geom_segment(
    aes(x = Cluster,y = 0,xend = Cluster, yend = 3000),linetype = "dashed",color = "gray12") +
  scale_fill_gradientn(
    colours = c("#6C5B7B", "#C06C84", "#F67280", "#F8B195")) +
  guides(
    fill = guide_colorsteps(
      barwidth = 15, barheight = 0.5, title.position = "top", title.hjust = 0.5)) +
  geom_text(aes(x=Cluster, y=Control,label=paste0(as.numeric(ctrl_cluster),' (', round(ctrl_prop,2),'%)')),color="gray12",size=3)+
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = -45 ,color = "gray12", size = 7),
    legend.position = "bottom")+
  coord_polar()


p2 <- ggplot(cluster_df.2) + 
  geom_col(
    aes(x = Cluster,y = Treatment,fill = Treatment),
    position = "dodge2",show.legend = TRUE,alpha = 0.9)+
  geom_hline(
    aes(yintercept = y), data.frame(y = c(0:7) * 500),color = "lightgrey") +
  scale_y_continuous(limits = c(-1000, 3500)) +
  geom_segment(
    aes(x = Cluster,y = 0,xend = Cluster, yend = 3500),linetype = "dashed",color = "gray12") +
  scale_fill_gradientn(
    colours = c("#6C5B7B", "#C06C84", "#F67280", "#F8B195")) +
  guides(
    fill = guide_colorsteps(
      barwidth = 15, barheight = 0.5, title.position = "top", title.hjust = 0.5)) +
  geom_text(aes(x=Cluster, y=Treatment,label=paste0(as.numeric(trt_cluster),' (', round(trt_prop,2),'%)')),color="gray12",size=3)+
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = -45 ,color = "gray12", size = 7),
    legend.position = "bottom")+
  coord_polar()




p4<-ggplot(cluster_df.3)+
  geom_bar(aes(Sample,Number,fill=Type),
           stat='identity',width = 0.5,position = 'fill')+
  scale_fill_manual(values=c("#E7E1EF", "#C994C7", "#C06C84",
                                      "#EFF3FF", "#BDD7E7", "#6BAED6", "#2171B5",
                                      "#A6611A", "#DFC27D", "#80CDC1", "#018571"))+
                                        labs(y='Proportion')+
  theme(panel.grid = element_blank())




cluster_df.3<- gather(cluster_df.2,Sample,Number,-Cluster)
cluster_df.3$Type <- NA
cluster_df.3$Type <- ifelse(cluster_df.3$Cluster %in% 
                              c('Luminal HS','Luminal AV', 'Basal'),
                            'Epithelial cells',
                            'Stromal cells')
cluster_df.3$Type[cluster_df.3$Cluster %in% 
                    c("Macrophage A", "Dendritic cell", "Macrophage B", "T cell")] <- 'Immune cells'


p1<-ggplot(cluster_df.3[cluster_df.3$Type=='Epithelial cells',])+
  geom_bar(aes(Sample,Number,fill=Cluster),
           stat='identity',width = 0.5,position = 'fill')+
  facet_grid(~Type)+
  scale_fill_manual(values=c("#FEE8C8", "#FDBB84", "#E34A33",
                                      "#EFF3FF", "#BDD7E7", "#6BAED6", "#2171B5",
                                      "#A6611A", "#DFC27D", "#80CDC1", "#018571"))+
                                        labs(y='Proportion')+
  theme(panel.grid = element_blank())


p2<-ggplot(cluster_df.3[cluster_df.3$Type=='Stromal cells',])+
  geom_bar(aes(Sample,Number,fill=Cluster),
           stat='identity',width = 0.5,position = 'fill')+
  facet_grid(~Type)+
  scale_fill_manual(values=c("#A6611A", "#DFC27D", "#80CDC1", "#018571"))+
  labs(y='Proportion')+
  theme(panel.grid = element_blank())


p3<-ggplot(cluster_df.3[cluster_df.3$Type=='Immune cells',])+
  geom_bar(aes(Sample,Number,fill=Cluster),
           stat='identity',width = 0.5,position = 'fill')+
  facet_grid(~Type)+
  scale_fill_manual(values=c("#EFF3FF", "#BDD7E7", "#6BAED6", "#2171B5",
                                      "#A6611A", "#DFC27D", "#80CDC1", "#018571"))+
                                        labs(y='Proportion')+
  theme(panel.grid = element_blank())


ggarrange(p4,p1,p2,p3,
          labels = c('C','D','E','F'))




prop_df <- data.frame('Cluster'=colnames(cluster_df),
                      'Control'=ctrl_prop,
                      'Treatment'=trt_prop,
                      'Difference'=(trt_prop-ctrl_prop))
#23.82 11.16
#62.64 67.72
#13.54 21.11


p1<-ggplot(prop_df) +
  geom_histogram(aes(x = Difference, y = Cluster, fill = ifelse(Difference > 0, "Greater than 0", "Less than or equal to 0")), stat = 'identity',show.legend=F) +
  xlab("Absolute Proportion Difference (%)") + ylab("Cluster") + #ggtitle("Distribution of Proportion Difference by Cluster") +
  scale_fill_manual(values = c("Greater than 0" = "#D73027", "Less than or equal to 0" = "#313695"))




cell_type<-colnames(cluster_df)
reserved_markers<-c()
for (i in 1:11){
  markers<-FindConservedMarkers(mg.combined, ident.1 = cell_type[i], grouping.var = "sample")
  reserved_markers<-c(reserved_markers,
                      rownames(markers %>% 
                                 arrange(desc(control_avg_log2FC)) %>% 
                                 head(2)))
}




DimPlot(mg.combined,split.by = 'sample',
        cells.highlight=list(
          Luminal_HS=WhichCells(mg.combined, idents='Luminal HS'),
          Luminal_AV=WhichCells(mg.combined, idents='Luminal AV'),
          Basal=WhichCells(mg.combined, idents='Basal'),
          Adipocyte=WhichCells(mg.combined, idents='Adipocyte')),
        cols.highlight=c("#2B83BA", "#ABDDA4", "#FDAE61","#D7191C"),label = T)




mg.combined$cell_type_sample<- paste(Idents(mg.combined), mg.combined$sample, sep = "_")
Idents(mg.combined) <- "cell_type_sample"


mg.combined$cell_type <- NA
mg.combined$cell_type[mg.combined$integrated_snn_res.0.3 %in% c(4,12)] <- "Luminal_AV"
mg.combined$cell_type[mg.combined$integrated_snn_res.0.3 %in% c(11)] <- "Luminal_HS"
mg.combined$cell_type[mg.combined$integrated_snn_res.0.3 %in% c(6)] <- "Basal"
mg.combined$cell_type[mg.combined$integrated_snn_res.0.3 %in% c(0,2,8)] <- "Adipocyte"




luminal_AVs <- subset(mg.combined, cell_type == "Luminal_AV")
Idents(luminal_AVs) <- "sample"
avg.luminal_AVs <- log1p(AverageExpression(luminal_AVs)$RNA)

p1<-ggplot(as.data.frame(avg.luminal_AVs), aes(control, treatment)) + geom_point() + ggtitle("Luminal AV")

luminal_AV_markers <- FindMarkers(mg.combined, ident.1 = "Luminal AV_treatment", ident.2 = "Luminal AV_control")

luminal_AV_top <- luminal_AV_markers %>% 
  filter(abs(avg_log2FC) > 1 & p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC))

p1<-LabelPoints(plot = p1, 
                points = rownames(luminal_AV_top),repel = T)


luminal_HSs <- subset(mg.combined, cell_type == "Luminal_HS")
Idents(luminal_HSs) <- "sample"
avg.luminal_HSs <- log1p(AverageExpression(luminal_HSs)$RNA)

p2<-ggplot(as.data.frame(avg.luminal_HSs), aes(control, treatment)) + geom_point() + ggtitle("Luminal HS")

luminal_HS_markers <- FindMarkers(mg.combined, ident.1 = "Luminal HS_treatment", ident.2 = "Luminal HS_control")

top5_max <- luminal_HS_markers %>% 
  filter(abs(avg_log2FC) > 1 & p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC)) %>% 
  head(5)

top5_min <- luminal_HS_markers %>% 
  filter(abs(avg_log2FC) > 1 & p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC)) %>% 
  tail(5)

luminal_HS_top<-rbind(top5_max,top5_min)

p2<-LabelPoints(plot = p2, 
                points = rownames(luminal_HS_top),repel = T)


basals <- subset(mg.combined, cell_type == "Basal")
Idents(basals) <- "sample"
avg.basals <- log1p(AverageExpression(basals)$RNA)

p3<-ggplot(as.data.frame(avg.basals), aes(control, treatment)) + geom_point() + ggtitle("Basal")

basal_markers <- FindMarkers(mg.combined, ident.1 = "Basal_treatment", ident.2 = "Basal_control")

basal_top <- basal_markers %>% 
  filter(abs(avg_log2FC) > 1 & p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC))

p3<-LabelPoints(plot = p3, 
                points = rownames(basal_top),repel = T)


adipocytes <- subset(mg.combined, cell_type == "Adipocyte")
Idents(adipocytes) <- "sample"
avg.adipocytes <- log1p(AverageExpression(adipocytes)$RNA)

p4<-ggplot(as.data.frame(avg.adipocytes), aes(control, treatment)) + geom_point() + ggtitle("Adipocyte")

adipocyte_markers <- FindMarkers(mg.combined, ident.1 = "Adipocyte_treatment", ident.2 = "Adipocyte_control")

top5_max <- adipocyte_markers %>% 
  filter(abs(avg_log2FC) > 1 & p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC)) %>% 
  head(4)

top5_min <- adipocyte_markers %>% 
  filter(abs(avg_log2FC) > 1 & p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC)) %>% 
  rownames()
tail(5)

adipocyte_top <- rbind(top5_max,top5_min)

p4<-LabelPoints(plot = p4, 
                points = rownames(adipocyte_top),repel = T)




Idents(luminal_AVs)<-''
p1<-VlnPlot(luminal_AVs, features = rownames(luminal_AV_top),
            stack=T,pt.size=0,
            flip = T,
            split.by = 'sample',
            split.plot = T)+
  theme(legend.position = 'none')


Idents(luminal_HSs)<-''
p2<-VlnPlot(luminal_HSs, features = rownames(luminal_HS_top),
            stack=T,pt.size=0,
            flip = T,
            split.by = 'sample',
            split.plot = T)+
  theme(legend.position = 'none')


Idents(basals)<-''
p3<-VlnPlot(basals, features = rownames(basal_top),
            stack=T,pt.size=0,
            flip = T,
            split.by = 'sample',
            split.plot = T)+
  theme(legend.position = 'none')



Idents(adipocytes)<-''
p4<-VlnPlot(adipocytes, features = rownames(adipocyte_top),
            stack=T,pt.size=0,
            flip = T,
            split.by = 'sample',
            split.plot = T)+
  theme(legend.position = 'none')




hall_markers <-
  msigdbr(species = 'Mus musculus',category = 'C2')
fgsea_sets<- hall_markers %>% 
  split(x = .$gene_symbol, f = .$gs_name)




luminal_AV_markers$genes<-rownames(luminal_AV_markers)
luminal_AV_genes <- luminal_AV_markers %>% 
  arrange(desc(avg_log2FC)) %>% 
  filter(p_val_adj < 0.05) %>%
  select(genes,avg_log2FC)
luminal_AV_ranks <- deframe(luminal_AV_genes)


luminal_HS_markers$genes<-rownames(luminal_HS_markers)
luminal_HS_genes <- luminal_HS_markers %>% 
  arrange(desc(avg_log2FC)) %>% 
  filter(p_val_adj < 0.05) %>%
  select(genes,avg_log2FC)
luminal_HS_ranks <- deframe(luminal_HS_genes)


fgseaRes1 <- fgsea(pathways = fgsea_sets, 
                   stats = luminal_HS_ranks,
                   minSize=10,
                   maxSize=500)


fgseaRes2 <- fgsea(pathways = fgsea_sets, 
                   stats = luminal_AV_ranks,
                   minSize=10,
                   maxSize=500)




basal_markers$genes<-rownames(basal_markers)
basal_genes <- basal_markers %>% 
  arrange(desc(avg_log2FC)) %>% 
  filter(p_val_adj < 0.05) %>%
  select(genes,avg_log2FC)
basal_ranks <- deframe(basal_genes)


fgseaRes3 <- fgsea(pathways = fgsea_sets, 
                   stats = basal_ranks,
                   minSize=10,
                   maxSize=500)


adipocyte_markers$genes<-rownames(adipocyte_markers)
adipocyte_genes <- adipocyte_markers %>% 
  arrange(desc(avg_log2FC)) %>% 
  filter(p_val_adj < 0.05) %>%
  select(genes,avg_log2FC)
adipocyte_ranks <- deframe(adipocyte_genes)


fgseaRes4 <- fgsea(pathways = fgsea_sets, 
                   stats = adipocyte_ranks,
                   minSize=10,
                   maxSize=500)




fgseaRes1 <- fgseaRes1 %>%
  filter(padj < 0.05)
fgseaRes2 <- fgseaRes2 %>%
  filter(padj < 0.05)
fgseaRes3 <- fgseaRes3 %>%
  filter(padj < 0.05)
fgseaRes4 <- fgseaRes4 %>%
  filter(padj < 0.05) %>%
  arrange(desc(NES))




p2<-plotGseaTable(fgsea_sets[c(topup$pathway,rev(fgseaRes1[NES<0,pathway]))],
                  luminal_HS_ranks,fgseaRes1,gseaParam = 0.5)




topup<-fgseaRes2 %>%
  arrange(desc(NES)) %>%
  head(10)
topdown<-fgseaRes2 %>%
  arrange(desc(NES)) %>%
  tail(10)


p1<-plotGseaTable(fgsea_sets[c(topup$pathway,topdown$pathway)],
                  luminal_AV_ranks,fgseaRes2,gseaParam = 0.5)


p3<-plotGseaTable(fgsea_sets[fgseaRes4$pathway],
                  adipocyte_ranks,fgseaRes4,gseaParam = 0.5)
