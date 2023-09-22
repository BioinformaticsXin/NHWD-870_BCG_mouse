#use CodeLib
options(encoding = "UTF-8")      ## for chines
## use getOption("encoding") to see if things were changed

project_dir <- "/data/home/lixin/Project/870_Mouse"

library(Seurat)
library(ggplot2)
library(ggthemes)

M_870 <- readRDS(paste0(project_dir,"/2_cell_annotation/6_annotation/M_870.rds"))

T_cell <- subset(M_870, idents = c('T cell/NK'));
all.genes <- rownames(T_cell)
T_cell <- ScaleData(T_cell, features = all.genes)
T_cell <- RunPCA(T_cell, features = VariableFeatures(object = T_cell))
T_cell <- FindNeighbors(T_cell)
T_cell <- FindClusters(T_cell, resolution = 0.3)
T_cell <- RunTSNE(T_cell, dims = 1:50)
T_cell <- RunUMAP(T_cell, dims = 1:50, label = T)
p <- FeaturePlot(T_cell, reduction="umap",features = c("Cd4","Cd8a","Nkg7","Ncr1"),min.cutoff = 0, max.cutoff = 4,ncol=2,cols = c("#DCDCDC","red")) + coord_fixed()
ggsave(p, file="/data/home/lixin/Project/870_Mouse/5_T_cell/T_feature.pdf")

p <- DimPlot(T_cell, reduction = "umap",label = TRUE, pt.size = 1,cols=c("#387FB9","#F91924","#4EB14B","#9D4CAA","#FC7F03","#A75728","#92caeb")) + coord_fixed()
ggsave(p, file="/data/home/lixin/Project/870_Mouse/5_T_cell/T_dim.pdf")

T_cell <- Seurat_marker(seurat_object = T_cell, resolution = 0.4, cores = 20, top_n = 45, out_dir = paste0(project_dir,"/5_T_cell"))


T_cell <- subset(T_cell, idents=c("0","1","2","4","7"))

cluster_to_cell <- list("0" = "CD8 T cell", "1" = "NK", "2" = "Cycling_CD8 T cell", "4" = "CD4 T cell", "7" = "NK")
cell_color <- c("#4EB14B","#c06ace","#fc952e","#92caeb")

marker_gene <- c("Cd8a","Nkg7","Ncr1","Klrd1","Xcl1","Top2a","Mki67","Cd4")
marker_gene <- unlist(marker_gene)
new.cluster.ids <- unlist(as.character(cluster_to_cell))
names(new.cluster.ids) <- names(cluster_to_cell)
T_cell <- RenameIdents(T_cell, new.cluster.ids)
T_cell@meta.data$cell.type <- Idents(T_cell)

levels(T_cell) <- c("CD4 T cell", "CD8 T cell", "Cycling_CD8 T cell", "NK")
p <- DimPlot(T_cell, reduction = "umap",label = TRUE, pt.size = 0.8,group.by="cell.type",cols=c("CD4 T cell"="#4EB14B","CD8 T cell"="#c06ace","Cycling_CD8 T cell"="#fc952e","NK"="#92caeb"))
ggsave(p, file="/data/home/lixin/Project/870_Mouse/5_T_cell/omit_double/T_dim.pdf", width=5.6, height=4)

p <- StackedVlnPlot(obj = T_cell, features = marker_gene,plot.margin = unit(c(-0.75, 0, -0.75, 0)),cols=c("CD4 T cell"="#4EB14B","CD8 T cell"="#c06ace","Cycling_CD8 T cell"="#fc952e","NK"="#92caeb"),ncol=1)
ggsave(p, file="/data/home/lixin/Project/870_Mouse/5_T_cell/omit_double/VlnPlot.pdf", width=3, height=6)

saveRDS(T_cell, file = paste0(project_dir, "/5_T_cell/omit_double/T_cell.rds"))


T_cell@meta.data$cell.type <- Idents(T_cell)
samples <- unique(T_cell@meta.data$SampleID)
Cell_sta <- c()
for(i in 1:length(samples)){
   this.meta <- T_cell@meta.data[which(T_cell@meta.data$SampleID == samples[i]),]
   this.sta <- data.frame(table(this.meta$cell.type), as.numeric(table(this.meta$cell.type)/sum(table(this.meta$cell.type))), SampleID = rep(samples[i], length(table(this.meta$cell.type))), stringsAsFactors=FALSE)
   Cell_sta <- rbind(Cell_sta, this.sta)
}
colnames(Cell_sta) <- c("Cell", "Number", "Proportion", "SampleID")

library(ggplot2)
library(ggthemes)
Cell_sta$SampleID <- factor(Cell_sta$SampleID, levels=rev(samples))
p1 <- ggplot(Cell_sta, aes(x=SampleID, y=Proportion,fill=Cell))+
geom_bar(stat="identity")+scale_fill_manual(values=c("CD4 T cell"="#4EB14B","CD8 T cell"="#c06ace","Cycling_CD8 T cell"="#fc952e","NK"="#92caeb"))+theme_classic();
ggsave(p1, file = paste0(project_dir, "/5_T_cell/omit_double/bar_cell_pro.pdf"), width=5,height=8)


T_cell <- readRDS(paste0(project_dir, "/5_T_cell/omit_double/T_cell.rds"))
all.genes <- rownames(T_cell)
T_cell <- ScaleData(T_cell, features = all.genes)
T_cell <- RunPCA(T_cell, features = VariableFeatures(object = T_cell))
T_cell <- FindNeighbors(T_cell)
T_cell <- FindClusters(T_cell, resolution = 0.3)
T_cell <- RunTSNE(T_cell, dims = 1:50)
T_cell <- RunUMAP(T_cell, dims = 1:50, label = T)

p <- DimPlot(T_cell, reduction = "umap",label = TRUE, pt.size = 0.8,group.by="cell.type",cols=c("CD4 T cell"="#4EB14B","CD8 T cell"="#c06ace","Cycling_CD8 T cell"="#fc952e","NK"="#92caeb"))
ggsave(p, file="/data/home/lixin/Project/870_Mouse/5_T_cell/omit_double/T_dim.pdf", width=5.6, height=4)
