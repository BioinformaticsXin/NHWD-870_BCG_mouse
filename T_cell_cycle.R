#use CodeLib
options(encoding = "UTF-8")      ## for chines
## use getOption("encoding") to see if things were changed

project_dir <- "/data/home/lixin/Project/870_Mouse"

T_cell <- readRDS(paste0(project_dir,"/5_T_cell/omit_double/T_cell_update.rds"))

out_dir <- paste0(project_dir, "/5_T_cell/cell_cycle")
library(RColorBrewer)
p1 <- DimPlot(T_cell, reduction = "umap", label = TRUE, pt.size = .8, cols =brewer.pal(5,"Set2"), group.by="Phase")
ggsave(p1, file = paste0(out_dir, "/umap_Phase.pdf"),height=4,width=4.3)

p1 <- DimPlot(T_cell, reduction = "umap", label = TRUE, pt.size = .8, cols =brewer.pal(5,"Set2"), group.by="Phase", split.by="SampleType")
ggsave(p1, file = paste0(out_dir, "/umap_Phase_SampleType.pdf"),height=5,width=14)


cells <- as.character(unique(T_cell@meta.data$cell.type))
samples <- unique(T_cell@meta.data$SampleID)
cycling_meta <- T_cell@meta.data
Phase_sta <- c()
for(j in 1:length(cells)){
 for(i in 1:length(samples)){
    this.meta <- cycling_meta[which(cycling_meta$SampleID == samples[i] & cycling_meta$cell.type == cells[j]),]
    if(dim(this.meta)[1] > 0){
       this.sta <- data.frame(table(this.meta$Phase), as.numeric(table(this.meta$Phase)/sum(table(this.meta$Phase))), SampleID = rep(samples[i], length(table(this.meta$Phase))), Cell = rep(cells[j], length(table(this.meta$Phase))), stringsAsFactors=FALSE)
       Phase_sta <- rbind(Phase_sta, this.sta)
    }
    
 }
}

colnames(Phase_sta) <- c("Phase", "Number", "Proportion", "SampleID", "Cell")


library(ggplot2)
library(ggthemes)
Phase_sta$SampleID <- factor(Phase_sta$SampleID, levels=rev(samples))
p1 <- ggplot(Phase_sta, aes(x=SampleID, y=Proportion,fill=Phase))+
geom_bar(stat="identity")+scale_fill_manual(values=c("G1"="#66C2A5","G2M"="#FC8D62","S"="#8DA0CB"))+theme_classic()+facet_wrap(~ Cell, ncol=4);
ggsave(p1, file = paste0(out_dir, "/bar_Phase_pro.pdf"), width=10,height=5)
