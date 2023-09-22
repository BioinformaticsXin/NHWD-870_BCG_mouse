#use CodeLib
options(encoding = "UTF-8")      ## for chines
## use getOption("encoding") to see if things were changed

project_dir <- "/hengya/home/lixin/Project/870_Mouse"
T_cell <- readRDS(paste0(project_dir,"/5_T_cell/omit_double/T_cell_update.rds"))
CD4_T_cell <- subset(T_cell,ident=c("CD4 T cell"));
Idents(object = CD4_T_cell) <- "SampleType";
levels(CD4_T_cell) <- c("CON", "870", "870_BCG")


dif_BCG_CON <- Dif_cluster_enrich(seurat_object = CD4_T_cell,clustera = "870_BCG",clusterb = "CON",group.by = NULL, out_dir = "/hengya/home/lixin/Project/870_Mouse/5_T_cell/function_anno/CD4/BCG_CON", OrgDb = "org.Mm.eg.db", sig_gene_p = 0.05, sig_gene_fd = 1, organism = "mmu", cut_p = "p_val")


library(dplyr)
dif_BCG_CON.GO_up <- select(dif_BCG_CON$GO_up, c("ONTOLOGY","ID","Description","pvalue","Count","GeneRatio"))
dif_BCG_CON.GO_up[,1] <- paste0("GO.",dif_BCG_CON.GO_up[,1])

dif_BCG_CON.KEGG_up <- select(dif_BCG_CON$KEGG_up, c("ID","Description","pvalue","Count","GeneRatio"))
dif_BCG_CON.KEGG_up <- dif_BCG_CON.KEGG_up[1:6,]
dif_BCG_CON.KEGG_up <- data.frame(ONTOLOGY=rep("KEGG",nrow(dif_BCG_CON.KEGG_up)), dif_BCG_CON.KEGG_up, stringsAsFactors=F)

dif_BCG_CON.GO_up <- rbind(dif_BCG_CON.GO_up,dif_BCG_CON.KEGG_up)


GeneRatio <- function(x){
   a <- as.numeric(unlist(strsplit(x, split="\\/")))[1]
   b <- as.numeric(unlist(strsplit(x, split="\\/")))[2]
   return(a/b)
}
dif_BCG_CON.GO_up$GeneRatio <- GeneRatio(dif_BCG_CON.GO_up$GeneRatio)

library(ggplot2)
library(ggthemes)
dif_BCG_CON.GO_up$Count <- as.character(dif_BCG_CON.GO_up$Count)
p <- ggplot(dif_BCG_CON.GO_up, aes(x=Description,y=Count))+geom_point(aes(color=-log10(pvalue),size=GeneRatio))+
coord_flip()+scale_color_gradientn(colors=c("#5D2265","#0095A2","#FCE557"))+theme_bw()+
facet_grid(ONTOLOGY ~ ., scales = "free", space='free') 
ggsave(p, file="/hengya/home/lixin/Project/870_Mouse/5_T_cell/function_anno/CD4/BCG_CON/up_bubble_plot.pdf", width=7, height=6)

#KEGG down
dif_BCG_CON.KEGG_down <- select(dif_BCG_CON$KEGG_down, c("ID","Description","pvalue","Count","GeneRatio"))
dif_BCG_CON.KEGG_down <- data.frame(ONTOLOGY=rep("KEGG",nrow(dif_BCG_CON.KEGG_down)), dif_BCG_CON.KEGG_down, stringsAsFactors=F)
dif_BCG_CON.KEGG_down$GeneRatio <- GeneRatio(dif_BCG_CON.KEGG_down$GeneRatio)

library(ggplot2)
library(ggthemes)
dif_BCG_CON.KEGG_down$Count <- as.character(dif_BCG_CON.KEGG_down$Count)
p <- ggplot(dif_BCG_CON.KEGG_down, aes(x=Description,y=Count))+geom_point(aes(color=-log10(pvalue),size=GeneRatio))+
coord_flip()+scale_color_gradientn(colors=c("#5D2265","#0095A2","#FCE557"))+theme_bw()+
facet_grid(ONTOLOGY ~ ., scales = "free", space='free') 
ggsave(p, file="/hengya/home/lixin/Project/870_Mouse/5_T_cell/function_anno/CD4/BCG_CON/down_bubble_plot.pdf", width=5, height=3)


library(openxlsx)
dif_BCG_CON_gene_up <- read.xlsx("/hengya/home/lixin/Project/870_Mouse/5_T_cell/function_anno/BCG_CON/dif_gene.xlsx", sheet=1)
dif_BCG_CON_gene_down <- read.xlsx("/hengya/home/lixin/Project/870_Mouse/5_T_cell/function_anno/BCG_CON/dif_gene.xlsx", sheet=2)

library(Seurat)
library(RColorBrewer)
p <- DoHeatmap(subset(CD4_T_cell, ident=c("CON","870_BCG")), features = c(dif_BCG_CON_gene_up$SYMBOL, dif_BCG_CON_gene_down$SYMBOL), size = 3.5,group.colors = brewer.pal(3,"Paired")[2:3]) + 
scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
midpoint = 0, guide = "colourbar", aesthetics = "fill") + theme(axis.text.y = element_text(size = 4))
ggsave(p, file="/hengya/home/lixin/Project/870_Mouse/5_T_cell/function_anno/BCG_CON/heatmap.pdf", width=4, height=4)


