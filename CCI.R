#use CodeLib
options(encoding = "UTF-8")      ## for chines
## use getOption("encoding") to see if things were changed


project_dir <- "/hengya/home/lixin/Project/870_Mouse"
M_870 <- readRDS("/hengya/home/lixin/Project/870_Mouse/2_cell_annotation/6_annotation/update/M_870.rds")
YUMM1.7 <- subset(M_870, idents = c('YUMM1.7'))
T_cell <- readRDS(paste0(project_dir,"/5_T_cell/omit_double/T_cell.rds"))
T_cell <- subset(T_cell, idents = c('CD4 T cell', 'CD8 T cell', 'Cycling_CD8 T cell'))
YUMM1.7_T <- Seurat_merge(seurat_object_list <- list(YUMM1.7=YUMM1.7,T_cell=T_cell))

Idents(object = YUMM1.7_T) <- "SampleType";
YUMM1.7_T_CON <- subset(YUMM1.7_T, idents = c('CON'));
YUMM1.7_T_870 <- subset(YUMM1.7_T, idents = c('870'));
YUMM1.7_T_870_BCG <- subset(YUMM1.7_T, idents = c('870_BCG'));

cellchat_CON <- CellChat_analysis(seurat_object = YUMM1.7_T_CON, seurat_annotations = "cell.type", CellChatDB = "CellChatDB.mouse",
cores = 20, search_subsetDB = NULL, PPIdata = "PPI.mouse", raw.use = FALSE,
out_dir = "/hengya/home/lixin/Project/870_Mouse/6_YUMM1.7/CellChat/CON")

CellChat_pathway_plot(cellchat = cellchat_CON, pathway = "CXCL", vertex.receiver = seq(1,2), out_dir = "/hengya/home/lixin/Project/870_Mouse/6_YUMM1.7/CellChat/CON")
CellChat_pathway_plot(cellchat = cellchat_CON, pathway = "CCL", vertex.receiver = seq(1,2), out_dir = "/hengya/home/lixin/Project/870_Mouse/6_YUMM1.7/CellChat/CON")
p <- plotGeneExpression(cellchat_CON, features = c("Ccl3", "Ccl4", "Ccl5", "Ccl8", "Ccr5"), type="dot")
ggsave(p, file=paste0(project_dir, "/6_YUMM1.7/CellChat/CON", "/VlnPlot_CCL_CON.pdf"), width=6, height=3)


cellchat_870 <- CellChat_analysis(seurat_object = YUMM1.7_T_870, seurat_annotations = "cell.type", CellChatDB = "CellChatDB.mouse",
cores = 20, search_subsetDB = NULL, PPIdata = "PPI.mouse", raw.use = FALSE,
out_dir = "/hengya/home/lixin/Project/870_Mouse/6_YUMM1.7/CellChat/870")

CellChat_pathway_plot(cellchat = cellchat_870, pathway = "CXCL", vertex.receiver = seq(1,2), out_dir = "/hengya/home/lixin/Project/870_Mouse/6_YUMM1.7/CellChat/870")
CellChat_pathway_plot(cellchat = cellchat_870, pathway = "CCL", vertex.receiver = seq(1,2), out_dir = "/hengya/home/lixin/Project/870_Mouse/6_YUMM1.7/CellChat/870")
p <- plotGeneExpression(cellchat_870, features = c("Ccl3", "Ccl4", "Ccl5", "Ccl8", "Ccr5"), type="dot")
ggsave(p, file=paste0(project_dir, "/6_YUMM1.7/CellChat/870", "/VlnPlot_CCL_870.pdf"), width=2, height=4)



cellchat_870_BCG <- CellChat_analysis(seurat_object = YUMM1.7_T_870_BCG, seurat_annotations = "cell.type", CellChatDB = "CellChatDB.mouse",
cores = 20, search_subsetDB = NULL, PPIdata = "PPI.mouse", raw.use = FALSE,
out_dir = "/hengya/home/lixin/Project/870_Mouse/6_YUMM1.7/CellChat/870_BCG")

CellChat_pathway_plot(cellchat = cellchat_870_BCG, pathway = "CXCL", vertex.receiver = seq(1,2), out_dir = "/hengya/home/lixin/Project/870_Mouse/6_YUMM1.7/CellChat/870_BCG")
CellChat_pathway_plot(cellchat = cellchat_870_BCG, pathway = "CCL", vertex.receiver = seq(1,2), out_dir = "/hengya/home/lixin/Project/870_Mouse/6_YUMM1.7/CellChat/870_BCG")


library(RColorBrewer)
object.list <- list(cellchat_CON = cellchat_CON, cellchat_870 = cellchat_870, cellchat_870_BCG = cellchat_870_BCG)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("cellchat_CON", "cellchat_870", "cellchat_870_BCG")) # set factor level
p <- plotGeneExpression(cellchat, features = c("Ccl3", "Ccl4", "Ccl5", "Ccl8", "Ccr5", "Cxcl10", "Cxcr3"), split.by = "datasets", color.use = brewer.pal(5,"Paired"))
ggsave(p, file=paste0(project_dir, "/6_YUMM1.7/CellChat", "/VlnPlot_CCL_CXCL_merge.pdf"), width=5.5, height=4.5)


#############bubble plot############
pdf(paste0(project_dir, "/CCI/merge_cell/MK_RC.pdf"), width = 10, height = 3)
netVisual_bubble(cellchat.two, sources.use = 16, targets.use = 1:15, signaling = c("MK"), remove.isolate = FALSE)#> Comparing communications on a single object
dev.off()

