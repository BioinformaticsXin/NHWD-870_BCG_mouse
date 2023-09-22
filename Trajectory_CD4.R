#use CodeLib
options(encoding = "UTF-8")      ## for chines
## use getOption("encoding") to see if things were changed


project_dir <- "/data/home/lixin/Project/870_Mouse"
T_cell <- readRDS(paste0(project_dir,"/5_T_cell/omit_double/T_cell.rds"))
T_cell <- subset(T_cell,ident=c("CD4 T cell"));
Idents(object = T_cell) <- "SampleType";


library(RColorBrewer)
#T_cell_cds <- trajectory_analysis(seurat_object = T_cell, expressionFamily = negbinomial.size(), gene_method = "cluster", out_dir = paste0(project_dir, "/5_T_cell/Trajectory"))
T_cell_cds <- trajectory_analysis(seurat_object = T_cell, expressionFamily = negbinomial.size(), gene_method = "variable", cell_cluster_name = "SampleType", nfeature = 1500, top_n = 100, q_cutoff = 0.05, out_dir = paste0(project_dir, "/5_T_cell/Trajectory/CD4"))
T_cell_cds$cell.type <- factor(T_cell_cds$cell.type, levels=c("CD4 T cell","CD8 T cell","Cycling_CD8 T cell"))
T_cell_cds$Phase <- factor(T_cell_cds$Phase, levels=c("G1","S","G2M"))
T_cell_cds$SampleType <- factor(T_cell_cds$SampleType, levels=c("CON", "870", "870_BCG"))


plot1 <- plot_cell_trajectory(T_cell_cds, color_by = "State", cell_size = 1, show_branch_points=FALSE)+theme(legend.position = "​right") 

plot2 <- plot_cell_trajectory(T_cell_cds, color_by = "SampleType", cell_size = 1, show_branch_points=FALSE)+scale_colour_manual(values = brewer.pal(5,"Paired"))+theme(legend.position = "​right")

plotc <- plot1|plot2
ggsave(paste0(project_dir, "/5_T_cell/Trajectory/CD4", "/T_cell_Combination.pdf"), plot = plotc, width = 10, height = 5)

p1 <- plot_cell_trajectory(T_cell_cds, color_by = "State", cell_size = 1,show_branch_points=TRUE) + facet_wrap(~State, nrow = 1)
p2 <- plot_cell_trajectory(T_cell_cds, color_by = "SampleType", cell_size = 1,show_branch_points=FALSE) + facet_wrap(~SampleType, nrow = 1)+scale_colour_manual(values = brewer.pal(5,"Paired"))
p3 <- plot_cell_trajectory(T_cell_cds, color_by = "Phase", cell_size = 1,show_branch_points=FALSE) + facet_wrap(~Phase, nrow = 1)+scale_colour_manual(values = c("G1"="#66C2A5","G2M"="#FC8D62","S"="#8DA0CB"))
p4 <- plot_cell_trajectory(T_cell_cds, color_by = "cell.type", cell_size = 1,show_branch_points=FALSE) + facet_wrap(~cell.type, nrow = 1)+scale_colour_manual(values = c("CD4 T cell"="#4EB14B","CD8 T cell"="#c06ace","Cycling_CD8 T cell"="#fc952e","NK"="#92caeb"))
p5 <- plot_cell_trajectory(T_cell_cds, color_by = "Pseudotime") + scale_color_gradientn(colours=c("darkblue", "darkred", "orange"))
p6 <- plot_cell_trajectory(T_cell_cds, color_by = "SampleType", cell_size = 1,show_branch_points=FALSE)+scale_colour_manual(values = brewer.pal(5,"Set1"))

layout <- matrix(c(rep(c(1,2,3), each = 3), c(4, 5, 6)), nrow = 4, byrow = TRUE)
pdf(paste0(project_dir, "/5_T_cell/Trajectory/CD4", "/T_cell_trajectory_facet.pdf"), width = 10, height = 16)
multiplot(plotlist = list(p1, p2, p3, p5, p6, p4), layout = layout)
#ggsave(paste0(project_dir, "/5_T_cell/Trajectory/CD8", "/T_cell_trajectory_facet.pdf"), plot = plotc, width = 10, height = 16)
dev.off()

GM_state <- function(cds){
    if (length(unique(pData(cds)$State)) > 1){
        T0_counts <- table(pData(cds)$State, pData(cds)$SampleType)[,"0"]
        return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
    } else {
        return (1)
    }
}
T_cell_cds <- orderCells(T_cell_cds, root_state = 5)

p1 <- plot_cell_trajectory(T_cell_cds, color_by = "State", cell_size = 1,show_branch_points=TRUE) + facet_wrap(~State, nrow = 1)
p2 <- plot_cell_trajectory(T_cell_cds, color_by = "SampleType", cell_size = 1,show_branch_points=FALSE) + facet_wrap(~SampleType, nrow = 1)+scale_colour_manual(values = brewer.pal(5,"Set1"))
p3 <- plot_cell_trajectory(T_cell_cds, color_by = "Phase", cell_size = 1,show_branch_points=FALSE) + facet_wrap(~Phase, nrow = 1)+scale_colour_manual(values = c("G1"="#66C2A5","G2M"="#FC8D62","S"="#8DA0CB"))
p4 <- plot_cell_trajectory(T_cell_cds, color_by = "cell.type", cell_size = 1,show_branch_points=FALSE) + facet_wrap(~cell.type, nrow = 1)+scale_colour_manual(values = c("CD4 T cell"="#4EB14B","CD8 T cell"="#c06ace","Cycling_CD8 T cell"="#fc952e","NK"="#92caeb"))
p5 <- plot_cell_trajectory(T_cell_cds, color_by = "Pseudotime") + scale_color_gradientn(colours=c("darkblue", "darkred", "orange"))
p6 <- plot_cell_trajectory(T_cell_cds, color_by = "SampleType", cell_size = 1,show_branch_points=FALSE)+scale_colour_manual(values = brewer.pal(5,"Set1"))
#plotc <- p1/p2/p3/(p5|p4)
layout <- matrix(c(rep(c(1,2,3), each = 3), c(4, 5, 6)), nrow = 4, byrow = TRUE)
pdf(paste0(project_dir, "/5_T_cell/Trajectory/CD4", "/T_cell_trajectory_facet.pdf"), width = 10, height = 16)
multiplot(plotlist = list(p1, p2, p3, p5, p6, p4), layout = layout)
#ggsave(paste0(project_dir, "/5_T_cell/Trajectory/CD8", "/T_cell_trajectory_facet.pdf"), plot = plotc, width = 10, height = 16)
dev.off()


p <- plot_cell_trajectory(T_cell_cds, markers = c("Top2a"), use_color_gradient = TRUE, cell_size = 1) + facet_wrap(~cell.type, nrow = 1) +scale_color_gradient2 (low="#3581bc", mid="#f6b1a5", high="#db1916")
ggsave(paste0(project_dir, "/5_T_cell/Trajectory/CD4", "/Top2a.pdf"), plot = p, width = 3, height = 4)

p <- plot_cell_trajectory(T_cell_cds, markers = c("Mki67"), use_color_gradient = TRUE, cell_size = 1) + facet_wrap(~cell.type, nrow = 1) +scale_color_gradient2 (low="#3581bc", mid="#f6b1a5", high="#db1916")
ggsave(paste0(project_dir, "/5_T_cell/Trajectory/CD4", "/Mki67.pdf"), plot = p, width = 3, height = 4)

