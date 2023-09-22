library(Seurat)
library(ggplot2)
library(ggthemes)

project_dir <- "/data/home/lixin/Project/870_Mouse"
M_870 <- readRDS(paste0(project_dir, "/2_cell_annotation/6_annotation/update/M_870.rds"))
#count
counts <- as.matrix(GetAssayData(object = M_870, slot = "counts"));

#gene
gene <- read.table("/data/apps/infercnv-master/annotation_file/Mus_musculus.GRCm38.99_gene_pos.txt",header=F,sep="\t",stringsAsFactors=F);
gene <- unique(gene,);
gene <- gene[gene[,1] %in% rownames(counts),]
counts <- counts[gene[,1],];

#annotation
annotation <- data.frame(id=rownames(M_870@meta.data),SampleType=M_870@meta.data$cell.type,stringsAsFactors=F)
rownames(annotation) <- annotation[,1];
annotation <- annotation[colnames(counts),];
write.table(annotation, paste0(project_dir,'/infercnv/annotation.txt'), sep = '\t', row.names = F, col.names = F, quote = F)
write.table(gene, paste0(project_dir,'/infercnv/gene.txt'), sep = '\t', row.names = F, col.names = F, quote = F)
write.table(counts, paste0(project_dir,'/infercnv/counts.txt'), sep = '\t', row.names = T, col.names = T, quote = F)



#run infercnv
library(infercnv);
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste0(project_dir,"/infercnv/counts.txt"), 
annotations_file=paste0(project_dir,"/infercnv/annotation.txt"), delim="\t", 
gene_order_file=paste0(project_dir,"/infercnv/gene.txt"), ref_group_names=c("Monocyte/Macrophage","T cell/NK", "Fibroblast", "Endothelial"))
for(i in 2:5){
    infercnv_re = infercnv::run(infercnv_obj, cutoff=0.1, out_dir=paste0(project_dir,"/infercnv/cluster",i), 
    cluster_by_groups=F, k_obs_groups=i,denoise=F, HMM=T,output_format="pdf", num_threads=256)
}

