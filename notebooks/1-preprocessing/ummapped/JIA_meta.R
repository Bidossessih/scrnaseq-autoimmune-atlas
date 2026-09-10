
path="/nfs/team205/bh14/Datasets/Remapped/raw_adata/ummapped/"


jia_path="/lustre/scratch117/casm/team274/my4/oldScratch/Temp/Arthritis/"

jia = readRDS(paste0(jia_path, "Paediatric_arthritis_UCL_all_channel_EmptyDrops_CellRangers_pass_cells_count_matrix_filt_NONscaled.RDS"))

jia = Seurat::UpdateSeuratObject(object = jia)

jia_meta = jia@meta.data

table(jia_meta$tissue_type)

jia_meta = jia_meta[,c(5:8, 14)]

colnames(jia_meta) = c("Sample.ID", "donor_id", "tissue", "cell_source", "sex")
write.csv(jia_meta, file = paste0(path, "JIA_Martin.csv" ))
