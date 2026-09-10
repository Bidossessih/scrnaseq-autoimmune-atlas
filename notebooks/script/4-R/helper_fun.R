
#' Identify DE genes from two model
#' DE genes from fit2 (background, control for example) were substract from DE genes obtained in fit1 (disease for example)
#' @param alpha is the threshold of p-value
#' @param max_lfc is the maximum log FC. This parameter used to filter out genes that are low expressed
#' and that tend to have extreme log FC
#' @param min_lfc is the min log FC
#

filter_res_2fit = function(fit1, fit2, contrast1 = `nhoo_test_conditiontest_nhood`, 
                           contrast2 = `nhoo_test_conditiontest_nhood`, alpha=0.05, max_lfc = 15, min_lfc=1) {
  de_res_ <- test_de(fit1, contrast = eval(parse(text=contrast1)))
  
  des_res2 = test_de(fit2, contrast = eval(parse(text=contrast2)))
  
  de_res_$lfc2 = des_res2$lfc
  de_res_$adj_pval2 = des_res2$adj_pval
  
  de_res_ = de_res_ %>% filter(adj_pval < alpha) %>%
    filter((lfc >= min_lfc &
              lfc2 < 0) | (lfc <= -(min_lfc) & lfc2 > 0)) %>%
    arrange(desc(lfc)) %>% dplyr::select(1:7)
  
  de_res_$lfc <- ifelse(abs(de_res_$lfc) > max_lfc, sign(de_res_$lfc) * Inf, de_res_$lfc)
  
  return(de_res_)
  
}

#' Identify DE genes from substracting two contrast based on the same model 
#'

filter_res = function(fit, contr_test= "`nhoo_test_conditiontest_nhood`", contr_ctrl = "`disease_shortHealthy`",
                      alpha=0.05, max_lfc = 15){
  de_res <- test_de(fit, contrast = eval(parse(text=contr_test)))
  
  des_res2 = test_de(fit, contrast = eval(parse(text=contr_ctrl)))
  de_res$lfc2 = des_res2$lfc
  de_res$adj_pval2 = des_res2$adj_pval
  
  de_res = de_res %>% filter(adj_pval < alpha) %>%
    filter((lfc >= 1 &
              lfc2 < 0) | (lfc <= -1 & lfc2 > 0)) %>%
    arrange(desc(lfc)) %>% dplyr::select(1:7)
  
  de_res$lfc <- ifelse(abs(de_res$lfc) > max_lfc, sign(de_res$lfc) * Inf, de_res$lfc)
  
  return(de_res)
}


#' Test DE in cells in neighborhoods 
#' 
test_de_nhoo <- function(test_adata=test_adata, contrast='nhoo_test_conditionenriched',
                         nhoo_ref = 'depleted', alpha_level = 0.1, min_lfc = log2(1.5)){
  
  
  # Prepare disease data
  sce_dis = test_adata[, test_adata$disease_short!="Healthy"]
  
  # Some donors appears in differents nhoods with opposite annotations (depleted in one and enriched in  other)
  # We remove these donors from the dataset
  
  if(any(duplicated(colnames(sce_dis)))){
    
    dup = colnames(sce_dis)[duplicated(colnames(sce_dis))]
    
    print(str_glue('{length(dup)} duplited samples have been detected and removed from sce_dis.'))
    
    sce_dis=sce_dis[,!(colnames(sce_dis) %in% dup)]
  }
  
  sce_dis$disease_short = factor(sce_dis$disease_short)
  sce_dis$Cell_type_scANVI = factor(sce_dis$Cell_type_scANVI)
  sce_dis$nhoo_test_condition = factor(sce_dis$nhoo_test_condition)
  sce_dis$donor_id = factor(sce_dis$donor_id)
  sce_dis$dataset_id = factor(sce_dis$dataset_id)
  
  print("")
  print("sce_dis")
  print(sce_dis)
  
  #Prepare ctrl data
  sce_ctrl = test_adata[, test_adata$disease_short == "Healthy"]
  
  # Some donors appears in differents nhoods with opposite annotations (depleted in one and enriched in  other)
  # We remove these donors from the dataset
  
  if(any(duplicated(colnames(sce_ctrl)))){
    
    dup2 = colnames(sce_ctrl)[duplicated(colnames(sce_ctrl))]
    
    print(str_glue('{length(dup2)} duplited samples have been detected and removed from sce_ctrl.'))
    
    sce_ctrl = sce_ctrl[,!duplicated(colnames(sce_ctrl))]
    sce_ctrl=sce_ctrl[,!(colnames(sce_ctrl) %in% dup)]
  }
  sce_ctrl$disease_short = factor(sce_ctrl$disease_short)
  sce_ctrl$Cell_type_scANVI = factor(sce_ctrl$Cell_type_scANVI)
  sce_ctrl$nhoo_test_condition = factor(sce_ctrl$nhoo_test_condition)
  sce_ctrl$donor_id = factor(sce_ctrl$donor_id)
  sce_ctrl$dataset_id = factor(sce_ctrl$dataset_id)
  
  
  cat("#")
  cat("sce_ctrl")
  print(sce_ctrl)
  
  if(any(c(length(unique(sce_ctrl$nhoo_test_condition))==1, 
           length(unique(sce_dis$nhoo_test_condition))==1))){
    cat("lengh of nhoo_test_condition: \n")
    cat("\t -sce_ctrl \n")
    print(length(unique(sce_ctrl$nhoo_test_condition)))
    
    cat("\t -sce_dis \n")
    print(length(unique(sce_dis$nhoo_test_condition)))
    
    stop("nhoo_test_condition must be factors with 2 or more levels")
  }
  
  
  ## Fit
  my_design = ~ 1 + dataset_id + disease_short + nhoo_test_condition + nhoo_test_condition:disease_short
  
  # ctrl
  cat("#")
  cat("fitting ctrl....")
  # Start the clock!
  ptm <- proc.time()
  #matDS = model.matrix(my_design, data = colData(sce_ctrl))
  #my_design = ~ 0 + donor_id + nhoo_test_condition
  
  fit <- glm_gp(sce_ctrl, design = my_design, reference_level = nhoo_ref)
  de_res_ctrl <- test_de(fit, contrast = nhoo_test_conditionenriched)
  cat("#")
  print(proc.time() - ptm)
  
  # disease
  
  cat("#")
  cat("fitting disease....")
  # Start the clock!
  ptm <- proc.time()
  
  fit1 <- glm_gp(test_adata, design = my_design, reference_level = nhoo_ref)
  
  de_res_disease <- test_de(fit1, contrast = nhoo_test_conditionenriched)
  cat("#")
  print(proc.time() - ptm)
  
  
  ## Test result 
  de_res_disease$lfc_ctrl = de_res_ctrl$lfc
  de_res_disease$adj_pvalctrl = de_res_ctrl$adj_pval
  
  de_res_disease = de_res_disease %>% filter(adj_pval < alpha_level)
  de_res_disease = de_res_disease %>% filter((lfc>=min_lfc & lfc_ctrl < 0) | (lfc <= -min_lfc & lfc_ctrl > 0))
  
  return(de_res_disease)
}




#' This function prepare the test data for DE analysis and have as input a depleted 
#' data and an enriched data. It will be used to compare an enriched cell types 
#' a group of depleted cell types.

prep_de_data = function(depl_="", enrch_="", test_adata=NULL, nhood = FALSE, 
                        hvgs=TRUE, cel_of_int_vs_rest = FALSE, cell_type = NULL){
  #combine data
  if(depl_!="" & enrch_!="" & is.null(test_adata)){
    test_adata = cbind(depl_, enrch_)
  }
  
  counts(test_adata) <- assay(test_adata)
  
  test_adata = logNormCounts(test_adata)
  
  # genes where at least some counts are more than 10 count
  non_empty_rows <- which(rowSums2(assay(test_adata)) > 10)
  test_adata <- test_adata[non_empty_rows, ]
  ## Feature selection with scran. 10000 hvgs will be selected
  dec <- modelGeneVar(test_adata)
  
  if(hvgs){
    hvgs <- getTopHVGs(dec, n = 10000)
    test_adata <- test_adata[hvgs,]
  }
  
  
  counts(test_adata) <- as.matrix(counts(test_adata))
  
  assay(test_adata, "X") = as.matrix(assay(test_adata, "X"))
  
  test_adata$disease_short = factor(test_adata$disease_short) %>% droplevels()
  test_adata$Cell_type_scANVI = factor(test_adata$Cell_type_scANVI) %>% droplevels()
  test_adata$nhoo_test_condition = factor(test_adata$nhoo_test_condition) %>% droplevels()
  test_adata$donor_id = factor(test_adata$donor_id) %>% droplevels()
  
  test_adata$disease_bin = ifelse(test_adata$disease_short == "Healthy", "Healthy", "Disease")
  
  test_adata$disease_bin = factor(test_adata$disease_bin, levels = c("Healthy", "Disease"))
  
  if(nhood){
    test_adata$nhoo_test_condition = factor(test_adata$nhoo_test_condition, 
                                            levels = c("depleted", "enriched"))
    test_adata$nhoo_test_condition = droplevels(test_adata$nhoo_test_condition)
    
    print(table(test_adata$disease_short, test_adata$nhoo_test_condition))
  }
  
 
  
  test_adata$disease_short = droplevels(test_adata$disease_short)
  test_adata$disease_bin = droplevels(test_adata$disease_bin)
  test_adata$Cell_type_scANVI = droplevels(test_adata$Cell_type_scANVI)
  test_adata$donor_id = droplevels(test_adata$donor_id)
  
  
  #Remove duplicate
  dup = colnames(test_adata)[duplicated(colnames(test_adata))]
  test_adata = test_adata[, !(colnames(test_adata) %in% dup)]
  
  table(test_adata$disease_short, test_adata$nhoo_test_condition)
  
  test_adata$disease_short = droplevels(test_adata$disease_short)
  test_adata$donor_id = droplevels(test_adata$donor_id)
  test_adata$dataset_id = droplevels(test_adata$dataset_id)
  
  
  # update nhood
  # to show a specific cell type vs the rest
  
  if(cel_of_int_vs_rest & !is.null(cell_type)) {
    test_adata$nhoo_test_condition = ifelse(test_adata$Cell_type_scANVI == cell_type,
                                            "test_nhood",
                                            "rest")
  }
  
  
  return(test_adata)
}



#' This second function (prep_de_data_2) is used when we would like to compare a group 
#' of cell with all other cell. For example to compare a depleted nhood with cells 
#' within all other neighborhoods.
#' 
#' @param enr_vs_depl bool indicate whether the comparison will be performed between enriched and depleted.
#

prep_de_data_2 = function(depl_=test_adata_depleted, enrch_=test_adata_enriched, cell_type=NA,
                          dir_fc, enr_vs_depl, hvgs_=TRUE){
  
  dir_ = c(depleted="depleted", enriched="enriched")
  
  depl_$nhoo_test_condition_ori = depl_$nhoo_test_condition
  depl_$nhoo_test_condition = "depleted"
  
  enrch_$nhoo_test_condition_ori = enrch_$nhoo_test_condition
  enrch_$nhoo_test_condition = "enriched"
  
  #combine data
  test_adata = cbind(depl_, enrch_)
  
  print(test_adata)
  # Reannotate nhoo_test_condition
  # nhood of interest will be renamed 'test_nhood' and all other nhood will be renamed 'rest'
  if(enr_vs_depl & !is.na(cell_type)) {
    test_adata$nhoo_test_condition = ifelse(
      test_adata$Cell_type_scANVI == cell_type &
        test_adata$nhoo_test_condition_ori == "shared" &
        test_adata$nhoo_test_condition == dir_fc,
      "test_nhood",
      ifelse(
        test_adata$Cell_type_scANVI == cell_type &
          test_adata$nhoo_test_condition_ori == "shared" &
          test_adata$nhoo_test_condition == dir_[dir_ != dir_fc] %>% as.character(),
        "rest",
        "exclude"
      )
    )
    
    test_adata = test_adata[test_adata$nhoo_test_condition != "exclude"]
    
  } else if (is.na(cell_type) & enr_vs_depl) {
    test_adata = test_adata[, test_adata$nhoo_test_condition_ori == "shared"]
    print(test_adata)
  } else{
    test_adata$nhoo_test_condition = ifelse(
      test_adata$Cell_type_scANVI == cell_type &
        test_adata$nhoo_test_condition_ori == "shared" &
        test_adata$nhoo_test_condition == dir_fc,
      "test_nhood",
      "rest"
    )
  }
  
  test_adata$nhoo_test_condition = factor(test_adata$nhoo_test_condition)
  
  counts(test_adata) <- assay(test_adata)
  
  test_adata = logNormCounts(test_adata)
  
  
  if(hvgs_) {
    # genes where at least some counts are more than 10 count
    non_empty_rows <- which(rowSums2(assay(test_adata)) > 10)
    test_adata <- test_adata[non_empty_rows, ]
    ## Feature selection with scran. 10000 hvgs will be selected
    dec <- modelGeneVar(test_adata)
    hvgs <- getTopHVGs(dec, n = 5000)
    test_adata <- test_adata[hvgs, ]
  }
  
  counts(test_adata) <- as.matrix(counts(test_adata))
  
  assay(test_adata, "X") = as.matrix(assay(test_adata, "X"))
  
  ### drop empty unused levels, overwise glmGampoi will complain
  
  test_adata$disease_short = factor(test_adata$disease_short) %>% droplevels()
  test_adata$Cell_type_scANVI = factor(test_adata$Cell_type_scANVI) %>% droplevels()
  test_adata$nhoo_test_condition = factor(test_adata$nhoo_test_condition) %>% droplevels()
  test_adata$donor_id = factor(test_adata$donor_id) %>% droplevels()
  
  test_adata$disease_bin = ifelse(test_adata$disease_short == "Healthy", "Healthy", "Disease")
  
  test_adata$disease_bin = factor(test_adata$disease_bin, levels = c("Healthy", "Disease"))
  
  
  test_adata$disease_short = droplevels(test_adata$disease_short)
  test_adata$disease_bin = droplevels(test_adata$disease_bin)
  test_adata$Cell_type_scANVI = droplevels(test_adata$Cell_type_scANVI)
  test_adata$donor_id = droplevels(test_adata$donor_id)
  test_adata$nhoo_test_condition = droplevels(test_adata$nhoo_test_condition)
  
  print(table(test_adata$disease_short, test_adata$nhoo_test_condition))
  
  print("")
  
  print(table(test_adata$Cell_type_scANVI, test_adata$nhoo_test_condition))
  
  
  test_adata$disease_short = droplevels(test_adata$disease_short)
  test_adata$donor_id = droplevels(test_adata$donor_id)
  test_adata$dataset_id = droplevels(test_adata$dataset_id)
  
  
  return(test_adata)
}


## Prepare data for Heatmap plotting

prep_heatmap = function(depl_, enrch_, cell_type, dir_fc){
  
  depl_ = test_adata_depleted
  enrch_ = test_adata_enriched
  depl_$nhoo_test_condition_ori = depl_$nhoo_test_condition
  depl_$nhoo_test_condition = "depleted"
  
  enrch_$nhoo_test_condition_ori = enrch_$nhoo_test_condition
  enrch_$nhoo_test_condition = "enriched"
  
  #combine data
  test_adata = cbind(depl_, enrch_)
  
  # Reannotate nhoo_test_condition
  # nhood of interest will be renamed 'test_nhood' and all other nhood will be renamed 'rest'
  
  
  
  test_adata$nhoo_test_condition = ifelse(
    test_adata$Cell_type_scANVI == cell_type &
      test_adata$nhoo_test_condition_ori == "shared" &
      test_adata$nhoo_test_condition == dir_fc,
    "test_nhood",
    "rest"
  )
  
  test_adata$nhoo_test_condition = factor(test_adata$nhoo_test_condition)
  
  counts(test_adata) <- assay(test_adata)
  
  test_adata = logNormCounts(test_adata)
  
  test_adata
}






#' Function  DE from loop
#' 


run_de_2 = function(depl_, enrch_, cell_type, tissue_source, dir_fc, sce_path, enr_vs_depl=FALSE) {
  # merge test data
  test_adata = prep_de_data_2(depl_=depl_, enrch_=enrch_, cell_type=cell_type,
                              dir_fc=dir_fc, enr_vs_depl=enr_vs_depl)
  cat("\n")
  
  print("detecting DE in disease...")
  
  
  my_design = ~ 1 + disease_short + nhoo_test_condition
  
  
  test_adata$disease_short = factor(test_adata$disease_short)
  test_adata$disease_bin = factor(test_adata$disease_bin)
  test_adata$donor_id = factor(test_adata$donor_id)
  test_adata$dataset_id = factor(test_adata$dataset_id)
  test_adata$Cell_type_scANVI = str_replace(test_adata$Cell_type_scANVI, "/", "_") %>% factor()
  
  fit1 <- glm_gp(test_adata, design = my_design, reference_level = "rest")
  
  print(summary(fit1))
  cat("\n")
  print("Result...")
  
  
  de_res = filter_res(fit=fit1, contr_test= "`nhoo_test_conditiontest_nhood`", 
                      contr_ctrl = "`disease_shortHealthy`", alpha=0.1, max_lfc = 15)
  print("Top_10 up-regulated")
  print(de_res$name[1:10])
  
  print(str_glue("Up: {length(de_res$lfc[de_res$lfc > 0])}"))
  print(str_glue("Down: {length(de_res$lfc[de_res$lfc < 0])}"))
  
  dir_ = c(depleted="depleted", enriched="enriched")
  
  if(enr_vs_depl){
    ctrl_dir=str_glue("_{dir_[dir_ != dir_fc] %>% as.character()}")
  }else{
    ctrl_dir=""
  }
  
  write.csv(de_res, file = str_glue("{sce_path}Result_{tissue_source}_de_{cell_type}_{dir_fc}_vs_rest{ctrl_dir}.csv"))
  cat("\n")
  print(str_glue("Saved to {sce_path}Result_{tissue_source}_de_{cell_type}_{dir_fc}_vs_rest{ctrl_dir}.csv"))
  cat("\n")
}



#' Function to plot heatmap of DE
#'
# height 

plot_heatmap = function(test_adata,
                        deg=NA,
                        cel.selected,
                        selected_gene = NULL,
                        n_top = 15,
                        aggregate = TRUE,
                        aggr_only_nhoo = FALSE,
                        show_column_names = FALSE,
                        aggr_fun = "mean",
                        col_pal = vega_20,
                        test.column="Cell_type_scANVI",
                        column.label2 = "nhoo_test_condition",
                        col_test=c("grey", "grey50"),
                        scale = TRUE,
                        remove_inf = TRUE,
                        cluster.col = FALSE,
                        labels_annot_bar = TRUE,
                        label.only.disease = FALSE,
                        figure_Path,
                        save=FALSE,
                        fig.width=10,
                        fig.height=8,
                        grid_width=5) {
  
  if(is.data.frame(deg)) {
    # Remove Inf from deg dataset
    if (remove_inf) {
      deg = deg[!(deg$lfc %in% c(Inf,-Inf)), ]
    }
    # extract top genes
    up = deg %>% filter(lfc > 0) %>% arrange(desc(lfc))
    down = deg %>% filter(lfc < 0) %>% arrange(lfc)
    top_ = na.omit(rbind(up[1:n_top,], down[1:n_top,]))
  }
  
  #We subset test_data with the top DE genes
  #or a list of genes of interest (selected_gene)
  if(!is.null(selected_gene)) {
    data = test_adata[selected_gene,]
  } else{
    data = test_adata[top_$name,]
  }
  
  
  if(!aggr_only_nhoo){
    data$sample = paste0(data$disease_short, "_", data[[test.column]])
    data$sample = factor(data$sample)
    levels(data$sample)
  } else {
    data$sample = paste0(data$disease_bin, "_", data[[test.column]])
    data$sample = factor(data$sample)
    levels(data$sample)
  }
  
  
  if(aggregate){
    data = scuttle::aggregateAcrossCells(data,
                                         data$sample,
                                         use.assay.type = "logcounts",
                                         statistics = aggr_fun)
  }
  
  
  data$nhoo_test_condition = ifelse(data$nhoo_test_condition == "test_nhood",
                                    cel.selected,
                                    "Rest")
  
  # pheno data from sce
  pheno = as.data.frame(colData(data)) %>% arrange(eval(parse(text=test.column)))
  
  top_data = assay(data, "logcounts")
  
  top_data = top_data[,rownames(pheno)]
  
  #
  if (!all.equal(colnames(top_data), rownames(pheno))) {
    stop("Colnames of expression data and rownames of pheno don't match.")
  }
  
  
  data$disease_short = droplevels(data$disease_short)
  
  data$nhoo_test_condition = factor(data$nhoo_test_condition)
  
  # column annotation
  #disease
  dis_color = c()
  Disease = factor(test_adata$disease_short)
  
  for (i in 1:length(levels(Disease))) {
    dis = Disease[i]
    dis_color = c(dis_color, col_pal[i])
  }
  
  names(dis_color) = levels(Disease)
  # second annotation bar (por example test group)
  
  da_color = c()
  
  for (i in 1:length(unique(data[[column.label2]]))) {
    da = unique(data[[column.label2]])[i]
    da_color = c(da_color, col_test[i])
  }
  
  names(da_color) = unique(data[[column.label2]])
  
  # col annotation
  
  if(labels_annot_bar){
    labels_test_groups = as.character(1:length(pheno[[column.label2]]))#unique(pheno[[column.label2]]) %>% str_replace("Rest", "Other")
    
    
    ha_col = HeatmapAnnotation(foo = anno_block(align_to = 1:length(pheno[[column.label2]]), labels = labels_test_groups, labels_gp = gpar(col = "white")))
  }else{
    labels_test_groups = labels_disease= NULL
  }
  
  if(!label.only.disease) {
    ha_col = HeatmapAnnotation(
      df = data.frame(
        Test_group = factor(pheno[[column.label2]]),
        Disease = factor(pheno$disease_short)
      ),
      labels = labels_test_groups, 
      #labels_gp = gpar(col = "white", fontsize = 10),
      show_legend = c(TRUE, TRUE),
      show_annotation_name = F,
      col = list(
        Test_group = da_color,
        Disease = dis_color
      ),
      simple_anno_size = unit(0.25, "cm"),
      simple_anno_size_adjust = TRUE,
      annotation_legend_param = list(
        Test_group = list(
          title_gp = gpar(fontsize = 12, fontface = "bold"),
          labels_gp = gpar(fontsize = 12),
          grid_width = unit(5, "mm"),
          grid_height = unit(5, "mm")
        ),
        Disease = list(
          title_gp = gpar(fontsize = 12, fontface = "bold"),
          labels_gp = gpar(fontsize = 12),
          grid_width = unit(5, "mm"),
          grid_height = unit(5, "mm")
        )
      )
    )
  } else {
    ha_col = HeatmapAnnotation(
      df = data.frame(
        Disease = factor(pheno$disease_short)
      ),
      labels = labels_test_groups, 
      #labels_gp = gpar(col = "white", fontsize = 10),
      show_legend = c(TRUE, TRUE),
      show_annotation_name = F,
      col = list(
        Disease = dis_color
      ),
      simple_anno_size = unit(0.25, "cm"),
      simple_anno_size_adjust = TRUE,
      annotation_legend_param = list(
        Disease = list(
          title_gp = gpar(fontsize = 12, fontface = "bold"),
          labels_gp = gpar(fontsize = 12),
          grid_width = unit(5, "mm"),
          grid_height = unit(5, "mm")
        )
      )
    )
  }
  
  
  
  
  if (scale) {
    top_data = t(scale(t(top_data)))
  }
  
  
  # Included FC direction annotation if selected_gene is not provided
  if (is.null(selected_gene)) {
    top_$FC_dir = ifelse(top_$lfc > 0, "Up", "Down")
    
    top_$FC_dir = factor(top_$FC_dir)
    # Row annotation
    
    ha_row = HeatmapAnnotation(
      df = data.frame(FC_dir = top_$FC_dir),
      show_annotation_name = F,
      col = list(FC_dir = c(Up = "darkred", Down = "darkgreen")),
      simple_anno_size = unit(0.25, "cm"),
      simple_anno_size_adjust = TRUE,
      which = "row"
    )
    
   # monocle3::(cds[top_gene,], show_rownames = T)
    
    
    heat_m = Heatmap(
      top_data,
      cluster_rows = T,
      name = "Z-score",
      #col = col_fun,
      cluster_columns = cluster.col,
      show_row_dend = FALSE,
      rect_gp = gpar(col = "white"),
      show_column_names = F,
      show_row_names = T,
      top_annotation = ha_col,
      left_annotation = ha_row,
      row_names_side = "left",
      width = ncol(top_data) * unit(grid_width, "mm"),
      row_names_gp = gpar(fontsize = 10),
      column_labels = pheno[[test.column]],
      #column_dend_height = unit(10, "mm"),
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = 15, fontface = "bold"),
        labels_gp = gpar(fontsize = 15),
        legend_height = unit(5, "cm"),
        legend_width = unit(1, "cm"),
        grid_width = unit(5, "mm"),
        title_gap = unit(3, "mm")
      ),
      height = nrow(top_data) * unit(grid_width, "mm")
    ) 
  } else{
    heat_m = Heatmap(
      top_data,
      cluster_rows = F,
      name = "Z-score",
      #col = col_fun,
      cluster_columns = cluster.col,
      show_row_dend = FALSE,
      rect_gp = gpar(col = "white"),
      show_column_names = show_column_names,
      show_row_names = T,
      top_annotation = ha_col,
      #left_annotation = ha_row,
      row_names_side = "left",
      width = ncol(top_data) * unit(grid_width, "mm"),
      column_labels = pheno[[test.column]],
      row_names_gp = gpar(fontsize = 10),
      #column_dend_height = unit(10, "mm"),
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = 15, fontface = "bold"),
        labels_gp = gpar(fontsize = 15),
        legend_height = unit(5, "cm"),
        legend_width = unit(1, "cm"),
        grid_width = unit(5, "mm"),
        title_gap = unit(3, "mm")
      ),
      height = nrow(top_data) * unit(grid_width, "mm")
    )
  }
  
  heat_m = ggplotify::as.ggplot(heat_m)
  
  plot(heat_m)
  
  if(save) {
    ggsave(
      str_glue("{figure_Path}.svg"),
      width = fig.width,
      height = fig.height,
      bg = 'transparent'
    )
    
    ggsave(
      str_glue("{figure_Path}.png"),
      width = fig.width,
      height = fig.height,
      bg = 'transparent'
    )
    
    annot = t(data.frame(Test_group=pheno[[column.label2]], Disease=pheno$disease_short))
    colnames(annot) = colnames(top_data)
    write.csv(rbind(annot, top_data), str_glue("{figure_Path}.csv"))
  }
  
  #return(heat_m)
  
}


#' Heatmap
#' 
#' @param gsea_net_heat df containing gsea

plot_heatmap_wout_bar = function(test_adata,
                               cel.selected,
                               selected_gene = NULL,
                               n_top = 15,
                               aggregate = TRUE,
                               gsea_net_heat = NULL,
                               show_column_names = FALSE,
                               aggr_fun = "mean",
                               col_pal = vega_20,
                               broadcell = "\nmacrophages",
                               ref.vs.other = TRUE,
                               test.column = "Cell_type_scANVI",
                               column.title.bg = "grey",
                               column.label = "disease_short",
                               path_title = "Pathways",
                               col_test = c("grey", "grey50"),
                               scale = TRUE,
                               remove_inf = TRUE,
                               cluster.col = FALSE,
                               cluster.row = FALSE,
                               labels_annot_bar = TRUE,
                               label.only.disease = FALSE,
                               figure_Path,
                               save = FALSE,
                               fig.width = 10,
                               fig.height = 8,
                               font = 15,
                               path.col.font = 15,
                               title.font = 15,
                               title.face = "bold",
                               fontface = "bold",
                               column.title.col = "black",
                               column.name.col = "grey40",
                               title_g = 5,
                               grid_width1 = 5,
                               grid_height1 = 5,
                               legend_h = 10,
                               legend_w = 4,
                               grid_width = 5,
                               grid_width_row_annot,
                               rowannot = TRUE,
                               ...) {
 
  
  # row annotation
  colors = structure(c("grey90", "grey35"), names = c(0,1))
  
  # col_par = list()
  # 
  # legend_par = list()
  # 
  # for (column in colnames(gsea_net_heat)) {
  #   col_par[[column]] = colors
  #   legend_par[[column]] = list(
  #     title_gp = gpar(fontsize = 28, fontface = "bold"),
  #     labels_gp = gpar(fontsize = 28),
  #     grid_width = unit(1.5, "cm"),
  #     grid_height = unit(1.5, "cm"),
  #     at = c(0, 1),
  #     labels = c("NO", "YES"),
  #     title = "Associated Term"
  #     
  #   )
  # }
  # 
  # if(rowannot) {
  #   annotation_row = HeatmapAnnotation(
  #     df = gsea_net_heat,
  #     show_annotation_name = T,
  #     show_legend = rep(c(TRUE, FALSE), c(1, ncol(gsea_net_heat) - 1)),
  #     col = col_par,
  #     which = "row",
  #     simple_anno_size = unit(0.6, "cm"),
  #     simple_anno_size_adjust = TRUE,
  #     annotation_legend_param = legend_par,
  #     gp = gpar(col = "white"),
  #     annotation_width = unit(grid_width_row_annot, "cm"),
  #     annotation_name_gp = gpar(fontsize = font)
  #   )
  # } else{
  #   annotation_row = NULL
  # }
  
  ### Heatmap 1
  
  # annotation_hm1 = HeatmapAnnotation(
  #   df = gsea_net_heat,
  #   show_annotation_name = T,
  #   show_legend = rep(c(TRUE, FALSE), c(1, ncol(gsea_net_heat) - 1)),
  #   col = col_par,
  #   which = "row",
  #   simple_anno_size = unit(0.6, "cm"),
  #   simple_anno_size_adjust = TRUE,
  #   annotation_legend_param = legend_par,
  #   gp = gpar(col = "white"),
  #   annotation_width = unit(grid_width_row_annot, "cm"),
  #   annotation_name_gp = gpar(fontsize = font)
  # )
  
  gsea_net_heat = gsea_net_heat[selected_gene,]
    
  heat_m1 = Heatmap(
    gsea_net_heat,
    cluster_rows = F,
    name = path_title,
    col = colors,
    cluster_columns = F,
    show_row_dend = F,
    show_column_dend = F,
    rect_gp = gpar(col = "white"),
    show_column_names = TRUE,
    show_row_names = T,
    column_title = path_title,
    column_title_gp = gpar(
      fill = column.title.bg,
      col = column.title.col,
      border = "transparent",
      fontsize = title.font,
      fontface = title.face
    ),
    ##top_annotation = ha_col,
    row_names_side = "left",
    width = ncol(gsea_net_heat) * unit(grid_width, "mm"),
    row_names_gp = gpar(fontsize = font),
    column_names_gp = gpar(fontsize = path.col.font, col = column.name.col),
    #column_dend_height = unit(10, "mm"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = font, fontface = fontface),
      labels_gp = gpar(fontsize = font),
      #legend_height = unit(legend_h1, "cm"),
      #legend_width = unit(legend_w1, "cm"),
      grid_width = unit(grid_width1, "cm"),
      grid_height = unit(grid_height1, "cm"),
      title_gap = unit(title_g, "cm"),
      at = c(0, 1),
      labels = c("NO", "YES")
    ), 
    height = nrow(gsea_net_heat) * unit(grid_width, "mm")
  )
  
  heat_m1
  
  #Heatmap m2
  
  #We subset test_data with the top DE genes
  #or a list of genes of interest (selected_gene)
  test_adata$disease_short = droplevels(test_adata$disease_short)
  
  test_adata = test_adata[selected_gene,]
  
  
  test_adata$sample = paste0(test_adata$disease_short, "_", test_adata[[test.column]])
  test_adata$sample = factor(test_adata$sample)
  #levels(test_adata$sample)
  
  
  
  if (aggregate) {
    test_adata = scuttle::aggregateAcrossCells(test_adata,
                                               test_adata$sample,
                                               use.assay.type = "logcounts",
                                               statistics = aggr_fun)
  }
  
  
  # pheno data from sce
  pheno = as.data.frame(colData(test_adata)) %>% arrange(eval(parse(text = test.column)))
  
  if(ref.vs.other) {
    pheno[[test.column]] = str_replace(pheno[[test.column]], "Other", paste0("All other ", broadcell)) %>% as.factor()
    pheno[[test.column]] = relevel(pheno[[test.column]], ref = paste0("All other ", broadcell))
  }
  
 
  top_data = assay(test_adata, "logcounts")
  
  top_data = top_data[, rownames(pheno)]
  
  #
  if (!all.equal(colnames(top_data), rownames(pheno))) {
    stop("Colnames of expression data and rownames of pheno don't match.")
  }
  
  
  if (scale) {
    top_data = t(scale(t(top_data)))
  }
  
  #column split
  split = factor(pheno[[test.column]], levels = levels(pheno[[test.column]]))
  
  top_data = top_data[rownames(gsea_net_heat),]
  
  if(!all.equal(rownames(gsea_net_heat), rownames(top_data))) {
    stop("rows of heatmap don't match")
  }
 
  heat_m2  = Heatmap(
    top_data,
    cluster_rows = FALSE,
    name = "Z-score",
    #col = col_fun,
    cluster_columns = cluster.col,
    show_row_dend = F,
    show_column_dend = F,
    rect_gp = gpar(col = "white"),
    show_column_names = show_column_names,
    show_row_names = T,
    column_split = split,
    column_title_gp = gpar(
      fill = column.title.bg,
      col = column.title.col,
      border = "transparent",
      fontsize = title.font,
      fontface = title.face
    ),
    ##top_annotation = ha_col,
    #left_annotation = annotation_row,
    row_names_side = "left",
    width = ncol(top_data) * unit(grid_width, "mm"),
    column_labels = pheno[[column.label]],
    row_names_gp = gpar(fontsize = font),
    column_names_gp = gpar(fontsize = font),
    #column_dend_height = unit(10, "mm"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = font, fontface = fontface),
      labels_gp = gpar(fontsize = font),
      legend_height = unit(legend_h, "cm"),
      legend_width = unit(legend_w, "cm"),
      grid_width = unit(legend_w, "cm"),
      title_gap = unit(title_g, "cm")
    ),
    height = nrow(top_data) * unit(grid_width, "mm")
  )
  
  
  heat_m2
  # 
  # 
  heat_m = heat_m1 + heat_m2
  # 
  # heat_m = ggplotify::as.ggplot(heat_m)
  
  draw(heat_m)
  
  if (save) {
    # ggsave(
    #   str_glue("{figure_Path}.svg"),
    #   width = fig.width,
    #   height = fig.height,
    #   bg = 'transparent', dpi = 320
    # )
    # 
    # ggsave(
    #   str_glue("{figure_Path}.png"),
    #   width = fig.width,
    #   height = fig.height,
    #   bg = 'transparent', dpi = 320
    # )
    
    
    svg(filename = str_glue("{figure_Path}.svg"),
      width = fig.width,
      height = fig.height,
      bg = 'transparent'
    )
    print(heat_m)
    dev.off()
    
    
    svg(filename = str_glue("{figure_Path}.png"),
        width = fig.width,
        height = fig.height,
        bg = 'transparent'
    )
    print(heat_m)
    dev.off()
    
    # 
    # png(filename = str_glue("{figure_Path}.png"),
    #   width = fig.width,
    #   height = fig.height,
    #   bg = 'transparent', res = 320
    # )
    # 
    # print(heat_m)
    # dev.off()
    
    
    annot = t(data.frame(
      Test_group = pheno[[column.label]],
      Disease = pheno$disease_short
    ))
    colnames(annot) = colnames(top_data)
    
    df = rbind(annot, top_data)
    
    write.csv(df, str_glue("{figure_Path}.csv"))
  }
  
  
 return(heat_m)
}


plot_heatmap_2 = function(test_adata=test_adata_prep,
                        deg=int_macr_gene_program_Naïve_treat_RA,
                        cel.selected="Intermediate_macrophage",
                        selected_gene = NA,
                        n_top = 15,
                        aggregate = TRUE,
                        aggr_fun = "mean",
                        col_pal = vega_20,
                        test.column="nhoo_test_condition",
                        col_test=c("grey", "grey50"),
                        scale = TRUE,
                        remove_inf = TRUE,
                        cluster.col = FALSE,
                        grid_width=5) {
  
  # Remove Inf from deg dataset
  if(remove_inf) {
    deg = deg[!(deg$lfc %in% c(Inf, -Inf)),]
  }
  # extract top genes
  up = deg %>% filter(lfc > 0) %>% arrange(desc(lfc))
  down = deg %>% filter(lfc < 0) %>% arrange(lfc)
  top_ = na.omit(rbind(up[1:n_top, ], down[1:n_top, ]))
  
  #We subset test_data with the top DE genes
  #or a list of genes of interest (selected_gene)
  if(!is.na(selected_gene)) {
    data = test_adata[selected_gene,]
  } else{
    data = test_adata[top_$name,]
  }
  
  data$sample = paste0(data$disease_bin, "_", data$nhoo_test_condition)
  data$sample = factor(data$sample)
  levels(data$sample)
  
  
  
  if(aggregate){
    data = scuttle::aggregateAcrossCells(data,
                                         data$sample,
                                         use.assay.type = "logcounts",
                                         statistics = aggr_fun)
  }
  
  
  data$nhoo_test_condition = ifelse(data$nhoo_test_condition == "test_nhood",
                                    cel.selected,
                                    "Rest")
  
  # pheno data from sce
  pheno = as.data.frame(colData(data)) %>% arrange(eval(parse(text=test.column)))
  
  top_data = assay(data, "logcounts")
  
  top_data = top_data[,rownames(pheno)]
  
  #
  if (!all.equal(colnames(top_data), rownames(pheno))) {
    stop("Colnames of expression data and rownames of pheno don't match.")
  }
  
  
  data$disease_short = droplevels(data$disease_short)
  
  data$nhoo_test_condition = factor(data$nhoo_test_condition)
  
  # column annotation
  #disease
  dis_color = c()
  Disease = factor(test_adata$disease_short)
  
  for (i in 1:length(levels(Disease))) {
    dis = Disease[i]
    dis_color = c(dis_color, col_pal[i])
  }
  
  names(dis_color) = levels(Disease)
  # second annotation bar (por example test group)
  
  da_color = c()
  
  for (i in 1:length(unique(data[[test.column]]))) {
    da = unique(data[[test.column]])[i]
    da_color = c(da_color, col_test[i])
  }
  
  names(da_color) = unique(data[[test.column]])
  
  # col annotation
  ha_col = HeatmapAnnotation(
    df = data.frame(
      Test_group = factor(pheno[[test.column]]),
      Disease = factor(pheno$disease_short)
    ),
    show_legend = c(TRUE, TRUE),
    show_annotation_name = F,
    col = list(
      Test_group = da_color,
      Disease = dis_color
    ),
    simple_anno_size = unit(0.25, "cm"),
    simple_anno_size_adjust = TRUE,
    annotation_legend_param = list(
      Test_group = list(
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 12),
        grid_width = unit(5, "mm"),
        grid_height = unit(5, "mm")
      ),
      Disease = list(
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 12),
        grid_width = unit(5, "mm"),
        grid_height = unit(5, "mm")
      )
    )
  )
  
  
  if (scale) {
    top_data = t(scale(t(top_data)))
  }
  
  
  # Included FC direction annotation if selected_gene is not provided
  if (is.na(selected_gene)) {
    top_$FC_dir = ifelse(top_$lfc > 0, "Up", "Down")
    
    top_$FC_dir = factor(top_$FC_dir)
    # Row annotation
    
    ha_row = HeatmapAnnotation(
      df = data.frame(FC_dir = top_$FC_dir),
      show_annotation_name = F,
      col = list(FC_dir = c(Up = "darkred", Down = "darkgreen")),
      simple_anno_size = unit(0.25, "cm"),
      simple_anno_size_adjust = TRUE,
      which = "row"
    )
    
    
    
    
    Heatmap(
      top_data,
      cluster_rows = F,
      name = "Z-score",
      #col = col_fun,
      cluster_columns = cluster.col,
      show_row_dend = FALSE,
      rect_gp = gpar(col = "white"),
      show_column_names = F,
      show_row_names = T,
      top_annotation = ha_col,
      left_annotation = ha_row,
      row_names_side = "left",
      width = ncol(top_data) * unit(grid_width, "mm"),
      row_names_gp = gpar(fontsize = 10),
      #column_dend_height = unit(10, "mm"),
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = 15, fontface = "bold"),
        labels_gp = gpar(fontsize = 15),
        legend_height = unit(5, "cm"),
        legend_width = unit(1, "cm"),
        grid_width = unit(5, "mm"),
        title_gap = unit(3, "mm")
      ),
      height = nrow(top_data) * unit(grid_width, "mm")
    )
  } else{
    Heatmap(
      top_data,
      cluster_rows = F,
      name = "Z-score",
      #col = col_fun,
      cluster_columns = cluster.col,
      show_row_dend = FALSE,
      rect_gp = gpar(col = "white"),
      show_column_names = F,
      show_row_names = T,
      top_annotation = ha_col,
      #left_annotation = ha_row,
      row_names_side = "left",
      width = ncol(top_data) * unit(grid_width, "mm"),
      row_names_gp = gpar(fontsize = 10),
      #column_dend_height = unit(10, "mm"),
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = 15, fontface = "bold"),
        labels_gp = gpar(fontsize = 15),
        legend_height = unit(5, "cm"),
        legend_width = unit(1, "cm"),
        grid_width = unit(5, "mm"),
        title_gap = unit(3, "mm")
      ),
      height = nrow(top_data) * unit(grid_width, "mm")
    )
  }
  
}


#' Function to prepare gsea data for network visualization
#' 
#' @import RColorBrewer
#

gsea_2_net = function(gsea_res,
                      pval.threshold = 0.05,
                      pval = "Adjusted P-value",
                      Combined.Score = "Combined.Score",
                      Number_genes = 2,
                      gsea_Path = NULL,
                      openTarget = NULL) {
  
  
  
  for (i in 1:nrow(gsea_res)) {
    gsea_res$Overlap_2[i] = eval(parse(text = gsea_res$Overlap[i])) %>% as.numeric()
    gsea_res$Number_genes[i] = str_remove(gsea_res$Overlap[i], "/.*") %>% as.numeric()
  }
  
  # select gene sets with at least 2 genes
  gsea_res = gsea_res %>% filter(.[[pval]] < pval.threshold & Number_genes>=Number_genes)
  
  
  tab = na.omit(
    data.frame(
      Term = NA,
      symbol = NA,
      color = NA,
      combined_score = NA,
      overlap = NA,
      number_of_genes = NA,
      adj.pval = NA
    )
  )
  
  # Create color palette
  colourCount = length(unique(gsea_res$Term))
  
  getPalette = colorRampPalette(vega_20)
  
  Palette = getPalette(colourCount)

  
  #Some gsea failed the p-value filtering and have 0 nrow
  #If statement is used avoid errors
  if(nrow(gsea_res)>0){
    
    for (i in 1:nrow(gsea_res)) {
      ge = str_split(gsea_res$Genes[i], pattern = ";")
      
      
      tab0 = data.frame(
        Term = gsea_res$Term[i],
        symbol = ge[[1]],
        color = Palette[i],
        combined_score = gsea_res[[Combined.Score]][i],
        overlap = gsea_res$Overlap_2[i],
        number_of_genes = gsea_res$Number_genes[i],
        adj.pval = gsea_res$`Adjusted P-value`[i]
      )
      tab = rbind(tab, tab0)
    }
    
    if(!is.null(openTarget)){
      tab = merge(tab, openTarget, by = "symbol", all.x=TRUE)
      
      tab$targeted = NA
      
      for (i in 1:nrow(tab)) {
        if (!is.na(tab$drugName[i])) {
          tab$targeted[i] = "Yes"
        } else{
          tab$targeted[i] = "No"
        }
      }
    }
    
    
    print(head(tab))
  }
  
  if(!is.null(gsea_Path)) {
    write.csv(tab, str_glue("{gsea_Path}.csv"), row.names = FALSE)
    #write_delim(tab, str_glue("{gsea_Path}.tsv"))
  }
  
  return(tab)
  
}


#' Fuction toplot normalized cell count per donor (Cell_per_1000_per_donor)
#' 
#' @param data data with at least the following columns: "disease_tissue", "Cell_type_scANVI", "total_per_donor_per_celltype", "Cell_per_1000"
#' @param xlab label of x axis
#' @param ylab label of y axis
#' @param color_dis_tis legend map
#' @param disease_column column containing X variable
#' @param facet_var String. Column to be used for facet_wrap if facet is TRUE
#' @param group_sorting_col col name to be used to sort groups in facet_grid

plot_Cell_per_1000_per_donor = function(data, color_dis_tis, xlab = "Tissue origin", ylab = "", plot.title = "", legend_title="Groups", disease_column="disease_tissue", 
                                        Cell_type = "Cell_type_scANVI", show.legend=TRUE, axis.size=18, legend.size = 18,
                                        title.size = 18, legend.keys.size=0.6, facet = FALSE, facet_var = "Tissue", facet_nrow = 1,
                                        group_sorting_col = "Groups", xfactor.level=NA, same.x.order=FALSE, point.size=1) {
  
  if(!same.x.order){
    if (!facet) {
      if (show.legend) {
        ggplot(data = data,
               aes(
                 y = Cell_per_1000,
                 x = reorder_within(
                   x = .data[[disease_column]],
                   by = Cell_per_1000,
                   within = .data[[Cell_type]],
                   fun = mean
                 ),
                 fill = get(disease_column)
               )) +
          geom_boxplot(show.legend = show.legend) +
          geom_beeswarm(stat = "identity",
                        shape = 21,
                        show.legend = F,
                        size = point.size) +
          labs(
            x = xlab,
            y = ylab,
            title = plot.title,
            fill = legend_title
          ) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
          theme_classic() + #coord_flip() +
          scale_color_manual(values = color_dis_tis) +
          scale_fill_manual(values = color_dis_tis) +
          scale_x_reordered() +
          theme(
            axis.text.x = element_blank(),
            
            axis.text.y = element_text(size = axis.size),
            axis.title = element_text(size = axis.size + 1),
            legend.key.size = unit(legend.keys.size, 'cm'),
            legend.title = element_text(size = legend.size + 1),
            legend.text = element_text(size = legend.size),
            plot.title = element_text(size = title.size, face = "plain"),
            
            axis.ticks.x = element_blank()
          )
      } else {
        ggplot(data = data,
               aes(
                 y = Cell_per_1000,
                 x = reorder_within(
                   x = .data[[disease_column]],
                   by = Cell_per_1000,
                   within = .data[[Cell_type]],
                   fun = mean
                 ),
                 fill = get(disease_column)
               )) +
          geom_boxplot(show.legend = show.legend) +
          geom_beeswarm(
            stat = "identity",
            shape = 21,
            aes(alpha = 0.5),
            show.legend = F,
            size = point.size
          ) +
          labs(
            x = xlab,
            y = ylab,
            title = plot.title,
            fill = legend_title
          ) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
          theme_bw() + #coord_flip() +
          scale_color_manual(values = color_dis_tis) +
          scale_fill_manual(values = color_dis_tis) +
          scale_x_reordered() +
          theme(
            axis.text.x = element_text(size = axis.size, angle = 90, vjust = 0.5, hjust = 1),
            axis.title = element_text(size = axis.size + 1),
            axis.text.y = element_text(size = axis.size),
            legend.key.size = unit(legend.keys.size, 'cm'),
            legend.text = element_text(size = legend.size),
            panel.grid = element_blank(),
            plot.title = element_text(size = title.size, face = "plain"),
            legend.title = element_text(size = legend.size + 1)
            #axis.ticks.x = element_blank()
          )
      }
    } else {
      
      ##rename axis.x.tick text
      
      x_tick = rename_x_label(data=data, col1=group_sorting_col, col2=disease_column)
      
      ## with facet
      
      ggplot(data = data,
             aes(
               y = Cell_per_1000,
               x = reorder(
                 .data[[disease_column]],
                 Cell_per_1000
               ),
               fill = get(group_sorting_col)
             )) +
        geom_boxplot(show.legend = show.legend) +
        geom_beeswarm(stat = "identity",
                      shape = 21,
                      show.legend = F,
                      size = point.size) +
        labs(
          x = xlab,
          y = ylab,
          title = plot.title,
          fill = legend_title
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        theme_bw() + #coord_flip() +
        scale_color_manual(values = color_dis_tis) +
        scale_fill_manual(values = color_dis_tis) +
        scale_x_discrete(labels=x_tick) +
        ggforce::facet_row(vars(.data[[facet_var]]), scales = 'free', space = 'free') +
        theme(
          axis.text.x = element_text(size = axis.size, angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = axis.size),
          axis.title = element_text(size = axis.size + 1),
          legend.position = "none",
          plot.title = element_text(size = title.size, face = "plain"),
          strip.text.x = element_text(size = title.size)
        )
      
      
    }
    
    ######### Maintain the same order of x-axis
  }else {
    if (!facet) {
      if (show.legend) {
        ggplot(data = data,
               aes(
                 y = Cell_per_1000,
                 x = factor(
                   x = .data[[disease_column]],
                   xfactor.level
                 ),
                 fill = get(disease_column)
               )) +
          geom_boxplot(show.legend = show.legend) +
          geom_beeswarm(stat = "identity",
                        shape = 21,
                        show.legend = F,
                        size = point.size) +
          labs(
            x = xlab,
            y = ylab,
            title = plot.title,
            fill = legend_title
          ) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
          theme_classic() + #coord_flip() +
          scale_color_manual(values = color_dis_tis) +
          scale_fill_manual(values = color_dis_tis) +
          theme(
            axis.text.x = element_blank(),
            
            axis.text.y = element_text(size = axis.size),
            axis.title = element_text(size = axis.size + 1),
            legend.key.size = unit(legend.keys.size, 'cm'),
            legend.title = element_text(size = legend.size + 1),
            legend.text = element_text(size = legend.size),
            plot.title = element_text(size = title.size, face = "plain"),
            axis.ticks.x = element_blank()
          )
      } else {
        ggplot(data = data,
               aes(
                 y = Cell_per_1000,
                 x = factor(
                   x = .data[[disease_column]],
                   xfactor.level
                 ),
                 fill = get(disease_column)
               )) +
          geom_boxplot(show.legend = show.legend) +
          geom_beeswarm(
            stat = "identity",
            shape = 21,
            aes(alpha = 0.5),
            show.legend = F,
            size = point.size
          ) +
          labs(
            x = xlab,
            y = ylab,
            title = plot.title,
            fill = legend_title
          ) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
          theme_bw() + #coord_flip() +
          scale_color_manual(values = color_dis_tis) +
          scale_fill_manual(values = color_dis_tis) +
          theme(
            axis.text.x = element_text(size = axis.size, angle = 90, vjust = 0.5, hjust = 1),
            axis.title = element_text(size = axis.size + 1),
            axis.text.y = element_text(size = axis.size),
            legend.key.size = unit(legend.keys.size, 'cm'),
            legend.text = element_text(size = legend.size),
            panel.grid = element_blank(),
            plot.title = element_text(size = title.size, face = "plain"),
            legend.title = element_text(size = legend.size + 1)
            #axis.ticks.x = element_blank()
          )
      }
    } else {
      
      ##rename axis.x.tick text
      
      x_tick = rename_x_label(data=data, col1=group_sorting_col, col2=disease_column)
      
      ## with facet
      
      ggplot(data = data,
             aes(
               y = Cell_per_1000,
               x = factor(
                 x = .data[[disease_column]],
                 xfactor.level
               ),
               fill = get(group_sorting_col)
             )) +
        geom_boxplot(show.legend = show.legend) +
        geom_beeswarm(stat = "identity",
                      shape = 21,
                      show.legend = F,
                      size = point.size) +
        labs(
          x = xlab,
          y = ylab,
          title = plot.title,
          fill = legend_title
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        theme_bw() + #coord_flip() +
        scale_color_manual(values = color_dis_tis) +
        scale_fill_manual(values = color_dis_tis) +
        scale_x_discrete(labels=x_tick) +
        ggforce::facet_row(vars(.data[[facet_var]]), scales = 'free', space = 'free') +
        theme(
          axis.text.x = element_text(size = axis.size, angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = axis.size),
          axis.title = element_text(size = axis.size + 1),
          legend.position = "none",
          plot.title = element_text(size = title.size, face = "plain"),
          strip.text.x = element_text(size = title.size)
        )
      
      
    }
  }
  
  
}


color_disease_tis = c(
  IBD_colon = "#1f77b4",
  IBD_ileum = "#8c564b",
  Naïve_treat_RA_syno_membr = "#ff7f0e",
  Remission_RA_syno_membr = "#2ca02c",
  Oligo_JIA_syno_fld = "#d62728",
  Resistant_treat_RA_syno_membr = "#9467bd",
  Uninfl_IBD_colon = "#bcbd22",
  #"#e377c2",
  Healthy_colon = "#7f7f7f",
  Healthy_ileum = "#f3f3f3",
  Healthy_syno_membr = "#7b7c9c",
  IBD = "#1f77b4",
  Naïve_treat_RA = "#ff7f0e",
  Remission_RA = "#2ca02c",
  Oligo_JIA = "#d62728",
  Resistant_treat_RA = "#9467bd",
  Uninfl_IBD = "#bcbd22",
  Healthy_dermis = "#7f7f7f",
  Healthy_epidermis = "#7f7f7f",
  Healthy_lung = "#7f7f7f",
  Eczema_dermis  = "#c1019d", 
  Eczema_epidermis= "#1ee7bd", 
  Psoriasis_dermis= "#e99999",
  Psoriasis_epidermis= "#1ee7ed",
  Asthma_lung= "#98df8a",
  hashimoto_thyroid= "#9edae5",
  Paed_IBD = "#9177f2",
  Eczema  = "#a1011d",
  Psoriasis = "#e999b9",
  Asthma = "#98df8a",
  Hashimoto_thyroid= "#9edae5",
  hashimoto = "#9edae5",
  Paediatric_IBD = "#9177f2"
)


#"#1f77b4" "#aec7e8" "#ff7f0e" "#ffbb78" "#2ca02c" "#98df8a" "#d62728" "#ff9896" "#9467bd" "#c5b0d5" "#8c564b" "#c49c94" "#e377c2" "#f7b6d2" "#7f7f7f" "#c7c7c7" "#bcbd22"
#"#dbdb8d" "#17becf" "#9edae5"

#Healthy all_CD all_SLE all_hashimoto all_sjo_syndr all_spondoarthritis

color_disease_pbmc = c(
  IBD = "#1f77b4",
  SLE = "#ff7f0e",
  hashimoto = "#2ca02c",
  `Hashimoto's thyroiditis` = "#2ca02c",
  sjo_syndr = "#d62728",
  spondyloarthritis = "#9467bd",
  Healthy = "#7f7f7f",
  Oligo_JIA = "#dbdb8d",
  "Hashimoto thyroiditis" = "#2ca02c",
  IBD_and_spondoarthritis = "#f1c1c2",
  "Oligoarticular JIA" = "#dbdb8d",
  "Psoriatic arthritis" = "#e377c2",
  "Sjögren's syndrome" = "#d62728",
  Spondyloarthritis = "#9467bd",
  IBD_and_Spondoarthritis = "#0e77f2"
)

#"cMon_1"    "CCR7_CDE3" "cMon_2"    "DC2"       "Non-cMon"  "Mon_SIRPA" "DC1"       "pDC"


color_cell_pbmc0 = c(
  cMon_1 = "#1f77b4",
  CCR7_CDE3 = "#ff7f0e",
  DC2 = "#2ca02c",
  cMon_2 = "#d62728",
  `Activated cMon` = "#1f77b4",
  cMon = "#d62728",
  `Non-cMon` = "#9467bd",
  Mon_SIRPA = "#1bbbbb",
  DC1 = "#f1c1c2",
  pDC = "#dbdb8d",
  "Naive B cells" = "#8c564b",
  ABC = "#bcbd22", 
  "Memory B cells_1" = "#e377c2",
  "B cells_CR1" = "#ad494a",
  "Memory B cells_2" = "#ff9896",
  "Plasma cells" = "#c49c94",
  "Cycling B cells" = "#9edae5",
  CD16__NK_CD56_plus = "#1f77b4",
  CD16_plus_IFNG_plus_NK = "#ff7f0e",
  CD16_plus_NK = "#2ca02c",
  Cycling_NK = "#d62728",
  Cycling_T_cell = "#9467bd",
  DN = "#7e1e7f",
  GZMK_producing_CD4 = "#f1c1c2",
  MAIT_act = "#dbdb8d",
  MAIT_cells = "#8c564b",
  NKT = "#bcbd22",
  TCD4_CD8 = "#e377c2",
  "Tcm/Naive_cytotoxic" = "#ad494a",
  "Tcm/Naive_cytotoxic_1" = "#ff9896",
  "Tcm/Naive_cytotoxic_2" = "#c49c94",
  "Tcm/Naive_h_activated" = "#9edae5",
  "Tcm/Naive_helper" = "#ffbb78",
  Tcm_CCR6 = "#98df8a",
  "Tem/Temra_cytotoxic" = "#aec7e8",
  "Tem/Temra_cytotoxic_T_cells" = "#8c6d31",
  "Tem/Trm_cytotoxic_T_cells" = "#999999",
  Tem_cytotoxic = "grey90",
  Tregs = "#1e77f2",
  γδ_T = "#7f77f2",
  γδ_T_ZNF683_plus= "#1bbbbb"
)

color_cell_pbmc1 = c(
  cMon_1 = "#1f77b4",
  CCR7_CDE3 = "#ff7f0e",
  DC2 = "#2ca02c",
  cMon_2 = "#d62728",
  `Activated cMon` = "#1f77b4",
  cMon = "#d62728",
  `Non-cMon` = "#9467bd",
  Mon_SIRPA = "#1bbbbb",
  DC1 = "#f1c1c2",
  pDC = "#dbdb8d",
  `Naive B cells` = "#8c564b",
  ABC = "#bcbd22",
  `Memory B cells_1` = "#e377c2",
  `B cells_CR1` = "#ad494a",
  `Memory B cells_2` = "#ff9896",
  `Plasma cells` = "#c49c94",
  `Cycling B cells` = "#9edae5",
  
  CD16__NK_CD56_plus = "#1f77b4",
  `CD16- NK CD56+` = "#1f77b4",
  
  CD16_plus_IFNG_plus_NK = "#ff7f0e",
  `CD16+ IFNG+ NK` = "#ff7f0e",
  
  CD16_plus_NK = "#2ca02c",
  `CD16+ NK` = "#2ca02c",
  
  Cycling_NK = "#d62728",
  `Cycling T cell` = "#9467bd",
  
  DN = "#7e1e7f",
  
  GZMK_producing_CD4 = "#f1c1c2",
  `GZMK producing CD4` = "#f1c1c2",
  
  MAIT_act = "#dbdb8d",
  MAIT_cells = "#8c564b",
  `MAIT cells` = "#8c564b",
  NKT = "#bcbd22",
  TCD4_CD8 = "#e377c2",
  
  Tcm_CCR6 = "#98df8a",
  `Tcm CCR6` = "#98df8a",
  
  `Tcm/Naive_cytotoxic` = "#ad494a",
  `Tcm/Naive cytotoxic` = "#ad494a",
  `Tcm/Naive_h_activated` = "#9edae5",
  `Tcm/Naive h activated` = "#9edae5",
  
  `Tcm/Naive_helper` = "#ffbb78",
  `Tcm/Naive helper` = "#ffbb78",
  
  Tem_cytotoxic = "grey90",
  `Tem cytotoxic` = "grey90",
  
  `Tem/Temra_cytotoxic` = "#aec7e8",
  `Tem/Temra cytotoxic` = "#aec7e8",
  
  `Tem/Trm_cytotoxic_T_cells` = "#8c6d31",
  `Tem/Trm cytotoxic T cells` = "#8c6d31",
  
  Tregs = "#1e77f2",
  
  `γδ_T` = "#7f77f2",
  `γδ T` = "#7f77f2",
  
  γδ_T_ZNF683_plus = "#1bbbbb",
  `γδ T ZNF683+` = "#1bbbbb"
)

color_cell_pbmc = c(
  `CD16- NK CD56+` = "#1f77b4",
  `CD16+ IFNG+ NK` = "#ff7f0e",
  `CD16+ NK` = "#2ca02c",
  `Cycling T cell` = "#9467bd",
  DN = "#7e1e7f",
  `GZMK producing CD4` = "#f1c1c2",
  `MAIT cells` = "#8c564b",
  MAIT_act = "#dbdb8d",
  NKT = "#bcbd22",
  TCD4_CD8 = "#e377c2",
  `Tcm CCR6` = "#98df8a",
  `Tcm/Naive cytotoxic` = "#ad494a",
  `Tcm/Naive h activated` = "#9edae5",
  `Tcm/Naive helper` = "#ffbb78",
  `Tem cytotoxic` = "grey90",
  `Tem/Temra cytotoxic` = "#aec7e8",
  `Tem/Trm cytotoxic T cells` = "#8c6d31",
  Tregs = "#1e77f2",
  `γδ T` = "#7f77f2",
  `γδ T ZNF683+` = "#1bbbbb"
)

#"#1f77b4" "#ff7f0e" "#2ca02c" "#d62728" "#9467bd" "#8c564b" "#e377c2" "#bcbd22" "#17becf" "#aec7e8" "#ffbb78" "#98df8a" "#ff9896" "#c5b0d5" "#c49c94" "#f7b6d2" "#dbdb8d" "#9edae5" "#ad494a" "#8c6d31"

color_cell_tissue = c(
  "Intermediate macrophage" = "#A6CEE3",
  "Non_cMon" = "#d62728",
  "CD14+CD16+ Monocytes" = "#d62728",
  "CD14+CD16+ Mon" = "#d62728",
  "IL1B+ macrophage" = "#A6CEE3",
  "IL1B high macrophage" = "#A6CEE3",
  "Mast cell" = "#7f7f7f",
  "Mac CCL2+LYVE1+" = "#98D277",
  "Mac LYVE1+" = "#6F9E4C",
  "Mac MERTK_hi_TREM2_low" = "#ad494a",
  "Mac NUPR1_TREM2_hi" = "#F06C45",
  "Mac S100A8 prod" = "#8c554a",
  "Mac TREM2_hi" = "#D9A295",
  Mac_CX3CR1 =  "#7D54A5",
  "Mac CX3CR1+" =  "#7D54A5",
  Mac1 = "#e373c2",
  Mac2 = "#f7b6d2",
  Cycling = "#1f77b4",
  DC2 = "#c5b0d5", 
  DC1 = "#ffbb77",
  "Mig DC" = "#c7c7c7", 
  pDC = "#bcbc22",
  moDC = "grey90"
) 


color_cell_trajectory = c(
  cMon_1 = "#1f77b4",
  cMon_2 = "#d62728",
  Mon_SIRPA = "#1bbbbb",
  "Mac CX3CR1+" =  "#7D54A5",
  "IL1B high macrophage" = "#A6CEE3",
  "Mac S100A8 prod" = "#8c554a"
) 

color_cell_trajectory0 = c(
  cMon_1 = "grey90",
  cMon_2 = "#d62728",
  Mon_SIRPA = "grey90",
  "Mac CX3CR1+" =  "#7D54A5",
  "IL1B high macrophage" = "grey90",
  "Mac S100A8 prod" = "#8c554a"
) 


lineage3 = c(
  "T cells" = "#E9B5A8",
  "B cells" = "#DEC98F",
  "Myeloid cells" = "#B2C0A0",
  "Mast cells" = "#D9A9B2",
  "MNP and T/B Doublets" = "#D0CDC5",
  "NK cells"= "#7f7f7f"
)

lineage2 = c(
  "T cells" = "#9CCFE0",
  "B cells" = "#A98BC0",
  "Myeloid cells" = "#79AFC8",
  "Mast cells" = "#72CACA",
  "MNP and T/B Doublets" = "#C8C8C8",
  "NK cells" = "#999999"
)


lineage = c(
  "T cells" = "#D98C7A",
  "B cells" = "#C7A15A",
  "Myeloid cells" = "#8E9A6B",
  "Mast cells" = "#C17C8A",
  "MNP and T/B Doublets" = "#B7B1A3",
  "NK cells" = "#999999"
)



#"#1f77b4" "#ff7f0e" "#2ca02c" "#d62728" "#9467bd" "#8c564b" "#e377c2" "#bcbd22" "#17becf" "#aec7e8" "#ffbb78" "#98df8a" "#ff9896" "#c5b0d5" "#c49c94" "#f7b6d2" "#dbdb8d"
#"#9edae5" "#ad494a" "#8c6d31"

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}

rename_x_label = function(data, col1, col2){
  x = as.character(data[[col1]])
  names(x) = as.character(data[[col2]])
  x[!duplicated(names(x))]
}


#' Function to plot umap
#' 
#' @param umap_data data frame containing umap and cell meta data
#' @param cluster string indicating the column containing cluster factors
#' @param selected_clusters clusters to be included. Default is all for all clusters. A vector of clusters can be provided.
#' @param umap_1 colname of umap1 or tsne1 column 
#' @param umap_2 colname of umap2 or tsne2 column
#' @param arrow.length int length of umap axis arrow
#' @param repel "text" or "label"

plot_umap = function(umap_data = umap,
                     cluster = "Cell_type_scANVI",
                     umap_1 = "umap1",
                     umap_2 = "umap2",
                     UMAP_1_lab = "UMAP_1",
                     UMAP_2_lab = "UMAP_2",
                     selected_clusters = "all",
                     shape.size = 0.5,
                     repel_text_size = 2,
                     cluster_palette = vega_20,
                     arrow.length=3,
                     font="bold",
                     axis.lab.size=3,
                    legend.position = "none",
                     repel = "text", cust.alpha = 0.4, ...) {
  
  theme_umap_2 <- theme_bw(base_size = 25) + theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(size = 30),
    legend.position = legend.position,
    axis.title = element_blank(),
    panel.background = element_rect(fill = 'transparent'),
    #transparent panel bg
    plot.background = element_rect(fill = 'transparent', color = NA),
    #transparent plot bg
    legend.background = element_rect(fill = 'transparent'),
    #transparent legend bg
    legend.box.background = element_rect(fill = 'transparent')
  )
  
  # Label data
  if (all(selected_clusters != "all")) {
    umap_data = umap_data[umap_data[[cluster]] %in% selected_clusters,]
  }
  
  
  umap_rep = umap_data %>% group_by(.data[[cluster]]) %>%
    summarise(umap1_center = median(.data[[umap_1]]),
              umap2_center = median(.data[[umap_2]]))
  
  
  
  #Plot
  if(repel == "label") {
    umap_data %>%
      ggplot(aes(
        x = .data[[umap_1]],
        y = .data[[umap_2]],
        color = .data[[cluster]],
        label = .data[[cluster]]
      )) +
      geom_point(show.legend = FALSE, size = shape.size) +
      scale_colour_manual(values = cluster_palette) +
      geom_label_repel(
      data = umap_rep,
      aes(label = .data[[cluster]], x = umap1_center, y = umap2_center, fontface = font),
      point.padding = unit(0.5, "lines"), size = repel_text_size, color = "grey20", fill = alpha("white", cust.alpha)
    ) +
      
      geom_segment(
        data = umap_data,
        aes(
          x = min(.data[[umap_1]]),
          y = min(.data[[umap_2]]),
          xend = min(.data[[umap_1]]) + arrow.length,
          yend = min(.data[[umap_2]])
        ),
        colour = "black",
        size = 1,
        arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
        linejoin = "mitre"
      ) +
      
      geom_segment(
        data = umap_data,
        aes(
          x = min(.data[[umap_1]]),
          y = min(.data[[umap_2]]),
          xend = min(.data[[umap_1]]) ,
          yend = min(.data[[umap_2]]) + arrow.length
        ),
        colour = "black",
        size = 1,
        arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
        linejoin = "mitre"
      ) +
      
      annotate(
        "text",
        x = (min(umap_data[[umap_1]])+ (min(umap_data[[umap_1]]) + arrow.length))/2,
        y = min(umap_data[[umap_2]]),
        label = UMAP_1_lab,
        color = "black",
        size = axis.lab.size,
        fontface = "bold",
        angle = 0,
        hjust = 0.5,
        vjust = 1.5
      ) +
      
      annotate(
        "text",
        x = min(umap_data[[umap_1]]),
        y = (min(umap_data[[umap_2]]) + (min(umap_data[[umap_2]]) + arrow.length))/2,
        label = UMAP_2_lab,
        color = "black",
        size = axis.lab.size,
        fontface = "bold" ,
        angle = 90,
        hjust = 0.5,
        vjust = -0.5
      ) +
      #coord_fixed(ratio = 1) #+
      
      theme_umap_2
  }else if(repel == "text"){
    umap_data %>%
      ggplot(aes(
        x = .data[[umap_1]],
        y = .data[[umap_2]],
        color = .data[[cluster]],
        label = .data[[cluster]]
      )) +
      geom_point(show.legend = FALSE, size = shape.size) +
      scale_colour_manual(values = cluster_palette) +
      geom_text_repel(data = umap_rep,
                      aes(label = .data[[cluster]], x = umap1_center, 
                          y = umap2_center, fontface = font),
                      seed = 42, colour = "black",
                      segment.color = 'black', size = repel_text_size, ...) +
      
      geom_segment(
        data = umap_data,
        aes(
          x = min(.data[[umap_1]]),
          y = min(.data[[umap_2]]),
          xend = min(.data[[umap_1]]) + arrow.length,
          yend = min(.data[[umap_2]])
        ),
        colour = "black",
        size = 1,
        arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
        linejoin = "mitre"
      ) +
      
      geom_segment(
        data = umap_data,
        aes(
          x = min(.data[[umap_1]]),
          y = min(.data[[umap_2]]),
          xend = min(.data[[umap_1]]) ,
          yend = min(.data[[umap_2]]) + arrow.length
        ),
        colour = "black",
        size = 1,
        arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
        linejoin = "mitre"
      ) +
      
      annotate(
        "text",
        x = (min(umap_data[[umap_1]])+ (min(umap_data[[umap_1]]) + arrow.length))/2,
        y = min(umap_data[[umap_2]]),
        label = UMAP_1_lab,
        color = "black",
        size = axis.lab.size,
        fontface = "bold",
        angle = 0,
        hjust = 0.5,
        vjust = 1.5
      ) +
      
      annotate(
        "text",
        x = min(umap_data[[umap_1]]),
        y = (min(umap_data[[umap_2]]) + (min(umap_data[[umap_2]]) + arrow.length))/2,
        label = UMAP_2_lab,
        color = "black",
        size = axis.lab.size,
        fontface = "bold" ,
        angle = 90,
        hjust = 0.5,
        vjust = -0.5
      ) +
      #coord_fixed(ratio = 1) #+
      
      theme_umap_2
  }
  
  }
  



#' Function to plot nhood graph
#' 
#' @param nhood_data data frame containing nhood graph data
#' @param cluster string indicating the column containing nhood annotation
#' @param selected_clusters clusters to be included. Default is all for all clusters. A vector of clusters can be provided.
#' @param nhood_1 colname of umap1 or tsne1 column 
#' @param nhood_2 colname of umap2 or tsne2 column
#' @param arrow.length int length of umap axis arrow
#' @param FC_dir direction of Fc enriched_depleted, enriched ordepleted
#' 

plot_nhood_graph = function(nhood_data = Nhood_Graph_df_2,
                            cluster = "nhood_annotation",
                            nhood_1 = "nhood_graph1",
                            nhood_2 = "nhood_graph2",
                            NHOOD_1_lab = "NHOOD_1",
                            NHOOD_2_lab = "NHOOD_2",
                            selected_clusters = "all",
                            shape.size = 2,
                            repel_text_size = 2,
                            cluster_palette = "OrRd",
                            arrow.length=3,
                            font="bold",
                            FC = TRUE,
                            logFC = "logFC",
                            FC_dir = "enriched_depleted",
                            ncol, ...) {
  
  theme_nhood <- theme_bw(base_size = 25) + theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    #plot.title = element_text(size = 30),
    legend.position = "none",
    axis.title = element_blank(),
    ...
  )
  
  
  theme_nhood_2 <- theme_bw(base_size = 25) + theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    #plot.title = element_text(size = 30),
    legend.position = "none",
    axis.title = element_blank(),
    ...
  )
  
  
  # Label data
  if (all(selected_clusters != "all")) {
    nhood_data = nhood_data[nhood_data[[cluster]] %in% selected_clusters,]
  }
  
  
  nhood_rep = nhood_data %>% group_by(.data[[cluster]]) %>%
    summarise(umap1_center = median(.data[[nhood_1]]),
              umap2_center = median(.data[[nhood_2]]))
  
  nhood_data$logFC_cleaned = ifelse(abs(nhood_data$logFC) > 0 & nhood_data$SpatialFDR < 0.1, nhood_data$logFC, NA)
  
  # if(FC_dir == "enriched_depleted") {
  #   nhood_data$logFC_cleaned = ifelse(abs(nhood_data$logFC) > 0 & nhood_data$SpatialFDR < 0.1, nhood_data$logFC, NA)
  #   #remove non significant
  #   nhood_data$Nhood_size = ifelse(abs(nhood_data$logFC) > 0 & nhood_data$SpatialFDR < 0.1, nhood_data$Nhood_size, 1)
  # } else if(FC_dir == "enriched") {
  #   nhood_data$logFC_cleaned = ifelse(nhood_data$logFC > 0 & nhood_data$SpatialFDR < 0.1, nhood_data$logFC, NA)
  # 
  #   nhood_data$Nhood_size = ifelse(nhood_data$logFC > 0 & nhood_data$SpatialFDR < 0.1, nhood_data$Nhood_size, 1)
  # } else if(FC_dir == "depleted"){
  #   nhood_data$logFC_cleaned = ifelse(nhood_data$logFC < 0 & nhood_data$SpatialFDR < 0.1, nhood_data$logFC, NA)
  #   #remove non significant
  #   nhood_data$Nhood_size = ifelse(nhood_data$logFC < 0 & nhood_data$SpatialFDR < 0.1, nhood_data$Nhood_size, 1)
  # }
  
  #Plot
  if(!FC) {
    nhood_plot=nhood_data %>%
      ggplot(aes(
        x = .data[[nhood_1]],
        y = .data[[nhood_2]],
        color = .data[[cluster]],
        label = .data[[cluster]]
      )) +
      geom_point(show.legend = FALSE, size = shape.size) +
      scale_colour_manual(values = cluster_palette) +
      
      geom_label_repel(
        data = nhood_rep,
        aes(label = .data[[cluster]], x = umap1_center, y = umap2_center, fontface = font),
        point.padding = unit(0.5, "lines"), size = repel_text_size
      ) +
      
      geom_segment(
        data = nhood_data,
        aes(
          x = min(.data[[nhood_1]]),
          y = min(.data[[nhood_2]]),
          xend = min(.data[[nhood_1]]) + arrow.length,
          yend = min(.data[[nhood_2]])
        ),
        colour = "black",
        size = 1,
        arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
        linejoin = "mitre"
      ) +
      
      geom_segment(
        data = nhood_data,
        aes(
          x = min(.data[[nhood_1]]),
          y = min(.data[[nhood_2]]),
          xend = min(.data[[nhood_1]]) ,
          yend = min(.data[[nhood_2]]) + arrow.length
        ),
        colour = "black",
        size = 1,
        arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
        linejoin = "mitre"
      ) +
      
      annotate(
        "text",
        x = min(nhood_data[[nhood_1]]) + 1.5,
        y = min(nhood_data[[nhood_2]]) - 1,
        label = NHOOD_1_lab,
        color = "black",
        size = 3,
        fontface = "bold"
      ) +
      
      annotate(
        "text",
        x = min(nhood_data[[nhood_1]]) - 1,
        y = min(nhood_data[[nhood_2]]) + 1.5,
        label = NHOOD_2_lab,
        color = "black",
        size = 3,
        fontface = "bold" ,
        angle = 90
      ) +
      
      theme_nhood
   
    
  } else {
    lapply(unique(nhood_data$disease_short), function(dis){
      
      gg <- ggplot(nhood_data[nhood_data$disease_short == dis,], aes(
          x = .data[[nhood_1]],
          y = .data[[nhood_2]],
          color = logFC_cleaned,
          size = Nhood_size
        )) +
        geom_point() +
        labs(color = "logFC", size = "Nhood \nsize", title = dis) +
      #scale_fill_brewer(palette = cluster_palette) +
      scale_colour_gradient2(low = 'blue', mid = 'white', high = "#dc1c13", na.value = "grey") +
        #facet_wrap(~ disease_short, ncol = ncol, scales = "free") +
        theme_void() + theme(...)
      
      gg
      
      
    }) -> nhood_plot
    
    nhood_plot = cowplot::plot_grid(plotlist = nhood_plot, ncol = ncol)
    
    #nhood_plot[["ncol"]] <- ncol
    
    #do.call(gridExtra::grid.arrange, nhood_plot)
    
  }
  
  nhood_plot
}




plot_ggplot_gene2 = function(Normdata, gene, title, subtitle = "", xlab = "", ylab, pheno,
                             size = 2, shape = 20) {
  
  pheno = pheno
  
  rownames(pheno) = NULL
  rownames(pheno) = pheno[,1]
  
  colnames(pheno) = c("Sample", "Groups")
  
  pheno[,2] = factor(pheno[,2])
  
  pheno = pheno[match(rownames(pheno), colnames(Normdata)), ]
  
  data = data.frame(Gene = Normdata[gene, ], Groups = pheno[,2])
  
  min(data[,1]) 
  max(data[,1])
  
  yRange1 = 3.5
  
  yRange2 = 6
  
  ggplot(data, aes(Groups, Gene, color = Groups, fill = "white"))  +
    stat_summary(
      fun = median,
      fun.min = median,
      fun.max = median,
      geom = "crossbar",
      width = 0.6,
      color = "black",
      size = 0.2
    ) +
    scale_y_log10() +
    geom_quasirandom(
      size = size,
      shape = shape,
      bandwidth = 0.4,
      width = 0.2,
      dodge.width = 1
    ) +
    theme_classic() +
    labs(x = xlab, y = ylab, title = title) +
    scale_colour_brewer(palette = "Set1") +
    # theme_prism(palette = "shades_of_gray", base_size = 10) +
    theme(plot.title.position = "plot",
          legend.position = "none")
}


#' Function to plot GSEA barplot or lollipot plot
#' 
#' 
#' 

gsea_plot = function(gsea, n_top = NULL, x = "Term", y = "combined_score", log_y = TRUE, bar_plot = TRUE,
                     ylab="log2(Combined score)", xlab = "Hallmark GSEA",
                     figure_Path,
                     save=FALSE,
                     fig.width=10,
                     fig.height=8,
                     ...){
  
  if(y == "Adjusted P-value"){
    gsea = gsea %>% arrange(.data[[y]])
  }
  
  
  if(!is.null(n_top)) {
    gsea = na.omit(gsea[1:n_top,])
  }
  
  
  if(log_y & y != "Adjusted P-value"){
    gsea[[y]] = log2(gsea[[y]])
  } else if(y == "Adjusted P-value"){
    gsea[[y]] = -log10(gsea[[y]])
  }
  
  #clean term removing reactome accession
  gsea$Term = str_remove(gsea$Term, " R-HSA.*$")
  
  if(bar_plot){
    bplot=ggplot(gsea, aes(y = .data[[y]], x = reorder(.data[[x]], .data[[y]]), fill = "red")) +
      geom_bar(stat="identity") + 
      scale_x_discrete(expand = c(0, 0.5)) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(x=xlab, y=ylab) +
      coord_flip() +
      theme_light() +  
      theme(#axis_text_x = element_text(hjust=1, angle = 45),
        legend.position = "none",
        legend.text = element_text(size=10),
        panel.grid.minor = element_blank(),
        ...)
  }else{
    bplot=ggplot(gsea, aes(y = .data[[y]], x = reorder(.data[[x]], .data[[y]]))) +
      geom_point(color = "red", fill = "red", shape = 21, stroke = 1) +  
      coord_flip() +
      geom_segment(aes(xend = .data[[x]], x=.data[[x]], y=0, yend=.data[[y]]), colour="red") +
      #scale_x_discrete(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0.5)) +
      labs(x=xlab, y=ylab) +
      theme_bw() +
      theme(
        # panel.grid.major.x = element_blank(),
        # panel.grid.minor.x = element_blank(),
        #axis.title=element_text(size=14),
        #axis.text.y = element_text(size = 14)
        #panel.grid.major.y = element_line(colour = "grey60", linetype = "dashed")
        panel.background = element_rect(fill = 'transparent'),
        #transparent panel bg
        plot.background = element_rect(fill = 'transparent', color = NA),
        #transparent plot bg
        legend.background = element_rect(fill = 'transparent', colour = 'transparent', linewidth = 0),
        #transparent legend bg
        legend.box.background = element_rect(fill = 'transparent', colour = 'transparent', linewidth = 0),
        ...
      )
    
  }
  
  
  if(bar_plot){
    fig.type = "bar_plot"
  }else{
    fig.type = "lolli_plot"
  }
  
  bplot
  
  if(save) {
    ggsave(
      str_glue("{figure_Path}_{fig.type}.svg"),
      width = fig.width,
      height = fig.height,
      bg = 'transparent'
    )
    
    ggsave(
      str_glue("{figure_Path}_{fig.type}.png"),
      width = fig.width,
      height = fig.height,
      bg = 'transparent'
    )
    
    coln = c({{x}}, {{y}}, "Adjusted P-value")
    gsea = gsea[,colnames(gsea) %in% coln]
    
    write.csv(gsea, str_glue("{figure_Path}.csv"))
  }
  
  plot(bplot)
}


#' Get vector of genes with non-zero expression 
#' 
#' @param scedata SCE format
#' @param cell_type cell type to be used to subset the data. String or vector
#' @param tissue tissue to be used to subset the data. String or vector
#' @param cell_col column containing cell_type
#' @param tissue_col column containing tissue
#' @param cell_threshold threshold indicating the percentage of cell
#' @param dir direction of the comparison. "low" for lower than the indicated cell_threshold
#' "high" for higher or equal to the indicated cell_threshold


non_zero_genes_exprSCE = function(scedata, cell_type, cell_col, tissue_col, tissue, cell_threshold,
                                  dir = "high") {
  test_adata = scedata[, scedata[[cell_col]] %in% cell_type]
  table(test_adata[[tissue_col]])
  
  test_adata = test_adata[, test_adata[[tissue_col]] %in% tissue]
  
  counts(test_adata) = as.matrix(counts(test_adata))
  
  test_adata = logNormCounts(test_adata)
  
  test_adata = assay(test_adata, "logcounts")
  
  perc = data.frame(Perc = apply(test_adata, 1, function(x)
    sum(x > 0) / length(x)))
  
  if(dir == "high"){
    test_adata = test_adata[perc$Perc >= cell_threshold, ]
  }else if(dir == "low"){
    test_adata = test_adata[perc$Perc < cell_threshold, ]
  }
  
  
  return(rownames(test_adata))
}




#' ROC Curve
#' 
#' @param scedata sce containing log Normalized data
#' @param id Gene or ensembl ID
#' @param pheno_col vector with two string, columns. col 1: Sample; col 2: Groups
#' @param groups_col colnames of groups column
#' @param levels levels of included groups
#' @param cell_col column containing cell_types
#' @import pROC
#' #id = Results.RA_vs_Control$Gene.stable.ID
#' @export

#create ROC curve

roc_auc <- function(scedata, id, pheno_col, cell_type, cell_col, levels, groups_col, groups_included, main_groups_level) {
  scedata = scedata[,scedata[[cell_col]] %in% cell_type & scedata[[groups_col]] %in% groups_included]
  counts(scedata) = as.matrix(counts(scedata))
  
  scedata = logNormCounts(scedata)
  pheno = data.frame(Sample = colnames(scedata), Groups = scedata[[groups_col]])
  
  Normdata = assay(scedata, "logcounts")
  
  Normdata = Normdata[id, ]
  roc_data = data.frame(Gene=id, AUC=NA)
  rownames(roc_data) = roc_data$Gene
  
  pheno = pheno[colnames(Normdata),]
  
  roc_data = roc_data[match(rownames(roc_data), rownames(Normdata)),]
  
  
  pheno$Groups = ifelse(pheno$Groups %in% main_groups_level, 1L, 0L) %>% as.factor()
  
  
  pheno = pheno[match(rownames(pheno), colnames(Normdata)), ]
  
  for (i in 1:nrow(roc_data)) {
    roc_wt <- roc(pheno$Groups, as.numeric(t(Normdata)[,i]), quiet = TRUE)
    roc_data$AUC[i] = auc(roc_wt)
    
  }
  
  roc_data = arrange(roc_data, desc(AUC))
  return(roc_data)
  
}



#' Correlation matrix
#' 
#' 
#' @export

cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  
  z = ifelse(z < 0.05, "*", "ns")
  #dimnames(z) <- list(colnames(x), colnames(x))
  z
}


#' Reorder the correlation matrix
#' 
#' @export

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}



#' Cor matrix plot
#' 
#' @param ... parameter of pheatmap
#' @inheritParams pheatmap::pheatmap
#' @export
cor.mat.plot = function(mat, transpose = FALSE, title, labelmat=TRUE, ...) {
  
  if(transpose) {
    mat = t(mat)
  }
  # mat = reorder_cormat(mat)
  
  ComplexHeatmap::pheatmap(mat, 
                           display_numbers = labelmat,
                           cluster_cols = FALSE,
                           cluster_rows = FALSE,
                           main = title,
                           ...)
}



#'
#'
#'
#' @export


plot_ggplot_gene_col = function(sce = test_adata_, pheno = NULL, normalized.data = FALSE, gene, title, subtitle = "", xlab = "", ylab = "",
                                size = 2, shape = 20, cell_type, gg_x="tissue_disease", angle = 90, colmap,
                                base_size = 14, base_family = "sans", base_fontface = "bold", 
                                plot.title.position = "plot", legend.position = "right", show.x.text = FALSE) {
  
  base_line_size = base_size/14
  
  if(!normalized.data) {
    counts(sce) <- assay(sce)
    sce = logNormCounts(sce)
    sce = sce[, sce$Cell_type_scANVI %in% cell_type]
    data = data.frame(Gene = assay(sce, "logcounts")[gene,], as.data.frame(colData(sce)))
  } else {
    data = data.frame(Gene = sce[gene,], pheno)
  }
  
  
  if(angle==45){
    h.just=v.just=1
  } else {
    h.just = 1
    v.just = 0
  }
  
  
  plot = ggplot(data, aes(reorder(.data[[gg_x]], Gene, median), Gene, colour = .data[[gg_x]]))  +
    #geom_beeswarm() + 
    geom_point(
      position = position_quasirandom(), size=size) +
    # scale_y_log10() +
    stat_summary(
      fun = median,
      fun.min = median,
      fun.max = median,
      geom = "crossbar",
      width = size/5,
      color = "grey30",
      size = size/10
    ) +
    theme_classic() +
    labs(x = xlab, y = ylab, title = title, subtitle = subtitle) +
    scale_colour_manual(values = colmap)
    #ggprism::theme_prism(palette = "shades_of_gray", base_size = 10) +
  
  if(show.x.text){
    plot = plot + theme(plot.title.position = plot.title.position,
          legend.position = legend.position,
          axis.text.x = element_text(angle = angle, vjust = v.just, hjust = h.just, size = base_size, family = base_family, face = base_fontface,colour = "black"),
          axis.text.y = element_text(size = base_size, family = base_family, face = base_fontface, colour = "black"),
          legend.title = element_text(size = base_size, family = base_family, face = base_fontface, colour = "black"),
          legend.text = element_text(size = base_size, family = base_family, colour = "black"),
          axis.line = element_line(linewidth = base_line_size, colour = "black"),
          axis.ticks = element_line(linewidth = base_line_size, colour = "black")
    )
  }else{
    plot = plot + theme(plot.title.position = plot.title.position,
                        legend.position = legend.position,
                        axis.text.x = element_blank(),
                        axis.text.y = element_text(size = base_size, family = base_family, face = base_fontface, colour = "black"),
                        legend.title = element_text(size = base_size, family = base_family, face = base_fontface, colour = "black"),
                        legend.text = element_text(size = base_size, family = base_family, colour = "black"),
                        axis.line = element_line(linewidth = base_line_size, colour = "black"),
                        axis.ticks = element_blank()
    )
  }
  
  plot
  
}




plot_violin_facet = function(sce = test_adata_,
                             gene = c("NR1H3", "GBP2"),
                             pheno = NULL, normalized.data = FALSE,
                             title = "",
                             xlab = "",
                             ylab = "",
                             size = 2,
                             shape = 20,
                             cell_type = c("Mac_CX3CR1"),
                             gg_x = "disease_short",
                             angle = 90,
                             colmap = vega_20,
                             base_size = 14,
                             base_line_size = base_size/14,
                             base_family = "sans",
                             base_fontface = "bold",
                             legend.title = "Diseases",
                             plot.title.position = "plot",
                             plot.title ="",
                             legend.position = "right",
                             strip.colour="black",
                             strip.fill = "grey",
                             show.x.text = FALSE,
                             show.legend = TRUE,
                             ncol=6,
                             figure_Path,
                             save = FALSE,
                             fig.width = 10,
                             fig.height = 8,
                             stat=FALSE,
                             stat.ref="Healthy",
                             stat.method="t.test",
                             stat.paired = FALSE,
                             label.y = 0.9,
                             scale.data = FALSE,
                             violin = F,
                             include_violin_boxplot = FALSE,
                             only.point = FALSE,
                             crossbar.color = "black",
                             crossbar.width = 1,
                             crossbar.size = 1,
                             violin_boxplot_w = 0.1,
                             signif.size = 8,
                             label.y.npc = 1, 
                             y_expand_1 = 0,
                             y_expand_2 =0.1, ...) {
  
  
  if(!normalized.data) {
    counts(sce) <- assay(sce)
    
    sce = logNormCounts(sce)
    sce = sce[, sce$Cell_type_scANVI %in% cell_type]
    
    if (scale.data) {
      gene_data = data.frame(t(assay(sce, "logcounts")[gene,])) %>% t() %>%
        scale() %>% t() %>% as.data.frame() %>%
        mutate(ID = rownames(.)) %>%
        pivot_longer(.,
                     cols = gene,
                     values_to = "Expr",
                     names_to = "genes")
    } else {
      gene_data = data.frame(t(assay(sce, "logcounts")[gene, ])) %>% mutate(ID = rownames(.)) %>% pivot_longer(.,
                                                                                                               cols = gene,
                                                                                                               values_to = "Expr",
                                                                                                               names_to = "genes")
    }
    
    pheno = as.data.frame(colData(sce)) %>% mutate(ID = rownames(.))
  } else {
    gene_data = as.data.frame(t(sce[gene,])) %>% mutate(ID = rownames(.)) %>% pivot_longer(.,
                                                                           cols = gene,
                                                                           values_to = "Expr",
                                                                           names_to = "genes")
    pheno = pheno %>% mutate(ID = rownames(.))
  }
  
  data = inner_join(gene_data, pheno, by="ID")
  
  if(angle==45){
    h.just=v.just=1
  } else {
    h.just = 1
    v.just = 0
  }
  
  
  if(violin){
    plot = ggplot(data, aes(reorder(.data[[gg_x]], Expr, median), Expr, fill = .data[[gg_x]])) +
      geom_violin(trim=FALSE, show.legend = show.legend) +theme_classic()
    
    if(include_violin_boxplot){
      plot = plot + geom_boxplot(width=violin_boxplot_w, fill="white")
    }
      
    
  } else {
    if(only.point) {
      plot = ggplot(data, aes(reorder(.data[[gg_x]], Expr, median), Expr, fill = .data[[gg_x]]))  +
        geom_point(
          position = position_quasirandom(...),
          show.legend = F,
          shape = 21,
          size = size
        ) + 
        stat_summary(
          fun = median,
          fun.min = median,
          fun.max = median,
          geom = "crossbar",
          width = crossbar.width,
          color = crossbar.color,
          size = crossbar.size
        ) + theme_bw()
    } else {
      plot = ggplot(data, aes(reorder(.data[[gg_x]], Expr, median), Expr, fill = .data[[gg_x]]))  +
        geom_point(
          position = position_quasirandom(...),
          show.legend = F,
          shape = 21,
          size = size
        ) + geom_boxplot(show.legend = show.legend) + theme_bw()
    }
      
  }
  
  
  
  plot = plot +
    scale_y_continuous(expand = expansion(mult = c(y_expand_1, y_expand_2))) +
    #geom_point(position = position_quasirandom(), show.legend = F, shape = 21, size=size) +
    # scale_y_log10() +
    facet_wrap(~genes, scales = "free", ncol = ncol) +
    labs(x = xlab, y = ylab, fill = legend.title) +
    scale_colour_manual(values = colmap) +
    scale_fill_manual(values = colmap) + 
    ggtitle(label = plot.title)
  
  
  # ### Manual statistics
  # stat.test <- data %>%
  #   group_by(genes) %>%
  #   t_test(Expr ~ Cluster) %>%
  #   adjust_pvalue(method = "holm") %>%
  #   add_significance() %>% 
  #   add_y_position()
  # 
  # 
  # stat.test1 = lapply(gene, function(i){
  #   stat = data.frame(
  #     genes = i,
  #     ggpubr::compare_means(
  #       data = data[data$genes == i, ],
  #       formula = Expr ~ .data[[gg_x]],
  #       method = stat.method,
  #       paired = stat.paired,
  #       ref.group = stat.ref
  #     )
  #   )
  #   stat
  # }) %>% do.call(rbind, .)
  # 
  # # Make facet and add p-values
  # y_position =  data %>% group_by(genes) %>% 
  #   get_y_position(., formula = Expr ~ Cluster,
  #   ref.group = stat.ref,
  #   scales = "free"
  # )
  # 
  # stat.test1
  # 
  # stat.test1 = inner_join(stat.test, y_position[,c(1,4,5)], by="genes")
  
  if(stat){
    
    # plot + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01)
    plot = plot + ggpubr::stat_compare_means(label = "p.signif", method = stat.method,
                         ref.group = stat.ref, paired = stat.paired, label.y.npc = label.y.npc,
                         size = signif.size)
  }
  
  #ggprism::theme_prism(palette = "shades_of_gray", base_size = 10) +
  
  if(show.x.text){
    plot = plot + theme(plot.title.position = plot.title.position,
                        plot.title = element_text(size = base_size, family = base_family, face = base_fontface, colour = "black"),
                        legend.position = legend.position,
                        axis.text.x = element_text(angle = angle, vjust = v.just, hjust = h.just, size = base_size, family = base_family, face = base_fontface,colour = "black"),
                        axis.text.y = element_text(size = base_size, family = base_family, face = base_fontface, colour = "black"),
                        legend.title = element_text(size = base_size, family = base_family, face = base_fontface, colour = "black"),
                        legend.text = element_text(size = base_size, family = base_family, colour = "black"),
                        axis.line = element_line(linewidth = base_line_size, colour = "black"),
                        axis.ticks = element_line(linewidth = base_line_size, colour = "black"),
                        strip.text = element_text(size = base_size, family = base_family, face = base_fontface, colour = "black"),
                        strip.background=element_rect(colour=strip.colour, fill = strip.fill),
                        panel.background = element_rect(fill = 'transparent'),
                        #transparent panel bg
                        plot.background = element_rect(fill = 'transparent', color = NA),
                        #transparent plot bg
                        legend.background = element_rect(fill = 'transparent', colour = 'transparent', linewidth = 0),
                        #transparent legend bg
                        legend.box.background = element_rect(fill = 'transparent', colour = 'transparent', linewidth = 0)
    )
  }else{
    plot = plot + theme(plot.title.position = plot.title.position,
                        plot.title = element_text(size = base_size, family = base_family, face = base_fontface, colour = "black"),
                        legend.position = legend.position,
                        axis.text.x = element_blank(),
                        axis.text.y = element_text(size = base_size, family = base_family, face = base_fontface, colour = "black"),
                        legend.title = element_text(size = base_size, family = base_family, face = base_fontface, colour = "black"),
                        legend.text = element_text(size = base_size, family = base_family, colour = "black"),
                        axis.line = element_line(linewidth = base_line_size, colour = "black"),
                        axis.ticks.y = element_line(linewidth = base_line_size, colour = "black"),
                        axis.ticks.x = element_blank(),
                        strip.text = element_text(size = base_size, family = base_family, face = base_fontface, colour = "black"),
                        strip.background=element_rect(colour=strip.colour, fill = strip.fill),
                        panel.background = element_rect(fill = 'transparent'),
                        #transparent panel bg
                        plot.background = element_rect(fill = 'transparent', color = NA),
                        #transparent plot bg
                        legend.background = element_rect(fill = 'transparent', colour = 'transparent', linewidth = 0),
                        #transparent legend bg
                        legend.box.background = element_rect(fill = 'transparent', colour = 'transparent', linewidth = 0)
                        
    )
  }
  
  plot
  
  if(save) {
   
    ggsave(
      str_glue("{figure_Path}.png"),
      width = fig.width,
      height = fig.height,
      bg = 'transparent',
      dpi = 300
    )
    
  }
  
  plot
  
}


#' Abundance plot
#' 
#' 
#' 

plot_abundance_bar = function(data,
                              count.on.top = TRUE,
                              text.size = 5,
                              axis.text.size = 20,
                              key.size = 0.7,
                              fill_column = "Cell_type_scANVI2",
                              y_column = "Cell_count",
                              x_column = "disease_short_2",
                              y_label = "Composition (%)",
                              fill_label = "Cellular compartment",
                              margin.bottom = 0,
                              margin.top = 50,
                              margin.right = 0,
                              margin.left = 20,
                              Palette) {
  
  
  temp = data %>%
    group_by(.data[[x_column]]) %>% 
    mutate(total=sum(.data[[y_column]])) %>% 
    dplyr::select(all_of(x_column), total) %>% unique()
  
  #colnames(temp)[1] = "X_text"
  
  if(count.on.top) {
    ggplot(data = data,
           aes(y = .data[[y_column]], x = as.factor(.data[[x_column]]))) +
      geom_text(
        data = temp,
        aes(x = .data[[x_column]], y = Inf, label = .data[["total"]]),
        vjust = 0,
        hjust = 0,
        angle = 45,
        size = text.size
      ) +
      coord_cartesian(clip = "off") +
      geom_bar(aes(fill = .data[[fill_column]]), position = "fill", stat = "identity") +
      #geom_col(width = 0.4) +
      labs(x = "", y = y_label, fill = fill_label) +
      scale_y_continuous(expand = c(0, 0.005)) +
      scale_fill_manual(values = Palette) +
      theme_classic() + #coord_flip() +
      theme(
        axis.text.x = element_text(
          hjust = 1,
          vjust = 1,
          size = 20,
          angle = 45
        ),
        axis.text.y = element_text(size = axis.text.size),
        plot.margin = margin(margin.top, margin.right, margin.bottom, margin.left, unit = "pt"),
        axis.title = element_text(size = axis.text.size),
        legend.key.size = unit(key.size, 'cm'),
        legend.text = element_text(size = axis.text.size),
        axis.line = element_line(colour = "black", linewidth = 0.6),
        legend.key.width = unit(key.size, "cm"),
        legend.key.height = unit(key.size, "cm"),
        legend.title = element_text(size = axis.text.size),
        panel.background = element_rect(fill = 'transparent'),
        #transparent panel bg
        plot.background = element_rect(fill = 'transparent', color = NA),
        #transparent plot bg
        legend.background = element_rect(fill = 'transparent'),
        #transparent legend bg
        #legend.box.background = element_rect(fill = 'transparent'),
        #legend.box = element_blank()
      ) +
      guides(fill = guide_legend(title.position = "top"))
  } else {
    ggplot(data = data,
           aes(fill = .data[[fill_column]], y = .data[[y_column]], x = .data[[x_column]])) +
      geom_bar(position = "fill", stat = "identity") +
      #geom_col(width = 0.4) +
      labs(x = "", y = y_label, fill = fill_label) +
      scale_y_continuous(expand = c(0, 0.005)) +
      scale_fill_manual(values = Palette) +
      theme_classic() + #coord_flip() +
      theme(
        axis.text.x = element_text(
          hjust = 1,
          vjust = 1,
          size = 20,
          angle = 45
        ),
        axis.text.y = element_text(size = axis.text.size),
        axis.title = element_text(size = axis.text.size),
        legend.key.size = unit(key.size, 'cm'),
        legend.text = element_text(size = axis.text.size),
        axis.line = element_line(colour = "black", linewidth = 0.6),
        legend.key.width = unit(key.size, "cm"),
        legend.key.height = unit(key.size, "cm"),
        legend.title = element_text(size = axis.text.size),
        panel.background = element_rect(fill = 'transparent'),
        #transparent panel bg
        plot.background = element_rect(fill = 'transparent', color = NA),
        #transparent plot bg
        legend.background = element_rect(fill = 'transparent'),
        #transparent legend bg
        #legend.box.background = element_rect(fill = 'transparent'),
        #legend.box = element_blank()
      ) +
      guides(fill = guide_legend(title.position = "top"))
  }
  
  
}



# =============================================================================
# study_composition_ch()
#
# ComplexHeatmap version of the study-composition panel:
#   heatmap : disease x tissue, fill = log10(cells), value printed in cell
#   right   : barplot of total cells (or donors) per disease + printed total
#
# Self-contained: does not depend on the patchwork/ggplot version.
#
# Input: a per-cell metadata table (e.g. `umap2`) containing
#        tissue_harmo, disease_label_trans, donor_id, dataset_id
# =============================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(scales)
library(grid)
library(circlize)
library(ComplexHeatmap)


# -- helpers ------------------------------------------------------------------

# scales::cut_short_scale() errors on NA input, so format only the real values
.fmt_raw <- scales::label_number(accuracy = 1, scale_cut = scales::cut_short_scale())

.fmt_n <- function(x) {
  out <- rep(NA_character_, length(x))
  ok  <- !is.na(x)
  if (any(ok)) out[ok] <- .fmt_raw(x[ok])
  out
}

# long counts table, shared by both the matrix and the annotations
prep_counts <- function(data,
                        disease_col = "disease_label_trans",
                        tissue_col  = "tissue_harmo",
                        donor_col   = "donor_id",
                        dataset_col = "dataset_id",
                        min_cells   = 0) {
  
  needed <- c(disease_col, tissue_col, donor_col, dataset_col)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols)) {
    stop("Column(s) not found in `data`: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }
  
  # transmute, not mutate: `data` may already have columns named disease/tissue
  meta <- data %>%
    transmute(
      disease    = as.character(.data[[disease_col]]),
      tissue     = as.character(.data[[tissue_col]]),
      donor_id   = .data[[donor_col]],
      dataset_id = .data[[dataset_col]]
    ) %>%
    filter(!is.na(disease), !is.na(tissue))
  
  counts <- meta %>%
    group_by(disease, tissue) %>%
    summarise(
      n_cells   = n(),
      n_donors  = n_distinct(donor_id),
      n_studies = n_distinct(dataset_id),
      .groups   = "drop"
    ) %>%
    filter(n_cells >= min_cells) %>%
    complete(disease, tissue)          # NA cells = tissue not profiled
  
  # totals restricted to the surviving combinations, so the bar always agrees
  # with the sum of the printed values. Donors counted from `meta` directly --
  # summing n_donors across tissues would double count multi-tissue donors.
  keep <- counts %>% filter(!is.na(n_cells)) %>% distinct(disease, tissue)
  
  totals <- meta %>%
    semi_join(keep, by = c("disease", "tissue")) %>%
    group_by(disease) %>%
    summarise(
      total_cells  = n(),
      total_donors = n_distinct(donor_id),
      .groups      = "drop"
    )
  
  list(counts = counts, totals = totals, meta = meta)
}


# -- main function ------------------------------------------------------------

# base-R label wrapper: no stringr dependency. strwrap() breaks on whitespace
# only, so a single long token (e.g. "hepatosplenomegaly") stays on one line.
.wrap_lab <- function(x, width) {
  if (is.null(width) || !is.finite(width)) return(x)
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
}

study_composition_ch <- function(
    data           = umap2,
    disease_col    = "disease_harmo",
    tissue_col     = "tissue_harmo",
    donor_col      = "donor_id",
    dataset_col    = "dataset_id",
    min_cells      = 0,
    show_donors    = TRUE,
    format_values  = FALSE,    # TRUE -> "420k", FALSE -> "420318"
    cluster_rows   = FALSE,
    cluster_columns = FALSE,
    legend_at      = c(1e3, 1e4, 1e5, 5e5),
    low_col        = "#E2E2E2",
    mid_col        = "#BDBDBD",
    high_col       = "#6A6A6A",
    cell_fontsize  = 7.5,
    cell_width     = 1,
    cell_height    = 1,
    lwd_size       = 1.2,
    name_fontsize  = 9,
    wrap_rows      = NULL,     # chars per line for row names; NULL = no wrapping
    wrap_cols      = NULL,     # chars per line for column names; NULL = no wrapping
    col_names_rot  = 45,       # set to 0 when wrap_cols is used
    na_col         = "white",
    row_title      = NULL,
    column_title   = NULL,
    row_name_alpha = TRUE,
    flip           = FALSE     # FALSE: disease = rows | TRUE: disease = columns
) {
  
  prepped <- prep_counts(data, disease_col, tissue_col,
                         donor_col, dataset_col, min_cells)
  counts <- prepped$counts
  totals <- prepped$totals
  n_donors_tot <- dplyr::n_distinct(prepped$meta$donor_id)
  n_cells_tot  <- nrow(prepped$meta)
  
  # ---- wide matrices -------------------------------------------------------
  mat_cells <- counts %>%
    dplyr::select(disease, tissue, n_cells) %>%
    pivot_wider(names_from = tissue, values_from = n_cells) %>%
    column_to_rownames("disease") %>%
    as.matrix()
  
  mat_donors <- counts %>%
    dplyr::select(disease, tissue, n_donors) %>%
    pivot_wider(names_from = tissue, values_from = n_donors) %>%
    column_to_rownames("disease") %>%
    as.matrix()
  mat_donors <- mat_donors[rownames(mat_cells), colnames(mat_cells), drop = FALSE]
  
  mat_log <- log10(mat_cells)
  
  # ---- ordering ------------------------------------------------------------
  # ComplexHeatmap draws the first row at the TOP, the opposite of ggplot,
  # so sort descending to get the largest disease at the top.
  if (row_name_alpha) {
    row_ord <- sort(unique(counts$disease))
  } else {
    row_ord <- totals %>% arrange(desc(total_cells)) %>% pull(disease)
  }
  row_ord <- row_ord[row_ord %in% rownames(mat_cells)]
  
  col_ord <- counts %>%
    group_by(tissue) %>%
    summarise(tot = sum(n_cells, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(tot)) %>%
    pull(tissue)
  
  mat_cells  <- mat_cells[row_ord, col_ord, drop = FALSE]
  mat_donors <- mat_donors[row_ord, col_ord, drop = FALSE]
  mat_log    <- mat_log[row_ord, col_ord, drop = FALSE]
  
  totals <- totals[match(row_ord, totals$disease), ]
  
  # keep the true dimensions before any transpose, for the return value
  n_disease <- length(row_ord)
  n_tissue  <- length(col_ord)
  
  if (is.null(column_title)) {
    column_title <- sprintf("%s cells from %s donors \n across %s tissues",
                            n_cells_tot, n_donors_tot, n_tissue)
  }
  
  # ---- flip ----------------------------------------------------------------
  # transpose all three matrices together so cell_fun's (i, j) stay aligned
  # with the matrix that Heatmap() is actually drawing.
  if (flip) {
    mat_cells  <- t(mat_cells)
    mat_donors <- t(mat_donors)
    mat_log    <- t(mat_log)
  }
  
  # ---- wrap displayed labels ----------------------------------------------
  # applied AFTER the flip, so wrap_rows always refers to the drawn rows.
  # cell_fun indexes mat_cells positionally, so its labels are unaffected.
  rownames(mat_log) <- .wrap_lab(rownames(mat_log), wrap_rows)
  colnames(mat_log) <- .wrap_lab(colnames(mat_log), wrap_cols)
  
  # dist() cannot handle NA; cluster on a zero-filled copy if asked
  mat_clust <- mat_log
  mat_clust[is.na(mat_clust)] <- 0
  
  # ---- colour scale --------------------------------------------------------
  rng <- range(mat_log, na.rm = TRUE)
  col_fun <- colorRamp2(
    c(rng[1], mean(rng), rng[2]),
    c(low_col, mid_col, high_col)
  )
  ink_cut <- rng[1] + 0.62 * diff(rng)     # above this, white text
  
  legend_at <- legend_at[legend_at >= 10^rng[1] & legend_at <= 10^rng[2]]
  if (!length(legend_at)) legend_at <- 10^pretty(rng, 3)
  
  # ---- in-cell labels ------------------------------------------------------
  fmt_cell <- if (format_values) .fmt_n else function(x) format(x, big.mark = "")
  
  cell_fun <- function(j, i, x, y, w, h, fill) {
    #grid.rect(x, y, w, h, gp = gpar(col = "black", lwd = 1.2, fill = NA))
    v <- mat_cells[i, j]
    if (is.na(v)) return(invisible(NULL))
    txt <- if (show_donors) paste0(fmt_cell(v), "\nn=", mat_donors[i, j]) else fmt_cell(v)
    grid.text(txt, x, y, gp = gpar(fontsize = cell_fontsize,
                                   col = if (log10(v) > ink_cut) "white" else "grey10",
                                   lineheight = 0.95))
  }
  
  # ---- heatmap -------------------------------------------------------------
  # cell_width / cell_height always refer to the DRAWN cell, so they read from
  # the post-transpose matrix and need no swapping.
  ht <- Heatmap(
    if (cluster_rows || cluster_columns) mat_clust else mat_log,
    width  = unit(cell_width,  "cm") * ncol(mat_log),
    height = unit(cell_height, "cm") * nrow(mat_log),
    name            = "Cells",
    col             = col_fun,
    na_col          = na_col,
    rect_gp         = gpar(col = "black", lwd = lwd_size),
    cluster_rows    = cluster_rows,
    cluster_columns = cluster_columns,
    row_order       = if (cluster_rows)    NULL else seq_len(nrow(mat_log)),
    column_order    = if (cluster_columns) NULL else seq_len(ncol(mat_log)),
    row_names_side        = "left",
    column_names_side     = "top",
    column_names_rot      = col_names_rot,
    column_names_centered = (col_names_rot == 0),
    # rotated single-line labels need WIDTH reserved; wrapped horizontal
    # labels need HEIGHT. Measuring the wrong one clips the text.
    column_names_max_height = if (col_names_rot != 0) {
      max_text_width(colnames(mat_log), gp = gpar(fontsize = name_fontsize))
    } else {
      max_text_height(colnames(mat_log), gp = gpar(fontsize = name_fontsize)) * 1.15
    },
    row_names_max_width = max_text_width(rownames(mat_log),
                                         gp = gpar(fontsize = name_fontsize)) * 1.05,
    row_names_gp    = gpar(fontsize = name_fontsize),
    column_names_gp = gpar(fontsize = name_fontsize),
    row_title       = row_title,
    column_title    = column_title,
    column_title_gp = gpar(fontsize = name_fontsize),
    cell_fun        = cell_fun,
    heatmap_legend_param = list(
      at        = log10(legend_at),
      labels    = .fmt_n(legend_at),
      title     = "Cells",
      direction = "horizontal",
      legend_width = unit(4, "cm"),
      title_position = "topcenter",
      labels_gp = gpar(fontsize = name_fontsize - 1),
      title_gp  = gpar(fontsize = name_fontsize)
    )
  )
  
  invisible(list(
    heatmap = ht,
    counts  = counts,
    totals  = totals,
    mat     = mat_cells,
    flip    = flip,
    dims    = c(disease = n_disease, tissue = n_tissue)
  ))
}





# -- draw + export ------------------------------------------------------------

draw_composition_ch <- function(res, ...) {
  draw(res$heatmap,
       heatmap_legend_side    = "bottom",
       annotation_legend_side = "bottom",
       merge_legend           = TRUE,
       ...)
}

save_composition_ch <- function(res, file, width = NULL, height = NULL,
                                res_dpi = 300, editable = TRUE, background) {
  n_dis <- res$dims[["disease"]]
  n_tis <- res$dims[["tissue"]]
  w <- width  %||% (3.2 + 0.85 * n_tis)
  h <- height %||% (2.0 + 0.42 * n_dis)
  
  if (grepl("\\.svg$", file)) {
    svglite::svglite(
      file,
      width  = w,
      height = h,
      # TRUE locks glyph spacing with textLength; FALSE keeps text editable
      fix_text_size = !editable,
      system_fonts  = list(sans = "Helvetica"),
      bg = background
    )
  } else if (grepl("\\.pdf$", file)) {
    cairo_pdf(file, width = w, height = h)
  } else if (grepl("\\.(tif|tiff)$", file)) {
    tiff(file, width = w, height = h, units = "in", res = res_dpi, compression = "lzw")
  } else {
    png(file, width = w, height = h, units = "in", res = res_dpi)
  }
  on.exit(dev.off(), add = TRUE)
  draw_composition_ch(res)
  invisible(file)
}


