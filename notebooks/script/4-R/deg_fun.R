#' Function to determine cell gene program


run_cell_gene_program = function(test_adata,
                                 
                                 dir_fc,
                                 
                                 lfc = 1,
                                 
                                 p.val = 0.05,
                                 
                                 to_all = FALSE,
                                 
                                 cel_type = "",
                                 
                                 tissue_disease = c()) {
  test_adata = test_adata[, test_adata$tissue_disease %in% tissue_disease]
  
  counts(test_adata) <- assay(test_adata)
  
  test_adata = logNormCounts(test_adata)
  
  # genes where at least some counts are more than 10 count
  
  non_empty_rows <- which(rowSums2(assay(test_adata)) > 10)
  
  test_adata <- test_adata[non_empty_rows, ]
  
  ## Feature selection with scran. 10000 hvgs will be selected
  
  dec <- modelGeneVar(test_adata)
  
  hvgs <- getTopHVGs(dec, n = 5000)
  
  test_adata <- test_adata[hvgs,]
  
  counts(test_adata) <- as.matrix(counts(test_adata))
  
  assay(test_adata, "X") = as.matrix(assay(test_adata, "X"))
  
  test_adata$disease_short = factor(test_adata$disease_short) %>% droplevels()
  
  test_adata$nhoo_test_condition = droplevels(test_adata$nhoo_test_condition)
  
  test_adata$tissue_disease = factor(test_adata$tissue_disease) %>% droplevels()
  
  test_adata$donor_id = factor(test_adata$donor_id) %>% droplevels()
  
  test_adata$dataset_id = factor(test_adata$dataset_id) %>% droplevels()
  
  test_adata$Cell_type_scANVI = str_replace(test_adata$Cell_type_scANVI, "/", "_") %>% factor()
  
  test_adata$Cell_type_scANVI = factor(test_adata$Cell_type_scANVI) %>% droplevels()
  
  print(table(test_adata$disease_short, test_adata$nhoo_test_condition))
  
  # Comparison will be done only with the cell of interest
  
  if (!to_all) {
    test_adata = test_adata[, test_adata$Cell_type_scANVI == cel_type]
    
    test_adata$disease_nhoo = ifelse(test_adata$nhoo_test_condition == dir_fc,
                                     
                                     dir_fc,
                                     
                                     "Rest") %>% factor()
    
    # Prepare data for comparison with all other cells
    
  } else {
    test_adata$disease_nhoo = ifelse(
      test_adata$nhoo_test_condition == dir_fc &
        test_adata$Cell_type_scANVI == cel_type,
      
      dir_fc,
      "Rest"
    ) %>% factor()
    
  }
  
  my_design = ~ 1 + disease_nhoo
  
  fit1 <-
    glm_gp(test_adata, design = my_design, reference_level = "Rest")
  
  de_res1 <-
    test_de(fit1, contrast = str_glue("disease_nhoo{dir_fc}"))
  
  de_res1$dirFC = ifelse(
    de_res1$lfc >= lfc & de_res1$adj_pval < p.val,
    "Up",
    ifelse(de_res1$lfc <= -lfc & de_res1$adj_pval < p.val,
           "Down",
           "NS")
  )
  
  de_res = de_res1 %>% filter(dirFC != "NS") %>% arrange(desc(lfc))
  de_res
  
}

#' Disease gene program
#' 


disease_gene_program = function(test_adata,
                                
                                lfc = 1,
                                
                                dir_fc,
                                
                                p.val = 0.05,
                                
                                cel_type="",
                                var_column = "tissue_disease",
                                var.included = c()) {
  test_adata = test_adata[, test_adata[[var_column]] %in% var.included]
  
  
  # Select only one cell type of interest
  test_adata = test_adata[, test_adata$Cell_type_scANVI %in% cel_type]
  
  
  counts(test_adata) <- assay(test_adata)
  
  
  
  test_adata = logNormCounts(test_adata)
  
  
  
  # genes where at least some counts are more than 10 count
  
  non_empty_rows <- which(rowSums2(assay(test_adata)) > 10)
  
  test_adata <- test_adata[non_empty_rows,]
  
  ## Feature selection with scran. 10000 hvgs will be selected
  
  dec <- modelGeneVar(test_adata)
  
  hvgs <- getTopHVGs(dec, n = 5000)
  
  test_adata <- test_adata[hvgs, ]
  
  
  
  counts(test_adata) <- as.matrix(counts(test_adata))
  
  
  
  assay(test_adata, "X") = as.matrix(assay(test_adata, "X"))
  
  
  
  test_adata$disease_short = factor(test_adata$disease_short) %>% droplevels()
  
  
  
  test_adata$nhoo_test_condition = droplevels(test_adata$nhoo_test_condition)
  
  test_adata[[var_column]] = factor(test_adata[[var_column]]) %>% droplevels()
  
  test_adata$donor_id = factor(test_adata$donor_id) %>% droplevels()
  
  test_adata$dataset_id = factor(test_adata$dataset_id) %>% droplevels()
  
  test_adata$Cell_type_scANVI = str_replace(test_adata$Cell_type_scANVI, "/", "_") %>% factor()
  
  
  
  test_adata$Cell_type_scANVI = factor(test_adata$Cell_type_scANVI) %>% droplevels()
  
  print(table(test_adata[[var_column]], test_adata$nhoo_test_condition))
  
  test_adata$disease_test = ifelse(
    test_adata$disease_short == "Healthy",
    "Healthy",
    ifelse(
      test_adata$disease_short != "Healthy" &
        test_adata$nhoo_test_condition == dir_fc,
      
      "Disease_nhoo",
      "Disease_rest"
    )
  ) %>% factor()
  
  
  my_design = ~ 1 + disease_test
  
  fit1 <- glm_gp(test_adata, design = my_design,
                 reference_level = "Healthy")
  
  
  de_res1 <- test_de(fit1, contrast = 'disease_testDisease_nhoo')
  
  de_res1$dirFC = ifelse(
    de_res1$lfc >= lfc & de_res1$adj_pval < p.val,
    
    "Up",
    
    ifelse(de_res1$lfc <= -lfc & de_res1$adj_pval < p.val,
           
           "Down",
           
           "NS")
    
  )
  
  de_res2 <- test_de(fit1, contrast = ('disease_testDisease_rest'))
  
  de_res2$dirFC = ifelse(
    de_res2$lfc >= lfc & de_res2$adj_pval < p.val,
    
    "Up",
    
    ifelse(de_res2$lfc <= -lfc & de_res2$adj_pval < p.val,
           
           "Down",
           
           "NS")
    
  )
  
  
  de_res1$dirFC2 = de_res2$dirFC
  
  de_res = de_res1 %>% filter(dirFC != "NS") %>% filter(dirFC != dirFC2) %>% arrange(desc(lfc))
  
  de_res
  
}

#' Cell gene program vs all


cell_gene_program_vs_all = function(test_adata,
                                    lfc = 1,
                                    p.val = 0.05,
                                    cel_type="",
                                    tissue_disease = c(),
                                    covariable.dataset = FALSE,
                                    hvgs = TRUE,
                                    remove_empty = TRUE) {
  
  test_adata = test_adata[, test_adata$tissue_disease %in% tissue_disease]
  
  counts(test_adata) <- assay(test_adata)
  
  test_adata = logNormCounts(test_adata)
  
  # genes where at least some counts are more than 10 count
  
  if(remove_empty){
    non_empty_rows <- which(rowSums2(assay(test_adata)) > 10)
    test_adata <- test_adata[non_empty_rows, ]
  }
  
  
  ## Feature selection with scran. 10000 hvgs will be selected
  if(hvgs){
    dec <- modelGeneVar(test_adata)
    
    hvgs <- getTopHVGs(dec, n = 5000)
    
    test_adata <- test_adata[hvgs,]
  }
  
  
  counts(test_adata) <- as.matrix(counts(test_adata))
  
  assay(test_adata, "X") = as.matrix(assay(test_adata, "X"))
  
  
  test_adata$disease_short = factor(test_adata$disease_short) %>% droplevels()
  
  test_adata$nhoo_test_condition = droplevels(test_adata$nhoo_test_condition)
  
  test_adata$tissue_disease = factor(test_adata$tissue_disease) %>% droplevels()
  
  test_adata$donor_id = factor(test_adata$donor_id) %>% droplevels()
  
  test_adata$dataset_id = factor(test_adata$dataset_id) %>% droplevels()
  
  test_adata$Cell_type_scANVI = str_replace(test_adata$Cell_type_scANVI, "/", "_") %>% factor()
  
  test_adata$Cell_type_scANVI = factor(test_adata$Cell_type_scANVI) %>% droplevels()
  
  print(table(test_adata$disease_short, test_adata$nhoo_test_condition))
  
  # Comparison will be done only with the cell of interest
  
  test_adata$disease_cell = ifelse(test_adata$Cell_type_scANVI == cel_type, "DiseaseCell", "Rest") %>% factor()
  
  
  if(covariable.dataset){
    my_design = ~ 1 + disease_cell + dataset_id
  } else {
    my_design = ~ 1 + disease_cell
  }
  
  fit1 <- glm_gp(test_adata,
                 
                 design = my_design,
                 
                 reference_level = "Rest", size_factors = test_adata$sizeFactor)
  
  de_res1 <- test_de(fit1, contrast = 'disease_cellDiseaseCell')
  
  de_res1$dirFC = ifelse(
    
    de_res1$lfc >= lfc & de_res1$adj_pval < p.val,
    
    "Up",
    
    ifelse(
      
      de_res1$lfc <= -lfc & de_res1$adj_pval < p.val,
      
      "Down",
      
      "NS"
      
    )
    
  )
  
  
  de_res = de_res1 %>% filter(dirFC != "NS") %>% arrange(desc(lfc))
  de_res
  
}




#' ROC Curve
#' 
#' @param Normdata Normalized data
#' @param id Gene or ensembl ID
#' @param pheno data frame with two colonnes. col 1: Sample; col 2: Groups
#' 
#' @import dplyr
#' @import pROC
#' #id = Results.RA_vs_Control$Gene.stable.ID
#' @export

#create ROC curve

roc_auc_AIID <- function(Normdata, id, pheno, levels = c("Rest", "GOI")) {
  Normdata = Normdata[id, ]
  roc_data = data.frame(Gene.stable.ID=id, AUC=NA)
  rownames(roc_data) = roc_data$Gene.stable.ID
  roc_data = roc_data[match(rownames(roc_data), rownames(Normdata)),]
  
  pheno = pheno
  
  rownames(pheno) = NULL
  rownames(pheno) = pheno[,1]
  
  colnames(pheno) = c("Sample", "Groups")
  
  pheno[,2] = factor(pheno[,2], levels = levels)
  
  pheno = pheno[match(rownames(pheno), colnames(Normdata)), ]
  
  for (i in 1:nrow(roc_data)) {
    roc_wt <- pROC::roc(as.factor(pheno[,2]), as.numeric(t(Normdata)[,i]), levels = levels, quiet = TRUE)
    roc_data$AUC[i] = pROC::auc(roc_wt)
    #cat("Processed", i, nrow(roc_data), "\n")
  }
  
  roc_data = arrange(roc_data, desc(AUC))
  return(roc_data)
  
}




#' Function to compare two groups within a cell type
#' 
#' @param test_adata SingleCellExperiment format
#' 
#' @param lfc log2 fold change threshold
#' @param var.included category to be included in the comparison. For example tissue disease combination to be included
#' @param var_column column names of column containing var.included
#' @param  cel_type specific cell types
#' @param ctrl control groups
#' @param test_col column containing comparison factor


run_deg = function(test_adata=test_adata_,
                   lfc = 1,
                   p.val = 0.05,
                   cel_type = "",
                   test_col,
                   ctrl,
                   var_column = "tissue_disease",
                   var.included = c(),
                   covariable.dataset = FALSE,
                   levels = c("Rest", "GOI"),
                   covar = "disease_short",
                   hvg =TRUE,
                   auc = TRUE) {
  
  test_adata = test_adata[, test_adata[[var_column]] %in% var.included & test_adata$Cell_type_scANVI %in% cel_type]
  
  counts(test_adata) <- assay(test_adata)
  
  test_adata = logNormCounts(test_adata)
  
  # genes where at least some counts are more than 10 count
  
  non_empty_rows <- which(rowSums2(assay(test_adata)) > 10)
  
  test_adata <- test_adata[non_empty_rows,]
  
  ## Feature selection with scran. 10000 hvgs will be selected
  if(hvg) {
    dec <- modelGeneVar(test_adata)
    
    hvgs <- getTopHVGs(dec, n = 5000)
    
    test_adata <- test_adata[hvgs,]
  }
  
  counts(test_adata) <- as.matrix(counts(test_adata))
  
  assay(test_adata, "X") = as.matrix(assay(test_adata, "X"))
  
  test_adata$disease_short = factor(test_adata$disease_short) %>% droplevels()
  
  test_adata$nhoo_test_condition = droplevels(test_adata$nhoo_test_condition)
  
  test_adata[[var_column]] = factor(test_adata[[var_column]]) %>% droplevels()
  
  test_adata$donor_id = factor(test_adata$donor_id) %>% droplevels()
  
  test_adata$dataset_id = factor(test_adata$dataset_id) %>% droplevels()
  
  test_adata$Cell_type_scANVI = str_replace(test_adata$Cell_type_scANVI, "/", "_") %>% factor()
  
  test_adata$Cell_type_scANVI = factor(test_adata$Cell_type_scANVI) %>% droplevels()
  
  test_adata$disease_test=test_adata[[test_col]]
  
  print(table(test_adata$disease_test, test_adata$nhoo_test_condition))
  
  # Consider "colon_Uninfl_CD" Control
  if("Uninfl_CD" %in% levels(test_adata$disease_test)){
    test_adata$disease_test = str_replace(test_adata$disease_test, "Uninfl_CD", "Healthy")
    print("Uninfl_CD have been transformed to Healthy")
    print(table(test_adata$disease_test, test_adata$nhoo_test_condition))
    
  }
  test_adata$disease_bin = ifelse(test_adata$disease_short=="Healthy", "Healthy", "Disease") %>% factor(levels = c("Healthy", "Disease"))
  
  
  test_adata$Cell_disease_bin = paste0(test_adata$Cell_type_scANVI, "_", test_adata$disease_bin)
  
  test_adata$disease_test = ifelse(test_adata$disease_test%in%ctrl, "Rest", "GOI") #GOI ==> group of interest
  
  
  test_adata$disease_test = factor(test_adata$disease_test) %>% droplevels()
  
  print(table(test_adata$disease_test, test_adata$nhoo_test_condition))
  
  # Comparison between two groups
  
  if(covariable.dataset){
    test_adata$covar = test_adata[[covar]]
    my_design = ~ disease_test + covar
  } else {
    my_design = ~ disease_test
  }
  
  #my_design = ~ 1 + disease_test
  
  # comparison Enriched vs depleted
  
  
  fit1 <- glm_gp(test_adata, design = my_design, reference_level = "Rest", size_factors = test_adata$sizeFactor)
  
  #colnames(fit1$Beta)
  
  de_res1 <- test_de(fit1, contrast = disease_testGOI)
  
  de_res1$dirFC = ifelse(
    de_res1$lfc >= lfc & de_res1$adj_pval < p.val,
    "Up",
    ifelse(de_res1$lfc <= -lfc & de_res1$adj_pval < p.val,
           "Down",
           "NS")
  )
  
  de_res = de_res1 %>% filter(dirFC != "NS") %>% arrange(desc(lfc))
  
  #### ROC
  
  normdata = assay(test_adata, "logcounts")[de_res$name,]
  
  pheno=colData(test_adata) %>% as.data.frame()
  pheno$sample = rownames(pheno)
  pheno = dplyr::select(pheno, sample, disease_test)
  
  #all.equal(rownames(pheno), colnames(normdata))
  if(auc){
    de_roc = roc_auc_AIID(Normdata = normdata, id=rownames(normdata), pheno=pheno, levels=factor(unique(pheno$disease_test), levels = levels))
    
    colnames(de_roc)[1]="name"
    
    de_res = inner_join(de_res, de_roc, by = "name")
  }
  
  
  
  # tmp <- data.frame(gene = rep(de_res$name[1:6], times = ncol(test_adata)),
  #                   expression = c(assay(test_adata, "logcounts")[de_res$name[1:6], ]),
  #                   celltype = rep(paste0(test_adata$Cell_type_scANVI, "_", test_adata$disease_bin), each = 6))
  # 
  # ggplot(tmp, aes(x = celltype, y = expression)) +
  #   geom_jitter(height = 0.1) +
  #   stat_summary(geom = "crossbar", fun = "mean", color = "red") +
  #   facet_wrap(~ gene, scales = "free_y") +
  #   ggtitle("cMon_1 vs. cMon_2")
  
  
  de_res
  
}





#' Disease gene program
#' 


disease_nhoo_program = function(test_adata,
                                
                                lfc = 1,
                                
                                dir_fc,
                                
                                p.val = 0.05,
                                
                                cel_type="",
                                var_column = "tissue_disease",
                                comparison,
                                covariable.dataset = TRUE,
                                hvgs=TRUE) {
 # test_adata = test_adata[, test_adata[[var_column]] %in% var.included]
  
  
  # Select only one cell type of interest
  test_adata = test_adata[, test_adata$Cell_type_scANVI %in% cel_type]
  
  
  counts(test_adata) <- assay(test_adata)
  
  
  
  test_adata = logNormCounts(test_adata)
  
  
  
  # genes where at least some counts are more than 10 count
  
  non_empty_rows <- which(rowSums2(assay(test_adata)) > 10)
  
  test_adata <- test_adata[non_empty_rows,]
  
  ## Feature selection with scran. 5000 hvgs will be selected
  if(hvgs) {
    dec <- modelGeneVar(test_adata)
    
    hvgs <- getTopHVGs(dec, n = 5000)
    
    test_adata <- test_adata[hvgs,]
  }
  
  
  counts(test_adata) <- as.matrix(counts(test_adata))
  
  
  
  assay(test_adata, "X") = as.matrix(assay(test_adata, "X"))
  
  
  
  test_adata$disease_short = factor(test_adata$disease_short) %>% droplevels()
  
  
  
  test_adata$nhoo_test_condition = droplevels(test_adata$nhoo_test_condition)
  
  test_adata[[var_column]] = factor(test_adata[[var_column]]) %>% droplevels()
  
  test_adata$donor_id = factor(test_adata$donor_id) %>% droplevels()
  
  test_adata$dataset_id = factor(test_adata$dataset_id) %>% droplevels()
  
  test_adata$Cell_type_scANVI = str_replace(test_adata$Cell_type_scANVI, "/", "_") %>% factor()
  
  
  
  test_adata$Cell_type_scANVI = factor(test_adata$Cell_type_scANVI) %>% droplevels()
  
  print(table(test_adata[[var_column]], test_adata$nhoo_test_condition))
  
  test_adata$disease_bin = ifelse(test_adata$disease_short=="Healthy", "Healthy", "Disease")
  
  test_adata$nhoo = case_when(test_adata$nhoo_test_condition==dir_fc ~ "TestNhoo",
                               TRUE ~ "Rest") %>% as.factor()
  
  test_adata$disease_test = paste0(test_adata$disease_bin, "_", test_adata$nhoo) %>% as.factor()
  
  print(table(test_adata$disease_test))
  
  if(covariable.dataset){
    my_design = ~ 0 + disease_test + dataset_id
  } else {
    my_design = ~ 0 + disease_test
  }
  
  
  fit1 <- glm_gp(test_adata, design = my_design, size_factors = test_adata$sizeFactor)
  
  
  
  colnames(fit1$Beta)
  
  if(comparison == "selected_nhoo"){
    de_res <- test_de(fit1, contrast = `disease_testDisease_TestNhoo` - (`disease_testDisease_Rest` + `disease_testHealthy_TestNhoo` + `disease_testHealthy_Rest`)/3)
  } else if(comparison == "dis_vs_ctrl"){
    de_res <- test_de(fit1, contrast = (`disease_testDisease_TestNhoo` + `disease_testDisease_Rest`) - (`disease_testHealthy_TestNhoo` + `disease_testHealthy_Rest`))
  } else if(comparison == "between_nhoo"){
    de_res <- test_de(fit1, contrast = (`disease_testHealthy_TestNhoo` + `disease_testDisease_TestNhoo`) - (`disease_testDisease_Rest` + `disease_testHealthy_Rest`))
  } else if(comparison == "enriched_vs_depleted"){
    
    test_adata$disease_test = paste0(test_adata$disease_bin, "_", test_adata$nhoo_test_condition) %>% as.factor()
    my_design = ~ 0 + disease_test
    fit1 <- glm_gp(test_adata, design = my_design, size_factors = test_adata$sizeFactor)
    
    de_res <- test_de(fit1, contrast = `disease_testDisease_Depleted` - `disease_testHealthy_Depleted`)
    de_res2 <- test_de(fit1, contrast = `disease_testDisease_Enriched`- `disease_testHealthy_Enriched`)
   }
  
  
  
  if(comparison != "enriched_vs_depleted"){
    
    if(nrow(de_res) > 0) {
      de_res$dirFC = ifelse(
        de_res$lfc >= lfc & de_res$adj_pval < p.val,
        
        "Up",
        
        ifelse(de_res$lfc <= -lfc & de_res$adj_pval < p.val,
               
               "Down",
               
               "NS")
        
      )
      
      
      de_res = de_res %>% filter(dirFC != "NS") %>% arrange(desc(lfc))
      
    }
    
  } else {
    
    de_res[["dirFC_rest_dis_depleted"]] = ifelse(
      de_res$lfc >= lfc & de_res$adj_pval < p.val, "Up",
      
      ifelse(de_res$lfc <= -lfc & de_res$adj_pval < p.val,
             
             "Down", "NS")
      
    )
    
    de_res[["dirFC_dis_enriched"]] = ifelse(
      de_res2$lfc >= lfc & de_res2$adj_pval < p.val, "Up",
      
      ifelse(de_res2$lfc <= -lfc & de_res2$adj_pval < p.val,
             
             "Down", "NS")
      
    )
    
    de_res[["lfc_dis_enriched"]] = de_res2$lfc
    de_res[["adj_pval_dis_enriched"]] = de_res2$adj_pval
    
  }
  
  
  
  de_res
  
}




