library(tidyverse)
library(RColorBrewer)
rm(list = ls())
library(ggtext)
library(cowplot)


##### Milo plot function ####

plot_milo = function(data, condition, yintercept=0){
  colourCount = length(unique(data$nhood_annotation))
  
  getPalette = colorRampPalette(brewer.pal(12, "Paired"))
  
  Palette = getPalette(colourCount)
  
  names(Palette) = unique(data$nhood_annotation)
  
  DE_col = c("Non sig" = "gray70",  "Enriched" = "darkgreen", Depleted = "red")
  
  #Prepare data to include intercept of v_line
  data = data[data$condition == condition,] 
  
  #vline intercepts
  v_int = unique(data$X_2_max)
  #v_int = v_int[-(length(v_int))]
  
  data %>% 
    ggplot() +
    geom_segment(aes(x = X_2, y = 0, xend=X_2, yend=LogFC,
                     colour = direction_FC)) + 
    geom_hline(yintercept = yintercept, linewidth=0.2, color="grey40") +
    geom_vline(xintercept = v_int, linewidth=0.5, color="grey40", linetype =2) + 
    #geom_segment(aes(x = X_3, y = 0, xend=X_3, yend=max(LogFC),
    #                 colour = "grey40"), linewidth = 0.5,
    #             linetype =2) +
    #scale_y_log10() + 
    scale_colour_manual(values = DE_col) +
    #scale_alpha_identity() +
    scale_x_continuous(expand = c(0,0)) +
    theme_void() +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.text = element_blank(),
    ) 
  
}


# Function to plot horizontal plot with annotation bar
# y are setted to 0 to get only the annotation

plot_milo_with_annot = function(data, 
                                condition, legend.position = "none"){
  DE_col = c("Non sig" = "gray60",  "Enriched" = "darkgreen", Depleted = "red")
  
  colourCount = length(unique(data$nhood_annotation))
  
  getPalette = colorRampPalette(brewer.pal(12, "Paired"))
  
  Palette = getPalette(colourCount)
  
  names(Palette) = unique(data$nhood_annotation)
  #print(Palette)
  #Palette = c(Palette, nhoo_col)
  
  #Prepare data to include intercept of v_line
  data = data[data$condition == condition,]
  
  
  
  data %>%  
    ggplot() +
    geom_segment(aes(x = X_2, y = 0, xend=X_2, yend=0), show.legend=FALSE) + 
    #geom_hline(yintercept = yintercept, linewidth=0.2, color="grey40") +
    
    #geom_vline(xintercept = v_int, linewidth=0.5, color="grey40", linetype =2) + 
    geom_ribbon(aes(x = X_2, ymin = 0, ymax = 0.5, 
                    fill = nhood_annotation)) +
    #geom_segment(aes(y=max(LogFC)+4, yend=max(LogFC)+10, xend=X_2, x=X_2, colour = nhood_annotation), 
    #linetype=1, linewidth=60) +
    #geom_segment(aes(x = X_2, y = 0, xend=X_2, yend=LogFC,
    #  colour = direction_FC, alpha = my_alpha)) +
    #scale_y_log10() +
    scale_color_manual(values = DE_col) +
    scale_fill_manual(values = Palette) +
    #scale_alpha_identity() +
    ylim(c(0, 1)) +
    scale_x_continuous(expand = c(0,0)) +
    theme_void() +
    theme(
      legend.position = legend.position,
      legend.key.size = unit(2, "cm"),
      legend.text = element_text(size = 25),
      legend.title = element_text(size = 35),
      axis.title = element_blank(),
      axis.text = element_blank(),
    ) 
  
}


#' Plotting of horizontal differential abundance 
#' 
#' @param nhoo_prep data frame prepared by 
#' @param rel_widths and rel_heights proportion of the main plot 
#' and the legend. See cowplot
#' 

milo_plot_list = function(nhoo_prep, rel_widths = c(1, 0.25),
                          rel_heights = c(3, 2)){
  plotlist = list()
  
  n_nhoo = nhoo_prep$index_cell %>% unique() %>% length()
  n_condition = nhoo_prep$condition %>% unique() %>% length()
  
  #### cluster conditions 
  dummie_data = nhoo_prep %>% mutate(dumFC = case_when(
    direction_FC == "Depleted" ~ -1,
    direction_FC == "Enriched" ~ 1,
    TRUE ~ 0
  )) %>%
    select(index_cell, dumFC, condition) %>%
    pivot_wider(names_from = condition, values_from = dumFC) %>%
    tibble::column_to_rownames(var = "index_cell")
  
  clusters_cond <- hclust(dist(t(dummie_data)))
  
  #plot(clusters_cond)
  
  disease_clustered = clusters_cond$labels[clusters_cond$order]
  
  ### Legend
  
  legend = ggdraw(get_legend(
    plot_milo_with_annot(
      data = nhoo_prep,
      condition = disease_clustered[1],
      legend.position = "right"
    )
  ))
  
  #Plotting annotation bar
  # the first disease in the plotlist will be use
  plotlist[["anot"]] = plot_milo_with_annot(data = nhoo_prep,
                                            condition = disease_clustered[1],
                                            legend.position = "none")
  
  # annotation and condition plots list
  for (i in 1:length(disease_clustered)) {
    
    dis = disease_clustered[i]
    plotlist[[dis]] = plot_milo(data = nhoo_prep,
                                condition = dis)
    
  }
  
  #plot grid with cowplot
  
  milo_hori_plot = ggdraw(plot_grid(
    plot_grid(
      plotlist = plotlist,
      labels = c("", names(plotlist)[-1]),# label may have the order of plot in plotlist
      label_size = 35,
      ncol = 1,
      nrow = length(plotlist),
      align = "hv",
      #label_y = 1,
      hjust = 0, vjust = c(0.5, 1.2, rep(0.5, length(plotlist)-2))
      
    ),
    plot_grid(legend, ncol = 1),
    rel_widths = rel_widths,
    rel_heights = rel_heights
  ))
  
  
  return(milo_hori_plot)
}



#' Prepare neighbourhood abundance data for horizontal plotting
#' 
#' This function take as input the result of neighbourhood abundance 
#' test and returns a data frame formated to be used in the plotting across condition.
#' It returns a data frame
#' @param nhood_res data frame containing the results of neighbourhood (nhood) differential abundance.
#' The data frame must have the following columns label:
#' - logFC log of fold change
#' - nhood_annotation label of the column containing nhood annotation (eg: cell type names)
#' - index_cell column with nhood labels
#' - Sig label of column 
#' - SpatialFDR 
#' @param Sig (valor bool) is the label  of significance column
#' if Sig is not in the data frame, provide the level of significance thres. The default is 0.1
#' @param condition string indicating the label of the tested conditions column. 
#' Default is "disease".
#' @param index_cell string label of nhoo column
#' @param n_shared interger indicating whether to filter out nhood that are not shared
#' between n conditions. The default is 1 for at least one condition. 
#' In this case all nhoods are keeped in the data for visualization.


prep_nhoo = function(nhood_res, condition = "disease", Sig = NULL, thres = 0.1,
                     dis_pref = "all_", index_cell = "index_cell", 
                     nhood_annotation="nhood_annotation", n_shared=1,
                     SpatialFDR = "SpatialFDR", logFC = "logFC"){
  #nhood_res = Myel_tissue #to be removed
  
  # harmonize colnames
  if(is.character(Sig)){
    colnames(nhood_res)[which(colnames(nhood_res) %in% Sig)] = "Sig"
  }
  
  colnames(nhood_res)[which(colnames(nhood_res) %in% index_cell)] = "index_cell"
  
  colnames(nhood_res)[which(colnames(nhood_res) %in% nhood_annotation)] = "nhood_annotation"
  
  colnames(nhood_res)[which(colnames(nhood_res) %in% condition)] = "condition"
  
  
  ### Set the direction of the abundance: depleted, enriched or non significant
  
  if(!is.null(Sig) & is.character(Sig)) {
    
    nhood_res = nhood_res %>% 
      mutate(direction_FC = case_when(logFC < 0 & Sig ~ "Depleted",
                                      logFC > 0 & Sig ~ "Enriched",
                                      TRUE ~ "Non sig"))
  } else if(is.null(Sig)) {
    nhood_res = nhood_res %>% 
      mutate(direction_FC = case_when(logFC < 0 & Sig ~ "Depleted",
                                      logFC > 0 & Sig ~ "Enriched",
                                      TRUE ~ "Non sig"))
  } else {
    stop("Please provide a valide label of significance column")
  }
  
  # remove disease profix if present. optional
  
  nhood_res$condition = str_remove(nhood_res$condition, dis_pref) %>% as.factor()
  
  
  nhood_res$direction_FC = as.factor(nhood_res$direction_FC)
  
  condition = unique(nhood_res$condition) %>% as.character()
  
  nhoo_name = unique(nhood_res$index_cell)
  
  nhoo = unique(nhood_res$nhood_annotation) %>% as.character() %>% sort()
  
  
  for(i in 1:length(nhoo)){
    #subset data of each nhood_annotation
    tab = nhood_res[nhood_res$nhood_annotation==nhoo[i], ]
    
    # Transform enriched nhoo in 1, depleted in -1 and non signicant in 0
    # The magnitude of the FC is not included. The most important 
    # for the comparison between conditions is the 
    # direction of the FC.
    #This way, nhoo with the same FC direction will tend to be clustered together
    
    dummie_data = tab %>% mutate(dumFC = ifelse(
      direction_FC == "Depleted",
      -1,
      ifelse(direction_FC == "Enriched", 1,
             0)
    )) %>%
      select(index_cell, dumFC, condition) %>%
      pivot_wider(names_from = condition, values_from = dumFC) %>%
      tibble::column_to_rownames(var = "index_cell")
    
    # Select nhoo shared by at least n_shared diseases (in the same direction)
    
    #count number of in which each nhood is depleted
    dummie_data$count1 = apply(dummie_data, 1, function(x) length(which(x==-1)))
    #count number of in which each nhood is enriched. the column of count1 has been removed.
    dummie_data$count2 = apply(dummie_data[,-c(ncol(dummie_data))], 1, function(x) length(which(x==1)))
    
    #Filter out nhood that are not shared in at least n_shared condition
    dummie_data = dummie_data %>% 
      mutate(resCount = pmax(count1, count2)) %>% 
      filter(resCount >= n_shared) %>% select(-c(count1,count2, resCount))
    
    #cluster nhood using euclidean distance.
    clusters <- hclust(dist(dummie_data))
    # order nhood based on hclust
    cell_index = clusters$labels[clusters$order]
    
    # Set the order of nhoo
    if(i==1){
      nhoo_oder = data.frame(index_cell=cell_index,
                             X_2 = 1:length(cell_index),
                             X_2_max = max(1:length(cell_index)),
                             X_2_median = median(1:length(cell_index)))
    }else{
      nhoo_oder0 = data.frame(index_cell = cell_index,
                              X_2 = (nrow(nhoo_oder) + 1):(nrow(nhoo_oder) +
                                                             length(cell_index)))
      nhoo_oder0$X_2_max = max(nhoo_oder0$X_2)#Will be used to draw vline
      nhoo_oder0$X_2_median = median(nhoo_oder0$X_2) # Can be used to annotate the annotation bar
      
      nhoo_oder = rbind(nhoo_oder, nhoo_oder0)
    }
    
  }
  
  
  data_ordered = merge(nhoo_oder, nhood_res, by="index_cell",
                       all.x=TRUE)
  
  colnames(data_ordered)[which(colnames(data_ordered) %in% c(logFC, SpatialFDR))] = c("LogFC", "Spa_FDR")
  
  data_ordered$direction_FC = as.factor(data_ordered$direction_FC)
  return(data_ordered)
}

