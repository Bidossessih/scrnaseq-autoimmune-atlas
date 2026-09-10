library(tidyverse)
library(RColorBrewer)
library(ggtext)
library(cowplot)
library(patchwork)
library(ggdendro)


#' Prepare neighbourhood abundance data for horizontal plotting
#'
#' Takes the result of a neighbourhood abundance test and returns a data frame
#' formatted for cross-condition plotting.
#'
#' @param nhood_res Data frame of neighbourhood (nhood) differential abundance results.
#'   Must contain columns for logFC, nhood annotation, nhood index, significance and SpatialFDR.
#' @param condition String; label of the tested-condition column. Default "disease".
#' @param Sig String; label of the (boolean) significance column. If NULL, provide `thres`.
#' @param thres Numeric significance threshold used when `Sig` is NULL. Default 0.1.
#' @param dis_pref String prefix to strip from condition names. Optional.
#' @param index_cell String; label of the nhood index column.
#' @param nhood_annotation String; label of the nhood annotation column (e.g. cell type names).
#' @param n_shared Integer; keep only nhoods shared (same direction) by at least this many
#'   conditions. Default 0 keeps all nhoods.
#' @param SpatialFDR String; label of the SpatialFDR column.
#' @param logFC String; label of the logFC column.
#' @param order_by_kth Logical. If TRUE, order nhoods by `kth_distance` (a column that must
#'   then be present in `nhood_res`); otherwise order by hierarchical clustering.
milo_prep_data = function(nhood_res, condition = "disease", Sig = NULL, thres = 0.1,
                          dis_pref = "all_", index_cell = "index_cell",
                          nhood_annotation = "nhood_annotation", n_shared = 0,
                          SpatialFDR = "SpatialFDR", logFC = "logFC", order_by_kth = TRUE){
  
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
    stop("Please provide a valid label of significance column")
  }
  
  # remove disease prefix if present. optional
  
  nhood_res$condition = str_remove(nhood_res$condition, dis_pref) %>% as.factor()
  
  
  nhood_res$direction_FC = as.factor(nhood_res$direction_FC)
  
  condition = unique(nhood_res$condition) %>% as.character()
  
  nhoo_name = unique(nhood_res$index_cell)
  
  nhoo = unique(nhood_res$nhood_annotation) %>% as.character() %>% sort()
  
  for(i in 1:length(nhoo)){
    # subset data of each nhood_annotation
    tab = nhood_res[nhood_res$nhood_annotation == nhoo[i], ]
    
    # Transform enriched nhoo in 1, depleted in -1 and non significant in 0.
    # The magnitude of the FC is not included; what matters for the comparison
    # between conditions is the direction of the FC, so nhoos with the same FC
    # direction tend to cluster together.
    
    dummie_data = tab %>% mutate(dumFC = ifelse(
      direction_FC == "Depleted",
      -1,
      ifelse(direction_FC == "Enriched", 1,
             0)
    )) %>%
      dplyr::select(index_cell, dumFC, condition) %>%
      pivot_wider(names_from = condition, values_from = dumFC) %>%
      tibble::column_to_rownames(var = "index_cell")
    
    # Select nhoos shared by at least n_shared conditions (in the same direction)
    
    # count number of conditions in which each nhood is depleted
    dummie_data$count1 = apply(dummie_data, 1, function(x) length(which(x == -1)))
    # count number of conditions in which each nhood is enriched (count1 column excluded)
    dummie_data$count2 = apply(dummie_data[, -c(ncol(dummie_data))], 1, function(x) length(which(x == 1)))
    
    # Filter out nhoods that are not shared in at least n_shared conditions
    dummie_data = dummie_data %>%
      mutate(resCount = pmax(count1, count2)) %>%
      filter(resCount >= n_shared) %>% dplyr::select(-c(count1, count2, resCount))
    
    # cluster nhoods using euclidean distance
    clusters <- hclust(dist(dummie_data))
    # order nhoods based on hclust
    cell_index = clusters$labels[clusters$order]
    
    # Set the order of nhoos
    if(i == 1){
      # Order based on kth distance
      if(order_by_kth) {
        kth_data = tab %>% select(index_cell, kth_distance) %>%
          unique() %>% arrange(kth_distance)
        
        cell_index = kth_data$index_cell
      }
      
      nhoo_oder = data.frame(
        index_cell = cell_index,
        X_2 = 1:length(cell_index),
        X_2_max = max(1:length(cell_index)),
        X_2_median = median(1:length(cell_index))
      )
      
    }else{
      # Order based on kth distance
      if(order_by_kth) {
        kth_data = tab %>% select(index_cell, kth_distance) %>%
          unique() %>% arrange(kth_distance)
        
        cell_index = kth_data$index_cell
      }
      
      nhoo_oder0 = data.frame(index_cell = cell_index,
                              X_2 = (nrow(nhoo_oder) + 1):(nrow(nhoo_oder) +
                                                             length(cell_index))
      )
      
      nhoo_oder0$X_2_max = max(nhoo_oder0$X_2)          # Will be used to draw vline
      nhoo_oder0$X_2_median = median(nhoo_oder0$X_2)    # Can be used to annotate the annotation bar
      
      nhoo_oder = rbind(nhoo_oder, nhoo_oder0)
    }
    
  }
  
  
  data_ordered = merge(nhoo_oder, nhood_res, by = "index_cell",
                       all.x = TRUE)
  
  colnames(data_ordered)[which(colnames(data_ordered) %in% c(logFC, SpatialFDR))] = c("LogFC", "Spa_FDR")
  
  data_ordered$direction_FC = as.factor(data_ordered$direction_FC)
  return(data_ordered)
}


#' Build the main Milo plot
#'
#' @param data Data frame. Output of `milo_prep_data`, with `LogFC_transf` and `y0` added.
#' @param col.pal Named vector of colours for cell/nhood labels.
#' @param ns_col Colour of non-significant neighbourhoods.
#' @param enriched_col Colour of enriched neighbourhoods.
#' @param depleted_col Colour of depleted neighbourhoods.
#' @param annot.bar.height Height of the top annotation bar.
#' @param annot_size Size of cell labels.
#' @param fc_annot_size Size of the fold-change legend.
#' @param annot_col Colour of cell labels.
#' @param dendrogram Logical; whether to draw the dendrogram.
#' @param bar.title Title placed on the annotation bar.
#' @param y.axis.label Character; label drawn along the y axis. "" (default) draws nothing.
#' @param y.axis.label.size Size of the y-axis label.
milo_build_main2 = function(data,
                           col.pal,
                           ns_col = "gray70",
                           enriched_col = "darkgreen",
                           depleted_col = "red",
                           annot.bar.height = 0.5,
                           condition_tile,
                           condition_tile_custom,
                           bar.title,
                           left.annot,
                           annot_size = 25,
                           legend.margin.custom = -75,
                           fc_annot_size = 35,
                           legend.spacing = 2,
                           fc_legend.size = 2,
                           annot_col = "black",
                           ddata_seg,
                           dendrogram,
                           dendro.width = dendro.width,
                           y.axis.label = "",
                           y.axis.label.size = 35,
                           ...) {
  
  #####
  colourCount = length(unique(data$nhood_annotation))
  
  getPalette = colorRampPalette(brewer.pal(12, "Paired"))
  
  Palette = getPalette(colourCount)
  
  names(Palette) = sort(unique(data$nhood_annotation))
  
  if(all(col.pal == "")){
    Palette = Palette
  } else if(all(col.pal != "") & length(col.pal) < colourCount) {
    Palette = Palette
    Palette = Palette[order(names(Palette))]
  } else if(length(col.pal[unique(data$nhood_annotation)]) == colourCount) {
    Palette = col.pal[sort(unique(data$nhood_annotation))]
  } else {
    Palette = col.pal[1:colourCount]
    
    names(Palette) = unique(data$nhood_annotation)
    Palette = Palette[order(names(Palette))]
    
  }
  
  
  DE_col = c("Non sig" = ns_col,  Enriched = enriched_col, Depleted = depleted_col)
  
  # vline intercepts
  v_int = unique(data$X_2_max)
  yintercept = unique(data$y0)
  
  
  ### Annotation bar
  
  annot_data = data[data$condition == unique(data$condition)[1], ] %>% arrange(X_2)
  
  ### define plot compartments
  
  #### Milo data ####
  
  milo = geom_segment(data = data,
                      aes(
                        x = X_2,
                        y = y0,
                        xend = X_2,
                        yend = LogFC_transf,
                        colour = direction_FC
                      ), show.legend = TRUE)
  
  #### Dendro data ####
  
  if(dendrogram){
    dendro = geom_segment(data = ddata_seg, aes(
      y = x ,
      x = -5 * y,
      yend = xend,
      xend = -5 * yend,
      colour = "grey40"
    ), linewidth = dendro.width)
  } else{
    dendro = geom_blank()
  }
  
  
  #### Text annotation: condition label ####
  ymax = max(yintercept) + 0.8 + annot.bar.height
  
  annot_bar = geom_rect(aes(
    ymin = max(yintercept) + 0.8,
    ymax = ymax,
    xmin = c(1, sort(v_int)[1:(length(v_int) - 1)]),
    xmax = c(sort(v_int)),
    fill = names(Palette)
  ))
  
  #### Title annotation: general information about the plot ####
  annot.bar.title = annotate(
    "text",
    x = Inf,
    y = ymax + 1,
    label = bar.title,
    size = annot_size,
    colour = "black",
    fontface = "bold",
    hjust = 1,
    vjust = 1
  )
  
  #### Y-axis label ####
  if(!is.null(y.axis.label) && y.axis.label != ""){
    y_axis_annot = annotate(
      "text",
      x = -Inf,
      y = (min(yintercept) + max(yintercept)) / 2,
      label = y.axis.label,
      size = y.axis.label.size,
      colour = "black",
      fontface = "bold",
      angle = 90,
      hjust = 0.5,
      vjust = 1
    )
  } else {
    y_axis_annot = geom_blank()
  }
  
  #### second Text annotation: tissue of each row or another label ####
  
  if(!is.null(left.annot)){
    
    if(!all(names(condition_tile) %in% names(left.annot))){
      print("Condition_title:")
      print(names(condition_tile))
      
      print("left.annot:")
      print(names(left.annot))
      
      stop("Condition_title names don't match the names of left.annot")
    }
    
    print("Names of Condition_title and left.annot match")
    
    
    left.annot = left.annot[names(condition_tile)]
    
    condition_tile2 = condition_tile
    
    if(is.null(condition_tile_custom)) {
      stop("Provide valid condition_tile_custom. condition_tile_custom may be a named vector.")
    }
    names(condition_tile) = as.character(condition_tile_custom[names(condition_tile)])
    
    
    
  } else {
    left.annot = rep("", length(condition_tile))
  }
  
  
  #### Text annotation: condition label (left) ####
  annot_text = annotate(
    "text",
    x = -Inf,
    y = as.numeric(condition_tile) + 0.5,
    label = names(condition_tile),
    size = annot_size,
    colour = annot_col,
    fontface = "bold",
    hjust = 0,
    vjust = 0.5
  )
  
  
  #### Text annotation: left.annot / tissue label (right) ####
  annot_text2 = annotate(
    "text",
    x = Inf,
    y = as.numeric(condition_tile) + 0.5,
    label = left.annot,
    size = annot_size,
    colour = "black",
    fontface = "plain",
    hjust = 1,
    vjust = 0.5
  )
  
  ggplot() +
    milo +
    # Horizontal lines
    geom_segment(aes(
      x = 1,
      y = yintercept,
      xend = Inf,
      yend = yintercept,
      colour = "grey40"
    )) +
    # vertical lines
    geom_segment(aes(
      x = v_int,
      y = 1,
      xend = v_int,
      yend = max(yintercept) + 0.8,
      colour = "grey40"
    ),
    linetype = 2) +
    scale_x_continuous(expand = c(0.05, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_colour_manual(
      values = DE_col,
      guide = guide_legend(
        title = "Fold change direction",
        ncol = 1,
        title.position = "top",
        title.theme = element_text(size = fc_annot_size + 1,  face = "bold"),
        label.theme = element_text(size = fc_annot_size),
        override.aes = list(linewidth = fc_legend.size)
      )
    ) +
    theme_void() +
    annot_text +
    annot_text2 +
    annot_bar  +
    annot.bar.title +
    y_axis_annot +
    scale_fill_manual('nhood_annotation',
                      values = Palette,
                      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    dendro +
    ylim(0, ymax + 1) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      legend.spacing = unit(legend.spacing, "cm"),
      legend.box.margin = margin(t = legend.margin.custom),
      ...
    )
}





#' Build the main Milo plot
#'
#' @param data Data frame. Output of `milo_prep_data`, with `LogFC_transf` and `y0` added.
#' @param col.pal Named vector of colours for cell/nhood labels.
#' @param ns_col Colour of non-significant neighbourhoods.
#' @param enriched_col Colour of enriched neighbourhoods.
#' @param depleted_col Colour of depleted neighbourhoods.
#' @param annot.bar.height Height of the top annotation bar.
#' @param annot_size Size of cell labels.
#' @param fc_annot_size Size of the fold-change legend.
#' @param annot_col Colour of cell labels.
#' @param dendrogram Logical; whether to draw the dendrogram.
#' @param bar.title Title placed on the annotation bar.
#' @param y.axis.label Character; y-axis title drawn in the left margin. "" (default) draws nothing.
#' @param y.axis.label.size Size of the y-axis title.
#' @param margin.left,margin.right,margin.top,margin.bottom Plot margins in points.
#'   The right margin holds the right-side labels; the left margin holds the y-axis title.
milo_build_main = function(data,
                           col.pal,
                           ns_col = "gray70",
                           enriched_col = "darkgreen",
                           depleted_col = "red",
                           annot.bar.height = 0.5,
                           condition_tile,
                           condition_tile_custom,
                           bar.title,
                           left.annot,
                           annot_size = 25,
                           legend.margin.custom = -75,
                           fc_annot_size = 35,
                           legend.spacing = 2,
                           fc_legend.size = 2,
                           annot_col = "black",
                           ddata_seg,
                           dendrogram,
                           dendro.width = dendro.width,
                           y.axis.label = "",
                           y.axis.label.size = 35,
                           margin.left = 60,
                           margin.right = 160,
                           margin.top = 10,
                           margin.bottom = 10,
                           ...) {
  
  #####
  colourCount = length(unique(data$nhood_annotation))
  
  getPalette = colorRampPalette(brewer.pal(12, "Paired"))
  
  Palette = getPalette(colourCount)
  
  names(Palette) = sort(unique(data$nhood_annotation))
  
  if(all(col.pal == "")){
    Palette = Palette
  } else if(all(col.pal != "") & length(col.pal) < colourCount) {
    Palette = Palette
    Palette = Palette[order(names(Palette))]
  } else if(length(col.pal[unique(data$nhood_annotation)]) == colourCount) {
    Palette = col.pal[sort(unique(data$nhood_annotation))]
  } else {
    Palette = col.pal[1:colourCount]
    
    names(Palette) = unique(data$nhood_annotation)
    Palette = Palette[order(names(Palette))]
    
  }
  
  
  DE_col = c("Non sig" = ns_col,  Enriched = enriched_col, Depleted = depleted_col)
  
  # vline intercepts
  v_int = unique(data$X_2_max)
  yintercept = unique(data$y0)
  
  
  ### Annotation bar
  
  annot_data = data[data$condition == unique(data$condition)[1], ] %>% arrange(X_2)
  
  ### define plot compartments
  
  #### Milo data ####
  
  milo = geom_segment(data = data,
                      aes(
                        x = X_2,
                        y = y0,
                        xend = X_2,
                        yend = LogFC_transf,
                        colour = direction_FC
                      ), show.legend = TRUE)
  
  #### Dendro data ####
  
  if(dendrogram){
    dendro = geom_segment(data = ddata_seg, aes(
      y = x ,
      x = -5 * y,
      yend = xend,
      xend = -5 * yend,
      colour = "grey40"
    ), linewidth = dendro.width)
  } else{
    dendro = geom_blank()
  }
  
  
  #### Text annotation: condition label ####
  ymax = max(yintercept) + 0.8 + annot.bar.height
  
  annot_bar = geom_rect(aes(
    ymin = max(yintercept) + 0.8,
    ymax = ymax,
    xmin = c(1, sort(v_int)[1:(length(v_int) - 1)]),
    xmax = c(sort(v_int)),
    fill = names(Palette)
  ))
  
  #### Title annotation: general information about the plot ####
  annot.bar.title = annotate(
    "text",
    x = Inf,
    y = ymax + 1,
    label = bar.title,
    size = annot_size,
    colour = "black",
    fontface = "bold",
    hjust = 1,
    vjust = 1
  )
  
  #### second Text annotation: tissue of each row or another label ####
  
  if(!is.null(left.annot)){
    
    if(!all(names(condition_tile) %in% names(left.annot))){
      print("Condition_title:")
      print(names(condition_tile))
      
      print("left.annot:")
      print(names(left.annot))
      
      stop("Condition_title names don't match the names of left.annot")
    }
    
    print("Names of Condition_title and left.annot match")
    
    
    left.annot = left.annot[names(condition_tile)]
    
    condition_tile2 = condition_tile
    
    if(is.null(condition_tile_custom)) {
      stop("Provide valid condition_tile_custom. condition_tile_custom may be a named vector.")
    }
    names(condition_tile) = as.character(condition_tile_custom[names(condition_tile)])
    
    
    
  } else {
    left.annot = rep("", length(condition_tile))
  }
  
  
  #### Text annotation: condition label (left, in the dendrogram lane) ####
  annot_text = annotate(
    "text",
    x = 0,
    y = as.numeric(condition_tile) + 0.5,
    label = names(condition_tile),
    size = annot_size,
    colour = annot_col,
    fontface = "bold",
    hjust = 0,
    vjust = 0.5
  )
  
  
  #### Text annotation: left.annot / tissue label (right, in the right margin) ####
  annot_text2 = annotate(
    "text",
    x = Inf,
    y = as.numeric(condition_tile) + 0.5,
    label = left.annot,
    size = annot_size,
    colour = "black",
    fontface = "plain",
    hjust = 1,
    vjust = 0.5
  )
  
  
  #### Y-axis title: drawn in the left margin (no in-panel placement) ####
  draw_y_axis = !is.null(y.axis.label) && y.axis.label != ""
  y_axis_title = if(draw_y_axis) {
    element_text(size = y.axis.label.size, face = "bold", angle = 90,
                 margin = margin(r = 5))
  } else {
    element_blank()
  }
  
  
  ggplot() +
    milo +
    # Horizontal lines
    geom_segment(aes(
      x = 1,
      y = yintercept,
      xend = Inf,
      yend = yintercept,
      colour = "grey40"
    )) +
    # vertical lines
    geom_segment(aes(
      x = v_int,
      y = 1,
      xend = v_int,
      yend = max(yintercept) + 0.8,
      colour = "grey40"
    ),
    linetype = 2) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_colour_manual(
      values = DE_col,
      guide = guide_legend(
        title = "Fold change direction",
        ncol = 1,
        title.position = "top",
        title.theme = element_text(size = fc_annot_size + 1,  face = "bold"),
        label.theme = element_text(size = fc_annot_size),
        override.aes = list(linewidth = fc_legend.size)
      )
    ) +
    coord_cartesian(clip = "off") +
    theme_void() +
    labs(y = if(draw_y_axis) y.axis.label else NULL) +
    annot_text +
    annot_text2 +
    annot_bar  +
    annot.bar.title +
    scale_fill_manual('nhood_annotation',
                      values = Palette,
                      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    dendro +
    ylim(0, ymax + 1) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = y_axis_title,
      axis.text = element_blank(),
      legend.spacing = unit(legend.spacing, "cm"),
      legend.box.margin = margin(t = legend.margin.custom),
      plot.margin = margin(t = margin.top, r = margin.right,
                           b = margin.bottom, l = margin.left, unit = "pt"),
      ...
    )
}
#' Return the horizontal Milo plot
#'
#' @param left.annot Left label of each row (level). Named vector with plot label names
#'   (`condition_tile`) as names, e.g. `c(disease1 = "tissue1", disease2 = "tissue2")`.
#' @param y.axis.label Character; label drawn along the y axis. Passed to `milo_build_main`.
#' @param y.axis.label.size Size of the y-axis label.
milo_plot_horizontal = function(nhoo_prep,
                                rel_widths = c(0.15, 1, 0.25),
                                col.pal = "",
                                label_size = 35,
                                annot_size = 25,
                                fc_annot_size = 35,
                                fc_legend.size,
                                legend.margin.custom = -75,
                                annot_col = "black",
                                left.annot = NULL,
                                expand_y_sc = 0.1,
                                legend.spacing,
                                anot_bar_height = 0.05,
                                dendrogram = TRUE,
                                dendro.width = 1,
                                condition_tile_custom = NULL,
                                bar.title = "",
                                y.axis.label = "",
                                y.axis.label.size = 35,
                                ...) {
  plotlist = list()
  n_nhoo = nhoo_prep$index_cell %>% unique() %>% length()
  n_condition = nhoo_prep$condition %>% unique() %>% length()
  
  #### cluster conditions
  dummie_data = nhoo_prep %>% mutate(dumFC = case_when(
    direction_FC == "Depleted" ~ -1,
    direction_FC == "Enriched" ~ 1,
    TRUE ~ 0
  )) %>%
    dplyr::select(index_cell, dumFC, condition) %>%
    pivot_wider(names_from = condition, values_from = dumFC) %>%
    tibble::column_to_rownames(var = "index_cell")
  clusters_cond <- hclust(dist(t(dummie_data)))
  dhc <- as.dendrogram(clusters_cond)
  
  # Rectangular lines
  ddata <- dendro_data(dhc, type = "rectangle")
  ddata_seg <- segment(ddata)
  
  # Make sure that the dendrogram and the plot have the same order
  disease_clustered = clusters_cond$labels[clusters_cond$order]
  
  
  # transform
  norm_factor = nhoo_prep %>% group_by(condition) %>%
    summarise(max_val = 3 * max(abs(LogFC))) %>% as.data.frame()
  
  rownames(norm_factor) = as.character(norm_factor$condition)
  norm_factor = norm_factor[rev(disease_clustered), ]
  
  
  nhoo_prep$LogFC_transf = nhoo_prep$y0 = NA
  
  for (i in 1:nrow(norm_factor)) {
    ind = which(nhoo_prep$condition == norm_factor$condition[i])
    nhoo_prep$LogFC_transf[ind] = i + (nhoo_prep$LogFC[ind] / norm_factor$max_val[i])
    nhoo_prep$y0[ind] = i
  }
  
  
  dis_label = length(disease_clustered):1
  
  names(dis_label) = disease_clustered
  
  main_plot = milo_build_main(
    data = nhoo_prep,
    col.pal = col.pal,
    fc_annot_size = fc_annot_size,
    annot_size = annot_size,
    annot_col = annot_col,
    fc_legend.size = fc_legend.size,
    legend.margin.custom = legend.margin.custom,
    condition_tile = dis_label,
    legend.spacing = legend.spacing,
    ddata_seg = ddata_seg,
    dendrogram = dendrogram,
    dendro.width = dendro.width,
    bar.title = bar.title,
    left.annot = left.annot,
    condition_tile_custom = condition_tile_custom,
    y.axis.label = y.axis.label,
    y.axis.label.size = y.axis.label.size,
    ...
  )
  
  return(main_plot)
}