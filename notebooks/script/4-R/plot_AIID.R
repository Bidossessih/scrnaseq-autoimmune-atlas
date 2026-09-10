library(tidyverse)
library(RColorBrewer)
library(ggtext)
library(cowplot)
library(patchwork)
library(ggdendro)

##### Milo plot function ####

plot_milo = function(data, condition, yintercept=0, col.pal,
                     ns_col = "gray70",  enriched_col = "darkgreen", depleted_col = "red", 
                     plot_tile, aspect.ratio, ...){
  
  
  #####
  colourCount = length(unique(data$nhood_annotation))
  
  getPalette = colorRampPalette(brewer.pal(12, "Paired"))
  
  Palette = getPalette(colourCount)
  
  names(Palette) = unique(data$nhood_annotation)
  #print(Palette)
  #Palette = c(Palette, nhoo_col)
  
  if(col.pal==""){
    Palette = Palette
  } else if(col.pal!="" & length(col.pal) < colourCount) {
    Palette = Palette
  } else {
    Palette = col.pal[1:colourCount]
    
    names(Palette) = unique(data$nhood_annotation)
    
  }
  
  
  DE_col = c("Non sig" = ns_col,  "Enriched" = enriched_col, Depleted = depleted_col)
  
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
    #coord_flip() +
    scale_colour_manual(values = DE_col) +
    #scale_alpha_identity() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_void() +
    ggtitle(plot_tile) +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.text = element_blank(),
      plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
      ...
    ) 
  
}


# Function to plot horizontal plot with annotation bar
# y are setted to 0 to get only the annotation

plot_milo_with_annot = function(data, col.pal, condition, legend.position = "none", 
                                aspect.ratio){
  DE_col = c("Non sig" = "gray60",  "Enriched" = "darkgreen", Depleted = "red")
  
  colourCount = length(unique(data$nhood_annotation))
  
  getPalette = colorRampPalette(brewer.pal(12, "Paired"))
  
  Palette = getPalette(colourCount)
  
  names(Palette) = unique(data$nhood_annotation)
  #print(Palette)
  #Palette = c(Palette, nhoo_col)
  
  if(col.pal==""){
    Palette = Palette
  } else if(col.pal!="" & length(col.pal) < colourCount) {
    Palette = Palette
  } else {
    Palette = col.pal[1:colourCount]
    
    names(Palette) = unique(data$nhood_annotation)
    
  }
  
  
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
    scale_y_continuous(expand = c(0,0)) +
    theme_void() +
    theme(
      legend.position = legend.position,
      legend.key.size = unit(2, "cm"),
      legend.text = element_text(size = 25),
      legend.title = element_text(size = 35),
      axis.title = element_blank(),
      axis.text = element_blank()
    ) 
  
}


#' Plotting of horizontal differential abundance 
#' 
#' @param nhoo_prep data frame prepared by 
#' @param rel_widths and rel_heights proportion of the main plot 
#' and the legend. See cowplot
#' 

milo_plot_list = function(nhoo_prep, rel_widths = c(1, 0.25),
                          rel_heights = c(3, 2), col.pal=""){
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
  
  as.dendrogram(clusters_cond)
  
  dhc <- as.dendrogram(clusters_cond)
  # Rectangular lines
  ddata <- dendro_data(dhc, type = "rectangle")
  p <- ggplot(segment(ddata)) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    coord_flip() + 
    scale_y_reverse(expand = c(0.2, 0)) + 
    scale_x_continuous(expand = c(0,0)) +
    theme_void() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
    ) 
  p
  
  
  disease_clustered = clusters_cond$labels[clusters_cond$order]
  
  ### Legend
  
  legend = ggdraw(get_legend(
    plot_milo_with_annot(
      data = nhoo_prep,
      condition = disease_clustered[1],
      legend.position = "right",
      col.pal = col.pal
    )
  ))
  
  #Plotting annotation bar
  # the first disease in the plotlist will be use
  plotlist[["anot"]] = plot_milo_with_annot(data = nhoo_prep,
                                            condition = disease_clustered[1],
                                            legend.position = "none",
                                            col.pal = col.pal)
  
  # annotation and condition plots list
  for (i in 1:length(disease_clustered)) {
    
    dis = disease_clustered[i]
    plotlist[[dis]] = plot_milo(data = nhoo_prep,
                                condition = dis,
                                col.pal = "")
    
  }
  
  #plot grid with cowplot
  
  milo_hori_plot = ggdraw(plot_grid(
    plot_grid(plot_grid(NULL, p, ncol = 1, rel_heights = c(1,length(plotlist)-1)), ncol = 2, nrow = 1, rel_widths = c(0.15,1),
    plot_grid(
      plotlist = plotlist,
      labels = c("", names(plotlist)[-1]),# label may have the order of plot in plotlist
      label_size = 35,
      ncol = 1,
      nrow = length(plotlist),
      align = "hv",
      #label_y = 1,
      hjust = 0, vjust = c(0.5, 1.2, rep(0.5, length(plotlist)-2))
      
    )),
    plot_grid(legend, ncol = 1),
    rel_widths = rel_widths,
    rel_heights = rel_heights
  ))
  
  
  return(list(milo_hori_plot=milo_hori_plot, clusters_cond=clusters_cond))
}


#' Plotting of horizontal differential abundance with dendrogram
#' 
#' @param nhoo_prep data frame prepared by 
#' @param rel_widths and rel_heights proportion of the main plot 
#' and the legend. See cowplot
#' 


milo_horz_plot = function(nhoo_prep,
                          rel_widths = c(0.15, 1, 0.25),
                          col.pal = "",
                          label_size = 35,
                          expand_y_sc = 0.1,
                          anot_bar_height = 0.05,
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
    select(index_cell, dumFC, condition) %>%
    pivot_wider(names_from = condition, values_from = dumFC) %>%
    tibble::column_to_rownames(var = "index_cell")
  clusters_cond <- hclust(dist(t(dummie_data)))
  as.dendrogram(clusters_cond)
  dhc <- as.dendrogram(clusters_cond)
  # Rectangular lines
  ddata <<- dendro_data(dhc, type = "rectangle")
  #ddata = segment(ddata)
  #ddata$x=ddata$x+x1
  #ddata$xend=ddata$xend+x2
  dendroplot <- ggplot(segment(ddata)) +
    geom_segment(aes(
      x = x ,
      y = y,
      xend = xend,
      yend = yend
    )) +
    coord_flip() +
    scale_y_reverse(expand = rep(expand_y_sc, 2)) +
    #scale_x_continuous(expand=c(0,0)) +
    #xlim(c(0, length(clusters_cond$order))) +
    theme_void() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),)
  # Make sure that the dendrogram and the plot have the same order
  disease_clustered = clusters_cond$labels[clusters_cond$order]
  ### Legend
  bar_legend = ggdraw(get_legend(
    plot_milo_with_annot(
      data = nhoo_prep,
      condition = disease_clustered[1],
      legend.position = "right",
      col.pal = col.pal
    )
  ))
  #Plotting annotation bar
  # the first disease in the plotlist will be use
  anot = plot_milo_with_annot(
    data = nhoo_prep,
    condition = disease_clustered[1],
    legend.position = "none",
    col.pal = col.pal
  )
  
  
  # annotation and condition plots list
  for (i in 1:length(disease_clustered)) {
    dis = disease_clustered[i]
    plotlist[[dis]] = plot_milo(
      data = nhoo_prep,
      condition = dis,
      col.pal = "",
      plot_tile = dis,
      ...
    )
  }
  # Here is the annotation bar
  # Two fake plots (NULL) have been used to align the annotation with
  # the horizontal plot.
  anot_bar = (plot_spacer() + anot + plot_spacer() + plot_layout(
    ncol = 3,
    nrow = 1,
    widths = rel_widths
  ))
  # Horizontal plot composed by:
  # - the dendrogram (dendroplot)
  # - Milo plots
  # - the legend
  main_plot = wrap_plots(plotlist, nrow = length(plotlist), ncol = 1)
  
  main_horiz_plot = (dendroplot + main_plot + bar_legend +
                       plot_layout(ncol = 3, widths = rel_widths))
  
  if (anot_bar_height > 0.15) {
    anot_bar_height = 1 / (2 * length(plotlist)) # half of the space occupied by each horizontal plots
  }
  
  horiz_plot = anot_bar / main_horiz_plot +  plot_layout(
    ncol = 1,
    nrow = 2,
    heights = c(anot_bar_height, 1)
  )
  
  #plot(horiz_plot)
  
  blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
    theme_void()
  
  dum_plot = list()
  
  dum_plot[["b1"]] = blankPlot
  
  dum_plot[["anot"]] = anot
  
  dum_plot[["b2"]] = blankPlot
  
  dum_plot[["dendroplot"]] = dendroplot
  
  for (i in 1:length(disease_clustered)) {
    dis = disease_clustered[i]
    dum_plot[[dis]] = plot_milo(
      data = nhoo_prep,
      condition = dis,
      col.pal = "",
      plot_tile = dis,
      ...
    )
  }
  
  dum_plot[["legend"]] = legend
  
  
  #return(list(anot, dendroplot, plotlist, legend))
  return(horiz_plot)
}



plot_milo_simple = function(data,
                            condition,
                            col.pal,
                            ns_col = "gray70",
                            enriched_col = "darkgreen",
                            depleted_col = "red",
                            annot.bar.height = 0.5,
                            condition_tile,
                            annot_size = 25,
                            annot_col = "black",
                            aspect.ratio,
                            ddata_seg,
                            dendrogram,
                            ...) {
  
  
  
  #####
  colourCount = length(unique(data$nhood_annotation))
  
  getPalette = colorRampPalette(brewer.pal(12, "Paired"))
  
  Palette = getPalette(colourCount)
  
  names(Palette) = sort(unique(data$nhood_annotation))
  #Palette = Palette[order(names(Palette))]
  #print(Palette)
  #Palette = c(Palette, nhoo_col)
  
  if(col.pal==""){
    Palette = Palette
  } else if(col.pal!="" & length(col.pal) < colourCount) {
    Palette = Palette
    Palette = Palette[order(names(Palette))]
  } else {
    Palette = col.pal[1:colourCount]
    
    names(Palette) = unique(data$nhood_annotation)
    Palette = Palette[order(names(Palette))]
    
  }
  
  
  DE_col = c("Non sig" = ns_col,  "Enriched" = enriched_col, Depleted = depleted_col)
  
  #Prepare data to include intercept of v_line
  #data = data[data$condition == condition,] 
  
  #vline intercepts
  v_int = unique(data$X_2_max)
  #v_int = v_int[-(length(v_int))]
  yintercept=unique(data$y0)
  
  
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
                      ))
  
  #### Dendro data ####
  
  if(dendrogram){
    dendro = geom_segment(data = ddata_seg, aes(
      y = x ,
      x = -5 * y,
      yend = xend,
      xend = -5 * yend
    ))
  } else{
    dendro = geom_blank()
  }
  
  #### Text annotaion: condition label ####
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
  
  #### Text annotaion: condition label ####
  ymax = max(yintercept) + 0.8 + annot.bar.height
  
  annot_bar = geom_rect(aes(
    ymin = max(yintercept) + 0.8,
    ymax = ymax,
    xmin = c(1, sort(v_int)[1:(length(v_int) - 1)]),
    xmax = c(sort(v_int)),
    fill = names(Palette)
  ))
  
  
  #data %>%
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
    scale_colour_manual(values = DE_col, guide = "none") + #guide = "none" to remove legend
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    annot_text +
    annot_bar  +
    scale_fill_manual('nhood_annotation',
                      values = Palette,
                      guide = guide_legend(override.aes = list(alpha = 1))) +
    dendro+
    ylim(0, ymax) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      #plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
      ...
    ) 
}


milo_horz_plot_simple = function(nhoo_prep,
                                 rel_widths = c(0.15, 1, 0.25),
                                 col.pal = "",
                                 label_size = 35,
                                 annot_size = 25, 
                                 annot_col = "black",
                                 expand_y_sc = 0.1,
                                 anot_bar_height = 0.05,
                                 dendrogram = TRUE,
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
    select(index_cell, dumFC, condition) %>%
    pivot_wider(names_from = condition, values_from = dumFC) %>%
    tibble::column_to_rownames(var = "index_cell")
  clusters_cond <- hclust(dist(t(dummie_data)))
  as.dendrogram(clusters_cond)
  dhc <- as.dendrogram(clusters_cond)
  # Rectangular lines
  ddata <<- dendro_data(dhc, type = "rectangle")
  ddata_seg <<- segment(ddata)
  #ddata$x=ddata$x+x1
  ddata_seg <- segment(ddata)
  #ddata$xend=ddata$xend+x2
  dendroplot <- ggplot(segment(ddata)) +
    geom_segment(aes(
      x = x ,
      y = y,
      xend = xend,
      yend = yend
    )) +
    coord_flip() +
    scale_y_reverse(expand = rep(expand_y_sc, 2)) +
    #scale_x_continuous(expand=c(0,0)) +
    #xlim(c(0, length(clusters_cond$order))) +
    theme_void() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),)
  # Make sure that the dendrogram and the plot have the same order
  disease_clustered = clusters_cond$labels[clusters_cond$order]
  ### Legend
  bar_legend = ggdraw(get_legend(
    plot_milo_with_annot(
      data = nhoo_prep,
      condition = disease_clustered[1],
      legend.position = "right",
      col.pal = col.pal
    )
  ))
  #Plotting annotation bar
  # the first disease in the plotlist will be use
  anot = plot_milo_with_annot(
    data = nhoo_prep,
    condition = disease_clustered[1],
    legend.position = "none",
    col.pal = col.pal
  )
  
  #transform
  norm_factor = nhoo_prep %>% group_by(condition) %>%
    summarise(max_val = 3*max(abs(LogFC))) %>% as.data.frame()
  
  rownames(norm_factor) = as.character(norm_factor$condition)
  norm_factor = norm_factor[rev(disease_clustered),]
  
  
  nhoo_prep$LogFC_transf = nhoo_prep$y0 = NA
  
  for (i in 1:nrow(norm_factor)) {
    ind = which(nhoo_prep$condition == norm_factor$condition[i])
    nhoo_prep$LogFC_transf[ind] = i + (nhoo_prep$LogFC[ind] / norm_factor$max_val[i])
    nhoo_prep$y0[ind] = i
  }
  
  
  # annotation and condition plots list
  ##for (i in 1:length(disease_clustered)) {
  ## dis = disease_clustered[i]
  ## plotlist[[dis]] = plot_milo(
  ##data = nhoo_prep,
  ##condition = dis,
  ##col.pal = "",
  ##plot_tile = dis,
  ## ...
  ##)
  ##}
  
  
  
  # Here is the annotation bar
  # Two fake plots (NULL) have been used to align the annotation with
  # the horizontal plot.
  anot_bar = (plot_spacer() + anot + plot_spacer() + plot_layout(
    ncol = 3,
    nrow = 1,
    widths = rel_widths
  ))
  # Horizontal plot composed by:
  # - the dendrogram (dendroplot)
  # - Milo plots
  # - the legend
  dis_label = length(disease_clustered):1
  
  names(dis_label) = disease_clustered
  #data <<- nhoo_prep 
  
  main_plot = plot_milo_simple(
    data = nhoo_prep,
    col.pal = "",
    annot_size = annot_size, annot_col = annot_col,
    condition_tile = dis_label,
    ddata_seg = ddata_seg,
    dendrogram = dendrogram,
    ...
  )
  #wrap_plots(plotlist, nrow = length(plotlist), ncol = 1)
  
  main_horiz_plot = (main_plot + bar_legend +
                       plot_layout(ncol = 2, widths = rel_widths))
  
  if (anot_bar_height > 0.15) {
    anot_bar_height = 1 / (2 * length(plotlist)) # half of the space occupied by each horizontal plots
  }
  
  horiz_plot = anot_bar / main_horiz_plot +  plot_layout(
    ncol = 1,
    nrow = 2,
    heights = c(anot_bar_height, 1)
  )
  
  #plot(horiz_plot)
  
  
  #return(list(anot, dendroplot, plotlist, legend))
  #return(horiz_plot)
  return(main_plot)
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
                     nhood_annotation="nhood_annotation", n_shared=0,
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
  
  #print(head(nhood_res))
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
    
    #print(head(dummie_data))
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

gsea_plot = function(gsea, n_top = 10, x = "Term", y = "combined_score", log_y = TRUE, bar_plot = TRUE,
                     ylab="log2(Combined score)", xlab = "Hallmark GSEA"){
  
  gsea = na.omit(gsea[1:n_top,])
  
  if(log_y){
    gsea[[y]] = log2(gsea[[y]])
  }
  
  if(bar_plot){
    ggplot(gsea, aes(y = .data[[y]], x = reorder(.data[[x]], .data[[y]]), fill = "red")) +
      geom_bar(stat="identity") + 
      #scale_x_continuous(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      labs(x=xlab, y=ylab) +
      coord_flip() +
      theme_light() + #coord_flip() + 
      theme(#axis_text_x = element_text(hjust=1, angle = 45),
        legend.position = "none",
        legend.text = element_text(size=10),
        panel.grid.minor = element_blank())
  }else{
    ggplot(gsea, aes(y = .data[[y]], x = reorder(.data[[x]], .data[[y]]))) +
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
        axis.title=element_text(size=14),
        axis.text.y = element_text(size = 14)
        #panel.grid.major.y = element_line(colour = "grey60", linetype = "dashed")
      )
    
  }
  
  
}

