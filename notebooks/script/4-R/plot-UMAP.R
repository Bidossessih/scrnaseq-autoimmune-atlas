#' Create UMAP plots side by side.
#' The left plot is colored by marker.
#' The right plot is colored by cluster.
#' @param dat A dataframe with columns T1, T2, marker, cluster


#source('C:/Wilfried/scRNAseq_AR/Reumato_single-cell_rnaseq/Dados/RNA_sequencing_result/Data_R/meta-colors.R')

#function to print disease


plot_umap_disease <- function(dat, umap_x = "UMAP_1", umap_y = "UMAP_2", title = NULL, myplot) {
  n_nonzero  <- sum(dat$marker > 0)
  # tsne_title <- bquote("tSNE of PCA on Log"[2]~"(CPM + 1)")
  point_size <- 1.0
  if (nrow(dat) < 5000) {
    point_size <- 2.0
  }
  # fill_values <- quantile_breaks(dat$marker, n = 9)
  # fill_values <- fill_values / max(fill_values)
  fill_palette <- RColorBrewer::brewer.pal(9, "Greens")
  theme_umap_1 <- theme_bw(base_size = 25) + theme(
    legend.position = "bottom",
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title      = element_text(size = 30),
    legend.text    = element_text(size = 18)
  )
 
  # Umap with clusters
  
  plot <- ggplot() +
    geom_point(data = dat,
      #data    = subset(dat, Erosion == "Abscence"),
      mapping = aes_string(x = umap_x, y = umap_y, fill = "Erosion"),
      size    = point_size,
      shape   = 21,
      
    ) +
    guides(
      fill = guide_legend(
        nrow = 1,
        override.aes = list(size = 3)
      ),
      alpha = "none"
    ) + scale_fill_manual(values = c("#e0e0e0", "#e41a1c"), name = "Erosion") + #, "#1b9e77"
    #scale_fill_brewer(palette = "Set3", name = "Erosion") +
    #labs(x = NULL, y = NULL, title = substitute(italic(x), list(x = title))) +
    
    theme_umap_1 + labs(x = NULL, y = NULL, title = "")
  
  
  
  plot #+ plot_annotation(
    #caption = bottom_text,
    #theme = theme(
     # plot.caption = element_text(size = 25, hjust = 0.5)
   # )
  #)
}



# Function to plot cluster
plot_umap_cluster <- function(dat, umap_x = "UMAP_1", umap_y = "UMAP_2", title = NULL, group = NULL, shape = NULL) {
  stopifnot("Please define group" = !is.null(group) )
  
  stopifnot("Shape must be numeric" = is.numeric(shape) )
  
  n_nonzero  <- sum(dat$marker > 0)
  # tsne_title <- bquote("tSNE of PCA on Log"[2]~"(CPM + 1)")
  point_size <- 1.0
  if (nrow(dat) < 5000) {
    point_size <- 2.0
  }
  # fill_values <- quantile_breaks(dat$marker, n = 9)
  # fill_values <- fill_values / max(fill_values)
  fill_palette <- RColorBrewer::brewer.pal(9, "Greens")
 
  theme_umap_2 <- theme_bw(base_size = 25) + theme(
    legend.position = "bottom",
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title      = element_text(size = 30),
    legend.text    = element_text(size = 18)
  )
  
  
  
  # Make a plot showing the clustering results.
  dat$cluster <- factor(dat$cluster)
  
  plot <- ggplot() +
    geom_point(
      data    = dat,
      # mapping = aes(x = T1, y = T2, fill = cluster),
      mapping = aes_string(x = umap_x, y = umap_y, fill = group),
      size    = point_size,
      shape   = shape,
      stroke  = 0.1
    ) +
    scale_fill_brewer(palette = "Paired", name = group) +
    #scale_fill_manual(values = meta_colors$cluster, name = "Cluster") +
    guides(
      fill = guide_legend(
        nrow = 6,
        override.aes = list(size = 3)
      )
    ) +
    labs(x = NULL, y = NULL, title = title) +
    
    theme_umap_2
  
  
  plot #+ plot_annotation(
  #caption = bottom_text,
  #theme = theme(
  # plot.caption = element_text(size = 25, hjust = 0.5)
  # )
  #)
}
#library(ggplot2)

#plot_umap(umap_ggplot)


# marker <- "THY1"
# gene_ix <- which(gene_symbols == marker)
# meta$marker <- as.numeric(log2cpm[gene_ix,])
# cell_ix <- seq(nrow(meta))
# tsne_x <- "T1_all"
# tsne_y <- "T2_all"
# plot_tsne(meta[cell_ix,], tsne_x, tsne_y, title = marker)
