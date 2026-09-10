#' Function to get association data from OpenTarget
#' 
#' @param diseaseIfo string indicating the IFO code of the disease of interest
#' @param index pagination is set to 0. 
#' @param threshold float indicating the the threshold of genetic association score. 
#' Genes with L2G scores higher than 0.5 are expected to be causal for the respective
#' trait association in 50% of cases.
#' @import httr


disease2gwas = function(diseaseIfo, index=0, threshold = 0.5) {
  # Build query string to get general information about AR and genetic constraint and tractability assessments
  query_string = "query disease($efoId: String!, $index: Int!, $size: Int!) {
   disease(efoId: $efoId) {
    id
    name
    associatedTargets(page: {index: $index, size: $size}) {
      count
      rows {
        target {
          id
          approvedSymbol
          
        }
        
        score
        
        datatypeScores {
          componentId: id
          score
        }
        
      }
    }
    
   }
}"
  
  
  # Set base URL of GraphQL API endpoint
  base_url <- "https://api.platform.opentargets.org/api/v4/graphql"
  
  
  # Check row count
  # get only one line
  v_dum = list("efoId" = diseaseIfo, "index" = 0, "size" = 1)
  
  r_dum <- httr::POST(url = base_url,
                   body = list(query = query_string, variables = v_dum),
                   encode = 'json')
  #get the total row length. This count will be used in the query to retrieve all the rows
  
  count = httr::content(r_dum)$data$disease$associatedTargets$count
  
  # Set variables object of arguments to be passed to endpoint
  
  variables <- list("efoId" = diseaseIfo, "index" = index, "size" = count)
  # Construct POST request body object with query string and variables
  post_body <- list(query = query_string, variables = variables)
  
  
  # Perform POST request
  
  response <- httr::POST(url = base_url,
            body = post_body,
            encode = 'json')
  
  
  
  op_t = response$content %>% rawToChar() %>% jsonlite::fromJSON()
  
  op_t2 = op_t$data$disease$associatedTargets$rows$datatypeScores
  
 # op_t2 = op_t[["disease"]][["associatedTargets"]][["rows"]]
  
  data = op_t$data$disease$associatedTargets$rows$target
  
  data$over.all.score = op_t$data$disease$associatedTargets$rows$score
  
  data$genetic_association = NA
  
  
  for (i in 1:nrow(data)) {
    tab = op_t$data$disease$associatedTargets$rows$datatypeScores[[i]]
    
    if("genetic_association" %in% tab$componentId) {
      data$genetic_association[i] = round(tab$score[tab$componentId == "genetic_association"], 2)
    }
  }
    
  data = data[(!is.na(data$genetic_association) & data$genetic_association >= threshold),]
  rownames(data) = NULL
  return(data)
  
}



#' Function to expand network based on a vector of protein of interest (seed)
#' 
#' @param graph reference interactome from which the network will be expanded
#' @param poi vector of protein of interest used as seed
#' @param exclude genes to be excluded. Generally down-regulated gene or non expressed genes
#' @param include genes to be included. This param is used to remove non-expressed intermediate genes


xpandhelper = function(graph, poi, exclude=NULL, include = NULL, returnIntProt=TRUE){
  
  interProt = c()
  
  for (p in poi) {
    
    # incoming
    inc = neighbors(graph, v = p, mode = "in") %>% as_ids()
    
    # outgoing
    outg = neighbors(graph, v = p, mode = "out") %>% as_ids()
    
    # filter genes
    if(!is.null(exclude)){
      
      outg = outg[!(outg %in% exclude)]
      
      inc = inc[!(inc %in% exclude)]
      
    } 
    
    if(!is.null(include)){
      outg = outg[outg %in% include]
      
      inc = inc[inc %in% include]
    }
    
    
    tab = data.frame(source = c(inc, rep(p, length(outg))),
                     target = c(rep(p, length(inc)), outg))
    
    if(!exists("data_e")){
      data_e = tab
    } else {
      data_e = rbind(data_e, tab)
    }
    
    interProt = c(interProt, inc, outg) %>% unique()
  }
  
  #return also vector of intermediate proteins
  if(returnIntProt){
    response = list(data_e, interProt)
  }else{
    response = data_e
  }
  
  return(response)
  
}


poi2netXpansion = function(graph = OPI_g, poi, exclude = NULL, intermediate=FALSE, include = NULL){
  
  # check if all poi are present in the graph
  p_int = intersect(V(graph) %>% as_ids(), poi)
  
  if(length(poi)!=length(p_int)){
    p_dif = setdiff(poi, V(graph) %>% as_ids())
    
    cat(paste(str_glue("{length(p_dif)} are not present in the network"), "\n"))
    
    cat(p_dif)
    print("#")
    
    poi = p_int
  }
  
  # main network expansion
  dataXpanded = xpandhelper(graph=graph, poi=poi, exclude=exclude, returnIntProt = TRUE, include = include)[[1]]
  print("dim of dataXpanded")
  
  cat(dim(dataXpanded))
  print("#")
  print("Number of vertices:")
  
  print(as.character(length(unique(c(dataXpanded$source, dataXpanded$target)))))
  
  intProt = xpandhelper(graph=graph, poi=poi, exclude=exclude, returnIntProt = TRUE, include = include)[[2]]
  if(!is.null(include)){
    intProt = c(intProt, include)
  }
  # Check link between intermediate proteins
  if(intermediate) {
    intProtXpanded = xpandhelper(
      graph = graph,
      poi = intProt,
      exclude = NULL,
      include = intProt,
      returnIntProt = FALSE
    )
    
    # Combine table
    print("dim of intProtXpanded")
    
    cat(dim(intProtXpanded))
    print("#")
    print("Number of intermediate protein:")
    print(length(intProt))
    print("Number of vertices:")
    print(as.character(length(unique(c(intProtXpanded$source, intProtXpanded$target)))))
    
    data_e = rbind(dataXpanded, intProtXpanded)
  } else{
    data_e = dataXpanded
  }
  return(data_e)
}



#' Function to annotate vertex
#' 
#' @param expanded_g data.frame, expanded network
#' @param poi.meta data.frame, metadata of protein of interest. the first column contains gene names

annotateVertex = function(expanded_g, poi.meta){
  tab = data.frame(unique(c(expanded_g$source,expanded_g$target)))
  colnames(tab) = colnames(poi.meta)[1]
  tab = unique(tab)
  anot = merge(tab, poi.meta, by=colnames(tab), all.x = TRUE)
  
  anot
}


#' transform FC to score ranging from 0 to 1 using min/max normalization
#' nlogfc (y-min(X))/(max(X)-min(X))
#' y is each value of logFC

nlogFC = function(df=int_macr_RA, col="meanLfc", ratio=TRUE){
  minx = min(df[,col])
  
  max_x = max(df[,col])
  
  if(ratio) {
    df = df %>% mutate(lfc_score = (.data[[col]]) / max_x)
  }else{
    df = df %>% mutate(lfc_score = (.data[[col]] - minx) / (max_x - minx))
  }
}



#' Function to construct network based on single disease information
#' 
#' @param graph reference interactome from which the network will be expanded
#' @param poi.meta data.frame containing metadata of protein of interest. The first column contains the official names of genes
#' @param exclude genes to be excluded. Generally down-regulated gene or non expressed genes
#' @param intermediate bool indicating whether edge between intermediate proteins must be included
#' @param walk.step The length of the random walks to perform


disease2net = function(graph = OPI_g, poi.meta, walk.step = 6, exclude = NULL, include = NULL, intermediate = FALSE, min.gene.cluster=5) {
  #poi==protein of interest
  poi = poi.meta[,1] %>% as.factor()
  
  dis_data_exp = poi2netXpansion(graph = graph, poi=poi, include = include, exclude = exclude, intermediate = intermediate)
  
  dis_v_meta = annotateVertex(expanded_g = dis_data_exp, poi.meta = poi.meta)
  
  dis_v_meta$Origin = ifelse(is.na(dis_v_meta$Origin), "Intermediate", dis_v_meta$Origin)
  
  #  set score of intermediate proteins to 0
  dis_v_meta$Score = ifelse(is.na(dis_v_meta$Score), 0, dis_v_meta$Score)
  
  # Set the edge weight to sum of their vertices' score
  dis_data_exp$weight = NA
  
  for (i in 1:nrow(dis_data_exp)) {
    w = sum(c(dis_v_meta$Score[dis_v_meta$name %in% dis_data_exp$source[i]],
              dis_v_meta$Score[dis_v_meta$name %in% dis_data_exp$target[i]]),
            na.rm = TRUE)
    
    dis_data_exp$weight[i] = w
  }
  
  net = graph_from_data_frame(d=dis_data_exp,vertices=dis_v_meta,directed=T)
  
  page.rank.net=page_rank(net, personalized = as.numeric(V(net)$Score), weights = as.numeric(E(net)$weight))
  
  V(net)$page.rank = page.rank.net$vector
  
  degree=igraph::degree(net)
  V(net)$degree = degree
  
  net.del=igraph::simplify(net,
                           remove.loops = T,
                           remove.multiple = T ,
                           edge.attr.comb = c(weight="max","ignore"))
  
  
  # Remove nodes connected to only one node
  net.del = delete_vertices(net.del, which(as.numeric(V(net.del)$degree) <= 1))
  
  cwt.net=cluster_walktrap(net.del, 
                           weights = E(net.del)$weight, 
                           steps = walk.step,
                           merges = F, 
                           modularity = TRUE, 
                           membership = TRUE)
  
  
  V(net.del)$cwt = cwt.net$membership %>% as.factor()
  
  
  print(table(cwt.net$membership))
  # Remove clusters with less than min.gene.cluster
  #cwtTab.net = data.frame(cwt.net$membership %>% table())
  #colnames(cwtTab.net)[1] = "cluster"
  
  cwtnames = table(cwt.net$membership)[table(cwt.net$membership) < min.gene.cluster] %>% names() %>% as.numeric()
  
  print(cwtnames)
  
  #vtx = V(net.del)[as.numeric(V(net.del)$cwt) %in% cwtnames]
  
  net.del = delete_vertices(net.del, V(net.del)$cwt %in% cwtnames)
  
  # Then remove nodes that are not connected and cluster
  degree=igraph::degree(net.del)
  V(net.del)$degree = degree
  
 # net.del = delete_vertices(net.del, V(net.del)[as.numeric(V(net.del)$degree) == 0])
  
  
  # Reset cluster levels to get consecutive numbers after removing some clusters
  x = V(net.del)$cwt
  x = droplevels(x)
  
  levels(x) = 1:length(levels(x))
  
  V(net.del)$cwt = x
  
  
  print("Clusters summary:")
  
  print(table(V(net.del)$cwt))
  
  gplot = ggraph(net.del, layout = "kk") +
    geom_edge_link0(width=0.2, colour = "gray90") +
    #scale_color_brewer(palette = "Paired") +
    geom_edge_fan(aes(alpha = after_stat(index)), show.legend = FALSE) +
    geom_node_point(aes(color = cwt, size = page.rank)) +
    theme_graph()
  
  plot(gplot)
  
  net.del
  
}


#' Transform graph to list based on community membership
#' 


community2list = function(graph){
  comt = V(graph)$cwt %>% levels() %>% as.numeric()
  net_list = list()
  
  for (i in comt) {
    net_list[[i]] = V(graph)[as.numeric(V(graph)$cwt) == i] %>% as_ids()
  }
  
  net_list
}

#' Jaccard
#' 

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


net2similarity = function(g1, g2){
  tab = matrix(nrow = length(g1), ncol = length(g2)) %>% as.data.frame()
  
  overlapTab = matrix(nrow = length(g1), ncol = length(g2)) %>% as.data.frame()
  
  l1 = g1 %>% unlist()
  
  l2 = g2 %>% unlist()
  
  totalgene1 = length(l1)
  totalgene2 = length(l2)
  
  
  rownames(tab) = paste0("disease1_", 1:length(g1))
  
  colnames(tab) = paste0("disease2_", 1:length(g2))
  
  for (i in 1:length(g1)) {
    
    for (j in 1:length(g2)) {
      siml = jaccard(g1[[i]], g2[[j]])
      
      tab[i,j] = siml
      # mean of overlap
      overlapTab[i,j] = (length(intersect(g1[[i]], g2[[j]]))/totalgene1 + length(intersect(g1[[i]], g2[[j]]))/totalgene2)/2
    }
    
  }
  
  overlapTab = tibble::rownames_to_column(overlapTab, var = "disease1") %>% 
    pivot_longer(!disease1, names_to = "disease2", values_to = "percent_overlap")
  
  tab = tibble::rownames_to_column(tab, var = "disease1") %>% 
    pivot_longer(!disease1, names_to = "disease2", values_to = "jaccardSim") 
  
  tab$percent_overlap = overlapTab$percent_overlap
  
  tab = filter(tab, jaccardSim > 0)
  
  # Rank by jaccardSim and percent_overlap
  # Then calculate the rank product (geometric mean of ranks)
  
  tab = tab %>% arrange(desc(jaccardSim)) %>%
    mutate(Rank1=1:nrow(tab)) %>% arrange(desc(percent_overlap)) %>% 
    mutate(Rank2=1:nrow(tab))
  # Rank product
  tab$Rank = (tab$Rank1*tab$Rank2)**0.5
  
  tab = tab %>% arrange(Rank) %>% mutate(Rank=1:nrow(tab))
  
  return(tab)
}


#' Select expressed genes in at least a fraction of sample
#' 
#' @param scedata pseudobulk data in sce format
#' @param min.count minimum count per pseudobulk. Default 10 counts
#' @param sample.size.percent percentage of sample
#' 

psdobkl2expressed = function(scedata, min.count=10, sample.size.percent = 0.1){
  
  scedata$Cell_type_scANVI = gsub(" ", "_", scedata$Cell_type_scANVI) %>% factor()
  
  scedata$Cell_type_scANVI = gsub("[:+:]", "_plus", scedata$Cell_type_scANVI) %>% factor()
  
  scedata$Cell_type_scANVI = gsub("-", "_", scedata$Cell_type_scANVI) %>% factor()
  
  scedata$donor_id = factor(gsub(" ", "_", scedata$donor_id))
  
  scedata$donor_id = factor(gsub("-", "_", scedata$donor_id))
  
  scedata$disease_short = droplevels(scedata$disease_short)
  
  scedata$tissue = gsub( " ", "_", scedata$tissue)
  
  scedata$tissue_disease = paste0(scedata$tissue, "_", scedata$disease_short) %>% factor()
  
  counts(scedata) <- assay(scedata)
  
  scedata = scedata[, scedata$disease_short != "Healthy"]
  # genes where at least some counts are more than min.count count
  n.sample = round(nrow(colData(scedata))*sample.size.percent, 0)
  non_empty_rows <- rowSums(assay(scedata) > min.count) >= n.sample
  #non_empty_rows <- which(rowSums2(assay(scedata)) > min.count)
  scedata <- scedata[non_empty_rows, ]
  
  gene = dimnames(scedata)[[1]] %>% unlist()
  
  gene
  
}


#' Function to get association data from OpenTarget
#' 
#' @param diseaseIfo string indicating the IFO code of the disease of interest
#' @param index pagination is set to 0. 
#' @param threshold float indicating the the threshold of genetic association score. 
#' Genes with L2G scores higher than 0.5 are expected to be causal for the respective
#' trait association in 50% of cases.
#' @import httr


disease2gwas = function(diseaseIfo, index=0, threshold = 0.5) {
  # Build query string to get general information about AR and genetic constraint and tractability assessments
  query_string = "query disease($efoId: String!, $index: Int!, $size: Int!) {
   disease(efoId: $efoId) {
    id
    name
    associatedTargets(page: {index: $index, size: $size}) {
      count
      rows {
        target {
          id
          approvedSymbol
          
        }
        
        score
        
        datatypeScores {
          componentId: id
          score
        }
        
      }
    }
    
   }
}"
  
  
  # Set base URL of GraphQL API endpoint
  base_url <- "https://api.platform.opentargets.org/api/v4/graphql"
  
  
  # Check row count
  # get only one line
  v_dum = list("efoId" = diseaseIfo, "index" = 0, "size" = 1)
  
  r_dum <- httr::POST(url = base_url,
                   body = list(query = query_string, variables = v_dum),
                   encode = 'json')
  #get the total row length. This count will be used in the query to retrieve all the rows
  
  count = httr::content(r_dum)$data$disease$associatedTargets$count
  
  # Set variables object of arguments to be passed to endpoint
  
  variables <- list("efoId" = diseaseIfo, "index" = index, "size" = count)
  # Construct POST request body object with query string and variables
  post_body <- list(query = query_string, variables = variables)
  
  
  # Perform POST request
  
  response <- httr::POST(url = base_url,
            body = post_body,
            encode = 'json')
  
  
  
  op_t = response$content %>% rawToChar() %>% jsonlite::fromJSON()
  
  op_t2 = op_t$data$disease$associatedTargets$rows$datatypeScores
  
 # op_t2 = op_t[["disease"]][["associatedTargets"]][["rows"]]
  
  data = op_t$data$disease$associatedTargets$rows$target
  
  data$over.all.score = op_t$data$disease$associatedTargets$rows$score
  
  data$genetic_association = NA
  
  
  for (i in 1:nrow(data)) {
    tab = op_t$data$disease$associatedTargets$rows$datatypeScores[[i]]
    
    if("genetic_association" %in% tab$componentId) {
      data$genetic_association[i] = round(tab$score[tab$componentId == "genetic_association"], 2)
    }
  }
    
  data = data[(!is.na(data$genetic_association) & data$genetic_association >= threshold),]
  rownames(data) = NULL
  return(data)
  
}



#' Function to expand network based on a vector of protein of interest (seed)
#' 
#' @param graph reference interactome from which the network will be expanded
#' @param poi vector of protein of interest used as seed
#' @param exclude genes to be excluded. Generally down-regulated gene or non expressed genes
#' @param include genes to be included. This param is used to remove non-expressed intermediate genes


xpandhelper = function(graph, poi, exclude=NULL, include = NULL, returnIntProt=TRUE){
  
  interProt = c()
  
  for (p in poi) {
    
    # incoming
    inc = neighbors(graph, v = p, mode = "in") %>% as_ids()
    
    # outgoing
    outg = neighbors(graph, v = p, mode = "out") %>% as_ids()
    
    # filter genes
    if(!is.null(exclude)){
      
      outg = outg[!(outg %in% exclude)]
      
      inc = inc[!(inc %in% exclude)]
      
    } 
    
    if(!is.null(include)){
      outg = outg[outg %in% include]
      
      inc = inc[inc %in% include]
    }
    
    
    tab = data.frame(source = c(inc, rep(p, length(outg))),
                     target = c(rep(p, length(inc)), outg))
    
    if(!exists("data_e")){
      data_e = tab
    } else {
      data_e = rbind(data_e, tab)
    }
    
    interProt = c(interProt, inc, outg) %>% unique()
  }
  
  #return also vector of intermediate proteins
  if(returnIntProt){
    response = list(data_e, interProt)
  }else{
    response = data_e
  }
  
  return(response)
  
}


poi2netXpansion = function(graph = OPI_g, poi, exclude = NULL, intermediate=FALSE, include = NULL){
  
  # check if all poi are present in the graph
  p_int = intersect(V(graph) %>% as_ids(), poi)
  
  if(length(poi)!=length(p_int)){
    p_dif = setdiff(poi, V(graph) %>% as_ids())
    
    cat(paste(str_glue("{length(p_dif)} are not present in the network"), "\n"))
    
    cat(p_dif)
    print("#")
    
    poi = p_int
  }
  
  # main network expansion
  dataXpanded = xpandhelper(graph=graph, poi=poi, exclude=exclude, returnIntProt = TRUE, include = include)[[1]]
  print("dim of dataXpanded")
  
  cat(dim(dataXpanded))
  print("#")
  print("Number of vertices:")
  
  print(as.character(length(unique(c(dataXpanded$source, dataXpanded$target)))))
  
  intProt = xpandhelper(graph=graph, poi=poi, exclude=exclude, returnIntProt = TRUE, include = include)[[2]]
  if(!is.null(include)){
    intProt = c(intProt, include)
  }
  # Check link between intermediate proteins
  if(intermediate) {
    intProtXpanded = xpandhelper(
      graph = graph,
      poi = intProt,
      exclude = NULL,
      include = intProt,
      returnIntProt = FALSE
    )
    
    # Combine table
    print("dim of intProtXpanded")
    
    cat(dim(intProtXpanded))
    print("#")
    print("Number of intermediate protein:")
    print(length(intProt))
    print("Number of vertices:")
    print(as.character(length(unique(c(intProtXpanded$source, intProtXpanded$target)))))
    
    data_e = rbind(dataXpanded, intProtXpanded)
  } else{
    data_e = dataXpanded
  }
  return(data_e)
}



#' Function to annotate vertex
#' 
#' @param expanded_g data.frame, expanded network
#' @param poi.meta data.frame, metadata of protein of interest. the first column contains gene names

annotateVertex = function(expanded_g, poi.meta){
  tab = data.frame(unique(c(expanded_g$source,expanded_g$target)))
  colnames(tab) = colnames(poi.meta)[1]
  tab = unique(tab)
  anot = merge(tab, poi.meta, by=colnames(tab), all.x = TRUE)
  
  anot
}


#' transform FC to score ranging from 0 to 1 using min/max normalization
#' nlogfc (y-min(X))/(max(X)-min(X))
#' y is each value of logFC

nlogFC = function(df=int_macr_RA, col="meanLfc", ratio=TRUE){
  minx = min(df[,col])
  
  max_x = max(df[,col])
  
  if(ratio) {
    df = df %>% mutate(lfc_score = (.data[[col]]) / max_x)
  }else{
    df = df %>% mutate(lfc_score = (.data[[col]] - minx) / (max_x - minx))
  }
}



#' Function to construct network based on single disease information
#' 
#' @param graph reference interactome from which the network will be expanded
#' @param poi.meta data.frame containing metadata of protein of interest. The first column contains the official names of genes
#' @param exclude genes to be excluded. Generally down-regulated gene or non expressed genes
#' @param intermediate bool indicating whether edge between intermediate proteins must be included
#' @param walk.step The length of the random walks to perform


disease2net = function(graph = OPI_g, poi.meta, walk.step = 6, exclude = NULL, include = NULL, intermediate = FALSE, min.gene.cluster=5) {
  #poi==protein of interest
  poi = poi.meta[,1] %>% as.factor()
  
  dis_data_exp = poi2netXpansion(graph = graph, poi=poi, include = include, exclude = exclude, intermediate = intermediate)
  
  dis_v_meta = annotateVertex(expanded_g = dis_data_exp, poi.meta = poi.meta)
  
  dis_v_meta$Origin = ifelse(is.na(dis_v_meta$Origin), "Intermediate", dis_v_meta$Origin)
  
  #  set score of intermediate proteins to 0
  dis_v_meta$Score = ifelse(is.na(dis_v_meta$Score), 0, dis_v_meta$Score)
  
  # Set the edge weight to sum of their vertices' score
  dis_data_exp$weight = NA
  
  for (i in 1:nrow(dis_data_exp)) {
    w = sum(c(dis_v_meta$Score[dis_v_meta$name %in% dis_data_exp$source[i]],
              dis_v_meta$Score[dis_v_meta$name %in% dis_data_exp$target[i]]),
            na.rm = TRUE)
    
    dis_data_exp$weight[i] = w
  }
  
  net = graph_from_data_frame(d=dis_data_exp,vertices=dis_v_meta,directed=T)
  
  page.rank.net=page_rank(net, personalized = as.numeric(V(net)$Score), weights = as.numeric(E(net)$weight))
  
  V(net)$page.rank = page.rank.net$vector
  
  degree=igraph::degree(net)
  V(net)$degree = degree
  
  net.del=igraph::simplify(net,
                           remove.loops = T,
                           remove.multiple = T ,
                           edge.attr.comb = c(weight="max","ignore"))
  
  
  # Remove nodes connected to only one node
  net.del = delete_vertices(net.del, which(as.numeric(V(net.del)$degree) <= 1))
  
  cwt.net=cluster_walktrap(net.del, 
                           weights = E(net.del)$weight, 
                           steps = walk.step,
                           merges = F, 
                           modularity = TRUE, 
                           membership = TRUE)
  
  
  V(net.del)$cwt = cwt.net$membership %>% as.factor()
  
  
  print(table(cwt.net$membership))
  # Remove clusters with less than min.gene.cluster
  #cwtTab.net = data.frame(cwt.net$membership %>% table())
  #colnames(cwtTab.net)[1] = "cluster"
  
  cwtnames = table(cwt.net$membership)[table(cwt.net$membership) < min.gene.cluster] %>% names() %>% as.numeric()
  
  print(cwtnames)
  
  #vtx = V(net.del)[as.numeric(V(net.del)$cwt) %in% cwtnames]
  
  net.del = delete_vertices(net.del, V(net.del)$cwt %in% cwtnames)
  
  # Then remove nodes that are not connected and cluster
  degree=igraph::degree(net.del)
  V(net.del)$degree = degree
  
 # net.del = delete_vertices(net.del, V(net.del)[as.numeric(V(net.del)$degree) == 0])
  
  
  # Reset cluster levels to get consecutive numbers after removing some clusters
  x = V(net.del)$cwt
  x = droplevels(x)
  
  levels(x) = 1:length(levels(x))
  
  V(net.del)$cwt = x
  
  
  print("Clusters summary:")
  
  print(table(V(net.del)$cwt))
  
  gplot = ggraph(net.del, layout = "kk") +
    geom_edge_link0(width=0.2, colour = "gray90") +
    #scale_color_brewer(palette = "Paired") +
    geom_edge_fan(aes(alpha = after_stat(index)), show.legend = FALSE) +
    geom_node_point(aes(color = cwt, size = page.rank)) +
    theme_graph()
  
  plot(gplot)
  
  net.del
  
}


#' Transform graph to list based on community membership
#' 


community2list = function(graph){
  comt = V(graph)$cwt %>% levels() %>% as.numeric()
  net_list = list()
  
  for (i in comt) {
    net_list[[i]] = V(graph)[as.numeric(V(graph)$cwt) == i] %>% as_ids()
  }
  
  net_list
}

#' Jaccard
#' 

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


net2similarity = function(g1, g2){
  tab = matrix(nrow = length(g1), ncol = length(g2)) %>% as.data.frame()
  
  overlapTab = matrix(nrow = length(g1), ncol = length(g2)) %>% as.data.frame()
  
  l1 = g1 %>% unlist()
  
  l2 = g2 %>% unlist()
  
  totalgene1 = length(l1)
  totalgene2 = length(l2)
  
  
  rownames(tab) = paste0("disease1_", 1:length(g1))
  
  colnames(tab) = paste0("disease2_", 1:length(g2))
  
  for (i in 1:length(g1)) {
    
    for (j in 1:length(g2)) {
      siml = jaccard(g1[[i]], g2[[j]])
      
      tab[i,j] = siml
      # mean of overlap
      overlapTab[i,j] = (length(intersect(g1[[i]], g2[[j]]))/totalgene1 + length(intersect(g1[[i]], g2[[j]]))/totalgene2)/2
    }
    
  }
  
  overlapTab = tibble::rownames_to_column(overlapTab, var = "disease1") %>% 
    pivot_longer(!disease1, names_to = "disease2", values_to = "percent_overlap")
  
  tab = tibble::rownames_to_column(tab, var = "disease1") %>% 
    pivot_longer(!disease1, names_to = "disease2", values_to = "jaccardSim") 
  
  tab$percent_overlap = overlapTab$percent_overlap
  
  tab = filter(tab, jaccardSim > 0)
  
  # Rank by jaccardSim and percent_overlap
  # Then calculate the rank product (geometric mean of ranks)
  
  tab = tab %>% arrange(desc(jaccardSim)) %>%
    mutate(Rank1=1:nrow(tab)) %>% arrange(desc(percent_overlap)) %>% 
    mutate(Rank2=1:nrow(tab))
  # Rank product
  tab$Rank = (tab$Rank1*tab$Rank2)**0.5
  
  tab = tab %>% arrange(Rank) %>% mutate(Rank=1:nrow(tab))
  
  return(tab)
}


#' Select expressed genes in at least a fraction of sample
#' 
#' @param scedata pseudobulk data in sce format
#' @param min.count minimum count per pseudobulk. Default 10 counts
#' @param sample.size.percent percentage of sample
#' 

psdobkl2expressed = function(scedata, min.count=10, sample.size.percent = 0.1){
  
  scedata$Cell_type_scANVI = gsub(" ", "_", scedata$Cell_type_scANVI) %>% factor()
  
  scedata$Cell_type_scANVI = gsub("[:+:]", "_plus", scedata$Cell_type_scANVI) %>% factor()
  
  scedata$Cell_type_scANVI = gsub("-", "_", scedata$Cell_type_scANVI) %>% factor()
  
  scedata$donor_id = factor(gsub(" ", "_", scedata$donor_id))
  
  scedata$donor_id = factor(gsub("-", "_", scedata$donor_id))
  
  scedata$disease_short = droplevels(scedata$disease_short)
  
  scedata$tissue = gsub( " ", "_", scedata$tissue)
  
  scedata$tissue_disease = paste0(scedata$tissue, "_", scedata$disease_short) %>% factor()
  
  counts(scedata) <- assay(scedata)
  
  scedata = scedata[, scedata$disease_short != "Healthy"]
  # genes where at least some counts are more than min.count count
  n.sample = round(nrow(colData(scedata))*sample.size.percent, 0)
  non_empty_rows <- rowSums(assay(scedata) > min.count) >= n.sample
  #non_empty_rows <- which(rowSums2(assay(scedata)) > min.count)
  scedata <- scedata[non_empty_rows, ]
  
  gene = dimnames(scedata)[[1]] %>% unlist()
  
  gene
  
}


#' Barplot

barplot4netannot = function(data,
                            X,
                            Y,
                            colmap,
                            label = FALSE,
                            xlab = "",
                            ylab = "",
                            name = "Groups",
                            key.size = 0.4,
                            label_extra_position = 0,
                            flip = FALSE) {
  
  gplot = ggplot(data=data, aes(x=reorder(.data[[X]], .data[[Y]]), y=.data[[Y]], fill = .data[[X]])) +
    geom_bar(stat="identity")+
    scale_fill_manual(values = colmap) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    #coord_flip() +
    labs(x=xlab, y=ylab, fill = name) +
    theme_bw() + 
    theme(#axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      panel.grid = element_blank(),
      legend.key.width = unit(key.size, "cm"),
      legend.key.height = unit(key.size, "cm"))
  
  if(label){
    gplot = gplot + geom_text(aes(label=.data[[Y]], y = .data[[Y]] + label_extra_position), size=3.5)
  } 
  
  if(flip){
    gplot = gplot + coord_flip() + theme(
      axis.text.x = element_text(color = "black")
    )
  } else {
    gplot = gplot + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
      )
  }
  
  gplot
}


#' Function to check hue colors
plotcolor = function(col){
  scales::show_col(col)
}




#' Function to retrieve known drugs of a vector of targets genes
#' 
#' @param geneid a vector of gene symbol of ensembl id
#' @param ensId bool indicating if geneid is ensembl id
#' @import biomaRt, dplyr, ghql, jsonlite
#' 


## Set up to query Open Targets Platform API
otp_cli <- GraphqlClient$new(url = 'https://api.platform.opentargets.org/api/v4/graphql')
otp_qry <- Query$new()


gene2knowndrugs = function(geneid, ensId = FALSE) {
  
  
  otp_qry$query('query_string', 'query targetAnnotation($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    knownDrugs {
            uniqueDrugs
            rows {
              drug {
                id
                name
                maximumClinicalTrialPhase
                mechanismsOfAction {
                  uniqueActionTypes
                  uniqueTargetTypes
                }
              }
            }
          }
    }
  }'
  )
  
  
  
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  
  if(!ensId) {
    ensemblId = getBM(
      attributes = c('hgnc_symbol', 'ensembl_gene_id', 'chromosome_name'),
      filters = 'hgnc_symbol',
      values = geneid,
      mart = ensembl
    ) %>%
      filter(chromosome_name %in% as.character(1:23))
    
  } else {
    ensemblId = getBM(
      attributes = c('hgnc_symbol', 'ensembl_gene_id'),
      filters = 'ensembl_gene_id',
      values = geneid,
      mart = ensembl
    )
    
  }
  
  result = data.frame(gene.name = NA, ensemblId = NA, drug_id=NA, drug_name=NA, ClinicalTrialPhase=NA, ActionTypes=NA, TargetTypes=NA) %>% na.omit()
  
  
  
  for (i in 1:nrow(ensemblId)) {
    
    ## Execute the query
    variables <- list(ensemblId = ensemblId$ensembl_gene_id[i])
    
    
    tab <- fromJSON(otp_cli$exec(otp_qry$queries$query_string, variables, flatten = TRUE))$data$target$knownDrugs$rows$drug %>% 
      as.data.frame()
    
    
    if(nrow(tab) > 0){
      tab = data.frame(tab[,1:3], data.frame(tab$mechanismsOfAction))
      colnames(tab) = c("drug_id", "drug_name", "ClinicalTrialPhase", "ActionTypes", "TargetTypes")
      tab$ActionTypes = as.character(tab$ActionTypes)
      tab$TargetTypes = as.character(tab$TargetTypes)
      
      tab = data.frame(gene.name = ensemblId$hgnc_symbol[i], ensemblId = ensemblId$ensembl_gene_id[i], tab)
    } else {
      
      tab = data.frame(
        gene.name = ensemblId$hgnc_symbol[i],
        ensemblId = ensemblId$ensembl_gene_id[i],
        drug_id = NA,
        drug_name = NA,
        ClinicalTrialPhase = NA,
        ActionTypes=NA, 
        TargetTypes=NA
      )
    }
    
    
    result = rbind(result, tab)
    
  }
  
  
  return(result)
  
}





create_tripartite_ggplot <- function(data,
                                     drug_col = "Drug",
                                     gene_col = "Gene",
                                     indication_col = "Indication",
                                     weight_col = "Weight",
                                     text_face = "bold",
                                     node_size = 4,
                                     label_size = 3,
                                     line_alpha = 0.5,
                                     pathwidth = 0.8,
                                     color_by = "drug",
                                     drug_color = "Set2",
                                     indication_color = "Spectral",
                                     gene_color = "Paired",
                                     use_curves = TRUE,
                                     text_size = 16,
                                     show_size = TRUE,
                                     trial_dot_size = 14,
                                     size_trial_lim = 0.05,
                                     size_gene_lim = 0.05,
                                     trial_dot_color = "red",
                                     title = "Drug → Gene → Indication Tripartite Network") {
  
  
  # Rename columns for consistency
  df <- data %>%
    rename(
      Drug = !!sym(drug_col),
      Gene = !!sym(gene_col),
      Indication = !!sym(indication_col),
      Weight = !!sym(weight_col)
    )
  
  # Get unique entities
  
  drugs <- sort(unique(df$Drug))
  genes <- sort(unique(df$Gene))
  indications <- sort(unique(df$Indication))
  
  # Assign y-positions (evenly spaced)
  drug_y <- setNames(seq(1, 0, length.out = length(drugs)), drugs)
  gene_y <- setNames(seq(1, 0, length.out = length(genes)), genes)
  indication_y <- setNames(seq(1, 0, length.out = length(indications)), indications)
  
  size_trial <- df %>% select(Gene, Indication) %>% 
    unique() %>% 
    group_by(Gene) %>% mutate(size = n()) %>%
    ungroup() %>% select(Gene, size) %>% 
    unique() %>% arrange(Gene)
  
  # X positions for the three columns
  x_drug <- 0
  x_gene <- 1
  x_indication <- 2
  
  # Create node coordinates
  nodes_drug <- data.frame(
    x = x_drug,
    y = drug_y,
    label = names(drug_y),
    type = "Drug"
  )
  
  nodes_gene <- data.frame(
    x = x_gene,
    y = gene_y,
    label = names(gene_y),
    type = "Gene"
  )
  
  nodes_trial_size <- data.frame(
    x = x_gene,
    y = gene_y,
    label = as.character(size_trial$size),
    type = "Size_trial"
  )
  
  nodes_indication <- data.frame(
    x = x_indication,
    y = indication_y,
    label = names(indication_y),
    type = "Indication"
  )
  
 
  nodes <- bind_rows(nodes_drug, nodes_gene, nodes_indication, nodes_trial_size)
  nodes$type <- factor(nodes$type, levels = c("Drug", "Gene", "Size_trial", "Indication"))
  
  # Create edge data for Drug → Gene
  edges_drug_gene <- df %>%
    group_by(Drug, Gene) %>%
    summarise(Weight = sum(Weight), .groups = "drop") %>%
    mutate(
      x_start = x_drug,
      y_start = drug_y[Drug],
      x_end = x_gene,
      y_end = gene_y[Gene],
      from_type = "Drug",
      color_group = Drug
    )
  
  # Create edge data for Gene → Indication
  edges_gene_indication <- df %>%
    group_by(Gene, Indication) %>%
    summarise(Weight = sum(Weight), .groups = "drop") %>%
    mutate(
      x_start = x_gene,
      y_start = gene_y[Gene],
      x_end = x_indication,
      y_end = indication_y[Indication],
      from_type = "Gene",
      color_group = Gene
    )
  
  # Combine edges
  edges <- bind_rows(edges_drug_gene, edges_gene_indication)
  
  # Set color grouping based on user preference
  if (color_by == "drug") {
    # For gene-indication edges, map back to drugs
    drug_gene_map <- df %>%
      select(Drug, Gene) %>%
      distinct()
    
    edges_gene_indication_colored <- edges_gene_indication %>%
      left_join(drug_gene_map, by = "Gene", relationship = "many-to-many") %>%
      mutate(color_group = Drug) %>%
      select(-Drug)
    
    edges <- bind_rows(edges_drug_gene, edges_gene_indication_colored)
    color_values <- setNames(
      colorRampPalette(brewer.pal(8, drug_color))(length(drugs)),
      drugs
    )
    legend_title <- "Drug"
    
  } else if (color_by == "indication") {
    # For drug-gene edges, map forward to indications
    gene_indication_map <- df %>%
      select(Gene, Indication) %>%
      distinct()
    
    edges_drug_gene_colored <- edges_drug_gene %>%
      left_join(gene_indication_map, by = "Gene", relationship = "many-to-many") %>%
      mutate(color_group = Indication) %>%
      select(-Indication)
    
    edges_gene_indication$color_group <- edges_gene_indication$Indication
    
    edges <- bind_rows(edges_drug_gene_colored, edges_gene_indication)
    color_values <- setNames(
      colorRampPalette(brewer.pal(11, indication_color))(length(indications)),
      indications
    )
    legend_title <- "Indication"
    
  } else {
    # Color by gene
    edges_drug_gene$color_group <- edges_drug_gene$Gene
    edges_gene_indication$color_group <- edges_gene_indication$Gene
    edges <- bind_rows(edges_drug_gene, edges_gene_indication)
    color_values <- setNames(
      colorRampPalette(brewer.pal(8, gene_color))(length(genes)),
      genes
    )
    legend_title <- "Gene"
  }
  
  # Normalize weights for line thickness
  edges$line_width <- scales::rescale(edges$Weight, to = c(0.3, 2.5))
  
  # Build the plot
  p <- ggplot()
  
  # Add edges (curves or straight lines)
  if (use_curves) {
    # Create bezier curve data for smooth connections
    create_bezier <- function(x0, y0, x1, y1, n = 50) {
      t <- seq(0, 1, length.out = n)
      # Control points for cubic bezier
      cx1 <- x0 + (x1 - x0) * 0.4
      cx2 <- x0 + (x1 - x0) * 0.6
      
      x <- (1-t)^3 * x0 + 3*(1-t)^2*t * cx1 + 3*(1-t)*t^2 * cx2 + t^3 * x1
      y <- (1-t)^3 * y0 + 3*(1-t)^2*t * y0 + 3*(1-t)*t^2 * y1 + t^3 * y1
      
      data.frame(x = x, y = y)
    }
    
    # Generate curve data for each edge
    curve_data <- edges %>%
      rowwise() %>%
      mutate(curve = list(create_bezier(x_start, y_start, x_end, y_end))) %>%
      unnest(curve) %>%
      mutate(edge_id = paste(x_start, y_start, x_end, y_end, color_group))
    
    p <- p +
      geom_path(
        data = curve_data,
        aes(x = x, y = y, group = edge_id, color = color_group),
        linewidth = pathwidth,
        alpha = line_alpha
      )
    
  } else {
    # Straight lines
    p <- p +
      geom_segment(
        data = edges,
        aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
            color = color_group, linewidth = line_width),
        alpha = line_alpha
      ) +
      scale_linewidth_identity()
  }
  
  # Add nodes
  p <- p +
    geom_point(
      data = nodes,
      aes(x = x, y = y, shape = type),
      size = node_size,
      fill = "white",
      color = "black",
      stroke = 1.2
    ) +
    scale_shape_manual(
      values = c("Drug" = 16, "Gene" = 15, "Indication" = 18),
      name = "Node Type"
    )
  
  # Add labels
  p <- p +
    # Drug labels (left side)
    geom_text(
      data = filter(nodes, type == "Drug"),
      aes(x = x - 0.05, y = y, label = label),
      hjust = 1, size = label_size, fontface = text_face
    ) +
    # Gene labels (center)
    geom_text(
      data = filter(nodes, type == "Gene"),
      aes(x = x, y = y + size_gene_lim, label = label),
      hjust = 0.5, size = label_size, fontface = text_face
    ) +
    
    # Indication labels (right side)
    geom_text(
      data = filter(nodes, type == "Indication"),
      aes(x = x + 0.05, y = y, label = label),
      hjust = 0, size = label_size, fontface = text_face
    )
  
  # Number of trial
  if(show_size) {
    p <- p +
      geom_point(
        data = filter(nodes, type == "Size_trial"),
        aes(x = x, y = y - size_trial_lim), # Use the same x and adjusted y as the text
        shape = 19,              # Solid circle shape
        color = trial_dot_color,           # Set the color to red
        size = trial_dot_size                 # Adjust size to be slightly larger than text
      ) +
      geom_text(
        data = filter(nodes, type == "Size_trial"),
        aes(x = x, y = y - size_trial_lim, label = label),
        hjust = 0.5,
        vjust = 0.5,
        size = label_size,
        fontface = "bold",
        color = "white",
      )
  }
  
  
  # Apply colors and theme
  p <- p +
    scale_color_manual(values = color_values, name = legend_title) +
    scale_x_continuous(
      breaks = c(0, 1, 2),
      labels = c("Drugs", "Genes", "Indications"),
      limits = c(-0.5, 2.8)
    ) +
    labs(title = title, x = "", y = "") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = text_size+5),
      axis.text.x = element_text(size = text_size + 2, colour = "black", face = "bold"),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = text_size + 1),
      legend.text = element_text(size = text_size, face = text_face)
      #plot.margin = margin(20, 20, 20, 20)
    ) +
    guides(
      color = guide_legend(override.aes = list(linewidth = pathwidth + 1, alpha = 1)),
      shape = guide_legend(override.aes = list(size = node_size+1))
    )
  
  return(p)
}



create_tripartite_ggplot_2 <- function(data,
                                       drug_col = "Drug",
                                       gene_col = "Gene",
                                       indication_col = "Indication",
                                       weight_col = "Weight",
                                       text_face = "bold",
                                       node_size = 4,
                                       label_size = 3,
                                       line_alpha = 0.5,
                                       pathwidth = 0.8,
                                       color_by = "drug",
                                       color_name = NULL,
                                       drug_color = "Set2",
                                       indication_color = "Spectral",
                                       gene_color = "Paired",
                                       custom_color = c(
                                         "#1f77b4",
                                         "#2ca02c",
                                         "#d62728",
                                         "#9467bd",
                                         "#8c564b",
                                         "#e377c2",
                                         "#bcbd22",
                                         "#17becf",
                                         "#aec7e8",
                                         "#ff7f0e"
                                       ),
                                       use_curves = TRUE,
                                       text_size = 16,
                                       title = "Drug → Gene → Indication Tripartite Network"
) {
  
  
  # Rename columns for consistency
  df <- data %>%
    rename(
      Drug = !!sym(drug_col),
      Gene = !!sym(gene_col),
      Indication = !!sym(indication_col),
      Weight = !!sym(weight_col)
    )
  
  # Get unique entities
  
  drugs <- sort(unique(df$Drug))
  genes <- sort(unique(df$Gene))
  indications <- sort(unique(df$Indication))
  
  # Assign y-positions (evenly spaced)
  drug_y <- setNames(seq(1, 0, length.out = length(drugs)), drugs)
  gene_y <- setNames(seq(1, 0, length.out = length(genes)), genes)
  indication_y <- setNames(seq(1, 0, length.out = length(indications)), indications)
  
  # X positions for the three columns
  x_drug <- 0
  x_gene <- 1
  x_indication <- 2
  
  # Create node coordinates
  nodes_drug <- data.frame(
    x = x_drug,
    y = drug_y,
    label = names(drug_y),
    type = "Drug"
  )
  
  nodes_gene <- data.frame(
    x = x_gene,
    y = gene_y,
    label = names(gene_y),
    type = "Gene"
  )
  
  nodes_indication <- data.frame(
    x = x_indication,
    y = indication_y,
    label = names(indication_y),
    type = "Indication"
  )
  
  nodes <- bind_rows(nodes_drug, nodes_gene, nodes_indication)
  nodes$type <- factor(nodes$type, levels = c("Drug", "Gene", "Indication"))
  
  # Create edge data for Drug → Gene
  edges_drug_gene <- df %>%
    group_by(Drug, Gene) %>%
    summarise(Weight = sum(Weight), .groups = "drop") %>%
    mutate(
      x_start = x_drug,
      y_start = drug_y[Drug],
      x_end = x_gene,
      y_end = gene_y[Gene],
      from_type = "Drug",
      color_group = Drug
    )
  
  # Create edge data for Gene → Indication
  edges_gene_indication <- df %>%
    group_by(Gene, Indication) %>%
    summarise(Weight = sum(Weight), .groups = "drop") %>%
    mutate(
      x_start = x_gene,
      y_start = gene_y[Gene],
      x_end = x_indication,
      y_end = indication_y[Indication],
      from_type = "Gene",
      color_group = Gene
    )
  
  # Combine edges
  edges <- bind_rows(edges_drug_gene, edges_gene_indication)
  
  # Set color grouping based on user preference
  if (tolower(color_by) == "drug") {
    # For gene-indication edges, map back to drugs
    drug_gene_map <- df %>%
      select(Drug, Gene) %>%
      distinct()
    
    edges_gene_indication_colored <- edges_gene_indication %>%
      left_join(drug_gene_map, by = "Gene", relationship = "many-to-many") %>%
      mutate(color_group = Drug) %>%
      select(-Drug)
    
    edges <- bind_rows(edges_drug_gene, edges_gene_indication_colored)
    color_values <- setNames(
      colorRampPalette(brewer.pal(8, drug_color))(length(drugs)),
      drugs
    )
    legend_title <- "Drug"
    
  } else if (tolower(color_by) == "indication") {
    # For drug-gene edges, map forward to indications
    gene_indication_map <- df %>%
      select(Gene, Indication) %>%
      distinct()
    
    edges_drug_gene_colored <- edges_drug_gene %>%
      left_join(gene_indication_map, by = "Gene", relationship = "many-to-many") %>%
      mutate(color_group = Indication) %>%
      select(-Indication)
    
    edges_gene_indication$color_group <- edges_gene_indication$Indication
    
    edges <- bind_rows(edges_drug_gene_colored, edges_gene_indication)
    color_values <- setNames(
      colorRampPalette(brewer.pal(11, indication_color))(length(indications)),
      indications
    )
    legend_title <- "Indication"
    
  } else if (tolower(color_by) == "gene") {
    # Color by gene
    edges_drug_gene$color_group <- edges_drug_gene$Gene
    edges_gene_indication$color_group <- edges_gene_indication$Gene
    edges <- bind_rows(edges_drug_gene, edges_gene_indication)
    color_values <- setNames(
      colorRampPalette(brewer.pal(8, gene_color))(length(genes)),
      genes
    )
    legend_title <- "Gene"
    
  } else if (color_by %in% names(df)) {
    # Color by custom column
    edges_drug_gene <- df %>%
      group_by(Drug, Gene, !!sym(color_by)) %>%
      summarise(Weight = sum(Weight), .groups = "drop") %>%
      mutate(
        x_start = x_drug,
        y_start = drug_y[Drug],
        x_end = x_gene,
        y_end = gene_y[Gene],
        from_type = "Drug",
        color_group = !!sym(color_by)
      )
    
    edges_gene_indication <- df %>%
      group_by(Gene, Indication, !!sym(color_by)) %>%
      summarise(Weight = sum(Weight), .groups = "drop") %>%
      mutate(
        x_start = x_gene,
        y_start = gene_y[Gene],
        x_end = x_indication,
        y_end = indication_y[Indication],
        from_type = "Gene",
        color_group = !!sym(color_by)
      )
    
    edges <- bind_rows(edges_drug_gene, edges_gene_indication)
    
    unique_groups <- sort(unique(as.character(edges$color_group)))
    n_groups <- length(unique_groups)
    
    set.seed(123)
    
    color_values <- setNames(
      sample(custom_color, n_groups),
      unique_groups
    )
    
    if(!is.null(color_name)) {
      legend_title <- color_name
    } else {
      legend_title <- color_by
    }
    
    
  } else {
    stop(paste0("color_by column '", color_by, "' not found in data.\n",
                "Available columns: ", paste(names(df), collapse = ", "), "\n",
                "Or use: 'drug', 'gene', 'indication'"))
  }
  
  # Normalize weights for line thickness
  edges$line_width <- scales::rescale(edges$Weight, to = c(0.3, 2.5))
  
  # Build the plot
  p <- ggplot()
  
  # Add edges (curves or straight lines)
  if (use_curves) {
    # Create bezier curve data for smooth connections
    create_bezier <- function(x0, y0, x1, y1, n = 50) {
      t <- seq(0, 1, length.out = n)
      # Control points for cubic bezier
      cx1 <- x0 + (x1 - x0) * 0.4
      cx2 <- x0 + (x1 - x0) * 0.6
      
      x <- (1-t)^3 * x0 + 3*(1-t)^2*t * cx1 + 3*(1-t)*t^2 * cx2 + t^3 * x1
      y <- (1-t)^3 * y0 + 3*(1-t)^2*t * y0 + 3*(1-t)*t^2 * y1 + t^3 * y1
      
      data.frame(x = x, y = y)
    }
    
    # Generate curve data for each edge
    curve_data <- edges %>%
      rowwise() %>%
      mutate(curve = list(create_bezier(x_start, y_start, x_end, y_end))) %>%
      unnest(curve) %>%
      mutate(edge_id = paste(x_start, y_start, x_end, y_end, color_group))
    
    p <- p +
      geom_path(
        data = curve_data,
        aes(x = x, y = y, group = edge_id, color = color_group),
        linewidth = pathwidth,
        alpha = line_alpha
      )
    
  } else {
    # Straight lines
    p <- p +
      geom_segment(
        data = edges,
        aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
            color = color_group, linewidth = line_width),
        alpha = line_alpha
      ) +
      scale_linewidth_identity()
  }
  
  # Add nodes
  p <- p +
    geom_point(
      data = nodes,
      aes(x = x, y = y, shape = type),
      size = node_size,
      fill = "white",
      color = "black",
      stroke = 1.2
    ) +
    scale_shape_manual(
      values = c("Drug" = 16, "Gene" = 15, "Indication" = 18),
      name = "Node Type"
    )
  
  # Add labels
  p <- p +
    # Drug labels (left side)
    geom_text(
      data = filter(nodes, type == "Drug"),
      aes(x = x - 0.05, y = y, label = label),
      hjust = 1, size = label_size, fontface = text_face
    ) +
    # Gene labels (center)
    geom_text(
      data = filter(nodes, type == "Gene"),
      aes(x = x, y = y + 0.03, label = label),
      hjust = 0.5, size = label_size, fontface = text_face
    ) +
    # Indication labels (right side)
    geom_text(
      data = filter(nodes, type == "Indication"),
      aes(x = x + 0.05, y = y, label = label),
      hjust = 0, size = label_size, fontface = text_face
    )
  
  # Apply colors and theme
  p <- p +
    scale_color_manual(values = color_values, name = legend_title) +
    scale_x_continuous(
      breaks = c(0, 1, 2),
      labels = c("Drugs", "Genes", "Indications"),
      limits = c(-0.5, 2.8)
    ) +
    labs(title = title, x = "", y = "") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = text_size+2),
      axis.text.x = element_text(size = text_size + 2, colour = "black", face = "bold"),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = text_size + 1),
      legend.text = element_text(size = text_size, face = text_face)
    ) +
    guides(
      color = guide_legend(override.aes = list(linewidth = pathwidth + 1, alpha = 1)),
      shape = guide_legend(override.aes = list(size = node_size+1))
    )
  
  return(p)
}



#' Function to check hue colors
plotcolor = function(col){
  scales::show_col(col)
}




#' Function to retrieve known drugs of a vector of targets genes
#' 
#' @param geneid a vector of gene symbol of ensembl id
#' @param ensId bool indicating if geneid is ensembl id
#' @import biomaRt, dplyr, ghql, jsonlite
#' 


## Set up to query Open Targets Platform API
otp_cli <- GraphqlClient$new(url = 'https://api.platform.opentargets.org/api/v4/graphql')
otp_qry <- Query$new()


gene2knowndrugs = function(geneid, ensId = FALSE) {
  
  
  otp_qry$query('query_string', 'query targetAnnotation($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    knownDrugs {
            uniqueDrugs
            rows {
              drug {
                id
                name
                maximumClinicalTrialPhase
                mechanismsOfAction {
                  uniqueActionTypes
                  uniqueTargetTypes
                }
              }
            }
          }
    }
  }'
  )
  
  
  
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  
  if(!ensId) {
    ensemblId = getBM(
      attributes = c('hgnc_symbol', 'ensembl_gene_id', 'chromosome_name'),
      filters = 'hgnc_symbol',
      values = geneid,
      mart = ensembl
    ) %>%
      filter(chromosome_name %in% as.character(1:23))
    
  } else {
    ensemblId = getBM(
      attributes = c('hgnc_symbol', 'ensembl_gene_id'),
      filters = 'ensembl_gene_id',
      values = geneid,
      mart = ensembl
    )
    
  }
  
  result = data.frame(gene.name = NA, ensemblId = NA, drug_id=NA, drug_name=NA, ClinicalTrialPhase=NA, ActionTypes=NA, TargetTypes=NA) %>% na.omit()
  
  
  
  for (i in 1:nrow(ensemblId)) {
    
    ## Execute the query
    variables <- list(ensemblId = ensemblId$ensembl_gene_id[i])
    
    
    tab <- fromJSON(otp_cli$exec(otp_qry$queries$query_string, variables, flatten = TRUE))$data$target$knownDrugs$rows$drug %>% 
      as.data.frame()
    
    
    if(nrow(tab) > 0){
      tab = data.frame(tab[,1:3], data.frame(tab$mechanismsOfAction))
      colnames(tab) = c("drug_id", "drug_name", "ClinicalTrialPhase", "ActionTypes", "TargetTypes")
      tab$ActionTypes = as.character(tab$ActionTypes)
      tab$TargetTypes = as.character(tab$TargetTypes)
      
      tab = data.frame(gene.name = ensemblId$hgnc_symbol[i], ensemblId = ensemblId$ensembl_gene_id[i], tab)
    } else {
      
      tab = data.frame(
        gene.name = ensemblId$hgnc_symbol[i],
        ensemblId = ensemblId$ensembl_gene_id[i],
        drug_id = NA,
        drug_name = NA,
        ClinicalTrialPhase = NA,
        ActionTypes=NA, 
        TargetTypes=NA
      )
    }
    
    
    result = rbind(result, tab)
    
  }
  
  
  return(result)
  
}





create_tripartite_ggplot <- function(data,
                                     drug_col = "Drug",
                                     gene_col = "Gene",
                                     indication_col = "Indication",
                                     weight_col = "Weight",
                                     text_face = "bold",
                                     node_size = 4,
                                     label_size = 3,
                                     line_alpha = 0.5,
                                     pathwidth = 0.8,
                                     color_by = "drug",
                                     drug_color = "Set2",
                                     indication_color = "Spectral",
                                     gene_color = "Paired",
                                     use_curves = TRUE,
                                     text_size = 16,
                                     show_size = TRUE,
                                     trial_dot_size = 14,
                                     size_trial_lim = 0.05,
                                     size_gene_lim = 0.05,
                                     trial_dot_color = "red",
                                     title = "Drug → Gene → Indication Tripartite Network") {
  
  
  # Rename columns for consistency
  df <- data %>%
    rename(
      Drug = !!sym(drug_col),
      Gene = !!sym(gene_col),
      Indication = !!sym(indication_col),
      Weight = !!sym(weight_col)
    )
  
  # Get unique entities
  
  drugs <- sort(unique(df$Drug))
  genes <- sort(unique(df$Gene))
  indications <- sort(unique(df$Indication))
  
  # Assign y-positions (evenly spaced)
  drug_y <- setNames(seq(1, 0, length.out = length(drugs)), drugs)
  gene_y <- setNames(seq(1, 0, length.out = length(genes)), genes)
  indication_y <- setNames(seq(1, 0, length.out = length(indications)), indications)
  
  size_trial <- df %>% select(Gene, Indication) %>% 
    unique() %>% 
    group_by(Gene) %>% mutate(size = n()) %>%
    ungroup() %>% select(Gene, size) %>% 
    unique() %>% arrange(Gene)
  
  # X positions for the three columns
  x_drug <- 0
  x_gene <- 1
  x_indication <- 2
  
  # Create node coordinates
  nodes_drug <- data.frame(
    x = x_drug,
    y = drug_y,
    label = names(drug_y),
    type = "Drug"
  )
  
  nodes_gene <- data.frame(
    x = x_gene,
    y = gene_y,
    label = names(gene_y),
    type = "Gene"
  )
  
  nodes_trial_size <- data.frame(
    x = x_gene,
    y = gene_y,
    label = as.character(size_trial$size),
    type = "Size_trial"
  )
  
  nodes_indication <- data.frame(
    x = x_indication,
    y = indication_y,
    label = names(indication_y),
    type = "Indication"
  )
  
 
  nodes <- bind_rows(nodes_drug, nodes_gene, nodes_indication, nodes_trial_size)
  nodes$type <- factor(nodes$type, levels = c("Drug", "Gene", "Size_trial", "Indication"))
  
  # Create edge data for Drug → Gene
  edges_drug_gene <- df %>%
    group_by(Drug, Gene) %>%
    summarise(Weight = sum(Weight), .groups = "drop") %>%
    mutate(
      x_start = x_drug,
      y_start = drug_y[Drug],
      x_end = x_gene,
      y_end = gene_y[Gene],
      from_type = "Drug",
      color_group = Drug
    )
  
  # Create edge data for Gene → Indication
  edges_gene_indication <- df %>%
    group_by(Gene, Indication) %>%
    summarise(Weight = sum(Weight), .groups = "drop") %>%
    mutate(
      x_start = x_gene,
      y_start = gene_y[Gene],
      x_end = x_indication,
      y_end = indication_y[Indication],
      from_type = "Gene",
      color_group = Gene
    )
  
  # Combine edges
  edges <- bind_rows(edges_drug_gene, edges_gene_indication)
  
  # Set color grouping based on user preference
  if (color_by == "drug") {
    # For gene-indication edges, map back to drugs
    drug_gene_map <- df %>%
      select(Drug, Gene) %>%
      distinct()
    
    edges_gene_indication_colored <- edges_gene_indication %>%
      left_join(drug_gene_map, by = "Gene", relationship = "many-to-many") %>%
      mutate(color_group = Drug) %>%
      select(-Drug)
    
    edges <- bind_rows(edges_drug_gene, edges_gene_indication_colored)
    color_values <- setNames(
      colorRampPalette(brewer.pal(8, drug_color))(length(drugs)),
      drugs
    )
    legend_title <- "Drug"
    
  } else if (color_by == "indication") {
    # For drug-gene edges, map forward to indications
    gene_indication_map <- df %>%
      select(Gene, Indication) %>%
      distinct()
    
    edges_drug_gene_colored <- edges_drug_gene %>%
      left_join(gene_indication_map, by = "Gene", relationship = "many-to-many") %>%
      mutate(color_group = Indication) %>%
      select(-Indication)
    
    edges_gene_indication$color_group <- edges_gene_indication$Indication
    
    edges <- bind_rows(edges_drug_gene_colored, edges_gene_indication)
    color_values <- setNames(
      colorRampPalette(brewer.pal(11, indication_color))(length(indications)),
      indications
    )
    legend_title <- "Indication"
    
  } else {
    # Color by gene
    edges_drug_gene$color_group <- edges_drug_gene$Gene
    edges_gene_indication$color_group <- edges_gene_indication$Gene
    edges <- bind_rows(edges_drug_gene, edges_gene_indication)
    color_values <- setNames(
      colorRampPalette(brewer.pal(8, gene_color))(length(genes)),
      genes
    )
    legend_title <- "Gene"
  }
  
  # Normalize weights for line thickness
  edges$line_width <- scales::rescale(edges$Weight, to = c(0.3, 2.5))
  
  # Build the plot
  p <- ggplot()
  
  # Add edges (curves or straight lines)
  if (use_curves) {
    # Create bezier curve data for smooth connections
    create_bezier <- function(x0, y0, x1, y1, n = 50) {
      t <- seq(0, 1, length.out = n)
      # Control points for cubic bezier
      cx1 <- x0 + (x1 - x0) * 0.4
      cx2 <- x0 + (x1 - x0) * 0.6
      
      x <- (1-t)^3 * x0 + 3*(1-t)^2*t * cx1 + 3*(1-t)*t^2 * cx2 + t^3 * x1
      y <- (1-t)^3 * y0 + 3*(1-t)^2*t * y0 + 3*(1-t)*t^2 * y1 + t^3 * y1
      
      data.frame(x = x, y = y)
    }
    
    # Generate curve data for each edge
    curve_data <- edges %>%
      rowwise() %>%
      mutate(curve = list(create_bezier(x_start, y_start, x_end, y_end))) %>%
      unnest(curve) %>%
      mutate(edge_id = paste(x_start, y_start, x_end, y_end, color_group))
    
    p <- p +
      geom_path(
        data = curve_data,
        aes(x = x, y = y, group = edge_id, color = color_group),
        linewidth = pathwidth,
        alpha = line_alpha
      )
    
  } else {
    # Straight lines
    p <- p +
      geom_segment(
        data = edges,
        aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
            color = color_group, linewidth = line_width),
        alpha = line_alpha
      ) +
      scale_linewidth_identity()
  }
  
  # Add nodes
  p <- p +
    geom_point(
      data = nodes,
      aes(x = x, y = y, shape = type),
      size = node_size,
      fill = "white",
      color = "black",
      stroke = 1.2
    ) +
    scale_shape_manual(
      values = c("Drug" = 16, "Gene" = 15, "Indication" = 18),
      name = "Node Type"
    )
  
  # Add labels
  p <- p +
    # Drug labels (left side)
    geom_text(
      data = filter(nodes, type == "Drug"),
      aes(x = x - 0.05, y = y, label = label),
      hjust = 1, size = label_size, fontface = text_face
    ) +
    # Gene labels (center)
    geom_text(
      data = filter(nodes, type == "Gene"),
      aes(x = x, y = y + size_gene_lim, label = label),
      hjust = 0.5, size = label_size, fontface = text_face
    ) +
    
    # Indication labels (right side)
    geom_text(
      data = filter(nodes, type == "Indication"),
      aes(x = x + 0.05, y = y, label = label),
      hjust = 0, size = label_size, fontface = text_face
    )
  
  # Number of trial
  if(show_size) {
    p <- p +
      geom_point(
        data = filter(nodes, type == "Size_trial"),
        aes(x = x, y = y - size_trial_lim), # Use the same x and adjusted y as the text
        shape = 19,              # Solid circle shape
        color = trial_dot_color,           # Set the color to red
        size = trial_dot_size                 # Adjust size to be slightly larger than text
      ) +
      geom_text(
        data = filter(nodes, type == "Size_trial"),
        aes(x = x, y = y - size_trial_lim, label = label),
        hjust = 0.5,
        vjust = 0.5,
        size = label_size,
        fontface = "bold",
        color = "white",
      )
  }
  
  
  # Apply colors and theme
  p <- p +
    scale_color_manual(values = color_values, name = legend_title) +
    scale_x_continuous(
      breaks = c(0, 1, 2),
      labels = c("Drugs", "Genes", "Indications"),
      limits = c(-0.5, 2.8)
    ) +
    labs(title = title, x = "", y = "") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = text_size+5),
      axis.text.x = element_text(size = text_size + 2, colour = "black", face = "bold"),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = text_size + 1),
      legend.text = element_text(size = text_size, face = text_face)
      #plot.margin = margin(20, 20, 20, 20)
    ) +
    guides(
      color = guide_legend(override.aes = list(linewidth = pathwidth + 1, alpha = 1)),
      shape = guide_legend(override.aes = list(size = node_size+1))
    )
  
  return(p)
}



create_tripartite_ggplot_2 <- function(data,
                                       drug_col = "Drug",
                                       gene_col = "Gene",
                                       indication_col = "Indication",
                                       weight_col = "Weight",
                                       text_face = "bold",
                                       node_size = 4,
                                       label_size = 3,
                                       line_alpha = 0.5,
                                       pathwidth = 0.8,
                                       color_by = "drug",
                                       color_name = NULL,
                                       drug_color = "Set2",
                                       indication_color = "Spectral",
                                       gene_color = "Paired",
                                       custom_color = c(
                                         "#1f77b4",
                                         "#2ca02c",
                                         "#d62728",
                                         "#9467bd",
                                         "#8c564b",
                                         "#e377c2",
                                         "#bcbd22",
                                         "#17becf",
                                         "#aec7e8",
                                         "#ff7f0e"
                                       ),
                                       use_curves = TRUE,
                                       text_size = 16,
                                       title = "Drug → Gene → Indication Tripartite Network"
) {
  
  
  # Rename columns for consistency
  df <- data %>%
    rename(
      Drug = !!sym(drug_col),
      Gene = !!sym(gene_col),
      Indication = !!sym(indication_col),
      Weight = !!sym(weight_col)
    )
  
  # Get unique entities
  
  drugs <- sort(unique(df$Drug))
  genes <- sort(unique(df$Gene))
  indications <- sort(unique(df$Indication))
  
  # Assign y-positions (evenly spaced)
  drug_y <- setNames(seq(1, 0, length.out = length(drugs)), drugs)
  gene_y <- setNames(seq(1, 0, length.out = length(genes)), genes)
  indication_y <- setNames(seq(1, 0, length.out = length(indications)), indications)
  
  # X positions for the three columns
  x_drug <- 0
  x_gene <- 1
  x_indication <- 2
  
  # Create node coordinates
  nodes_drug <- data.frame(
    x = x_drug,
    y = drug_y,
    label = names(drug_y),
    type = "Drug"
  )
  
  nodes_gene <- data.frame(
    x = x_gene,
    y = gene_y,
    label = names(gene_y),
    type = "Gene"
  )
  
  nodes_indication <- data.frame(
    x = x_indication,
    y = indication_y,
    label = names(indication_y),
    type = "Indication"
  )
  
  nodes <- bind_rows(nodes_drug, nodes_gene, nodes_indication)
  nodes$type <- factor(nodes$type, levels = c("Drug", "Gene", "Indication"))
  
  # Create edge data for Drug → Gene
  edges_drug_gene <- df %>%
    group_by(Drug, Gene) %>%
    summarise(Weight = sum(Weight), .groups = "drop") %>%
    mutate(
      x_start = x_drug,
      y_start = drug_y[Drug],
      x_end = x_gene,
      y_end = gene_y[Gene],
      from_type = "Drug",
      color_group = Drug
    )
  
  # Create edge data for Gene → Indication
  edges_gene_indication <- df %>%
    group_by(Gene, Indication) %>%
    summarise(Weight = sum(Weight), .groups = "drop") %>%
    mutate(
      x_start = x_gene,
      y_start = gene_y[Gene],
      x_end = x_indication,
      y_end = indication_y[Indication],
      from_type = "Gene",
      color_group = Gene
    )
  
  # Combine edges
  edges <- bind_rows(edges_drug_gene, edges_gene_indication)
  
  # Set color grouping based on user preference
  if (tolower(color_by) == "drug") {
    # For gene-indication edges, map back to drugs
    drug_gene_map <- df %>%
      select(Drug, Gene) %>%
      distinct()
    
    edges_gene_indication_colored <- edges_gene_indication %>%
      left_join(drug_gene_map, by = "Gene", relationship = "many-to-many") %>%
      mutate(color_group = Drug) %>%
      select(-Drug)
    
    edges <- bind_rows(edges_drug_gene, edges_gene_indication_colored)
    color_values <- setNames(
      colorRampPalette(brewer.pal(8, drug_color))(length(drugs)),
      drugs
    )
    legend_title <- "Drug"
    
  } else if (tolower(color_by) == "indication") {
    # For drug-gene edges, map forward to indications
    gene_indication_map <- df %>%
      select(Gene, Indication) %>%
      distinct()
    
    edges_drug_gene_colored <- edges_drug_gene %>%
      left_join(gene_indication_map, by = "Gene", relationship = "many-to-many") %>%
      mutate(color_group = Indication) %>%
      select(-Indication)
    
    edges_gene_indication$color_group <- edges_gene_indication$Indication
    
    edges <- bind_rows(edges_drug_gene_colored, edges_gene_indication)
    color_values <- setNames(
      colorRampPalette(brewer.pal(11, indication_color))(length(indications)),
      indications
    )
    legend_title <- "Indication"
    
  } else if (tolower(color_by) == "gene") {
    # Color by gene
    edges_drug_gene$color_group <- edges_drug_gene$Gene
    edges_gene_indication$color_group <- edges_gene_indication$Gene
    edges <- bind_rows(edges_drug_gene, edges_gene_indication)
    color_values <- setNames(
      colorRampPalette(brewer.pal(8, gene_color))(length(genes)),
      genes
    )
    legend_title <- "Gene"
    
  } else if (color_by %in% names(df)) {
    # Color by custom column
    edges_drug_gene <- df %>%
      group_by(Drug, Gene, !!sym(color_by)) %>%
      summarise(Weight = sum(Weight), .groups = "drop") %>%
      mutate(
        x_start = x_drug,
        y_start = drug_y[Drug],
        x_end = x_gene,
        y_end = gene_y[Gene],
        from_type = "Drug",
        color_group = !!sym(color_by)
      )
    
    edges_gene_indication <- df %>%
      group_by(Gene, Indication, !!sym(color_by)) %>%
      summarise(Weight = sum(Weight), .groups = "drop") %>%
      mutate(
        x_start = x_gene,
        y_start = gene_y[Gene],
        x_end = x_indication,
        y_end = indication_y[Indication],
        from_type = "Gene",
        color_group = !!sym(color_by)
      )
    
    edges <- bind_rows(edges_drug_gene, edges_gene_indication)
    
    unique_groups <- sort(unique(as.character(edges$color_group)))
    n_groups <- length(unique_groups)
    
    set.seed(123)
    
    color_values <- setNames(
      sample(custom_color, n_groups),
      unique_groups
    )
    
    if(!is.null(color_name)) {
      legend_title <- color_name
    } else {
      legend_title <- color_by
    }
    
    
  } else {
    stop(paste0("color_by column '", color_by, "' not found in data.\n",
                "Available columns: ", paste(names(df), collapse = ", "), "\n",
                "Or use: 'drug', 'gene', 'indication'"))
  }
  
  # Normalize weights for line thickness
  edges$line_width <- scales::rescale(edges$Weight, to = c(0.3, 2.5))
  
  # Build the plot
  p <- ggplot()
  
  # Add edges (curves or straight lines)
  if (use_curves) {
    # Create bezier curve data for smooth connections
    create_bezier <- function(x0, y0, x1, y1, n = 50) {
      t <- seq(0, 1, length.out = n)
      # Control points for cubic bezier
      cx1 <- x0 + (x1 - x0) * 0.4
      cx2 <- x0 + (x1 - x0) * 0.6
      
      x <- (1-t)^3 * x0 + 3*(1-t)^2*t * cx1 + 3*(1-t)*t^2 * cx2 + t^3 * x1
      y <- (1-t)^3 * y0 + 3*(1-t)^2*t * y0 + 3*(1-t)*t^2 * y1 + t^3 * y1
      
      data.frame(x = x, y = y)
    }
    
    # Generate curve data for each edge
    curve_data <- edges %>%
      rowwise() %>%
      mutate(curve = list(create_bezier(x_start, y_start, x_end, y_end))) %>%
      unnest(curve) %>%
      mutate(edge_id = paste(x_start, y_start, x_end, y_end, color_group))
    
    p <- p +
      geom_path(
        data = curve_data,
        aes(x = x, y = y, group = edge_id, color = color_group),
        linewidth = pathwidth,
        alpha = line_alpha
      )
    
  } else {
    # Straight lines
    p <- p +
      geom_segment(
        data = edges,
        aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
            color = color_group, linewidth = line_width),
        alpha = line_alpha
      ) +
      scale_linewidth_identity()
  }
  
  # Add nodes
  p <- p +
    geom_point(
      data = nodes,
      aes(x = x, y = y, shape = type),
      size = node_size,
      fill = "white",
      color = "black",
      stroke = 1.2
    ) +
    scale_shape_manual(
      values = c("Drug" = 16, "Gene" = 15, "Indication" = 18),
      name = "Node Type"
    )
  
  # Add labels
  p <- p +
    # Drug labels (left side)
    geom_text(
      data = filter(nodes, type == "Drug"),
      aes(x = x - 0.05, y = y, label = label),
      hjust = 1, size = label_size, fontface = text_face
    ) +
    # Gene labels (center)
    geom_text(
      data = filter(nodes, type == "Gene"),
      aes(x = x, y = y + 0.03, label = label),
      hjust = 0.5, size = label_size, fontface = text_face
    ) +
    # Indication labels (right side)
    geom_text(
      data = filter(nodes, type == "Indication"),
      aes(x = x + 0.05, y = y, label = label),
      hjust = 0, size = label_size, fontface = text_face
    )
  
  # Apply colors and theme
  p <- p +
    scale_color_manual(values = color_values, name = legend_title) +
    scale_x_continuous(
      breaks = c(0, 1, 2),
      labels = c("Drugs", "Genes", "Indications"),
      limits = c(-0.5, 2.8)
    ) +
    labs(title = title, x = "", y = "") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = text_size+2),
      axis.text.x = element_text(size = text_size + 2, colour = "black", face = "bold"),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = text_size + 1),
      legend.text = element_text(size = text_size, face = text_face)
    ) +
    guides(
      color = guide_legend(override.aes = list(linewidth = pathwidth + 1, alpha = 1)),
      shape = guide_legend(override.aes = list(size = node_size+1))
    )
  
  return(p)
}
