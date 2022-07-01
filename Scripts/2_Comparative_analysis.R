library(tidyverse)
library(tidygraph)
library(ggraph)
library(visNetwork)
library(igraph)
library(stringr)
library(beepr)
library(rlang)
library(rstatix)
library(ggpubr)
library(ggdist)
library(scales)
library(apcluster)

####################
# Сравнение графов #
####################

compare_graphs_with_formula <- function(graph_1, graph_2, names) {
  
  # FIXME to remove after
  # graph_1 <- `graph_list_w-plots_3136192901`[[1]]
  # graph_2 <- `graph_list_w-plots_3136217149`[[1]]
  
  comparison_graphs_list <- list()
  
  nodes_1 <- graph_1 %>%
    activate(nodes) %>%
    as_tibble()
  
  nodes_2 <- graph_2 %>%
    activate(nodes) %>%
    as_tibble() %>%
    mutate(id = id + nrow(nodes_1))
  
  diff <- nodes_1 %>%
    full_join(nodes_2, by = "label") %>%
    as_tibble() %>%
    mutate(diff = abs(coef.x - coef.y)) %>%
    pull(diff)
  
  nodes <- nodes_1 %>%
    bind_rows(nodes_2, .id = "graph") %>%
    mutate(graph = ifelse(graph == 1, names[1], names[2])) %>%
    mutate(diff = rep(diff, 2))
  
  edges_1 <- graph_1 %>%
    activate(edges) %>%
    as_tibble()
  
  edges_2 <- graph_2 %>%
    activate(edges) %>%
    as_tibble()
  
  # edges_joined <- edges_1 %>%
  #   full_join(edges_2, by = c("from", "to")) %>%
  #   mutate(join_weight = ifelse(!is.na(new_weight.x)&!is.na(new_weight.y),
  #                               new_weight.x + new_weight.y,
  #                                ifelse(is.na(new_weight.x), new_weight.y, new_weight.x))) %>%
  #   distinct() %>%
  #   select(from, to, new_weight.x, new_weight.y, join_weight) %>%
  #   mutate(graph_type = ifelse(is.na(new_weight.x)|is.na(new_weight.y), TRUE, FALSE)) %>%
  #   select(-new_weight.x, -new_weight.y)
  
  edges_joined <- edges_1 %>%
    full_join(edges_2, by = c("from", "to")) %>%
    mutate(join_weight = ifelse(!is.na(new_weight.x)&!is.na(new_weight.y),
                                abs(new_weight.x - new_weight.y), 0)) %>%
    distinct() %>%
    select(from, to, new_weight.x, new_weight.y, join_weight) %>%
    mutate(graph_type = ifelse(is.na(new_weight.x)|is.na(new_weight.y), TRUE, FALSE)) %>%
    select(-new_weight.x, -new_weight.y)
  
  edges_1 <- edges_1 %>%
    left_join(edges_joined, by = c("from", "to")) %>%
    distinct() 
  edges_2 <- edges_2 %>%
    left_join(edges_joined, by = c("from", "to")) %>%
    mutate(from = from + nrow(nodes_1), to = to + nrow(nodes_1)) %>%
    distinct()
  edges <- edges_1 %>%
    bind_rows(edges_2, .id = "graph") %>%
    mutate(graph = ifelse(graph == 1, names[1], names[2]))
  
  graph <- tbl_graph(
    nodes = nodes, 
    edges = edges,
    directed = F
  )
  
  comparison_graphs_list[[1]] <- graph
  
  nodes_xy <- graph %>%
    activate(nodes) %>%
    as_tibble() %>%
    select(x,y) 
  names(nodes_xy) <- c("x1", "y1")
  
  graph_1a <- tbl_graph(
    nodes = nodes_1, 
    edges = edges_1,
    directed = F
  )
  
  edges_2 <- edges_2 %>%
    mutate(from = from - nrow(nodes_1), to = to - nrow(nodes_1))
  
  graph_2a <- tbl_graph(
    nodes = nodes_2, 
    edges = edges_2,
    directed = F
  )
  
  node1 <- V(graph_1a)$label
  node2 <- V(graph_2a)$label
  node_index <- c(node1,node2)
  
  matrix_1 <- as.matrix(as_adj(graph_1a, names = T)) * as.numeric(nodes_1$coef)
  matrix_2 <- as.matrix(as_adj(graph_2a, names = T)) * as.numeric(nodes_2$coef)
  
  all_net_matrix <- abs(matrix_1 - matrix_2)
  
  comparison_graphs_list[[2]] <- all_net_matrix
  
  cluster_scores <- list()
  sim_scores <- c()

  for (i in 1:10) {

    clusters <- apcluster(all_net_matrix)

    # calc_cluster_score <- function(cluster){
    #   LS1 <- sapply(cluster,function(x)median(matrix_1[x, cluster])) %>% sum
    #   LS2 <- sapply(cluster,function(x)median(matrix_2[x, cluster])) %>% sum
    #   LS <- median(LS1, LS2)
    #   return(LS)
    # }

    calc_cluster_score <- function(cluster){
      LS <- sapply(cluster,function(x) max(all_net_matrix[x, cluster])) %>% mean()
      return(LS)
    }

    score <- sapply(1:length(clusters),function(x)1-calc_cluster_score(clusters[[x]]))
    cluster_score <- data.frame(score = score, cluster = sapply(clusters@clusters, function(x) paste0(x, collapse = " ")))

    cluster_scores[[i]] <- cluster_score

    sim_score <- median(cluster_score$score)
    sim_scores[i] <- sim_score

  }

  sim_score_m <- median(sim_scores, na.rm = T)

  comparison_graphs_list[[3]] <- cluster_scores
  comparison_graphs_list[[4]] <- sim_scores
  comparison_graphs_list[[5]] <- sim_score_m
  
  diff_95 <- graph %>%
    activate(nodes) %>%
    as_tibble() %>%
    summarise(q95 = quantile(diff, .95)) %>%
    pull()
  
  graph_plot_g <- graph %>%
    activate(nodes) %>%
    mutate(label = ifelse(diff > diff_95, label, NA)) %>%
    mutate(coef_2 = ifelse(diff > diff_95, coef*2, coef),
           diff = ifelse(diff > diff_95, 1, 1))
  
  graph_plot <- graph_plot_g %>%
    ggraph(layout = nodes_xy) + 
    geom_edge_link(aes(color = graph, alpha = new_weight, linetype = graph_type), width = 1.5) +
    geom_node_point(aes(color = graph, size = coef_2, alpha = diff)) +
    geom_node_text(aes(label = label, size = coef_2), repel = T) +
    scale_edge_linetype_manual(values=c("solid", "dotdash")) +
    scale_color_manual(values = c("yellow", "blue"), aesthetics = "color") +
    scale_size_continuous(trans = "log10") +
    facet_nodes(~ graph) +
    theme_graph() +
    theme(plot.title = element_text(hjust = 0.5))
  
  comparison_graphs_list[[3]] <- graph_plot
  # comparison_graphs_list[[6]] <- graph_plot
  
  ggsave(filename = str_c("comparison_graphs_list", str_c(names, collapse = ", "), "_2.png", sep = "_"), plot = graph_plot, scale = 3)
  
  assign(x = str_c("comparison_graphs_list", str_c(names, collapse = ", "), sep = "_"), value = comparison_graphs_list, envir = .GlobalEnv)
  
}  
