files_to_import <- c("1890760379","1846921335","1717818438")
files_to_import_or <- str_c(files_to_import, collapse = "|")
system_bonds <- c("328_350","303_328","no_bonds")

system_info <- data.frame(system_id = files_to_import, bonds = system_bonds)

nodes_files <- list.files(path = "Data/NaPi2b/str_data/")
edges_files <- list.files(path = "Data/NaPi2b/aminoacids_interactions/")
nodes_to_import <- nodes_files[str_detect(nodes_files, files_to_import_or)]
edges_to_import <- edges_files[str_detect(edges_files, files_to_import_or)]

# Graph - control

for (i in 1:length(files_to_import)) {
  
  name <- system_info[i,1]
  system_name <- system_info[i,2]
  filenames_c <- find_file(name)
  nodes <- read_csv(filenames_c[1]) %>%
    dplyr::mutate(resid = str_c(resno, resid, sep = "_"))
  names(nodes)[1:2] <- c('id', 'label')
  edges <- read_csv(filenames_c[2]) %>%
    dplyr::filter(persent_ring != 0) %>%
    dplyr::mutate(five_aa = ifelse(abs(number.x - number.y) >= 5, 2, 1))
  names(edges) <- c('to', 'from', 'weight', 'five_aa')
  
  names_not_parameters <- str_detect(names(nodes), 
                                     str_c(c("id", "label", "x", "y", "z", "type", "topology",
                                             "amino", "importance_ring"), collapse = "|"))
  parameters <- c("default", names(nodes)[!names_not_parameters])
  
  get_networks(nodes, edges, parameters)
  get_network_parameters(graph_list)
  compare_node_parameters(parameter_list)
  stating_parameter_formula(comparison_node_parameters)
  graph_with_formula(nodes, edges, formula_coef, system_name)
  
}

beep()

################################
# Compare control with 328-350 #
################################

# compare_graphs_with_formula(`graph_list_w-plots_no_bonds`[[1]],
#                             `graph_list_w-plots_328_350`[[1]],
#                             system_bonds[c(3,1)])

quant_graph_c_328350 <- `comparison_graphs_list_no_bonds, 328_350`[[1]] %>%
  activate(nodes) %>%
  as_tibble() %>%
  summarise(quan_75_coef = quantile(diff, probs = 0.95)) %>%
  pull() %>%
  unname()

nodes_c_328350_1 <- `comparison_graphs_list_no_bonds, 328_350`[[1]] %>%
  activate(nodes) %>%
  as_tibble() %>%
  head(690) 

nodes_c_328350_2 <- `comparison_graphs_list_no_bonds, 328_350`[[1]] %>%
  activate(nodes) %>%
  as_tibble() %>%
  tail(690)

imp_aa_c_328350 <- nodes_c_328350_1 %>%
  full_join(nodes_c_328350_2, by = c("label", "diff")) %>%
  as_tibble() %>%
  filter(diff >= quant_graph_c_328350) %>%
  arrange(desc(diff)) %>%
  group_by(type.x) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

# write.csv(x = as.data.frame(imp_aa_c_328350), file = "Results/napi2b_c_328350_imp_aa.csv")

################################
# Compare control with 303-328 #
################################

compare_graphs_with_formula(`graph_list_w-plots_no_bonds`[[1]],
                            `graph_list_w-plots_303_328`[[1]],
                            system_bonds[c(3,2)])

quant_graph_c_303328 <- `comparison_graphs_list_no_bonds, 303_328`[[1]] %>%
  activate(nodes) %>%
  as_tibble() %>%
  summarise(quan_75_coef = quantile(diff, probs = 0.95)) %>%
  pull() %>%
  unname()

nodes_c_303328_1 <- `comparison_graphs_list_no_bonds, 303_328`[[1]] %>%
  activate(nodes) %>%
  as_tibble() %>%
  head(690)

nodes_c_303328_2 <- `comparison_graphs_list_no_bonds, 303_328`[[1]] %>%
  activate(nodes) %>%
  as_tibble() %>%
  tail(690)

imp_aa_c_303328 <- nodes_c_303328_1 %>%
  full_join(nodes_c_303328_2, by = c("label", "diff")) %>%
  as_tibble() %>%
  filter(diff >= quant_graph_c_328350) %>%
  arrange(desc(diff)) %>%
  group_by(type.x) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

# write.csv(x = as.data.frame(imp_aa_c_303328), file = "Results/napi2b_c_303328_imp_aa.csv")

################################
# Compare 328-350 with 303-328 #
################################

# compare_graphs_with_formula(`graph_list_w-plots_328_350`[[1]],
#                             `graph_list_w-plots_303_328`[[1]],
#                             system_bonds[c(1,2)])
# beep()

quant_graph_328350_303328 <- `comparison_graphs_list_328_350, 303_328`[[1]] %>%
  activate(nodes) %>%
  as_tibble() %>%
  summarise(quan_75_coef = quantile(diff, probs = 0.95)) %>%
  pull() %>%
  unname()

nodes_328350_303328_1 <- `comparison_graphs_list_328_350, 303_328`[[1]] %>%
  activate(nodes) %>%
  as_tibble() %>%
  head(690)

nodes_328350_303328_2 <- `comparison_graphs_list_328_350, 303_328`[[1]] %>%
  activate(nodes) %>%
  as_tibble() %>%
  tail(690)

imp_aa_328350_303328 <- nodes_328350_303328_1 %>%
  full_join(nodes_328350_303328_2, by = c("label", "diff")) %>%
  as_tibble() %>%
  filter(diff >= quant_graph_c_328350) %>%
  arrange(desc(diff)) %>%
  group_by(type.x) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

# write.csv(x = as.data.frame(imp_aa_328350_303328), file = "Results/napi2b_328350_303328_imp_aa.csv")

full_domain_diff <- imp_aa_c_328350 %>%
  full_join(imp_aa_c_303328, "type.x") %>%
  full_join(imp_aa_328350_303328, "type.x")


ggsave(filename = str_c("Plots_pre/comparison_graphs_list_nb_303_328.png", sep = "_"), plot = `comparison_graphs_list_no_bonds, 303_328`[[3]], scale = 5)
ggsave(filename = str_c("Plots_pre/comparison_graphs_list_nb_328_350.png", sep = "_"), plot = `comparison_graphs_list_no_bonds, 328_350`[[3]], scale = 5)
ggsave(filename = str_c("Plots_pre/comparison_graphs_list_328_350_303_328.png", sep = "_"), plot = `comparison_graphs_list_328_350, 303_328`[[3]], scale = 5)


