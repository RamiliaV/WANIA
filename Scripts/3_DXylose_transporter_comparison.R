#############################
# IMPORT AND GRAPH CREATION #
#############################

# Graph 1
nodes <- read_csv("Data/charmm-gui-3136192901.csv") %>%
  dplyr::mutate(resid = str_c(resno, resid, sep = "_"))
names(nodes)[1:2] <- c('id', 'label')
edges <- read_csv("Data/charmm-gui-3136192901 (1).csv") %>%
  dplyr::filter(persent_ring != 0) %>%
  dplyr::mutate(five_aa = ifelse(abs(number.x - number.y) >= 5, 2, 1))
names(edges) <- c('to', 'from', 'weight', 'five_aa')

parameters <- c("default", names(nodes)[c(8:10,12:15)])

get_networks(nodes, edges, parameters)
get_network_parameters(graph_list)
compare_node_parameters(parameter_list)
stating_parameter_formula(comparison_node_parameters)
graph_with_formula(nodes, edges, formula_coef, names[1])

# Graph 2
nodes <- read_csv("Data/charmm-gui-3136217149.csv") %>%
  dplyr::mutate(resid = str_c(resno, resid, sep = "_"))
names(nodes)[1:2] <- c('id', 'label')
edges <- read_csv("Data/charmm-gui-3136217149 (1).csv") %>%
  dplyr::filter(persent_ring != 0) %>%
  dplyr::mutate(five_aa = ifelse(abs(number.x - number.y) >= 5, 2, 1))
names(edges) <- c('to', 'from', 'weight', 'five_aa')

parameters <- c("default", names(nodes)[c(8:10,12:15)])

get_networks(nodes, edges, parameters)
get_network_parameters(graph_list)
compare_node_parameters(parameter_list)
stating_parameter_formula(comparison_node_parameters)
graph_with_formula(nodes, edges, formula_coef, names[2])

compare_graphs_with_formula(`graph_list_w-plots_3136192901`[[1]],
                            `graph_list_w-plots_3136217149`[[1]],
                            names)
beep()