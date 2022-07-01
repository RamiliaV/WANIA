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
library(readxl)

muts <- read_excel("Results_pre/napi2b_diff_3sys.xlsx", sheet=4)
muts <- muts %>%
  select(1,2,5)

nodes_1 <- `graph_list_w-plots_no_bonds`[[1]] %>%
  activate(nodes) %>%
  as_tibble()
edges_1 <- `graph_list_w-plots_no_bonds`[[1]] %>%
  activate(edges) %>%
  as_tibble() 

nodes_2 <- `graph_list_w-plots_303_328`[[1]] %>%
  activate(nodes) %>%
  mutate(id = id + nrow(nodes_1))%>%
  as_tibble() 
edges_2 <- `graph_list_w-plots_303_328`[[1]] %>%
  activate(edges) %>% 
  as_tibble() %>%
  mutate(from = from + nrow(nodes_1), to = to + nrow(nodes_1))

nodes_3 <- `graph_list_w-plots_328_350`[[1]] %>%
  activate(nodes) %>%
  mutate(id = id + nrow(nodes_1)*2)%>%
  as_tibble() 
edges_3 <- `graph_list_w-plots_328_350`[[1]] %>%
  activate(edges) %>%
  as_tibble() %>%
  mutate(from = from + nrow(nodes_1)*2, to = to + nrow(nodes_1)*2)

nodes_d3 <- nodes_1 %>%
  bind_rows(nodes_2, nodes_3, .id = "graph") %>%
  mutate(graph = ifelse(graph == 1, "no bonds", 
                        ifelse(graph == 2,"303-328", "328-350")))

edges_d3 <- edges_1 %>%
  bind_rows(edges_2, edges_3, .id = "graph") %>%
  mutate(graph = ifelse(graph == 1, "no bonds", 
                        ifelse(graph == 2,"303-328", "328-350")))

nodes_d3 <- mutate(nodes_d3, id = id - 1)
edges_d3 <- mutate(edges_d3, from = from - 1, to = to - 1)

nodes_d3 <- nodes_d3 %>%
  mutate(label_num = as.numeric(str_extract(label, "(\\d)+"))) %>%
  left_join(muts, by = c("label_num"="Aanumber")) %>%
  mutate(`Cancer type` = ifelse(!is.na(`Cancer type`), `Cancer type`, "No mutation")) %>%
  mutate(Mutation = ifelse(!is.na(Mutation), Mutation, "No mutation")) %>%
  mutate(Domain = ifelse(Mutation=="No mutation", NA, str_c(type)))

nodes_d3 <- nodes_d3 %>%
  mutate(title = paste0("<p><b>", Mutation,"</b><br/>", type, "<b><br/>", `Cancer type`, "<b><br/>", round(coef, 3), "</b></p>")) %>%
  mutate(label = str_c(str_c(label, round(coef, 3), sep=" - "), graph, sep = "\n"))

names(nodes_d3)[1] <- "group"
names(nodes_d3)[15] <- "value"

# visNetwork(nodes_d3, edges_d3, height = "1200px", width = "100%") %>%
#   visIgraphLayout() %>%
#   visOptions(highlightNearest = TRUE,  selectedBy = "Domain") %>%
#   visNodes(font = list(size = 50, strokeColor = "white", strokeWidth=1)) %>%
#   visGroups(groupname = "328_350", color ="red") %>%
#   visGroups(groupname = "303_328", color ="blue") %>%
#   visLegend(position = "right", main = "") %>%
#   visSave(file = "Results/interactive_network_mutation_by_domain.html", background = "white")

color <- c("#43c29e",
                   "#cb3f7d",
                   "#5bbc73",
                   "#be50a8",
                   "#9bb03f",
                   "#6f71d9",
                   "#d49d3b",
                   "#553687",
                   "#698d3c",
                   "#c282d3",
                   "#a67e3a",
                   "#5e87d3",
                   "#cc6837",
                   "#36dee6",
                   "#d84d56",
                   "#842b60",
                   "#d67663",
                   "#d7719d",
                   "#893018",
                   "#d45e75",
                   "#902639")
domain_color <- data.frame(type = unique(nodes_d3$type), color = color)
quan_75 <- quantile(nodes_d3$value, probs = 0.95)
nodes_d3 <- nodes_d3 %>%
  mutate(coef_75 = ifelse(value>quan_75, "Yes", NA)) %>%
  # mutate(Domain = ifelse(coef_75 == "Yes"&Mutation!="No mutation", type, NA))
  mutate(Domain = ifelse(coef_75 == "Yes", "group", NA)) 
  # %>%
  # left_join(domain_color)
  
visNetwork(nodes_d3, edges_d3, height = "600px", width = "700px") %>%
  visIgraphLayout() %>%
  visOptions(highlightNearest = TRUE,  selectedBy = "Domain") %>%
  visNodes(font = list(size = 80, strokeColor = "white", strokeWidth=1)) %>%
  visLegend(position = "right", main = "") %>%
  visGroups(groupname = "328_350", color ="red") %>%
  visGroups(groupname = "303_328", color ="blue") %>%
  visSave(file = "Results/interactive_network_mutation_upper25_3.html", background = "white")

imp_aa <- nodes_d3 %>% 
  filter(!is.na(Domain)) %>%
  select(group, label, type, Mutation, `Cancer type`) %>%
  separate(label, into = c("label", "score"), sep = " - ")

write_csv2(imp_aa, "Results_pre/imp_aa_3sys.csv")

imp_aa_2 <- imp_aa %>% 
  group_by(type, group) %>%
  summarise(aa = str_c(label, collapse = ", ")) %>%
  pivot_wider(names_from = group, values_from = aa)
write_csv2(imp_aa_2, "Results_pre/imp_aa_3sys_all.csv")

imp_aa_3 <- imp_aa %>% 
  group_by(label) %>%
  summarise(n = n())
write_csv2(imp_aa_3, "Results_pre/imp_aa_3sys_list.csv")

write.csv2(edges_1, "edges_nb.csv")
write.csv2(edges_2, "edges_303.csv")
write.csv2(edges_3, "edges_328.csv")
