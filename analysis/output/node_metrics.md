node-level metrics
================
Elena Quintero
2024-02-07

``` r
library(here)
library(tidyverse)
library(magrittr)
library(reshape2)
library(ggrepel)
library(patchwork)

theme_set(theme_bw())
```

Read networks:

``` r
#Read list of names
nets_names <- list.files(path = here("networks/standardized/"), pattern = "_int")

#Create empty list
nets <- list()

#Add all read files to the list
for (i in 1:length(nets_names)){
  net_file <- paste0("networks/standardized/", nets_names[[i]])
  net <- read.csv(here(net_file))
  
  if(ncol(net)<2){
    net <- read.csv(here(net_file), sep = ";") #For those datasets saved with ; separation
  } else
  {}  

  net %<>% column_to_rownames("ind")
  nets[[i]] <- net #add to list
  names(nets)[i] <- nets_names[[i]] #name the list

  # assign(nets_names[[i]], net)
}
```

Read nets characteristics (name, site etc…)

``` r
nets_char <- readxl::read_xlsx(here("data/selected_individual_networks.xlsx"))

nets_ID <- nets_char %>%
  select(net_id = code_ID,
         ref,
         plant_sp = `focus species`,
         pop = pop_site) %>%
  mutate(net_name = paste0(plant_sp, "_", pop)) %>%
  group_by(net_id) %>%
  mutate(net_n = cur_group_id())
```

Individual level metrics:

``` r
library(bipartite)

bipartite.ind.level <- data.frame()

for (i in 1:length(nets)){

  ind.level <- specieslevel(nets[[i]], level = "lower")

  ind.level %<>% rownames_to_column("ind") %>%
  mutate(net_id = names(nets[i]),
         net_id = str_sub(net_id, 3, 7))
  
  bipartite.ind.level <- rbind(bipartite.ind.level, ind.level)
}

bipartite.ind.level %<>% group_by(net_id) %>% 
  mutate(net_n = cur_group_id()) %>% 
  mutate(ind_ID = paste0(net_id, "_", ind))
  
glimpse(bipartite.ind.level)
```

    ## Rows: 1,022
    ## Columns: 24
    ## Groups: net_id [48]
    ## $ ind                           <chr> "301", "302", "303", "304", "305", "30…
    ## $ degree                        <dbl> 11, 15, 9, 19, 5, 11, 7, 7, 9, 7, 8, 7…
    ## $ normalised.degree             <dbl> 0.4074074, 0.5555556, 0.3333333, 0.703…
    ## $ species.strength              <dbl> 0.85051099, 2.68245083, 0.13992168, 3.…
    ## $ interaction.push.pull         <dbl> -0.013589910, 0.112163389, -0.09556425…
    ## $ nestedrank                    <dbl> 0.20512821, 0.07692308, 0.58974359, 0.…
    ## $ PDI                           <dbl> 0.9523748, 0.9401682, 0.9577912, 0.927…
    ## $ species.specificity.index     <dbl> 0.5295019, 0.4518455, 0.5550037, 0.436…
    ## $ resource.range                <dbl> 0.6153846, 0.4615385, 0.6923077, 0.307…
    ## $ PSI                           <dbl> 0.85051099, 2.68245083, 0.13992168, 3.…
    ## $ node.specialisation.index.NSI <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
    ## $ betweenness                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
    ## $ weighted.betweenness          <dbl> 0.00000000, 0.00000000, 0.00000000, 0.…
    ## $ closeness                     <dbl> 0.025, 0.025, 0.025, 0.025, 0.025, 0.0…
    ## $ weighted.closeness            <dbl> 0.031765756, 0.030524084, 0.014375754,…
    ## $ Fisher.alpha                  <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ partner.diversity             <dbl> 1.429485, 1.753430, 1.322201, 1.784777…
    ## $ effective.partners            <dbl> 4.176546, 5.774377, 3.751668, 5.958250…
    ## $ proportional.similarity       <dbl> 0.7557782, 0.7745749, 0.7720482, 0.801…
    ## $ proportional.generality       <dbl> 0.7998651, 1.1058714, 0.7184952, 1.141…
    ## $ d                             <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
    ## $ net_id                        <chr> "01_01", "01_01", "01_01", "01_01", "0…
    ## $ net_n                         <int> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
    ## $ ind_ID                        <chr> "01_01_301", "01_01_302", "01_01_303",…

Mean niche overlap:

Using bray curtis dissmilarity index

``` r
niche_overlap <- data.frame()

for (i in 1:length(nets)) {
  
  net <- as.data.frame(nets[[i]])
  net_id = names(nets[i]) %>% str_sub(3, 7)
  
  bray.overlap <- as.matrix(vegdist(net, method="bray")) %>% 
    as.data.frame() %>%
    rownames_to_column("ind") %>%
    summarise(ind = ind,
              mean.bray.overlap = rowMeans(select(., -1))) %>%
    mutate(net_id = net_id)
  
  jaccard.overlap <- as.matrix(vegdist(net, method="jaccard")) %>% 
    as.data.frame() %>%
    rownames_to_column("ind") %>%
    summarise(ind = ind,
              mean.jaccard.overlap = rowMeans(select(., -1))) %>%
    mutate(net_id = net_id)
  
  overlap <- left_join(bray.overlap, jaccard.overlap)
  
  niche_overlap <- rbind(niche_overlap, overlap)
}

niche_overlap %<>% mutate(ind_ID = paste0(net_id, "_", ind)) %>%
  mutate(mean.bray.overlap = 1 - mean.bray.overlap,
         mean.jaccard.overlap = 1 - mean.jaccard.overlap) #inverse, where 1 is full overlap, 0 no overlap

glimpse(niche_overlap)
```

    ## Rows: 1,022
    ## Columns: 5
    ## $ ind                  <chr> "301", "302", "303", "304", "305", "306", "307"…
    ## $ mean.bray.overlap    <dbl> 0.5754007, 0.5782921, 0.5968529, 0.4159960, 0.5…
    ## $ net_id               <chr> "01_01", "01_01", "01_01", "01_01", "01_01", "0…
    ## $ mean.jaccard.overlap <dbl> 0.4182767, 0.4215563, 0.4405808, 0.2837609, 0.3…
    ## $ ind_ID               <chr> "01_01_301", "01_01_302", "01_01_303", "01_01_3…

Eigenvalue centrality - igraph:

``` r
library(igraph)

metrics.igraph.ind <- data.frame()

for (i in seq_along(nets)){
  
  net <- as.data.frame(nets[[i]])

  inet <- graph_from_incidence_matrix(net, weighted= T,  add.names=NULL)
  
  eigen.centrality = eigen_centrality(inet)$vector
  net_id = str_sub(names(nets[i]), 3, 7)
  as_calc <- cbind(net_id, eigen.centrality) %>%
    as.data.frame() %>%
    rownames_to_column("ind") %>%
    mutate(ind_ID = paste0(net_id, "_", ind))
  
  metrics.igraph.ind <- rbind(metrics.igraph.ind, as_calc) 
}

metrics.igraph.ind %<>% 
  mutate(eigen.centrality = as.numeric(eigen.centrality))

glimpse(metrics.igraph.ind)
```

    ## Rows: 1,601
    ## Columns: 4
    ## $ ind              <chr> "301", "302", "303", "304", "305", "306", "307", "3…
    ## $ net_id           <chr> "01_01", "01_01", "01_01", "01_01", "01_01", "01_01…
    ## $ eigen.centrality <dbl> 0.24815105, 0.22843692, 0.11912347, 0.51744591, 0.0…
    ## $ ind_ID           <chr> "01_01_301", "01_01_302", "01_01_303", "01_01_304",…

Closeness and betweenness - tnet:

``` r
library(tnet)

mynets = nets[-c(18, 30, 31)] # remove nets "08_01", "12_06", "12_07"

metrics.tnet.ind <- data.frame()

for (i in seq_along(mynets)){
  
  net <- as.data.frame(mynets[[i]])
  
  inet <- graph_from_incidence_matrix(net, weighted= T,  add.names=NULL)
  NodeLabels <- V(inet)$name
  NodeType <- V(inet)$type
  nodes_id <- cbind(NodeLabels, NodeType) %>% as.data.frame() %>% 
    filter(NodeType == "FALSE") %>%
    select(ind = NodeLabels) %>%
    rownames_to_column("node")

  #w.tm <- tnet::as.tnet(net) #problem with Paco's net, so I do it manually
  tm <- get.edgelist(inet, names=FALSE)
  w <- E(inet)$weight
  w.tm <- cbind(tm, w)
  
  om.tm <- projecting_tm(tm, "Newman") #One mode projection
  om.w.tm <- projecting_tm(w.tm, "Newman") #One mode projection
  
  closeness <- tnet::closeness_w(om.tm) %>% as.data.frame() %>% 
    select(node, tm.closeness = closeness)
  w.closeness <- tnet::closeness_w(om.w.tm) %>% as.data.frame()  %>% 
    select(node, tm.w.closeness = closeness)
  betweenness <- tnet::betweenness_w(om.tm) %>% as.data.frame() %>% 
    select(node, tm.betweenness = betweenness)
  w.betweenness <- tnet::betweenness_w(om.w.tm) %>% as.data.frame() %>% 
    select(node, tm.w.betweenness = betweenness)
  w.mets <- full_join(closeness, w.closeness) %>%  # don't do cbind because in some nets closeness is not calculated for all nodes (bec of singletons)
    full_join(betweenness) %>% full_join(w.betweenness) %>% 
    mutate(node = as.character(node)) %>%
    left_join(nodes_id, by = "node") %>%
    mutate(tm.closeness = ifelse(is.na(tm.closeness), 0, tm.closeness),
           tm.w.closeness = ifelse(is.na(tm.w.closeness), 0, tm.w.closeness)) # assign zero values to nodes where no closeness.

  net_id = str_sub(names(mynets[i]), 3, 7)
  
  tnet_metrics <- cbind(net_id, w.mets[, -1])
  
  metrics.tnet.ind <- rbind(metrics.tnet.ind, tnet_metrics)
}

metrics.tnet.ind %<>% mutate(ind_ID = paste0(net_id, "_", ind))

glimpse(metrics.tnet.ind)
```

    ## Rows: 1,003
    ## Columns: 7
    ## $ net_id           <chr> "01_01", "01_01", "01_01", "01_01", "01_01", "01_01…
    ## $ tm.closeness     <dbl> 0.03039829, 0.03329066, 0.02345089, 0.03965148, 0.0…
    ## $ tm.w.closeness   <dbl> 0.031202628, 0.031920244, 0.014090767, 0.074394795,…
    ## $ tm.betweenness   <dbl> 0.0, 0.5, 0.0, 149.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,…
    ## $ tm.w.betweenness <dbl> 0, 0, 0, 89, 0, 0, 0, 0, 0, 0, 0, 0, 0, 74, 0, 0, 0…
    ## $ ind              <chr> "301", "302", "303", "304", "305", "306", "307", "3…
    ## $ ind_ID           <chr> "01_01_301", "01_01_302", "01_01_303", "01_01_304",…

``` r
ind.level.df <- bipartite.ind.level %<>% 
  left_join(niche_overlap, by = c("ind", "net_id", "ind_ID")) %>%
  left_join(metrics.igraph.ind, by = c("ind", "net_id", "ind_ID")) %>%
  left_join(metrics.tnet.ind, by = c("ind", "net_id", "ind_ID"))

glimpse(ind.level.df)
```

    ## Rows: 1,022
    ## Columns: 31
    ## Groups: net_id [48]
    ## $ ind                           <chr> "301", "302", "303", "304", "305", "30…
    ## $ degree                        <dbl> 11, 15, 9, 19, 5, 11, 7, 7, 9, 7, 8, 7…
    ## $ normalised.degree             <dbl> 0.4074074, 0.5555556, 0.3333333, 0.703…
    ## $ species.strength              <dbl> 0.85051099, 2.68245083, 0.13992168, 3.…
    ## $ interaction.push.pull         <dbl> -0.013589910, 0.112163389, -0.09556425…
    ## $ nestedrank                    <dbl> 0.20512821, 0.07692308, 0.58974359, 0.…
    ## $ PDI                           <dbl> 0.9523748, 0.9401682, 0.9577912, 0.927…
    ## $ species.specificity.index     <dbl> 0.5295019, 0.4518455, 0.5550037, 0.436…
    ## $ resource.range                <dbl> 0.6153846, 0.4615385, 0.6923077, 0.307…
    ## $ PSI                           <dbl> 0.85051099, 2.68245083, 0.13992168, 3.…
    ## $ node.specialisation.index.NSI <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
    ## $ betweenness                   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
    ## $ weighted.betweenness          <dbl> 0.00000000, 0.00000000, 0.00000000, 0.…
    ## $ closeness                     <dbl> 0.025, 0.025, 0.025, 0.025, 0.025, 0.0…
    ## $ weighted.closeness            <dbl> 0.031765756, 0.030524084, 0.014375754,…
    ## $ Fisher.alpha                  <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ partner.diversity             <dbl> 1.429485, 1.753430, 1.322201, 1.784777…
    ## $ effective.partners            <dbl> 4.176546, 5.774377, 3.751668, 5.958250…
    ## $ proportional.similarity       <dbl> 0.7557782, 0.7745749, 0.7720482, 0.801…
    ## $ proportional.generality       <dbl> 0.7998651, 1.1058714, 0.7184952, 1.141…
    ## $ d                             <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
    ## $ net_id                        <chr> "01_01", "01_01", "01_01", "01_01", "0…
    ## $ net_n                         <int> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
    ## $ ind_ID                        <chr> "01_01_301", "01_01_302", "01_01_303",…
    ## $ mean.bray.overlap             <dbl> 0.5754007, 0.5782921, 0.5968529, 0.415…
    ## $ mean.jaccard.overlap          <dbl> 0.4182767, 0.4215563, 0.4405808, 0.283…
    ## $ eigen.centrality              <dbl> 0.24815105, 0.22843692, 0.11912347, 0.…
    ## $ tm.closeness                  <dbl> 0.03039829, 0.03329066, 0.02345089, 0.…
    ## $ tm.w.closeness                <dbl> 0.031202628, 0.031920244, 0.014090767,…
    ## $ tm.betweenness                <dbl> 0.0, 0.5, 0.0, 149.5, 0.0, 0.0, 0.0, 0…
    ## $ tm.w.betweenness              <dbl> 0, 0, 0, 89, 0, 0, 0, 0, 0, 0, 0, 0, 0…

``` r
write_csv(ind.level.df, here("data/ind.level.df.csv"))
```
