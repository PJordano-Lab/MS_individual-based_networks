Web-level metrics
================
Elena Quintero
2024-02-07

``` r
library(here)
library(tidyverse)
library(magrittr)
library(reshape2)
library(bipartite)
library(ggrepel)
library(patchwork)


theme_set(theme_bw())
```

Read all networks:

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

``` r
#Read list of names
nets_names_sp <- list.files(path = here("networks/standardized/"), pattern = "stsp_")

#Create empty list
nets_sp <- list()

#Add all read files to the list
for (i in 1:length(nets_names_sp)){
  net_file_sp <- paste0("networks/standardized/", nets_names_sp[[i]])
  net_sp <- read.csv(here(net_file_sp))
  
  if(ncol(net_sp)<2){
    net_sp <- read.csv(here(net_file_sp), sep = ";") #For those datasets saved with ; separation
  } else
  {}  

  net_sp %<>% column_to_rownames("sp")
  nets_sp[[i]] <- net_sp #add to list
  names(nets_sp)[i] <- nets_names_sp[[i]] #name the list

  # assign(nets_names_sp[[i]], net_sp)
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

``` r
nets_char_sp <- readxl::read_xlsx(here("data/selected_networks.xlsx"))

nets_ID_sp <- nets_char_sp %>% 
  filter (type == "sp") %>%
  select(net_id = code_ID,
         ref) %>%
  mutate(net_name = paste0(net_id, "_", ref)) %>%
  group_by(net_id) %>%
  mutate(net_n = cur_group_id())
```

Network level metrics:

``` r
indices <- as.list(c("number of species", "connectance", "weighted connectance", 
                     "links per species", "NODF", "weighted NODF",
                     "cluster coefficient", "weighted cluster coefficient", 
                     "generality", "linkage density", 
                     "Shannon diversity", "interaction evenness", "Alatalo interaction evenness", "H2", 
                     "ISA", "SA", "web asymmetry",
                     "niche overlap", "extinction slope", "robustness"
                     ))

# METRICS NOT USED
# "Fisher alpha", "number of compartments", "compartment diversity", "mean number of links",
# "mean number of shared partners", "degree distribution", "togetherness", "C score", "V ratio"
# "discrepancy", "vulnerability", "fc"
```

``` r
net.level.df <- data.frame()

for (i in 1:length(nets)){
  
  net.level <- networklevel(nets[[i]], indices)
  
  net.overlap2 <- networklevel(nets[[i]], "niche overlap", dist = "bray") %>% 
    t() %>% as.data.frame() %>%
    rename(niche.overlap.bray.HL = niche.overlap.HL ) %>%
    rename(niche.overlap.bray.LL = niche.overlap.LL)
  
  net.overlap3 <- networklevel(nets[[i]], "niche overlap", dist = "jaccard") %>%
    t() %>% as.data.frame() %>%
    rename(niche.overlap.jaccard.HL = niche.overlap.HL ) %>%
    rename(niche.overlap.jaccard.LL = niche.overlap.LL)
  
  links = sum(nets[[i]])
  
  net.level <- cbind(as.data.frame(t(net.level)), net.overlap2, net.overlap3) %>%
    mutate(net_size = number.of.species.HL * number.of.species.LL) %>%
    mutate(net_id = names(nets[i]),
         net_id = str_sub(net_id, 3, 7)) %>%
    mutate(links = links)
  
  net.level.df <- rbind(net.level.df, net.level)
}

net.level.df %<>% mutate(net_n = order(net_id))

glimpse(net.level.df)
```

    ## Rows: 48
    ## Columns: 36
    ## $ connectance                      <dbl> 0.3629630, 0.2093750, 0.3914286, 0.4400000, 0.3844156, 0.2846154, 0.…
    ## $ `web asymmetry`                  <dbl> -0.19402985, -0.42857143, -0.55555556, -0.55555556, -0.52173913, -0.…
    ## $ `links per species`              <dbl> 5.8507463, 2.3928571, 3.0444444, 3.4222222, 3.2173913, 1.6086957, 1.…
    ## $ `cluster coefficient`            <dbl> 0.22500000, 0.06250000, 0.25714286, 0.27142857, 0.17142857, 0.192307…
    ## $ NODF                             <dbl> 67.21158, 63.51606, 72.48350, 71.68956, 74.34611, 42.77003, 27.02206…
    ## $ `weighted NODF`                  <dbl> 52.79227, 43.93032, 37.73835, 38.51104, 41.01261, 18.95083, 16.29902…
    ## $ `interaction strength asymmetry` <dbl> 0.098096343, 0.303002612, 1.194491723, 0.604743584, 0.442346250, 0.2…
    ## $ `specialisation asymmetry`       <dbl> 1.00000000, -0.48846320, -0.82692306, 0.05225481, -0.76459563, -0.22…
    ## $ `linkage density`                <dbl> 16.046511, 9.152386, 10.990321, 13.940699, 13.173851, 3.928852, 3.22…
    ## $ `weighted connectance`           <dbl> 0.2395002, 0.1634355, 0.2442293, 0.3097933, 0.2863881, 0.1708196, 0.…
    ## $ `Shannon diversity`              <dbl> 4.9335047, 3.9591261, 3.6537360, 4.1713244, 4.2982020, 3.2323745, 3.…
    ## $ `interaction evenness`           <dbl> 0.7063286, 0.6127286, 0.6237244, 0.7120812, 0.7219933, 0.6640681, 0.…
    ## $ `Alatalo interaction evenness`   <dbl> 0.6550250, 0.5924512, 0.5490386, 0.6632521, 0.6502579, 0.7120396, 0.…
    ## $ H2                               <dbl> 0.09933567, 0.24470407, 0.17791608, 0.13305221, 0.22670003, 0.496503…
    ## $ number.of.species.HL             <dbl> 27, 16, 10, 10, 11, 10, 10, 13, 12, 11, 7, 17, 20, 9, 7, 8, 37, 3, 1…
    ## $ number.of.species.LL             <dbl> 40, 40, 35, 35, 35, 13, 14, 14, 13, 12, 11, 18, 19, 17, 15, 30, 27, …
    ## $ cluster.coefficient.HL           <dbl> 0.9360087, 0.7633714, 0.9610723, 0.9495592, 0.8898621, 0.5728938, 0.…
    ## $ cluster.coefficient.LL           <dbl> 0.4297705, 0.2429504, 0.3870170, 0.4480171, 0.4216733, 0.2866667, 0.…
    ## $ weighted.cluster.coefficient.HL  <dbl> 1.00000000, 0.75836703, 0.82419380, 0.93258022, 0.87484441, 0.380436…
    ## $ weighted.cluster.coefficient.LL  <dbl> 0.9675836, 0.8203872, 0.9500581, 0.9613604, 0.9671022, 0.5767548, 0.…
    ## $ niche.overlap.HL                 <dbl> 0.28227633, 0.09229992, 0.14771050, 0.22153444, 0.18881442, 0.152051…
    ## $ niche.overlap.LL                 <dbl> 0.8377570, 0.6984493, 0.9032446, 0.9060923, 0.8430284, 0.4520093, 0.…
    ## $ extinction.slope.HL              <dbl> 10.7307404, 3.7217600, 4.2930684, 5.5611144, 4.8479593, 3.1093811, 2…
    ## $ extinction.slope.LL              <dbl> 7.8229402, 2.5961339, 5.6495242, 6.7967055, 7.1907273, 2.4643146, 2.…
    ## $ robustness.HL                    <dbl> 0.9145234, 0.7824115, 0.8078489, 0.8408830, 0.8250309, 0.7547800, 0.…
    ## $ robustness.LL                    <dbl> 0.8599260, 0.7138990, 0.8181298, 0.8268341, 0.8678092, 0.7013677, 0.…
    ## $ generality.HL                    <dbl> 27.560125, 15.435135, 20.212817, 25.501566, 23.612408, 5.218376, 4.0…
    ## $ vulnerability.LL                 <dbl> 4.532896, 2.869637, 1.767824, 2.379832, 2.735294, 2.639328, 2.420494…
    ## $ niche.overlap.bray.HL            <dbl> 0.09729529, 0.03883194, 0.07298345, 0.09792473, 0.11689444, 0.127150…
    ## $ niche.overlap.bray.LL            <dbl> 0.5567506, 0.3390878, 0.4458126, 0.5419888, 0.5085545, 0.2842427, 0.…
    ## $ niche.overlap.jaccard.HL         <dbl> 0.05728874, 0.02443975, 0.04183463, 0.05846665, 0.06721856, 0.085077…
    ## $ niche.overlap.jaccard.LL         <dbl> 0.40377235, 0.22904359, 0.31947895, 0.40320436, 0.36587822, 0.179518…
    ## $ net_size                         <dbl> 1080, 640, 350, 350, 385, 130, 140, 182, 156, 132, 77, 306, 380, 153…
    ## $ net_id                           <chr> "01_01", "01_02", "02_01", "02_02", "02_03", "03_01", "03_02", "03_0…
    ## $ links                            <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
    ## $ net_n                            <int> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 2…

``` r
write_csv(net.level.df, here("data/net.level.df.csv"))
```

``` r
net.level.df.sp <- data.frame()

for (i in 1:length(nets_sp)){
  
  net.level.sp <- networklevel(nets_sp[[i]], indices)
  
  net.overlap2 <- networklevel(nets_sp[[i]], "niche overlap", dist = "bray") %>% 
    t() %>% as.data.frame() %>%
    rename(niche.overlap.bray.HL = niche.overlap.HL ) %>%
    rename(niche.overlap.bray.LL = niche.overlap.LL)
  
  net.overlap3 <- networklevel(nets_sp[[i]], "niche overlap", dist = "jaccard") %>%
    t() %>% as.data.frame() %>%
    rename(niche.overlap.jaccard.HL = niche.overlap.HL ) %>%
    rename(niche.overlap.jaccard.LL = niche.overlap.LL)
  
  links = sum(nets_sp[[i]])
  
  net.level.sp <- cbind(as.data.frame(t(net.level.sp)), net.overlap2, net.overlap3)  %>%
    mutate(net_size = number.of.species.HL * number.of.species.LL) %>%
    mutate(net_id = names(nets_sp[i]),
         net_id = str_sub(net_id, 1, 10)) %>%
    mutate(links = links)
  
  net.level.df.sp <- rbind(net.level.df.sp, net.level.sp)
}

net.level.df.sp %<>% mutate(net_n = order(net_id))

glimpse(net.level.df.sp)
```

    ## Rows: 68
    ## Columns: 36
    ## $ connectance                      <dbl> 0.2533333, 0.4411765, 0.2731092, 0.5333333, 0.4000000, 0.4333333, 0.…
    ## $ `web asymmetry`                  <dbl> 0.18032787, 0.03030303, 0.24444444, 0.25000000, 0.09090909, 0.090909…
    ## $ `links per species`              <dbl> 3.737705, 3.636364, 2.888889, 1.000000, 1.090909, 1.181818, 1.928571…
    ## $ `cluster coefficient`            <dbl> 0.20000000, 0.31250000, 0.23529412, 0.33333333, 0.30000000, 0.400000…
    ## $ NODF                             <dbl> 65.78425, 79.13246, 63.17203, 57.69231, 72.00000, 74.66667, 93.02326…
    ## $ `weighted NODF`                  <dbl> 44.06047, 49.74284, 43.96037, 46.15385, 62.00000, 69.33333, 66.86047…
    ## $ `interaction strength asymmetry` <dbl> -0.102223525, -0.104346451, -0.050103350, -5.707375801, -8.985732463…
    ## $ `specialisation asymmetry`       <dbl> 0.94419400, 0.29044024, -0.92564302, 0.07279023, -0.20313124, -0.357…
    ## $ `linkage density`                <dbl> 6.956794, 4.669495, 5.398074, 1.719472, 1.933557, 2.416312, 2.926579…
    ## $ `weighted connectance`           <dbl> 0.11404581, 0.14149985, 0.11995719, 0.21493401, 0.17577793, 0.219664…
    ## $ `Shannon diversity`              <dbl> 3.930055, 3.199001, 3.823299, 1.053506, 1.181889, 1.903583, 2.483697…
    ## $ `interaction evenness`           <dbl> 0.5777458, 0.5706589, 0.6201200, 0.3890277, 0.3474921, 0.5596802, 0.…
    ## $ `Alatalo interaction evenness`   <dbl> 0.5009230, 0.3861832, 0.5661735, 0.5253060, 0.5974077, 0.6089427, 0.…
    ## $ H2                               <dbl> 0.2552644, 0.1663608, 0.2895236, 0.4166363, 0.1721675, 0.1994522, 0.…
    ## $ number.of.species.HL             <dbl> 36, 17, 28, 5, 6, 6, 6, 9, 6, 10, 7, 29, 27, 14, 8, 7, 14, 21, 30, 3…
    ## $ number.of.species.LL             <dbl> 25, 16, 17, 3, 5, 5, 8, 31, 7, 16, 12, 35, 5, 11, 15, 71, 50, 7, 8, …
    ## $ cluster.coefficient.HL           <dbl> 0.4609968, 0.7811052, 0.5530137, 0.8613333, 0.5002278, 0.6767442, 0.…
    ## $ cluster.coefficient.LL           <dbl> 0.6337912, 0.8263240, 0.5193060, 0.7776000, 0.9745634, 0.7790698, 0.…
    ## $ weighted.cluster.coefficient.HL  <dbl> 0.92170934, 0.98287291, 0.88756451, 0.00000000, 0.00000000, 0.000000…
    ## $ weighted.cluster.coefficient.LL  <dbl> 0.84405594, 0.92836212, 0.78648457, 0.00000000, 0.00000000, 0.000000…
    ## $ niche.overlap.HL                 <dbl> 0.4865477, 0.8117259, 0.4264740, 0.5900290, 0.9973706, 0.6646970, 0.…
    ## $ niche.overlap.LL                 <dbl> 0.2791500, 0.5524205, 0.3208485, 0.4279841, 0.5465940, 0.5163291, 0.…
    ## $ extinction.slope.HL              <dbl> 5.539376, 7.508850, 5.382114, 3.379500, 1.583034, 2.008398, 3.697882…
    ## $ extinction.slope.LL              <dbl> 4.489972, 6.245220, 3.309007, 1.628800, 1.418173, 2.238060, 4.691325…
    ## $ robustness.HL                    <dbl> 0.8199792, 0.8754105, 0.8254271, 0.7696667, 0.6106077, 0.6682441, 0.…
    ## $ robustness.LL                    <dbl> 0.8002365, 0.8574808, 0.7556634, 0.5865000, 0.5991066, 0.6810434, 0.…
    ## $ generality.HL                    <dbl> 3.752654, 3.755649, 4.924750, 1.199885, 1.222159, 2.222414, 3.172118…
    ## $ vulnerability.LL                 <dbl> 10.160935, 5.583341, 5.871397, 2.239059, 2.644955, 2.610211, 2.68104…
    ## $ niche.overlap.bray.HL            <dbl> 0.18347485, 0.32788005, 0.17002209, 0.18106342, 0.25084925, 0.259552…
    ## $ niche.overlap.bray.LL            <dbl> 0.12607116, 0.24763884, 0.15276247, 0.12789392, 0.27078096, 0.199363…
    ## $ niche.overlap.jaccard.HL         <dbl> 0.11586431, 0.21762683, 0.10950731, 0.11822761, 0.17878896, 0.172641…
    ## $ niche.overlap.jaccard.LL         <dbl> 0.07820296, 0.15759876, 0.09188458, 0.07248097, 0.20690476, 0.136471…
    ## $ net_size                         <dbl> 900, 272, 476, 15, 30, 30, 48, 279, 42, 160, 84, 1015, 135, 154, 120…
    ## $ net_id                           <chr> "stsp_01_01", "stsp_01_02", "stsp_02_01", "stsp_03_01", "stsp_03_02"…
    ## $ links                            <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
    ## $ net_n                            <int> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 2…

``` r
write_csv(net.level.df.sp, here("data/net.level.df.sp.csv"))
```

Assortativity & centrality:

``` r
library(igraph)

metrics.igraph.ind <- data.frame()

for (i in seq_along(nets)){
  
  net <- as.data.frame(nets[[i]])

  inet <- graph_from_incidence_matrix(net, weighted= T,  add.names=NULL)
  
  assortativity = assortativity_degree(inet)
  centr_binary = centr_eigen(inet, scale = TRUE, normalized = TRUE)$centralization
  centr.weighted.obs <- sum(max(eigen_centrality(inet)$vector) - eigen_centrality(inet)$vector)
  centralization.w <- centr.weighted.obs/ centr_eigen(inet)$theoretical_max
  net_id = str_sub(names(nets[i]), 3, 7)
  as_calc <- cbind(net_id, assortativity, centr_binary, centralization.w)
  
  metrics.igraph.ind <- rbind(metrics.igraph.ind, as_calc)
}
```

``` r
metrics.igraph.sp <- data.frame()

for (i in seq_along(nets_sp)){
  
  net <- as.data.frame(nets_sp[[i]])

  inet <- graph_from_incidence_matrix(net, weighted=T,  add.names=NULL)
  
  assortativity = assortativity_degree(inet)
  centr_binary = centr_eigen(inet, scale = TRUE, normalized = TRUE)$centralization
  centr.weighted.obs <- sum(max(eigen_centrality(inet)$vector) - eigen_centrality(inet)$vector)
  centralization.w <- centr.weighted.obs/ centr_eigen(inet)$theoretical_max
  net_id = str_sub(names(nets_sp[i]), 6, 10)
  as_calc <- cbind(net_id, assortativity, centr_binary, centralization.w)
  
  metrics.igraph.sp <- rbind(metrics.igraph.sp, as_calc)
}
```

``` r
metrics.igraph <- metrics.igraph.ind %>%
  mutate(type = "ind") %>%
  full_join(metrics.igraph.sp) %>%
  mutate(type = ifelse(is.na(type), "sp", type)) %>%
  mutate(assortativity = as.numeric(assortativity),
         centr_binary = as.numeric(centr_binary),
         centralization.w = as.numeric(centralization.w))

glimpse(metrics.igraph)
```

    ## Rows: 116
    ## Columns: 5
    ## $ net_id           <chr> "01_01", "01_02", "02_01", "02_02", "02_03", "03_01", "03_02", "03_03", "03_04", "03…
    ## $ assortativity    <dbl> -0.5071021, -0.6033629, -0.7374228, -0.7297445, -0.6540469, -0.4644660, -0.2413978, …
    ## $ centr_binary     <dbl> 0.6210147122, 0.7887528765, 0.7122377610, 0.6821950187, 0.7050468979, 0.6976499423, …
    ## $ centralization.w <dbl> 0.8549594, 0.9221565, 0.9284492, 0.8970345, 0.8961652, 0.8593742, 0.9149339, 0.82652…
    ## $ type             <chr> "ind", "ind", "ind", "ind", "ind", "ind", "ind", "ind", "ind", "ind", "ind", "ind", …

``` r
write.csv(metrics.igraph, here("data/igraph_metrics.csv"))
```
