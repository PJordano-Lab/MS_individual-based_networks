Plot networks
================
Elena Quintero
2024-07-10

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

### Read all networks:

``` r
#Read list of names
nets_names <- list.files(path = here("networks/nets_std/"))

#Create empty list
nets <- list()

#Add all read files to the list
for (i in 1:length(nets_names)) {
  
  net_file <- paste0("networks/nets_std/", nets_names[[i]])
  net <- read.csv(here(net_file))
  
  if (str_detect(nets_names[i], "sp_")) {
    net %<>% column_to_rownames("sp")
  } else
  {
    net %<>% column_to_rownames("ind")
    }  

  nets[[i]] <- net #add to list
  names(nets)[i] <- nets_names[[i]] #name the list
  
}
```

### Read nets characteristics

(name, site etc…)

``` r
nets_char_ind <- readxl::read_xlsx(here("data/selected_individual_networks.xlsx"))

nets_ID_ind <- nets_char_ind %>%
  select(net_id = code_ID,
         ref,
         plant_sp = `focus species`,
         pop = pop_site) %>%
  filter(! net_id %in%  c("08_01", "12_01", "12_06", "12_07")) %>%
  mutate(type = "ind-based",
         net_name = paste0(plant_sp, "_", pop))

nets_char_sp <- readxl::read_xlsx(here("data/selected_networks.xlsx"))

nets_ID_sp <- nets_char_sp %>% 
  filter (type == "sp") %>%
  select(net_id = code_ID,
         ref, site) %>%
  filter(! net_id %in% c("03_01", "03_02", "03_03", "03_04",
                       "05_01", "09_01", "29_12", "31_01", "44_01")) %>%
  mutate(type = "sp-based",
         net_name = paste0(ref, "_", site))

nets_ID <- full_join(nets_ID_ind, nets_ID_sp) %>%
  group_by(type, net_id) %>%
  mutate(net_n = cur_group_id(),
         net_name = paste0("Net: ", net_n, " - ", net_name))
```

### Plot webs:

``` r
for (i in 1:length(nets)){
  
  if (str_detect(names(nets)[i], "sp_")) {
    this_net <- str_sub(names(nets)[i], 4, 8)
    net_info <- nets_ID %>% 
      filter(type == "sp-based", net_id == this_net)
  }
  else{
    this_net <- str_sub(names(nets)[i], 1, 5)
    net_info <- nets_ID %>% 
      filter(type == "ind-based", net_id == this_net)
  }
  
  col.animals = MetBrewer::met.brewer("Lakota", length(nets[[i]]))
  
  plotweb(nets[[i]],
        col.high = col.animals, bor.col.high = col.animals,
        col.low = "grey50", bor.col.low = "grey50",
        col.interaction = col.animals, 
        bor.col.interaction = col.animals,
        text.rot = 90)
  
  title(main = net_info$net_name)
}
```

![](net_plots_files/figure-gfm/network-1.png)<!-- -->![](net_plots_files/figure-gfm/network-2.png)<!-- -->![](net_plots_files/figure-gfm/network-3.png)<!-- -->![](net_plots_files/figure-gfm/network-4.png)<!-- -->![](net_plots_files/figure-gfm/network-5.png)<!-- -->![](net_plots_files/figure-gfm/network-6.png)<!-- -->![](net_plots_files/figure-gfm/network-7.png)<!-- -->![](net_plots_files/figure-gfm/network-8.png)<!-- -->![](net_plots_files/figure-gfm/network-9.png)<!-- -->![](net_plots_files/figure-gfm/network-10.png)<!-- -->![](net_plots_files/figure-gfm/network-11.png)<!-- -->![](net_plots_files/figure-gfm/network-12.png)<!-- -->![](net_plots_files/figure-gfm/network-13.png)<!-- -->![](net_plots_files/figure-gfm/network-14.png)<!-- -->![](net_plots_files/figure-gfm/network-15.png)<!-- -->![](net_plots_files/figure-gfm/network-16.png)<!-- -->![](net_plots_files/figure-gfm/network-17.png)<!-- -->![](net_plots_files/figure-gfm/network-18.png)<!-- -->![](net_plots_files/figure-gfm/network-19.png)<!-- -->![](net_plots_files/figure-gfm/network-20.png)<!-- -->![](net_plots_files/figure-gfm/network-21.png)<!-- -->![](net_plots_files/figure-gfm/network-22.png)<!-- -->![](net_plots_files/figure-gfm/network-23.png)<!-- -->![](net_plots_files/figure-gfm/network-24.png)<!-- -->![](net_plots_files/figure-gfm/network-25.png)<!-- -->![](net_plots_files/figure-gfm/network-26.png)<!-- -->![](net_plots_files/figure-gfm/network-27.png)<!-- -->![](net_plots_files/figure-gfm/network-28.png)<!-- -->![](net_plots_files/figure-gfm/network-29.png)<!-- -->![](net_plots_files/figure-gfm/network-30.png)<!-- -->![](net_plots_files/figure-gfm/network-31.png)<!-- -->![](net_plots_files/figure-gfm/network-32.png)<!-- -->![](net_plots_files/figure-gfm/network-33.png)<!-- -->![](net_plots_files/figure-gfm/network-34.png)<!-- -->![](net_plots_files/figure-gfm/network-35.png)<!-- -->![](net_plots_files/figure-gfm/network-36.png)<!-- -->![](net_plots_files/figure-gfm/network-37.png)<!-- -->![](net_plots_files/figure-gfm/network-38.png)<!-- -->![](net_plots_files/figure-gfm/network-39.png)<!-- -->![](net_plots_files/figure-gfm/network-40.png)<!-- -->![](net_plots_files/figure-gfm/network-41.png)<!-- -->![](net_plots_files/figure-gfm/network-42.png)<!-- -->![](net_plots_files/figure-gfm/network-43.png)<!-- -->![](net_plots_files/figure-gfm/network-44.png)<!-- -->![](net_plots_files/figure-gfm/network-45.png)<!-- -->![](net_plots_files/figure-gfm/network-46.png)<!-- -->![](net_plots_files/figure-gfm/network-47.png)<!-- -->![](net_plots_files/figure-gfm/network-48.png)<!-- -->![](net_plots_files/figure-gfm/network-49.png)<!-- -->![](net_plots_files/figure-gfm/network-50.png)<!-- -->![](net_plots_files/figure-gfm/network-51.png)<!-- -->![](net_plots_files/figure-gfm/network-52.png)<!-- -->![](net_plots_files/figure-gfm/network-53.png)<!-- -->![](net_plots_files/figure-gfm/network-54.png)<!-- -->![](net_plots_files/figure-gfm/network-55.png)<!-- -->![](net_plots_files/figure-gfm/network-56.png)<!-- -->![](net_plots_files/figure-gfm/network-57.png)<!-- -->![](net_plots_files/figure-gfm/network-58.png)<!-- -->![](net_plots_files/figure-gfm/network-59.png)<!-- -->![](net_plots_files/figure-gfm/network-60.png)<!-- -->![](net_plots_files/figure-gfm/network-61.png)<!-- -->![](net_plots_files/figure-gfm/network-62.png)<!-- -->![](net_plots_files/figure-gfm/network-63.png)<!-- -->![](net_plots_files/figure-gfm/network-64.png)<!-- -->![](net_plots_files/figure-gfm/network-65.png)<!-- -->![](net_plots_files/figure-gfm/network-66.png)<!-- -->![](net_plots_files/figure-gfm/network-67.png)<!-- -->![](net_plots_files/figure-gfm/network-68.png)<!-- -->![](net_plots_files/figure-gfm/network-69.png)<!-- -->![](net_plots_files/figure-gfm/network-70.png)<!-- -->![](net_plots_files/figure-gfm/network-71.png)<!-- -->![](net_plots_files/figure-gfm/network-72.png)<!-- -->![](net_plots_files/figure-gfm/network-73.png)<!-- -->![](net_plots_files/figure-gfm/network-74.png)<!-- -->![](net_plots_files/figure-gfm/network-75.png)<!-- -->![](net_plots_files/figure-gfm/network-76.png)<!-- -->![](net_plots_files/figure-gfm/network-77.png)<!-- -->![](net_plots_files/figure-gfm/network-78.png)<!-- -->![](net_plots_files/figure-gfm/network-79.png)<!-- -->![](net_plots_files/figure-gfm/network-80.png)<!-- -->![](net_plots_files/figure-gfm/network-81.png)<!-- -->![](net_plots_files/figure-gfm/network-82.png)<!-- -->![](net_plots_files/figure-gfm/network-83.png)<!-- -->![](net_plots_files/figure-gfm/network-84.png)<!-- -->![](net_plots_files/figure-gfm/network-85.png)<!-- -->![](net_plots_files/figure-gfm/network-86.png)<!-- -->![](net_plots_files/figure-gfm/network-87.png)<!-- -->![](net_plots_files/figure-gfm/network-88.png)<!-- -->![](net_plots_files/figure-gfm/network-89.png)<!-- -->![](net_plots_files/figure-gfm/network-90.png)<!-- -->![](net_plots_files/figure-gfm/network-91.png)<!-- -->![](net_plots_files/figure-gfm/network-92.png)<!-- -->![](net_plots_files/figure-gfm/network-93.png)<!-- -->![](net_plots_files/figure-gfm/network-94.png)<!-- -->![](net_plots_files/figure-gfm/network-95.png)<!-- -->![](net_plots_files/figure-gfm/network-96.png)<!-- -->![](net_plots_files/figure-gfm/network-97.png)<!-- -->![](net_plots_files/figure-gfm/network-98.png)<!-- -->![](net_plots_files/figure-gfm/network-99.png)<!-- -->![](net_plots_files/figure-gfm/network-100.png)<!-- -->![](net_plots_files/figure-gfm/network-101.png)<!-- -->![](net_plots_files/figure-gfm/network-102.png)<!-- -->![](net_plots_files/figure-gfm/network-103.png)<!-- -->![](net_plots_files/figure-gfm/network-104.png)<!-- -->![](net_plots_files/figure-gfm/network-105.png)<!-- -->
