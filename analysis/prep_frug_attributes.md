Generate frugivore attributes
================
Elena Quintero
2024-02-07

``` r
library(here)
library(tidyverse)
library(tidylog)
library(magrittr)
```

Load networks:

``` r
nets_names <- list.files(path = here("networks/standardized/"), pattern = "_int")

nets <- list()

for (i in 1:length(nets_names)){
  net_file <- paste0("networks/standardized/", nets_names[[i]])
  net <- read.csv(here(net_file))
  
  if(ncol(net)<2){
    net <- read.csv(here(net_file), sep = ";") 
  } else
  {}  

  net %<>% column_to_rownames("ind")
  nets[[i]] <- net
  names(nets)[i] <- nets_names[[i]]
}
```

Calculate aggregated interactions by bird species and define core birds:

``` r
frug.cons <- data.frame()

for (i in 1:length(nets_names)){
  tmp <- nets[[i]] %>% 
    rownames_to_column("ind") %>% 
    melt() %>%
    mutate(net_id = names(nets[i]),
           net_id = str_sub(net_id, 3, 7)) %>%
    filter(value!= 0) %>%
    group_by(net_id) %>% 
    mutate(n_ind = n_distinct(ind)) %>%
    group_by(net_id, variable) %>% 
    summarise(ints = sum(value),
              n_ind = first(n_ind),
              n_plants_int = n_distinct(ind)) %>%
    mutate(per.plants = n_plants_int/n_ind)
  
  frug.cons <- rbind(frug.cons, tmp)
}

glimpse(frug.cons)
```

    ## Rows: 579
    ## Columns: 6
    ## Groups: net_id [48]
    ## $ net_id       <chr> "01_01", "01_01", "01_01", "01_01", "01_01", "01_01", "…
    ## $ variable     <fct> Chloris_chloris, Cyanistes_caeruleus, Erithacus_rubecul…
    ## $ ints         <dbl> 2.942993e-01, 3.153531e-04, 2.317130e-01, 7.216522e-05,…
    ## $ n_ind        <int> 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,…
    ## $ n_plants_int <int> 38, 6, 40, 4, 10, 29, 25, 40, 29, 29, 9, 9, 17, 10, 15,…
    ## $ per.plants   <dbl> 0.950, 0.150, 1.000, 0.100, 0.250, 0.725, 0.625, 1.000,…

Add network information

``` r
ind.level.df <- read.csv(here("data/node.level.selection.csv"))

net_info <- ind.level.df %>% 
  distinct(net_id, net_name, plant_sp, pop, continent, country, plant_plot_rank)

head(net_info)
```

    ##   net_id            plant_sp                          pop   country
    ## 1  01_01  Pistacia lentiscus                    El Puntal     Spain
    ## 2  01_02  Pistacia lentiscus       Laguna de las Madroñas     Spain
    ## 3  02_01 Juniperus phoenicea                 Colonizacion     Spain
    ## 4  02_02 Juniperus phoenicea                       Ojillo     Spain
    ## 5  02_03 Juniperus phoenicea                   El Marqués     Spain
    ## 6  03_01 Lithraea molleoides Los Hornillos (low invasion) Argentina
    ##                                           net_name continent plant_plot_rank
    ## 1                     Pistacia lentiscus_El Puntal    Europe              27
    ## 2        Pistacia lentiscus_Laguna de las Madroñas    Europe              27
    ## 3                 Juniperus phoenicea_Colonizacion    Europe              23
    ## 4                       Juniperus phoenicea_Ojillo    Europe              23
    ## 5                   Juniperus phoenicea_El Marqués    Europe              23
    ## 6 Lithraea molleoides_Los Hornillos (low invasion)   America               1

``` r
frug.cons %<>% left_join(net_info) %>%
  rename(frug_sp = variable)

glimpse(frug.cons)
```

    ## Rows: 579
    ## Columns: 12
    ## Groups: net_id [48]
    ## $ net_id          <chr> "01_01", "01_01", "01_01", "01_01", "01_01", "01_01"…
    ## $ frug_sp         <fct> Chloris_chloris, Cyanistes_caeruleus, Erithacus_rube…
    ## $ ints            <dbl> 2.942993e-01, 3.153531e-04, 2.317130e-01, 7.216522e-…
    ## $ n_ind           <int> 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, …
    ## $ n_plants_int    <int> 38, 6, 40, 4, 10, 29, 25, 40, 29, 29, 9, 9, 17, 10, …
    ## $ per.plants      <dbl> 0.950, 0.150, 1.000, 0.100, 0.250, 0.725, 0.625, 1.0…
    ## $ plant_sp        <chr> "Pistacia lentiscus", "Pistacia lentiscus", "Pistaci…
    ## $ pop             <chr> "El Puntal", "El Puntal", "El Puntal", "El Puntal", …
    ## $ country         <chr> "Spain", "Spain", "Spain", "Spain", "Spain", "Spain"…
    ## $ net_name        <chr> "Pistacia lentiscus_El Puntal", "Pistacia lentiscus_…
    ## $ continent       <chr> "Europe", "Europe", "Europe", "Europe", "Europe", "E…
    ## $ plant_plot_rank <int> 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, …

Number of unique animal species:

``` r
length(unique(frug.cons$frug_sp)) #246 frugivore species
```

    ## [1] 246

Get frugivores body mass:

``` r
elton.bird <- readr::read_delim(here("data/animal_elton_traits/BirdFuncDat.txt"), delim = "\t") %>% 
  dplyr::select(frug_sp_elton = Scientific, family = BLFamilyLatin, bmass.elton = `BodyMass-Value`) %>%
  mutate(class = "bird")

elton.mam <- readr::read_delim(here("data/animal_elton_traits/MamFuncDat.txt"), delim = "\t") %>% 
  dplyr::select(frug_sp_elton = Scientific, family = MSWFamilyLatin, bmass.elton = `BodyMass-Value`) %>%
  mutate(class = "mammal")

elton <- rbind(elton.bird, elton.mam)

head(elton)
```

    ## # A tibble: 6 × 4
    ##   frug_sp_elton              family        bmass.elton class
    ##   <chr>                      <chr>               <dbl> <chr>
    ## 1 Struthio camelus           Struthionidae     111000  bird 
    ## 2 Rhea americana             Rheidae            23000  bird 
    ## 3 Rhea pennata               Rheidae            23900  bird 
    ## 4 Casuarius casuarius        Casuariidae        44000  bird 
    ## 5 Casuarius bennetti         Casuariidae        35000. bird 
    ## 6 Casuarius unappendiculatus Casuariidae        46074. bird

``` r
frugs <- frug.cons %>% 
  mutate(frug_sp_original = str_replace_all(frug_sp, "_", " ")) %>%
  ungroup() %>% distinct(frug_sp_original) %>%
  mutate(frug_sp_elton = frug_sp_original) %>%
  mutate(frug_sp_elton = ifelse(frug_sp_elton == "Cyanistes caeruleus", "Parus caeruleus", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Cyanopica cooki", "Cyanopica cyanus", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Saxicola rubicola", "Saxicola torquatus", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Chloris chloris", "Carduelis chloris", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Aburria jacutinga", "Pipile jacutinga", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Curruca melanocephala", "Sylvia melanocephala", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Curruca communis", "Sylvia communis", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Curruca conspicillata", "Sylvia conspicillata", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Curruca hortensis", "Sylvia hortensis", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Curruca iberiae", "Sylvia cantillans", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Curruca undata", "Sylvia undata", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Regulus ignicapillus", "Regulus ignicapilla", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Dixiphia pipra", "Pipra pipra", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Dixiphia rubrocapilla", "Pipra rubrocapilla", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Ceratopipra rubrocapilla", "Pipra rubrocapilla", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Machaeroptus regulus", "Machaeropterus regulus", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Pipraidea bonariensis", "Thraupis bonariensis", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Gennetta genetta", "Genetta genetta", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Microspingus melanoleucus", "Poospiza melanoleuca", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Thectocercus acuticaudata", "Aratinga acuticaudata", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Cyclaris gujanensis", "Cyclarhis gujanensis", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Rubigula gularis", "Pycnonotus melanicterus", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Psilopogon viridis", "Megalaima viridis", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Gracula indica", "Gracula religiosa", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Tangara brasiliensis", "Tangara mexicana", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Pteroglossus erythropygius", "Pteroglossus torquatus", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Lanio cristatus", "Tachyphonus cristatus", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Setophaga pitiayumi", "Parula pitiayumi", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Spilopelia senegalensis", "Stigmatopelia senegalensis", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Spilopelia chinensis", "Stigmatopelia chinensis", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Cinnyris asiaticus", "Nectarinia asiatica", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Leptocoma zeylonica", "Nectarinia zeylonica", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Saltatricula atricollis", "Saltator atricollis", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Semnopithecus johnii", "Trachypithecus johnii", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Turdus simillimus", "Turdus merula", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Oriolus kundoo", "Oriolus oriolus", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Geokichla citrina", "Zoothera citrina", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Hypsipetes ganeesa", "Hypsipetes leucocephalus", frug_sp_elton),
         frug_sp_elton = ifelse(frug_sp_elton == "Sturnia blythii", "Sturnus malabaricus", frug_sp_elton)) %>% 
  mutate(frug_sp_elton = ifelse(frug_sp_elton == "Rattus sp", "Rattus rattus", frug_sp_elton), #Doñana most common
         frug_sp_elton = ifelse(frug_sp_elton == "Ramphastos sp", "Ramphastos brevis", frug_sp_elton)) %>% #following Lamperty most observed species
  ungroup() %>% 
  left_join(elton, by = "frug_sp_elton") %>%
  mutate(bmass = ifelse(frug_sp_elton == "Tupinambis rufescens", 10000, bmass.elton),
         bmass = ifelse(frug_sp_elton == "Leontocebus nigrifrons", 420, bmass),
         bmass = ifelse(frug_sp_elton == "Psammodromus algirus", 11, bmass), #doi:10.1016/j.anbehav.2003.06.008
         #elton %>% filter(grepl("Amazilia", frug_sp_elton)) %>% summarise(mean(bmass.elton))
         bmass = ifelse(frug_sp_elton == "Amazilia sp", 4.66, bmass),
         #elton %>% filter(grepl("Euphonia", frug_sp_elton)) %>% summarise(mean(bmass.elton))
         bmass = ifelse(frug_sp_elton == "Euphonia sp", 13.3, bmass),
         bmass = ifelse(frug_sp_elton == "Pipreola sp", 50, bmass),
         #bmass = ifelse(frug_sp_elton == "Rodentia sp", NA, bmass),
         bmass = ifelse(frug_sp_elton == "Dacnis sp.", 13.5, bmass)
         ) %>% #mean D. cayana and D. nigripes
  mutate(class = ifelse(frug_sp_elton == "Tupinambis rufescens", "reptile", class),
         class = ifelse(frug_sp_elton == "Leontocebus nigrifrons", "mammal", class),
         class = ifelse(frug_sp_elton == "Psammodromus algirus", "reptile", class),
         class = ifelse(frug_sp_elton == "Dacnis sp.", "bird", class),
         class = ifelse(frug_sp_elton == "Amazilia sp", "bird", class),
         class = ifelse(frug_sp_elton == "Euphonia sp", "bird", class),
         class = ifelse(frug_sp_elton == "Pipreola sp", "bird", class),
         class = ifelse(frug_sp_elton == "Rodentia sp", "mammal", class),
         class = ifelse(frug_sp_elton == "Rattus sp", "mammal", class)) %>%
  arrange(family)

glimpse(frugs)
```

    ## Rows: 246
    ## Columns: 6
    ## $ frug_sp_original <chr> "Aegithalos caudatus", "Aegithina tiphia", "Galerid…
    ## $ frug_sp_elton    <chr> "Aegithalos caudatus", "Aegithina tiphia", "Galerid…
    ## $ family           <chr> "Aegithalidae", "Aegithinidae", "Alaudidae", "Bovid…
    ## $ bmass.elton      <dbl> 8.60, 12.00, 42.68, 900000.00, 292.00, 2790.79, 547…
    ## $ class            <chr> "bird", "bird", "bird", "mammal", "bird", "bird", "…
    ## $ bmass            <dbl> 8.60, 12.00, 42.68, 900000.00, 292.00, 2790.79, 547…

``` r
#frugivores with no body mass info:
frugs %>% 
  filter(is.na(bmass)) %>%
  distinct(frug_sp_original)
```

    ## # A tibble: 1 × 1
    ##   frug_sp_original
    ##   <chr>           
    ## 1 Rodentia sp

Save more complete dataset with interaction data:

``` r
frug.cons.full <- frug.cons %>%
  mutate(frug_sp_original = str_replace_all(frug_sp, "_", " ")) %>%
  left_join(frugs, by = "frug_sp_original")

glimpse(frug.cons.full)
```

    ## Rows: 579
    ## Columns: 18
    ## Groups: net_id [48]
    ## $ net_id           <chr> "01_01", "01_01", "01_01", "01_01", "01_01", "01_01…
    ## $ frug_sp          <fct> Chloris_chloris, Cyanistes_caeruleus, Erithacus_rub…
    ## $ ints             <dbl> 2.942993e-01, 3.153531e-04, 2.317130e-01, 7.216522e…
    ## $ n_ind            <int> 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,…
    ## $ n_plants_int     <int> 38, 6, 40, 4, 10, 29, 25, 40, 29, 29, 9, 9, 17, 10,…
    ## $ per.plants       <dbl> 0.950, 0.150, 1.000, 0.100, 0.250, 0.725, 0.625, 1.…
    ## $ plant_sp         <chr> "Pistacia lentiscus", "Pistacia lentiscus", "Pistac…
    ## $ pop              <chr> "El Puntal", "El Puntal", "El Puntal", "El Puntal",…
    ## $ country          <chr> "Spain", "Spain", "Spain", "Spain", "Spain", "Spain…
    ## $ net_name         <chr> "Pistacia lentiscus_El Puntal", "Pistacia lentiscus…
    ## $ continent        <chr> "Europe", "Europe", "Europe", "Europe", "Europe", "…
    ## $ plant_plot_rank  <int> 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27,…
    ## $ frug_sp_original <chr> "Chloris chloris", "Cyanistes caeruleus", "Erithacu…
    ## $ frug_sp_elton    <chr> "Carduelis chloris", "Parus caeruleus", "Erithacus …
    ## $ family           <chr> "Fringillidae", "Paridae", "Muscicapidae", "Sylviid…
    ## $ bmass.elton      <dbl> 26.00, 13.30, 17.70, 11.00, 19.60, 14.59, 15.10, 11…
    ## $ class            <chr> "bird", "bird", "bird", "bird", "bird", "bird", "bi…
    ## $ bmass            <dbl> 26.00, 13.30, 17.70, 11.00, 19.60, 14.59, 15.10, 11…

``` r
write_csv(frug.cons.full, here("data/frugivores_interactions.csv"))
```
