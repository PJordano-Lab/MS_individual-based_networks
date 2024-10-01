Null model network comparisons
================
Elena Quintero
2024-09-30

``` r
library(here)
library(tidyverse)
library(magrittr)
library(bipartite)
library(igraph)
library(tnet)
```

## NULL MODELS

Read matrices:

First I create a function to calculate for each matrix, the scale factor
necessary so that every non-zero value has at least number different
that zero before the decimals. Then I do a ceiling (function to round
interaction to integers). By using a scale factor decimals are not lost.

``` r
#Read list of names

net_names <- list.files(path = here("networks/nets_std"), pattern = "_int")

#Create empty list

nets <- list()

#Add all read files to the list

for (i in 1:length(net_names)){
  
  net_file <- paste0("networks/nets_std/", net_names[[i]])
  net <- read.csv(here(net_file))
  
  if (str_detect(net_names[i], "sp_")) {
    
    net %<>% column_to_rownames("sp")
    
  }
  
  else {
    
    net %<>% column_to_rownames("ind")
    
  }
  
  nets[[i]] <- net #add to list
  names(nets)[i] <- net_names[[i]] #name the list
  
}
```

First I create a function to calculate for each matrix, the scale factor
necessary so that every non-zero value has at least number different
that zero before the decimals. Then I do a ceiling (function to round
interaction to integers). By using a scale factor small decimals are not
lost.

``` r
calc_order_magnitude <- function(matrix) {
  # Excluir ceros para evitar división por cero
  non_zero_values <- matrix[matrix != 0]
  
  # Obtener el valor mínimo absoluto de los valores no cero
  min_value <- min(abs(non_zero_values))
  
  # Calcular el orden de magnitud necesario para convertirlo en un número mayor o igual a 1
  order_magnitude <- 10^ceiling(-log10(min_value))
  
  return(order_magnitude)
}

# Calcular el factor de escala para todas las redes
scale_factors <- lapply(nets, calc_order_magnitude)
```

``` r
nets_multiplied <- list()

for (i in 1:length(nets)){
  
  nets_multiplied[[i]] <- ceiling(nets[[i]] * scale_factors[[i]])
  
  names(nets_multiplied)[i] <- net_names[[i]] #name the list

}
```

Generate 1000 random nets following *Patefield*’s algorithm (keeps
constant marginal rows and columns)

``` r
null.list <- list()
#null.list <- vector("list", 5)

for (i in 1:length(nets_multiplied)) {
  
  print(i)
  
  nulls <- nullmodel(nets_multiplied[[i]], N = 1000, method = "r2dtable") 
  
  scale_back <- function(list){
    list / scale_factors[[i]]
    }
  
  nulls.fix <- lapply(nulls, scale_back)
  
  null.list <- append(null.list, list(nulls.fix))
  #null.list[[i]] <- nulls.fix
  
  names(null.list)[[i]] <- names(nets_multiplied)[i]
  
}
```

    ## [1] 1
    ## [1] 2
    ## [1] 3
    ## [1] 4
    ## [1] 5
    ## [1] 6
    ## [1] 7
    ## [1] 8
    ## [1] 9
    ## [1] 10
    ## [1] 11
    ## [1] 12
    ## [1] 13
    ## [1] 14
    ## [1] 15
    ## [1] 16
    ## [1] 17
    ## [1] 18
    ## [1] 19
    ## [1] 20
    ## [1] 21
    ## [1] 22
    ## [1] 23
    ## [1] 24
    ## [1] 25
    ## [1] 26
    ## [1] 27
    ## [1] 28
    ## [1] 29
    ## [1] 30
    ## [1] 31
    ## [1] 32
    ## [1] 33
    ## [1] 34
    ## [1] 35
    ## [1] 36
    ## [1] 37
    ## [1] 38
    ## [1] 39
    ## [1] 40
    ## [1] 41
    ## [1] 42
    ## [1] 43
    ## [1] 44
    ## [1] 45
    ## [1] 46
    ## [1] 47
    ## [1] 48
    ## [1] 49
    ## [1] 50
    ## [1] 51
    ## [1] 52
    ## [1] 53
    ## [1] 54
    ## [1] 55
    ## [1] 56
    ## [1] 57
    ## [1] 58
    ## [1] 59
    ## [1] 60
    ## [1] 61
    ## [1] 62
    ## [1] 63
    ## [1] 64
    ## [1] 65
    ## [1] 66
    ## [1] 67
    ## [1] 68
    ## [1] 69
    ## [1] 70
    ## [1] 71
    ## [1] 72
    ## [1] 73
    ## [1] 74
    ## [1] 75
    ## [1] 76
    ## [1] 77
    ## [1] 78
    ## [1] 79
    ## [1] 80
    ## [1] 81
    ## [1] 82
    ## [1] 83
    ## [1] 84
    ## [1] 85
    ## [1] 86
    ## [1] 87
    ## [1] 88
    ## [1] 89
    ## [1] 90
    ## [1] 91
    ## [1] 92
    ## [1] 93
    ## [1] 94
    ## [1] 95
    ## [1] 96
    ## [1] 97
    ## [1] 98
    ## [1] 99
    ## [1] 100
    ## [1] 101
    ## [1] 102
    ## [1] 103
    ## [1] 104
    ## [1] 105

Generate 1000 random nets following *Vazquez* algorithm (keeps constant
marginal rows and columns)

``` r
#vaznull.list <- list()
vaznull.list <- vector("list", 5)

for (i in 1:length(nets_multiplied)) {
  
  print(i)
  
  nulls <- nullmodel(nets_multiplied[[i]], N = 1000, method = "vaznull")
  
  scale_back <- function(list){
    list / scale_factors[[i]]
    }
  
  nulls.fix <- lapply(nulls, scale_back)
  
  #vaznull.list <- append(vaznull.list, list(nulls.fix))
  vaznull.list[[i]] <- nulls.fix
  
  names(vaznull.list)[[i]] <- names(nets_multiplied)[i]
  
}
```

    ## [1] 1
    ## [1] 2
    ## [1] 3
    ## [1] 4
    ## [1] 5
    ## [1] 6
    ## [1] 7
    ## [1] 8
    ## [1] 9
    ## [1] 10
    ## [1] 11
    ## [1] 12
    ## [1] 13
    ## [1] 14
    ## [1] 15
    ## [1] 16
    ## [1] 17
    ## [1] 18
    ## [1] 19
    ## [1] 20
    ## [1] 21
    ## [1] 22
    ## [1] 23
    ## [1] 24
    ## [1] 25
    ## [1] 26
    ## [1] 27
    ## [1] 28
    ## [1] 29
    ## [1] 30
    ## [1] 31
    ## [1] 32
    ## [1] 33
    ## [1] 34
    ## [1] 35
    ## [1] 36
    ## [1] 37
    ## [1] 38
    ## [1] 39
    ## [1] 40
    ## [1] 41
    ## [1] 42
    ## [1] 43
    ## [1] 44
    ## [1] 45
    ## [1] 46
    ## [1] 47
    ## [1] 48
    ## [1] 49
    ## [1] 50
    ## [1] 51
    ## [1] 52
    ## [1] 53
    ## [1] 54
    ## [1] 55
    ## [1] 56
    ## [1] 57
    ## [1] 58
    ## [1] 59
    ## [1] 60
    ## [1] 61
    ## [1] 62
    ## [1] 63
    ## [1] 64
    ## [1] 65
    ## [1] 66
    ## [1] 67
    ## [1] 68
    ## [1] 69
    ## [1] 70
    ## [1] 71
    ## [1] 72
    ## [1] 73
    ## [1] 74
    ## [1] 75
    ## [1] 76
    ## [1] 77
    ## [1] 78
    ## [1] 79
    ## [1] 80
    ## [1] 81
    ## [1] 82
    ## [1] 83
    ## [1] 84
    ## [1] 85
    ## [1] 86
    ## [1] 87
    ## [1] 88
    ## [1] 89
    ## [1] 90
    ## [1] 91
    ## [1] 92
    ## [1] 93
    ## [1] 94
    ## [1] 95
    ## [1] 96
    ## [1] 97
    ## [1] 98
    ## [1] 99
    ## [1] 100
    ## [1] 101
    ## [1] 102
    ## [1] 103
    ## [1] 104
    ## [1] 105

Check Vazquez null models sum up to 1:

``` r
sum(vaznull.list[[1]][[1]])
```

    ## [1] 1.000199

``` r
sum(vaznull.list[[81]][[981]])
```

    ## [1] 1.019

## Calculate network level metrics:

``` r
indices <- as.list(c("number of species",
                     "connectance", 
                     "weighted NODF",
                     "interaction evenness", 
                     "Alatalo interaction evenness"))
```

``` r
# Progress bar initialization
pb <- txtProgressBar(min = 0, max = length(nets), style = 3)
```

    ##   |                                                                              |                                                                      |   0%

``` r
# Initialize web.metrics list
web.metrics <- lapply(1:length(vaznull.list), function(x) vector("list", 1000))

# Function to compute bipartite metrics
compute_bipartite_metrics <- function(net) {
  m.bipart <- networklevel(net, indices)
  return(as.data.frame(t(m.bipart)))
}

# Function to compute modularity
compute_modularity <- function(net) {
  mod <- computeModules(net)
  return(mod@likelihood)
}

# Function to compute igraph metrics
compute_igraph_metrics <- function(net) {
  inet <- graph_from_biadjacency_matrix(net, weighted = TRUE, add.names = NULL)
  assortativity <- assortativity_degree(inet)
  eigen.centrality <- eigen_centrality(inet)$value
  centr_binary <- centr_eigen(inet, scale = TRUE, normalized = TRUE)$centralization
  centr_weighted_obs <- sum(max(eigen_centrality(inet)$vector) - eigen_centrality(inet)$vector)
  centralization_w <- centr_weighted_obs / centr_eigen(inet)$theoretical_max
  return(cbind(assortativity, eigen.centrality, centr_binary, centralization_w))
}

# Umbrella function to process each network and iteration
process_network <- function(i, null_list) {
  
  net_nulls <- null_list[[i]]
  
  metrics_list <- lapply(1:1000, function(j) {
    
    m.bipart <- compute_bipartite_metrics(net_nulls[[j]])
    
    M <- compute_modularity(net_nulls[[j]])
    
    igraph.metrics <- compute_igraph_metrics(net_nulls[[j]])
    
    # Merge all metrics
    metrics.nulls <- cbind(m.bipart, M, igraph.metrics) %>%
      select_all(~gsub("\\s+", ".", .)) %>%
      mutate(net_size = number.of.species.HL * number.of.species.LL,
             iter = j)
    
    # Add type and net_id
    if (str_detect(names(null_list)[i], "sp")) {
      metrics.nulls <- metrics.nulls %>%
        mutate(type = "sp", net_id = str_sub(names(null_list)[i], 4, 8))
    } else {
      metrics.nulls <- metrics.nulls %>%
        mutate(type = "ind", net_id = str_sub(names(null_list)[i], 1, 5))
    }
    
    return(metrics.nulls)
  })
  
  return(metrics_list)
}
```

Run umbrella function for each *Vazquez* null model network:

``` r
# Main lapply for each network
web.metrics <- lapply(1:length(vaznull.list), function(i) {
  metrics <- process_network(i, vaznull.list)
  
  # Update progress bar
  setTxtProgressBar(pb, i)
  
  return(metrics)
})
```

    ##   |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%

Makes one data.table from a list of many

``` r
web.metrics.df <- data.table::rbindlist(
  lapply(1:length(web.metrics), function(i) {
  data.table::rbindlist(web.metrics[[i]])
  })
  )

web.metrics.df <- web.metrics.df |> 
  relocate(type, .before = everything()) |> 
  relocate(net_id, .after = type) |> 
  relocate(iter, .after = net_id) |> 
  mutate(net_code = paste0(type, "_", net_id)) |> 
  group_by(iter) |> 
  mutate(net_n = order(net_code))

glimpse(web.metrics.df)
```

    ## Rows: 105,000
    ## Columns: 17
    ## Groups: iter [1,000]
    ## $ type                         <chr> "ind", "ind", "ind", "ind", "ind", "ind",…
    ## $ net_id                       <chr> "01_01", "01_01", "01_01", "01_01", "01_0…
    ## $ iter                         <int> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13…
    ## $ connectance                  <dbl> 0.362963, 0.362963, 0.362963, 0.362963, 0…
    ## $ weighted.NODF                <dbl> 65.84884, 62.68206, 66.91053, 64.04740, 6…
    ## $ interaction.evenness         <dbl> 0.6949892, 0.6962069, 0.6971196, 0.696918…
    ## $ Alatalo.interaction.evenness <dbl> 0.6206832, 0.6188404, 0.6158226, 0.617120…
    ## $ number.of.species.HL         <dbl> 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 2…
    ## $ number.of.species.LL         <dbl> 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 4…
    ## $ M                            <dbl> 0.013992342, 0.013195380, 0.012084882, 0.…
    ## $ assortativity                <dbl> -0.5815894, -0.5672870, -0.5823415, -0.56…
    ## $ eigen.centrality             <dbl> 0.1117928, 0.1114919, 0.1114095, 0.111374…
    ## $ centr_binary                 <dbl> 0.6122276, 0.6112887, 0.6092668, 0.611419…
    ## $ centralization_w             <dbl> 0.8795333, 0.8790228, 0.8790518, 0.879315…
    ## $ net_size                     <dbl> 1080, 1080, 1080, 1080, 1080, 1080, 1080,…
    ## $ net_code                     <chr> "ind_01_01", "ind_01_01", "ind_01_01", "i…
    ## $ net_n                        <int> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…

``` r
saveRDS(web.metrics.df, here("data/net_level_nulls_vazquez.rds"))
```

Run umbrella function for each *Patefield* null model network:

``` r
web.metrics.pat <- lapply(1:length(null.list), function(i) {
  metrics <- process_network(i, null.list)

  # Update progress bar
  setTxtProgressBar(pb, i)

  return(metrics)
})
```

    ##   |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%

``` r
web.metrics.pat.df <- data.table::rbindlist(
  lapply(1:length(web.metrics.pat), function(i) {
  data.table::rbindlist(web.metrics.pat[[i]])
  })
  )

web.metrics.pat.df <- web.metrics.pat.df |>
  relocate(type, .before = everything()) |>
  relocate(net_id, .after = type) |>
  relocate(iter, .after = net_id) |>
  mutate(net_code = paste0(type, "_", net_id)) |>
  group_by(iter) |>
  mutate(net_n = order(net_code))

glimpse(web.metrics.pat.df)
```

    ## Rows: 105,000
    ## Columns: 17
    ## Groups: iter [1,000]
    ## $ type                         <chr> "ind", "ind", "ind", "ind", "ind", "ind",…
    ## $ net_id                       <chr> "01_01", "01_01", "01_01", "01_01", "01_0…
    ## $ iter                         <int> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13…
    ## $ connectance                  <dbl> 0.9101852, 0.9101852, 0.9157407, 0.907407…
    ## $ weighted.NODF                <dbl> 61.22583, 59.09732, 56.77112, 56.01830, 6…
    ## $ interaction.evenness         <dbl> 0.7145744, 0.7145695, 0.7145706, 0.714569…
    ## $ Alatalo.interaction.evenness <dbl> 0.5703252, 0.5708836, 0.5712427, 0.570515…
    ## $ number.of.species.HL         <dbl> 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 2…
    ## $ number.of.species.LL         <dbl> 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 4…
    ## $ M                            <dbl> 0.003141923, 0.003141632, 0.003284954, 0.…
    ## $ assortativity                <dbl> -0.7973251, -0.7975083, -0.8114309, -0.78…
    ## $ eigen.centrality             <dbl> 0.1089190, 0.1088673, 0.1088336, 0.108902…
    ## $ centr_binary                 <dbl> 0.1794218, 0.1797270, 0.1751363, 0.181908…
    ## $ centralization_w             <dbl> 0.8786098, 0.8784934, 0.8782771, 0.878311…
    ## $ net_size                     <dbl> 1080, 1080, 1080, 1080, 1080, 1080, 1080,…
    ## $ net_code                     <chr> "ind_01_01", "ind_01_01", "ind_01_01", "i…
    ## $ net_n                        <int> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…

``` r
saveRDS(web.metrics.pat.df, here("data/net_level_nulls.rds"))
```

Compare null models for same net using different scale factors
(e.g. 1000, 10000 and 100000):
