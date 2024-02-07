Webs standarisation
================
Elena Quintero
2024-02-07

``` r
library(here)
library(dplyr)
library(tidyverse)
```

Get list of net names in web folder:

``` r
nets_names_ind <- list.files(path = here("networks/"), pattern = "_int") # For individual-based networks
nets_names_sp <- list.files(path = here("networks/"), pattern = "sp_")   # For species-based networks
```

Read webs and apply standardization:

``` r
# Individual-based networks standardization
for (i in 1:length(nets_names_ind)){
  net_file <- paste0("networks/", nets_names_ind[[i]])
  net <- read.csv(here(net_file))
  if(ncol(net)<2){
    net <- read.csv(here(net_file), sep = ";")
  } else
  {}  
  
  net_mat <- net %>% column_to_rownames("ind") %>% as.matrix()
  
  mat_gts <- net_mat/sum(net_mat)
  
  net_gts <- as.data.frame(mat_gts) %>% rownames_to_column("ind")
  
  assign(paste0("st", nets_names_ind[[i]]), net_gts)
  
  write_csv(net_gts, here(paste0("networks/standardized/", "st", nets_names_ind[[i]])))
}

# Species-based networks standardization
for (i in 1:length(nets_names_sp)){
  net_file <- paste0("networks/", nets_names_sp[[i]])
  net <- read.csv(here(net_file))
  if(ncol(net)<2){
    net <- read.csv(here(net_file), sep = ";")
  } else
  {}  
  
  net_mat <- net %>% column_to_rownames("sp") %>% as.matrix()
  
  mat_gts <- net_mat/sum(net_mat)
  
  net_gts <- as.data.frame(mat_gts) %>% rownames_to_column("sp")
  
  assign(paste0("st", nets_names_sp[[i]]), net_gts)
  
  write_csv(net_gts, here(paste0("networks/standardized/", "st", nets_names_sp[[i]])))
}
```

    ## Warning in read.table(file = file, header = header, sep = sep, quote = quote, : incomplete final line found by
    ## readTableHeader on '/Users/Elena/Documents/PROJECTS/MS_individual-based_networks/networks/sp_03_01.csv'

    ## Warning in read.table(file = file, header = header, sep = sep, quote = quote, : incomplete final line found by
    ## readTableHeader on '/Users/Elena/Documents/PROJECTS/MS_individual-based_networks/networks/sp_03_01.csv'

Check all new nets sum up to 1:

``` r
nets_names_st_ind <- list.files(path = here("networks/standardized/"), pattern = "_int") # For individual-based networks
nets_names_st_sp <- list.files(path = here("networks/standardized/"), pattern = "stsp_") # For species-based networks

# For individual-based networks
for (i in 1:length(nets_names_st_ind)){
  net_file <- paste0("networks/standardized/", nets_names_st_ind[[i]])
  net <- read.csv(here(net_file))
  
  print(nets_names_st_ind[i])
  print(sum(net[,-1]))
}
```

    ## [1] "st01_01_int.csv"
    ## [1] 1
    ## [1] "st01_02_int.csv"
    ## [1] 1
    ## [1] "st02_01_int.csv"
    ## [1] 1
    ## [1] "st02_02_int.csv"
    ## [1] 1
    ## [1] "st02_03_int.csv"
    ## [1] 1
    ## [1] "st03_01_int.csv"
    ## [1] 1
    ## [1] "st03_02_int.csv"
    ## [1] 1
    ## [1] "st03_03_int.csv"
    ## [1] 1
    ## [1] "st03_04_int.csv"
    ## [1] 1
    ## [1] "st03_05_int.csv"
    ## [1] 1
    ## [1] "st03_06_int.csv"
    ## [1] 1
    ## [1] "st04_01_int.csv"
    ## [1] 1
    ## [1] "st05_01_int.csv"
    ## [1] 1
    ## [1] "st06_01_int.csv"
    ## [1] 1
    ## [1] "st06_02_int.csv"
    ## [1] 1
    ## [1] "st06_03_int.csv"
    ## [1] 1
    ## [1] "st07_01_int.csv"
    ## [1] 1
    ## [1] "st08_01_int.csv"
    ## [1] 1
    ## [1] "st08_02_int.csv"
    ## [1] 1
    ## [1] "st08_03_int.csv"
    ## [1] 1
    ## [1] "st08_04_int.csv"
    ## [1] 1
    ## [1] "st09_01_int.csv"
    ## [1] 1
    ## [1] "st10_01_int.csv"
    ## [1] 1
    ## [1] "st11_01_int.csv"
    ## [1] 1
    ## [1] "st12_01_int.csv"
    ## [1] 1
    ## [1] "st12_02_int.csv"
    ## [1] 1
    ## [1] "st12_03_int.csv"
    ## [1] 1
    ## [1] "st12_04_int.csv"
    ## [1] 1
    ## [1] "st12_05_int.csv"
    ## [1] 1
    ## [1] "st12_06_int.csv"
    ## [1] 1
    ## [1] "st12_07_int.csv"
    ## [1] 1
    ## [1] "st12_08_int.csv"
    ## [1] 1
    ## [1] "st12_09_int.csv"
    ## [1] 1
    ## [1] "st12_10_int.csv"
    ## [1] 1
    ## [1] "st13_01_int.csv"
    ## [1] 1
    ## [1] "st13_02_int.csv"
    ## [1] 1
    ## [1] "st14_01_int.csv"
    ## [1] 1
    ## [1] "st15_01_int.csv"
    ## [1] 1
    ## [1] "st16_01_int.csv"
    ## [1] 1
    ## [1] "st16_02_int.csv"
    ## [1] 1
    ## [1] "st16_03_int.csv"
    ## [1] 1
    ## [1] "st16_04_int.csv"
    ## [1] 1
    ## [1] "st16_05_int.csv"
    ## [1] 1
    ## [1] "st17_01_int.csv"
    ## [1] 1
    ## [1] "st18_01_int.csv"
    ## [1] 1
    ## [1] "st18_02_int.csv"
    ## [1] 1
    ## [1] "st19_01_int.csv"
    ## [1] 1
    ## [1] "st20_01_int.csv"
    ## [1] 1

``` r
# For species-based networks
for (i in 1:length(nets_names_st_sp)){
  net_file <- paste0("networks/standardized/", nets_names_st_sp[[i]])
  net <- read.csv(here(net_file))
  
  print(nets_names_st_sp[i])
  print(sum(net[,-1]))
}
```

    ## [1] "stsp_01_01.csv"
    ## [1] 1
    ## [1] "stsp_01_02.csv"
    ## [1] 1
    ## [1] "stsp_02_01.csv"
    ## [1] 1
    ## [1] "stsp_03_01.csv"
    ## [1] 1
    ## [1] "stsp_03_02.csv"
    ## [1] 1
    ## [1] "stsp_03_03.csv"
    ## [1] 1
    ## [1] "stsp_03_04.csv"
    ## [1] 1
    ## [1] "stsp_04_01.csv"
    ## [1] 1
    ## [1] "stsp_05_01.csv"
    ## [1] 1
    ## [1] "stsp_06_01.csv"
    ## [1] 1
    ## [1] "stsp_07_01.csv"
    ## [1] 1
    ## [1] "stsp_08_01.csv"
    ## [1] 1
    ## [1] "stsp_09_01.csv"
    ## [1] 1
    ## [1] "stsp_10_01.csv"
    ## [1] 1
    ## [1] "stsp_11_01.csv"
    ## [1] 1
    ## [1] "stsp_12_01.csv"
    ## [1] 1
    ## [1] "stsp_13_01.csv"
    ## [1] 1
    ## [1] "stsp_14_01.csv"
    ## [1] 1
    ## [1] "stsp_15_01.csv"
    ## [1] 1
    ## [1] "stsp_15_02.csv"
    ## [1] 1
    ## [1] "stsp_15_03.csv"
    ## [1] 1
    ## [1] "stsp_15_04.csv"
    ## [1] 1
    ## [1] "stsp_16_01.csv"
    ## [1] 1
    ## [1] "stsp_17_01.csv"
    ## [1] 1
    ## [1] "stsp_18_01.csv"
    ## [1] 1
    ## [1] "stsp_19_01.csv"
    ## [1] 1
    ## [1] "stsp_20_01.csv"
    ## [1] 1
    ## [1] "stsp_21_01.csv"
    ## [1] 1
    ## [1] "stsp_22_01.csv"
    ## [1] 1
    ## [1] "stsp_23_01.csv"
    ## [1] 1
    ## [1] "stsp_24_01.csv"
    ## [1] 1
    ## [1] "stsp_25_01.csv"
    ## [1] 1
    ## [1] "stsp_25_02.csv"
    ## [1] 1
    ## [1] "stsp_25_03.csv"
    ## [1] 1
    ## [1] "stsp_26_01.csv"
    ## [1] 1
    ## [1] "stsp_27_01.csv"
    ## [1] 1
    ## [1] "stsp_28_01.csv"
    ## [1] 1
    ## [1] "stsp_29_01.csv"
    ## [1] 1
    ## [1] "stsp_29_02.csv"
    ## [1] 1
    ## [1] "stsp_29_03.csv"
    ## [1] 1
    ## [1] "stsp_29_04.csv"
    ## [1] 1
    ## [1] "stsp_29_05.csv"
    ## [1] 1
    ## [1] "stsp_29_06.csv"
    ## [1] 1
    ## [1] "stsp_29_07.csv"
    ## [1] 1
    ## [1] "stsp_29_08.csv"
    ## [1] 1
    ## [1] "stsp_29_09.csv"
    ## [1] 1
    ## [1] "stsp_29_10.csv"
    ## [1] 1
    ## [1] "stsp_29_11.csv"
    ## [1] 1
    ## [1] "stsp_29_12.csv"
    ## [1] 1
    ## [1] "stsp_29_13.csv"
    ## [1] 1
    ## [1] "stsp_30_01.csv"
    ## [1] 1
    ## [1] "stsp_31_01.csv"
    ## [1] 1
    ## [1] "stsp_32_01.csv"
    ## [1] 1
    ## [1] "stsp_33_01.csv"
    ## [1] 1
    ## [1] "stsp_34_01.csv"
    ## [1] 1
    ## [1] "stsp_35_01.csv"
    ## [1] 1
    ## [1] "stsp_36_01.csv"
    ## [1] 1
    ## [1] "stsp_37_01.csv"
    ## [1] 1
    ## [1] "stsp_38_01.csv"
    ## [1] 1
    ## [1] "stsp_38_02.csv"
    ## [1] 1
    ## [1] "stsp_39_01.csv"
    ## [1] 1
    ## [1] "stsp_40_01.csv"
    ## [1] 1
    ## [1] "stsp_41_01.csv"
    ## [1] 1
    ## [1] "stsp_42_01.csv"
    ## [1] 1
    ## [1] "stsp_43_01.csv"
    ## [1] 1
    ## [1] "stsp_44_01.csv"
    ## [1] 1
    ## [1] "stsp_45_01.csv"
    ## [1] 1
    ## [1] "stsp_45_02.csv"
    ## [1] 1

# 
