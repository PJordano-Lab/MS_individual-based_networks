Network selection and dataset preparation
================
Elena Quintero
2024-02-07

``` r
library(here)
library(tidyverse)
library(magrittr)
library(tidylog)

theme_set(theme_minimal())
```

# 1. Network level

Load ind-sp data:

``` r
net.level.ind <- read.csv(here("data/net.level.df.csv")) %>%
  rename(niche.overlap.horn.HL = niche.overlap.HL) %>%
  rename(niche.overlap.horn.LL = niche.overlap.LL)
```

    ## rename: renamed one variable (niche.overlap.horn.HL)

    ## rename: renamed one variable (niche.overlap.horn.LL)

``` r
mod.ind <- read.csv(here("data/modularity_estimation_gts.csv")) %>%
  rename(net_id = network) %>%
  mutate(net_id = str_sub(net_id, 3, 7)) %>% 
  dplyr::select(-X)
```

    ## rename: renamed one variable (net_id)

    ## mutate: changed 48 values (100%) of 'net_id' (0 new NA)

``` r
net.level.ind %<>% left_join(mod.ind, by = "net_id") %>%
  relocate(M, .after = weighted.NODF) %>%
  mutate(type = "ind")
```

    ## left_join: added one column (M)

    ##            > rows only in x    0

    ##            > rows only in y  ( 0)

    ##            > matched rows     48

    ##            >                 ====

    ##            > rows total       48

    ## relocate: columns reordered (connectance, web.asymmetry, links.per.species, cluster.coefficient, NODF, …)

    ## mutate: new variable 'type' (character) with one unique value and 0% NA

Load sp-sp data:

``` r
net.level.sp <- read.csv(here("data/net.level.df.sp.csv")) %>%
  mutate(net_id = str_sub(net_id, 6, 10)) %>%
  rename(niche.overlap.horn.HL = niche.overlap.HL) %>%
  rename(niche.overlap.horn.LL = niche.overlap.LL)
```

    ## mutate: changed 68 values (100%) of 'net_id' (0 new NA)

    ## rename: renamed one variable (niche.overlap.horn.HL)

    ## rename: renamed one variable (niche.overlap.horn.LL)

``` r
mod.sp <- read.csv(here("data/modularity_estimation_gts_sp.csv")) %>%
  rename(net_id = network) %>%
  mutate(net_id = str_sub(net_id, 6, 10)) %>% 
  dplyr::select(-X)
```

    ## rename: renamed one variable (net_id)

    ## mutate: changed 68 values (100%) of 'net_id' (0 new NA)

``` r
net.level.sp %<>% left_join(mod.sp, by = "net_id") %>%
  relocate(M, .after = weighted.NODF) %>%
  mutate(type = "sp")
```

    ## left_join: added one column (M)

    ##            > rows only in x    0

    ##            > rows only in y  ( 0)

    ##            > matched rows     68

    ##            >                 ====

    ##            > rows total       68

    ## relocate: columns reordered (connectance, web.asymmetry, links.per.species, cluster.coefficient, NODF, …)

    ## mutate: new variable 'type' (character) with one unique value and 0% NA

Load dissasortativity and centrality metrics:

``` r
igraph.metrics <- read.csv(here("data/igraph_metrics.csv")) %>%
  mutate(net_code = paste0(type, "_", net_id)) %>%
  dplyr::select(-c(X, net_id))
```

    ## mutate: new variable 'net_code' (character) with 116 unique values and 0% NA

Load net characteristics:

``` r
all_nets_char <- readxl::read_xlsx(here("data/selected_networks.xlsx")) %>%
  mutate(net_code = paste0(type, "_", code_ID))
```

    ## mutate: new variable 'net_code' (character) with 116 unique values and 0% NA

``` r
nets_char <- all_nets_char %>%
  dplyr::select(net_code, type, code_ID, ref, bioregion, country, plant_sp)

#glimpse(nets_char)
```

Merge all datasets:

``` r
net.level <- full_join(net.level.ind, net.level.sp) %>%
  mutate(net_code = paste0(type, "_", net_id)) %>%
  left_join(igraph.metrics) %>%
  mutate(unique.ints = connectance * net_size) %>% 
  left_join(nets_char)
```

    ## Joining, by = c("connectance", "web.asymmetry", "links.per.species",
    ## "cluster.coefficient", "NODF", "weighted.NODF", "M",
    ## "interaction.strength.asymmetry", "specialisation.asymmetry",
    ## "linkage.density", "weighted.connectance", "Shannon.diversity",
    ## "interaction.evenness", "Alatalo.interaction.evenness", "H2",
    ## "number.of.species.HL", "number.of.species.LL", "cluster.coefficient.HL",
    ## "cluster.coefficient.LL", "weighted.cluster.coefficient.HL",
    ## "weighted.cluster.coefficient.LL", "niche.overlap.horn.HL",
    ## "niche.overlap.horn.LL", "extinction.slope.HL", "extinction.slope.LL",
    ## "robustness.HL", "robustness.LL", "generality.HL", "vulnerability.LL",
    ## "niche.overlap.bray.HL", "niche.overlap.bray.LL", "niche.overlap.jaccard.HL",
    ## "niche.overlap.jaccard.LL", "net_size", "net_id", "links", "net_n", "type")
    ## full_join: added no columns
    ## > rows only in x 48
    ## > rows only in y 68
    ## > matched rows 0
    ## > =====
    ## > rows total 116
    ## mutate: new variable 'net_code' (character) with 116 unique values and 0% NA
    ## Joining, by = c("type", "net_code")
    ## left_join: added 3 columns (assortativity, centr_binary, centralization.w)
    ## > rows only in x 0
    ## > rows only in y ( 0)
    ## > matched rows 116
    ## > =====
    ## > rows total 116
    ## mutate: new variable 'unique.ints' (double) with 85 unique values and 0% NA
    ## Joining, by = c("type", "net_code")
    ## left_join: added 5 columns (code_ID, ref, bioregion, country, plant_sp)
    ## > rows only in x 0
    ## > rows only in y ( 0)
    ## > matched rows 116
    ## > =====
    ## > rows total 116

``` r
#glimpse(net.level)
```

``` r
net.level %>% 
  #filter(type == "sp") %>%
  select(type, plant_sp, net_id, plant_sp, 
         number.of.species.HL, number.of.species.LL,
         net_size, unique.ints) %>%
  mutate(n_nodes = number.of.species.HL + number.of.species.LL) %>%
  arrange(n_nodes) %>%
  filter(n_nodes > 15) %>% 
  #count(as.factor(type)) %>%
  summary()
```

    ## select: dropped 41 variables (connectance, web.asymmetry, links.per.species, cluster.coefficient, NODF, …)

    ## mutate: new variable 'n_nodes' (integer) with 51 unique values and 0% NA

    ## filter: removed 11 rows (9%), 105 rows remaining

    ##      type             plant_sp            net_id         
    ##  Length:105         Length:105         Length:105        
    ##  Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character  
    ##                                                          
    ##                                                          
    ##                                                          
    ##  number.of.species.HL number.of.species.LL    net_size       unique.ints   
    ##  Min.   : 5.0         Min.   : 5.00        Min.   :  55.0   Min.   : 21.0  
    ##  1st Qu.:10.0         1st Qu.:11.00        1st Qu.: 135.0   1st Qu.: 41.0  
    ##  Median :14.0         Median :17.00        Median : 261.0   Median : 72.0  
    ##  Mean   :17.6         Mean   :21.37        Mean   : 380.9   Mean   : 89.5  
    ##  3rd Qu.:21.0         3rd Qu.:29.00        3rd Qu.: 476.0   3rd Qu.:119.0  
    ##  Max.   :88.0         Max.   :77.00        Max.   :2904.0   Max.   :419.0  
    ##     n_nodes      
    ##  Min.   : 16.00  
    ##  1st Qu.: 25.00  
    ##  Median : 35.00  
    ##  Mean   : 38.97  
    ##  3rd Qu.: 47.00  
    ##  Max.   :121.00

Remove small nets:

``` r
net.level.selection <- net.level %>%
  mutate(n_nodes = number.of.species.HL + number.of.species.LL) %>%
  filter(n_nodes > 15)
```

    ## mutate: new variable 'n_nodes' (integer) with 51 unique values and 0% NA

    ## filter: removed 11 rows (9%), 105 rows remaining

``` r
removed_nets <- setdiff(net.level$net_code, net.level.selection$net_code)
removed_nets
```

    ##  [1] "ind_08_01" "ind_12_01" "ind_12_06" "ind_12_07" "sp_03_01"  "sp_03_02" 
    ##  [7] "sp_03_03"  "sp_03_04"  "sp_05_01"  "sp_31_01"  "sp_44_01"

11 nets removed (4 ind-sp, 7 sp-sp) “ind_08_01” “ind_12_01” “ind_12_06”
“ind_12_07” “sp_03_01” “sp_03_02” “sp_03_03” “sp_03_04” “sp_05_01”
“sp_31_01” “sp_44_01”

``` r
net.level.selection %<>% 
  arrange(type, net_id) %>%
  mutate(net_n = order(net_code)) %>%
  mutate(study = str_sub(net_id, 1,2),
         study = paste0(type, "_", study))
```

    ## mutate: changed 88 values (84%) of 'net_n' (0 new NA)

    ## mutate: new variable 'study' (character) with 61 unique values and 0% NA

``` r
glimpse(net.level.selection)
```

    ## Rows: 105
    ## Columns: 50
    ## $ connectance                     <dbl> 0.3629630, 0.2093750, 0.3914286, 0.4…
    ## $ web.asymmetry                   <dbl> -0.19402985, -0.42857143, -0.5555555…
    ## $ links.per.species               <dbl> 5.850746, 2.392857, 3.044444, 3.4222…
    ## $ cluster.coefficient             <dbl> 0.22500000, 0.06250000, 0.25714286, …
    ## $ NODF                            <dbl> 67.21158, 63.51606, 72.48350, 71.689…
    ## $ weighted.NODF                   <dbl> 52.79227, 43.93032, 37.73835, 38.511…
    ## $ M                               <dbl> 0.15799305, 0.23921401, 0.09040128, …
    ## $ interaction.strength.asymmetry  <dbl> 0.098096343, 0.303002612, 1.19449172…
    ## $ specialisation.asymmetry        <dbl> 1.00000000, -0.48846320, -0.82692306…
    ## $ linkage.density                 <dbl> 16.046511, 9.152386, 10.990321, 13.9…
    ## $ weighted.connectance            <dbl> 0.2395002, 0.1634355, 0.2442293, 0.3…
    ## $ Shannon.diversity               <dbl> 4.933505, 3.959126, 3.653736, 4.1713…
    ## $ interaction.evenness            <dbl> 0.7063286, 0.6127286, 0.6237244, 0.7…
    ## $ Alatalo.interaction.evenness    <dbl> 0.6550250, 0.5924512, 0.5490386, 0.6…
    ## $ H2                              <dbl> 0.09933567, 0.24470407, 0.17791608, …
    ## $ number.of.species.HL            <int> 27, 16, 10, 10, 11, 10, 10, 13, 12, …
    ## $ number.of.species.LL            <int> 40, 40, 35, 35, 35, 13, 14, 14, 13, …
    ## $ cluster.coefficient.HL          <dbl> 0.9360087, 0.7633714, 0.9610723, 0.9…
    ## $ cluster.coefficient.LL          <dbl> 0.4297705, 0.2429504, 0.3870170, 0.4…
    ## $ weighted.cluster.coefficient.HL <dbl> 1.00000000, 0.75836703, 0.82419380, …
    ## $ weighted.cluster.coefficient.LL <dbl> 0.9675836, 0.8203872, 0.9500581, 0.9…
    ## $ niche.overlap.horn.HL           <dbl> 0.28227633, 0.09229992, 0.14771050, …
    ## $ niche.overlap.horn.LL           <dbl> 0.8377570, 0.6984493, 0.9032446, 0.9…
    ## $ extinction.slope.HL             <dbl> 10.730740, 3.721760, 4.293068, 5.561…
    ## $ extinction.slope.LL             <dbl> 7.822940, 2.596134, 5.649524, 6.7967…
    ## $ robustness.HL                   <dbl> 0.9145234, 0.7824115, 0.8078489, 0.8…
    ## $ robustness.LL                   <dbl> 0.8599260, 0.7138990, 0.8181298, 0.8…
    ## $ generality.HL                   <dbl> 27.560125, 15.435135, 20.212817, 25.…
    ## $ vulnerability.LL                <dbl> 4.532896, 2.869637, 1.767824, 2.3798…
    ## $ niche.overlap.bray.HL           <dbl> 0.09729529, 0.03883194, 0.07298345, …
    ## $ niche.overlap.bray.LL           <dbl> 0.5567506, 0.3390878, 0.4458126, 0.5…
    ## $ niche.overlap.jaccard.HL        <dbl> 0.05728874, 0.02443975, 0.04183463, …
    ## $ niche.overlap.jaccard.LL        <dbl> 0.40377235, 0.22904359, 0.31947895, …
    ## $ net_size                        <int> 1080, 640, 350, 350, 385, 130, 140, …
    ## $ net_id                          <chr> "01_01", "01_02", "02_01", "02_02", …
    ## $ links                           <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, …
    ## $ net_n                           <int> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 1…
    ## $ type                            <chr> "ind", "ind", "ind", "ind", "ind", "…
    ## $ net_code                        <chr> "ind_01_01", "ind_01_02", "ind_02_01…
    ## $ assortativity                   <dbl> -0.5071021, -0.6033629, -0.7374228, …
    ## $ centr_binary                    <dbl> 0.6210147, 0.7887529, 0.7122378, 0.6…
    ## $ centralization.w                <dbl> 0.8549594, 0.9221565, 0.9284492, 0.8…
    ## $ unique.ints                     <dbl> 392, 134, 137, 154, 148, 37, 33, 46,…
    ## $ code_ID                         <chr> "01_01", "01_02", "02_01", "02_02", …
    ## $ ref                             <chr> "Quintero et al 2023 Ecology letters…
    ## $ bioregion                       <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA, …
    ## $ country                         <chr> "Spain", "Spain", "Spain", "Spain", …
    ## $ plant_sp                        <chr> "Pistacia lentiscus", "Pistacia lent…
    ## $ n_nodes                         <int> 67, 56, 45, 45, 46, 23, 24, 27, 25, …
    ## $ study                           <chr> "ind_01", "ind_01", "ind_02", "ind_0…

``` r
write_csv(net.level.selection, here("data/net.level.selection.csv"))
```

Tables:

``` r
library(knitr)
library(kableExtra)
```

    ## 
    ## Attaching package: 'kableExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     group_rows

``` r
table <- net.level.selection %>%
  dplyr::select(net_n, type, net_code, ref, net_size,
                n_plants = number.of.species.LL,
                n_frugivores = number.of.species.HL) %>%
  kable() %>% 
  kable_classic(font_size = 12, full_width = F)

table
```

<table class=" lightable-classic" style="font-size: 12px; font-family: &quot;Arial Narrow&quot;, &quot;Source Sans Pro&quot;, sans-serif; width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
net_n
</th>
<th style="text-align:left;">
type
</th>
<th style="text-align:left;">
net_code
</th>
<th style="text-align:left;">
ref
</th>
<th style="text-align:right;">
net_size
</th>
<th style="text-align:right;">
n_plants
</th>
<th style="text-align:right;">
n_frugivores
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_01_01
</td>
<td style="text-align:left;">
Quintero et al 2023 Ecology letters
</td>
<td style="text-align:right;">
1080
</td>
<td style="text-align:right;">
40
</td>
<td style="text-align:right;">
27
</td>
</tr>
<tr>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_01_02
</td>
<td style="text-align:left;">
Quintero et al 2023 Ecology letters
</td>
<td style="text-align:right;">
640
</td>
<td style="text-align:right;">
40
</td>
<td style="text-align:right;">
16
</td>
</tr>
<tr>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_02_01
</td>
<td style="text-align:left;">
Isla et al 2023 PRSB
</td>
<td style="text-align:right;">
350
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:right;">
4
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_02_02
</td>
<td style="text-align:left;">
Isla et al 2023 PRSB
</td>
<td style="text-align:right;">
350
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:right;">
5
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_02_03
</td>
<td style="text-align:left;">
Isla et al 2023 PRSB
</td>
<td style="text-align:right;">
385
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
11
</td>
</tr>
<tr>
<td style="text-align:right;">
6
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_03_01
</td>
<td style="text-align:left;">
Vergara-Tabares et al 2022 Oikos
</td>
<td style="text-align:right;">
130
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:right;">
7
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_03_02
</td>
<td style="text-align:left;">
Vergara-Tabares et al 2022 Oikos
</td>
<td style="text-align:right;">
140
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:right;">
8
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_03_03
</td>
<td style="text-align:left;">
Vergara-Tabares et al 2022 Oikos
</td>
<td style="text-align:right;">
182
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
13
</td>
</tr>
<tr>
<td style="text-align:right;">
9
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_03_04
</td>
<td style="text-align:left;">
Vergara-Tabares et al 2022 Oikos
</td>
<td style="text-align:right;">
156
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
12
</td>
</tr>
<tr>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_03_05
</td>
<td style="text-align:left;">
Vergara-Tabares et al 2022 Oikos
</td>
<td style="text-align:right;">
132
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
11
</td>
</tr>
<tr>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_03_06
</td>
<td style="text-align:left;">
Vergara-Tabares et al 2022 Oikos
</td>
<td style="text-align:right;">
77
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:right;">
12
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_04_01
</td>
<td style="text-align:left;">
Rodriguez-Sanchez 2010 PhD thesis
</td>
<td style="text-align:right;">
306
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
17
</td>
</tr>
<tr>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_05_01
</td>
<td style="text-align:left;">
Jordano 1995 Ecology
</td>
<td style="text-align:right;">
380
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
20
</td>
</tr>
<tr>
<td style="text-align:right;">
14
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_06_01
</td>
<td style="text-align:left;">
Friedemann et al 2022 Oikos
</td>
<td style="text-align:right;">
153
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:right;">
15
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_06_02
</td>
<td style="text-align:left;">
Friedemann et al 2022 Oikos
</td>
<td style="text-align:right;">
105
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:right;">
16
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_06_03
</td>
<td style="text-align:left;">
Friedemann et al 2022 Oikos
</td>
<td style="text-align:right;">
240
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:right;">
17
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_07_01
</td>
<td style="text-align:left;">
Cecropia Frugivory course
</td>
<td style="text-align:right;">
999
</td>
<td style="text-align:right;">
27
</td>
<td style="text-align:right;">
37
</td>
</tr>
<tr>
<td style="text-align:right;">
18
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_08_02
</td>
<td style="text-align:left;">
Gopal et al 2020 Biotropica
</td>
<td style="text-align:right;">
264
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
11
</td>
</tr>
<tr>
<td style="text-align:right;">
19
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_08_03
</td>
<td style="text-align:left;">
Gopal et al 2020 Biotropica
</td>
<td style="text-align:right;">
175
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:right;">
20
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_08_04
</td>
<td style="text-align:left;">
Gopal et al 2020 Biotropica
</td>
<td style="text-align:right;">
672
</td>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
21
</td>
</tr>
<tr>
<td style="text-align:right;">
21
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_09_01
</td>
<td style="text-align:left;">
Crestani et al 2018 Acta Oecologica
</td>
<td style="text-align:right;">
396
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
22
</td>
</tr>
<tr>
<td style="text-align:right;">
22
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_10_01
</td>
<td style="text-align:left;">
Lamperty et al 2021 Biotropica
</td>
<td style="text-align:right;">
279
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:right;">
23
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_11_01
</td>
<td style="text-align:left;">
Corema album - unpublished
</td>
<td style="text-align:right;">
360
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
15
</td>
</tr>
<tr>
<td style="text-align:right;">
24
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_12_02
</td>
<td style="text-align:left;">
Ramaswami et al 2017 Plant Ecology
</td>
<td style="text-align:right;">
165
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
11
</td>
</tr>
<tr>
<td style="text-align:right;">
25
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_12_03
</td>
<td style="text-align:left;">
Ramaswami et al 2017 Plant Ecology
</td>
<td style="text-align:right;">
72
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:right;">
26
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_12_04
</td>
<td style="text-align:left;">
Ramaswami et al 2017 Plant Ecology
</td>
<td style="text-align:right;">
65
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:right;">
27
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_12_05
</td>
<td style="text-align:left;">
Ramaswami et al 2017 Plant Ecology
</td>
<td style="text-align:right;">
80
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:right;">
28
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_12_08
</td>
<td style="text-align:left;">
Ramaswami et al 2017 Plant Ecology
</td>
<td style="text-align:right;">
80
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:right;">
29
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_12_09
</td>
<td style="text-align:left;">
Ramaswami et al 2017 Plant Ecology
</td>
<td style="text-align:right;">
140
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:right;">
30
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_12_10
</td>
<td style="text-align:left;">
Ramaswami et al 2017 Plant Ecology
</td>
<td style="text-align:right;">
195
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
13
</td>
</tr>
<tr>
<td style="text-align:right;">
31
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_13_01
</td>
<td style="text-align:left;">
Jácome-Flores et al 2019 Oikos
</td>
<td style="text-align:right;">
234
</td>
<td style="text-align:right;">
39
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:right;">
32
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_13_02
</td>
<td style="text-align:left;">
Jácome-Flores et al 2019 Oikos
</td>
<td style="text-align:right;">
144
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:right;">
33
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_14_01
</td>
<td style="text-align:left;">
Guerra et al 2017 Oecología
</td>
<td style="text-align:right;">
135
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:right;">
34
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_15_01
</td>
<td style="text-align:left;">
Juniperus macrocarpa - unpublished
</td>
<td style="text-align:right;">
286
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
11
</td>
</tr>
<tr>
<td style="text-align:right;">
35
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_16_01
</td>
<td style="text-align:left;">
Miguel et al 2018 Oikos
</td>
<td style="text-align:right;">
234
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:right;">
36
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_16_02
</td>
<td style="text-align:left;">
Miguel et al 2018 Oikos
</td>
<td style="text-align:right;">
280
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:right;">
37
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_16_03
</td>
<td style="text-align:left;">
Miguel et al 2018 Oikos
</td>
<td style="text-align:right;">
540
</td>
<td style="text-align:right;">
54
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:right;">
38
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_16_04
</td>
<td style="text-align:left;">
Miguel et al 2018 Oikos
</td>
<td style="text-align:right;">
280
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:right;">
39
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_16_05
</td>
<td style="text-align:left;">
Miguel et al 2018 Oikos
</td>
<td style="text-align:right;">
261
</td>
<td style="text-align:right;">
29
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:right;">
40
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_17_01
</td>
<td style="text-align:left;">
Visotto et al 2022 Oecologia
</td>
<td style="text-align:right;">
416
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
16
</td>
</tr>
<tr>
<td style="text-align:right;">
41
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_18_01
</td>
<td style="text-align:left;">
Phyllirea angustifolia - unpublished
</td>
<td style="text-align:right;">
160
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
16
</td>
</tr>
<tr>
<td style="text-align:right;">
42
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_18_02
</td>
<td style="text-align:left;">
Phyllirea angustifolia - unpublished
</td>
<td style="text-align:right;">
108
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
12
</td>
</tr>
<tr>
<td style="text-align:right;">
43
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_19_01
</td>
<td style="text-align:left;">
Thiel et al 2023 Biotropica
</td>
<td style="text-align:right;">
1032
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
43
</td>
</tr>
<tr>
<td style="text-align:right;">
44
</td>
<td style="text-align:left;">
ind
</td>
<td style="text-align:left;">
ind_20_01
</td>
<td style="text-align:left;">
Osyris lanceolata - unpublished
</td>
<td style="text-align:right;">
266
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
14
</td>
</tr>
<tr>
<td style="text-align:right;">
45
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_01_01
</td>
<td style="text-align:left;">
Olesen et al 2011 PRSB
</td>
<td style="text-align:right;">
900
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
36
</td>
</tr>
<tr>
<td style="text-align:right;">
46
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_01_02
</td>
<td style="text-align:left;">
Olesen et al 2011 PRSB
</td>
<td style="text-align:right;">
272
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
17
</td>
</tr>
<tr>
<td style="text-align:right;">
47
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_02_01
</td>
<td style="text-align:left;">
Garcia-Castaño et al 2011 PhD thesis
</td>
<td style="text-align:right;">
476
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
28
</td>
</tr>
<tr>
<td style="text-align:right;">
48
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_04_01
</td>
<td style="text-align:left;">
Beehler 1983 Auk
</td>
<td style="text-align:right;">
279
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:right;">
49
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_06_01
</td>
<td style="text-align:left;">
Frost 1980 Acta XVII Cong. Internat. Ornitho.
</td>
<td style="text-align:right;">
160
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:right;">
50
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_07_01
</td>
<td style="text-align:left;">
Guitián 1983 PhD thesis
</td>
<td style="text-align:right;">
84
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:right;">
51
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_08_01
</td>
<td style="text-align:left;">
Galetti & Pizo 1996 Ararajuba
</td>
<td style="text-align:right;">
1015
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
29
</td>
</tr>
<tr>
<td style="text-align:right;">
52
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_09_01
</td>
<td style="text-align:left;">
Kantak 1979 Auk
</td>
<td style="text-align:right;">
135
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
27
</td>
</tr>
<tr>
<td style="text-align:right;">
53
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_10_01
</td>
<td style="text-align:left;">
Snow & Snow 1988 Bird and berries
</td>
<td style="text-align:right;">
154
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
14
</td>
</tr>
<tr>
<td style="text-align:right;">
54
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_11_01
</td>
<td style="text-align:left;">
Noma 1997 Ecological Research
</td>
<td style="text-align:right;">
120
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:right;">
55
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_12_01
</td>
<td style="text-align:left;">
Crome 1975 Aust Wildl Res
</td>
<td style="text-align:right;">
497
</td>
<td style="text-align:right;">
71
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:right;">
56
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_13_01
</td>
<td style="text-align:left;">
Snow & Snow 1971 Auk
</td>
<td style="text-align:right;">
700
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
14
</td>
</tr>
<tr>
<td style="text-align:right;">
57
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_14_01
</td>
<td style="text-align:left;">
Baird 1980 Wilson Bulletin
</td>
<td style="text-align:right;">
147
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
21
</td>
</tr>
<tr>
<td style="text-align:right;">
58
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_15_01
</td>
<td style="text-align:left;">
Menke et al 2012 Oikos
</td>
<td style="text-align:right;">
240
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
30
</td>
</tr>
<tr>
<td style="text-align:right;">
59
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_15_02
</td>
<td style="text-align:left;">
Menke et al 2012 Oikos
</td>
<td style="text-align:right;">
266
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
38
</td>
</tr>
<tr>
<td style="text-align:right;">
60
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_15_03
</td>
<td style="text-align:left;">
Menke et al 2012 Oikos
</td>
<td style="text-align:right;">
272
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
34
</td>
</tr>
<tr>
<td style="text-align:right;">
61
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_15_04
</td>
<td style="text-align:left;">
Menke et al 2012 Oikos
</td>
<td style="text-align:right;">
312
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
39
</td>
</tr>
<tr>
<td style="text-align:right;">
62
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_16_01
</td>
<td style="text-align:left;">
Pizo 2004 Ornitologia Neotropical
</td>
<td style="text-align:right;">
735
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
49
</td>
</tr>
<tr>
<td style="text-align:right;">
63
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_17_01
</td>
<td style="text-align:left;">
Schleuning et al 2011 Ecology
</td>
<td style="text-align:right;">
2904
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
88
</td>
</tr>
<tr>
<td style="text-align:right;">
64
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_18_01
</td>
<td style="text-align:left;">
Castro 2007 PhD thesis
</td>
<td style="text-align:right;">
784
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
16
</td>
</tr>
<tr>
<td style="text-align:right;">
65
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_19_01
</td>
<td style="text-align:left;">
Correia 1997 PhD thesis
</td>
<td style="text-align:right;">
585
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
45
</td>
</tr>
<tr>
<td style="text-align:right;">
66
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_20_01
</td>
<td style="text-align:left;">
Alves 2008 PhD thesis
</td>
<td style="text-align:right;">
390
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
30
</td>
</tr>
<tr>
<td style="text-align:right;">
67
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_21_01
</td>
<td style="text-align:left;">
Athie 2009 PhD thesis
</td>
<td style="text-align:right;">
270
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
30
</td>
</tr>
<tr>
<td style="text-align:right;">
68
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_22_01
</td>
<td style="text-align:left;">
Fadini & de Marco 2004 Ararajuba
</td>
<td style="text-align:right;">
700
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
28
</td>
</tr>
<tr>
<td style="text-align:right;">
69
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_23_01
</td>
<td style="text-align:left;">
Hasui 1994 MSc thesis
</td>
<td style="text-align:right;">
572
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
22
</td>
</tr>
<tr>
<td style="text-align:right;">
70
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_24_01
</td>
<td style="text-align:left;">
Silva 2011 MSc thesis
</td>
<td style="text-align:right;">
440
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
20
</td>
</tr>
<tr>
<td style="text-align:right;">
71
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_25_01
</td>
<td style="text-align:left;">
Ribeiro da Silva et al 2015 Restoration Ecology
</td>
<td style="text-align:right;">
180
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
15
</td>
</tr>
<tr>
<td style="text-align:right;">
72
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_25_02
</td>
<td style="text-align:left;">
Ribeiro da Silva et al 2015 Restoration Ecology
</td>
<td style="text-align:right;">
667
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
29
</td>
</tr>
<tr>
<td style="text-align:right;">
73
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_25_03
</td>
<td style="text-align:left;">
Ribeiro da Silva et al 2015 Restoration Ecology
</td>
<td style="text-align:right;">
196
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
14
</td>
</tr>
<tr>
<td style="text-align:right;">
74
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_26_01
</td>
<td style="text-align:left;">
Robinson 2015 PhD thesis
</td>
<td style="text-align:right;">
168
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
28
</td>
</tr>
<tr>
<td style="text-align:right;">
75
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_27_01
</td>
<td style="text-align:left;">
Rodrigues 2015 PhD thesis
</td>
<td style="text-align:right;">
1740
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
58
</td>
</tr>
<tr>
<td style="text-align:right;">
76
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_28_01
</td>
<td style="text-align:left;">
Burns 2013 Ecology
</td>
<td style="text-align:right;">
144
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:right;">
77
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_29_01
</td>
<td style="text-align:left;">
Albrecht et al 2015 J of Ecology
</td>
<td style="text-align:right;">
96
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
12
</td>
</tr>
<tr>
<td style="text-align:right;">
78
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_29_02
</td>
<td style="text-align:left;">
Albrecht et al 2015 J of Ecology
</td>
<td style="text-align:right;">
77
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
11
</td>
</tr>
<tr>
<td style="text-align:right;">
79
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_29_03
</td>
<td style="text-align:left;">
Albrecht et al 2015 J of Ecology
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
13
</td>
</tr>
<tr>
<td style="text-align:right;">
80
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_29_04
</td>
<td style="text-align:left;">
Albrecht et al 2015 J of Ecology
</td>
<td style="text-align:right;">
104
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
13
</td>
</tr>
<tr>
<td style="text-align:right;">
81
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_29_05
</td>
<td style="text-align:left;">
Albrecht et al 2015 J of Ecology
</td>
<td style="text-align:right;">
80
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:right;">
82
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_29_06
</td>
<td style="text-align:left;">
Albrecht et al 2015 J of Ecology
</td>
<td style="text-align:right;">
80
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:right;">
83
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_29_07
</td>
<td style="text-align:left;">
Albrecht et al 2015 J of Ecology
</td>
<td style="text-align:right;">
88
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
11
</td>
</tr>
<tr>
<td style="text-align:right;">
84
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_29_08
</td>
<td style="text-align:left;">
Albrecht et al 2015 J of Ecology
</td>
<td style="text-align:right;">
130
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
13
</td>
</tr>
<tr>
<td style="text-align:right;">
85
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_29_09
</td>
<td style="text-align:left;">
Albrecht et al 2015 J of Ecology
</td>
<td style="text-align:right;">
120
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
15
</td>
</tr>
<tr>
<td style="text-align:right;">
86
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_29_10
</td>
<td style="text-align:left;">
Albrecht et al 2015 J of Ecology
</td>
<td style="text-align:right;">
144
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
16
</td>
</tr>
<tr>
<td style="text-align:right;">
87
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_29_11
</td>
<td style="text-align:left;">
Albrecht et al 2015 J of Ecology
</td>
<td style="text-align:right;">
120
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
20
</td>
</tr>
<tr>
<td style="text-align:right;">
88
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_29_12
</td>
<td style="text-align:left;">
Albrecht et al 2015 J of Ecology
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
11
</td>
</tr>
<tr>
<td style="text-align:right;">
89
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_29_13
</td>
<td style="text-align:left;">
Albrecht et al 2015 J of Ecology
</td>
<td style="text-align:right;">
152
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
19
</td>
</tr>
<tr>
<td style="text-align:right;">
90
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_30_01
</td>
<td style="text-align:left;">
Andrade et al 2011 Rev. Brasileira de Ornitologia
</td>
<td style="text-align:right;">
374
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
17
</td>
</tr>
<tr>
<td style="text-align:right;">
91
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_32_01
</td>
<td style="text-align:left;">
Yang et al 2013 Ecosphere
</td>
<td style="text-align:right;">
680
</td>
<td style="text-align:right;">
34
</td>
<td style="text-align:right;">
20
</td>
</tr>
<tr>
<td style="text-align:right;">
92
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_33_01
</td>
<td style="text-align:left;">
Garcia et al 2000 Rev. Biología Tropical
</td>
<td style="text-align:right;">
84
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:right;">
93
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_34_01
</td>
<td style="text-align:left;">
Gorchov et al 1995 Oikos
</td>
<td style="text-align:right;">
1386
</td>
<td style="text-align:right;">
77
</td>
<td style="text-align:right;">
18
</td>
</tr>
<tr>
<td style="text-align:right;">
94
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_35_01
</td>
<td style="text-align:left;">
Palmeirim et al 1989 Oecologia
</td>
<td style="text-align:right;">
490
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
14
</td>
</tr>
<tr>
<td style="text-align:right;">
95
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_36_01
</td>
<td style="text-align:left;">
Lopez & Vaughan 2004 Acta Chiropterologica
</td>
<td style="text-align:right;">
490
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
14
</td>
</tr>
<tr>
<td style="text-align:right;">
96
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_37_01
</td>
<td style="text-align:left;">
Heleno et al 2013 Proc Biol Sci
</td>
<td style="text-align:right;">
645
</td>
<td style="text-align:right;">
43
</td>
<td style="text-align:right;">
15
</td>
</tr>
<tr>
<td style="text-align:right;">
97
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_38_01
</td>
<td style="text-align:left;">
Hernandez-Montero et al 2015 PloS ONE
</td>
<td style="text-align:right;">
154
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:right;">
98
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_38_02
</td>
<td style="text-align:left;">
Hernandez-Montero et al 2015 PloS ONE
</td>
<td style="text-align:right;">
114
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:right;">
99
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_39_01
</td>
<td style="text-align:left;">
Passos et al 2003 Rev. Brasileira de Zoologia
</td>
<td style="text-align:right;">
168
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:right;">
100
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_40_01
</td>
<td style="text-align:left;">
Pedro 1992 PhD thesis
</td>
<td style="text-align:right;">
91
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:right;">
101
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_41_01
</td>
<td style="text-align:left;">
Poulin et al 1999 J Tropical Ecol
</td>
<td style="text-align:right;">
340
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
20
</td>
</tr>
<tr>
<td style="text-align:right;">
102
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_42_01
</td>
<td style="text-align:left;">
Sarmento et al 2014 Zoologia
</td>
<td style="text-align:right;">
1120
</td>
<td style="text-align:right;">
56
</td>
<td style="text-align:right;">
20
</td>
</tr>
<tr>
<td style="text-align:right;">
103
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_43_01
</td>
<td style="text-align:left;">
Stiebel & Barlein 2008 Vogelwarte
</td>
<td style="text-align:right;">
930
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
31
</td>
</tr>
<tr>
<td style="text-align:right;">
104
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_45_01
</td>
<td style="text-align:left;">
Saavedra et al 2014 Oecologia
</td>
<td style="text-align:right;">
1476
</td>
<td style="text-align:right;">
36
</td>
<td style="text-align:right;">
41
</td>
</tr>
<tr>
<td style="text-align:right;">
105
</td>
<td style="text-align:left;">
sp
</td>
<td style="text-align:left;">
sp_45_02
</td>
<td style="text-align:left;">
Saavedra et al 2014 Oecologia
</td>
<td style="text-align:right;">
460
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
23
</td>
</tr>
</tbody>
</table>

``` r
# save_kable(table, "list_nets.html")
# webshot::webshot("list_nets.html", "list_nets.pdf")
```

# 2. Node level

``` r
ind.level.df <- read.csv(here("data/ind.level.df.csv"))
```

Load ind networks characteristics:

``` r
ind_nets_char <- readxl::read_xlsx(here("data/selected_individual_networks.xlsx"))

ind_nets_ID <- ind_nets_char %>%
  dplyr::select(net_id = code_ID,
         ref,
         plant_sp = `focus species`,
         pop = pop_site,
         country) %>%
  mutate(net_name = paste0(plant_sp, "_", pop))
```

    ## mutate: new variable 'net_name' (character) with 48 unique values and 0% NA

Load net colors:

``` r
net_cols <- read_csv(here("data/net_colors.csv")) %>%
  dplyr::select(continent, country, plant_sp, plant_plot_rank, 
                cols_continent3, cols_continent4)
```

Add network characteristics, colors and filter small nets as above:

``` r
ind.level.selec <- ind.level.df %>%
  mutate(study = str_sub(net_id, 1,2)) %>%
  left_join(ind_nets_ID) %>%
  left_join(net_cols) %>%
  filter(!net_id %in% c("08_01", "12_01", "12_06", "12_07")) %>%
  group_by(net_id) %>%
  mutate(net_n = cur_group_id()) %>%
  mutate(ind_ID = paste0(net_id, "_", ind)) %>% #create a plant id unique name
  relocate(net_n, .before = everything()) %>%
  relocate(net_id, .after = net_n) %>%
  relocate(study, .after = net_id) %>%
  relocate(ind_ID, .after = study)
```

    ## mutate: new variable 'study' (character) with 20 unique values and 0% NA

    ## Joining, by = "net_id"
    ## left_join: added 5 columns (ref, plant_sp, pop, country, net_name)
    ## > rows only in x 0
    ## > rows only in y ( 0)
    ## > matched rows 1,022
    ## > =======
    ## > rows total 1,022
    ## Joining, by = c("plant_sp", "country")
    ## left_join: added 4 columns (continent, plant_plot_rank, cols_continent3,
    ## cols_continent4)
    ## > rows only in x 27
    ## > rows only in y ( 0)
    ## > matched rows 995
    ## > =======
    ## > rows total 1,022
    ## filter: removed 27 rows (3%), 995 rows remaining
    ## group_by: one grouping variable (net_id)
    ## mutate (grouped): changed 607 values (61%) of 'net_n' (0 new NA)
    ## mutate (grouped): no changes
    ## relocate: columns reordered (net_n, ind, degree, normalised.degree,
    ## species.strength, …)
    ## relocate: columns reordered (net_n, net_id, ind, degree, normalised.degree,
    ## …)
    ## relocate: columns reordered (net_n, net_id, study, ind, degree, …)
    ## relocate: columns reordered (net_n, net_id, study, ind_ID, ind, …)

``` r
glimpse(ind.level.selec)
```

    ## Rows: 995
    ## Columns: 41
    ## Groups: net_id [44]
    ## $ net_n                         <int> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
    ## $ net_id                        <chr> "01_01", "01_01", "01_01", "01_01", "0…
    ## $ study                         <chr> "01", "01", "01", "01", "01", "01", "0…
    ## $ ind_ID                        <chr> "01_01_301", "01_01_302", "01_01_303",…
    ## $ ind                           <chr> "301", "302", "303", "304", "305", "30…
    ## $ degree                        <int> 11, 15, 9, 19, 5, 11, 7, 7, 9, 7, 8, 7…
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
    ## $ mean.bray.overlap             <dbl> 0.5754007, 0.5782921, 0.5968529, 0.415…
    ## $ mean.jaccard.overlap          <dbl> 0.4182767, 0.4215563, 0.4405808, 0.283…
    ## $ eigen.centrality              <dbl> 0.24815105, 0.22843692, 0.11912347, 0.…
    ## $ tm.closeness                  <dbl> 0.03039829, 0.03329066, 0.02345089, 0.…
    ## $ tm.w.closeness                <dbl> 0.031202628, 0.031920244, 0.014090767,…
    ## $ tm.betweenness                <dbl> 0.0, 0.5, 0.0, 149.5, 0.0, 0.0, 0.0, 0…
    ## $ tm.w.betweenness              <dbl> 0, 0, 0, 89, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ ref                           <chr> "Quintero et al 2023 Ecology letters",…
    ## $ plant_sp                      <chr> "Pistacia lentiscus", "Pistacia lentis…
    ## $ pop                           <chr> "El Puntal", "El Puntal", "El Puntal",…
    ## $ country                       <chr> "Spain", "Spain", "Spain", "Spain", "S…
    ## $ net_name                      <chr> "Pistacia lentiscus_El Puntal", "Pista…
    ## $ continent                     <chr> "Europe", "Europe", "Europe", "Europe"…
    ## $ plant_plot_rank               <dbl> 27, 27, 27, 27, 27, 27, 27, 27, 27, 27…
    ## $ cols_continent3               <chr> "#FF9B00", "#FF9B00", "#FF9B00", "#FF9…
    ## $ cols_continent4               <chr> "#FFCA28FF", "#FFCA28FF", "#FFCA28FF",…

``` r
write_csv(ind.level.selec, here("data/node.level.selection.csv"))
```
