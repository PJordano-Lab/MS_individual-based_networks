library(dplyr)

Shannon <- function(x) {x*log(x)}  
# Note this doesn't calculate Shannon index, just part of it


## TNW

tnw <- function(net) {
  net |> 
    summarise(across(everything(), sum)) |> 
    mutate(across(everything(), Shannon)) |> 
    mutate(TNW = -rowSums(across(everything()), na.rm = TRUE)) |> 
    select(TNW)
}



## WIC

wic <- function(net) {
  net |> 
    rowwise() |> 
    mutate(plant.sum = rowSums(across(everything()))) |> 
    mutate(across(!matches("plant.sum"), function(x, y = plant.sum) {x/y})) |> 
    mutate(across(!matches("plant.sum"), Shannon)) |> 
    mutate(shannon = -rowSums(across(!matches("plant.sum")), na.rm = TRUE)) |> 
    mutate(shannon.prop = plant.sum * shannon) |> 
    ungroup() |> 
    summarise(wic = sum(shannon.prop))
}

