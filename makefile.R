#### Master script for the project ####

library(rmarkdown)

#### DATA PREPARATION #####

# Control networks by sampling effort and standardize networks
# (both for species-based and individual-based)

render("analysis/nets_standarization.Rmd")

# FOR INDIVIDUAL-BASED NETWORKS:

# Bayesian modelling of interaction counts

nets <- list.files(here::here("networks/nets_raw"), pattern = "int.csv")
nets <- sort(gsub("_int.csv", "", nets))
params <- list(beta = 0.01, iter = 10000, model = "varying_preferences")
for (net in nets) {
  params$net <- net
  if (net %in% c("11_01", "15_01", "16_03", "16_04", "19_01")) {
    params$iter = 50000
  }
  if (net %in% c("02_01")) {
    params$iter = 100000
  }
  render(input = "analysis/Bayes_varying_pref.Rmd",
         output_file = net,
         output_dir = "analysis/output/Modelling/",
         knit_root_dir = here::here(),
         params = params)
}

# Control networks posteriors by sampling effort and standardize networks

render("analysis/nets_standarization_Bayes.Rmd")


#####

# Asses sampling completeness

render("analysis/sampling_completeness.Rmd")

# Obtain network-level metrics for species-based and individual-based nets

render("analysis/net_level_metrics.Rmd")

# Null models for network-level metrics for species-based and individual-based nets

render("analysis/net_level_null_models.Rmd")

# Obtain web metrics for individual-based network posteriors

render("analysis/net_metrics_posteriors.Rmd")

# Selection of frugivore attributes for analysis

render("analysis/prep_frug_attributes.Rmd")

# Calculate frugivore contribution to overall interactions (Fig. 4)

render("analysis/frugivore_contribution.Rmd")

# Selection of network metrics

render("analysis/selection_net_level.Rmd") # Explore net-level metrics
render("analysis/selection_node_level.Rmd") # Explore node-level metrics

# Compare individual-based and sp-based networks (Fig. 2)

render("analysis/net_level_comparison.Rmd")

# Calculate Individual specialization (Fig. 3)

render("analysis/WIC_TNW_calc.Rmd")

# Node-level analyses for individual-based networks (Fig. 5)

render("analysis/PCA_cluster_analysis.Rmd")
render("analysis/PCA_plant_individuals.Rmd")
