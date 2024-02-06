#### Master script for the project ####

library(rmarkdown)

#### DATA PREPARATION #####

# Each network has been standardized dividing each interaction by the total number
# of interactions (Grand Total Standardization). All standardized interactions matrices
# are in ("networks/standardized").
render("analysis/webs_standarization.Rmd")

#Explore sampling effort per individual:
render("analysis/sampling_effort.Rmd")

# To obtain main web metrics:
render("analysis/mod_estimation.Rmd")
render("analysis/web_metrics.Rmd")
render("analysis/node_metrics.Rmd")

# Selection of webs for analysis:
render("analysis/web_selection.Rmd") #filter by non-representative nets and include network characteristics

# To asses sampling completeness:
render("analysis/sampling_completeness.Rmd")

# Selection of frugivore attributes for analysis:
render("analysis/prep_frug_attributes.Rmd")

# Calculate frugivore contribution to overall interactions:
render("analysis/frugivore_contribution.Rmd")

# Metrics selection:
render("analysis/selection_net_level.Rmd") # Explore net-level metrics
render("analysis/selection_node_level.Rmd") # Explore node-level metrics

# PCA analysis of NET-LEVEL and NODE-LEVEL metrics:
render("analysis/net_level_comparison.Rmd")
render("analysis/node_level_analysis.Rmd")
render("analysis/PCA_3d_plot.Rmd")

# Analysis individual specialization:
render("analysis/ind_specialization.Rmd")



