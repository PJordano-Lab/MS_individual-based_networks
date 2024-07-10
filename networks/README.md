# Plant-frugivore interaction networks compilation

For this study, we compiled studies on frugivory ecological networks
with publicly-available data, both at the species and the individual
plant level. **Species-based networks** were gathered from published
studies at community scale. **Individual-based networks**, which are
scarcer, were compiled from phyto-centric studies (plant-based). We
combined published studies with unpublished data-sets, gathering data
for 21 different study systems, including data-sets from our own field
studies with different Mediterranean species (n = 9).

Please note that some of the studies selected present more than one
network from different communities (in species-based studies) or
populations (in individual-based studies), see file
[list_nets_full_ref.pdf](list_nets_full_ref.pdf) for the origin and
basic information of each network. Data is entered as adjacency
matrices, where rows represent plant species (or individuals) and
columns represent animal species. When possible we referred the
interaction value to the coarsest level, that is, frugivore visitation
events, otherwise number of fruits consumed. For analytic purposes, we
discarded networks in which the number of interacting nodes (plants and
frugivore species) was less than 15 or plants were less than six (n = 13
networks). The final data-set consists of 105 networks, of which 46 are
individual-based networks and 59 are species-based networks

### Folder estructure

This folder contains six sub-folders with the following information:

1)  **nets_raw** - contains raw interaction counts from the original
    individual-based network studies.

2)  **nets_se** - contains sampling effort information on each plant
    individual for each individual-based network studies.

3)  **nets_sp_raw** - contains interaction counts from original
    species-based network studies.

4)  **nets_std** - contains both individual-based and species-based
    networks controlled by sampling effort when needed (in the case of
    some ind-based nets, n = 21) and standardized to proportions, so
    that each interaction represents a proportion from the total
    interaction count, therefore the sum of all interactions equals 1.
    These nets are generated as output from the code
    “nets_standarization.Rmd”.

5)  **nets_post** - contains 1000 posterior distribution counts for each
    individual-based network reconstructed using a Bayesian-framework
    that accounts for differential sampling effort and interaction
    uncertainty. These nets are generated as output from the code
    “Bayes_varying_pref.Rmd”.

6)  **nets_post_std** - contains individual-based networks posteriors
    controlled by sampling effort and standardized to proportions. These
    nets are generated as output from the code
    “nets_standarization_Bayes.Rmd”.

The name of each file comprises **two numbers** separated by a dash
“\_”, first number indicates the study and second number the population
(in case there is more than one). Study and population numbers
identifiers are provided in “selected_individual_networks.xlsx” and in
[list_nets_full_ref.pdf](list_nets_full_ref.pdf) .
