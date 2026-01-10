# Plant traits decouple ecosystem carbon fluxes from temperature along tropical, temperate, and Arctic elevation gradients 

This repository contains all data and code used for the manuscript "Plant traits decouple ecosystem carbon fluxes from temperature along tropical, temperate, and Arctic elevation gradients" by Trepel, Niittynen, et al (under review)

The code is organized on four main folders: cleaning and preparation (R/data_cleaning_and_prep), getting the environemental data (R/environmental_data), performing the analysis (R/analysis) and create the visuals (R/viz). The data folder contains, well, the data and the the builds folder contains the plots. 

The R/environmental_data folder contains:
1. ....
2. ....

R/data_cleaning_and_prep:
1. get_south_africa_temp_and_par.R: Get temperature and PAR for South African data
2. make_vcg_3D_species_lists.R: Get the species lists for the Norway data
3. clean_data.R: script to clean the raw trait, flux and community data
4. traits_trap.R: calculate community weighted mean trait values
5. get_plant_fun_div.R: Calculate plant functional diversity (not actually used in the manuscript)
6. combine_all.R: Coombine analysis ready dataset

R/analysis: 
1. multiple_tiers_psem.R: Run piecewise structural equation models for the different model tiers (1-4, see paper for details)
2. get_model_predictions.R: get model predicted vs observed values for all six gradients
3. get_model_predictions_all_traits.R: get model predicted vs observed values for the four gradients with complete trait data
4. relative_variable_importance.R: Get relative variable importance
5. sensitivity/vpd_instead_growing_season_temp_models.R: Sensitivity analysis replacing growing season temperature with Vapor Pressure Deficit to test if moisture is more limiting
6. sensitivity/log_transform_fluxes.R: sensitivity analysis to test if temperature becomes significant when using the strict Arrhenius term: ln(flux) ~ -1/kT, with k = Boltzmann constant and T = temperature in Kelvin
7. sensitivity/loo_gradient_psem.R: repeating model tier 4, itertively leaving out one gradient

R/viz: 
1. map.R: create the map and data distribution for Fig. 1
2. fluxes_along_gradient.R: code for Fig. S1, showing the (lack of) response of fluxes along gradients of growing season temperature
3. all_traits_pca_along_gradient.R: Show how the trait principal component axes vary along gradients of growing season temperature 
4. correlations.R: visualize correlations
5. arrhenius_plot_supplement.R: plot ln(flux) ~ -1/kT, with k = Boltzmann constant and T = temperature in Kelvin for both growing season and instantaneous temperature 

*most manuscript figures have been combined in inkscape later on

We recommend running the scripts in this order. 

The data can be found in the datafolder and is decribed in the manuscript and in the associated datapapers: 
- V. Vandvik, et al., Plant traits and vegetation data from climate warming experiments along an 1100 m elevation gradient in Gongga Mountains, China. Sci Data 7, 189 (2020)
- V. Vandvik, et al., Plant traits and associated data from a warming experiment, a seabird colony, and along elevation in Svalbard. Sci Data 10, 578 (2023)
- A. H. Halbritter, et al., Plant trait and vegetation data along a 1314 m elevation gradient with fire history in Puna grasslands, Perú. Sci Data 11, 225 (2024)
- V. Vandvik, et al., Plant traits and associated ecological data from global change experiments and climate gradients in Norway. Sci Data 12, 1477 (2025)
- A. H. Halbritter, et al., Plant traits and associated ecological data from Afromontane grasslands of Maloti-Drakensberg, South Africa. Sci Data 12, 1778 (2025)

Please don't hesitate to contact us at jonas.trepel[at]bio.au.dk or pekka.niittynen[at]oulu.fi if you have any questions!  
