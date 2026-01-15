library(data.table)
library(traitstrap)
library(tidyverse)
library(RColorBrewer)

source("R/functions/traitstrap_own.R")

trait <- fread("data/processed_data/preliminary_data/prelim_traits.csv") %>% 
  rename(PlotID = plot_id, 
         Site = site, 
         Tier = tier, 
         Gradient = gradient, 
         Taxon = species,
         Trait = trait_name, 
         Value = trait_value) %>% 
  dplyr::select(PlotID, Site, Tier, Trait, Value, Taxon, Gradient) %>% 
  filter(!is.na(Value)) %>% 
  filter(!Taxon == "" ) %>% 
  filter(!(Trait == "sla_cm2_g" & Value > 1000))

unique(trait$Trait)


trait %>% 
  filter(Trait == "sla_cm2_g")
# Filter out trait outliers by comparing to plot level medians
trait <- trait %>% 
  arrange(Tier, PlotID, Trait, Taxon) %>% 
  # filter(Trait == "sla_cm2_g") %>% 
  group_by(PlotID, Trait, Taxon) %>% 
  mutate(n = n(),
         Tmd = median(Value, na.rm = TRUE),
         sd = sd(Value, na.rm = TRUE)) %>% 
  mutate(er = abs(Tmd - Value)/sd) %>% 
  filter(!(er > 2.25 & n >= 4 & !is.na(sd))) %>% 
  select(-c(n:er)) %>% ungroup()

unique(trait$Trait)


# Filter out trait outliers by comparing them to country level medians if less than 4 measurements in a plot
trait <- trait %>% 
  arrange(Tier, Trait, Taxon) %>% 
  # filter(Trait == "sla_cm2_g") %>%
  group_by(PlotID, Trait, Taxon) %>% 
  mutate(Tmd = median(Value, na.rm = TRUE),
         sd = sd(Value, na.rm = TRUE)) %>% 
  arrange(Tier, PlotID, Trait, Taxon) %>% 
  mutate(n = n()) %>% ungroup %>% 
  mutate(er = abs(Tmd - Value)/sd) %>% 
  filter(!(er > 2.25 & n < 4 & !is.na(sd))) %>% 
  select(-c(Tmd:er)) %>% ungroup()

table(trait[trait$Tier == "South_Africa_2023", ]$Trait)

community <- fread("data/processed_data/preliminary_data/prelim_cover.csv") %>% 
  mutate(gradient = case_when(
    grepl("South_Africa", tier) ~ "Drakensberg", 
    grepl("Peru", tier) ~ "Central Andes", 
    grepl("China", tier) ~ "Eastern Himalayas", 
    grepl("Colorado", tier) ~ "Rocky Mountains", 
    grepl("Norway", tier) ~ "Southern Scandes", 
    grepl("Svalbard", tier) ~ "Svalbard" 
  )) %>% 
  rename(PlotID = plot_id, 
         Site = site, 
         Tier = tier, 
         Gradient = gradient, 
         Taxon = species,
         Cover = cover) %>% dplyr::select(PlotID, Site, Tier, Cover, Taxon, Gradient) %>%
  filter(!Taxon == "" )

table(community$Gradient)
unique(community$PlotID)
unique(trait[grepl("South_Africa", trait$Tier),]$Taxon)
unique(community[grepl("South_Africa", community$Tier),]$Taxon)

unique(community[grepl("US", community$PlotID),]$Taxon)

trait_taxa <- unique(trait[grepl("South_Africa", trait$Tier), ]$Taxon)
community_taxa <- unique(community[grepl("South_Africa", community$Tier), ]$Taxon)
(common_taxa <- intersect(trait_taxa, community_taxa))


table(trait[grepl("South_Africa", trait$Tier), ]$Value)
# Hirarchy: 

#"plot_id" nested in "site" which is nested in "tier" which is nested in gradient (somewhat)


### start with trait filling

trait_filling <- trait_fill(
  # input data (mandatory)
  comm = community,
  traits = trait,
  
  # specifies columns in your data (mandatory)
  abundance_col = "Cover",
  taxon_col = "Taxon",
  trait_col = "Trait",
  value_col = "Value",
  
  # specifies sampling hierarchy
  scale_hierarchy = c("Gradient", "Tier", "Site", "PlotID"),
  
  # min number of samples
  min_n_in_sample = 5, 
  
  #don't calculate traits on global scale
  global = FALSE
)

trait_missing(
  filled_trait = trait_filling,
  comm = community
)

# run nonparametric bootstrapping
np_bootstrapped_moments <- trait_np_bootstrap(
  trait_filling, 
  nrep = 100
)

np_bootstrapped_moments
np_bootstrapped_moments[np_bootstrapped_moments$Tier == "South_Africa_2023" &
                          np_bootstrapped_moments$Trait == "sla_cm2_g", ]

# summarizes bootstrapping output
sum_boot_moment <- trait_summarise_boot_moments(
  np_bootstrapped_moments
)

sum_boot_moment[sum_boot_moment$Tier == "South_Africa_2023" & sum_boot_moment$Trait == "sla_cm2_g", ]
sum_boot_moment[sum_boot_moment$Tier == "Peru_2018" & sum_boot_moment$Trait == "p_percent", ]

sum_boot_moment
unique(sum_boot_moment[grepl("US", sum_boot_moment$PlotID),]$mean)

sum_boot_moment <- sum_boot_moment 

fwrite(sum_boot_moment, "data/processed_data/preliminary_data/prelim_traitstrap.csv")

unique(sum_boot_moment$Trait)
unique(sum_boot_moment$Tier)

sum_boot_moment %>% 
  filter(Trait == "sla_cm2_g") %>% 
  ggplot(aes(x = Tier, y = mean)) + 
  geom_boxplot() +
  geom_jitter()

sum_boot_moment %>% 
  filter(Trait == "leaf_area_cm2") %>% 
  mutate(mean = log(mean+1)) %>%
  ggplot(aes(x = Tier, y = mean)) + 
  geom_boxplot() +
  geom_jitter()

sum_boot_moment %>% 
  filter(Trait == "n_percent") %>% 
  #mutate(mean = log(mean+1)) %>%
  ggplot(aes(x = Tier, y = mean)) + 
  geom_boxplot() +
  geom_jitter()

