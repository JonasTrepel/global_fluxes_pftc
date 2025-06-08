library(data.table)
library(traitstrap)
library(tidyverse)
library(RColorBrewer)

source("R/functions/traitstrap_own.R")

trait <- fread("data/processed_data/preliminary_data/prelim_traits.csv") %>% 
  rename(PlotID = plot_id, 
         Site = site, 
         Tier = tier, 
         Taxon = species,
         Trait = trait_name, 
         Value = trait_value) %>% 
  dplyr::select(PlotID, Site, Tier, Trait, Value, Taxon) %>% 
  filter(!is.na(Value)) %>% 
  filter(!Taxon == "" ) %>% 
  filter(!(Trait == "sla_cm2_g" & Value > 1000))

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

community <- fread("data/processed_data/preliminary_data/prelim_cover.csv") %>% 
  rename(PlotID = plot_id, 
         Site = site, 
         Tier = tier, 
         Taxon = species,
         Cover = cover) %>% dplyr::select(PlotID, Site, Tier, Cover, Taxon) %>% filter(!Taxon == "" )

unique(community$PlotID)
unique(trait[grepl("US", trait$PlotID),]$PlotID)
unique(community[grepl("US", community$PlotID),]$Taxon)


# Hirarchy: 

#"plot_id" nested in "site" which is nested in "tier" 


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
  scale_hierarchy = c("Tier", "Site", "PlotID"),
  
  # min number of samples
  min_n_in_sample = 5
)

trait_missing(
  filled_trait = trait_filling,
  comm = community
)

# run nonparametric bootstrapping
np_bootstrapped_moments <- trait_np_bootstrap_own(
  trait_filling, 
  nrep = 100
)

np_bootstrapped_moments

# summarizes bootstrapping output
sum_boot_moment <- trait_summarise_boot_moments_own(
  np_bootstrapped_moments
)
sum_boot_moment
unique(sum_boot_moment[grepl("US", sum_boot_moment$PlotID),]$mean)

sum_boot_moment <- sum_boot_moment 

fwrite(sum_boot_moment, "data/processed_data/preliminary_data/prelim_traitstrap.csv")

unique(sum_boot_moment$Trait)
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
