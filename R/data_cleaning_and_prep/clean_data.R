## GLOBAL FLUXES 
## April 23, 2024

## INITIAL DATA CLEANING AND COMBINATION 

## We have data from China, Colorado, Peru, Norway, South Africa, Svalbard 
## Variables of interest: 
 # 1. NEE
 # 2. Temperature
 # 3. Elevation
 # 4. Biomass
 # 5. Plot locations and metadata (i.e., plot ID etc)
 # 6. Datetime
 # 7. Traits: 
 #   - SLA
 #   - Leaf size
 #   - Leaf area
 #   - Leaf N (if available for all)
 #   - Leaf thickness
 #   - Wet mass
 #   - Dry mass
 #   - Plant height



## Suggested workflow: load data locationwise and clean step by step. 
## Standardize column names to be lowercase only
## plot_id should always contain Country_Site_Plot
## Combine in the end.  

library(data.table)
library(tidyverse)
library(tidylog)


# China ---------------------

ch_meta <- fread("data/raw_data/china/metach.csv") %>% 
  rename(site = Site, 
         elevation = Elevation, 
         latitude = Latitude, 
         longitude = Longitude) %>% 
  dplyr::select(-c(Gradient, Country))
#unfortunately only site metadata, not plotwise


## Flux data --------------------

#### 2016 ---------
# https://www.nature.com/articles/s41597-020-0529-0#Sec14
ch_flux_2016_raw <- fread("data/raw_data/china/PFTC_CO2flux_all_limits.txt") # Data from Inge
glimpse(ch_flux_2016_raw)

quantile(ch_flux_2016_raw$PAR, na.rm = T)

unique(ch_flux_2016_raw$type)

ch_flux_2016 <- ch_flux_2016_raw %>% 
  mutate(flux_best = NEE_lm,
         tier = "China_2016", 
         elevation = as.numeric(gsub("elev", "", site))) %>%
  dplyr::select(c(datetime, type, elevation, treatment, block, temp.C, flux_best, tier, PAR)) %>% 
  left_join(ch_meta) %>% 
  mutate(plot_id = paste0("CH_",site, block), 
         type = ifelse(type == "photo", "nee", "reco")) %>% 
  rename(temperature = temp.C, 
         plot = block, 
         par = PAR) %>% 
  filter(treatment %in% c("c", "otc"))
  
table(ch_flux_2016$block)
table(ch_flux_2016$plot_id)

table(ch_flux_2016$treatment)

# c = control, 
# otc = open top chamber 
# tt0 = local control 
# tt1 = warming
# tt2 = cooling
# tt3: extreme warming 
# tt4 = extreme cooling
# otc = open top chamber 

## for now I'd only move forward with the control and the open top chamber 

#### 



ggplot(data = ch_flux_2016) + 
  geom_boxplot(aes(x = type, y = flux_best)) +
  geom_jitter(aes(x = type, y = flux_best)) #This is weird. why all the negative fluxes in the photo part?
#also, I guess they followed the convention that negative fluxes represent outgonig CO2 while positive ones is uptake


## Traits ---------------------------------

#chemical 
ch_cht_raw <- fread("data/raw_data/china/PFTC1.2_China_2015_2016_ChemicalTraits.csv") 

ch_cht <- ch_cht_raw %>% 
  mutate(plot_id = paste0("CH_", destBlockID), 
         year = year(Date),
         tier = "China_2016", 
         leaf_number = NA) %>% 
  filter(year == 2016) %>%
  rename(site = Site,
         elevation = Elevation, 
         treatment = Treatment, 
         species = Taxon, 
         date = Date, 
         p_percent = P_percent, 
         c_percent = C_percent, 
         n_percent = N_percent, 
         cn_ratio = CN_ratio) %>% 
  dplyr::select(c(species, site, elevation, treatment, date, tier, p_percent, c_percent, n_percent, cn_ratio, plot_id, leaf_number)) %>% 
  pivot_longer(cols = c("p_percent", "c_percent", "n_percent", "cn_ratio"), names_to = "trait_name", values_to = "trait_value") %>% 
  left_join(ch_meta)


#leaves 
ch_lt_raw <- fread("data/raw_data/china/PFTC1.2_China_2015_2016_LeafTraits.csv")
names(ch_lt_raw)
#leaf number 
ch_lt <- ch_lt_raw %>% 
  mutate(plot_id = paste0("CH_", destBlockID), 
         year = year(Date),
         tier = "China_2016") %>% 
  filter(year == 2016) %>%
  rename(site = Site,
         elevation = Elevation, 
         treatment = Treatment, 
         species = Taxon, 
         wet_mass_g = Wet_Mass_g, 
         dry_mass_g = Dry_Mass_g,
         mean_leaf_thickness_mm = Leaf_Thickness_Ave_mm,
         leaf_area_cm2 = Leaf_Area_cm2,
         sla_cm2_g = SLA_cm2_g,
         ldmc = LDMC, 
         date = Date, 
         leaf_number = Leaf_number) %>% 
  dplyr::select(c(species, site, elevation, treatment, wet_mass_g, dry_mass_g, mean_leaf_thickness_mm, leaf_area_cm2, 
           sla_cm2_g, ldmc, date, plot_id, leaf_number, tier)) %>% 
  pivot_longer(cols = c("wet_mass_g", "dry_mass_g", "mean_leaf_thickness_mm", "leaf_area_cm2", 
           "sla_cm2_g", "ldmc"), names_to = "trait_name", values_to = "trait_value") %>% 
  left_join(ch_meta)


names(ch_cht)
names(ch_lt)

quantile(ch_lt_raw$SLA_cm2_g, na.rm = T )
#bind
ch_traits_raw <- rbind(ch_lt, ch_cht)

ch_traits <- ch_traits_raw %>% 
  mutate(treatment = tolower(treatment)) %>% 
  filter(treatment %in% c("c", "otc"))
unique(ch_traits$treatment)

## Cover and Biomass ----------------------------------------
ch_bio_raw <- fread("data/raw_data/china/China_2016_Biomass_cleaned.csv")
table(ch_bio_raw$site)
table(ch_bio_raw$plot)
table(ch_flux_2016$plot_id) #cool - we have plot level data here, but plot naming is not consistent between traits and fluxes
table(ch_traits$plot_id)

unique(ch_bio_raw$plot)
unique(ch_traits$plot_id)

ch_bio <- ch_bio_raw %>% 
  mutate(plot_gsub = gsub("I-", "", plot), 
         plot_id = paste0("CH_", site, plot_gsub), 
         tier = "China_2016") %>% 
  rename(species = speciesName) %>% 
  dplyr::select(-plot_gsub)


ch_tr <- ch_traits %>% 
  dplyr::select(plot_id, treatment) %>% unique()

ch_height_trait <- ch_bio %>% 
  dplyr::select(species, plot_id, site, tier,
                trait_value = height) %>% 
  mutate(trait_name = "plant_height_cm", 
         leaf_number = NA, 
         treatment = NA, date = NA) %>% 
  left_join(ch_meta)
  
  
## Summary -------------------------------

ch_flux_2016
setdiff(names(ch_traits), names(ch_height_trait))
ch_traits <- rbind(ch_traits, ch_height_trait) %>% unique()
ch_bio #no reasonable plot ID either... 

unique(ch_flux_2016$plot_id)
unique(ch_bio$plot_id)
unique(ch_traits$plot_id)




# Peru  ------------------------------
#https://www.nature.com/articles/s41597-024-02980-3

pe_meta_raw <- fread("data/raw_data/peru/PU.10_PFTC3.10_2020_Peru_Coordinates.csv")

names(pe_meta_raw) <- tolower(names(pe_meta_raw))

pe_meta <- pe_meta_raw %>% 
  rename(plot = plotid) %>% 
  filter(!is.na(plot)) %>%
  mutate(plot_id = paste0("PE_", site, "_", treatment ,"_",  plot)) %>% dplyr::select(-comment) %>% 
  dplyr::select(plot_id, elevation, longitude, latitude, burn_year, plot)

# PAR data -------------------

pe_par_2018_nd <- fread("data/raw_data/peru/peru_ecosystem_flux_output_2018.csv") %>%
  mutate( Treatment = ifelse(Site == "TRE" & Treatment == "B", "NB", Treatment),
          plot_id = paste0("PE_", Site, "_", Treatment ,"_",  Plot), 
         year = 2018, 
         type = case_when(Flux == "P" ~ "nee", Flux == "R" ~ "reco")) %>%
  dplyr::select(plot_id, par = PAR, type) %>%
  unique() %>% 
  mutate(tier = "Peru_2018")
  
pe_par_2019 <- fread("data/raw_data/peru/Flux_tent_environment_data_2019.csv") %>%
  mutate(Treatment = ifelse(Site == "TRE" & Treatment == "B", "NB", Treatment),
         plot_id = paste0("PE_", Site, "_", Treatment ,"_",  Plot), 
         mean_par = rowMeans(select(., contains("PAR")), na.rm = TRUE),
         date = dmy(Date), 
         type = case_when(Flux_type == "p" ~ "nee", Flux_type == "r" ~ "reco")) %>%
  filter(!is.na(type)) %>% 
  dplyr::select(plot_id, date, par = mean_par, type) %>%
  unique() %>% 
  mutate(tier = "Peru_2019")

pe_par_2020 <- fread("data/raw_data/peru/Tent_flux_env_data_all_plots.xlsx - tent_flux_env_data_control_2020.csv") %>%
  mutate(Treatment = ifelse(Site == "TRE" & Treatment == "B", "NB", Treatment),
         plot_id = paste0("PE_", Site, "_", Treatment ,"_",  Plot), 
         date = dmy(Date), 
         type = case_when(Flux == "P" ~ "nee", Flux == "R" ~ "reco")) %>%
  filter(!is.na(type)) %>% 
  dplyr::select(plot_id, date, par = PAR, type) %>%
  unique() %>% 
  mutate(tier = "Peru_2020")

## Flux data ----------------------------

pe_flux_raw <- fread("data/raw_data/peru/PFTC3_Puna_PFTC5_Peru_2018_2020_Cflux.csv")
names(pe_flux_raw)
unique(pe_flux_raw[,.(site, treatment)])


pe_flux_2018_dates <- pe_flux_raw %>% 
  filter(year == 2018) %>% 
  mutate(date = paste0(day, "-", month, "-", year),
         plot_id = paste0("PE_", site, "_", treatment ,"_",  plot_id)) %>%
  unique() %>% 
  dplyr::select(date, plot_id) %>% 
  unique()

unique(pe_flux_2018_dates$date)
pe_par_2018 <- pe_par_2018_nd %>% left_join(pe_flux_2018_dates) %>% 
  mutate(date = dmy(date))

pe_par <- rbind(pe_par_2018, pe_par_2019, pe_par_2020) %>% 
  unique() %>% 
  filter(!is.na(date))


n_distinct(pe_flux_2018_dates$plot_id) #thank god, only one date per plot! 


pe_flux <- pe_flux_raw %>% 
  rename(plot = plot_id, 
         type = flux, 
         temperature = t_ave) %>%  #some temperatures (if this column even is temperature) are oddly low: quantile(pe_flux_raw$t_ave, c(0.05, 0.95)) 
                                    # no way it was -60°C there. fix below
  mutate(
    tier = paste0("Peru_", year), 
    flux_best = linear_model,
    treatment = ifelse(site == "TRE" & treatment == "B", "NB", treatment), ###--> to be able to use the metadata as there seem to be a discrepancy between the treatment_ It's NB in metadata and B in the flux data 
    plot_id = paste0("PE_", site, "_", treatment ,"_",  plot),
    temperature = ifelse(temperature < -10, NA, temperature),
    date = paste0(day, "-", month, "-", year), 
    type = tolower(type), 
    type = ifelse(grepl("nee", type), "nee", type), 
    treatment = tolower(treatment), 
    par = NA, 
    date = as.Date(date, format = "%d-%B-%Y")) %>% 
  dplyr::select(c(type, temperature, plot_id, flux_best, tier, date, site, treatment)) %>% 
  left_join(pe_meta)  %>% 
  as.data.table() %>%
  unique() %>% 
  left_join(pe_par)

## Trait data ----------------------------

pe_t_raw <- fread("data/raw_data/peru/PFTC3-Puna-PFTC5_Peru_2018-2020_FunctionalTraits_clean.csv")

unique(pe_t_raw$trait)
quantile(pe_t_raw[pe_t_raw$trait == "sla_cm2_g" & course == "PFTC3", ]$value , na.rm = T )
quantile(pe_t_raw[pe_t_raw$trait == "c_percent" & course == "PFTC3", ]$value , na.rm = T )
unique(pe_t_raw[pe_t_raw$trait == "c_percent" & year == 2019, ]$value)
sum(is.na(unique(pe_t_raw[pe_t_raw$trait == "n_percent" & course == "PFTC5", ]$value)))


pe_traits <- pe_t_raw %>% 
  rename(trait_name = trait, 
         trait_value = value, 
         plot = plot_id, 
         leaf_number = leaf_nr, 
         species = taxon
         ) %>% 
  filter(trait_name %in% c("plant_height_cm", "wet_mass_g", "dry_mass_g", "mean_leaf_thickness_mm", "leaf_area_cm2", "sla_cm2_g",
                           "ldmc", "p_percent", "c_percent", "n_percent", "cn_ratio")) %>% 
  mutate( tier = paste0("Peru_", year), 
    plot_id = paste0("PE_", site, "_", treatment ,"_",  plot)) %>% 
  dplyr::select(c(species, site, elevation, treatment, plot_id, leaf_number, trait_name, trait_value, burn_year, year)) %>% 
  left_join(pe_meta)  %>% dplyr::select(-plot)

unique(pe_traits$trait_name)
unique(pe_traits[,.(site, treatment)])


## Biomass/Height -----------------------------

pe_height <- fread("data/raw_data/peru/PFTC3-Puna-Peru_2018-2019_CommunityStructure_clean.csv") %>% 
  filter(variable == "median_height") %>% 
  mutate(
    tier = paste0("Peru_", year)) %>% 
  dplyr::select(site, tier, treatment, value, plot_id) %>% 
  rename(height = value) %>% 
  mutate(plot_id = paste0("PE_", site, "_", treatment ,"_",  plot_id)) 


## Cover -----------------------------

pe_cover_raw <- fread("data/raw_data/peru/PFTC3-Puna-PFTC5_Peru_2018-2020_CommunityCover_clean.csv")

pe_cover <- pe_cover_raw %>% 
  mutate(tier = paste0("Peru_", year),
         plot_id = paste0("PE_", site, "_", treatment ,"_",  plot_id)) %>% 
  rename(species = taxon)


## Summary -----------------------------

pe_flux
pe_traits
#pe_bio #no plot ID... # but site - prob. better than nothing 
pe_cover
unique(pe_flux$plot_id)
unique(pe_traits$plot_id)
unique(pe_cover$plot_id)

# Svalbard  ------------------------------

sv_meta_raw <- fread("data/raw_data/svalbard/PFTC4_Svalbard_2018_metaItex.csv")
sv_coords_raw <- fread("data/raw_data/svalbard/PFTC4_Svalbard_Coordinates_ITEX.csv")
sv_coords.grad <- fread("data/raw_data/svalbard/PFTC4_Svalbard_Coordinates_Gradient.csv") %>% 
  mutate(plot_id = paste0("SV_", Gradient, Site, PlotID)) %>% 
  rename(elevation = Elevation_m,
         longitude = Longitude_E,
         latitude = Latitude_N) %>% 
  select(plot_id, elevation, latitude, longitude)

sv_meta <- sv_meta_raw %>% 
  left_join(sv_coords_raw) %>% 
  dplyr::select(-c(V8, V9, V10, V11, V12, V13, V14, Project, New_Site_name)) %>% 
  rename(Latitude = Latitude_N, 
         Longitude = Longitude_E, 
         Elevation = Elevation_m) %>% 
  mutate(plotid = PlotID, 
         plot_id = paste0("SV_ITEX", plotid), 
         plot_id = gsub("-", "_L", plot_id),
         plot = plotid, 
         Treatment = ifelse(Treatment == "CTL", "c", "otc")) %>% 
  dplyr::select(-c(plotid, plot, Treatment, PlotID, Site))

names(sv_meta) <- tolower(names(sv_meta))


## Flux data ----------------------------
#https://www.nature.com/articles/s41597-023-02467-7

#sv_flux_raw <- fread("data/raw_data/svalbard/Cflux_SV_ITEX_2018.csv") 
sv_flux_raw <- get(load("data/raw_data/svalbard/ITEX_all.Rdata"))
glimpse(sv_flux_raw)
names(sv_flux_raw)  <- tolower(names(sv_flux_raw))
#a first glimpse suggests that the weather was shit all the time and the PAR was unreasonable low. 
#we may have to exclude Svalbard
glimpse(sv_flux_raw)
sv_flux <- sv_flux_raw %>% 
  pivot_longer(cols = c(nee_ln, er_ln), names_to = "type", values_to = "flux_best") %>% 
  rename(plot = plotid) %>% 
  mutate(
    tier = paste0("Svalbard_", year), 
    plot_id = paste0("SV_", site, plot),
    plot_id = gsub("-", "_L", plot_id), 
    treatment = ifelse(treatment == "CTL", "c", "otc"), 
    par = par, 
    temperature = ifelse(!is.na(cantemp_light), (cantemp_light + cantemp_dark)/2, cantemp_dark), 
    type = ifelse(type == "nee_ln", "nee", "reco")) %>% 
  dplyr::select(c(type, plot_id, flux_best, tier, date, site, treatment, plot, par, temperature)) %>% 
  group_by(site, plot_id, plot, type) %>% 
  slice_max(par) %>% 
  ungroup() %>% 
  left_join(sv_meta, by = c("plot_id")) %>% as.data.table() 

setdiff(sv_meta$plot_id, sv_flux$plot_id)
setdiff(sv_flux$plot_id, sv_meta$plot_id)

# Svalbard gradients

sv_flux2 <- fread("data/raw_data/svalbard/Cflux_SV_Gradient_2018.csv") %>% 
  setNames(tolower(names(.))) %>% 
  rename(flux_best = nee,
         type = cover,
         par = par_mean,
         temperature = ir_mean) %>% 
  select(date, site, plotid, gradient, starttime, type, par, temperature, flux_best, rsqd, treatment) %>% 
  mutate(type = ifelse(type == "L", "nee", "reco")) %>% 
  filter(treatment == "Control") %>% 
  filter(!is.na(flux_best)) %>% 
  mutate(flux_best = flux_best*(-1)) %>% 
  # mutate(flux_best = ifelse(type == "reco", flux_best*(-1), flux_best)) %>% 
  group_by(date, site, plotid, gradient) %>% 
  mutate(temperature = zoo::na.approx(temperature)) %>% 
  ungroup() %>% 
  rename(plot = plotid) %>% 
  mutate(
    tier = paste0("Svalbard_", year(dmy(date))), 
    plot_id = paste0("SV_", gradient, site, plot),
    treatment = "c", 
    site = paste0(gradient, site)) %>% 
  dplyr::select(c(type, plot_id, flux_best, tier, date, site, treatment, plot, par, temperature)) %>% 
  group_by(site, plot_id, plot, type) %>% 
  slice_max(par) %>% 
  ungroup() %>% 
  left_join(., sv_coords.grad)

sv_flux <- bind_rows(sv_flux, sv_flux2 %>% mutate(date = dmy(date),
                                                  site = as.character(site)))

## Trait data ----------------------------

### ITEX 
sv_t_raw <- fread("data/raw_data/svalbard/PFTC4_Svalbard_2018_ITEX_Traits.csv")
names(sv_t_raw)  <- tolower(names(sv_t_raw))
table(sv_t_raw$trait)


sv_traits_itex <- sv_t_raw %>% 
  rename(trait_name = trait, 
       trait_value = value, 
       plot = plotid, 
       species = taxon, 
       elevation = elevation_m, 
       latitude = latitude_n, 
       longitude = longitude_e
) %>% 
  mutate( 
    plot = gsub("CH", "CAS", plot), 
    plot = gsub("SB", "BIS", plot), 
    plot = gsub("DH", "DRY", plot), 
    plot = gsub("-", "_L", plot), 
    site = "ITEX",
    treatment = ifelse(treatment == "CTL", "c", "otc"), 
    tier = paste0("Svalbard_", year), 
          plot_id = paste0("SV_", site, plot), 
          trait_name = tolower(trait_name)) %>% 
  filter(trait_name %in% c("plant_height_cm", "wet_mass_g", "dry_mass_g", "mean_leaf_thickness_mm", "leaf_area_cm2", "sla_cm2_g",
                           "ldmc", "p_percent", "c_percent", "n_percent", "cn_ratio")) %>% 
  dplyr::select(c(species, site, elevation, treatment, plot_id, trait_name, trait_value, longitude, latitude, tier, year))
sv_traits_itex

## Gradient 

sv_t_grad_raw <- fread("data/raw_data/svalbard/PFTC4_Svalbard_2018_Gradient_Traits.csv")

sv_gradient_t <- sv_t_grad_raw %>% 
  rename(trait_name = Trait, 
         trait_value = Value, 
         plot = PlotID, 
         elevation = Elevation_m, 
         latitude = Latitude_N, 
         longitude = Longitude_E, 
         species = Taxon, 
         site = Site, 
         gradient = Gradient, 
         year = Year
  ) %>% 
  mutate(
    treatment = "c", 
    tier = paste0("Svalbard_", year), 
    plot_id = paste0("SV_", gradient, site, plot),
    treatment = "c", 
    site = paste0(gradient, site), 
    species = str_to_sentence(species), 
    trait_name = tolower(trait_name)
  ) %>% 
  filter(trait_name %in% c("plant_height_cm", "wet_mass_g", "dry_mass_g", "mean_leaf_thickness_mm", "leaf_area_cm2", "sla_cm2_g",
                           "ldmc", "p_percent", "c_percent", "n_percent", "cn_ratio")) %>% 
  dplyr::select(c(species, site, elevation, treatment, plot_id, trait_name, trait_value, longitude, latitude, tier, year))

sv_traits <- rbind(sv_traits_itex, sv_gradient_t)

## Biomass/Height -----------------------------

# can't find any biomass for svalbard (which i guess makes sense)

sv_height_raw <- fread("data/raw_data/svalbard/PFTC4_Svalbard_2003_2015_ITEX_Community_Structure.csv") %>% 
  filter(Year == 2015)
names(sv_height_raw)  <- tolower(names(sv_height_raw))
sv_height_itex <- sv_height_raw %>% 
  mutate(plot = gsub("CH", "CAS", plotid), 
         plot = gsub("SB", "BIS", plot), 
         plot = gsub("DH", "DRY", plot), 
         plot = gsub("-", "_L", plot), 
         site = "ITEX",
         treatment = ifelse(treatment == "CTL", "c", "otc"), 
         tier = paste0("Svalbard_", year), 
         plot_id = paste0("SV_", site, plot)) %>% dplyr::select(plot_id, height)

sv_height_raw.gr <- fread("data/raw_data/svalbard/PFTC4_Svalbard_2018_Community_Structure_Gradient.csv") 
names(sv_height_raw.gr)  <- tolower(names(sv_height_raw.gr))

sv_gradient_height <- sv_height_raw.gr %>% 
  filter(variable == "MedianHeight_cm") %>% 
  rename(
         plot = plotid, 
         elevation = elevation_m, 
         latitude = latitude_n, 
         longitude = longitude_e, 
         height = value, 
  ) %>% 
  mutate(
    plot_id = paste0("SV_", gradient, site, plot),
  ) %>% 
  dplyr::select(height, plot_id)

sv_height <- rbind(sv_height_itex, sv_gradient_height)

## Cover -----------------------

sv_cover_raw <- fread("data/raw_data/svalbard/PFTC4_Svalbard_2003_2015_ITEX_Community.csv")
names(sv_cover_raw)  <- tolower(names(sv_cover_raw))
sv_cover_itex <- sv_cover_raw %>% 
  mutate(plot = gsub("CH", "CAS", plotid), 
         plot = gsub("SB", "BIS", plot), 
         plot = gsub("DH", "DRY", plot), 
         plot = gsub("-", "_L", plot), 
         site = "ITEX",
         treatment = ifelse(treatment == "CTL", "c", "otc"), 
         tier = paste0("Svalbard_", year), 
         plot_id = paste0("SV_", site, plot), 
         taxon = str_to_sentence(taxon)) %>% 
  filter(treatment == "c") %>% 
  rename(cover = abundance, 
         species = taxon) %>% 
  dplyr::select(site, plot_id, cover, species, tier)

sv_cover_g_raw <- fread("data/raw_data/svalbard/PFTC4_Svalbard_2018_Community_Gradient.csv")
names(sv_cover_g_raw)  <- tolower(names(sv_cover_g_raw))

sv_cover_gradient <- sv_cover_g_raw %>% 
  mutate(
    tier = paste0("Svalbard_", year), 
    plot_id = paste0("SV_", gradient, site, plotid),
    treatment = "c", 
    site = paste0(gradient, site), 
    species = str_to_sentence(taxon)) %>% 
  dplyr::select(site, plot_id, cover, species, tier)

sv_cover <- rbind(sv_cover_gradient, sv_cover_itex)

## Summary -----------------------------

sv_flux
sv_traits
sv_cover

unique(sv_flux$plot_id)
unique(sv_traits$plot_id)
unique(sv_cover$plot_id)


# Norway ------------------------------
#https://docs.google.com/document/d/1nXt4sljpExC_fSIVvGJ6-bDw9cAKANp6cwCWNklOHUU/edit
## meta data missing 

## Flux data ----------------------------

no_flux_raw <- fread("data/raw_data/norway/PFTC6_24h_cflux_allsites_2022.csv") 
quantile(no_flux_raw$PARavg, na.rm = T)

plot(no_flux_raw$flux, no_flux_raw$flux_corrected)
no_flux_raw2 <- no_flux_raw %>% 
  mutate(rowID = paste0("Row", 1:nrow(.)), 
         Minute = minute(datetime), 
         Hour = hour(datetime), 
         Day = day(datetime), 
         Month = month(datetime), 
         Year = year(datetime))  %>% 
  filter(flag %in% c("ok") & warming == "A") %>% 
  group_by(destSiteID, turfID, type) %>% 
  slice_max(PARavg, n = 3) %>% 
  group_by(destSiteID, turfID, type, warming) %>% 
  summarize(
    Minute = round(mean(Minute, na.rm = T), 0), 
    Hour = round(mean(Hour, na.rm = T), 0), 
    Day = round(mean(Day, na.rm = T),0), 
    Month = median(Month, na.rm = T), 
    Year = median(Year, na.rm = T), 
    temp_airavg = mean(temp_airavg, na.rm = T), 
    flux_corrected = mean(flux_corrected, na.rm = T),
    PARavg = mean(PARavg, na.rm = T)
  ) %>% 
  mutate(
    rawDateT = paste0(Year, "-", Month, "-", Day, " ", Hour, ":", Minute), 
    datetime = ymd_hm(rawDateT)
  ) %>%  dplyr::select(destSiteID, turfID, type, temp_airavg, flux_corrected, datetime, PARavg, warming)


glimpse(no_flux_raw)
table(no_flux_raw$flag)
hist(no_flux_raw$PARavg)

no_flux <- no_flux_raw2 %>% 
  rename(temperature = temp_airavg, 
         flux_best = flux_corrected) %>% 
  mutate(hour = hour(datetime), 
         flux_best = (flux_best*(-1))*0.2778, #x0.2778 is to convert from mmol m⁻² hr⁻¹ to µmol m⁻² s⁻¹
         tier = paste0("Norway_2022"), 
         plot_id = paste0("NO_", destSiteID, "_", turfID), 
         type = case_when(
           type == "ER" ~ "reco", 
           type == "NEE" ~ "nee", 
           type == "GPP" ~ "gpp"), 
         par = PARavg,  
         treatment = ifelse(warming == "A", "c", "transplanted"), 
         latitude = case_when(
           destSiteID == "Vik" ~ 60.8802, 
           destSiteID == "Lia" ~ 60.8599, 
           destSiteID == "Joa" ~ 60.8618, 
           destSiteID == "Hog" ~ 60.8760, 
         ), 
         longitude = case_when(
           destSiteID == "Vik" ~ 7.1699, 
           destSiteID == "Lia" ~ 7.1950, 
           destSiteID == "Joa" ~ 7.1680, 
           destSiteID == "Hog" ~ 7.17666, 
         ), 
         elevation = case_when(
           destSiteID == "Vik" ~ 469, 
           destSiteID == "Lia" ~ 1290, 
           destSiteID == "Joa" ~ 920, 
           destSiteID == "Hog" ~ 700, 
         ), 
         ) %>% 
  filter(type != "gpp" & treatment == "c") %>% #only peak time fluxes # alternatively, we could do a mean between 10 and 16 or so. 
  rename(site = destSiteID) %>% 
  dplyr::select(c(turfID, type, temperature, datetime, flux_best, site, tier, plot_id, par, treatment, elevation, longitude, latitude)) %>% 
  as.data.table()#no idea what's going on here w plot ID etc
         
unique(no_flux[,c("site", "turfID")])

ggplot(data = no_flux) + 
  geom_boxplot(aes(x = type, y = flux_best)) +
  geom_hline(yintercept = 0) +
  geom_jitter(aes(x = type, y = flux_best)) #CO2 uptake = positive, release = negative


ggplot(data = no_flux[grepl("nee", type), ]) + 
  geom_jitter(aes(x = temperature, y = flux_best)) +
  geom_smooth(aes(x = temperature, y = flux_best), method = "lm")

## Trait data ----------------------------


no_t_threed_raw <- fread("data/raw_data/norway/PFTC6_ThreeD_clean_leaf_traits_2022.csv") 
names(no_t_threed_raw)
unique(no_t_threed_raw$siteID)
unique(no_t_threed_raw$Namount_kg_ha_y)

unique(no_t_threed_raw$destSiteID)

unique(no_t_threed_raw[siteID %in% c("Vikesland") & warming == "A", ]$blockID)

unique(no_t_threed_raw[siteID %in% c("Vikesland") & warming == "W", .(blockID, turfID)])


unique(no_flux[site %in% c("Vik"), ]$turfID)

quantile(no_t_threed_raw[no_t_threed_raw$trait == "sla_cm2_g", ]$value , na.rm = T )


no_traits_raw <- no_t_threed_raw %>% 
  rename(trait_name = trait,
         trait_value = value) %>% 
  filter(trait_name %in% c("plant_height_cm", "wet_mass_g", "dry_mass_g", "mean_leaf_thickness_mm", "leaf_area_cm2", "sla_cm2_g",
                           "ldmc", "p_percent", "c_percent", "n_percent", "cn_ratio")) %>% 
  filter(gradient == "gradient") %>% 
  mutate(tier = paste0("Norway_2022"), 
         burn_year = NA, 
         aspect = NA, 
         longitude  = NA, 
         latitude = NA, 
         elevation = NA, 
         site = siteID, 
         date = as_date(date),
         year = year(date), 
         datetime = as_datetime(NA), 
         treatment = ifelse(warming == "A" & Namount_kg_ha_y == 0, "c", "transplanted or fertilized"), 
         site = case_when(
           siteID == "Hogsete" ~ "Hog",
           siteID == "Vikesland" ~ "Vik",
           siteID == "Liahovden" ~ "Lia",
           siteID == "Joasete" ~ "Joa",
         ), 
         plot_id = paste0("NO_", site, "_", turfID), 
         latitude = case_when(
           site == "Vik" ~ 60.8802, 
           site == "Lia" ~ 60.8599, 
           site == "Joa" ~ 60.8618, 
           site == "Hog" ~ 60.8760 
         ), 
         longitude = case_when(
           site == "Vik" ~ 7.1699, 
           site == "Lia" ~ 7.1950, 
           site == "Joa" ~ 7.1680, 
           site == "Hog" ~ 7.17666, 
         ), 
         elevation = case_when(
           site == "Vik" ~ 469, 
           site == "Lia" ~ 1290, 
           site == "Joa" ~ 920, 
           site == "Hog" ~ 700, 
         ), 
         
         plot_id = ifelse(is.na(turfID), NA, plot_id)

  ) %>% 
  filter(treatment == "c" & site %in% c("Vik", "Lia", "Joa", "Hog")) %>%
  dplyr::select(c(blockID, plot_id, site, elevation, aspect, trait_name, trait_value, species, tier, treatment,
           longitude, latitude, year, date, datetime))
  
unique(no_traits_raw$site)
unique(no_traits_raw$elevation)
unique(no_flux$plot_id)



#no_traits %>% dplyr::select(site, treatment, turfID) %>% unique()
#unique(no_traits$trait_name)



## Biomass/Height -----------------------------

## nothing available on OSF: https://osf.io/fcbw4/
## get height the funky way (from traits)
## will do below because I need the cover 

## Cover --------------------------------

### Code from Aud Halbritter
library("RSQLite")
con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "data/raw_data/norway/seedclim.sqlite")

dbListTables(con)


vcg_community_raw <- tbl(con, "turf_community")  |>
  dplyr::select(-cf, -flag) |>
  left_join(tbl(con, "turfs"), by = "turfID") |>
  ### adding to Auds code 
  mutate(tier = paste0("Norway_2022")) %>%
  # only control plots
  filter(TTtreat %in% c("TTC")) |>
  dplyr::select(-RTtreat, -GRtreat, -destinationPlotID) |>
  
  # join plot, block and site IDs
  left_join(tbl(con, "plots"), by = c("originPlotID" = "plotID")) |>
  rename("plotID" = originPlotID) |>
  dplyr::select(-aspect, -slope) |>
  left_join(tbl(con, "blocks"), by = c("blockID")) |>
  dplyr::select(-aspect, -slope) |>
  left_join(tbl(con, "sites"), by = c("siteID")) |>
  dplyr::select(-comment, -norwegian_name, -site_code, -c(biogeographic_zone:precipitation_level)) |>
  
  # filter 2 sites, and last year
  filter(siteID %in% c("Hogsete", "Vikesland"),
         year == 2019) |>
  mutate(site = case_when(
    siteID == "Hogsete" ~ "Hog", 
    siteID == "Vikesland" ~ "Vik"),
    turfID = paste("TTC", plotID),
    plot_id = paste0("NO_",site, "_", turfID) 
    ) %>%
  left_join(tbl(con, "taxon"), by = "species") |>
  group_by(site, plot_id, blockID, species_name) |>
  summarise(cover = mean(cover)) |>
  rename(species = species_name) |>
  collect() %>% 
  filter(plot_id %in% c("NO_Vik_TTC 141", "NO_Vik_TTC 146", "NO_Hog_TTC 110", "NO_Hog_TTC 115","NO_Hog_TTC 101")) %>% 
  as.data.table()

vcg_community_raw

turfBlockLeg <- vcg_community_raw %>% 
  dplyr::select(plot_id, blockID, site) %>% 
  mutate(blockID = gsub("Hog", "", blockID), 
         blockID = gsub("Vik", "", blockID), 
         blockID = as.numeric(blockID)) %>%
  unique() %>% 
  rename(plot_id_vcg = plot_id)

vcg_community <- vcg_community_raw %>% dplyr::select(-c(blockID))


# 3D
threeD_community <- read_csv("data/raw_data/norway/Three-D_clean_cover_2019-2022.csv") |>
  filter(grazing == "C",
         Nlevel %in% c(1, 2, 3),
         warming == "A", 
         year == 2022) |>
  mutate(destBlockID = as.character(destBlockID),
         elevation = if_else(destSiteID == "Lia", 1290, 920),
         siteID = destSiteID, 
         plot_id = paste0("NO_", destSiteID, "_", turfID)) |>
  group_by(siteID, species, plot_id) |>
  summarise(cover = mean(cover)) %>% 
  filter(plot_id %in% unique(no_traits_raw$plot_id)) %>% 
  dplyr::select(siteID, species, plot_id, cover) %>% 
  rename(site = siteID) 
threeD_community
vcg_community


no_cover <- rbind(threeD_community, vcg_community) %>% 
  mutate(plot_id = ifelse(site %in% c("Hog", "Vik"), NA, plot_id), 
         tier = paste0("Norway_2022"))

no_cover_ch <- rbind(threeD_community, vcg_community) %>% 
  mutate(tier = paste0("Norway_2022"))

summary(no_cover)
unique(no_cover$site)


## Summary -----------------------------

no_flux
no_traits <- no_traits_raw %>% 
  left_join(turfBlockLeg) %>% 
  mutate(plot_id = ifelse(is.na(plot_id), plot_id_vcg, plot_id)) %>% #
  dplyr::select(-plot_id_vcg)
no_traits     
no_cover

## Height --------
no_median.cover <- no_cover_ch %>% 
  group_by(plot_id) %>% 
  summarize(medianCov = median(cover, na.rm = T))

no_height <- no_traits %>% 
  filter(trait_name == "plant_height_cm") %>% 
  left_join(no_cover_ch[, c("plot_id", "cover", "species")]) %>% 
  left_join(no_median.cover) %>% 
  mutate(cover = ifelse(is.na(cover), medianCov, cover)) %>% 
  group_by(plot_id) %>% 
  summarize(
    meanHeight = sum(trait_value * cover, na.rm = TRUE) / sum(cover, na.rm = TRUE)
  )



# Colorado ------------------------------

us_meta_raw <- fread("data/raw_data/colorado/rmbl_site_info.csv")
us_meta <- us_meta_raw %>% 
  rename(
    latitude = lat, 
    longitude = long, 
    elevation = elev
  ) %>% 
  mutate(site = tolower(site)) %>% 
  dplyr::select(-lat_long)

## Flux data ----------------------------

us_flux_raw <- fread("data/raw_data/colorado/rmbl_gradient_flux_data_12042023.csv")

unique(us_flux_raw$time)
unique(us_flux_raw$plot)
unique(us_flux_raw$season)


unique(us_flux_raw$site)
us_flux <- us_flux_raw %>% 
  rename(flux_best = linear, 
         temperature = x7500_amb_temp) %>% #or better the measured one?
  mutate(
    site = tolower(site),
    plot = gsub("RESP", "", plot), 
    plot = gsub("A1II", "1", plot), 
    plot = gsub("A1I", "1", plot), 
    plot = gsub("A2I", "2", plot), 
    plot = gsub("A2II", "2", plot), 
    plot = gsub("B1II", "1", plot), 
    plot = gsub("B1I", "1", plot), 
    plot = gsub("B2II", "2", plot), 
    plot = gsub("B2I", "2", plot), 
    plot = gsub("C1II", "1", plot), 
    plot = gsub("C1I", "1", plot), 
    plot = gsub("C2II", "2", plot), 
    plot = gsub("C2I", "2", plot), 
    plot = gsub("2I", "2", plot), 
    plot = tolower(plot),
    temperature = ifelse(temperature < -10, NA, temperature),
    tier = paste0("Colorado_", year),
    plot_id = paste0("US_", site, plot),
    type = case_when(
      time == "DAY" ~ "nee", 
      time == "DAY RESP" ~ "reco",
      time == "NIGHT" ~ "night_reco"),
    plot_id = paste0("US_", site, plot)) %>% 
  dplyr::select(c(flux_best, temperature, site, date, season, type, plot, plot_id, tier)) %>% 
  filter(!is.na(type)) %>% 
  filter(site %in% c("almont", "cbt", "road", "pfeiler" ,"pbm" , "cinnamon") & type %in% c("nee", "reco")) %>% 
  left_join(us_meta) %>% filter(!is.na(flux_best)) %>% 
  filter(season == "Peak")
us_flux
us_flux_raw[time == "DAY RESP"]


ggplot(data = us_flux) + 
  geom_boxplot(aes(x = type, y = flux_best)) +
  geom_hline(yintercept = 0) +
  geom_jitter(aes(x = type, y = flux_best)) #CO2 uptake = positive, release = negative.

quantile(us_flux$temperature, na.rm = TRUE)
ggplot(data = us_flux[grepl("nee", type), ]) + 
  geom_jitter(aes(x = temperature, y = flux_best)) +
  geom_smooth(aes(x = temperature, y = flux_best), method = "lm")
  
unique(us_flux$site)
table(us_flux$tier)


## Trait data ----------------------------

us_t_raw <- fread("data/raw_data/colorado/rmbl_trait_data_master.csv")
glimpse(us_t_raw) 


#get plant height since this does not seem to be available for plots with numbers otherwise
us_pl_h <- us_t_raw %>%
  filter(!is.na(height)) %>% 
  dplyr::select(height, species) %>% 
  group_by(species) %>% 
  summarize(height = mean(height, na.rm = T))


us_traits <- us_t_raw %>% 
  filter(!is.na(plot)) %>% 
  dplyr::select(-height) %>% 
  left_join(us_pl_h) %>% 
  rename(
    elevation = elev, 
    latitude = lat, 
    longitude = long, 
    wet_mass_g = leaf_mass_fresh, 
    dry_mass_g = leaf_mass_dry,
    sla_cm2_g = SLA,
    p_percent = percent_P, 
    c_percent = percent_C,
    n_percent = percent_N,
    cn_ratio = CN_ratio,
    mean_leaf_thickness_mm = leaf_thickness,
    leaf_number = leaf_num, 
    plant_height_cm = height
  ) %>% 
  mutate(
    site = tolower(site),
    plot = tolower(plot), 
    plot_id = paste("US_", site, plot), 
    plot_id = gsub(" ", "", plot_id),
    tier = paste0("Colorado_", year),
    ldmc = LDMC/1000, ## assuming colorado gives LDMC in mg instead of g 
    leaf_area_cm2 = leaf_area/100 ## assuming colorado gives LDMC in mm2 instead of cm2
  ) %>% 
  filter(site %in% c("almont", "cbt", "road", "pfeiler", "pbm", "cinnamon")) %>% #same sites as the flux data
  pivot_longer(cols = c("plant_height_cm", "wet_mass_g", "dry_mass_g", "leaf_area_cm2", "sla_cm2_g", 
                       "ldmc", "p_percent", "c_percent", "n_percent", "cn_ratio", "mean_leaf_thickness_mm"), 
               names_to = "trait_name", values_to = "trait_value") %>% 
  dplyr::select(c(plot_id, site, plot, trait_name, trait_value, species, year, tier)) %>% 
  left_join(us_meta) %>% #vast majority (like 60,000) is NA for the plot_ 
  #filter(!grepl("NA", plot_id)) %>% 
  mutate(trait_value = ifelse(year == 2016 & trait_name == "leaf_area_cm2", trait_value*10, trait_value))

us_traits[us_traits$site == "pbm" & us_traits$trait_name == "p_percent", ]$trait_value

us_traits
unique(us_traits[us_traits$site == "cinnamon", ]$plot)

## Biomass/Height -----------------------------
#cover
us_coverH_raw <- fread("data/raw_data/colorado/veg_cover_data_rmbl.csv") %>% 
  mutate(site = tolower(site),
         plot = tolower(plot), 
         plot_id = paste("US_", site, plot), 
         plot_id = gsub(" ", "", plot_id),
         tier = paste0("Colorado_", year), 
         vegCover = (herb + shrub + grass + cactus + forb)) %>% 
  filter(tier %in% c("Colorado_2016", "Colorado_2018")) %>% 
  filter(plot_id %in% unique(us_flux$plot_id)) %>% 
  dplyr::select(plot_id, tier, vegCover, site) %>% 
  group_by(tier, site, plot_id) %>% 
  summarize(vegCover = mean(vegCover, na.rm = T))


## Height -------------
us_height_raw <- fread("data/raw_data/colorado/veg_height_data_rmbl.csv") %>% 
  mutate(site = tolower(site),
         plot = tolower(plot), 
         plot_id = paste("US_", site, plot), 
         plot_id = gsub(" ", "", plot_id),
         tier = paste0("Colorado_", year)) %>% 
  filter(tier %in% c("Colorado_2016", "Colorado_2018")) %>% 
  filter(plot_id %in% unique(us_flux$plot_id)) %>% 
  rename(vegHeight = mean) %>% 
  dplyr::select(plot_id, tier, vegHeight, site) %>% unique() %>% 
  group_by(tier, site, plot_id) %>% 
  summarize(vegHeight = mean(vegHeight, na.rm = T))

us_height <- us_height_raw %>% 
  left_join(us_coverH_raw) %>% unique()


# don't have anything yet

## Cover --------------------------

us_cover <- fread("data/raw_data/colorado/rmbl_plot_data_2016.csv") %>% 
  mutate(site = tolower(site),
         plot = tolower(plot), 
         plot_id = paste("US_", site, plot), 
         plot_id = gsub(" ", "", plot_id),
         tier = paste0("Colorado_", year), 
         taxon_std = case_when(
           .default = taxon_std, 
           taxon_std == "Cirsium sp." & plot_id %in% c("US_road4", "US_road5") ~ "Cirsium undulatum", 
           taxon_std == "Antennaria parvifolia" & plot_id %in% c("US_almont2") ~ "Antennaria media", 
           taxon_std == "Poa pratensis" & plot_id %in% c("US_cbt2", "US_cbt3", "US_cbt5") ~ "Poa reflexa", 
           taxon_std == "Erigeron coulteri" & plot_id %in% c("US_road1", "US_road3", "US_pfeiler4") ~ "Erigeron sp.", 
           taxon_std == "Erigeron coulteri" & plot_id %in% c("US_road1", "US_road3", "US_pfeiler4") ~ "Erigeron sp.", 
           taxon_std == "Senecio crassulus" & plot_id %in% c("US_pfeiler3") ~ "Senecio sp.", 
           taxon_std == "Hydrophyllum fendleri" & plot_id %in% c("US_pfeiler3") ~ "Hydrophyllum capitatum")) %>% 
  filter(abundance > 0 & abundance < 100) %>% 
  rename(species = taxon_std, 
         cover = abundance) %>% 
  dplyr::select(site, plot_id, cover, species, tier) %>% 
  rename(site_id = site) 
  

#unique(us_cover_raw$treatment)
## Summary -----------------------------

us_flux
us_traits

unique(us_traits$plot_id)

unique(us_flux$plot_id)

# South Africa ------------------------------


##metadata 
library(sf) 
sa_wp_raw <- read_sf("data/raw_data/south_africa/PFCT7_site_waypoints.gpx") %>%
  dplyr::select(name) %>%
  filter(name != "PFTC5 E4")

sa_wp_coords <- st_coordinates(sa_wp_raw)

sa_wp <- sa_wp_raw %>% 
  cbind(sa_wp_coords) %>% 
  rename(latitude = Y, 
         longitude = X) %>%
  mutate(
    aspect = case_when(
      grepl("E", name) ~ "east", 
      grepl("W", name) ~ "west", 
    ),
    elevation = case_when(
      grepl("C1", name) ~ 2000, 
      grepl("C2", name) ~ 2200, 
      grepl("C3", name) ~ 2400, 
      grepl("C4", name) ~ 2600, 
      grepl("C5", name) ~ 2800, 
      grepl("C6", name) ~ 3000
    ), 
    site = case_when(
      grepl("C1", name) ~ 1, 
      grepl("C2", name) ~ 2, 
      grepl("C3", name) ~ 3, 
      grepl("C4", name) ~ 4, 
      grepl("C5", name) ~ 5, 
      grepl("C6", name) ~ 6
    ),
    plot = case_when(
      grepl("E1", name)| grepl("W1", name) ~ 1, 
      grepl("E2", name)| grepl("W2", name) ~ 2, 
      grepl("E3", name)| grepl("W3", name) ~ 3, 
      grepl("E4", name)| grepl("W4", name) ~ 4, 
      grepl("E5", name)| grepl("W5", name) ~ 5
    ), 
    name = NULL, 
    geometry = NULL, 
    unique_site = paste0(site, "_", aspect)
  ) %>% as.data.table()

# Function to infer plots 2-4
calc_intermediate_plots <- function(site_data) {
  lat1 <- site_data$latitude[site_data$plot == 1]
  lon1 <- site_data$longitude[site_data$plot == 1]
  lat5 <- site_data$latitude[site_data$plot == 5]
  lon5 <- site_data$longitude[site_data$plot == 5]
  
  data.frame(
    unique_site = unique(site_data$unique_site),
    plot = 2:4,
    latitude = c(
      (3 * lat1 + lat5) / 4,
      (lat1 + lat5) / 2,
      (lat1 + 3 * lat5) / 4
    ),
    longitude = c(
      (3 * lon1 + lon5) / 4,
      (lon1 + lon5) / 2,
      (lon1 + 3 * lon5) / 4
    ), 
    aspect = unique(site_data$aspect), 
    elevation = unique(site_data$elevation), 
    site = unique(site_data$site)
  )
}

# Calculate intermediate plots for each site
intermediate_plots <- sa_wp %>%
  group_by(unique_site) %>%
  do(calc_intermediate_plots(.))

# Combine original and intermediate plots
sa_wp2 <- bind_rows(sa_wp, intermediate_plots) %>%
  arrange(unique_site, plot) %>% filter(!is.na(aspect)) %>% mutate(site = paste0(elevation, aspect))

sa_wp_sf <- st_as_sf(sa_wp2, 
                     coords = c("longitude", "latitude"), 
                     crs = 4326)
mapview::mapview(sa_wp_sf)

fwrite(sa_wp2, "data/raw_data/south_africa/PFCT7_plot_locations_clean.csv")


### get dates 

sa_meta_raw <- fread("data/raw_data/south_africa/x_PFTC7_clean_ecosystem_fluxes_2023.csv")

sa_meta <- sa_meta_raw %>% 
  filter(day_night == "day")  %>% 
  filter(flag %in% c("okay", "manual_flux_time_selection"))%>% 
  mutate(unique_site = paste0(site_id, "_", aspect)) %>%
  dplyr::select(unique_site, date ) %>% 
  unique() %>% 
  mutate(discard = case_when(
    .default = "keep", 
    unique_site == "3_west" & date == "2023-12-09" ~ "discard", 
    date == "2023-12-03" ~ "discard"
  )) %>% 
  filter(discard == "keep") %>% 
  dplyr::select(-discard)

sa_meta

sa_wp3 <- sa_wp2 %>% left_join(sa_meta)


## Flux data ----------------------------
#temperature first_.. 

sa_temp <- fread("data/raw_data/south_africa/pftc7_flux_par_and_tmp.csv") %>% 
  dplyr::select(file, date, temperature, par) %>% 
  unique() %>% 
  mutate(par = case_when(
    .default = par, #update NAs using field notebook
    file == "1_2000_east_1_day_photo.txt" & date == as.Date("2023-12-03") ~ 838.1,
    file == "1_2000_east_1_day_resp.txt" & date == as.Date("2023-12-03") ~ 5.1,
    file == "1_2000_east_2_day_photo.txt" & date == as.Date("2023-12-03") ~ 944.3,
    file == "1_2000_east_2_day_redo_resp.txt" & date == as.Date("2023-12-03") ~ 5.9,
    file == "1_2000_east_2_day_resp.txt" & date == as.Date("2023-12-03") ~ 5.9,
    file == "1_2000_east_3_day_photo.txt" & date == as.Date("2023-12-03") ~ 884.2,
    file == "1_2000_east_3_day_resp.txt" & date == as.Date("2023-12-03") ~ 5.9,
    file == "1_2000_east_4_day_photo.txt" & date == as.Date("2023-12-03") ~ 719,
    file == "1_2000_east_4_day_resp.txt" & date == as.Date("2023-12-03") ~ 5.9,
    file == "1_2000_east_5_day_photo.txt" & date == as.Date("2023-12-03") ~ 609.2,
    file == "1_2000_west_1_day_photo.txt" & date == as.Date("2023-12-03") ~ 979.5,
    file == "1_2000_west_1_day_redo2_resp.txt" & date == as.Date("2023-12-03") ~ 5.7,
    file == "1_2000_west_1_day_redo_resp.txt" & date == as.Date("2023-12-03") ~ 5.6,
    file == "1_2000_west_2_day_photo.txt" & date == as.Date("2023-12-03") ~ 1123,
    file == "1_2000_west_2_day_resp.txt" & date == as.Date("2023-12-03") ~ 6.8,
    file == "1_2000_west_3_day_photo.txt" & date == as.Date("2023-12-03") ~ 970,
    file == "1_2000_west_3_day_resp.txt" & date == as.Date("2023-12-03") ~ 6.3,
    file == "1_2000_west_4_day_photo.txt" & date == as.Date("2023-12-03") ~ 1007,
    file == "1_2000_west_4_day_resp.txt" & date == as.Date("2023-12-03") ~ 6.8,
    file == "2_2200_east_1_day_photo.txt" & date == as.Date("2023-12-05") ~ 1390,
    file == "2_2200_east_1_day_resp.txt" & date == as.Date("2023-12-05") ~ 2.7,
    file == "2_2200_east_2_day_photo.txt" & date == as.Date("2023-12-05") ~ 1261,
    file == "2_2200_east_2_day_resp.txt" & date == as.Date("2023-12-05") ~ 3.3,
    file == "2_2200_east_3_day_photo.txt" & date == as.Date("2023-12-05") ~ 1200,
    file == "2_2200_east_3_day_resp.txt" & date == as.Date("2023-12-05") ~ 2.9,
    file == "2_2200_east_4_day_photo.txt" & date == as.Date("2023-12-05") ~ 1128,
    file == "2_2200_east_4_day_resp.txt" & date == as.Date("2023-12-05") ~ 3.5,
    file == "2_2200_east_5_day_photo.txt" & date == as.Date("2023-12-05") ~ 1152,
    file == "2_2200_east_5_day_redo2_photo.txt" & date == as.Date("2023-12-05") ~ 449,
    file == "2_2200_east_5_day_redo_photo.txt" & date == as.Date("2023-12-05") ~ 670,
    file == "2_2200_east_5_day_resp.txt" & date == as.Date("2023-12-05") ~ 3.2,
    file == "2_2200_west_1_day_photo.txt" & date == as.Date("2023-12-05") ~ 1023,
    file == "2_2200_west_1_day_redo_photo.txt" & date == as.Date("2023-12-05") ~ 1068,
    file == "2_2200_west_1_day_resp.txt" & date == as.Date("2023-12-05") ~ 2.5,
    file == "2_2200_west_2_day_photo.txt" & date == as.Date("2023-12-05") ~ 1172,
    file == "2_2200_west_2_day_resp.txt" & date == as.Date("2023-12-05") ~ 2.3,
    file == "2_2200_west_3_day_photo.txt" & date == as.Date("2023-12-05") ~ 1174,
    file == "2_2200_west_3_day_resp.txt" & date == as.Date("2023-12-05") ~ 6.4,
    file == "2_2200_west_4_day_photo.txt" & date == as.Date("2023-12-05") ~ 677,
    file == "2_2200_west_4_day_resp.txt" & date == as.Date("2023-12-05") ~ 6.4,
    file == "2_2200_west_5_day_photo.txt" & date == as.Date("2023-12-05") ~ 1071,
    file == "2_2200_west_5_day_resp.txt" & date == as.Date("2023-12-05") ~ 7))

#https://docs.google.com/document/d/1P2X-3IIQE6IQwvgvDQe9YJLqVWd2srqzxJZWWoM3mig/edit
sa_flux_raw <- fread("data/raw_data/south_africa/x_PFTC7_clean_ecosystem_fluxes_2023_with_file.csv")
names(sa_flux_raw)
sa_flux <- sa_flux_raw %>% 
  unique() %>% 
  left_join(sa_temp) %>% 
  filter(flag %in% c("okay", "manual_flux_time_selection")) %>% 
  rename(flux_best = flux_value, 
         elevation = elevation_m_asl, 
         plot = plot_id) %>%
  filter(flux_type %in% c("nee", "resp_day", "resp_night")) %>% 
  mutate(tier = "South_Africa_2023", 
         plot_id = paste0("SA_", elevation, aspect, plot), 
         site = paste0(elevation, aspect),
         type = case_when(
           flux_type == "nee" ~ "nee", 
           flux_type == "resp_day" ~ "reco", 
           flux_type == "resp_night" ~ "night_reco", 
         )) %>% 
  dplyr::select(c(plot_id, tier, flux_best, type, site, plot, aspect, temperature, par, elevation)) %>% 
  left_join(sa_wp3) %>% 
  filter(type != "night_reco")

range(sa_flux[type == "nee", flux_best], na.rm = T)
range(sa_flux[type == "reco", flux_best], na.rm = T)


sa_flux %>% 
  filter(is.na(par)) %>% 
  select(plot_id, type) %>% 
  unique()

## Trait data ----------------------------

sa_t_raw <- fread("data/raw_data/south_africa/iv_PFTC7_clean_elevationgradient_traits_2023.csv")
table(sa_t_raw$problem_flag)
table(sa_t_raw$traits)

sa_traits <- sa_t_raw %>% 
  rename(elevation = elevation_m_asl, 
         site = site_id, 
         plot = plot_id, 
         trait_name = traits, 
         trait_value = value) %>% 
  mutate(tier = "South_Africa_2023", 
         plot_id = paste0("SA_", elevation, aspect, plot), 
         year = 2023, 
         trait_name = ifelse(trait_name == "veg_height_cm", "plant_height_cm", trait_name), 
         species = str_to_sentence(species), 
         site = paste0(elevation, aspect), 
         trait_name = case_when(
           .default = trait_name,
           trait_name == "sla" ~ "sla_cm2_g",
           trait_name == "dry_mass" ~ "dry_mass_g",
           trait_name == "veg_height" ~ "plant_height_cm",
           trait_name == "wet_mass" ~ "wet_mass_g",
           trait_name == "leaf_area" ~ "leaf_area_cm2",
           trait_name == "leaf_thickness" ~ "leaf_thickness_mm"
           )) %>% 
  filter(trait_name %in% c("plant_height_cm", "wet_mass_g", "dry_mass_g", "leaf_area_cm2", "leaf_thickness_mm", "sla_cm2_g",
                        "ldmc")) %>% 
  dplyr::select(c(plot_id, site, plot, elevation, aspect, trait_name, trait_value, species, tier, year)) %>% 
  left_join(sa_wp2)
    
unique(sa_traits$species)
## Biomass -----------------------------

## I don't think we have clean biomass data yet 

## Cover --------------------------------
sa_cover_raw <- fread("data/raw_data/south_africa/i_PFTC7_clean_elevationgradient_community_2023.csv")
sa_cover <- sa_cover_raw %>% 
  rename(site = site_id, 
         plot = plot_id) %>% 
  mutate(tier = "South_Africa_2023", 
         plot = round(plot, 0), 
         plot_id = paste0("SA_", elevation_m_asl, aspect, plot), 
         year = 2023, 
         species = gsub("_", " ", species), 
         species = str_to_sentence(species)) %>% 
  filter(!site == 6)

## Height ------------

sa_height_raw <- fread("data/raw_data/south_africa/ii_PFTC7_clean_elevationgradient_community_structure_2023.csv")

sa_height <- sa_height_raw %>% 
  mutate(tier = "South_Africa_2023",
         plot = round(plot_id, 0),
         plot_id = paste0("SA_", elevation_m_asl, aspect, plot)) %>% 
  filter(variable == "vegetation_height" & site_id != 6) %>% 
    rename(meanHeight = value, 
           elevation = elevation_m_asl, 
           site = site_id)  %>% 
    dplyr::select(tier, site, plot_id, meanHeight, elevation, aspect) %>% 
  mutate(site = paste0(elevation, aspect)) %>% dplyr::select(-elevation)


## Summary -----------------------------
sa_flux 
sa_traits

# combine everything ---------------------------------------------------

## fluxes ------------

names(ch_flux_2016)# tier, "datetime", "type", "elevation", "treatment", "plot", "temperature", "flux_best", "site", "latitude", "longitude", "plot_id" 
ch_flux_2016.fin <- ch_flux_2016 %>% 
  mutate(burn_year = NA, 
         datetime = dmy_hm(datetime),
         date = date(datetime),
         date = as_date(date, format = '%Y.%m.%d'),
         
         year = year(date), 
         
         aspect = NA, 
         season = "not_part_of_the_dataset")

names(pe_flux) #extra = burn_year
pe_flux_fin <- pe_flux %>% 
  rename(date_old = date) %>% 
  mutate(
       date = gsub("March", "03", date_old),
       date = gsub("July", "07", date),
       date = gsub("April", "04", date),
       date = gsub("November", "11", date),
       datetime= as_datetime(NA), 
     #  date = as_date(date, format = '%d-%m-%Y'),
       date = as_date(date, format = '%Y-%m-%d'),
       year = year(date), 
       aspect = NA, 
       par = NA,
       season = "not_part_of_the_dataset") %>% 
  dplyr::select(-date_old)


names(sv_flux) 
sv_flux_fin <- sv_flux %>% 
  mutate(burn_year = NA, 
         datetime = as_datetime(NA),
         aspect = NA, 
         date = as_date(date, format = '%d.%m.%Y'),
         date = as_date(date, format = '%Y.%m.%d'),
         year = year(date), 
         season = "not_part_of_the_dataset")

names(no_flux) #
no_flux_fin <- no_flux %>%
  mutate(burn_year = NA,
         date = date(datetime),
         date = as_date(date, format = '%Y.%m.%d'),
         year = year(date), 
         aspect = NA,
         season = "not_part_of_the_dataset") %>% rename(plot = turfID)

names(us_flux)
us_flux_fin <- us_flux %>% 
  mutate(burn_year = NA, 
         datetime = as_datetime(NA), 
         date = as_date(date), 
         date = as_date(date, format = '%Y.%m.%d'),
         year = year(date), 
         aspect = NA, 
         treatment = NA,
         par = NA)

setdiff(names(no_flux_fin), names(pe_flux_fin))
setdiff(names(pe_flux_fin), names(no_flux_fin))


names(sa_flux)
sa_flux_fin <- sa_flux %>% 
  mutate(burn_year = NA, 
         datetime = as_datetime(NA), 
         date = as_date(date, format = '%Y.%m.%d'),
         treatment = NA,
         year = year(date), 
         season = "not_part_of_the_dataset") %>% 
  filter(!type == "night_reco") %>% 
  dplyr::select(-unique_site)

fluxes.combined_raw <- rbind(ch_flux_2016.fin, 
                         pe_flux_fin, 
                         sv_flux_fin, 
                         no_flux_fin,
                         us_flux_fin, 
                         sa_flux_fin)

setdiff(names(ch_flux_2016.fin), names(sa_flux_fin))
setdiff(names(sa_flux_fin), names(ch_flux_2016.fin))


unique(fluxes.combined_raw[grepl("Peru", tier), ]$treatment)

fluxes.combined <- fluxes.combined_raw %>% 
  mutate(treatment = ifelse(is.na(treatment), "c", treatment),
         location = case_when(
    grepl("China", tier) ~ "China", 
    grepl("Colorado", tier) ~ "USA", 
    grepl("Norway", tier) ~ "Norway", 
    grepl("Peru", tier) ~ "Peru", 
    grepl("South_Africa", tier) ~ "South Africa", 
    grepl("Svalbard", tier) ~ "Svalbard"
  ), 
  location = as.factor(location)) %>% 
  filter(!tier %in% c("Colorado_2009", "Colorado_2010")) %>% 
  filter(treatment %in% c("c", "nb", "b"))


unique(fluxes.combined[is.na(fluxes.combined$latitude), ])
table(fluxes.combined$type)
unique(fluxes.combined[tier == "South_Africa_2023", ]$par)


fwrite(fluxes.combined, "data/processed_data/preliminary_data/prelim_fluxes.csv")

nrow(fluxes.combined[type == "nee" & !grepl("Colo", tier), ])
nrow(fluxes.combined[type == "reco" & !grepl("Colo", tier), ])



nrow(fluxes.combined[type == "nee" & treatment != "otc", ])
ggplot(data = fluxes.combined[type == "nee" & treatment != "otc", ]) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_jitter(aes(x = temperature, y = flux_best, color = tier), size = 2, alpha = 0.75) +
  scale_color_viridis_d() +
  labs(y = "NEE", y = "Temperature") +
  theme_bw() +
  geom_smooth(aes(x = temperature, y = flux_best), method = "lm")

ggplot(data = fluxes.combined[treatment != "otc", ]) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_jitter(aes(x = type, y = flux_best, color = location), size = 1.5, alpha = 0.75) +
  geom_boxplot(aes(x = type, y = flux_best), size = 1, alpha = 0.75) +
  scale_color_viridis_d() +
  labs(y = "Flux (µmol C02/m2/s)", x = "Flux type", color = "Location") +
  theme_bw()

fwrite(fluxes.combined, "data/processed_data/preliminary_data/prelim_fluxes.csv")



unique(fluxes.combined[is.na(par),]$tier)

m0 <- lm(flux_best ~ 1, data = fluxes.combined[type == "nee" & treatment != "otc", ])
summary(m0)

m1 <- lm(flux_best ~ temperature, data = fluxes.combined[type == "nee" & treatment != "otc", ])
summary(m1)



fluxes.combined2 <- fluxes.combined %>%
  dplyr::select(plot_id, tier, date, longitude, latitude) %>% 
  unique()

flux_sf <- st_as_sf(fluxes.combined2, 
                    coords = c("longitude", "latitude"), 
                    crs = 4326)

mapview::mapview(flux_sf)

flux_sf2 <- flux_sf %>% 
  dplyr::select(plot_id, tier, date) %>% 
  mutate(date = paste0(date), 
         flux_id = paste0(plot_id, 1:nrow(.))) %>% st_make_valid()
str(flux_sf2)
mapview::mapview(flux_sf2)
n_distinct(flux_sf2$plot_id)
n_distinct(flux_sf2$flux_id)

st_write(flux_sf2, "data/processed_data/preliminary_data/prelim_flux_loc.gpkg", append = FALSE)
st_write(flux_sf2, "data/processed_data/preliminary_data/prelim_flux_loc.shp", append = FALSE)

mapview::mapview(flux_sf2 %>% filter(plot_id == "PE_ACJ_C_1"))


## traits ---------------------------------

names(ch_traits) #get rid of date and replace w year because who cares

ch_traits.fin <- ch_traits %>% 
  mutate(burn_year = NA, 
         aspect = NA, 
         date = as_date(date),
         year = year(date)) %>% dplyr::select(-date, -leaf_number)

names(pe_traits)
pe_traits.fin <- pe_traits %>%
  mutate(aspect = NA,
         tier = paste0("Peru_", year)) %>% 
  dplyr::select(-c(leaf_number))

unique(pe_traits.fin$tier)
pe_traits.fin[tier == "Peru_2020" & trait_name == "p_percent" & treatment == "C", ]

names(sv_traits) #
sv_traits.fin <- sv_traits %>%
  mutate(aspect = NA, 
         burn_year = NA)
unique(sv_traits$tier)
# names(no_traits) #not yet finished 
no_traits.fin <- no_traits %>%
  mutate(aspect = NA, 
         burn_year = NA) %>% dplyr::select(-datetime, -date, -blockID)

unique(no_traits.fin$trait_name)

names(us_traits) # get rd of 
us_traits.fin <- us_traits %>%
     mutate(aspect = NA, 
            burn_year = NA, 
            treatment = NA) %>% dplyr::select(-plot)
unique(us_traits.fin$site)
us_traits.fin[us_traits.fin$trait_name == "p_percent" & 
                us_traits.fin$site == "road", ]$trait_value


names(sa_traits) # get rd of 

sa_traits.fin <- sa_traits %>% 
  mutate(burn_year = NA, 
         treatment = NA) %>% dplyr::select(-plot, - unique_site)

names(ch_traits.fin)
      names(pe_traits.fin) 
      names(sv_traits.fin)  
      names(no_traits.fin) 
      names(us_traits.fin) 
      names(sa_traits.fin) 

traits.combined_raw <- rbind(ch_traits.fin, 
                         pe_traits.fin, 
                         sv_traits.fin, 
                         no_traits.fin,
                         us_traits.fin,
                         sa_traits.fin) %>% mutate(
                           species = gsub("_", " ", species), 
                           species = str_to_sentence(species),
                           trait_name = ifelse(trait_name == "mean_leaf_thickness_mm", "leaf_thickness_mm", trait_name)
                         ) %>% 
  mutate(gradient = case_when(
    grepl("South_Africa", tier) ~ "Drakensberg", 
    grepl("Peru", tier) ~ "Central Andes", 
    grepl("China", tier) ~ "Eastern Himalayas", 
    grepl("Colorado", tier) ~ "Rocky Mountains", 
    grepl("Norway", tier) ~ "Southern Scandes", 
    grepl("Svalbard", tier) ~ "Svalbard" 
  ))


mean_narm <- function(x){m <- mean(x, na.rm = T); return(m)}

traits.combined_raw[traits.combined_raw$tier == "Colorado_2016" &
    traits.combined_raw$trait_name == "n_percent", ]$trait_value



traits.stoich <- traits.combined_raw %>% 
  filter(trait_name %in% c("p_percent", "c_percent", "n_percent")) %>% 
  pivot_wider(names_from = "trait_name",
              values_from = "trait_value",
              id_cols = c("species", "plot_id", "site", "tier", "elevation",
                          "treatment", "latitude", "longitude", "burn_year", 
                          "aspect", "year"), 
              values_fn = mean_narm) %>% 
  mutate(cp_ratio = c_percent/p_percent, 
         np_ratio = n_percent/p_percent) %>% 
  dplyr::select(-c("c_percent", "n_percent", "p_percent")) %>% 
  pivot_longer(cols = c("cp_ratio", "np_ratio"), 
               names_to = "trait_name", values_to = "trait_value") %>% 
  filter(!is.na(trait_value)) %>% 
  mutate(gradient = case_when(
    grepl("South_Africa", tier) ~ "Drakensberg", 
    grepl("Peru", tier) ~ "Central Andes", 
    grepl("China", tier) ~ "Eastern Himalayas", 
    grepl("Colorado", tier) ~ "Rocky Mountains", 
    grepl("Norway", tier) ~ "Southern Scandes", 
    grepl("Svalbard", tier) ~ "Svalbard" 
  ))

summary(traits.stoich)

traits.stoich_co <- traits.combined_raw %>% 
  filter(tier == "Colorado_2016") %>% 
  filter(trait_name %in% c("p_percent", "c_percent", "n_percent")) %>% 
  pivot_wider(names_from = "trait_name",
              values_from = "trait_value",
              id_cols = c("species", "plot_id", "site", "tier", "elevation",
                          "treatment", "latitude", "longitude", "burn_year", 
                          "aspect", "year"), 
              values_fn = mean_narm) %>% 
  mutate(cn_ratio = c_percent/n_percent) %>% 
  dplyr::select(-c("c_percent", "n_percent", "p_percent")) %>% 
  pivot_longer(cols = c("cn_ratio"), 
               names_to = "trait_name", values_to = "trait_value") %>% 
  filter(!is.na(trait_value)) %>% 
  mutate(gradient = case_when(
    grepl("South_Africa", tier) ~ "Drakensberg", 
    grepl("Peru", tier) ~ "Central Andes", 
    grepl("China", tier) ~ "Eastern Himalayas", 
    grepl("Colorado", tier) ~ "Rocky Mountains", 
    grepl("Norway", tier) ~ "Southern Scandes", 
    grepl("Svalbard", tier) ~ "Svalbard" 
  ))


traits.combined <- rbind(traits.combined_raw, traits.stoich, traits.stoich_co)

setdiff(names(traits.combined_raw), names(traits.stoich))
setdiff(names(sv_traits.fin), names(no_traits.fin))

unique(traits.combined$trait_name)
table(traits.combined$gradient)

fwrite(traits.combined, "data/processed_data/preliminary_data/prelim_traits.csv")
## cover --------------------------

#china
ch_cover_final <- ch_bio %>% dplyr::select(site, plot_id, cover, species, tier) %>% 
  rename(site_id = site) 

#peru
pe_cover_final <- pe_cover %>% dplyr::select(site, plot_id, cover, species, tier) %>% 
  rename(site_id = site)

#svalbard
###  taking the mean cover of all the previous years 
sv_cover_final <- sv_cover %>% 
  rename(site_id = site) %>% 
  group_by(site_id, plot_id, species) %>% 
  summarize(cover = mean(cover, na.rm = T)) %>% 
  mutate(tier = "Svalbard_2018", 
         species = str_to_sentence(species))
table(sv_cover_final$tier)

# norway  
no_cover_final <- no_cover %>% dplyr::select(site, plot_id, cover, species, tier) %>% 
  rename(site_id = site)

# colorado 
us_cover_final <- us_cover 

#south africa
sa_cover_final <- sa_cover %>% 
  mutate(site = paste0(elevation_m_asl, aspect)) %>% dplyr::select(site, plot_id, cover, species, tier) %>% 
  rename(site_id = site) 


# Combine multiple data frames (ch_cover_final, pe_cover_final, sv_cover_final, sa_cover_final)
# into one data frame 'cover_all_raw'. Then, rename 'site_id' column to 'site'.
# Group the data by 'site', 'plot_id', 'tier', and 'species', and calculate
# the mean cover value for each group
cover_all_raw <- rbind(ch_cover_final, pe_cover_final, sv_cover_final, no_cover_final, sa_cover_final) %>% 
  rename(site = site_id) %>% 
  group_by(tier, site, plot_id, species) %>% 
  summarize(cover = mean(cover, na.rm = T))

# Calculate the mean cover for each tier and store it in 'tier.means'.
tier.means <- cover_all_raw %>% 
  group_by(tier) %>% 
  summarize(cover_tm = median(cover, na.rm = T))

# Select unique combinations of 'site', 'plot_id', 'tier', and 'species' from 'traits.combined' dataset_
all.trait_species <- traits.combined %>% 
  dplyr::select(site, plot_id, tier, species) %>% unique()

# Merge 'all.trait_species' with 'cover_all_raw' by matching columns, keeping all unique rows.
# Then merge the result with 'tier.means'. Replace NA values in 'cover' with 'cover_tm' (mean cover per tier)
# or with 1 if 'cover_tm' is also NA.
all.trait_species.comb <- all.trait_species %>%
  left_join(cover_all_raw) %>%
  unique() %>% 
  left_join(tier.means) %>% 
  mutate(cover = ifelse(is.na(cover), cover_tm, cover)) %>% 
  mutate(cover = ifelse(is.na(cover), 1, cover),
         species = gsub("_", " ", species), 
         species = str_to_sentence(species))

# Output summary statistics of the final data frame 'all.trait_species.comb'.
summary(all.trait_species.comb)

fwrite(all.trait_species.comb, "data/processed_data/preliminary_data/prelim_cover.csv")

## species richness 

dt_sp <- all.trait_species.comb %>% 
  group_by(tier, plot_id) %>% 
  summarize(SpeciesRichness = n()) %>% 
  mutate(VegPlotSizeM2 = case_when(
    grepl("China", tier) ~ 0.25*0.25, 
    grepl("Colorado", tier) ~ 0.50*0.50, 
    grepl("Norway", tier) ~ 0.50*0.50, 
    grepl("Peru", tier) ~ 1.2*1.2, 
    grepl("South_Africa", tier) ~ 1.2*1.2, 
    grepl("Svalbard", tier) ~ 0.75*0.75, 
  )) 

fwrite(dt_sp, "data/processed_data/preliminary_data/prelim_species_richness.csv")


## height n cover ---------

#china
ch_ch <- ch_bio %>% 
  group_by(tier, site, plot_id) %>% 
  summarize(
    vegHeight = sum(height * cover, na.rm = TRUE) / sum(cover, na.rm = TRUE), 
    coverSum = sum(cover, na.rm = T))

#peru
pe_height

pe_cover2 <- pe_cover %>% 
  left_join(pe_height) %>% 
  group_by(tier, site, plot_id) %>%
  summarize(coverSum = sum(cover))
  
  
pe_ch <- pe_cover2 %>% 
  left_join(pe_height) %>% 
  rename(vegHeight = height) %>% 
  dplyr::select(-treatment)

#svalbard 

sv_ch <- sv_cover %>% 
  left_join(sv_height) %>%
  group_by(tier, site, plot_id) %>% 
  summarize(coverSum = sum(cover, na.rm = T), 
            vegHeight = mean(height, na.rm = T)) %>% 
  filter(tier %in% c("Svalbard_2015", "Svalbard_2018")) %>% 
  mutate(tier = "Svalbard_2018")
unique(sv_ch$tier)
## norway

no_ch <- no_cover_ch %>% 
  left_join(no_height) %>% 
  rename(height = meanHeight) %>%
  group_by(tier, site, plot_id) %>% 
  summarize(coverSum = sum(cover, na.rm = T), 
            vegHeight = mean(height, na.rm = T))

## colorado 

us_ch <- us_height %>% rename(
  coverSum = vegCover,
)

## South Africa 

sa_ch <- sa_cover %>% 
  mutate(site = paste0(elevation_m_asl, aspect)) %>% 
  left_join(sa_height) %>% 
  rename(height = meanHeight) %>%
  group_by(tier, site, plot_id) %>% 
  summarize(coverSum = sum(cover, na.rm = T), 
            vegHeight = mean(height, na.rm = T)) %>% 
  mutate(site = as.character(site))
 
## combine
allCovHeight <- rbind(ch_ch, pe_ch, sv_ch, no_ch, us_ch, sa_ch) %>% 
  mutate(HeightXCover = vegHeight*coverSum)

fwrite(allCovHeight, "data/processed_data/preliminary_data/prelim_coverXheight.csv")


# ## spatial autocorrelation --------------------------------
# 
# sac <- flux_sf2 %>% dplyr::select(plot_id) %>% group_by(plot_id) %>% summarize() #%>% as.data.table() %>% mutate(geometry = NULL) %>% unique()
# 
# ### spatial predictors: ------
# library(spatialRF)
# 
# #The euclidian distance matrix:
# sf_use_s2(FALSE)
# c1 <- st_centroid(sac)
# coords <- st_coordinates(c1)
# c2 <- cbind(c1, coords)
# 
# distance.matrix <-as.matrix(dist(cbind(c2$X, c2$Y)))
# diag(distance.matrix) <- 0 #ged rid of diagonal
# distance.thresholds <- unname(round(quantile(distance.matrix, c(seq(0.05, .95, 0.05))), 1))
# 
# #several distances
# mems <- spatialRF::mem_multithreshold(
#   distance.matrix = distance.matrix,
#   distance.thresholds = distance.thresholds
# )
# 
# 
# # rank by moran's I
# 
# mem.rank <- spatialRF::rank_spatial_predictors(
#   distance.matrix = distance.matrix,
#   spatial.predictors.df = mems,
#   ranking.method = "moran"
# )
# 
# #order the data frame
# mems2 <- mems[, mem.rank$ranking]
# head(mems2)
# 
# ## add spatial predictors
# sac$spatial_predictor1 <- mems2[,1]
# sac$spatial_predictor2 <- mems2[,2]
# sac$spatial_predictor3 <- mems2[,3]
# sac$spatial_predictor4 <- mems2[,4]
# sac$spatial_predictor5 <- mems2[,5]
# 
# 
# sps <- sac %>%
#   as.data.table() %>%
#   mutate(geometry = NULL) %>%
#   dplyr::select(contains("spatial_predictor"), plot_id) %>%
#   unique()
# 
# unique(sps$plot_id)
# 
# fwrite(sps, "data/processed_data/preliminary_data/spatial_predictors.csv")
