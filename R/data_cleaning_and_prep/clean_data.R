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

ch.meta <- fread("data/raw_data/china/metaCH.csv") %>% 
  rename(site = Site, 
         elevation = Elevation, 
         latitude = Latitude, 
         longitude = Longitude) %>% 
  dplyr::select(-c(Gradient, Country))
#unfortunately only site metadata, not plotwise


## Flux data --------------------

#Data that I've found on OSF which shouldn't be there though. 
ch.flux.2016.2 <- get(load("data/raw_data/china/standardControlFluxCH_2016.Rdata"))
glimpse(ch.flux.2016.2)
# odd. let's ignore that for now and continue with the data we got from Inge

# Data directly from Inge

#### 2015 --------------
ch.flux.2015 <- fread("data/raw_data/china/PFTC_CO2flux_2015.txt") # Data from Inge

ggplot(data = ch.flux.2015) +
  geom_point(aes(x=NEE_exp, y = NEE_lm)) +
  geom_abline() #interesting. It looks like NEE exp is constantly higher than Lm. perhaps more conservative to go for lm?

table(ch.flux.2015$type) # ah, only respiration anyways. 
#Given that we don't have an objective way to select the better flux (R2 is of course always better for the non-linear), 
#we should probably leave it out for now 

#### 2016 ---------
# https://www.nature.com/articles/s41597-020-0529-0#Sec14
ch.flux.2016.raw <- fread("data/raw_data/china/PFTC_CO2flux_all_limits.txt") # Data from Inge
glimpse(ch.flux.2016.raw)
table(ch.flux.2016.raw$site)
table(ch.flux.2016.raw$type)
table(ch.flux.2016.raw$block)
table(ch.flux.2016.raw$datetime)

quantile(ch.flux.2016.raw$PAR, na.rm = T)

names(ch.flux.2016.raw)
unique(ch.flux.2016.raw$type)
ch.flux.2016 <- ch.flux.2016.raw %>% 
  mutate(flux_best = NEE_lm,
         tier = "China_2016", 
         elevation = as.numeric(gsub("elev", "", site))) %>%
  dplyr::select(c(datetime, type, elevation, treatment, block, temp.C, flux_best, tier, PAR)) %>% 
  left_join(ch.meta) %>% 
  mutate(plot_id = paste0("CH_",site, block), 
         type = ifelse(type == "photo", "nee", "reco")) %>% 
  rename(temperature = temp.C, 
         plot = block, 
         par = PAR) %>% 
  filter(treatment %in% c("c", "otc"))
  
table(ch.flux.2016$block)
table(ch.flux.2016$plot_id)

table(ch.flux.2016$treatment)

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



ggplot(data = ch.flux.2016) + 
  geom_boxplot(aes(x = type, y = flux_best)) +
  geom_jitter(aes(x = type, y = flux_best)) #This is weird. why all the negative fluxes in the photo part?
#also, I guess they followed the convention that negative fluxes represent outgonig CO2 while positive ones is uptake

ggplot(data = ch.flux.2016[type == "nee"]) + 
  geom_jitter(aes(x = temperature, y = flux_best)) +
  ylim(-20, 20) +
  geom_smooth(aes(x = temperature, y = flux_best), method = "lm") # no strong relationship 
##


## Traits ---------------------------------

#chemical 
ch.cht.raw <- fread("data/raw_data/china/PFTC1.2_China_2015_2016_ChemicalTraits.csv") 

ch.cht <- ch.cht.raw %>% 
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
  left_join(ch.meta)


#leaves 
ch.lt.raw <- fread("data/raw_data/china/PFTC1.2_China_2015_2016_LeafTraits.csv")
names(ch.lt.raw)
#leaf number 
ch.lt <- ch.lt.raw %>% 
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
  left_join(ch.meta)


names(ch.cht)
names(ch.lt)

quantile(ch.lt.raw$SLA_cm2_g, na.rm = T )
#bind
ch.traits.raw <- rbind(ch.lt, ch.cht)

ch.traits <- ch.traits.raw %>% 
  mutate(treatment = tolower(treatment)) %>% 
  filter(treatment %in% c("c", "otc"))
unique(ch.traits$treatment)

## Cover and Biomass ----------------------------------------
ch.bio.raw <- fread("data/raw_data/china/China_2016_Biomass_cleaned.csv")
table(ch.bio.raw$site)
table(ch.bio.raw$plot)
table(ch.flux.2016$plot_id) #cool - we have plot level data here, but plot naming is not consistent between traits and fluxes
table(ch.traits$plot_id)

unique(ch.bio.raw$plot)
unique(ch.traits$plot_id)

ch.bio <- ch.bio.raw %>% 
  mutate(plot_gsub = gsub("I-", "", plot), 
         plot_id = paste0("CH_", site, plot_gsub), 
         tier = "China_2016") %>% 
  rename(species = speciesName) %>% 
  dplyr::select(-plot_gsub)


ch.tr <- ch.traits %>% 
  dplyr::select(plot_id, treatment) %>% unique()

ch.height.trait <- ch.bio %>% 
  dplyr::select(species, plot_id, site, tier,
                trait_value = height) %>% 
  mutate(trait_name = "plant_height_cm", 
         leaf_number = NA, 
         treatment = NA, date = NA) %>% 
  left_join(ch.meta)
  
  
## Summary -------------------------------

ch.flux.2016
setdiff(names(ch.traits), names(ch.height.trait))
ch.traits <- rbind(ch.traits, ch.height.trait) %>% unique()
ch.bio #no reasonable plot ID either... 

unique(ch.flux.2016$plot_id)
unique(ch.bio$plot_id)
unique(ch.traits$plot_id)




# Peru  ------------------------------
#https://www.nature.com/articles/s41597-024-02980-3

pe.meta.raw <- fread("data/raw_data/peru/PU.10_PFTC3.10_2020_Peru_Coordinates.csv")

names(pe.meta.raw) <- tolower(names(pe.meta.raw))

pe.meta <- pe.meta.raw %>% 
  rename(plot = plotid) %>% 
  filter(!is.na(plot)) %>%
  mutate(plot_id = paste0("PE_", site, "_", treatment ,"_",  plot)) %>% dplyr::select(-comment) %>% 
  dplyr::select(plot_id, elevation, longitude, latitude, burn_year, plot)


## Flux data ----------------------------

pe.flux.raw <- fread("data/raw_data/peru/PFTC3_Puna_PFTC5_Peru_2018_2020_Cflux.csv")
names(pe.flux.raw)
unique(pe.flux.raw[,.(site, treatment)])


pe.flux <- pe.flux.raw %>% 
  rename(plot = plot_id, 
         type = flux, 
         temperature = t_ave) %>%  #some temperatures (if this column even is temperature) are oddly low: quantile(pe.flux.raw$t_ave, c(0.05, 0.95)) 
                                    # no way it was -60°C there. fix below
  mutate(
    tier = paste0("Peru_", year), 
    flux_best = linear_model,
    treatment = ifelse(site == "TRE" & treatment == "B", "NB", treatment), ###--> to be able to use the metadata as there seem to be a discrepancy between the treatment. It's NB in metadata and B in the flux data 
    plot_id = paste0("PE_", site, "_", treatment ,"_",  plot),
    temperature = ifelse(temperature < -10, NA, temperature),
    date = paste0(day, "-", month, "-", year), 
    type = tolower(type), 
    type = ifelse(grepl("nee", type), "nee", type), 
    treatment = tolower(treatment), 
    par = NA, 
    date = as.Date(date, format = "%d-%B-%Y"), 
    date = format(date)) %>% 
  dplyr::select(c(type, temperature, plot_id, flux_best, tier, date, site, treatment)) %>% 
  left_join(pe.meta) 


ggplot(data = pe.flux) + 
  geom_boxplot(aes(x = type, y = flux_best)) +
  geom_hline(yintercept = 0) +
  geom_jitter(aes(x = type, y = flux_best)) #CO2 uptake = positive, release = negative. But wtf is NEE1, NEE2, NEE3?


ggplot(data = pe.flux[grepl("nee", type), ]) + 
  geom_jitter(aes(x = temperature, y = flux_best)) +
  geom_smooth(aes(x = temperature, y = flux_best), method = "lm") #nothing going on here?


## Trait data ----------------------------

pe.t.raw <- fread("data/raw_data/peru/PFTC3-Puna-PFTC5_Peru_2018-2020_FunctionalTraits_clean.csv")

unique(pe.t.raw$trait)
quantile(pe.t.raw[pe.t.raw$trait == "sla_cm2_g", ]$value , na.rm = T )

pe.traits <- pe.t.raw %>% 
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
  left_join(pe.meta)  %>% dplyr::select(-plot)

unique(pe.traits$trait_name)
unique(pe.traits[,.(site, treatment)])


## Biomass/Height -----------------------------

pe.height <- fread("data/raw_data/peru/PFTC3-Puna-Peru_2018-2019_CommunityStructure_clean.csv") %>% 
  filter(variable == "median_height") %>% 
  mutate(
    tier = paste0("Peru_", year)) %>% 
  dplyr::select(site, tier, treatment, value, plot_id) %>% 
  rename(height = value) %>% 
  mutate(plot_id = paste0("PE_", site, "_", treatment ,"_",  plot_id)) 


## Cover -----------------------------

pe.cover.raw <- fread("data/raw_data/peru/PFTC3-Puna-PFTC5_Peru_2018-2020_CommunityCover_clean.csv")

pe.cover <- pe.cover.raw %>% 
  mutate(tier = paste0("Peru_", year),
         plot_id = paste0("PE_", site, "_", treatment ,"_",  plot_id)) %>% 
  rename(species = taxon)


## Summary -----------------------------

pe.flux
pe.traits
#pe.bio #no plot ID... # but site - prob. better than nothing 
pe.cover
unique(pe.flux$plot_id)
unique(pe.traits$plot_id)
unique(pe.cover$plot_id)

# Svalbard  ------------------------------

sv.meta.raw <- fread("data/raw_data/svalbard/PFTC4_Svalbard_2018_metaItex.csv")
sv.coords.raw <- fread("data/raw_data/svalbard/PFTC4_Svalbard_Coordinates_ITEX.csv")
sv.coords.grad <- fread("data/raw_data/svalbard/PFTC4_Svalbard_Coordinates_Gradient.csv") %>% 
  mutate(plot_id = paste0("SV_", Gradient, Site, PlotID)) %>% 
  rename(elevation = Elevation_m,
         longitude = Longitude_E,
         latitude = Latitude_N) %>% 
  select(plot_id, elevation, latitude, longitude)

sv.meta <- sv.meta.raw %>% 
  left_join(sv.coords.raw) %>% 
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

names(sv.meta) <- tolower(names(sv.meta))


## Flux data ----------------------------
#https://www.nature.com/articles/s41597-023-02467-7

#sv.flux.raw <- fread("data/raw_data/svalbard/Cflux_SV_ITEX_2018.csv") 
sv.flux.raw <- get(load("data/raw_data/svalbard/ITEX_all.Rdata"))
glimpse(sv.flux.raw)
names(sv.flux.raw)  <- tolower(names(sv.flux.raw))
#a first glimpse suggests that the weather was shit all the time and the PAR was unreasonable low. 
#we may have to exclude Svalbard
glimpse(sv.flux.raw)
sv.flux <- sv.flux.raw %>% 
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
  left_join(sv.meta, by = c("plot_id")) %>% as.data.table() 

setdiff(sv.meta$plot_id, sv.flux$plot_id)
setdiff(sv.flux$plot_id, sv.meta$plot_id)

# Svalbard gradients

sv.flux2 <- fread("data/raw_data/svalbard/Cflux_SV_Gradient_2018.csv") %>% 
  setNames(tolower(names(.))) %>% 
  rename(flux_best = nee,
         type = cover,
         par = par_mean,
         temperature = ir_mean) %>% 
  select(date, site, plotid, gradient, starttime, type, par, temperature, flux_best, rsqd) %>% 
  mutate(type = ifelse(type == "L", "nee", "reco")) %>% 
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
  left_join(., sv.coords.grad)

sv.flux <- bind_rows(sv.flux, sv.flux2 %>% mutate(date = dmy(date),
                                                  site = as.character(site)))

ggplot(data = sv.flux) + 
  geom_boxplot(aes(x = type, y = flux_best)) +
  geom_hline(yintercept = 0) +
  geom_jitter(aes(x = type, y = flux_best)) #CO2 uptake = positive, release = negative.

ggplot(data = sv.flux) + 
  geom_boxplot(aes(x = treatment, y = temperature)) +
  geom_hline(yintercept = 0) +
  geom_jitter(aes(x = treatment, y = temperature)) #CO2 uptake = positive, release = negative.


## Trait data ----------------------------

### ITEX 
sv.t.raw <- fread("data/raw_data/svalbard/PFTC4_Svalbard_2018_ITEX_Traits.csv")
names(sv.t.raw)  <- tolower(names(sv.t.raw))
table(sv.t.raw$trait)


sv.traits.itex <- sv.t.raw %>% 
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
sv.traits.itex

## Gradient 

sv.t.grad.raw <- fread("data/raw_data/svalbard/PFTC4_Svalbard_2018_Gradient_Traits.csv")

sv.gradient.t <- sv.t.grad.raw %>% 
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

sv.traits <- rbind(sv.traits.itex, sv.gradient.t)

## Biomass/Height -----------------------------

# can't find any biomass for svalbard (which i guess makes sense)

sv.height.raw <- fread("data/raw_data/svalbard/PFTC4_Svalbard_2003_2015_ITEX_Community_Structure.csv") %>% 
  filter(Year == 2015)
names(sv.height.raw)  <- tolower(names(sv.height.raw))
sv.height.itex <- sv.height.raw %>% 
  mutate(plot = gsub("CH", "CAS", plotid), 
         plot = gsub("SB", "BIS", plot), 
         plot = gsub("DH", "DRY", plot), 
         plot = gsub("-", "_L", plot), 
         site = "ITEX",
         treatment = ifelse(treatment == "CTL", "c", "otc"), 
         tier = paste0("Svalbard_", year), 
         plot_id = paste0("SV_", site, plot)) %>% dplyr::select(plot_id, height)

sv.height.raw.gr <- fread("data/raw_data/svalbard/PFTC4_Svalbard_2018_Community_Structure_Gradient.csv") 
names(sv.height.raw.gr)  <- tolower(names(sv.height.raw.gr))

sv.gradient.height <- sv.height.raw.gr %>% 
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

sv.height <- rbind(sv.height.itex, sv.gradient.height)

## Cover -----------------------

sv.cover.raw <- fread("data/raw_data/svalbard/PFTC4_Svalbard_2003_2015_ITEX_Community.csv")
names(sv.cover.raw)  <- tolower(names(sv.cover.raw))
sv.cover.itex <- sv.cover.raw %>% 
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

sv.cover.g.raw <- fread("data/raw_data/svalbard/PFTC4_Svalbard_2018_Community_Gradient.csv")
names(sv.cover.g.raw)  <- tolower(names(sv.cover.g.raw))

sv.cover.gradient <- sv.cover.g.raw %>% 
  mutate(
    tier = paste0("Svalbard_", year), 
    plot_id = paste0("SV_", gradient, site, plotid),
    treatment = "c", 
    site = paste0(gradient, site), 
    species = str_to_sentence(taxon)) %>% 
  dplyr::select(site, plot_id, cover, species, tier)

sv.cover <- rbind(sv.cover.gradient, sv.cover.itex)

## Summary -----------------------------

sv.flux
sv.traits
sv.cover

unique(sv.flux$plot_id)
unique(sv.traits$plot_id)
unique(sv.cover$plot_id)


# Norway ------------------------------
#https://docs.google.com/document/d/1nXt4sljpExC_fSIVvGJ6-bDw9cAKANp6cwCWNklOHUU/edit
## meta data missing 

## Flux data ----------------------------

no.flux.raw <- fread("data/raw_data/norway/PFTC6_24h_cflux_allsites_2022.csv") 
quantile(no.flux.raw$PARavg, na.rm = T)
no.flux.raw2 <- no.flux.raw %>% 
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


glimpse(no.flux.raw)
table(no.flux.raw$flag)
hist(no.flux.raw$PARavg)

no.flux <- no.flux.raw2 %>% 
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
         
unique(no.flux[,c("site", "turfID")])

ggplot(data = no.flux) + 
  geom_boxplot(aes(x = type, y = flux_best)) +
  geom_hline(yintercept = 0) +
  geom_jitter(aes(x = type, y = flux_best)) #CO2 uptake = positive, release = negative


ggplot(data = no.flux[grepl("nee", type), ]) + 
  geom_jitter(aes(x = temperature, y = flux_best)) +
  geom_smooth(aes(x = temperature, y = flux_best), method = "lm")

## Trait data ----------------------------


no.t.threed.raw <- fread("data/raw_data/norway/PFTC6_ThreeD_clean_leaf_traits_2022.csv") 
names(no.t.threed.raw)
unique(no.t.threed.raw$siteID)
unique(no.t.threed.raw$Namount_kg_ha_y)

unique(no.t.threed.raw$destSiteID)

unique(no.t.threed.raw[siteID %in% c("Vikesland") & warming == "A", ]$blockID)

unique(no.t.threed.raw[siteID %in% c("Vikesland") & warming == "W", .(blockID, turfID)])


unique(no.flux[site %in% c("Vik"), ]$turfID)

quantile(no.t.threed.raw[no.t.threed.raw$trait == "sla_cm2_g", ]$value , na.rm = T )


no.traits.raw <- no.t.threed.raw %>% 
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
  
unique(no.traits.raw$site)
unique(no.traits.raw$elevation)
unique(no.flux$plot_id)



#no.traits %>% dplyr::select(site, treatment, turfID) %>% unique()
#unique(no.traits$trait_name)



## Biomass/Height -----------------------------

## nothing available on OSF: https://osf.io/fcbw4/
## get height the dodgy way (from traits)
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
  filter(plot_id %in% unique(no.traits.raw$plot_id)) %>% 
  dplyr::select(siteID, species, plot_id, cover) %>% 
  rename(site = siteID) 
threeD_community
vcg_community


no.cover <- rbind(threeD_community, vcg_community) %>% 
  mutate(plot_id = ifelse(site %in% c("Hog", "Vik"), NA, plot_id), 
         tier = paste0("Norway_2022"))

no.cover.ch <- rbind(threeD_community, vcg_community) %>% 
  mutate(tier = paste0("Norway_2022"))

summary(no.cover)
unique(no.cover$site)


## Summary -----------------------------

no.flux
no.traits <- no.traits.raw %>% 
  left_join(turfBlockLeg) %>% 
  mutate(plot_id = ifelse(is.na(plot_id), plot_id_vcg, plot_id)) %>% #
  dplyr::select(-plot_id_vcg)
no.traits     
no.cover

## Height --------
no.median.cover <- no.cover.ch %>% 
  group_by(plot_id) %>% 
  summarize(medianCov = median(cover, na.rm = T))

no.height <- no.traits %>% 
  filter(trait_name == "plant_height_cm") %>% 
  left_join(no.cover.ch[, c("plot_id", "cover", "species")]) %>% 
  left_join(no.median.cover) %>% 
  mutate(cover = ifelse(is.na(cover), medianCov, cover)) %>% 
  group_by(plot_id) %>% 
  summarize(
    meanHeight = sum(trait_value * cover, na.rm = TRUE) / sum(cover, na.rm = TRUE)
  )



# Colorado ------------------------------

us.meta.raw <- fread("data/raw_data/colorado/rmbl_site_info.csv")
us.meta <- us.meta.raw %>% 
  rename(
    latitude = lat, 
    longitude = long, 
    elevation = elev
  ) %>% 
  mutate(site = tolower(site)) %>% 
  dplyr::select(-lat_long)

## Flux data ----------------------------

us.flux.raw <- fread("data/raw_data/colorado/rmbl_gradient_flux_data_12042023.csv")

unique(us.flux.raw$time)
unique(us.flux.raw$plot)
unique(us.flux.raw$season)


unique(us.flux.raw$site)
us.flux <- us.flux.raw %>% 
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
  left_join(us.meta) %>% filter(!is.na(flux_best)) %>% 
  filter(season == "Peak")
us.flux
us.flux.raw[time == "DAY RESP"]


ggplot(data = us.flux) + 
  geom_boxplot(aes(x = type, y = flux_best)) +
  geom_hline(yintercept = 0) +
  geom_jitter(aes(x = type, y = flux_best)) #CO2 uptake = positive, release = negative.

quantile(us.flux$temperature, na.rm = TRUE)
ggplot(data = us.flux[grepl("nee", type), ]) + 
  geom_jitter(aes(x = temperature, y = flux_best)) +
  geom_smooth(aes(x = temperature, y = flux_best), method = "lm")
  
unique(us.flux$site)

## Trait data ----------------------------

us.t.raw <- fread("data/raw_data/colorado/rmbl_trait_data_master.csv")
glimpse(us.t.raw) 

unique(us.t.raw$site)
unique(us.t.raw[us.t.raw$site == "Cinnamon", ]$plot)
unique(us.t.raw[us.t.raw$site == "Cinnamon", ]$plot_name)
unique(us.t.raw[us.t.raw$site == "Cinnamon", ]$plot_code)
unique(us.t.raw[us.t.raw$site == "Almont", ]$plot)
unique(us.t.raw[us.t.raw$site == "CBT" & !is.na(plot), ]$year)
unique(us.t.raw[us.t.raw$site == "Almont" & !is.na(plot), ]$year)
unique(us.t.raw[us.t.raw$site == "Cinnamon" & !is.na(plot), ]$year)
unique(us.t.raw[us.t.raw$site == "Road" & !is.na(plot), ]$year)
unique(us.t.raw[us.t.raw$site == "PBM" & !is.na(plot), ]$year)
unique(us.t.raw[us.t.raw$site == "Pfeiler" & !is.na(plot), ]$year)


unique(us.t.raw[!is.na(us.t.raw$height), ]$plot)

#get plant height since this does not seem to be available for plots with numbers otherwise
us_pl_h <- us.t.raw %>%
  filter(!is.na(height)) %>% 
  dplyr::select(height, species) %>% 
  group_by(species) %>% 
  summarize(height = mean(height, na.rm = T))


us.traits <- us.t.raw %>% 
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
  left_join(us.meta) %>% #vast majority (like 60,000) is NA for the plot. 
  #filter(!grepl("NA", plot_id)) %>% 
  mutate(trait_value = ifelse(year == 2016 & trait_name == "leaf_area_cm2", trait_value*10, trait_value))

us.traits[us.traits$trait_name == "plant_height_cm", ]$trait_value

us.traits
unique(us.traits[us.traits$site == "cinnamon", ]$plot)

## Biomass/Height -----------------------------
#cover
us.coverH.raw <- fread("data/raw_data/colorado/veg_cover_data_rmbl.csv") %>% 
  mutate(site = tolower(site),
         plot = tolower(plot), 
         plot_id = paste("US_", site, plot), 
         plot_id = gsub(" ", "", plot_id),
         tier = paste0("Colorado_", year), 
         vegCover = (herb + shrub + grass + cactus + forb)) %>% 
  filter(tier %in% c("Colorado_2016", "Colorado_2018")) %>% 
  filter(plot_id %in% unique(us.flux$plot_id)) %>% 
  dplyr::select(plot_id, tier, vegCover, site) %>% 
  group_by(tier, site, plot_id) %>% 
  summarize(vegCover = mean(vegCover, na.rm = T))


## Height -------------
us.height.raw <- fread("data/raw_data/colorado/veg_height_data_rmbl.csv") %>% 
  mutate(site = tolower(site),
         plot = tolower(plot), 
         plot_id = paste("US_", site, plot), 
         plot_id = gsub(" ", "", plot_id),
         tier = paste0("Colorado_", year)) %>% 
  filter(tier %in% c("Colorado_2016", "Colorado_2018")) %>% 
  filter(plot_id %in% unique(us.flux$plot_id)) %>% 
  rename(vegHeight = mean) %>% 
  dplyr::select(plot_id, tier, vegHeight, site) %>% unique() %>% 
  group_by(tier, site, plot_id) %>% 
  summarize(vegHeight = mean(vegHeight, na.rm = T))

us.height <- us.height.raw %>% 
  left_join(us.coverH.raw) %>% unique()


# don't have anything yet

## Cover --------------------------

us.cover <- fread("data/raw_data/colorado/rmbl_plot_data_2016.csv") %>% 
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
  

#unique(us.cover.raw$treatment)
## Summary -----------------------------

us.flux
us.traits

unique(us.traits$plot_id)

unique(us.flux$plot_id)

# South Africa ------------------------------


##metadata 
library(sf) 
sa.wp.raw <- read_sf("data/raw_data/south_africa/PFCT7_site_waypoints.gpx") %>%
  dplyr::select(name) %>%
  filter(name != "PFTC5 E4")

sa.wp.coords <- st_coordinates(sa.wp.raw)

sa.wp <- sa.wp.raw %>% 
  cbind(sa.wp.coords) %>% 
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
calc.intermediate.plots <- function(site.data) {
  lat1 <- site.data$latitude[site.data$plot == 1]
  lon1 <- site.data$longitude[site.data$plot == 1]
  lat5 <- site.data$latitude[site.data$plot == 5]
  lon5 <- site.data$longitude[site.data$plot == 5]
  
  data.frame(
    unique_site = unique(site.data$unique_site),
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
    aspect = unique(site.data$aspect), 
    elevation = unique(site.data$elevation), 
    site = unique(site.data$site)
  )
}

# Calculate intermediate plots for each site
intermediate.plots <- sa.wp %>%
  group_by(unique_site) %>%
  do(calc.intermediate.plots(.))

# Combine original and intermediate plots
sa.wp2 <- bind_rows(sa.wp, intermediate.plots) %>%
  arrange(unique_site, plot) %>% filter(!is.na(aspect)) %>% mutate(site = paste0(elevation, aspect))

sa.wp.sf <- st_as_sf(sa.wp2, 
                     coords = c("longitude", "latitude"), 
                     crs = 4326)
mapview::mapview(sa.wp.sf)

fwrite(sa.wp2, "data/raw_data/south_africa/PFCT7_plot_locations_clean.csv")


### get dates 

sa.meta.raw <- fread("data/raw_data/south_africa/PFTC7_SA_raw_fluxes_2023 - Gradient.csv")

sa.meta <- sa.meta.raw %>% 
  rename(sampling_day = `day (NOT DATE!!!)`) %>% 
  filter(day.night == "day") %>% 
  mutate(unique_site = paste0(siteID, "_", aspect), 
         keep = case_when(
           .default = "keep",
           unique_site == "5_east" & Attempt == 1 ~ "discard", 
           unique_site == "5_west" & Attempt == 1 ~ "discard", 
           unique_site == "3_east" & Attempt == 1 ~ "discard", 
           unique_site == "3_west" & Attempt == 1 ~ "discard", 
           unique_site == "1_west" & Attempt == 1 ~ "discard", 
           unique_site == "1_east" & sampling_day == 3 ~ "discard")) %>%
  filter(keep == "keep") %>% 
  dplyr::select(unique_site, siteID, sampling_day, Attempt) %>% 
  arrange(siteID) %>% 
  unique() %>% 
  mutate(day = sampling_day + 2, 
         month = 12, 
         year = 2023, 
         date = as_date(paste0(year, "-", month, "-", day))
         ) %>% 
  dplyr::select(unique_site, date )
           
sa.meta

sa.wp3 <- sa.wp2 %>% left_join(sa.meta)


## Flux data ----------------------------
#https://docs.google.com/document/d/1P2X-3IIQE6IQwvgvDQe9YJLqVWd2srqzxJZWWoM3mig/edit
sa.flux.raw <- fread("data/raw_data/south_africa/PFTC7_licor_nee_flagged.csv")
names(sa.flux.raw)
sa.flux <- sa.flux.raw %>% 
  filter(flag %in% c("okay", "manual_flux_time_selection")) %>% 
  rename(flux_best = flux_value, 
         temperature = tav) %>%
  mutate(tier = "South_Africa_2023", 
         plot_id = paste0("SA_", elevation, aspect, plot), 
         site = paste0(elevation, aspect),
         type = case_when(
           measurement == "photo" ~ "nee", 
           measurement == "resp" & day.night == "day" ~ "reco", 
           measurement == "resp" & day.night == "night" ~ "night_reco", 
         )) %>% 
  dplyr::select(c(plot_id, tier, flux_best, type, temperature, site, plot, aspect, elevation)) %>% 
  left_join(sa.wp3)

range(sa.flux[type == "nee", flux_best], na.rm = T)
range(sa.flux[type == "nee", flux_best], na.rm = T)
range(sa.flux[type == "night_reco", flux_best], na.rm = T)

sum(!is.na(sa.flux$flux_best))

ggplot(data = sa.flux) + 
  geom_boxplot(aes(x = type, y = flux_best)) +
  geom_hline(yintercept = 0) +
  geom_jitter(aes(x = type, y = flux_best)) #CO2 uptake = positive, release = negative. But wtf is NEE1, NEE2, NEE3?

ggplot(data = sa.flux[grepl("nee", type), ]) + 
  geom_jitter(aes(x = temperature, y = flux_best)) +
  geom_smooth(aes(x = temperature, y = flux_best), method = "lm")


## Trait data ----------------------------

sa.t.raw <- fread("data/raw_data/south_africa/PFTC7_SA_clean_traits_2023.csv")
table(sa.t.raw$problem_flag)
table(sa.t.raw$traits)

sa.traits <- sa.t.raw %>% 
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
         site = paste0(elevation, aspect)) %>% 
  filter(trait_name %in% c("plant_height_cm", "wet_mass_g", "dry_mass_g", "leaf_area_cm2", "leaf_thickness_mm", "sla_cm2_g",
                        "ldmc")) %>% 
  dplyr::select(c(plot_id, site, plot, elevation, aspect, trait_name, trait_value, species, tier, year)) %>% 
  left_join(sa.wp2)
    
unique(sa.traits$species)
## Biomass -----------------------------

## I don't think we have clean biomass data yet 

## Cover --------------------------------
sa.cover.raw <- fread("data/raw_data/south_africa/PFTC7_SA_clean_community_19Apr2024.csv")
sa.cover <- sa.cover.raw %>% 
  rename(site = site_id, 
         plot = plot_id) %>% 
  mutate(tier = "South_Africa_2023", 
         aspect = ifelse(aspect == "W", "west", "east"), 
         plot = round(plot, 0), 
         plot_id = paste0("SA_", elevation, aspect, plot), 
         year = 2023, 
         species = gsub("_", " ", species), 
         species = str_to_sentence(species)) %>% 
  filter(!site == 6)

## Height ------------

sa.height.raw <- fread("data/raw_data/south_africa/PFTC7_SA_clean_community_structure_2023.csv")

sa.height <- sa.height.raw %>% 
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
sa.flux 
sa.traits

# combine everything ---------------------------------------------------

## fluxes ------------

names(ch.flux.2016)# tier, "datetime", "type", "elevation", "treatment", "plot", "temperature", "flux_best", "site", "latitude", "longitude", "plot_id" 
ch.flux.2016.fin <- ch.flux.2016 %>% 
  mutate(burn_year = NA, 
         datetime = dmy_hm(datetime),
         date = date(datetime),
         date = as_date(date, format = '%Y.%m.%d'),
         
         year = year(date), 
         
         aspect = NA, 
         season = "not_part_of_the_dataset")

names(pe.flux) #extra = burn_year
pe.flux.fin <- pe.flux %>% 
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


names(sv.flux) 
sv.flux.fin <- sv.flux %>% 
  mutate(burn_year = NA, 
         datetime = as_datetime(NA),
         aspect = NA, 
         date = as_date(date, format = '%d.%m.%Y'),
         date = as_date(date, format = '%Y.%m.%d'),
         year = year(date), 
         season = "not_part_of_the_dataset")

names(no.flux) #big mess, have to clear that up before we go on
no.flux.fin <- no.flux %>%
  mutate(burn_year = NA,
         date = date(datetime),
         date = as_date(date, format = '%Y.%m.%d'),
         year = year(date), 
         aspect = NA,
         season = "not_part_of_the_dataset") %>% rename(plot = turfID)

names(us.flux)
us.flux.fin <- us.flux %>% 
  mutate(burn_year = NA, 
         datetime = as_datetime(NA), 
         date = as_date(date), 
         date = as_date(date, format = '%Y.%m.%d'),
         year = year(date), 
         aspect = NA, 
         treatment = NA,
         par = NA)

setdiff(names(no.flux.fin), names(pe.flux.fin))
setdiff(names(pe.flux.fin), names(no.flux.fin))


names(sa.flux)
sa.flux.fin <- sa.flux %>% 
  mutate(burn_year = NA, 
         datetime = as_datetime(NA), 
         date = as_date(date, format = '%Y.%m.%d'),
         treatment = NA,
         year = year(date), 
         par = NA,
         season = "not_part_of_the_dataset") %>% 
  filter(!type == "night_reco") %>% 
  dplyr::select(-unique_site)

fluxes.combined.raw <- rbind(ch.flux.2016.fin, 
                         pe.flux.fin, 
                         sv.flux.fin, 
                         no.flux.fin,
                         us.flux.fin, 
                         sa.flux.fin)

setdiff(names(ch.flux.2016.fin), names(sa.flux.fin))
setdiff(names(sa.flux.fin), names(ch.flux.2016.fin))


unique(fluxes.combined.raw[grepl("Peru", tier), ]$treatment)

fluxes.combined <- fluxes.combined.raw %>% 
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

flux.sf <- st_as_sf(fluxes.combined2, 
                    coords = c("longitude", "latitude"), 
                    crs = 4326)

mapview::mapview(flux.sf)

flux.sf2 <- flux.sf %>% 
  dplyr::select(plot_id, tier, date) %>% 
  mutate(date = paste0(date), 
         flux_id = paste0(plot_id, 1:nrow(.))) %>% st_make_valid()
str(flux.sf2)
mapview::mapview(flux.sf2)
n_distinct(flux.sf2$plot_id)
n_distinct(flux.sf2$flux_id)

st_write(flux.sf2, "data/processed_data/preliminary_data/prelim_flux_loc.gpkg", append = FALSE)
st_write(flux.sf2, "data/processed_data/preliminary_data/prelim_flux_loc.shp", append = FALSE)

mapview::mapview(flux.sf2 %>% filter(plot_id == "PE_ACJ_C_1"))


## traits ---------------------------------

names(ch.traits) #get rid of date and replace w year because who cares

ch.traits.fin <- ch.traits %>% 
  mutate(burn_year = NA, 
         aspect = NA, 
         date = as_date(date),
         year = year(date)) %>% dplyr::select(-date, -leaf_number)

names(pe.traits)
pe.traits.fin <- pe.traits %>%
  mutate(aspect = NA,
         tier = paste0("Peru_", year)) %>% 
  dplyr::select(-c(leaf_number))

names(sv.traits) #
sv.traits.fin <- sv.traits %>%
  mutate(aspect = NA, 
         burn_year = NA)
unique(sv.traits$tier)
# names(no.traits) #not yet finished 
no.traits.fin <- no.traits %>%
  mutate(aspect = NA, 
         burn_year = NA) %>% dplyr::select(-datetime, -date, -blockID)

names(us.traits) # get rd of 
us.traits.fin <- us.traits %>%
     mutate(aspect = NA, 
            burn_year = NA, 
            treatment = NA) %>% dplyr::select(-plot)

names(sa.traits) # get rd of 

sa.traits.fin <- sa.traits %>% 
  mutate(burn_year = NA, 
         treatment = NA) %>% dplyr::select(-plot, - unique_site)

names(ch.traits.fin)
      names(pe.traits.fin) 
      names(sv.traits.fin)  
      names(no.traits.fin) 
      names(us.traits.fin) 
      names(sa.traits.fin) 

traits.combined.raw <- rbind(ch.traits.fin, 
                         pe.traits.fin, 
                         sv.traits.fin, 
                         no.traits.fin,
                         us.traits.fin,
                         sa.traits.fin) %>% mutate(
                           species = gsub("_", " ", species), 
                           species = str_to_sentence(species),
                           trait_name = ifelse(trait_name == "mean_leaf_thickness_mm", "leaf_thickness_mm", trait_name)
                         )

mean_narm <- function(x){m <- mean(x, na.rm = T); return(m)}

traits.stoich <- traits.combined.raw %>% 
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
  filter(!is.na(trait_value))

summary(traits.stoich)
traits.combined <- rbind(traits.combined.raw, traits.stoich)

setdiff(names(traits.combined.raw), names(traits.stoich))
setdiff(names(sv.traits.fin), names(no.traits.fin))

unique(traits.combined$trait_name)
fwrite(traits.combined, "data/processed_data/preliminary_data/prelim_traits.csv")
## cover --------------------------

#china
ch.cover.final <- ch.bio %>% dplyr::select(site, plot_id, cover, species, tier) %>% 
  rename(site_id = site) 

#peru
pe.cover.final <- pe.cover %>% dplyr::select(site, plot_id, cover, species, tier) %>% 
  rename(site_id = site)

#svalbard
###  taking the mean cover of all the previous years 
sv.cover.final <- sv.cover %>% 
  rename(site_id = site) %>% 
  group_by(site_id, plot_id, species) %>% 
  summarize(cover = mean(cover, na.rm = T)) %>% 
  mutate(tier = "Svalbard_2018", 
         species = str_to_sentence(species))
table(sv.cover.final$tier)

# norway  
no.cover.final <- no.cover %>% dplyr::select(site, plot_id, cover, species, tier) %>% 
  rename(site_id = site)

# colorado 
us.cover.final <- us.cover 

#south africa
sa.cover.final <- sa.cover %>% 
  mutate(site = paste0(elevation, aspect)) %>% dplyr::select(site, plot_id, cover, species, tier) %>% 
  rename(site_id = site) 


# Combine multiple data frames (ch.cover.final, pe.cover.final, sv.cover.final, sa.cover.final)
# into one data frame 'cover.all.raw'. Then, rename 'site_id' column to 'site'.
# Group the data by 'site', 'plot_id', 'tier', and 'species', and calculate
# the mean cover value for each group
cover.all.raw <- rbind(ch.cover.final, pe.cover.final, sv.cover.final, no.cover.final, sa.cover.final) %>% 
  rename(site = site_id) %>% 
  group_by(tier, site, plot_id, species) %>% 
  summarize(cover = mean(cover, na.rm = T))

# Calculate the mean cover for each tier and store it in 'tier.means'.
tier.means <- cover.all.raw %>% 
  group_by(tier) %>% 
  summarize(cover_tm = median(cover, na.rm = T))

# Select unique combinations of 'site', 'plot_id', 'tier', and 'species' from 'traits.combined' dataset.
all.trait.species <- traits.combined %>% 
  dplyr::select(site, plot_id, tier, species) %>% unique()

# Merge 'all.trait.species' with 'cover.all.raw' by matching columns, keeping all unique rows.
# Then merge the result with 'tier.means'. Replace NA values in 'cover' with 'cover_tm' (mean cover per tier)
# or with 1 if 'cover_tm' is also NA.
all.trait.species.comb <- all.trait.species %>%
  left_join(cover.all.raw) %>%
  unique() %>% 
  left_join(tier.means) %>% 
  mutate(cover = ifelse(is.na(cover), cover_tm, cover)) %>% 
  mutate(cover = ifelse(is.na(cover), 1, cover),
         species = gsub("_", " ", species), 
         species = str_to_sentence(species))

# Output summary statistics of the final data frame 'all.trait.species.comb'.
summary(all.trait.species.comb)

fwrite(all.trait.species.comb, "data/processed_data/preliminary_data/prelim_cover.csv")

## species richness 

dt.sp <- all.trait.species.comb %>% 
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

fwrite(dt.sp, "data/processed_data/preliminary_data/prelim_species_richness.csv")


## height n cover ---------

#china
ch.ch <- ch.bio %>% 
  group_by(tier, site, plot_id) %>% 
  summarize(
    vegHeight = sum(height * cover, na.rm = TRUE) / sum(cover, na.rm = TRUE), 
    coverSum = sum(cover, na.rm = T))

#peru
pe.height

pe.cover2 <- pe.cover %>% 
  left_join(pe.height) %>% 
  group_by(tier, site, plot_id) %>%
  summarize(coverSum = sum(cover))
  
  
pe.ch <- pe.cover2 %>% 
  left_join(pe.height) %>% 
  rename(vegHeight = height) %>% 
  dplyr::select(-treatment)

#svalbard 

sv.ch <- sv.cover %>% 
  left_join(sv.height) %>%
  group_by(tier, site, plot_id) %>% 
  summarize(coverSum = sum(cover, na.rm = T), 
            vegHeight = mean(height, na.rm = T)) %>% 
  filter(tier %in% c("Svalbard_2015", "Svalbard_2018")) %>% 
  mutate(tier = "Svalbard_2018")
unique(sv.ch$tier)
## norway

no.ch <- no.cover.ch %>% 
  left_join(no.height) %>% 
  rename(height = meanHeight) %>%
  group_by(tier, site, plot_id) %>% 
  summarize(coverSum = sum(cover, na.rm = T), 
            vegHeight = mean(height, na.rm = T))

## colorado 

us.ch <- us.height %>% rename(
  coverSum = vegCover,
)

## South Africa 

sa.ch <- sa.cover %>% 
  mutate(site = paste0(elevation, aspect)) %>% 
  left_join(sa.height) %>% 
  rename(height = meanHeight) %>%
  group_by(tier, site, plot_id) %>% 
  summarize(coverSum = sum(cover, na.rm = T), 
            vegHeight = mean(height, na.rm = T)) %>% 
  mutate(site = as.character(site))
 
## combine
allCovHeight <- rbind(ch.ch, pe.ch, sv.ch, no.ch, us.ch, sa.ch) %>% 
  mutate(HeightXCover = vegHeight*coverSum)

fwrite(allCovHeight, "data/processed_data/preliminary_data/prelim_coverXheight.csv")


# ## spatial autocorrelation --------------------------------
# 
# sac <- flux.sf2 %>% dplyr::select(plot_id) %>% group_by(plot_id) %>% summarize() #%>% as.data.table() %>% mutate(geometry = NULL) %>% unique()
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
