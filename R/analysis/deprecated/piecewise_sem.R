library("lavaan")
library("piecewiseSEM")
library(data.table)
library(tidyverse)
library(lme4)
library(MuMIn)
library(glmmTMB)
library("semPlot")
library(performance)
library(gridExtra)


dt_raw <- fread("data/processed_data/clean_data/global_fluxes_main_data.csv") %>% 
  dplyr::select(
    # identifiers 
    country, site, plot_id,# treatment,
    
    # fluxes
    nee, reco, gpp,
    
    # environmental 
    elevation, map, mat,
    temperature_nee, temperature_reco, temperature_gpp,
    
    # trait means 
    sla_cm2_g, ldmc, leaf_area_cm2, dry_mass_g, plant_height_cm,
    n_percent, cn_ratio, c_percent, p_percent,
    
    # others 
    height_x_cover, species_richness,
    functional_diversity_q1, lat,
    
    #pca axis 
    all_traits_pc1,  all_traits_pc2,
    chem_traits_pc1, chem_traits_pc2,
    morph_traits_pc1,  morph_traits_pc2,
    
    # Anomalies
    gpp_anomaly_country,
    nee_anomaly_country,
    reco_anomaly_country,
    temperature_gpp_anomaly_country,
    temperature_nee_anomaly_country,
    temperature_reco_anomaly_country,
    mat_anomaly_country,
    elevation_anomaly_country,
    
    leaf_area_anomaly_country, 
    sla_anomaly_country, 
    height_x_cover_anomaly_country, 
    all_traits_pc1_anomaly_country,  all_traits_pc2_anomaly_country,
    chem_traits_pc1_anomaly_country, chem_traits_pc2_anomaly_country,
    morph_traits_pc1_anomaly_country,  morph_traits_pc2_anomaly_country,
    
    # functional_diversity_q1_anomaly_country,
    # species_richness_anomaly_country,
    # # n_percent_anomaly_country,
    # # p_percent_anomaly_country,
    # # cn_ratio_anomaly_country,
    # # c_percent_anomaly_country,
    # par_anomaly_country,
    # soil_moisture_anomaly_country,
    # woodiness_anomaly_country,
    # grassiness_anomaly_country
  ) #%>% filter(complete.cases(.))

dt <- dt_raw %>% 
  dplyr::select(-c(nee, reco, gpp, lat)) %>% 
  mutate(leaf_area_cm2 = log(leaf_area_cm2),
         species_richness = log(species_richness),
         functional_diversity_q1 = log(functional_diversity_q1),
         height_x_cover = log(height_x_cover)
         
  ) %>% 
  mutate(across(where(is.numeric), ~as.numeric(scale(.x)))) %>% 
  left_join(dt_raw[, c("plot_id", "nee", "reco", "gpp", "lat")]) %>% 
  distinct()

# 1. NEE ------------------------------------------
## Moprhological Trait PCA --------------------
m_nee_1 <- psem(
  
  # model for temperature
  glmmTMB(temperature_nee_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for sla / ldmc axis 
  glmmTMB(morph_traits_pc1_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for lead area / dry mass axis 
  glmmTMB(morph_traits_pc2_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  #model for flux
  glmmTMB(nee_anomaly_country ~
            temperature_nee_anomaly_country +
            height_x_cover_anomaly_country +
            morph_traits_pc1_anomaly_country + 
            morph_traits_pc2_anomaly_country +
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  # Correlated errors
  height_x_cover_anomaly_country %~~% morph_traits_pc1_anomaly_country,
  height_x_cover_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  morph_traits_pc1_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  temperature_nee_anomaly_country %~~% height_x_cover_anomaly_country,
  temperature_nee_anomaly_country %~~% morph_traits_pc1_anomaly_country,
  temperature_nee_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  
  data = dt
)
ss_nee_1 <- summary(m_nee_1)
ss_nee_1
dSep(m_nee_1)
LLchisq(m_nee_1)
AIC(m_nee_1)
plot(m_nee_1)
anova(m_nee_1)


## Traits separately --------------------
m_nee_2 <- psem(
  
  # model for temperature
  glmmTMB(temperature_nee_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for sla 
  glmmTMB(sla_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for lead area
  glmmTMB(leaf_area_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  #model for flux
  glmmTMB(nee_anomaly_country ~
            temperature_nee_anomaly_country +
            height_x_cover_anomaly_country +
            sla_anomaly_country + 
            leaf_area_anomaly_country +
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  # Correlated errors
  height_x_cover_anomaly_country %~~% sla_anomaly_country,
  height_x_cover_anomaly_country %~~% leaf_area_anomaly_country,
  sla_anomaly_country %~~% leaf_area_anomaly_country,
  temperature_nee_anomaly_country %~~% height_x_cover_anomaly_country,
  temperature_nee_anomaly_country %~~% sla_anomaly_country,
  temperature_nee_anomaly_country %~~% leaf_area_anomaly_country,
  
  data = dt
)
ss_nee_2 <- summary(m_nee_2)
ss_nee_2
dSep(m_nee_2)
LLchisq(m_nee_2)
AIC(m_nee_2)
plot(m_nee_2)
anova(m_nee_2)

## All Traits PCA --------------------
dt_at <- dt %>% filter(!is.na(all_traits_pc1_anomaly_country))

m_nee_3 <- psem(
  
  # model for temperature
  glmmTMB(temperature_nee_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_at),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_at),
  
  # model for sla / ldmc axis 
  glmmTMB(all_traits_pc1_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_at),
  
  # model for lead area / dry mass axis 
  glmmTMB(all_traits_pc2_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_at),
  
  #model for flux
  glmmTMB(nee_anomaly_country ~
            temperature_nee_anomaly_country +
            height_x_cover_anomaly_country +
            all_traits_pc1_anomaly_country + 
            all_traits_pc2_anomaly_country +
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_at),
  
  # Correlated errors
  height_x_cover_anomaly_country %~~% all_traits_pc1_anomaly_country,
  height_x_cover_anomaly_country %~~% all_traits_pc2_anomaly_country,
  all_traits_pc1_anomaly_country %~~% all_traits_pc2_anomaly_country,
  temperature_nee_anomaly_country %~~% height_x_cover_anomaly_country,
  temperature_nee_anomaly_country %~~% all_traits_pc1_anomaly_country,
  temperature_nee_anomaly_country %~~% all_traits_pc2_anomaly_country,
  
  data = dt_at 
)
ss_nee_3 <- summary(m_nee_3)
ss_nee_3
dSep(m_nee_3)
LLchisq(m_nee_3)
AIC(m_nee_3)
plot(m_nee_3)
anova(m_nee_3)

plot(dt$nee_anomaly_country, dt$mat_anomaly_country)
plot(dt_at$nee_anomaly_country, dt_at$mat_anomaly_country)

## chem Traits PCA --------------------
dt_ct <- dt %>% filter(!is.na(chem_traits_pc1_anomaly_country))

m_nee_4 <- psem(
  
  # model for temperature
  glmmTMB(temperature_nee_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_ct),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_ct),
  
  # model for sla / ldmc axis 
  glmmTMB(chem_traits_pc1_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_ct),
  
  # model for lead area / dry mass axis 
  glmmTMB(chem_traits_pc2_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_ct),
  
  #model for flux
  glmmTMB(nee_anomaly_country ~
            temperature_nee_anomaly_country +
            height_x_cover_anomaly_country +
            chem_traits_pc1_anomaly_country + 
            chem_traits_pc2_anomaly_country +
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_ct),
  
  # Correlated errors
  height_x_cover_anomaly_country %~~% chem_traits_pc1_anomaly_country,
  height_x_cover_anomaly_country %~~% chem_traits_pc2_anomaly_country,
  chem_traits_pc1_anomaly_country %~~% chem_traits_pc2_anomaly_country,
  temperature_nee_anomaly_country %~~% height_x_cover_anomaly_country,
  temperature_nee_anomaly_country %~~% chem_traits_pc1_anomaly_country,
  temperature_nee_anomaly_country %~~% chem_traits_pc2_anomaly_country,
  
  data = dt_ct 
)
ss_nee_4 <- summary(m_nee_4)
ss_nee_4
dSep(m_nee_4)
LLchisq(m_nee_4)
AIC(m_nee_4)
plot(m_nee_4)
anova(m_nee_4)



# 2. Reco ------------------------------------------
## Moprhological Trait PCA --------------------
m_reco_1 <- psem(
  
  # model for temperature
  glmmTMB(temperature_reco_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for sla / ldmc axis 
  glmmTMB(morph_traits_pc1_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for lead area / dry mass axis 
  glmmTMB(morph_traits_pc2_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  #model for flux
  glmmTMB(reco_anomaly_country ~
            temperature_reco_anomaly_country +
            height_x_cover_anomaly_country +
            morph_traits_pc1_anomaly_country + 
            morph_traits_pc2_anomaly_country +
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  # Correlated errors
  height_x_cover_anomaly_country %~~% morph_traits_pc1_anomaly_country,
  height_x_cover_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  morph_traits_pc1_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  temperature_reco_anomaly_country %~~% height_x_cover_anomaly_country,
  temperature_reco_anomaly_country %~~% morph_traits_pc1_anomaly_country,
  temperature_reco_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  
  data = dt
)
ss_reco_1 <- summary(m_reco_1)
ss_reco_1
dSep(m_reco_1)
LLchisq(m_reco_1)
AIC(m_reco_1)
plot(m_reco_1)
anova(m_reco_1)


## Traits separately --------------------
m_reco_2 <- psem(
  
  # model for temperature
  glmmTMB(temperature_reco_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for sla 
  glmmTMB(sla_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for lead area
  glmmTMB(leaf_area_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  #model for flux
  glmmTMB(reco_anomaly_country ~
            temperature_reco_anomaly_country +
            height_x_cover_anomaly_country +
            sla_anomaly_country + 
            leaf_area_anomaly_country +
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  # Correlated errors
  height_x_cover_anomaly_country %~~% sla_anomaly_country,
  height_x_cover_anomaly_country %~~% leaf_area_anomaly_country,
  sla_anomaly_country %~~% leaf_area_anomaly_country,
  temperature_reco_anomaly_country %~~% height_x_cover_anomaly_country,
  temperature_reco_anomaly_country %~~% sla_anomaly_country,
  temperature_reco_anomaly_country %~~% leaf_area_anomaly_country,
  
  data = dt
)
ss_reco_2 <- summary(m_reco_2)
ss_reco_2
dSep(m_reco_2)
LLchisq(m_reco_2)
AIC(m_reco_2)
plot(m_reco_2)
anova(m_reco_2)

## All Traits PCA --------------------
dt_at <- dt %>% filter(!is.na(all_traits_pc1_anomaly_country))

m_reco_3 <- psem(
  
  # model for temperature
  glmmTMB(temperature_reco_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_at),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_at),
  
  # model for sla / ldmc axis 
  glmmTMB(all_traits_pc1_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_at),
  
  # model for lead area / dry mass axis 
  glmmTMB(all_traits_pc2_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_at),
  
  #model for flux
  glmmTMB(reco_anomaly_country ~
            temperature_reco_anomaly_country +
            height_x_cover_anomaly_country +
            all_traits_pc1_anomaly_country + 
            all_traits_pc2_anomaly_country +
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_at),
  
  # Correlated errors
  height_x_cover_anomaly_country %~~% all_traits_pc1_anomaly_country,
  height_x_cover_anomaly_country %~~% all_traits_pc2_anomaly_country,
  all_traits_pc1_anomaly_country %~~% all_traits_pc2_anomaly_country,
  temperature_reco_anomaly_country %~~% height_x_cover_anomaly_country,
  temperature_reco_anomaly_country %~~% all_traits_pc1_anomaly_country,
  temperature_reco_anomaly_country %~~% all_traits_pc2_anomaly_country,
  
  data = dt_at 
)
ss_reco_3 <- summary(m_reco_3)
ss_reco_3
dSep(m_reco_3)
LLchisq(m_reco_3)
AIC(m_reco_3)
plot(m_reco_3)
anova(m_reco_3)

plot(dt$reco_anomaly_country, dt$mat_anomaly_country)
plot(dt_at$reco_anomaly_country, dt_at$mat_anomaly_country)

## chem Traits PCA --------------------
dt_ct <- dt %>% filter(!is.na(chem_traits_pc1_anomaly_country))

m_reco_4 <- psem(
  
  # model for temperature
  glmmTMB(temperature_reco_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_ct),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_ct),
  
  # model for sla / ldmc axis 
  glmmTMB(chem_traits_pc1_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_ct),
  
  # model for lead area / dry mass axis 
  glmmTMB(chem_traits_pc2_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_ct),
  
  #model for flux
  glmmTMB(reco_anomaly_country ~
            temperature_reco_anomaly_country +
            height_x_cover_anomaly_country +
            chem_traits_pc1_anomaly_country + 
            chem_traits_pc2_anomaly_country +
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_ct),
  
  # Correlated errors
  height_x_cover_anomaly_country %~~% chem_traits_pc1_anomaly_country,
  height_x_cover_anomaly_country %~~% chem_traits_pc2_anomaly_country,
  chem_traits_pc1_anomaly_country %~~% chem_traits_pc2_anomaly_country,
  temperature_reco_anomaly_country %~~% height_x_cover_anomaly_country,
  temperature_reco_anomaly_country %~~% chem_traits_pc1_anomaly_country,
  temperature_reco_anomaly_country %~~% chem_traits_pc2_anomaly_country,
  
  data = dt_ct 
)
ss_reco_4 <- summary(m_reco_4)
ss_reco_4
dSep(m_reco_4)
LLchisq(m_reco_4)
AIC(m_reco_4)
plot(m_reco_4)
anova(m_reco_4)


# 3. GPP ------------------------------------------
## Moprhological Trait PCA --------------------
m_gpp_1 <- psem(
  
  # model for temperature
  glmmTMB(temperature_gpp_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for sla / ldmc axis 
  glmmTMB(morph_traits_pc1_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for lead area / dry mass axis 
  glmmTMB(morph_traits_pc2_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  #model for flux
  glmmTMB(gpp_anomaly_country ~
            temperature_gpp_anomaly_country +
            height_x_cover_anomaly_country +
            morph_traits_pc1_anomaly_country + 
            morph_traits_pc2_anomaly_country +
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  # Correlated errors
  height_x_cover_anomaly_country %~~% morph_traits_pc1_anomaly_country,
  height_x_cover_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  morph_traits_pc1_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  temperature_gpp_anomaly_country %~~% height_x_cover_anomaly_country,
  temperature_gpp_anomaly_country %~~% morph_traits_pc1_anomaly_country,
  temperature_gpp_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  
  data = dt
)
ss_gpp_1 <- summary(m_gpp_1)
ss_gpp_1
dSep(m_gpp_1)
LLchisq(m_gpp_1)
AIC(m_gpp_1)
plot(m_gpp_1)
anova(m_gpp_1)


## Traits separately --------------------
m_gpp_2 <- psem(
  
  # model for temperature
  glmmTMB(temperature_gpp_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for sla 
  glmmTMB(sla_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for lead area
  glmmTMB(leaf_area_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  #model for flux
  glmmTMB(gpp_anomaly_country ~
            temperature_gpp_anomaly_country +
            height_x_cover_anomaly_country +
            sla_anomaly_country + 
            leaf_area_anomaly_country +
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  # Correlated errors
  height_x_cover_anomaly_country %~~% sla_anomaly_country,
  height_x_cover_anomaly_country %~~% leaf_area_anomaly_country,
  sla_anomaly_country %~~% leaf_area_anomaly_country,
  temperature_gpp_anomaly_country %~~% height_x_cover_anomaly_country,
  temperature_gpp_anomaly_country %~~% sla_anomaly_country,
  temperature_gpp_anomaly_country %~~% leaf_area_anomaly_country,
  
  data = dt
)
ss_gpp_2 <- summary(m_gpp_2)
ss_gpp_2
dSep(m_gpp_2)
LLchisq(m_gpp_2)
AIC(m_gpp_2)
plot(m_gpp_2)
anova(m_gpp_2)

## All Traits PCA --------------------
dt_at <- dt %>% filter(!is.na(all_traits_pc1_anomaly_country))

m_gpp_3 <- psem(
  
  # model for temperature
  glmmTMB(temperature_gpp_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_at),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_at),
  
  # model for sla / ldmc axis 
  glmmTMB(all_traits_pc1_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_at),
  
  # model for lead area / dry mass axis 
  glmmTMB(all_traits_pc2_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_at),
  
  #model for flux
  glmmTMB(gpp_anomaly_country ~
            temperature_gpp_anomaly_country +
            height_x_cover_anomaly_country +
            all_traits_pc1_anomaly_country + 
            all_traits_pc2_anomaly_country +
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_at),
  
  # Correlated errors
  height_x_cover_anomaly_country %~~% all_traits_pc1_anomaly_country,
  height_x_cover_anomaly_country %~~% all_traits_pc2_anomaly_country,
  all_traits_pc1_anomaly_country %~~% all_traits_pc2_anomaly_country,
  temperature_gpp_anomaly_country %~~% height_x_cover_anomaly_country,
  temperature_gpp_anomaly_country %~~% all_traits_pc1_anomaly_country,
  temperature_gpp_anomaly_country %~~% all_traits_pc2_anomaly_country,
  
  data = dt_at 
)
ss_gpp_3 <- summary(m_gpp_3)
ss_gpp_3
dSep(m_gpp_3)
LLchisq(m_gpp_3)
AIC(m_gpp_3)
plot(m_gpp_3)
anova(m_gpp_3)

plot(dt$gpp_anomaly_country, dt$mat_anomaly_country)
plot(dt_at$gpp_anomaly_country, dt_at$mat_anomaly_country)

## chem Traits PCA --------------------
dt_ct <- dt %>% filter(!is.na(chem_traits_pc1_anomaly_country))

m_gpp_4 <- psem(
  
  # model for temperature
  glmmTMB(temperature_gpp_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_ct),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_ct),
  
  # model for sla / ldmc axis 
  glmmTMB(chem_traits_pc1_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_ct),
  
  # model for lead area / dry mass axis 
  glmmTMB(chem_traits_pc2_anomaly_country ~ 
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_ct),
  
  #model for flux
  glmmTMB(gpp_anomaly_country ~
            temperature_gpp_anomaly_country +
            height_x_cover_anomaly_country +
            chem_traits_pc1_anomaly_country + 
            chem_traits_pc2_anomaly_country +
            mat_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_ct),
  
  # Correlated errors
  height_x_cover_anomaly_country %~~% chem_traits_pc1_anomaly_country,
  height_x_cover_anomaly_country %~~% chem_traits_pc2_anomaly_country,
  chem_traits_pc1_anomaly_country %~~% chem_traits_pc2_anomaly_country,
  temperature_gpp_anomaly_country %~~% height_x_cover_anomaly_country,
  temperature_gpp_anomaly_country %~~% chem_traits_pc1_anomaly_country,
  temperature_gpp_anomaly_country %~~% chem_traits_pc2_anomaly_country,
  
  data = dt_ct 
)
ss_gpp_4 <- summary(m_gpp_4)
ss_gpp_4
dSep(m_gpp_4)
LLchisq(m_gpp_4)
AIC(m_gpp_4)
plot(m_gpp_4)
anova(m_gpp_4)



## Extract model estimates ---------------------------------------

model_list <- list(
  m_nee_1 = m_nee_1,
  m_nee_2 = m_nee_2,
  m_nee_3 = m_nee_3,
  m_nee_4 = m_nee_4,
  m_reco_1 = m_reco_1,
  m_reco_2 = m_reco_2,
  m_reco_3 = m_reco_3,
  m_reco_4 = m_reco_4,
  m_gpp_1 = m_gpp_1,
  m_gpp_2 = m_gpp_2,
  m_gpp_3 = m_gpp_3,
  m_gpp_4 = m_gpp_4
)

dt_res <- data.frame()

for(i in 1:length(model_list)){
  
  
  
  m <- model_list[[i]]
  
  m_name <- names(model_list)[i]
  
  m_s <- summary(m)
  
  tmp_covs <- m_s$coefficients %>% 
    as.data.table()
  names(tmp_covs)[9] <- "sig_symbol"  
  
  tmp_r2 <- m_s$R2
  
  tmp_dt <- tmp_covs %>% 
    left_join(tmp_r2) %>% 
    filter(!grepl("~~", Response)) %>% 
    mutate(model_name = m_name, 
           aic = as.numeric(AIC_psem(m)[1]),
           aicc = as.numeric(AIC_psem(m)[2]))
  
  dt_res <- rbind(dt_res, tmp_dt)
  
}


dt_aic <- dt_res %>% 
  mutate(m_group = ifelse(grepl("1", model_name) | grepl("2", model_name), "full", "nsa"),
         response_tier = case_when(
           .default = "veg",
           grepl("reco", Response) & !grepl("temperature", Response) ~ "flux",
           grepl("nee", Response) & !grepl("temperature", Response) ~ "flux",
           grepl("gpp", Response) & !grepl("temperature", Response) ~ "flux")) %>% 
  filter(response_tier == "flux") %>% 
  arrange(desc(Response), aicc) %>% 
  dplyr::select(model_name, Response, aic, aicc) %>% unique()
dt_aic

dt_est <- dt_res %>% 
  rename(p_value = P.Value, 
         std_error = Std.Error, 
         rsq_c = Conditional, 
         rsq_m = Marginal, 
         predictor = Predictor, 
         response = Response, 
         estimate = Estimate) %>% 
  mutate(
    std_error = as.numeric(std_error),
    predictor_tier = case_when(
    grepl("1", model_name) ~ "morph_traits_pca",
    grepl("2", model_name) ~ "traits_separately",
    grepl("3", model_name) ~ "all_traits_pca", 
    grepl("4", model_name) ~ "chem_traits_pca"),
    response_tier = case_when(
      .default = "veg",
      grepl("reco", response) & !grepl("temperature", response) ~ "flux",
      grepl("nee", response) & !grepl("temperature", response) ~ "flux",
      grepl("gpp", response) & !grepl("temperature", response) ~ "flux"),
    rsq_label = paste0("(R2m = ", round(rsq_m, 2), "; R2c = ", round(rsq_c,2), ")"), 
    clean_response = case_when(
      response == "nee_anomaly_country" ~ paste0("NEE\n", rsq_label),
      response == "gpp_anomaly_country" ~ paste0("GPP\n", rsq_label),
      response == "reco_anomaly_country" ~ paste0("Reco\n", rsq_label),
      response == "height_x_cover_anomaly_country" ~ paste0("'Biomass'\n", rsq_label), 
      response == "morph_traits_pc1_anomaly_country" ~ paste0("Morph. Traits PCA 1\n", rsq_label),
      response == "morph_traits_pc2_anomaly_country" ~ paste0("Morph. Traits PCA 2\n", rsq_label),
      response == "all_traits_pc1_anomaly_country" ~ paste0("All Traits PCA 1\n", rsq_label),
      response == "all_traits_pc2_anomaly_country" ~ paste0("All Traits PCA 2\n", rsq_label),
      response == "chem_traits_pc1_anomaly_country" ~ paste0("Chem. Traits PCA 1\n", rsq_label),
      response == "chem_traits_pc2_anomaly_country" ~ paste0("Chem. Traits PCA 2\n", rsq_label),
      grepl("temperature", response) ~ paste0("Local Temperature\n", rsq_label),
      response == "sla_anomaly_country" ~ paste0("SLA\n", rsq_label),
      response == "leaf_area_anomaly_country" ~ paste0("Leaf Area\n", rsq_label),
    ), 
    ci_lb = estimate - 1.96*std_error,
    ci_ub = estimate + 1.96*std_error, 
    clean_term =  case_when(
      .default = predictor,
      grepl("temperature_", predictor) ~ "Local Temperature",
      predictor == "mat_anomaly_country" ~ "MAT",
      predictor == "height_x_cover_anomaly_country"      ~ "'Biomass'",
      predictor == "morph_traits_pc1_anomaly_country"          ~ "Morph. Traits PC1",
      predictor == "morph_traits_pc2_anomaly_country"          ~ "Morph. Traits PC2",
      predictor == "all_traits_pc1_anomaly_country"          ~ "All Traits PC1",
      predictor == "all_traits_pc2_anomaly_country"          ~ "All Traits PC2",
      predictor == "chem_traits_pc1_anomaly_country"          ~ "Chem. Traits PC1",
      predictor == "chem_traits_pc2_anomaly_country"          ~ "Chem. Traits PC2",
      predictor == "sla_anomaly_country"                 ~ "SLA",
      predictor == "leaf_area_anomaly_country"           ~ "Leaf Area"), 
    significance = case_when(
      .default = "Non significant", 
      ci_lb > 0 ~ "Significantly positive", 
      ci_ub < 0 ~ "Significantly negative", 
    ))


# Plots -------------
theme_est <-   theme(legend.position = "none", 
                     legend.box="vertical",
                     plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                     panel.grid = element_line(color = "seashell"), 
                     #axis.title.x = element_blank(), 
                     axis.text = element_text(size = 12), 
                     axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
                     panel.border = element_rect(color = NA), 
                     panel.background = element_rect(fill = "snow2"), 
                     strip.text.x = element_text(size = 12), 
                     strip.text.y = element_text(size = 12, face = "bold"), 
                     strip.background = element_rect(fill = "seashell", color = "seashell") )

# traits separately -------
p_f_ts <- dt_est %>% 
  filter(!grepl("ntercept", predictor)) %>% 
  filter(predictor_tier == "traits_separately" & response_tier == "flux") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dotted", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = significance),
                  linewidth = 1.2, size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non significant" = "grey", "Significantly positive" = "#fab255","Significantly negative" = "#0f7ba2")) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Estimate", y = "",) +
  theme_est
p_f_ts

p_v_ts <- dt_est %>% 
  filter(!grepl("ntercept", predictor) & grepl("gpp", model_name)) %>% 
  filter(predictor_tier == "traits_separately" & response_tier == "veg") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dotted", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = significance),
                  linewidth = 1.2, size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non significant" = "grey", "Significantly positive" = "#fab255","Significantly negative" = "#0f7ba2")) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Estimate", y = "",) +
  theme_est
p_v_ts


p_ts <- grid.arrange(p_v_ts, p_f_ts, heights = c(1, 2))
ggsave(plot = p_ts,
       "builds/plots/psem_estimates_traits_separately.png",
       dpi = 600, 
       height = 5, width = 10)

#morph. traits PCA -------
p_f_mtpca <- dt_est %>% 
  filter(!grepl("ntercept", predictor)) %>% 
  filter(predictor_tier == "morph_traits_pca" & response_tier == "flux") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dotted", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = significance),
                  linewidth = 1.2, size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non significant" = "grey", "Significantly positive" = "#fab255","Significantly negative" = "#0f7ba2")) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Estimate", y = "",) +
  theme_est
p_f_mtpca

p_v_mtpca <- dt_est %>% 
  filter(!grepl("ntercept", predictor) & grepl("gpp", model_name)) %>% 
  filter(predictor_tier == "morph_traits_pca" & response_tier == "veg") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dotted", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = significance),
                  linewidth = 1.2, size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non significant" = "grey", "Significantly positive" = "#fab255","Significantly negative" = "#0f7ba2")) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Estimate", y = "",) +
  theme_est
p_v_mtpca


p_mtpca <- grid.arrange(p_v_mtpca,p_f_mtpca, heights = c(1, 2))
ggsave(plot = p_mtpca,
       "builds/plots/psem_estimates_morph_traits_pca.png",
       dpi = 600, 
       height = 5, width = 10)

#All traits PCA 
p_f_atpca <- dt_est %>% 
  filter(!grepl("ntercept", predictor)) %>% 
  filter(predictor_tier == "all_traits_pca" & response_tier == "flux") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dotted", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = significance),
                  linewidth = 1.2, size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non significant" = "grey", "Significantly positive" = "#fab255","Significantly negative" = "#0f7ba2")) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Estimate", y = "",) +
  theme_est
p_f_atpca

p_v_atpca <- dt_est %>% 
  filter(!grepl("ntercept", predictor) & grepl("gpp", model_name)) %>% 
  filter(predictor_tier == "all_traits_pca" & response_tier == "veg") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dotted", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = significance),
                  linewidth = 1.2, size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non significant" = "grey", "Significantly positive" = "#fab255","Significantly negative" = "#0f7ba2")) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Estimate", y = "",) +
  theme_est
p_v_atpca


p_atpca <- grid.arrange(p_v_atpca,p_f_atpca, heights = c(1, 2))
ggsave(plot = p_atpca,
       "builds/plots/psem_estimates_all_traits_pca.png",
       dpi = 600, 
       height = 5, width = 10)

# Chemical traits PCA
p_f_ctpca <- dt_est %>% 
  filter(!grepl("ntercept", predictor)) %>% 
  filter(predictor_tier == "chem_traits_pca" & response_tier == "flux") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dotted", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = significance),
                  linewidth = 1.2, size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non significant" = "grey", "Significantly positive" = "#fab255","Significantly negative" = "#0f7ba2")) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Estimate", y = "",) +
  theme_est
p_f_ctpca

p_v_ctpca <- dt_est %>% 
  filter(!grepl("ntercept", predictor) & grepl("gpp", model_name)) %>% 
  filter(predictor_tier == "chem_traits_pca" & response_tier == "veg") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dotted", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = significance),
                  linewidth = 1.2, size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non significant" = "grey", "Significantly positive" = "#fab255","Significantly negative" = "#0f7ba2")) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Estimate", y = "",) +
  theme_est
p_v_ctpca


p_ctpca <- grid.arrange(p_v_ctpca,p_f_ctpca, heights = c(1, 2))
ggsave(plot = p_ctpca,
       "builds/plots/psem_estimates_chem_traits_pca.png",
       dpi = 600, 
       height = 5, width = 10)



