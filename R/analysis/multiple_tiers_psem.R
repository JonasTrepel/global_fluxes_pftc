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
    sla_cm2_g, ldmc, leaf_area_cm2, dry_mass_t1_g, plant_height_cm,
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
    elevation_anomaly_country,
    elevation_anomaly_country,
    
    leaf_area_anomaly_country, 
    sla_anomaly_country, 
    height_x_cover_anomaly_country, 
    all_traits_pc1_anomaly_country,  all_traits_pc2_anomaly_country,
    chem_traits_pc1_anomaly_country, chem_traits_pc2_anomaly_country,
    morph_traits_pc1_anomaly_country,  morph_traits_pc2_anomaly_country,
    
    # functional_diversity_q1_anomaly_country,
    # species_richness_t1_anomaly_country,
    # # n_percent_anomaly_country,
    # # p_percent_anomaly_country,
    # # cn_ratio_anomaly_country,
    # # c_percent_anomaly_country,
    # par_anomaly_country,
    # soil_moisture_anomaly_country,
    # woodiness_t1_anomaly_country,
    # grassiness_t1_anomaly_country
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


#########################################################################
###########################    FIRST TIER     ###########################
#########################################################################

## Moprhological Trait PCA --------------------
m_t1_nee <- psem(
  
  #model for flux
  glmmTMB(nee_anomaly_country ~
            elevation_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  data = dt
)
ss_t1_nee <- summary(m_t1_nee)
ss_t1_nee
dSep(m_t1_nee)
LLchisq(m_t1_nee)
AIC(m_t1_nee)
plot(m_t1_nee)
anova(m_t1_nee)


# 2. Reco ------------------------------------------

m_t1_reco <- psem(
  
  #model for flux
  glmmTMB(reco_anomaly_country ~
            elevation_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  data = dt
)
ss_t1_reco <- summary(m_t1_reco)
ss_t1_reco
dSep(m_t1_reco)
LLchisq(m_t1_reco)
AIC(m_t1_reco)
plot(m_t1_reco)
anova(m_t1_reco)

# 3. GPP ------------------------------------------
m_t1_gpp <- psem(

  
  #model for flux
  glmmTMB(gpp_anomaly_country ~
            elevation_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  data = dt
)
ss_t1_gpp <- summary(m_t1_gpp)
ss_t1_gpp
dSep(m_t1_gpp)
LLchisq(m_t1_gpp)
AIC(m_t1_gpp)
plot(m_t1_gpp)
anova(m_t1_gpp)


#########################################################################
###########################    SECOND TIER     ###########################
#########################################################################

# 1. NEE ------------------------------------------

## Moprhological Trait PCA --------------------
m_t2_nee <- psem(
  
  # model for temperature
  glmmTMB(temperature_nee_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  #model for flux
  glmmTMB(nee_anomaly_country ~
            temperature_nee_anomaly_country +
            elevation_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  # Correlated errors
  # height_x_cover_anomaly_country %~~% morph_traits_pc1_anomaly_country,
  # height_x_cover_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  # morph_traits_pc1_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  # temperature_nee_anomaly_country %~~% height_x_cover_anomaly_country,
  # temperature_nee_anomaly_country %~~% morph_traits_pc1_anomaly_country,
  # temperature_nee_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  
  data = dt
)
ss_t2_nee <- summary(m_t2_nee)
ss_t2_nee
dSep(m_t2_nee)
LLchisq(m_t2_nee)
AIC(m_t2_nee)
plot(m_t2_nee)
anova(m_t2_nee)


# 2. Reco ------------------------------------------

m_t2_reco <- psem(
  
  # model for temperature
  glmmTMB(temperature_reco_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),

  #model for flux
  glmmTMB(reco_anomaly_country ~
            temperature_reco_anomaly_country +
            elevation_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  # Correlated errors
  # height_x_cover_anomaly_country %~~% morph_traits_pc1_anomaly_country,
  # height_x_cover_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  # morph_traits_pc1_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  # temperature_reco_anomaly_country %~~% height_x_cover_anomaly_country,
  # temperature_reco_anomaly_country %~~% morph_traits_pc1_anomaly_country,
  # temperature_reco_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  
  data = dt
)
ss_t2_reco <- summary(m_t2_reco)
ss_t2_reco
dSep(m_t2_reco)
LLchisq(m_t2_reco)
AIC(m_t2_reco)
plot(m_t2_reco)
anova(m_t2_reco)

# 3. GPP ------------------------------------------
m_t2_gpp <- psem(
  
  # model for temperature
  glmmTMB(temperature_gpp_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),

  #model for flux
  glmmTMB(gpp_anomaly_country ~
            temperature_gpp_anomaly_country +
            elevation_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  # Correlated errors

  data = dt
)
ss_t2_gpp <- summary(m_t2_gpp)
ss_t2_gpp
dSep(m_t2_gpp)
LLchisq(m_t2_gpp)
AIC(m_t2_gpp)
plot(m_t2_gpp)
anova(m_t2_gpp)

#########################################################################
###########################    THIRD TIER     ###########################
#########################################################################

# 1. NEE ------------------------------------------

## Moprhological Trait PCA --------------------
m_t3_nee <- psem(
  
  # model for temperature
  glmmTMB(temperature_nee_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for biomass
  glmmTMB(height_x_cover_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  
  #model for flux
  glmmTMB(nee_anomaly_country ~
            temperature_nee_anomaly_country +
            elevation_anomaly_country +
            height_x_cover_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  # Correlated errors
    height_x_cover_anomaly_country %~~% temperature_nee_anomaly_country,
  # height_x_cover_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  # morph_traits_pc1_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  # temperature_nee_anomaly_country %~~% height_x_cover_anomaly_country,
  # temperature_nee_anomaly_country %~~% morph_traits_pc1_anomaly_country,
  # temperature_nee_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  
  data = dt
)
ss_t3_nee <- summary(m_t3_nee)
ss_t3_nee
dSep(m_t3_nee)
LLchisq(m_t3_nee)
AIC(m_t3_nee)
plot(m_t3_nee)
anova(m_t3_nee)


# 2. Reco ------------------------------------------

m_t3_reco <- psem(
  
  # model for temperature
  glmmTMB(temperature_reco_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for biomass
  glmmTMB(height_x_cover_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  
  #model for flux
  glmmTMB(reco_anomaly_country ~
            temperature_reco_anomaly_country +
            elevation_anomaly_country +
            height_x_cover_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  # Correlated errors
  height_x_cover_anomaly_country %~~% temperature_reco_anomaly_country,
  # height_x_cover_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  # morph_traits_pc1_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  # temperature_reco_anomaly_country %~~% height_x_cover_anomaly_country,
  # temperature_reco_anomaly_country %~~% morph_traits_pc1_anomaly_country,
  # temperature_reco_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  
  data = dt
)
ss_t3_reco <- summary(m_t3_reco)
ss_t3_reco
dSep(m_t3_reco)
LLchisq(m_t3_reco)
AIC(m_t3_reco)
plot(m_t3_reco)
anova(m_t3_reco)

# 3. GPP ------------------------------------------
m_t3_gpp <- psem(
  
  # model for temperature
  glmmTMB(temperature_gpp_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for biomass
  glmmTMB(height_x_cover_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  
  #model for flux
  glmmTMB(gpp_anomaly_country ~
            temperature_gpp_anomaly_country +
            elevation_anomaly_country +
            height_x_cover_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  # Correlated errors
  height_x_cover_anomaly_country %~~% temperature_gpp_anomaly_country,
  # height_x_cover_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  # morph_traits_pc1_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  # temperature_reco_anomaly_country %~~% height_x_cover_anomaly_country,
  # temperature_reco_anomaly_country %~~% morph_traits_pc1_anomaly_country,
  # temperature_reco_anomaly_country %~~% morph_traits_pc2_anomaly_country,
  
  data = dt
)
ss_t3_gpp <- summary(m_t3_gpp)
ss_t3_gpp
dSep(m_t3_gpp)
LLchisq(m_t3_gpp)
AIC(m_t3_gpp)
plot(m_t3_gpp)
anova(m_t3_gpp)


# 1. NEE ------------------------------------------
## Moprhological Trait PCA --------------------
m_t4_nee_1 <- psem(
  
  # model for temperature
  glmmTMB(temperature_nee_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for sla / ldmc axis 
  glmmTMB(morph_traits_pc1_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for lead area / dry mass axis 
  glmmTMB(morph_traits_pc2_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  #model for flux
  glmmTMB(nee_anomaly_country ~
            temperature_nee_anomaly_country +
            height_x_cover_anomaly_country +
            morph_traits_pc1_anomaly_country + 
            morph_traits_pc2_anomaly_country +
            elevation_anomaly_country +
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
ss_t4_nee_1 <- summary(m_t4_nee_1)
ss_t4_nee_1
dSep(m_t4_nee_1)
LLchisq(m_t4_nee_1)
AIC(m_t4_nee_1)
plot(m_t4_nee_1)
anova(m_t4_nee_1)


## Traits separately --------------------
m_t4_nee_2 <- psem(
  
  # model for temperature
  glmmTMB(temperature_nee_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for sla 
  glmmTMB(sla_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for lead area
  glmmTMB(leaf_area_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  #model for flux
  glmmTMB(nee_anomaly_country ~
            temperature_nee_anomaly_country +
            height_x_cover_anomaly_country +
            sla_anomaly_country + 
            leaf_area_anomaly_country +
            elevation_anomaly_country +
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
ss_t4_nee_2 <- summary(m_t4_nee_2)
ss_t4_nee_2
dSep(m_t4_nee_2)
LLchisq(m_t4_nee_2)
AIC(m_t4_nee_2)
plot(m_t4_nee_2)
anova(m_t4_nee_2)

## All Traits PCA --------------------
dt_at <- dt %>% filter(!is.na(all_traits_pc1_anomaly_country))

m_t4_nee_3 <- psem(
  
  # model for temperature
  glmmTMB(temperature_nee_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_at),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_at),
  
  # model for sla / ldmc axis 
  glmmTMB(all_traits_pc1_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_at),
  
  # model for lead area / dry mass axis 
  glmmTMB(all_traits_pc2_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_at),
  
  #model for flux
  glmmTMB(nee_anomaly_country ~
            temperature_nee_anomaly_country +
            height_x_cover_anomaly_country +
            all_traits_pc1_anomaly_country + 
            all_traits_pc2_anomaly_country +
            elevation_anomaly_country +
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
ss_t4_nee_3 <- summary(m_t4_nee_3)
ss_t4_nee_3
dSep(m_t4_nee_3)
LLchisq(m_t4_nee_3)
AIC(m_t4_nee_3)
plot(m_t4_nee_3)
anova(m_t4_nee_3)


# 2. Reco ------------------------------------------
## Moprhological Trait PCA --------------------
m_t4_reco_1 <- psem(
  
  # model for temperature
  glmmTMB(temperature_reco_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for sla / ldmc axis 
  glmmTMB(morph_traits_pc1_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for lead area / dry mass axis 
  glmmTMB(morph_traits_pc2_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  #model for flux
  glmmTMB(reco_anomaly_country ~
            temperature_reco_anomaly_country +
            height_x_cover_anomaly_country +
            morph_traits_pc1_anomaly_country + 
            morph_traits_pc2_anomaly_country +
            elevation_anomaly_country +
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
ss_t4_reco_1 <- summary(m_t4_reco_1)
ss_t4_reco_1
dSep(m_t4_reco_1)
LLchisq(m_t4_reco_1)
AIC(m_t4_reco_1)
plot(m_t4_reco_1)
anova(m_t4_reco_1)


## Traits separately --------------------
m_t4_reco_2 <- psem(
  
  # model for temperature
  glmmTMB(temperature_reco_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for sla 
  glmmTMB(sla_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for lead area
  glmmTMB(leaf_area_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  #model for flux
  glmmTMB(reco_anomaly_country ~
            temperature_reco_anomaly_country +
            height_x_cover_anomaly_country +
            sla_anomaly_country + 
            leaf_area_anomaly_country +
            elevation_anomaly_country +
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
ss_t4_reco_2 <- summary(m_t4_reco_2)
ss_t4_reco_2
dSep(m_t4_reco_2)
LLchisq(m_t4_reco_2)
AIC(m_t4_reco_2)
plot(m_t4_reco_2)
anova(m_t4_reco_2)

## All Traits PCA --------------------
dt_at <- dt %>% filter(!is.na(all_traits_pc1_anomaly_country))

m_t4_reco_3 <- psem(
  
  # model for temperature
  glmmTMB(temperature_reco_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_at),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_at),
  
  # model for sla / ldmc axis 
  glmmTMB(all_traits_pc1_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_at),
  
  # model for lead area / dry mass axis 
  glmmTMB(all_traits_pc2_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_at),
  
  #model for flux
  glmmTMB(reco_anomaly_country ~
            temperature_reco_anomaly_country +
            height_x_cover_anomaly_country +
            all_traits_pc1_anomaly_country + 
            all_traits_pc2_anomaly_country +
            elevation_anomaly_country +
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
ss_t4_reco_3 <- summary(m_t4_reco_3)
ss_t4_reco_3
dSep(m_t4_reco_3)
LLchisq(m_t4_reco_3)
AIC(m_t4_reco_3)
plot(m_t4_reco_3)
anova(m_t4_reco_3)

# 3. GPP ------------------------------------------
## Moprhological Trait PCA --------------------
m_t4_gpp_1 <- psem(
  
  # model for temperature
  glmmTMB(temperature_gpp_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for sla / ldmc axis 
  glmmTMB(morph_traits_pc1_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for lead area / dry mass axis 
  glmmTMB(morph_traits_pc2_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  #model for flux
  glmmTMB(gpp_anomaly_country ~
            temperature_gpp_anomaly_country +
            height_x_cover_anomaly_country +
            morph_traits_pc1_anomaly_country + 
            morph_traits_pc2_anomaly_country +
            elevation_anomaly_country +
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
ss_t4_gpp_1 <- summary(m_t4_gpp_1)
ss_t4_gpp_1
dSep(m_t4_gpp_1)
LLchisq(m_t4_gpp_1)
AIC(m_t4_gpp_1)
plot(m_t4_gpp_1)
anova(m_t4_gpp_1)


## Traits separately --------------------
m_t4_gpp_2 <- psem(
  
  # model for temperature
  glmmTMB(temperature_gpp_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for sla 
  glmmTMB(sla_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt),
  
  # model for lead area
  glmmTMB(leaf_area_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt),
  
  #model for flux
  glmmTMB(gpp_anomaly_country ~
            temperature_gpp_anomaly_country +
            height_x_cover_anomaly_country +
            sla_anomaly_country + 
            leaf_area_anomaly_country +
            elevation_anomaly_country +
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
ss_t4_gpp_2 <- summary(m_t4_gpp_2)
ss_t4_gpp_2
dSep(m_t4_gpp_2)
LLchisq(m_t4_gpp_2)
AIC(m_t4_gpp_2)
plot(m_t4_gpp_2)
anova(m_t4_gpp_2)

## All Traits PCA --------------------
dt_at <- dt %>% filter(!is.na(all_traits_pc1_anomaly_country))

m_t4_gpp_3 <- psem(
  
  # model for temperature
  glmmTMB(temperature_gpp_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_at),
  
  
  # model for veg volume / biomass 
  glmmTMB(height_x_cover_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_at),
  
  # model for sla / ldmc axis 
  glmmTMB(all_traits_pc1_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site),
          na.action = na.omit,
          data = dt_at),
  
  # model for lead area / dry mass axis 
  glmmTMB(all_traits_pc2_anomaly_country ~ 
            elevation_anomaly_country +
            ( 1 | site), 
          na.action = na.omit,
          data = dt_at),
  
  #model for flux
  glmmTMB(gpp_anomaly_country ~
            temperature_gpp_anomaly_country +
            height_x_cover_anomaly_country +
            all_traits_pc1_anomaly_country + 
            all_traits_pc2_anomaly_country +
            elevation_anomaly_country +
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
ss_t4_gpp_3 <- summary(m_t4_gpp_3)
ss_t4_gpp_3
dSep(m_t4_gpp_3)
LLchisq(m_t4_gpp_3)
AIC(m_t4_gpp_3)
plot(m_t4_gpp_3)
anova(m_t4_gpp_3)


## Extract model estimates ---------------------------------------

model_list <- list(
  m_t1_nee = m_t1_nee,
  m_t1_reco = m_t1_reco,
  m_t1_gpp = m_t1_gpp,

  m_t2_nee = m_t2_nee,
  m_t2_reco = m_t2_reco,
  m_t2_gpp = m_t2_gpp,
  
  m_t3_nee = m_t3_nee,
  m_t3_reco = m_t3_reco,
  m_t3_gpp = m_t3_gpp,
  
  m_t4_nee_1 = m_t4_nee_1,
  m_t4_nee_2 = m_t4_nee_2,
  m_t4_nee_3 = m_t4_nee_3,
  m_t4_reco_1 = m_t4_reco_1,
  m_t4_reco_2 = m_t4_reco_2,
  m_t4_reco_3 = m_t4_reco_3,
  m_t4_gpp_1 = m_t4_gpp_1,
  m_t4_gpp_2 = m_t4_gpp_2,
  m_t4_gpp_3 = m_t4_gpp_3
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
  filter(grepl("t4", model_name)) %>%
  mutate(m_group = ifelse(grepl("_1", model_name) | grepl("_2", model_name), "full", "nsa"),
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
    model_tier = case_when(
      grepl("t1", model_name) ~ "tier_1",
      grepl("t2", model_name) ~ "tier_2",
      grepl("t3", model_name) ~ "tier_3",
      grepl("t4", model_name) ~ "tier_4"),
    predictor_tier = case_when(
      grepl("1", model_name) ~ "morph_traits_pca",
      grepl("2", model_name) ~ "traits_separately",
      grepl("3", model_name) ~ "all_traits_pca"),
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
      predictor == "elevation_anomaly_country" ~ "Elevation",
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

# Tier 1 ----------------
p_t1_flux <- dt_est %>% 
  filter(!grepl("ntercept", predictor)) %>% 
  filter(model_tier == "tier_1" & response_tier == "flux") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dotted", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = significance),
                  linewidth = 1.2, size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non significant" = "grey", "Significantly positive" = "#fab255","Significantly negative" = "#0f7ba2")) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Estimate", y = "",) +
  theme_est
p_t1_flux

ggsave(plot = p_t1_flux,
       "builds/plots/psem_t1_estimates_all_traits_pca.png",
       dpi = 600, 
       height = 3, width = 8)


# Tier 2 ----------------------
p_t2_flux <- dt_est %>% 
  filter(!grepl("ntercept", predictor)) %>% 
  filter(model_tier == "tier_2" & response_tier == "flux") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dotted", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = significance),
                  linewidth = 1.2, size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non significant" = "grey", "Significantly positive" = "#fab255","Significantly negative" = "#0f7ba2")) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Estimate", y = "",) +
  theme_est
p_t2_flux

p_t2_pred <- dt_est %>% 
  filter(!grepl("ntercept", predictor)) %>% 
  filter(model_tier == "tier_2" & response_tier == "veg" & grepl("gpp", model_name)) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dotted", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = significance),
                  linewidth = 1.2, size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non significant" = "grey", "Significantly positive" = "#fab255","Significantly negative" = "#0f7ba2")) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Estimate", y = "",) +
  theme_est
p_t2_pred

p_t2 <- grid.arrange(p_t2_pred ,p_t2_flux, heights = c(1, 2))

ggsave(plot = p_t2,
       "builds/plots/psem_t2_estimates_all_traits_pca.png",
       dpi = 600, 
       height = 5, width = 8)


# Tier 3 ------------ 
p_t3_flux <- dt_est %>% 
  filter(!grepl("ntercept", predictor)) %>% 
  filter(model_tier == "tier_3" & response_tier == "flux") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dotted", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = significance),
                  linewidth = 1.2, size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non significant" = "grey", "Significantly positive" = "#fab255","Significantly negative" = "#0f7ba2")) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Estimate", y = "",) +
  theme_est
p_t3_flux

p_t3_pred <- dt_est %>% 
  filter(!grepl("ntercept", predictor)) %>% 
  filter(model_tier == "tier_3" & response_tier == "veg"  & grepl("gpp", model_name)) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dotted", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = significance),
                  linewidth = 1.2, size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non significant" = "grey", "Significantly positive" = "#fab255","Significantly negative" = "#0f7ba2")) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Estimate", y = "",) +
  theme_est
p_t3_pred

p_t3 <- grid.arrange(p_t3_pred ,p_t3_flux, heights = c(1, 2))

ggsave(plot = p_t3,
       "builds/plots/psem_t3_estimates_all_traits_pca.png",
       dpi = 600, 
       height = 5, width = 8)


# Tier 4--------------------------

# traits separately -------
p_t4_flux_ts <- dt_est %>% 
  filter(!grepl("ntercept", predictor)) %>% 
  filter(model_tier == "tier_4" & predictor_tier == "traits_separately" & response_tier == "flux") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dotted", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = significance),
                  linewidth = 1.2, size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non significant" = "grey", "Significantly positive" = "#fab255","Significantly negative" = "#0f7ba2")) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Estimate", y = "",) +
  theme_est
p_t4_flux_ts

p_t4_pred_ts <- dt_est %>% 
  filter(!grepl("ntercept", predictor) & grepl("gpp", model_name)) %>% 
  filter(model_tier == "tier_4" & predictor_tier == "traits_separately" & response_tier == "veg") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dotted", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = significance),
                  linewidth = 1.2, size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non significant" = "grey", "Significantly positive" = "#fab255","Significantly negative" = "#0f7ba2")) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Estimate", y = "",) +
  theme_est
p_t4_pred_ts


p_t4_ts <- grid.arrange(p_t4_pred_ts, p_t4_flux_ts, heights = c(1, 2))

ggsave(plot = p_t4_ts,
       "builds/plots/psem_t4_estimates_traits_separately.png",
       dpi = 600, 
       height = 5, width = 10)

#morph. traits PCA -------
p_t4_flux_mtpca <- dt_est %>% 
  filter(!grepl("ntercept", predictor)) %>% 
  filter(model_tier == "tier_4" & predictor_tier == "morph_traits_pca" & response_tier == "flux") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dotted", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = significance),
                  linewidth = 1.2, size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non significant" = "grey", "Significantly positive" = "#fab255","Significantly negative" = "#0f7ba2")) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Estimate", y = "",) +
  theme_est
p_t4_flux_mtpca

p_t4_pred_mtpca <- dt_est %>% 
  filter(!grepl("ntercept", predictor) & grepl("gpp", model_name)) %>% 
  filter(model_tier == "tier_4" & predictor_tier == "morph_traits_pca" & response_tier == "veg") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dotted", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = significance),
                  linewidth = 1.2, size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non significant" = "grey", "Significantly positive" = "#fab255","Significantly negative" = "#0f7ba2")) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Estimate", y = "",) +
  theme_est
p_t4_pred_mtpca


p_t4_mtpca <- grid.arrange(p_t4_pred_mtpca,p_t4_flux_mtpca, heights = c(1, 2))

ggsave(plot = p_t4_mtpca,
       "builds/plots/psem_t4_estimates_morph_traits_pca.png",
       dpi = 600, 
       height = 5, width = 10)

#All traits PCA 
p_t4_flux_atpca <- dt_est %>% 
  filter(!grepl("ntercept", predictor)) %>% 
  filter(model_tier == "tier_4" & predictor_tier == "all_traits_pca" & response_tier == "flux") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dotted", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = significance),
                  linewidth = 1.2, size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non significant" = "grey", "Significantly positive" = "#fab255","Significantly negative" = "#0f7ba2")) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Estimate", y = "",) +
  theme_est
p_t4_flux_atpca

p_t4_pred_atpca <- dt_est %>% 
  filter(!grepl("ntercept", predictor) & grepl("gpp", model_name)) %>% 
  filter(model_tier == "tier_4" & predictor_tier == "all_traits_pca" & response_tier == "veg") %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dotted", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = significance),
                  linewidth = 1.2, size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non significant" = "grey", "Significantly positive" = "#fab255","Significantly negative" = "#0f7ba2")) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 5)  +
  theme_bw() +
  labs(x = "Estimate", y = "",) +
  theme_est
p_t4_pred_atpca


p_t4_atpca <- grid.arrange(p_t4_pred_atpca,p_t4_flux_atpca, heights = c(1, 2))
ggsave(plot = p_t4_atpca,
       "builds/plots/psem_t4_estimates_all_traits_pca.png",
       dpi = 600, 
       height = 5, width = 10)

