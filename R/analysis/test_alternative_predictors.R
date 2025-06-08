# test alternative hypothesis 


library("lavaan")
library("piecewiseSEM")
library(data.table)
library(tidyverse)
library(lme4)
library(MuMIn)
library(glmmTMB)
library("sjPlot")
library(performance)
library(gridExtra)


dt_raw <- fread("data/processed_data/clean_data/global_fluxes_main_data.csv") %>% 
  dplyr::select(
    # identifiers 
    country, site, plot_id, treatment,
    
    # fluxes
    nee, reco, gpp,
    
    # environmental 
    elevation, map, mat,
    temperature_nee, temperature_reco, temperature_gpp,
    
    # trait means 
    sla_cm2_g, ldmc, leaf_area_cm2, dry_mass_g,
    
    # others 
    height_x_cover, species_richness,
    functional_diversity_q1, lat,
    
    #pca axis 
    all_traits_pc1, all_traits_pc2,
    chem_traits_pc1, chem_traits_pc2,
    morph_traits_pc1, morph_traits_pc2,
    
    # Country-level means
    gpp_country_mean,
    nee_country_mean,
    reco_country_mean, 
    temperature_gpp_country_mean,
    temperature_nee_country_mean,
    temperature_reco_country_mean, 
    mat_country_mean, 
    
    leaf_area_country_mean, 
    sla_country_mean, 
    height_x_cover_country_mean, 
    species_richness_country_mean,
    all_traits_pc1_country_mean, all_traits_pc2_country_mean,
    chem_traits_pc1_country_mean, chem_traits_pc2_country_mean,
    morph_traits_pc1_country_mean, morph_traits_pc2_country_mean,
    
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
    all_traits_pc1_anomaly_country, all_traits_pc2_anomaly_country,
    chem_traits_pc1_anomaly_country, chem_traits_pc2_anomaly_country,
    morph_traits_pc1_anomaly_country, morph_traits_pc2_anomaly_country,
    
    functional_diversity_q1_anomaly_country,
    species_richness_anomaly_country,
    n_percent_anomaly_country,
    p_percent_anomaly_country,
    cn_ratio_anomaly_country,
    c_percent_anomaly_country,
    par_anomaly_country,
    soil_moisture_anomaly_country,
    woodiness_anomaly_country,
    grassiness_anomaly_country
  ) %>% 
  filter(!is.na(sla_cm2_g))

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


### Test response of alternative hypotheses to climate -------
# PAR
m_alt_par <- glmmTMB(par_anomaly_country ~ 
                       mat_anomaly_country +
                       (1 | site), 
                     na.action = na.omit, data = dt)
summary(m_alt_par); check_collinearity(m_alt_par); r.squaredGLMM(m_alt_par)

# Woodiness
m_alt_wood <- glmmTMB(woodiness_anomaly_country ~ 
                        mat_anomaly_country +
                        (1 | site), 
                      na.action = na.omit, data = dt)
summary(m_alt_wood); check_collinearity(m_alt_wood); r.squaredGLMM(m_alt_wood)

# Grassiness
m_alt_grass <- glmmTMB(grassiness_anomaly_country ~ 
                         mat_anomaly_country +
                         (1 | site), 
                       na.action = na.omit, data = dt)
summary(m_alt_grass); check_collinearity(m_alt_grass); r.squaredGLMM(m_alt_grass)


# Species Richness
m_alt_sr <- glmmTMB(species_richness_anomaly_country ~ 
                        mat_anomaly_country +
                        (1 | site), 
                      na.action = na.omit, data = dt)
summary(m_alt_sr); check_collinearity(m_alt_sr); r.squaredGLMM(m_alt_sr)


# Functional Diversity
m_alt_fd <- glmmTMB(functional_diversity_q1_anomaly_country ~ 
                      mat_anomaly_country +
                      (1 | site), 
                    na.action = na.omit, data = dt)
summary(m_alt_fd); check_collinearity(m_alt_fd); r.squaredGLMM(m_alt_fd)



# Soil moisture
m_alt_soil <- glmmTMB(soil_moisture_anomaly_country ~ 
                        mat_anomaly_country +
                        (1 | site), 
                      na.action = na.omit, data = dt)
summary(m_alt_soil); check_collinearity(m_alt_soil); r.squaredGLMM(m_alt_soil)

#NEE ---------
nee_base_formula <- "nee_anomaly_country ~ temperature_nee_anomaly_country +
                  height_x_cover_anomaly_country +
                  sla_anomaly_country + 
                  leaf_area_anomaly_country +
                  mat_anomaly_country"

# PAR
m_nee_no_par <- glmmTMB(
  formula = as.formula(paste(nee_base_formula, "+ (1 | site)")),
  data = dt[!is.na(dt$par_anomaly_country)],
  na.action = na.omit
)
summary(m_nee_no_par); check_collinearity(m_nee_no_par); r.squaredGLMM(m_nee_no_par); AICc(m_nee_no_par)

m_nee_with_par <- glmmTMB(
  formula = as.formula(paste(nee_base_formula, "+ par_anomaly_country + (1 | site)")),
  data = dt[!is.na(dt$par_anomaly_country)],
  na.action = na.omit
)
summary(m_nee_with_par); check_collinearity(m_nee_with_par); r.squaredGLMM(m_nee_with_par); AICc(m_nee_with_par)


# WOODINESS
m_nee_no_wood <- glmmTMB(
  formula = as.formula(paste(nee_base_formula, "+ (1 | site)")),
  data = dt[!is.na(dt$woodiness_anomaly_country)],
  na.action = na.omit
)
summary(m_nee_no_wood); check_collinearity(m_nee_no_wood); r.squaredGLMM(m_nee_no_wood); AICc(m_nee_no_wood)

m_nee_with_wood <- glmmTMB(
  formula = as.formula(paste(nee_base_formula, "+ woodiness_anomaly_country + (1 | site)")),
  data = dt[!is.na(dt$woodiness_anomaly_country)],
  na.action = na.omit
)
summary(m_nee_with_wood); check_collinearity(m_nee_with_wood); r.squaredGLMM(m_nee_with_wood); AICc(m_nee_with_wood)


# GRASSINESS
m_nee_no_grass <- glmmTMB(
  formula = as.formula(paste(nee_base_formula, "+ (1 | site)")),
  data = dt[!is.na(dt$grassiness_anomaly_country)],
  na.action = na.omit
)
summary(m_nee_no_grass); check_collinearity(m_nee_no_grass); r.squaredGLMM(m_nee_no_grass); AICc(m_nee_no_grass)

m_nee_with_grass <- glmmTMB(
  formula = as.formula(paste(nee_base_formula, "+ grassiness_anomaly_country + (1 | site)")),
  data = dt[!is.na(dt$grassiness_anomaly_country)],
  na.action = na.omit
)
summary(m_nee_with_grass); check_collinearity(m_nee_with_grass); r.squaredGLMM(m_nee_with_grass); AICc(m_nee_with_grass)


# SOIL MOISTURE
m_nee_no_soil <- glmmTMB(
  formula = as.formula(paste(nee_base_formula, "+ (1 | site)")),
  data = dt[!is.na(dt$soil_moisture_anomaly_country)],
  na.action = na.omit
)
summary(m_nee_no_soil); check_collinearity(m_nee_no_soil); r.squaredGLMM(m_nee_no_soil); AICc(m_nee_no_soil)

m_nee_with_soil <- glmmTMB(
  formula = as.formula(paste(nee_base_formula, "+ soil_moisture_anomaly_country + (1 | site)")),
  data = dt[!is.na(dt$soil_moisture_anomaly_country)],
  na.action = na.omit
)
summary(m_nee_with_soil); check_collinearity(m_nee_with_soil); r.squaredGLMM(m_nee_with_soil); AICc(m_nee_with_soil)

# SPECIES RICHNESS
m_nee_no_sr <- glmmTMB(
  formula = as.formula(paste(nee_base_formula, "+ (1 | site)")),
  data = dt[!is.na(dt$species_richness_anomaly_country)],
  na.action = na.omit
)
summary(m_nee_no_sr); check_collinearity(m_nee_no_sr); r.squaredGLMM(m_nee_no_sr); AICc(m_nee_no_sr)

m_nee_with_sr <- glmmTMB(
  formula = as.formula(paste(nee_base_formula, "+ species_richness_anomaly_country + (1 | site)")),
  data = dt[!is.na(dt$species_richness_anomaly_country)],
  na.action = na.omit
)
summary(m_nee_with_sr); check_collinearity(m_nee_with_sr); r.squaredGLMM(m_nee_with_sr); AICc(m_nee_with_sr)

# FUNCTIONAL DIVERSITY
  m_nee_no_fd <- glmmTMB(
    formula = as.formula(paste(nee_base_formula, "+ (1 | site)")),
    data = dt[!is.na(dt$functional_diversity_q1_anomaly_country)],
    na.action = na.omit
  )
  
summary(m_nee_no_fd); check_collinearity(m_nee_no_fd); r.squaredGLMM(m_nee_no_fd); AICc(m_nee_no_fd)

m_nee_with_fd <- glmmTMB(
  formula = as.formula(paste(nee_base_formula, "+ functional_diversity_q1_anomaly_country + (1 | site)")),
  data = dt[!is.na(dt$functional_diversity_q1_anomaly_country)],
  na.action = na.omit
)
summary(m_nee_with_fd); check_collinearity(m_nee_with_fd); r.squaredGLMM(m_nee_with_fd); AICc(m_nee_with_fd)


#Reco ---------
reco_base_formula <- "reco_anomaly_country ~ temperature_reco_anomaly_country +
                  height_x_cover_anomaly_country +
                  sla_anomaly_country + 
                  leaf_area_anomaly_country +
                  mat_anomaly_country"

# PAR
m_reco_no_par <- glmmTMB(
  formula = as.formula(paste(reco_base_formula, "+ (1 | site)")),
  data = dt[!is.na(dt$par_anomaly_country)],
  na.action = na.omit
)
summary(m_reco_no_par); check_collinearity(m_reco_no_par); r.squaredGLMM(m_reco_no_par); AICc(m_reco_no_par)

m_reco_with_par <- glmmTMB(
  formula = as.formula(paste(reco_base_formula, "+ par_anomaly_country + (1 | site)")),
  data = dt[!is.na(dt$par_anomaly_country)],
  na.action = na.omit
)
summary(m_reco_with_par); check_collinearity(m_reco_with_par); r.squaredGLMM(m_reco_with_par); AICc(m_reco_with_par)


# WOODINESS
m_reco_no_wood <- glmmTMB(
  formula = as.formula(paste(reco_base_formula, "+ (1 | site)")),
  data = dt[!is.na(dt$woodiness_anomaly_country)],
  na.action = na.omit
)
summary(m_reco_no_wood); check_collinearity(m_reco_no_wood); r.squaredGLMM(m_reco_no_wood); AICc(m_reco_no_wood)

m_reco_with_wood <- glmmTMB(
  formula = as.formula(paste(reco_base_formula, "+ woodiness_anomaly_country + (1 | site)")),
  data = dt[!is.na(dt$woodiness_anomaly_country)],
  na.action = na.omit
)
summary(m_reco_with_wood); check_collinearity(m_reco_with_wood); r.squaredGLMM(m_reco_with_wood); AICc(m_reco_with_wood)


# GRASSINESS
m_reco_no_grass <- glmmTMB(
  formula = as.formula(paste(reco_base_formula, "+ (1 | site)")),
  data = dt[!is.na(dt$grassiness_anomaly_country)],
  na.action = na.omit
)
summary(m_reco_no_grass); check_collinearity(m_reco_no_grass); r.squaredGLMM(m_reco_no_grass); AICc(m_reco_no_grass)

m_reco_with_grass <- glmmTMB(
  formula = as.formula(paste(reco_base_formula, "+ grassiness_anomaly_country + (1 | site)")),
  data = dt[!is.na(dt$grassiness_anomaly_country)],
  na.action = na.omit
)
summary(m_reco_with_grass); check_collinearity(m_reco_with_grass); r.squaredGLMM(m_reco_with_grass); AICc(m_reco_with_grass)


# SOIL MOISTURE
m_reco_no_soil <- glmmTMB(
  formula = as.formula(paste(reco_base_formula, "+ (1 | site)")),
  data = dt[!is.na(dt$soil_moisture_anomaly_country)],
  na.action = na.omit
)
summary(m_reco_no_soil); check_collinearity(m_reco_no_soil); r.squaredGLMM(m_reco_no_soil); AICc(m_reco_no_soil)

m_reco_with_soil <- glmmTMB(
  formula = as.formula(paste(reco_base_formula, "+ soil_moisture_anomaly_country + (1 | site)")),
  data = dt[!is.na(dt$soil_moisture_anomaly_country)],
  na.action = na.omit
)
summary(m_reco_with_soil); check_collinearity(m_reco_with_soil); r.squaredGLMM(m_reco_with_soil); AICc(m_reco_with_soil)

# SPECIES RICHNESS
m_reco_no_sr <- glmmTMB(
  formula = as.formula(paste(reco_base_formula, "+ (1 | site)")),
  data = dt[!is.na(dt$species_richness_anomaly_country)],
  na.action = na.omit
)
summary(m_reco_no_sr); check_collinearity(m_reco_no_sr); r.squaredGLMM(m_reco_no_sr); AICc(m_reco_no_sr)

m_reco_with_sr <- glmmTMB(
  formula = as.formula(paste(reco_base_formula, "+ species_richness_anomaly_country + (1 | site)")),
  data = dt[!is.na(dt$species_richness_anomaly_country)],
  na.action = na.omit
)
summary(m_reco_with_sr); check_collinearity(m_reco_with_sr); r.squaredGLMM(m_reco_with_sr); AICc(m_reco_with_sr)

# FUNCTIONAL DIVERSITY
m_reco_no_fd <- glmmTMB(
  formula = as.formula(paste(reco_base_formula, "+ (1 | site)")),
  data = dt[!is.na(dt$functional_diversity_q1_anomaly_country)],
  na.action = na.omit
)

summary(m_reco_no_fd); check_collinearity(m_reco_no_fd); r.squaredGLMM(m_reco_no_fd); AICc(m_reco_no_fd)

m_reco_with_fd <- glmmTMB(
  formula = as.formula(paste(reco_base_formula, "+ functional_diversity_q1_anomaly_country + (1 | site)")),
  data = dt[!is.na(dt$functional_diversity_q1_anomaly_country)],
  na.action = na.omit
)
summary(m_reco_with_fd); check_collinearity(m_reco_with_fd); r.squaredGLMM(m_reco_with_fd); AICc(m_reco_with_fd)




#GPP ---------
gpp_base_formula <- "gpp_anomaly_country ~ temperature_gpp_anomaly_country +
                  height_x_cover_anomaly_country +
                  sla_anomaly_country + 
                  leaf_area_anomaly_country +
                  mat_anomaly_country"

# PAR
m_gpp_no_par <- glmmTMB(
  formula = as.formula(paste(gpp_base_formula, "+ (1 | site)")),
  data = dt[!is.na(dt$par_anomaly_country)],
  na.action = na.omit
)
summary(m_gpp_no_par); check_collinearity(m_gpp_no_par); r.squaredGLMM(m_gpp_no_par); AICc(m_gpp_no_par)

m_gpp_with_par <- glmmTMB(
  formula = as.formula(paste(gpp_base_formula, "+ par_anomaly_country + (1 | site)")),
  data = dt[!is.na(dt$par_anomaly_country)],
  na.action = na.omit
)
summary(m_gpp_with_par); check_collinearity(m_gpp_with_par); r.squaredGLMM(m_gpp_with_par); AICc(m_gpp_with_par)


# WOODINESS
m_gpp_no_wood <- glmmTMB(
  formula = as.formula(paste(gpp_base_formula, "+ (1 | site)")),
  data = dt[!is.na(dt$woodiness_anomaly_country)],
  na.action = na.omit
)
summary(m_gpp_no_wood); check_collinearity(m_gpp_no_wood); r.squaredGLMM(m_gpp_no_wood); AICc(m_gpp_no_wood)

m_gpp_with_wood <- glmmTMB(
  formula = as.formula(paste(gpp_base_formula, "+ woodiness_anomaly_country + (1 | site)")),
  data = dt[!is.na(dt$woodiness_anomaly_country)],
  na.action = na.omit
)
summary(m_gpp_with_wood); check_collinearity(m_gpp_with_wood); r.squaredGLMM(m_gpp_with_wood); AICc(m_gpp_with_wood)


# GRASSINESS
m_gpp_no_grass <- glmmTMB(
  formula = as.formula(paste(gpp_base_formula, "+ (1 | site)")),
  data = dt[!is.na(dt$grassiness_anomaly_country)],
  na.action = na.omit
)
summary(m_gpp_no_grass); check_collinearity(m_gpp_no_grass); r.squaredGLMM(m_gpp_no_grass); AICc(m_gpp_no_grass)

m_gpp_with_grass <- glmmTMB(
  formula = as.formula(paste(gpp_base_formula, "+ grassiness_anomaly_country + (1 | site)")),
  data = dt[!is.na(dt$grassiness_anomaly_country)],
  na.action = na.omit
)
summary(m_gpp_with_grass); check_collinearity(m_gpp_with_grass); r.squaredGLMM(m_gpp_with_grass); AICc(m_gpp_with_grass)


# SOIL MOISTURE
m_gpp_no_soil <- glmmTMB(
  formula = as.formula(paste(gpp_base_formula, "+ (1 | site)")),
  data = dt[!is.na(dt$soil_moisture_anomaly_country)],
  na.action = na.omit
)
summary(m_gpp_no_soil); check_collinearity(m_gpp_no_soil); r.squaredGLMM(m_gpp_no_soil); AICc(m_gpp_no_soil)

m_gpp_with_soil <- glmmTMB(
  formula = as.formula(paste(gpp_base_formula, "+ soil_moisture_anomaly_country + (1 | site)")),
  data = dt[!is.na(dt$soil_moisture_anomaly_country)],
  na.action = na.omit
)
summary(m_gpp_with_soil); check_collinearity(m_gpp_with_soil); r.squaredGLMM(m_gpp_with_soil); AICc(m_gpp_with_soil)

# SPECIES RICHNESS
m_gpp_no_sr <- glmmTMB(
  formula = as.formula(paste(gpp_base_formula, "+ (1 | site)")),
  data = dt[!is.na(dt$species_richness_anomaly_country)],
  na.action = na.omit
)
summary(m_gpp_no_sr); check_collinearity(m_gpp_no_sr); r.squaredGLMM(m_gpp_no_sr); AICc(m_gpp_no_sr)

m_gpp_with_sr <- glmmTMB(
  formula = as.formula(paste(gpp_base_formula, "+ species_richness_anomaly_country + (1 | site)")),
  data = dt[!is.na(dt$species_richness_anomaly_country)],
  na.action = na.omit
)
summary(m_gpp_with_sr); check_collinearity(m_gpp_with_sr); r.squaredGLMM(m_gpp_with_sr); AICc(m_gpp_with_sr)

# FUNCTIONAL DIVERSITY
m_gpp_no_fd <- glmmTMB(
  formula = as.formula(paste(gpp_base_formula, "+ (1 | site)")),
  data = dt[!is.na(dt$functional_diversity_q1_anomaly_country)],
  na.action = na.omit
)

summary(m_gpp_no_fd); check_collinearity(m_gpp_no_fd); r.squaredGLMM(m_gpp_no_fd); AICc(m_gpp_no_fd)

m_gpp_with_fd <- glmmTMB(
  formula = as.formula(paste(gpp_base_formula, "+ functional_diversity_q1_anomaly_country + (1 | site)")),
  data = dt[!is.na(dt$functional_diversity_q1_anomaly_country)],
  na.action = na.omit
)
summary(m_gpp_with_fd); check_collinearity(m_gpp_with_fd); r.squaredGLMM(m_gpp_with_fd); AICc(m_gpp_with_fd)

# gotta loop through that thing 

resps <- c("nee_anomaly_country", "reco_anomaly_country", "gpp_anomaly_country")

vars <- c(#"n_percent_anomaly_country",
           # "p_percent_anomaly_country",
            # "cn_ratio_anomaly_country",
             # "c_percent_anomaly_country",
            "par_anomaly_country",
            "soil_moisture_anomaly_country",
            "woodiness_anomaly_country",
            "grassiness_anomaly_country", 
          "species_richness_anomaly_country", 
          "functional_diversity_q1_anomaly_country")

dt_guide <- CJ(variable = vars, 
               response = resps) %>% 
  mutate(tier = paste0(response, "_", variable))

dt_res <- data.frame()

for(i in 1:nrow(dt_guide)){
  
  var <- dt_guide[i, ]$variable
  resp <- dt_guide[i, ]$response
  tier <- dt_guide[i, ]$tier
  
  dt_sub <- dt %>% filter(!is.na(.[[var]]))
  
  #Formulas 
  form_without <- as.formula(paste0(resp, " ~ ", "temperature_nee_anomaly_country +
    height_x_cover_anomaly_country +
    sla_anomaly_country + 
    leaf_area_anomaly_country +
    mat_anomaly_country + (1 | site)"))
  
  form_with <- as.formula(paste0(resp, " ~ ", var, " + temperature_nee_anomaly_country +
    height_x_cover_anomaly_country +
    sla_anomaly_country + 
    leaf_area_anomaly_country +
    mat_anomaly_country + (1 | site)"))
  
  # Model with and without the alternative hypothesis 
  m_without <- glmmTMB(
    formula = form_without,
    data = dt_sub,
    na.action = na.omit
  )
  
  m_with <- glmmTMB(
    formula = form_with,
    data = dt_sub,
    na.action = na.omit
  )
  
  
  delta_aicc <- AICc(m_with) - AICc(m_without)
  delta_aic <- AIC(m_with) - AIC(m_without)
  delta_bic <- BIC(m_with) - BIC(m_without)

  rsq_m_with <-  as.numeric(r.squaredGLMM(m_with)[1])
  rsq_c_with <-  as.numeric(r.squaredGLMM(m_with)[2])
  
  rsq_m_without <-  as.numeric(r.squaredGLMM(m_without)[1])
  rsq_c_without <-  as.numeric(r.squaredGLMM(m_without)[2])
  
  delta_rsq_m <- rsq_m_with - rsq_m_without
  
  tmp <- data.table(tier = tier)
  
  
  tmp <- tmp %>% 
    mutate(
      response = resp,
      new_var = var,
      rsq_m_with = round(rsq_m_with, 3),
      rsq_c_with = round(rsq_c_with, 3),
      rsq_m_without = round(rsq_m_without, 3),
      rsq_c_without = round(rsq_c_without, 3),
      delta_rsq_m = round(delta_rsq_m, 2),
      delta_aicc = delta_aicc, 
      delta_aic = delta_aic, 
      delta_bic = delta_bic, 
      n = nrow(dt_sub))
  
  
  ## extract estimates 
  tidy_m <- broom.mixed::tidy(m_with)
  
  
  ## bring in good shape 
  tmp_est <- tidy_m %>%
    filter(effect == "fixed") %>% 
    mutate(
      ci_ub = estimate + (std.error*1.96),
      ci_lb = estimate - (std.error*1.96), 
      response = resp, 
      tier = tier, 
      new_var = var) %>% 
    filter(!effect == "ran_pars") %>% 
    mutate(group = NULL, effect = NULL, component = NULL) %>% 
    left_join(tmp)
  
  dt_res <- rbind(dt_res, tmp_est)
  
  print(paste0(i, "/", nrow(dt_guide)))
}

dt_clean <- dt_res %>% 
  filter(term == new_var) %>% 
  mutate(
  response_clean = case_when(
    response == "nee_anomaly_country"  ~ "NEE",
    response == "reco_anomaly_country" ~ "Reco",
    response == "gpp_anomaly_country"  ~ "GPP",
    TRUE ~ response
  ),
  
  variable_clean = case_when(
    new_var == "n_percent_anomaly_country"           ~ "Leaf Nitrogen",
    new_var == "p_percent_anomaly_country"           ~ "Leaf Phosphorus",
    new_var == "cn_ratio_anomaly_country"            ~ "Leaf C:N Ratio",
    new_var == "c_percent_anomaly_country"           ~ "Leaf Carbon",
    new_var == "par_anomaly_country"                 ~ "PAR",
    new_var == "soil_moisture_anomaly_country"       ~ "Soil Moisture",
    new_var == "woodiness_anomaly_country"           ~ "Woodiness",
    new_var == "grassiness_anomaly_country"          ~ "Grassiness",
    new_var == "species_richness_anomaly_country"    ~ "Species Richness",
    new_var == "functional_diversity_q1_anomaly_country" ~ "Functional Diversity",
    TRUE ~ new_var
  ), 
  est_ci = paste0(round(estimate, 2), " (", round(ci_lb,2), ";", round(ci_ub,2), ")"), 
  p.value = round(p.value, 3),
  delta_aicc = round(delta_aicc, 2)
) %>%
  dplyr::select(
    Response = response_clean,
    Variable = variable_clean,
    Estimate = est_ci,
    p = p.value,
    R2m = rsq_m_with,
    `Delta R2m` = delta_rsq_m,
    `Delta AICc` = delta_aicc,
    N = n
  )

dt_clean

fwrite(dt_clean, "builds/model_outputs/alternative_hypothesis.csv")

### Get predictions for alternative hypotheses vs climate --------

model_list <- list(
  m_alt_par = m_alt_par,
  m_alt_grass = m_alt_grass,
  m_alt_wood = m_alt_wood,
  # m_alt_n = m_alt_n,
  # m_alt_p = m_alt_p,
  # m_alt_c = m_alt_c,
  # m_alt_cn = m_alt_cn, 
  m_alt_sr = m_alt_sr,
  m_alt_fd = m_alt_fd, 
  m_alt_soil = m_alt_soil
)

dt_results <- data.frame()
dt_pred <- data.frame()
dt_op <- data.frame()
for(i in 1:length(model_list)){
  
  
  m <- model_list[[i]]
  
  m_name <- names(model_list)[i]
  
  response_name <- all.vars(formula(m))[1]
  
  rsq_m <-  as.numeric(r.squaredGLMM(m)[1])
  rsq_c <-  as.numeric(r.squaredGLMM(m)[2])
  
  tmp <- data.table(model_name = m_name)
  
  
  tmp <- tmp %>% 
    mutate(
      response = response_name, 
      rsq_m = round(rsq_m, 3),
      rsq_c = round(rsq_c, 3))
  
  
  ## extract estimates 
  tidy_m <- broom.mixed::tidy(m)
  
  
  ## bring in good shape 
  tmp_est <- tidy_m %>%
    filter(effect == "fixed") %>% 
    mutate(
      ci_ub = estimate + (std.error*1.96),
      ci_lb = estimate - (std.error*1.96), 
      model_name = m_name) %>% 
    filter(!effect == "ran_pars") %>% 
    mutate(group = NULL) %>% 
    left_join(tmp)
  
  dt_results <- rbind(dt_results, tmp_est)
  
  var_names <- tidy_m %>%
    dplyr::select(term) %>%
    filter(!grepl("ntercept", term) & term != "sd__Observation") %>%
    filter(!grepl("\\:", term))
  
  if(nrow(var_names) > 0){
    
    for(j in 1:nrow(var_names)){
      
      var <- var_names[j,] %>% pull()
      
      clean_var = paste0(gsub("_scaled", "", var))
      
      p <- plot_model(m, term = var, type = "pred", ci.lvl = 0.95, show.values = T)
      
      marg_tmp <- p$data %>% 
        as.data.table() %>% 
        rename(var_value = x) %>% 
        mutate(term = var) %>% 
        dplyr::select(-group, -group_col) 
      
      if(j==1){
        marg <- marg_tmp}else{
          marg <- rbind(marg, marg_tmp)}
    }
    
    tmp_pred <- marg %>% 
      mutate(
        model_name = m_name, 
        response = response_name)
    
    dt_pred <- rbind(dt_pred, tmp_pred)
    
  }
  

  print(paste0(m_name, " done (",i, "/", length(model_list), ")"))
}

dt_sig <- dt_results %>% 
  mutate(sig = ifelse(p.value < 0.05, "significant", "non-significant"))

dt_plot <- dt_pred %>% 
  left_join(dt_sig[, c("term", "model_name", "response", "sig")]) %>% 
  mutate(clean_response = case_when(
    model_name == "m_alt_par" ~ "PAR", 
    model_name == "m_alt_grass" ~ "Grassiness", 
    model_name == "m_alt_wood" ~ "Woodiness", 
    model_name == "m_alt_n" ~ "Leaf N", 
    model_name == "m_alt_p" ~ "Leaf P", 
    model_name == "m_alt_c" ~ "Leaf C", 
    model_name == "m_alt_cn" ~ "Leaf C:N", 
    model_name == "m_alt_sr" ~ "Species Richness", 
    model_name == "m_alt_fd" ~ "Functional Diversity", 
    model_name == "m_alt_soil" ~ "Soil Moisture" 
  ), 
  clean_term = "MAT")


dt_long <- dt %>% 
  pivot_longer(cols = unique(dt_pred$response), 
               names_to = "response", values_to = "response_value") %>% 
  mutate(clean_response = case_when(
    response == "par_anomaly_country" ~ "PAR", 
    response == "grassiness_anomaly_country" ~ "Grassiness", 
    response == "woodiness_anomaly_country" ~ "Woodiness", 
    response == "n_percent_anomaly_country" ~ "Leaf N", 
    response == "p_percent_anomaly_country" ~ "Leaf P", 
    response == "c_percent_anomaly_country" ~ "Leaf C", 
    response == "cn_ratio_anomaly_country" ~ "Leaf C:N", 
    response == "species_richness_anomaly_country" ~ "Species Richness", 
    response == "functional_diversity_q1_anomaly_country" ~ "Functional Diversity", 
    response == "soil_moisture_anomaly_country" ~ "Soil Moisture" 
  ))

p_pred <- dt_plot %>% 
  ggplot() +
  geom_point(data = dt_long, 
             aes(x = mat_anomaly_country, y = response_value), alpha = 0.25) +
  geom_ribbon(aes(x = var_value, ymin = conf.low, ymax = conf.high), alpha = 0.25) +
  geom_line(aes(x = var_value, y = predicted, linetype = sig), linewidth = 1.1) +
  scale_linetype_manual(values = c("significant" = "solid", "non-significant" = "dashed")) +
  labs(x = "Response Value", y = "MAT") +
  facet_wrap(~clean_response, scales = "free", ncol = 4) +
  theme(legend.position = "none", 
        legend.box="vertical",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.grid = element_line(color = "seashell"), 
        #axis.title.x = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        #panel.border = element_rect(color = NA), 
        panel.background = element_rect(fill = "snow2"), 
        strip.text.x = element_text(size = 14), 
        strip.text.y = element_text(size = 14, face = "bold"), 
        strip.background = element_rect(fill = "seashell", color = "seashell") )

p_pred
ggsave(plot = p_pred, "builds/plots/alternative_hypotheses_vs_mat.png", 
       dpi = 600, 
       width = 10, height = 5)
