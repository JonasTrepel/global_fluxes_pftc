# test alternative hypothesis 

# funcitons 

predict_with_ci <- function(model, data, response_name,
                            predictor_col = NULL,
                            country_mean = NULL, country_col = NULL) {
  
  preds <- predict(model, newdata = data %>% mutate(site = NA, country = NA), se.fit = TRUE)
  
  # Calculate cis 
  pred_df <- data.frame(
    response = response_name,
    predictor = predictor_col, 
    pred = preds$fit,
    se = preds$se.fit
  ) %>%
    mutate(
      ci_lb = pred - 1.96 * se,
      ci_ub = pred + 1.96 * se
    )
  
  # add offset back 
  if (!is.null(country_mean)) {
    pred_df <- pred_df %>%
      mutate(
        pred = pred + country_mean,
        ci_lb = ci_lb + country_mean,
        ci_ub = ci_ub + country_mean, 
        country = country_col
      )
  }
  
  # Rename columns appropriately
  #colnames(pred_df) <- c("response_name", paste(response_name, c("pred", "se", "ci_lb", "ci_ub"), sep = "_"), "country")
  
  return(pred_df)
}

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


dt <- fread("data/processed_data/clean_data/global_fluxes_main_data.csv") %>% 
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


dt %>% ggplot() +
  geom_point(aes(x = elevation, y = gpp_anomaly_country)) +
  facet_wrap(~country, scales = "free")


#### Fluxes predict by elevation ------
m_ele_gpp <- glmmTMB(gpp_anomaly_country ~
                    elevation_anomaly_country +
                    (1 | site), 
                  na.action = na.omit,
                  data = dt)
summary(m_ele_gpp); r.squaredGLMM(m_ele_gpp) #0.00

m_ele_nee <- glmmTMB(nee_anomaly_country ~
                   elevation_anomaly_country +
                   (1 | site), 
                 na.action = na.omit,
                 data = dt)
summary(m_ele_nee); r.squaredGLMM(m_ele_nee) #0.02

m_ele_reco <- glmmTMB(reco_anomaly_country ~
                   elevation_anomaly_country +
                   (1 | site), 
                 na.action = na.omit,
                 data = dt)
summary(m_ele_reco); r.squaredGLMM(m_ele_reco) #0.03



dt_ele_pred <- rbind(
                predict_with_ci(model = m_ele_gpp,
                          data = dt,
                          response_name = "gpp",
                          predictor_col = dt$elevation,
                          country_mean = dt$gpp_country_mean, 
                          country_col = dt$country), 
                predict_with_ci(model = m_ele_nee,
                                data = dt,
                                response_name = "nee",
                                predictor_col = dt$elevation,
                                country_mean = dt$nee_country_mean, 
                                country_col = dt$country), 
                predict_with_ci(model = m_ele_reco,
                                data = dt,
                                response_name = "reco",
                                predictor_col = dt$elevation,
                                country_mean = dt$reco_country_mean, 
                                country_col = dt$country)) %>% 
  rename(elevation = predictor, 
         flux_type = response)



p_elev <- dt %>% 
  pivot_longer(cols = c("gpp", "nee", "reco"), 
               names_to = "flux_type", values_to = "flux_value") %>% 
  ggplot(aes(x = elevation, y = flux_value)) +
  geom_point(alpha = .5, color = "black") +
  geom_ribbon(data = dt_ele_pred, aes(x = elevation, ymin = ci_lb, ymax = ci_ub), 
              alpha = .1, color = "grey", inherit.aes = FALSE) +
  geom_line(data = dt_ele_pred, aes(x = elevation, y = pred),
            alpha = .75, linewidth = 1.1, linetype = "dashed", color = "black") +
  facet_grid(rows = vars(flux_type), cols = vars(country), scales = "free") +
  labs(x = "Elevation (m)", y = "Flux Value", title = "Fluxes vs. Elevation") +
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

p_elev
ggsave(plot = p_elev, "builds/plots/supplement/fluxes_vs_elevation.png", dpi = 600, height = 5, width = 10)


#### Fluxes predicted by Temperature -----------------------

m_mat_gpp <- glmmTMB(gpp_anomaly_country ~
                       mat_anomaly_country +
                       (1 | site), 
                     na.action = na.omit,
                     data = dt)
summary(m_mat_gpp); r.squaredGLMM(m_mat_gpp) #R2c ~ 0

m_mat_nee <- glmmTMB(nee_anomaly_country ~
                       mat_anomaly_country +
                       (1 | site), 
                     na.action = na.omit,
                     data = dt)
summary(m_mat_nee); r.squaredGLMM(m_mat_nee) #R2c = 0.02

m_mat_reco <- glmmTMB(reco_anomaly_country ~
                        mat_anomaly_country +
                        (1 | site), 
                      na.action = na.omit,
                      data = dt)
summary(m_mat_reco); r.squaredGLMM(m_mat_reco) #R2c = 0.02



dt_mat_pred <- rbind(
  predict_with_ci(model = m_mat_gpp,
                  data = dt,
                  response_name = "gpp",
                  predictor_col = dt$mat,
                  country_mean = dt$gpp_country_mean, 
                  country_col = dt$country), 
  predict_with_ci(model = m_mat_nee,
                  data = dt,
                  response_name = "nee",
                  predictor_col = dt$mat,
                  country_mean = dt$nee_country_mean, 
                  country_col = dt$country), 
  predict_with_ci(model = m_mat_reco,
                  data = dt,
                  response_name = "reco",
                  predictor_col = dt$mat,
                  country_mean = dt$reco_country_mean, 
                  country_col = dt$country)) %>% 
  rename(mat = predictor, 
         flux_type = response)



p_mat <- dt %>% 
  pivot_longer(cols = c("gpp", "nee", "reco"), 
               names_to = "flux_type", values_to = "flux_value") %>% 
  ggplot(aes(x = mat, y = flux_value)) +
  geom_point(alpha = .5, color = "black") +
  geom_ribbon(data = dt_mat_pred, aes(x = mat, ymin = ci_lb, ymax = ci_ub), 
              alpha = .1, color = "grey", inherit.aes = FALSE) +
  geom_line(data = dt_mat_pred, aes(x = mat, y = pred),
            alpha = .75, linewidth = 1.1, linetype = "dashed", color = "black") +
  facet_grid(rows = vars(flux_type), cols = vars(country), scales = "free") +
  labs(x = "MAT (°C)", y = "Flux Value", title = "Fluxes vs. MAT") +
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

p_mat

#Fluxes predicted by local temperature ---------------------

m_temp_gpp <- glmmTMB(gpp_anomaly_country ~
                        temperature_gpp_anomaly_country +
                        (1 | site), 
                      na.action = na.omit,
                      data = dt)
summary(m_temp_gpp); r.squaredGLMM(m_temp_gpp) #R2c = 0.0319909

m_temp_nee <- glmmTMB(nee_anomaly_country ~
                        temperature_nee_anomaly_country +
                        (1 | site), 
                      na.action = na.omit,
                      data = dt)
summary(m_temp_nee); r.squaredGLMM(m_temp_nee) #R2c = 0.006. 

m_temp_reco <- glmmTMB(reco_anomaly_country ~
                         temperature_reco_anomaly_country +
                         (1 | site), 
                       na.action = na.omit,
                       data = dt)
summary(m_temp_reco); r.squaredGLMM(m_temp_reco) # 0.04

# Predictions and confidence intervals for plots
dt_temp_pred <- rbind(
  predict_with_ci(model = m_temp_gpp,
                  data = dt,
                  response_name = "gpp",
                  predictor_col = dt$temperature_gpp,
                  country_mean = dt$gpp_country_mean, 
                  country_col = dt$country), 
  predict_with_ci(model = m_temp_nee,
                  data = dt,
                  response_name = "nee",
                  predictor_col = dt$temperature_nee,
                  country_mean = dt$nee_country_mean, 
                  country_col = dt$country), 
  predict_with_ci(model = m_temp_reco,
                  data = dt,
                  response_name = "reco",
                  predictor_col = dt$temperature_reco,
                  country_mean = dt$reco_country_mean, 
                  country_col = dt$country)) %>% 
  rename(temperature_gpp = predictor, 
         flux_type = response) 

# Plotting the fluxes vs. temperature

flux_temp <- dt %>% 
  pivot_longer(cols = c("temperature_nee", "temperature_gpp", "temperature_reco"), 
               names_to = "flux_temp_type", values_to = "flux_temp") %>% 
  mutate(flux_type = gsub("temperature_", "", flux_temp_type)) %>% 
  select(flux_type, plot_id, flux_temp)

p_temp <- dt %>% 
  pivot_longer(cols = c("gpp", "nee", "reco"), 
               names_to = "flux_type", values_to = "flux_value") %>% 
  left_join(flux_temp) %>% 
  ggplot(aes(x = flux_temp, y = flux_value)) +
  geom_point(alpha = .5, color = "black") +
  geom_ribbon(data = dt_temp_pred, aes(x = temperature_gpp, ymin = ci_lb, ymax = ci_ub), 
              alpha = .1, color = "grey", inherit.aes = FALSE) +
  geom_line(data = dt_temp_pred %>% filter(flux_type %in% c("gpp", "reco")), aes(x = temperature_gpp, y = pred),
            alpha = .75, linewidth = 1.1, linetype = "solid", color = "black") +
  geom_line(data = dt_temp_pred %>% filter(flux_type %in% c("nee")), aes(x = temperature_gpp, y = pred),
            alpha = .75, linewidth = 1.1, linetype = "dashed", color = "black") +
  facet_grid(rows = vars(flux_type), cols = vars(country), scales = "free") +
  labs(x = "Local (instantanous) Temperature (°C)", y = "Flux Value", title = "Fluxes vs. Local Temperature") +
  theme(legend.position = "none", 
        legend.box="vertical",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.grid = element_line(color = "seashell"), 
        axis.text = element_text(size = 12), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        panel.background = element_rect(fill = "snow2"), 
        strip.text.x = element_text(size = 14), 
        strip.text.y = element_text(size = 14, face = "bold"), 
        strip.background = element_rect(fill = "seashell", color = "seashell") )

p_temp

dt %>% select(temperature_nee, temperature_reco, temperature_gpp) %>% 
  filter(complete.cases(.)) %>% cor() #cor temp GPP with both = 0.97

dt %>% select(mat_anomaly_country, elevation_anomaly_country, temperature_gpp_anomaly_country) %>% 
  filter(complete.cases(.)) %>% cor()


p_c <- grid.arrange(p_mat, p_temp)
ggsave(plot = p_c, "builds/plots/fluxes_vs_tmep_and_climate.png", dpi = 600, height = 10, width = 10)

###### Check fluxes vs climate worldwide 
m_mat_gpp_g <- glmmTMB(gpp ~
                       mat +
                       (1 | country/site), 
                     na.action = na.omit,
                     data = dt)
summary(m_mat_gpp_g)

m_mat_nee_g <- glmmTMB(nee ~
                       mat +
                       (1 | country/site), 
                     na.action = na.omit,
                     data = dt)
summary(m_mat_nee_g)

m_mat_reco_g <- glmmTMB(reco ~
                          mat +
                          (1 | country/site), 
                      na.action = na.omit,
                      data = dt)
summary(m_mat_reco_g)



dt_mat_pred_g <- rbind(
  predict_with_ci(model = m_mat_gpp_g,
                  data = dt,
                  response_name = "gpp",
                  predictor_col = dt$mat,
                  country_mean = 0, 
                  country_col = dt$country), 
  predict_with_ci(model = m_mat_nee_g,
                  data = dt,
                  response_name = "nee",
                  predictor_col = dt$mat,
                  country_mean = 0, 
                  country_col = dt$country), 
  predict_with_ci(model = m_mat_reco_g,
                  data = dt,
                  response_name = "reco",
                  predictor_col = dt$mat,
                  country_mean = 0, 
                  country_col = dt$country)) %>% 
  rename(mat = predictor, 
         flux_type = response)



p_mat_g <- dt %>% 
  pivot_longer(cols = c("gpp", "nee", "reco"), 
               names_to = "flux_type", values_to = "flux_value") %>% 
  ggplot(aes(x = mat, y = flux_value)) +
  geom_point(alpha = .5, color = "black") +
  geom_ribbon(data = dt_mat_pred_g, aes(x = mat, ymin = ci_lb, ymax = ci_ub), 
              alpha = .1, color = "grey", inherit.aes = FALSE) +
  geom_line(data = dt_mat_pred_g, aes(x = mat, y = pred),
            alpha = .75, linewidth = 1.1, linetype = "dashed", color = "black") +
  facet_wrap(~flux_type, scales = "free") +
  labs(x = "MAT (°C)", y = "Flux Value", title = "Fluxes vs. MAT (across gradients)", 
       subtitle = "flux ~ MAT + (1 | country/site)") +
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

p_mat_g


ggsave(plot = p_mat_g, "builds/plots/supplement/fluxes_vs_mat_across_gradients.png", dpi = 600, height = 3.75, width = 9)


######### Check models including both climate and temp ---------------------

m_b_gpp <- glmmTMB(gpp_anomaly_country ~
                     temperature_gpp_anomaly_country +
                     mat_anomaly_country + 
                     (1 | site), 
                   na.action = na.omit,
                   data = dt)
summary(m_b_gpp); r.squaredGLMM(m_b_gpp) #R2c = 0.032

m_b_nee <- glmmTMB(nee_anomaly_country ~
                        temperature_nee_anomaly_country +
                        mat_anomaly_country + 
                        (1 | site), 
                      na.action = na.omit,
                      data = dt)
summary(m_b_nee); r.squaredGLMM(m_b_nee) #R2c = 0.02. 

m_b_reco <- glmmTMB(reco_anomaly_country ~
                         temperature_reco_anomaly_country +
                         mat_anomaly_country + 
                         (1 | site), 
                       na.action = na.omit,
                       data = dt)
summary(m_b_reco); r.squaredGLMM(m_b_reco) #R2c = 0.08
#NOTE: results are the same when using elevation instead of MAT

##### EXCLUDE SOUTH AFRICA --------

dt_nsa <- dt %>% filter(country != "South Africa")

#### Fluxes predict by elevation ------
m_ele_gpp_nsa <- glmmTMB(gpp_anomaly_country ~
                           elevation_anomaly_country +
                           (1 | site), 
                         na.action = na.omit,
                         data = dt_nsa)
summary(m_ele_gpp_nsa); r.squaredGLMM(m_ele_gpp_nsa) #R2m: 0.02

m_ele_nee_nsa <- glmmTMB(nee_anomaly_country ~
                           elevation_anomaly_country +
                           (1 | site), 
                         na.action = na.omit,
                         data = dt_nsa) 
summary(m_ele_nee_nsa);  r.squaredGLMM(m_ele_nee_nsa) #R2m: 0.1289201

m_ele_reco_nsa <- glmmTMB(reco_anomaly_country ~
                            elevation_anomaly_country +
                            (1 | site), 
                          na.action = na.omit,
                          data = dt_nsa)
summary(m_ele_reco_nsa);   r.squaredGLMM(m_ele_reco_nsa) #R2m: 0.15

# Prediction dataframe for elevation-based models
dt_ele_pred_nsa <- rbind(
  predict_with_ci(model = m_ele_gpp_nsa,
                  data = dt_nsa,
                  response_name = "gpp",
                  predictor_col = dt_nsa$elevation,
                  country_mean = dt_nsa$gpp_country_mean, 
                  country_col = dt_nsa$country), 
  predict_with_ci(model = m_ele_nee_nsa,
                  data = dt_nsa,
                  response_name = "nee",
                  predictor_col = dt_nsa$elevation,
                  country_mean = dt_nsa$nee_country_mean, 
                  country_col = dt_nsa$country), 
  predict_with_ci(model = m_ele_reco_nsa,
                  data = dt_nsa,
                  response_name = "reco",
                  predictor_col = dt_nsa$elevation,
                  country_mean = dt_nsa$reco_country_mean, 
                  country_col = dt_nsa$country)) %>% 
  rename(elevation = predictor, 
         flux_type = response)

# Plot for elevation-based fluxes (NSA)
p_elev_nsa <- dt_nsa %>% 
  pivot_longer(cols = c("gpp", "nee", "reco"), 
               names_to = "flux_type", values_to = "flux_value") %>% 
  ggplot(aes(x = elevation, y = flux_value)) +
  geom_point(alpha = .5, color = "black") +
  geom_ribbon(data = dt_ele_pred_nsa, aes(x = elevation, ymin = ci_lb, ymax = ci_ub), 
              alpha = .1, color = "grey", inherit.aes = FALSE) +
  geom_line(data = dt_ele_pred_nsa %>% filter(!flux_type == "gpp"), aes(x = elevation, y = pred),
            alpha = .75, linewidth = 1.1, linetype = "solid", color = "black") +
  geom_line(data = dt_ele_pred_nsa %>% filter(flux_type == "gpp"), aes(x = elevation, y = pred),
            alpha = .75, linewidth = 1.1, linetype = "dashed", color = "black") +
  facet_grid(rows = vars(flux_type), cols = vars(country), scales = "free") +
  labs(x = "Elevation (m)", y = "Flux Value", title = "Fluxes vs. Elevation (without South Africa)") +
  theme(legend.position = "none", 
        legend.box="vertical",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.grid = element_line(color = "seashell"), 
        axis.text = element_text(size = 12), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        panel.background = element_rect(fill = "snow2"), 
        strip.text.x = element_text(size = 14), 
        strip.text.y = element_text(size = 14, face = "bold"), 
        strip.background = element_rect(fill = "seashell", color = "seashell"))

p_elev_nsa


#### Fluxes predicted by Temperature -----------------------
m_mat_gpp_nsa <- glmmTMB(gpp_anomaly_country ~
                           mat_anomaly_country +
                           (1 | site), 
                         na.action = na.omit,
                         data = dt_nsa)
summary(m_mat_gpp_nsa); r.squaredGLMM(m_mat_gpp_nsa) #R2m: 0.01094702


m_mat_nee_nsa <- glmmTMB(nee_anomaly_country ~
                           mat_anomaly_country +
                           (1 | site), 
                         na.action = na.omit,
                         data = dt_nsa)
summary(m_mat_nee_nsa); r.squaredGLMM(m_mat_nee_nsa) #R2m: 0.1465595

m_mat_reco_nsa <- glmmTMB(reco_anomaly_country ~
                            mat_anomaly_country +
                            (1 | site), 
                          na.action = na.omit,
                          data = dt_nsa)
summary(m_mat_reco_nsa); r.squaredGLMM(m_mat_reco_nsa) #R2m: 0.127445

dt_mat_pred_nsa <- rbind(
  predict_with_ci(model = m_mat_gpp_nsa,
                  data = dt_nsa,
                  response_name = "gpp",
                  predictor_col = dt_nsa$mat,
                  country_mean = dt_nsa$gpp_country_mean, 
                  country_col = dt_nsa$country), 
  predict_with_ci(model = m_mat_nee_nsa,
                  data = dt_nsa,
                  response_name = "nee",
                  predictor_col = dt_nsa$mat,
                  country_mean = dt_nsa$nee_country_mean, 
                  country_col = dt_nsa$country), 
  predict_with_ci(model = m_mat_reco_nsa,
                  data = dt_nsa,
                  response_name = "reco",
                  predictor_col = dt_nsa$mat,
                  country_mean = dt_nsa$reco_country_mean, 
                  country_col = dt_nsa$country)) %>% 
  rename(mat = predictor, 
         flux_type = response)

p_mat_nsa <- dt_nsa %>% 
  pivot_longer(cols = c("gpp", "nee", "reco"), 
               names_to = "flux_type", values_to = "flux_value") %>% 
  ggplot(aes(x = mat, y = flux_value)) +
  geom_point(alpha = .5, color = "black") +
  geom_ribbon(data = dt_mat_pred_nsa, aes(x = mat, ymin = ci_lb, ymax = ci_ub), 
              alpha = .1, color = "grey", inherit.aes = FALSE) +
  geom_line(data = dt_mat_pred_nsa %>% filter(!flux_type == "gpp"), aes(x = mat, y = pred),
            alpha = .75, linewidth = 1.1, linetype = "solid", color = "black") +
  geom_line(data = dt_mat_pred_nsa %>% filter(flux_type == "gpp"), aes(x = mat, y = pred),
            alpha = .75, linewidth = 1.1, linetype = "dashed", color = "black") +
  facet_grid(rows = vars(flux_type), cols = vars(country), scales = "free") +
  labs(x = "MAT (°C)", y = "Flux Value", title = "Fluxes vs. MAT (without South Africa)") +
  theme(legend.position = "none", 
        legend.box="vertical",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.grid = element_line(color = "seashell"), 
        axis.text = element_text(size = 12), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        panel.background = element_rect(fill = "snow2"), 
        strip.text.x = element_text(size = 14), 
        strip.text.y = element_text(size = 14, face = "bold"), 
        strip.background = element_rect(fill = "seashell", color = "seashell"))

p_mat_nsa

# Save plot
p_c_nsa <- grid.arrange(p_elev_nsa, p_mat_nsa)
ggsave(plot = p_c_nsa, "builds/plots/fluxes_vs_elevation_and_climate_without_south_africa.png", dpi = 600, height = 10, width = 10)


###### Check fluxes vs climate worldwide 
m_mat_gpp_g_nsa <- glmmTMB(gpp ~
                             mat +
                             (1 | country/site), 
                           na.action = na.omit,
                           data = dt_nsa)
summary(m_mat_gpp_g_nsa)

m_mat_nee_g_nsa <- glmmTMB(nee ~
                             mat +
                             (1 | country/site), 
                           na.action = na.omit,
                           data = dt_nsa)
summary(m_mat_nee_g_nsa)

m_mat_reco_g_nsa <- glmmTMB(reco ~
                              mat +
                              (1 | country/site), 
                            na.action = na.omit,
                            data = dt_nsa)
summary(m_mat_reco_g_nsa)

# Predictions and confidence intervals
dt_mat_pred_g_nsa <- rbind(
  predict_with_ci(model = m_mat_gpp_g_nsa,
                  data = dt_nsa,
                  response_name = "gpp",
                  predictor_col = dt_nsa$mat,
                  country_mean = 0, 
                  country_col = dt_nsa$country), 
  predict_with_ci(model = m_mat_nee_g_nsa,
                  data = dt_nsa,
                  response_name = "nee",
                  predictor_col = dt_nsa$mat,
                  country_mean = 0, 
                  country_col = dt_nsa$country), 
  predict_with_ci(model = m_mat_reco_g_nsa,
                  data = dt_nsa,
                  response_name = "reco",
                  predictor_col = dt_nsa$mat,
                  country_mean = 0, 
                  country_col = dt_nsa$country)) %>% 
  rename(mat = predictor, 
         flux_type = response)

# Plotting
p_mat_g_nsa <- dt_nsa %>% 
  pivot_longer(cols = c("gpp", "nee", "reco"), 
               names_to = "flux_type", values_to = "flux_value") %>% 
  ggplot(aes(x = mat, y = flux_value)) +
  geom_point(alpha = .5, color = "black") +
  geom_ribbon(data = dt_mat_pred_g_nsa, aes(x = mat, ymin = ci_lb, ymax = ci_ub), 
              alpha = .1, color = "grey", inherit.aes = FALSE) +
  geom_line(data = dt_mat_pred_g_nsa, aes(x = mat, y = pred),
            alpha = .75, linewidth = 1.1, linetype = "dashed", color = "black") +
  geom_line(data = dt_mat_pred_g_nsa %>% filter(flux_type == "reco"), aes(x = mat, y = pred),
            alpha = .75, linewidth = 1.1, linetype = "solid", color = "black") +
  facet_wrap(~flux_type, scales = "free") +
  labs(x = "MAT (°C)", y = "Flux Value", 
       title = "Fluxes vs. MAT (across gradients, without South Africa)", 
       subtitle = "flux ~ MAT + (1 | country/site)") +
  theme(legend.position = "none", 
        legend.box = "vertical",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.grid = element_line(color = "seashell"), 
        axis.text = element_text(size = 12), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        panel.background = element_rect(fill = "snow2"), 
        strip.text.x = element_text(size = 14), 
        strip.text.y = element_text(size = 14, face = "bold"), 
        strip.background = element_rect(fill = "seashell", color = "seashell"))

# Plot object
p_mat_g_nsa


ggsave(plot = p_mat_g_nsa, "builds/plots/supplement/fluxes_vs_mat_across_gradients_without_south_africa.png", dpi = 600, height = 3.75, width = 9)





##### Random slopes -------
m_glmm <- glmmTMB(gpp ~
          temperature_gpp +
          height_x_cover +
          sla_cm2_g + 
          leaf_area_cm2 +
          mat +
          (temperature_gpp + height_x_cover + sla_cm2_g + leaf_area_cm2 + mat | country) +
          (1 | site), 
        na.action = na.omit,
        data = dt)
summary(m_glmm); r.squaredGLMM(m_glmm)


m_b <- brms::brm(gpp ~
                    temperature_gpp +
                    height_x_cover +
                    sla_cm2_g + 
                    leaf_area_cm2 +
                    mat +
                    (temperature_gpp + height_x_cover + sla_cm2_g + leaf_area_cm2 + mat | country) +
                    (1 | site), 
                  data = dt)
summary(m_b); brms::bayes_R2(m_b)

