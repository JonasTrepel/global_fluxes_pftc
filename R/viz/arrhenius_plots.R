# test alternative hypothesis 

# funcitons 

predict_with_ci <- function(model, data, response_name,
                            predictor_col = NULL,
                            country_mean = NULL, country_col = NULL) {
  
  preds <- predict(model, newdata = data, se.fit = TRUE)
  
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
  
  # Rename cols
  #colnames(pred_df) <- c("response_name", paste(response_name, c("pred", "se", "ci_lb", "ci_ub"), sep = "_"), "country")
  
  return(pred_df)
}

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
    country, site, plot_id, treatment, gradient,lat,
    
    # fluxes
    nee, reco, gpp,
    
    # environmental 
    downscaled_temp, map, mat,
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
    downscaled_temp_anomaly_country,
    
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
    par_nee_anomaly_country,
    soil_moisture_anomaly_country,
    woodiness_anomaly_country,
    grassiness_anomaly_country
  ) %>% 
  filter(!is.na(sla_cm2_g)) 


dt_mod <- dt %>% 
  # Create absolute temperature (K) and transformed temperature for NEE, RECO, GPP
  mutate(
    temp_k_nee  = temperature_nee  + 273.15,
    temp_k_reco = temperature_reco + 273.15,
    temp_k_gpp  = temperature_gpp  + 273.15,
    
    downscaled_temp_k = downscaled_temp + 273.15, 
    
    temp_k_nee_trans  = -1 / ((8.617e-5)*temp_k_nee),
    temp_k_reco_trans = -1 / ((8.617e-5)*temp_k_reco),
    temp_k_gpp_trans  = -1 / ((8.617e-5)*temp_k_gpp),
    
    downscaled_temp_k_trans = -1 / ((8.617e-5)*downscaled_temp_k)
  ) %>%
  group_by(gradient) %>% 
  mutate(
    temp_k_nee_trans_mean  = mean(temp_k_nee_trans,  na.rm = TRUE),
    temp_k_reco_trans_mean = mean(temp_k_reco_trans, na.rm = TRUE),
    temp_k_gpp_trans_mean  = mean(temp_k_gpp_trans,  na.rm = TRUE),
    
    downscaled_temp_k_trans_mean  = mean(downscaled_temp_k_trans,  na.rm = TRUE),
    
    temp_k_nee_trans_anomaly  = temp_k_nee_trans  - temp_k_nee_trans_mean,
    temp_k_reco_trans_anomaly = temp_k_reco_trans - temp_k_reco_trans_mean,
    temp_k_gpp_trans_anomaly  = temp_k_gpp_trans  - temp_k_gpp_trans_mean,
    
    downscaled_temp_k_trans_anomaly  = downscaled_temp_k_trans  - downscaled_temp_k_trans_mean,
    
    positive_nee_anomaly_country  = nee_anomaly_country  + abs(min(nee_anomaly_country,  na.rm = TRUE)) + 1,
    positive_reco_anomaly_country = reco_anomaly_country + abs(min(reco_anomaly_country, na.rm = TRUE)) + 1,
    positive_gpp_anomaly_country  = gpp_anomaly_country  + abs(min(gpp_anomaly_country,  na.rm = TRUE)) + 1,
    
    log_positive_nee_anomaly_country = log(positive_nee_anomaly_country),
    log_positive_gpp_anomaly_country = log(positive_gpp_anomaly_country),
    log_positive_reco_anomaly_country = log(positive_reco_anomaly_country)
  ) %>% 
  ungroup()

# models 
#### Fluxes predict by temperature ------
m_gpp <- glmmTMB(log_positive_gpp_anomaly_country ~
                       scale(downscaled_temp_k_trans_anomaly) +
                       scale(temp_k_gpp_trans_anomaly) +
                       scale(par_nee_anomaly_country) +
                       (1 | site), 
                     na.action = na.omit,
                     data = dt_mod)
summary(m_gpp); r.squaredGLMM(m_gpp) #0.00

m_nee <- glmmTMB(log_positive_nee_anomaly_country ~
                       scale(downscaled_temp_k_trans_anomaly) +
                       scale(temp_k_nee_trans_anomaly) +
                       scale(par_nee_anomaly_country) +
                       (1 | site), 
                     na.action = na.omit,
                     data = dt_mod)
summary(m_nee); r.squaredGLMM(m_nee) #0.02

m_reco <- glmmTMB(log_positive_reco_anomaly_country ~
                        scale(downscaled_temp_k_trans_anomaly) +
                        scale(temp_k_reco_trans_anomaly) +
                        (1 | site), 
                      na.action = na.omit,
                      data = dt_mod)
summary(m_reco); r.squaredGLMM(m_reco) #0.03


# growing season temperature 
dt_pred_gst <- rbind(
  predict_with_ci(model = m_gpp,
                  data = dt_mod %>% mutate(site = NA, country = NA, 
                                           temp_k_gpp_trans_anomaly = 0, 
                                           par_nee_anomaly_country = 0),
                  response_name = "ln(GPP)",
                  predictor_col = dt_mod$downscaled_temp_k_trans,
                  country_mean = 0, 
                  country_col = dt_mod$country), 
  predict_with_ci(model = m_nee,
                  data = dt_mod %>% mutate(site = NA, country = NA, 
                                           temp_k_nee_trans_anomaly = 0,
                                           par_nee_anomaly_country = 0),
                  response_name = "ln(NEE)",
                  predictor_col = dt_mod$downscaled_temp_k_trans,
                  country_mean = 0, 
                  country_col = dt_mod$country), 
  predict_with_ci(model = m_reco,
                  data = dt_mod %>% mutate(site = NA, country = NA, 
                                          temp_k_reco_trans_anomaly = 0, 
                                          par_nee_anomaly_country = 0),
                  response_name = "ln(Reco)",
                  predictor_col = dt_mod$downscaled_temp_k_trans,
                  country_mean = 0, 
                  country_col = dt_mod$country)) %>% 
  rename(predictor = predictor, 
         flux_type = response) %>% 
  left_join(unique(dt[, c("country", "gradient")]))


glimpse(dt_pred_gst)


p_gst = dt_mod %>% 
  rename(`ln(GPP)` = log_positive_gpp_anomaly_country, 
         `ln(NEE)` = log_positive_nee_anomaly_country, 
         `ln(Reco)` = log_positive_reco_anomaly_country) %>%
  pivot_longer(cols = c("ln(GPP)", "ln(NEE)", "ln(Reco)"), 
               names_to = "flux_type", values_to = "flux_value") %>% 
  mutate(gradient = fct_reorder(gradient, lat), 
         temp = case_when(
           flux_type == "ln(GPP)" ~ temp_k_gpp_trans_anomaly, 
           flux_type == "ln(NEE)" ~ temp_k_nee_trans_anomaly, 
           flux_type == "ln(Reco)" ~ temp_k_reco_trans_anomaly 
         )) %>% 
  ggplot(aes(x = downscaled_temp_k_trans, y = flux_value)) +
  geom_point(alpha = .5, aes(color = gradient)) +
  geom_ribbon(data = dt_pred_gst, aes(x = predictor, ymin = ci_lb, ymax = ci_ub), 
              alpha = .5, fill = "snow3", inherit.aes = FALSE) +
  geom_line(data = dt_pred_gst, aes(x = predictor, y = pred),
            alpha = .75, linetype = "dashed", color = "black", linewidth = 0.75) +
  facet_grid(rows = vars(flux_type), cols = vars(gradient), scales = "free") +
  MetBrewer::scale_color_met_d(name = "Archambault") +
  labs(x = "-1/kT", y = "ln (flux + offset)", subtitle = "Growing Season Temperature", title = "A") +
  theme(legend.position = "none", 
        legend.box="vertical",
        plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.grid = element_line(color = "snow"), 
        #axis.title.x = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1), 
        #panel.border = element_rect(color = NA), 
        panel.background = element_rect(fill = "snow"), 
        strip.text.x = element_text(size = 12), 
        strip.text.y = element_text(size = 12, face = "bold"), 
        strip.background = element_rect(fill = "linen", color = "linen") )

p_gst 

#Instantaneous temperature 
dt_pred_it <- rbind(
  predict_with_ci(model = m_gpp,
                  data = dt_mod %>% mutate(site = NA, country = NA, 
                                           downscaled_temp_k_trans_anomaly = 0,
                                           par_nee_anomaly_country = 0),
                  response_name = "ln(GPP)",
                  predictor_col = dt_mod$temp_k_gpp_trans_anomaly,
                  country_mean = 0, 
                  country_col = dt_mod$country), 
  predict_with_ci(model = m_nee,
                  data = dt_mod %>% mutate(site = NA, country = NA, 
                                           downscaled_temp_k_trans_anomaly = 0,
                                           par_nee_anomaly_country = 0),
                  response_name = "ln(NEE)",
                  predictor_col = dt_mod$temp_k_nee_trans_anomaly,
                  country_mean = 0, 
                  country_col = dt_mod$country), 
  predict_with_ci(model = m_reco,
                  data = dt_mod %>% mutate(site = NA, country = NA, 
                                           downscaled_temp_k_trans_anomaly = 0,
                                           par_nee_anomaly_country = 0),
                  response_name = "ln(Reco)",
                  predictor_col = dt_mod$temp_k_reco_trans_anomaly,
                  country_mean = 0, 
                  country_col = dt_mod$country)) %>% 
  rename(predictor = predictor, 
         flux_type = response) %>% 
  left_join(unique(dt[, c("country", "gradient")]))


glimpse(dt_pred_it)


p_it = dt_mod %>% 
  rename(`ln(GPP)` = log_positive_gpp_anomaly_country, 
         `ln(NEE)` = log_positive_nee_anomaly_country, 
         `ln(Reco)` = log_positive_reco_anomaly_country) %>%
  pivot_longer(cols = c("ln(GPP)", "ln(NEE)", "ln(Reco)"), 
               names_to = "flux_type", values_to = "flux_value") %>% 
  mutate(gradient = fct_reorder(gradient, lat), 
         temp = case_when(
           flux_type == "ln(GPP)" ~ temp_k_gpp_trans_anomaly, 
           flux_type == "ln(NEE)" ~ temp_k_nee_trans_anomaly, 
           flux_type == "ln(Reco)" ~ temp_k_reco_trans_anomaly 
         )) %>% 
  ggplot(aes(x = temp, y = flux_value)) +
  geom_point(alpha = .5, aes(color = gradient)) +
  geom_ribbon(data = dt_pred_it, aes(x = predictor, ymin = ci_lb, ymax = ci_ub), 
              alpha = .5, fill = "snow3", inherit.aes = FALSE) +
  geom_line(data = dt_pred_it, aes(x = predictor, y = pred),
            alpha = .75, linetype = "dashed", color = "black", linewidth = 0.75) +
  facet_grid(rows = vars(flux_type), cols = vars(gradient), scales = "free") +
  MetBrewer::scale_color_met_d(name = "Archambault") +
  labs(x = "-1/kT", y = "ln (flux + offset)", subtitle = "Instantaneous Temperature", title = "B") +
  theme(legend.position = "none", 
        legend.box="vertical",
        plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.grid = element_line(color = "snow"), 
        #axis.title.x = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1), 
        #panel.border = element_rect(color = NA), 
        panel.background = element_rect(fill = "snow"), 
        strip.text.x = element_text(size = 12), 
        strip.text.y = element_text(size = 12, face = "bold"), 
        strip.background = element_rect(fill = "linen", color = "linen") )

p_it 

library(patchwork)

(p_comb <- p_gst / p_it)

ggsave(plot = p_comb, "builds/plots/supplement/arrhenius_plots_fluxes_vs_temp.png", dpi = 600, height = 11, width = 10)


##########################################################################################################
#### get the plots with anomaly ready ####################################################################
##########################################################################################################

# growing season temperature 
dt_pred_gst_a <- rbind(
  predict_with_ci(model = m_gpp,
                  data = dt_mod %>% mutate(site = NA, country = NA, 
                                           temp_k_gpp_trans_anomaly = 0, 
                                           par_nee_anomaly_country = 0),
                  response_name = "GPP",
                  predictor_col = dt_mod$downscaled_temp_k_trans_anomaly,
                  country_mean = 0, 
                  country_col = dt_mod$country), 
  predict_with_ci(model = m_nee,
                  data = dt_mod %>% mutate(site = NA, country = NA, 
                                           temp_k_nee_trans_anomaly = 0,
                                           par_nee_anomaly_country = 0),
                  response_name = "NEE",
                  predictor_col = dt_mod$downscaled_temp_k_trans_anomaly,
                  country_mean = 0, 
                  country_col = dt_mod$country), 
  predict_with_ci(model = m_reco,
                  data = dt_mod %>% mutate(site = NA, country = NA, 
                                           temp_k_reco_trans_anomaly = 0, 
                                           par_nee_anomaly_country = 0),
                  response_name = "Reco",
                  predictor_col = dt_mod$downscaled_temp_k_trans_anomaly,
                  country_mean = 0, 
                  country_col = dt_mod$country)) %>% 
  rename(predictor = predictor, 
         flux_type = response) %>% 
  left_join(unique(dt[, c("country", "gradient")]))


glimpse(dt_pred_gst_a)

c(MetBrewer::met.brewer(name = "Archambault", n = 6))
unique(dt[, c("country", "gradient")])

palette = c("Central Andes" = "#7c4b73", 
            "Drakensberg" = "#88a0dc", 
            "Southern Scandes" = "#e78429",
            "Svalbard" = "#f9d14a",
            "Rocky Mountains" = "#ab3329",
            "Eastern Himalayas"  = "#ed968c")

p_gst_a = dt_mod %>% 
  rename(`GPP` = log_positive_gpp_anomaly_country, 
         `NEE` = log_positive_nee_anomaly_country, 
         `Reco` = log_positive_reco_anomaly_country) %>%
  pivot_longer(cols = c("GPP", "NEE", "Reco"), 
               names_to = "flux_type", values_to = "flux_value") %>% 
  mutate(gradient = fct_reorder(gradient, lat), 
         temp = case_when(
           flux_type == "GPP" ~ temp_k_gpp_trans_anomaly, 
           flux_type == "NEE" ~ temp_k_nee_trans_anomaly, 
           flux_type == "Reco" ~ temp_k_reco_trans_anomaly 
         )) %>% 
  ggplot(aes(x = downscaled_temp_k_trans_anomaly, y = flux_value)) +
  geom_point(alpha = .5, aes(color = gradient)) +
  geom_ribbon(data = dt_pred_gst_a, aes(x = predictor, ymin = ci_lb, ymax = ci_ub), 
              alpha = .1, linetype = "dashed", inherit.aes = FALSE, color = "grey25", linewidth = 0.5) +
  geom_line(data = dt_pred_gst_a, aes(x = predictor, y = pred),
            alpha = 1, linetype = "dashed", linewidth = 0.75, color = "black") +
  facet_wrap(~flux_type, scales = "free") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(x = "-1/kT (deviation from gradient mean)", y = "ln (flux + offset)", subtitle = "Growing Season Temperature", title = "") +
  theme(legend.position = "none", 
        legend.box="vertical",
        plot.title = element_blank(),
        plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
        panel.grid = element_line(color = "snow"), 
        #axis.title.x = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.text.x = element_text(size = 10, angle = 0, hjust = 1), 
        #panel.border = element_rect(color = NA), 
        panel.background = element_rect(fill = "snow"), 
        strip.text.x = element_text(size = 12), 
        strip.text.y = element_text(size = 12, face = "bold"), 
        strip.background = element_rect(fill = "linen", color = "linen") )

p_gst_a 

#Instantaneous temperature 
dt_pred_it <- rbind(
  predict_with_ci(model = m_gpp,
                  data = dt_mod %>% mutate(site = NA, country = NA, 
                                           downscaled_temp_k_trans_anomaly = 0,
                                           par_nee_anomaly_country = 0),
                  response_name = "GPP",
                  predictor_col = dt_mod$temp_k_gpp_trans_anomaly,
                  country_mean = 0, 
                  country_col = dt_mod$country), 
  predict_with_ci(model = m_nee,
                  data = dt_mod %>% mutate(site = NA, country = NA, 
                                           downscaled_temp_k_trans_anomaly = 0,
                                           par_nee_anomaly_country = 0),
                  response_name = "NEE",
                  predictor_col = dt_mod$temp_k_nee_trans_anomaly,
                  country_mean = 0, 
                  country_col = dt_mod$country), 
  predict_with_ci(model = m_reco,
                  data = dt_mod %>% mutate(site = NA, country = NA, 
                                           downscaled_temp_k_trans_anomaly = 0,
                                           par_nee_anomaly_country = 0),
                  response_name = "Reco",
                  predictor_col = dt_mod$temp_k_reco_trans_anomaly,
                  country_mean = 0, 
                  country_col = dt_mod$country)) %>% 
  rename(predictor = predictor, 
         flux_type = response) %>% 
  left_join(unique(dt[, c("country", "gradient")]))


glimpse(dt_pred_it)


p_it_a = dt_mod %>% 
  rename(`GPP` = log_positive_gpp_anomaly_country, 
         `NEE` = log_positive_nee_anomaly_country, 
         `Reco` = log_positive_reco_anomaly_country) %>%
  pivot_longer(cols = c("GPP", "NEE", "Reco"), 
               names_to = "flux_type", values_to = "flux_value") %>% 
  mutate(gradient = fct_reorder(gradient, lat), 
         temp = case_when(
           flux_type == "GPP" ~ temp_k_gpp_trans_anomaly, 
           flux_type == "NEE" ~ temp_k_nee_trans_anomaly, 
           flux_type == "Reco" ~ temp_k_reco_trans_anomaly 
         )) %>% 
  ggplot(aes(x = temp, y = flux_value)) +
  geom_point(alpha = .5, aes(color = gradient)) +
  geom_ribbon(data = dt_pred_it, aes(x = predictor, ymin = ci_lb, ymax = ci_ub), 
              alpha = .1, linetype = "dashed", inherit.aes = FALSE, color = "grey25", linewidth = 0.5) +
  geom_line(data = dt_pred_it, aes(x = predictor, y = pred),
            alpha = 1, linetype = "dashed", linewidth = 0.75, color = "black") +
  facet_wrap(~flux_type, scales = "free") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(x = "-1/kT (deviation from gradient mean)", y = "ln (flux + offset)", subtitle = "Instantaneous Temperature") +
  theme(legend.position = "none", 
        legend.box="vertical",
        plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
        panel.grid = element_line(color = "snow"), 
        #axis.title.x = element_blank(), 
        axis.text = element_text(size = 10), 
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10, angle = 0, hjust = 1), 
        #panel.border = element_rect(color = NA), 
        panel.background = element_rect(fill = "snow"), 
        strip.text.x = element_text(size = 12), 
        strip.text.y = element_text(size = 12, face = "bold"), 
        strip.background = element_rect(fill = "linen", color = "linen") )

p_it_a 

library(patchwork)

(p_comb_a <- p_gst_a / p_it_a)

ggsave(plot = p_comb_a, "builds/plots/supplement/arrhenius_plots_fluxes_vs_temp_anomalies.png", dpi = 600, height = 6, width = 10)
