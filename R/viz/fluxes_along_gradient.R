# test alternative hypothesis 

# funcitons 

predict_with_ci <- function(model, data, response_name,
                            predictor_col = NULL,
                            country_mean = NULL, country_col = NULL) {
  
  preds <- predict(model, newdata = data %>% mutate(site = NA, country = NA, 
                                                    par_nee_anomaly_country = 0), se.fit = TRUE)
  
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
    country, site, plot_id, treatment, gradient,
    
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


dt %>% ggplot() +
  geom_point(aes(x = downscaled_temp, y = gpp_anomaly_country)) +
  facet_wrap(~country, scales = "free")


#### Fluxes predict by downscaled_temp ------
m_ele_gpp <- glmmTMB(gpp_anomaly_country ~
                    scale(downscaled_temp_anomaly_country) +
                      scale(par_nee_anomaly_country) +
                    (1 | site), 
                  na.action = na.omit,
                  data = dt)
summary(m_ele_gpp); r.squaredGLMM(m_ele_gpp) #0.00

m_ele_nee <- glmmTMB(nee_anomaly_country ~
                  scale(downscaled_temp_anomaly_country) +
                     scale(par_nee_anomaly_country) +
                   (1 | site), 
                 na.action = na.omit,
                 data = dt)
summary(m_ele_nee); r.squaredGLMM(m_ele_nee) #0.02

m_ele_reco <- glmmTMB(reco_anomaly_country ~
                   downscaled_temp_anomaly_country +
                   (1 | site), 
                 na.action = na.omit,
                 data = dt)
summary(m_ele_reco); r.squaredGLMM(m_ele_reco) #0.03



dt_ele_pred <- rbind(
                predict_with_ci(model = m_ele_gpp,
                          data = dt,
                          response_name = "GPP",
                          predictor_col = dt$downscaled_temp,
                          country_mean = dt$gpp_country_mean, 
                          country_col = dt$country), 
                predict_with_ci(model = m_ele_nee,
                                data = dt,
                                response_name = "NEE",
                                predictor_col = dt$downscaled_temp,
                                country_mean = dt$nee_country_mean, 
                                country_col = dt$country), 
                predict_with_ci(model = m_ele_reco,
                                data = dt,
                                response_name = "Reco",
                                predictor_col = dt$downscaled_temp,
                                country_mean = dt$reco_country_mean, 
                                country_col = dt$country)) %>% 
  rename(downscaled_temp = predictor, 
         flux_type = response) %>% 
  left_join(unique(dt[, c("country", "gradient")]))

#### MANUSCRIPT PLOT -----------------------

dt <- dt %>%
  mutate(gradient = fct_reorder(gradient, lat))

dt_ele_pred <- dt_ele_pred %>%
  mutate(gradient = factor(gradient, levels = levels(dt$gradient)))

p_elev <- dt %>% 
  rename(GPP = gpp, 
         NEE = nee, 
         Reco = reco) %>%
  pivot_longer(cols = c("GPP", "NEE", "Reco"), 
               names_to = "flux_type", values_to = "flux_value") %>% 
  mutate(gradient = fct_reorder(gradient, lat)) %>% 
  ggplot(aes(x = downscaled_temp, y = flux_value)) +
  geom_point(alpha = .5, aes(color = gradient)) +
  geom_ribbon(data = dt_ele_pred, aes(x = downscaled_temp, ymin = ci_lb, ymax = ci_ub), 
              alpha = .5, fill = "snow3", inherit.aes = FALSE) +
  geom_line(data = dt_ele_pred, aes(x = downscaled_temp, y = pred),
            alpha = .75, linetype = "dashed", color = "black", linewidth = 0.75) +
  facet_grid(rows = vars(flux_type), cols = vars(gradient), scales = "free") +
  MetBrewer::scale_color_met_d(name = "Archambault") +
  labs(x = "Growing Season Temperature (°C)", y = "µmol m⁻² s⁻¹", title = NULL) +
  theme(legend.position = "none", 
        legend.box="vertical",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.grid = element_line(color = "snow"), 
        #axis.title.x = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1), 
        #panel.border = element_rect(color = NA), 
        panel.background = element_rect(fill = "snow"), 
        strip.text.x = element_text(size = 12), 
        strip.text.y = element_text(size = 12, face = "bold"), 
        strip.background = element_rect(fill = "linen", color = "linen") )

p_elev
ggsave(plot = p_elev, "builds/plots/fluxes_vs_downscaled_temp.png", dpi = 600, height = 5, width = 10)

