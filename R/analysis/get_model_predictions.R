# get model predictions and plot observed vs predicted...

library(glmmTMB)
library(tidyverse)
library(data.table)
library(sjPlot)
library(MetBrewer)

dt_raw <- fread("data/processed_data/clean_data/global_fluxes_main_data.csv") %>% 
  dplyr::select(
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



dt %>% 
  ggplot() +
  geom_point(aes(x = leaf_area_anomaly_country, y = gpp_anomaly_country)) +
  geom_smooth(aes(x = leaf_area_anomaly_country, y = gpp_anomaly_country), method = "lm") +
  facet_wrap(~country, scales = "free_x")


### Models --------------


# model for temperature
m_temp <- glmmTMB(temperature_gpp_anomaly_country ~ 
          mat_anomaly_country +
          ( 1 | site),
        na.action = na.omit,
        data = dt)


# model for veg volume / biomass 
m_vv <- glmmTMB(height_x_cover_anomaly_country ~ 
          mat_anomaly_country +
          ( 1 | site),
        na.action = na.omit,
        data = dt)

# model for sla 
m_sla <- glmmTMB(sla_anomaly_country ~ 
          mat_anomaly_country +
          ( 1 | site),
        na.action = na.omit,
        data = dt)

# model for lead area
m_la <- glmmTMB(leaf_area_anomaly_country ~ 
          mat_anomaly_country +
          ( 1 | site), 
        na.action = na.omit,
        data = dt)

#model for NEE
m_nee <- glmmTMB(nee_anomaly_country ~
          temperature_nee_anomaly_country +
          height_x_cover_anomaly_country +
          sla_anomaly_country + 
          leaf_area_anomaly_country +
          mat_anomaly_country +
          ( 1 | site), 
        na.action = na.omit,
        data = dt)

#model for Reco
m_reco <- glmmTMB(reco_anomaly_country ~
                   temperature_reco_anomaly_country +
                   height_x_cover_anomaly_country +
                   sla_anomaly_country + 
                   leaf_area_anomaly_country +
                   mat_anomaly_country +
                   ( 1 | site), 
                 na.action = na.omit,
                 data = dt)

m_gpp <- glmmTMB(gpp_anomaly_country ~
                    temperature_reco_anomaly_country +
                    height_x_cover_anomaly_country +
                    sla_anomaly_country + 
                    leaf_area_anomaly_country +
                    mat_anomaly_country +
                    ( 1 | site), 
                  na.action = na.omit,
                  data = dt)


## Extract estimates and predict --------------- 

model_list <- list(
  m_temp = m_temp,
  m_vv = m_vv,
  m_sla = m_sla,
  m_la = m_la,
  m_nee = m_nee,
  m_reco = m_reco,
  m_gpp = m_gpp
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
  
  #observed vs predicted
  
  tmp_op <- data.frame(observed = dt[[response_name]], 
                       country = dt$country,
                       predicted_fixed_only = predict(m, newdata = dt %>% mutate(site = NA)),
                       predicted = predict(m),
                      response = response_name, 
                      model_name = m_name)
      
  dt_op <- rbind(dt_op, tmp_op)            
  
  print(paste0(m_name, " done (",i, "/", length(model_list), ")"))
}


#Plot observed vs predicted ------------------

p_op_fm <- dt_op %>% 
  mutate(clean_response = case_when(
    response == "temperature_gpp_anomaly_country" ~ "Local Temperature", 
    response == "height_x_cover_anomaly_country" ~ "'Biomass'", 
    response == "sla_anomaly_country" ~ "SLA", 
    response == "leaf_area_anomaly_country" ~ "Leaf Area", 
    response == "nee_anomaly_country" ~ "NEE", 
    response == "reco_anomaly_country" ~ "Reco", 
    response == "gpp_anomaly_country" ~ "GPP", 
  )) %>% 
  ggplot(aes(x = observed, y = predicted)) +
  geom_point(alpha = 0.3, size = 0.75) +
  geom_abline(linetype = "dashed", color = "red") +
  facet_wrap(~clean_response, ncol = 4) +
  theme_bw() +
  labs(title = "a)") +
  ylim(min(dt_op$observed), max(dt_op$observed))
p_op_fm

p_op_fo <- dt_op %>% 
  mutate(clean_response = case_when(
    response == "temperature_gpp_anomaly_country" ~ "Local Temperature", 
    response == "height_x_cover_anomaly_country" ~ "'Biomass'", 
    response == "sla_anomaly_country" ~ "SLA", 
    response == "leaf_area_anomaly_country" ~ "Leaf Area", 
    response == "nee_anomaly_country" ~ "NEE", 
    response == "reco_anomaly_country" ~ "Reco", 
    response == "gpp_anomaly_country" ~ "GPP", 
  )) %>% 
  ggplot(aes(x = observed, y = predicted_fixed_only)) +
  geom_point(alpha = 0.3, size = 0.75) +
  geom_abline(linetype = "dashed", color = "red") +
  facet_wrap(~clean_response, ncol = 4) +
  theme_bw() +
  labs(y = "predicted (fixed effects only)", title = "b)") +
  ylim(min(dt_op$observed), max(dt_op$observed))
p_op_fo

p_op <- grid.arrange(p_op_fm, p_op_fo)
ggsave(plot = p_op, "builds/plots/supplement/observed_vs_predicted.png", dpi = 600, height = 8, width = 7.25)


# Plot predictions ---------

dt_predictions <- dt_pred %>% 
  left_join(dt_results[, c("model_name", "term", "response", "p.value")]) %>% 
  rename(ci_ub = conf.high, 
         ci_lb = conf.low) %>% 
  mutate(sig_yn = ifelse(p.value >= 0.05, "no", "yes"), 
         clean_term =  case_when(
           .default = term,
           grepl("temperature_", term) ~ "Local Temperature",
           term == "mat_anomaly_country" ~ "MAT",
           term == "height_x_cover_anomaly_country"      ~ "'Biomass'",
           term == "morph_traits_pc1_anomaly_country"          ~ "Morph. Traits PC1",
           term == "morph_traits_pc2_anomaly_country"          ~ "Morph. Traits PC2",
           term == "all_traits_pc1_anomaly_country"          ~ "All Traits PC1",
           term == "all_traits_pc2_anomaly_country"          ~ "All Traits PC2",
           term == "chem_traits_pc1_anomaly_country"          ~ "Chem. Traits PC1",
           term == "chem_traits_pc2_anomaly_country"          ~ "Chem. Traits PC2",
           term == "sla_anomaly_country"                 ~ "SLA",
           term == "leaf_area_anomaly_country"           ~ "Leaf Area"), 
         clean_response = case_when(
           grepl("nee", response) ~ "NEE", 
           grepl("reco", response) ~ "Reco", 
           grepl("gpp", response) ~ "GPP", 
           
         ))


dt_long <- dt %>% 
  pivot_longer(cols = c(unique(dt_predictions$term)), 
               values_to = "var_value", names_to = "term") %>% 
  mutate(
         clean_term =  case_when(
           .default = term,
           grepl("temperature_", term) ~ "Local Temperature",
           term == "mat_anomaly_country" ~ "MAT",
           term == "height_x_cover_anomaly_country"      ~ "'Biomass'",
           term == "morph_traits_pc1_anomaly_country"          ~ "Morph. Traits PC1",
           term == "morph_traits_pc2_anomaly_country"          ~ "Morph. Traits PC2",
           term == "all_traits_pc1_anomaly_country"          ~ "All Traits PC1",
           term == "all_traits_pc2_anomaly_country"          ~ "All Traits PC2",
           term == "chem_traits_pc1_anomaly_country"          ~ "Chem. Traits PC1",
           term == "chem_traits_pc2_anomaly_country"          ~ "Chem. Traits PC2",
           term == "sla_anomaly_country"                 ~ "SLA",
           term == "leaf_area_anomaly_country"           ~ "Leaf Area"))

theme_pred <-   theme(legend.position = "none", 
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
# NEE 

#regular
dt_nee_flux <- dt_predictions %>%
  filter(response == "nee_anomaly_country")
p_nee <- dt_nee_flux %>% 
  ggplot() +
  geom_point(data = dt_long %>% filter(term %in% c(unique(dt_nee_flux$term))) %>% 
               mutate(country = fct_reorder(country, lat)),
             aes(x = var_value, y = nee_anomaly_country, color = country), alpha = 0.75) +
  geom_ribbon(aes(x = var_value, ymin = ci_lb, ymax = ci_ub), alpha = 0.25) +
  scale_color_met_d(name = "Archambault") +
  scale_fill_met_d(name = "Archambault") +
  geom_line(aes(x = var_value, y = predicted, linetype = sig_yn), linewidth = 1.01) +
  scale_linetype_manual(values = c(no = "dashed", yes = "solid")) +
  facet_wrap(~clean_term, scales = "free_x", ncol = 5) +
  labs(y = "NEE") +
  theme_bw() + 
  theme_pred
p_nee

dt_reco_flux <- dt_predictions %>%
  filter(response == "reco_anomaly_country")
p_reco <- dt_reco_flux %>% 
  ggplot() +
  geom_point(data = dt_long %>% filter(term %in% c(unique(dt_reco_flux$term))) %>% 
               mutate(country = fct_reorder(country, lat)),
             aes(x = var_value, y = reco_anomaly_country, color = country), alpha = 0.75) +
  geom_ribbon(aes(x = var_value, ymin = ci_lb, ymax = ci_ub), alpha = 0.25) +
  scale_color_met_d(name = "Archambault") +
  scale_fill_met_d(name = "Archambault") +
  geom_line(aes(x = var_value, y = predicted, linetype = sig_yn), linewidth = 1.01) +
  scale_linetype_manual(values = c(no = "dashed", yes = "solid")) +
  facet_wrap(~clean_term, scales = "free_x", ncol = 5) +
  labs(y = "Reco") +
  theme_bw() +
  theme_pred
p_reco

dt_gpp_flux <- dt_predictions %>%
  filter(response == "gpp_anomaly_country")
p_gpp <- dt_gpp_flux %>% 
  ggplot() +
  geom_point(data = dt_long %>% filter(term %in% c(unique(dt_gpp_flux$term))) %>% 
               mutate(country = fct_reorder(country, lat)),
             aes(x = var_value, y = gpp_anomaly_country, color = country), alpha = 0.75) +
  geom_ribbon(aes(x = var_value, ymin = ci_lb, ymax = ci_ub), alpha = 0.25) +
  scale_color_met_d(name = "Archambault") +
  scale_fill_met_d(name = "Archambault") +
  geom_line(aes(x = var_value, y = predicted, linetype = sig_yn), linewidth = 1.01) +
  scale_linetype_manual(values = c(no = "dashed", yes = "solid")) +
  facet_wrap(~clean_term, scales = "free_x", ncol = 5) +
  theme_bw() +
  labs(y = "GPP") +
  theme_pred 
p_gpp

p_pred_flux <- grid.arrange(p_gpp, p_nee, p_reco, ncol = 1, heights = c(1, 1, 1))
ggsave(plot = p_pred_flux, "builds/plots/flux_predictions_.png", dpi = 600, height = 7, width = 8)

# plot relationships of other responses 

dt_temp_flux <- dt_predictions %>%
  filter(response == "temperature_gpp_anomaly_country")
p_temp <- dt_temp_flux %>% 
  ggplot() +
  geom_point(data = dt %>% 
               mutate(country = fct_reorder(country, lat)),
             aes(x = mat_anomaly_country, y = temperature_gpp_anomaly_country, color = country), alpha = 0.75) +
  geom_ribbon(aes(x = var_value, ymin = ci_lb, ymax = ci_ub), alpha = 0.25) +
  scale_color_met_d(name = "Archambault") +
  scale_fill_met_d(name = "Archambault") +
  geom_line(aes(x = var_value, y = predicted, linetype = sig_yn), linewidth = 1.01) +
  scale_linetype_manual(values = c(no = "dashed", yes = "solid")) +
  theme_bw() +
  labs( x= "MAT", y = "Local Temperature") +
  theme_pred
p_temp


dt_vv_flux <- dt_predictions %>%
  filter(response == "height_x_cover_anomaly_country")
p_vv <- dt_vv_flux %>% 
  ggplot() +
  geom_point(data = dt %>% 
               mutate(country = fct_reorder(country, lat)),
             aes(x = mat_anomaly_country, y = height_x_cover_anomaly_country, color = country), alpha = 0.75) +
  geom_ribbon(aes(x = var_value, ymin = ci_lb, ymax = ci_ub), alpha = 0.25) +
  scale_color_met_d(name = "Archambault") +
  scale_fill_met_d(name = "Archambault") +
  geom_line(aes(x = var_value, y = predicted, linetype = sig_yn), linewidth = 1.01) +
  scale_linetype_manual(values = c(no = "dashed", yes = "solid")) +
  labs( x= "MAT", y = "'Biomass'") +
  theme_bw() +
  theme_pred
p_vv


dt_sla_flux <- dt_predictions %>%
  filter(response == "sla_anomaly_country")
p_sla <- dt_sla_flux %>% 
  ggplot() +
  geom_point(data = dt %>% 
               mutate(country = fct_reorder(country, lat)),
             aes(x = mat_anomaly_country, y = sla_anomaly_country, color = country), alpha = 0.75) +
  geom_ribbon(aes(x = var_value, ymin = ci_lb, ymax = ci_ub), alpha = 0.25) +
  scale_color_met_d(name = "Archambault") +
  scale_fill_met_d(name = "Archambault") +
  geom_line(aes(x = var_value, y = predicted, linetype = sig_yn), linewidth = 1.01) +
  scale_linetype_manual(values = c(no = "dashed", yes = "solid")) +
  labs( x= "MAT", y = "SLA") +
  theme_bw() +
  theme_pred
p_sla


dt_la_flux <- dt_predictions %>%
  filter(response == "leaf_area_anomaly_country")
p_la <- dt_la_flux %>% 
  ggplot() +
  geom_point(data = dt %>% 
               mutate(country = fct_reorder(country, lat)),
             aes(x = mat_anomaly_country, y = leaf_area_anomaly_country, color = country), alpha = 0.75) +
  geom_ribbon(aes(x = var_value, ymin = ci_lb, ymax = ci_ub), alpha = 0.25) +
  scale_color_met_d(name = "Archambault") +
  scale_fill_met_d(name = "Archambault") +
  geom_line(aes(x = var_value, y = predicted, linetype = sig_yn), linewidth = 1.01) +
  scale_linetype_manual(values = c(no = "dashed", yes = "solid")) +
  labs( x= "MAT", y = "Leaf Area") +
  theme_bw() +
  theme_pred
p_la


p_1 <- grid.arrange(p_temp, p_vv, p_sla, p_la, ncol = 4)


p_pred <- grid.arrange(p_1, p_pred_flux, heights = c(1, 3))  
ggsave(plot = p_pred, "builds/plots/model_output_bivariate_rels.png", dpi = 600, height = 10, width = 11)
  
  
  