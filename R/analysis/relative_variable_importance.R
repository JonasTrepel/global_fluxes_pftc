library(glmm.hp)
library(tidyverse)
library(data.table)
library(glmmTMB)





## Prepare models -------

fm_t4_nee_1 <- glmmTMB(nee_anomaly_country ~
                         temperature_nee_anomaly_country +
                         elevation_anomaly_country +
                         height_x_cover_anomaly_country +
                         morph_traits_pc1_anomaly_country +
                         morph_traits_pc2_anomaly_country +
                         ( 1 | site), na.action = na.omit,
                       data = m_dat)
summary(fm_t4_nee_1)

fm_t4_nee_2 <- glmmTMB(nee_anomaly_country ~
                         temperature_nee_anomaly_country +
                         elevation_anomaly_country +
                         height_x_cover_anomaly_country +
                         leaf_area_anomaly_country +
                         sla_anomaly_country +
                         ( 1 | site), na.action = na.omit,
                       data = m_dat)
summary(fm_t4_nee_2)

fm_t4_nee_3 <- glmmTMB(nee_anomaly_country ~
                         temperature_nee_anomaly_country +
                         elevation_anomaly_country +
                         height_x_cover_anomaly_country +
                         all_traits_pc1_anomaly_country +
                         all_traits_pc2_anomaly_country +
                         ( 1 | site), na.action = na.omit,
                       data = m_dat)
summary(fm_t4_nee_3)



#Reco

fm_t4_reco_1 <- glmmTMB(reco_anomaly_country ~
                         temperature_reco_anomaly_country +
                         elevation_anomaly_country +
                         height_x_cover_anomaly_country +
                         morph_traits_pc1_anomaly_country +
                         morph_traits_pc2_anomaly_country +
                         ( 1 | site), na.action = na.omit,
                       data = m_dat)
summary(fm_t4_reco_1)

fm_t4_reco_2 <- glmmTMB(reco_anomaly_country ~
                         temperature_reco_anomaly_country +
                         elevation_anomaly_country +
                         height_x_cover_anomaly_country +
                         leaf_area_anomaly_country +
                         sla_anomaly_country +
                         ( 1 | site), na.action = na.omit,
                       data = m_dat)
summary(fm_t4_reco_2)

fm_t4_reco_3 <- glmmTMB(reco_anomaly_country ~
                         temperature_reco_anomaly_country +
                         elevation_anomaly_country +
                         height_x_cover_anomaly_country +
                         all_traits_pc1_anomaly_country +
                         all_traits_pc2_anomaly_country +
                         ( 1 | site), na.action = na.omit,
                       data = m_dat)
summary(fm_t4_reco_3)

#Gpp
fm_t4_gpp_1 <- glmmTMB(gpp_anomaly_country ~
                         temperature_gpp_anomaly_country +
                         elevation_anomaly_country +
                         height_x_cover_anomaly_country +
                         morph_traits_pc1_anomaly_country +
                         morph_traits_pc2_anomaly_country +
                         ( 1 | site), na.action = na.omit,
                       data = m_dat)
summary(fm_t4_gpp_1)

fm_t4_gpp_2 <- glmmTMB(gpp_anomaly_country ~
                         temperature_gpp_anomaly_country +
                         elevation_anomaly_country +
                         height_x_cover_anomaly_country +
                         leaf_area_anomaly_country +
                         sla_anomaly_country +
                         ( 1 | site), na.action = na.omit,
                       data = m_dat)
summary(fm_t4_gpp_2)

fm_t4_gpp_3 <- glmmTMB(gpp_anomaly_country ~
                         temperature_gpp_anomaly_country +
                         elevation_anomaly_country +
                         height_x_cover_anomaly_country +
                         all_traits_pc1_anomaly_country +
                         all_traits_pc2_anomaly_country +
                         ( 1 | site), na.action = na.omit,
                       data = m_dat)
summary(fm_t4_gpp_3)

#extract AIC etc 
model_list_fm <- list(
  
  fm_t4_nee_1 = fm_t4_nee_1,
  fm_t4_nee_2 = fm_t4_nee_2,
  fm_t4_nee_3 = fm_t4_nee_3,
  
  fm_t4_reco_1 = fm_t4_reco_1,
  fm_t4_reco_2 = fm_t4_reco_2,
  fm_t4_reco_3 = fm_t4_reco_3,
  
  fm_t4_gpp_1 = fm_t4_gpp_1,
  fm_t4_gpp_2 = fm_t4_gpp_2,
  fm_t4_gpp_3 = fm_t4_gpp_3
  
)

dt_res_fm <- data.frame()

for(i in 1:length(model_list_fm)){
  
  
  m <- model_list_fm[[i]]
  
  m_name <- names(model_list_fm)[i]
  
  v_par <- glmm.hp(m)
  
  tmp_dt <- v_par$hierarchical.partitioning %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "predictor") %>% 
    mutate( 
    model_name = m_name,
    rsq_m = as.numeric(r.squaredGLMM(m)[1]), 
    rsq_c = as.numeric(r.squaredGLMM(m)[2]))
  
  dt_res_fm <- rbind(dt_res_fm, tmp_dt)
  
  print(paste0(i, " done"))
}

dt_vp <- dt_res_fm %>% 
  mutate(
    predictor_tier = case_when(
      .default = "not_applicable",
      grepl("_1", model_name) ~ "morph_traits_pca",
      grepl("_2", model_name) ~ "traits_separately", 
      grepl("_3", model_name) ~ "all_traits_pca"), 
    response = case_when(
      grepl("reco", model_name) ~ "Reco",
      grepl("nee", model_name) ~ "NEE",
      grepl("gpp", model_name)~ "GPP"), 
  clean_term =  case_when(
    .default = predictor,
    grepl("temperature_", predictor) ~ "Local Temperature",
    predictor == "elevation_anomaly_country" ~ "Elevation",
    predictor == "height_x_cover_anomaly_country"      ~ "'Biomass'",
    predictor == "morph_traits_pc1_anomaly_country"          ~ "Morph. Traits PC1",
    predictor == "morph_traits_pc2_anomaly_country"          ~ "Morph. Traits PC2",
    predictor == "all_traits_pc1_anomaly_country"          ~ "All Traits PC1",
    predictor == "all_traits_pc2_anomaly_country"          ~ "All Traits PC2",
    predictor == "sla_anomaly_country"                 ~ "SLA",
    predictor == "leaf_area_anomaly_country"           ~ "Leaf Area"), 
  label = paste0(response, "\n(R2m = ", round(rsq_m, 2), ")")) %>% 
  rename(var_imp = `I.perc(%)`) 



p_vp_ts <- dt_vp %>% 
  filter(predictor_tier == "traits_separately") %>% 
  ggplot() +
  geom_col(aes(x = var_imp, y = clean_term)) +
  facet_wrap(~label) +
  labs(x = "Percentage of Marginal R2", y = "") +
  theme_minimal() +
  theme(legend.position = "none",
        legend.box="vertical",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.grid = element_line(color = "seashell"), 
        #axis.title.x = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.text.x = element_text(size = 10, angle = 0, hjust = 1), 
        panel.background = element_rect(fill = "snow2", color = NA), 
        strip.text.x = element_text(size = 12), 
        strip.text.y = element_text(size = 12, face = "bold"), 
        strip.background = element_rect(fill = "seashell", color = "seashell") )

p_vp_ts

ggsave(plot = p_vp_ts, "builds/plots/variable_importance_traits_separately.png", dpi = 600, height = 3, width = 7)

p_vp_mt <- dt_vp %>% 
  filter(predictor_tier == "morph_traits_pca") %>% 
  ggplot() +
  geom_col(aes(x = var_imp, y = clean_term)) +
  facet_wrap(~label) +
  labs(x = "Percentage of Marginal R2", y = "") +
  theme_minimal() +
  theme(legend.position = "none",
        legend.box="vertical",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.grid = element_line(color = "seashell"), 
        #axis.title.x = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.text.x = element_text(size = 10, angle = 0, hjust = 1), 
        panel.background = element_rect(fill = "snow2", color = NA), 
        strip.text.x = element_text(size = 12), 
        strip.text.y = element_text(size = 12, face = "bold"), 
        strip.background = element_rect(fill = "seashell", color = "seashell") )

p_vp_mt


p_vp_at <- dt_vp %>% 
  filter(predictor_tier == "all_traits_pca") %>% 
  ggplot() +
  geom_col(aes(x = var_imp, y = clean_term)) +
  facet_wrap(~label) +
  labs(x = "Percentage of Marginal R2", y = "") +
  theme_minimal() +
  theme(legend.position = "none",
        legend.box="vertical",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.grid = element_line(color = "seashell"), 
        #axis.title.x = element_blank(), 
        axis.text = element_text(size = 10), 
        axis.text.x = element_text(size = 10, angle = 0, hjust = 1), 
        panel.background = element_rect(fill = "snow2", color = NA), 
        strip.text.x = element_text(size = 12), 
        strip.text.y = element_text(size = 12, face = "bold"), 
        strip.background = element_rect(fill = "seashell", color = "seashell") )

p_vp_at
