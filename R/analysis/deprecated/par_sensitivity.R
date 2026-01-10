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
    country, site, plot_id, treatment, gradient,
    
    # fluxes
    nee, reco, gpp,
    
    # environmental 
    elevation, map, mat, par,
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
    elevation_anomaly_country,
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


nee_base_formula <- "nee_anomaly_country ~ temperature_nee_anomaly_country +
                  height_x_cover_anomaly_country +
                  sla_anomaly_country + 
                  leaf_area_anomaly_country +
                  elevation_anomaly_country"

# NEE
m_nee_no_par <- glmmTMB(
  formula = nee_anomaly_country ~ 
    temperature_nee_anomaly_country +
    height_x_cover_anomaly_country +
    sla_anomaly_country + 
    leaf_area_anomaly_country +
    elevation_anomaly_country + (1 | site),
  data = dt[!is.na(dt$par_anomaly_country)],
  na.action = na.omit
)
summary(m_nee_no_par); check_collinearity(m_nee_no_par); r.squaredGLMM(m_nee_no_par); AICc(m_nee_no_par)

m_nee_with_par <- glmmTMB(
  formula = nee_anomaly_country ~ 
    temperature_nee_anomaly_country +
    height_x_cover_anomaly_country +
    sla_anomaly_country + 
    leaf_area_anomaly_country +
    elevation_anomaly_country +
    par_anomaly_country + (1 | site),
  data = dt[!is.na(dt$par_anomaly_country)],
  na.action = na.omit
)
summary(m_nee_with_par); check_collinearity(m_nee_with_par); r.squaredGLMM(m_nee_with_par); AICc(m_nee_with_par)

#GPP
m_gpp_no_par <- glmmTMB(
  formula = gpp_anomaly_country ~ 
    temperature_gpp_anomaly_country +
    height_x_cover_anomaly_country +
    sla_anomaly_country + 
    leaf_area_anomaly_country +
    elevation_anomaly_country + (1 | site),
  data = dt[!is.na(dt$par_anomaly_country)],
  na.action = na.omit
)
summary(m_gpp_no_par); check_collinearity(m_gpp_no_par); r.squaredGLMM(m_gpp_no_par); AICc(m_gpp_no_par)

m_gpp_with_par <- glmmTMB(
  formula = gpp_anomaly_country ~ 
    temperature_gpp_anomaly_country +
    height_x_cover_anomaly_country +
    sla_anomaly_country + 
    leaf_area_anomaly_country +
    elevation_anomaly_country +
    par_anomaly_country + (1 | site),
  data = dt[!is.na(dt$par_anomaly_country)],
  na.action = na.omit
)
summary(m_gpp_with_par); check_collinearity(m_gpp_with_par); r.squaredGLMM(m_gpp_with_par); AICc(m_gpp_with_par)

#Reco
m_reco_no_par <- glmmTMB(
  formula = reco_anomaly_country ~ 
    temperature_reco_anomaly_country +
    height_x_cover_anomaly_country +
    sla_anomaly_country + 
    leaf_area_anomaly_country +
    elevation_anomaly_country + (1 | site),
  data = dt[!is.na(dt$par_anomaly_country)],
  na.action = na.omit
)
summary(m_reco_no_par); check_collinearity(m_reco_no_par); r.squaredGLMM(m_reco_no_par); AICc(m_reco_no_par)

m_reco_with_par <- glmmTMB(
  formula = reco_anomaly_country ~ 
    temperature_reco_anomaly_country +
    height_x_cover_anomaly_country +
    sla_anomaly_country + 
    leaf_area_anomaly_country +
    elevation_anomaly_country +
    par_anomaly_country + (1 | site),
  data = dt[!is.na(dt$par_anomaly_country)],
  na.action = na.omit
)
summary(m_reco_with_par); check_collinearity(m_reco_with_par); r.squaredGLMM(m_reco_with_par); AICc(m_reco_with_par)



### Get predictions for alternative hypotheses vs climate --------

model_list <- list(
  m_reco_no_par = m_reco_no_par,
  m_reco_with_par = m_reco_with_par,
  
  m_nee_no_par = m_nee_no_par,
  m_nee_with_par = m_nee_with_par,
  
  m_gpp_no_par = m_gpp_no_par,
  m_gpp_with_par = m_gpp_with_par
)

dt_est <- data.frame()

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
  
  dt_est <- rbind(dt_est, tmp_est)
  
  var_names <- tidy_m %>%
    dplyr::select(term) %>%
    filter(!grepl("ntercept", term) & term != "sd__Observation") %>%
    filter(!grepl("\\:", term))
  
  print(paste0(m_name, " done (",i, "/", length(model_list), ")"))
}


dt_clean <- dt_est %>% 
  mutate(
    clean_response = case_when(
      response == "nee_anomaly_country"  ~ "NEE",
      response == "reco_anomaly_country" ~ "Reco",
      response == "gpp_anomaly_country"  ~ "GPP",
      TRUE ~ response
    ),
    
    clean_response = ifelse(grepl("no_par", model_name), 
                            paste0(clean_response, "\n(No PAR)"), 
                            paste0(clean_response, "\n(With PAR)")),
    
    clean_term =  case_when(
      .default = term,
      grepl("temperature_", term) ~ "Inst. Temperature",
      term == "(Intercept)" ~ "Intercept", 
      term == "elevation_anomaly_country" ~ "Elevation",
      term == "height_x_cover_anomaly_country"      ~ "'Biomass'",
      term == "morph_traits_pc1_anomaly_country"          ~ "Morph. Traits PC1",
      term == "morph_traits_pc2_anomaly_country"          ~ "Morph. Traits PC2",
      term == "all_traits_pc1_anomaly_country"          ~ "All Traits PC1",
      term == "all_traits_pc2_anomaly_country"          ~ "All Traits PC2",
      term == "chem_traits_pc1_anomaly_country"          ~ "Chem. Traits PC1",
      term == "chem_traits_pc2_anomaly_country"          ~ "Chem. Traits PC2",
      term == "sla_anomaly_country"                 ~ "SLA",
      term == "leaf_area_anomaly_country"           ~ "Leaf Area", 
      term == "par_anomaly_country"                 ~ "PAR",
    ), 
    significance = case_when(
      .default = "Non significant", 
      ci_lb > 0 ~ "Significantly positive", 
      ci_ub < 0 ~ "Significantly negative", 
    )) %>% 
  mutate(
    clean_term = factor(clean_term,
                        levels = c("PAR",
                                   "'Biomass'",
                                   "Leaf Area",
                                   "SLA",
                                   "Chem. Traits PC2", "Chem. Traits PC1",
                                   "All Traits PC2", "All Traits PC1",
                                   "Morph. Traits PC2", "Morph. Traits PC1",
                                   "Inst. Temperature",
                                   "Elevation")),
    response = factor(response,
                      levels = c("gpp_anomaly_country",
                                 "nee_anomaly_country",
                                 "reco_anomaly_country",
                                 "temperature_gpp_anomaly_country",
                                 "temperature_reco_anomaly_country",
                                 "temperature_nee_anomaly_country",
                                 "morph_traits_pc1_anomaly_country",
                                 "morph_traits_pc2_anomaly_country",
                                 "all_traits_pc1_anomaly_country",
                                 "all_traits_pc2_anomaly_country",
                                 "chem_traits_pc1_anomaly_country",
                                 "chem_traits_pc2_anomaly_country",
                                 "sla_anomaly_country",
                                 "leaf_area_anomaly_country",
                                 "height_x_cover_anomaly_country")
    )
  )


resp_clean_comb <- dt_clean %>% 
  select(response, clean_response) %>%
  unique() %>% 
  arrange(response) 

dt_clean$clean_response <- factor(dt_clean$clean_response, levels = c(
  resp_clean_comb$clean_response
))

# Plots -------------
scico::scico(9, palette = 'bam')
#"#001959" "#818231" "#F9CCF9"

theme_est <- theme(legend.position = "none", 
                     legend.box="vertical",
                     plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                     panel.grid = element_line(color = "snow"), 
                     #axis.title.x = element_blank(), 
                     axis.text = element_text(size = 12), 
                     axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
                     panel.border = element_rect(color = NA), 
                     panel.background = element_rect(fill = "snow"), 
                     strip.text.x = element_text(size = 12), 
                     strip.text.y = element_text(size = 12, face = "bold"), 
                     strip.background = element_rect(fill = "linen", color = "linen") )

# Tier 1 ----------------
p <- dt_clean %>% 
  filter(!grepl("ntercept", term)) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linewidth = 1, linetype = "dotted", color = "grey25") +
  geom_pointrange(aes(y = clean_term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = significance),
                  linewidth = 1.2, size = 0.9, alpha = 0.9) +
  scale_color_manual(values = c("Non significant" = "grey", "Significantly positive" = "#4B802E","Significantly negative" = "#A4428B")) +
  facet_wrap(~clean_response, scales = "free_x", ncol = 6)  +
  theme_bw() +
  labs(x = "Estimate", y = "",) +
  theme_est
p

ggsave(plot = p,
       "builds/plots/supplement/par_flux_estimates.png",
       dpi = 600, 
       height = 3, width = 9)


library(ggridges)
dt_raw %>% 
  filter(!is.na(par_anomaly_country)) %>% 
  ggplot() +
  geom_density_ridges(aes(x = par, y = gradient))
