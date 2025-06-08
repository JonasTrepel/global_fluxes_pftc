library(tidyverse)
library(data.table)
library(rnaturalearth)
library(sf)
library(mapview)
library(scico)
library(MetBrewer)
library(gridExtra)

dt <- fread("data/processed_data/clean_data/global_fluxes_main_data.csv") %>% 
  select(country, spatial_block, site, plot_id, mmp, mat, elevation, 
         temperature_nee, temperature_reco, temperature_gpp,
         height_x_cover, sla_cm2_g, ldmc, leaf_area_cm2, nee, reco, gpp,
         lat, lon, species_richness
  ) %>% 
  rename(local_temperature = temperature_gpp) %>% 
  filter(complete.cases(.))
table(dt$country)
sf_fluxes <- read_sf("data/processed_data/preliminary_data/prelim_flux_loc.gpkg") %>% 
  rename(plot_id = plot_id) 

dt_sf <- dt %>%
  left_join(sf_fluxes) %>%
  st_as_sf() %>%
  st_transform(crs = 'ESRI:54009') %>% 
  mutate(country = fct_reorder(country, lat)) 

mapview(dt_sf)

world <- rnaturalearth::ne_countries() %>% filter(!name_en == "Antarctica") %>% st_transform(crs = 'ESRI:54009')

world %>% st_transform(crs = 'ESRI:54009')

#c(MetBrewer::met.brewer(name = "Archambault", n = 6))

map <- ggplot() +
  geom_sf(data = world)+
  geom_sf(data = dt_sf, aes(color = country), size = 4) +
  #scale_color_scico_d(palette = "batlow") +
  scale_color_met_d(name = "Archambault") +
  theme_void() +
  theme(legend.position = "none")
map

p_empty <- ggplot() + theme_void()
map_e <- grid.arrange(p_empty, map, widths = c(0.5, 2), ncol = 2)

##### make boxplots ####

units <- c(
  elevation = "m",
  mat = "°C",
  mmp = "mm",
  nee = "µmol/m²/s",
  reco = "µmol/m²/s",
  gpp = "µmol/m²/s", 
  sla_cm2_g = "cm²/g",          
  leaf_area_cm2 = "cm²",       
  local_temperature = "°C", 
  species_richness = "dimensionless",
  height_x_cover = "dimensionless", 
  ldmc = "g"
)

dt_units <- data.table(
  var_name = names(units), 
  unit = units)

dt_long <- dt %>% 
  pivot_longer(cols = c(elevation, mat, mmp, nee, reco, gpp, ldmc, species_richness,
                        sla_cm2_g, leaf_area_cm2, local_temperature, height_x_cover), 
               names_to = "var_name", values_to = "var_value") %>% 
  left_join(dt_units) %>% 
  mutate(plotvar_name = ifelse(!var_name %in% c("height_x_cover", "species_richness"), paste0(var_name, " (", unit, ")"), var_name)) %>% 
  mutate(plotvar_name = gsub("local_temperature", "Local Temperature", plotvar_name),
         plotvar_name = gsub("leaf_area_cm2", "Leaf Area", plotvar_name),
         plotvar_name = gsub("height_x_cover", "'Biomass'", plotvar_name), 
         plotvar_name = gsub("ldmc", "LDMC", plotvar_name), 
         plotvar_name = gsub("sla_cm2_g", "SLA", plotvar_name), 
         plotvar_name = gsub("gpp", "GPP", plotvar_name), 
         plotvar_name = gsub("nee", "NEE", plotvar_name), 
         plotvar_name = gsub("reco", "Reco", plotvar_name), 
         plotvar_name = gsub("elevation", "Elevation", plotvar_name), 
         plotvar_name = gsub("mat", "MAT", plotvar_name), 
         plotvar_name = gsub("mmp", "MMP", plotvar_name), 
         plotvar_name = gsub("species_richness", "Species Richness", plotvar_name), 
         
         )

library(ggridges)
pba1 <- dt_long %>%
  filter(var_name %in% c("nee", "gpp", "reco", "mat", "elevation")) %>% 
  mutate(country = fct_reorder(country, lat), 
         plotvar_name = factor(plotvar_name, levels = c("GPP (µmol/m²/s)",
                                                        "Reco (µmol/m²/s)",
                                                        "NEE (µmol/m²/s)",
                                                        "Elevation (m)", 
                                                        "MAT (°C)"))) %>%
  ggplot() +
  geom_density_ridges(aes(y = country, x = var_value, color = country, fill = country), alpha = 0.75, size = 0.75) +
  scale_color_met_d(name = "Archambault") +
  scale_fill_met_d(name = "Archambault") +
  facet_wrap(~plotvar_name, scales = "free_x", ncol = 5) +
  theme_bw() +
  labs(y = "", x = "Variable Value") +
  theme(legend.position = "none", 
        legend.box="vertical",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #panel.grid = element_line(color = "seashell"), 
        #axis.title.x = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        panel.border = element_rect(color = NA), 
        panel.background = element_rect(fill = "snow2"), 
        strip.text.x = element_text(size = 14), 
        strip.text.y = element_text(size = 14, face = "bold"), 
        strip.background = element_rect(fill = "seashell", color = "seashell") )
pba1



pba3 <- dt_long %>%
  filter(var_name %in% c("sla_cm2_g", "leaf_area_cm2", "height_x_cover", "local_temperature")) %>% 
  mutate(country = fct_reorder(country, lat)) %>%
  ggplot() +
  geom_density_ridges(aes(y = country, x = var_value, color = country, fill = country), alpha = 0.75, size = 0.75) +
  scale_color_met_d(name = "Archambault") +
  scale_fill_met_d(name = "Archambault") +
  facet_wrap(~plotvar_name, scales = "free_x", ncol = 6) +
  theme_bw() +
  labs(y = "", x = "Variable Value") +
  theme(legend.position = "none", 
        legend.box="vertical",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        #panel.grid = element_line(color = "seashell"), 
        #axis.title.x = element_blank(), 
        axis.text = element_text(size = 12), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        panel.border = element_rect(color = NA), 
        panel.background = element_rect(fill = "snow2"), 
        strip.text.x = element_text(size = 14), 
        strip.text.y = element_text(size = 14, face = "bold"), 
        strip.background = element_rect(fill = "seashell", color = "seashell") )
pba3

pba_fus <- grid.arrange(pba1, pba2, ncol = 2)             

comb_plot_a <- grid.arrange(pba1, map_e, pba3, heights = c(1.1, 1.5, 1.2), 
                          padding = unit(0.5, "lines"))
ggsave(plot = comb_plot_a, "builds/plots/map.png", dpi = 600, height = 10, width = 10)

