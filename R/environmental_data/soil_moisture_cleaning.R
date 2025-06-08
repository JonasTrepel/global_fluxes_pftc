library(tidyverse)
library(sf)
library(zoo)
library(data.table)

d <- read_csv("data/processed_data/preliminary_data/prelim_fluxes.csv")

# Svalbard

dt <- d %>% filter(startsWith(tier, "Sval")) %>% 
  select(plot_id, date) %>% 
  distinct()

sm <- read_csv("data/raw_data/svalbard/Cflux_SV_ITEX_2018.csv") %>% 
  select(PlotID, Date, starts_with("SM")) %>% 
  distinct() %>% 
  mutate(PlotID = gsub("-","_",paste0("SV_ITEX", PlotID))) %>% 
  group_by(PlotID, Date) %>% 
  rowwise() %>% 
  summarise(soil_moisture = mean(c(SM1,SM2,SM3,SM4))) %>% 
  ungroup() %>% 
  distinct() %>% 
  rename(plot_id = PlotID) %>% 
  mutate(date = dmy(Date)) %>% 
  select(-Date)

dt$plot_id[!dt$plot_id %in% sm$plot_id]

svalbard <- left_join(dt, sm)

# Norway

dt <- d %>% filter(startsWith(tier, "Nor")) %>% 
  select(plot_id, date) %>% 
  distinct()

sm <- read_csv("data/environmental_data/Soil_moisture_raw/PFTC6_microclimate_2022.csv") %>% 
  filter(warming == "A", climate_variable == "soil_moisture") %>% 
  select(datetime, turfID, origSiteID, value) %>% 
  mutate(date = as_date(datetime)) %>% 
  group_by(turfID, origSiteID, date) %>% 
  summarise(soil_moisture = mean(value, na.rm = TRUE)) %>% ungroup() %>% 
  mutate(origSiteID = substr(origSiteID, 1, 3)) %>% 
  mutate(plot_id = paste0("NO_", origSiteID, "_", turfID)) %>% 
  distinct() %>% 
  select(plot_id, date, soil_moisture) %>% 
  mutate(soil_moisture = soil_moisture*100)

dt$plot_id[!dt$plot_id %in% sm$plot_id]

norway <- left_join(dt, sm)

# South Africa

dt <- d %>% filter(startsWith(tier, "South")) %>% 
  select(plot_id, elevation, aspect, plot, date) %>% 
  rename(site = elevation) %>% 
  distinct() %>% 
  mutate(date = as_date(gsub("2024","2023",date)))
unique(dt$plot_id)
sm <- read_csv("data/environmental_data/Soil_moisture_raw/PFTC7_Tomst_Data.csv") %>% 
  select(datetime, site, aspect, plot, moist_vol) %>% 
  mutate(date = as_date(datetime)) %>% 
  group_by(site, aspect, plot, date) %>% 
  summarise(soil_moisture = mean(moist_vol, na.rm = TRUE)) %>% ungroup() %>% 
  mutate(site = case_match(
    site,
    1 ~ 2000,
    2 ~ 2200,
    3 ~ 2400,
    4 ~ 2600,
    5 ~ 2800
  )) %>% 
  mutate(plot_id = paste0("SA_", site, aspect, plot)) %>% 
  distinct()

dt$plot_id[!dt$plot_id %in% sm$plot_id]

south_africa <- left_join(dt, sm %>% mutate(plot = as.character(plot))) %>% 
  arrange(site, aspect, plot) %>% 
  group_by(site, aspect) %>% 
  mutate(soil_moisture = na.approx(soil_moisture, maxgap = 4)) %>% 
  ungroup() %>% 
  select(plot_id, date, soil_moisture)

# China

dt <- d %>% filter(startsWith(tier, "China")) %>% 
  select(plot_id, date) %>% 
  distinct()
unique(dt$plot_id)


# Peru

dt <- d %>% filter(startsWith(tier, "Colorado")) %>% 
  select(site, date) %>% 
  distinct()
unique(dt$site)

us.moist.raw <- fread("data/raw_data/colorado/soil_moisture_data.csv") 

us.moist.raw2 <- us.moist.raw %>% 
  mutate(site = tolower(site_clean),
         date = ymd(date), 
         tier = paste0("Colorado_", year)) %>% 
  filter(tier %in% c("Colorado_2016", "Colorado_2018")) %>% 
  mutate(
    percent_moisture = ifelse(percent_moisture < 0, 0, percent_moisture), 
    percent_moisture = ifelse(percent_moisture > 1, NA, percent_moisture)) %>% 
  filter(site %in% unique(dt$site) & date %in% unique(dt$date)) %>% 
  group_by(site, date) %>% 
  summarize(
    soil_moisture = mean(percent_moisture, na.rm = T)*100) %>% 
  mutate(unit = "percent", 
         date = as.Date(date))


us.moist <- d %>% 
  filter(startsWith(tier, "Colorado")) %>% 
  dplyr::select(plot_id, site, date) %>% 
  unique() %>% 
  left_join(us.moist.raw2) %>% 
  filter(!is.na(soil_moisture)) %>% 
  dplyr::select(plot_id, soil_moisture, date)



# Combine

all <- bind_rows(svalbard, norway, south_africa, us.moist)

write_csv(all, "data/environmental_data/soil_moist_temp.csv")
