library(sf)
library(tidyverse)
library(rgee)

source("R/environmental_data/FUNCTIONS_extract_from_GEE.R")
source("R/environmental_data/FUNCTIONS_clean_rs_data.R")

# Intialize the Earth Engine with rgee
ee_Initialize(drive = T)

tdir <- "C://MyTemp/RS"
if(!dir.exists(tdir)){
  dir.create(tdir)
}

############################################

p <- st_read("data/coordinates.gpkg") %>% 
  rename(geometry = geom) %>% 
  rowid_to_column("SITE_ID")

# Exctract Landsat
task_list <- Landsat_export_points(pixel_coords_sf = p,
                                   sample_id_from = "SITE_ID",
                                   start_doy = 1,
                                   end_doy = 366,
                                   start_date = "2010-01-01",
                                   end_date = as.character(today()))
# Monitor export progress, waiting for last export to have finished
map(task_list, ee_monitoring, max_attempts = 500, quiet = FALSE, task_time  = 20)

# Copy exported file(s) to tempfolder in R using ee_drive_to_local()
gls <- googledrive::drive_find(pattern = "landsat_export.*.csv")
if(nrow(gls) != length(task_list)){stop("Weird number of files in Google Drive")}
for(ii in 1:nrow(gls)){
  googledrive::drive_download(gls$id[ii], path = paste0(tdir, "/", gls$name[ii]), overwrite = T)
}
temp_files <- paste0(tdir, "/", gls$name)

ld <- lapply(temp_files, function(x){
  print(x)
  d <- read_csv(x)
  return(d)
}) %>% bind_rows()

ee_clean_container(name = "landsat_export")
unlink(unlist(temp_files))

# Clean landsat data
ld <- clean_landsat(ld, sample_id_from = "SITE_ID")

ld %>% write_csv("data/environmental_data/Landsat_data.csv")

################################################################
# SENTINEL 2 data

# Exctract Sentinel
task_list <- Sentinel_export_points(pixel_coords_sf = p,
                                   sample_id_from = "SITE_ID",
                                   start_doy = 1,
                                   end_doy = 366,
                                   start_date = "2010-01-01",
                                   end_date = as.character(today()), 
                                   scale = 10)
# Monitor export progress, waiting for last export to have finished
map(task_list, ee_monitoring, max_attempts = 500, quiet = FALSE, task_time  = 20)

# Copy exported file(s) to tempfolder in R using ee_drive_to_local()
gls <- googledrive::drive_find(pattern = "sentinel_export.*.csv")
if(nrow(gls) != length(task_list)){stop("Weird number of files in Google Drive")}
for(ii in 1:nrow(gls)){
  googledrive::drive_download(gls$id[ii], path = paste0(tdir, "/", gls$name[ii]), overwrite = T)
}
temp_files <- paste0(tdir, "/", gls$name)

s2d <- lapply(temp_files, function(x){
  print(x)
  d <- read_csv(x)
  return(d)
}) %>% bind_rows()

ee_clean_container(name = "sentinel_export")
unlink(unlist(temp_files))

# Clean landsat data
s2d <- clean_sentinel(s2d, sample_id_from = "SITE_ID")

s2d %>% write_csv("data/environmental_data/Sentinel2_data.csv")

#######################################################
# Combine

library(modelbased)
library(mgcv)

p <- st_read("data/coordinates.gpkg") %>% 
  rename(geometry = geom) %>% 
  rowid_to_column("SITE_ID")

ld <- read_csv("data/environmental_data/Landsat_data.csv")
s2d <- read_csv("data/environmental_data/Sentinel2_data.csv")

d <- bind_rows(ld, s2d) %>% 
  arrange(SITE_ID, date) %>% 
  filter(date > "2017-01-01") %>% 
  mutate(NDVI = ifelse(NDVI < 0, 0, NDVI))

unique(d$SITE_ID)
id <- 138
d %>% 
  filter(SITE_ID == id) %>% 
  ggplot(aes(x = date, y = NDVI, color = satellite, group = satellite)) +
  geom_point()
d %>% 
  filter(SITE_ID == id) %>% 
  ggplot(aes(x = yday(date), y = NDVI, color = satellite, group = satellite)) +
  geom_point()



d %>% 
  filter(SITE_ID == id) %>% 
  gam(NDVI ~ s(as.numeric(date)) + s(yday(date)), data = .) -> gmod
summary(gmod)

# data <- d
# id_col <- "SITE_ID"
# ndvi_col <- "NDVI"
# date_col <- "date"
# satellite_col <- "satellite"

calc_phenology <- function(data, id_col, ndvi_col, date_col, satellite_col){
  # i <- 1000
  data <- data %>% 
    rename(ID = any_of(id_col),
           ndvi = any_of(ndvi_col),
           date = any_of(date_col),
           satellite = any_of(satellite_col)) %>% 
    mutate(date = as.Date(date)) %>% 
    mutate(yd = yday(date),
           jd = as.numeric(date)) %>% 
    select(ID, date, yd, jd, ndvi, satellite) %>% 
    mutate(satellite = factor(satellite))
  
  df_all <- lapply(unique(data$ID), function(x){
    # x <- 6
    df <- data %>% filter(ID == x)
    
    m1 <- gam(ndvi ~ s(yd, bs="cc", k = 5) + satellite, data = df,
        knots = list(yd = c(1, 366)),
        family = betar(link = "logit"))
    
    df$rsds <- df$ndvi - predict(m1, df, "response")
    
    df <- df %>%
      filter(abs(rsds) < 0.2)
    
    m2 <- gam(ndvi ~ s(yd, bs="cc", k = 5) + satellite, data = df,
             knots = list(yd = c(1, 366)),
             family = betar(link = "logit"))
    
    es <- estimate_slopes(m2, trend = "yd", at = "yd", regrid = "response",
                          length = diff(range(m2$model$yd))+1)
    
    # plot(es)
    dft <- full_join(tibble(ID = rep(x, times = 365),
                            yd = 1:365,
                            ndvi = predict(m2, tibble(yd = 1:365, satellite = "LANDSAT_8"), type = "response") %>% as.numeric),
                     es %>% tibble() %>% 
                       select(yd, Coefficient, SE, CI_low, CI_high, p),
                     by = "yd") %>% 
      mutate(ID2 = x)
    
    # gg <- dft %>%
    #   ggplot(aes(y = ndvi, x = yd)) +
    #   geom_line() +
    #   ggtitle(x)
    # print(gg)
    return(dft)
    
  }) %>% 
    bind_rows()
  
  return(df_all)
}

res <- calc_phenology(data = d, 
                       id_col = "SITE_ID", ndvi_col = "NDVI", date_col = "date", satellite_col = "satellite")

res <- left_join(res, p %>% st_drop_geometry() %>% rename(ID = SITE_ID))

res <- res %>% 
  mutate(area = ifelse(area == "south_africa", "South_Africa", area))

res %>% 
  filter(ID %in% sample(1:222, 20)) %>% 
  ggplot(aes(y = ndvi, x = yd, color = area)) +
  geom_line() + 
  facet_wrap(vars(ID))

maxndvi <- res %>% 
  group_by(area, site, elev, treatment, plot) %>% 
  summarise(ndvi_max = max(ndvi),
            ndvi_min = min(ndvi),
            ndvi_avg = mean(ndvi)) %>% 
  ungroup

f <- read_csv("data/processed_data/preliminary_data/prelim_fluxes.csv") %>% 
  select(datetime, date, elevation, treatment, plot, tier, site, latitude, longitude, plot_id) %>% 
  mutate(area = str_replace_all(tier, "[:digit:]", "")) %>% 
  mutate(area = gsub('.{0,1}$', '', area)) %>% 
  mutate(area = ifelse(area == "Colorado", "USA", area))
unique(f$area)
table(res$area)

res2 <- res %>% 
  select(area, site, elev, treatment, plot, yd, ndvi, Coefficient) %>% 
  rename(ndvi_change = Coefficient)

# Combine china

china <- full_join(f %>% filter(area == "China"),
                   res2 %>% filter(area == "China") %>% 
                     mutate(site = substr(site, 1, 1)) %>% 
                     select(area, site, yd, ndvi, ndvi_change)) %>% 
  distinct() %>% 
  filter(yday(date) == yd) %>% 
  select(-yd)

# Combine Peru

peru <- full_join(f %>% filter(area == "Peru") %>% 
                    mutate(plot = as.character(extract_numeric(plot_id))),
                  res2 %>% filter(area == "Peru") %>% 
                    mutate(treatment = tolower(treatment)) %>% 
                    select(site, treatment, plot, yd, ndvi, ndvi_change)) %>% 
  distinct() %>% 
  filter(yday(date) == yd) %>% 
  select(-yd)

# Combine Svalbard

svalbard <- full_join(f %>% filter(area == "Svalbard") %>% 
                    mutate(site2 = substr(plot, 1, 3)),
                  res2 %>% filter(area == "Svalbard") %>% 
                    rename(site2 = site) %>% 
                    select(site2, yd, ndvi, ndvi_change)) %>% 
  distinct() %>% 
  select(-site2) %>% 
  filter(yday(date) == yd) %>% 
  select(-yd)

# Combine Norway

norway <- full_join(f %>% filter(area == "Norway"),
                  res2 %>% filter(area == "Norway") %>% 
                    rename(plot_id = plot) %>% 
                    select(plot_id, yd, ndvi, ndvi_change)) %>% 
  distinct() %>% 
  filter(yday(date) == yd) %>% 
  select(-yd)

# Combine USA

usa <- full_join(f %>% filter(area == "USA"),
                  res2 %>% filter(area == "USA") %>% 
                    mutate(site = tolower(site)) %>% 
                    select(site, yd, ndvi, ndvi_change)) %>% 
  distinct() %>% 
  filter(yday(date) == yd) %>% 
  select(-yd)

# Combine South_Africa

sa <- full_join(f %>% filter(area == "South_Africa"),
                res2 %>% filter(area == "South_Africa") %>% 
                  mutate(treatment = ifelse(treatment == "W", "west", "east")) %>% 
                  mutate(site = factor(site, levels = c("1","2","3","4","5"), labels = c("2000","2200","2400","2600","2800"))) %>% 
                  mutate(plot_id = paste0("SA_",site,treatment,plot)) %>% 
                  select(plot_id, yd, ndvi, ndvi_change)) %>% 
  distinct() %>% 
  filter(yday(date) == yd) %>% 
  select(-yd)

# Bring together

ndvimom <- bind_rows(china, peru, norway, svalbard, usa, sa)

ndvimom %>% distinct() %>% 
  write_csv("data/environmental_data/NDVI_moment.csv")


#################################################
# Polish

# Combine china

china <- full_join(f %>% filter(area == "China"),
                   maxndvi %>% filter(area == "China") %>% 
                     mutate(site = substr(site, 1, 1)) %>% 
                     select(area, site, ndvi_max, ndvi_min, ndvi_avg)) %>% 
  distinct()

# Combine Peru

peru <- full_join(f %>% filter(area == "Peru") %>% 
                    mutate(plot = as.character(extract_numeric(plot_id))),
                  maxndvi %>% filter(area == "Peru") %>% 
                    mutate(treatment = tolower(treatment)) %>% 
                    select(site, treatment, plot, ndvi_max, ndvi_min, ndvi_avg)) %>% 
  distinct()

# Combine Svalbard

svalbard <- full_join(f %>% filter(area == "Svalbard") %>% 
                        mutate(site2 = substr(plot, 1, 3)),
                      maxndvi %>% filter(area == "Svalbard") %>% 
                        rename(site2 = site) %>% 
                        select(site2, ndvi_max, ndvi_min, ndvi_avg)) %>% 
  distinct() %>% 
  select(-site2)

# Combine Norway

norway <- full_join(f %>% filter(area == "Norway"),
                    maxndvi %>% filter(area == "Norway") %>% 
                      rename(plot_id = plot) %>% 
                      select(plot_id, ndvi_max, ndvi_min, ndvi_avg)) %>% 
  distinct()

# Combine USA

usa <- full_join(f %>% filter(area == "USA"),
                 maxndvi %>% filter(area == "USA") %>% 
                   mutate(site = tolower(site)) %>% 
                   select(site, ndvi_max, ndvi_min, ndvi_avg)) %>% 
  distinct()

# Combine South_Africa

sa <- full_join(f %>% filter(area == "South_Africa"),
                maxndvi %>% filter(area == "South_Africa") %>% 
                  mutate(treatment = ifelse(treatment == "W", "west", "east")) %>% 
                  mutate(site = factor(site, levels = c("1","2","3","4","5"), labels = c("2000","2200","2400","2600","2800"))) %>% 
                  mutate(plot_id = paste0("SA_",site,treatment,plot)) %>% 
                  select(plot_id, ndvi_max, ndvi_min, ndvi_avg)) %>% 
  distinct()

# Bring together

ndvimax <- bind_rows(china, peru, norway, svalbard, usa, sa) %>% 
  select(-datetime, -date)

############################################
# Combine all

ndvimax %>% distinct() %>% 
  write_csv("data/environmental_data/NDVI_means.csv")


