#get SA temperature and PAR sorted
# PAR ------------------

source("R/functions/read_flux_files.R")
#source("/functions/fixFileNames.R")
source("R/functions/get_flux_df.R")
source("R/functions/read_and_add_par.R")
source("R/functions/fix_file_names.R")



## path with flux files: 
fix_file_names(path = "data/raw_data/south_africa/raw_li7500_measurements/")

list.files("data/raw_data/south_africa/raw_li7500_measurements/", recursive = TRUE)


flux_files <- read_flux_files(path = "data/raw_data/south_africa/raw_li7500_measurements/",
                           photo = "photo", ## specify the patterns in filenames to categorize fluxes
                           resp = "resp", 
                           ambient = "a",
                           recursive = TRUE)


## get Flux metadata (South Africa specific unfortunately)

flux_meta <- tibble(filename = unlist(flux_files),
                    file = basename(filename)) %>%
  mutate(site = unlist(lapply(file, function(x) str_split(x, "_")[[1]][1])),
         elevation = unlist(lapply(file, function(x) str_split(x, "_")[[1]][2])),
         aspect = unlist(lapply(file, function(x) str_split(x, "_")[[1]][3])),
         plot = unlist(lapply(file, function(x) str_split(x, "_")[[1]][4])),
         day_night = unlist(lapply(file, function(x) str_split(x, "_")[[1]][5])),
         measurement = unlist(lapply(file, function(x) gsub(".txt","",tail(str_split(x, "_")[[1]],1)))),
         redo = grepl("redo", file, ignore.case = T), 
         measurement = gsub("\\(1\\)", "", measurement), 
         plotID = paste0(elevation, aspect, plot))

flux_meta


## read fluxes into a dataframe 
dt_flux_raw <- get_flux_df(fluxfiles = flux_files,
                       skip = 3, #default
                       device = "LI7500" #default
) %>% 
  mutate(air_temperature = as.numeric(air_temperature)) %>% 
  left_join(flux_meta)


## Add PAR -------------

par_files <- list.files("data/raw_data/south_africa/raw_PAR", full.names = T)
par_skip = 7


dt_flux <- read_and_add_par(par_files = par_files, #default = NULL
                        flux_df = dt_flux_raw, #default = NULL 
                        par_time_zone = "Etc/GMT+8", #default = UTC
                        flux_time_zone = "Africa/Johannesburg", #default = UTC
                        flux_df_time_col = "date_time", #default
                        par_skip = 7 #default
)


## get metadata/field records 

dt_par_tmp <- dt_flux %>% 
  mutate(plot_id = paste0("SA_", elevation, aspect, plot), 
         date = date(date_time),
         flux_type = case_when(
           .default = "ambient", 
           grepl("_resp", filename) ~ "reco", 
           grepl("_photo", filename) ~ "nee", 
         )) %>% 
  group_by(date, file) %>% 
  summarize(temperature = mean(air_temperature, na.rm = T), 
            par = mean(par, na.rm = T)) %>%
  unique()


fwrite(dt_par_tmp, "data/raw_data/south_africa/pftc7_flux_par_and_tmp.csv")


#### PAR manually -------------------

par_files <- list.files("data/raw_data/south_africa/raw_PAR", full.names = T)
par_skip = 7
#Site 1 --------------------

dt_flux_s1 <- dt_flux_raw %>%
  filter(site == 1) %>%
  filter(!grepl("old", filename)) %>% 
  filter(day_night == "day")

par_files_s1 <- par_files[3]

par_list <- list()

if(!is.null(par_files_s1)){
  par_list[[1]] <- suppressMessages(suppressWarnings(read_delim(par_files_s1[1], delim="\t", skip = par_skip)))
  if(length(par_files_s1) > 1){
    for(f in c(2:length(par_files_s1))){
      par_list[[f]] <-  suppressMessages(suppressWarnings(read_delim(par_files_s1[f], delim="\t", skip = par_skip)))
    }
    par_comb_s1 <- bind_rows(par_list)
  }
  else{
    par_comb_s1 <- par_list[[1]]
  }
} 


dt_flux_s1$POSIXct <- as.POSIXct(dt_flux_s1$date_time, 
                                 format="%Y-%m-%d %H:%M:%S", 
                                 tz = "Africa/Johannesburg")

range(dt_flux_s1$POSIXct, na.rm = T)

par_comb_s1$POSIXct_uc <- as.POSIXct(paste(par_comb_s1$Date, par_comb_s1$Time),
                                     format="%Y-%m-%d %H:%M:%S", 
                                     tz = "Etc/GMT+8")

range(par_comb_s1$POSIXct_uc, na.rm = T)

common_tz = "Africa/Johannesburg"

par_comb_s1$POSIXct <- with_tz(par_comb_s1$POSIXct_uc, tzone = common_tz)

range(par_comb_s1$POSIXct, na.rm = T)
range(par_comb_s1$POSIXct_uc, na.rm = T)
range(dt_flux_s1$POSIXct, na.rm = T)

dt_flux_s1 <- inner_join(par_comb_s1, dt_flux_s1, by="POSIXct") %>% 
  rename(par = INPUT1) %>% 
  dplyr::select(c(par, filename, date_time))

summary(dt_flux_s1)


# Site 2 -------

#no measurements, sadly 

#Site 3 --------------------

# first attempts 
dt_flux_s3.1 <- dt_flux_raw %>%
  filter(site == 3) %>%
  filter(!grepl("old", filename)) %>% 
  filter(day_night == "day") %>% 
  mutate(date = date(date_time)) %>% 
  filter(date == "2023-12-09")

unique(dt_flux_s3.1$date)
par_files_s3.1 <- par_files[1:2]

par_list <- list()

if(!is.null(par_files_s3.1)){
  par_list[[1]] <- suppressMessages(suppressWarnings(read_delim(par_files_s3.1[1], delim="\t", skip = par_skip)))
  if(length(par_files_s3.1) > 1){
    for(f in c(2:length(par_files_s3.1))){
      par_list[[f]] <-  suppressMessages(suppressWarnings(read_delim(par_files_s3.1[f], delim="\t", skip = par_skip)))
    }
    par_comb_s3.1 <- bind_rows(par_list)
  }
  else{
    par_comb_s3.1 <- par_list[[1]]
  }
} 

# Redo
dt_flux_s3.1$POSIXct <- as.POSIXct(dt_flux_s3.1$date_time, 
                                   format="%Y-%m-%d %H:%M:%S", 
                                   tz = "Africa/Johannesburg")

range(dt_flux_s3.1$POSIXct, na.rm = T)

par_comb_s3.1$POSIXct_uc <- as.POSIXct(paste(par_comb_s3.1$Date, par_comb_s3.1$Time),
                                       format="%Y-%m-%d %H:%M:%S", 
                                       tz = "Etc/GMT+8")

range(par_comb_s3.1$POSIXct_uc, na.rm = T)

common_tz = "Africa/Johannesburg"

par_comb_s3.1$POSIXct <- with_tz(par_comb_s3.1$POSIXct_uc, tzone = common_tz)

range(par_comb_s3.1$POSIXct, na.rm = T)
range(par_comb_s3.1$POSIXct_uc, na.rm = T)
range(dt_flux_s3.1$POSIXct, na.rm = T)

dt_flux_s3.1 <- inner_join(par_comb_s3.1, dt_flux_s3.1, by="POSIXct") %>% 
  rename(par = INPUT1) %>% 
  dplyr::select(c(par, filename, date_time))

summary(dt_flux_s3.1)

# Redo 

dt_flux_s3.2 <- dt_flux_raw %>%
  filter(site == 3) %>%
  filter(!grepl("old", filename)) %>% 
  filter(day_night == "day") %>% 
  mutate(date = date(date_time)) %>% 
  filter(date == "2023-12-14")

unique(dt_flux_s3.2$date)
par_files_s3.2 <- par_files[4]

par_list <- list()

if(!is.null(par_files_s3.2)){
  par_list[[1]] <- suppressMessages(suppressWarnings(read_delim(par_files_s3.2[1], delim="\t", skip = par_skip)))
  if(length(par_files_s3.2) > 1){
    for(f in c(2:length(par_files_s3.2))){
      par_list[[f]] <-  suppressMessages(suppressWarnings(read_delim(par_files_s3.2[f], delim="\t", skip = par_skip)))
    }
    par_comb_s3.2 <- bind_rows(par_list)
  }
  else{
    par_comb_s3.2 <- par_list[[1]]
  }
} 


dt_flux_s3.2$POSIXct <- as.POSIXct(dt_flux_s3.2$date_time, 
                                   format="%Y-%m-%d %H:%M:%S", 
                                   tz = "Africa/Johannesburg")

range(dt_flux_s3.2$POSIXct, na.rm = T)

par_comb_s3.2$POSIXct_uc <- as.POSIXct(paste(par_comb_s3.2$Date, par_comb_s3.2$Time),
                                       format="%Y-%m-%d %H:%M:%S", 
                                       tz = "Etc/GMT+8")

range(par_comb_s3.2$POSIXct_uc, na.rm = T)

common_tz = "Africa/Johannesburg"

par_comb_s3.2$POSIXct <- with_tz(par_comb_s3.2$POSIXct_uc, tzone = common_tz)

range(par_comb_s3.2$POSIXct, na.rm = T)
range(par_comb_s3.2$POSIXct_uc, na.rm = T)
range(dt_flux_s3.2$POSIXct, na.rm = T)

dt_flux_s3.2 <- inner_join(par_comb_s3.2, dt_flux_s3.2, by="POSIXct") %>% 
  rename(par = INPUT1) %>% 
  dplyr::select(c(par, filename, date_time))

summary(dt_flux_s3.2)

# Site 4 ------------------
dt_flux_s4 <- dt_flux_raw %>%
  filter(site == 4) %>%
  filter(!grepl("old", filename)) %>% 
  filter(day_night == "day") %>% 
  mutate(date = date(date_time)) 

unique(dt_flux_s4$date)
par_files_s4 <- par_files[5]

par_list <- list()

if(!is.null(par_files_s4)){
  par_list[[1]] <- suppressMessages(suppressWarnings(read_delim(par_files_s4[1], delim="\t", skip = par_skip)))
  if(length(par_files_s4) > 1){
    for(f in c(2:length(par_files_s4))){
      par_list[[f]] <-  suppressMessages(suppressWarnings(read_delim(par_files_s4[f], delim="\t", skip = par_skip)))
    }
    par_comb_s4 <- bind_rows(par_list)
  }
  else{
    par_comb_s4 <- par_list[[1]]
  }
} 


dt_flux_s4$POSIXct <- as.POSIXct(dt_flux_s4$date_time, 
                                 format="%Y-%m-%d %H:%M:%S", 
                                 tz = "Africa/Johannesburg")

range(dt_flux_s4$POSIXct, na.rm = T)

par_comb_s4$POSIXct_uc <- as.POSIXct(paste(par_comb_s4$Date, par_comb_s4$Time),
                                     format="%Y-%m-%d %H:%M:%S", 
                                     tz = "Etc/GMT+8")

range(par_comb_s4$POSIXct_uc, na.rm = T)

common_tz = "Africa/Johannesburg"

par_comb_s4$POSIXct <- with_tz(par_comb_s4$POSIXct_uc, tzone = common_tz)

range(par_comb_s4$POSIXct, na.rm = T)
range(par_comb_s4$POSIXct_uc, na.rm = T)
range(dt_flux_s4$POSIXct, na.rm = T)

dt_flux_s4 <- inner_join(par_comb_s4, dt_flux_s4, by="POSIXct") %>% 
  rename(par = INPUT1) %>% 
  dplyr::select(c(par, filename, date_time))

summary(dt_flux_s4)

# Site 5 ------------------
dt_flux_s5 <- dt_flux_raw %>%
  filter(site == 5) %>%
  filter(!grepl("old", filename)) %>% 
  filter(day_night == "day") %>% 
  mutate(date = date(date_time)) 

unique(dt_flux_s5$date)
par_files_s5 <- par_files[6:7]

par_list <- list()

if(!is.null(par_files_s5)){
  par_list[[1]] <- suppressMessages(suppressWarnings(read_delim(par_files_s5[1], delim="\t", skip = par_skip)))
  if(length(par_files_s5) > 1){
    for(f in c(2:length(par_files_s5))){
      par_list[[f]] <-  suppressMessages(suppressWarnings(read_delim(par_files_s5[f], delim="\t", skip = par_skip)))
    }
    par_comb_s5 <- bind_rows(par_list)
  }
  else{
    par_comb_s5 <- par_list[[1]]
  }
} 


dt_flux_s5$POSIXct <- as.POSIXct(dt_flux_s5$date_time, 
                                 format="%Y-%m-%d %H:%M:%S", 
                                 tz = "Africa/Johannesburg")

range(dt_flux_s5$POSIXct, na.rm = T)

par_comb_s5$POSIXct_uc <- as.POSIXct(paste(par_comb_s5$Date, par_comb_s5$Time),
                                     format="%Y-%m-%d %H:%M:%S", 
                                     tz = "Etc/GMT+8")

range(par_comb_s5$POSIXct_uc, na.rm = T)

common_tz = "Africa/Johannesburg"

par_comb_s5$POSIXct <- with_tz(par_comb_s5$POSIXct_uc, tzone = common_tz)

range(par_comb_s5$POSIXct, na.rm = T)
range(par_comb_s5$POSIXct_uc, na.rm = T)
range(dt_flux_s5$POSIXct, na.rm = T)

dt_flux_s5 <- inner_join(par_comb_s5, dt_flux_s5, by="POSIXct") %>% 
  rename(par = INPUT1) %>% 
  dplyr::select(c(par, filename, date_time))

summary(dt_flux_s5)



# Temperature ------------------

# Extract temperatures from the LI7500 data
# devtools::install_github("PaulESantos/co2fluxtent")
library(co2fluxtent)
library(tidyverse)
library(data.table)

# Read in all the ambient files ----
# Remember to fix the file names if you've just downloaded from OSF
source("R/functions/fix_file_names.R")
source("R/functions/test_flux_files.R")
fix_file_names(path = "data/raw_data/south_africa/raw_li7500_measurements/")

list.files("data/raw_data/south_africa/raw_li7500_measurements/", recursive = TRUE)
# Look for flux files in a folder
licor_files <- Map(c, co2fluxtent::read_files("data/raw_data/south_africa/raw_li7500_measurements/LI7500_Site 1"),
                   co2fluxtent::read_files("data/raw_data/south_africa/raw_li7500_measurements/LI7500_Site 2"),
                   co2fluxtent::read_files("data/raw_data/south_africa/raw_li7500_measurements/LI7500_Site 3"),
                   co2fluxtent::read_files("data/raw_data/south_africa/raw_li7500_measurements/LI7500_Site 4"),
                   co2fluxtent::read_files("data/raw_data/south_africa/raw_li7500_measurements/LI7500_Site 5"))

## clean file names
# Check if the files are ok
licor_files <- test_flux_files(licor_files, skip = 3, min_rows = 50) ##removed three files

#All
#licor_files_ambient <- licor_files[["ambient_names"]]
files <- unlist(licor_files) |>
  # Fixed missing file names thanks to this post https://stackoverflow.com/questions/69357657/purrrmap-dfr-gives-number-of-list-element-as-id-argument-not-value-of-list-e
  purrr::set_names()

#Read them in
temp <- map_df(set_names(files), function(file) {
  file %>%
    set_names() %>%
    map_df(~ read_delim(file = file, skip = 3, delim = "\t")) #important! reading in American format
}, .id = "File")

# Gather site, plot etc. information from the filenames
t_ambient = temp |>
  # extract metadata
  mutate(File = basename(File),
         File = str_remove(File, ".txt")
  ) |>
  # Separate into relevant info
  separate(File, into = c("siteID", "elevation", "aspect", "plot", "day.night"), remove = FALSE) |>
  # Get just the relevant data
  select(File:day.night, Time:`Sequence Number`, `Temperature (C)`, `CO2 Signal Strength`) |>
  # Filter to just data with a decent signal strength
  filter(`CO2 Signal Strength` >= 90) %>% 
  rename(temperature = `Temperature (C)`, 
         date = Date) %>% 
  mutate(plot_id = paste0("SA_", elevation, aspect, plot), 
         flux_type = case_when(
           .default = "ambient", 
           grepl("_resp", File) ~ "reco", 
           grepl("_photo", File) ~ "nee", 
         )) %>% 
  group_by(day.night, date, plot_id, flux_type) %>% 
  summarize(temperature = mean(temperature, na.rm = T))
  

fwrite(t_ambient, "data/raw_data/south_africa/LI7500_plot_flux_temperature.csv")
