
require(dplyr)
require(lubridate)

# par_files <- par_files
# flux_df <- flux_dfRaw
# par_time_zone <- "America/Mexico_City"
# flux_time_zone = "Africa/Johannesburg"

# adapted from a funciton from Michael Mustri

read_and_add_par <- function(par_files = NULL, 
                          flux_df = NULL,
                          par_time_zone = "UTC",
                          flux_time_zone = "UTC",
                          flux_df_time_col = "DateTime",
                          par_skip = 7){
  
  
  if(is.null(par_files)){
    print("Please provide PAR file names")
  }
  
  if(is.null(flux_df)){
    print("Please provide flux data")
  }
  
  parDT <- list()
  
  if(!is.null(par_files)){
    parDT[[1]] <- suppressMessages(suppressWarnings(read_delim(par_files[1], delim="\t", skip = par_skip)))
    if(length(par_files) > 1){
      for(f in c(2:length(par_files))){
        parDT[[f]] <-  suppressMessages(suppressWarnings(read_delim(par_files[f], delim="\t", skip = par_skip)))
      }
      parComb <- bind_rows(parDT)
    }
    else{
      parComb <- parDT[[1]]
    }
  }
  else{
    parComb = NULL
  }
  
  
  if(!is.null(parComb)){
    flux_df$POSIXct <- as.POSIXct(flux_df[[flux_df_time_col]], 
                                 format="%Y-%m-%d %H:%M:%S", 
                                 tz = flux_time_zone)
    

    parComb$POSIXct_uc <- as.POSIXct(paste(parComb$Date, parComb$Time),
                                     format="%Y-%m-%d %H:%M:%S", 
                                     tz = par_time_zone)
    
    range(parComb$POSIXct_uc, na.rm = T)
    range(flux_df$POSIXct, na.rm = T)
  

    commonTZ = flux_time_zone
    
    parComb$POSIXct <- with_tz(parComb$POSIXct_uc, tzone = commonTZ)
    
    range(parComb$POSIXct, na.rm = T)
    range(parComb$POSIXct_uc, na.rm = T)
    range(flux_df$POSIXct, na.rm = T)
    
    dt_par_fluxRaw <- inner_join(parComb, flux_df, by="POSIXct") %>% 
      rename(par = INPUT1) %>% 
      dplyr::select(c(par, filename, date_time))
    
    print(paste0("Succesfully joined PAR and fluxes, resturning flux_df with PAR column"))
  }
  else{
    
    flux_df$POSIXct <- as.POSIXct(flux_df[[flux_df_time_col]], 
                                 format="%Y-%m-%d %H:%M:%S", 
                                 tz = flux_time_zone)
    dt_par_fluxRaw <- flux_df %>% 
      dplyr::select(c(Filename, DateTime)) %>% mutate(PAR = NA)
    
    print(paste0("Something went wrong with the PAR files - returning flux_df without PAR column"))
  }
  
  dt_par_flux <- suppressMessages(flux_df %>% left_join(dt_par_fluxRaw))
  
  if(sum(is.na(dt_par_flux$par)) > 0 &  sum(is.na(dt_par_flux$par)) < nrow(dt_par_flux)){
    print(paste0("Uhoh, some PAR values are NAs (", round(sum(is.na(dt_par_flux$par))/nrow(dt_par_flux)*100, 1), "%)"))
  }
  return(dt_par_flux)
}
