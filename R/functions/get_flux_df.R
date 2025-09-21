get_flux_df <- function(fluxfiles, skip = 3, device = "LI7500"){  
  # Function to read multiple flux data files and process them into a combined data table.
  # Arguments:
  # fluxfiles: A list of file paths to flux data files.
  # skip: The number of lines to skip at the start of each file (default is 3 which makes at least sense for the Licor7500).
  # device: The device used to record the flux data (default is "Licor7500"), it could be useful to include more devices later on.
  
  if(device == "LI7500"){ 
    
    # Convert the list of flux files to a vector and remove names
    fluxFileV <- unlist(fluxfiles) %>% unname()  
    
    # Create an empty data table to store the combined flux data
    fluxDT <- data.table::data.table()  
    
    # Loop through each unique file in the fluxFileV vector
    for(file in unique(fluxFileV)){  
      
      # Read the flux data file, skipping the first few lines as specified by `skip`.
      # The delimiter is tab-separated ("\t").
      rawFlux <- suppressMessages(readr::read_delim(file, skip = skip, delim = "\t"))  
      
      # Process the raw flux data and select specific columns for analysis
      tempFlux <- rawFlux %>%
        dplyr::select(Time, Date, `CO2 (umol/mol)`, `H2O (mmol/mol)`, `Temperature (C)`, `Pressure (kPa)`, `CO2 Signal Strength`) %>%
        
        # Rename the columns for easier use later
        rename(
          conc_co2 = `CO2 (umol/mol)`,
          conc_h2o = `H2O (mmol/mol)`, 
          air_temperature = `Temperature (C)`, 
          pressure_kpa = `Pressure (kPa)`,  
          signal_strength = `CO2 Signal Strength`,  
          licor_time = Time, 
          licor_date = Date   
        ) %>% 
        
        # Perform data transformations
        mutate(
          licor_time = gsub("\\:000", "", licor_time),  # Clean time column by removing trailing ":000" - no idea why that's there 
          licor_date = lubridate::ymd(licor_date),
          licor_time = lubridate::hms(licor_time),
          date_time = lubridate::ymd_hms(paste0(licor_date, " ", licor_time)),  
          start_time = min(date_time), 
          filename = file, 
          conc_co2 = as.numeric(conc_co2),
          conc_h2o = as.numeric(conc_h2o)
        ) %>% 
        dplyr::select(-c(licor_date, licor_time))  
      
      # Append the processed data (tempFlux) to the combined data table (fluxDT)
      fluxDT <- rbind(fluxDT, tempFlux)  
      
      # Print a message indicating that the file has been processed
      print(paste0("File done: ", file))  
    }  
    
  }else{
    print(paste0("Sorry, can do only Licor7500 at this point :("))
  }  
  
  return(fluxDT)  
}
