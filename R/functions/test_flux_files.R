test_flux_files <- function(fluxfiles, skip = 3, min_rows = 20){

  fluxfiles$photo_names <- unlist(lapply(fluxfiles$photo_names, function(filename){
    e <- try(suppressMessages(suppressWarnings(input <- read_delim(filename, skip = skip, delim = "\t"))), silent = TRUE)

    if(class(e)[[1]] == "try-error"){
      print(paste(filename, "error in reading the file, likely not a LI-7500 data file"))
      return(NULL)
    } else {
      if(!"CO2 (mmol/m^3)" %in% names(input)){
        print(paste(filename, "is likely not a LI-7500 data file"))
        return(NULL)
      } else {
        if(nrow(input) < min_rows){
          print(paste(filename, "has too few observations"))
          return(NULL)
        } else {
          return(filename)
        }
      }
    }
  }))

  fluxfiles$ambient_names <- unlist(lapply(fluxfiles$ambient_names, function(filename){
    e <- try(suppressMessages(suppressWarnings(input <- read_delim(filename, skip = skip, delim = "\t"))), silent = TRUE)

    if(class(e)[[1]] == "try-error"){
      print(paste(filename, "error in reading the file, likely not a LI-7500 data file"))
      return(NULL)
    } else {
      if(!"CO2 (mmol/m^3)" %in% names(input)){
        print(paste(filename, "is likely not a LI-7500 data file"))
        return(NULL)
      } else {
        if(nrow(input) < min_rows){
          print(paste(filename, "has too few observations"))
          return(NULL)
        } else {
          return(filename)
        }
      }
    }
  }))

  fluxfiles$resp_names <- unlist(lapply(fluxfiles$resp_names, function(filename){
    e <- try(suppressMessages(suppressWarnings(input <- read_delim(filename, skip = skip, delim = "\t"))), silent = TRUE)

    if(class(e)[[1]] == "try-error"){
      print(paste(filename, "error in reading the file, likely not a LI-7500 data file"))
      return(NULL)
    } else {
      if(!"CO2 (mmol/m^3)" %in% names(input)){
        print(paste(filename, "is likely not a LI-7500 data file"))
        return(NULL)
      } else {
        if(nrow(input) < min_rows){
          print(paste(filename, "has too few observations"))
          return(NULL)
        } else {
          return(filename)
        }
      }
    }
  }))

  return(fluxfiles)
}

filter_flagged <- function(fluxfiles, flux_df){
  x1 <- list(photo_names = fluxfiles$photo_names, resp_names = fluxfiles$resp_names) %>%
    lapply(., function(x){ x[x %in% licor_nee$filename[licor_nee$flagged]]})

  return(c(x1, fluxfiles["ambient_names"]))

}
