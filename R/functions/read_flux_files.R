## function to read LiCOR 75000 files
## Adapted from Pekka Niittyinen 

read_flux_files <- function (path, photo = "photo", resp = "resp", ambient = "a", recursive = TRUE) 
{
  if (missing(path)) {
    readline("Please set the working directory to the folder\n that contains the LiCOR files to be analyzed.\n Do so with the upcoming prompt. Note that you must choose a file in the folder that\n you want to set as the working directory. \nPlease press 'return' to continue.")
    
    if(recursive == TRUE){
      files <- list.files(dirname(file.choose()), full.names = TRUE, recursive = TRUE)}else{
        files <- list.files(dirname(file.choose()), full.names = TRUE, recursive = FALSE)
      }
  }
  else {
    if(recursive == TRUE){
      files <- list.files(path, full.names = TRUE, recursive = TRUE)}else{
        files <- list.files(path, full.names = TRUE, recursive = FALSE)
      }
  }
  photo.names <- grep(paste0("[^", resp, "].txt$"), grep(paste0("[^_", 
                                                                ambient, "]\\.txt$"), files, value = TRUE), value = TRUE)
  ambient_pattern <- paste0(ambient, ".txt$")
  ambient_names <- grep(ambient_pattern, files, value = TRUE)
  respiration_pattern <- paste0(resp, ".txt$")
  resp.names <- grep(respiration_pattern, files, value = TRUE)
  if (length(photo.names) == 0) {
    message("No matching photo files found.")
  }
  if (length(resp.names) == 0) {
    message("No matching resp files found.")
  }
  if (length(ambient_names) == 0) {
    message("No matching ambient files found. First second of photo files will be used.")
  }
  licor_files <- list(photo_names = photo.names, ambient_names = ambient_names, 
                      resp_names = resp.names)
  invisible(licor_files)
}
