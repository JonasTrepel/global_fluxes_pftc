#####################
library(tidyverse)
library(googlesheets4)
library(stringr)

# Read the author table from Drive
ss <- "https://docs.google.com/spreadsheets/d/1q7ny-lVg2VgSbdDFWq4dXM1YcIw_v_alZYbH-iGqlOI/edit?gid=0#gid=0"
d <- read_sheet(ss, col_types = "c") 

# Extract components of author contributions
conce <- as.character(na.omit(d$conceptualization))
datacur <- as.character(na.omit(d$data_curation))
anal <- as.character(na.omit(d$analysis))
fund <- as.character(na.omit(d$funding))
inve <- as.character(na.omit(d$investigation))
meth <- as.character(na.omit(d$methodology))
adm <- as.character(na.omit(d$administration))
reso <- as.character(na.omit(d$resources))
supe <- as.character(na.omit(d$supervision))
visu <- as.character(na.omit(d$visualization))
writ <- as.character(na.omit(d$writing))
edi <- as.character(na.omit(d$editing))

# Combine to a single string
paste(paste0(paste(conce[-length(conce)], collapse = ", "), ", and ", conce[length(conce)], " conceptualized the study."),
      paste0("All authors participated in data collection."),
      paste0(paste(datacur[-length(datacur)], collapse = ", "), ", and ", datacur[length(datacur)], " curated the data."),
      paste0(paste(meth[-length(meth)], collapse = ", "), ", and ", meth[length(meth)], " created the methodology."),
      paste0(paste(anal[-length(anal)], collapse = ", "), " and ", anal[length(anal)], " conducted the analyses."),
      paste0("JT did the visualisations."),
      paste0(paste(writ[-length(writ)], collapse = ", "), ", and ", writ[length(writ)], " wrote the first draft."),
      paste0("All authors edited the manuscript."),
      paste0(paste(fund[-length(fund)], collapse = ", "), ", and ", fund[length(fund)], " acquired funding and provided resources."),
      paste0(paste(adm[-length(adm)], collapse = ", "), ", and ", adm[length(adm)], " were responsible for administration."),
      paste0(paste(supe[-length(supe)], collapse = ", "), ", and ", supe[length(supe)], " conducted supervision."), sep = " "
)

###################
# Author list & affiliations

# Transform to long format
d2 <- d %>% 
  select(given_name, family_name, affiliation_1:affiliation_3) %>% 
  pivot_longer(cols = affiliation_1:affiliation_3, values_to = "Affiliation") %>%
  drop_na(Affiliation)

# Create a running number for the affiliations in order of the author list
d2 <- d2 %>% 
  mutate(Affiliation2 = factor(Affiliation, levels = unique(Affiliation))) %>%
  group_by(Affiliation2) %>%
  mutate(Aff_number = cur_group_id()) %>%
  ungroup() %>% 
  select(-Affiliation2)

# Helper function for superscript digits
superscript <- function(n) {
  superscript_digits <- c("⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹")
  digits <- as.integer(unlist(strsplit(as.character(n), "")))
  paste0(superscript_digits[digits + 1], collapse = "")
}

# Combine the author names and affiliation indicators
lapply(1:nrow(d), function(x){
  
  author_name <- d %>% slice(x) %>% select(given_name, family_name) %>% as.list() %>% paste(collapse = " ")
  author_affiliations <- d2 %>% 
    filter(given_name == d$given_name[[x]],
           family_name == d$family_name[[x]])
  
  paste0(author_name, " ",  paste(lapply(author_affiliations$Aff_number, superscript), collapse = ","))
  
}) %>% 
  paste0(collapse = ", ")
# The output can be just copied from the console to any text editor

# This will construct the list of affiliations
lapply(unique(d2$Aff_number), function(x){
  
  paste0(superscript(x), (d2 %>% filter(Aff_number == x) %>% pull(Affiliation) %>% unique))
  
  
}) %>% 
  paste(collapse = ", ")

out <- lapply(unique(d2$Aff_number), function(x){
  paste0(
    superscript(x),
    (d2 %>% filter(Aff_number == x) %>% pull(Affiliation) %>% unique)
  )
}) %>% paste(collapse = ", ")

writeLines(out, pipe("pbcopy"))



mails = unique(d$email) %>% paste(collapse = ", ")
writeLines(mails, pipe("pbcopy"))
