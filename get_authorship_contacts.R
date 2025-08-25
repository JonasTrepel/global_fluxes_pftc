library(data.table)
library(tidyverse)


dt <- fread("/Users/jonas/Downloads/PFTCs participants - PFTC participants.csv")
unique(dt$Role)
names(dt)

dt_inc <- dt %>% 
  mutate(ask = case_when(
    .default = "no", 
    Role %in% c("Instructor", "group leader", "student / group leader", 
                "Instructor / outreach") ~ "yes", 
    grepl("(4)", `PFTC2 (China) 2016`) ~ "yes", 
    grepl("(1)", `PFTC2 (China) 2016`) ~ "yes", 
    grepl("(4)", `PFTC3 (Peru) 2018`) ~ "maybe",
    grepl("(1)", `PFTC3 (Peru) 2018`) ~ "maybe",
    grepl("(4)", `PFTC4 (Svalbard) 2018`) ~ "yes",
    grepl("(3)", `PFTC4 (Svalbard) 2018`) ~ "yes",
    grepl("(4)", `PFTC5 (Peru) 2020`) ~ "maybe", 
    grepl("(1)", `PFTC5 (Peru) 2020`) ~ "maybe",
    grepl("(3)", `PFTC7 (South Africa) 2023`) ~ "yes",
    grepl("(1)", `PFTC7 (South Africa) 2023`) ~ "yes"
  ))
table(dt_inc$ask)
