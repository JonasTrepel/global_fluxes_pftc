library(data.table)
library(tidyverse)


dt <- fread("R/data_cleaning_and_prep/organization/pftc_participants_old.csv") %>% 
  filter(!Participant == "")
unique(dt$Role)
names(dt)

dt_inc <- dt %>% 
  mutate(ask = case_when(
    .default = "no", 
    Role %in% c("Instructor", "group leader", "student / group leader", 
                "Instructor / outreach", "Site PI") ~ "yes", 
    grepl("(4)", `PFTC2 (China) 2016`) ~ "yes", 
    grepl("(1)", `PFTC2 (China) 2016`) ~ "yes", 
    grepl("(4)", `PFTC3 (Peru) 2018`) ~ "maybe",
    grepl("(1)", `PFTC3 (Peru) 2018`) ~ "maybe",
    grepl("(4)", `PFTC4 (Svalbard) 2018`) ~ "yes",
    grepl("(3)", `PFTC4 (Svalbard) 2018`) ~ "yes",
    grepl("(4)", `PFTC5 (Peru) 2020`) ~ "maybe", 
    grepl("(1)", `PFTC5 (Peru) 2020`) ~ "maybe",
    #grepl("student", `PFTC6 (Norway) 2022`) ~ "maybe",
    grepl("(3)", `PFTC7 (South Africa) 2023`) ~ "yes",
    grepl("(1)", `PFTC7 (South Africa) 2023`) ~ "yes", 
  )) 
table(dt_inc$ask)

c(unique(dt_inc[dt_inc$ask == "yes" | dt_inc$ask == "maybe", ]$Email))

c(unique(dt_inc[dt_inc$ask == "yes" | dt_inc$ask == "maybe", ]$Participant))

mails <- paste(c(unique(dt_inc[dt_inc$ask == "yes" | dt_inc$ask == "maybe", ]$Email)), collapse = "; ")
head(mails, all = T)
 ## plus all the norway students 

#original mail went to: 

sent_to = c(
  "pengahui@yahoo.com",
  "akuonanizakeyo@gmail.com",
  "brummera@email.arizona.edu",
  "alexander-vagenes@hotmail.no",
  "anna.gowera@students.wits.ac.za",
  "aud.halbritter@uib.no",
  "agualdo@uni-goettingen.de",
  "b.seaman.e@gmail.com",
  "by9102@gmail.com",
  "bismark.ofosu-bamfo@uenr.edu.gh",
  "benquist@email.arizona.edu",
  "bmaitner@gmail.com",
  "u16112017@tuks.co.za",
  "carmen.vzri@gmail.com",
  "christiansen.ct@gmail.com",
  "christien@sasenvgroup.co.za",
  "claire.ponsac@hotmail.fr",
  "dagmar.egelkraut@uib.no",
  "e.nh94@hotmail.com",
  "ranfei@imde.ac.cn",
  "francisco.navarrorosales@wolfson.ox.ac.uk",
  "greciarivas2399@outlook.com",
  "hanna.lee@ntnu.no",
  "hrosedawson@gmail.com",
  "u19053721@tuks.co.za",
  "imma.oliveras@ouce.ox.ac.uk",
  "ialt@norceresearch.no",
  "sallobravo@gmail.com",
  "j.white2@kew.org",
  "jonas.trepel@gmail.com",
  "jonathan.henn@colorado.edu",
  "j.lee8@westernsydney.edu.au",
  "julia.kemppinen@helsinki.fi",
  "kagoff@ncsu.edu",
  "kari.klanderud@nmbu.no",
  "karopank@gmail.com",
  "lesego0981@gmail.com",
  "liyenne.hagenberg@umu.se",
  "lcaviere@udec.cl",
  "lorahbeth@gmail.com",
  "loreleipatrick@gmail.com",
  "mm2809@cam.ac.uk",
  "marcella.cross@botany.ubc.ca",
  "marcus.spiegel@ouce.ox.ac.uk",
  "mary.linabury@gmail.com",
  "matiss@email.arizona.edu",
  "sullivanmks@gmail.com",
  "mustrim@arizona.edu",
  "nadine.arzt@uib.no",
  "mohamm88@msu.edu",
  "nina.roth@natgeo.su.se",
  "onalennag37@gmail.com",
  "pekka.niittynen@helsinki.fi",
  "p.b.eidesen@ibv.uio.no",
  "peter.leroux@up.ac.za",
  "glaciobotanist@gmail.com",
  "hansda.priya95@gmail.com",
  "ragnhild.gya@uib.no",
  "ragnhild.s.stokka@hotmail.com",
  "becca.harris@colostate.edu",
  "richard.telford@uib.no",
  "ruben.erik.roos@nmbu.no",
  "smduranm@gmail.com",
  "nasirsara41@gmail.com",
  "sean.michaletz@botany.ubc.ca",
  "sehoya@umn.edu",
  "silje.ostman@uib.no",
  "siri.haugum@muho.no",
  "tin.satriawan@u.nus.edu",
  "vbuzzard@email.arizona.edu",
  "vigdis.vandvik@uib.no",
  "wfarfan@gmail.com",
  "sunxiangyang@imde.ac.cn",
  "xiaoxiangzhao94@yahoo.com",
  "ymalhi1@gmail.com",
  "holle.mukhlish@gmail.com",
  "marta.baumane@bio.ku.dk",
  "b.g.oberholzer@outlook.com",
  "sorrelh@umich.edu",
  "correio.marta@gmail.com",
  "sam.ahler@colorado.edu",
  "j.atkinson@unsw.edu.au",
  "erkelenz@mail.muni.cz",
  "emil.andersen@umu.se",
  "pia_marie.bradler@tu-dresden.de",
  "coraena.ls@gmail.com",
  "a.d.elsy@stir.ac.uk",
  "celeste.mare@bio.au.dk",
  "m.s.e.eshelman@sms.ed.ac.uk",
  "uckydic@gmail.com",
  "coskun.guclu@gmail.com",
  "eckbergj@umich.edu",
  "jonathan.vonoppen@unibas.ch",
  "josef.garen@botany.ubc.ca",
  "juliachacon@gmail.com",
  "marc.maciasfauria@ouce.ox.ac.uk",
  "sonya.geange@uib.no",
  "sehoya.cotner@uib.no",
  "paulefrens@gmail.com"
)

# Groups 
# 
# China: 
#   Fluxes: 4
# Traits: 1
# Peru 2018: 
#   Fluxes: 4
# Traits: 1
# Svalbard 2018: 
#   Fluxes: 4
# Traits: 3
# Peru 2020
# Fluxes: 4
# Traits: 1
# Norway 2022
# Fully unclear 
# South Africa 2023
# Fluxes: 3
# Traits: 1


dt_inc[dt_inc$ask == "yes" | dt_inc$ask == "maybe", ] %>% 
  select(Participant, Email) %>% 
  unique()

#old list 
dt_inc[dt_inc$ask == "yes" | dt_inc$ask == "maybe", ] %>% 
  select(Participant, Email) %>%
  mutate(Email = tolower(Email)) %>% 
  filter(!Email %in% sent_to) %>% 
  unique()

dt_inc[dt_inc$ask == "yes" | dt_inc$ask == "maybe", ] %>% 
  select(Participant, Email) %>%
  mutate(Email = tolower(Email)) %>% 
  filter(Email %in% sent_to) %>% 
  unique()


should_mails = dt_inc[dt_inc$ask == "yes" | dt_inc$ask == "maybe", ] %>% 
  select(Email) %>%
  unique() %>%
  mutate(Email = tolower(Email)) %>% pull()

#folks who got the mail but are not in the list... 
setdiff(sent_to, should_mails)

