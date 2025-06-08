### calculate functional diversity

library(mFD)
library(tidyverse)
library(data.table)
library(traitstrap)

traits <- fread("data/processed_data/preliminary_data/prelim_traits.csv") %>%
  dplyr::select(plot_id, site, tier, trait_name, trait_value, species) %>%
  filter(!species == "") %>%
  filter(!(trait_name == "sla_cm2_g" & trait_value > 1000)) %>%
  filter(trait_name %in% c("sla_cm2_g", "ldmc", "leaf_area_cm2", "dry_mass_g", "plant_height_cm"))

community <- fread("data/processed_data/preliminary_data/prelim_cover.csv") %>%
  dplyr::select(plot_id, site, tier, cover, species) %>%
  filter(!species == "")


trait_filling <- trait_fill(
  # input data (mandatory)
  comm = community,
  traits = traits,

  # specifies columns in your data (mandatory)
  abundance_col = "cover",
  taxon_col = "species",
  trait_col = "trait_name",
  value_col = "trait_value",

  # specifies sampling hierarchy
  scale_hierarchy = c("tier", "site", "plot_id"),

  # min number of samples
  min_n_in_sample = 1
)

incomplete_trait_names <- trait_missing(
  filled_trait = trait_filling,
  comm = community
)

quantile(incomplete_trait_names[incomplete_trait_names$n_traits < 4, ]$max_abun, c(0, .05, .1, .2, .8, .9, .95, 1))



sp_tr_raw <- trait_filling %>%
  as.data.table() %>%
  filter(!plot_id == "") %>%
  filter(!grepl("NA", plot_id)) %>%
  pivot_wider(names_from = trait_name, values_from = trait_value, values_fn = mean) %>%
  group_by(species) %>%
  summarize(
    dry_mass_g = mean(dry_mass_g, na.rm = TRUE),
    leaf_area_cm2 = mean(leaf_area_cm2, na.rm = TRUE),
    sla_cm2_g = mean(sla_cm2_g, na.rm = TRUE),
    ldmc = mean(ldmc, na.rm = TRUE),
    plant_height_cm = mean(plant_height_cm, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(complete.cases(.))

sp_tr <- sp_tr_raw %>%
  remove_rownames() %>%
  column_to_rownames(var = "species")

# get cover fraction for each plot
dt_cover_frac <- data.table()

for (plot in unique(community$plot_id)) {
  sp_tr_raw_cov_i <- trait_filling %>%
    as.data.table() %>%
    filter(plot_id %in% c(plot))

  if (n_distinct(sp_tr_raw_cov_i$trait_name) < 4) {
    next
  }

  sp_tr_raw_cov <- sp_tr_raw_cov_i %>%
    pivot_wider(names_from = trait_name, values_from = trait_value, values_fn = mean) %>%
    group_by(plot_id, species) %>%
    summarize(
      dry_mass_g = mean(dry_mass_g, na.rm = TRUE),
      leaf_area_cm2 = mean(leaf_area_cm2, na.rm = TRUE),
      sla_cm2_g = mean(sla_cm2_g, na.rm = TRUE),
      ldmc = mean(ldmc, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    filter(complete.cases(.))

  unique_sp <- sp_tr_raw_cov %>%
    filter(plot_id == plot) %>%
    dplyr::select(species) %>%
    unique() %>%
    pull()

  cover_frac <- community %>%
    filter(plot_id == plot) %>%
    mutate(cover_sum_total = sum(cover, na.rm = T)) %>%
    filter(species %in% unique_sp) %>%
    mutate(
      cover_sum_traits = sum(cover, na.rm = T),
      cover_frac = cover_sum_traits / cover_sum_total
    ) %>%
    dplyr::select(plot_id, cover_frac) %>%
    unique()

  dt_cover_frac <- rbind(cover_frac, dt_cover_frac)
  quantile(dt_cover_frac$cover_frac)
}



tr_cat <- data.table(
  trait_name = c("dry_mass_g", "leaf_area_cm2", "sla_cm2_g", "ldmc", "plant_height_cm"),
  trait_type = c("Q", "Q", "Q", "Q", "Q"),
  trait_weight = c(1, 1, 1, 1, 1)
)



fspace <- mFD::tr.cont.fspace(
  sp_tr        = sp_tr,
  pca          = TRUE,
  nb_dim       = 4,
  scaling      = "scale_center",
  compute_corr = "pearson"
)

dist_mat <- as.matrix(fspace$sp_dist_init)
dist_mat[1:5, 1:5]

fdist <- mFD::funct.dist(
  sp_tr = sp_tr,
  tr_cat = tr_cat,
  metric = "euclidean",
  scale_euclid = "scale_center",
  weight_type = "equal"
)

fspaces_quality <- mFD::quality.fspaces(
  sp_dist = fdist,
)

round(fspaces_quality$"quality_fspaces", 3) %>% arrange(mad) # 6 dimensional space it is


p <- quality.fspaces.plot(
  fspaces_quality = fspaces_quality,
  quality_metric = "mad",
  fspaces_plot = c("pcoa_1d", "pcoa_2d", "pcoa_3d", "pcoa_4d", "pcoa_5d")
)
p


### build matrix for species weights (i.e., species as columns and rows contain their biomass/cover). Assemblages should be row names

weight_mat <- community %>%
  filter(!plot_id == "") %>%
  filter(!grepl("NA", plot_id)) %>%
  dplyr::select(species, cover, plot_id) %>%
  ## assign mean cover in case cover is missing
  dplyr::filter(species %in% unique(sp_tr_raw$species)) %>%
  pivot_wider(names_from = "species", values_from = "cover", values_fn = mean) %>%
  as.data.table() %>%
  mutate(
    across(where(is.list), ~ sapply(., toString)),
    across(where(is.character) & !all_of("plot_id"), ~ as.numeric(.)),
    across(everything(), ~ ifelse(is.na(.), 0, .))
  ) %>%
  remove_rownames() %>%
  column_to_rownames(var = "plot_id") %>%
  as.matrix()

### Get the occurrence dataframe:
asb_sp_summ <- mFD::asb.sp.summary(asb_sp_w = weight_mat)
asb_sp_occ <- asb_sp_summ$"asb_sp_occ"



#### get matrix of species coordinates

sp_faxes_coord <- fspaces_quality$"details_fspaces"$"sp_pc_coord"

#### distance based functional diversity metric (Hill numbers )
fd.hill.res <- alpha.fd.hill(asb_sp_w = weight_mat, sp_dist = fdist)

dt_dbfd <- fd.hill.res$asb_FD_Hill %>%
  as.data.frame() %>%
  rownames_to_column(var = "plot_id") %>%
  rename(
    functional_diversity_q1 = FD_q1,
    functional_diversity_q2 = FD_q2,
    functional_diversity_q0 = FD_q0,
  )


### compute distance based functional diversity metrics (needs at least 2 species per plot)
dt_sr <- community %>%
  group_by(plot_id) %>%
  mutate(species_richness = n_distinct(species)) %>%
  ungroup()

nono_plots <- unique(dt_sr[dt_sr$species_richness < 2, "plot_id"]) %>% pull()
weight_mat_multi <- weight_mat[!(rownames(weight_mat) %in% c(nono_plots)), ]

quantile(dt_sr$species_richness)

alpha_fd_indices <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord[, c("PC1", "PC2", "PC3", "PC4", "PC5")],
  asb_sp_w         = weight_mat_multi,
  ind_vect         = c("fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE
)

dt_afdi <- alpha_fd_indices$functional_diversity_indices %>%
  rownames_to_column(var = "plot_id") %>%
  as.data.table() %>%
  dplyr::select("plot_id", "fdis") %>%
  rename(
    functional_dispersion = fdis,
  )


alpha_fd_indices$functional_diversity_indices
hist(log(alpha_fd_indices$functional_diversity_indices$fdis))

hist(log(dt_dbfd$functional_diversity_q1))

dt_fun_div <- dt_dbfd %>%
  left_join(dt_afdi) %>%
  left_join(dt_cover_frac)

fwrite(dt_fun_div, "data/processed_data/preliminary_data/plant_functional_diversity.csv")
