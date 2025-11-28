suppressPackageStartupMessages({
library(tidyverse)
library(rstan)
library(brms)
library(tidybayes)
library(Matrix)
})

source("plot_support.R")

# database preset
data_preset <- "500km-snapshot-latest@v1.0"
humanEcologyGrid::download_preset(data_preset)

# create outputs
fs::dir_create("plots")
fs::dir_create("models")

cli::cli_h1("Case study 2: Marginalization and Language Endangerment")

# load the data
grid <- humanEcologyGrid::get_grid_cells(data_preset, split_at_meridian = -30) |> filter(!cell_is_ocean_ice)
data_case_study2_full <- humanEcologyGrid::get_data(data_preset) |>
  select(cell_id, epr_groups_num, epr_groups_rel_unpriviledged_num, languages_num, languages_endangered_num) |>
  mutate(p_endangered = languages_endangered_num/languages_num) |>
  mutate(p_marginalized = epr_groups_rel_unpriviledged_num/epr_groups_num) |>
  collect() |>
  drop_na() |>
  left_join(as_tibble(grid) |> select(cell_id, geometry), join_by(cell_id))

# for the model, only take the data with sufficient coverage
data_case_study2 <- data_case_study2_full |>
  filter(epr_groups_num >= 5) |>
  filter(languages_num >= 5)


# discard disconnected cells
W <- humanEcologyGrid::get_grid_adjacency(
  data_preset,
  output = "matrix",
  cell_ids = data_case_study2$cell_id,
  isolated_cell_policy = "allow"
)
while (any(rowSums(W) == 1)) {
  masked_cell_ids <- rownames(W)[rowSums(W) > 1]
  W <- W[masked_cell_ids, masked_cell_ids]
}
data_case_study2 <- filter(data_case_study2, cell_id %in% rownames(W))

data_case_study2_full <- mutate(
  data_case_study2_full,
  included = cell_id %in% data_case_study2$cell_id,
  p_endangered = if_else(included, p_endangered, NA),
  p_marginalized = if_else(included, p_marginalized, NA),
  alpha = if_else(included, 0.9, 0.5)
)


# run the model
case_study2_model <- local({
  # attempt to load the cached model
  # model_uri <- "models/case_study2_model_beta_binomial_dispersion.rds"
  model_uri <- "models/case_study2_model.rds"
  tryCatch(return(suppressWarnings(readRDS(model_uri))), error = force)

  # adjacency matrix
  W <- humanEcologyGrid::get_grid_adjacency(data_preset, output = "matrix", cell_ids = data_case_study2$cell_id)

  cli::cli_progress_step("[Case Study 2] Running the model")

  # Binomial model with exact sparse CAR spatial autoregressive component
  model <- brm(
    languages_endangered_num | trials(languages_num) ~ 1 +  log1p(epr_groups_rel_unpriviledged_num) + car(W, gr = cell_id, type = "escar"),
    data    = data_case_study2,
    data2   = list(W = W),
    family  = beta_binomial(),
    backend = "cmdstanr",
    chains = 4L,
    cores = 16L,
    threads = threading(4L),
    #iter   = 20000L,
    #warmup   = 15000L,
    iter = 20000L,
    #warmup  = 15000L,
    # control=list(adapt_delta=0.99, max_treedepth = 15),
    refresh = 100L
  )
  cli::cli_progress_done()

  saveRDS(model, model_uri)

  model
})

print(case_study2_model)

# calculate the posterior bands
case_study2_posterior_bands <- local({
  draws <- spread_draws(
    case_study2_model,
    b_Intercept,
    b_log1pepr_groups_rel_unpriviledged_num
  )

  # predictor range
  pred_range <- seq.int(0, max(data_case_study2$epr_groups_rel_unpriviledged_num), by = 0.01)

  # calculate the posterior bands
  draws |>
    crossing(epr_groups_rel_unpriviledged_num = pred_range) |>
    mutate(p_endangered = plogis(b_Intercept + b_log1pepr_groups_rel_unpriviledged_num * log1p(epr_groups_rel_unpriviledged_num))) |>
    group_by(epr_groups_rel_unpriviledged_num) |>
    median_hdi(p_endangered, width = 0.99)
})

print(case_study2_posterior_bands)


# plot case study 2 maps and data
save_plot("plots/map_marginalization.png", width = 1, plot = {
  ggplot(sf::st_as_sf(data_case_study2_full)) +
  geom_sf(data = landmass_sf, fill = "lavender") +
  geom_sf(aes(fill = p_marginalized, alpha = alpha)) +
  colorspace::scale_fill_continuous_sequential(
    "Heat 2",
    name = "",
    na.value = "gray90"
  ) +
  scale_alpha_identity() +
  coord_sf(crs = 8859, expand = FALSE) +
  map_theme
})

save_plot("plots/map_endangerment.png", width = 1, plot = {
  ggplot(sf::st_as_sf(data_case_study2_full)) +
  geom_sf(data = landmass_sf, fill = "lavender", alpha = 0.5) +
  geom_sf(aes(fill = p_endangered, alpha = alpha)) +
  colorspace::scale_fill_continuous_sequential(
    "Heat 2",
    name = "",
    na.value = "gray90"
  ) +
  scale_alpha_identity() +
  coord_sf(crs = 8859, expand = FALSE) +
  map_theme
})


save_plot("plots/marginalization_endangerment_model_plot.png", width = 1,  plot = {
  ggplot(case_study2_posterior_bands, aes(x = epr_groups_rel_unpriviledged_num, y = p_endangered)) +
    # plot the data
    geom_point(data = data_case_study2, alpha = 0.25, size = 1.5, shape = 16, stroke = 1, color = "#56B4E9") +
    # plot the 95 % HDI
    geom_ribbon(aes(ymin = p_endangered.lower, ymax = p_endangered.upper), fill = "#E69F00", color = NA, alpha = 0.2) +
    geom_line(linewidth = 1.2, alpha = 0.9, color = "#E69F00") +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = scales::percent_format(accuracy = 1)
    ) +
    annotate(
      "text",
      x = 0.5, y = 0.2,
      hjust = 0, vjust = 0,
      size = 2,
      angle = 7.5,
      label = "Expectation of the posterior predictive (99% HDI)\nspatial effect marginalized"
    ) +
    # axis labels
    labs(
      x = "Number of Politically Marginalized Groups in Grid Cell",
      y = "Proportion of Endangered Languages"
    ) +
    plot_theme
})


cat("Done\n")
