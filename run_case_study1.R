suppressPackageStartupMessages({
library(tidyverse)
library(rstan)
library(brms)
library(tidybayes)
library(Matrix)
})

source("plot_support.R")

# database preset
data_preset <- "500km-yearly@v1.0"
humanEcologyGrid::download_preset(data_preset)

# create outputs
fs::dir_create("plots")
fs::dir_create("models")

cli::cli_h1("Case study 1: NTL and GDP")

# Pacific-centered grid
grid <- humanEcologyGrid::get_grid_cells(data_preset, split_at_meridian = -30) |> filter(!cell_is_ocean_ice)

# data for the year 2005
data_case_study1 <- humanEcologyGrid::get_data(data_preset) |>
  select(cell_id, year, gdp_ppp_usd, ntl_radiance_mean = ntl_radiance_p99capped_mean) |>
  filter(year == 2005, ntl_radiance_mean > 0, gdp_ppp_usd > 0) |>
  collect() |>
  drop_na() |>
  left_join(as_tibble(grid) |> select(cell_id, geometry), join_by(cell_id))


# run the model
case_study1_model <- local({
  # attempt to load the cached model
  model_uri <- "models/case_study1_model.rds"
  tryCatch(return(suppressWarnings(readRDS(model_uri))), error = force)

  # remove disconnected cells
  W <- humanEcologyGrid::get_grid_adjacency(
    data_preset,
    output = "matrix",
    cell_ids = data_case_study1$cell_id,
    isolated_cell_policy = "allow"
  )
  while (any(rowSums(W) == 1)) {
    masked_cell_ids <- rownames(W)[rowSums(W) > 1]
    W <- W[masked_cell_ids, masked_cell_ids]
  }
  data_case_study1 <- filter(data_case_study1, cell_id %in% rownames(W))
  W <- humanEcologyGrid::get_grid_adjacency(data_preset, output = "matrix", cell_ids = data_case_study1$cell_id)

  cli::cli_progress_step("[Case Study 1] Running the model")
  # log-log Gausssian model with an exact sparse CAR spatial autoregressive component and explicitly modeled sd
  model <- brm(
    bf(
      log(gdp_ppp_usd) ~ 1 + log(ntl_radiance_mean) + car(W, gr = cell_id, type = "escar"),
      sigma ~ log(ntl_radiance_mean)
    ),
    data  = data_case_study1,
    data2   = list(W = W),
    family  = gaussian(),
    backend = "cmdstanr",
    chains = 4L,
    cores = 12L,
    threads = threading(6L),
    iter   = 26000L,
    warmup = 20000L,
    control=list(adapt_delta=0.99),
    refresh = 50L
  )
  cli::cli_progress_done()

  saveRDS(model, model_uri)

  model
})

print(case_study1_model)

# calculate the posterior bands
case_study1_posterior_bands <- local({
  draws <- spread_draws(
    case_study1_model,
    b_Intercept,
    b_logntl_radiance_mean,
    b_sigma_Intercept,
    b_sigma_logntl_radiance_mean
  ) |> ungroup()

  # predictor range (equal spacing on log scale)
  radiance_range <- range(data_case_study1$ntl_radiance_mean)
  radiance_range <- exp(seq(log(radiance_range[[1]]), log(radiance_range[[2]]), length.out = 50L))

  # calculate the posterior bands
  draws |>
    crossing(ntl_radiance_mean = radiance_range) |>
    mutate(sigma = b_sigma_Intercept + b_sigma_logntl_radiance_mean*log(ntl_radiance_mean)) |>
    mutate(log_gdp_ppp_usd = b_Intercept + b_logntl_radiance_mean * log(ntl_radiance_mean)) |>
    mutate(gdp_ppp_usd = exp(log_gdp_ppp_usd + 0.5 * sigma^2))|>
    group_by(ntl_radiance_mean) |>
    point_interval(gdp_ppp_usd, .width = 0.99, .point = median, .interval = function(x, ...) hdi(x, ..., n = 128))
})

print(case_study1_posterior_bands)


save_plot("plots/radiance_gdp_2005_model_plot.png", width = 1, plot = {
  ggplot(case_study1_posterior_bands, aes(x = ntl_radiance_mean, y = gdp_ppp_usd)) +
    # plot the data
    geom_point(data = data_case_study1, alpha = 0.25, size = 1.5, shape = 16, stroke = 1, color = "#56B4E9") +
    # plot the 95 % HDI
    #geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = "#E69F00", color = NA, alpha = 0.2) +
    #geom_line(linewidth = 1.2, alpha = 0.9, color = "#E69F00") +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = "#E69F00", color = NA, alpha = 0.2) +
    # geom_ribbon(aes(ymin = .lower_50, ymax = .upper_50), fill = "#E69F00", color = NA, alpha = 0.25) +
    geom_line(linewidth = 1.2, alpha = 0.9, color = "#E69F00") +
    annotate(
      "text",
      x = 1e-12, y = 0.8e3,
      hjust = 0, vjust = 0,
      angle = 10,
      size = 2,
      label = "Expectation of the posterior predictive (99% HDI)\nspatial effect marginalized"
    ) +
    # log scale
    scale_y_log10(breaks = scales::log_breaks(n = 6), labels = scales::label_math(10^.x, format = log10)) +
    scale_x_log10(breaks = scales::log_breaks(n = 7), labels = scales::label_math(10^.x, format = log10)) +
    # axis labels
    labs(
      x = "Mean NTL Radiance",
      y = "GDP (USD 2005)"
    ) +
    plot_theme
})


# plot case study 1 maps and data
save_plot("plots/map_ntl_radiance_2005.png", width = 1, plot = {
  ggplot(sf::st_as_sf(data_case_study1)) +
  geom_sf(data = landmass_sf, fill = "lavender") +
  geom_sf(aes(fill = ntl_radiance_mean), alpha = 0.9) +
  colorspace::scale_fill_continuous_sequential(
    "Heat 2",
    name = latex2exp::TeX("$W sr^{-1} m^{-2}$"),
    transform = "log10",
    breaks = scales::log_breaks(n = 7),
    labels = scales::label_math(10^.x, format = log10),
    na.value = "gray90"
  ) +
  coord_sf(crs = 8859, expand = FALSE) +
  map_theme
})

save_plot("plots/map_gdp_2005.png", width = 1, plot = {
  ggplot(sf::st_as_sf(data_case_study1)) +
  geom_sf(data = landmass_sf, fill = "lavender") +
  geom_sf(aes(fill = gdp_ppp_usd), alpha = 0.9) +
  colorspace::scale_fill_continuous_sequential(
    "Heat 2",
    name = "USD (2005)",
    transform = "log10",
    breaks = scales::log_breaks(n = 6),
    labels = scales::label_math(10^.x, format = log10),
    na.value = "gray90"
  ) +
  coord_sf(crs = 8859, expand = FALSE) +
  map_theme
})

cat("Done\n")
