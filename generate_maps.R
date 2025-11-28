suppressPackageStartupMessages({
library(tidyverse)
library(rstan)
library(brms)
library(tidybayes)
})

source("plot_support.R")

# database preset
data_preset <- "500km-snapshot-latest@v1.0"
humanEcologyGrid::download_preset(data_preset)

# create outputs
fs::dir_create("plots")

grid <- humanEcologyGrid::get_grid_cells(data_preset, split_at_meridian = -30) |> filter(!cell_is_ocean_ice)
data <- humanEcologyGrid::get_data(data_preset) |>
  select(
    cell_id,
    railways_density,
    languages_diversity_index,
    population_urban_num,
    population_total_num,
    gdp_ppp_usd,
    land_roughness_sd,
    msa_mean
  ) |>
  mutate(
    population_urban_prop = population_urban_num/population_total_num
  ) |>
  collect() |>
  drop_na() |>
  left_join(as_tibble(grid) |> select(cell_id, geometry), join_by(cell_id))


  save_plot("plots/railways.png", width = 1, plot = {
    ggplot(sf::st_as_sf(data)) +
    geom_sf(data = landmass_sf, fill = "lavender") +
    geom_sf(aes(fill = railways_density), alpha = 0.9) +
    colorspace::scale_fill_continuous_sequential(
      "Heat 2",
      name = "",
      na.value = "gray90",
      breaks = scales::pretty_breaks(n = 5)
    ) +
    coord_sf(crs = 8859, expand = FALSE) +
    map_theme
  })


  save_plot("plots/gdp.png", width = 1, plot = {
    ggplot(sf::st_as_sf(data)) +
    geom_sf(data = landmass_sf, fill = "lavender") +
    geom_sf(aes(fill = gdp_ppp_usd), alpha = 0.9) +
    colorspace::scale_fill_continuous_sequential(
      "Heat 2",
      name = "USD",
      transform = "log10",
      breaks = scales::log_breaks(n = 6),
      labels = scales::label_math(10^.x, format = log10),
      na.value = "gray90",
    ) +
    coord_sf(crs = 8859, expand = FALSE) +
    map_theme
  })


  save_plot("plots/urban_population.png", width = 1, plot = {
    ggplot(sf::st_as_sf(data)) +
    geom_sf(data = landmass_sf, fill = "lavender") +
    geom_sf(aes(fill = population_urban_prop), alpha = 0.9) +
    colorspace::scale_fill_continuous_sequential(
      "Heat 2",
      name = "",
      na.value = "gray90",
      limits = c(0, 1),
      breaks = scales::pretty_breaks(n = 5)
    ) +
    coord_sf(crs = 8859, expand = FALSE) +
    map_theme
  })

  save_plot("plots/terrain_rri_sd.png", width = 1, plot = {
    ggplot(sf::st_as_sf(data)) +
    geom_sf(data = landmass_sf, fill = "lavender") +
    geom_sf(aes(fill = land_roughness_sd), alpha = 0.9) +
    colorspace::scale_fill_continuous_sequential(
      "Heat 2",
      name = "m",
      na.value = "gray90",
      breaks = scales::pretty_breaks(n = 5)
    ) +
    coord_sf(crs = 8859, expand = FALSE) +
    map_theme
  })

  save_plot("plots/msa.png", width = 1, plot = {
    ggplot(sf::st_as_sf(data)) +
    geom_sf(data = landmass_sf, fill = "lavender") +
    geom_sf(aes(fill = msa_mean), alpha = 0.9) +
    colorspace::scale_fill_continuous_sequential(
      "Heat 2",
      name = "",
      na.value = "gray90",
      limits = c(0, 1),
      breaks = scales::pretty_breaks(n = 5)
    ) +
    coord_sf(crs = 8859, expand = FALSE) +
    map_theme
  })

  save_plot("plots/languages_diversity.png", width = 1, plot = {
    ggplot(sf::st_as_sf(data)) +
    geom_sf(data = landmass_sf, fill = "lavender") +
    geom_sf(aes(fill = languages_diversity_index), alpha = 0.9) +
    colorspace::scale_fill_continuous_sequential(
      "Heat 2",
      name = "",
      na.value = "gray90",
      limits = c(0, 1),
      breaks = scales::pretty_breaks(n = 5)
    ) +
    coord_sf(crs = 8859, expand = FALSE) +
    map_theme
  })


  cli::cli_alert_success("Done")
