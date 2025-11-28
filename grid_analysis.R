library(tidyverse)
library(vctrs)
library(s2)
#library(sf)
options(width = 140)

# grid building
make_grid <- function(level, ..., model = c("IVEA3H", "ISEA3H", "FULLER3H"), engine = c("DGGAL", "DGGRIDR")) {
  cache_uri <- fs::path("grid_cache", sprintf("%s_%s_%s.qs2", engine, model, level))

  if (fs::file_exists(cache_uri)) {
    grid <- qs2::qs_read(cache_uri)
  } else {
    cli::cli_progress_step("Running {engine}.{model}.{level}")
    fun <- get(paste0("make_grid.", engine), mode = "function")
    grid <- fun(level, ..., model = model)
    qs2::qs_save(grid, cache_uri)
  }

  grid
}


make_grid.DGGRIDR <- function(level, model = c("ISEA3H", "FULLER3H")) {
  model <- rlang::arg_match(model)
  proj <- switch(model,
    FULLER3H = "FULLER",
    ISEA3H = "ISEA"
  )
  wsg84 <- sf::st_crs(4326L)
  empty_line <- sf::st_linestring() |> sf::st_sfc(crs = wsg84)

  dggridR::dgconstruct(res = level, aperture = 3, projection = proj, topology = "HEXAGON") |>
    dggridR::dgearthgrid() |>
    tibble::as_tibble() |>
    dplyr::mutate(
      id = seq_len(dplyr::n()),
      shape = sapply(geometry, function(xy) {
        nrow(xy[[1]]) - 1L
      }),
      centroid = s2::s2_centroid(geometry),
      area = s2::s2_area(geometry)/1e6,
      # dggridR does not provide neighborhood information, so we calculate it manually
      # we take advantage of the fact that the cells are hexagons to find all centroids
      # within  hexagon circumradius + error
      neighbors = {
        # look for neighbours within 120% of the circumradius
        R <- sqrt(max(area)*2/3/sqrt(3))*1e3
        max_dist <- (R*2)*1.2

        neighbors <- parallel::mclapply(id, mc.cores = 8L, function(i) {
          # brute-force neighbor search (s2_prepared_dwithin indexes using the first argument)
          nn <- which(s2_prepared_dwithin(centroid, centroid[[i]], max_dist))
          nn <- nn[nn != i]

          if (length(nn) != shape[[i]]) return(integer())
          nn
        })

        if (any(vctrs::list_sizes(neighbors) == 0)) stop("neighbor detection failed")
        neighbors
      },
      # distances from the grid cell center to its corners
      radial_distances = {
        xy <- sf::st_coordinates(geometry) |> as_tibble() |> slice(-1, .by = L2)
        id <- xy$L2
        xy <- s2::s2_geog_point(xy$X, xy$Y)



        # centroid to corner distance
        d <- s2::s2_distance(centroid[id], xy)/1e3
        vec_chop(d, sizes = field(vec_group_rle(id), "length"))
      },
      # cell edge lengths
      edge_lengths = {
        xy <- sf::st_coordinates(geometry) |> as_tibble()
        id <- xy$L2
        xy <- s2::s2_geog_point(xy$X, xy$Y)

        purrr::map(vec_chop(xy, sizes = field(vec_group_rle(id), "length")), function(xy) {
          s2_distance(xy[-1], xy[-length(xy)])/1e3
        })
      }
    ) |>
    dplyr::select(
      id,
      shape,
      centroid,
      area,
      neighbors,
      radial_distances,
      edge_lengths,
      geometry
    )
}

make_grid.DGGAL <- function(level, ..., model = c("IVEA3H", "ISEA3H")) {
  model <- rlang::arg_match(model)

  # run the grid tool
  json_file <- withr::local_tempfile(fileext = ".jsonl")
  json <- system2("docker", c(
    "run",
    "--rm",
    # "--tty",
    "ghcr.io/human-ecology-grid/pipeline",
    "ivea-grid",
    paste("-t ", model),
    paste("--level ", level),
    "-s 0"
  ), stdout = json_file)

  yyjsonr::read_ndjson_file(json_file) |>
    as_tibble() |>
    mutate(
      centroid = {
        xy <- unlist(centroid)
        s2::s2_geog_point(xy[seq(1, length(xy), by = 2)], xy[seq(2, length(xy), by = 2)])
      },
      # we use base vertices to compute radial distances and edge lengths
      # these can be recovered from the neighbour edge data
      radial_distances = map2(centroid, neighbors, function(cc, X) {
        xy <- unlist(X$neighbor_edge)
        n <- length(xy)
        xy <- tibble(lon = xy[seq.int(1L, n, by = 2L)], lat = xy[seq.int(2L, n, by = 2L)]) |> distinct()
        xy <- s2::s2_geog_point(xy$lon, xy$lat)


        s2::s2_distance(cc, xy)/1e3
      }),
      edge_lengths = map(neighbors, function(X) {
        xy <- unlist(X$neighbor_edge)
        n <- length(xy)

        edges <- s2::s2_make_line(
          xy[seq.int(1L, n, by = 2L)],
          xy[seq.int(2L, n, by = 2L)],
          feature_id = rep(seq_len(nrow(X)), each = 2L)
        )

        s2::s2_length(edges)/1e3
      }),
      # only need neighbour ids
      neighbors = map(neighbors, function(X) {
        X$neighbor_id
      }),
      # validate and process cell geometry
      geometry = {
        sizes <- list_sizes(geometry)
        xy <- list_unchop(geometry)

        geometry <- s2::s2_make_polygon(xy$lon, xy$lat, feature_id = rep.int(seq_along(sizes), sizes), check = FALSE)

        # repair polygons
        invalid_at <- which(!s2::s2_is_valid(geometry))
        if (length(invalid_at) > 0) {
          geometry[invalid_at] <- s2::s2_rebuild(geometry[invalid_at], options = s2::s2_options(
            split_crossing_edges = TRUE,
            dimensions = "polygon"
          ))
          stopifnot(all(s2::s2_is_valid(geometry)))
        }

        geometry
      },
      area = s2::s2_area(geometry)/1e6
    ) |>
    select(
      id,
      shape,
      centroid,
      area,
      neighbors,
      radial_distances,
      edge_lengths,
      geometry
    )
}


# calculate grid statistics
calc_grid_stats <- function(grid) {
  # only care about hexagonal cells
  mutate(grid, is_hex = (shape == 6)) |>
  mutate(edge_variability = map_dbl(edge_lengths, function(x) max(x)/min(x) - 1)) |>
  summarize(
    n = nrow(grid),
    area_variability = max(area[is_hex])/min(area[is_hex]) - 1,
    # spacing
    {
      nn <- pick(id, neighbors) |> unnest(neighbors) |> filter(id < neighbors)
      nn <- filter(nn, is_hex[match(nn$id, id)] & is_hex[match(nn$neighbors, id)])
      spacing <- s2::s2_distance(centroid[match(nn$id, id)], centroid[match(nn$neighbors, id)])/1e3

      xx <- tibble(
        #avg_spacing = mean(spacing),
        # cv_spacing = sd(spacing)/mean(spacing),
        spacing_variability = max(spacing)/min(spacing) - 1
      )
    },
    avg_edge_variability = mean(edge_variability[is_hex]),
    max_edge_variability = max(edge_variability[is_hex]),
    #sd_edge_variability = sd(edge_variability[is_hex]),
  ) |>
  mutate(across(ends_with("_variability") , function(x)
    round(x*100, 1)
  ))
}

# run grid stats calculation
grid_parameters <- expand_grid(
  level = seq.int(5, 9),
  tribble(
    ~ model,    ~ engine,
    "FULLER3H",   "DGGRIDR",
    "ISEA3H",     "DGGRIDR",
    "ISEA3H",     "DGGAL",
    "IVEA3H",     "DGGAL"
  )
) |> arrange(engine, model)


source("plot_support.R")

plot_grid_edge_variability <- function(grid, file) {



  grid_sf <- grid |>
    mutate(edge_variability = map_dbl(edge_lengths, function(x) max(x)/min(x) - 1)) |>
    mutate(geometry = s2::s2_intersection(geometry, pacific_centered_clip_plate)) |>
    sf::st_as_sf(sf_column_name = "geometry")


  gradient_steps <- c(
    # zero (neutral)
    "#F5F0E1",
    # positives (blues):
    "#BDD7E7",
    "#6BAED6",
    "#3182BD",
    "#08519C"
  )

  save_plot(file, width = 1, plot = {
    ggplot(grid_sf) +
    geom_sf(data = landmass_sf, fill = "orange") +
    geom_sf(aes(fill = edge_variability), alpha = 0.75) +
    scale_fill_gradientn(
      colors = gradient_steps,
      limits = c(0, 0.25),
      space = "Lab",
      name = "",
      oob = scales::squish,
      labels = function(x) paste0(x*100, "%")
    ) +
    # colorspace::scale_fill_continuous_sequential(
    #   "TealGrn",
    #   name = "",
    #   limits = c(0, 0.3),
    #   na.value = "gray90"
    # ) +
    coord_sf(crs = 8859, expand = FALSE) +
    map_theme +
    theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 9),
      legend.key.size = ggplot2::unit(1, 'cm'),
      legend.key.height = ggplot2::unit(6, "pt"),
      legend.box.spacing = ggplot2::unit(-0.25, "cm"),
      legend.box.margin = margin(t = 8)
    )
  })
}

plot_grid_edge_variability(
  make_grid(level = 6, engine = "DGGAL", model = "IVEA3H"),
  "plots/ivea_edge_length_variability.png"
)

plot_grid_edge_variability(
  make_grid(level = 6, engine = "DGGRIDR", model = "ISEA3H"),
  "plots/isea_edge_length_variability.png"
)


stop()


stats <- pmap_dfr(grid_parameters, function(level, engine, model) {
  grid <- make_grid(level = level, engine = engine, model = model)
  stats <- calc_grid_stats(grid)
  tibble(
    model = model,
    engine = engine,
    level = level,
    calc_grid_stats(grid)
  )
})

write_csv(stats, "grid_stats.csv")

stats
