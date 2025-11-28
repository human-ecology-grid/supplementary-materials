# landcover map for pacific-centered plotting
pacific_centered_clip_plate <- local({
  m <- -30.0
  tol <- 0.01

  s2::s2_make_polygon(
    c(m - tol, -180.0, m + tol, m + tol, -180.0, m - tol),
    c(90.0, 90.0, 90.0, -90.0, -90.0, -90.0),
    oriented = TRUE
  )
})

landmass_sf <- local({
  rnaturalearth::countries110 |>
  tibble::as_tibble() |>
  dplyr::filter(CONTINENT != "Antarctica") |>
  dplyr::pull(geometry) |>
  s2::as_s2_geography(check = FALSE) |>
  s2::s2_union_agg(options = s2::s2_options(dimensions = "polygon")) |>
  s2::s2_intersection(pacific_centered_clip_plate) |>
  sf::st_as_sf()
})

page_width <- 14

save_plot <- function(file, width = 1.0, asp = 0.5, dpi = 300, plot) {
  ggsave(file, width = width*page_width, height = asp*width*page_width, dpi = dpi, units = "cm", plot = plot)
}

# common map theme
map_theme <- ggplot2::theme_minimal() + ggplot2::theme(
  legend.position = "bottom",
  legend.title = ggplot2::element_text(size = 11),
  legend.text = ggplot2::element_text(size = 7),
  legend.key.size = ggplot2::unit(1, 'cm'),
  legend.key.height = ggplot2::unit(6, "pt"),
  legend.box.spacing = ggplot2::unit(-0.25, "cm"),
  legend.box.margin = ggplot2::margin(0, 0, 0, -20),
  axis.title = ggplot2::element_blank(),
  axis.text = ggplot2::element_blank(),
  axis.ticks = ggplot2::element_blank(),
  axis.line = ggplot2::element_blank(),
  panel.border = ggplot2::element_blank(),
  plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm"),
  panel.grid = ggplot2::element_blank(),
  panel.grid.major  = ggplot2::element_line(colour = "grey80", linewidth   = 0.2, linetype = "dotted"),
  panel.background = ggplot2::element_blank()
)

# common plot theme
plot_theme <- ggplot2::theme_minimal() + ggplot2::theme(
  axis.title  = element_text(size = 9, colour = "grey20"),
  axis.text = element_text(size = 6, colour = "grey20"),
  axis.title.x = element_text(margin = margin(t = 10, b = 6)),
  axis.title.y = element_text(margin = margin(r = 10, l = 6)),
  legend.position = "bottom",
  legend.title = ggplot2::element_text(size = 11),
  legend.text = ggplot2::element_text(size = 7),
  legend.key.size = ggplot2::unit(1, 'cm'),
  legend.key.height = ggplot2::unit(6, "pt"),
  legend.box.spacing = ggplot2::unit(-0.25, "cm"),
  legend.box.margin = ggplot2::margin(0, 0, 0, -20),
  plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm"),
  panel.grid = ggplot2::element_blank(),
  panel.grid.major  = ggplot2::element_line(colour = "grey80", linewidth   = 0.2, linetype = "dotted"),
  panel.background = ggplot2::element_blank()
)
