# utils/build_completeness_correlation_plot.R
# Extracted from AmAni.summarize.R. Sourced by AmAni.summarize.R — relies on
# the packages (tidyverse, sf, patchwork, ggrepel) loaded there.

build_completeness_correlation_plot <- function(cfr) {
  d <- sf::st_drop_geometry(cfr$delta) |>
    dplyr::filter(!is.na(plant_completeness), !is.na(animal_completeness))
  if (nrow(d) < 3) return(NULL)

  pearson  <- suppressWarnings(cor(d$plant_completeness, d$animal_completeness))
  spearman <- suppressWarnings(cor(d$plant_completeness, d$animal_completeness,
                                   method = "spearman"))

  ggplot(d, aes(x = plant_completeness, y = animal_completeness)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "black", fill = "grey80", linewidth = 0.6) +
    geom_point(aes(colour = retained), size = 2.5, alpha = 0.8) +
    scale_colour_manual(values = c(`TRUE` = "#2c7fb8", `FALSE` = "grey60"),
                        labels = c(`TRUE` = "Retained", `FALSE` = "Filtered out"),
                        name = "Animal cell") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4, size = 3.5,
             label = sprintf("Pearson r = %.2f\nSpearman rho = %.2f\nn = %d cells",
                             pearson, spearman, nrow(d))) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = "Plant vs animal completeness in overlapping grid cells",
         subtitle = "Chao1 completeness per cell; dashed line = 1:1",
         x = "Plant completeness", y = "Animal completeness") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
}
