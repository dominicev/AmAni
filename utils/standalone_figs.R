# utils/standalone_figs.R — Figures that can be run on their own
#
#   plot_spatial_extent_overlap()  animal vs plant spatial-extent comparison
#                                  (also called by summarize_all() /
#                                  summarize_organism())
#   compile_workflow_diagram()     annotated pipeline flowchart
#
# Sourced by AmAni.summarize.R, which loads the packages; uses
# organism_label() and ORGANISM_GROUPS from utils/func_core.R.

# =============================================================================
# --- 7. Spatial extent overlap comparison ------------------------------------
# =============================================================================
#
# plot_spatial_extent_overlap() loads per-organism spatial extent distributions
# (animal + plant) and produces a faceted comparison plot for the organisms
# with the smallest animal spatial extents.
#
# Requires plotdata_spatial_extent_distributions.rds (saved by the pipeline
# after the next re-run) or falls back to plotdata_plant_extent_model.rds
# (plant data only — available now for most organisms).

plot_spatial_extent_overlap <- function(include = NULL, group = NULL,
                                        max_cells = 10,
                                        output_dir = NULL, save_pdf = TRUE) {

  if (!exists("project_root")) source("AmAni.manager.R")
  if (is.null(output_dir)) output_dir <- file.path(project_root, "output")

  all_dirs <- list.dirs(output_dir, full.names = FALSE, recursive = FALSE)
  all_dirs <- all_dirs[!all_dirs %in% c("cross_organism_summary", "shared")]

  if (!is.null(group) && group %in% names(ORGANISM_GROUPS)) {
    include <- union(include, ORGANISM_GROUPS[[group]])
  }
  if (!is.null(include)) {
    all_dirs <- all_dirs[all_dirs %in% include]
  }

  # Load combined distributions (preferred) or plant-only fallback
  combined_data <- purrr::map_dfr(all_dirs, function(org) {
    pd_dir <- file.path(output_dir, org, "animal_model", "plot_data")

    # Preferred: saved combined distribution (animal + plant)
    combined_path <- file.path(pd_dir, "plotdata_spatial_extent_distributions.rds")
    if (file.exists(combined_path)) {
      d <- readRDS(combined_path)
      return(dplyr::mutate(d, organism = org))
    }

    # Fallback: plant-only from extent model
    ext_path <- file.path(pd_dir, "plotdata_plant_extent_model.rds")
    if (file.exists(ext_path)) {
      d <- readRDS(ext_path)
      return(dplyr::mutate(d$plant_df, organism = org))
    }
    NULL
  })

  if (nrow(combined_data) == 0) {
    message("No spatial extent data found.")
    return(invisible(NULL))
  }

  # Show all matched organisms, filtered to 1..max_cells on x-axis
  focus_orgs <- unique(combined_data$organism)

  plot_data <- combined_data |>
    dplyr::filter(num_grid_cells <= max_cells) |>
    dplyr::mutate(
      organism = organism_label(organism),
      num_grid_cells = factor(num_grid_cells, levels = seq_len(max_cells))
    )

  n_orgs <- dplyr::n_distinct(plot_data$organism)

  p <- ggplot(plot_data, aes(x = num_grid_cells, y = n, fill = group)) +
    geom_col(position = "dodge", color = "grey30", linewidth = 0.2) +
    scale_fill_manual(
      values = c("Plant species" = "#2d6a4f", "Animal species" = "#457b9d"),
      na.value = "grey60", name = NULL
    ) +
    facet_wrap(~organism, scales = "free_y",
               ncol = min(n_orgs, 4)) +
    labs(
      title = paste0("Spatial Extent Distribution (1\u2013", max_cells, " grid cells)"),
      subtitle = if (!is.null(group)) paste0("Group: ", group) else NULL,
      x = "Number of grid cells occupied",
      y = "Number of species"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 7),
      legend.position = "top"
    )

  if (save_pdf) {
    summary_dir <- file.path(output_dir, "cross_organism_summary")
    dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    out_file <- file.path(summary_dir,
                          paste0("spatial_extent_overlap_", ts, ".pdf"))
    pdf(out_file, width = 12, height = 8)
    print(p)
    dev.off()
    message("Wrote: ", out_file)
  }

  p
}

# =============================================================================
# --- 8. Workflow compilation diagram -----------------------------------------
# =============================================================================
#
# compile_workflow_diagram() produces a ggplot-based flowchart illustrating
# the major input/output stages of the pipeline, annotated with example
# values from a chosen organism.

compile_workflow_diagram <- function(organism = "butterflies",
                                     output_dir = NULL,
                                     save_pdf = TRUE) {

  if (!exists("project_root")) source("AmAni.manager.R")
  if (is.null(output_dir)) output_dir <- file.path(project_root, "output")

  # --- Helper: extract first number after a colon from text lines ---
  .ev <- function(lines, pattern) {
    l <- grep(pattern, lines, value = TRUE, ignore.case = TRUE)
    if (!length(l)) return("?")
    raw <- sub(".*:\\s*", "", l[1])
    m <- regmatches(raw, regexpr("[0-9][0-9,]*\\.?[0-9]*", raw))
    if (!length(m)) "?" else m
  }

  org_am  <- file.path(output_dir, organism, "animal_model")
  org_ext <- file.path(output_dir, organism, "extinction_model")

  # --- Load example values ---
  pre_cells <- pre_spp <- post_cells <- post_spp <- "?"
  ds_files <- list.files(org_am, "^dataset_summary.*\\.txt$", full.names = TRUE)
  if (length(ds_files)) {
    ds <- readLines(ds_files[which.max(file.mtime(ds_files))])
    pre_cells <- .ev(ds, "^\\s*Grid cells:")
    pre_spp   <- .ev(ds, "^\\s*Total species:")
    ai <- grep("After filtering", ds)
    if (length(ai)) {
      da <- ds[ai[1]:length(ds)]
      post_cells <- .ev(da, "^\\s*Grid cells:")
      post_spp   <- .ev(da, "^\\s*Total species:")
    }
  }

  gamma_best <- "?"
  gm_files <- list.files(org_am, "^gamma_richness.*\\.txt$", full.names = TRUE)
  if (length(gm_files)) {
    gm <- readLines(gm_files[which.max(file.mtime(gm_files))])
    gamma_best <- .ev(gm, "Best estimate:")
  }

  med_lost <- "?"
  ext_file <- file.path(org_ext, "summary_all_models.csv")
  if (file.exists(ext_file)) {
    ext <- read.csv(ext_file)
    r <- ext[ext$scenario == "FreqBoundary" & ext$model == "Clustered", ]
    if (nrow(r)) med_lost <- as.character(round(r$median_spp_lost[1]))
  }

  org_display <- organism_label(organism)

  # --- Box definitions: cx, cy, label, fill ---
  bw <- 2.0   # box half-width
  bh <- 0.65  # box half-height

  boxes <- dplyr::tribble(
    ~cx,   ~cy,  ~label,                                                       ~fill,
    # Plant lane (top)
    1.5,    3.0,  "GBIF Plants\n(Amazon, ~9k spp)",                            "#2d6a4f",
    5.0,    3.0,  "Stage 0\u20131\nPlant grid + endemism\n(22 km cells)",      "#52b788",

    # Animal lane (middle)
    1.5,    1.5,  paste0(org_display, "\n",
                         pre_spp, " spp / ", pre_cells, " cells"),             "#1d3557",
    5.0,    1.5,  paste0("Stage 2\n",
                         post_spp, " spp / ", post_cells, " cells\n",
                         "Extent model + scaling"),                            "#457b9d",
    8.5,    1.5,  paste0("Simulated Occ.\n\u2248",
                         gamma_best, " spp\n(gamma est.)"),                    "#264653",

    # Deforestation lane (bottom)
    1.5,    0.0,  "RAISG Raster\n(2001\u20132020)",                            "#6d4c41",
    5.0,    0.0,  "Stage 3\nPatch extraction\nLognormal + power-law",          "#a1887f",

    # Extinction + results (right)
    11.5,   0.75, paste0("Stage 4: Extinction\nNeutral / Clustered\n",
                         "/ Shifted (100 reps)"),                              "#c1121f",
    14.5,   0.75, paste0("Results\n\u2248", med_lost,
                         " spp lost\n(FreqBound., Clustered)"),                "#e76f51"
  )

  boxes <- boxes |>
    dplyr::mutate(
      xmin = cx - bw / 2, xmax = cx + bw / 2,
      ymin = cy - bh,      ymax = cy + bh
    )

  # --- Arrow definitions ---
  arrows <- dplyr::tribble(
    ~x1,    ~y1,   ~x2,    ~y2,
    # Plant GBIF → Stage 0-1
    2.5,    3.0,   3.0,    3.0,
    # Stage 0-1 → Stage 2 (plant endemism reference)
    5.0,    2.35,  5.0,    2.15,
    # Animal data → Stage 2
    2.5,    1.5,   3.0,    1.5,
    # Stage 2 → Simulated occ.
    6.0,    1.5,   6.5,    1.5,
    # Simulated occ → Stage 4
    9.5,    1.5,   10.5,   1.15,
    # Deforestation → Stage 3
    2.5,    0.0,   3.0,    0.0,
    # Stage 3 → Stage 4
    6.0,    0.0,   10.5,   0.4,
    # Stage 4 → Results
    12.5,   0.75,  13.5,   0.75
  )

  # --- Build plot ---
  p <- ggplot() +
    geom_rect(data = boxes,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = boxes$fill, color = "white", linewidth = 0.9,
              alpha = 0.92) +
    geom_text(data = boxes,
              aes(x = cx, y = cy, label = label),
              size = 2.5, color = "white", fontface = "bold",
              lineheight = 1.0) +
    geom_segment(data = arrows,
                 aes(x = x1, y = y1, xend = x2, yend = y2),
                 arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
                 color = "grey40", linewidth = 0.55) +
    # Lane labels
    annotate("text", x = 0.0, y = 3.0, label = "Plants",
             angle = 90, size = 3.2, color = "#2d6a4f", fontface = "bold") +
    annotate("text", x = 0.0, y = 1.5, label = "Animals",
             angle = 90, size = 3.2, color = "#1d3557", fontface = "bold") +
    annotate("text", x = 0.0, y = 0.0, label = "Deforestation",
             angle = 90, size = 3.2, color = "#6d4c41", fontface = "bold") +
    coord_cartesian(xlim = c(-0.5, 16.5), ylim = c(-1.0, 3.9),
                    expand = FALSE) +
    labs(title = paste0("AmAni Pipeline Overview: ", org_display),
         subtitle = "Major input/output stages (minor filtering steps omitted)") +
    theme_void() +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 10, 10)
    )

  if (save_pdf) {
    summary_dir <- file.path(output_dir, "cross_organism_summary")
    dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    out_file <- file.path(summary_dir,
                          paste0("workflow_diagram_", organism, "_", ts, ".pdf"))
    pdf(out_file, width = 18, height = 6)
    print(p)
    dev.off()
    message("Wrote: ", out_file)
  }

  p
}
