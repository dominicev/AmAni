
# computeCellStats.R
# Helper: compute per-site completeness statistics from occurrence data.
#
# Source this file wherever you need site-level completeness diagnostics —
# it has no side effects and does not modify the global environment.

# ─────────────────────────────────────────────────────────────────────────────
#' Compute per-site completeness statistics
#'
#' Given a table of individual occurrence records (one row per record), this
#' function aggregates abundance by site × species, then computes species
#' richness estimates and sampling completeness for each site using one of
#' several methods.
#'
#' @section Input format:
#' \code{occ_data} must be a \code{data.frame}, \code{data.table}, or \code{sf}
#' object where:
#'
#' \itemize{
#'   \item Each \strong{row} is one occurrence record (a single individual or
#'         detection event). Multiple rows for the same species in the same site
#'         are treated as repeated detections and count toward abundance.
#'   \item \strong{site_col} — name of the column holding the site / grid-cell
#'         identifier (e.g. \code{"grid_id"}). Values may be integer, character,
#'         or factor; they are used as grouping keys and are not modified.
#'   \item \strong{species_col} — name of the column holding the species /
#'         morphotype identifier (e.g. \code{"Morphotype"}). Each unique value
#'         is treated as a distinct taxon within a site.
#'   \item Geometry columns (for \code{sf} inputs) are dropped automatically
#'         before computation — you do not need to strip them beforehand.
#' }
#'
#' There is no requirement on row ordering; all other columns are ignored.
#'
#' @section Completeness methods:
#' All methods return a value in [0, 1] quantifying sampling completeness.
#'
#' \strong{Chao1} (default): ratio-based, \eqn{S_{obs} / \hat{S}_{Chao1}}.
#' Uses singletons (\eqn{f_1}) and doubletons (\eqn{f_2}).
#'
#' \strong{ACE}: Abundance-based Coverage Estimator.
#' \eqn{S_{obs} / \hat{S}_{ACE}}. Partitions species into rare (abundance
#' \eqn{\le 10}) and abundant (\eqn{> 10}); uses the full frequency
#' distribution of rare species for a more stable estimate when
#' \eqn{f_2} is small.
#'
#' \strong{Coverage}: bias-corrected Good-Turing sample coverage
#' (Chao & Jost 2012). Estimates the probability that a new individual
#' belongs to an already-observed species. Not a richness ratio but a
#' widely-used completeness proxy that accounts for abundance structure.
#'
#' @param occ_data   A \code{data.frame}, \code{data.table}, or \code{sf}
#'                   object of individual occurrence records. See Input format.
#' @param site_col   Character. Name of the site / grid-cell ID column.
#'                   Default \code{"grid_id"}.
#' @param species_col Character. Name of the species / morphotype column.
#'                   Default \code{"Morphotype"}.
#' @param method     Character. Completeness method: \code{"Chao1"} (default),
#'                   \code{"ACE"}, or \code{"Coverage"}.
#' @param export     Logical. If \code{TRUE}, save the result list to an RDS
#'                   file. Default \code{FALSE}.
#' @param export_path Character. File path for the RDS output (including
#'                   filename and \code{.rds} extension). Required when
#'                   \code{export = TRUE}. If the directory does not exist it
#'                   will be created.
#' @param verbose    Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{\code{spp_abund}}{Long-format \code{data.frame}: one row per
#'     site × species combination, with columns \code{<site_col>},
#'     \code{<species_col>}, and \code{abundance} (detection count).}
#'   \item{\code{cell_stats}}{Wide-format \code{data.frame}: one row per site,
#'     with columns \code{<site_col>}, \code{n_occurrences} (total records),
#'     \code{n_species} (observed richness), \code{f1} (singletons),
#'     \code{f2} (doubletons), \code{richness_est} (estimated richness;
#'     NA for Coverage method), and
#'     \code{completeness} (sampling completeness, 0–1).}
#' }
#'
#' @examples
#' \dontrun{
#' # From an sf object (geometry stripped automatically):
#' out <- compute_cell_stats(animal_grid)
#' out$cell_stats
#'
#' # Custom column names + export:
#' out <- compute_cell_stats(
#'   occ_data     = my_df,
#'   site_col     = "plot_id",
#'   species_col  = "species",
#'   export       = TRUE,
#'   export_path  = "outputs/cell_stats.rds"
#' )
#' }
# ─────────────────────────────────────────────────────────────────────────────
compute_cell_stats <- function(
    occ_data,
    site_col    = "grid_id",
    species_col = "Morphotype",
    method      = "Chao1",
    export      = FALSE,
    export_path = NULL,
    verbose     = TRUE
) {

  method <- match.arg(method, c("Chao1", "ACE", "Coverage"))

  # ── Input validation ────────────────────────────────────────────────────────
  if (!inherits(occ_data, c("data.frame", "data.table", "sf"))) {
    stop("occ_data must be a data.frame, data.table, or sf object.")
  }
  if (!site_col    %in% names(occ_data)) stop("site_col '",    site_col,    "' not found in occ_data.")
  if (!species_col %in% names(occ_data)) stop("species_col '", species_col, "' not found in occ_data.")

  if (export && (is.null(export_path) || !nzchar(export_path))) {
    stop("export_path must be provided when export = TRUE.")
  }

  # ── Drop sf geometry if present ─────────────────────────────────────────────
  if (inherits(occ_data, "sf")) {
    occ_data <- sf::st_drop_geometry(occ_data)
  }

  if (verbose) message("[compute_cell_stats] Computing per-site abundance...")

  # ── Step 1: site × species abundance counts ──────────────────────────────
  # Count how many times each species was recorded in each site.
  # Result: long-format table with one row per site-species pair.
  spp_abund <- occ_data[, c(site_col, species_col)]
  spp_abund <- aggregate(
    list(abundance = rep(1L, nrow(spp_abund))),
    by  = spp_abund[, c(site_col, species_col), drop = FALSE],
    FUN = sum
  )

  if (verbose) message("[compute_cell_stats] Aggregating to site-level statistics (method: ", method, ")...")

  # ── Step 2: per-site summary statistics ──────────────────────────────────
  # For each site: total records, observed richness, singleton count (f1),
  # doubleton count (f2). These feed the Chao1 estimator.
  cell_stats <- aggregate(
    cbind(
      n_occurrences = spp_abund$abundance,
      n_species     = 1L,             # will sum to observed richness
      f1            = as.integer(spp_abund$abundance == 1),
      f2            = as.integer(spp_abund$abundance == 2)
    ),
    by  = spp_abund[, site_col, drop = FALSE],
    FUN = sum
  )

  # ── Step 3: Richness estimate and completeness ────────────────────────────
  if (method == "Chao1") {
    # Chao1 bias-corrected: f1^2 / (2*f2) when f2 > 0, else f1*(f1-1)/2
    cell_stats$richness_est <- cell_stats$n_species + ifelse(
      cell_stats$f2 > 0,
      cell_stats$f1^2 / (2 * cell_stats$f2),
      cell_stats$f1 * (cell_stats$f1 - 1) / 2
    )
    cell_stats$completeness <- cell_stats$n_species / cell_stats$richness_est

  } else if (method == "ACE") {
    # Abundance-based Coverage Estimator.
    # Partition species into rare (abundance <= 10) and abundant (> 10).
    # Uses the full frequency distribution of rare species.
    ace_by_site <- tapply(
      spp_abund$abundance,
      spp_abund[[site_col]],
      function(abunds) {
        s_obs   <- length(abunds)
        rare    <- abunds[abunds <= 10]
        s_rare  <- length(rare)
        s_abund <- s_obs - s_rare
        if (s_rare == 0) return(c(richness_est = s_obs, completeness = 1))
        n_rare  <- sum(rare)
        f1      <- sum(rare == 1)
        if (n_rare == 0) return(c(richness_est = s_obs, completeness = 1))
        c_rare  <- 1 - f1 / n_rare
        if (c_rare == 0) return(c(richness_est = NA_real_, completeness = NA_real_))
        # Coefficient of variation (gamma^2) of rare species
        freq_counts <- tabulate(rare, nbins = 10)
        k_vec       <- seq_len(10)
        num         <- sum(k_vec * (k_vec - 1) * freq_counts)
        gamma_sq    <- max(s_rare / c_rare * num / (n_rare * (n_rare - 1)) - 1, 0)
        s_ace       <- s_abund + s_rare / c_rare + f1 * gamma_sq / c_rare
        c(richness_est = s_ace, completeness = s_obs / s_ace)
      }
    )
    ace_df <- do.call(rbind, ace_by_site)
    # Match to cell_stats row order
    idx <- match(cell_stats[[site_col]], rownames(ace_df))
    cell_stats$richness_est <- ace_df[idx, "richness_est"]
    cell_stats$completeness <- ace_df[idx, "completeness"]

  } else if (method == "Coverage") {
    # Bias-corrected Good-Turing sample coverage (Chao & Jost 2012).
    # C_hat = 1 - (f1/n) * ((n-1)*f1) / ((n-1)*f1 + 2*f2)
    n  <- cell_stats$n_occurrences
    f1 <- cell_stats$f1
    f2 <- cell_stats$f2
    cell_stats$richness_est <- NA_real_
    cell_stats$completeness <- ifelse(
      n == 0, NA_real_,
      ifelse(f1 == 0, 1,
             1 - (f1 / n) * ((n - 1) * f1) / ((n - 1) * f1 + 2 * f2))
    )
  }

  # Clamp to [0, 1] for safety
  cell_stats$completeness <- pmin(pmax(cell_stats$completeness, 0), 1)

  if (verbose) {
    message("[compute_cell_stats] Done (", method, "). Sites: ", nrow(cell_stats),
            " | Median completeness: ", round(median(cell_stats$completeness, na.rm = TRUE), 3))
  }

  # ── Assemble return object ───────────────────────────────────────────────
  out <- list(
    spp_abund  = spp_abund,
    cell_stats = cell_stats
  )

  # ── Optional export ──────────────────────────────────────────────────────
  if (export) {
    export_dir <- dirname(export_path)
    if (!dir.exists(export_dir)) {
      dir.create(export_dir, recursive = TRUE)
      if (verbose) message("[compute_cell_stats] Created directory: ", export_dir)
    }
    saveRDS(out, file = export_path)
    if (verbose) message("[compute_cell_stats] Results saved to: ", export_path)
  }

  return(out)
}

# ─────────────────────────────────────────────────────────────────────────────
#' Compute Chao2 (incidence-based) sampling completeness across sites
#'
#' Chao1/ACE/Coverage (see \code{compute_cell_stats()}) estimate completeness
#' \emph{within a single site} from that site's abundance distribution. Chao2
#' instead estimates a \strong{single, pooled (regional/gamma-scale)}
#' completeness value from incidence data \emph{across} a set of sites —
#' here, grid cells. Each grid cell is treated as one presence/absence
#' sampling unit; species detected in only 1 or 2 grid cells (\eqn{Q_1},
#' \eqn{Q_2} — the incidence analogues of singletons and doubletons) drive
#' the estimate of how many species remain undetected region-wide.
#'
#' Use this for an overall "proportion of species sampled" figure for the
#' whole dataset (e.g. \code{prop_species} in SimulateScaledRichness.R).
#' It does \strong{not} produce a per-cell value and should not be used to
#' filter individual grid cells — see \code{compute_cell_stats()} for that.
#'
#' @section Chao2 estimator (bias-corrected):
#' \deqn{
#'   \hat{S}_{Chao2} = S_{obs} + \frac{m-1}{m} \times
#'     \begin{cases}
#'       Q_1^2 / (2 Q_2) & \text{if } Q_2 > 0 \\
#'       Q_1 (Q_1 - 1) / 2 & \text{if } Q_2 = 0
#'     \end{cases}
#' }
#' where \eqn{m} is the number of grid cells (sampling units). Completeness
#' is \eqn{S_{obs} / \hat{S}_{Chao2}}, clamped to [0, 1].
#'
#' @param spp_abund   Long-format data.frame as returned by
#'                     \code{compute_cell_stats()$spp_abund} — one row per
#'                     site × species combination (abundance column is not
#'                     used; only presence per site matters here).
#' @param site_col    Character. Name of the site / grid-cell column.
#'                     Default \code{"grid_id"}.
#' @param species_col Character. Name of the species / morphotype column.
#'                     Default \code{"Morphotype"}.
#'
#' @return A named list: \code{S_obs}, \code{Q1}, \code{Q2}, \code{m}
#'   (number of grid cells), \code{chao2} (estimated total richness), and
#'   \code{completeness} (\eqn{S_{obs} / \hat{S}_{Chao2}}, in [0, 1]).
# ─────────────────────────────────────────────────────────────────────────────
compute_chao2_completeness <- function(
    spp_abund,
    site_col    = "grid_id",
    species_col = "Morphotype"
) {
  if (!site_col    %in% names(spp_abund)) stop("site_col '",    site_col,    "' not found in spp_abund.")
  if (!species_col %in% names(spp_abund)) stop("species_col '", species_col, "' not found in spp_abund.")

  # Number of distinct grid cells each species occurs in (incidence).
  incidence <- table(spp_abund[[species_col]])

  S_obs <- length(incidence)
  Q1    <- sum(incidence == 1)
  Q2    <- sum(incidence == 2)
  m     <- length(unique(spp_abund[[site_col]]))

  rare_term <- if (Q2 > 0) Q1^2 / (2 * Q2) else Q1 * (Q1 - 1) / 2
  chao2 <- S_obs + ((m - 1) / m) * rare_term
  completeness <- min(max(S_obs / chao2, 0), 1)

  list(S_obs = S_obs, Q1 = Q1, Q2 = Q2, m = m, chao2 = chao2, completeness = completeness)
}
