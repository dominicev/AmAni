# classify_iucn_organism.R — per-organism IUCN threat classification
#
# Wraps compute_iucn_category_counts() + plot_iucn_threat() for a single
# organism and writes, to its extinction_model/ directory:
#   <organism>_iucn_threat_classification.pdf
#   <organism>_iucn_threat_classification_data.rds
#
# The *_data.rds file is what summarize_organism(), summarize_all() and the
# publication figures read, so this must run after every extinction model run.
# Called at the end of run_extinction_model() and at the start of
# summarize_organism(). Requires project_root (set by AmAni.manager.R).

source(file.path(project_root, "4-ExtModeller", "scripts", "utils",
                 "iucn_threat_classification.R"))

# Args:
#   organism:        organism name (directory under output/)
#   extinction_dir:  the organism's extinction_model/ directory; defaults to
#                    <project_root>/output/<organism>/extinction_model
#   organism_label:  plot title (defaults to the organism name)
# Returns (invisibly) list(classification, plot), or NULL when the extinction
# model summaries are missing.
classify_iucn_organism <- function(organism,
                                   extinction_dir = NULL,
                                   organism_label = organism) {

  if (is.null(extinction_dir)) {
    extinction_dir <- file.path(project_root, "output", organism, "extinction_model")
  }

  thresholds_path <- file.path(extinction_dir, "summary_range_reduction_thresholds.csv")
  if (!file.exists(thresholds_path)) {
    message("Skipping IUCN classification for '", organism,
            "' — file not found: ", thresholds_path)
    return(invisible(NULL))
  }

  message("Classifying IUCN threat categories for: ", organism)
  classification <- compute_iucn_category_counts(extinction_dir)
  p <- plot_iucn_threat(
    classification,
    output_dir      = extinction_dir,
    filename_prefix = paste0(organism, "_iucn_threat_classification"),
    organism_label  = organism_label
  )

  invisible(list(classification = classification, plot = p))
}
