# utils/collect_pub_data.R
# Sourced by AmAni.summarize.R — relies on the packages (tidyverse, sf,
# patchwork, ggrepel) loaded there, on organism_label() and ORGANISM_GROUPS
# from utils/organism_label.R / utils/organism_groups.R, and on the four
# collectors (collect_extinction_summary, collect_organism_summary,
# collect_dataset_distributions, collect_iucn).

# Bundle everything the publication figures need into one named list:
#   $organisms  organism_dir + display label + coarse taxon group
#   $extinction species/population loss summaries (all models x scenarios)
#   $summary    per-organism filtering, gamma and sampling-extension stats
#   $dists      per-cell five-number summaries before/after filtering
#   $iucn       IUCN category counts (all models x scenarios)
#
# `include` filters to specific organism directory names; `group` adds a
# named group from ORGANISM_GROUPS (e.g. group = "vertebrates"). Like the
# rest of the project, the working directory is assumed to be the project
# root, so output_dir defaults to "output".

collect_pub_data <- function(include = NULL, group = NULL,
                             output_dir = "output") {
  if (!is.null(group)) {
    if (!group %in% names(ORGANISM_GROUPS))
      stop("Unknown group '", group, "'. Available: ",
           paste(names(ORGANISM_GROUPS), collapse = ", "), call. = FALSE)
    include <- union(include, ORGANISM_GROUPS[[group]])
  }

  dirs <- .organism_dirs(output_dir)
  if (!is.null(include)) {
    missing <- setdiff(include, dirs)
    if (length(missing) > 0)
      warning("No output directory for: ", paste(missing, collapse = ", "),
              call. = FALSE)
    dirs <- intersect(dirs, include)
  }
  if (length(dirs) == 0)
    stop("No organism output directories found in '", output_dir, "'.",
         call. = FALSE)

  # Coarse taxonomic grouping for colouring (Vertebrates vs everything else)
  is_vert <- dirs %in% ORGANISM_GROUPS$vertebrates

  list(
    organisms = tibble::tibble(
      organism_dir = dirs,
      organism     = organism_label(dirs),
      taxon_group  = ifelse(is_vert, "Vertebrates", "Invertebrates")
    ),
    extinction = collect_extinction_summary(dirs, output_dir),
    summary    = collect_organism_summary(dirs, output_dir),
    dists      = collect_dataset_distributions(dirs, output_dir),
    iucn       = collect_iucn(dirs, output_dir)
  )
}
