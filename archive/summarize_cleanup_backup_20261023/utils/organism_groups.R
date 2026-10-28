# utils/organism_groups.R
# Extracted from AmAni.summarize.R. Sourced by AmAni.summarize.R — relies on
# the packages (tidyverse, sf, patchwork, ggrepel) loaded there.

# Predefined organism groups for easy filtering in summarize_all().
# Usage: summarize_all(group = "insects")
#        summarize_all(group = "insects", results_name = "insects_summary")
# By default each run archives any existing output into
# <results_name>/_archive/<timestamp>/ before writing new files (archive = FALSE
# to disable).
ORGANISM_GROUPS <- list(
  insects     = c("cockroaches", "allRoaches", "beetles", "lepidoptera",
                  "hemiptera", "odonata", "butterfl-curated", "odonata-curated", "hymenoptera",
                  "bees", "flies", "DungBeetles"),
  comparisons = c("butterflies", "butterfl-curated", "odonata", "odonata-curated", "odonata-combined"),
  vertebrates = c("antbirds", "finches", "frogs", "mammals", "tyrannidae",
                  "squamates", "squamata"),
  birds       = c("antbirds", "finches", "tyrannidae"),
  herps       = c("frogs", "squamates", "squamata")
)
