# retained_cell_site_names.R — Parks, reserves and field stations at the
# grid cells highlighted in the P3 retained-cell map (fig_retained_cells_combined).
#
# One row per grid_id (EPSG:3857 cells, ~22 km wide). Where a biological
# station sits inside a park/reserve the label names both ("Park - Station").
# Used by label_retained_cells() (utils/label_retained_cells.R); cells without
# a row there fall back to "<country> (cell <id>)".
#
#   set         why the cell is logged
#                 richness_top  top plant completeness x richness cells
#                 country_top   highest plant completeness x richness cell in a country
#                 retained_top  most-retained cells (>= 5 organisms)
#                 low_plant     retained by >= 3 organisms but poor plant sampling
#   confidence  how the label was established
#                 user          confirmed by the user on satellite imagery
#                 high/medium/low  matched from coordinates + web/own knowledge, unchecked
#                 unidentified  no site identified yet
#
# Labels are cell-level best matches, not exhaustive lists of what a cell holds.

RETAINED_CELL_SITES <- tibble::tribble(
  ~grid_id, ~set,            ~site_name,                                                              ~confidence,
  # --- Top plant completeness x richness ---------------------------------------
  21604,    "richness_top",  "Saul, Parc amazonien de Guyane",                                        "medium",
  17793,    "richness_top",  "Yasuni National Park, Napo River / Pompeya road",                       "medium",
  17617,    "richness_top",  "Yasuni National Park - PUCE Yasuni Scientific Station",                 "high",
  6916,     "richness_top",  "Tambopata National Reserve - Explorer's Inn / Posada Amazonas",         "medium",
  15008,    "richness_top",  "Amacayacu National Park",                                               "medium",
  18322,    "richness_top",  "Cuyabeno Wildlife Reserve",                                             "medium",
  16292,    "richness_top",  "Rio Urubu State Forest",                                                "user",
  17442,    "richness_top",  "Waorani territory / Yasuni National Park",                              "user",
  17792,    "richness_top",  "Western Yasuni, Napo River lodges",                                     "user",
  7264,     "richness_top",  "Los Amigos Conservation Concession - CICRA Biological Station",         "high",
  16290,    "richness_top",  "Western BDFFP, BR-174 / Ramal do Rio Branquinho",                       "user",
  6724,     "richness_top",  "Machu Picchu Historic Sanctuary",                                       "user",
  15347,    "richness_top",  "Confluence of Amazon and Napo Rivers",                                  "user",
  # --- Top cell in each remaining country (not yet identified) ------------------------
  2191,     "country_top",   "Bosque de Helechos Gigantes",                                           "user",
  19640,    "country_top",   "Waiwoi Village, Rupununi",                                              "user",
  21766,    "country_top",   "Central Suriname Nature Reserve",                                       "user",
  19077,    "country_top",   "Cerro de la Neblina",                                                   "user",
  # --- Most-retained cells ---------------------------------------------------------
  9975,     "retained_top",  "Cristalino State Park - Cristalino Jungle Lodge / RPPN",                "high",
  22489,    "retained_top",  "Kaw-Roura Marshland Reserve - Tresor Regional Nature Reserve",          "high",
  22664,    "retained_top",  "Montsinery, D5 road",                                                   "user",
  22839,    "retained_top",  "Montagne des Singes, Kourou",                                           "medium",
  7786,     "retained_top",  "Manu National Park - Cocha Cashu Biological Station",                   "high",
  14465,    "retained_top",  "Tamshiyacu-Tahuayo Regional Conservation Area",                         "medium",
  17435,    "retained_top",  "Jatun Sacha Biological Station, Napo River",                            "medium",
  17968,    "retained_top",  "Limoncocha Biological Reserve",                                         "medium",
  21959,    "retained_top",  "Nouragues Nature Reserve - Nouragues Research Station",                 "high",
  # --- Retained but poorly sampled plants ----------------------------------------------
  21678,    "low_plant",     "Rio Guamal",                                                            "user",
  7259,     "low_plant",     "Amarakaeri Communal Reserve",                                           "user",
  7083,     "low_plant",     "Amarakaeri Communal Reserve (west)",                                    "user",
  7437,     "low_plant",     "Madre de Dios River lodges near Manu",                                  "user",
  7445,     "low_plant",     "Finca Las Piedras Research Station",                                    "user",
  9588,     "low_plant",     "Balneario Batistao",                                                    "user",
  19730,    "low_plant",     "Finca Bella Vista Eduardo Bravo, La Drymophila, El Cedral, Los Ariscos, Hostal de la Montana", "user",
  12862,    "low_plant",     "Rural area near Casa Alto Rioja",                                       "user",
  18673,    "low_plant",     "Colombia-Ecuador border (unidentified)",                                "unidentified"
)
