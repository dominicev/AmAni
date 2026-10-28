##set the working directory
library(knitr)
library(rstudioapi)
knitr::opts_knit$set(root.dir = dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##load other dependencies

library(data.table)
setDTthreads(0)   # 0 = use all available cores
sessionInfo()

library(sf)
library(utils)
library(lwgeom)
library(ggplot2)
library(terra)
library(tidyr)  # For replace_na
library(readr)
library(tidyr)
library(stringr)
library(parzer)
library(MASS)
library(dplyr)
library(readxl)
library(vegan)
library(scales)
library(WorldFlora)
library(readr)
library(svglite)
library(rnaturalearth)

#library(taxonlookup)
#library(installr)
#updateR()


###set directory for loading data
dataDir<-paste0(dirname(getwd()),"/0-CuratingPlantData")

###load species occurences
species_counts<-readRDS(paste0(dataDir,"/species_counts3.rds"))

###load weighted endemism
species_endemism<-readRDS(paste0(dataDir,"/species_endemism2.rds"))

#make spatial grids

grid_sf<-species_counts
#names(grid_sf)[1]<-"grid_id"



#grid_sf <- st_join(grid_sf, species_endemism, by = "grid_id") #slower version of above

#hist(grid_sf$weighted_endemism) #this takes a few minutes
grid_sf_SE <- grid_sf |>
  left_join(
    st_drop_geometry(species_endemism),
    by = "grid_id"
  )

#check for outliers in weighted endemism
# 1) Compute cutoff (extremely high outliers only)
q95  <- quantile(grid_sf_SE$weighted_endemism, 0.95, na.rm = TRUE)
iqr <- IQR(grid_sf_SE$weighted_endemism, na.rm = TRUE)
upper_cutoff <- q95 + 3 * iqr ##remove all values above 367

# 2) Create a filtered column
grid_sf_SE <- grid_sf_SE |>
  dplyr::mutate(weight_endemism_filt = dplyr::if_else(
    weighted_endemism > upper_cutoff, NA_real_, as.numeric(weighted_endemism)))

##save temp output
saveRDS(grid_sf_SE, "filteredWeighted_species_endemism.rds")
weighted_endemism_filt<-readRDS("filteredWeighted_species_endemism.rds")


# Raw endemism: count of species that are endemic (i.e., occur in only one grid cell)
# First, identify endemic species
endemics <- species_counts %>%
  group_by(scientificName) %>%
  filter(n_distinct(grid_id) == 1) %>%
  ungroup()

# Then, count number of endemic species per grid cell
raw_endemism <- endemics %>%
  group_by(grid_id) %>%
  summarise(raw_endemism = n_distinct(scientificName), .groups = "drop")

saveRDS(raw_endemism, "raw_species_endemism2.rds")
raw_endemism<-readRDS("raw_species_endemism2.rds")

#names(grid_sf)[1]<-"grid_id"

grid_sf_RE <- grid_sf |>
  left_join(
    st_drop_geometry(raw_endemism),
    by = "grid_id"
  )

#grid_sf_RE <- st_join(grid_sf, raw_endemism, by = "grid_id") #slow version of above

#Check for outliers in the endemism data
# 1) Compute cutoff (extremely high outliers only)
q95  <- quantile(grid_sf_RE$raw_endemism, 0.95, na.rm = TRUE)
iqr <- IQR(grid_sf_RE$raw_endemism, na.rm = TRUE)
upper_cutoff <- q95 + 3 * iqr ##remove all values above 172

# 2) Create a filtered column
grid_sf_RE <- grid_sf_RE |>
  dplyr::mutate(raw_endemism_filt = dplyr::if_else(
    raw_endemism > upper_cutoff, NA_real_, as.numeric(raw_endemism)
  ))

#save temporary output
saveRDS(grid_sf_RE, "filteredRaw_species_endemism.rds")
raw_endemism_filt<-readRDS("filteredRaw_species_endemism.rds")

##Now, remove the grid cells in the occurence data that we filtered from the endemism data
filtered_ids1 <- grid_sf_RE |>
  dplyr::filter(raw_endemism > upper_cutoff) |>
  dplyr::pull(grid_id)
union(filtered_ids1,filtered_ids1)

filtered_ids2 <- grid_sf_SE |>
  dplyr::filter(weighted_endemism > upper_cutoff) |>
  dplyr::pull(grid_id)
union(filtered_ids1,filtered_ids2)


###also, i should remove all grid cells that correspond to cities, as many of the plants there are likely cultivated or, in the worst case scenario, botanical gardens.
###now find the grid cells that correspond to cities with a population >10,000 people

##get bounding box for my data

grid_sf_clean <- st_as_sf(
  grid_sf_SE,
  sf_column_name = attr(grid_sf_SE, "sf_column"),
  crs = st_crs(grid_sf_SE)
)
grid_sf_clean <- st_cast(grid_sf_SE, "POINT", warn = FALSE)
grid_ll <- st_transform(grid_sf_clean, 4326)

bbox <- st_bbox(grid_ll)
min_lon <- bbox["xmin"]
max_lon <- bbox["xmax"]
min_lat <- bbox["ymin"]
max_lat <- bbox["ymax"]

###import city data

cities_sf <- ne_download(
  scale = 10,               # highest resolution
  type = "populated_places",
  category = "cultural",
  returnclass = "sf"
)

sort(colnames(cities_sf))
#filter cities list by min max lat longs
coords <- st_coordinates(cities_sf)

cities_crop <- cities_sf |>
  dplyr::filter(
    coords[,1] >= min_lon &
      coords[,1] <= max_lon &
      coords[,2] >= min_lat &
      coords[,2] <= max_lat
  )

sort(colnames(cities_crop))
cities_crop$ADM1NAME

#now get a list of only cities with pops >10,000
cities_SA_large <- cities_crop |>
  dplyr::filter(
    POP_MAX >= 10000)

#take these cities, find the corresponding grid ID in my data, and add those to the list of bad grid IDs
grid_sf <- st_as_sf(
  grid_sf,
  sf_column_name = "geometry",
  crs = 3857
)

cities_big <- st_transform(cities_SA_large, st_crs(grid_sf))

#get the radius of each city
cities_big <- cities_big |>
  dplyr::mutate(
    buffer_radius_m = sqrt(MAX_AREAKM / pi) * 1000
  )
summary(cities_big$buffer_radius_m)

#create a polygon using the radius of the city and its point location
cities_poly <- st_buffer(
  cities_big,
  dist = cities_big$buffer_radius_m
)

cities_poly <- st_simplify(cities_poly, dTolerance = 500)

#now find the grifd_ids in our data whose point falls within the city polygon
grid_ids_cities <- grid_sf |>
  st_join(cities_poly, join = st_intersects) |>
  dplyr::filter(!is.na(NAME)) |>
  dplyr::pull(grid_id) |>
  unique()


length(grid_ids_cities)
length(cities_poly)




length(union(filtered_ids1, filtered_ids2))
length(union(union(filtered_ids1, filtered_ids2), cities_poly))
length(grid_ids_cities)+length(union(filtered_ids1, filtered_ids2))

#based on he above lengths, it seems that all the outliers were already cities! Good that we're removing the rest of them.

bad_grid_ids <- union(union(filtered_ids1, filtered_ids2), grid_ids_cities)

#this lists all the grid cells we want to remove

setDT(species_counts)

species_counts_filt <- species_counts[!grid_id %in% bad_grid_ids]
nrow(species_counts) - nrow(species_counts_filt)

inherits(species_counts_filt, "sf")
species_counts_filt <- st_as_sf(species_counts_filt)

#count how many species records were in the grid cells we removed
species_counts |>
  dplyr::filter(grid_id %in% bad_grid_ids) |>
  dplyr::count(grid_id, sort = TRUE)

saveRDS(species_counts_filt,"species_counts4_filtered.rds")


##since the same grid cells weren't removed from the weighted endemism and the raw endemism set, we should do the same thing as for species occurences, so all the datasets have exactly the same grid cells


#grid_sf_RE
setDT(grid_sf_RE)
grid_sf_RE_filteredOuliers <- grid_sf_RE[!grid_id %in% bad_grid_ids]
inherits(grid_sf_RE_filteredOuliers, "sf")
grid_sf_RE_filteredOuliers <- st_as_sf(grid_sf_RE_filteredOuliers)
nrow(grid_sf_RE)
nrow(grid_sf_RE_filteredOuliers)

geom <- st_geometry(grid_sf_RE_filteredOuliers)
attrs <- as.data.frame(grid_sf_RE_filteredOuliers)
attrs$geometry <- NULL

grid_sf_RE_filteredOuliers <- st_sf(
  attrs,
  geometry = geom,
  crs = 3857
)




#grid_sf_SE
setDT(grid_sf_SE)
grid_sf_SE_filteredOuliers <- grid_sf_SE[!grid_id %in% bad_grid_ids]
inherits(grid_sf_SE_filteredOuliers, "sf")
grid_sf_SE_filteredOuliers <- st_as_sf(grid_sf_SE_filteredOuliers)

geom <- st_geometry(grid_sf_SE_filteredOuliers)
attrs <- as.data.frame(grid_sf_SE_filteredOuliers)
attrs$geometry <- NULL

grid_sf_SE_filteredOuliers <- st_sf(
  attrs,
  geometry = geom,
  crs = 3857
)


saveRDS(grid_sf_RE_filteredOuliers, "filteredRaw_species_endemism_noCities.rds")
raw_endemism_filt_BEST<-readRDS("filteredRaw_species_endemism_noCities.rds")

saveRDS(grid_sf_SE_filteredOuliers, "filteredWeighted_species_endemism_noCities.rds")
weighted_endemism_filt_BEST<-readRDS("filteredWeighted_species_endemism_noCities.rds")



temp<-weighted_endemism_filt_BEST |>
  dplyr::group_by(grid_id) |>
  dplyr::slice(1) |>
  dplyr::ungroup()

all(weighted_endemism_filt_BEST$weight_endemism_filt == weighted_endemism_filt_BEST$weight_endemism_filt)

new<-weighted_endemism_filt_BEST[,!(names(weighted_endemism_filt_BEST) %in% c("weight_endemism_filt", "scientificName", "count"))]
new2<-raw_endemism_filt_BEST[,!(names(raw_endemism_filt_BEST) %in% c("raw_endemism_filt",  "scientificName",  "count"))]

new_unique <- new |>
  group_by(grid_id) |>
  slice(1) |>
  ungroup()

new_unique2 <- new2 |>
  group_by(grid_id) |>
  slice(1) |>
  ungroup()


allEndemismAndOccurenceFiltered <- new_unique |>
  left_join(
    st_drop_geometry(new_unique2),
    by = "grid_id"
  )

saveRDS(allEndemismAndOccurenceFiltered, "FINAL_allEndemism_FullyFiltered.rds")
allEndemismAndOccurenceFiltered<-readRDS("FINAL_allEndemism_FullyFiltered.rds")

####plot each endemism dataset


###weighted plot
class(weighted_endemism_filt_BEST)

weighted_endemism_filt<-c()
raw_endemism_filt<-c()

#remove multi-point geometries
weighted_endemism_filt_BEST <- st_cast(
  weighted_endemism_filt_BEST,
  "POINT",
  warn = FALSE
)

{
  
#remove redundant points
grid_plot <- weighted_endemism_filt_BEST |>
  dplyr::group_by(grid_id) |>
  dplyr::slice(1) |>
  dplyr::ungroup()

#plot
pWeight<-ggplot(grid_plot) +
  geom_sf(
    aes(color = weight_endemism_filt),
    size = 0.6
  ) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal()


ggsave(
  filename = "weighted_endemism.svg",
  plot = pWeight,
  device = "svg",
  width = 8,
  height = 6,
  units = "in"
)
}


###raw plot
class(raw_endemism_filt_BEST)


#remove multi-point geometries
raw_endemism_filt_BEST <- st_cast(
  raw_endemism_filt_BEST,
  "POINT",
  warn = FALSE
)

{
  #remove redundant points
  grid_plot <- raw_endemism_filt_BEST |>
    dplyr::group_by(grid_id) |>
    dplyr::slice(1) |>
    dplyr::ungroup()
  
  #plot
  pRaw<-ggplot(grid_plot) +
    geom_sf(
      aes(color = raw_endemism_filt),
      size = 0.6
    ) +
    scale_color_viridis_c(option = "plasma") +
    theme_minimal()
  
  
  ggsave(
    filename = "raw_endemism.svg",
    plot = pRaw,
    device = "svg",
    width = 8,
    height = 6,
    units = "in"
  )
}










