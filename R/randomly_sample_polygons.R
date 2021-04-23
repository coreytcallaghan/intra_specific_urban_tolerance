## This scripts is used to randomly sample
## A large buffer
## and then get all eBird checklists within that buffer

# packages
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(rnaturalearth)
library(urbnmapr)

# read in checklists with viirs
checklists <- read_csv("Data/checklists_viirs_scores/ebird_samples_viirs_scores.csv") %>%
  dplyr::select(2, 3) %>%
  rename(viirs=first) %>%
  rename(SAMPLING_EVENT_IDENTIFIER=SAMPLIN)

# read in eBird dataset
ebird_data <- readRDS("Data/ebird_data_raw_May.RDS") %>% 
  bind_rows(readRDS("Data/ebird_data_raw_Jun.RDS")) %>%
  bind_rows(readRDS("Data/ebird_data_raw_Jul.RDS")) %>%
  bind_rows(readRDS("Data/ebird_data_raw_Aug.RDS")) %>%
  left_join(., read_csv("Data/checklists_mod_scores/ebird_samples_mod_scores.csv") %>%
              dplyr::select(first, SAMPLIN) %>%
              rename(ghm=first) %>%
              rename(SAMPLING_EVENT_IDENTIFIER=SAMPLIN), by="SAMPLING_EVENT_IDENTIFIER") %>%
  left_join(., read_csv("Data/Clements-Checklist-v2019-August-2019.csv") %>%
              dplyr::filter(category=="species") %>%
              dplyr::select(category, `English name`, `scientific name`, order, family) %>%
              rename(COMMON_NAME=`English name`,
                     SCIENTIFIC_NAME=`scientific name`)) %>%
  dplyr::filter(!family %in% c("Strigidae (Owls)", "Tytonidae (Barn-Owls)",
                               "Stercorariidae (Skuas and Jaegers)", "Alcidae (Auks, Murres, and Puffins)",
                               "Sulidae (Boobies and Gannets)", "Procellariidae (Shearwaters and Petrels)",
                               "Hydrobatidae (Northern Storm-Petrels)", "Oceanitidae (Southern Storm-Petrels)")) %>%
  dplyr::filter(complete.cases(BCR_CODE))

# create a checklists sf object for spatial stuff below
checklists_sf <- checklists %>%
  left_join(ebird_data %>%
              dplyr::select(SAMPLING_EVENT_IDENTIFIER, LONGITUDE, LATITUDE) %>%
              distinct()) %>%
  dplyr::filter(complete.cases(.)) %>%
  st_as_sf(coords=c("LONGITUDE", "LATITUDE"), crs=4326)

# quick plot of viirs spatial data to check it makes sense
ggplot()+
  geom_sf(data=checklists_sf, aes(color=log10(viirs)))+
  scale_color_viridis_c()

# read in US shapefile
states_sf <- get_urbn_map("states", sf = TRUE) %>%
  dplyr::filter(state_name != "Alaska") %>%
  dplyr::filter(state_name != "Hawaii")

states_sf %>% 
  ggplot(aes()) +
  geom_sf(fill = "grey", color = "#ffffff")

st_crs(states_sf)

# Now my goal is write a function
# that randomly creates a polygon
# and then selects all eBird checklists within that polygon
# correlates this with the data
# summarizes the median scores and adjusted scores of all species
# within that polygon
# and also summarizes the number of observations for each species
# then writes out this summary somehow

# first create a random sample of 10,000 points in the us
potential_points <- st_sample(states_sf, 10000) %>%
  st_as_sf() %>%
  mutate(ID=1:nrow(.))

states_sf %>% 
  ggplot(aes()) +
  geom_sf(fill = "grey", color = "#ffffff")+
  geom_sf(data=potential_points)

random_sample_function <- function(id_number, grain_size){
  
  message(paste0("Analyzing poin id number: ", id_number))
  
  # filter to point
  point <- potential_points %>%
    dplyr::filter(ID==id_number)
  
  point2 <- point %>%
    st_transform(crs=4326)
  
  # create a random buffer that is specified by the grain size
  buff <- point %>%
    st_buffer(grain_size)
  
  # plot test
  states_sf %>% 
    ggplot(aes()) +
    geom_sf(fill = "grey", color = "#ffffff")+
    geom_sf(data=buff)+
    geom_sf(data=point)
  
  buff2 <- buff %>%
    st_transform(crs=st_crs(checklists_sf))
  
  states_sf %>% 
    ggplot(aes()) +
    geom_sf(fill = "grey", color = "#ffffff")+
    geom_sf(data=buff2)+
    geom_sf(data=point)
  
  # now get all eBird checklists that fall within that buffer
  lists <- checklists_sf %>%
    st_intersects(buff2) %>%
    as.data.frame()
  
  buff_ebird_dat <- checklists_sf %>%
    as.data.frame() %>%
    dplyr::select(SAMPLING_EVENT_IDENTIFIER, viirs) %>%
    mutate(row.id=1:nrow(.)) %>%
    right_join(., lists, by="row.id") %>%
    dplyr::select(-row.id, -col.id) %>%
    left_join(ebird_data, by="SAMPLING_EVENT_IDENTIFIER")
  
  species_summary <- buff_ebird_dat %>%
    dplyr::filter(complete.cases(.)) %>%
    group_by(COMMON_NAME, SCIENTIFIC_NAME) %>%
    summarize(obs=n(),
              mean_viirs=mean(viirs),
              median_viirs=median(viirs),
              sd_viirs=sd(viirs),
              max_viirs=max(viirs)) %>%
    dplyr::filter(obs>=100)
  
  all_samples_summary <- buff_ebird_dat %>%
    dplyr::select(SAMPLING_EVENT_IDENTIFIER, viirs) %>%
    distinct() %>%
    dplyr::filter(complete.cases(.)) %>%
    summarize(total_obs=n(),
              total_mean_viirs=mean(viirs),
              total_median_viirs=median(viirs),
              total_sd_viirs=sd(viirs),
              total_max_viirs=max(viirs))
  
  final_summary <- bind_cols(species_summary, all_samples_summary) %>%
    mutate(point_ID=id_number) %>%
    mutate(LATITUDE=st_coordinates(point2)[2]) %>%
    mutate(LONGITUDE=st_coordinates(point2[1]))
  
  size=grain_size/1000
  
  saveRDS(final_summary, paste0("random_polygon_results/", size, "_km/random_sample_", id_number, ".RDS"))
  
}

lapply(c(1:1000), function(x){random_sample_function(x, 500000)})
lapply(c(1001:2000), function(x){random_sample_function(x, 500000)})
lapply(c(1:1000), function(x){random_sample_function(x, 100000)})
lapply(c(1001:2000), function(x){random_sample_function(x, 100000)})
