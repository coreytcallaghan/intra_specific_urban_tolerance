# packages
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(rnaturalearth)
library(urbnmapr)


sf_use_s2(FALSE)

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
# potential_points <- st_sample(states_sf, 10000) %>%
#   st_as_sf() %>%
#   mutate(ID=1:nrow(.))
# 
# saveRDS(potential_points, "Data/potential_sample_points.RDS")

potential_points <- readRDS("Data/potential_sample_points.RDS")

# Get some of the most well-sampled points
# summarize data function
aggregate_data_500 <- function(file_name){
  
  dat <- readRDS(paste0("random_polygon_results/500_km/", file_name)) %>%
    mutate(sample=gsub("random_sample_", "", file_name)) %>%
    mutate(sample=gsub(".RDS", "", sample))
  
  return(dat)
  
}

data_500_km_all <- bind_rows(lapply(files, aggregate_data_500))

temp <- data_500_km_all %>% 
  group_by(point_ID) %>% 
  summarize(N=sum(obs)) %>%
  arrange(desc(N)) %>%
  slice(1:200) %>%
  sample_n(10)

# to match our other script
temp <- data.frame(point_ID=c(527, 557, 634, 898, 1234, 1414, 1453, 1661, 1894, 2000))

assess_buffer_size_function <- function(id_number, grain_size){
  
  message(paste0("Analyzing point id number: ", id_number))
  
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
  
  # now for every species
  # get the area of a concave hull
  species_area <- function(species_name){
    
    species_obs <- buff_ebird_dat %>%
      dplyr::filter(COMMON_NAME==species_name)
    
    hull <- species_obs %>%
      st_as_sf(coords=c("LONGITUDE", "LATITUDE"), crs=4326) %>%
      concaveman::concaveman(.) %>%
      st_area()
    
    # summary df
    summary_df <- data.frame(area=as.numeric(hull),
                             COMMON_NAME=species_name)
    
    return(summary_df)
  }
  
  # do this for every species with >100 observations
  species_list <- buff_ebird_dat %>%
    group_by(COMMON_NAME) %>%
    summarize(N=n()) %>%
    dplyr::filter(N>=100)
  
  species_area_results <- bind_rows(lapply(species_list$COMMON_NAME, species_area)) %>%
    mutate(point_ID=id_number) %>%
    mutate(buffer_area=as.numeric(st_area(buff2)))
  
  return(species_area_results)
  
}

buffer_size_area_results_500000 <- bind_rows(lapply(temp$point_ID, function(x){assess_buffer_size_function(x, 500000)}))

buffer_size_area_results_v100000 <- bind_rows(lapply(temp$point_ID, function(x){assess_buffer_size_function(x, 100000)}))

buffer_size_area_results_v1000000 <- bind_rows(lapply(temp$point_ID, function(x){assess_buffer_size_function(x, 1000000)}))

buffer_size_area_results_v250000 <- bind_rows(lapply(temp$point_ID, function(x){assess_buffer_size_function(x, 250000)}))

buffer_size_assessment_all <- buffer_size_area_results_500000 %>%
  mutate(buffer_size_km=500) %>%
  bind_rows(buffer_size_area_results_v100000 %>%
              mutate(buffer_size_km=100)) %>%
  bind_rows(buffer_size_area_results_v1000000 %>%
              mutate(buffer_size_km=1000)) %>%
  bind_rows(buffer_size_area_results_v250000 %>%
              mutate(buffer_size_km=250))

buffer_size_assessment_all %>%
  group_by(buffer_size_km) %>%
  summarize(total_species_points=n(),
            total_unique_species=length(unique(COMMON_NAME)))

buffer_size_assessment_all %>%
  mutate(area_percent=(area/buffer_area)*100) %>%
  mutate(buffer_size_km=as.character(buffer_size_km)) %>%
  mutate(buffer_size_km=factor(buffer_size_km, levels=c("100", "250", "500", "1000"))) %>%
  ggplot(., aes(x=fct_inorder(buffer_size_km), y=area_percent))+
  geom_boxplot()+
  coord_flip()+
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Buffer radius size (km)")+
  ylab("Percent area overlap of a concave hull")+
  scale_x_discrete(breaks=c("100", "250", "500", "1000"), limits=c("100", "250", "500", "1000"), labels=c("100", "250", "500", "1000"))

