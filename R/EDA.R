# prelim EDA
# packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(readr)
library(sf)
library(ggridges)
library(forcats)
library(patchwork)

files <- list.files("random_polygon_results/500_km")




# summarize data function
aggregate_data_500 <- function(file_name){
  
  dat <- readRDS(paste0("random_polygon_results/500_km/", file_name)) %>%
    mutate(sample=gsub("random_sample_", "", file_name)) %>%
    mutate(sample=gsub(".RDS", "", sample))
  
  return(dat)
  
}

data_500_km_all <- bind_rows(lapply(files, aggregate_data_500))

length(unique(data_500_km_all$sample))

species_per_sample <- data_500_km_all %>%
  group_by(sample) %>%
  summarize(number_species=length(unique(COMMON_NAME)))

# Add lat/lng to the dataframe
sample_lat_lng <- readRDS("Data/potential_sample_points.RDS") %>%
  mutate(lat=st_coordinates(.)[,2]) %>%
  mutate(lng=st_coordinates(.)[,1]) %>%
  st_set_geometry(NULL) %>%
  rename(sample=ID) %>%
  mutate(sample=as.character(as.integer(sample)))

data_500_km_all <- data_500_km_all %>%
  left_join(., sample_lat_lng)

# Make a couple of big plots
ggplot(data_500_km_all, aes(x=total_median_viirs, y=UT_median_mean))+
  geom_point()+
  theme_bw()


data_500_km_all %>%
  group_by(sample) %>%
  summarize(number_species=n(),
            total_median_viirs=mean(total_median_viirs),
            positive_species=sum(UT_median>0),
            negative_species=sum(UT_median<0)) %>%
  mutate(proportion_negative=negative_species/number_species) %>%
  ggplot(., aes(x=total_median_viirs, y=proportion_negative))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))

ggplot(data_500_km_all, aes(x=UT_median, y=UT_median_mean))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()


species_sd <- data_500_km_all %>%
  group_by(COMMON_NAME) %>%
  summarize(N=n(),
            sd_urban=sd(UT_median),
            mean_urban=mean(UT_median),
            sd_total=sd(total_median_viirs),
            range_urban=max(UT_median)-min(UT_median))

# do species that are, on average across all buffers, selecting towards urban habitat
# versus those that are selecting against it have more variability in their urban scores among
# regions?
species_sd %>%
  mutate(association=ifelse(mean_urban>0, "positive", "negative")) %>%
  ggplot(., aes(x=association, y=sd_urban, fill=association))+
  geom_violin()+
  coord_flip()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  scale_fill_brewer(palette="Dark2")

# now look at the sd of the buffers
# a species is found in
# does that correlate with the sd of the species
sd_buffers <- species_sd %>%
  mutate(testing=sd_urban-sd_total)

ggplot(species_sd, aes(y=sd_urban, x=sd_total))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  scale_x_log10()+
  scale_y_log10()+
  ylab("Species-specific SD of urban scores")+
  xlab("Buffer SD of urban sampling")+
  geom_smooth(method="lm")+
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed")

ggplot(species_sd, aes(y=sd_urban, x=sd_total, size=N))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  scale_x_log10()+
  scale_y_log10()+
  ylab("Species-specific SD of urban scores")+
  xlab("Buffer SD of urban sampling")+
  geom_smooth(method="lm")+
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed")


mod <- lm(log10(sd_urban) ~ log10(sd_total) +log10(N), data=species_sd)
summary(mod)

resids <- species_sd %>%
  mutate(residual=resid(mod))


ggplot(species_sd, aes(x=sd_urban, y=mean_urban))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  #scale_x_log10()+
  #scale_y_log10()+
  xlab("Species-specific SD of urban scores")+
  ylab("Mean urban score of all buffers")+
  geom_smooth(method="lm")

ggplot(species_sd, aes(x=N, y=sd_urban))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  scale_x_log10()+
  scale_y_log10()+
  xlab("Number of buffers")+
  ylab("Species-specific SD of urban scores")+
  geom_smooth(method="lm")

# Try getting their rank for each sample for every species
test <- data_500_km_all %>%
  mutate(adjusted_viirs=median_viirs-total_median_viirs) %>%
  ungroup() %>%
  group_by(sample) %>%
  arrange(sample, desc(adjusted_viirs)) %>%
  mutate(rank=)

# try CV
# look at range of differences
# look at whether positive and negative mixes a lot



# Pick some species
# and make some ggridges figures where the ridges are high
# medium
# low ranked as urbanization

# read in data necessary to do the sampling for example species
potential_points <- readRDS("Data/potential_sample_points.RDS")

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



species_ex_function <- function(species_name){
  
  summary <- data_500_km_all %>%
    dplyr::filter(COMMON_NAME==species_name)
  
  # now read in data for the species
  # from those samples where it occurs
  get_data_function <- function(sample_number, grain_size=500000){
    
    # filter to point
    point <- potential_points %>%
      dplyr::filter(ID==sample_number)
    
    point2 <- point %>%
      st_transform(crs=4326)
    
    # create a random buffer that is specified by the grain size
    buff <- point %>%
      st_buffer(grain_size)
    
    buff2 <- buff %>%
      st_transform(crs=st_crs(checklists_sf))
    
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
      left_join(ebird_data, by="SAMPLING_EVENT_IDENTIFIER") %>%
      dplyr::filter(COMMON_NAME==species_name) %>%
      mutate(sample=sample_number)
    
    }
  
  dat <- bind_rows(lapply(unique(summary$sample), get_data_function))
  
  final_sp_dat <- dat %>%
    left_join(., summary %>%
                ungroup() %>%
                dplyr::select(sample, total_median_viirs)) %>%
    mutate(urban_category=case_when(total_median_viirs >= quantile(summary$total_median_viirs, 0.66) ~ "High",
                                    total_median_viirs < quantile(summary$total_median_viirs, 0.66) & total_median_viirs >= quantile(summary$total_median_viirs, 0.33) ~ "Medium",
                                    total_median_viirs < quantile(summary$total_median_viirs, 0.33) ~ "Low"))
  
  return(final_sp_dat)
  
}


# do this for four species for now
ex_sp_list <- c("Canada Warbler", "White-winged Dove",
                "Northern Cardinal", "Roseate Spoonbill")


example_sp_dat <- bind_rows(lapply(ex_sp_list, species_ex_function))

cawa <- example_sp_dat %>%
  dplyr::filter(COMMON_NAME=="Canada Warbler") %>%
  ggplot(., aes(y=factor(urban_category, levels=c("Low", "Medium", "High")), x=viirs))+
  geom_density_ridges()+
  scale_x_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Urbanization (VIIRS night-time lights)")+
  ylab("Categorical urbanization of buffer")+
  facet_wrap(~COMMON_NAME)

cawa

rosp <- example_sp_dat %>%
  dplyr::filter(COMMON_NAME=="Roseate Spoonbill") %>%
  ggplot(., aes(y=factor(urban_category, levels=c("Low", "Medium", "High")), x=viirs))+
  geom_density_ridges()+
  scale_x_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("Urbanization (VIIRS night-time lights)")+
  ylab("")+
  facet_wrap(~COMMON_NAME)

rosp

wwdo <- example_sp_dat %>%
  dplyr::filter(COMMON_NAME=="White-winged Dove") %>%
  ggplot(., aes(y=factor(urban_category, levels=c("Low", "Medium", "High")), x=viirs))+
  geom_density_ridges()+
  scale_x_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("")+
  ylab("")+
  facet_wrap(~COMMON_NAME)

wwdo

noca <- example_sp_dat %>%
  dplyr::filter(COMMON_NAME=="Northern Cardinal") %>%
  ggplot(., aes(y=factor(urban_category, levels=c("Low", "Medium", "High")), x=viirs))+
  geom_density_ridges()+
  scale_x_log10()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  xlab("")+
  ylab("Categorical urbanization of buffer")+
  facet_wrap(~COMMON_NAME)

noca



noca + wwdo + cawa + rosp + plot_layout(ncol=2)





# Some other exploration
data_500_km_all %>%
  ungroup() %>%
  dplyr::filter(COMMON_NAME %in% ex_sp_list) %>%
  ggplot(., aes(y=UT_median, x=total_median_viirs))+
  geom_point()+
  facet_wrap(~COMMON_NAME, scales="free")+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  geom_smooth(method="lm")



# run a simple lm and get the parameter estimate
# between UT median and the total median VIIRS for the buffers that a species
# occurs in

species_model_function <- function(species_name){
  
  tmp <- data_500_km_all %>%
    dplyr::filter(COMMON_NAME==species_name)
  
  
  mod <- lm(UT_median ~ total_median_viirs, data=tmp)
  mod2 <- mgcv::gam(UT_median ~ total_median_viirs + s(lng, lat, k=4), data=tmp)
  
  summary <- broom::tidy(mod2, parametric=TRUE) %>%
    mutate(N=nrow(tmp)) %>%
    mutate(prop_negative=sum(tmp$UT_median<0)/nrow(tmp)) %>%
    mutate(COMMON_NAME=species_name)
   
  return(summary)
  
}

species_results <- bind_rows(lapply(species_sd %>%
                                      dplyr::filter(N>10) %>%
                                      .$COMMON_NAME, species_model_function))



species_results %>%
  dplyr::filter(term=="total_median_viirs") %>%
  ggplot(., aes(y=estimate, x=prop_negative))+
  geom_point()+
  geom_smooth(method="lm")




























###############################################################
###############################################################
# OLD STUFF
# summarize data function
aggregate_data_100 <- function(file_name){
  
  dat <- readRDS(paste0("random_polygon_results/100_km/", file_name)) %>%
    mutate(sample=gsub("random_sample_", "", file_name)) %>%
    mutate(sample=gsub(".RDS", "", sample))
  
  return(dat)
  
}

data_100_km_all <- bind_rows(lapply(files, aggregate_data_100))

length(unique(data_100_km_all$sample))

species_per_sample <- data_100_km_all %>%
  group_by(sample) %>%
  summarize(number_species=length(unique(COMMON_NAME)))

species_sd <- data_100_km_all %>%
  mutate(adjusted_viirs=median_viirs-total_median_viirs) %>%
  group_by(COMMON_NAME) %>%
  summarize(N=n(),
            sd_urban=sd(adjusted_viirs),
            mean_urban=mean(adjusted_viirs),
            sd_total=sd(total_median_viirs))

ggplot(species_sd, aes(x=sd_urban, y=sd_total))+
  geom_point()+
  theme_bw()+
  theme(axis.text=element_text(color="black"))+
  scale_x_log10()+
  scale_y_log10()+
  xlab("Species-specific SD of urban scores")+
  ylab("Buffer SD of urban sampling")+
  geom_smooth(method="lm")


# any samples with at least 10 species
samples_to_test <- species_per_sample %>%
  dplyr::filter(number_species>=10)

# what about pairwise R2 of every sample
# with more than 10 species (10 is arbitrary for now)
pairwise_r2_function <- function(sample_number){
  
  xxx <- data_100_km_all %>%
    dplyr::filter(sample==sample_number) %>%
    mutate(adjusted_viirs=median_viirs-total_median_viirs)
  
  # get R2
  get_r2 <- function(sample_comparison){
    
    zzz <- data_100_km_all %>%
      dplyr::filter(sample==sample_comparison) %>%
      mutate(adjusted_viirs=median_viirs-total_median_viirs) %>%
      dplyr::select(COMMON_NAME, adjusted_viirs) %>%
      left_join(., xxx %>%
                  dplyr::select(COMMON_NAME, adjusted_viirs), by="COMMON_NAME") %>%
      ungroup() %>%
      dplyr::filter(complete.cases(.))
    
    summary_df <- data.frame(r2=summary(lm(adjusted_viirs.x ~ adjusted_viirs.y, data=zzz))$r.squared) %>%
      mutate(number_species=nrow(zzz)) %>%
      mutate(comparison=paste0(sample_number, "_vs_", sample_comparison))
    
  }
  
  r2s <- bind_rows(lapply(samples_to_test %>%
                            dplyr::filter(sample != sample_number) %>%
                            .$sample, get_r2))
  
}

# lots of comparisons!
lapply_with_error <- function(X,FUN,...){    
  lapply(X, function(x, ...) tryCatch(FUN(x, ...),
                                      error=function(e) NULL))
}

r2_comparisons <- bind_rows(lapply_with_error(samples_to_test$sample, pairwise_r2_function))

temp <- r2_comparisons %>%
  dplyr::filter(number_species>=5)

hist(temp$r2)



