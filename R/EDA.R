# prelim EDA


# packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggplot2)

files <- list.files("random_polygon_results/100_km")


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

species_sd <- data_500_km_all %>%
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

