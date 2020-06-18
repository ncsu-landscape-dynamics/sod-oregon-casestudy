### to make automanagement R data file 

library(PoPS)
library(raster)
library(doParallel)
library(iterators)
library(foreach)
library(parallel)
library(lubridate)

source('/rpops/R/helpers.R')
source('/rpops/R/checks.R')
source('/rpops/R/uncertainty_propogation.R')

species <- c("SOD_EU1", "SOD_NA1")
infected_files <- c("/input data/end_inf_2019eu.tif", "/input data/end_inf_2019.tif")
host_file <- "/input data/lide_100m_median_2019.tif"
total_plants_file <- "/input data/lemma_max100m.tif"
temperature_file <- ""
temperature_coefficient_file <- #"/average_weather_2020_2024.tif"
  ## NOTE: weather_coef_2020_2024 file too large to attach on github
  ## If you are a member of the Center for Geospatial Analytics, you can access the file here:
  ## https://drive.google.com/drive/u/0/folders/1Ho1C1vFjechtY5zZVrfDrDEtPlVBiGYO
  ## Otherwise, please reach out to Devon Gaydos at devon.gaydos@usda.gov or to Chris Jones at cmjone25@ncsu.edu
  
precipitation_coefficient_file <-""
use_lethal_temperature <- FALSE
temp <- TRUE
precip <- FALSE
season_month_start <- 1
season_month_end <- 12
time_step <- "week"
start_date = "2020-01-01"
end_date = "2024-12-31"
lethal_temperature <- -35
lethal_temperature_month <- as.integer(1)
treatments_file <- ""
treatment_years <- c(0)
treatment_dates <- "2024-12-01"
treatment_method <- "all infected"
management <- FALSE
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
reproductive_rate <- c(1.6, 0.6)
percent_natural_dispersal <- c(1.0, 1.0)
natural_kernel_type <- c("exponential", "exponential")
anthropogenic_kernel_type <- c("cauchy", "cauchy")
natural_distance_scale <- c(242, 354)
anthropogenic_distance_scale <- c(0.0, 0.0)
natural_dir <- c("N", "N")
natural_kappa <- c(2, 2)
anthropogenic_dir <- c("NONE", "NONE")
anthropogenic_kappa <- c(0, 0)
pesticide_duration <- c(0)
pesticide_efficacy <- 1.0
random_seed <- 42
output_frequency = "year"
movements_file = ""
use_movements = FALSE



treatment_metric_check <- treatment_metric_checks(treatment_method)
if (!treatment_metric_check$checks_passed) {
  return(treatment_metric_check$failed_check)
}

time_check <- time_checks(end_date, start_date, time_step, output_frequency)
if(time_check$checks_passed) {
  number_of_time_steps <- time_check$number_of_time_steps
  number_of_years <- time_check$number_of_years
  number_of_outputs <- time_check$number_of_outputs
} else {
  return(time_check$failed_check)
}

percent_check <- percent_checks(percent_natural_dispersal)
if (percent_check$checks_passed){
  use_anthropogenic_kernel <- percent_check$use_anthropogenic_kernel
} else {
  return(percent_check$failed_check)
}

multispecies_check <- multispecies_checks(species, infected_files, reproductive_rate, percent_natural_dispersal, natural_kernel_type, anthropogenic_kernel_type, 
                                          natural_distance_scale, anthropogenic_distance_scale, natural_dir, natural_kappa, anthropogenic_dir, anthropogenic_kappa)
if (!multispecies_check$checks_passed){
  return(percent_check$failed_check)
}

infected_check <- initial_raster_checks(infected_files)
if (infected_check$checks_passed) {
  infected <- infected_check$raster
  # if (raster::nlayers(infected) > 1) {
  #   infected <- output_from_raster_mean_and_sd(infected)
  # }
} else {
  return(infected_check$failed_check)
}

host_check <- secondary_raster_checks(host_file, infected)
if (host_check$checks_passed) {
  host <- host_check$raster
  if (raster::nlayers(host) > 1) {
    host <- output_from_raster_mean_and_sd(host)
  }
} else {
  return(host_check$failed_check)
}

total_plants_check <- secondary_raster_checks(total_plants_file, infected)
if (total_plants_check$checks_passed) {
  total_plants <- total_plants_check$raster
  if (raster::nlayers(total_plants) > 1) {
    total_plants <- output_from_raster_mean_and_sd(total_plants)
  }
} else {
  return(total_plants_check$failed_check)
}

infected_speci <- infected
susceptible_speci <- stack()
for (r in 1:length(infected_files)) {
  infected_name <-  paste("infected", r, sep = "")
  susceptible <- host - infected[[r]]
  susceptible[susceptible < 0] <- 0
  
  susceptible_speci <- stack(susceptible_speci, susceptible)
}

if (use_movements) {
  movements_check <- movement_checks(movements_file, infected, start_date, end_date)
  if (movements_check$checks_passed) {
    movements <- movements_check$movements
    movements_dates <- movements_check$movements_dates
    movements_r <- movements_check$movements_r
  } else {
    return(movements_check$failed_check)
  }
} else {
  movements <- list(0,0,0,0,0)
  movements_dates <- start_date
}

infected_species <- list(as.matrix(infected_speci[[1]]))
susceptible_species <- list(as.matrix(susceptible_speci[[1]]))
for(u in 2:nlayers(infected_speci)) {
  infected_species[[u]] <- as.matrix(infected_speci[[u]])
  susceptible_species[[u]] <- as.matrix(susceptible_speci[[u]])
}

if (use_lethal_temperature == TRUE) {
  temperature_check <- secondary_raster_checks(temperature_file, infected)
  if (temperature_check$checks_passed) {
    temperature_stack <- temperature_check$raster
  } else {
    return(temperature_check$failed_check)
  }
  
  temperature <- list(raster::as.matrix(temperature_stack[[1]]))
  for(o in 2:number_of_years) {
    temperature[[o]] <- raster::as.matrix(temperature_stack[[o]])
  }
} else {
  temperature <- host
  raster::values(temperature) <- 1
  temperature <- list(raster::as.matrix(temperature))
}

weather <- FALSE
if (temp == TRUE) {
  temperature_coefficient_check <- secondary_raster_checks(temperature_coefficient_file, infected)
  if (temperature_coefficient_check$checks_passed) {
    temperature_coefficient <- temperature_coefficient_check$raster
  } else {
    return(temperature_coefficient_check$failed_check)
  }
  
  weather <- TRUE
  weather_coefficient_stack <- temperature_coefficient
  if (precip ==TRUE){
    precipitation_coefficient_check <- secondary_raster_checks(precipitation_coefficient_file, infected)
    if (precipitation_coefficient_check$checks_passed) {
      precipitation_coefficient <- precipitation_coefficient_check$raster
    } else {
      return(precipitation_coefficient_check$failed_check)
    }
    
    weather_coefficient_stack <- weather_coefficient_stack * precipitation_coefficient
  }
} else if(precip == TRUE){
  precipitation_coefficient_check <- secondary_raster_checks(precipitation_coefficient_file, infected)
  if (precipitation_coefficient_check$checks_passed) {
    precipitation_coefficient <- precipitation_coefficient_check$raster
  } else {
    return(precipitation_coefficient_check$failed_check)
  }
  
  weather <- TRUE
  weather_coefficient_stack <- precipitation_coefficient
}

if (weather == TRUE){
  # weather_coefficient_stack <- raster::reclassify(weather_coefficient_stack, matrix(c(NA,0), ncol = 2, byrow = TRUE), right = NA)
  weather_coefficient <- list(raster::as.matrix(weather_coefficient_stack[[1]]))
  for(h in 2:number_of_time_steps) {
    weather_coefficient[[h]] <- raster::as.matrix(weather_coefficient_stack[[h]])
  }
} else {
  weather_coefficient <- host
  raster::values(weather_coefficient) <- 1
  weather_coefficient <- list(raster::as.matrix(weather_coefficient))
}

if (management == TRUE) {
  treatments_check <- secondary_raster_checks(treatments_file, infected)
  if (treatments_check$checks_passed) {
    treatment_stack <- treatments_check$raster
  } else {
    return(treatments_check$failed_check)
  }
  
  treatment_check <- treatment_checks(treatment_stack, treatments_file, pesticide_duration, treatment_dates, pesticide_efficacy)
  if (treatment_check$checks_passed) {
    treatment_maps <- treatment_check$treatment_maps
  } else {
    return(treatment_check$failed_check)
  }
} else {
  treatment_map <- host
  raster::values(treatment_map) <- 0
  treatment_maps <- list(raster::as.matrix(treatment_map))
}

ew_res <- raster::xres(susceptible)
ns_res <- raster::yres(susceptible)
num_cols <- raster::ncol(susceptible)
num_rows <- raster::nrow(susceptible)

mortality_tracker <- infected[[1]]
raster::values(mortality_tracker) <- 0

total_plants <- raster::as.matrix(total_plants)
mortality_tracker <- raster::as.matrix(mortality_tracker)
mortality <- mortality_tracker
resistant <- mortality_tracker

years <- seq(year(start_date), year(end_date), 1)



save.image(file = "/automanagement.RData")

