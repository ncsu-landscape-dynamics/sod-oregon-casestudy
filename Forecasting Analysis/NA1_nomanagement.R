
library(PoPS)
library(raster)
#invisible(utils::memory.limit(73728)) 

source('/rpops/R/uncertainty_propogation.R')
source('/rpops/R/helpers.R')
source('/rpops/R/checks.R')

## Read in data
infected_file = "/input data/end_inf_2019.tif"
host_file = "/input data/lide_100m_median_2019.tif"
total_plants_file = "/input data/lemma_max100m.tif"
temperature_coefficient_file = "/input data/average_weather_2020_2024.tif"
treatments_file = ""
#treatment_years = c(2016,2017)
treatment_month = 12
treatment_method = "all infected"
treatment_dates <- "2024-12-31"

use_lethal_temperature = FALSE
management = FALSE
temp = TRUE
precip = FALSE
season_month_start = 1
season_month_end = 12
time_step = "week"
start_date = "2020-01-01"
end_date = "2024-12-31"
lethal_temperature = NA
lethal_temperature_month = NA
random_seed = 42
reproductive_rate = .6
mortality_on = FALSE
mortality_rate = 0.05
mortality_time_lag = 2
percent_natural_dispersal <- 1.0
natural_kernel_type <- "exponential"
anthropogenic_kernel_type <- "cauchy"
natural_distance_scale <- 354
anthropogenic_distance_scale <- 0.0
natural_dir <- "N"
natural_kappa <- 2
anthropogenic_dir <- "NONE"
anthropogenic_kappa <- 0
output_frequency = "year"
pesticide_duration = c(0)
pesticide_efficiacy <- 0.0
mask <- NULL
output_frequency <- "year"
movements = ""
use_movements = FALSE
movements_dates = ""



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

if (is.null(random_seed)) {
  random_seed = round(stats::runif(1, 1, 1000000))
}

infected_check <- initial_raster_checks(infected_file)
if (infected_check$checks_passed) {
  infected <- infected_check$raster
  if (raster::nlayers(infected) > 1) {
    infected <- output_from_raster_mean_and_sd(infected)
  }
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

susceptible <- host - infected
susceptible[susceptible < 0] <- 0

if (use_lethal_temperature == TRUE) {
  temperature_check <- secondary_raster_checks(temperature_file, infected)
  if (temperature_check$checks_passed) {
    temperature_stack <- temperature_check$raster
  } else {
    return(temperature_check$failed_check)
  }
  
  temperature <- list(raster::as.matrix(temperature_stack[[1]]))
  for(i in 2:number_of_years) {
    temperature[[i]] <- raster::as.matrix(temperature_stack[[i]])
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
  for(i in 2:number_of_time_steps) {
    weather_coefficient[[i]] <- raster::as.matrix(weather_coefficient_stack[[i]])
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

mortality_tracker <- infected
raster::values(mortality_tracker) <- 0

infected <- raster::as.matrix(infected)
susceptible <- raster::as.matrix(susceptible)
total_plants <- raster::as.matrix(total_plants)
mortality_tracker <- raster::as.matrix(mortality_tracker)
mortality <- mortality_tracker
resistant <- mortality_tracker

reproductive_rate_check <- uncertainty_check(reproductive_rate, round_to = 1, n = 1)
if (reproductive_rate_check$checks_passed) {
  reproductive_rate <- reproductive_rate_check$value
} else {
  return(reproductive_rate_check$failed_check)
}

natural_distance_scale_check <- uncertainty_check(natural_distance_scale, round_to = 0, n = 1)
if (natural_distance_scale_check$checks_passed) {
  natural_distance_scale <- natural_distance_scale_check$value
} else {
  return(natural_distance_scale_check$failed_check)
}

anthropogenic_distance_scale_check <- uncertainty_check(anthropogenic_distance_scale, round_to = 0, n = 1)
if (anthropogenic_distance_scale_check$checks_passed) {
  anthropogenic_distance_scale <- anthropogenic_distance_scale_check$value
} else {
  return(anthropogenic_distance_scale_check$failed_check)
}

percent_natural_dispersal_check <- uncertainty_check(percent_natural_dispersal, round_to = 3, n = 1)
if (percent_natural_dispersal_check$checks_passed) {
  percent_natural_dispersal <- percent_natural_dispersal_check$value
} else {
  return(percent_natural_dispersal_check$failed_check)
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


data_na <- pops_model(random_seed = random_seed, 
                   use_lethal_temperature = use_lethal_temperature, 
                   lethal_temperature = lethal_temperature, lethal_temperature_month = lethal_temperature_month,
                   infected = infected,
                   susceptible = susceptible,
                   total_plants = total_plants,
                   mortality_on = mortality_on,
                   mortality_tracker = mortality_tracker,
                   mortality = mortality,
                   treatment_maps = treatment_maps,
                   treatment_dates = treatment_dates,
                   pesticide_duration = pesticide_duration,
                   resistant = resistant,
                   use_movements = use_movements,
                   movements = movements,
                   movements_dates = movements_dates,
                   weather = weather,
                   temperature = temperature,
                   weather_coefficient = weather_coefficient,
                   ew_res = ew_res, ns_res = ns_res, num_rows = num_rows, num_cols = num_cols,
                   time_step = time_step, reproductive_rate = reproductive_rate,
                   mortality_rate = mortality_rate, mortality_time_lag = mortality_time_lag,
                   season_month_start = season_month_start, season_month_end = season_month_end,
                   start_date = start_date, end_date = end_date,
                   treatment_method = treatment_method,
                   natural_kernel_type = natural_kernel_type, anthropogenic_kernel_type = anthropogenic_kernel_type, 
                   use_anthropogenic_kernel = use_anthropogenic_kernel, percent_natural_dispersal = percent_natural_dispersal,
                   natural_distance_scale = natural_distance_scale, anthropogenic_distance_scale = anthropogenic_distance_scale, 
                   natural_dir = natural_dir, natural_kappa = natural_kappa,
                   anthropogenic_dir = anthropogenic_dir, anthropogenic_kappa = anthropogenic_kappa, 
                   output_frequency = output_frequency
)


data2_na <- list(data_na)
for (seed in 1:500) {
  random_seed <- round(stats::runif(1, 1, 1000000))
  data2_na[[seed]] <- pops_model(random_seed = random_seed, 
                              use_lethal_temperature = use_lethal_temperature, 
                              lethal_temperature = lethal_temperature, lethal_temperature_month = lethal_temperature_month,
                              infected = infected,
                              susceptible = susceptible,
                              total_plants = total_plants,
                              mortality_on = mortality_on,
                              mortality_tracker = mortality_tracker,
                              mortality = mortality,
                              treatment_maps = treatment_maps,
                              treatment_dates = treatment_dates,
                              pesticide_duration = pesticide_duration,
                              resistant = resistant,
                              use_movements = use_movements,
                              movements = movements,
                              movements_dates = movements_dates,
                              weather = weather,
                              temperature = temperature,
                              weather_coefficient = weather_coefficient,
                              ew_res = ew_res, ns_res = ns_res, num_rows = num_rows, num_cols = num_cols,
                              time_step = time_step, reproductive_rate = reproductive_rate,
                              mortality_rate = mortality_rate, mortality_time_lag = mortality_time_lag,
                              season_month_start = season_month_start, season_month_end = season_month_end,
                              start_date = start_date, end_date = end_date,
                              treatment_method = treatment_method,
                              natural_kernel_type = natural_kernel_type, anthropogenic_kernel_type = anthropogenic_kernel_type, 
                              use_anthropogenic_kernel = use_anthropogenic_kernel, percent_natural_dispersal = percent_natural_dispersal,
                              natural_distance_scale = natural_distance_scale, anthropogenic_distance_scale = anthropogenic_distance_scale, 
                              natural_dir = natural_dir, natural_kappa = natural_kappa,
                              anthropogenic_dir = anthropogenic_dir, anthropogenic_kappa = anthropogenic_kappa, 
                              output_frequency = output_frequency)
}


# num_iterations = 100
# number_cores = NA
# infected_years_file = "H:/Shared Drives/APHIS  Projects/PoPS/Case Studies/sudden_oak_death/Oregon/Infection Data - Proportion/300mLEMMA/Cumulative Infections/cum_inf_2018eu.tif"
# 
# test_validate <- validate(infected_years_file, num_iterations, number_cores,
#                           infected_file, host_file, total_plants_file, reproductive_rate = 1.0,
#                           use_lethal_temperature = FALSE, temp = TRUE, precip = FALSE, management = TRUE, 
#                           mortality_on = TRUE, temperature_file = "", temperature_coefficient_file, 
#                           precipitation_coefficient_file ="", treatments_file,
#                           season_month_start = 1, season_month_end = 12, time_step = time_step,
#                           start_time = 2018, end_time = 2018, treatment_years = c(2018),
#                           dispersal_kern = "exponential", percent_short_distance_dispersal = 1.0,
#                           natural_distance_scale = natural_distance_scale, long_distance_scale = 0.0,
#                           lethal_temperature = NA, lethal_temperature_month = 1,
#                           mortality_rate = 0.05, mortality_time_lag = 2,
#                           natural_dir = wind_dir, natural_kappa = natural_kappa)



infected = raster(infected_file)
infected[is.na(infected)] <- 0
infected_stack <- stack()
infected_list <- list()
for (j in 1:length(data2_na)) {
  infected_stack <- stack()
  for (q in 1:length(data2_na[[1]]$infected)) {
    infected[] <- data2_na[[j]]$infected[[q]]
    infected_stack <- stack(infected_stack, infected)
  }
  infected_list[[j]] <- infected_stack
}

m <- c(0, .99, 0,
       1, Inf,1)
rcl_mat <- matrix(m, ncol=3, byrow=TRUE)
prediction_list <- list()
probability <- prediction <- reclassify(infected_list[[1]], rcl_mat)
probability[probability > 0] <- 0
for (i in 1:length(infected_list)) {
  prediction <- reclassify(infected_list[[i]], rcl_mat)
  prediction_list[[i]] <- prediction
  probability <- probability + prediction
}

all_2020_inf = stack()
all_2021_inf = stack()
all_2022_inf = stack()
all_2023_inf = stack()
all_2024_inf = stack()

for (k in 1:length(infected_list)) {
  run = infected_list[[k]]
  run_stack = stack(run)
  inf_2020 = subset(run_stack, 1:1)
  inf_2021 = subset(run_stack, 2:2)
  inf_2022 = subset(run_stack, 3:3)
  inf_2023 = subset(run_stack, 4:4)
  inf_2024 = subset(run_stack, 5:5)
  all_2020_inf = stack(all_2020_inf, inf_2020)
  all_2021_inf = stack(all_2021_inf, inf_2021)
  all_2022_inf = stack(all_2022_inf, inf_2022)
  all_2023_inf = stack(all_2023_inf, inf_2023)
  all_2024_inf = stack(all_2024_inf, inf_2024)
}

probability_before_treat_2020 = subset(probability, 1:1)
probability_before_treat_2021 = subset(probability, 2:2)
probability_before_treat_2022 = subset(probability, 3:3)
probability_before_treat_2023 = subset(probability, 4:4)
probability_before_treat_2024 = subset(probability, 5:5)

# to write to file
CRS = crs(host)
crs(probability_before_treat_2020) <- CRS
crs(probability_before_treat_2021) <- CRS
crs(probability_before_treat_2022) <- CRS
crs(probability_before_treat_2023) <- CRS
crs(probability_before_treat_2024) <- CRS

north_rates_2020 = data.frame()
east_rates_2020 = data.frame()
west_rates_2020 = data.frame()
south_rates_2020 = data.frame()
for (k in 1:length(infected_list)) {
  ratelist = data2_na[[k]]$rates[[1]]
  north_rate = data.frame(ratelist[[1]])
  south_rate = data.frame(ratelist[[2]])
  east_rate = data.frame(ratelist[[3]])
  west_rate = data.frame(ratelist[[4]])
  north_rates_2020 <- rbind(north_rates_2020, north_rate)
  south_rates_2020 <- rbind(south_rates_2020, south_rate)
  west_rates_2020 <- rbind(west_rates_2020, west_rate)
  east_rates_2020 <- rbind(east_rates_2020, east_rate)
}

north_rates_2021 = data.frame()
east_rates_2021 = data.frame()
west_rates_2021 = data.frame()
south_rates_2021 = data.frame()
for (k in 1:length(infected_list)) {
  ratelist = data2_na[[k]]$rates[[2]]
  north_rate = data.frame(ratelist[[1]])
  south_rate = data.frame(ratelist[[2]])
  east_rate = data.frame(ratelist[[3]])
  west_rate = data.frame(ratelist[[4]])
  north_rates_2021 <- rbind(north_rates_2021, north_rate)
  south_rates_2021 <- rbind(south_rates_2021, south_rate)
  west_rates_2021 <- rbind(west_rates_2021, west_rate)
  east_rates_2021 <- rbind(east_rates_2021, east_rate)
}

north_rates_2022 = data.frame()
east_rates_2022 = data.frame()
west_rates_2022 = data.frame()
south_rates_2022 = data.frame()
for (k in 1:length(infected_list)) {
  ratelist = data2_na[[k]]$rates[[3]]
  north_rate = data.frame(ratelist[[1]])
  south_rate = data.frame(ratelist[[2]])
  east_rate = data.frame(ratelist[[3]])
  west_rate = data.frame(ratelist[[4]])
  north_rates_2022 <- rbind(north_rates_2022, north_rate)
  south_rates_2022 <- rbind(south_rates_2022, south_rate)
  west_rates_2022 <- rbind(west_rates_2022, west_rate)
  east_rates_2022 <- rbind(east_rates_2022, east_rate)
}

north_rates_2023 = data.frame()
east_rates_2023 = data.frame()
west_rates_2023 = data.frame()
south_rates_2023 = data.frame()
for (k in 1:length(infected_list)) {
  ratelist = data2_na[[k]]$rates[[4]]
  north_rate = data.frame(ratelist[[1]])
  south_rate = data.frame(ratelist[[2]])
  east_rate = data.frame(ratelist[[3]])
  west_rate = data.frame(ratelist[[4]])
  north_rates_2023 <- rbind(north_rates_2023, north_rate)
  south_rates_2023 <- rbind(south_rates_2023, south_rate)
  west_rates_2023 <- rbind(west_rates_2023, west_rate)
  east_rates_2023 <- rbind(east_rates_2023, east_rate)
}

north_rates_2024 = data.frame()
east_rates_2024 = data.frame()
west_rates_2024 = data.frame()
south_rates_2024 = data.frame()
for (k in 1:length(infected_list)) {
  ratelist = data2_na[[k]]$rates[[5]]
  north_rate = data.frame(ratelist[[1]])
  south_rate = data.frame(ratelist[[2]])
  east_rate = data.frame(ratelist[[3]])
  west_rate = data.frame(ratelist[[4]])
  north_rates_2024 <- rbind(north_rates_2024, north_rate)
  south_rates_2024 <- rbind(south_rates_2024, south_rate)
  west_rates_2024 <- rbind(west_rates_2024, west_rate)
  east_rates_2024 <- rbind(east_rates_2024, east_rate)
}



writeRaster(probability_before_treat_2020, "/No Management Results/na_probability_2020.tif", format='GTiff')
writeRaster(probability_before_treat_2021, "/No Management Results/na_probability_2021.tif", format='GTiff')
writeRaster(probability_before_treat_2022, "/No Management Results/na_probability_2022.tif", format='GTiff')
writeRaster(probability_before_treat_2023, "/No Management Results/na_probability_2023.tif", format='GTiff')
writeRaster(probability_before_treat_2024, "/No Management Results/na_probability_2024.tif", format='GTiff')

writeRaster(all_2020_inf, "/No Management Results/na_inf_2020.tif", format='GTiff')
writeRaster(all_2021_inf, "/No Management Results/na_inf_2021.tif", format='GTiff')
writeRaster(all_2022_inf, "/No Management Results/na_inf_2022.tif", format='GTiff')
writeRaster(all_2023_inf, "/No Management Results/na_inf_2023.tif", format='GTiff')
writeRaster(all_2024_inf, "/No Management Results/na_inf_2024.tif", format='GTiff')


