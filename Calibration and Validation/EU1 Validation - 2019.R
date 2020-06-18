
## To validate spore rate and distance parameters for EU1 Sudden Oak Death disease strain in Oregon

## Uses 100m tanoak proportion data processed from LEMMA OSU Group
## Works with PoPS Package
## Uses 2018 infection data, runs for 2019, compares to 2019 actual infections

## If PoPS isn't already installed
# library(devtools)
# devtools::install_github("ncsu-landscape-dynamics/rpops", ref="feature/disagreements")

library(PoPS)
library(raster)
invisible(utils::memory.limit(73728)) 

## Read in data
infected_file = "/validation input data/inf_2018_eu.tif"
host_file = "/validation input data/tanoak_2018.tif"
total_plants_file = "/validation input data/max_density.tif"
weather_coefficient_file = "/validation input data/weather_coef_2017_2018.tif"
treatments_file = ""
treatment_years = c(2018)
treatment_month = 12
treatment_method = "all infected"

actual_inf_2019 = raster("/validation input data/inf_2019_eu.tif")

use_lethal_temperature = FALSE
management = FALSE
weather = TRUE
season_month_start = 1
season_month_end = 12
time_step = "week"
start_time = 2019
end_time = 2019
lethal_temperature = NA
lethal_temperature_month = NA
random_seed = 42
reproductive_rate = 1.6
mortality_on = TRUE
mortality_rate = 0.05
mortality_time_lag = 2
percent_natural_dispersal <- 1.0
natural_kernel_type <- "exponential"
anthropogenic_kernel_type <- "cauchy"
natural_distance_scale <- 242
anthropogenic_distance_scale <- 0.0
natural_dir <- "N"
natural_kappa <- 3
anthropogenic_dir <- "NONE"
anthropogenic_kappa <- 0


number_of_years = end_time-start_time+1

if (time_step == "week") {
  number_of_time_steps = (end_time-start_time+1)*52
} else if (time_step == "month") {
  number_of_time_steps = (end_time-start_time+1)*12
} else if (time_step == "day") {
  number_of_time_steps = (end_time-start_time+1)*365
}

infected = raster(infected_file)
infected[is.na(infected)] <- 0
dataType(infected) <- "INT2U"
host = raster(host_file)
host[is.na(host)] <- 0
dataType(host) <- "INT2U"
susceptible = host - infected
susceptible[is.na(susceptible)] <- 0
dataType(susceptible) <- "INT2U"
ew_res = xres(susceptible)
ns_res = yres(susceptible)
num_cols <- raster::ncol(susceptible)
num_rows <- raster::nrow(susceptible)
total_plants = raster(total_plants_file)
total_plants[is.na(total_plants)] <- 0
dataType(total_plants) <- "INT2U"
weather_coefficient_allyears = stack(weather_coefficient_file)
weather_coefficient_stack = subset(weather_coefficient_allyears, 53:105)


if (!(raster::extent(infected) == raster::extent(host) && raster::extent(infected) == raster::extent(total_plants))) {
  print("Extents of input rasters do not match. Ensure that all of your input rasters have the same extent") 
} else {
  print("Extents okay") 
}

if (!(raster::xres(infected) == raster::xres(host) && raster::xres(infected) == raster::xres(total_plants) && raster::yres(infected) == raster::yres(host) && raster::yres(infected) == raster::yres(total_plants))) {
  print("Resolution of input rasters do not match. Ensure that all of your input rasters have the same resolution")
} else {
  print("Resolutions okay")
}

if (!(raster::compareCRS(host,infected) && raster::compareCRS(host, total_plants))) {
  print("Coordinate reference system (crs) of input rasters do not match. Ensure that all of your input rasters have the same crs")
} else {
  print("coordinates okay")
}


if (use_lethal_temperature == TRUE) {
  temperature_stack = stack(temperature_file)
  temperature_stack[is.na(temperature_stack)] <- 0
  temperature = list(as.matrix(temperature_stack[[1]]))
  for(i in 2:number_of_years) {
    temperature[[i]] <- as.matrix(temperature_stack[[i]])
  }
} else {
  temperature <- host
  values(temperature) <- 1
  temperature <- list(as.matrix(temperature))
}

# weather = FALSE
# if (temp == TRUE) {
#   temperature_coefficient_all = stack(temperature_coefficient_file)
#   temperature_coefficient = subset(temperature_coefficient_all, 1873:1976)
#   weather = TRUE
#   weather_coefficient_stack = temperature_coefficient
#   if (precip ==TRUE){
#     precipitation_coefficient_all = stack(precipitation_coefficient_file)
#     precipitation_coefficient = subset(precipitation_coefficient_all, 1873:1976)
#     weather_coefficient_stack = weather_coefficient_stack * precipitation_coefficient
#   }
# } else if(precip == TRUE){
#   precipitation_coefficient = stack(precipitation_coefficient_file)
#   weather = TRUE
#   weather_coefficient_stack = precipitation_coefficient
# }

if (weather == TRUE){
  weather_coefficient_stack[is.na(weather_coefficient_stack)] <- 0
  weather_coefficient <- list(as.matrix(weather_coefficient_stack[[1]]))
  for(i in 2:number_of_time_steps) {
    weather_coefficient[[i]] <- as.matrix(weather_coefficient_stack[[i]])
  }
} else {
  weather_coefficient <- host
  values(weather_coefficient) <- 1
  weather_coefficient <- list(as.matrix(weather_coefficient))
}

if (management == TRUE) {
  treatment_stack <- stack(treatments_file)
  treatment_stack[is.na(treatment_stack)] <- 0
  treatment_maps <- list(as.matrix(treatment_stack[[1]]))
  if (nlayers(treatment_stack) >= 2) {
    for(i in 2:nlayers(treatment_stack)) {
      treatment_maps[[i]] <- as.matrix(treatment_stack[[i]])
    }
  }
  treatment_years = treatment_years
} else {
  treatment_map <- host
  values(treatment_map) <- 0
  treatment_maps = list(as.matrix(treatment_map))
  treatment_years = c(0)
}



mortality_tracker = infected
values(mortality_tracker) <- 0
mortality <- mortality_tracker


if(percent_natural_dispersal == 1.0) {
  use_anthropogenic_kernel = FALSE
} else if (percent_natural_dispersal < 1.0  && percent_natural_dispersal >= 0.0) {
  use_anthropogenic_kernel = TRUE
} else {
  return("Percent natural dispersal must be between 0.0 and 1.0")
}


infected = as.matrix(infected)
susceptible = as.matrix(susceptible)
total_plants = as.matrix(total_plants)
mortality_tracker = as.matrix(mortality_tracker)
mortality <- as.matrix(mortality)

data <- pops_model(random_seed = random_seed, 
                   use_lethal_temperature = use_lethal_temperature, 
                   lethal_temperature = lethal_temperature, lethal_temperature_month = lethal_temperature_month,
                   infected = infected,
                   susceptible = susceptible,
                   total_plants = total_plants,
                   mortality_on = mortality_on,
                   mortality_tracker = mortality_tracker,
                   mortality = mortality,
                   treatment_maps = treatment_maps,
                   treatment_years = treatment_years,
                   weather = weather,
                   temperature = temperature,
                   weather_coefficient = weather_coefficient,
                   ew_res = ew_res, ns_res = ns_res, num_rows = num_rows, num_cols = num_cols,
                   time_step = time_step, reproductive_rate = reproductive_rate,
                   mortality_rate = mortality_rate, mortality_time_lag = mortality_time_lag,
                   season_month_start = season_month_start, season_month_end = season_month_end,
                   start_time = start_time, end_time = end_time,
                   treatment_month = treatment_month, treatment_method = treatment_method,
                   natural_kernel_type = natural_kernel_type, anthropogenic_kernel_type = anthropogenic_kernel_type, 
                   use_anthropogenic_kernel = use_anthropogenic_kernel, percent_natural_dispersal = percent_natural_dispersal,
                   natural_distance_scale = natural_distance_scale, anthropogenic_distance_scale = anthropogenic_distance_scale, 
                   natural_dir = natural_dir, natural_kappa = natural_kappa,
                   anthropogenic_dir = anthropogenic_dir, anthropogenic_kappa = anthropogenic_kappa
)



data2 <- list(data)
for (seed in 1:1000) {
  random_seed <- round(stats::runif(1, 1, 1000000))
  data2[[seed]] <- pops_model(random_seed = random_seed, 
                              use_lethal_temperature = use_lethal_temperature, 
                              lethal_temperature = lethal_temperature, lethal_temperature_month = lethal_temperature_month,
                              infected = infected,
                              susceptible = susceptible,
                              total_plants = total_plants,
                              mortality_on = mortality_on,
                              mortality_tracker = mortality_tracker,
                              mortality = mortality,
                              treatment_maps = treatment_maps,
                              treatment_years = treatment_years,
                              weather = weather,
                              temperature = temperature,
                              weather_coefficient = weather_coefficient,
                              ew_res = ew_res, ns_res = ns_res, num_rows = num_rows, num_cols = num_cols,
                              time_step = time_step, reproductive_rate = reproductive_rate,
                              mortality_rate = mortality_rate, mortality_time_lag = mortality_time_lag,
                              season_month_start = season_month_start, season_month_end = season_month_end,
                              start_time = start_time, end_time = end_time,
                              treatment_month = treatment_month, treatment_method = treatment_method,
                              natural_kernel_type = natural_kernel_type, anthropogenic_kernel_type = anthropogenic_kernel_type, 
                              use_anthropogenic_kernel = use_anthropogenic_kernel, percent_natural_dispersal = percent_natural_dispersal,
                              natural_distance_scale = natural_distance_scale, anthropogenic_distance_scale = anthropogenic_distance_scale, 
                              natural_dir = natural_dir, natural_kappa = natural_kappa,
                              anthropogenic_dir = anthropogenic_dir, anthropogenic_kappa = anthropogenic_kappa)
}



infected = raster(infected_file)
infected[is.na(infected)] <- 0
infected_stack <- stack()
infected_list <- list()
for (j in 1:length(data2)) {
  infected_stack <- stack()
  for (q in 1:length(data2[[1]]$infected_before_treatment)) {
    infected[] <- data2[[j]]$infected_before_treatment[[q]]
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

all_2019_inf = stack()

for (k in 1:length(infected_list)) {
  run = infected_list[[k]]
  run_stack = stack(run)
  inf_2019 = subset(run_stack, 1:1)
  all_2019_inf = stack(all_2019_inf, inf_2019)
}

probability_before_treat_2019 = subset(probability, 1:1)

# to write to file
CRS = crs(host)
crs(probability_before_treat_2019) <- CRS


writeRaster(probability_before_treat_2019, "/EU1 Validation - 2019/probability_2019.tif", format='GTiff')
writeRaster(all_2019_inf, "/EU1 Validation - 2019/inf_2019.tif", format='GTiff')



#### VALIDATION COMPONENT####


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

quantity_allocation_disagreement <- function(reference, comparison){
  # test that the comparison raster is the same extent, resolution, and crs as the reference (if not end)
  raster::compareRaster(reference, comparison)
  compare <- reference - comparison
  compare3 <- reference + comparison
  # extent = extent(reference)
  # num_of_cells = max(cellsFromExtent(reference, extent(reference)))
  
  # calculate number of infected patches
  NP_ref <- landscapemetrics::lsm_c_np(reference, directions = 8)$value[2]
  if (sum(comparison[comparison >0]) == 0) {
    NP_comp <- 0
    ENN_MN_comp <- 0
    PARA_MN_comp <- 0
    LPI_comp <- 0
  } else {
    NP_comp <- landscapemetrics::lsm_c_np(comparison, directions = 8)$value[2]
  }
  
  change_NP <- abs((NP_comp - NP_ref)/NP_ref)
  if (change_NP > 1) {change_NP <- 1}
  if (change_NP >= 1) {change_NP <- 1}
  
  # calculate the mean euclidean distance between patches
  if (NP_ref > 1) {
    ENN_MN_ref <- landscapemetrics::lsm_c_enn_mn(reference, directions = 8, verbose = TRUE)$value[2]
  } else  if (NP_ref == 1) {
    ENN_MN_ref <- 0
  }
  
  if (sum(comparison[comparison > 0]) != 0 && NP_comp > 1) {
    ENN_MN_comp <- landscapemetrics::lsm_c_enn_mn(comparison, directions = 8, verbose = TRUE)$value[2]
  } else  if (sum(comparison[comparison > 0]) != 0 && NP_comp <= 1) {
    ENN_MN_comp <- 0
  }
  
  if (ENN_MN_ref != 0) {
    change_ENN_MN <- abs((ENN_MN_comp - ENN_MN_ref)/ENN_MN_ref)
    if (change_ENN_MN > 1) {change_ENN_MN <- 1}
  } else if (ENN_MN_comp == 0 && ENN_MN_ref == 0) {
    change_ENN_MN <- 0
  } else {
    change_ENN_MN <- 1
  }
  
  # calculate the mean perimeter-area ratio of patches and the difference
  PARA_MN_ref <- landscapemetrics::lsm_c_para_mn(reference, directions = 8)$value[2]
  if (sum(comparison[comparison >0]) == 0) {
    PARA_MN_comp <- 0
  } else if (sum(comparison[comparison >0]) != 0) {
    PARA_MN_comp <- landscapemetrics::lsm_c_para_mn(comparison, directions = 8)$value[2]
  }
  
  change_PARA_MN <- abs((PARA_MN_comp - PARA_MN_ref)/PARA_MN_ref)
  if (change_PARA_MN > 1) {change_PARA_MN <- 1}
  
  # calculate the largest patch index and difference
  LPI_ref <- landscapemetrics::lsm_c_lpi(reference, directions = 8)$value[2]
  if (sum(comparison[comparison >0]) == 0) {
    LPI_comp <- 0
  } else if (sum(comparison[comparison >0]) != 0) {
    LPI_comp <- landscapemetrics::lsm_c_lpi(comparison, directions = 8)$value[2]
  }
  
  change_LPI <- abs((LPI_comp - LPI_ref)/(LPI_ref))
  if (change_LPI > 1) {change_LPI <- 1}
  
  # calculate configuration disagreement between reference and comparison
  configuration_disagreement <- ((change_NP + change_ENN_MN + change_PARA_MN + change_LPI) / 4)
  
  
  ## create data frame for comparison
  output <- data.frame(quantity_disagreement = 0, allocation_disagreement = 0, total_disagreement = 0, omission = 0, commission = 0 , number_of_infected_comp = 0, directional_disagreement = 0, landscape_similarity = 0)
  output$total_disagreement <- sum(compare[compare == 1]) + abs(sum(compare[compare == -1]))
  output$quantity_disagreement <- abs(sum(compare[compare == 1]) + sum(compare[compare == -1]))
  output$allocation_disagreement <- output$total_disagreement - output$quantity_disagreement
  output$allocation_disagreement_num_pixels <- output$allocation_disagreement/2
  output$omission <- abs(sum(compare[compare == 1]))
  output$commission <- abs(sum(compare[compare == -1]))
  output$number_of_infected_comp <- sum(comparison[comparison == 1])
  output$directional_disagreement <- sum(compare[compare == 1]) + sum(compare[compare == -1])
  output$configuration_disagreement <- configuration_disagreement
  output$true_positives <- abs(sum(compare3[compare3 ==2]))/2
  # output$true_negatives <- num_of_cells - output$omission - output$commission - output$true_positives
  # output$odds_ratio = (output$true_positives*output$true_negatives)/(output$omission*output$commission)
  
  
  return(output)
}

to_write_csv_2019 = "/EU1 Validation - 2019/disagreements_2019.csv"



## 2017:

disagreement_df_2019 <- data.frame()
rcl <- c(1, 3000, 1, 0, 0.99, NA)
rclmat <- matrix(rcl, ncol=3, byrow=TRUE)
for (i in 1:nlayers(all_2019_inf)){
  reference = reclassify(actual_inf_2019, rclmat)
  comparison = subset(all_2019_inf, i:i)
  comparison = reclassify(comparison, rclmat)
  output = data.frame(quantity_allocation_disagreement(reference, comparison))
  disagreement_df_2019 <- rbind(disagreement_df_2019, output)
}

write.csv(disagreement_df_2019, to_write_csv_2019)

average_2019_quantity_disagreement = mean(disagreement_df_2019$quantity_disagreement)
print(paste("Average 2019 Quantity Disagreement: ", average_2019_quantity_disagreement))
min_2019_quantity_disagreement = min(disagreement_df_2019$quantity_disagreement)
print(paste("Min 2019 Quantity Disagreement: ", min_2019_quantity_disagreement))
max_2019_quantity_disagreement = max(disagreement_df_2019$quantity_disagreement)
print(paste("Max 2019 quantity disagreement: ", max_2019_quantity_disagreement))
mode_2019_quantity_disagreement = getmode(disagreement_df_2019$quantity_disagreement)
print(paste("Mode 2019 Quantity Disagreement: ", mode_2019_quantity_disagreement))

average_2019_allocation_disagreement = mean(disagreement_df_2019$allocation_disagreement)
print(paste("Average 2019 Allocation Disagreement: ", average_2019_allocation_disagreement))
min_2019_allocation_disagreement = min(disagreement_df_2019$allocation_disagreement)
print(paste("Min 2019 Allocation Disagreement: ", min_2019_allocation_disagreement))
max_2019_allocation_disagreement = max(disagreement_df_2019$allocation_disagreement)
print(paste("Max 2019 Allocation disagreement: ", max_2019_allocation_disagreement))
mode_2019_allocation_disagreement = getmode(disagreement_df_2019$allocation_disagreement)
print(paste("Mode 2019 Allocation Disagreement: ", mode_2019_allocation_disagreement))

average_2019_true_positives = mean(disagreement_df_2019$true_positives)
print(paste("Average 2019 true positives: ", average_2019_true_positives))
min_2019_true_positives = min(disagreement_df_2019$true_positives)
print(paste("Min 2019 true positives: ", min_2019_true_positives))
max_2019_true_positives = max(disagreement_df_2019$true_positives)
print(paste("Max 2019 true positives: ", max_2019_true_positives))
mode_2019_true_positives = getmode(disagreement_df_2019$true_positives)
print(paste("Mode 2019 true positives: ",mode_2019_true_positives))

average_2019_configuration_disagreement = mean(disagreement_df_2019$configuration_disagreement)
print(paste("Average 2019 landscape: ", average_2019_configuration_disagreement))
min_2019_configuration_disagreement = min(disagreement_df_2019$configuration_disagreement)
print(paste("Min 2019 landscape: ", min_2019_configuration_disagreement))
max_2019_configuration_disagreement = max(disagreement_df_2019$configuration_disagreement)
print(paste("Max 2019 landscape: ", max_2019_configuration_disagreement))
mode_2019_configuration_disagreement = getmode(disagreement_df_2019$configuration_disagreement)
print(paste("Mode 2019 landscape: ",mode_2019_configuration_disagreement))

