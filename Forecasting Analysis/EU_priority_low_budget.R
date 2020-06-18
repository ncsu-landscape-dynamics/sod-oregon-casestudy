
## Testing Scenarios of Sudden Oak Death Management in Oregon
## In this code, we focus on the EU1 strain of the disease first.
## If all EU1 are treated and resources remain, NA1 will be removed. 

## Budget: 2.5 million dollars per year
## 300ft buffer around treatments

library(PoPS)
library(raster)
library(doParallel)
library(iterators)
library(foreach)
library(parallel)
library(lubridate)

load('/automanagement.Rdata')
source("/rpops/R/helpers.R")


#budget:
num_iterations = 500
number_of_cores = 10
cost_per_meter_sq = 1.37
budget = 2500000
buffer <- 91
points <- data.frame(i = as.integer(1), j = ncol(infected[[1]]))
direction_first <- TRUE
treatment_priority = "ranked"
treatment_rank = c(1,0)
selection_method = "Points"
selection_priority = 'group size'
priority = "group size"
treatment_efficacy = 1
direction_first = TRUE




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

## management module information
num_cells <- round((budget/cost_per_meter_sq)/(ew_res*ns_res))
buffer_cells <- buffer/ew_res
years_simulated <- length(years)

random_seeds <- round(stats::runif(num_iterations, 1, 1000000))

treatment_speci <- raster()

## preallocate data frames and lists
s <- paste("sim_", years,sep = "")
s2 <- c("management_year", s)
infected_areas <- data.frame(matrix(ncol = length(years) + 1, nrow = length(years)))
names(infected_areas) <- s2
infected_areas$management_year <- years
number_infecteds <- infected_areas
west_rates <- infected_areas
east_rates <- infected_areas
south_rates <- infected_areas
north_rates <- infected_areas

## maintain starting state of the system for each simulation
start_date2 <- start_date
end_date2 <- end_date
years2 <- years
infected_species_2 <- infected_species
susceptible_species_2 <- susceptible_species
weather_coefficient_2 <- weather_coefficient
infected_speci_2 <- infected_speci
susceptible_speci_2 <- susceptible_speci

year_names <- paste(years)

points <- data.frame(i = as.integer(1), j = ncol(infected[[1]]))

if (is.na(number_of_cores) || number_of_cores > parallel::detectCores()) {
  core_count <- parallel::detectCores() - 1
} else {
  core_count <- number_of_cores
}
cl <- makeCluster(core_count)
registerDoParallel(cl)

sim_start_time <- Sys.time()

runs <- foreach(p = 1:num_iterations, .combine = rbind, .packages = c("raster", "PoPS", "foreach", "lubridate", "rlist")) %dopar% {
  infected_speci <- infected_speci_2
  infected_species <- infected_species_2
  susceptible_speci <- susceptible_speci_2
  susceptible_species <- susceptible_species_2
  
  run_years <-   foreach(y = 1:years_simulated, .combine = rbind, .packages = c("raster", "PoPS", "foreach", "lubridate", "rlist")) %do% {
    
    
    print("start")
    treatment <- treatmentAuto(rasts = infected_speci, rasts2 = susceptible_speci, 
                               method = selection_method, priority = selection_priority,
                               number_of_locations = num_cells, points = points, 
                               treatment_efficacy = treatment_efficacy, 
                               buffer_cells = buffer_cells, direction_first = direction_first, 
                               treatment_rank = treatment_rank, treatment_priority = treatment_priority)
    
    if (y == 1) {
      treatment_dates <- c(paste(years[y], "-12", "-01", sep = ""))
      treatment_maps <- list(as.matrix(treatment))
    } else if (y > 1) {
      treatment_dates <- c(treatment_dates, paste(years[y], "-12", "-01", sep = ""))
      treatment_maps <- list.append(treatment_maps, as.matrix(treatment))
      pesticide_duration <- c(pesticide_duration, 0)
    }
    
    management <- TRUE
    print("end_treatment")
    
    species_run <-   foreach(i = 1:length(infected_files), .combine = rbind, .packages = c("raster", "PoPS")) %do% {
      data <- pops_model(random_seed = random_seeds[p], 
                         use_lethal_temperature = use_lethal_temperature, 
                         lethal_temperature = lethal_temperature, 
                         lethal_temperature_month = lethal_temperature_month,
                         infected = infected_species[[i]],
                         susceptible = susceptible_species[[i]],
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
                         ew_res = ew_res, 
                         ns_res = ns_res, 
                         num_rows = num_rows, 
                         num_cols = num_cols,
                         time_step = time_step, 
                         reproductive_rate = reproductive_rate[[i]],
                         mortality_rate = mortality_rate, 
                         mortality_time_lag = mortality_time_lag,
                         season_month_start = season_month_start, 
                         season_month_end = season_month_end,
                         start_date = start_date,
                         end_date = end_date,
                         treatment_method = treatment_method,
                         natural_kernel_type = natural_kernel_type[[i]],
                         anthropogenic_kernel_type = anthropogenic_kernel_type[[i]], 
                         use_anthropogenic_kernel = use_anthropogenic_kernel,
                         percent_natural_dispersal = percent_natural_dispersal[[i]],
                         natural_distance_scale = natural_distance_scale[[i]],
                         anthropogenic_distance_scale = anthropogenic_distance_scale[[i]], 
                         natural_dir = natural_dir[[i]],
                         natural_kappa = natural_kappa[[i]],
                         anthropogenic_dir = anthropogenic_dir[[i]],
                         anthropogenic_kappa = anthropogenic_kappa[[i]],
                         output_frequency = output_frequency
      )
      
      infected_run <- raster::stack(lapply(1:length(data$infected), function(x) host))
      susceptible_run <- raster::stack(lapply(1:length(data$infected), function(x) host))
      
      for (q in 1:raster::nlayers(infected_run)) {
        infected_run[[q]] <- data$infected[[q]]
        susceptible_run[[q]] <- data$susceptible[[q]]
      }
      
      number_infected <- data$number_infected
      spread_rate <- data$rates
      infected_area <- data$area_infected
      
      west_rate <- c()
      east_rate <- c()
      south_rate <- c()
      north_rate <- c()
      
      for (k in 1:length(spread_rate)) {
        west_rate[k] <- spread_rate[[k]][4]
        east_rate[k] <- spread_rate[[k]][3]
        south_rate[k] <- spread_rate[[k]][2]
        north_rate[k] <- spread_rate[[k]][1]
      }
      
      to.species_run <- list(infected_run, susceptible_run, number_infected, infected_area, west_rate, east_rate, north_rate, south_rate)
    }
    
    infections_out <- c()
    susceptibles_out <- c()
    number_infecteds_out <- c()
    infected_areas_out <- c()
    west_rate_out <- c()
    east_rate_out <- c()
    north_rate_out <- c()
    south_rate_out <- c()
    infected_speci <- stack()
    susceptible_speci <- stack()
    
    for (t in 1:length(infected_species)) {
      infected_speci <- stack(infected_speci, species_run[[t]][[y]])
      susceptible_speci <- stack(susceptible_speci, species_run[[t+length(infected_species)]][[y]])
      
      infections_out[[t]] <- c(species_run[[t]])
      susceptibles_out[[t]] <- c(species_run[[t+length(infected_species)]])
      number_infecteds_out[[t]] <- c(species_run[[t+(2*length(infected_species))]])
      infected_areas_out[[t]] <- c(species_run[[t+(3*length(infected_species))]])
      west_rate_out[[t]] <- c(species_run[[t+(4*length(infected_species))]])
      east_rate_out[[t]] <- c(species_run[[t+(5*length(infected_species))]])
      north_rate_out[[t]] <- c(species_run[[t+(6*length(infected_species))]])
      south_rate_out[[t]] <- c(species_run[[t+(7*length(infected_species))]])
    }
    print("outer loop")
    
    outputs <- list(infections_out, susceptibles_out, number_infecteds_out, infected_areas_out, west_rate_out, east_rate_out, north_rate_out, south_rate_out, treatment)
  }
  
  infecteds <- list()
  susceptibles <- list()
  probabilities <- list()
  number_infecteds <- list()
  infected_areas <- list()
  west_rates <- list()
  east_rates <- list()
  north_rates <- list()
  south_rates <- list()
  treatments <- list()
  
  for (sp in 1:length(infected_files)) {
    infecs <- list()
    susceps <- list()
    probs <- list()
    number_infects <- list()
    infected_ars <- list()
    west_rats <- list()
    east_rats <- list()
    north_rats <- list()
    south_rats <- list()
    for (l in 1:years_simulated) {
      
      infecs[[l]] <- run_years[[l]][[sp]][[1]]
      susceps[[l]] <- run_years[[l + years_simulated]][[sp]][[1]]
      number_infects[[l]] <- run_years[[l + years_simulated * 2]][[sp]]
      infected_ars[[l]] <- run_years[[l + years_simulated * 3]][[sp]]
      west_rats[[l]] <- run_years[[l + years_simulated * 4]][[sp]]
      east_rats[[l]] <- run_years[[l + years_simulated * 5]][[sp]]
      north_rats[[l]] <- run_years[[l + years_simulated * 6]][[sp]]
      south_rats[[l]] <- run_years[[l + years_simulated * 7]][[sp]]
      treatments[[l]] <- run_years[[l + years_simulated * 8]][[1]]
    }
    names(infecs) <- year_names
    names(susceps) <- year_names
    names(number_infects) <- year_names
    names(infected_ars) <- year_names
    names(west_rats) <- year_names
    names(east_rats) <- year_names
    names(north_rats) <- year_names
    names(south_rats) <- year_names
    names(treatments) <- year_names
    infecteds[[sp]] <- infecs
    susceptibles[[sp]] <- susceps
    probabilities[[sp]] <- probs
    number_infecteds[[sp]] <- number_infects
    infected_areas[[sp]] <- infected_ars
    west_rates[[sp]] <- west_rats
    east_rates[[sp]] <- east_rats
    north_rates[[sp]] <- north_rats
    south_rates[[sp]] <- south_rats
  }
  
  names(infecteds) <- species
  names(susceptibles) <- species
  names(number_infecteds) <- species
  names(infected_areas) <- species
  names(west_rates) <- species
  names(east_rates) <- species
  names(north_rates) <- species
  names(south_rates) <- species
  
  outputs <- list(infecteds, susceptibles, number_infecteds, infected_areas, west_rates, east_rates, north_rates, south_rates, treatments)
  names(outputs) <- c('infecteds', 'susceptibles', 'number_infecteds', 'infected_areas', 'west_rates', 'east_rates', 'north_rates', 'south_rates', 'treatments')
  
  outputs
}

stopCluster(cl)




treatments_2020_stack = stack()
treatments_2021_stack = stack()
treatments_2022_stack = stack()
treatments_2023_stack = stack()
treatments_2024_stack = stack()

EU1_2020_stack = stack()
EU1_2021_stack = stack()
EU1_2022_stack = stack()
EU1_2023_stack = stack()
EU1_2024_stack = stack()

NA1_2020_stack = stack()
NA1_2021_stack = stack()
NA1_2022_stack = stack()
NA1_2023_stack = stack()
NA1_2024_stack = stack()


for(i in 1:500){
  both = runs[[i]]
  eu_list = both$SOD_EU1
  na_list = both$SOD_NA1
  eu = eu_list$`2024`
  eu_2020 = subset(eu, 1:1)
  eu_2021 = subset(eu, 2:2)
  eu_2022 = subset(eu, 3:3)
  eu_2023 = subset(eu, 4:4)
  eu_2024 = subset(eu, 5:5)
  na = na_list$`2024`
  na_2020 = subset(na, 1:1)
  na_2021 = subset(na, 2:2)
  na_2022 = subset(na, 3:3)
  na_2023 = subset(na, 4:4)
  na_2024 = subset(na, 5:5)
  EU1_2020_stack = stack(EU1_2020_stack, eu_2020)
  EU1_2021_stack = stack(EU1_2021_stack, eu_2021)
  EU1_2022_stack = stack(EU1_2022_stack, eu_2022)
  EU1_2023_stack = stack(EU1_2023_stack, eu_2023)
  EU1_2024_stack = stack(EU1_2024_stack, eu_2024)
  
  NA1_2020_stack = stack(NA1_2020_stack, na_2020)
  NA1_2021_stack = stack(NA1_2021_stack, na_2021)
  NA1_2022_stack = stack(NA1_2022_stack, na_2022)
  NA1_2023_stack = stack(NA1_2023_stack, na_2023)
  NA1_2024_stack = stack(NA1_2024_stack, na_2024)
}

for(t in 4001:4500){
  treat = runs[[t]]
  treat_2020 = treat$`2020`
  treat_2021 = treat$`2021`
  treat_2022 = treat$`2022`
  treat_2023 = treat$`2023`
  treat_2024 = treat$`2024`
  treatments_2020_stack = stack(treatments_2020_stack, treat_2020)
  treatments_2021_stack = stack(treatments_2021_stack, treat_2021)
  treatments_2022_stack = stack(treatments_2022_stack, treat_2022)
  treatments_2023_stack = stack(treatments_2023_stack, treat_2023)
  treatments_2024_stack = stack(treatments_2024_stack, treat_2024)
}


writeRaster(treatments_2020_stack, "/EU1 Priority Results/2.5 mill 300ft buff/treatments_2020.tif", format="GTiff")
writeRaster(treatments_2021_stack, "/EU1 Priority Results/2.5 mill 300ft buff/treatments_2021.tif", format="GTiff")
writeRaster(treatments_2022_stack, "/EU1 Priority Results/2.5 mill 300ft buff/treatments_2022.tif", format="GTiff")
writeRaster(treatments_2023_stack, "/EU1 Priority Results/2.5 mill 300ft buff/treatments_2023.tif", format="GTiff")
writeRaster(treatments_2024_stack, "/EU1 Priority Results/2.5 mill 300ft buff/treatments_2024.tif", format="GTiff")

writeRaster(EU1_2020_stack, "/EU1 Priority Results/2.5 mill 300ft buff/EU1_2020.tif", format="GTiff")
writeRaster(EU1_2021_stack, "/EU1 Priority Results/2.5 mill 300ft buff/EU1_2021.tif", format="GTiff")
writeRaster(EU1_2022_stack, "/EU1 Priority Results/2.5 mill 300ft buff/EU1_2022.tif", format="GTiff")
writeRaster(EU1_2023_stack, "/EU1 Priority Results/2.5 mill 300ft buff/EU1_2023.tif", format="GTiff")
writeRaster(EU1_2024_stack, "/EU1 Priority Results/2.5 mill 300ft buff/EU1_2024.tif", format="GTiff")

writeRaster(NA1_2020_stack, "/EU1 Priority Results/2.5 mill 300ft buff/NA1_2020.tif", format="GTiff")
writeRaster(NA1_2021_stack, "/EU1 Priority Results/2.5 mill 300ft buff/NA1_2021.tif", format="GTiff")
writeRaster(NA1_2022_stack, "/EU1 Priority Results/2.5 mill 300ft buff/NA1_2022.tif", format="GTiff")
writeRaster(NA1_2023_stack, "/EU1 Priority Results/2.5 mill 300ft buff/NA1_2023.tif", format="GTiff")
writeRaster(NA1_2024_stack, "/EU1 Priority Results/2.5 mill 300ft buff/NA1_2024.tif", format="GTiff")

sim_end_time <- Sys.time()
time_taken <- sim_end_time - sim_start_time
print(paste0("SIMULATION TIME:  ", time_taken))
