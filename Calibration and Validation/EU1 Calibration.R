## To estimate parameters for EU1 Sudden Oak Death disease strain in Oregon

## Uses 100m tanoak proportion data processed from LEMMA OSU Group
## Works with PoPS Package

## If PoPS isn't already installed
# library(devtools)
# devtools::install_github("ncsu-landscape-dynamics/rpops", ref="feature/disagreements")

library(PoPS)
library(raster)

## For zenith:
infected_file = "/calibration input data/inf_2016_eu.tif"
host_file = "/calibration input data/tanoak_2015.tif"
total_plants_file = "/calibration input data/max_density.tif"
infected_years_file = "/calibration input data/inf_stack_2017_2018_eu.tif"
treatments_file = "/calibration input data/management_2015_2017.tif"
weather_coefficient_file = "/calibration input data/weather_coef_2016_2017.tif"
treatment_years = c(2016,2017)
management = TRUE


## FOR MCMC:
num_iterations = 100000
start_reproductive_rate <- 3.0
start_natural_distance_scale <- 300
sd_reproductive_rate <- 0.2
sd_natural_distance_scale <- 10
number_of_cores <- 11
success_metric = "quantity"

start.time <- Sys.time()


EU1_lemma_params_start2015_100m_expon <- PoPS::calibrate(infected_years_file, num_iterations, start_reproductive_rate, number_of_cores,
                                              start_natural_distance_scale, sd_reproductive_rate, sd_natural_distance_scale,
                                              infected_file, host_file, total_plants_file, 
                                              temp = TRUE, temperature_coefficient_file = weather_coefficient_file, 
                                              precip = FALSE, precipitation_coefficient_file, 
                                              time_step = "week", reproductive_rate,
                                              season_month_start = 1, season_month_end = 12, 
                                              start_time = 2016, end_time = 2017, 
                                              use_lethal_temperature = FALSE, temperature_file,
                                              lethal_temperature = NA, lethal_temperature_month = 1,
                                              mortality_on = TRUE, mortality_rate = .05, mortality_time_lag = 2, 
                                              management, treatment_years, treatments_file,
                                              treatment_method = "all infected", treatment_month = 12,
                                              percent_natural_dispersal = 1.0,
                                              natural_kernel_type = "exponential", anthropogenic_kernel_type = "cauchy",
                                              natural_distance_scale = 20, anthropogenic_distance_scale = 0.0,
                                              natural_dir = "N", natural_kappa = 3, 
                                              anthropogenic_dir = "NONE", anthropogenic_kappa = 0,
                                              mask = NULL, success_metric = "quantity")

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

sp_rate <- data.frame(table(EU1_lemma_params_start2015_100m_expon$reproductive_rate))
sp_rate <- sp_rate[sp_rate$Freq >0,]

dist <- data.frame(table(EU1_lemma_params_start2015_100m_expon$natural_distance_scale))
dist <- dist[dist$Freq >0,]

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

most_common_sp_rate <- getmode(EU1_lemma_params_start2015_100m_expon$reproductive_rate)
print(most_common_sp_rate)

most_common_dist <- getmode(EU1_lemma_params_start2015_100m_expon$natural_distance_scale)
print(most_common_dist)

## to find best combo of params:
parameter_sets <- as.data.frame(table(EU1_lemma_params_start2015_100m_expon[,9:10]))