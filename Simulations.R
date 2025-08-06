# Simulations
library(doParallel)
library(foreach)
library(doSNOW)
library(phenometrics)
library(dplyr)

# source('TOE.R')

# Samples
samples = seq(5, 255, 5)

# Spread
spread = seq(4.42, 25, 0.42)

# Trend
trend = seq(-0.01, -1.01, -0.02)

# Years
years = seq(1, 50, 1)

# baseline
baseline = 150

# annual spread as a ratio of mean spread
sfact = 1.5

# Number of tests to run
ntests = 100

# test values
tspread = spread[round(length(spread)/2)]
tsamples = samples[round(length(samples)/2)]
ttrend = trend[round(length(trend)/2)]

# Simulation parameters
sim_params = expand.grid(samples, trend, spread); colnames(sim_params) = c('samples', 'trend', 'spread')
# sim_params$year = NA
# sim_params$year2 = NA
sim_params$row = 1:nrow(sim_params)

# Filter out rows where spread, samples, or trend are not test vales
sim_params = filter(sim_params, (spread == tspread) | (samples == tsamples) | (trend == ttrend))







# Generate simulated datasets
sim_dataset = function(row){
    
  # Generate test dataset
  means = rnorm(length(years), baseline, row$spread) + row$trend*(1:length(years))
    
  simdata = NULL
  for(y in 1:length(years)){
      
      # simdata = rbind(simdata, cbind(event = rnorm(row$samples, baseline+(years[y]*row$trend) #+ rnorm(1, 0, jiggle*row$spread)
      #                                              ,row$spread), year = years[y]))
      
    simdata = rbind(simdata, cbind(event = rnorm(row$samples, means[y], row$spread*sfact), year = years[y]))
      
  }
    
    # Change negatives to 0
  simdata = ifelse(simdata < 0, 0, simdata)
    
  # return
  return(simdata)

}


# calculate emergence test
sim_emergence = function(row){
  
  # Container objects
  eyear_stat = rep(NA, ntests)
  eyear_emp = rep(NA, ntests)
  results = data.frame(eyear_stat = NA, eyear_emp = NA)
  
  # Loop through tests
  for(i in 1:ntests){
    
    # Simulate data
    simdata = sim_dataset(row)
    
    # Run emergence tests
    em_stat = ks_tope(as.data.frame(simdata), plot = F)
    em_emp = emp_tope(as.data.frame(simdata), plot = F)
    
    # calculate minimum emergence year
    if(1 %in% em_stat$emerged) {eyear_stat[i] = min(em_stat[em_stat$emerged == 1,'year'], na.rm = T)} else {eyear_stat[i] = NA}
    if(1 %in% em_emp$emerged) {eyear_emp[i] = min(em_emp[em_emp$emerged == 1,'year'], na.rm = T)} else {eyear_emp[i] = NA}
  
  } # End loop through tests

  # if more than half are NA, set NA, else set to mean
  if(sum(is.na(eyear_stat)) > (ntests/2)){results$eyear_stat = NA} else {results$eyear_stat = mean(eyear_stat, na.rm = T)}
  if(sum(is.na(eyear_emp)) > (ntests/2)){results$eyear_emp = NA} else {results$eyear_emp = mean(eyear_emp, na.rm = T)}
  
  # return
  return(results)
  
}

# Set up cluster
cl = makeCluster(detectCores()-1, outfile = '')
registerDoSNOW(cl)

# Set up progress bar
pb <- txtProgressBar(min = 1, max = nrow(sim_params), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Loop through sim_params
start = Sys.time()

# set seed
set.seed(123)

# Loop through sim_params rows (efficiently)
sim_results = NULL
sim_results = foreach(i = 1:nrow(sim_params), .options.snow = opts, .combine = rbind) %dopar% {
  
  # Progress bar
  setTxtProgressBar(pb, i)
  
  # Library
  library(phenometrics)
  
  # Calculate mean emergences
  sim_emergence(sim_params[i,])

  
}

Sys.time() - start

stopCluster(cl)

# Combine
sim_params = cbind(sim_params, sim_results)

save.image('sim_results_comb.RData')


