# Simulations
library(doParallel)
library(foreach)
library(doSNOW)
library(phenometrics)
library(tidyverse)

# Samples
samples = 130
# samples = 20

# Spread
spread = 14.710

# Trend
trend = -0.51

# Years
years = seq(1, 100, 1)

# Baseline length
base_l = seq(5, 50, 1)

# baseline
baseline = 150

# annual spread as a ratio of mean spread
sfact = 1.5

# Number of tests to run
ntests = 1000

# Emergence thresholds
emt = c(5,10,15)

# Simulation parameters
sim_params_base = expand.grid(samples, trend, spread); colnames(sim_params_base) = c('samples', 'trend', 'spread')
# sim_params_base = expand.grid(samples, trend, spread, base_l); colnames(sim_params_base) = c('samples', 'trend', 'spread', 'base_l')
# sim_params_base$year = NA
# sim_params_base$year2 = NA

# # Reduce combinations
# sim_params_base = sim_params_base %>% 
#   filter(spread != max(spread) & spread != min(spread)) %>%
#   filter(trend != max(trend) & trend != min(trend))
  
# Add row column
sim_params_base$row = 1:nrow(sim_params_base)










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
sim_emergence = function(row, base_l){
  
  # Container objects
  eyear_stat = matrix(NA, ncol = ntests, nrow = length(base_l))
  eyear_stat_static = matrix(NA, ncol = ntests, nrow = length(base_l))
  eyear_emp = matrix(NA, ncol = ntests, nrow = length(base_l))
  eyear_emp_static = matrix(NA, ncol = ntests, nrow = length(base_l))
  # results = data.frame(row, eyear_stat = NA, eyear_emp = NA, eyear_emp_static = NA, prop_stat = NA, prop_emp = NA, prop_emp_static = NA, base_l = base_l)
  
  # Loop through tests
  for(i in 1:ntests){
    
    # Simulate data
    simdata = as.data.frame(sim_dataset(row))
    
    # Loop through base_l
    for(j in 1:length(base_l)){
    
      # Run emergence tests
      em_stat = ks_tope(simdata[simdata$year <= base_l[j],], plot = F, emt = 5)
      em_stat_static = ks_tope(simdata, plot = F, subset = which(simdata$year <= base_l[j]), emt = 5)
      em_emp = emp_tope(simdata[simdata$year <= base_l[j],], plot = F, emt = 5)
      em_emp_static = emp_tope(simdata, plot = F, max_y = base_l[j], emt = 5)
      
      # calculate minimum emergence year
      if(1 %in% em_stat$emerged) {eyear_stat[j,i] = min(em_stat[em_stat$emerged == 1,'year'], na.rm = T)}
      if(1 %in% em_stat_static$emerged) {eyear_stat_static[j,i] = min(em_stat_static[em_stat_static$emerged == 1,'year'], na.rm = T)}
      if(1 %in% em_emp$emerged) {eyear_emp[j,i] = min(em_emp[em_emp$emerged == 1,'year'], na.rm = T)}
      if(1 %in% em_emp_static$emerged) {eyear_emp_static[j,i] = min(em_emp_static[em_emp_static$emerged == 1,'year'], na.rm = T)}
      
    }
    
  } # End loop through tests
  
  # # Calculate mean emergence
  # results$eyear_stat = rowMeans(eyear_stat, na.rm = T)
  # results$eyear_stat_static = rowMeans(eyear_stat_static, na.rm = T)
  # results$eyear_emp = rowMeans(eyear_emp, na.rm = T)
  # results$eyear_emp_static = rowMeans(eyear_emp_static, na.rm = T)
  # 
  # # Calculate proportion emerged
  # results$prop_stat = 1-((rowSums(is.na(eyear_stat)))/ntests)
  # results$prop_stat_static = 1-((rowSums(is.na(eyear_stat_static)))/ntests)
  # results$prop_emp = 1-((rowSums(is.na(eyear_emp)))/ntests)
  # results$prop_emp_static = 1-((rowSums(is.na(eyear_emp_static)))/ntests)

  eyear_stat = cbind(base_l, method = 'statistical', eyear_stat)
  eyear_stat_static = cbind(base_l, method = 'statistical_static', eyear_stat_static)
  eyear_emp = cbind(base_l, method = 'empirical', eyear_emp)
  eyear_emp_static = cbind(base_l, method = 'empirical_static', eyear_emp_static)
  outputs = rbind(eyear_stat, eyear_stat_static, eyear_emp, eyear_emp_static)
  
  results = data.frame(row, outputs)
  
  # return
  return(results)
  
}


# calculate emergence test
sim_emergence_par = function(row, base_l, method, emt = 5){
  
  # Set up progress bar
  pb <- txtProgressBar(min = 1, max = ntests, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Loop through tests
  em_out = foreach(i = 1:ntests, .packages = 'phenometrics', .combine = cbind, .options.snow = opts,
                   .export = c('ntests', 'years', 'sim_dataset', 'baseline', 'sfact')) %dopar% {
                     
    # Progress bar
    setTxtProgressBar(pb, i)
    
    # Simulate data
    simdata = as.data.frame(sim_dataset(row))
    
    # Container objects
    eyear = rep(NA, length(base_l)*length(method))
    
    # Loop through methods
    for(k in 1:length(method)){
      
      # loop through base_l
      for(j in 1:length(base_l)){
      
        # Run emergence tests
        if(method[k] == 'statistical'){
          
          # statistical method using only base_l years
          em_out = ks_tope(simdata[simdata$year <= base_l[j],], plot = F, emt = emt)
          
        } else if (method[k] == 'empirical'){
          
          # empirical method using only base_l years
          em_out = emp_tope(simdata[simdata$year <= base_l[j],], plot = F, emt = emt)
          
        } else if (method[k] == 'statistical_static'){
          
          # statistical method using all years but only base_l for baseline
          em_out = ks_tope(simdata, plot = F, subset = which(simdata$year <= base_l[j]), emt = emt)
          
        } else if (method[k] == 'empirical_static'){
          
          # empirical method using all years but only base_l for baseline
          em_out = emp_tope(simdata, plot = F, max_y = base_l[j], emt = emt)
          
        } else {stop('method not recognized')}
  
        # calculate minimum emergence year
        if(1 %in% em_out$emerged) {eyear[((k-1)*length(base_l)) + j] = min(em_out[em_out$emerged == 1,'year'], na.rm = T)}
        
      }
      
    }
    
    # Print for foreach
    eyear
    
  } # End foreach loop through tests
  
  # Initialize reformat
  method_name = NULL
  # Loop through k
  for(k in 1:length(method)){

    # Add method column
    method_name = c(method_name, rep(method[k], length(base_l)))

  }
  
  # rename em_out columns
  colnames(em_out) = 1:ntests
  
  # Print results
  results = data.frame(row, base_l = rep(base_l, length(method)), emt, method = method_name, em_out)
  
  # close progress
  close(pb)
  
  # return
  return(results)
  
} # End sim_emergence_par function







# Set up cluster
cl = makeCluster(detectCores()-1, outfile = '')
registerDoSNOW(cl)



# start time
start = Sys.time()

# Initialize loop
sim_results = NULL

# Run simulations
for(i in 1:length(emt)){
  
  # set seed
  set.seed(123)
  
  # Loop through emts 
  sim_results = rbind(sim_results, sim_emergence_par(sim_params_base, base_l[base_l >= emt[i]], method = c('statistical', 'empirical', 'statistical_static', 'empirical_static'), emt = emt[i]))
  
}

# End time
Sys.time() - start

stopCluster(cl)

# # Set up progress bar
# pb <- txtProgressBar(min = 1, max = nrow(sim_params_base), style = 3)
# progress <- function(n) setTxtProgressBar(pb, n)
# opts <- list(progress = progress)
# 
# # Loop through sim_params_base
# start = Sys.time()
# 

# 
# # Loop through sim_params_base rows (efficiently)
# sim_results = NULL
# sim_results = foreach(i = 1:nrow(sim_params_base), .options.snow = opts, .combine = rbind) %dopar% {
# 
#   # Progress bar
#   setTxtProgressBar(pb, i)
# 
#   # Library
#   library(phenometrics)
# 
#   # Calculate mean emergences
#   sim_emergence(sim_params_base[i,], base_l)
# 
# 
# }
# 
# # If only one row
# # sim_results = sim_emergence(sim_params_base, base_l)
# 
# Sys.time() - start



# save.image('sim_results_baselines_static.RData')

# Combine
# sim_params_base = cbind(sim_params_base, sim_results)

# Convert to long format, calculate relative emergence
sim_results_l = pivot_longer(sim_results, cols = paste0("X", 1:ntests), names_to = 'test', values_to = 'ToE') %>%
  mutate(test = as.numeric(gsub('X', '', test)), ToE = as.numeric(ToE), base_l = as.numeric(base_l)) %>% 
  group_by(row, test, method, emt) %>%
  mutate(add_em = ifelse((is.na(lag(ToE))) & (!is.na(ToE)), 1, 0)) %>%
  mutate(sub_em = ifelse((!is.na(lag(ToE))) & (is.na(ToE)), 1, 0)) %>%
  ungroup()

# # Calculate relevant counts
# sim_results_counts = sim_results_l %>% group_by(spread, samples, trend, row, method, base_l) %>%
#   summarize(prop_em = sum(!is.na(ToE)), add_em = sum(add_em), sub_em = -sum(sub_em)) %>%
#   mutate(add_em_cum = cumsum(add_em), sub_em_cum = cumsum(sub_em))
# 
# # Plot
# ggplot(sim_results_counts, aes(x = base_l, color = method))  + 
#   geom_col(aes(y = add_em), fill = '#F8766D', color = 'black') + 
#   geom_col(aes(y = sub_em), fill = '#00BFC4', color = 'black') +
#   geom_line(aes(y = prop_em), lwd = 2) +
#   ylim(c(-15, 100)) + theme_bw() +
#   facet_wrap(~method)

# Sim results curve
sim_results_curve = sim_results_l %>%
  group_by(method, row, base_l, emt) %>%
  summarize(ToE = mean(ToE, na.rm = T)) %>%
  filter(method %in% c('statistical_static', 'empirical_static')) %>%
  mutate(group = paste(method, emt))
  
ggplot(data = sim_results_curve, aes(x = base_l, y = ToE, color = method, alpha = as.factor(emt), group = group)) + 
  geom_line(linewidth = 2) + theme_classic() +
  scale_alpha_manual(name = 'Emergence Threshold', values = c(0.2, 0.6, 1)) +
  scale_color_discrete(name = 'Method', labels = c('Empirical', 'Statistical')) +
  labs(x = 'Baseline Length', y = 'Time of Emergence (ToE)')







