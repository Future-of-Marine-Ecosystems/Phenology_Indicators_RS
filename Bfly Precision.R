# Butterfly response power testing
# Reid Steele, 11/27/2024

# Example species:
# Shifting: Vanessa atalanta
# Decoupling: Coenonympha pamphilus
# Combination: Thymelicus lineola

# Libraries
library(tidyverse)
library(survival)
library(lubridate)
library(plotly)
library(data.table)
library(runner)

# Run intialization
source('Bfly Initial.R')

# Prepped data for testing
data = data_prep(bfly_f, 'Vanessa atalanta')

########################################################################################

# Generate a predictability curve
pred_curve = function(data, min_y = 5){
  
  # Gather final ("correct") classification
  class_f = classify(data)$class
  
  # Container object
  class_prec = NULL
  
  # Cycle through truncated time series length
  for(i in min_y:length(unique(data$YEAR))){
    
    # Generate test windows
    windows = runner(1:length(unique(data$YEAR)), k = i)[i:length(unique(data$YEAR))]

    # Container object
    class_j = NULL
    
    # Loop through windows
    for(j in 1:length(windows)){
      
      # Reformat window of interest
      win = unlist(windows[j])
      
      # Filter data
      data_f = filter(data, YEAR %in% unique(data$YEAR)[win])
      
      # Run classification and contain
      class_j = c(class_j, classify(data_f)$class)
      
    } # End loop through windows
    
    # Calculate correctness 
    class_i = cbind(tsl = i, correct = sum(class_j == class_f)/length(class_j))
    
    # Concatenate result
    class_prec = rbind(class_prec, class_i)
    
  } # End loop through truncated time series
  
  # return
  return(class_prec)
  
} # End pred_curve function


# Generate a predictability curve
pred_curve_ks = function(data, min_y = 5){
  
  # Gather final ("correct") classification
  
  # Run emergence tests
  em = max(ks_emergence(data, plot = F)$emerged, na.rm = T)
  dc = max(ks_decouple(data, plot = F)$emerged, na.rm = T)
  
  # Assign class
  class_f = em + dc * 2
  if (class_f == 0) {class_f = 'no change'}
  if (class_f == 1) {class_f = 'shifting'}
  if (class_f == 2) {class_f = 'decoupling'}
  if (class_f == 3) {class_f = 'combination'}
  
  # Container object
  class_prec = NULL
  
  # Cycle through truncated time series length
  for(i in min_y:length(unique(data$YEAR))){
    
    # Generate test windows
    windows = runner(1:length(unique(data$YEAR)), k = i)[i:length(unique(data$YEAR))]
    
    # Container object
    class_j = NULL
    
    # Loop through windows
    for(j in 1:length(windows)){
      
      # Reformat window of interest
      win = unlist(windows[j])
      
      # Filter data
      data_f = filter(data, YEAR %in% unique(data$YEAR)[win])
      
      # Run emergence tests
      em = max(ks_emergence(data_f, plot = F)$emerged, na.rm = T)
      dc = max(ks_decouple(data_f, plot = F)$emerged, na.rm = T)
      
      # Assign class
      class = em + dc*2
      if(class == 0){class = 'no change'}
      if(class == 1){class = 'shifting'}
      if(class == 2){class = 'decoupling'}
      if(class == 3){class = 'combination'}
      
      # Run classification and contain
      class_j = c(class_j, class)
      
    } # End loop through windows
    
    # Calculate correctness 
    class_i = cbind(tsl = i, correct = sum(class_j == class_f)/length(class_j))
    
    # Concatenate result
    class_prec = rbind(class_prec, class_i)
    
  } # End loop through truncated time series
  
  # return
  return(class_prec)
  
} # End pred_curve function


tt1 = Sys.time()
t1 = pred_curve_ks(data)
Sys.time() - tt1

tt1 = Sys.time()
t3 = pred_curve(data)
Sys.time() - tt1


# Loop through species

# Container object
precision = NULL
tot = length(unique(bfly_f$SPECIES_NAME))

# Loop
for(s in unique(bfly_f$SPECIES_NAME)){
  
  # Run data prep
  sp = data_prep(bfly_f, s)
  
  # Run precision tests
  lin = cbind(pred_curve(sp), species = s, test = 'linear trend')
  kst = cbind(pred_curve_ks(sp), species = s, test = 'ks test')
  
  # Gather result
  precision = rbind(precision, lin, kst)
  
  # Paste Progress
  cur = match(s, unique(bfly_f$SPECIES_NAME))
  cat(paste(cur, 'of', tot))
  
} # End loop through species

# backup = precision

# Set to data framem
precision = as.data.frame(precision)

# Fix non-numerics
precision$tsl = as.numeric(precision$tsl)
precision$correct = as.numeric(precision$correct)

# Calculate mean curves across all species
precision_summ = precision %>% group_by(tsl, test) %>%
  summarize(prop = mean(correct), prop_sd = sd(correct))

# Plot curves
ggplot(precision_summ, aes(x = tsl, y = prop, color = test)) + geom_line(linewidth = 1) +
  geom_ribbon(aes(y = prop, fill = test, ymin = prop - prop_sd, ymax = prop + prop_sd), alpha = 0.1, colour = NA) +
  theme_classic() + coord_cartesian(ylim = c(0,1)) +
  labs(x = 'Time Series Length', y = 'Proportion Matching Final Classification', 
       color = 'Classification Test', fill = 'Classification Test')



