# Predictability testing

source('Bfly Initial.R')
load('cpr_results.RData')

# # Experimenting with first-differencing
# va = filter(bfly_f, species == 'Pieris brassicae')
# va_f = filter(va, SITENAME %in% va[va$year == 1979, 'SITENAME'])
# 
# ggplot(va, aes(x = year, y = event)) + geom_point() + stat_smooth(method = 'lm')
# ggplot(va_f, aes(x = year, y = event)) + geom_point() + stat_smooth(method = 'lm')
# group_by(va_f, year) %>% summarize(event = sd(event)) %>%
#   ggplot(aes(x = year, y = event)) + geom_point() + stat_smooth(method = 'lm')
# 
# cf = filter(cpr_f, species == 4)
# 
# ggplot(cf, aes(x = year, y = event)) + geom_point() + stat_smooth(method = 'lm')


# Minimum number of years
minyear = 10

# Output
output = NULL

# Loop through species
for(i in 1:length(unique(bfly_f$species))){
  
  # Species data
  species_data = bfly_f[bfly_f$species == unique(bfly_f$species)[i],]
  
  # Time series years
  ts_yrs = (minyear + min(species_data$year) - 1):max(species_data$year)
  
  # Loop through end years
  for(end_yr in ts_yrs){
    
    # Filter to appropriate ts length
    species_data_trunc = species_data[species_data$year <= end_yr,]
    
    # Run emergence test
    res = try(ks_tope(species_data_trunc, plot = F))
    
    # Go next if error
    if(class(res) == "try-error"){next} else {
    
      # Check if any emergences
      if(any(res$emerged == T, na.rm = T)){
        
        # If there is an emergence, retern emergence year
        eyear = min(res[res$emerged == T,'year'], na.rm = T)
        
      } else {
        
        # If not, set NA
        eyear = NA
        
      } # End if/else emergence testing
      
      # Construct data frame row to add
      add = data.frame(species = unique(bfly_f$species)[i], # Species name
                       max_year = end_yr, # Last year of time series
                       ts_length = end_yr - min(species_data$year) + 1, # Time series length
                       eyear = eyear) # Emergence year
      
      # Add and concatenate output
      output = rbind(output, add)
    
    } # End try if/else
    
  } # End end year for loop
  
} # End species for loop

# Plot output
ggplot(output, aes(x = ts_length, y = eyear, color = as.factor(species))) + geom_line()



# Predictability function
predictability = function(species_data, minyears, method = 'empirical', ...){
  
  # Time series years
  ts_yrs = (minyears + min(species_data$year) - 1):max(species_data$year)
  
  # Initialize loop
  output = NULL
  
  # Loop through end years
  for(end_yr in ts_yrs){
    
    # Filter to appropriate ts length
    species_data_trunc = species_data[species_data$year <= end_yr,]
    
    # Run emergence test
    if(method == 'empirical') {
      
      res = try(emp_tope(species_data_trunc, plot = F, ...))
      
    } else if(method == 'statistical'){
      
      res = try(ks_tope(species_data_trunc, plot = F, ...))
      
    } else {stop()}
    
    
    # Go next if error
    if(class(res) == "try-error"){next} else {
      
      # Check if any emergences
      if(any(res$emerged == T, na.rm = T)){
        
        # If there is an emergence, return emergence year
        eyear = min(res[res$emerged == T,'year'], na.rm = T)
        
      } else {
        
        # If not, set NA
        eyear = NA
        
      } # End if/else emergence testing
      
      # Construct data frame row to add
      add = data.frame(max_year = end_yr, # Last year of time series
                       ts_length = end_yr - min(species_data$year) + 1, # Time series length
                       eyear = eyear) # Emergence year
      
      # Add and concatenate output
      output = rbind(output, add)
      
    } # End try if/else
    
  } # End end year for loop
  
  return(output)
  
} # End function
  





output_cpr
output_ks_ally = output
output_ks_window = output
output_emp_ally = output
output_emp_window = output

output_emp_window$emerged = as.numeric(!is.na(output_emp_window$eyear))

ggplot(output_emp_window, aes(x = max_year, y = species, fill = as.factor(emerged))) + geom_tile()








