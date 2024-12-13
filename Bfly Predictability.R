# Butterfly response predictability testing
# Reid Steele, 11/28/2024

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
library(spatstat.utils)

# Run intialization
source('Bfly Initial.R')

# Prepped data for testing
data = data_prep(bfly_f, 'Vanessa atalanta')

# # First 10 years
# data = filter(data, YEAR >= max(YEAR - 6))

# # Calculate mean temperature
# temp_mean = temp %>% group_by(YEAR) %>%
#   summarize(Temp_m = mean(Temp)) %>% mutate(Temp_a = scale(Temp_m))
# 
# # Mean data
# data_mean = data %>% group_by(YEAR) %>%
#   summarize(Arr_m = mean(event)) %>% mutate(Arr_a = scale(Arr_m))
# 
# # Plot
# plot(diff(Temp_a) ~ YEAR[-1], data = temp_mean, type = 'l', lwd = 2)
# arr_mean = data %>% group_by(YEAR) %>%
#   summarize(Arr_m = mean(event)) %>% mutate(Arr_a = scale(Arr_m))
# lines(diff(Arr_a) ~ YEAR[-1], data = data_mean, type = 'l', col = 'red', lwd = 2)
# 
# cor(diff(temp_mean$Temp_a)[-1:-6], diff(data_mean$Arr_a))

# # Plot arrival over time
# ggplot(data_f, (aes(y = event, x = YEAR))) + geom_point() +
#   labs(x = 'Year', y = 'Day of Arrival')  +
#   stat_smooth(method = 'glm') +
#   theme(text = element_text(size = 20))
# 
# # Plot temperature over time
# ggplot(data_f, (aes(y = Temp, x = YEAR))) + geom_point() +
#   labs(x = 'Year', y = 'Temperature at Arrival')  +
#   stat_smooth(method = 'glm') +
#   theme(text = element_text(size = 20))

########################################################################################

# Function to calculate correlation of stochasticity with temperature
stoch_corr = function(data, temp){

  # # Calculate range of arrival times
  # arr_range = range(data$event)
  # 
  # # Filter temperature to arrival times
  # arr_temp_env = filter(temp, inside.range(event, arr_range))
  
  # Dont filter
  arr_temp_env = temp
  
  # Calculate mean temperature,
  arr_temp_env_mean = group_by(arr_temp_env, YEAR) %>% filter(YEAR %in% data$YEAR) %>%
  summarize(Temp_m = mean(Temp)) %>% mutate(Temp_a = scale(Temp_m))# Filter to data years
  
  # Calculate mean arrival
  data_mean = data %>% group_by(YEAR) %>% summarize(Arr_m = mean(event)) %>%
    mutate(Arr_a = scale(Arr_m))
  
  # Join into single data frame
  data_cor = left_join(data_mean, arr_temp_env_mean)
  
  # # Plot
  # plot(diff(Temp_a*-1) ~ YEAR[-1], data = data_cor, type = 'l', lwd = 2)
  # lines(diff(Arr_a) ~ YEAR[-1], data = data_cor, type = 'l', lwd = 2, col = 'red')
  
  # # Detrend data
  # data_d = data.frame(year = data_cor$YEAR[-1], Temp_d = diff(data_cor$Temp_a), Arr_d = diff(data_cor$event))
  
  # # Run lm perm test
  # perm = lmp(Arr_d ~ Temp_d, data = data_d)
  # 
  # # Pull p
  # p = summary(perm)$coefficients['Temp_d','Pr(Prob)']
  
  # Calculate correlation between detrended temperature and arrival
  cor = cor(diff(data_cor$Temp_a), diff(data_cor$Arr_a))
  
  # Return
  return(cor)

} # End stoch_corr function



# Function to calculate correlation of stochasticity with temperature
stoch_corr_perm = function(data, temp, nsim = 1000){
  
  # # Calculate range of arrival times
  # arr_range = range(data$event)
  # 
  # # Filter temperature to arrival times
  # arr_temp_env = filter(temp, inside.range(event, arr_range))
  
  # Dont filter
  arr_temp_env = temp
  
  # Calculate mean temperature,
  arr_temp_env_mean = group_by(arr_temp_env, YEAR) %>% filter(YEAR %in% data$YEAR) %>%
    summarize(Temp_m = mean(Temp)) %>% mutate(Temp_a = scale(Temp_m))# Filter to data years
  
  # Calculate mean arrival
  data_mean = data %>% group_by(YEAR) %>% summarize(Arr_m = mean(event)) %>%
    mutate(Arr_a = scale(Arr_m))
  
  # Join into single data frame
  data_cor = left_join(data_mean, arr_temp_env_mean)
  
  # Loop through permutations
  res <- numeric(nsim) # container for residuals
  for (i in 1:nsim) {
    
    ## standard approach: scramble response value
    perm <- sample(nrow(data_mean))
    bdat <- transform(data_mean,Arr_a=Arr_a[perm])
    
    ## compute & store difference in means; store the value
    res[i] <- cor(diff(data_cor$Temp_a), diff(bdat$Arr_a))^2
  }
  
  obs <- cor(diff(data_cor$Temp_a), diff(data_cor$Arr_a))^2
  
  
  # # Plot
  # plot(diff(Temp_a*-1) ~ YEAR[-1], data = data_cor, type = 'l', lwd = 2)
  # lines(diff(Arr_a) ~ YEAR[-1], data = data_cor, type = 'l', lwd = 2, col = 'red')
  
  # # Detrend data
  # data_d = data.frame(year = data_cor$YEAR[-1], Temp_d = diff(data_cor$Temp_a), Arr_d = diff(data_cor$event))
  
  # # Run lm perm test
  # perm = lmp(Arr_d ~ Temp_d, data = data_d)
  # 
  # # Pull p
  # p = summary(perm)$coefficients['Temp_d','Pr(Prob)']
  
  # Calculate correlation between detrended temperature and arrival
  percentile = ecdf(res)
  cor_quant = percentile(obs)
  
  # Return
  return(cor_quant)
  
} # End stoch_corr_perm function








# Container Object
stoch_species = NULL

# Loop through species
for(i in 1:length(unique(bfly_f$SPECIES_NAME))){
  
  # Run data_prep
  data = data_prep(bfly_f, unique(bfly_f$SPECIES_NAME)[i])
  
  # Gather classification
  a = classify(data)
  
  # Gather stock_species correlation coefficient
  b = stoch_corr(data, temp)
  
  # Combine
  c = cbind(a, b)
  
  # Concatenate
  stoch_species = rbind(stoch_species, c)
  
} # End loop through species


head(stoch_species)

ggplot(stoch_species, aes(x = class, y = b, fill = class)) + geom_boxplot() +
  labs(ylab = 'Correlation')

stoch_species$r2 = stoch_species$b^2

ggplot(stoch_species, aes(x = class, y = r2, fill = class)) + geom_boxplot()

test = aov(b ~ class, data = stoch_species)
summary(test)
TukeyHSD(test)

