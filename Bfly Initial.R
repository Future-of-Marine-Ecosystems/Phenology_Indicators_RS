# Butterfly functions and initialization
# Reid Steele, 11/21/2024

# Libraries
library(tidyverse)
library(survival)
library(lubridate)
library(plotly)

# Working Directory
# setwd('C:/Users/R3686/OneDrive - Dalhousie University/Documents/Acoustic Tracking/Data/Butterflies')

# Read in data
bfly = read.csv('./Butterflies/data/ukbmsphenology2022.csv')
temp = read.csv('./Butterflies/data/uktemps.csv')

# Rename temperature column
colnames(temp) = c('Date', 'Temp')

# Format dates
temp$Date = strptime(temp$Date, format = '%Y-%m-%d', tz = 'UTC')

# Add year
temp$YEAR = year(temp$Date)

# Filter to years in butterfly data
temp = filter(temp, temp$YEAR %in% bfly$YEAR)

# Add julianday
temp$yday = yday(temp$Date)

# Correct to butterfly dates
temp$FIRSTDAY = ifelse(leap_year(temp$YEAR), temp$yday-91, temp$yday-90)

# Filter temp
temp = temp[,c('Temp', 'YEAR', 'FIRSTDAY')]

# Calculate mean temperature
temp_mean = group_by(temp, YEAR) %>% summarize(temp = mean(Temp)) %>% mutate(temp_norm = (temp - mean(temp))/sd(temp))

# Run linear model
temp_lm = lm(data = temp_mean, temp_norm ~ YEAR)
summary(temp_lm)

# Plot normalize temperature
ggplot(temp_mean, aes(x = YEAR, y = temp_norm)) + geom_line() + geom_point() + stat_smooth(method = 'lm')

# Join temperature
bfly = left_join(bfly, temp, by = c('YEAR', 'FIRSTDAY'))

# Plot butterflies

# Filter first brood
bfly_0 = filter(bfly, NUMBER_OF_BROODS == 1)
bfly_1 = filter(bfly, BROOD == 1)
bfly_f = rbind(bfly_1, bfly_0)
bfly_f = arrange(bfly_f, YEAR) %>% filter(YEAR >= 1979)

# # plot
# ggplot(bfly_f, aes(x = YEAR, y = FIRSTDAY)) + geom_point() +
#   stat_smooth(method = 'glm') + facet_wrap(~SPECIES_NAME)

# Function for prepping data
data_prep = function(input, species){
  
  # Filter data
  data = filter(input, SPECIES_NAME == species) # data
  
  # # Remove single record years
  # data = group_by(data, YEAR)
  # data = data %>% filter(!(YEAR %in% group_keys(data)[[1]][which(group_size(data) == 1)]))
  
  # Calculate survival in-year
  df = NULL
  for(i in 1:length(unique(data$YEAR))){
    
    dy = filter(data, YEAR == unique(data$YEAR)[i])
    
    s = survfit(formula = Surv(FIRSTDAY) ~ 1, data = dy)
    
    df = rbind(df, data.frame(YEAR = unique(data$YEAR)[i], FIRSTDAY = s$time, surv_y = s$surv))
    
  } # end year loop
  
  # Join survival in-year
  data = left_join(data, df)
  
  
  # Rate of change
  
  # Survival model by year
  surv_by = survfit(formula = Surv(FIRSTDAY) ~ YEAR, data = data)
  
  # Calculate quantiles on all data
  quant_all = quantile(surv_by, abs(1-data$surv_y))$quantile
  
  # Fix row names
  rownames(quant_all) = gsub('year=', '', rownames(quant_all))
  
  # Calculate average quantile
  quant_avg = colMeans(quant_all)
  
  # Join average quantile to original data
  data_quants = cbind(data, quantile = quant_avg)
  
  # Calculate residual
  data_quants$residual = data_quants$quantile - data_quants$FIRSTDAY
  
  # Return
  return(data_quants)
  
} # End data_prep function


# Classify phenological responses
classify = function(data){
  
  # Slope storage object
  slopes = data.frame(species = NULL, event = NULL, temp = NULL)
  
  # Slope loop
  for(i in 1:length(unique(data$SPECIES_NAME))){
    
    # Prep data
    data_s = data_prep(data, unique(data$SPECIES_NAME)[i])
    
    # Normalize data
    data_s$FIRSTDAY = (data_s$FIRSTDAY - mean(data_s$FIRSTDAY))/sd(data_s$FIRSTDAY)
    data_s$Temp = (data_s$Temp - mean(data_s$Temp))/sd(data_s$Temp)
    
    # Run linear models, extract slopes
    e_summ = summary(lm(data = data_s, FIRSTDAY ~ YEAR))
    t_summ = summary(lm(data = data_s, Temp ~ YEAR))
    
    # Run linear models, extract slopes
    e_slope = e_summ$coefficients['YEAR','Estimate']
    t_slope = t_summ$coefficients['YEAR','Estimate']
    
    # Run linear models, extract slopes
    e_p = e_summ$coefficients['YEAR','Pr(>|t|)']
    t_p = t_summ$coefficients['YEAR','Pr(>|t|)']
    
    # Object to add to data frame
    add = data.frame(species = unique(data$SPECIES_NAME)[i], 
                     event = e_slope, temp = t_slope, event_p = e_p, temp_p = t_p)
    
    # Add to data frame
    slopes = rbind(slopes, add)
    
  } # End slope loop
  
  # # Calculate difference in slopes
  # slopes$diff = slopes$event+slopes$temp
  # slopes$diff_n = slopes$event+temp_lm$coefficients['YEAR']
  
  # Create significance column
  slopes$event_sig = ifelse(slopes$event_p < 0.05, TRUE, FALSE)
  slopes$temp_sig = ifelse(slopes$temp_p < 0.05, TRUE, FALSE)
  
  # Create classification column
  slopes$class = ifelse((slopes$event_sig == TRUE) & (slopes$temp_sig == TRUE), 'combination', NA)
  slopes$class = ifelse((slopes$event_sig == FALSE) & (slopes$temp_sig == FALSE), 'no change', slopes$class)
  slopes$class = ifelse((slopes$event_sig == TRUE) & (slopes$temp_sig == FALSE), 'shifting', slopes$class)
  slopes$class = ifelse((slopes$event_sig == FALSE) & (slopes$temp_sig == TRUE), 'decoupling', slopes$class)
  
  # Classify overshoots
  slopes$class = ifelse((slopes$class == 'combination') & (slopes$temp < 0), 'shifting', slopes$class)
  slopes$class = ifelse((slopes$class == 'combination') & (slopes$event > 0), 'decoupling', slopes$class)
  
  # Return
  return(slopes)
  
} # End classify function
