# Butterfly response classification and comparison
# Reid Steele, 11/21/2024

# Example species:
# Shifting: Vanessa atalanta
# Decoupling: Coenonympha pamphilus
# Combination: Thymelicus lineola

# Libraries
library(tidyverse)
library(survival)
library(lubridate)
library(plotly)

# Run intialization
source('Bfly Initial.R')

# Prepped data for testing
data_quants = data_prep(bfly_f, 'Vanessa atalanta')

########################################################################################

# # Plot survival curves
# ggplot(data_quants, aes(x = quantile, y = surv_y)) + geom_line(linewidth = 2) +
#   labs(color = 'Year', x = 'Day of Year', y = 'Proportion to Arrive') +
#   geom_point(data = data_quants, aes(x = FIRSTDAY, y = abs(surv_y), color = as.factor(YEAR)), size = 3) + 
#   theme(text = element_text(size = 20))
# 
# # Plot survival residuals over time
# ggplot(data_quants, (aes(y = residual, x = YEAR))) + geom_point() +
#   labs(x = 'Year', y = 'Difference from Mean Arrival')  + 
#   stat_smooth(method = 'glm') + 
#   theme(text = element_text(size = 20))
# 
# # Plot arrival over time
# ggplot(data_quants, (aes(y = FIRSTDAY, x = YEAR))) + geom_point() +
#   labs(x = 'Year', y = 'Day of Arrival')  + 
#   stat_smooth(method = 'glm') + 
#   theme(text = element_text(size = 20))
# 
# # Plot temerature over time
# ggplot(data_quants, (aes(y = Temp, x = YEAR))) + geom_point() +
#   labs(x = 'Year', y = 'Temperature at Arrival')  + 
#   stat_smooth(method = 'glm') + 
#   theme(text = element_text(size = 20))
# 
# # Plot temerature at arrival
# ggplot(data_quants, (aes(y = Temp, x = FIRSTDAY, color = YEAR))) + geom_point() +
#   labs(x = 'Day of Arrival', y = 'Temperature at Arrival')  + 
#   stat_smooth(method = 'glm') + 
#   theme(text = element_text(size = 20))
# 
# 
# # Run linear models, extract slopes
# e_slope = summary(lm(data = data_quants, FIRSTDAY ~ YEAR))$coefficients['YEAR',]
# t_slope = summary(lm(data = data_quants, Temp ~ YEAR))$coefficients['YEAR',]


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
  
# Run function
slopes = classify(bfly_f)

# Load in trait data
traits = read.csv('./Butterflies/bfly_traits/data/ecological_traits_2022.csv', skip = 1)

# Fix text encoding
traits$scientific_name = iconv(traits$scientific_name, from = "ISO-8859-1", to = "UTF-8")
traits$scientific_name = gsub('\xa0', ' ', traits$scientific_name, useBytes = T)

# Join classifications to traits
class = left_join(slopes, traits, by = join_by(x$species == y$scientific_name)) %>%
  select(species, class, long_term_distribution_trend_gb) %>% 
  filter((is.na(long_term_distribution_trend_gb) == F) & (nchar(long_term_distribution_trend_gb) > 0)) %>%
  mutate(long_term_distribution_trend_gb = as.numeric(gsub('%', '', long_term_distribution_trend_gb))) %>%
  filter(class != 'no change')

# Plot abundance by class
ggplot(class, aes(x = class, y = long_term_distribution_trend_gb, group = class, fill = class)) +
  geom_boxplot()

# Filter out NAs in abundance trend
test = aov(long_term_distribution_trend_gb ~ class, data = class)
summary(test)
tukey = TukeyHSD(test)
plot(tukey)

# Reformat to long format
slopes_l = pivot_longer(slopes, cols = colnames(slopes)[c(-1,-length(colnames(slopes)))])

# plot slopes
int = ggplot(slopes, aes(x = event, y = temp, color = class, text = species)) +
  theme_classic() + geom_hline(yintercept = 0, linetype = 'dashed') + 
  geom_vline(xintercept = 0, linetype = 'dashed') +  
  geom_abline(slope = 1 , intercept = 0) +
  geom_point(size = 2)
ggplotly(int)

