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
# # Plot temperature over time
# ggplot(data_quants, (aes(y = Temp, x = YEAR))) + geom_point() +
#   labs(x = 'Year', y = 'Temperature at Arrival')  + 
#   stat_smooth(method = 'glm') + 
#   theme(text = element_text(size = 20))
# 
# # Plot temperature at arrival
# ggplot(data_quants, (aes(y = Temp, x = FIRSTDAY, color = YEAR))) + geom_point() +
#   labs(x = 'Day of Arrival', y = 'Temperature at Arrival')  + 
#   stat_smooth(method = 'glm') + 
#   theme(text = element_text(size = 20))
# 
# 
# # Run linear models, extract slopes
# e_slope = summary(lm(data = data_quants, FIRSTDAY ~ YEAR))$coefficients['YEAR',]
# t_slope = summary(lm(data = data_quants, Temp ~ YEAR))$coefficients['YEAR',]


# Loop through species
slopes = NULL
for(i in 1:length(unique(bfly_f$SPECIES_NAME))){
  
  # Run data_prep
  data = data_prep(bfly_f, unique(bfly_f$SPECIES_NAME)[i])
  
  # Classify
  slopes = rbind(slopes, classify(data))
  
} # End classification loop

# Load in trait data
traits = read.csv('./Butterflies/bfly_traits/data/ecological_traits_2022.csv', skip = 1)

# Fix text encoding
traits$scientific_name = iconv(traits$scientific_name, from = "ISO-8859-1", to = "UTF-8")
traits$scientific_name = gsub('\xa0', ' ', traits$scientific_name, useBytes = T)

# Join classifications to traits
class = left_join(slopes, traits, by = join_by(x$species == y$scientific_name)) %>%
  select(species, class, forewing_maximum) %>% 
  filter((is.na(forewing_maximum) == F) & (nchar(forewing_maximum) > 0)) %>%
  mutate(forewing_maximum = as.numeric(gsub('%', '', forewing_maximum))) %>%
  filter(class != 'no change')

# Join classifications to traits (categorical)
class = left_join(slopes, traits, by = join_by(x$species == y$scientific_name)) %>%
  select(species, class, obligate_univoltine) %>% 
  mutate(cat = as.factor(ifelse(is.na(obligate_univoltine), 0 , 1))) %>%
  filter(class != 'no change')

# Plot abundance by class
ggplot(class, aes(x = class, y = cat, group = class, fill = class)) +
  geom_boxplot()

# Filter out NAs in abundance trend
test = aov(obligate_univoltine ~ class, data = class)
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

