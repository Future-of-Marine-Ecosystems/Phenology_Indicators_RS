# UK butterfly phenology indicator testing
# Reid Steele 06/10/2024

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

# plot
ggplot(bfly_f, aes(x = YEAR, y = FIRSTDAY)) + geom_point() +
  stat_smooth(method = 'glm') + facet_wrap(~SPECIES_NAME)

# Moving average quantiles

# Setup
ma = 3 # Moving average length
toe_data = NULL# output object

# Example species:
# Shifting: Vanessa atalanta
# Decoupling: Coenonympha pamphilus
# Combination: Thymelicus lineola

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







# Prepped data for testing
data = data_prep(bfly_f, 'Coenonympha pamphilus')






# Loop through years
for(i in ma:length(unique(data$YEAR))){
  
  # Select years
  years_ma = seq(unique(data$YEAR)[i]-(ma+1), unique(data$YEAR)[i], 1)
  
  # Filter data
  data_f = filter(data, YEAR %in% years_ma)
  
  # Run survival fit
  sfit = survfit(Surv(FIRSTDAY) ~ YEAR, data = data_f)
  
  # Calculate quantiles
  quants = quantile(sfit, abs(1-data_f$surv_y))$quantile
  
  # Average quantiles
  quants_ma = colMeans(quants)
  
  # Output data
  out = cbind(year = max(years_ma), quants_ma, surv_y = data_f$surv_y)
  
  # bind
  toe_data = rbind(toe_data, out)
  
} # End moving average loop

# Set data frame
toe_data = as.data.frame(toe_data)

# Plot moving average lines
ggplot(toe_data, aes(x = quants_ma, y = surv_y, color = year, group = year)) + geom_line() +
  scale_colour_gradient(low = "blue", high = "red")

# Run hypothesis test
surv_diff <- survdiff(Surv(quants_ma) ~ year, data = toe_data_f)
surv_diff

ggplot(data, aes(x = FIRSTDAY, y = surv_y, color = as.factor(YEAR), group = YEAR)) + geom_line()

# No moving average

# Significance matrix - visualize TOE ideas
p_matrix = matrix(NA, length(unique(toe_data$year)), length(unique(toe_data$year)))

# Rename rows/columns
rownames(p_matrix) = (unique(toe_data$year))
colnames(p_matrix) = (unique(toe_data$year))

# Loop through rows
for(i in 1:nrow(p_matrix)){
  
  # Gather test years
  y1 = as.numeric(rownames(p_matrix)[i])
  y2 = (unique(toe_data$year))
    
    # Loop through columns
    for(j in 1:ncol(p_matrix)){
      
      if(i == j){
        
        p_matrix[i,j] = NA
        
      } else {
      
        # filter data
        df = filter(toe_data, year %in% c(y1, y2[j]))
        
        # Run significance test
        st = survdiff(Surv(quants_ma) ~ year, data = df)
        
        # extract p value
        p = st$pvalue
        
        # Enter into matrix
        p_matrix[i,j] = p
      
      } # End i==j else
        
    } # End matrix column loop
  
} # End matrix row loop



# Reformat to long format
p_matrix_plot = reshape2::melt(p_matrix)

# Plot p-value matrix
ggplot(p_matrix_plot, aes(x = Var1, y = Var2)) +
  geom_raster(aes(fill = value)) +
  scale_fill_gradient(low = "blue", high = "red") + 
  theme_bw() + labs(x = 'Year', y = 'Year', fill = 'p-value')

# Create significance column
p_matrix_plot$sig = ifelse(p_matrix_plot$value <= 0.05, 1, 0)

# Plot p-value matrix
ggplot(p_matrix_plot, aes(x = Var1, y = Var2)) +
  geom_raster(aes(fill = sig)) +
  scale_fill_gradient(low = "red", high = "blue") + 
  theme_bw() + labs(x = 'Year', y = 'Year', fill = 'significance')

# Emergence storage object
emergence = NULL
e_y = 5

# Gather emergence
for(i in 1:length(unique(toe_data$year))){
  
  # Filter year, Filter out years below the data year
  em_data = filter(p_matrix_plot, Var1 == unique(toe_data$year)[i]) %>%
    filter(Var1 <= Var2)
  
  # Change NA to 0
  em_data$sig[which(is.na(em_data$sig))] = 0
  
  # Calculate first year that starts run of e_y years that are significantly different
  em_calc = rle(em_data$sig) %>%
    unclass() %>%
    as.data.frame() %>%
    mutate(end = cumsum(lengths),
           start = c(1, dplyr::lag(end)[-1] + 1)) %>%
    filter((lengths >= e_y) & (values == 1)) %>%
    filter(start == min(start))
  
  # Create add object
  if(nrow(em_calc) > 0){add = (em_data[em_calc$start, 'Var2'])} else {add = NA}
  
  # Grab emergence
  emergence = c(emergence, add)

} # End emergence loop

# Emergence box
em_df = data.frame(year = unique(toe_data$year), emergence)

# Plot year of emergence
ggplot(em_df, aes(x = year, y = emergence)) + geom_line()

# Add own year
em_df_line = rbind(em_df, data.frame(year = unique(toe_data$year), emergence = unique(toe_data$year)))

# Plot p-value matrix
ggplot(p_matrix_plot, aes(x = Var1, y = Var2)) +
  geom_raster(aes(fill = sig), alpha = 0.75) +
  scale_fill_gradient(low = "red", high = "blue") + 
  theme_bw() + labs(x = 'Year', y = 'Year', fill = 'significance') +
  geom_line(data = em_df_line, aes(y = year, x = emergence, group = year), linewidth = 2) + 
  geom_point(data = em_df, aes(y = year, x = emergence), size = 3)

# Reduce to one side
p_matrix_plot_red = filter(p_matrix_plot, Var1 >= Var2)

# Plot p-value matrix
ggplot(p_matrix_plot_red, aes(x = Var1, y = Var2)) +
  # geom_raster(aes(fill = (sig)), alpha = 0.75) +
  geom_raster(aes(fill = as.factor(sig))) +
  scale_fill_brewer(palette = 'Set1', na.value = 'gray') + 
  theme_bw() + labs(x = 'Year', y = 'Year', fill = 'significance') +
  geom_line(data = em_df_line, aes(y = year, x = emergence, group = year), linewidth = 2) + 
  geom_point(data = em_df, aes(y = year, x = emergence), size = 3) +
  geom_text(data = filter(p_matrix_plot_red, Var1 == Var2), aes(x = Var1, y = Var2, label = Var2),
             nudge_y = 2, angle = 270)


ggplot(em_df, aes(x = year, y = emergence)) + geom_line()



# # Interactive plot
# test = ggplot(toe_data, aes(x = quants_ma, y = surv_y, color = year, group = year)) + geom_line() +
#   scale_colour_gradient(low = "blue", high = "red")



# Generate data
# Good example species
# Shift = Vanessa atalanta
# Decouple = Coenonympha pamphilus
# Combination = Thymelicus lineola
data_quants = data_prep(bfly_f, 'Vanessa atalanta')



# Plot survival curves
ggplot(data_quants, aes(x = quantile, y = surv_y)) + geom_line(linewidth = 2) +
  labs(color = 'Year', x = 'Day of Year', y = 'Proportion to Arrive') +
  geom_point(data = data_quants, aes(x = FIRSTDAY, y = abs(surv_y), color = as.factor(YEAR)), size = 3) + 
  theme(text = element_text(size = 20))

# Plot survival residuals over time
ggplot(data_quants, (aes(y = residual, x = YEAR))) + geom_point() +
  labs(x = 'Year', y = 'Difference from Mean Arrival')  + 
  stat_smooth(method = 'glm') + 
  theme(text = element_text(size = 20))

# Plot arrival over time
ggplot(data_quants, (aes(y = FIRSTDAY, x = YEAR))) + geom_point() +
  labs(x = 'Year', y = 'Day of Arrival')  + 
  stat_smooth(method = 'glm') + 
  theme(text = element_text(size = 20))

# Plot temerature over time
ggplot(data_quants, (aes(y = Temp, x = YEAR))) + geom_point() +
  labs(x = 'Year', y = 'Temperature at Arrival')  + 
  stat_smooth(method = 'glm') + 
  theme(text = element_text(size = 20))

# Plot temerature at arrival
ggplot(data_quants, (aes(y = Temp, x = FIRSTDAY, color = YEAR))) + geom_point() +
  labs(x = 'Day of Arrival', y = 'Temperature at Arrival')  + 
  stat_smooth(method = 'glm') + 
  theme(text = element_text(size = 20))


# Run linear models, extract slopes
e_slope = summary(lm(data = data_quants, FIRSTDAY ~ YEAR))$coefficients['YEAR',]
t_slope = summary(lm(data = data_quants, Temp ~ YEAR))$coefficients['YEAR',]





# Slope storage object
slopes = data.frame(species = NULL, event = NULL, temp = NULL)

# Slope loop
for(i in 1:length(unique(bfly_f$SPECIES_NAME))){
  
  # Prep data
  data_s = data_prep(bfly_f, unique(bfly_f$SPECIES_NAME)[i])
  
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
  add = data.frame(species = unique(bfly_f$SPECIES_NAME)[i], 
                   event = e_slope, temp = t_slope, event_p = e_p, temp_p = t_p)
  
  # Add to data frame
  slopes = rbind(slopes, add)
  
} # End slope loop

# Calculate difference in slopes
slopes$diff = slopes$event+slopes$temp
slopes$diff_n = slopes$event+temp_lm$coefficients['YEAR']

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

# Reformat to long format
slopes_l = pivot_longer(slopes, cols = colnames(slopes)[c(-1,-length(colnames(slopes)))])

# plot slope difference
ggplot(slopes, aes(x = species, y = diff, fill = class)) + geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(x = 'Species', y = 'Slope Difference')

# plot slope difference from overall temp
ggplot(slopes, aes(x = species, y = diff_n, fill = class)) + geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(x = 'Species', y = 'Slope Difference')

# plot event slopes
ggplot(slopes, aes(x = species, y = event, fill = event_sig)) + geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(x = 'Species', y = 'Arrival Slope')

# plot temperature slopes
ggplot(slopes, aes(x = species, y = temp, fill = temp_sig)) + geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(x = 'Species', y = 'Temperature at Arrival Slope')


# plot temperature slopes
ggplot(filter(slopes_l, name %in% c('event', 'temp')), aes(x = species, y = value, fill = name)) + 
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(x = 'Species', y = 'Slope')


# plot slopes
test = ggplot(slopes, aes(x = event, y = temp, color = class, text = species)) +
  theme_classic() + geom_hline(yintercept = 0, linetype = 'dashed') + 
  geom_vline(xintercept = 0, linetype = 'dashed') +  
  geom_abline(slope = 1 , intercept = 0) +
  geom_point(size = 2)
ggplotly(test)

# Plot all linear models - events
lm_events = ggplot(bfly_f, aes(x = YEAR, y = FIRSTDAY, color = SPECIES_NAME)) +
  stat_smooth(method = 'lm', se = F)
ggplotly(lm_events)

# Plot all linear models - temperature
lm_temps = ggplot(bfly_f, aes(x = YEAR, y = Temp, color = SPECIES_NAME)) +
  stat_smooth(method = 'lm', se = F)
ggplotly(lm_temps)













