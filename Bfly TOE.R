
# UK butterfly phenology Time of Emergence
# Reid Steele 06/10/2024


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

# Run intialization
source('Bfly Initial.R')

# Prepped data for testing
data = data_prep(bfly_f, 'Vanessa atalanta')
data2 = data_prep(bfly_f, 'Coenonympha pamphilus')
data3 = data_prep(bfly_f, 'Thymelicus lineola')

########################################################################################

# # ToE Method 1: Survival Curve Differences
# 
# # Moving average quantiles
# 
# # Setup
# ma = 3 # Moving average length
# toe_data = NULL# output object
# 
# # Loop through years
# for(i in ma:length(unique(data$YEAR))){
#   
#   # Select years
#   years_ma = seq(unique(data$YEAR)[i]-(ma-1), unique(data$YEAR)[i], 1)
#   
#   # Filter data
#   data_f = filter(data, YEAR %in% years_ma)
#   
#   # Run survival fit
#   sfit = survfit(Surv(FIRSTDAY) ~ YEAR, data = data_f)
#   
#   # Calculate quantiles
#   quants = quantile(sfit, abs(1-data_f$surv_y))$quantile
#   
#   # Average quantiles
#   quants_ma = apply(quants, 2, mean)
#   
#   # Output data
#   out = cbind(year = max(years_ma), quants_ma, surv_y = data_f$surv_y)
#   rownames(out) = NULL
#   
#   # bind
#   toe_data = rbind(toe_data, out)
#   
# } # End moving average loop
# 
# # Set data frame
# toe_data = as.data.frame(toe_data); rownames(toe_data) = NULL
# 
# # Sort
# toe_data = arrange(toe_data, year, quants_ma)
# 
# # Plot moving average lines
# ggplot(toe_data, aes(x = quants_ma, y = surv_y, color = year, group = year)) + geom_line() +
#   scale_colour_gradient(low = "blue", high = "red")
# 
# # Run hypothesis test
# surv_diff <- survdiff(Surv(quants_ma) ~ year, data = toe_data_f)
# surv_diff
# 
# ggplot(data, aes(x = FIRSTDAY, y = surv_y, color = as.factor(YEAR), group = YEAR)) + geom_line()
# 
# # No moving average
# 
# # Significance matrix - visualize TOE ideas
# p_matrix = matrix(NA, length(unique(toe_data$year)), length(unique(toe_data$year)))
# 
# # Rename rows/columns
# rownames(p_matrix) = (unique(toe_data$year))
# colnames(p_matrix) = (unique(toe_data$year))
# 
# # Loop through rows
# for(i in 1:nrow(p_matrix)){
#   
#   # Gather test years
#   y1 = as.numeric(rownames(p_matrix)[i])
#   y2 = (unique(toe_data$year))
#     
#     # Loop through columns
#     for(j in 1:ncol(p_matrix)){
#       
#       if(i == j){
#         
#         p_matrix[i,j] = NA
#         
#       } else {
#       
#         # filter data
#         df = filter(toe_data, year %in% c(y1, y2[j]))
#         
#         # Run significance test
#         st = survdiff(Surv(quants_ma) ~ year, data = df)
#         
#         # extract p value
#         p = st$pvalue
#         
#         # Enter into matrix
#         p_matrix[i,j] = p
#       
#       } # End i==j else
#         
#     } # End matrix column loop
#   
# } # End matrix row loop
# 
# 
# 
# # Reformat to long format
# p_matrix_plot = reshape2::melt(p_matrix)
# 
# # Plot p-value matrix
# ggplot(p_matrix_plot, aes(x = Var1, y = Var2)) +
#   geom_raster(aes(fill = value)) +
#   scale_fill_gradient(low = "blue", high = "red") + 
#   theme_bw() + labs(x = 'Year', y = 'Year', fill = 'p-value')
# 
# # Create significance column
# p_matrix_plot$sig = ifelse(p_matrix_plot$value <= 0.05, 1, 0)
# 
# # Plot p-value matrix
# ggplot(p_matrix_plot, aes(x = Var1, y = Var2)) +
#   geom_raster(aes(fill = sig)) +
#   scale_fill_gradient(low = "red", high = "blue") + 
#   theme_bw() + labs(x = 'Year', y = 'Year', fill = 'significance')
# 
# # Emergence storage object
# emergence = NULL
# e_y = 5
# 
# # Gather emergence
# for(i in 1:length(unique(toe_data$year))){
#   
#   # Filter year, Filter out years below the data year
#   em_data = filter(p_matrix_plot, Var1 == unique(toe_data$year)[i]) %>%
#     filter(Var1 <= Var2)
#   
#   # Change NA to 0
#   em_data$sig[which(is.na(em_data$sig))] = 0
#   
#   # Calculate first year that starts run of e_y years that are significantly different
#   em_calc = rle(em_data$sig) %>%
#     unclass() %>%
#     as.data.frame() %>%
#     mutate(end = cumsum(lengths),
#            start = c(1, dplyr::lag(end)[-1] + 1)) %>%
#     filter((lengths >= e_y) & (values == 1)) %>%
#     filter(start == min(start))
#   
#   # Create add object
#   if(nrow(em_calc) > 0){add = (em_data[em_calc$start, 'Var2'])} else {add = NA}
#   
#   # Grab emergence
#   emergence = c(emergence, add)
# 
# } # End emergence loop
# 
# # Emergence box
# em_df = data.frame(year = unique(toe_data$year), emergence)
# 
# # Plot year of emergence
# ggplot(em_df, aes(x = year, y = emergence)) + geom_line()
# 
# # Add own year
# em_df_line = rbind(em_df, data.frame(year = unique(toe_data$year), emergence = unique(toe_data$year)))
# 
# # Plot p-value matrix
# ggplot(p_matrix_plot, aes(x = Var1, y = Var2)) +
#   geom_raster(aes(fill = sig), alpha = 0.75) +
#   scale_fill_gradient(low = "red", high = "blue") + 
#   theme_bw() + labs(x = 'Year', y = 'Year', fill = 'significance') +
#   geom_line(data = em_df_line, aes(y = year, x = emergence, group = year), linewidth = 2) + 
#   geom_point(data = em_df, aes(y = year, x = emergence), size = 3)
# 
# # Reduce to one side
# p_matrix_plot_red = filter(p_matrix_plot, Var1 >= Var2)
# 
# # Plot p-value matrix
# ggplot(p_matrix_plot_red, aes(x = Var1, y = Var2)) +
#   # geom_raster(aes(fill = (sig)), alpha = 0.75) +
#   geom_raster(aes(fill = as.factor(sig))) +
#   scale_fill_brewer(palette = 'Set2', na.value = 'gray') + 
#   theme_bw() + labs(x = 'Year', y = 'Year', fill = 'significance') +
#   geom_line(data = em_df_line, aes(y = year, x = emergence, group = year), linewidth = 2) + 
#   geom_point(data = em_df, aes(y = year, x = emergence), size = 3) +
#   geom_text(data = filter(p_matrix_plot_red, Var1 == Var2), aes(x = Var1, y = Var2, label = Var2),
#              nudge_y = 2, angle = 270)
# 
# # Plot emergence year over time
# ggplot(em_df, aes(x = year, y = emergence)) + geom_line()


########################################################################################

# ToE Method 2: KS Test on Detrended Data

# Calculate time of emergence using KS Test
ks_emergence = function(data,    # Input data
                        emt = 5, # Emergence threshold (number of years for emergence)
                        plot = T # Plot results?
                        ){

  # Fit linear model (arrival)
  lm_ev = lm(FIRSTDAY ~ YEAR, data = data)
  
  # Pull out slope
  lm_ev_summ = summary(lm_ev)
  lm_ev_b = lm_ev_summ$coefficients['YEAR','Estimate']
  
  # # Plot LM
  # plot(FIRSTDAY ~ YEAR, data = data, pch = 16)
  # abline(lm_ev, col = 'blue', lwd = 2)
  
  # Gather years
  years = unique(data$YEAR)
  
  # Adjust year to order
  years_0 = years-min(years)
  
  # Generate match index
  ind = match(data$YEAR, years)
  
  # Calculate adjustment
  adj = ind*lm_ev_b
  
  # Detrend
  data_d = data
  data_d$FIRSTDAY = data_d$FIRSTDAY-adj
  
  # Fit linear model (arrival)
  lm_ev_d = lm(FIRSTDAY ~ YEAR, data = data_d)
  
  # Pull out slope
  lm_ev_summ_d = summary(lm_ev_d)
  lm_ev_b_d = lm_ev_summ_d$coefficients['YEAR','Estimate']
  
  # # Plot LM
  # plot(FIRSTDAY ~ YEAR, data = data_d, pch = 16)
  # abline(lm_ev_d, col = 'blue', lwd = 2)
  
  # p value container
  ks_p = rep(NA, length(unique(data$YEAR)))
  
  # Cycle years and perform KS test
  for(i in 1:length(unique(data$YEAR))){
    
    # Run KS test
    ks = ks.test(data_d[data_d$YEAR==unique(data$YEAR)[i], 'FIRSTDAY'], 
                 data[data$YEAR==unique(data$YEAR)[i], 'FIRSTDAY'], 
                 alternative = 'less')
    
    # Fill ks_p
    ks_p[i] = ks$p.value
    
  } # End KS test loop
  
  # Calculate emergence
  emerged = frollapply(ks_p, n = emt, function(x){all(x<0.05)}, align = 'left')
  
  # Set all to 1 after emergence
  
  # Loop through emereged
  for(i in 1:length(emerged)){
    
    # If emerged is NA, leave it
    if(is.na(emerged[i])){
      
      emerged[i] = NA
      
    } else {
      
      # If previous value of emerged is 1, make all subsequent values 1
      if((is.na(shift(emerged)[i]) == F) & (shift(emerged)[i]==1)){
        
        emerged[i] = 1
        
      } # End second if
      
    } # End else
    
  } # End loop
  
  # Create data frame
  ks_p_df = cbind.data.frame(year = unique(data$YEAR), p = ks_p, emerged)
  
  # Plot if requested
  if(plot == T){
  
    # Plot p-value curve
    plot(ks_p ~ unique(data$YEAR), type = 'l', pch = 16, lwd = 2,
         main = unique(data$SPECIES_NAME),
         ylab = 'KS test p value', xlab = 'Year')
    abline(h = 0.05, lty = 'dashed', col = 'blue') # Add significance threshold line
    abline(v = ks_p_df[min(which(emerged == 1)),'year'], lwd = 2, col = 'red') # Add Emerged Year
    text(y = 0.5, x = ks_p_df[min(which(emerged == 1)),'year'] + 4, label = paste('TOE:', ks_p_df[min(which(emerged == 1)),'year']), col = 'red')

  } # End plot if
  
  # Return data frame
  return(ks_p_df)
  
} # end KS emergence function

# Test time of emergence function
ks_emergence(data)
ks_emergence(data2)
ks_emergence(data3)









# Calculate time of decoupling using KS Test
ks_decouple = function(data,    # Input data
                       emt = 5, # Emergence threshold (number of years for emergence)
                       plot = T # Plot results?
){
  
  # Fit linear model (arrival)
  lm_ev = lm(Temp ~ YEAR, data = data)
  
  # Pull out slope
  lm_ev_summ = summary(lm_ev)
  lm_ev_b = lm_ev_summ$coefficients['YEAR','Estimate']
  
  # # Plot LM
  # plot(Temp ~ YEAR, data = data, pch = 16)
  # abline(lm_ev, col = 'blue', lwd = 2)
  
  # Gather years
  years = unique(data$YEAR)
  
  # Adjust year to order
  years_0 = years-min(years)
  
  # Generate match index
  ind = match(data$YEAR, years)
  
  # Calculate adjustment
  adj = ind*lm_ev_b
  
  # Detrend
  data_d = data
  data_d$Temp = data_d$Temp-adj
  
  # Fit linear model (arrival)
  lm_ev_d = lm(Temp ~ YEAR, data = data_d)
  
  # Pull out slope
  lm_ev_summ_d = summary(lm_ev_d)
  lm_ev_b_d = lm_ev_summ_d$coefficients['YEAR','Estimate']
  
  # # Plot LM
  # plot(Temp ~ YEAR, data = data_d, pch = 16)
  # abline(lm_ev_d, col = 'blue', lwd = 2)
  
  # p value container
  ks_p = rep(NA, length(unique(data$YEAR)))
  
  # Cycle years and perform KS test
  for(i in 1:length(unique(data$YEAR))){
    
    # Run KS test
    ks = ks.test(data_d[data_d$YEAR==unique(data$YEAR)[i], 'Temp'], 
                 data[data$YEAR==unique(data$YEAR)[i], 'Temp'], 
                 alternative = 'greater')
    
    # Fill ks_p
    ks_p[i] = ks$p.value
    
  } # End KS test loop
  
  # Calculate emergence
  emerged = frollapply(ks_p, n = emt, function(x){all(x<0.05)}, align = 'left')
  
  # Set all to 1 after emergence
  
  # Loop through emereged
  for(i in 1:length(emerged)){
    
    # If emerged is NA, leave it
    if(is.na(emerged[i])){
      
      emerged[i] = NA
      
    } else {
      
      # If previous value of emerged is 1, make all subsequent values 1
      if((is.na(shift(emerged)[i]) == F) & (shift(emerged)[i]==1)){
        
        emerged[i] = 1
        
      } # End second if
      
    } # End else
    
  } # End loop
  
  # Create data frame
  ks_p_df = cbind.data.frame(year = unique(data$YEAR), p = ks_p, emerged)
  
  # Plot if requested
  if(plot == T){
    
    # Plot p-value curve
    plot(ks_p ~ unique(data$YEAR), type = 'l', pch = 16, lwd = 2,
         main = unique(data$SPECIES_NAME),
         ylab = 'KS test p value', xlab = 'Year')
    abline(h = 0.05, lty = 'dashed', col = 'blue') # Add significance threshold line
    abline(v = ks_p_df[min(which(emerged == 1)),'year'], lwd = 2, col = 'red') # Add Emerged Year
    text(y = 0.5, x = ks_p_df[min(which(emerged == 1)),'year'] + 4, label = paste('TOD:', ks_p_df[min(which(emerged == 1)),'year']), col = 'red')
    
  } # End plot if
  
  # Return data frame
  return(ks_p_df)
  
} # end KS emergence function

# Test time of decoupling function
ks_decouple(data)
ks_decouple(data2)
ks_decouple(data3)



# Calculate time of mismatch using KS Test
ks_mismatch = function(sp1,     # Input data for species 1
                       sp2,     # Input data for species 2
                       alt = 'two.sided',     # Alternative hypothesis for ks test: greater, less, 2 sided
                       emt = 5, # Emergence threshold (number of years for emergence)
                       plot = T # Plot results?
){
  
  # Standardize
  sp1$FIRSTDAY = scale(sp1$FIRSTDAY)
  sp2$FIRSTDAY = scale(sp2$FIRSTDAY)
  
  # Filter to common sites
  sp1 = filter(sp1, SITENAME %in% sp2$SITENAME)
  sp2 = filter(sp2, SITENAME %in% sp1$SITENAME)
  
  # Fit linear model (arrival)
  lm_ev = lm(FIRSTDAY ~ YEAR, data = sp1)
  
  # Pull out slope
  lm_ev_summ = summary(lm_ev)
  lm_ev_b = lm_ev_summ$coefficients['YEAR','Estimate']
  
  # # Plot LM
  # plot(FIRSTDAY ~ YEAR, data = sp1, pch = 16)
  # abline(lm_ev, col = 'blue', lwd = 2)
  
  # Fit linear model (arrival)
  lm_ev_d = lm(FIRSTDAY ~ YEAR, data = sp2)
  
  # Pull out slope
  lm_ev_summ_d = summary(lm_ev_d)
  lm_ev_b_d = lm_ev_summ_d$coefficients['YEAR','Estimate']
  
  # # Plot LM
  # plot(FIRSTDAY ~ YEAR, data = sp2, pch = 16)
  # abline(lm_ev_d, col = 'blue', lwd = 2)
  
  # p value container
  ks_p = rep(NA, length(unique(sp1$YEAR)))
  
  # Cycle years and perform KS test
  for(i in 1:length(unique(sp1$YEAR))){
    
    # Run KS test
    ks = ks.test(sp2[sp2$YEAR==unique(sp2$YEAR)[i], 'FIRSTDAY'], 
                 sp1[sp1$YEAR==unique(sp1$YEAR)[i], 'FIRSTDAY'], 
                 alternative = alt)
    
    # Fill ks_p
    ks_p[i] = ks$p.value
    
  } # End KS test loop
  
  # Calculate emergence
  emerged = frollapply(ks_p, n = emt, function(x){all(x<0.05)}, align = 'left')
  
  # Set all to 1 after emergence
  
  # Loop through emereged
  for(i in 1:length(emerged)){
    
    # If emerged is NA, leave it
    if(is.na(emerged[i])){
      
      emerged[i] = NA
      
    } else {
      
      # If previous value of emerged is 1, make all subsequent values 1
      if((is.na(shift(emerged)[i]) == F) & (shift(emerged)[i]==1)){
        
        emerged[i] = 1
        
      } # End second if
      
    } # End else
    
  } # End loop
  
  # Create data frame
  ks_p_df = cbind.data.frame(year = unique(sp1$YEAR), p = ks_p, emerged)
  
  # Plot if requested
  if(plot == T){
    
    # Plot p-value curve
    plot(ks_p ~ unique(sp1$YEAR), type = 'l', pch = 16, lwd = 2,
         main = paste(unique(sp1$SPECIES_NAME), 'x', unique(sp2$SPECIES_NAME)),
         ylab = 'KS test p value', xlab = 'Year')
    abline(h = 0.05, lty = 'dashed', col = 'blue') # Add significance threshold line
    abline(v = ks_p_df[min(which(emerged == 1)),'year'], lwd = 2, col = 'red') # Add Emerged Year
    text(y = 0.5, x = ks_p_df[min(which(emerged == 1)),'year'] + 4, label = paste('TOM:', ks_p_df[min(which(emerged == 1)),'year']), col = 'red')
    
  } # End plot if
  
  # Return data frame
  return(ks_p_df)
  
} # end KS emergence function

# Test ks_mismatch
ks_mismatch(data, data2, alt = 'two.sided') # Check if species 2 difference is less than species 1
ks_mismatch(data, data3, alt = 'greater') # Check if species 2 difference is less than species 1
ks_mismatch(data, data3, alt = 'less') # Check if species 2 difference is less than species 1


# Emergence - species curves
em_curve = NULL
dc_curve = NULL

# Loop through species
for(i in 1:length(unique(bfly_f$SPECIES_NAME))){
  
  # Load in data
  sp = data_prep(bfly_f, unique(bfly_f$SPECIES_NAME)[i])
  
  # Run curves
  em_curve = rbind(em_curve, cbind(species = unique(bfly_f$SPECIES_NAME)[i], ks_emergence(sp, plot = F)))
  dc_curve = rbind(dc_curve, cbind(species = unique(bfly_f$SPECIES_NAME)[i], ks_decouple(sp, plot = F)))
  
} # End species loop

# Rename emerged coulumn to decoupled in dc_curve column
colnames(dc_curve)[which(colnames(dc_curve) =='emerged')] = 'decoupled'

# Join data frames
all_curves = left_join(em_curve, dc_curve, by = c('species', 'year'))

# Add combination curve
all_curves = all_curves %>% mutate(combination = ifelse((emerged == 1) & (decoupled == 1), 1, 0)) %>%
  mutate(shift = ifelse(combination == 1, 0, emerged), # Remove combined shifts
         decouple = ifelse(combination == 1, 0, decoupled), # Remove combined decouples
         unaffected = ifelse(shift+decouple+combination == 0, 1, 0)) # Add no climate change column

# Calculate final metrics
conc_curve = all_curves %>% group_by(year) %>% # Group by year
  summarize(n_em = sum(emerged), # Number of emerged (shift) species
            n_dc = sum(decoupled), # Number of decoupled species
            n_cb = sum(combination), # Number of combined species
            n_shift = sum(shift), # Number of shifted, non combined species
            n_decouple = sum(decouple), # Number of decoupled, non combined species
            n_ncc = sum(unaffected), # Number of species with no climate change signal
            sp = length(unique(species))) # Total number of species

# Remove NA years
conc_curve = conc_curve[is.na(conc_curve$n_em) == F,]

# Plot
plot(n_ncc ~ year, data = conc_curve, col = 'black', type = 'l', lwd = 2,
     xlab = 'Year', ylab = 'Number of Species Emerged', ylim = c(0, 60))
lines(n_em ~ year, data = conc_curve, col = 'chartreuse2', type = 'l', lwd = 2, lty = 'dashed')
lines(n_dc ~ year, data = conc_curve, col = 'firebrick1', type = 'l', lwd = 2, lty = 'dashed')
lines(n_shift ~ year, data = conc_curve, col = 'chartreuse2', type = 'l', lwd = 2)
lines(n_decouple ~ year, data = conc_curve, col = 'firebrick1', type = 'l', lwd = 2)
lines(n_cb ~ year, data = conc_curve, col = 'darkorange', type = 'l', lwd = 2)
legend('topright', col = c('chartreuse2', 'firebrick1', 'darkorange'), lwd = 2,
       legend = c('Shifted', 'Decoupled', 'Combination'), bty = 'n', lty = 'solid')
legend('top', lty = c('solid', 'dashed'), legend = c('Classification','Effect'), lwd = 2, bty = 'n')

# Try ggplot
conc_curve_l = select(conc_curve, year, n_shift, n_decouple, n_cb, n_ncc)

# Rename columns
colnames(conc_curve_l) = c('year', 'shift', 'decouple', 'combination', 'no signal')

# Pivot longer
conc_curve_l = pivot_longer(conc_curve_l, cols = c('shift', 'decouple', 'combination', 'no signal'), 
                            names_to = 'classification', values_to = 'n_species')

# Plot curve
ggplot(conc_curve_l, aes(x = year, y = n_species, color = classification)) +
  geom_line(linewidth = 1) + theme_classic() + labs(x = 'Year', y = 'Number of Species')

# Matrix plot
curve_matrix = matrix(NA, length(unique(all_curves$species)), length(unique(all_curves$year)))
rownames(curve_matrix) = unique(all_curves$species)
colnames(curve_matrix) = unique(all_curves$year)

# Condense all curves
all_curves$class = all_curves$unaffected
all_curves$class = ifelse(all_curves$shift == 1, 2, all_curves$class)
all_curves$class = ifelse(all_curves$decouple == 1, 3, all_curves$class)
all_curves$class = ifelse(all_curves$combination == 1, 4, all_curves$class)

# Reformat to matrix
curve_matrix = all_curves %>% select(species, year, class)

# Plot p-value matrix
ggplot(curve_matrix, aes(x = year, y = species)) +
  # geom_raster(aes(fill = (sig)), alpha = 0.75) +
  geom_raster(aes(fill = as.factor(class))) +
  scale_fill_brewer(palette = 'Set2', na.value = 'gray', labels = c('No Signal', 'Shift', 'Decouple', 'Combination')) +
  theme_bw() + labs(x = 'Species', y = 'Year', fill = 'Classification')
  
  
  
  
  
  
  