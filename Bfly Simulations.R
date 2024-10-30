# UK butterfly simulations
# Reid Steele 06/10/2024

# Libraries
library(tidyverse)
library(survival)

# Working Directory
setwd('C:/Users/R3686/OneDrive - Dalhousie University/Documents/Acoustic Tracking/Data/Butterflies')

# Read in data
bfly = read.csv('./data/ukbmsphenology2022.csv')

# Filter species of interest
ata = filter(bfly, SPECIES_NAME == 'Vanessa atalanta')
napi = filter(bfly, SPECIES_NAME == 'Pieris napi')

# Filter brood 1 for napi
napi = filter(napi, BROOD == 1)

# Remove years with one record
ata = group_by(ata, YEAR)
ata = ata %>% filter(!(YEAR %in% group_keys(ata)[[1]][which(group_size(ata) == 1)]))
napi = group_by(napi, YEAR)
napi = napi %>% filter(!(YEAR %in% group_keys(napi)[[1]][which(group_size(napi) == 1)]))

# Arrival date survival analysis
surv_ata_a = survfit(formula = Surv(FIRSTDAY) ~ YEAR, data = ata)
autoplot(surv_ata_a)

surv_napi_a = survfit(formula = Surv(FIRSTDAY) ~ YEAR, data = napi)
autoplot(surv_napi_a)






# Calculate quantiles - atalanta
surv_ata_q = quantile(surv_ata_a, seq(0,1,0.01))

# Reformat quantiles to data frame
surv_ata_quants = data.frame(rownames(surv_ata_q$quantile), surv_ata_q$quantile); colnames(surv_ata_quants) = c('year', seq(0,100,1))
surv_ata_quants = pivot_longer(surv_ata_quants, cols = -1, names_to = 'surv', values_to = 'first_yday')
surv_ata_quants$year = gsub('YEAR=', '', as.character(surv_ata_quants$year))
surv_ata_quants$surv = as.numeric(surv_ata_quants$surv); surv_ata_quants$year = as.numeric(surv_ata_quants$year) 

# Calculate mean, sd, upper, and lower appearance - filter first 10 years
surv_ata_qavg = group_by(surv_ata_quants, surv) %>% filter(year <= (min(year) + 10)) %>%
  summarize(first_yday_mean = mean(first_yday), sd = sd(first_yday), upper = first_yday_mean + 2*sd, lower = first_yday_mean - 2*sd)

# # Plot quantiles - all years
# ggplot(surv_ata_quants, aes(x = first_yday, y = surv, color = year)) + geom_line()

# Plot quantiles - historical years
ggplot(filter(surv_ata_quants, year <= (min(year) + 10)), aes(x = first_yday, y = surv, color = as.factor(year))) + geom_line()

# Plot quantiles with historical average (first 10 years)
ggplot(surv_ata_qavg, aes(x = first_yday_mean, y = surv)) +
  geom_line(data = surv_ata_quants, aes(x = first_yday, y = surv, color = as.factor(year))) + geom_line(linewidth = 2) +
  geom_line(data = surv_ata_qavg, aes(x = lower, y = surv), linetype = 'dashed') + 
  geom_line(data = surv_ata_qavg, aes(x = upper, y = surv), linetype = 'dashed')

# Gather 50% quantile
surv_ata_50 = filter(surv_ata_quants, surv == 50)

# Identify 50% quantiles outside 95% CI
surv_ata_50$sig = ifelse(surv_ata_50$first_yday < as.numeric(surv_ata_qavg[surv_ata_qavg$surv == 50, 'lower']), 1, 0)

which(surv_ata_50$sig == 1)






# Calculate quantiles - napi
surv_napi_q = quantile(surv_napi_a, seq(0,1,0.01))

# Reformat quantiles to data frame
surv_napi_quants = data.frame(rownames(surv_napi_q$quantile), surv_napi_q$quantile); colnames(surv_napi_quants) = c('year', seq(0,100,1))
surv_napi_quants = pivot_longer(surv_napi_quants, cols = -1, names_to = 'surv', values_to = 'first_yday')
surv_napi_quants$year = gsub('YEAR=', '', as.character(surv_napi_quants$year))
surv_napi_quants$surv = as.numeric(surv_napi_quants$surv); surv_napi_quants$year = as.numeric(surv_napi_quants$year) 

# Calculate mean, sd, upper, and lower appearance - filter first 10 years
surv_napi_qavg = group_by(surv_napi_quants, surv) %>% filter(year <= (min(year) + 10)) %>%
  summarize(first_yday_mean = mean(first_yday), sd = sd(first_yday), upper = first_yday_mean + 2*sd, lower = first_yday_mean - 2*sd)

# Plot quantiles - all years
ggplot(surv_napi_quants, aes(x = first_yday, y = surv, color = as.factor(year))) + geom_line()

# Plot quantiles - historical years
ggplot(filter(surv_napi_quants, year <= (min(year) + 10)), aes(x = first_yday, y = surv, color = as.factor(year))) + geom_line()

# Plot quantiles with historical average (first 10 years)
ggplot(surv_napi_qavg, aes(x = first_yday_mean, y = surv)) +
  geom_line(data = surv_napi_quants, aes(x = first_yday, y = surv, color = as.factor(year))) + geom_line(linewidth = 2) +
  geom_line(data = surv_napi_qavg, aes(x = lower, y = surv), linetype = 'dashed') + 
  geom_line(data = surv_napi_qavg, aes(x = upper, y = surv), linetype = 'dashed')

# Gather 50% quantile
surv_napi_50 = filter(surv_napi_quants, surv == 50)

# Identify 50% quantiles outside 95% CI
surv_napi_50$sig = ifelse(surv_napi_50$first_yday < as.numeric(surv_napi_qavg[surv_napi_qavg$surv == 50, 'lower']), 1, 0)

which(surv_napi_50$sig == 1)




data = surv_ata_quants

# Moving average
mal = 25:75 # survival range
years = unique(data$year) # years
win = 3 # moving average window
may = matrix(nrow = length(years), ncol = length(mal), NA)# container for moving average years

# Gather past values
# Loop through data years
for(i in 1:nrow(may)){
  
  # rows = (i-(mal-1)):i
  
  # Loop through moving average years
  for(j in 1:length(mal)){
    
    # # Check if row is NA
    # if(rows[j] <= 0){
    #   
    #   # Set NA if row is negative or 0
    #   may[i,j] = NA
    #   
    # } else {
    #   
    #   # Otherwise set data
    #   may[i,j] = data[rows[j]]
    #   
    # } # End ifelse
    
    may[i,j] = data$first_yday[which((data$year == years[i]) & (data$surv == mal[j]))]
    
  } # end moving average loop
  
} # end data year loop

# Historical average
havg_ata = filter(surv_ata_quants, (year <= (min(year) + 9))) %>%
  filter(surv >= 25) %>% filter(surv <= 75)

pvals = NA
for(i in win:nrow(may)){
  
  pvals[i] = t.test(havg_ata$first_yday, may[(i-win+1):i,], alternate = 'l')$p.value
  
}

#plot(rowMeans(may) ~ surv_ata_50$year)
plot(pvals ~ surv_ata_50$year)







data = surv_napi_quants

# Moving average
mal = 25:75 # survival range
years = unique(data$year) # years
win = 3 # moving average window
may = matrix(nrow = length(years), ncol = length(mal), NA)# container for moving average years

# Gather past values
# Loop through data years
for(i in 1:nrow(may)){
  
  # rows = (i-(mal-1)):i
  
  # Loop through moving average years
  for(j in 1:length(mal)){
    
    # # Check if row is NA
    # if(rows[j] <= 0){
    #   
    #   # Set NA if row is negative or 0
    #   may[i,j] = NA
    #   
    # } else {
    #   
    #   # Otherwise set data
    #   may[i,j] = data[rows[j]]
    #   
    # } # End ifelse
    
    may[i,j] = data$first_yday[which((data$year == years[i]) & (data$surv == mal[j]))]
    
  } # end moving average loop
  
} # end data year loop

# Historical average
havg_napi = filter(surv_napi_quants, (year <= (min(year) + 9))) %>%
  filter(surv >= 25) %>% filter(surv <= 75)

pvals = NA
for(i in win:nrow(may)){
  
  pvals[i] = t.test(havg_napi$first_yday, may[(i-win+1):i,], alternate = 'l')$p.value
  
}

#plot(rowMeans(may) ~ surv_napi_50$year)
plot(pvals ~ surv_napi_50$year)






# Checking probability of going below 95%

# Gather test window quantiles
surv_ata_win = filter(surv_ata_quants, surv %in% mal)

# Summarize historical data
havg_ata_summ = havg_ata %>% group_by(surv) %>% summarize(h_mean = mean(first_yday), h_sd = sd(first_yday),
                                                  h_lower = h_mean - h_sd*2, h_upper = h_mean + h_sd*2)

# Join historical summaries
surv_ata_win = left_join(surv_ata_win, havg_ata_summ) %>%
  # Check if first_yday is below lower bound or above upper bound
  mutate(b_l = first_yday < h_lower, a_u = first_yday > h_upper) 

# Calculate proportion above or below for each year
surv_ata_prob = surv_ata_win %>% group_by(year) %>% summarize(prob_b = sum(b_l)/n(), prob_a = sum(a_u)/n())

# Pivot long
surv_ata_prob_l = pivot_longer(surv_ata_prob, cols = -1)

# plot
ggplot(surv_ata_prob_l, aes(x = year, y = value, color = name)) + geom_point()


# Plot all years middle
ggplot(surv_ata_win) + 
  geom_point(aes(x = first_yday, y = surv, color = b_l)) + 
  geom_line(aes(x = h_lower, y = surv, linetype = 'dashed')) + 
  geom_line(aes(x = h_upper, y = surv, linetype = 'dashed')) + 
  facet_wrap(~year)


# Gather past values
# Loop through data years
for(i in 1:nrow(may)){
  
  # rows = (i-(mal-1)):i
  
  # Loop through moving average years
  for(j in 1:length(mal)){
    
    # # Check if row is NA
    # if(rows[j] <= 0){
    #   
    #   # Set NA if row is negative or 0
    #   may[i,j] = NA
    #   
    # } else {
    #   
    #   # Otherwise set data
    #   may[i,j] = data[rows[j]]
    #   
    # } # End ifelse
    
    may[i,j] = data$first_yday[which((data$year == years[i]) & (data$surv == mal[j]))]
    
  } # end moving average loop
  
} # end data year loop

save.image('Bfly_Initial.RData')




df = NULL
for(i in 1:length(unique(ata$YEAR))){
  
  data = filter(ata, YEAR == unique(ata$YEAR)[i])
  
  s = survfit(formula = Surv(FIRSTDAY) ~ 1, data = data)
  
  df = rbind(df, data.frame(YEAR = unique(ata$YEAR)[i], FIRSTDAY = s$time, surv_y = s$surv))
  
} # end year loop

ata = left_join(ata, df)



df = NULL
for(i in 1:length(unique(napi$YEAR))){
  
  data = filter(napi, YEAR == unique(napi$YEAR)[i])
  
  s = survfit(formula = Surv(FIRSTDAY) ~ 1, data = data)
  
  df = rbind(df, data.frame(YEAR = unique(napi$YEAR)[i], FIRSTDAY = s$time, surv_y = s$surv))
  
} # end year loop

napi = left_join(napi, df)



# Join average quantile to original data
test = quantile(surv_ata_a, abs(1-ata$surv_y))$quantile
rownames(test) = gsub('year=', '', rownames(test))
# test = test[as.numeric(rownames(test)) %in% 2011:2018,]
test2 = colMeans(test)

test3 = cbind(ata, quantile = test2)

ggplot(test3, aes(x = quantile, y = surv_y)) + geom_line() +
  geom_point(data = test3, aes(x = FIRSTDAY, y = abs(surv_y), color = as.factor(YEAR)))

ggplot(test3, (aes(y = quantile, x = FIRSTDAY))) + geom_point()

ggplot(test3, (aes(y = quantile - FIRSTDAY, group = YEAR, x = YEAR))) + geom_boxplot()

ggplot(test3, (aes(y = quantile - FIRSTDAY, x = YEAR, color = surv_y))) + geom_point()

ggplot(test3, (aes(y = quantile - FIRSTDAY, x = YEAR))) + geom_point() + 
  stat_smooth(method = 'glm') + labs(x = 'Year', y = 'Difference from Mean Arrival') + 
  theme(text = element_text(size = 20))

ggplot(test3, (aes(y = quantile - FIRSTDAY, x = YEAR, color = surv_y))) + geom_point()  + 
  stat_summary(fun=mean, aes(group=1), geom="line", colour="blue")


# ggplot(test3, (aes(y = sst, x = YEAR, color = surv_y))) + geom_point() +
#   stat_smooth(method = 'glm')
# 
# ggplot(test3, (aes(y = quantile - FIRSTDAY, x = sst, color = as.factor(YEAR)))) + geom_point()
# 
# ggplot(test3, (aes(y = quantile - FIRSTDAY, x = sst, color = FIRSTDAY))) + geom_point() +
#   stat_smooth(method = 'glm')


q = group_by(ata, YEAR) %>% summarise(arr = mean(FIRSTDAY))
ggplot(q, aes(x = YEAR, y = arr)) + geom_point()
# ggplot(q, aes(x = YEAR, y = sst)) + geom_point()


ggplot(ata, aes(x = YEAR, y = FIRSTDAY)) + geom_point() + stat_smooth(method = 'glm')
ggplot(test3, (aes(y = quantile - FIRSTDAY, x = YEAR, color = surv_y))) + geom_point() + scale_y_reverse()

x11()











# Join average quantile to original data
test = quantile(surv_napi_a, abs(1-napi$surv_y))$quantile
rownames(test) = gsub('year=', '', rownames(test))
# test = test[as.numeric(rownames(test)) %in% 2011:2018,]
test2 = colMeans(test)

test3 = cbind(napi, quantile = test2)

ggplot(test3, aes(x = quantile, y = surv_y)) + geom_line() +
  geom_point(data = test3, aes(x = FIRSTDAY, y = abs(surv_y), color = as.factor(YEAR)))

ggplot(test3, (aes(y = quantile, x = FIRSTDAY))) + geom_point()

ggplot(test3, (aes(y = quantile - FIRSTDAY, group = YEAR, x = YEAR))) + geom_boxplot()

ggplot(test3, (aes(y = quantile - FIRSTDAY, x = YEAR, color = surv_y))) + geom_point()

ggplot(test3, (aes(y = quantile - FIRSTDAY, x = YEAR, color = surv_y))) + geom_point() + 
  stat_smooth(method = 'glm')

# ggplot(test3, (aes(y = sst, x = YEAR, color = surv_y))) + geom_point() +
#   stat_smooth(method = 'glm')
# 
# ggplot(test3, (aes(y = quantile - FIRSTDAY, x = sst, color = as.factor(YEAR)))) + geom_point()
# 
# ggplot(test3, (aes(y = quantile - FIRSTDAY, x = sst, color = FIRSTDAY))) + geom_point() +
#   stat_smooth(method = 'glm')

q = group_by(napi, YEAR) %>% summarise(arr = mean(FIRSTDAY))
ggplot(q, aes(x = YEAR, y = arr)) + geom_point()
# ggplot(q, aes(x = YEAR, y = sst)) + geom_point()


ggplot(napi, aes(x = YEAR, y = FIRSTDAY)) + geom_point() + stat_smooth(method = 'glm')
ggplot(test3, (aes(y = quantile - FIRSTDAY, x = YEAR, color = surv_y))) + geom_point() + scale_y_reverse() + stat_smooth(method = 'glm')




surv_napi_a = survfit(formula = Surv(FIRSTDAY) ~ YEAR, data = napi)
autoplot(surv_napi_a)


