# Butterfly functions and initialization

# Libraries
library(tidyverse)
library(survival)
library(lubridate)
library(plotly)
library(phenometrics)

# source('TOE.R')

# Working Directory
# setwd('C:/Users/R3686/OneDrive - Dalhousie University/Documents/Acoustic Tracking/Data/Butterflies')

# Read in data
bfly = read.csv('./Butterflies/data/ukbmsphenology2022.csv')
temp = read.csv('./Butterflies/data/uktemps.csv')

# Select event
bfly$event = bfly$MEAN_FLIGHT_DATE

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
temp$event_r = ifelse(leap_year(temp$YEAR), temp$yday-91, temp$yday-90)

# Filter temp
temp = temp[,c('Temp', 'YEAR', 'event_r')]

# Calculate mean temperature
temp_mean = group_by(temp, YEAR) %>% summarize(temp = mean(Temp)) %>% mutate(temp_norm = (temp - mean(temp))/sd(temp))

# Run linear model
temp_lm = lm(data = temp_mean, temp_norm ~ YEAR)
summary(temp_lm)

# Plot normalize temperature
ggplot(temp_mean, aes(x = YEAR, y = temp_norm)) + geom_line() + geom_point() + stat_smooth(method = 'lm')

# Round event
bfly$event_r = round(bfly$event)

# Join temperature
bfly = left_join(bfly, temp, by = c('YEAR', 'event_r'))

# Plot butterflies

# Overwintering species
overwinter = c('Aglais io', 'Gonepteryx rhamni')

# # Filter first brood
bfly_0 = filter(bfly, NUMBER_OF_BROODS == 1)
bfly_1 = filter(bfly, BROOD == 1) %>% filter(!(SPECIES_NAME %in% overwinter))
bfly_2 = filter(bfly, SPECIES_NAME %in% overwinter) %>% filter(BROOD == 2)
bfly_f = rbind(bfly_1, bfly_0, bfly_2) %>% filter(SPECIES_NAME != 'Thymelicus lineola/sylvestris') %>%
  filter(!((SPECIES_NAME == 'Phengaris arion')&(YEAR==1987))) %>% # Remove single observation of Phengaris arion in 1987
  filter(FIRSTDAY != LASTDAY) %>% # Remove single observations
  filter(YEAR >= 1979) %>%
  arrange(YEAR) # Order on year

# bfly_f = arrange(bfly, YEAR) %>% filter((YEAR >= 1979) & (BROOD == 0)) %>%
#   filter(FIRSTDAY != LASTDAY)

# # plot all butterfly arrival time series
# ggplot(bfly_f, aes(x = YEAR, y = event)) + geom_point() +
#   stat_smooth(method = 'glm') + facet_wrap(~SPECIES_NAME)
# 
# # Gather all trends
# test = group_by(bfly_f, SPECIES_NAME) %>% summarize(trend = summary(lm(event ~ YEAR))$coefficients[2])



# Add common column names
bfly_f$species = bfly_f$SPECIES_NAME
bfly_f$year = bfly_f$YEAR
bfly_f$env = bfly_f$Temp











