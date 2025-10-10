# CPR Data Pre-Processing

# Libraries
library(tidyverse)
library(data.table)
library(rerddap)
library(rerddapXtracto)
library(phenometrics)

# Load in data
data = read.csv('./CPR/CPR_DataRequest_FONIDalhousie_03Feb25/CPR_FONIDalhousie_Data_LargeZooplankton_03022025.csv')
szoo = read.csv('./CPR/CPR_DataRequest_FONIDalhousie_03Feb25/CPR_FONIDalhousie_Data_SmallZooplankton_03022025.csv')
phyto = read.csv('./CPR/CPR_DataRequest_FONIDalhousie_03Feb25/CPR_FONIDalhousie_Data_Phytoplankton_03022025.csv')

# Remove redundant columns

# load in species lists
llzoo = read.csv('./CPR/CPR_DataRequest_FONIDalhousie_03Feb25/CPR_FONIDalhousie_List_LargeZooplankton_03022025.csv');llzoo$class = 'lzoo'
lszoo = read.csv('./CPR/CPR_DataRequest_FONIDalhousie_03Feb25/CPR_FONIDalhousie_List_SmallZooplankton_03022025.csv');lszoo$class = 'szoo'
lphyt = read.csv('./CPR/CPR_DataRequest_FONIDalhousie_03Feb25/CPR_FONIDalhousie_List_Phytoplankton_03022025.csv');lphyt$class = 'phyto'
slist = rbind(llzoo, lszoo, lphyt)

# Join
data = left_join(data, szoo) %>% left_join(phyto)

# Filter taxon list to species
slist = slist[grepl(' ', slist$name_worms),]

# # Filter latitude
# data = data[data$Latitude <= 58,]; data = data[data$Latitude >= 55,]

# Left join data

# Grid size
gsize = 4
grad = gsize/2

# Create centerpoints of grid
grid = expand.grid(lon = seq(floor(min(data$Longitude)) - grad, ceiling(max(data$Longitude)) + grad, gsize), 
                   lat = seq(floor(min(data$Latitude)) - grad, ceiling(max(data$Latitude)) + grad, gsize))
# Set to data table
setDT(grid)

# Add identifier
grid$cell = seq(1, nrow(grid), 1)

# Create bounds
bounds <- grid[, .(xl = lon - grad, xu = lon + grad
                   , yl = lat - grad, yu = lat + grad
                   , cell )]



# Zero-pad
data$Month = ifelse(nchar(data$Month) == 1, paste0('0', data$Month), data$Month)
data$Day = ifelse(nchar(data$Day) == 1, paste0('0', data$Day), data$Day)
data$Hour = ifelse(nchar(data$Hour) == 1, paste0('0', data$Hour), data$Hour)
data$Minute = ifelse(nchar(data$Minute) == 1, paste0('0', data$Minute), data$Minute)

# Set to data table
setDT(data)

# Assign grid cell number
out <- bounds[data, .(SampleId, cell), 
              on = .(xl <= Longitude, xu >= Longitude,  yl <= Latitude, yu >= Latitude)]

# Join
data = left_join(data, out) %>%
  mutate(Date = as.POSIXct(paste0(Year, '-', Month, '-', Day), tz = 'UTC')) %>% # Create date column
  left_join(grid)


# Pivot longer
data_long = pivot_longer(data, grep('X', colnames(data)), names_to = 'species', values_to = 'abundance') %>%
  filter(is.na(abundance) == F) %>% # Filter out Nans
  mutate(julianday = yday(Date), month = month(Date)) %>% # Create Julianday 
  arrange(Date) # Order by date

# Fix species code
data_long$species = as.numeric(gsub('X', '', data_long$species))

# Join on species info
data_long = left_join(data_long, slist, by = join_by(species == accepted_id))

# Filter out NA species 
data_long = filter(data_long, is.na(name_worms) == F)

# Rename year
data_long$year = data_long$Year

# Check for appropriate sampling effort
echeck = group_by(data_long, species, month) %>% filter(abundance > 0) %>%
  summarize(nmean = (n()))

# Average monthly samples, remove any under 200
mcheck = group_by(echeck, species) %>% summarize(mmean = mean(nmean), max = max(nmean))

# Filter out monthly samples under 200
data_filt = filter(data_long, species %in% mcheck$species[which(mcheck$max > 200)])

# Filter down mcheck and echeck
echeck = filter(echeck, species %in% mcheck$species[which(mcheck$max > 200)])
mcheck = filter(mcheck, max > 200)

# Check for binomial species
binom = data.frame(species = unique(data_filt$species), binomial = NA, nrow = NA)
binom$class = slist[match(binom$species, slist$accepted_id),'class']

# Loop through species and test for binomial distribution
for(i in 1:nrow(binom)){
  
  # Grab species
  tspec = binom$species[i]
  
  # Filter data to species
  tdata = filter(data_filt, species == tspec)
  
  # round out near zeros
  tdata$abundance = round(tdata$abundance)
  
  # Divide abundance by 1000 for phytoplankton
  if(binom$class[i] == 'phyto'){tdata$abundance = tdata$abundance/1000}
  
  # generate frequency distribution
  tfreq = rep(tdata$julianday, tdata$abundance)
  
  # Check for bimodality
  hist(tfreq, main = tspec)
  
  # Check number of rows
  binom$nrow[i] = nrow(filter(tdata, abundance > 0))
  
  # Progress
  print(i)
  
} # end binomial check loop

# Identify problem species
pspec = c(41, 55, 101, 107, 114, 117, 177, 194, 10624, 10680, 168, 10615, 157, 155, 102, 61, 112, 113, 750, 976, 84, 160, 165, 10641, 10623, 10541, 63)
qspec = c(61, 108, 112, 113, 125, 126, 155, 168, 750, 10615, 10623, 10684, 133, 148) #questionable species I am not entirely sure about

# remove bad species
data_filt = filter(data_filt, !(species %in% pspec))

# Calculate metrics
lz_metrics = group_by(data_filt, species, year, cell) %>%
  mutate(cumab = cumsum(abundance), prop = cumab/sum(abundance), tot = sum(abundance), day = yday(Date)) %>% 
  filter(tot > 0.1)

# Set threshold
threshold = 0.25

# Calculate event (annual threshold)
lz_event = group_by(lz_metrics, species, year, cell) %>%
  filter(tot >= 100) %>%
  summarize(event = julianday[min(which(prop > threshold))], nsamps = n(), totab = sum(tot), ltotab = log(totab)) %>%
  filter(is.na(event) == F) %>%
  mutate(Date = as.POSIXct(paste(event, year), format = '%j %Y', tz = 'UTC')) %>%# Create date column
  left_join(grid) %>% filter(year >= 1985)

# # Remote sense temperature
# dataInfo <- rerddap::info('noaacrwsstDaily', url = 'https://coastwatch.noaa.gov/erddap/') #daily, 0.025 degrees, gaps
# # Option 1 is to directly input tour variables into the Rctracto script
# start = Sys.time()
# event_temp <- rxtracto(dataInfo, parameter = 'analysed_sst', xcoord = lz_event$lon, ycoord = lz_event$lat, tcoord = lz_event$Date,
#                        xlen = gsize, ylen = gsize, progress_bar = T) #zcoord = rep(0., length(fitted.crw$lon)), xlen = .1, ylen = .1,  progress_bar = TRUE)
# end = Sys.time()
# end - start
# 
# # Save
# save(event_temp, file = './CPR/northseaL4.RData')

# # Load in temperature
load('./CPR/northseaL4.RData')




# Join temperature to data
lz_event$env = event_temp$`mean analysed_sst`

lz_event = filter(lz_event, nsamps >= 25)

# Calculate number of observations per species/grid cell
check = group_by(lz_metrics, species, cell, year) %>%
  summarize(n = n()) # %>%
# group_by(cell, species) %>% 
# summarize(mean_n = mean(n))

# Plot histograms
for(i in 1:length(unique(lz_event$species))){
  
  s = unique(lz_event$species)[i]
  
  hist(filter(lz_event, species == s)$ltotab, main = s)
  
} # End histogram loop

# Calculate 5% quantiles
quants = group_by(lz_event, species) %>% summarize(five = quantile(ltotab, 0.05))

# Calculate mean abundance
lz_event = group_by(lz_event, species) %>% mutate(mean = mean(ltotab), median = median(ltotab), stab = (ltotab - mean)/sd(ltotab))

# Calculate outlier ranges
lz_out = group_by(lz_event, species) %>% summarize(iqr = IQR(event), 
                                                   LQ = quantile(event, probs=c(.25), na.rm = T),
                                                   UQ = quantile(event, probs=c(.75), na.rm = T),
                                                   lower = LQ-1.5*iqr,
                                                   upper = UQ+1.5*iqr)
# Join to lz_event
lz_event = left_join(lz_event, lz_out) %>% mutate(outlier = ifelse((event <= upper) & (event >= lower), 0, 1))

# Filter out outliers
cpr_f = lz_event %>% filter(event <= upper) %>% filter(event >= lower)# %>% filter(ltotab >= quantile(ltotab, 0.05))

# Save results
save.image('cpr_results.RData')


# # # Class by year
# # class_by_year(planktys)
# 
ggplot(lz_event, aes(x = year, y = event)) + 
  facet_wrap(~species) + geom_point(aes(color = as.factor(outlier))) + geom_smooth(method = 'lm')
# 
ggplot(lz_event, aes(x = year, y = event, color = stab)) + scale_color_continuous(type = 'viridis') +
  facet_wrap(~species) + geom_point() + geom_smooth(method = 'lm')
# 
ggplot(lz_event, aes(x = year, y = env, color = stab)) + scale_color_continuous(type = 'viridis') +
  facet_wrap(~species) + geom_point() + geom_smooth(method = 'lm')
# 
ggplot(cpr_f, aes(x = year, y = event, color = stab)) + scale_color_continuous(type = 'viridis') +
  facet_wrap(~species) + geom_point() + geom_smooth(method = 'lm')
#
ggplot(cpr_f, aes(x = year, y = env, color = stab)) + scale_color_continuous(type = 'viridis') +
  facet_wrap(~species) + geom_point() + geom_smooth(method = 'lm')
