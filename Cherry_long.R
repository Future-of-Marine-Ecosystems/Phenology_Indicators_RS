# Cherry Blossom Bloom Time Analyses (long term)
# Reid Steele, April 28 2025

# Libraries
library(phenometrics)
library(tidyverse)

# Load in Data
# data = read.csv('./Cherry/cherry.csv')

data <- read.csv("https://ourworldindata.org/grapher/date-of-the-peak-cherry-tree-blossom-in-kyoto.csv?v=1&csvType=full&useColumnShortNames=true")

# Rename columns
colnames(data) = c('entity', 'code', 'year', 'ma', 'yr')

# remove NAs
data = arrange(data, year)
 
# Filter to needed columns
data_filt = select(data, year, yr)


# Calculate expanded dataset
data_exp = NULL
data_ma = NULL
for(i in 1:nrow(data_filt)){

  # Calculate minimum year
  curyear = data_filt$year[i]
  minyear = data_filt$year[i]-19

  # Filter to last 20 years
  ma_data = filter(data_filt, (year <= curyear) & (year >= minyear)) %>%
    filter(is.na(yr) == F)  # Remove NAs

  # Check for min years, calculate moving average if so
  if(nrow(ma_data) >= 5){

    # Calculate
    data_ma = rbind(data_ma, cbind(year = curyear, event = mean(ma_data$yr, na.rm = T)))
    data_exp = rbind(data_exp, cbind(year = curyear, event = ma_data$yr, eyear = ma_data$year))

  } # End if checking for min years

} # End expanded dataset loop

# Set data frames
data_ma = as.data.frame(data_ma); data_exp = as.data.frame(data_exp)

# # Remove very early years
data = left_join(data, data_ma) %>% filter(year >=1000) %>% filter(!is.na(event))

# Test emergence
ks_win = emp_tope(data, alt = 'two.sided', max_y = 50, unemergence = T, emt = 10, quants = c(0.25, 0.75))

# Gather emergences
em = ks_win[ks_win$emerged == 1,]

# Remove NAs
em = filter(em, is.na(em$p) == F)

# Remove multi-year runs
em = filter(em, !((em$year - 1) %in% em$year))

# Plot time series
# data_plot = filter(data_ma, year >= 895)
ggplot(data, aes(x = year, y = event)) + geom_line() +
  geom_vline(xintercept = em$year, color = 'red') + theme_classic() +
  geom_point(inherit.aes = F, aes(x = year, y = yr))

# Apply unemergence
ks_win_unem = ks_win

# Loop through rows
for(i in 2:(nrow(ks_win_unem))){
  
  # check if previous row is emerged
  if(ks_win_unem[i-1,'emerged'] == 1){
    
    # Check next 5 values
    crows = ks_win_unem[i:(i+4),]
    
    # Change to 1 unless all 5 values below threshold
    if(all(crows$p < 0.6)){ks_win_unem[i, 'emerged'] = 0} else {ks_win_unem[i, 'emerged'] = 1}
    
  }
  
}



# Calculate consecutive years of positive test results
con = rep(NA, nrow(ks_win))
for(i in 1:length(con)){
  
  # Dont check i-1 if i=1
  if(i == 1){con[i] = ks_win$p[i]} else {
    
    # Set to 0 if test result is 0, else add 1 to previous value
    if(ks_win$p[i] == 0){con[i] = 0} else{con[i] = ks_win$p[i] + con[i-1]}
    
  }
  
}

# Add to ks_win_unem
ks_win_unem = cbind(ks_win_unem, con)



