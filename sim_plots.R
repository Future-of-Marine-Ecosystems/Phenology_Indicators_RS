# Plot simulation results
library(plotly)
library(tidyverse)
library(RColorBrewer)
library(ggtext)

# load in results
load(('sim_results_comb.RData'))

# Add points from bfly and cpr data
source('Bfly Initial.R')
load('cpr_results.RData')

# Gather butterfly metrics
bfly_spread = group_by(bfly_f, species, year) %>%
  summarize(event = mean(event)) %>%
  group_by(species) %>%
  summarize(spread = sd(event))

bfly_sd = group_by(bfly_f, species, year) %>%
  summarize(sd = sd(event)) %>%
  group_by(species) %>%
  summarize(mean_sd = mean(sd, na.rm = T))

bfly_samples = group_by(bfly_f, species, year) %>%
  summarize(n = n()) %>%
  group_by(species) %>%
  summarize(samples = mean(n))

bfly_trend = group_by(bfly_f, species) %>%
  summarize(trend = lm(event ~ year)$coefficients[2])

# Check relationship between spread and SD
summary(lm(bfly_sd$mean_sd ~ bfly_spread$spread))
# sd ~ 1.5x spread

# join together
bfly_params = left_join(bfly_spread, bfly_sd) %>%
  left_join(bfly_samples) %>%
  left_join(bfly_trend) %>%
  mutate(dataset = 'UKBMS')





# Gather cpr metrics
cpr_spread = group_by(cpr_f, species, year) %>%
  summarize(event = mean(event)) %>%
  group_by(species) %>%
  summarize(spread = sd(event))

cpr_sd = group_by(cpr_f, species, year) %>%
  summarize(sd = sd(event)) %>%
  group_by(species) %>%
  summarize(mean_sd = mean(sd, na.rm = T))

cpr_samples = group_by(cpr_f, species, year) %>%
  summarize(n = n()) %>%
  group_by(species) %>%
  summarize(samples = mean(n))

cpr_trend = group_by(cpr_f, species) %>%
  summarize(trend = lm(event ~ year)$coefficients[2])

# Check relationship between spread and SD
summary(lm(cpr_sd$mean_sd ~ cpr_spread$spread))
# sd ~ 1.3x spread

# join together
cpr_params = left_join(cpr_spread, cpr_sd) %>%
  left_join(cpr_samples) %>%
  left_join(cpr_trend) %>%
  mutate(dataset = 'CPR')



# Combine into one dataset
data_params = rbind(bfly_params, cpr_params)



# Limits
range(data_params$samples)
range(data_params$spread)
range(data_params$trend)

# Filter to simulation parameter range
data_params = filter(data_params, spread <= max(sim_params$spread)) %>%
  filter(samples <= max(sim_params$samples)) %>%
  filter(trend <= max(sim_params$trend)) %>%
  filter(spread >= min(sim_params$spread)) %>%
  filter(samples >= min(sim_params$samples)) %>%
  filter(trend >= min(sim_params$trend))



# test values
tspread = spread[round(length(spread)/2)]
tsamples = samples[round(length(samples)/2)]
ttrend = trend[round(length(trend)/2)]

# filter to test spread
t_sim_params_spread = filter(sim_params, spread == tspread) %>%
  select(samples,trend,eyear_stat,eyear_emp)

# Plot statistical
samp_x_trend_stat = ggplot(t_sim_params_spread, aes(x = samples, y = trend, fill = eyear_stat)) + geom_raster() +
  geom_point(data = data_params, inherit.aes = F, aes(x = samples, y = trend, color = dataset), size = 5, show.legend = F) + scale_color_hue(h = c(15, 300)) +
  scale_fill_continuous(type = 'viridis', limits = c(0,40)) + theme_classic() + labs(x = 'Replicates (*r*)', y = 'Trend (*t*)', fill = 'Years to Emergence\n ') +
  theme(legend.position = 'right', legend.key.height= unit(2.5, 'cm'), legend.key.width= unit(2, 'cm'), 
        legend.direction = 'vertical',
        legend.text = element_text(size=35), legend.title = element_text(size = 40, angle = 270, vjust = 0.5),
        axis.text = element_markdown(size = 35), axis.title = element_markdown(size = 35), 
        plot.title = element_markdown(size = 35, hjust = 0.5, vjust = 0)) + 
  guides(fill = guide_colourbar(title.position="right", title.hjust = 0.5)) + ggtitle(paste0('Spread (*s*) = ', tspread))

# Plot empirical
samp_x_trend_emp = ggplot(t_sim_params_spread, aes(x = samples, y = trend, fill = eyear_emp)) + geom_raster(show.legend = F) +
  geom_point(data = data_params, inherit.aes = F, aes(x = samples, y = trend, color = dataset), size = 5) + scale_color_hue(h = c(15, 300)) +
  scale_fill_continuous(type = 'viridis', limits = c(0,40)) + theme_classic() + labs(x = 'Replicates (*r*)', y = 'Trend (*t*)', color = 'Dataset') +
  theme(legend.position = 'right', legend.key.height= unit(2.5, 'cm'), legend.key.width= unit(2, 'cm'), 
        legend.direction = 'vertical',
        legend.text = element_markdown(size=35), legend.title = element_text(size = 40),
        axis.text = element_markdown(size = 35), axis.title = element_markdown(size = 35), 
        plot.title = element_markdown(size = 35, hjust = 0.5)) + 
  guides(fill = guide_colourbar(title.position="right", title.hjust = 0.5))

# # Make into matrix
# t_sim_mat = spread(t_sim_params, trend, eyear_stat)
# rownames(t_sim_mat) = t_sim_mat[,1]; t_sim_mat = t_sim_mat[,-1]

# filter to test samples
t_sim_params_samples = filter(sim_params, samples == tsamples) %>%
  select(spread,trend,eyear_stat,eyear_emp)

# Plot statistical
spread_x_trend_stat = ggplot(t_sim_params_samples, aes(x = spread, y = trend, fill = eyear_stat)) + geom_raster() +
  geom_point(data = data_params, inherit.aes = F, aes(x = spread, y = trend, color = dataset), size = 5, show.legend = F) + scale_color_hue(h = c(15, 300)) +
  scale_fill_continuous(type = 'viridis', limits = c(0,40)) + theme_classic() + labs(x = 'Spread (*s*)', y = 'Trend (*t*)', fill = 'Years to Emergence') + 
  theme(axis.text = element_markdown(size = 35), axis.title = element_markdown(size = 35), 
        plot.title = element_markdown(size = 35, hjust = 0.5)) + ggtitle(paste0('Replicates (*r*) = ', tsamples))

# Plot empirical
spread_x_trend_emp = ggplot(t_sim_params_samples, aes(x = spread, y = trend, fill = eyear_emp)) + geom_raster() +
  geom_point(data = data_params, inherit.aes = F, aes(x = spread, y = trend, color = dataset), size = 5, show.legend = F) + scale_color_hue(h = c(15, 300)) +
  scale_fill_continuous(type = 'viridis', limits = c(0,40)) + theme_classic() + labs(x = 'Spread (*s*)', y = 'Trend (*t*)', fill = 'Years to Emergence') + 
  theme(axis.text = element_markdown(size = 35), axis.title = element_markdown(size = 35))

# # Make into matrix
# t_sim_mat = spread(t_sim_params, trend, eyear_stat)
# rownames(t_sim_mat) = t_sim_mat[,1]; t_sim_mat = t_sim_mat[,-1]

# filter to test trend
t_sim_params_trend = filter(sim_params, trend == ttrend) %>%
 select(samples,spread,eyear_stat,eyear_emp)

# Plot statistical
samp_x_spread_stat = ggplot(t_sim_params_trend, aes(x = samples, y = spread, fill = eyear_stat)) + geom_raster() +
  geom_point(data = data_params, inherit.aes = F, aes(x = samples, y = spread, color = dataset), size = 5, show.legend = F) + scale_color_hue(h = c(15, 300)) +
  scale_fill_continuous(type = 'viridis', limits = c(0,40)) + theme_classic()  + labs(x = 'Replicates (*r*)', y = 'Spread (*s*)', fill = 'Years to Emergence') + 
  theme(axis.text = element_markdown(size = 35), axis.title = element_markdown(size = 35), 
        plot.title = element_markdown(size = 35, hjust = 0.5, vjust = 0)) + ggtitle(paste0('Trend (*t*) = ', ttrend))

# Plot empirical
samp_x_spread_emp = ggplot(t_sim_params_trend, aes(x = samples, y = spread, fill = eyear_emp)) + geom_raster() +
  geom_point(data = data_params, inherit.aes = F, aes(x = samples, y = spread, color = dataset), size = 5, show.legend = F) + scale_color_hue(h = c(15, 300)) +
  scale_fill_continuous(type = 'viridis', limits = c(0,40)) + theme_classic()  + labs(x = 'Replicates (*r*)', y = 'Spread (*s*)', fill = 'Years to Emergence') + 
  theme(axis.text = element_markdown(size = 35), axis.title = element_markdown(size = 35))


# # Make into matrix
# t_sim_mat = spread(t_sim_params, spread, eyear_stat)
# rownames(t_sim_mat) = t_sim_mat[,1]; t_sim_mat = t_sim_mat[,-1]

# plot statistcal vs empirical results
sim_params_comp = sim_params

# Replace NAs 

# Calculate difference in mean emergence time
sim_params_comp$diff = sim_params_comp$eyear_stat - sim_params_comp$eyear_emp

# filter out NAs
sim_params_comp_na = filter(sim_params_comp, is.na(diff)) %>%
  filter((!is.na(eyear_stat)) | (!is.na(eyear_emp))) %>%
  mutate(eyear_stat_cor = ifelse(is.na(eyear_stat), 50, eyear_stat), eyear_emp_cor = ifelse(is.na(eyear_emp), 50, eyear_emp)) %>%
  mutate(diff = eyear_stat_cor - eyear_emp_cor)


# filter to test spread
t_sim_comp_spread = filter(sim_params_comp, spread == tspread) %>%
  select(samples,trend,eyear_stat,eyear_emp,diff)

# filter to test spread
t_sim_comp_spread_na = filter(sim_params_comp_na, spread == tspread) %>%
  select(samples,trend,eyear_stat,eyear_emp,diff)

# Plot difference
samp_x_trend_comp = ggplot(filter(t_sim_comp_spread, !is.na(diff)), aes(x = samples, y = trend, fill = diff)) + geom_raster() +
  geom_point(data = t_sim_comp_spread_na, size = 4, pch = 21) + scale_color_distiller(palette = 'PRGn') +
  scale_fill_distiller(palette = 'PRGn', limits = c(-35, 35)) + theme_classic() + labs(x = 'Replicates (*r*)', y = 'Trend (*t*)', fill = 'Statistical - Empirical\n ') +
  theme(legend.position = 'right', legend.key.height= unit(2.5, 'cm'), legend.key.width= unit(2, 'cm'), 
        legend.direction = 'vertical',
        legend.text = element_text(size=35), legend.title = element_text(size = 40, angle = 270, vjust = 0.5),
        axis.text = element_markdown(size = 35), axis.title = element_markdown(size = 35)) + 
  guides(fill = guide_colourbar(title.position="right", title.hjust = 0.5))






# filter to test samples
t_sim_comp_samples = filter(sim_params_comp, samples == tsamples) %>%
  select(spread,trend,eyear_stat,eyear_emp,diff)

# filter to test samples
t_sim_comp_samples_na = filter(sim_params_comp_na, samples == tsamples) %>%
  select(spread,trend,eyear_stat,eyear_emp,diff)

# Plot difference
spread_x_trend_comp = ggplot(filter(t_sim_comp_samples, !is.na(diff)), aes(x = spread, y = trend, fill = diff)) + geom_raster() +
  geom_point(data = t_sim_comp_samples_na, size = 4, pch = 21) + scale_color_distiller(palette = 'PRGn') +
  scale_fill_distiller(palette = 'PRGn', limits = c(-35, 35)) + theme_classic() + labs(x = 'Spread (*s*)', y = 'Trend (*t*)', fill = 'Statistical - Emprical') +
  theme(axis.text = element_markdown(size = 35), axis.title = element_markdown(size = 35))





# filter to test trend
t_sim_comp_trend = filter(sim_params_comp, trend == ttrend) %>%
  select(samples,spread,eyear_stat,eyear_emp,diff)

# filter to test trend
t_sim_comp_trend_na = filter(sim_params_comp_na, trend == ttrend) %>%
  select(samples,spread,eyear_stat,eyear_emp,diff)

# Plot difference
samp_x_spread_comp = ggplot(filter(t_sim_comp_trend, !is.na(diff)), aes(x = samples, y = spread, fill = diff)) + geom_raster() +
  geom_point(data = t_sim_comp_trend_na, size = 4, pch = 21) + scale_color_distiller(palette = 'PRGn') +
  scale_fill_distiller(palette = 'PRGn', limits = c(-35, 35)) + theme_classic() + labs(x = 'Replicates (*r*)', y = 'Spread (*s*)', fill = 'Statistical - Emprical') +
  theme(axis.text = element_markdown(size = 35), axis.title = element_markdown(size = 35))









