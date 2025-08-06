# Figure output for chapter 1
# June 30, 2025

# Libraries
library(ggpubr)
library(ggbreak)
library(data.table)
library(phenometrics)

# Time of emergence simulations
source('sim_plots.R')

# Gather legend s
leg = get_legend(samp_x_trend_stat)
leg2 = get_legend(samp_x_trend_emp)
leg3 = get_legend(samp_x_trend_comp)

# Blank space filler
space = ggplot() + theme_classic()

# export figure
tiff(filename = './Figures/simulations.tif', width = 2000, height = 1500, units = 'px')

# Arrange simulation plots
simfig = ggarrange(samp_x_trend_stat + theme(legend.position = "none"), 
              spread_x_trend_stat + theme(legend.position = "none"), 
              samp_x_spread_stat + theme(legend.position = "none"), 
              leg,
              samp_x_trend_emp + theme(legend.position = "none"), 
              spread_x_trend_emp + theme(legend.position = "none"), 
              samp_x_spread_emp + theme(legend.position = "none"), 
              leg2,
              samp_x_trend_comp + theme(legend.position = "none"), 
              spread_x_trend_comp + theme(legend.position = "none"), 
              samp_x_spread_comp + theme(legend.position = "none"), 
              leg3,
              ncol = 4, nrow = 3, widths = c(5,5,5,2))

annotate_figure(simfig, left = text_grob('       Delta                           Empirical                      Statistical', size = 45, rot = 90, vjust = 0.35))

dev.off()


# Cherry Blossom
source('Cherry_long.R')

# Join emergence to data_plot
data_plot = left_join(data, ks_win_unem)

# Gather emerged years
emerged_plot = filter(data_plot, emerged == 1)

# Calculate block starts/ends
max_ind = which((data_plot$emerged == 1) & ((shift(data_plot$emerged, type = 'lead') == 0) | (is.na(shift(data_plot$emerged, type = 'lead')))))
min_ind = which((data_plot$emerged == 1) & ((shift(data_plot$emerged, type = 'lag') == 0) | (is.na(shift(data_plot$emerged, type = 'lag')))))

# Create emergence blocks
emerged_block = tibble(xmin = data_plot[min_ind, 'year'], xmax = data_plot[max_ind, 'year'])
# data_block = tibble(xmin = emerged_plot[min_ind, 'year']-19, xmax = emerged_plot[max_ind, 'year'])
# center_block = tibble(xmin = emerged_plot[min_ind, 'year']-25, xmax = emerged_plot[max_ind, 'year']-25)
# line_block = tibble(x = c(emerged_plot[min_ind, 'year']-49, emerged_plot[min_ind, 'year']),
#                     y = rep(seq(97, 95, -0.5), 2), block = rep(seq(1,5,1), 2))
# y_block = tibble(x = rep(emerged_plot[min_ind, 'year']-49, 2),
#                  y = c(seq(97, 95, -0.5) - 0.25, seq(97, 95, -0.5) + 0.25), block = rep(seq(1,5,1), 2))

# Plot times series
ts = ggplot(data_plot) + geom_point(aes(x = year, y = event), size = 3) + labs(x = 'Year', y = 'Bloom Timing (Julian Day)') +
  theme_classic() + geom_rect(data = emerged_block, aes(ymin = -Inf, ymax = Inf, xmin = xmin, xmax = xmax), fill = 'red', alpha = 0.3) +
  geom_line(data = data_ma, aes(x = year, y = event), lwd = 2) +
  # geom_line(data = line_block, aes(x = x, y = y, group = block), color = 'red', lwd = 1) +
  # geom_line(data = y_block, aes(x = x, y = y, group = block), color = 'red', lwd = 1) +
  theme(axis.text = element_text(size = 35), axis.title = element_text(size = 35))
  #geom_errorbar(data = line_block, aes(xmin = xmin, xmax = xmax, y = y, group = block), color = 'red')#+
  #geom_rect(data = data_block, aes(ymin = -Inf, ymax = Inf, xmin = xmin, xmax = xmax), fill = 'red', alpha = 0.3) 

# Plot ks test power
power = ggplot(data_plot) + geom_line(aes(x = year, y = con), lwd = 1) + theme_classic() + 
  geom_hline(yintercept = 5, linetype = 'dashed', alpha = 0.3, lwd = 1) + labs(x = 'Year', y = 'Consecutive Positive Test Results') +
  geom_rect(data = emerged_block, aes(ymin = -Inf, ymax = Inf, xmin = xmin, xmax = xmax), fill = 'red', alpha = 0.3) +
  theme(axis.text = element_text(size = 35), axis.title = element_text(size = 35))



# export figure
tiff(filename = './Figures/cherry.tif', width = 1500, height = 1500, units = 'px')

# ggarrange
ggarrange(ts, power, ncol = 1, nrow = 2, align = 'v')

dev.off()



# Butterflies
source('Bfly Initial.R')

# colour scale
cscale = c(hcl(h=240, 0, 95, alpha = 0.8), # No signal
           hcl(h=60, 100, 100, alpha = 0.8), # Combination
           hcl(h=12, 200, 100, alpha = 0.8), # decouple
           hcl(h=120, 100, 100, alpha = 0.8))

# Load in trait data
bfly_traits = read.csv('./Butterflies/bfly_traits/data/ecological_traits_2022.csv', skip=1)

# Multivoltine species
species_multi = filter(bfly_traits, obligate_multivoltine == 1) %>%
  filter(scientific_name %in% bfly_f$species)

# Univoltine species
species_uni = filter(bfly_traits, scientific_name %in% bfly_f$species) %>%
  filter(!(scientific_name %in% species_multi$scientific_name))

# Run metrics
comm_all = community(bfly_f, method = 'statistical')
cbs_all = class_by_species(comm_all, plot = F)
cby_all = class_by_year(comm_all, plot = T)

# add voltine column
cbs_all$voltine = ifelse(cbs_all$species %in% species_multi$scientific_name, 'multivoltine', 'univoltine')

# Rearrange, including voltine
cbs_reorder = cbs_all %>% arrange(voltine, final, year_fc, species) %>% 
  mutate(order = match(species, unique(species)))



# Univoltine/multivoltine border
mv_uv = c(max(cbs_reorder[cbs_reorder$voltine == 'multivoltine', 'order'])+0.5,
          max(cbs_reorder[cbs_reorder$voltine == 'multivoltine', 'order'])+0.5)

# Univoltine/Bivoltine label positions
voltine_labs = data.frame(x = min(cbs_reorder$year)-2, 
                          y = c(median(cbs_reorder[cbs_reorder$voltine == 'multivoltine', 'order']),
                                median(cbs_reorder[cbs_reorder$voltine == 'univoltine', 'order'])),
                          label = c('Multivoltine', 'Univoltine'))

# Make NA into 0
cbs_reorder$class_c = ifelse(is.na(cbs_reorder$class_c), 0, cbs_reorder$class_c)

# Plot class by species
cbs_plot = ggplot(cbs_reorder, aes(x = year, y = order)) +
  geom_raster(aes(fill = as.character(class_c))) +
  scale_fill_manual(values = c(cscale[c(1,4,3,2)], 'grey50'), 
                    breaks = c('1' ,'2', '3', '4'),
                    labels = c('No Signal', 'Shift', 'Decouple', 'Undershift')) +
  theme_classic() + labs(x = 'Year', y = 'Species', fill = 'Classification') +
  scale_y_continuous(breaks = unique(cbs_reorder$order), labels = unique(cbs_reorder$species), position = 'right') + 
  theme(axis.text.y.left = element_blank(),
        axis.line.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 25),
        axis.text.y.right = element_text(size = 20), axis.title.y.right = element_text(size = 20),
        legend.text = element_text(size = 20), legend.title = element_text(size = 25)) +
  theme(legend.position = 'left') + 
  geom_line(data = filter(cbs_reorder, (class_c == final) & (class != 'no signal')), aes(x = year, y = order, group = species), inherit.aes = F) +
  geom_point(data = filter(cbs_reorder, class != 'no signal'), aes(x = year_fc, y = order, group = species), inherit.aes = F) +
  geom_text(data = voltine_labs, aes(x = x, y = y, label = label), inherit.aes = F, angle = 90, size = 12) + 
  scale_y_break(mv_uv, scales = 'fixed', expand = F, space = 1)  

# export figure
tiff(filename = './Figures/bfly_cbs.tif', width = 1500, height = 1500, units = 'px')

cbs_plot
  
# end figure
dev.off()



# CPR
load('cpr_results.RData')

# Run community
cpr_comm = community(cpr_f, quants = c(0.25, 0.75))

# Run class by species
planktys_sim_cbs = class_by_species(cpr_comm)

# Insert species name
planktys_sim_cbs$species = slist[match(planktys_sim_cbs$species, slist$accepted_id), 'name_worms']

# Gather class
planktys_sim_cbs$type = slist[match(planktys_sim_cbs$species, slist$name_worms), 'class']

# Change to phytoplankton and zooplankton
planktys_sim_cbs$type = ifelse(planktys_sim_cbs$type == 'phyto', 'Phytoplankton', 'Zooplankton')

# Rearrange, including type
planktys_sim_cbs_reorder = planktys_sim_cbs %>% arrange(type, final, year_fc, species) %>% 
  mutate(order = match(species, unique(species)))

# Set data frame
planktys_sim_cbs_reorder = as.data.frame(planktys_sim_cbs_reorder)

# phytoplankton/zooplanlton border
phy_zoo = c(max(planktys_sim_cbs_reorder[planktys_sim_cbs_reorder$type == 'Phytoplankton', 'order']+0.5),
          max(planktys_sim_cbs_reorder[planktys_sim_cbs_reorder$type == 'Phytoplankton', 'order'])++0.5)

# phytoplankton/zooplanlton label positions
plankton_labs = data.frame(x = min(planktys_sim_cbs_reorder$year)-5, 
                          y = c(median(planktys_sim_cbs_reorder[planktys_sim_cbs_reorder$type == 'Phytoplankton', 'order']),
                                median(planktys_sim_cbs_reorder[planktys_sim_cbs_reorder$type == 'Zooplankton', 'order'])),
                          label = c('Phytoplankton', 'Zooplankton'))

# Make NA into 0
planktys_sim_cbs_reorder$class_c = ifelse(is.na(planktys_sim_cbs_reorder$class_c), 0, planktys_sim_cbs_reorder$class_c)

# Plot class by species
plankton_cbs_plot = ggplot(planktys_sim_cbs_reorder, aes(x = year, y = order)) +
  geom_raster(aes(fill = as.character(class_c))) +
  scale_fill_manual(values = c(cscale[c(1,4,3,2)], 'grey50'), 
                    breaks = c('1' ,'2', '3', '4'),
                    labels = c('No Signal', 'Shift', 'Decouple', 'Undershift')) +
  theme_classic() + labs(x = 'Year', y = 'Species', fill = 'Classification') +
  scale_y_continuous(breaks = unique(planktys_sim_cbs_reorder$order), labels = unique(planktys_sim_cbs_reorder$species), position = 'right') + 
  theme(axis.text.y.left = element_blank(),
        axis.line.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 25),
        axis.text.y.right = element_text(size = 20), axis.title.y.right = element_text(size = 20),
        legend.text = element_text(size = 20), legend.title = element_text(size = 25)) +
  theme(legend.position = 'left') + 
  geom_line(data = filter(planktys_sim_cbs_reorder, (class_c == final) & (class != 'no signal')), aes(x = year, y = order, group = species), inherit.aes = F) +
  geom_point(data = filter(planktys_sim_cbs_reorder, class != 'no signal'), aes(x = year_fc, y = order, group = species), inherit.aes = F) +
  geom_text(data = plankton_labs, aes(x = x, y = y, label = label), inherit.aes = F, angle = 90, size = 12) + 
  scale_y_break(phy_zoo, scales = 'fixed', expand = F, space = 1)  

# export figure
tiff(filename = './Figures/cpr_cbs.tif', width = 1500, height = 1500, units = 'px')

plankton_cbs_plot

# end figure
dev.off()




# Class by year figures
planktys_cby = class_by_year(cpr_comm)

# Select and pivot
planktys_cby = select(planktys_cby, year, n_em, n_dc) %>%
  pivot_longer(cols = -1, names_to = 'emergence', values_to = 'proportion')
cby_all = select(cby_all, year, n_em, n_dc) %>%
  pivot_longer(cols = -1, names_to = 'emergence', values_to = 'proportion')

# Tie in identifiers
planktys_cby$name = 'CPR'
cby_all$name = 'UKBMS'

# Combine name and emergence into group for plotting
cby_comb = rbind(planktys_cby, cby_all)
cby_comb$group = paste(cby_comb$emergence, cby_comb$name, sep = '_')

# Reorder emergence
cby_comb$emergence = factor(cby_comb$emergence, levels = c('n_em', 'n_dc'))

# Plot
cby_plot = ggplot(cby_comb, aes(x = year, y = proportion, group = group, color = emergence)) +
  geom_line(aes(linetype = name), linewidth = 3) + theme_classic() +
  scale_color_manual(values = c(cscale[c(4,3)]), labels = c('POPE', 'POTE')) +
  theme(legend.title = element_blank(), axis.text = element_text(size = 35), 
        axis.title = element_text(size = 35), legend.text = element_text(size = 35), 
        legend.key.width = unit(3, 'cm'), legend.position = 'inside', legend.position.inside = c(0.9, 0.9),
        plot.margin = unit(c(1,1,1,1), "cm")) +
  labs(y = 'Proportion of Species Emerged', x = 'Year') + ylim(c(0,1))

# export figure
tiff(filename = './Figures/cby.tif', width = 1500, height = 1500, units = 'px')

cby_plot

# end figure
dev.off()




# Mismatch Figures
bfly_cm = comm_mismatch(bfly_f, method = 'statistical', base_y = 10)

# run by species
bfly_cmbs = comm_mm_by_spp(bfly_cm)

# add voltine column
bfly_cmbs$voltine = ifelse(bfly_cmbs$sp1 %in% species_multi$scientific_name, 'multivoltine', 'univoltine')

# Order species by mismatch ratio in final year, descending
sp_order = filter(bfly_cmbs, year == max(year)) %>% 
  arrange(voltine, desc(mismatch))

# Calculate species order for plotting
bfly_cmbs_reorder = bfly_cmbs %>% mutate(order = match(sp1, unique(sp_order$sp1))) %>%
  arrange(order, year)

# set data frame
bfly_cmbs_reorder = as.data.frame(bfly_cmbs_reorder)

# Univoltine/multivoltine border
mv_uv = c(max(bfly_cmbs_reorder[bfly_cmbs_reorder$voltine == 'multivoltine', 'order'])+0.5,
          max(bfly_cmbs_reorder[bfly_cmbs_reorder$voltine == 'multivoltine', 'order'])+0.5)

# Univoltine/Bivoltine label positions
voltine_labs = data.frame(x = min(bfly_cmbs_reorder$year)-2, 
                          y = c(median(bfly_cmbs_reorder[bfly_cmbs_reorder$voltine == 'multivoltine', 'order']),
                                median(bfly_cmbs_reorder[bfly_cmbs_reorder$voltine == 'univoltine', 'order'])),
                          label = c('Multivoltine', 'Univoltine'))

# Plot mismatch by species matrix
bfly_cmbs_plot = ggplot(bfly_cmbs_reorder, aes(x = year, y = order)) +
  geom_raster(aes(fill = mismatch)) +
  scale_fill_continuous(type = 'viridis') +
  theme_classic() + labs(x = 'Year', y = 'Species', fill = 'Proportion\nMismatched') +
  scale_y_continuous(breaks = unique(bfly_cmbs_reorder$order), labels = unique(bfly_cmbs_reorder$sp1), position = 'right') + 
  theme(axis.text.y.left = element_blank(),
        axis.line.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 25),
        axis.text.y.right = element_text(size = 20), axis.title.y.right = element_text(size = 20),
        legend.text = element_text(size = 20), legend.title = element_text(size = 25, margin = margin(0,0,1,0,'cm')), legend.position = 'left',
        legend.key.height= unit(3, 'cm'), legend.key.width = unit(2, 'cm')) + 
  geom_text(data = voltine_labs, aes(x = x, y = y, label = label), inherit.aes = F, angle = 90, size = 12) + 
  scale_y_break(mv_uv, scales = 'fixed', expand = F, space = 1)  

# export figure
tiff(filename = './Figures/bfly_cmbs.tif', width = 1500, height = 1500, units = 'px')

bfly_cmbs_plot

# end figure
dev.off()



# Mismatch Figures
planktys_sim_mm = comm_mismatch(cpr_f, quants = c(0.25, 0.75), base_y = 10)

# run by species
planktys_cmbs = comm_mm_by_spp(planktys_sim_mm)

# Insert species name
planktys_cmbs$sp1 = slist[match(planktys_cmbs$sp1, slist$accepted_id), 'name_worms']

# add type column
planktys_cmbs$type = slist[match(planktys_cmbs$sp1, slist$name_worms), 'class']

# Change to phytoplankton and zooplankton
planktys_cmbs$type = ifelse(planktys_cmbs$type == 'phyto', 'Phytoplankton', 'Zooplankton')

# Order species by mismatch ratio in final year, descending
sp_order = filter(planktys_cmbs, year == max(year)) %>% 
  arrange(type, desc(mismatch))

# Calculate species order for plotting
planktys_cmbs_reorder = planktys_cmbs %>% mutate(order = match(sp1, unique(sp_order$sp1))) %>%
  arrange(order, year)

# set data frame
planktys_cmbs_reorder = as.data.frame(planktys_cmbs_reorder)

# phytoplankton/zooplanlton border
phy_zoo = c(max(planktys_cmbs_reorder[planktys_cmbs_reorder$type == 'Phytoplankton', 'order']+0.5),
            max(planktys_cmbs_reorder[planktys_cmbs_reorder$type == 'Phytoplankton', 'order'])+0.5)

# phytoplankton/zooplanlton label positions
plankton_labs = data.frame(x = min(planktys_cmbs_reorder$year)-5, 
                           y = c(median(planktys_cmbs_reorder[planktys_cmbs_reorder$type == 'Phytoplankton', 'order']),
                                 median(planktys_cmbs_reorder[planktys_cmbs_reorder$type == 'Zooplankton', 'order'])),
                           label = c('Phytoplankton', 'Zooplankton'))

# Plot mismatch by species matrix
plankton_cmbs_plot = ggplot(planktys_cmbs_reorder, aes(x = year, y = order)) +
  geom_raster(aes(fill = mismatch)) +
  scale_fill_continuous(type = 'viridis') +
  theme_classic() + labs(x = 'Year', y = 'Species', fill = 'Proportion\nMismatched') +
  scale_y_continuous(breaks = unique(planktys_cmbs_reorder$order), labels = unique(planktys_cmbs_reorder$sp1), position = 'right') + 
  theme(axis.text.y.left = element_blank(),
        axis.line.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 25),
        axis.text.y.right = element_text(size = 20), axis.title.y.right = element_text(size = 20),
        legend.text = element_text(size = 20), legend.title = element_text(size = 25, margin = margin(0,0,1,0,'cm')), legend.position = 'left',
        legend.key.height= unit(3, 'cm'), legend.key.width = unit(2, 'cm')) + 
  geom_text(data = plankton_labs, aes(x = x, y = y, label = label), inherit.aes = F, angle = 90, size = 12) + 
  scale_y_break(phy_zoo, scales = 'fixed', expand = F, space = 1)  

# export figure
tiff(filename = './Figures/cpr_cmbs.tif', width = 1500, height = 1500, units = 'px')

plankton_cmbs_plot

# end figure
dev.off()






# Mismatch by year figures
bfly_cmby = comm_mm_by_year(bfly_cm)
plankton_cmby = comm_mm_by_year(planktys_sim_mm)

# Add group
bfly_cmby$group = 'UKBMS'
plankton_cmby$group = 'CPR'

# Bind
cmby_comb = rbind(bfly_cmby, plankton_cmby)

# Plot
cmby_plot = ggplot(cmby_comb, aes(x = year, y = mismatch, linetype = group)) + 
  geom_line(lwd = 2) + theme_classic() +
  labs(x = 'Year', y = 'Proportion of Species Mismatched (POMS)') + 
  theme(legend.title = element_blank(), axis.text = element_text(size = 35), 
        axis.title = element_text(size = 35), legend.text = element_text(size = 35), 
        legend.key.width = unit(3, 'cm'), legend.position = 'inside', legend.position.inside = c(0.9, 0.9),
        plot.margin = unit(c(1,1,1,1), "cm")) +
  labs(y = 'Proportion of Species Pairs Mismatched', x = 'Year') + ylim(c(0,1))

# export figure
tiff(filename = './Figures/cmby.tif', width = 1500, height = 1500, units = 'px')

cmby_plot

# end figure
dev.off()





# Example time series plot
ex_data = filter(bfly_f, species == 'Vanessa atalanta')

# summarize
ex_data_summ = group_by(ex_data, year) %>% 
  summarize(mean_event = mean(event), sd_event = sd(event)) %>%
  mutate(lower = mean_event - sd_event, upper = mean_event + sd_event)

# Plot 
ex_plot = ggplot(filter(ex_data_summ, year >= max(year)-15), aes(x = year, y = mean_event)) + geom_point(size = 5) +
  geom_errorbar(inherit.aes = F, aes(x = year, ymin = lower, ymax = upper), linewidth = 2, width = 0.2) +
  stat_smooth(method = 'lm', se = F, linewidth = 3) +
  theme_classic() + labs(x = 'Year', y = 'Event Timing (Julian Day)') + 
  theme(legend.title = element_blank(), axis.text = element_text(size = 35), 
        axis.title = element_text(size = 35))

# export figure
tiff(filename = './Figures/test_ts.tif', width = 1500, height = 1000, units = 'px')

ex_plot

# end figure
dev.off()










