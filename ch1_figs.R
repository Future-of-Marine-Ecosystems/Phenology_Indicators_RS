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

annotate_figure(simfig, 
                top = text_grob('    Trend vs. Replicates               Trend vs. Spread            Spread vs. Replicates', size = 45, hjust = 0.55, vjust = 0.4),
                left = text_grob('       Delta                           Empirical                      Statistical   ', size = 45, rot = 90, vjust = 0.35))

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
ts = ggplot(data_plot) + geom_point(aes(x = year, y = yr), size = 3) + labs(x = 'Year', y = 'Bloom Timing (Julian Day)') +
  theme_classic() + geom_rect(data = emerged_block, aes(ymin = -Inf, ymax = Inf, xmin = xmin, xmax = xmax), fill = 'red', alpha = 0.3) +
  geom_line(aes(x = year, y = event), lwd = 2) +
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
tiff(filename = './Figures/cherry_ma.tif', width = 1500, height = 1500, units = 'px')

# ggarrange
ggarrange(ts, power, ncol = 1, nrow = 2, align = 'v')

dev.off()



# Butterflies
source('Bfly Initial.R')

# colour scale
cscale = c(hcl(h=240, 0, 95, alpha = 0.8), # No signalhttp://127.0.0.1:29295/graphics/32674a51-291a-4d39-bcf3-872dcac5926c.png
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
  scale_color_manual(values = c(cscale[c(4,3)]), labels = c('POPES', 'POEES')) +
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
        plot.margin = unit(c(1,1,1,1), "cm")) + ylim(c(0,1))

# export figure
tiff(filename = './Figures/cmby.tif', width = 1500, height = 1500, units = 'px')

cmby_plot

# end figure
dev.off()




# Panel class by species plots

# Adjust plots for panelling

# Match reorders
bfly_cmbs_reorder_p = bfly_cmbs %>% mutate(order = match(sp1, unique(cbs_reorder$species))) %>%
  arrange(order, year)

planktys_cmbs_reorder_p = planktys_cmbs %>% mutate(order = match(sp1, unique(planktys_sim_cbs_reorder$species))) %>%
  arrange(order, year)

# Adjust labs
voltine_labs_p = voltine_labs
voltine_labs_p$x = 1979

plankton_labs_p = plankton_labs
plankton_labs_p$x = 1979
plankton_labs_p$y = c(6, 14)
plankton_labs_p$label = c('Phyto', 'Zoo')

# add spaces to cpr y axis
max(nchar(unique(cbs_reorder$species)))
max(nchar(unique(planktys_cmbs_reorder_p$sp1)))
planktys_sim_cbs_reorder_p = planktys_sim_cbs_reorder
planktys_sim_cbs_reorder_p$species = gsub('Rhizosolenia semispina', 'Rhizosolenia semispina       ', planktys_sim_cbs_reorder_p$species)


# Plot class by species
cbs_plot_p = ggplot(cbs_reorder, aes(x = year, y = order)) +
  geom_tile(aes(fill = as.character(class_c))) +
  scale_fill_manual(values = c(cscale[c(1,4,3,2)], 'grey50'), 
                    breaks = c('1' ,'2', '3', '4'),
                    labels = c('No Signal', 'Shift', 'Decouple', 'Undershift')) +
  theme_classic() + labs(x = 'Year', y = 'Species', fill = 'Classification') +
  scale_y_continuous(breaks = unique(cbs_reorder$order), labels = unique(cbs_reorder$species), position = 'right') + 
  ggtitle('Classification') +
  theme(plot.title = element_text(size = 50, hjust = 0.5),
        axis.text.y.left = element_blank(),
        axis.line.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 25),
        axis.text.y.right = element_blank(), axis.title.y.right = element_blank(),
        legend.text = element_text(size = 20), legend.title = element_text(size = 25)) +
  theme(legend.position = 'left') + 
  geom_line(data = filter(cbs_reorder, (class_c == final) & (class != 'no signal')), aes(x = year, y = order, group = species), inherit.aes = F) +
  geom_point(data = filter(cbs_reorder, class != 'no signal'), aes(x = year_fc, y = order, group = species), inherit.aes = F) +
  geom_text(data = voltine_labs_p, aes(x = x, y = y, label = label), inherit.aes = F, angle = 90, size = 10) + 
  geom_hline(yintercept = mv_uv) 


# Plot class by species
plankton_cbs_plot_p = ggplot(planktys_sim_cbs_reorder, aes(x = year, y = order)) +
  geom_tile(aes(fill = as.character(class_c))) +
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
        axis.text.y.right = element_blank(), axis.title.y.right = element_blank(),
        legend.text = element_text(size = 20), legend.title = element_text(size = 25)) +
  theme(legend.position = 'left') + 
  geom_line(data = filter(planktys_sim_cbs_reorder, (class_c == final) & (class != 'no signal')), aes(x = year, y = order, group = species), inherit.aes = F) +
  geom_point(data = filter(planktys_sim_cbs_reorder, class != 'no signal'), aes(x = year_fc, y = order, group = species), inherit.aes = F) +
  geom_text(data = plankton_labs_p, aes(x = x, y = y, label = label), inherit.aes = F, angle = 90, size = 10) + 
  geom_hline(yintercept = phy_zoo) 



# Plot mismatch by species matrix
bfly_cmbs_plot_p = ggplot(bfly_cmbs_reorder_p, aes(x = year, y = order)) +
  geom_tile(aes(fill = mismatch)) +
  scale_fill_continuous(type = 'viridis') +
  theme_classic() + labs(x = 'Year', y = 'UKBMS', fill = 'Proportion\nMismatched') +
  scale_y_continuous(breaks = unique(bfly_cmbs_reorder$order), labels = unique(bfly_cmbs_reorder$sp1), position = 'right') + 
  ggtitle('Mismatch') +
  theme(plot.title = element_text(size = 50, hjust = 0.5),
        axis.line.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 25),
        axis.text.y.right = element_text(size = 20), axis.title.y.right = element_text(size = 50),
        legend.text = element_text(size = 20), legend.title = element_text(size = 25, margin = margin(0,0,1,0,'cm')), legend.position = 'left',
        legend.key.height= unit(3, 'cm'), legend.key.width = unit(2, 'cm')) + 
  geom_hline(yintercept = mv_uv) 



# Plot mismatch by species matrix
plankton_cmbs_plot_p = ggplot(planktys_cmbs_reorder_p, aes(x = year, y = order)) +
  geom_tile(aes(fill = mismatch)) +
  scale_fill_continuous(type = 'viridis') +
  theme_classic() + labs(x = 'Year', y = 'CPR', fill = 'Proportion\nMismatched') +
  scale_y_continuous(breaks = unique(planktys_sim_cbs_reorder_p$order), labels = unique(planktys_sim_cbs_reorder_p$species), position = 'right') + 
  theme(axis.text.y.left = element_blank(),
        axis.line.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 25),
        axis.text.y.right = element_text(size = 20), axis.title.y.right = element_text(size = 50),
        legend.text = element_text(size = 20), legend.title = element_text(size = 25, margin = margin(0,0,1,0,'cm')), legend.position = 'left',
        legend.key.height= unit(3, 'cm'), legend.key.width = unit(2, 'cm')) + 
  geom_hline(yintercept = phy_zoo) 



# Get legends
class_leg = get_legend(cbs_plot)
mm_leg = get_legend(bfly_cmbs_plot)

# export figure
tiff(filename = './Figures/cby_panelled.tif', width = 1800, height = 1800, units = 'px')

ggarrange(cbs_plot_p + theme(legend.position = 'none') + xlim(c(1979, 2019)), 
          bfly_cmbs_plot_p + theme(legend.position = 'none') + xlim(c(1979, 2019)), 
          mm_leg,
          plankton_cbs_plot_p + theme(legend.position = 'none') + xlim(c(1979, 2019)), 
          plankton_cmbs_plot_p + theme(legend.position = 'none') + xlim(c(1979, 2019)), 
          class_leg,
          nrow = 2, ncol = 3, heights = c(3.6, 1), widths = c(2,3,1))

# end figure
dev.off()

# Panel indicators

# Create new mismatch data frame to rbind to cby_comb
cmby_add = cmby_comb

# Rename columns
colnames(cmby_add) = c('year', 'proportion', 'name')

# Add emergence
cmby_add$emergence = 'n_mm'

# Add group
cmby_add$group = paste(cmby_add$emergence, cmby_add$name, sep = '_')

# Reorder
cmby_add = select(cmby_add, year, emergence, proportion, name, group)

# Combine
inds = rbind(cmby_add, cby_comb)


# Option 1: plot all indicators together
ind_plot = ggplot(inds, aes(x = year, y = proportion, group = group, color = emergence)) +
  geom_line(aes(linetype = name), linewidth = 3) + theme_classic() +
  scale_color_manual(values = c(cscale[c(3,4)], 'black'), labels = c('PDI', 'PSI', 'CDI')) +
  theme(legend.title = element_blank(), axis.text = element_text(size = 35), 
        axis.title = element_text(size = 35), legend.text = element_text(size = 35), 
        legend.key.width = unit(3, 'cm'), legend.position = 'inside', legend.position.inside = c(0.2, 0.8),
        plot.margin = unit(c(1,1,1,1), "cm")) +
  labs(y = 'Proportion of Emergences', x = 'Year') + ylim(c(0,1))


# export figure
tiff(filename = './Figures/inds_comb.tif', width = 1500, height = 1500, units = 'px')

ind_plot

# end figure
dev.off()


# Option 2: Panel all indicators together

# steal legend
ind_plot2 = ggplot(inds, aes(x = year, y = proportion, group = group)) +
  geom_line(aes(linetype = name), linewidth = 3, color = 'grey50') + theme_classic() +
  # scale_color_manual(values = c(cscale[c(4,3)], 'black'), labels = c('PSI', 'PDI', 'CDI')) +
  theme(plot.title = element_text(size = 50, hjust = 0.5),legend.title = element_blank(), 
        axis.text = element_text(size = 35), 
        axis.title = element_text(size = 35), legend.text = element_text(size = 35), 
        legend.key.width = unit(3, 'cm'), legend.position = 'inside', legend.position.inside = c(0.5, 0.5),
        plot.margin = unit(c(1,1,1,1), "cm")) +
  labs(y = 'Proportion of Emergences', x = 'Year') + ylim(c(0,1))


inds_leg = get_legend(ind_plot2)

# Plot PSI
psi_plot = ggplot(filter(inds, emergence == 'n_em'), aes(x = year, y = proportion, group = group)) +
  geom_line(aes(linetype = name), linewidth = 3, color = cscale[c(4)]) + theme_classic() +
  ggtitle('Phenological Shift Index (PSI)') +
  theme(plot.title = element_text(size = 35, hjust = 0.5, face = 'bold'),
        legend.title = element_blank(), axis.text = element_text(size = 35), 
        axis.title = element_text(size = 35), legend.text = element_text(size = 35), 
        legend.key.width = unit(3, 'cm'), legend.position = 'none', legend.position.inside = c(0.2, 0.8),
        plot.margin = unit(c(1,1,1,1), "cm")) +
  labs(y = 'Proportion of Emergences', x = 'Year') + ylim(c(0,1))



# Plot PDI
pdi_plot = ggplot(filter(inds, emergence == 'n_dc'), aes(x = year, y = proportion, group = group)) +
  geom_line(aes(linetype = name), linewidth = 3, color = cscale[c(3)]) + theme_classic() +
  ggtitle('Phenological Decoupling Index (PDI)') +
  theme(plot.title = element_text(size = 35, hjust = 0.5, face = 'bold'),
        legend.title = element_blank(), axis.text = element_text(size = 35), 
        axis.title = element_text(size = 35), legend.text = element_text(size = 35), 
        legend.key.width = unit(3, 'cm'), legend.position = 'none', legend.position.inside = c(0.2, 0.8),
        plot.margin = unit(c(1,1,1,1), "cm")) +
  labs(y = 'Proportion of Emergences', x = 'Year') + ylim(c(0,1))



# Plot CDI
cdi_plot = ggplot(filter(inds, emergence == 'n_mm'), aes(x = year, y = proportion, group = group)) +
  geom_line(aes(linetype = name), linewidth = 3, color = 'black') + theme_classic() +
  ggtitle('Community Desynchrony Index (CDI)') +
  theme(plot.title = element_text(size = 35, hjust = 0.5, face = 'bold'),
        legend.title = element_blank(), axis.text = element_text(size = 35), 
        axis.title = element_text(size = 35), legend.text = element_text(size = 35), 
        legend.key.width = unit(3, 'cm'), legend.position = 'none', legend.position.inside = c(0.2, 0.8),
        plot.margin = unit(c(1,1,1,1), "cm")) +
  labs(y = 'Proportion of Emergences', x = 'Year') + ylim(c(0,1))


# export figure
tiff(filename = './Figures/inds_sep.tif', width = 1000, height = 1500, units = 'px')

# # Arrange
# ggarrange(psi_plot, pdi_plot, cdi_plot, inds_leg, nrow = 2, ncol = 2)
ggarrange(psi_plot, ggplot() + theme_classic(), pdi_plot, inds_leg, cdi_plot, ggplot() + theme_classic(),
          nrow = 3, ncol = 2, widths = c(3.75, 1))
# ggarrange(psi_plot, pdi_plot, cdi_plot, 
#           nrow = 3, ncol = 1)

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













# SI tables and figures

# Match species names in cpr_params
cpr_params$species = slist$name_worms[match(cpr_params$species, slist$accepted_id)]

# Bind to butterfly params
plot_params = rbind(bfly_params, cpr_params)

# Save as output table for SI
write.csv(plot_params, './Figures/sim_params_data.csv')

# export figure
tiff(filename = './Figures/spread_vs_sigma.tif', width = 1500, height = 1000, units = 'px')

# Correlate sigma and s
ggplot(plot_params, aes(x = spread, y = mean_sd, color = dataset, fill = dataset)) + 
  geom_point(size = 3) + geom_smooth(method = 'lm', linewidth = 2, formula = y~x) +  theme_classic() + 
  labs(x = 's', y = expression(sigma)) + 
  theme(legend.title = element_blank(), axis.text = element_text(size = 35), 
        legend.text = element_text(size = 35), 
        legend.key.width = unit(3, 'cm'),
        axis.title = element_text(size = 35)) + 
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")), size = 10, formula = y ~ x, show.legend = F) + 
  stat_regline_equation(label.x = 15, label.y = c(7.5,5), size = 10, show.legend = F, formula = y ~ x)

# end figure
dev.off()



# Example figures

# Remove negative
ex_data = filter(ex_data, event > 0)

# Create detrended data container 
ex_dt = ex_data

# Model
ex_lm = lm(data = ex_data, event ~ year)

# Get slope
ex_b = coef(ex_lm)[2]

# Subtract out slope
ex_dt$event = ex_dt$event - ((ex_dt$year - min(ex_dt$year)) * ex_b)

# summarize
ex_dt_summ = group_by(ex_dt, year) %>% 
  summarize(mean_event = mean(event), sd_event = sd(event)) %>%
  mutate(lower = mean_event - sd_event, upper = mean_event + sd_event)

# Quantiles
ex_quants = quantile(ex_dt_summ$mean_event, c(0.25, 0.75))

# Pull range and create y limits
ranges = c(range(ex_data$event), range(ex_dt$event))
range = c(floor(min(ranges)), ceiling(max(ranges)))

# Pull range and create y limits for summarized set
ranges_summ = rbind(ex_dt_summ, ex_data_summ) %>% summarize(lower = min(lower), upper = max(upper))
range_summ = c(floor(min(ranges_summ)), ceiling(max(ranges_summ)))

# Plot 
ex_eplot = ggplot(ex_data, aes(x = year, y = event)) + geom_point(size = 5) +
  # geom_errorbar(inherit.aes = F, aes(x = year, ymin = lower, ymax = upper), linewidth = 2) +
  stat_smooth(method = 'lm', se = F, linewidth = 3) +
  theme_classic() + labs(x = 'Year', y = 'Event Timing (Julian Day)') + 
  theme(legend.title = element_blank(), axis.text = element_text(size = 35), 
        axis.title = element_text(size = 35)) + ylim(range)



# Plot 
ex_dplot = ggplot(ex_dt, aes(x = year, y = event)) + geom_point(size = 5) +
  # geom_errorbar(inherit.aes = F, aes(x = year, ymin = lower, ymax = upper), linewidth = 2, width = 0.2) +
  stat_smooth(method = 'lm', se = F, linewidth = 3) +
  theme_classic() + labs(x = 'Year', y = 'Event Timing (Julian Day)') + 
  theme(legend.title = element_blank(), axis.text = element_text(size = 35), 
        axis.title = element_text(size = 35)) + ylim(range)

# Plot 
ex_dplot_m = ggplot(ex_dt_summ, aes(x = year, y = mean_event)) + geom_point(size = 5) +
  geom_errorbar(inherit.aes = F, aes(x = year, ymin = lower, ymax = upper), linewidth = 2) +
  stat_smooth(method = 'lm', se = F, linewidth = 3) +
  geom_hline(yintercept = min(ex_quants), linetype = 'dashed', color = 'red', linewidth = 2) +
  theme_classic() + labs(x = 'Year', y = 'Event Timing (Julian Day)') + 
  theme(legend.title = element_blank(), axis.text = element_text(size = 35), 
        axis.title = element_text(size = 35)) + ylim(range_summ)

# Emergence block
ex_tope = emp_tope(ex_data, plot = F)
ex_em_block = data.frame(xmin = min(ex_tope[which(ex_tope$emerged == 1), 'year']), xmax = max(ex_tope[which(ex_tope$emerged == 1), 'year']))


# Plot test results
ex_eplot_em = ggplot(ex_data_summ, aes(x = year, y = mean_event)) + 
  geom_errorbar(inherit.aes = F, aes(x = year, ymin = lower, ymax = upper), linewidth = 2) +
  geom_point(size = 5, aes(color = ifelse(mean_event < min(ex_quants), 'red', 'black'))) +
  stat_smooth(method = 'lm', se = F, linewidth = 3) +
  geom_rect(inherit.aes = F, data = ex_em_block, aes(ymin = -Inf, ymax = Inf, xmin = xmin, xmax = Inf), fill = 'red', alpha = 0.3) +
  scale_color_manual(values = c('black', 'red'), name = 'Test Result', labels = c('Negative', 'Positive')) +
  geom_hline(yintercept = min(ex_quants), linetype = 'dashed', color = 'red', linewidth = 2) +
  theme_classic() + labs(x = 'Year', y = 'Event Timing (Julian Day)') + 
  theme(legend.title = element_text(size = 30), axis.text = element_text(size = 35), 
        legend.text = element_text(size = 30), legend.background = element_blank(),
        legend.position = 'inside', legend.position.inside = c(0.14,0.14), 
        axis.title = element_text(size = 35)) + ylim(range_summ)

# Export figures
tiff(filename = './Figures/SI/emp_concept.tif', width = 1500, height = 1000, units = 'px')

# Arrange into empirical test conceptual figure
ggarrange(ex_eplot, ex_dplot, ex_eplot_em, ex_dplot_m, nrow = 2, ncol = 2)

# End figure
dev.off()


# Export figures
tiff(filename = './Figures/SI/original_data.tif', width = 750, height = 500, units = 'px')

ex_eplot

# End figure
dev.off()


# Export figures
tiff(filename = './Figures/SI/detrended_data.tif', width = 750, height = 500, units = 'px')

ex_dplot

# End figure
dev.off()


# Export figures
tiff(filename = './Figures/SI/detrended_quantile.tif', width = 750, height = 500, units = 'px')

ex_dplot_m

# End figure
dev.off()


# Export figures
tiff(filename = './Figures/SI/original_tests.tif', width = 750, height = 500, units = 'px')

ex_eplot_em

# End figure
dev.off()




# Pull a single year of data
tyear = 1990
ex_yr = filter(ex_dt, year == tyear) %>% select(species, year, event, env)
ex_yr_o = filter(ex_data, year == tyear) %>% select(species, year, event, env)


# Container
ex_dtboot = NULL

# Run sampling
for(i in 1:100){
  
  # Calculate sample
  ex_boot = sample(ex_yr$event, nrow(ex_yr), replace = T)
  
  # Add identifier
  ex_boot = data.frame(event = ex_boot, sample = i)
  
  # Run KS test
  ex_boot$p = ks.test(ex_yr_o$event, ex_boot$event, alternative = 'greater')$p.value
  
  # Concatenate
  ex_dtboot = rbind(ex_dtboot, ex_boot)
  
}

# Check for significance
ex_dtboot$sig = ifelse(ex_dtboot$p < 0.05, 1, 0)

# # Boxplot bootstrap samples
# ggplot(ex_dtboot, aes(y = event, x = sample, group = sample)) + geom_boxplot() + 
#   theme_classic() + labs(x = 'Sample', y = 'Event Timing (Julian Day)') + 
#   theme(legend.title = element_blank(), axis.text = element_text(size = 35), 
#         axis.title = element_text(size = 35))

# Add real data to ex_dtboot
ex_real = filter(ex_data, year == tyear) %>% mutate(sample = 0, p = 0, sig = 2) %>% select(event, sample, p, sig)
ex_dtboot = rbind(ex_real, ex_dtboot) %>%
  mutate(real = (ifelse(sample == 0, 'original', 'detrended bootstrap')))


# ECDF Bootstrap samples
ex_ecdf = ggplot(ex_dtboot, aes(x = event, group = sample, color = as.factor(sig), alpha = real, linewidth = real)) + stat_ecdf() +
  scale_alpha_manual(values = c(0.1,1), guide = 'none') +
  scale_linewidth_manual(values = c(0.5, 2), guide = 'none') +
  # stat_ecdf(data = filter(ex_data, year == 2015), aes(x = event), inherit.aes = F, color = 'red', linewidth = 2) +   
  guides(color = guide_legend(reverse = TRUE, nrow = 3, byrow = TRUE, override.aes = list(linewidth = c(2,1,1), alpha = c(1,0.2,0.2)))) +
  theme_classic() + labs(y = 'ECDF', x = 'Event Timing (Julian Day)') + ggtitle('Year = 1990') +
  theme(legend.title = element_blank(), axis.text = element_text(size = 35), 
        plot.title = element_text(hjust = 0.5, size = 30),
        legend.text = element_text(size = 30), legend.position = 'inside',
        legend.position.inside = c(0.3, 0.75), legend.spacing.y = unit(6, 'cm'),
        legend.key.width = unit(3, 'cm'), 
        axis.title = element_text(size = 35)) + 
  scale_color_manual(values = c('blue', 'red', 'black'), 
                       labels = c('Detrended Bootstrap\nKS Test p > 0.05', 'Detrended Bootstrap\nKS Test p < 0.05', 'Original'))
# Calculate boostrap means
ex_dtboot_hist = filter(ex_dtboot, real == 'detrended bootstrap') %>% group_by(sample) %>% summarize(bs_mean = mean(event))

# Plot bootstrap means
ex_bs_hist = ggplot(ex_dtboot_hist) + geom_histogram(color = 'black', aes(x = bs_mean), alpha = 0.5) +
  geom_vline(xintercept = mean(ex_yr$event), linewidth = 3, color = 'blue') +
  # geom_vline(xintercept = mean(ex_yr_o$event), linewidth = 3, color = 'black') +
  theme_classic() + scale_y_continuous(limits = c(0, NA), expand = c(0,0)) +
  labs(x = 'Bootstrap Mean Event Timing', y = 'Count') + ggtitle('Year = 1990') +
  theme(legend.title = element_blank(), axis.text = element_text(size = 35), 
        plot.title = element_text(hjust = 0.5, size = 30),
        legend.text = element_text(size = 35), legend.position = 'inside',
        legend.position.inside = c(0.3, 0.75),
        legend.key.width = unit(3, 'cm'),
        axis.title = element_text(size = 35))
  

# Calculate proportion of significant KS tests
ex_ks_tope = ks_tope(ex_data)
ex_ks_tope$test = ifelse(ex_ks_tope$p >= 0.6, 1, 0)

# Emergence block
ex_ks_block = data.frame(xmin = min(ex_ks_tope[which(ex_ks_tope$emerged == 1), 'year']), xmax = max(ex_ks_tope[which(ex_ks_tope$emerged == 1), 'year']))

# Plot proportion of significant KS tests
ex_ks_testres = ggplot(ex_ks_tope, aes(x = year, y = p)) + 
  geom_hline(yintercept = 0.6, color = 'red', linetype = 'dashed', linewidth = 2) +
  geom_line(linewidth = 2) + geom_point(aes(color = as.factor(test)), size = 5) + 
  geom_rect(inherit.aes = F, data = ex_ks_block, aes(ymin = -Inf, ymax = Inf, xmin = xmin, xmax = Inf), fill = 'red', alpha = 0.3) +
  scale_color_manual(values = c('black', 'red'), name = 'Test Result', labels = c('Negative', 'Positive')) +   
  theme_classic() + labs(x = 'Year', y = 'Proportion of KS Tests\nSignificant') + 
  theme(legend.title = element_text(size = 35), axis.text = element_text(size = 35), 
        legend.text = element_text(size = 35), legend.position = 'bottom',
        legend.key.width = unit(3, 'cm'),
        axis.title = element_text(size = 35))

# Export figures
tiff(filename = './Figures/SI/ecdf.tif', width = 750, height = 500, units = 'px')

ex_ecdf

# End figure
dev.off()


# Export figures
tiff(filename = './Figures/SI/boostrap_means.tif', width = 750, height = 500, units = 'px')

ex_bs_hist

# End figure
dev.off()


# Export figures
tiff(filename = './Figures/SI/stat_tests.tif', width = 900, height = 500, units = 'px')

ex_ks_testres

# End figure
dev.off()
