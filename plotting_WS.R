### Code to generate plots for Bronnvik et al., 2024
### Hester Bronnvik
### 2023-08-29
### hbronnvik@ab.mpg.de

library(cowplot)
library(tidyverse)

# the processed data used to look at relationships within the data
a_data <- readRDS("/home/hbronnvik/Documents/storkSSFs/a_data_2023-07-26.rds")%>% 
  mutate(season = ifelse(grepl("fall", track), "post", "pre"))

# the model results and predictions used to visualize results
pre_mod <- readRDS("/home/hbronnvik/Documents/storkSSFs/glmm_TMB_M_pre_2023-08-28.rds")
post_mod <- readRDS("/home/hbronnvik/Documents/storkSSFs/glmm_TMB_M_post_2023-08-28.rds")

pre_preds <- readRDS("/home/hbronnvik/Documents/storkSSFs/glmm_preds_pre_2023-08-28.rds")
post_preds <- readRDS("/home/hbronnvik/Documents/storkSSFs/glmm_preds_post_2023-08-28.rds")

# potential color schemes
# pansy, wisteria, apricot, tangerine
colfunc <- colorRampPalette(c("#9672D5", "#C5B0E8", "#FDC4AF", "#FB8F67"))
# flame, carrot, xanthous, yellow green, apple, avocado
colfunc <- colorRampPalette(c("#6A8532", "#87A330", "#A1C349", "#F3C053", "#F9A03F", "#EB5E28"))

# a universal palette (cool to warm)
# cerulean, Munsell, verdigris, Tiffany, light orange, melon, salmon, bittersweet
colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))
# universal facet labels
fac_labs <- c("Fall", "Spring")
names(fac_labs) <- c("post", "pre")

### plots for data exploration

# look at uplift distributions
png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/uplift_dists.png",
    height = 18, width = 20, units = "in", res = 300)
ggplot(a_data, aes(x = as.factor(migrations), y = w_star, fill = as.factor(used))) +
  geom_boxplot() +
  scale_fill_manual(values = c("#0096CC", "#FE6D5D")) +
  labs(x = "Migrations", y = "Uplift (m/s)", fill = "Use") +
  theme_classic() +
  facet_wrap(~season, labeller = labeller(season = fac_labs))
dev.off()

# look at wind support distributions
support <- ggplot(a_data, aes(as.factor(migrations), wind_support, fill = as.logical(used))) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot(linewidth = 1.5) +
  theme_classic() +
  labs(x = "Migrations", y = "Wind support (m/s)", fill = "Use case") +
  scale_color_manual(values = c("#0081A7", "#F07167")) +
  theme_classic() +
  theme(text = element_text(size = 25),
        axis.line = element_line(color = "black", linewidth = 1.5),
        axis.text = element_text(color = "black"),
        legend.key.size = unit(2, 'cm')) +
  facet_wrap(~season, labeller = labeller(season = fac_labs))
# look at wind speed distributions
speeds <- ggplot(a_data, aes(as.factor(migrations), wind_speed, fill = as.logical(used))) + 
  geom_boxplot(linewidth = 1.5) +
  theme_classic() +
  labs(x = "Migrations", y = "Wind speed (m/s)", fill = "Use case") +
  scale_color_manual(values = c("#0081A7", "#F07167")) +
  theme_classic() +
  theme(text = element_text(size = 25),
        axis.line = element_line(color = "black", linewidth = 1.5),
        axis.text = element_text(color = "black"),
        legend.key.size = unit(2, 'cm')) +
  facet_wrap(~season, labeller = labeller(season = fac_labs))

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/wind.png",
    height = 18, width = 20, units = "in", res = 300)
plot_grid(support, speeds, labels = c("A", "B"), ncol = 1,
          align = 'v', axis = 'l')
dev.off()

a_data$datestamp <- a_data$timestamp
year(a_data$datestamp) <- 2024

# look at the number of recordings of migration on each day for each age group
png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/migration_timing_sp23.png",
    height = 18, width = 20, units = "in", res = 300)
a_data %>% 
  mutate(yd = date(alignment)) %>% 
  group_by(yd, migrations) %>% 
  mutate(count = n()) %>% 
  slice(1) %>%
  ungroup() %>% 
  arrange(migrations) %>% 
  ggplot(aes(x = yd, y = count, group = as.factor(migrations), color = as.factor(migrations))) +
  geom_segment(aes(alpha = 0.7, x=yd, xend=yd, y=0, yend=count), linewidth = 1.5) +
  geom_point(size = 3) +
  labs(x = "Day", y = "Observations", color = "Migration") +
  scale_color_manual(values = c("#0081A7", "#0098B0", "#00AFB9", "#7FC4B8", "#F7A58F", "#F27E71", "#F07167", "#ED5145", "#EE5E53")) +
  guides(alpha = "none") +
  scale_x_date(date_breaks = "1 month", date_labels =  "%m-%d") +
  theme_classic() +
  theme(text = element_text(size = 25),
        axis.line = element_line(color = "black", linewidth = 1.5),
        axis.text = element_text(color = "black"),
        legend.key.size = unit(2, 'cm'))
dev.off()

# look at uplift distributions per age group
png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/seasonal_uplifts.png",
    height = 18, width = 20, units = "in", res = 300)
ggplot(a_data %>% mutate(season = ifelse(season == "pre", "Spring", "Fall")), aes(as.factor(migrations), w_star, fill = season)) +
  geom_boxplot(lwd = 1.5, aes(alpha = forcats::fct_rev(as.factor(migrations)))) +
  scale_fill_manual(values = c("#FE6D5D", "#0096CC"))  +
  guides(alpha = "none") +
  labs(x = "Migrations", y = "Uplift (m/s)", fill = "Season") +
  theme_classic() +
  theme(text = element_text(size = 25),
        axis.line = element_line(color = "black", linewidth = 1.5),
        axis.text = element_text(color = "black"))
dev.off()

# look at conspecific density distributions per age group
png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/density_dists.png",
    height = 18, width = 20, units = "in", res = 300)
ggplot(a_data, aes(x = as.factor(migrations), y = ud_pdf, fill = as.factor(used))) +
  geom_boxplot(lwd = 1.5, aes(alpha = forcats::fct_rev(as.factor(migrations)))) +
  scale_fill_manual(values = c("#0096CC", "#FE6D5D")) + 
  guides(alpha = "none") +
  labs(x = "Migrations", y = "Conspecific density", fill = "Use") +
  theme_classic() +
  theme(text = element_text(size = 25),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black")) +
  facet_wrap(~season, labeller = labeller(season = fac_labs))
dev.off()

# look at uplift distributions per season
res_aov <- aov(w_star ~ factor(season), data = a_data)
hist(res_aov$residuals)
# not at all normal
car::qqPlot(res_aov$residuals)
wilcox.test(w_star ~ factor(season), data = a_data)
png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/seasons_uplift.png",
    height = 18, width = 20, units = "in", res = 300)
print(ggplot(a_data %>% mutate(season = ifelse(season == "post", "Fall", "Spring")), aes(season, w_star, fill = season)) +
        geom_boxplot(linewidth = 2) +
        scale_color_manual(values = c("#007BA7", "#FE6F5E")) +
        labs(x = "Season", y = "Uplift potential (m/s)", fill = "Season") +
        theme_classic() +
        theme(axis.text = element_text(color = "black"),
              axis.line = element_line(linewidth = 1.5),
              text = element_text(size = 25)))
dev.off()

# look at overall distributions
fac_labs <- c("Fall", "Spring")
names(fac_labs) <- c("post", "pre")
w_star_var <- ggplot(a_data, aes(as.factor(migrations), w_star, group = as.factor(migrations), 
                                 fill = season, alpha = forcats::fct_rev(as.factor(migrations))))+
  geom_boxplot() +
  scale_fill_manual(values = c("#FE6F5E", "#007BA7")) +
  scale_alpha_discrete(
    range = c(0.3, 0.9),
    guide = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "", y = "Uplift potential (m/s)", fill = "Season") +
  theme_classic() +
  theme(text = element_text(size = 20, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        legend.position = "none") +
  facet_wrap(~season, labeller = labeller(season = fac_labs))
dens_var <- ggplot(a_data, aes(as.factor(migrations), ud_pdf, group = as.factor(migrations), 
                               fill = season, alpha = forcats::fct_rev(as.factor(migrations)))) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FE6F5E", "#007BA7")) +
  scale_alpha_discrete(
    range = c(0.3, 0.9),
    guide = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Migration", y = "Conspecific density", fill = "Season") +
  theme_classic() +
  theme(text = element_text(size = 20, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        legend.position = "none") +
  facet_wrap(~season, labeller = labeller(season = fac_labs))

vars_plot <-  plot_grid(w_star_var, dens_var, labels = c("A", "B"), ncol = 1,
                        align = 'v', axis = 'l')

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/distributions.png",
    height = 200, width = 300, units = "mm", res = 500)
vars_plot
dev.off()

# look at the distributions of variances

k_data <- a_data %>% 
  group_by(stratum) %>% 
  mutate(obs_diff = w_star[which(used == 1)]-mean(w_star[which(used == 0)]),
         strat_var_w = var(w_star),
         strat_var_s = var(ud_pdf)) %>% 
  slice(1) %>% 
  ungroup()
w_dens <- ggplot(k_data, aes(log(strat_var_w), fill = as.factor(migrations))) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = colfunc(9)) +
  labs(x = "log per-hour variance in uplift", y = "Density", fill = "Migrations") +
  theme_classic()  +
  theme(axis.text = element_text(color = "black", size = 35),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 45),
        legend.text = element_text(size = 35),
        strip.text.x = element_text(size = 30),
        legend.key.size = unit(1.5, 'cm')) +
  facet_wrap(~season, labeller = labeller(season = fac_labs))
s_dens <- ggplot(k_data, aes(log(strat_var_s), fill = as.factor(migrations))) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = colfunc(9)) +
  labs(x = "log per-hour variance in conspecific density", y = "Density", fill = "Migrations") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 35),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 45),
        legend.text = element_text(size = 35),
        strip.text.x = element_text(size = 30),
        legend.key.size = unit(1.5, 'cm')) +
  facet_wrap(~season, labeller = labeller(season = fac_labs))
png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/stratum_densities.png",
    height = 18, width = 20, units = "in", res = 300)
plot_grid(s_dens, w_dens, ncol = 1, align = 'v', axis = 'l')
dev.off()


### plots for the models

# look at model predictions for the 3-term interaction in spring
my_labels <- c("Low", "Mean", "High")
png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/preds_glmm_spring_3_09-04.png",
    height = 18, width = 20, units = "in", res = 300)
ggplot(pre_preds %>%  
         filter(interaction == "uplift_migration_ud") %>% 
         mutate(label = paste0("Spring ", migrations))) +
  geom_raster(aes(x = w_star, y = sqrt_ud, fill = probs, group = probs), interpolate = F) +
  # geom_contour(aes(x = z, y = y, z = value), color = "black") +
  scale_x_continuous(expand = c(0, 0), n.breaks = 4) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0,8.867e-07,1.755e-06), labels = my_labels) +
  labs(x = "Uplift (m/s)", y = "Conspecific density") +
  scale_fill_gradientn("Selection \nprobability", colours = colfunc(135)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 35),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 45),
        legend.text = element_text(size = 35),
        strip.text.x = element_text(size = 30),
        legend.key.size = unit(1.5, 'cm')) +
  facet_wrap(~as.factor(label))
dev.off()

# look at model predictions for the 3-term interaction in fall
png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/preds_glmm_fall_3_09-04.png",
    height = 18, width = 20, units = "in", res = 300)
ggplot(post_preds %>%  
         filter(interaction == "uplift_migration_ud") %>% 
         mutate(label = paste0("Fall ", migrations))) +
  geom_raster(aes(x = w_star, y = sqrt_ud, fill = probs, group = probs), interpolate = F) +
  # geom_contour(aes(x = z, y = y, z = value), color = "black") +
  scale_x_continuous(expand = c(0, 0), n.breaks = 4) +
  scale_y_continuous(expand = c(0, 0), breaks = c(1.489e-07,1.121e-06,2.092e-06), labels = my_labels) +
  labs(x = "Uplift (m/s)", y = "Conspecific density") +
  scale_fill_gradientn("Selection \nprobability", colours = colfunc(10)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 35),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 45),
        legend.text = element_text(size = 35),
        strip.text.x = element_text(size = 30),
        legend.key.size = unit(1.5, 'cm')) +
  facet_wrap(~as.factor(label))
dev.off()

# produce expectations
ex <- expand.grid(x = 1:2, y = 1:2)
ex$probs <- ifelse(ex$x == 2 | ex$y == 2, 0.5, 1)
ex$probs <- ifelse(ex$y == 2 & ex$x == 2, 0, ex$probs)

  yPosxNon <-
  ggplot(ex) +
  geom_raster(aes(x = x, y = y, fill = probs, group = probs), interpolate = F) +
  geom_hline(yintercept = 1.5, lwd = 1.5) +
  geom_vline(xintercept = 1.5, lwd = 1.5) +
  scale_x_continuous(expand = c(0, 0), n.breaks = 9) +
  scale_y_continuous(expand = c(0, 0), n.breaks = 6) +
  labs(x = "X", y = "Y", title = "Selection for Y and against X") +
  scale_fill_gradientn("Selection \nprobability", colours = colfunc(135)) +
  theme_classic() +
  theme(axis.text = element_text(color = "white", size = 0),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 45),
        legend.text = element_text(size = 35),
        strip.text.x = element_text(size = 30),
        legend.key.size = unit(1.5, 'cm'),
        axis.ticks = element_blank(),
        legend.position = "none")
  
  xyNeg <- xyNeg + theme(text = element_text(size = 25))

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/preds_template.png",
    height = 18, width = 20, units = "in", res = 600)
cowplot::plot_grid(xyNon, xyPos, yPos, yPosxNon, xNeg, xyNeg, yNeg, xPosyNon, xPos, ncol = 3,
          align = 'v', axis = 'l', label_size = 10)
dev.off()

# look at the 3-term interaction in 3D
i_data <- pre_preds %>%  
  filter(interaction == "uplift_migration_ud") %>% 
  dplyr::select("migrations", "sqrt_ud", "w_star", "probs")
plot3D::scatter3D(x = i_data$migrations, z = i_data$w_star, y = i_data$sqrt_ud, 
                  colvar = i_data$probs, col = colfunc(135), pch = 16, alpha = .4, 
                  cex = 1.5, xlab = "Migrations", zlab = "Uplift (m/s)", ylab = "Conspecific density",
                  ticktype = "detailed")

# look at model predictions for the 2-term interactions with the 3rd term held at its mean
inter_plots <- function(data, x_id, y_id, xlab, ylab, xmin, xmax){
  ggplot(data, aes_string(x_id, y_id, fill = "probs")) +
    geom_tile(color = "white", lwd = 0, linetype = 1) +
    scale_fill_gradientn(colors = colfunc(135)) +
    scale_y_continuous(expand=c(0, 0))+
    scale_x_continuous(expand=c(0, 0),
                       breaks=round(seq(from = xmin, to = xmax, length.out = 9), digits = 1))+
    labs(x = xlab, y = ylab, fill = "Selection probability") +
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 15),
          axis.line = element_line(linewidth = 1.2),
          text = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position = "none")
}

# spring interactions
w_age_s <- inter_plots(pre_preds %>% filter(interaction == "uplift_migration"),
                       "migrations", "w_star",
                       "Number of spring migrations", "Uplift (m/s)", 
                       1, 9) + labs(title = " ")
ud_age_s <- inter_plots(pre_preds %>% filter(interaction == "ud_migration"),
                        "migrations", "sqrt_ud",
                        "Number of spring migrations", "Social density", 
                        1, 9) + labs(title = "Spring")
ud_up_s <- inter_plots(pre_preds %>% filter(interaction == "ud_up"),
                       "w_star", "sqrt_ud",
                       "Uplift (m/s)", "Social density", 
                       -1.9, 3.6) + labs(title = " ")
# fall interactions
w_age_f <- inter_plots(post_preds %>% filter(interaction == "uplift_migration"),
                       "migrations", "w_star",
                       "Number of fall migrations", "Uplift (m/s)", 
                       1, 9) + labs(title = " ")
ud_age_f <- inter_plots(post_preds %>% filter(interaction == "ud_migration"),
                        "migrations", "sqrt_ud",
                        "Number of fall migrations", "Social density", 
                        1, 9) + labs(title = "Fall")
ud_up_f <- inter_plots(post_preds %>% filter(interaction == "ud_up"),
                       "w_star", "sqrt_ud",
                       "Uplift (m/s)", "Social density", 
                       -2.3, 3.9) + labs(title = " ")
# combine the plots
blank <- ggplot(post_preds %>% filter(interaction == "ud_up"), aes(w_star, sqrt_ud, fill = probs)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  scale_fill_gradientn(colors = colfunc(135)) +
  labs(x = "", y = "", fill = "Selection probability") +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 15),
        plot.margin = margin(c(0,1,0,0)))

legend <- get_legend(
  # create some space to the left of the legend
  blank + theme(legend.box.margin = margin(0, 0, 0, 10))
)

top_row <- plot_grid(ud_age_f, w_age_f, ud_up_f, labels = c("1", "2", "3"), 
                     label_x = 0.1, label_y = 0.9,
                     ncol = 4,align = 'v', axis = 'l')
bottom_row <- plot_grid(ud_age_s, w_age_s, ud_up_s, labels = c("1", "2", "3"), 
                        label_x = 0.1, label_y = 0.9, 
                        ncol = 4, align = 'v', axis = 'l')
full_plot <-  plot_grid(top_row, bottom_row, labels = c("A", "B"), ncol = 1,
                        align = 'v', axis = 'l')
full_plot + draw_grob(legend, scale = 0.7, x = 0.4)

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/preds_glmm_sqrt_ud_08-23.png",
    height = 18, width = 20, units = "in", res = 300)
print(full_plot + draw_grob(legend, scale = 0.7, x = 0.4))
dev.off()

library(rayshader)
i_data <- pre_preds %>%  
  filter(interaction == "uplift_migration_ud") %>% 
  dplyr::select("migrations", "sqrt_ud", "w_star", "probs") %>% 
  rename(x = migrations,
         y = sqrt_ud,
         z = w_star,
         value = probs)
heatmap_fill <- ggplot(i_data %>% group_by(x, y) %>% slice_sample(n = 1) %>% ungroup()) +
  geom_tile(aes(x, y, fill = value, color = z)) +
  scale_fill_gradientn(colors = colfunc(135)) +
  theme_classic() +
  theme(legend.position = "none")
heatmap_height <- ggplot(i_data %>% group_by(x, y) %>% slice_sample(n = 1) %>% ungroup()) +
  geom_tile(aes(x, y, fill = z, color = z)) +
  scale_fill_gradientn(colors = colfunc(135)) +
  theme_classic() +
  theme(legend.position = "none")
# rayshader::plot_gg(plot, width = 7, height = 4, raytrace = FALSE, preview = TRUE)
rayshader::plot_gg(heatmap_fill, heatmap_height, multicore = TRUE, raytrace = TRUE, width = 7, height = 4, 
                   scale = 300, windowsize = c(1400, 866), zoom = 0.6, phi = 30, theta = 30);
render_snapshot(clear = TRUE)

# look at the actual model results
pre_graph <- confint(pre_mod) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Factor") %>% 
  filter(!grepl("id", Factor)) %>% 
  mutate(Variable = c("Conspecific density", "Uplift", "Spring migrations", #"Step length", "Turning angle", 
                      "Conspecifics X Uplift", "Conspecifics X Migrations", "Uplift X Migrations", "Conspecifics X Uplift X Migrations"))
colnames(pre_graph)[2:3] <- c("Lower", "Upper") 

pre_graph$p <- as.data.frame(summary(pre_mod)[[6]][1])[,4]
pre_graph$significance <- ifelse(pre_graph$p < 0.001, "***", 
                                 ifelse(pre_graph$p > 0.001 & pre_graph$p < 0.01, "**",
                                        ifelse(pre_graph$p > 0.01 & pre_graph$p < 0.05, "*",
                                               "")))
pre_graph <- pre_graph %>% 
  mutate(Variable = factor(Variable, levels = c("Conspecifics X Uplift X Migrations", "Uplift X Migrations",
                                                "Conspecifics X Uplift", "Conspecifics X Migrations", "Spring migrations", 
                                                "Uplift", "Conspecific density"))) %>% 
  arrange(Variable)

post_graph <- confint(post_mod) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Factor") %>% 
  filter(!grepl("id", Factor)) %>% 
  mutate(Variable = c("Conspecific density", "Uplift", "Fall migrations", #"Step length", "Turning angle", 
                      "Conspecifics X Uplift", "Conspecifics X Migrations", "Uplift X Migrations", "Conspecifics X Uplift X Migrations"))
colnames(post_graph)[2:3] <- c("Lower", "Upper") 

post_graph$p <- as.data.frame(summary(post_mod)[[6]][1])[,4]
post_graph$significance <- ifelse(post_graph$p < 0.001, "***", 
                                  ifelse(post_graph$p > 0.001 & post_graph$p < 0.01, "**",
                                         ifelse(post_graph$p > 0.01 & post_graph$p < 0.05, "*",
                                                "")))
post_graph <- post_graph %>% 
  mutate(Variable = factor(Variable, levels = c("Conspecifics X Uplift X Migrations", "Uplift X Migrations",
                                                "Conspecifics X Uplift", "Conspecifics X Migrations", "Fall migrations", 
                                                "Uplift", "Conspecific density"))) %>% 
  arrange(Variable)

post_coefs <- ggplot(post_graph %>% filter(Variable != "Fall migrations"), aes(Estimate, Variable)) +
  geom_vline(xintercept = 0, lty = 2, linewidth = 1) +
  geom_pointrange(aes(xmin = Lower, xmax = Upper), #fill = significance_05, color = significance_05), 
                  linewidth = 1, size = 1) +
  labs(y = "", fill = "Significant", color = "Significant", title = "Fall") +
  theme_classic() +
  ggplot2::annotate(geom="text", x=16, y=c(1:6), label=post_graph$significance[post_graph$Variable != "Fall migrations"], color="black", size = 6)  +
  theme(axis.text = element_text(color = "black", size = 18),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 25),
        legend.text = element_text(size = 35),
        legend.key.size = unit(1.5, 'cm')) 
pre_coefs <- ggplot(pre_graph %>% filter(Variable != "Spring migrations"), aes(Estimate, Variable)) +
  geom_vline(xintercept = 0, lty = 2, linewidth = 1) +
  geom_pointrange(aes(xmin = Lower, xmax = Upper), #fill = significance_05, color = significance_05), 
                  linewidth = 1, size = 1) +
  labs(y = "", fill = "Significant", color = "Significant", title = "Spring") +
  theme_classic() +
  ggplot2::annotate(geom="text", x=5, y=c(1:6), label=pre_graph$significance[pre_graph$Variable != "Spring migrations"], color="black", size = 6)  +
  theme(axis.text = element_text(color = "black", size = 18),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 25),
        legend.text = element_text(size = 35),
        legend.key.size = unit(1.5, 'cm'))

coefs_plot <-  plot_grid(post_coefs, pre_coefs, labels = c("A", "B"), nrow = 2,
                         align = 'v', axis = 'l')
coefs_plot

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/coefs_glmm_05-09.png",
    height = 200, width = 300, units = "mm", res = 500)
coefs_plot
dev.off()


