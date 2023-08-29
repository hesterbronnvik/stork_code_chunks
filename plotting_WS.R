### Code to generate plots for Bronnvik et al., 2024
### Hester Bronnvik
### 2023-08-29
### hbronnvik@ab.mpg.de

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


fac_labs <- c("Fall", "Spring")
names(fac_labs) <- c("post", "pre")
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

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/migration_timing.png",
    height = 18, width = 20, units = "in", res = 300)
a_data %>% 
  mutate(yd = date(datestamp)) %>% 
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


### plots for the models

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/preds_glmm_spring_3_08-26.png",
    height = 18, width = 20, units = "in", res = 300)
ggplot(pre_preds %>%  
         filter(interaction == "uplift_migration_ud") %>% 
         mutate(label = paste0("Spring ", migrations))) +
  geom_raster(aes(x = w_star, y = sqrt_ud, fill = probs, group = probs), interpolate = F) +
  # geom_contour(aes(x = z, y = y, z = value), color = "black") +
  scale_x_continuous(expand = c(0, 0), n.breaks = 11) +
  scale_y_continuous(expand = c(0, 0), n.breaks = 6) +
  labs(x = "Uplift (m/s)", y = "Conspecific density") +
  scale_fill_gradientn("Selection", colours = colfunc(135)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 20),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 25),
        legend.text = element_text(size = 25),
        strip.text.x = element_text(size = 30)) +
  facet_wrap(~as.factor(label))
dev.off()

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/preds_glmm_fall_3_08-26.png",
    height = 18, width = 20, units = "in", res = 300)
ggplot(post_preds %>%  
         filter(interaction == "uplift_migration_ud") %>% 
         mutate(label = paste0("Fall ", migrations))) +
  geom_raster(aes(x = w_star, y = sqrt_ud, fill = probs, group = probs), interpolate = F) +
  # geom_contour(aes(x = z, y = y, z = value), color = "black") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Uplift (m/s)", y = "Conspecific density") +
  scale_fill_gradientn("Selection", colours = colfunc(10)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 20),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 25),
        legend.text = element_text(size = 25),
        strip.text.x = element_text(size = 30)) +
  facet_wrap(~as.factor(label))
dev.off()

inter_plots <- function(data, x_id, y_id, xlab, ylab, xmin, xmax){
  ggplot(data, aes_string(x_id, y_id, fill = "probs")) +
    geom_tile(color = "white", lwd = 1.5, linetype = 1) +
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

library(plotly)
plot_ly(z = as.matrix(a_data[, c("migrations", "w_star", "sqrt_ud_z")])) %>% add_surface()

# interpolate data onto grid
test_plotly <- with(a_data[a_data$season == "post",], akima::interp(migrations, w_star, sqrt_ud,
                                                                    duplicate = "mean"))
# plot surface over grid
axx <- list(title = "Migrations")
axy <- list(title = "Uplift (m/s)")
axz <- list(title = "Conspecific density")
plot_ly(x = test_plotly$x, y = test_plotly$y, z = test_plotly$z,
        type = "surface") %>% 
  layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))

p <- inter_preds %>%  
  filter(interaction == "uplift_migration_ud") %>% 
  ggplot(aes(x = "migrations", y = "w_star", z = "sqrt_ud")) +  
  geom_contour() +  
  scale_fill_distiller(palette = "Spectral", direction = -1)
ggplotly(p)

grd_3 <- expand.grid(x = (1:max(prep_d$migrations)),
                     y = seq(from = min(prep_d$w_star, na.rm = T), to = max(prep_d$w_star, na.rm = T), length.out = 15),
                     z = seq(from = min(prep_d$sqrt_ud, na.rm = T), to = max(prep_d$sqrt_ud, na.rm = T), length.out = 15)) %>% # quantile(prep_d$w_star, .9, na.rm = T)
  rename(migrations = x,
         w_star = y,
         sqrt_ud = z) %>% 
  mutate(step_length = mean(prep_d$step_length, na.rm = T),
         turning_angle = mean(prep_d$turning_angle, na.rm = T),
         interaction = "uplift_migration_ud")

plot3D::scatter3D(x = i_data$migrations, y = i_data$w_star, z = i_data$sqrt_ud, colvar = i_data$probs)
rgl::surface3d(x = i_data$migrations, y = i_data$w_star, z = i_data$sqrt_ud)

i_data <- inter_preds %>%  
  filter(interaction == "uplift_migration_ud") %>% 
  dplyr::select("migrations", "sqrt_ud", "w_star", "probs") %>% 
  rename(x = "migrations", 
         y = "sqrt_ud",
         z = "w_star",
         value = "probs")

p <- ggplot(i_data, aes(z, y, fill = value)) +
  geom_tile(color = "white", lwd = 0, linetype = 1) +
  scale_fill_gradientn(colors = colfunc(135)) +
  # scale_y_continuous(expand=c(0, 0))+
  # scale_x_continuous(expand=c(0, 0),
  # breaks=round(seq(from = min(y), to = max(y), length.out = 9), digits = 1))+
  labs(x = "Uplift (m/s)", y = "Conspecific density", fill = "Selection probability") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 15),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = "none") +
  facet_wrap(~x, strip.position="left") #+
# gganimate::transition_states(x,
#                   transition_length = 4,
#                   state_length = 0.2,
#                   wrap = F) +
# gganimate::exit_fade() +
# labs(title = "Migration {closest_state}") 

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/facet_spring_3d.png",
    height = 18, width = 20, units = "in", res = 500)
plot_gg(p,multicore = TRUE,width=10,height=8,scale=250,windowsize=c(1400,866),
        zoom = 0.8, phi = 60, theta = -30)
render_snapshot(clear = TRUE)
dev.off()

library(raster)
xyz <- rasterFromXYZ(i_data[1:8, 1:3], res=c(NA,NA), crs="", digits=5)
xyz <- raster(xmn = min(i_data$x), xmx = max(i_data$x),
              ymn = min(i_data$y), ymx = max(i_data$y),
              ncols = 8, nrows = 15)
raster::setValues(xyz, i_data$z)
rasterVis::plot3D(xyz)
# install.packages("rayshader")
library(rayshader)
library(ggplot2)
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



ex <- expand.grid(x = 1:8,
                  y = seq(from = -1.85, to = 3.58, length.out = 15),
                  z = seq(from = 0, to = 2.41, length.out = 15)) %>% 
  mutate(value = sample(seq(from = 0, to = 1, length.out = 3000), 1800))

vol <- simplify2array(by(i_data, i_data$x, as.matrix))

library(misc3d)
con <- misc3d::computeContour3d(vol, max(vol), 1)
misc3d::drawScene(makeTriangles(con))
contour3d(vol, 1)



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
  scale_color_manual(values = c("#bd68ee", "#f48c3d")) +
  labs(y = "", fill = "Significant", color = "Significant", title = "Fall") +
  theme_classic() +
  annotate(geom="text", x=16, y=c(1:6), label=post_graph$significance[post_graph$Variable != "Fall migrations"], color="black", size = 6) +
  theme(axis.text = element_text(color = "black", size = 15),
        axis.line = element_line(linewidth = 1),
        text = element_text(size = 15),
        legend.text = element_text(size = 10))
pre_coefs <- ggplot(pre_graph %>% filter(Variable != "Spring migrations"), aes(Estimate, Variable)) +
  geom_vline(xintercept = 0, lty = 2, linewidth = 1) +
  geom_pointrange(aes(xmin = Lower, xmax = Upper), #fill = significance_05, color = significance_05), 
                  linewidth = 1, size = 1) +
  scale_color_manual(values = c("#bd68ee", "#f48c3d")) +
  labs(y = "", fill = "Significant", color = "Significant", title = "Spring") +
  theme_classic() +
  annotate(geom="text", x=5, y=c(1:6), label=pre_graph$significance[pre_graph$Variable != "Spring migrations"], color="black", size = 6) +
  theme(axis.text = element_text(color = "black", size = 15),
        axis.line = element_line(linewidth = 1),
        text = element_text(size = 15),
        legend.text = element_text(size = 10))

coefs_plot <-  plot_grid(post_coefs, pre_coefs, labels = c("A", "B"), nrow = 2,
                         align = 'v', axis = 'l')
coefs_plot

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/coefs_glmm.png",
    height = 150, width = 300, units = "mm", res = 500)
coefs_plot
dev.off()


# look at homogeneity
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
