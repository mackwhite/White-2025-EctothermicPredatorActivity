###project: EctothermPredatorActivityDrivers
###author(s): MW
###goal(s): 
###date(s): February 2025
###note(s): 
# Housekeeping ------------------------------------------------------------

### define custom functions ----
### function to scale marsh stage to the range of temperature ---
scale_stage <- function(x) {
      (x - min_stage) / (max_stage - min_stage) * (max_temp - min_temp) + min_temp
}

### inverse of above function to rescale seconday y axis with real values ---
inv_scale_stage <- function(x) {
      (x - min_temp) / (max_temp - min_temp) * (max_stage - min_stage) + min_stage
}


### load necessary libraries ----
# install.packages("librarian")
librarian::shelf(tidyverse, readr, lme4, ggpubr, performance, suncalc, ggpubr,
                 scales, ggstats, ggeffects, visreg, mgcv, MuMIn, glmmTMB, corrplot)

### load necessary data ----
all <- read_csv('local-data/snook-acc-model-data.csv') |>
      mutate(time_og = time,
             time = hour(hms(time_og)),
             y = mean_acceleration) |>
      rename(weight = weight_vbf,
             temp = temp_c) |> 
      mutate(id = as.factor(id), station = as.factor(station)) |>
      dplyr::select(y, time, weight, temp, stage, id, station)
glimpse(all)

### rms conversion coefficients ---
slope <- 0.01922
intercept <- 0

og <- read_rds('local-data/archive/accelerometer-model-data-012025.RDS') |> 
      filter(species == "common snook") |> 
      rename(acceleration_raw = acceleration) |> 
      mutate(acceleration = slope*acceleration_raw + intercept) |> 
      rename(stage = marsh_stage,
             y = acceleration) |> 
      mutate(id = as.factor(id),
             station = as.factor(station),
             month = as.factor(month))

# individual acc boxplot --------------------------------------------------

summ <- og |> group_by(id) |> summarize(median = median(y)) |> arrange(median)
print(summ)
order <- c('5362', '5039', '5363', '5361', '5042', '5043', '5040', '5367')

a <- all |> 
      mutate(x = id) |> 
      mutate(x = factor(x, levels = order)) |> 
      ggplot(aes(x = x, y = y)) + 
      geom_boxplot(fill = "white", outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      labs(y = expression(bold("Activity (m/s"^2*")"))) +
      theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
            axis.title = element_blank(),
            plot.title = element_text(size = 16, face = "bold", colour = "black"),
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = "none",
            legend.text = element_text(face = 'bold'),
            legend.title = element_text(face = 'bold')) +
      geom_text(aes(x = max(as.numeric(factor(all$id))), 
                    y = max(all$y, na.rm = TRUE), 
                    label = 'Data Binned Hourly'),
                vjust = 0.25, hjust = 1.0, size = 5.0, fontface = "bold", inherit.aes = FALSE)
a

b <- og |> 
      mutate(x = id) |> 
      mutate(x = factor(x, levels = order)) |> 
      ggplot(aes(x = x, y = y)) + 
      geom_boxplot(fill = "white", outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      labs(y = expression(bold("Activity (m/s"^2*")"))) +
      # scale_y_continuous(breaks = c(0,50,100,150,200,250), limits = c(0,260)) +
      theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
            axis.title = element_blank(),
            plot.title = element_text(size = 16, face = "bold", colour = "black"),
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = "none",
            legend.text = element_text(face = 'bold'),
            legend.title = element_text(face = 'bold')) +
      geom_text(aes(x = max(as.numeric(factor(og$id))), 
                    y = max(og$y, na.rm = TRUE), 
                    label = 'Raw Data'),
                vjust = 0.25, hjust = 1.0, size = 5.0, fontface = "bold", inherit.aes = FALSE)

b

### combine a and b plots
s1 <- ggarrange(a, b, ncol = 1, align = "v", labels = c("a", "b"))

### add figure axis titles
s1 |> 
      annotate_figure(
            left = text_grob(label = expression(bold("Activity (m/s"^2*")")),
                             just = 'centre',
                             color = "black",
                             face = "bold",
                             size = 18,
                             rot = 90)
      )

# temperature and hydrology time series -----------------------------------
temp_hr <- read_csv("local-data/botcreek_temp_15mins.csv")
temphr1 <- temp_hr |> 
      mutate(datetime_est = as.POSIXct(datetime, format = "%m/%d/%y %H:%M")) |> 
      filter(datetime_est >= as.Date("2024-01-14") & datetime_est <= as.Date("2024-12-14"))

temp <- temphr1 |> 
      mutate(date = as.Date(datetime_est)) |> 
      group_by(date) |> 
      summarise(
            temp = mean(temp_c, na.rm = TRUE),
      )

stage <- readxl::read_xlsx('../MAP/data/hydrology/mo215_current.xlsx') |> 
      mutate(date = as.Date(Date)) |> 
      filter(date >= as.Date("2024-01-14") & date <= as.Date("2024-12-14")) |> 
      rename(stage = `Stage (cm)`) |> 
      select(date, stage) |> 
      arrange(date)

env <- temp |> left_join(stage)
glimpse(env)

env_long <- env |> 
      pivot_longer(cols = c('temp', 'stage'),
                   names_to = 'metric',
                   values_to = 'value')

keep <- c("env_long", "scale_stage", "inv_scale_stage")
rm(list = setdiff(ls(), keep))

env_wide <- env_long |>
      filter(metric %in% c("temp", "stage")) |>
      select(date, metric, value) |>
      pivot_wider(names_from = metric, values_from = value)

### defining range of data for transformations ----
min_temp <- min(env_wide$temp, na.rm = TRUE)
max_temp <- max(env_wide$temp, na.rm = TRUE)
min_stage <- min(env_wide$stage, na.rm = TRUE)
max_stage <- max(env_wide$stage, na.rm = TRUE)

### scale data for secondary y axis using function from above ----
env_wide <- env_wide |> 
      mutate(stage_scaled = scale_stage(stage))

### grab monthly breaks for below plot ----
month_breaks <- env_wide |> 
      mutate(first_of_month = floor_date(date, "month")) |>   # Get first day of month
      pull(first_of_month) |> 
      unique() |> 
      sort()

### set y axis breaks - manually because of the two axes ----
y_primary_limits <- c(16, 33)
y_primary_breaks <- seq(16, 32, by = 4)

# Plot with second axis
env_wide |> 
      ggplot(aes(x = date)) +
      geom_line(aes(y = temp, color = "Temperature"), linetype = 'solid', linewidth = 1.25) +
      geom_line(aes(y = stage_scaled, color = "Marsh Stage"), linetype = 'solid', linewidth = 1.25) +
      scale_y_continuous(
            name = "Temperature (°C)",
            limits = y_primary_limits,
            breaks = y_primary_breaks,
            sec.axis = sec_axis(
                  trans = ~ inv_scale_stage(.),
                  name = "Marsh Stage (cm)"
            )
      ) +
      scale_color_manual(values = c("Temperature" = "#e41a1c", 
                                    "Marsh Stage" = "#377eb8")) +
      scale_x_date(
            breaks = month_breaks,  
            labels = format(month_breaks, "%b")
      ) +
      labs(x = "Month", color = NULL) +
      theme(
            axis.text.x = element_text(size = 14, face = "bold", colour = "black", angle = 45, vjust = 0.8), 
            axis.text.y = element_text(size = 14, face = "bold", colour = "black"), 
            axis.title.x = element_text(size = 14, face = "bold", colour = "black"),
            axis.title.y.left = element_text(size = 15, face = "bold", colour = "#e41a1c"),
            axis.title.y.right = element_text(size = 15, face = "bold", colour = "#377eb8"),
            plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position = 'none',
            legend.text = element_text(size = 10, color = "black", face = 'bold'),
            legend.title = element_text(size = 10, color = "black", face = 'bold')
      )

# ggsave('output/env-ts-supp.png',
#        dpi = 600, units= 'in', bg = 'white', height = 5, width = 6)

# summary tables for tagged individuals -----------------------------------

all <- read_csv('local-data/snook-acc-model-data.csv') |>
      mutate(time_og = time,
             time = hour(hms(time_og)),
             y = mean_acceleration) |>
      rename(weight = weight_vbf,
             temp = temp_c) |> 
      mutate(id = as.factor(id), station = as.factor(station)) |>
      dplyr::select(y, time, weight, temp, stage, id, station,
                    fl_vbf, date)
glimpse(all)

total_obs <- tibble(
      id = c(5039,5040,5042,5043,5361,5362,5363,5367),
      obs = c(5313,18422,3565,9355,13696,2756,8876,6177)
) |> 
      mutate(id = as.character(id))

all |> 
      mutate(id = as.character(id)) |> 
      group_by(id) |> 
      summarize(
            weight = mean(weight*.001),
            fork_length = mean(fl_vbf/10),
            median_acc = median(y),
            mean_acc = mean(y),
            sd_acc = sd(y),
            cv_acc = sd(y)/mean(y),
            days_obs = n_distinct(date)) |> 
      ungroup() |> 
      left_join(total_obs, by = c("id")) #|> 
      # capture.output(file = "output/tagger-info.csv")

all |> 
      mutate(id = as.character(id)) |> 
      summarize(
            weight = mean(weight*.001),
            fork_length = mean(fl_vbf/10),
            median_acc = median(y),
            mean_acc = mean(y),
            sd_acc = sd(y),
            cv_acc = sd(y)/mean(y),
            days_obs = n_distinct(date)) |> 
      ungroup() #|> 
      # capture.output(file = "output/alltagger-summary-info.csv")
