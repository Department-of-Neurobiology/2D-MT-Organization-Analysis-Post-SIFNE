# --------------------------------------------------------------
# Script Name: plot_filaments_junctions.R
# Purpose: Make plots and save tables for collected data from SIFNE analysis
# Author: Nataliya Trushina
# Date: 2022-06-10
# Last Modified: 2023-05-04
# --------------------------------------------------------------

# Requirements: 
# A "concavity_*_collected" folder with "all_info_for_plot_{condition}.csv" files.
# config_2.yaml file in the directory of this R script.

# This script performs statistical testing on the data depending on the number of levels in the 'Condition' variable.
# If there are two levels, a two-sample t-test is performed using the 't.test' method, with the 'control_condition' as the reference group.
# If there are more than two levels, a one-way ANOVA is performed using the 'anova' method.
# The p-values obtained from these tests are displayed on the plot with the 'p.signif' label. 
# This label can be changed to 'p.format' to show the calculated values.
#
# For further analysis, the generated CSV files can be used to load the data into statistical software like GraphPad, R, or Python.
# This allows for more advanced statistical tes ts or creation of additional visualizations.

# Load required libraries
library(tidyr)
library(ggplot2)
library(Rmisc)
library(stringr)
library(plotly)
library(dplyr)
library(ggpubr)
library(scales)
library(gridExtra)
library(yaml)

# Set the working directory to the directory of the active R script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read the YAML configuration file
config <- yaml.load_file("config_2.yaml")

# Set working path from the YAML file
setwd(config$folder_path)

# Get the parameters from the YAML file
control_condition <- config$control_condition
conditions <- config$conditions
xintercept_stat <- config$distribution_plot$xintercept
dist_width <- config$distribution_plot$width
dist_height <- config$distribution_plot$height

# Define a custom color gradient for plots
col_gradient <- colorRampPalette(config$colors)

# Find all files matching the specified pattern and create a list of their paths
all_filenames_plot_together <- Sys.glob(file.path("all_info_*"))

# Initialize an empty data frame to store merged data
data_merge <- data.frame()

# Define a function for reading and merging data files
read_and_merge_data <- function(file_paths) {
  data_list <- lapply(file_paths, function(file) {
    read.table(file, sep = ",", header = TRUE)
  })
  merged_data <- bind_rows(data_list)
  return(merged_data)
}

# Read and merge data files
data_merge <- read_and_merge_data(all_filenames_plot_together)

# Rename the columns of the data_merge data frame
colnames(data_merge) <- c(
  "Region", "Filament ID",
  "Total_Length_um", "End_to_End_Distance_um", "Straightness",
  "Junction count", "Junction_per_total_length", 
  "Area.um.2", "Number.of.MTs", "MTs.per.um.2", 
  "Junctions.per.MT", "MT.density", "MT.mass", "MT.mass.Eike",
  "Region.length", "MT.density.by.region.legth", 
  "MT.mass.by.region.length", "MT.mass.by.region.length.Eike",
  "Name", "Condition"
)

# Update the condition names
table(data_merge$Condition)

# Reorder the conditions based on the order specified in the YAML file
data_merge$Condition <- factor(data_merge$Condition, levels = conditions)

# Define the custom color gradient based on the number of conditions
cbp <- col_gradient(nlevels(data_merge$Condition))

# Get the control condition from the YAML file
control_condition <- config$control_condition

# Save the merged data to a CSV file
write.table(data_merge, "all_track_info.csv", sep = ";", dec = '.', row.names = FALSE, col.names = TRUE)

# ----------------------------------------

# Create distribution plots
create_density_plot <- function(data, variable, x_label_title, output_file) {
  mu <- aggregate(data[[variable]] ~ Condition, data = data, FUN = xintercept_stat)
  colnames(mu) <- c("Condition", "grp.stat")
  
  # Calculate the density of each condition
  densities <- lapply(split(data[[variable]], data$Condition), density)
  
  # Find the maximum x value and density value for placing the label
  max_x <- max(sapply(densities, function(x) max(x$x)))
  max_y <- max(sapply(densities, function(x) max(x$y)))
  
  dens_plot <- ggplot(data, aes(x = !!sym(variable), color = Condition, fill = Condition)) + 
    geom_density(alpha = 0.4) +
    theme_classic() +
    geom_vline(data = mu, aes(xintercept = grp.stat, color = Condition), linetype = "dashed") +
    geom_text(data = mu, aes(x = max_x, y = max_y, label = paste(xintercept_stat, ": ", round(grp.stat, 3), sep = ""), color = Condition), hjust = 1, vjust = 1, position = position_nudge(y = rep(c(0, -0.1*max_y, 0.1*max_y), length.out = nrow(mu)))) +
    theme(
      plot.title = element_text(color = "black", size = 12),
      axis.ticks = element_line(colour = "black"),
      axis.text.x = element_text(color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 10)
    ) +
    scale_color_manual(values = cbp) +
    scale_fill_manual(values = cbp) +
    xlab(x_label_title) + 
    ylab("") +
    scale_x_continuous(breaks = pretty(data[[variable]], n = 10)) +
    scale_y_continuous(breaks = pretty_breaks(n = 5))
  ggsave(output_file, dens_plot, width = dist_width, height = dist_height)
  return(dens_plot)
}

# Create distribution histogram plots
create_hist_plot <- function(data, variable, x_label_title, output_file) {
  mu <- aggregate(data[[variable]] ~ Condition, data = data, FUN = xintercept_stat)
  colnames(mu) <- c("Condition", "grp.stat")

  hist_plot <- ggplot(data, aes(x = !!sym(variable), color = Condition, fill = Condition)) + 
    geom_histogram(binwidth=0.5, alpha = 0.1) + #aes(y=..density..), 
    theme_classic() +
    geom_vline(data = mu, aes(xintercept = grp.stat, color = Condition), linetype = "dashed") +
    theme(
      plot.title = element_text(color = "black", size = 12),
      axis.ticks = element_line(colour = "black"),
      axis.text.x = element_text(color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 10)
    ) +
    scale_color_manual(values = cbp) +
    scale_fill_manual(values = cbp) +
    xlab(x_label_title) + 
    ylab("") +
    scale_x_continuous(breaks = pretty_breaks(n = 10))
  ggsave(output_file, hist_plot, width = dist_width, height = dist_height)
  return(hist_plot)
}

# Call the function to create the plots with x_label_title added
dens_filaments <- create_density_plot(data_merge, "Total_Length_um", "Length (\u03BCm)", "distr_total_length.svg")
hist_filaments <- create_hist_plot(data_merge, "Total_Length_um", "Length (\u03BCm)", "distr_hist_total_length_true_num.svg")
dens_end_to_end <- create_density_plot(data_merge, "End_to_End_Distance_um", "End-to-End Distance (\u03BCm)", "distr_end_to_end.svg")
dens_straightness <- create_density_plot(data_merge, "Straightness", "MT straightness", "distr_straightness.svg")
dens_junctions <- create_density_plot(data_merge, "Junction_per_total_length", "Junctions per MT length", "distr_junct_per_length.svg")

all_distr_plots <- ggarrange(dens_filaments, dens_end_to_end, dens_junctions, dens_straightness,
                      labels = c("A", "B", "C", "D"),
                      ncol = 2, nrow = 2,
                      common.legend = TRUE, legend = "bottom")
ggsave("distr_all_plots.png", all_distr_plots, width = 10, height = 6)

# ----------------------------------------

# Round up the maximum y value to the next custom threshold value
round_y_max <- function(y_max) {
  # Define custom rounding thresholds. These thresholds are chosen to match the desired rounding behavior
  # for different ranges of numbers. For example, for numbers between 0.1 and 0.2, we want to round up to 0.2.
  rounding_thresholds <- c(0.1, 0.2, 0.4, 0.6, 1, 2, 4, 6, 8, 10, 20, 40, 60, 80, 100, 200, 400, 600, 800, 1000)
  # Find the interval in which y_max falls in the rounding_thresholds vector.
  # The findInterval function returns the index of the threshold that is immediately less than y_max, add one for the next one.
  round_up <- findInterval(y_max, rounding_thresholds) + 1
  return(rounding_thresholds[round_up])
}

# Custom function to create the plot
make_plot <- function(df, y_variable, y_label, control_condition, y_limits = NULL) {
  plot <- ggplot(df, aes(x = Condition, y = !!sym(y_variable))) +
    geom_boxplot(width = 0.5, color = "#555555") +
    geom_jitter(size = 2, width = 0.3, aes(color = Condition)) +
    scale_color_manual(values = cbp) +
    theme_classic() +
    theme(legend.position="none",
          axis.title.x = element_blank()) +
    expand_limits(y = 0) +
    ylab(y_label)

  if (is.null(y_limits)) {
    y_max <- max(df[[y_variable]], na.rm = TRUE)
    y_max_rounded <- round_y_max(y_max)
    plot <- plot + scale_y_continuous(limits = c(0, y_max_rounded), breaks = seq(0, y_max_rounded, y_max_rounded / 4))
  } else {
    y_limits <- c(y_limits[1], y_limits[2])
    plot <- plot + scale_y_continuous(limits = y_limits, breaks = seq(y_limits[1], y_limits[2], (y_limits[2] - y_limits[1]) / 4))
  }

  if (nlevels(data_merge$Condition) == 2) {
    plot <- plot + stat_compare_means(label = "p.signif", method = "t.test", ref.group = control_condition)
  } else {
    plot <- plot + stat_compare_means(label = "p.signif", method = "anova")
  }

  return(plot)
}

# List of variables, labels and y-axis limits
# Limits can be added for better representation of comparable parameters like total length and end-to-end length
# list(var = "Total_Length_um", label = "Length (um)", limits = c(0, 4))
variables <- list(
  list(var = "MT.mass", label = "MT mass"),
  # list(var = "MT.mass.Eike", label = "MT Mass Eike"),
  list(var = "MT.density", label = "MT density"),
  list(var = "MTs.per.um.2", label = "MTs per \u03BCm^2"),
  list(var = "Junctions.per.MT", label = "Junctions per MT"),
  list(var = "Total_Length_um", label = "MT length (\u03BCm)"),
  list(var = "End_to_End_Distance_um", label = "End-to-end distance (\u03BCm)"),
  list(var = "Straightness", label = "MT straightness", limits = c(0.8, 1)),
  list(var = "Junction_per_total_length", label = "Junctions per MT length"),
  list(var = "Area.um.2", label = "Area (\u03BCm^2)"),
  list(var = "Number.of.MTs", label = "Number of MTs"),
  list(var = "MTs.per.um.2", label = "MTs per \u03BCm^2"),
  list(var = "Region.length", label = "Region length"),
  list(var = "MT.mass.by.region.length", label = "MT mass by region length"),
  # list(var = "MT.mass.by.region.length.Eike", label = "MT Mass by Region Length Eike"),
  list(var = "MT.density.by.region.legth", label = "MT density by region length")
)

# ----------------------------------------

# Aggregate by cell
data_merge_to_average <- data_merge[c("Name", "Total_Length_um","End_to_End_Distance_um", 
                                      "Straightness", "Junction_per_total_length", 
                                      "Area.um.2", "Number.of.MTs", "MTs.per.um.2", 
                                      "Junctions.per.MT", "MT.density", "MT.mass", "MT.mass.Eike",
                                      "Region.length", "MT.density.by.region.legth", 
                                      "MT.mass.by.region.length", "MT.mass.by.region.length.Eike")]
df_for_plot_mean <- data_merge_to_average %>% group_by(`Name`) %>% dplyr::summarise(across(everything(), list(mean))) 
colnames(df_for_plot_mean) <- c("Name", "Total_Length_um","End_to_End_Distance_um", 
                                "Straightness", "Junction_per_total_length", 
                                "Area.um.2", "Number.of.MTs", "MTs.per.um.2", 
                                "Junctions.per.MT", "MT.density", "MT.mass", "MT.mass.Eike",
                                "Region.length", "MT.density.by.region.legth", 
                                "MT.mass.by.region.length", "MT.mass.by.region.length.Eike")
df_for_plot_mean_by_cell <- merge(df_for_plot_mean, unique(data_merge[c("Name", "Condition")]), by ="Name", all.x = TRUE)

# Sort the table by conditions
df_for_plot_mean_by_cell <- df_for_plot_mean_by_cell %>%
  arrange(Condition)

# ---------------------------------------------------------------------------- #

var_list <- setdiff(names(df_for_plot_mean_by_cell), "Condition")

data_list <- purrr::map(var_list, function(col_name) {
  temp <- df_for_plot_mean_by_cell %>%
    dplyr::select(Condition, all_of(col_name))
  
  # Flatten list columns by collapsing into a comma-separated string
  if (any(purrr::map_lgl(temp[[col_name]], is.list))) {
    temp[[col_name]] <- purrr::map_chr(temp[[col_name]], ~ paste(.x, collapse = ","))
  } else {
    if (all(purrr::map_lgl(temp[[col_name]], ~ stringr::str_detect(.x, "^[-+]?[0-9]*\\.?[0-9]+$")))) {
      temp[[col_name]] <- as.numeric(temp[[col_name]])
    } else {
      temp[[col_name]] <- as.character(temp[[col_name]])
    }
  }
  
  temp <- temp %>%
    dplyr::group_by(Condition) %>%
    dplyr::mutate(row_id = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  temp %>%
    tidyr::pivot_wider(names_from = Condition, values_from = all_of(col_name)) %>%
    dplyr::select(-row_id)
})

names(data_list) <- var_list

writexl::write_xlsx(data_list, path = "cell_means_by_cell.xlsx")

# ---------------------------------------------------------------------------- #

# writexl::write_xlsx(df_for_plot_mean_by_cell, "means_by_cell.xlsx")

# Create plots for each variable
cell_plots <- lapply(variables, function(x) {
  make_plot(df_for_plot_mean_by_cell, x$var, x$label, control_condition, x$limits)
})

# ---------------------------------------------------------------------------- #

# Aggregate by region
data_merge$`Cell region` <- paste(data_merge$Name, data_merge$Region)
data_merge_to_average <- data_merge[c("Cell region", "Total_Length_um",
                                      "End_to_End_Distance_um", "Straightness", 
                                      "Junction_per_total_length", 
                                      "Area.um.2", "Number.of.MTs", "MTs.per.um.2", 
                                      "Junctions.per.MT", "MT.density", 
                                      "MT.mass", "MT.mass.Eike",
                                      "Region.length", "MT.density.by.region.legth", 
                                      "MT.mass.by.region.length", 
                                      "MT.mass.by.region.length.Eike")]
df_for_plot_mean <- data_merge_to_average %>% group_by(`Cell region`) %>% dplyr::summarise(across(everything(), list(mean)))
colnames(df_for_plot_mean) <- c("Cell region", "Total_Length_um",
                                "End_to_End_Distance_um", "Straightness", 
                                "Junction_per_total_length",
                                "Area.um.2", "Number.of.MTs", 
                                "MTs.per.um.2", "Junctions.per.MT", 
                                "MT.density", "MT.mass", "MT.mass.Eike",
                                "Region.length", "MT.density.by.region.legth", 
                                "MT.mass.by.region.length", 
                                "MT.mass.by.region.length.Eike")
df_for_plot_mean_by_region <- merge(df_for_plot_mean, unique(data_merge[c("Cell region", "Condition", "Name", "Region")]), by ="Cell region", all.x = TRUE)

# Sort the table by conditions
df_for_plot_mean_by_region <- df_for_plot_mean_by_region %>%
  arrange(Condition)

# ---------------------------------------------------------------------------- #

var_list <- setdiff(names(df_for_plot_mean_by_region), "Condition")

data_list <- purrr::map(var_list, function(col_name) {
  temp <- df_for_plot_mean_by_region %>%
    dplyr::select(Condition, all_of(col_name))
  
  # Flatten list columns by collapsing into a comma-separated string
  if (any(purrr::map_lgl(temp[[col_name]], is.list))) {
    temp[[col_name]] <- purrr::map_chr(temp[[col_name]], ~ paste(.x, collapse = ","))
  } else {
    if (all(purrr::map_lgl(temp[[col_name]], ~ stringr::str_detect(.x, "^[-+]?[0-9]*\\.?[0-9]+$")))) {
      temp[[col_name]] <- as.numeric(temp[[col_name]])
    } else {
      temp[[col_name]] <- as.character(temp[[col_name]])
    }
  }
  
  temp <- temp %>%
    dplyr::group_by(Condition) %>%
    dplyr::mutate(row_id = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  temp %>%
    tidyr::pivot_wider(names_from = Condition, values_from = all_of(col_name)) %>%
    dplyr::select(-row_id)
})

names(data_list) <- var_list

writexl::write_xlsx(data_list, path = "cell_means_by_region.xlsx")

# ---------------------------------------------------------------------------- #

# writexl::write_xlsx(df_for_plot_mean_by_region, "mean_by_region.xlsx")

# Create plots for each variable
region_plots <- lapply(variables, function(x) {
  make_plot(df_for_plot_mean_by_region, x$var, x$label, control_condition, x$limits)
})

# ---------------------------------------------------------------------------- #

# Display and save plots in a grid (for both cell and region)
ggsave("means_by_cell.png", grid.arrange(grobs = cell_plots, ncol = 4), width = 10, height = 12)
ggsave("means_by_region.png", grid.arrange(grobs = region_plots, ncol = 4), width = 10, height = 12)

