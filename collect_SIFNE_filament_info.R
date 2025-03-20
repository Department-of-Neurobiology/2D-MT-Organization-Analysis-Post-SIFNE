# --------------------------------------------------------------
# Script Name: collect_SIFNE_filament_info.R
# Purpose: Collect data from SIFNE analysis output
# Author: Nataliya Trushina
# Date: 2022-06-10
# Last Modified: 2023-05-03
# --------------------------------------------------------------

# Requirements: 
# config_1.yaml file in the directory of this R script.

# Excel files exported from SIFNE. 
# One cell per folder with region folders inside:
# main_directory/
#   ├── cell_1/
#   │   ├── region_1/
#   │   │   ├── IntegratedInfo.xlsx
#   │   │   ├── (generated_plots_and_files)
#   │   ├── region_2/
#   │   │   ├── IntegratedInfo.xlsx
#   │   │   ├── (generated_plots_and_files)
#   │   └── ...
#   ├── cell_2/
#   │   ├── region_1/
#   │   │   ├── IntegratedInfo.xlsx
#   │   │   ├── (generated_plots_and_files)
#   │   ├── region_2/
#   │   │   ├── IntegratedInfo.xlsx
#   │   │   ├── (generated_plots_and_files)
#   │   └── ...
#   └── ...

# install.packages("svglite")
# devtools::install_github("marcusvolz/mathart")

# One can clear out previous analysis by deleting all JPG and PNG files in the current directory and its subdirectories.
# Run the following commands in Ubuntu terminal:
# find . -name "*.jpg" -type f -delete
# find . -name "*.png" -type f -delete

# Convert all TIFF files in the current directory and its subdirectories to SVG format using ImageMagick's 'convert' command.
# Run the following command in Ubuntu terminal:
# find . -type f -name '*.tif' -exec bash -c 'dir=$(dirname "{}"); base=$(basename "{}" .tif); output_svg="$dir/$base.svg"; convert "{}" "$output_svg"' \;
################################################

# Load libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mathart) 
library(sp)
library(sf)
library(yaml)

# Create a custom ggplot2 theme with a black background and white text/lines
theme_black <- function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Axis options
      axis.line = element_blank(),
      axis.text.x = element_text(size = base_size * 0.8, color = "white", lineheight = 0.9),
      axis.text.y = element_text(size = base_size * 0.8, color = "white", lineheight = 0.9),
      axis.ticks = element_line(color = "white", size = 0.2),
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),
      axis.ticks.length = unit(0.3, "lines"),
      
      # Legend options
      legend.background = element_rect(color = NA, fill = "black"),
      legend.key = element_rect(color = "white", fill = "black"),
      legend.key.size = unit(1.2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = element_text(size = base_size * 0.8, color = "white"),
      legend.title = element_text(size = base_size * 0.8, face = "bold", hjust = 0, color = "white"),
      legend.position = "right",
      legend.text.align = NULL,
      legend.title.align = NULL,
      legend.direction = "vertical",
      legend.box = NULL,
      
      # Panel options
      panel.background = element_rect(fill = "black", color = NA),
      panel.border = element_rect(fill = NA, color = "white"),
      panel.grid.major = element_line(color = "grey35"),
      panel.grid.minor = element_line(color = "grey20"),
      panel.spacing = unit(0.5, "lines"),
      
      # Facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),
      strip.text.x = element_text(size = base_size * 0.8, color = "white"),
      strip.text.y = element_text(size = base_size * 0.8, color = "white", angle = -90),
      
      # Plot options
      plot.background = element_rect(color = "black", fill = "black"),
      plot.title = element_text(size = base_size * 1.2, color = "white"),
      plot.margin = unit(rep(1, 4), "lines")
    )
}

# Find the indices of the largest element in a matrix
findmax <- function(m) {
  v = which.max(m) - 1
  c(v %% nrow(m) + 1, v %/% nrow(m) + 1)
}

# Calculate the maximum length chord across all pairs of vertex points in a polygon
max_chord <- function(polygon) {
  # Get polygon coordinates
  xy = st_coordinates(polygon)[, 1:2]
  
  # Compute the distance matrix and find the largest element
  df_dist = as.matrix(dist(xy))
  maxij = findmax(df_dist)
  
  # Define the largest chord using the elements with the largest distance
  chord = rbind(
    xy[maxij[1],],
    xy[maxij[2],]
  )
  chord
}

save_plot <- function(plot, file_prefix, cell, region, width, height) {
  # Save the plot as a PNG
  ggsave(plot, filename = paste("../", cell, "_", region, "_", file_prefix, ".png", sep=""), 
         type = "cairo", width = width, height = height)
  
  # Save the plot as an SVG
  ggsave(plot, filename = paste("../", cell, "_", region, "_", file_prefix, ".svg", sep=""), 
         width = width, height = height)
}

# Set the working directory to the directory of the active R script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read the configuration file
config <- yaml.load_file("config_1.yaml")

# Access the parameters from the config object
analysis_date <- config$analysis_date
genotype <- config$genotype
folder_path <- config$folder_path
concavity <- config$concavity
folder_regexp <- config$folder_regexp

# Set working paths
setwd(folder_path)

# IMPORTANT
# files_to_remove <- list.files(path = folder_path, pattern = "\\.(png|svg)$", recursive = TRUE, full.names = TRUE)
# file.remove(files_to_remove)

output_dir <- paste("../concavity_", concavity, "_collected", sep="")
dir.create(output_dir)
cell_folders <- Sys.glob(file.path(folder_regexp))

all_info <- {}
for (cell in cell_folders) {
  print(cell)
  # cell <- "20221108_B6_pos3_w529_h639_MAX_colored_15_18"
  setwd(paste(folder_path, cell, sep=""))
  region_folders <- Sys.glob(file.path("R*"))
  
  for (region in region_folders) {
    # Print the current region being processed
    print(region)
    
    # Set the working directory to the folder path, cell, and region
    setwd(paste(folder_path, cell, "/", region, sep=""))
    
    # Read in the data from the "IntegratedInfo.xlsx" file
    file <- "IntegratedInfo.xlsx"
    full_df <- read_excel(file, sheet = "Ultimate Filaments")
    
    # Filter out odd and even rows
    toDelete_1 <- seq(1, nrow(full_df), 2)
    toDelete_2 <- seq(2, nrow(full_df), 2)
    filament_df <- full_df[toDelete_1, ]
    
    # Add region and unit conversion columns to the data frame
    filament_df$Region <- region
    filament_df$`Total Length um` <- filament_df$`Total Length` * 0.026
    filament_df$`End-to-End Distance um` <- filament_df$`End-to-End Distance` * 0.026
    filament_df$Straightness <- filament_df$`End-to-End Distance` / filament_df$`Total Length`
    
    # Create a data frame with only coordinate information
    filament_df_coord <- select(filament_df, -c(
      "Straightness", "Orientation", "Total Length",
      "Region", "End-to-End Distance um", "Total Length um",
      "End-to-End Distance", "Centroid X", "Centroid Y",
      "1st X pos", "last X pos", "1st Y pos", "last Y pos",
    ))
    
    # Separate the X and Y coordinate columns
    toDelete_3 <- seq(1, ncol(filament_df_coord), 2)
    toDelete_4 <- seq(2, ncol(filament_df_coord), 2)
    filament_df_coord_x <- filament_df_coord[, toDelete_4]
    filament_df_coord_y <- filament_df_coord[, toDelete_3]
    
    # Add Filament ID to the X and Y coordinate data frames
    filament_df_coord_x$`Filament ID` <- filament_df_coord$`Filament ID`
    filament_df_coord_y$`Filament ID` <- filament_df_coord$`Filament ID`
    
    # Convert the X and Y coordinate data frames to long format
    filament_df_coord_x_long <- gather(filament_df_coord_x, Coordin.x, Coord.x, -"Filament ID")
    filament_df_coord_y_long <- gather(filament_df_coord_y, Coordin.y, Coord.y, -"Filament ID")
    
    # Combine the X and Y coordinate data frames
    filament_df_long_all <- cbind(filament_df_coord_x_long, select(filament_df_coord_y_long, -c("Filament ID")))
    
    # Remove rows with 0 values and incomplete cases
    filament_df_long_all[filament_df_long_all == 0] <- NA
    filament_df_long_all <- filament_df_long_all[complete.cases(filament_df_long_all), ]
    
    # Convert the Filament ID column to a factor
    filament_df_long_all$`Filament ID` <- as.factor(filament_df_long_all$`Filament ID`)
    
    # Set plot dimensions based on the range of the coordinates
    p_width <- (max(filament_df_long_all$Coord.x) - min(filament_df_long_all$Coord.x))/100
    p_height <- (max(filament_df_long_all$Coord.y) - min(filament_df_long_all$Coord.y))/100
    
    # Set a random seed and generate random colors for each filament
    set.seed(analysis_date)
    cols = rainbow(nlevels(filament_df_long_all$`Filament ID`), s=.6, v=.9)[sample(1:nlevels(filament_df_long_all$`Filament ID`), nlevels(filament_df_long_all$`Filament ID`))]
    
    # ---------------------------------------
    
    # Plot 1: Visualize filaments in random colors
    # p_filaments <- ggplot(filament_df_long_all, aes(x = Coord.y, y = Coord.x, group = `Filament ID`)) +
    #   geom_path(aes(color=`Filament ID`), size=1.2) +
    #   theme_black() +
    #   scale_color_manual(values=cols) +
    #   scale_y_reverse() +
    #   theme(
    #     axis.text.x = element_blank(),
    #     axis.text.y = element_blank(),
    #     axis.ticks = element_blank(),
    #     axis.title.x = element_blank(),
    #     axis.title.y = element_blank(),
    #     panel.grid.major = element_blank(),
    #     panel.grid.minor = element_blank(),
    #     legend.position = "none"
    #   )
    # 
    # save_plot(p_filaments, "filaments", cell, region, p_height, p_width)
    
    p_filaments <- ggplot(filament_df_long_all, aes(x = Coord.y, y = Coord.x, group = `Filament ID`)) +
      geom_path(aes(color = `Filament ID`), size = 1.2, lineend = "round", linejoin = "round") +
      scale_color_manual(values = cols) +
      theme_void() +
      scale_y_reverse() +
      theme(
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    save_plot(p_filaments, "filaments", cell, region, p_height, p_width)
    
    # ---------------------------------------
    
    # Area calculation
    dat <- as.matrix(filament_df_long_all[c("Coord.x", "Coord.y")])
    
    # Calculate the concave hull using the given coordinates
    ch <- concave_hull(dat, concavity = 3, length_threshold = 0)
    coords <- ch
    colnames(coords) <- c("Coord.x", "Coord.y")
    
    # Create a polygon object from the concave hull coordinates
    chull.poly <- Polygon(coords, hole = F)
    chull.area <- chull.poly@area
    
    # Convert the coordinates to a polygon object for further processing
    polygon <- st_polygon(list(as.matrix(rbind(coords, coords[1,]))))
    
    # Calculate the maximum chord of the polygon
    chord <- max_chord(polygon)
    
    # Plot 2: Visualize the concave hull and maximum chord on the image
    # p_concave_hull <- ggplot(filament_df_long_all, aes(x = Coord.y, y = Coord.x)) +
    #   geom_path(aes(x = Coord.y, y = Coord.x, group = `Filament ID`), color="white", size=1.2) +
    #   geom_path(data = as.data.frame(coords), aes(x = Coord.y, y = Coord.x), color="red", size=1.2) +
    #   geom_path(data = as.data.frame(chord), aes(x = Y, y = X), color="green", size=1.2) +
    #   theme_black() +
    #   scale_y_reverse() +
    #   theme(
    #     axis.text.x = element_blank(),
    #     axis.text.y = element_blank(),
    #     axis.ticks = element_blank(),
    #     axis.title.x = element_blank(),
    #     axis.title.y = element_blank(),
    #     panel.grid.major = element_blank(),
    #     panel.grid.minor = element_blank()
    #   )
    # 
    # save_plot(p_concave_hull, "area_calc", cell, region, p_height, p_width)
    
    p_concave_hull <- ggplot(filament_df_long_all, aes(x = Coord.y, y = Coord.x)) +
      geom_path(aes(group = `Filament ID`), color = "white", size = 1.2, lineend = "round", linejoin = "round") +
      geom_path(data = as.data.frame(coords), aes(x = Coord.y, y = Coord.x), color = "red", size = 1.2, lineend = "round", linejoin = "round") +
      geom_path(data = as.data.frame(chord), aes(x = Y, y = X), color = "green", size = 1.2, lineend = "round", linejoin = "round") +
      scale_color_manual(values = c("white", "red", "green")) +
      theme_void() +
      scale_y_reverse() +
      theme(
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    save_plot(p_concave_hull, "area_calc", cell, region, p_height, p_width)
    
    # ---------------------------------------
    
    # Junction count calculation
    junction_df <- full_df[toDelete_2, ]
    junction_df$`Junction count` <- NA
    
    # Iterate through each row of junction_df and calculate the junction count
    for (i in 1:nrow(junction_df)) {
      junction_df$`Junction count`[i] <- (sum(junction_df[i, 11:(ncol(junction_df) - 1)] != 0)) / 2
    }
    junction_df$`Cell Name` <- region
    
    # Extract and organize junction coordinates from junction_df
    junction_df_coord <- select(junction_df, -c("Cell Name", "Junction count", "Orientation", "Total Length", "End-to-End Distance", "Total Length", "End-to-End Distance", "Centroid X", "Centroid Y", "1st X pos", "last X pos", "1st Y pos", "last Y pos"))
    
    # Split the data frame into X and Y coordinate columns
    toDelete_5 <- seq(1, ncol(junction_df_coord), 2)
    toDelete_6 <- seq(2, ncol(junction_df_coord), 2)
    junct_df_coord_x <- junction_df_coord[, toDelete_5]
    junct_df_coord_y <- junction_df_coord[, toDelete_6]
    junct_df_coord_x$`Filament ID` <- junction_df_coord$`Filament ID`
    junct_df_coord_y$`Filament ID` <- junction_df_coord$`Filament ID`
    
    # Convert the data frames from wide to long format
    junct_df_coord_x_long <- gather(junct_df_coord_x, Coordin.x, Coord.x, -"Filament ID")
    junct_df_coord_y_long <- gather(junct_df_coord_y, Coordin.y, Coord.y, -"Filament ID")
    
    # Combine the X and Y coordinate data frames
    junct_df_long_all <- cbind(junct_df_coord_x_long, select(junct_df_coord_y_long, -c("Filament ID")))
    junct_df_long_all[junct_df_long_all == 0] <- NA
    junct_df_long_all <- junct_df_long_all[complete.cases(junct_df_long_all), ]
    junct_df_long_all$`Filament ID` <- as.factor(junct_df_long_all$`Filament ID`)
    
    # Remove duplicates and unnecessary columns
    junct_df_long_all <- select(junct_df_long_all, -c("Coordin.x", "Coordin.y"))
    junct_df_long_all <- junct_df_long_all[!duplicated(junct_df_long_all[c("Coord.x","Coord.y")]),]
    
    # Plot 3: Visualize filaments and junctions
    # p_junc <- ggplot(filament_df_long_all, aes(x = Coord.y, y = Coord.x, group = `Filament ID`)) +
    #   geom_path(color="white", size=1.2) +
    #   geom_point(data = junct_df_long_all, aes(x = Coord.x, y = Coord.y, group = `Filament ID`), size=1.5, color="red") +
    #   theme_black() +
    #   scale_y_reverse() +
    #   theme(
    #     axis.text.x = element_blank(),
    #     axis.text.y = element_blank(),
    #     axis.ticks = element_blank(),
    #     axis.title.x = element_blank(),
    #     axis.title.y = element_blank(),
    #     panel.grid.major = element_blank(),
    #     panel.grid.minor = element_blank()
    #   )
    # 
    # save_plot(p_junc, "junctions", cell, region, p_height, p_width)
    
    p_junc <- ggplot(filament_df_long_all, aes(x = Coord.y, y = Coord.x, group = `Filament ID`)) +
      geom_path(color = "white", size = 1.2, lineend = "round", linejoin = "round") +
      geom_point(data = junct_df_long_all, aes(x = Coord.x, y = Coord.y, group = `Filament ID`), size = 1.5, color = "red") +
      scale_color_manual(values = c("white", "red")) +
      theme_void() +
      scale_y_reverse() +
      theme(
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    save_plot(p_junc, "junctions", cell, region, p_height, p_width)
    
    # ---------------------------------------
    
    # Organizing data for export
    filament_df_w_jun <- merge(filament_df, junction_df[c("Filament ID", "Junction count")], by = "Filament ID")
    new_all_info <- as.data.frame(filament_df_w_jun[, c(
      "Region", "Filament ID",
      "Total Length um", "End-to-End Distance um", "Straightness",
      "Junction count"
    )])
    
    # Calculate additional metrics
    new_all_info$`Junction per total length` <- new_all_info$`Junction count`/new_all_info$`Total Length um`
    new_all_info$`Area um^2` <- chull.area*(0.026^2)
    new_all_info$`Number of MTs` <- nrow(new_all_info)
    new_all_info$`MTs per um^2` <- new_all_info$`Number of MTs`/new_all_info$`Area um^2`
    new_all_info$`Junctions per MT` <- new_all_info$`Junction count`/new_all_info$`Number of MTs`
    new_all_info$`MT density` <- new_all_info$`Number of MTs`/new_all_info$`Area um^2`
    new_all_info$`MT mass` <- sum(new_all_info$`Total Length um`)/new_all_info$`Area um^2`
    new_all_info$`MT mass Eike` <- new_all_info$`Number of MTs`*mean(new_all_info$`Total Length um`)/new_all_info$`Area um^2`
    new_all_info$`Region length` <- as.numeric(dist(chord)) * 0.026
    new_all_info$`MT density by region length` <- new_all_info$`Number of MTs`/as.numeric(dist(chord))
    new_all_info$`MT mass by region length` <- sum(new_all_info$`Total Length um`)/as.numeric(dist(chord))
    new_all_info$`MT mass by region length Eike` <- new_all_info$`Number of MTs`*mean(new_all_info$`Total Length um`)/as.numeric(dist(chord))
    new_all_info$`Cell name` <- cell
    new_all_info$Genotype <- genotype
    
    # Combine the new data with existing data
    all_info <- rbind(all_info, new_all_info)
  }
}

# Set the working directory to the analysis directory.
setwd(folder_path)

# Save the 'all_info' data frame as a CSV file in the specified output directory
write.table(all_info, paste(output_dir, "/all_info_for_plot_", genotype, ".csv", sep = ""), sep = ",", row.names = FALSE)
