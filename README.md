# 2D-MT-Organization-Analysis-Post-SIFNE
The scripts are used to perform data analysis and visualization of microtubule (MT) data. They extract and visualize filament information from SIFNE analysis output, creating custom plots for various MT parameters, such as MT mass, density, straightness, and junction counet, at both the cell and cell region levels.

Modify the configuration files (YAML) for the use of parameters and settings required for data analysis and visualization. 

"collect_SIFNE_filament_info.R" processes and visualizes output data from SIFNE filament analysis. It reads the Excel files generated by SIFNE, organizes and filters the data, and creates multiple plots to represent filament information.

"plot_filaments_junctions.R" calculates summary statistics for each experimental condition. It generates box plots and density plots for various MT parameters, such as MT mass, density, straightness, and junction count, at both the cell and region levels. The script also saves the aggregated data and plots to separate CSV and image files, respectively.
