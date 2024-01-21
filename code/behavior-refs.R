# miscellaneous code required for startup
require(MASS)
require(tidyverse)
require(colorRamps)
require(RColorBrewer)
require(spatstat)
require(viridis)
require(cowplot)
require(grid)
require(gridExtra)
require(colorspace)
require(scales)     #contains pretty_breaks() which allows for publication-quality x axes
library(zoo)
library(gmodels)
library(car)
library(multcomp)
library(rstatix)
library(emmeans)
library(ggpubr)
library(knitr)

# for 4quad
start_min <- 12
end_min <- 24

# constants for use in pipelines
topography_labels <- c("aplCoA eYFP", "aplCoA ChR2",
                       "mplCoA eYFP",  "mplCoA ChR2",
                       "pplCoA eYFP", "pplCoA ChR2")

projection_stim_labels <- c("MeA eYFP", "MeA ChR2",
                            "NAc eYFP", "NAc ChR2")
projection_labels <- c("MeA Gi", "MeA mCh", "NAc Gi", "NAc mCh")
projection_labels_tmt <- projection_labels_2PE <- projection_labels

genes_labels <- c("VGlut1 eYFP", "VGlut1 ChR2", "VGlut2 eYFP", "VGlut2 ChR2")
gene_labels <- c("VGlut1 Gi", "VGlut1 mCh", "VGlut2 Gi", "VGlut2 mCh")
gene_labels_2PE <- gene_labels_tmt <- gene_labels

# Linear plot constants
## 4-quad
### Performance Index
ap_linear_pi_axes <- c(-1.2, 1.2)
ap_linear_pi_label <-"Change in Performance Index [ON - OFF]"
### Port Distance
ap_linear_pd_axes <- c(-15, 15)
ap_linear_pd_label <- "Change in Mean Port Distance (cm) [ON - OFF]"

## EPM
### Open Time
ap_linear_opentime_axes <- c(-0.35, 0.35)
ap_linear_opentime_label <-"Change in Open Arm Time [ON - OFF]"
### Open Entry
ap_linear_entry_axes <- c(-15, 15)
ap_linear_entry_label <- "Change in Open Arm Entries [ON - OFF]"
### Distance
ap_linear_epm_distance_axes <- c(-8, 8)
ap_linear_epm_distance_label <- "Change in Distance Traveled (m) [ON - OFF]"

## OFT
### Center Time
ap_linear_centertime_axes <- c(-0.2, 0.2)
ap_linear_centertime_label <- "Change in Center Time [ON - OFF]"
### Center Time
ap_linear_cornertime_axes <- c(-0.25, 0.25)
ap_linear_cornertime_label <- "Change in Corner Time [ON - OFF]"
### Freezing
ap_linear_freezing_axes <- c(-0.2, 0.2)
ap_linear_freezing_label <- "Change in Freezing [ON - OFF]"
### Distance
ap_linear_oft_distance_axes <- c(-5, 5)
ap_linear_oft_distance_label <- "Change in Distance Traveled (m) [ON - OFF]"

# Shared labels

## 4quad
### Performance Index
rm_pi_label <- "Performance Index" # For repeated measures performance index
diff_pi_label <- "\u0394 Performance Index [ON - OFF]" # For looking at group differences in performance index
diff_pi_label_silence <- "\u0394 Performance Index" # For looking at group differences in performance index
### Port Distance
rm_pd_label <- "Mean Port Distance (cm)" # For repeated measures port distance
diff_pd_label <- "\u0394 Mean Port Distance (cm) [ON - OFF]" # For looking at group differences in port distance
diff_pd_label_silence <- "\u0394 Mean Port Distance (cm)" # For looking at group differences in port distance

## Elevated Plus Maze
epm_opentime_label <- "Proportion Time in Open Arms" # Open Arm Time
epm_entry_label <- "Open Arm Entries" # Open Arm Entries
epm_dist_label <- "Distance Traveled (m)" # Distance

## Open Field Test
oft_center_label <- "Proportion Time in Center" # Center Time
oft_corner_label <- "Proportion Time in Corners" # Corner Time
oft_freeze_label <- "Proportion Time Freezing" # Freezing
oft_dist_label <- "Distance Traveled (m)"
silence_diff_label <- "CNO - SAL"

# axis values for bar graphs
## 4quad
### Performance Index
ap_pi_axes <- c(-0.85, 0.85) # topography
genes_pi_axes <- c(-0.9, 0.9) # genes
proj_pi_axes <- c(-0.85, 0.85) # projections
### Port Distance
ap_pd_axes <- c(30, 46) # topography
ap_pd_axes_diff <- c(-7.5, 7.5) # topography difference
genes_pd_axes <- c(32, 44) # genes
genes_pd_axes_diff <- c(-10, 10) # genes difference
proj_pd_axes <- c(34, 46) # projections

## Elevated Plus Maze
### Open Arm Time
topo_ot_axes <- c(0, 0.75) # topographical pre-post
topo_ot_axes_min <- c(0, 1) # topographical min-wise
proj_ot_axes <- c(0, 0.75) # projection pre-post
proj_ot_axes_min <- c(0, 1) # projection min-wise
### Open Arm Entries
topo_entry_axes <- c(0, 25) # topographical pre-post
proj_entry_axes <- c(0, 25) # projection pre-post
proj_entry_axes_silence <- c(0, 50) # projection pre-post
topo_entry_axes_min <- c(0, 8) # topographical min-wise
proj_entry_axes_min <- c(0, 12) # projection min-wise
### Distance
topo_epm_dist_axes <- c(0, 25) # topographical pre-post
proj_epm_dist_axes <- c(0, 25) # projection pre-post
topo_epm_dist_axes_min <- c(0, 5) # topographical min-wise
proj_epm_dist_axes_min <- c(0, 6) # projection min-wise

## Open Field Test
### Center Time
topo_center_axes <- c(0, 0.5) # topographical pre-post
proj_center_axes <- c(0, 0.5) # projection pre-post
topo_center_axes_min <- c(0, 1) # topographical min-wise
proj_center_axes_min <- c(0, 1) # projection min-wise
### Corner Time
topo_corner_axes <- c(0, 0.6) # topographical pre-post
proj_corner_axes <- c(0, 0.6) # projection pre-post
topo_corner_axes_min <- c(0, 1) # topographical min-wise
proj_corner_axes_min <- c(0, 1) # projection min-wise
### Freezing
topo_freeze_axes <- c(0, 0.7) # topographical pre-post
proj_freeze_axes <- c(0, 0.7) # projection pre-post
topo_freeze_axes_min <- c(0, 1) # topographical min-wise
proj_freeze_axes_min <- c(0, 1) # projection min-wise
### Distance
topo_oft_dist_axes <- c(0, 25) # topographical pre-post
proj_oft_dist_axes <- c(0, 25) # projection pre-post
topo_oft_dist_axes_min <- c(0, 7) # topographical min-wise
proj_oft_dist_axes_min <- c(0, 7) # projection min-wise

# color scales for behavior and other plot constants

# metadata colors
# colors for plCoA regions
region_colors <- c("#E87D72", # anterior plCoA (aplCoA)
                   "#4D86C6") # posterior plCoA (pplCoA)

anterior_color <- "#E87D72"
posterior_color <- "#4D86C6"
mid_color <- "#C29DBC"

mea_color <- "#FA86C4" # pink
nac_color <- "#51087E" # indigo

VGlut1_color <- "#53868B" # slate blue
VGlut2_color <- "#F7931D" # orange

chr2_blue <- "#00B7FF" # ChR2 excitation wavelength (473nm)
eyfp_green <- "#52FF00" # eYFP excitation spectrum green (527nm)
chr2_eyfp <- c(chr2_blue, eyfp_green)

control_color <- "dimgrey"

tmt_colors <- "maroon"
pe_colors <- "lightseagreen"

# axis values for bar graphs

diff_of_axes <- c(-0.25, 0.25)
diff_cd_axes <- c(-5, 5)

pi_axes_tmt <- c(-1.2, 0.2)
pi_axes_tmt_genes <- c(-1, 0.3)
pi_axes_2PE <- c(-0.2, 0.8)
pi_axes_2PE_genes <- c(-0.2, 0.6)
diff_pi_axes_2PE <- c(-0.8, 0.8)
diff_pi_axes_2PE_genes <- c(-0.8, 0.8)

pd_axes_tmt <- c(-2, 15)
pd_axes_tmt_diff <- c(-7, 7)

pd_axes_2PE <- c(-8, 2)
pd_axes_2PE_genes <- c(-6, 2)
pd_axes_2PE_diff <- c(-7, 7)

epm_axes_entry_genes <- c(0, 50)

openfield_axes_dist_silence <- c(0, 50)
openfield_axes_dist_silence_epm <- c(0, 50)

openfield_center_axes_genes <- c(0, 0.4)

openfield_corner_axes_genes <- c(0, 0.9)

openfield_freeze_axes_genes <- c(0, 1)
