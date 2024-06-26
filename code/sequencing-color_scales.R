# generic colors to use for divergent and sequential continuous color scales
library(viridis)
library(colorjam)
library(colorspace)

# generate ggplot colors
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}


# colors for heatmaps, as well as the FDR on the original knee plot for libraries
heatmap_colors <- viridis(1000)

# colors for scatter intensity plots (count depth, marker gene expression, etc.)
scatter_intensity_colors <- viridis(1000)

# colors for dot plots
dotplot_colors <- viridis(1000)

# colors for scatter plots of feature intensity
feature_intensity_colors <- c("grey90", "blue")

# need discrete colors that are divergent but distinct from viridis
# colors for cell and gene filters
filter_colors <- c("#F8766D", # fail (red)
                   "dodgerblue") # pass (blue)

# metadata colors
# colors for plCoA regions
region_colors <- c("#E87D72", # anterior plCoA (aplCoA)
                   "#4D86C6") # posterior plCoA (pplCoA)

# colors for individual library batches
# Batch Colors
## colors for individual library batches
batch_colors_rna_aplCoA <- c("#B2182B", "#67001F", "#F4A582")  #  aplCoA batches (hot/red)
batch_colors_rna_pplCoA <- c("#4393C3","#053061", "#92C5DE") # pplCoA batches (cool/blue)
batch_colors_rna <- c(batch_colors_rna_aplCoA, batch_colors_rna_pplCoA)

batch_colors_atac <- c("#B2182B", "#67001F", #  aplCoA batches (hot/red)
                       "#4393C3","#053061") # pplCoA batches (cool/blue)

# DEG intensity colors
deg_intensity_colors_aplCoA <- c("#CD5C5C", # indian red, 1.5
                                 "#B22222", # firebrick, 2.5
                                 "#800000") # maroon, 5

deg_intensity_colors_pplCoA <- c("#98BAE3", # light blue, 1.5
                                 "#5981B1", # mid blue, 2.5
                                 "#15273C") # dark blue, 5

# cell type colors
# neurons vs glia
neuron_color <- "dodgerblue"
glia_color <- "forestgreen"

tissue_colors <- c(neuron_color, # neurons
                   glia_color) # glia


spatial_colors_initial <- qualitative_hcl(11, palette = "dynamic")
spatial_colors_section <- qualitative_hcl(21, palette = "dynamic")

spatial_colors <- c("#B2182B", # aplCoA.1
                    "#15273C", # pplCoA.1
                    "#67001F", # aplCoA.2
                    "goldenrod1", # L1.1
                    "darkgoldenrod3", # L1.2
                    "#5981B1", # pplCoA.2
                    "#F4A582", # aplCoA.3
                    "darkorange", # OLG
                    "#008944", # NLOT
                    "violetred3", # aplCoA.4
                    "#FF0800") # aplCoA.5

glut_color <- "#ADD8E6" # Glut, salmon, separated into VG1 and VG2 blue and orange
glut_type_colors <- c("#B5DFCF", # Glut.1/Zfp536
                      "#7179BA", # Glut.2/Satb2
                      "#9B9EC7", # Glut.3/Gulp1
                      "#59788E", # Glut.4/Fign
                      "#82EEFD", # Glut.5/Reln
                      "#757C88", # Glut.6/Gpc5
                      "#251E5D", # Glut.7/Smoc1
                      "#0492C2", # Glut.8/Ntng1
                      "#0A1772", # Glut.9/Chrm2
                      "#151E3D", # Glut.10/Ebf1
                      "#3944BC", # Glut.11/Vwc2
                      "#48AAAD", # Glut.12/Ebf2
                      "#1520A6", # Glut.13/Etv1
                      "#7285A5")

glut_type_colors2 <- c("#B5DFCF", # Glut.1/Zfp536
                       "#7179BA", # Glut.2/Satb2
                       "#E0115F", # Glut.3/Gulp1
                       "#59788E", # Glut.4/Fign
                       "#82EEFD", # Glut.5/Reln
                       "#48AAAD", # Glut.6/Gpc5
                       "#FA8072", # Glut.7/Smoc1
                       "#0492C2", # Glut.8/Ntng1
                       "#0A1772", # Glut.9/Chrm2
                       "#FF0800", # Glut.10/Ebf1
                       "#3944BC", # Glut.11/Vwc2
                       "#800020", # Glut.12/Ebf2
                       "#89ffc4", # Glut.13/Etv1
                       "#008944")

glut_type_colors2_display <- c("#008944",
                               "#FA8072",
                               "#FF0800",
                               "#E0115F",
                               "#800020",
                               "#89ffc4",
                               "#3944BC",
                               "#7179BA",
                               "#0492C2",
                               "#82EEFD",
                               "#B5DFCF",
                               "#59788E",
                               "#48AAAD",
                               "#0A1772")

gaba_color <- "#FA8072" # GABA, shades of red, this specifically salmon
gaba_type_colors <- c("#FA8072", # GABA.1/Lamp5 Salmon
                      "#ED2939", # GABA.2/Nxph2 Imperial red
                      "#E0115F", # GABA.3/Vip pinkish
                      "#800000", # GABA.4/Nell2 Maroon
                      "#FF0800", # GABA.5/Npas1 Candy Apple
                      "#CD5C5C", # GABA.6/Sst Indian Red
                      "#B22222", # GABA.7/Tacr1 Firebrick
                      "#EA3C53", # GABA.8/Vcan Desire
                      "#B43757", # GABA.9/Npsr1 Hibiscus
                      "#933A16", # GABA.10/Fign Rust
                      "#5E1914", # GABA.11/Rai14 Sangria
                      "#CA3433",
                      "#800020") # GABA.12 Persian

gaba_type_colors2 <- qualitative_hcl(13, palette = "dynamic")
gaba_color_index <- c(2, 8, 4, 7, 10, 9, 5, 12, 3, 1, 11, 6, 13)
gaba_type_colors2_display <- gaba_type_colors2[order(gaba_color_index)]

neuron_class_colors <- c(glut_color, gaba_color)

nonneuron_class_colors <- c("springgreen", # astrocytes
                            "violetred1", # microglia
                            "violetred3", # macrophages
                            "darkorange", # OPC
                            "chocolate2", # NFOL
                            "chocolate4", # MOL
                            "darkgoldenrod1", # endo
                            "gold1", # mural
                            "goldenrod1", # ABC
                            "lightgoldenrod1") # VLMC

class_colors <- c(neuron_class_colors,
                  nonneuron_class_colors)

astro_type_colors <- c("#89ffc4", # Astro.1/Luzp2
                       "#008944") # Astro.2/Myoc

macro_type_colors <- c("#663047", # Macro.1/F13a1
                       "#4B0082") # Macro.2/Cd74

nfol_type_colors <- c("#E05A00", # NFOL.1/Frmd4a
                      "#A43100") # NFOL.2/Sgcd

mol_type_colors <- c("#844421", # MOL.1/Mast4
                     "#553C2B") # MOL.2/Prom1

vlmc_type_colors <- c("#FBC901", # ABC.1/Slc4a10
                     "#FFB302",
                     "#8B8000") # ABC.2/Ranbp3l

abc_type_colors <- c("#FAFA0F", # VLMC.1/Mgat4c
                      "#FDFD96", # VLMC.2/Bnc2
                      "#F0E68C") # VLMC.3/Hecw2

type_colors <- c(glut_type_colors,
                 gaba_type_colors,
                 astro_type_colors,
                 "violetred1", # microglia
                 macro_type_colors,
                 "darkorange", # OPC
                 nfol_type_colors,
                 mol_type_colors,
                 "darkgoldenrod1", # endo
                 "lightgoldenrod1", # mural
                 abc_type_colors,
                 vlmc_type_colors)

immune_colors <- c("violetred1", macro_type_colors)

region_tissue_colors <- NULL
for (i in 1:length(tissue_colors)){
  for (j in 1:length(region_colors)){
    region_tissue_colors <- c(region_tissue_colors, substr(blend_colors(c(tissue_colors[i], region_colors[j])), 1, 7))
  }
}

region_class_colors <- NULL
for (i in 1:length(class_colors)){
  for (j in 1:length(region_colors)){
    region_class_colors <- c(region_class_colors, substr(blend_colors(c(class_colors[i], region_colors[j])), 1, 7))
  }
}

region_type_colors <- NULL
for (i in 1:length(type_colors)){
  for (j in 1:length(region_colors)){
    region_type_colors <- c(region_type_colors, substr(blend_colors(c(type_colors[i], region_colors[j])), 1, 7))
  }
}

region_neuron_class_colors <- NULL
for (i in 1:length(neuron_class_colors)){
  for (j in 1:length(region_colors)){
    region_neuron_class_colors <- c(region_neuron_class_colors, substr(blend_colors(c(neuron_class_colors[i], region_colors[j])), 1, 7))
  }
}

region_nonneuron_class_colors <- NULL
for (i in 1:length(nonneuron_class_colors)){
  for (j in 1:length(region_colors)){
    region_nonneuron_class_colors <- c(region_nonneuron_class_colors, substr(blend_colors(c(nonneuron_class_colors[i], region_colors[j])), 1, 7))
  }
}

region_glut_type_colors <- NULL
for (i in 1:length(glut_type_colors)){
  for (j in 1:length(region_colors)){
    region_glut_type_colors <- c(region_glut_type_colors, substr(blend_colors(c(glut_type_colors[i], region_colors[j])), 1, 7))
  }
}

region_gaba_type_colors <- NULL
for (i in 1:length(gaba_type_colors)){
  for (j in 1:length(region_colors)){
    region_gaba_type_colors <- c(region_gaba_type_colors, substr(blend_colors(c(gaba_type_colors[i], region_colors[j])), 1, 7))
  }
}

region_astro_type_colors <- NULL
for (i in 1:length(astro_type_colors)){
  for (j in 1:length(region_colors)){
    region_astro_type_colors <- c(region_astro_type_colors, substr(blend_colors(c(astro_type_colors[i], region_colors[j])), 1, 7))
  }
}

region_macro_type_colors <- NULL
for (i in 1:length(macro_type_colors)){
  for (j in 1:length(region_colors)){
    region_macro_type_colors <- c(region_macro_type_colors, substr(blend_colors(c(macro_type_colors[i], region_colors[j])), 1, 7))
  }
}

region_nfol_type_colors <- NULL
for (i in 1:length(nfol_type_colors)){
  for (j in 1:length(region_colors)){
    region_nfol_type_colors <- c(region_nfol_type_colors, substr(blend_colors(c(nfol_type_colors[i], region_colors[j])), 1, 7))
  }
}

region_mol_type_colors <- NULL
for (i in 1:length(mol_type_colors)){
  for (j in 1:length(region_colors)){
    region_mol_type_colors <- c(region_mol_type_colors, substr(blend_colors(c(mol_type_colors[i], region_colors[j])), 1, 7))
  }
}

region_abc_type_colors <- NULL
for (i in 1:length(abc_type_colors)){
  for (j in 1:length(region_colors)){
    region_abc_type_colors <- c(region_abc_type_colors, substr(blend_colors(c(abc_type_colors[i], region_colors[j])), 1, 7))
  }
}

region_vlmc_type_colors <- NULL
for (i in 1:length(vlmc_type_colors)){
  for (j in 1:length(region_colors)){
    region_vlmc_type_colors <- c(region_vlmc_type_colors, substr(blend_colors(c(vlmc_type_colors[i], region_colors[j])), 1, 7))
  }
}

batch_color_list <- c("batch_colors_rna", "region_colors")
batch_color_list_atac <- c("batch_colors_atac", "region_colors")
tissue_color_list <- c("tissue_colors", "class_colors")
