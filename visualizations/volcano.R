library(community)
library(ggplot2)
library(ggrepel) # to add text labels on the mean weigth vs number of interactions plot
suppressPackageStartupMessages(library(ComplexHeatmap)) # to plot heatmaps
suppressPackageStartupMessages(library(circlize)) # for gragient colors
library(igraph) # to plot circus plots
library(gridExtra) 
library(plotly)

path_in <- './computed_results/'
suppressWarnings(load(paste0(path_in, "/interactions.RData")))
data("visualization_functions")

# volcano plot
options(repr.plot.height = 5, repr.plot.width = 5)
    
plots <- plot_vulcano(interactions)
ggsave(filename = paste0("./plots/Volcano_plot_", ".png"),
         plot = plots,
         width = 10,
         height = 8,
         dpi = 300)