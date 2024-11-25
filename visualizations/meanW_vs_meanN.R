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

interaction_cell_types <- unique(paste(interactions$anno_interactions$sending_cell_type,
                                 interactions$anno_interactions$receiving_cell_type,
                                 sep = " to "))

interaction_type <- sapply(interaction_cell_types,
                           function(i){
                               ifelse(grepl("Ery",i),
                                      "engages Ery",
                                      ifelse(grepl("HSPC",i),
                                                   "engages HSPC",
                                                   "among immune cells"))
                           })

colors_interaction_type <- c("engages HSPC" = "magenta",
                             "among immune cells" = "deepskyblue",
                             "engages Ery" = "darkgoldenrod3")

# plot mumber of interactions vs mean interaction weights per cell type to cell type interaction
options(repr.plot.width = 10, repr.plot.height = 10)

plots <- plot_nrInt_vs_meanW_perCellType(interactions,
                                interaction_type = interaction_type,
                                colors = colors_interaction_type,
                                ylim = c(-3.65, -0.7),
                                label_font_size = 10)

# save plots as images
for (i in 1:length(plots)) {
  ggsave(filename = paste0("./plots/Log10_vs_NrInteractions_plot_", i, ".png"),
         plot = plots[[i]],
         width = 10,
         height = 8,
         dpi = 300)}