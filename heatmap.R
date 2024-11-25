library(community)
library(ggplot2)
library(ggrepel) # to add text labels on the mean weigth vs number of interactions plot
suppressPackageStartupMessages(library(ComplexHeatmap)) # to plot heatmaps
suppressPackageStartupMessages(library(circlize)) # for gragient colors
library(igraph) # to plot circus plots
library(gridExtra) 
library(plotly)
library(dplyr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_strings <- unlist(strsplit(args[1], " "))

load("./computed_results/interactions.RData")
data("visualization_functions")

idx <- which(interactions$anno_interactions$interaction_ID %in% input_strings)

filtered_list <- list()
for (name in names(interactions)) {
if (is.data.frame(interactions[[name]])) {
    if (nrow(interactions[[name]]) == 151744) {
    filtered_list[[name]] <- interactions[[name]][idx, ]
    } else {
    # Keep other dataframes as is
    filtered_list[[name]] <- interactions[[name]]
    }
} else if (is.list(interactions[[name]]) && length(interactions[[name]]) == 151744) {
    # Handle lists of length 151744
    filtered_list[[name]] <- interactions[[name]][idx]
} else {
    # Keep other elements as is
    filtered_list[[name]] <- interactions[[name]]
}
}

rm(interactions)
interactions <- filtered_list

# heatmap of interactions weight of top differential interactions
# top adjusted p value
idx_topsign <- interactions$anno_interactions$p.adj %in% unique(sort(interactions$anno_interactions$p.adj))
idx_topdown <- interactions$anno_interactions$log2FC_weights < -1 & idx_topsign
idx_topdown <- interactions$anno_interactions$interaction_ID %in% interactions$anno_interactions$interaction_ID[idx_topdown][1:10]

# all upregulated interactions
idx_topsign <- interactions$anno_interactions$p.adj %in% unique(sort(interactions$anno_interactions$p.adj))
idx_topup <- interactions$anno_interactions$log2FC_weights > 1 & idx_topsign
idx_topup <- interactions$anno_interactions$interaction_ID %in% interactions$anno_interactions$interaction_ID[idx_topup][1:9]

idx <- idx_topup | idx_topdown

# centred Heatmap
set.seed(3)
png("./plots/heatmap.png", width = 9.5, height = 8.5, units = "in", res = 300)
plot_heatmap(interactions,
            which_interactions = idx,
            values_to_plot = "weights",
            row_font_size = 18,
            column_font_size = 20,
            centered = TRUE,
            color_values = circlize::colorRamp2(c(-1, 0, 1), c("gray90", "white", "red3")),
            legend_title_font_size = 22,
            labels_font_size = 18)

dev.off()