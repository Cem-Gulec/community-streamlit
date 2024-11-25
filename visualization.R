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

update_network_colors <- function(immune_color, ery_color, gran_color) {
    # Load required libraries
    library(igraph)
    
    # Load the saved data
    path_in <- './computed_results/'
    load(paste0(path_in, "/interactions.RData"))
    
    # Re-define the interaction types (these were in the original visualization.R)
    cell_types <- unique(interactions$anno_cells$cell_type)
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
    
    # Update color mapping
    colors_interaction_type <- c(
        "engages HSPC" = "magenta",
        "among immune cells" = immune_color,
        "engages Ery" = ery_color,
        "engages Gran" = gran_color
    )
    
    # Get colors for edges
    colors <- colors_interaction_type[interaction_type]
    names(colors) <- names(interaction_type)
    
    # Source the original plotting function
    source("visualization.R")
    
    # Recreate network plots
    plorCelltypeNetwork(interactions,
                       edge.color = colors,
                       title_cex = 2.5,
                       vertex.label.cex = 3,
                       edge.arrow.size = 0.2,
                       verbose = FALSE)
}



plorCelltypeNetwork <- function(my_interactions
                                  ,amplify_edgeWidth = 50 # magnification for the edge width
                                  ,amplify_colorResolution = 20 # if color_palette is used
                                  ,nr_colors = 10 # if color_palette is used
                                  ,color_palette = c("red", "red4", "black")
                                  ,edge.color = NULL # custom edge color (named vector of colors, names are edge IDs)
                                  ,vertex.label.cex = 2
                                  ,vertex.shape="none"
                                  ,vertex.size = 10
                                  ,edge.arrow.size = 1
                                  ,title_cex = 1
                                  ,verbose = FALSE
                                  ,...){
        
        # extract health status
        health_status <- unique(my_interactions$anno_samples$health_status)
        
        # prepare the data for plotting
        data <- lapply(health_status
                       ,function(hs){
                               
                               idx_hs <- my_interactions$anno_samples$health_status == hs
                               
                               # create a matrix with mean weights
                               mean_w_mat <- matrix(NA
                                                    ,nrow = length(cell_types)
                                                    ,ncol = length(cell_types))
                               rownames(mean_w_mat) <- cell_types
                               colnames(mean_w_mat) <- cell_types
                               
                               # populate the matrix with mean weights
                               for(i in cell_types){
                                       for(j in cell_types){
                                               idx_send <- my_interactions$anno_interactions$sending_cell_type == i
                                               idx_rec <- my_interactions$anno_interactions$receiving_cell_type == j
                                               
                                               # calculate the mean for each significant interaction between the cell types of interest in the cohort
                                               mean_weight <- rowMeans(my_interactions$weights[idx_sign & idx_send & idx_rec
                                                                                               ,idx_hs]
                                               )
                                               # calculate the mean of non-zero values
                                               my_mean <- mean(mean_weight[mean_weight != 0])
                                               
                                               # populate the matrix
                                               ifelse(is.na(my_mean)
                                                      ,mean_w_mat[i,j] <- 0.000000001 # if an edge is completely missing, the plotting function throws an error
                                                      ,mean_w_mat[i,j] <- my_mean)
                                               
                                       }
                               }
                               
                               # create a graph object
                               g <- igraph::graph(unlist(lapply(rownames(mean_w_mat)
                                                                , function(i){
                                                                        lapply(colnames(mean_w_mat)
                                                                               ,function(j){
                                                                                       c(i,j)
                                                                               })
                                                                        
                                                                })
                               ))
                               
                               # collapse the graph ogject to unique edges
                               g_simp <- igraph::simplify(g, remove.loops = F)
                               
                               
                               
                               
                               
                               
                               summary <- data.frame(edge_ID = {# extract edges from my_graph object
                                       edges <- as_edgelist(g, names = TRUE)
                                       edges <- as.data.frame(sapply(1:nrow(edges),function(i)paste(edges[i,])))
                                       
                                       ID <- sapply(1:ncol(edges)
                                                    ,function(j) paste(edges[1,j],edges[2,j],sep = " to "))
                                       ID
                               }
                               ,edge.width = {
                                       # define edge width as the desired magnification of the edge weight
                                       edge.width <- unlist(lapply(1:nrow(mean_w_mat),function(i) mean_w_mat[i,]))
                                       edge.width <- edge.width*amplify_edgeWidth
                                       edge.width
                               }
                               ,color_resolution = {
                                       # define resolution for the colors
                                       color_resolution <- round(edge.width*amplify_colorResolution
                                                                 ,digits = 0)
                               }
                               ) 
                               
                               # return the list
                               list(g_simp = g_simp
                                    ,mean_w_mat = mean_w_mat
                                    ,summary = summary)
                               
                       })
        names(data) <- health_status
        
        if(verbose){
                print(str(data))
                lapply(health_status
                       ,function(hs){
                               print(hs)
                               print(data[[hs]]$summary)
                       })
        }
        
        # plot
        for(hs in health_status){
                
                # define maximum number of colors if color_palette is used
                nr_colors <- max(unlist(sapply(health_status,function(h){data[[h]]$summary$color_resolution})))
                
                # add color scheme for the edges 
                if(is.null(edge.color)){
                        edge.col <- colorRampPalette(color_palette)
                        edge.color <- edge.col(nr_colors)[df$colorResolution+1]
                } else {
                        edge.color <- edge.color[data[[hs]]$summary$edge_ID] # make sure it is sorted correctly
                }
                
                # plot network
                png(paste0("./plots/network_plot_", hs, ".png"), width = 7, height = 7, units = "in", res = 300)
                plot(data[[hs]]$g_simp
                     ,edge.width = data[[hs]]$summary$edge.width 
                     ,edge.arrow.size=edge.arrow.size
                     ,edge.alpha=0.5
                     ,edge.curved=0.05
                     ,edge.color = edge.color
                     ,edge.attr.comb=c(weight="sum", type="ignore")
                     ,vertex.size=vertex.size
                     ,vertex.label.cex = vertex.label.cex
                     ,vertex.shape=vertex.shape
                     ,main = ""
                     ,layout = layout_in_circle(data[[hs]]$g_simp) 
                     ,...
                )
                title(hs, cex.main = title_cex)
                dev.off()
        }
}

forestplot_for_category <- function(IDs, filename){
    # revert
    IDs <- IDs[length(IDs):1]
    
    rownames(interactions$anno_interactions) <- interactions$anno_interactions$interaction_ID
    my_anno_interactions <- interactions$anno_interactions[IDs,]
    my_anno_interactions$interaction_ID <- factor(my_anno_interactions$interaction_ID
                                                 ,ordered = TRUE
                                                        )
    
    print(paste(nrow(my_anno_interactions), "interactions in this category"))
    
    png(filename, width = 17, height = 4, units = "in", res = 300)
    plot_all_forests(my_idx = rep(TRUE,nrow(my_anno_interactions))
                     ,my_anno_interactions = my_anno_interactions
                     ,keep_order = TRUE
                     ,threshold = 1
                     ,legend_title_size = 0
                     ,legend_text_size = 20
                     ,component_lim = component_lim
                    )
}

# threshold for log2FC of the weights
threshold_log2FC <- interactions$thresholds$threshold_log2FC

idx_up <- interactions$anno_interactions$log2FC_weights > threshold_log2FC
idx_unchanged <- abs(interactions$anno_interactions$log2FC_weights) <= threshold_log2FC
idx_down <- interactions$anno_interactions$log2FC_weights < -threshold_log2FC
idx_sign <- !is.na(interactions$anno_interactions$sign) & interactions$anno_interactions$sign
idx_good <- interactions$anno_interactions$passed_QC_filter

# extract cell types
cell_types <- unique(interactions$anno_cells$cell_type)

# define broad type of interactions by cell types
immune_cell_types <- cell_types[!(cell_types %in% c("Ery", "HSPC"))]

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

interaction_type_by_ID <- sapply(interactions$anno_interactions$interaction_ID,
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

# claclulate mean weight of good interactions for each cell type to cell type communication (i.e. "T to B", "DC to T")
mean_weights_goodInteractions <- mean_weights(interactions)

# claclulate number of good interactions for each cell type to cell type communication (i.e. "T to B", "DC to T")
number_goodInteractions <- number_interactions(interactions)

immune_cell_types <- cell_types[!(cell_types %in% c("Ery", "HSPC"))]

# plot mumber of interactions vs mean interaction weights per cell type to cell type interaction
options(repr.plot.width = 10, repr.plot.height = 10)

xlim <- c(-20, 700)

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

print(paste0("downregulated:", sum(idx_down & idx_sign)))
print(paste0("upregulated:", sum(idx_up & idx_sign)))

# volcano plot
options(repr.plot.height = 5, repr.plot.width = 5)
    
plots <- plot_vulcano(interactions)
ggsave(filename = paste0("./plots/Volcano_plot_", ".png"),
         plot = plots,
         width = 10,
         height = 8,
         dpi = 300)

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

# cell type network
options(repr.plot.width = 7, repr.plot.height = 7)

colors <- colors_interaction_type[interaction_type]
names(colors) <- names(interaction_type)

plorCelltypeNetwork(interactions,
                    edge.color = colors,
                    title_cex = 2.5,
                    vertex.label.cex = 3,
                    edge.arrow.size = 0.2,
                    verbose = FALSE)


## Forest plots
ordered_IDs <- order_IDs_by_interaction_categories(interactions)

# split by category
component_lim <- find_component_limits(interactions, ordered_IDs)

# forest maps
IDs <- as.character(c(ordered_IDs$ID_order_no_change))
forestplot_for_category(IDs, "./plots/forest/forestplot_no_change.png")

IDs <- as.character(c(ordered_IDs$ID_order_rho_s_only_down
                         ,ordered_IDs$ID_order_phi_s_only_down
                         ,ordered_IDs$ID_order_p_s_only_down
                         ,ordered_IDs$ID_order_rho_r_only_down
                         ,ordered_IDs$ID_order_phi_r_only_down
                         ,ordered_IDs$ID_order_p_r_only_down
                     ))

forestplot_for_category(IDs, "./plots/forest/forestplot_only_down.png")

IDs <- as.character(c(ordered_IDs$ID_order_rho_s_only_up
                          ,ordered_IDs$ID_order_phi_s_only_up
                          ,ordered_IDs$ID_order_p_s_only_up
                          ,ordered_IDs$ID_order_rho_r_only_up
                          ,ordered_IDs$ID_order_phi_r_only_up
                          ,ordered_IDs$ID_order_p_r_only_up
                     ))

forestplot_for_category(IDs, "./plots/forest/forestplot_only_up.png")

IDs <- as.character(c(ordered_IDs$ID_order_concordantDown_s
                          ,ordered_IDs$ID_order_concordantDown_r
                          ,ordered_IDs$ID_order_concordantDown_b_one_one
                          ,ordered_IDs$ID_order_concordantDown_b_one_several
                          ,ordered_IDs$ID_order_concordantDown_b_several_one
                          ,ordered_IDs$ID_order_concordantDown_b_several_several
                     ))

forestplot_for_category(IDs, "./plots/forest/forestplot_concordantDown_one_several.png")

IDs <- as.character(c(ordered_IDs$ID_order_concordantUp_s
                          ,ordered_IDs$ID_order_concordantUp_r
                          ,ordered_IDs$ID_order_concordantUp_b_one_one
                          ,ordered_IDs$ID_order_concordantUp_b_one_several
                          ,ordered_IDs$ID_order_concordantUp_b_several_one
                          ,ordered_IDs$ID_order_concordantUp_b_several_several
                     ))

forestplot_for_category(IDs, "./plots/forest/forestplot_concordantUp_one_several.png")

IDs <- as.character(c(ordered_IDs$ID_order_insuffDown_s
                          ,ordered_IDs$ID_order_insuffDown_r
                          ,ordered_IDs$ID_order_insuffDown_b
                     ))

forestplot_for_category(IDs, "./plots/forest/forestplot_insuffDown.png")

IDs <- as.character(c(ordered_IDs$ID_order_insuffUp_s
                          ,ordered_IDs$ID_order_insuffUp_r
                          ,ordered_IDs$ID_order_insuffUp_b
                     ))

forestplot_for_category(IDs, "./plots/forest/forestplot_insuffUp.png")

IDs <- as.character(c(ordered_IDs$ID_order_suffComp_s
                      ,ordered_IDs$ID_order_suffComp_r
                      ,ordered_IDs$ID_order_suffComp_b
                 ))

forestplot_for_category(IDs, "./plots/forest/forestplot_suffComp.png")

output_dir <- "each_component_values/"

write.csv(interactions$a_s, file = paste0(output_dir, "interactions_a_s.csv"), row.names = TRUE)
write.csv(interactions$a_r, file = paste0(output_dir, "interactions_a_r.csv"), row.names = TRUE)
write.csv(interactions$e_s_l, file = paste0(output_dir, "interactions_e_s_l.csv"), row.names = TRUE)
write.csv(interactions$e_r_r, file = paste0(output_dir, "interactions_e_r_r.csv"), row.names = TRUE)
write.csv(interactions$rho_r, file = paste0(output_dir, "interactions_rho_r.csv"), row.names = TRUE)
write.csv(interactions$rho_s, file = paste0(output_dir, "interactions_rho_s.csv"), row.names = TRUE)

write.csv(interactions$phi_r_r, file = paste0(output_dir, "interactions_phi_r_r.csv"), row.names = TRUE)
write.csv(interactions$phi_s_l, file = paste0(output_dir, "interactions_phi_s_l.csv"), row.names = TRUE)

write.csv(interactions$p_r_r, file = paste0(output_dir, "interactions_p_r_r.csv"), row.names = TRUE)
write.csv(interactions$p_s_l, file = paste0(output_dir, "interactions_p_s_l.csv"), row.names = TRUE)
write.csv(interactions$weights, file = paste0(output_dir, "interactions_weights.csv"), row.names = TRUE)