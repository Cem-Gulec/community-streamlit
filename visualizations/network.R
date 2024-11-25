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

idx_sign <- !is.na(interactions$anno_interactions$sign) & interactions$anno_interactions$sign

# extract cell types
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

colors_interaction_type <- c("engages HSPC" = "magenta",
                             "among immune cells" = "deepskyblue",
                             "engages Ery" = "darkgoldenrod3")

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