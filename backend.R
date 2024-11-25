library(community)
library(data.table) #to read gz file
library(tidyverse)
library(gridExtra)

filterInteractions <- function(comm_result, threshold_log10_cum_weight, threshold_frac_samples_per_condition,
                                threshold_log10_meanexpr_per_condition, verbose = TRUE) {
    # plot the distribution of log10 cumulative interactions weight over the
    # fraction of samples in which the interactions is expressed
    options(repr.plot.height = 10
       ,repr.plot.width = 16)

    # calculate log10 cumulative interactions weights
    comm_result$anno_interactions$log10_cum_weight <- log10(rowSums(comm_result$weights) + 1)
    
    # identify control samples
    idx_control <- comm_result$anno_samples$case_or_control == "control"
    idx_case <- comm_result$anno_samples$case_or_control == "case"
    
    # calculate the fraction of samples expressing the interactions
    comm_result$anno_interactions$frac_samples_controls <- rowSums(comm_result$weights[,idx_control] != 0) / sum(idx_control) 
    comm_result$anno_interactions$frac_samples_cases <- rowSums(comm_result$weights[,idx_case] != 0) / sum(idx_case) 
    
    # set thresholds
    comm_result$thresholds$threshold_log10_cum_weight <- threshold_log10_cum_weight
    comm_result$thresholds$threshold_frac_samples_per_condition = threshold_frac_samples_per_condition
    comm_result$thresholds$threshold_log10_meanexpr_per_condition = threshold_log10_meanexpr_per_condition
    
    
    cumW <- plot_cumW(df = comm_result$anno_interactions, threshold_log10_cum_weight = threshold_log10_cum_weight)
    fracSamp <- plot_fracSamples(df = comm_result$anno_interactions, threshold_frac_samples_per_condition = threshold_frac_samples_per_condition)

    p <- arrangeGrob(fracSamp$ydensity
                                    ,fracSamp$QC_plot
                                    ,fracSamp$blankPlot
                                    ,fracSamp$xdensity
                                    ,ncol=2
                                    ,nrow=2
                                    ,widths=c(2.5, 5.5)
                                    ,heights=c(6.5, 1.5)
                    )
    #ggsave("./plots/above_plot2.png", plot = p, width = 10, height = 8, dpi = 300)

    # arrange plots
    grid.arrange(cumW
                    , p
                    , ncol=2
                    , widths = c(3.5,4.5)
                )
    #ggsave("./plots/above_plot1.png", plot = cumW, width = 10, height = 8, dpi = 300)
    
    
    meanLig_vs_meanRec <- plot_meanLig_vs_meanRec(comm_result$anno_interactions, threshold_log10_meanexpr_per_condition = threshold_log10_meanexpr_per_condition)
    
    #ggsave("./plots/meanlig-vs-meanrec.png", plot = meanLig_vs_meanRec, width = 15, height = 8, dpi = 300)
    
    # filter interactions which did not pass the threshold in any sample
    comm_result$anno_interactions$passed_log10_cum_weight_filter <- comm_result$anno_interactions$log10_cum_weight >
            threshold_log10_cum_weight
    comm_result$anno_interactions$passed_frac_samples_filter <- (comm_result$anno_interactions$frac_samples_controls >
                                                                                threshold_frac_samples_per_condition) | (comm_result$anno_interactions$frac_samples_cases > threshold_frac_samples_per_condition)
    comm_result$anno_interactions$passed_log10_meanexpr_control_filter <- (log10(comm_result$anno_interactions$mean_e_s_l_control +
                                                                                1) > threshold_log10_meanexpr_per_condition) & (log10(comm_result$anno_interactions$mean_e_r_r_control +
                                                                                                                                                1) > threshold_log10_meanexpr_per_condition)
    
    
    comm_result$anno_interactions$passed_log10_meanexpr_case_filter <- (log10(comm_result$anno_interactions$mean_e_s_l_case +
                                                                            1) > threshold_log10_meanexpr_per_condition) & (log10(comm_result$anno_interactions$mean_e_r_r_case +
                                                                                                                                            1) > threshold_log10_meanexpr_per_condition)
    
    comm_result$anno_interactions$passed_log10_meanexpr_per_condition_filter <- comm_result$anno_interactions$passed_log10_meanexpr_control_filter |
            comm_result$anno_interactions$passed_log10_meanexpr_case_filter
    
    # filter anno_interactions
    comm_result$anno_interactions$passed_QC_filter <- (comm_result$anno_interactions$passed_log10_cum_weight_filter &
                                                                comm_result$anno_interactions$passed_frac_samples_filter & comm_result$anno_interactions$passed_log10_meanexpr_per_condition_filter)
    
    samples <- names(comm_result$per_sample_anno_interactions)
    
    if (verbose) {
            cat(
                sum(!(comm_result$anno_interactions$passed_log10_cum_weight_filter & comm_result$anno_interactions$passed_frac_samples_filter)),
                "out of", nrow(comm_result$weights), "interactions do not pass the thresholds for log10 cumulative interactions weight >",
                threshold_log10_cum_weight, "and fraction of expressing samples >", threshold_frac_samples_per_condition,
                ".\nAlso", sum(!comm_result$anno_interactions$passed_log10_meanexpr_per_condition_filter),
                "interactions didn't pass the discrepancy filter.\nIn total,", sum(!comm_result$anno_interactions$passed_QC_filter),
                "bad quality interactions will be removed and", sum(comm_result$anno_interactions$passed_QC_filter),
                "good quality interactions will remain."
            )
    }
    
    return(comm_result)
}

input_dir <- "input_data/"
output_dir <- "computed_results/"

data("LR_database")

# # load counts
#print("load counts")
#counts <- fread(paste0(input_dir,"counts_lognorm.csv.gz"), header = TRUE)
#counts <- as.data.frame(counts)
#rownames(counts) <- counts$gene_symbol
#counts <- counts[,-1]

# load cell annotation
#print("load cell annotation")
#anno_cells <- read.table(paste0(input_dir,"anno_cells_norm.txt")
#                         ,sep = "\t"
#                         # ,row.names = 1
#                         ,header = TRUE
#                         )

# load sample annotation
print("load sample annotation")
anno_samples <- read.table(paste0(input_dir,"anno_samples_norm.txt")
                           ,sep = "\t"
                           # ,row.names = 1
                           ,header = TRUE
                           )

#colnames(counts) <- anno_cells$cell_ID
#rownames(anno_cells) <- anno_cells$cell_ID

threshold_celltype_size <- 6
threshold_nr_active_cells <- 6
threshold_expr <- 0.1

# Renaming the cell_ID.1 column in anno_cells to "cell_ID"
#colnames(anno_cells)[colnames(anno_cells) == "cell_ID.1"] <- "cell_ID"

counts <- read.csv(paste0(input_dir, "toy_counts.csv"), row.names = 1, check.names = FALSE)
cell_annot <- read.csv(paste0(input_dir, "toy_cell_annot.csv"), row.names = 1, check.names = FALSE)


print("calculate communication")
interactions <- calculate_communication(counts = counts
                                       ,anno_samples = anno_samples
                                       ,anno_cells = cell_annot
                                       ,threshold_celltype_size = threshold_celltype_size
                                       ,threshold_nr_active_cells = threshold_nr_active_cells
                                       ,threshold_expr = threshold_expr
                                       ,lrp_database = LR_database
                                       )

print("calculate general statistics")

interactions <- general_stat(comm_result = interactions
                                   ,verbose = FALSE
)

threshold_log10_cum_weight <-  0.01
threshold_frac_samples_per_condition <-  0.6
threshold_log10_meanexpr_per_condition <- 0.02

interactions_df <- interactions$anno_interactions
#write.table(interactions_df, file = "interactions_QC.txt", sep = "\t", row.names = FALSE, quote = FALSE)

print("filter weak interactions")
options(repr.plot.height = 10
       ,repr.plot.width = 16)
interactions <- filterInteractions(comm_result = interactions
                             ,threshold_frac_samples_per_condition = threshold_frac_samples_per_condition
                             ,threshold_log10_cum_weight = threshold_log10_cum_weight
                             ,threshold_log10_meanexpr_per_condition = threshold_log10_meanexpr_per_condition
)

# differential communication
threshold_log2FC <- 1
threshold_fdr <- 0.1

print("calculate differential communication")
interactions <- test_diff(comm_result = interactions
                          ,threshold_fdr = threshold_fdr
                          ,which_test = "t-test"
                          ,threshold_log2FC = threshold_log2FC
                          
                         )

diff_exp_interactions <- subset(interactions$anno_interactions, sign == TRUE)
#write.table(diff_exp_interactions, file = "diff_exp_interactions.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# calculate interactions of the individual components
interactions <- interaction_classes(interactions
                   ,threshold = threshold_log2FC)

#dir.create(output_dir)

#write.csv(interactions$weights,paste0(output_dir,"community_weights.csv"))
#write.csv(interactions$anno_interactions,paste0(output_dir,"community_anno_interactions.csv"))

#print("save interactions.RData")
#save(interactions, file = paste0(output_dir,"interactions.RData"))

#dir.create("visualizations/")
#write.table(interactions$thresholds, file = "visualizations/thresholds.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(interactions$anno_interactions, file = "visualizations/anno_interactions.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(interactions$anno_cells, file = "visualizations/anno_cells.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(interactions$anno_samples, file = "visualizations/anno_samples.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(interactions$weights, file = "visualizations/weights.txt", sep = "\t", row.names = FALSE, quote = FALSE)