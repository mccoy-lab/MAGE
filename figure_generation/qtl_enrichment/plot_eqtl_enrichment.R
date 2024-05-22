#!/usr/bin/env R

###########################################################
#                                                         #
# cis-eQTLs functional annotation enrichment analysis     #
#                                                         #
########################################################### 

# Ensure the packages and libraries below are installed
library(data.table)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(pheatmap)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

plot_foldenrich_fullheatmap <- function(local_dir, infile, annotation_file, output_fileprefix) {
  # Load input file
  input_file <- file.path(local_dir, infile)

  read_df <- fread(input_file, header = TRUE) %>%
    as.data.frame() %>%
    separate(Bed_File, c("celltype", "chromState", "txt"), sep = "[.]")

  read_df <- read_df[complete.cases(read_df), ]

  # Compute Enrichment: observed/expected and take log2
  read_df <- read_df %>%
    mutate(
      inbed = as.numeric(read_df[, 4]),
      expected = as.numeric(read_df[, 5]),
      foldEnrich = log2(inbed / expected)
    )

  read_tf <- dcast(read_df, celltype ~ chromState, value.var = "foldEnrich")

  data <- read_tf[, -1]
  mat_data <- as.matrix(data)
  mat_data <- apply(mat_data, 2, function(x) as.numeric(as.character(x)))
  mat_data[is.na(mat_data)] <- 0

  # Remove numeric prefixes and underscores from column names
  cleaned_colnames <- gsub("^[0-9]+_", "", colnames(data))
  colnames(mat_data) <- cleaned_colnames
  rownames(mat_data) <- read_tf$celltype

  roadmap_df <- fread(annotation_file, sep = "\t")
  anno_df <- left_join(data.frame(EpigenomeID = rownames(mat_data)), roadmap_df, by = "EpigenomeID")

  # Sorting anno_df based on rownames of mat_data to ensure consistent ordering
  anno_df <- anno_df[match(rownames(mat_data), anno_df$EpigenomeID), ]

  # Define color mapping for groups
  group_colors <- unique(anno_df[, c("Group", "Color")])
  color_mapping <- setNames(group_colors$Color, group_colors$Group)

  # Creating the row annotation object
  row_anno <- rowAnnotation(CellType = anno_df$Group, 
                            col = list(CellType = color_mapping),
                            width = unit(1, "cm"),
                            #annotation_name_side = "left",
                            annotation_legend_param = list(title = "CellType", 
                                                           title_gp = gpar(fontsize = 10, fontface = "bold"),
                                                           labels_gp = gpar(fontsize = 8)))

  
  #output_file <- file.path(local_dir, "GNR_fullHeatmap_Log2FoldEnrich_FinalVersionCall.pdf")
  output_file <- file.path(local_dir, paste0(output_fileprefix, "_fullHeatmap_Log2FoldEnrich_FinalVersionCall.pdf"))
  pdf(output_file)

  # Set custom fix color
  my_colors <- colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
  
  ht <- Heatmap(mat_data, name = "foldEnrich", 
        row_names_side = "right",
        col = my_colors,  # Apply your color scale
        #row_names_gp = gpar(fontsize = 3),
        column_names_gp = gpar(fontsize = 9.5),
        show_row_names = FALSE,
        left_annotation = row_anno,
        heatmap_legend_param = list(color_bar = "continuous", 
                                    title = "log2(foldEnrich)",
                                    title_gp = gpar(fontsize = 9, fontface = "bold" ),  # Adjust title font size here
                                    labels_gp = gpar(fontsize = 8),  # Adjust label font size here
                                    at = c(-4, -2, 0, 2, 4), 
                                    labels = c("<=-4", "-2", "0", "2", ">=4")))

  draw(ht)

  # Capture the clustering order of rows and columns
  row_order <- row_order(ht)
  column_order <- column_order(ht)

  dev.off()

  # Return the clustering order of rows and columns
  return(list(row_order = row_order, column_order = column_order))

  #return(ht_list)

}

infile="StatisticSummaryFile_chromHMM_15state_tl_topsnps.txt"
annotation_file="roadmap_annotation.txt"


local_dir="gregor/fulluniq_redo"
output_fileprefix="MAGE"
foldEnrich_heatmap_mage <- plot_foldenrich_fullheatmap(local_dir, infile, annotation_file, output_fileprefix)

# Extract row and column order from the drawing details
row_order_mage <- foldEnrich_heatmap_mage$row_order
col_order_mage  <- foldEnrich_heatmap_mage$column_order

local_dir="gregor/fulluniq_mageqtl"
output_fileprefix="GNR"
foldEnrich_heatmap_gnr  <- plot_foldenrich_fullheatmap(local_dir, infile, annotation_file, output_fileprefix)

# Extract row and column order from the drawing details
row_order_gnr <- foldEnrich_heatmap_gnr$row_order
col_order_gnr  <- foldEnrich_heatmap_gnr$column_order



#########################################
#                                       #
# For Multiple Annotations              # 
#                                       #
#########################################


# Sorting anno_df based on rownames of mat_data to ensure consistent ordering
anno_df <- anno_df[match(rownames(mat_data), anno_df$EpigenomeID), ]

# Define color mapping for groups
group_colors <- unique(anno_df[, c("Group", "Color")])
color_mapping <- setNames(group_colors$Color, group_colors$Group)

# Creating the row annotation object
row_anno1 <- rowAnnotation(Group = anno_df$Group, 
                          col = list(Group = color_mapping),
                          name = "CellType",
                          width = unit(1, "cm"))

row_anno2 <- rowAnnotation(Group = anno_df$Group, 
                          col = list(Group = color_mapping),
                          name = "CellTypeV2",
                          width = unit(5, "cm"))

# Creating a spacer annotation
spacer_anno <- rowAnnotation(df = data.frame(Spacer = rep("", nrow(anno_df))),
                             width = unit(0.001, "cm"), 
                             show_legend = FALSE)

# Combine the annotation objects into a single list
combined_row_anno <- c(row_anno1,spacer_anno,row_anno2)


# Define your color scale
my_colors <- colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

Heatmap(mat_data, name = "foldEnrich", 
        row_names_side = "right",
        col = my_colors,  # Apply your color scale
        show_row_names = FALSE,
        #row_names_gp = gpar(fontsize = 3),
        left_annotation = combined_row_anno,
        #right_annotation = combined_row_anno,
        heatmap_legend_param = list(color_bar = "continuous", 
                                    title = "log2(foldEnrich)", 
                                    at = c(-4, -2, 0, 2, 4), 
                                    labels = c("-4", "-2", "0", "2", "4")))


Heatmap(mat_data, name = "foldEnrich", 
        row_names_side = "right",
        show_row_names = FALSE,
        #row_names_gp = gpar(fontsize = 3),
        left_annotation = combined_row_anno,
        #right_annotation = combined_row_anno,
        heatmap_legend_param = list(color_bar = "continuous", 
                                    title = "log2(foldEnrich)", 
                                    at = c(-4, -2, 0, 2, 4), 
                                    labels = c("-4", "-2", "0", "2", "4")))



###############################
#                             #
# BINNED HEATMAP for Deciles  #
#                             #
###############################


# Load libraries
library(data.table)
library(tidyr)
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(circlize)

generate_decile_heatmaps <- function(elements, local_dir, annotation_file) {

  # Custom color scaling:
  my_colors <- colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

  # Step 3: Generate Heatmap
  for (cell_id in elements) {
    print(paste("Processing ...", cell_id))
    infile <- paste0(cell_id, "-StatisticSummaryFile_chromHMM_15state_tl_topsnps.txt")
    
    input_file <- file.path(local_dir, infile)
    
    read_df <- fread(input_file, header = FALSE)
    names(read_df) <- c("Bed_File", "InBed_Index_SNP", "ExpectNum_of_InBed_SNP",  "ID")
    
    read_df <- read_df %>%
      as.data.frame() %>%
      separate(Bed_File, c("celltype", "chromState", "txt"), sep = "[.]") %>%
      separate(ID, c("PValue", "DecileID"), sep = "[ ]") %>%
      filter(complete.cases(.)) %>%
      mutate(
        inbed = as.numeric(InBed_Index_SNP),
        expected = as.numeric(ExpectNum_of_InBed_SNP),
        foldEnrich = log2(inbed / expected)
      )
    
    read_tf <- dcast(read_df, DecileID ~ chromState, value.var = "foldEnrich")
    data <- read_tf[, -1]
    mat_data <- apply(as.matrix.noquote(data), 2, as.numeric)
    mat_data[is.na(mat_data)] <- 0
    
    # Clean up column names
    cleaned_colnames <- gsub("^[0-9]+_", "", colnames(data))
    colnames(mat_data) <- cleaned_colnames
    rownames(mat_data) <- read_tf$DecileID
    
    order_vec_rows <- rev(c("decile1", "decile2", "decile3", "decile4", "decile5", 
                            "decile6", "decile7", "decile8", "decile9", "decile10"))
    
    # Order the matrix rows
    mat_ordered <- mat_data[order_vec_rows, ]
    
    # Within for loop and after generating mat_ordered
    # Define color scale
    # my_colors <- colorRamp2(c(min(mat_ordered), 0, max(mat_ordered)), c("blue", "white", "red"))

    output_file <- file.path(local_dir, paste0(cell_id, "_MAGE_Log2FoldEnrichFinalversion.pdf"))
    pdf(output_file, width = 10, height = 8)

    ht <- Heatmap(mat_ordered, 
            name = "foldEnrich", 
            col = my_colors,
            show_row_names = TRUE, 
            show_column_names = TRUE, 
            cluster_rows = FALSE, 
            cluster_columns = TRUE, 
            row_names_side = "left", 
            column_title = "Chromatin States", 
            row_title = "Deciles",
            #heatmap_legend_param = list(title = "log2(foldEnrich)", at = c(-5, -3, -2, 0, 2, 3, 5), labels = c("-5", "-3", "-2", "0", "2", "3", "5"))
            #heatmap_legend_param = list(title = "log2(foldEnrich)", at = c(-4, -2, 0, 2, 4), labels = c("≤-4", "-2", "0", "2", "≥4"))
            heatmap_legend_param = list(title = "log2(foldEnrich)", at = c(-4, -2, 0, 2, 4), labels = c("<=-4", "-2", "0", "2", ">=4"))
    )

    draw(ht)

    dev.off()

  }

}


# Define your elements vector
elements <- c("E116", "E034", "E032", "E051", "E046")

# Define the local directory
local_dir <- "gregor/decileuniq_redo"

# Define the annotation file path (change this to your actual path)
annotation_file <- "roadmap_annotation.txt"

# Call the function
generate_decile_heatmaps(elements, local_dir, annotation_file)



#########################################################
#                                                       #
# TFBS VOLCANO PLOT                                     #
#                                                       #
#########################################################

# Pvalue Vs Fold(O/E) Enrichment Scatter Plot

plot_tf_volcano <- function(local_dir, infile) {

  # Load required libraries
  library(data.table)
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)

  input_file <- file.path(local_dir, infile)
  read_df <- fread(input_file, header = TRUE) %>% 
    as.data.frame() %>%
    separate(Bed_File, c("TF", "txt"),  sep = "[.]")


  bonferroni_threshold <- -log10(0.0001/nrow(read_df))

  rownames(read_df) <- read_df$TF

  enrich <- read_df[,3]/read_df[,4]
  read_df <- read_df %>%
    mutate(
      enrich = log2(enrich),
      log10PValue = -log10(PValue)
    ) %>%
    arrange(desc(log10PValue)) %>%
    filter(!is.na(PValue))

  options(ggrepel.max.overlaps = 100)

  label_data1 <- head(read_df[order(-read_df$log10PValue),], 15)
  label_data2 <- head(read_df[order(-read_df$enrich),], 15)
  label_data <- bind_rows(label_data1, label_data2)

  data1 <- read_df[read_df$TF == "POLR2A",]
  label_data <- bind_rows(label_data, data1)

  output_file <- file.path(local_dir, "test_MAGE_TFBS_VolcanoPlot_Log2Enrich.pdf")

  read_df <- read_df %>%
    mutate(point_color = ifelse(log10PValue < bonferroni_threshold, "grey", "#FF4040")) ##FF6666 (light red)

  plot <- ggplot(read_df, aes(x = enrich, y = log10PValue)) +
      geom_point(aes(color = point_color), size = 1.7) +
      scale_color_identity() +
      geom_vline(colour = "dodgerblue", xintercept = log2(2), linetype = 3, size = 0.9) +
      geom_vline(colour = "black", xintercept = log2(1), linetype = 3) +
      geom_hline(colour = "dodgerblue", yintercept = bonferroni_threshold, linetype = 3, size = 0.9) +
      geom_hline(colour = "black", yintercept = -log10(1), linetype = 3) +
      labs(y = "-log10(PValue)", x = "log2(Fold-enrichment)", title = "") +
      theme_minimal() + 
      theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
      geom_text_repel(
        data = label_data,
        aes(label = TF),
        segment.linetype = 6, 
        box.padding = 1,
        show.legend = FALSE,
        size = 2.7
      )

  ggsave(output_file, plot, width=7, height=7)

}

local_dir <- "gregor/fulluniq_mageqtl"
local_dir="gregor/fulluniq_redo"

infile <- "StatisticSummaryFile_wgEncodeRegTfbsCluster_tl_topsnps.txt"
plot_tf_volcano(local_dir, infile)



#############################################################
#                                                           #
# DNASE PROM ENH DYADIC                                     #
#                                                           #
#############################################################

# VIOLINBOX PLOT effect size distribution vs chromatin state
# FULL Dataset and not decile data


generate_violinbox_fullplot_dnase <- function(file_path, save_dir) {

  # Load required libraries
  library(ggplot2)
  library(data.table)
  library(dplyr)


  # Extract EID from the file name
  eid <- gsub(".*(E[0-9]+).*", "\\1", file_path)
    
  # Read in the data
  data <- fread(file_path, sep="\t")

  # Subset the data based on your columns of interest and take absolute values
  data_subset <- data.frame(log2_aFC = abs(as.numeric(data$V8)),
                            ChromatinState = data$V18)
                            #ColorCode = data$V18)

  # Calculate median for each ChromatinState
  median_values <- data_subset %>% 
    group_by(ChromatinState) %>% 
    summarise(medianValue = median(log2_aFC)) %>%
    arrange(medianValue)

  # Factorize the ChromatinState column based on order of median values
  data_subset$ChromatinState <- factor(data_subset$ChromatinState, levels = median_values$ChromatinState)

  # Create a color mapping for the chromatin states
  # color_mapping <- data_subset %>% distinct(ChromatinState, ColorCode)

  chromHMM_colors <- data.frame(
    ChromatinState = c(
      "prom", "enh", "dyadic"
    ),
    ColorCode = c(
      "#FF0000",  "#FFFF00", "#8B008B"
    ),
    stringsAsFactors = FALSE
  )

  # Calculate number of data points for each ChromatinState
  data_count <- data_subset %>% 
    group_by(ChromatinState) %>% 
    summarise(count = n())

  print(data_count)

  output_file <- file.path(save_dir, paste0(eid, "_MAGE_DNase_ViolinBoxplot.pdf"))

  # Create the plot
  plot <- ggplot(data_subset, aes(x = ChromatinState, y = log2_aFC, fill = ChromatinState)) +
        geom_violin(scale = "width", trim = TRUE, width = 0.6) + 
        geom_boxplot(width = 0.05, fill = "black", outlier.shape = NA) + 
        stat_summary(fun=median, geom="point", shape=23, size=2.5, color="black", fill="white") +  # This line adds the dot for median
        geom_text(data=data_count, aes(label=paste0("(n=", count, ")"), y=Inf), vjust=2, size=3) + # Number of data points
        
        scale_fill_manual(values = chromHMM_colors$ColorCode, 
                        breaks = chromHMM_colors$ChromatinState,
                        labels = chromHMM_colors$ChromatinState) + 

        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10),
              legend.key = element_blank(),
              legend.key.size = unit(0.6, "cm"),  # Adjust this value to decrease or increase the legend size
              legend.text = element_text(size = 8)) +  # Adjusting text size of the legend) +  # This removes the default square background) +  # Hide the default legend key

        #scale_shape_manual(values = rep(21, 15)) +
        guides(fill = guide_legend(override.aes = list(shape = 23, size = 0.25))) + 

        labs(title="Distribution of aFC per Chromatin State",
             x="Chromatin State Category", y="|log2(aFC)|")

  ggsave(output_file, plot = plot)

}

# Specify the directory where the files are and where to save the plots
output_directory <- "gregor/fulluniq/DNase_promenhdyadic_annotations"
output_directory <- "gregor/fulluniq_redo/DNase_promenhdyadic_annotations"

# Ensure the output directory exists; if not, create it
if (!dir.exists(output_directory)) {
  dir.create(output_directory)
}

# Get list of all files in the directory
file_list <- list.files("gregor/fulluniq/DNase_promenhdyadic_annotations", pattern = "dnaseintersected.bed$", full.names = TRUE)

# Apply the generate_plot function to each file, passing in the output directory for saving plots
lapply(file_list, generate_violinbox_fullplot_dnase, save_dir = output_directory)



#################################################
#                                               #
# GTEx Non-replicating (GNR) credible sets      #
#                                               #
#################################################


plot_fullheatmap_gnr <- function(local_dir, infile, annotation_file) {

  # Ensure necessary libraries are loaded
  library(data.table)
  library(tidyr)
  library(dplyr)
  library(pheatmap)
  library(reshape2)

  input_file <- file.path(local_dir, infile)

  read_df <- fread(input_file, header = T) %>% 
    as.data.frame() %>%
    separate(Bed_File, c("celltype", "chromState", "txt"), sep = "[.]")

  read_df <- read_df[complete.cases(read_df), ]

  # Compute Enrichment : observed/expected and take log2
  read_df <- read_df %>% 
    mutate(
      inbed = as.numeric(read_df[,4]),
      expected = as.numeric(read_df[,5]),
      foldEnrich = log2(inbed / expected)
    )

  read_tf <- dcast(read_df, celltype ~ chromState, value.var = "foldEnrich")

  data <- read_tf[, 2:ncol(read_tf)]
  mat_data <- apply(as.matrix.noquote(data), 2, as.numeric)
  mat_data[is.na(mat_data)] <- 0


  # Remove numeric prefixes and underscores from column names
  cleaned_colnames <- gsub("^[0-9]+_", "", colnames(data))
  colnames(mat_data) <- cleaned_colnames
  rownames(mat_data) <- read_tf$celltype

  roadmap_df <- fread(annotation_file, sep = "\t")
  anno_df <- left_join(data.frame(EpigenomeID = rownames(mat_data)), roadmap_df, by = "EpigenomeID")

  annotation_row <- data.frame(
    CellType = anno_df$Group
  )

  rownames(annotation_row) = rownames(mat_data)
  mat_colors <- list(
    CellType = anno_df$Color
  )
  names(mat_colors$CellType) <- anno_df$Group

  output_file <- file.path(local_dir, "GNR_chromHMM_15state_tl_topsnps_Log2FoldEnrich.pdf")

  pdf(output_file)

  pheatmap_plot <- pheatmap(mat_data,
    fontsize_row = 0.0000000001,
    fontsize_col = 9,
    annotation_row = annotation_row,
    annotation_colors = mat_colors,
    annotation_legend = FALSE,
    legend_labels = "Log2(Fold-Enrich)",
    main = ""
  )

  dev.off()

  return(list(row_order = pheatmap_plot$tree_row$order, col_order = pheatmap_plot$tree_col$order))

}


infile="StatisticSummaryFile_chromHMM_15state_tl_topsnps.txt"
annotation_file="roadmap_annotation.txt"

local_dir="gregor/fulluniq_mageqtl"

plot_fullheatmap_gnr(local_dir, infile, annotation_file)

fullheatmap_result_gnr <- plot_fullheatmap_gnr(local_dir, infile, annotation_file)

# Extract row and column order
row_order_gnr <- fullheatmap_result_gnr$row_order
column_order_gnr <- fullheatmap_result_gnr$col_order

# # Extract row and column order
# row_order_mage <- fullheatmap_mage_result$row_order
# column_order_mage <- fullheatmap_mage_result$col_order


#########################################
#                                       #
# DNASE Hypersensitivity Sites (DHS)    #
#                                       #
#########################################

# DHS analysis with 53 epigenomes

plot_DNase_pvalenrichment <- function(local_dir, infile, annotation_file) {

  # Load libraries
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(tidyr)

  input_file <- file.path(local_dir, infile)
  read_df <- fread(input_file, header = TRUE) %>% 
    as.data.frame() %>%
    separate(Bed_File, c("celltype", "txt"), sep = "[-]")
  
  # Roadmap annotation
  roadmap_df <- fread(annotation_file, sep = "\t")
  
  # First sync up with the matrix data
  anno_df <- left_join(read_df, roadmap_df, by = c("celltype" = "EpigenomeID"))
  anno_df$log10PValue <- -log10(anno_df$PValue)
  
  # Assign colors based on factor levels of groups
  anno_df[anno_df$celltype == "Control", "Group"] <- "Control"
  anno_df[anno_df$celltype == "Control", "Color"] <- "black"
  
  df1 <- data.frame(Group_level = levels(factor(anno_df$Group)))
  df2 <- anno_df[, c("Group", "Color")]
  df3 <- inner_join(df1, df2, by = c("Group_level" = "Group"))
  color_df <- df3 %>% distinct(Group_level, .keep_all = TRUE)
  color_df[color_df$Group_level == "ENCODE2012", "Color"] <- "dodgerblue"
  
  label_data <- anno_df[order(-anno_df$log10PValue),][1:5,]
  
  output_file <- file.path(local_dir, "GNR_DNase_DHSn53_pvalEnrichment.pdf")
  #pdf(file.path(local_dir, output_file), width = 9, height = 9)
  
  anno_df <- anno_df[complete.cases(anno_df),]
  
  plot <- ggplot(anno_df, aes(x = reorder(celltype, log10PValue), y = log10PValue)) +
    geom_point(size = 2.5, aes(color = Group)) + 
    scale_color_manual(values = color_df$Color) +
    geom_segment(aes(xend = celltype, yend = 0), size = 0.3, linetype = 1) +
    labs(y = "-log10(PValue)", x = "CellType", title = "SNP Enrichment in DHS(n=53)") +
    theme_minimal() + 
    theme(axis.text.x = element_text(size = 5.5, angle = 90, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    geom_text_repel(data = label_data,
                    aes(label = EpigenomeMnemonic), 
                    box.padding = 1,
                    max.overlaps = 100,
                    show.legend = FALSE, # this removes the 'a' from the legend
                    size = 2
    )
  
  ggsave(output_file, plot, width=7, height=7)

}

# Usage
annotation_file <- "roadmap_annotation.txt"
local_dir <- "gregor/fulluniq_mageqtl"

infile <- "StatisticSummaryFile_DNase_tl_topsnps.txt"
# plot_DNase_log2foldEnrichment(local_dir, infile, annotation_file)


plot_DNase_log2foldEnrichment <- function(local_dir, infile, annotation_file, output_fileprefix) {

  # Load libraries
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(tidyr)

  input_file <- file.path(local_dir, infile)
  
  read_df <- fread(input_file, header = TRUE) %>% 
    as.data.frame() %>%
    separate(Bed_File, c("celltype", "txt"), sep = "[-]")

  read_df <- read_df %>% 
    filter(complete.cases(.)) %>%
    mutate(
      inbed = as.numeric(InBed_Index_SNP),
      expected = as.numeric(ExpectNum_of_InBed_SNP),
      log2foldenrich = log2(inbed / expected)
    )
  
  # Roadmap annotation
  roadmap_df <- fread(annotation_file, sep = "\t")
  
  # First sync up with the matrix data
  anno_df <- left_join(read_df, roadmap_df, by = c("celltype" = "EpigenomeID"))
  #anno_df$log2foldenrich <- -log10(anno_df$PValue)
  
  # Assign colors based on factor levels of groups
  anno_df[anno_df$celltype == "Control", "Group"] <- "Control"
  anno_df[anno_df$celltype == "Control", "Color"] <- "black"
  
  df1 <- data.frame(Group_level = levels(factor(anno_df$Group)))
  df2 <- anno_df[, c("Group", "Color")]
  df3 <- inner_join(df1, df2, by = c("Group_level" = "Group"))
  color_df <- df3 %>% distinct(Group_level, .keep_all = TRUE)
  color_df[color_df$Group_level == "ENCODE2012", "Color"] <- "dodgerblue"
  
  label_data <- anno_df[order(-anno_df$log2foldenrich),][1:5,]
  
  output_file <- file.path(local_dir, paste0(output_fileprefix, "_DNase_DHSn53_log2foldEnrich.pdf"))
  #pdf(file.path(local_dir, output_file), width = 9, height = 9)
  
  anno_df <- anno_df[complete.cases(anno_df),]
  
  plot <- ggplot(anno_df, aes(x = reorder(celltype, log2foldenrich), y = log2foldenrich)) +
    geom_point(size = 2.5, aes(color = Group)) + 
    scale_color_manual(values = color_df$Color) +
    geom_segment(aes(xend = celltype, yend = 0), size = 0.3, linetype = 1) +
    labs(y = "log2(Fold-enrichment)", x = "CellType", title = "SNP Enrichment in DHS(n=53)") +
    theme_minimal() + 
    theme(axis.text.x = element_text(size = 5.5, angle = 90, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    geom_text_repel(data = label_data,
                    aes(label = EpigenomeMnemonic), 
                    box.padding = 1,
                    segment.linetype = 6,
                    max.overlaps = 100,
                    show.legend = FALSE, # this removes the 'a' from the legend
                    size = 2.9
    )
  
  ggsave(output_file, plot, width=8, height=8)

}

# Usage
annotation_file <- "roadmap_annotation.txt"
local_dir <- "gregor/fulluniq_redo"

output_fileprefix = "MAGE"
infile <- "StatisticSummaryFile_DNase_tl_topsnps.txt"
plot_DNase_log2foldEnrichment(local_dir, infile, annotation_file, output_fileprefix)

# For GTEx non replicating crediblesets
local_dir <- "gregor/fulluniq_mageqtl"

infile <- "StatisticSummaryFile_DNase_tl_topsnps.txt"
output_fileprefix = "GNR"
plot_DNase_log2foldEnrichment(local_dir, infile, annotation_file, output_fileprefix)


#######################################################
#                                                     #  
# TFBS VOLCANO PLOT                                   #
#                                                     #
#######################################################

# ENCODE TFBS Clustered CATEGORY
# Pvalue Vs Fold(O/E) Enrichment Scatter Plot

plot_tf_volcano <- function(local_dir, infile) {
  library(data.table)
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)

  input_file <- file.path(local_dir, infile)
  read_df <- fread(input_file, header = TRUE) %>% 
    as.data.frame() %>%
    separate(Bed_File, c("TF", "txt"),  sep = "[.]")


  bonferroni_threshold <- -log10(0.05/nrow(read_df))

  rownames(read_df) <- read_df$TF

  enrich <- read_df[,3]/read_df[,4]
  read_df <- read_df %>%
    mutate(
      enrich = log2(enrich),
      log10PValue = -log10(PValue)
    ) %>%
    arrange(desc(log10PValue)) %>%
    filter(!is.na(PValue))

  options(ggrepel.max.overlaps = 100)

  label_data1 <- head(read_df[order(-read_df$log10PValue),], 15)
  label_data2 <- head(read_df[order(-read_df$enrich),], 18)
  
  # Remove duplicate labels after checking the figures
  label_data2 <- label_data2 %>%
    filter(!(TF %in% c("PAX5", "EED", "POLR2A")))

  label_data <- bind_rows(label_data1, label_data2)

  #data1 <- read_df[read_df$TF == "POLR2A",]
  # label_data <- bind_rows(label_data, data1)

  output_file <- file.path(local_dir, "GNR_TFBS_VolcanoPlot_Log2Enrich.pdf")

  read_df <- read_df %>%
    mutate(point_color = ifelse(log10PValue < bonferroni_threshold, "grey", "#FF4040")) ##FF6666 (light red)

  plot <- ggplot(read_df, aes(x = enrich, y = log10PValue)) +
      geom_point(aes(color = point_color), size = 1.7) +
      scale_color_identity() +
      geom_vline(colour = "dodgerblue", xintercept = log2(2), linetype = 3, size = 0.9) +
      geom_vline(colour = "black", xintercept = log2(1), linetype = 3) +
      geom_hline(colour = "dodgerblue", yintercept = bonferroni_threshold, linetype = 3, size = 0.9) +
      geom_hline(colour = "black", yintercept = -log10(1), linetype = 3) +
      labs(y = "-log10(PValue)", x = "log2(Fold-enrichment)", title = "") +
      theme_minimal() + 
      theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
      geom_text_repel(
        data = label_data,
        aes(label = TF),
        segment.linetype = 6, 
        box.padding = 1,
        show.legend = FALSE,
        size = 2.7
      )

  ggsave(output_file, plot, width=7, height=7)

}

local_dir="gregor/fulluniq_redo"
local_dir <- "gregor/fulluniq_mageqtl"

infile <- "StatisticSummaryFile_wgEncodeRegTfbsCluster_tl_topsnps.txt"
plot_tf_volcano(local_dir, infile)


#############################################################
#                                                           #
# VIOLINBOX PLOT aFC distribution vs Chromatin state        #
#                                                           #
#############################################################

# FULL Dataset and not decile

generate_violinbox_fullplot_chromstates <- function(file_path, save_dir) {

  # Load required libraries
  library(ggplot2)
  library(data.table)
  library(dplyr)

  # Extract EID from the file name
  eid <- gsub(".*(E[0-9]+).*", "\\1", file_path)
  
  # Read in the data
  data <- fread(file_path, sep="\t")
  
  # Remove numeric prefixes from ChromatinState
  data$V13 <- gsub("^[0-9]+_", "", data$V13)
  
  # Subset the data based on your columns of interest and take absolute values
  data_subset <- data.frame(log2_aFC = abs(as.numeric(data$V8)),
                            ChromatinState = data$V13,
                            ColorCode = data$V18)
  
  # Calculate median for each ChromatinState
  median_values <- data_subset %>%
    group_by(ChromatinState) %>%
    summarise(medianValue = median(log2_aFC), .groups = "drop") %>%
    arrange(medianValue)
  
  # Factorize the ChromatinState column based on order of median values
  data_subset$ChromatinState <- factor(data_subset$ChromatinState, levels = median_values$ChromatinState)
  
  # Create a color mapping for the chromatin states
  chromHMM_colors <- data.frame(
    ChromatinState = c(
      "TssA", "TssAFlnk", "TxFlnk", "Tx", "TxWk", "EnhG", "Enh",
      "ZNF-Rpts", "Het", "TssBiv", "BivFlnk", "EnhBiv", "ReprPC",
      "ReprPCWk", "Quies"
    ),
    ColorCode = c(
      "#FF0000", "#FF4500", "#32CD32", "#008000", "#006400", "#C2E105", "#FFFF00",
      "#66CDAA", "#8A91D0", "#C71585", "#8B008B", "#9400D3", "#A020F0", "#D02090", "#FFFFFF"
    ),
    stringsAsFactors = FALSE
  )
  
  # Calculate number of data points for each ChromatinState
  data_count <- data_subset %>%
    group_by(ChromatinState) %>%
    summarise(count = n(), .groups = "drop")
  
  print(data_count)
  
  output_file <- file.path(save_dir, paste0(eid, "_chromatinstate_ViolinBoxplot_MAGE.pdf"))
  
  # Create the plot
  plot <- ggplot(data_subset, aes(x = ChromatinState, y = log2_aFC, fill = ChromatinState)) +
    geom_violin(scale = "width", trim = TRUE, width = 0.6) +
    geom_boxplot(width = 0.05, fill = "white", outlier.shape = NA) +
    stat_summary(fun=median, geom="point", shape=23, size=1.5, color="black", fill="white") +
    geom_text(data=data_count, aes(label=paste0("(n=", count, ")"), y=Inf), vjust=2, size=2.75) +
    scale_fill_manual(values = chromHMM_colors$ColorCode, 
                      breaks = chromHMM_colors$ChromatinState,
                      labels = chromHMM_colors$ChromatinState) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.key = element_blank(),
          legend.key.size = unit(0.6, "cm"),
          legend.text = element_text(size = 8)) +
    guides(fill = guide_legend(override.aes = list(shape = 23, size = 0.25))) +
    labs(title = "Distribution of |log2(aFC)| per Chromatin State",
         x = "Chromatin State", y = "|log2(aFC)|")
  
  ggsave(output_file, plot = plot, width = 9, height = 9)

}

# Example usage
# Specify the directory where the files are and where to save the plots
output_directory <- "gregor/fulluniq_redo/chromHMM_annotation"

# Ensure the output directory exists; if not, create it
if (!dir.exists(output_directory)) {
  dir.create(output_directory)
}

# # Get list of all files in the directory
file_list <- list.files("gregor/fulluniq_redo/chromHMM_annotation", pattern = "_intersected.bed$", full.names = TRUE)

# Apply the generate_plot function to each file, passing in the output directory for saving plots
lapply(file_list, generate_violinbox_fullplot_chromstates, save_dir = output_directory)



plot_pvalue_fullheatmap <- function(local_dir, infile, annotation_file, output_fileprefix, column_order = NULL, row_order = NULL) {
 
  # Ensure necessary libraries are loaded
  library(data.table)
  library(tidyr)
  library(dplyr)
  library(pheatmap)
  library(reshape2)

  input_file <- file.path(local_dir, infile)

  read_df <- fread(input_file, header = T) %>% 
    as.data.frame() %>%
    separate(Bed_File, c("celltype", "chromState", "txt"), sep = "[.]")

  read_df <- read_df[complete.cases(read_df), ]
  read_tf <- read_df %>% dcast(celltype ~ chromState, value.var = "PValue")

  data <- read_tf[, 2:ncol(read_tf)]
  mat_data <- apply(as.matrix.noquote(data), 2, as.numeric)
  mat_data[mat_data == 0] <- 1e-323 
  mat_data <- -log10(mat_data)
  mat_data[is.na(mat_data)] <- 0

  # Clean up column names
  cleaned_colnames <- gsub("^[0-9]+_", "", colnames(data))
  colnames(mat_data) <- cleaned_colnames
  rownames(mat_data) <- read_tf$celltype

  roadmap_df <- fread(annotation_file, sep = "\t")
  anno_df <- left_join(data.frame(EpigenomeID = rownames(mat_data)), roadmap_df, by = "EpigenomeID")

  annotation_row <- data.frame(
    CellType = anno_df$Group
  )

  rownames(annotation_row) = rownames(mat_data)
  mat_colors <- list(
    CellType = anno_df$Color
  )
  names(mat_colors$CellType) <- anno_df$Group

  quantile_breaks <- function(xs, n = 10) {
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
    breaks[!duplicated(breaks)]
  }

  breaks <- quantile_breaks(mat_data, n = 100)

  #output_file <- file.path(local_dir, "chromHMM_15state_tl_topsnps_w-o.encode.pdf")
  output_file <- file.path(local_dir, paste0(output_fileprefix, "_fullchromHMM_15state_tl_topsnps.pdf"))
  pdf(output_file)

  # Order the matrix data if row_order and/or column_order are provided
  if (!is.null(row_order)) {
    mat_data <- mat_data[row_order, , drop = FALSE]
  }
  if (!is.null(column_order)) {
    mat_data <- mat_data[, column_order, drop = FALSE]
  }


  pheatmap(mat_data,
    #color = colors,
    breaks = breaks,  # Add this line to use quantile breaks
    #fontsize_row = 3,
    fontsize_row = 0.0000000001,
    fontsize_col = 9,
    annotation_row = annotation_row,
    annotation_colors = mat_colors,
    annotation_legend = FALSE,
    legend_labels = "-log10(PValue)",
    main = "",
    cluster_cols = is.null(column_order),  # cluster columns only if column_order is NULL
    cluster_rows = is.null(row_order)  # cluster columns only if column_order is NULL
  )

  dev.off()
}


# Usage

infile="StatisticSummaryFile_chromHMM_15state_tl_topsnps.txt"
annotation_file="roadmap_annotation.txt"

local_dir="gregor/fulluniq_redo"
output_fileprefix = "MAGE_Final"
plot_pvalue_fullheatmap(local_dir, infile, annotation_file, output_fileprefix)

output_fileprefix = "GNR_Final_"
local_dir="gregor/fulluniq_mageqtl"
plot_pvalue_fullheatmap(local_dir, infile, annotation_file, output_fileprefix)




#############################################################
#                                                           #
# Identify regulatory region the SNPs are located at        #
#                                                           #
#############################################################

# install.packages("BiocManager")
# BiocManager::install("GenomicRanges")
library(GenomicRanges)
library(BiocFileCache)

bfc <- BiocFileCache(ask=FALSE)
Gm12878.hmm.file <- bfcrpath(bfc, "wgEncodeBroadHmmGm12878HMM.bed.gz")
Gm12878.hmm <- toGRanges(Gm12878.hmm.file)
Gm12878.hmm

# Creating GRanges object for the genomic ranges
reg_regions <- GRanges(seqnames = seqnames(Gm12878.hmm), 
                       ranges = ranges(Gm12878.hmm), 
                       strand = strand(Gm12878.hmm),
                       V4 = Gm12878.hmm$V4)

# Defining the SNPs
snps <- data.frame(
  rsid = c("rs115070172", "rs4930437", "rs7927381", "rs3082142"),
  chrom = c("chr11", "chr11", "chr11", "chr11"),
  variant_pos = c(67327106, 67329695, 67346743, 67348659)
)

# Creating a GRanges object for the SNPs
snp_ranges <- GRanges(seqnames = snps$chrom, 
                      ranges = IRanges(start = snps$variant_pos, end = snps$variant_pos))

# Finding overlaps
overlaps <- findOverlaps(snp_ranges, reg_regions)

# Getting the results
result <- data.frame(
  SNP = snps$rsid[queryHits(overlaps)],
  RegulatoryRegion = reg_regions$V4[subjectHits(overlaps)],
  Location = seqnames(reg_regions)[subjectHits(overlaps)],
  Range = ranges(reg_regions)[subjectHits(overlaps)]
)

print(result)

# > print(result)
#           SNP  RegulatoryRegion Location Range.start Range.end Range.width
# 1 rs115070172   7_Weak_Enhancer    chr11    67326424  67328824        2401
# 2   rs4930437   2_Weak_Promoter    chr11    67329624  67330224         601
# 3   rs7927381 3_Poised_Promoter    chr11    67346624  67347824        1201
# 4   rs3082142      12_Repressed    chr11    67347824  67350024        2201
#   Range.names
# 1      348001
# 2      348004
# 3      348013
# 4      348014


# Finding first nearest regulatory regions
nearest_hits <- nearest(snp_ranges, reg_regions)

# Calculating distances
distances <- distanceToNearest(snp_ranges, reg_regions)

# Getting the results
nearest_result <- data.frame(
  SNP = snps$rsid,
  NearestRegulatoryRegion = reg_regions$V4[nearest_hits],
  Location = seqnames(reg_regions)[nearest_hits],
  Range = ranges(reg_regions)[nearest_hits],
  Distance = mcols(distances)$distance
)

print(nearest_result)


# Finding second nearest regulatory regions (2nd nearest)
library(GenomicRanges)
library(BiocFileCache)

# Initialize BiocFileCache
bfc <- BiocFileCache(ask = FALSE)

# Download the data
Gm12878.hmm.file <- bfcrpath(bfc, "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmGm12878HMM.bed.gz")

# Convert it to a GRanges object
Gm12878.hmm <- toGRanges(Gm12878.hmm.file)

# Create a GRanges object for the regulatory regions
reg_regions <- GRanges(seqnames = seqnames(Gm12878.hmm), 
                       ranges = ranges(Gm12878.hmm), 
                       strand = strand(Gm12878.hmm),
                       V4 = Gm12878.hmm$V4)

# Define the SNPs
snps <- data.frame(
  rsid = c("rs115070172", "rs4930437", "rs7927381", "rs3082142"),
  chrom = c("chr11", "chr11", "chr11", "chr11"),
  variant_pos = c(67327106, 67329695, 67346743, 67348659)
)

# Create a GRanges object for the SNPs
snp_ranges <- GRanges(seqnames = snps$chrom, 
                      ranges = IRanges(start = snps$variant_pos, end = snps$variant_pos))

# Initialize a data frame to store the results
results <- data.frame(SNP = character(0), NearestRegulatoryRegion = character(0), 
                      Location = character(0), Range = character(0), 
                      Distance = numeric(0), stringsAsFactors = FALSE)

# For each SNP
for (i in 1:length(snp_ranges)) {
  # Find the regulatory region in which the SNP is located
  own_reg <- findOverlaps(snp_ranges[i], reg_regions, type = "within")
  
  # If the SNP is in a regulatory region, proceed
  if (length(own_reg) > 0) {
    # Remove that regulatory region
    other_reg_regions <- reg_regions[-subjectHits(own_reg)]
    
    # Find the nearest regulatory region to the SNP
    nearest_hit <- nearest(snp_ranges[i], other_reg_regions)
    
    if (length(nearest_hit) > 0) {
      # Calculate distance
      distance <- min(start(snp_ranges[i]) - end(other_reg_regions[nearest_hit]), 
                      start(other_reg_regions[nearest_hit]) - end(snp_ranges[i]))
      
      # Add the results to the data frame
      results <- rbind(results, data.frame(SNP = snps$rsid[i], 
                                           NearestRegulatoryRegion = other_reg_regions$V4[nearest_hit], 
                                           Location = seqnames(other_reg_regions)[nearest_hit], 
                                           Range = paste(start(other_reg_regions[nearest_hit]), 
                                                         "-", 
                                                         end(other_reg_regions[nearest_hit]), sep = ""), 
                                           Distance = distance))
    }
  }
}

print(results)

# > print(results)
#           SNP NearestRegulatoryRegion Location             Range Distance
# 1 rs115070172       5_Strong_Enhancer    chr11 67326224-67326424     -882
# 2   rs4930437       4_Strong_Enhancer    chr11 67329224-67329624     -471
# 3   rs7927381            12_Repressed    chr11 67344624-67346624    -2119
# 4   rs3082142       3_Poised_Promoter    chr11 67346624-67347824    -2035


############################################################
#                                                          # 
# Generate legend keys for regulatory states               #
#                                                          #
############################################################

library(ggplot2)

# Set dataframe
data <- data.frame(
  reg_elements = c("15_Repetitive/CNV", "13_Heterochrom/lo", "8_Insulator", "11_Weak_Txn", "7_Weak_Enhancer", "10_Txn_Elongation", "9_Txn_Transition", "2_Weak_Promoter", "1_Active_Promoter", "3_Poised_Promoter", "12_Repressed", "6_Weak_Enhancer", "14_Repetitive/CNV", "5_Strong_Enhancer", "4_Strong_Enhancer"),
  color = c("#F5F5F5", "#F5F5F5", "#0ABEFE", "#99FF66", "#FFFC04", "#00B050", "#00B050", "#FF6969", "#FF0000", "#CF0BC6", "#7F7F7F", "#FFFC04", "#F5F5F5", "#FACA00", "#FACA00"),
  order = as.numeric(sub(".*?([0-9]+).*", "\\1", data$reg_elements))
)

# Extract numeric part and use it to order the data frame
data$numeric_part <- as.numeric(sub(".*?([0-9]+).*", "\\1", data$reg_elements))
data <- data[order(data$numeric_part), ]

# Convert reg_elements to factor with levels in the correct order
data$reg_elements <- factor(data$reg_elements, levels = data$reg_elements)

# Create the plot
plot <- ggplot(data, aes(x = reg_elements, y = 1, fill = color)) + 
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() + 
  theme_minimal() +
  scale_fill_identity() +
  labs(x = "Regulatory Elements", y = "") +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank()) +
  geom_text(aes(y = 1, label = reg_elements), hjust = 0, size = 3) + # Adjusted text size
  theme(plot.margin = unit(c(1, 5, 1, 0), "cm")) # Increased right margin

print(plot)

# Define file path for saving the plot
file_path <- "GM12878_colorkey1.pdf"

# Save the plot as a PDF with increased width
ggsave(file_path, plot, device = "pdf", width = 40, height = 5) # Adjusted plot width



###################################################
#                                                 #
# Code for FUll-Heatmap aligned to Pvalue         #
#                                                 #
###################################################


plot_foldenrich_fullheatmap <- function(local_dir, infile, annotation_file, output_fileprefix) {
  
  # Ensure necessary libraries are loaded
  library(data.table)
  library(tidyr)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)

  input_file <- file.path(local_dir, infile)

  read_df <- fread(input_file, header = TRUE) %>%
    as.data.frame() %>%
    separate(Bed_File, c("celltype", "chromState", "txt"), sep = "[.]")

  read_df <- read_df[complete.cases(read_df), ]

  # Compute Enrichment: observed/expected and take log2
  read_df <- read_df %>%
    mutate(
      inbed = as.numeric(read_df[, 4]),
      expected = as.numeric(read_df[, 5]),
      foldEnrich = log2(inbed / expected)
    )

  read_tf <- dcast(read_df, celltype ~ chromState, value.var = "foldEnrich")

  data <- read_tf[, -1]
  mat_data <- as.matrix(data)
  mat_data <- apply(mat_data, 2, function(x) as.numeric(as.character(x)))
  mat_data[is.na(mat_data)] <- 0

  # Remove numeric prefixes and underscores from column names
  cleaned_colnames <- gsub("^[0-9]+_", "", colnames(data))
  colnames(mat_data) <- cleaned_colnames
  rownames(mat_data) <- read_tf$celltype

  roadmap_df <- fread(annotation_file, sep = "\t")
  anno_df <- left_join(data.frame(EpigenomeID = rownames(mat_data)), roadmap_df, by = "EpigenomeID")

  # Sorting anno_df based on rownames of mat_data to ensure consistent ordering
  anno_df <- anno_df[match(rownames(mat_data), anno_df$EpigenomeID), ]

  # Define color mapping for groups
  group_colors <- unique(anno_df[, c("Group", "Color")])
  color_mapping <- setNames(group_colors$Color, group_colors$Group)

  # Creating the row annotation object
  row_anno <- rowAnnotation(CellType = anno_df$Group, 
                            col = list(CellType = color_mapping),
                            width = unit(1, "cm"),
                            #annotation_name_side = "left",
                            annotation_legend_param = list(title = "CellType", 
                                                           title_gp = gpar(fontsize = 10, fontface = "bold"),
                                                           labels_gp = gpar(fontsize = 8)))

  
  #output_file <- file.path(local_dir, "GNR_fullHeatmap_Log2FoldEnrich_FinalVersionCall.pdf")
  output_file <- file.path(local_dir, paste0(output_fileprefix, "_fullHeatmap_Log2FoldEnrich_FinalVersionCall.pdf"))
  pdf(output_file)

  # Set custom fix color
  my_colors <- colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
  
  ht <- Heatmap(mat_data, name = "foldEnrich", 
        row_names_side = "right",
        col = my_colors,  # Apply your color scale
        #row_names_gp = gpar(fontsize = 3),
        column_names_gp = gpar(fontsize = 9.5),
        show_row_names = FALSE,
        left_annotation = row_anno,
        heatmap_legend_param = list(color_bar = "continuous", 
                                    title = "log2(foldEnrich)",
                                    title_gp = gpar(fontsize = 9, fontface = "bold" ),  # Adjust title font size here
                                    labels_gp = gpar(fontsize = 8),  # Adjust label font size here
                                    at = c(-4, -2, 0, 2, 4), 
                                    labels = c("<=-4", "-2", "0", "2", ">=4")))

  draw(ht)

  # Capture the clustering order of rows and columns
  row_order <- row_order(ht)
  column_order <- column_order(ht)

  dev.off()

  # Return the clustering order of rows and columns
  return(list(row_order = row_order, column_order = column_order))

  #return(ht_list)

}


infile="StatisticSummaryFile_chromHMM_15state_tl_topsnps.txt"
annotation_file="roadmap_annotation.txt"


local_dir="gregor/fulluniq_redo"
output_fileprefix="MAGE"
foldEnrich_heatmap_mage <- plot_foldenrich_fullheatmap(local_dir, infile, annotation_file, output_fileprefix)

# Extract row and column order from the drawing details
row_order_mage <- foldEnrich_heatmap_mage$row_order
col_order_mage  <- foldEnrich_heatmap_mage$column_order

local_dir="gregor/fulluniq_mageqtl"
output_fileprefix="GNR"
foldEnrich_heatmap_gnr  <- plot_foldenrich_fullheatmap(local_dir, infile, annotation_file, output_fileprefix)

# Extract row and column order from the drawing details
row_order_gnr <- foldEnrich_heatmap_gnr$row_order
col_order_gnr  <- foldEnrich_heatmap_gnr$column_order


plot_pvalue_fullheatmap <- function(local_dir, infile, annotation_file, output_fileprefix, column_order = NULL, row_order = NULL) {

  # Ensure necessary libraries are loaded
  library(data.table)
  library(tidyr)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(reshape2)

  input_file <- file.path(local_dir, infile)

  read_df <- fread(input_file, header = T) %>% 
    as.data.frame() %>%
    separate(Bed_File, c("celltype", "chromState", "txt"), sep = "[.]")

  read_df <- read_df[complete.cases(read_df), ]
  read_tf <- read_df %>% dcast(celltype ~ chromState, value.var = "PValue")

  data <- read_tf[, 2:ncol(read_tf)]
  mat_data <- apply(as.matrix.noquote(data), 2, as.numeric)
  mat_data[mat_data == 0] <- 1e-323 
  mat_data <- -log10(mat_data)
  mat_data[is.na(mat_data)] <- 0

  # Clean up column names
  cleaned_colnames <- gsub("^[0-9]+_", "", colnames(data))
  colnames(mat_data) <- cleaned_colnames
  rownames(mat_data) <- read_tf$celltype

  roadmap_df <- fread(annotation_file, sep = "\t")
  anno_df <- left_join(data.frame(EpigenomeID = rownames(mat_data)), roadmap_df, by = "EpigenomeID")

  # Sorting anno_df based on rownames of mat_data to ensure consistent ordering
  anno_df <- anno_df[match(rownames(mat_data), anno_df$EpigenomeID), ]

  # Define color mapping for groups
  group_colors <- unique(anno_df[, c("Group", "Color")])
  color_mapping <- setNames(group_colors$Color, group_colors$Group)

  row_anno <- rowAnnotation(CellType = anno_df$Group, 
                            col = list(CellType = color_mapping),
                            width = unit(1, "cm"),
                            annotation_legend_param = list(title = "CellType", 
                                                           title_gp = gpar(fontsize = 10, fontface = "bold"),
                                                           labels_gp = gpar(fontsize = 8)))


  # Order the matrix data if row_order and/or column_order are provided
  if (!is.null(row_order)) {
    mat_data <- mat_data[row_order, , drop = FALSE]
    anno_df <- anno_df[row_order, , drop = FALSE] # Make sure to reorder anno_df as well
    row_anno <- rowAnnotation(CellType = anno_df$Group, # Update row_anno based on the new order
                            col = list(CellType = color_mapping),
                            width = unit(1, "cm"),
                            annotation_legend_param = list(title = "CellType", 
                                                           title_gp = gpar(fontsize = 10, fontface = "bold"),
                                                           labels_gp = gpar(fontsize = 8)))
  }
  
  if (!is.null(column_order)) {
    mat_data <- mat_data[, column_order, drop = FALSE]
  }

  output_file <- file.path(local_dir, paste0(output_fileprefix, "_full_PvalueHeatmap_ClusterAligned1.pdf"))
  pdf(output_file)

  # CUSTOM COLORING OF PLOTS
  # Define the red monochromatic color palette
  redmono <- c("#EF3B2C", "#FB6A4A", "#FC9272", "#FCBBA1", "#FEE0D2", "#FFF5F0")
  redmono <- c("#EF3B2C", "#FB6A4A", "#FC9272", "#FCBBA1")
  bluemono = c("#9ECAE1", "#C6DBEF", "#DEEBF7", "#F7FBFF")

  # Create a color palette of length 100
  redmono_palette <- colorRampPalette(rev(redmono))(100)

  # Create breaks and corresponding colors
  breaks_blue <- seq(0, 5, length.out = length(bluemono) + 1)
  colors_blue <- bluemono

  breaks_red <- seq(6, 200, length.out = 100)  # Adjusted to have one less value
  colors_red <- redmono_palette

  # Combine breaks and colors
  breaks <- c(breaks_blue, breaks_red[-1])  # removing the overlapping 5
  colors <- c(colors_blue, colors_red)

  # Ensure that the lengths of breaks and colors are the same
  if(length(breaks) != length(colors)) {
    stop("Length of breaks should be equal to colors")
  }

  # Create the color mapping function
  my_colors <- colorRamp2(breaks, colors)

  #my_color <- my_colors(mat_data)
  ht <- Heatmap(mat_data, name = "foldEnrich", 
                col = my_colors,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = FALSE,
                show_column_names = TRUE,
                column_names_gp = gpar(fontsize = 9.5),
                left_annotation = row_anno,
                heatmap_legend_param = list(color_bar = "continuous", 
                                            title = "-log10(Pvalue)",
                                            title_gp = gpar(fontsize = 9, fontface = "bold"),
                                            labels_gp = gpar(fontsize = 8),
                                            at = c(0, 25, 50, 100, 150, 200), 
                                            labels = c("0", "25", "50", "100", "150", ">=200")))

  draw(ht)

  dev.off()
}

###################################################
#                                                 #    
# Column ordering based on original column orders #
#                                                 #
###################################################

local_dir="gregor/fulluniq_mageqtl"
output_fileprefix = "GNR_"

row_order_gnr <- foldEnrich_heatmap_gnr$row_order
column_order_gnr <- foldEnrich_heatmap_gnr$column_order
#plot_pvalue_fullheatmap(local_dir, infile, annotation_file, output_fileprefix, column_order_gnr)
plot_pvalue_fullheatmap(local_dir, infile, annotation_file, output_fileprefix, column_order_gnr, row_order_gnr)

# Extract row and column order
local_dir="gregor/fulluniq_redo"
output_fileprefix = "MAGE_"

row_order_mage <- foldEnrich_heatmap_mage$row_order
column_order_mage <- foldEnrich_heatmap_mage$column_order
#plot_pvalue_fullheatmap(local_dir, infile, annotation_file, output_fileprefix, column_order_mage)
plot_pvalue_fullheatmap(local_dir, infile, annotation_file, output_fileprefix, column_order_mage, row_order_mage)



###################################################
#                                                 #
# Code for DecileHeatmap aligned to Pvalue        #
#                                                 #
###################################################

# Aligns clustering with fold-enrich heatmaps

generate_decile_heatmaps <- function(elements, local_dir, annotation_file) {
  
  # Load libraries
  library(data.table)
  library(tidyr)
  library(dplyr)
  library(reshape2)
  library(ComplexHeatmap)
  library(circlize)

  # Custom color scaling:
  my_colors <- colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

  # Step 3: Generate Heatmap
  for (cell_id in elements) {
    print(paste("Processing ...", cell_id))
    infile <- paste0(cell_id, "-StatisticSummaryFile_chromHMM_15state_tl_topsnps.txt")
    
    input_file <- file.path(local_dir, infile)
    
    read_df <- fread(input_file, header = FALSE)
    names(read_df) <- c("Bed_File", "InBed_Index_SNP", "ExpectNum_of_InBed_SNP",  "ID")
    
    read_df <- read_df %>%
      as.data.frame() %>%
      separate(Bed_File, c("celltype", "chromState", "txt"), sep = "[.]") %>%
      separate(ID, c("PValue", "DecileID"), sep = "[ ]") %>%
      filter(complete.cases(.)) %>%
      mutate(
        inbed = as.numeric(InBed_Index_SNP),
        expected = as.numeric(ExpectNum_of_InBed_SNP),
        foldEnrich = log2(inbed / expected)
      )
    
    read_tf <- dcast(read_df, DecileID ~ chromState, value.var = "foldEnrich")
    data <- read_tf[, -1]
    mat_data <- apply(as.matrix.noquote(data), 2, as.numeric)
    mat_data[is.na(mat_data)] <- 0
    
    # Clean up column names
    cleaned_colnames <- gsub("^[0-9]+_", "", colnames(data))
    colnames(mat_data) <- cleaned_colnames
    rownames(mat_data) <- read_tf$DecileID
    
    order_vec_rows <- rev(c("decile1", "decile2", "decile3", "decile4", "decile5", 
                            "decile6", "decile7", "decile8", "decile9", "decile10"))
    
    # Order the matrix rows
    mat_ordered <- mat_data[order_vec_rows, ]
    
    # Inside your for loop and after generating mat_ordered
    # Define color scale
    # my_colors <- colorRamp2(c(min(mat_ordered), 0, max(mat_ordered)), c("blue", "white", "red"))

    output_file <- file.path(local_dir, paste0(cell_id, "_MAGE_Log2FoldEnrichFinalversionForCLUSTER.pdf"))
    pdf(output_file, width = 10, height = 8)

    ht <- Heatmap(mat_ordered, 
            name = "foldEnrich", 
            col = my_colors,
            show_row_names = TRUE, 
            show_column_names = TRUE, 
            cluster_rows = FALSE, 
            cluster_columns = TRUE, 
            row_names_side = "left", 
            column_title = "Chromatin States", 
            row_title = "Deciles",
            #heatmap_legend_param = list(title = "log2(foldEnrich)", at = c(-5, -3, -2, 0, 2, 3, 5), labels = c("-5", "-3", "-2", "0", "2", "3", "5"))
            #heatmap_legend_param = list(title = "log2(foldEnrich)", at = c(-4, -2, 0, 2, 4), labels = c("≤-4", "-2", "0", "2", "≥4"))
            heatmap_legend_param = list(title = "log2(foldEnrich)", at = c(-4, -2, 0, 2, 4), labels = c("<=-4", "-2", "0", "2", ">=4"))
    )

    draw(ht)

    # Capture the clustering order of rows and columns
    row_order <- row_order(ht)
    column_order <- column_order(ht)

    dev.off()
    
    # Return the clustering order of rows and columns
    return(list(row_order = row_order, column_order = column_order))

  }

}

# Define your elements vector
# elements <- c("E116", "E034", "E032", "E051", "E046")

# Define the local directory
local_dir <- "gregor/decileuniq_redo"
# Define the annotation file path (change this to your actual path)
annotation_file <- "roadmap_annotation.txt"

# Call the function
elements <- c("E116")
foldEnrich_heatmap_E116 <- generate_decile_heatmaps(elements, local_dir, annotation_file)
# Extract row and column order from the drawing details
row_order_E116 <- foldEnrich_heatmap_E116$row_order
col_order_E116  <- foldEnrich_heatmap_E116$column_order


elements <- c("E034")
foldEnrich_heatmap_E034 <- generate_decile_heatmaps(elements, local_dir, annotation_file)
# Extract row and column order from the drawing details
row_order_E034 <- foldEnrich_heatmap_E034$row_order
col_order_E034  <- foldEnrich_heatmap_E034$column_order

elements <- c("E032")
foldEnrich_heatmap_E032 <- generate_decile_heatmaps(elements, local_dir, annotation_file)
# Extract row and column order from the drawing details
row_order_E032 <- foldEnrich_heatmap_E032$row_order
col_order_E032  <- foldEnrich_heatmap_E032$column_order

elements <- c("E051")
foldEnrich_heatmap_E051 <- generate_decile_heatmaps(elements, local_dir, annotation_file)
# Extract row and column order from the drawing details
row_order_E051 <- foldEnrich_heatmap_E051$row_order
col_order_E051  <- foldEnrich_heatmap_E051$column_order


elements <- c("E046")
foldEnrich_heatmap_E046 <- generate_decile_heatmaps(elements, local_dir, annotation_file)
# Extract row and column order from the drawing details
row_order_E046 <- foldEnrich_heatmap_E046$row_order
col_order_E046  <- foldEnrich_heatmap_E046$column_order


generate_pvalue_decile_heatmap <- function(elements, local_dir, infile, annotation_file, output_fileprefix, column_order = NULL, row_order = NULL){

  # Load libraries
  library(data.table)
  library(tidyr)
  library(dplyr)
  library(pheatmap)
  library(reshape2)

  #suppressWarnings()

  for (cell_id in elements){
    print(paste("Processing ...", cell_id))
    infile <- paste0(cell_id, "-StatisticSummaryFile_chromHMM_15state_tl_topsnps.txt")
    annotation_file="roadmap_annotation.txt"

    input_file <- file.path(local_dir, infile)

    read_df <- fread(input_file, header = F)
    names(read_df) <- c("Bed_File", "InBed_Index_SNP", "ExpectNum_of_InBed_SNP",  "ID")
    
    read_df <-  read_df %>% as.data.frame() %>%
      separate(Bed_File, c("celltype", "chromState", "txt"), sep = "[.]") %>%
      separate(ID, c("PValue", "DecileID"), sep = "[ ]")

    read_df$PValue <- as.numeric(read_df$PValue)
    read_df <- read_df[complete.cases(read_df), ]
    read_tf <- read_df %>% dcast(DecileID ~ chromState, value.var = "PValue")

    data <- read_tf[, 2:ncol(read_tf)]
    mat_data <- apply(as.matrix.noquote(data), 2, as.numeric)

    mat_data <- -log10(mat_data)
    mat_data[is.na(mat_data)] <- 0

    # Clean up column names
    cleaned_colnames <- gsub("^[0-9]+_", "", colnames(data))
    colnames(mat_data) <- cleaned_colnames
    rownames(mat_data) <- read_tf$DecileID

    order_vec_rows <- rev(c("decile1", "decile2", "decile3", "decile4", "decile5", "decile6",
                        "decile7", "decile8", "decile9", "decile10"))

    # Order the matrix rows and columns
    mat_ordered <- mat_data[order_vec_rows, ]

    # SANITY CHECK:
    #row_order <- row_order_E116
    #col_order <- col_order_E116
    
    # Order the matrix data if row_order and/or column_order are provided
    if (!is.null(row_order)) {
      mat_ordered <- mat_ordered[row_order, , drop = FALSE]
    }
    if (!is.null(column_order)) {
      mat_ordered <- mat_ordered[, column_order, drop = FALSE]
    }


    output_file <- file.path(local_dir, paste0(cell_id, "_MAGE_PvalueHeatmap_clusterAligned1.pdf"))
    pdf(output_file, width = 10, height = 8)

    # CUSTOM COLORING OF PLOTS
    # Define the red monochromatic color palette
    redmono <- c("#EF3B2C", "#FB6A4A", "#FC9272", "#FCBBA1")
    bluemono = c("#9ECAE1", "#C6DBEF", "#DEEBF7", "#F7FBFF")

    # Create a color palette of length 100
    redmono_palette <- colorRampPalette(rev(redmono))(100)

    # Create breaks and corresponding colors
    breaks_blue <- seq(0, 5, length.out = length(bluemono) + 1)
    colors_blue <- bluemono

    breaks_red <- seq(6, 200, length.out = 100)  # Adjusted to have one less value
    colors_red <- redmono_palette

    # Combine breaks and colors
    breaks <- c(breaks_blue, breaks_red[-1])  # removing the overlapping 5
    colors <- c(colors_blue, colors_red)

    # Ensure that the lengths of breaks and colors are the same
    if(length(breaks) != length(colors)) {
      stop("Length of breaks should be equal to colors")
    }

    # Create the color mapping function
    my_colors <- colorRamp2(breaks, colors)

    #my_color <- my_colors(mat_ordered)
    ht <- Heatmap(mat_ordered, name = "PvalEnrich", 
                  col = my_colors,
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  show_row_names = TRUE,
                  show_column_names = TRUE,
                  row_names_side = "left", 
                  #column_title = "Chromatin States", 
                  row_title = "Deciles",
                  column_names_gp = gpar(fontsize = 11),
                  row_names_gp = gpar(fontsize = 11),
                  #left_annotation = row_anno,
                  heatmap_legend_param = list(color_bar = "continuous", 
                                              title = "-log10(Pvalue)",
                                              title_gp = gpar(fontsize = 9, fontface = "bold"),
                                              labels_gp = gpar(fontsize = 8),
                                              at = c(0, 25, 50, 100, 150, 200), 
                                              labels = c("0", "25", "50", "100", "150", ">=200")))

    draw(ht)
    dev.off()
  }
}


###################################################
#                                                 #    
# Column ordering based on original column orders #
#                                                 #
###################################################


# Define the local directory
local_dir <- "gregor/decileuniq_redo"
elements <- c("E116")
# row_order_E116 <- foldEnrich_heatmap_E116$row_order
# col_order_E116  <- foldEnrich_heatmap_E116$column_order
generate_pvalue_decile_heatmap(elements, local_dir, infile, annotation_file, output_fileprefix, col_order_E116, row_order_E116)

# Define the local directory
local_dir <- "gregor/decileuniq_redo"
elements <- c("E034")
# row_order_gnr <- foldEnrich_heatmap_gnr$row_order
# col_order_gnr <- foldEnrich_heatmap_gnr$col_order
generate_pvalue_decile_heatmap(elements, local_dir, infile, annotation_file, output_fileprefix, col_order_E034, row_order_E034)


# Define the local directory
local_dir <- "gregor/decileuniq_redo"
elements <- c("E032")
# row_order_mage <- foldEnrich_heatmap_mage$row_order
# col_order_mage <- foldEnrich_heatmap_mage$col_order
generate_pvalue_decile_heatmap(elements, local_dir, infile, annotation_file, output_fileprefix, col_order_E032, row_order_E032)


# Define the local directory
local_dir <- "gregor/decileuniq_redo"
elements <- c("E051")
# row_order_gnr <- foldEnrich_heatmap_gnr$row_order
# col_order_gnr <- foldEnrich_heatmap_gnr$col_order
generate_pvalue_decile_heatmap(elements, local_dir, infile, annotation_file, output_fileprefix, col_order_E051, row_order_E051)


# Define the local directory
local_dir <- "gregor/decileuniq_redo"
elements <- c("E046")
# row_order_gnr <- foldEnrich_heatmap_gnr$row_order
# col_order_gnr <- foldEnrich_heatmap_gnr$col_order
generate_pvalue_decile_heatmap(elements, local_dir, infile, annotation_file, output_fileprefix, col_order_E046, row_order_E046)



decile_pheatmap_based <- function(elements, local_dir, annotation_file) {
  
  # Load libraries
  library(data.table)
  library(tidyr)
  library(dplyr)
  library(pheatmap)
  library(reshape2)


  for (cell_id in elements) {
    print(paste("Processing ...", cell_id))
    infile <- paste0(cell_id, "-StatisticSummaryFile_chromHMM_15state_tl_topsnps.txt")
    
    input_file <- file.path(local_dir, infile)
    
    read_df <- fread(input_file, header = FALSE)
    names(read_df) <- c("Bed_File", "InBed_Index_SNP", "ExpectNum_of_InBed_SNP",  "ID")
    
    read_df <- read_df %>%
      as.data.frame() %>%
      separate(Bed_File, c("celltype", "chromState", "txt"), sep = "[.]") %>%
      separate(ID, c("PValue", "DecileID"), sep = "[ ]") %>%
      filter(complete.cases(.)) %>%
      mutate(
        inbed = as.numeric(InBed_Index_SNP),
        expected = as.numeric(ExpectNum_of_InBed_SNP),
        foldEnrich = log2(inbed / expected)
      )
    
    read_tf <- dcast(read_df, DecileID ~ chromState, value.var = "foldEnrich")
    data <- read_tf[, -1]
    mat_data <- apply(as.matrix.noquote(data), 2, as.numeric)
    mat_data[is.na(mat_data)] <- 0
    
    # Clean up column names
    cleaned_colnames <- gsub("^[0-9]+_", "", colnames(data))
    colnames(mat_data) <- cleaned_colnames
    rownames(mat_data) <- read_tf$DecileID
    
    order_vec_rows <- rev(c("decile1", "decile2", "decile3", "decile4", "decile5", 
                            "decile6", "decile7", "decile8", "decile9", "decile10"))
    
    # Order the matrix rows
    mat_ordered <- mat_data[order_vec_rows, ]
    
    output_file <- file.path(local_dir, paste0(cell_id, "_MAGE_chromHMM_15state_tl_topsnpsClusterOff_Log2FoldEnrich.pdf"))
    
    pdf(output_file)
    
    pheatmap(mat_ordered,
      cluster_rows = FALSE, 
      cluster_cols = TRUE,
      fontsize_row = 8,
      fontsize_col = 9,
      annotation_legend = FALSE
    )
    
    dev.off()
  }
}

# Define your elements vector
elements <- c("E116", "E034", "E032", "E051", "E046")

# Define the local directory
local_dir <- "gregor/decileuniq_redo"

# Define the annotation file path (change this to your actual path)
annotation_file <- "roadmap_annotation.txt"

# Call the function
decile_pheatmap_based(elements, local_dir, annotation_file)



#########################################
#                                       #
# Pheatmap based decile plots           #
#                                       #
#########################################


# Ensure necessary libraries are loaded
library(data.table)
library(tidyr)
library(dplyr)
library(pheatmap)
library(reshape2)

elements <- c("E116", "E034", "E032", "E051", "E046")

# Local dir
local_dir <- "gregor/decileuniq"

for (cell_id in elements){
    #cell_id <- "E046"
    infile <- paste0(cell_id, "-StatisticSummaryFile_chromHMM_15state_tl_topsnps.txt")
    annotation_file="roadmap_annotation.txt"

    input_file <- file.path(local_dir, infile)

    read_df <- fread(input_file, header = F)
    names(read_df) <- c("Bed_File", "InBed_Index_SNP", "ExpectNum_of_InBed_SNP",  "ID")
    
    read_df <-  read_df %>% as.data.frame() %>%
      separate(Bed_File, c("celltype", "chromState", "txt"), sep = "[.]") %>%
      separate(ID, c("PValue", "DecileID"), sep = "[ ]")

    read_df$PValue <- as.numeric(read_df$PValue)
    read_df <- read_df[complete.cases(read_df), ]
    read_tf <- read_df %>% dcast(DecileID ~ chromState, value.var = "PValue")

    data <- read_tf[, 2:ncol(read_tf)]
    mat_data <- apply(as.matrix.noquote(data), 2, as.numeric)

    mat_data <- -log10(mat_data)
    mat_data[is.na(mat_data)] <- 0

    colnames(mat_data) <- colnames(data)
    rownames(mat_data) <- read_tf$DecileID

    # Desired column and row order
    order_vec_cols <- c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh",
                 "8_ZNF-Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC",
                 "14_ReprPCWk", "15_Quies")

    order_vec_rows <- rev(c("decile1", "decile2", "decile3", "decile4", "decile5", "decile6",
                        "decile7", "decile8", "decile9", "decile10"))

    # Order the matrix rows and columns
    mat_ordered <- mat_data[order_vec_rows, order_vec_cols]


    output_file <- file.path(local_dir, paste0(cell_id, "chromHMM_15state_tl_topsnpsClusterOff.pdf"))

    pdf(output_file)

    pheatmap(mat_ordered,
      cluster_rows = FALSE, cluster_cols = FALSE,
      fontsize_row = 8,
      fontsize_col = 9,
      #annotation_row = annotation_row,
      #annotation_colors = mat_colors,
      annotation_legend = FALSE,
      legend_labels = "-log10(PValue)",
      main = paste(cell_id, "Decile-eQTL-Enrichment")
    )

    dev.off()

}
