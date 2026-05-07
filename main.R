library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
  counts <- read.delim(counts_csv, header = TRUE, stringsAsFactors = FALSE)
  meta <- read.csv(metafile_csv, stringsAsFactors = FALSE)
  
  # subset metadata
  meta <- meta[meta$timepoint %in% selected_times, ]
  
  # HARD FIX for DESeq2 directionality
  meta$timepoint <- factor(meta$timepoint)
  meta$timepoint <- relevel(meta$timepoint, ref = "vP0")
  
  # match sample columns
  sample_cols <- intersect(colnames(counts), meta$samplename)
  
  # subset counts
  counts_sub <- counts[, sample_cols, drop = FALSE]
  
  # reorder metadata to match counts column order
  meta <- meta[match(sample_cols, meta$samplename), ]
  
  # safety check
  stopifnot(all(meta$samplename == colnames(counts_sub)))
  
  # build SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = as.matrix(counts_sub)),
    colData = meta
  )
  
  return(se)
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  # build dataset
  dds <- DESeq2::DESeqDataSet(se, design = design)
  
  # 🔥 CRITICAL FIX: explicitly set reference level HERE
  dds$timepoint <- factor(dds$timepoint)
  dds$timepoint <- relevel(dds$timepoint, ref = "vP0")
  
  # run DESeq2
  dds <- DESeq2::DESeq(dds)
  
  # extract results
  res <- DESeq2::results(dds)
  
  return(list(
    dds = dds,
    results = as.data.frame(res)
  ))
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`. Ensure
#' that the column name for your rownames is called "genes". 
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
  # convert rownames to column
  deseq2_res$genes <- rownames(deseq2_res)
  
  # convert to tibble
  df <- tibble::as_tibble(deseq2_res)
  
  # create status column
  df$volc_plot_status <- dplyr::case_when(
    df$padj < padj_threshold & df$log2FoldChange > 0 ~ "UP",
    df$padj < padj_threshold & df$log2FoldChange < 0 ~ "DOWN",
    TRUE ~ "NS"
  )
  
  return(df)
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  p <- ggplot2::ggplot(labeled_results, ggplot2::aes(x = pvalue)) +
    ggplot2::geom_histogram(bins = 50) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Histogram of DESeq2 p-values",
      x = "p-value",
      y = "Frequency"
    )
  
  return(p)
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  # filter significant genes
  sig_data <- labeled_results[labeled_results$padj < padj_threshold, ]
  
  p <- ggplot2::ggplot(sig_data, ggplot2::aes(x = log2FoldChange)) +
    ggplot2::geom_histogram(bins = 50) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Log2 Fold Change (Significant Genes)",
      x = "log2 Fold Change",
      y = "Frequency"
    )
  
  return(p)
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
  # get top genes by smallest padj
  top_genes <- labeled_results[order(labeled_results$padj), "genes"]
  top_genes <- head(top_genes$genes, num_genes)
  
  # extract normalized counts
  norm_counts <- DESeq2::counts(dds_obj, normalized = TRUE)
  
  # subset to top genes
  norm_counts <- norm_counts[top_genes, , drop = FALSE]
  
  # reshape for plotting
  plot_df <- reshape2::melt(norm_counts)
  colnames(plot_df) <- c("gene", "sample", "count")
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = sample, y = count)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ gene, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(
      title = "Normalized Counts for Top Differentially Expressed Genes",
      x = "Sample",
      y = "Normalized Count"
    )
  
  return(p)
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  # avoid log issues with NA or 0 padj
  labeled_results$padj_adj <- ifelse(
    is.na(labeled_results$padj) | labeled_results$padj == 0,
    NA,
    labeled_results$padj
  )
  
  labeled_results$neg_log10_padj <- -log10(labeled_results$padj_adj)
  
  p <- ggplot2::ggplot(
    labeled_results,
    ggplot2::aes(
      x = log2FoldChange,
      y = neg_log10_padj,
      color = volc_plot_status
    )
  ) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Volcano Plot",
      x = "log2 Fold Change",
      y = "-log10 adjusted p-value",
      color = "Status"
    )
  
  return(p)
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  id_map <- read.delim(id2gene_path, stringsAsFactors = FALSE)
  
  colnames(id_map)[1] <- "genes"
  colnames(id_map)[2] <- "symbol"
  
  merged <- merge(labeled_results, id_map, by = "genes", all.x = TRUE)
  
  # fallback for missing symbols
  merged$symbol[is.na(merged$symbol) | merged$symbol == ""] <- merged$genes
  
  # 🔥 FIX: remove NA/Inf log2FC (fgsea requirement)
  merged <- merged[is.finite(merged$log2FoldChange), ]
  
  merged <- merged[order(merged$log2FoldChange, decreasing = TRUE), ]
  
  vec <- merged$log2FoldChange
  names(vec) <- merged$symbol
  
  return(vec)
}

#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  # read gene sets from GMT file
  pathways <- fgsea::gmtPathways(gmt_file_path)
  
  # run fgsea
  fgsea_res <- fgsea::fgsea(
    pathways = pathways,
    stats = rnk_list,
    minSize = min_size,
    maxSize = max_size
  )
  
  # convert to tibble and sort by p-value
  fgsea_df <- as.data.frame(fgsea_res)
  fgsea_df <- fgsea_df[order(fgsea_df$pval), ]
  
  return(tibble::as_tibble(fgsea_df))
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
  # top positive NES
  top_pos <- fgsea_results[order(fgsea_results$NES, decreasing = TRUE), ]
  top_pos <- head(top_pos, num_paths)
  
  # top negative NES
  top_neg <- fgsea_results[order(fgsea_results$NES, decreasing = FALSE), ]
  top_neg <- head(top_neg, num_paths)
  
  # combine
  plot_df <- rbind(top_pos, top_neg)
  
  # reorder pathways for plotting
  plot_df$pathway <- factor(plot_df$pathway, levels = plot_df$pathway[order(plot_df$NES)])
  
  # plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = NES, y = pathway, fill = NES)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Top Enriched Pathways (FGSEA)",
      x = "Normalized Enrichment Score (NES)",
      y = "Pathway"
    )
  
  return(p)
}

# debugging
head(counts)
colnames(counts)