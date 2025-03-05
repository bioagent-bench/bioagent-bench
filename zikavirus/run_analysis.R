# Load required libraries
library(DESeq2)
library(dplyr)
library(readr)
library(gprofiler2)

# Load expression matrix
expr_df <- read.csv("processing/3_normalized/repCpmMatrix_featureCounts.csv", row.names = 1)
expr_df <- expr_df[rowSums(expr_df) > 0, ]

# Filter out lowly expressed genes
mask_low_vals <- rowSums(expr_df > 0.3) > 2
expr_df <- expr_df[mask_low_vals, ]

# Load metadata
meta_df <- read.csv("data/metadata.csv", row.names = "Run")
meta_df <- meta_df[colnames(expr_df), ]
layout_platform <- c(PAIRED = "MiSeq", SINGLE = "NextSeq 500")
meta_df$platform <- layout_platform[meta_df$LibraryLayout]

meta_df$condition <- ifelse(meta_df$infection_status == "Zika infected", "infected", "control")
meta_df$condition <- factor(meta_df$condition, levels = c("control", "infected"))

platform_results <- list()
cd_results <- data.frame(row.names = rownames(expr_df))

# Perform differential expression analysis per platform
for (platform in unique(meta_df$platform)) {
  platform_meta <- meta_df[meta_df$platform == platform, ]
  platform_counts <- expr_df[, rownames(platform_meta)]
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(platform_counts),
    colData = platform_meta,
    design = ~ condition
  )
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("condition", "infected", "control"))
  res_df <- as.data.frame(res)
  res_df$abs_log2FC <- abs(res_df$log2FoldChange)
  platform_results[[platform]] <- res_df
  cd_results[[platform]] <- res_df$log2FoldChange
  
  sorted_res <- res_df[order(res_df$abs_log2FC, decreasing = TRUE), ]
  top_degs <- head(sorted_res, 600)

  # Separate upregulated and downregulated genes
  up_genes <- rownames(top_degs[top_degs$log2FoldChange > 0, ])
  dn_genes <- rownames(top_degs[top_degs$log2FoldChange < 0, ])
  
  platform_results[[paste0(platform, "-up")]] <- up_genes
  platform_results[[paste0(platform, "-dn")]] <- dn_genes

  write.csv(res_df, paste0("results/", platform, "_DESeq2_results.csv"))
  write.csv(up_genes, paste0("results/", platform, "_upregulated_genes.csv"), row.names = FALSE)
  write.csv(dn_genes, paste0("results/", platform, "_downregulated_genes.csv"), row.names = FALSE)
}

write.csv(cd_results, "results/combined_log2FC_results.csv")

# Function to run g:Profiler enrichment
run_gprofiler <- function(gene_list, platform_name, regulation_type) {
  enrichment <- gost(
    query = gene_list, 
    organism = "hsapiens",
    sources = c("KEGG", "GO:BP", "ChEA", "TF")
  )
  
  if (!is.null(enrichment$result)) {
    enrichment_clean <- enrichment$result
    enrichment_clean[] <- lapply(enrichment_clean, function(x) if (is.list(x)) sapply(x, paste, collapse = "; ") else x)
    
    output_file <- paste0("results/", platform_name, "_", regulation_type, "_enrichment.csv")
    write.csv(enrichment_clean, output_file, row.names = FALSE)
    return(output_file)
  } else {
    return(NULL)
  }
}

# Run enrichment for each platform separately
enrichment_results <- data.frame(Platform = character(), Size = integer(), File = character(), stringsAsFactors = FALSE)

for (platform in unique(meta_df$platform)) {
  up_genes <- platform_results[[paste0(platform, "-up")]]
  dn_genes <- platform_results[[paste0(platform, "-dn")]]
  
  up_file <- run_gprofiler(up_genes, platform, "up")
  dn_file <- run_gprofiler(dn_genes, platform, "dn")

  enrichment_results <- rbind(enrichment_results, data.frame(Platform = paste0(platform, "-up"), Size = length(up_genes), File = up_file))
  enrichment_results <- rbind(enrichment_results, data.frame(Platform = paste0(platform, "-dn"), Size = length(dn_genes), File = dn_file))
}

# Save enrichment metadata
write.csv(enrichment_results, "results/enrichment_summary.csv", row.names = FALSE)

# Load previously saved Zika-induced (upregulated) genes
miseq_up_genes <- read.csv("results/MiSeq_upregulated_genes.csv", stringsAsFactors = FALSE)[,1]
nextseq_up_genes <- read.csv("results/NextSeq 500_upregulated_genes.csv", stringsAsFactors = FALSE)[,1]

# Combine upregulated genes from both platforms
zika_up_genes <- unique(c(miseq_up_genes, nextseq_up_genes))

# Perform enrichment analysis using g:Profiler on MGI database
brain_morphology_enrichment <- gost(
  query = zika_up_genes, 
  organism = "hsapiens",
)

# Check if there are brain morphology-related terms in the results
if (!is.null(brain_morphology_enrichment$result)) {
  enriched_terms <- brain_morphology_enrichment$result
  
  # Extract terms related to brain morphology (filter by keywords)
  brain_related_terms <- enriched_terms %>%
    filter(grepl("brain|neural|cortex|microcephaly", term_name, ignore.case = TRUE))
  
  # Extract overlapping genes from enrichment results
  overlapping_genes <- unique(unlist(strsplit(as.character(brain_related_terms$intersection), ", ")))

  # Save results
  write.csv(brain_related_terms, "results/brain_morphology_enrichment.csv", row.names = FALSE)
  write.csv(overlapping_genes, "results/overlapping_genes_brain_morphology.csv", row.names = FALSE)

  # Print summary
  cat("Overlapping genes related to brain morphology identified and saved in 'results/overlapping_genes_brain_morphology.csv'.\n")
} else {
  cat("No significant enrichment found for brain morphology-related terms.\n")
}
