if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("phyloseq", quietly = TRUE)) {
    BiocManager::install("phyloseq")
}

if (!requireNamespace("RColorBrewer", quietly = TRUE) || 
    !requireNamespace("patchwork", quietly = TRUE) || 
    !requireNamespace("dplyr", quietly = TRUE) || 
    !requireNamespace("tidyr", quietly = TRUE)) {
    install.packages(c("RColorBrewer", "patchwork", "dplyr", "tidyr"), repos = "https://cloud.r-project.org")
}

setwd(".")

library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("patchwork")
library("dplyr")
library("tidyr")

merged_metagenomes <- import_biom("data/processing/5_biom/cuatroc.biom")

merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

merged_metagenomes <- subset_taxa(merged_metagenomes, Kingdom == "Bacteria")
merged_metagenomes <- subset_taxa(merged_metagenomes, Genus != "")
percentages <- transform_sample_counts(merged_metagenomes, function(x) x*100 / sum(x) )
percentages_glom <- tax_glom(percentages, taxrank = 'Phylum')
percentages_df <- psmelt(percentages_glom)
percentages_df$Phylum <- as.factor(percentages_df$Phylum)
phylum_colors_rel<- colorRampPalette(brewer.pal(8,"Dark2")) (length(levels(percentages_df$Phylum)))
relative_plot <- ggplot(data=percentages_df, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  scale_fill_manual(values = phylum_colors_rel)

wide_percentages_df <- percentages_df %>%
  select(OTU, Sample, Abundance, Kingdom, Phylum) %>%
  tidyr::pivot_wider(
    names_from = Sample,
    values_from = Abundance,
    values_fill = 0 
  )

write.csv(wide_percentages_df, "data/results/phylum_relative_abundances.csv", row.names = FALSE)
ggsave("data/results/phylum_relative_abundances.pdf", relative_plot, width = 10, height = 6)
